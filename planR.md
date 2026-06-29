# planR.md — Review-remediation: close the gaps found in the Im(A) / `MonomialPhaseLog` change

**Goal.** Correct the issues surfaced by the max-effort code review of the
`MonomialPhaseLog` / `LiftedMonoPhase` work (the Im(A) phase threaded through
`ProcessSector` / `ProcessSectorLifted` / `IBPProcessSector` / the C++ codegen,
`planAXpDIVv2.md`). The core lifted/IBP math, indexing, and codegen were verified
**correct**; what remains are (R1) a *missing clean refusal* on an unsupported
path that today returns a wrong number, (R2) a detection footgun, (R3) a test
that does not exercise its headline quantity and can pass vacuously, plus two
low-severity hardening items (R4, R5) and optional cleanup (R6).

**Discipline (inherited from `plan.md` §1.4 / §8.5 and `planAXpDIVv2.md`).**
Same rules apply: a known limitation must produce a **clean `$Failed`**, never a
wrong number; every numeric PASS needs ≥2 independent oracles; the real-exponent /
`Direct` path stays **byte-identical** (#25). No fix here may change real-A / real-B
output.

---

## Severity / ordering

| ID | Issue | Severity | Kind |
|----|-------|----------|------|
| **R1** | Unlifted complex-exponent **+ divergent** path drops the phase silently (no `splitliftdiv` refusal) | **High** — silent wrong number | correctness (refusal gap) |
| **R2** | `hasImagA` mis-fires on a symbolic non-`eps` monomial exponent | Medium | robustness |
| **R3** | New test never exercises `Const_A`/`logZ0`; can pass vacuously | Medium | test quality |
| **R4** | Im(A) pole-shift guard nested inside `If[lmf =!= None]` | Low (latent) | hardening |
| **R5** | Codegen emits a dead `exp(0)` phase term when the phase vanishes | Low | hygiene |
| **R6** | Duplication of the `imAList`/`hasImagA` triplet + realify idiom (4 sites) | Low | cleanup |

Do **R1 first** (it is the only one that can return a wrong physics number), then
R3 (so the suite actually guards the fix), then R2/R4/R5, then R6 if time permits.

---

## R1 — Mirror the `splitliftdiv` refusal into the unlifted branch (**High**)

### Diagnosis (verified against source)
- `cxSplit` has **no `isLifted` gate**: `tropical_eval.wl:3880`
  `cxSplit = (cxMode === "SplitRealImag") && (hasImagB || hasImagA)`. So `cxSplit`
  is reachable in the unlifted driver. (Pre-existing for `hasImagB`; this change
  widened it to `hasImagA`.)
- The clean refusal for `cxSplit && divergent` exists **only in the lifted
  branch**: `tropical_eval.wl:4072-4074`
  `If[cxSplit && Length[divergentSectors] > 0, Message[TropicalEval::splitliftdiv]; Return[$Failed]]`.
- The unlifted "Standard mode" branch (`tropical_eval.wl:4082-4111`) has **no such
  guard**. A divergent sector returns from `ProcessSector` via the divergent
  branch `tropical_eval.wl:451-478`, which **omits `MonomialPhaseLog`**; and
  `tropical_eval.wl:4107` reattaches `ImagPolyExponents` **only for convergent**
  sectors (`!TrueQ[sd["IsDivergent"]]`). Net: a divergent unlifted complex sector
  carries neither Im(A) nor Im(B) phase, and flows into the subtraction /
  `ProcessDivergentSector` machinery (no phase awareness) → **phase-less, wrong
  Laurent, no error**.
- Routing: with `Method->Automatic` the unlifted divergence probe
  (`tropical_eval.wl:3937-3941`) routes to the IBP driver, but `Method->"None"` /
  `"Subtraction"` (and `eps===None`) stay in this driver and hit the gap.

### Fix
Add the unlifted twin of the 4072 guard. In the **unlifted** post-`Table`
section (the `else` arm of `If[isLifted, …]`), after `allSectorData` is built and
the divergent sectors are partitioned, refuse cleanly:

```mathematica
(* planR.md R1: SplitRealImag (Im(A) and/or Im(B)) on an UNLIFTED divergent
   integral is not phase-aware on this inline subtraction path — same limitation
   as the lifted branch (line 4072 / planCXLIFTDIV.md C1).  Refuse cleanly rather
   than emit a phase-less (wrong) result. *)
If[cxSplit && AnyTrue[allSectorData,
       (AssociationQ[#] && TrueQ[#["IsDivergent"]]) &],
  Message[TropicalEval::splitliftdiv];  Return[$Failed]];
```

Place it in the unlifted branch right after the `allSectorData = Table[…]`
completes (mirror the lifted branch's partition-then-guard structure). If the
unlifted branch already partitions into convergent/divergent sets, reuse that
`Select`; otherwise the inline `AnyTrue` above is sufficient.

### Scope decision (confirm before coding)
"Unlifted + divergent + complex exponent" is **not a declared capability**
(`planAXpDIVv2.md` and the lift+divergence work are all lifted). So the correct
behavior is a **clean `$Failed`**, not a new feature. If we later *want* to support
it, that is a separate plan (the IBP route already handles lifted; an unlifted
analog would re-derive the per-sector phase through `ProcessDivergentSector`).
For now: refuse.

### Verification
- **Regression (must `$Failed`, not a number):** add a test feeding an unlifted
  `SplitRealImag` spec with complex `A` (or `B`) and a divergent sector on
  `Method->"None"`; assert the return is `$Failed` and that
  `TropicalEval::splitliftdiv` fired (`Quiet[…, TropicalEval::splitliftdiv]` +
  check via a message-capture wrapper).
- **No-regression:** the existing lifted complex+divergent corpus (`cc_45`,
  `cc_47`, the new `#43` selfgate test) must be **unchanged** — they route through
  the lifted branch / IBP and never reach the new unlifted guard.
- **Byte-identity (#25):** real-A / real-B unlifted divergent integrals have
  `cxSplit=False`, so the new `If` is never entered — confirm a real golden is
  byte-identical.

---

## R2 — Robust Im(A) detection: don't mis-read a symbolic non-`eps` exponent (Medium)

### Diagnosis
`tropical_eval.wl:3877` (and the duplicate at `:5703`):
```mathematica
imAList  = Im[integrandSpec["MonomialExponents"] /. If[eps =!= None, eps -> 0, {}]];
hasImagA = AnyTrue[imAList, (!TrueQ[PossibleZeroQ[#]]) &];
```
If a monomial exponent carries a symbol *other* than the regulator (an undeclared
kinematic/dimension symbol), `Im[sym]` stays **unevaluated**, so
`PossibleZeroQ[Im[sym]]` is `False` → `hasImagA` spuriously `True`. The subsequent
realify `# - I*imAList` then produces a genuinely non-real exponent →
`liftcomplex` `$Failed`, or `Im[sym]` leaks into the double-typed C++ domain
indicator. (Mirrors the pre-existing `imBList` pattern at `:3869`/`:5698`.)

### Fix (pick one; A preferred)
- **A (declare-real):** substitute the declared real symbols to real before `Im[]`,
  e.g. wrap with `Assuming[Element[symbols, Reals], …]` or
  `ComplexExpand[Im[…], targetSymbols]` so only *genuinely* complex literals
  survive. The natural "real symbols" set = `KinematicSymbols ∪ {eps}` (both are
  real). Concretely:
  ```mathematica
  realSyms = Join[Lookup[integrandSpec, "KinematicSymbols", {}],
                  If[eps === None, {}, {eps}]];
  imAList  = Assuming[Element[realSyms, Reals],
               Simplify[Im[integrandSpec["MonomialExponents"]]]];
  ```
  This also removes the `eps->0` hack (eps declared real ⇒ `Im[eps]=0`
  automatically) — but keep `eps->0` if any exponent is `a + b·eps` with complex
  `b` and you specifically want the eps-free imaginary part (it is — preserve the
  current `/. eps->0` semantics; combine: assume reals *then* `eps->0`).
- **B (defensive):** after computing `imAList`, assert every entry is an explicit
  real number; if any entry is still symbolic, `Message` a new
  `TropicalEval::cxmonosym` and `Return[$Failed]` (clean refusal rather than leak).

Apply identically at both call sites (`:3877` and `:5703`) — or fold into the R6
helper so it is written once.

### Verification
- A spec with a free symbol `s` in a monomial exponent (real-intended, no `Im`)
  must give `hasImagA=False` (option A) or a clean `$Failed` with `cxmonosym`
  (option B) — **not** a `liftcomplex` abort or a `cx`-poisoned indicator.
- The existing `#43` test (`A_2 = -1/2 - 3/2 I`, literal complex) must still set
  `hasImagA=True` and pass unchanged.

---

## R3 — Strengthen the `#43` selfgate test (Medium)

### Diagnosis (`TEST/phase2_selfgate.wl`, the new `#43` Module ~lines 9-52 of the diff)
1. **Headline quantity untested.** `payload` extracts `LiftedMonoPhase["Num"]` and
   `MonomialPhaseLog["Coeffs"]` but **never `["Const"]`**. `Const_A =
   (Im(eaug)_p/m_p)·logZ0` is the single most novel term of `planAXpDIVv2` — and
   in this fixture it is **identically 0** (the pivot slot's Im(A) is 0), so the
   `Log[z0]` path is never numerically exercised. The sibling complex-B test
   (~line 204) extracts the *whole* `MonoFactorLog`/`LiftedMonoFactor` association
   incl. `Const`; this test regressed that coverage.
2. **Vacuous pass.** `FreeQ[payload, _Real]` is `True` for `Missing[…]["Num"]` and
   for `{}`. If a regression drops the `LiftedMonoPhase` / `MonomialPhaseLog` key
   or yields empty `IBPTerms`, the assertion passes **green**. No
   `KeyExistsQ` / `=!=None` / `Length>0` guard exists.

### Fix
- **Add a fixture whose `Const_A ≠ 0`** so the `logZ0` term is exercised. Choose a
  lift rule / pivot where `Im(eaug)_p ≠ 0` (i.e. the pivot direction carries
  imaginary monomial weight). Verify `divSec["LiftedMonoPhase"]["Const"]` is a
  nonzero `rational·Log[radical]` and include it in `payload`.
- **Include `Const` in `payload`** for all three pieces:
  ```mathematica
  payload = {
    divSecs[[1]]["LiftedMonoPhase"]["Num"],
    divSecs[[1]]["LiftedMonoPhase"]["Const"],
    ibp["BoundaryData"]["MonomialPhaseLog"]["Coeffs"],
    ibp["BoundaryData"]["MonomialPhaseLog"]["Const"],
    #["MonomialPhaseLog"]["Coeffs"] & /@ ibp["IBPTerms"],
    #["MonomialPhaseLog"]["Const"] & /@ ibp["IBPTerms"]};
  ```
- **Guard against vacuous pass** — assert structure *before* the `FreeQ`:
  ```mathematica
  report["complex-A: LiftedMonoPhase is a populated Association",
         AssociationQ[divSecs[[1]]["LiftedMonoPhase"]] &&
         KeyExistsQ[divSecs[[1]]["LiftedMonoPhase"], "Num"]];
  report["complex-A: IBP produced ≥1 term", Length[ibp["IBPTerms"]] > 0];
  report["complex-A: every IBP term carries MonomialPhaseLog",
         AllTrue[ibp["IBPTerms"],
                 AssociationQ[#["MonomialPhaseLog"]] &]];
  ```
  Only then assert `FreeQ[payload, _Real]`.
- **(Optional, strong) independent oracle.** `planAXpDIVv2` M0 Test E matched
  NIntegrate to `1.2e-8`. Add a small numeric cross-check (pin `epsA`, integrate,
  compare to NIntegrate) so the test validates *correctness*, not just exactness +
  internal self-consistency. Keep it behind the existing M0 tolerance.

### Verification
- Mutating `Const -> 1.5` (a float) anywhere in the payload must now turn the test
  **red** (currently it stays green — confirm the regression first, then fix).
- Deleting the `MonomialPhaseLog` key from a term must turn the test red via the
  new `AllTrue` guard.

---

## R4 — Decouple the Im(A) pole-shift guard from `LiftedMonoFactor` (Low, latent)

### Diagnosis
`tropical_eval.wl:~4985` (`IBPProcessSector` `splitdivmono` guard): the new
`If[lmp =!= None, ckImag += lmp["Num"][[k]]]` sits **inside** the outer
`If[lmf =!= None, …]`. Safe **today** only because `liftedMonoFactorNum` (`lmf`)
is computed unconditionally (`tropical_eval.wl:1339`), so `lmf` always accompanies
`lmp` on divergent lifted sectors. Structurally fragile: any future sector that
sets `LiftedMonoPhase` without `LiftedMonoFactor` would skip the Im(A) Case-B
refusal → silently wrong complex pole.

### Fix
Hoist the check so it fires whenever **either** phase source is present:
```mathematica
If[lmf =!= None || lmp =!= None,
  Module[{imB = Lookup[sectorData, "ImagPolyExponents", None], ckImag = 0},
    If[ListQ[imB] && lmf =!= None,
      ckImag += Sum[imB[[j]] * lmf["DExp"][[j, k]], {j, Length[lmf["DExp"]]}]];
    If[lmp =!= None, ckImag += lmp["Num"][[k]]];
    If[!TrueQ[PossibleZeroQ[ckImag]],
      Message[TropicalEval::splitdivmono, sectorData["ConeIndex"], ckImag];
      Return[$Failed]]]];
```
(Note the `lmf =!= None` added to the `imB` sum guard so `lmf["DExp"]` is never
indexed when `lmf` is absent.)

### Verification
- No behavior change on the current corpus (lmf always present): `cc_45`, `cc_47`,
  `#43` unchanged.
- A synthetic divergent sector with `lmp` set, `lmf` absent, and
  `lmp["Num"][[k]] ≠ 0` must now `$Failed` with `splitdivmono` (previously it would
  have silently proceeded).

---

## R5 — Skip the dead phase term when the Im(A) phase vanishes (Low, hygiene)

### Diagnosis
`emitBaseFuncBody`, `tropical_eval.wl:~2436`: the `If[useMonoPhase, …]` block
appends `cx(0.0,1.0)*(const+coeff)` whenever `monomialPhaseLog =!= None`, with no
per-term `PossibleZeroQ` skip (the Im(B) loop guards each term). A complex-A sector
whose `Im(A).M` flattens to all-zero with `Const=0` still emits
`result *= std::exp(cx(0,1)*(0.))`. Real-exponent goldens are **unaffected** (real
path never enters the block, so #25 holds), but it is a dead op.

### Fix
Guard the append with a zero-check mirroring the Im(B) path:
```mathematica
If[useMonoPhase &&
   (!TrueQ[PossibleZeroQ[mpl["Const"]]] ||
     AnyTrue[mpl["Coeffs"], (!TrueQ[PossibleZeroQ[#]]) &]),
  … append the term …];
```
(Optional cosmetic: the emitted comment header always says "SplitRealImag
oscillatory phase (…)" even in the complex-A / real-B case where B is not split —
the inner `Which[]` label is already accurate, so this is purely a header-string
nicety, defer unless trivial.)

### Verification
- A complex-A sector with a zero flattened phase emits **no** `exp(phaseSum)` line
  for that term — confirm the generated C++ no longer contains `* (0.)`.
- A complex-A sector with a nonzero phase is unchanged (term still emitted).

---

## R6 — Factor the duplicated Im(A) detection + realify idioms (Low, cleanup)

### Diagnosis
- The detection triplet `imAList = Im[… /. eps->0]` + `hasImagA = AnyTrue[…]` is
  copy-pasted in `EvaluateTropicalMC` (`:3877`) and `EvaluateTropicalMCLifted`
  (`:5703`), each mirroring the `imBList`/`hasImagB` pair above it.
- The realify idiom `MapAt[# - I*imAList &, …, {Key["MonomialExponents"]}]` appears
  at **4 sites** (`:3926`, `:4030`, `:4093`, `:5713`), each re-explaining the
  subtract-not-`Re` rationale in a long comment.

### Fix
Introduce two small private helpers near the other internal utilities, and route
R2's hardening through them so the load-bearing `eps->0` / declare-real guard
lives in exactly one place:
```mathematica
(* {imA, hasImagA} for the monomial exponents (eps-free imaginary part). *)
imagMonoInfo[spec_, eps_] := Module[{imA},
  imA = Im[spec["MonomialExponents"] /. If[eps =!= None, eps -> 0, {}]];
  {imA, AnyTrue[imA, (!TrueQ[PossibleZeroQ[#]]) &]}];

(* Realify A by SUBTRACTING its eps-free imaginary part (NOT MapAt[Re], which
   would wrap the real regulator eps in Re[eps] and leak it into the C++). *)
realifyMonoA[spec_, imA_] := MapAt[# - I*imA &, spec, {Key["MonomialExponents"]}];
```
Replace the 2 detection sites and 4 realify sites. The single docstring on
`realifyMonoA` replaces the four repeated rationale comments.

### Verification
- Pure refactor: the lifted + unlifted complex corpus and a real golden must be
  **byte-identical** before/after. Diff the generated C++.

---

## Sequencing & test gate

1. **R1** (+ its `$Failed` regression test) — closes the one wrong-number path.
2. **R3** (strengthen `#43` + add a `Const_A ≠ 0` fixture, optional NIntegrate
   oracle) — so the suite actually guards R1's siblings and the headline term.
3. **R2** — fold the declare-real fix into R6's `imagMonoInfo` if R6 is done first.
4. **R4**, **R5** — local hardening/hygiene.
5. **R6** — cleanup last (or first, then R2/R1 build on the helpers).

**Gate before commit:** full `TEST/phase2_selfgate.wl` green; the lifted complex
corpus (`cc_45`, `cc_47`, `#43`) unchanged; ≥1 real golden byte-identical (#25);
the new R1 regression returns `$Failed` with `splitliftdiv`. No fix lands that
alters real-A / real-B output.

---

## Out of scope (record, do not silently absorb)
- **Supporting** unlifted complex + divergent (R1 refuses it). A real implementation
  (re-derive the per-sector phase through `ProcessDivergentSector`) is a separate
  plan if a physics case needs it.
- The pre-existing symbolic-eps inline Subtraction refusal (`splitliftdiv`,
  Method `"None"/"Subtraction"`) on the **lifted** route is intentional and stays.
