# planAXpDIVv2.md — Complex *monomial* exponents (A_i) × lifting: completing the `SplitRealImag` split

**Goal.** Make the auxiliary-variable **lift** work when the **monomial** exponents
`A_i` (the powers on the bare `x_i^{A_i}`) are complex — the bubble at imaginary
`μ` (principal series). Today the lift aborts on every cone with
`TropicalEval::liftcomplex`; this plan removes that abort by extending the
existing `ComplexExponentMode -> "SplitRealImag"` machinery to split the monomial
exponents on the **same footing** as the polynomial exponents `B_j`: decompose on
`Re(A)` (so the domain indicator is real) and carry `Im(A)` into the *same*
oscillatory phase that already transports `Im(B)`. Diagnosis and the one-character
confirmation are in `LIFTCOMPLISSUE/liftcomplex_issue.tex`.

**Thesis (no new math beyond planCXLIFTDIV — it is the symmetric twin).**
`planCXLIFTDIV.md` split the *polynomial* exponents `B_j = Re(B_j) + i·Im(B_j)`:
the real part builds the importance measure / domain indicator, the imaginary
part becomes `exp(i·Im(B_j)·(log Q_j + MonoFactorLog_j))`. The exact analog holds
for the *monomial* exponents. The post-δ effective exponent is **linear in both
families** —
```
ã_a = rawA_a + Σ_k B_k d̃_{k,a},   rawA_a = Σ_i (A_i+1) M_{ia}      (tropical_eval.wl:384, 422–423)
```
— so `Im(A)` contributes `Σ_i Im(A_i) M_{ia}` to `Im(ã_a)`, exactly as `Im(B)`
contributes `Σ_k Im(B_k) d̃_{k,a}`. Taking `Re` of **both** families makes `ã`
real and the half-space `Boole[log y_p* ≤ 0]` well-posed; the leftover `Im(A)`
rides an oscillatory phase on the real measure. The note proves the decomposition
half empirically: feeding `Re` of *both* `A` and `B` to `ProcessSectorLifted`
turns the `bubblepm` `μ=3i/2` point from **611/611 cones failing** (`liftcomplex`)
to **611/611 cones succeeding**. Flattening is untouched (it is geometric, fixed
by the cone rays / support — `LIFTCOMPLISSUE` §"What is and is not broken"). What
remains is wiring `Im(A)` through the phase codegen, the monomial analog of
planCXLIFTDIV §3.

Written to slot beside `plan.md` / `planAXpDIV.md` / `planCXLIFTDIV.md`: same
exactness discipline (§1.4 / §6.7), the same "≥2 independent oracles per numeric
PASS" rule (§8.5), and the same "a known limitation produces a clean `$Failed`,
never a wrong number".

---

## M0 — what the sandbox verified (DONE)

`SANDBOX/sandbox_M0_complexmono.wl` and `SANDBOX/sandbox_M0_lifted.wl` (run
2026-06-29) settle the plan's assumptions on self-contained complex-*monomial*
toys (real poly coeffs; `Im(A)={3/2,−3/2}`, bubble-like). Results:

| test | claim | result |
|---|---|---|
| **A** (lifted) | `Re(B)`-only ⇒ `liftcomplex` every cone; `Re(A)+Re(B)` ⇒ all succeed | **CONFIRMED** — 4/4 fail vs 4/4 succeed (reproduces `LIFTCOMPLISSUE` table; engine's 611/611, here at small scale) |
| **D'** (unlifted) | monomial phase `exp(i·Σ_a (Im(A)·M)_a/ã_a·log y'_a)` is correct & load-bearing | **CONFIRMED** — matches NIntegrate `<1%`; dropping it `>5%` off |
| **E** (lifted) | lifted phase **incl. the new `Const_A`** (the `log z0` term, §3) is correct & load-bearing | **CONFIRMED** — matches NIntegrate to `1.2×10⁻⁸`; dropping it 18× off |
| **C** (unlifted) | *(earlier draft)* unlifted `SplitRealImag` silently mis-evaluates complex `A` | **REFUTED** — agrees with `Direct`/NIntegrate to `5×10⁻⁴`; no bug (plan corrected) |

**Net:** the primary fix (realify `A` ⇒ clears `liftcomplex`) and the entire new
math (the monomial phase, **including the lifted `Const_A` closed form**) are
empirically verified before implementation. The one course-correction: there is
no unlifted correctness bug — the unlifted path is exact for any split (§1).

---

## 1. Where we are today (the half-split)

`SplitRealImag` re-implements only **one** of the two inputs to `ã`. Every site
that prepares a spec for decomposition takes `Re` of the **polynomial** exponents
and leaves the **monomial** exponents complex:

| site (`tropical_eval.wl`) | driver / purpose | current |
|---|---|---|
| `3814` | `EvaluateTropicalMC` — `routeToIBP` divergence detection (lifted) | `MapAt[Re, liftedSpec, {Key["PolynomialExponents"]}]` |
| `3912` | `EvaluateTropicalMC` — Step-2 **lifted** convergent processing | `MapAt[Re, specForProc, {Key["PolynomialExponents"]}]` |
| `3968` | `EvaluateTropicalMC` — Step-2 **unlifted** convergent processing | `MapAt[Re, specToUse, {Key["PolynomialExponents"]}]` |
| `5545` | `evaluateTropicalIBPDriver` — **lifted divergent** `specForProc` | `MapAt[Re, integrandSpec, {Key["PolynomialExponents"]}]` |

and the `cxSplit` gate that turns the mode on is keyed on `Im(B)` only, at **two**
sites: `cxSplit = (cxMode==="SplitRealImag") && hasImagB` in `EvaluateTropicalMC`
(`tropical_eval.wl:3774` — covers convergent **lifted and unlifted**), and
`cxSplit = isLifted && (cxMode==="SplitRealImag") && hasImagB` in
`evaluateTropicalIBPDriver` (`:5543` — the divergent route). `EvaluateTropicalMCLifted`
(`:4382`) is a thin wrapper that delegates to `EvaluateTropicalMC` (`:4463`), so it
inherits whatever the main driver does.

Consequences — **one** is a correctness failure, the other is not (both M0-verified, §M0):

- **Lifted path — loud failure (`liftcomplex`), the real bug.** With `A` left complex,
  `rawA = (A+1).M` is complex, so `ã = rawA + Σ B^R d̃` is complex from the `A`
  half. `ProcessSectorLifted`'s per-slot classifier `classDir` hits its first
  branch — `!TrueQ[PossibleZeroQ[Im[a0]]] → "complex"` (`tropical_eval.wl:1105`)
  — on every admissible pivot, and the guard fires
  (`tropical_eval.wl:1183`, message text at `:209`). The bubble at imaginary `μ`
  is the **first** integrand anyone has lifted whose *monomial* exponents are
  complex (`Im(A) = {−3/2, 3/2, −3/2, 3/2, 3}`, `n=5`, 611 cones); all prior
  complex tests (`cc_21b`, `cc_45`) had complex `B` with **real `A`**, so the
  half-split was invisible. The message itself betrays the blind spot — it says
  *"check that polynomial exponents B are real"*, never mentioning `A`.

- **Unlifted path — correct in value, only variance-suboptimal (NOT a bug).**
  The unlifted `ProcessSector` has **no** realness guard (`tropical_eval.wl:287`
  tests only `Re[a0]≤0`), so it never refuses complex `A`. An unlifted
  `SplitRealImag` call with complex `A` builds a *complex* `a_eff` (only `Re(B)`
  was realified) and attaches only `Im(B)` to the phase — yet the result is
  **still numerically correct**, because the tropical change of variables
  `y_i=(y'_i)^{1/a_eff_i}` is an **exact identity for any `a_eff`**, real or
  complex; the real/imag split only moves variance between "measure" and "phase",
  never the expected value. *(M0 Test C: unlifted `SplitRealImag` agrees with
  `Direct` and NIntegrate to `5×10⁻⁴`.)* This matches `LIFTCOMPLISSUE` §2 ("the
  unlifted scan completed normally"). An **earlier draft of this plan wrongly
  called this a silent-wrong bug; M0 refuted that.** Realifying `A` in the
  unlifted path is therefore an *optional* variance improvement (a real measure
  for VEGAS), not a correctness fix — and it is only safe **paired with** the
  `Im(A)` phase (§4.4): realifying `A` *without* the phase would turn a correct
  result into a wrong one (M0 Test D' without-phase: `>5%` off; lifted Test E:
  18× off). The **only** correctness failure is the lifted `liftcomplex` above.

(`Direct` mode is already correct for complex `A` *unlifted* — `exp(A·log x)` via
the branch-cut-free log-exp form, SUMMARY §2e. `Direct` **+ lift** is refused for
*any* complex exponent — a complex `ã` has no real domain indicator — so
`SplitRealImag` is the only mode that can lift a complex integrand at all; this
matches planCXLIFTDIV's conclusion for `B`.)

---

## 2. Why there is no fundamental obstruction (the math)

Write `A_i = A_i^R + i A_i^I` and `B_j = B_j^R + i B_j^I`. The flattened effective
exponent (`tropical_eval.wl:384`, `:422–423`) splits cleanly:

```
ã_a = ã_a^R + i·Im(ã_a),
  ã_a^R   = Σ_i (A_i^R + 1) M_{ia}  +  Σ_k B_k^R d̃_{k,a}      ← real measure / domain indicator
  Im(ã_a) = Σ_i A_i^I M_{ia}        +  Σ_k B_k^I d̃_{k,a}      ← oscillatory phase
              └─ NEW (monomial) ─┘      └─ existing (poly) ─┘
```

The sector integrand `∏_a y_a^{ã_a − 1} ∏_j Q_j^{B_j}` factors as **(real
measure)·(phase)**:

```
∏_a y_a^{ã_a^R − 1} ∏_j Q_j^{B_j^R}                         (flattened by ã^R — the existing real path)
  ×  exp( i Σ_a (Σ_i A_i^I M_{ia}) log y_a                  ← NEW: bare-monomial phase
        + i Σ_j B_j^I  log Q_j  +  (dropped tropical monomial factor) ).   ← existing MonoFactorLog phase
```

Two facts make this tractable:

1. **The domain indicator becomes real.** Because the measure is flattened by
   `ã^R` (real), the indicator coefficient `ic_a = mOther_a / ã_a^R`
   (`tropical_eval.wl:1169–1185`) is real and `Boole[log y_p* ≤ 0]` is
   well-defined. This is exactly the note's 611/611 result, and it is *all* of
   the decomposition-side fix: realify `A` as well as `B` before
   `ProcessSectorLifted`.

2. **The new phase term is the monomial twin of `MonoFactorLog`.** Re-express the
   bare-monomial phase in flattened coordinates with `log y_a = log y'_a / ã_a^R`:
   ```
   exp( i Σ_a [ (Im(A)·M)_a / ã_a^R ] log y'_a ),     (Im(A)·M)_a := Σ_i A_i^I M_{ia}.
   ```
   The numerator `(Im(A)·M)_a` uses `A` (not `A+1`: the `+1` is the **real**
   Jacobian/measure, not part of the imaginary phase). When `Im(A)=0` the term
   vanishes identically ⇒ the existing output is reproduced byte-for-byte.

The divergence machinery (IdentifyDivergences / IBPReduceSector / pole assembly)
operates on `Re(ã)` and is unchanged — the pole still lives in the
`Re(ã_k) = c_k ε → 0` direction, exactly as in the real-`A` planAXpDIV case.

---

## 3. The only new math — per-piece `MonomialPhaseLog`

Just like `MonoFactorLog` (planCXLIFTDIV §3), the new phase term must be
**re-derived per IBP piece**, because IBP shifts monomial exponents and
re-flattens each boundary/term with its own `α0^{(piece)}`. The structure that
makes this cheap is the same: the **numerator is shared, only the flattening
divisor differs**.

```
MonomialPhaseLog_a^{(piece)} = Const_A  +  (Im(A)·M)_a / α0_a^{(piece)} · log y'_a,
  (Im(A)·M)_a = Σ_i Im(A_i) M_{ia}            ← sector-level, shared across boundary + all terms
  α0_a^{(piece)}                              ← that piece's flattening (boundary drops slot k)
```

- `M` (the cone change-of-variables matrix) and `Im(A)` are sector-level, so the
  numerator `(Im(A)·M)_a` is computed **once** per sector; each IBP piece just
  re-divides it by its own `α0` — identical to how `MonoFactorLog`'s numerators
  `D_{j,i}` are re-divided (planCXLIFTDIV §4.1–4.2). At the boundary (`y_k = 1`),
  `log y_k = 0`, so slot `k` drops.

- **Lifted `Const_A` (the genuinely new bookkeeping) — DERIVED & VERIFIED (M0 Test E).**
  Apply the engine's pivot substitution `log Y_p = (logZ0 − Σ_jj mOther_jj log y_jj)/m_p`
  (the same one used for the `B`-phase, `tropical_eval.wl:1042–1090`, `:1286–1297`)
  to the bare monomial `Σ_i Im(A_i) log x_i = Σ_a Im(eaug)_a log Y_a`, where
  `Im(eaug) = Im(A_aug)·M_aug` (`A_aug = {Im A_1,…,Im A_n, 0}`, aux appended with 0;
  `M_aug = RayMatrix`, the augmented CoV). This gives, in flattened coords:
  ```
  Const_A     = ( Im(eaug)_p / m_p ) · logZ0
  Coeffs_jj   = ( Im(eaug)_{rIdx[jj]} − Im(eaug)_p · mOther_jj / m_p ) / ã_jj
  ```
  with `p = PivotIndex`, `m_p = ZRow[[p]]`, `rIdx = DeleteCases[Range[n+1], p]`,
  `mOther = ZRow[[rIdx]]`, `ã = NewExponents` — **all exposed sector keys**. This
  is exactly the lifted `MonoFactorLog` structure (`:1286–1297`) with the
  polynomial minima `d^aug` replaced by `Im(eaug)` and **no `rcMin` term** (the
  bare monomial is not a re-cleared polynomial). For **unlifted** sectors `m_p`/z0
  are absent and this reduces to `Const_A = 0`, `Coeffs_a = (Im(A)·M)_a/ã_a` (§2).
  *M0 Test E confirmed this lifted form to `1.2×10⁻⁸` vs NIntegrate (18× off if
  dropped), so the closed form above is the implementation target — not a TODO.*

Nothing else in the Laurent assembly changes: `Re(A)`, `Re(B)` drive `ã^R`,
`c_k`, `r_k`, `PrefactorBase`, `DLogPrefactor` (all real); `Im(A)` and `Im(B)`
ride the phase.

---

## 4. Design — data-flow changes

Mirrors the `Im(B)` plumbing already in place; **every change is gated on
`Im(A) ≠ 0` and reduces to the current code when `Im(A)=0`**, so real-`A` and
prior complex-`B` goldens stay byte-identical (#25).

### 4.1 Detection — gate `SplitRealImag` on `Im(A)` too
At **both** `cxSplit` sites (`:3774`, `:5543`) extend the `Im(B)` test so the mode
also engages when the monomial half is complex, **preserving each site's existing
structure** (the IBP-driver gate keeps its `isLifted &&` conjunct — unlifted
*divergent* integrals correctly fall back to `Direct`, out of scope here):
```
imAList  = Im[integrandSpec["MonomialExponents"]];
hasImagA = AnyTrue[imAList, (!TrueQ[PossibleZeroQ[#]]) &];
cxSplit  = … && (hasImagB || hasImagA);   (* :3774 as-is; :5543 keeps "isLifted &&" *)
```
Without this, a complex-`A`/real-`B` integrand would silently skip the split.

### 4.2 Decomposition side — take `Re(A)` wherever `Re(B)` is taken
At each of the four sites in §1 (`3814`, `3912`, `3968`, `5545`) add the monomial
realification next to the existing polynomial one:
```
specForProc = MapAt[Re, specForProc, {Key["PolynomialExponents"]}];
If[hasImagA, specForProc = MapAt[Re, specForProc, {Key["MonomialExponents"]}]];
```
**This single change clears `liftcomplex` (the note's 611/611).** Grep-audit for
any other `MapAt[Re, …, {Key["PolynomialExponents"]}]` that may be added later and
keep the two in lockstep.

> **Coupling (M0-critical): §4.2 must not ship without §4.4.** Realifying `A`
> *without* attaching the `Im(A)` phase produces a confidently **wrong** number
> (M0 Test D' / E without-phase: `>5%` / 18× off), because it drops the `Im(A)`
> oscillation that was previously carried (lifted: by failing loudly; unlifted: by
> staying in the complex measure). So M1 (realify) and M2 (phase) are **one
> atomic correctness unit** — M1's gates below assert only `0 liftcomplex` +
> byte-identity, *not* a correct value; correctness is gated at M2. Do not land
> M1 to `main` without M2.

### 4.3 Sector data — attach `Im(A)` and the monomial numerator
Symmetric to `ImagPolyExponents` / `MonoFactorLog`:
- Attach `"ImagMonoExponents" -> imAList` to **every** sector the driver builds
  wherever it currently attaches `ImagPolyExponents`: `tropical_eval.wl:3926`
  (lifted convergent) and `:3971` (unlifted convergent) in `EvaluateTropicalMC`,
  and `:5598` (lifted) in `evaluateTropicalIBPDriver`.
- In `ProcessSector` / `ProcessSectorLifted`, alongside `MonoFactorLog`
  (`:485–489` unlifted, `:1275–1339`/`:1379` lifted), store the **un-flattened**
  monomial-phase numerator and the lifted constant:
  ```
  "MonomialPhaseLog" -> <| "Num" -> (Im(A) . M),        (* length n, shared       *)
                           "Const" -> Const_A |>         (* 0 unlifted; §3 lifted  *)
  ```
  Absent on real-`A` sectors (codegen treats absent as "no monomial phase").
- IBP per-piece: in `IBPProcessSector` re-derive `MonomialPhaseLog^{(piece)}` for
  the boundary and each term by dividing `Num` by that piece's `α0^{(piece)}`
  (boundary: drop slot `k`), attaching it next to the per-piece `MonoFactorLog`
  (planCXLIFTDIV §4.2 sites).

### 4.4 Codegen — one more phase term in `emitBaseFuncBody`
`emitBaseFuncBody[…, monoFactorLogs_:None, imagExps_:None]`
(`tropical_eval.wl:2303–2377`) already builds the `Im(B)` phase block
(`useSplit` test `:2312`, `phaseSum` assembly `:2354–2377`). Add a parallel
optional pair `(monomialPhaseLog_:None, imagMonoExps_:None)` and, when
`Im(A) ≠ 0`, append one term to the **same** `phaseSum`:
```
// monomial oscillatory phase (Im(A); flattened bare powers + lift constant)
phaseSum += cx(0.0, 1.0) * ( Const_A + Σ_a (Num_a/α0_a) * log_y[a] );
```
Thread the two new arguments through **every** `emitBaseFuncBody` call site
(convergent `:2461–2468`; G0/G1 `:2486`,`:2509`; IBP boundary `:2641–2674`; IBP
terms `:2695–2722`; and the hand-written remainder integrand) using
`Lookup[piece,"MonomialPhaseLog",None]` / `Lookup[sd,"ImagMonoExponents",None]`.
With `Im(A)=None` the term is skipped ⇒ **byte-identical** C++ (#25). One emitter
now carries both phases, so the convergent and divergent routes are covered by
the same code path.

### 4.5 Both routes, both phases
- **Convergent route (e.g. `bubblepm`):** the convergent emitter (`:2461`) now
  passes both `MonoFactorLog` and `MonomialPhaseLog`; nothing else changes.
- **Divergent route (e.g. `bubblepp`):** the IBP boundary/term emitters
  (`:2641`,`:2695`) and `IBPProcessSector` already re-derive per-piece
  `MonoFactorLog` (planCXLIFTDIV C1); they now re-derive `MonomialPhaseLog` the
  same way. The pinned-ε `LaurentFromSubtraction` route inherits it for free
  (every sector is convergent at fixed ε, so the convergent emitter runs). The
  symbolic-ε inline Subtraction path stays refused by the existing
  `splitliftdiv` (`:215`) — unchanged scope.

### 4.6 `ValidateLiftedDecomposition` / `ValidateIBP` / `ValidateSubtraction`
The WL NIntegrate reconstructions that already multiply in the `Im(B)` phase must
also multiply in `exp(i·(Const_A + Σ_a (Num_a/α0_a) log y'_a))`, so the exact
validators remain a true oracle for the emitted C++.

---

## 5. Exactness discipline (plan.md §1.4 / §6.7)

- `Re(A)`, `Im(A)`, `Im(A).M`, `Const_A`, `ã^R` are built with `Re[…]`/`Im[…]`
  and **exact** arithmetic (`M` integer, `z0` exact radical) — **no `N[…]`**
  before the `MmaToC` boundary; numericization happens only inside
  `mmaToCInternal`, as for `Im(B)`/`MonoFactorLog`.
- The classification stays exact: `Re(A)` realification means
  `PossibleZeroQ[Im[ã]]` is now provably `True` on every slot in `SplitRealImag`
  mode, so `classDir` (`:1105`) never reaches the `"complex"` branch — `liftcomplex`
  becomes unreachable in this mode *by construction*, not by a tolerance test.
- Extend the #43 corpus with a **complex-`A`** lifted spec and assert
  `FreeQ[symbolicOutput, _Real]` on `Im(A).M`, `Const_A`, and the `ã^R`-derived
  fields before `MmaToC`.

---

## 6. Phased plan (with gates; ≥2 oracles each, none judged by its producer)

**Phase M0 — Sandbox the math. ✅ DONE (see the "M0 — verified" section above).**
`SANDBOX/sandbox_M0_complexmono.wl` + `sandbox_M0_lifted.wl` confirmed, on
complex-`A` toys: (a) `Re(A)+Re(B)` clears `liftcomplex` (Test A); (b) the lifted
`Const_A` closed form (§3) matches NIntegrate to `1.2×10⁻⁸` (Test E); (c) the
monomial phase is load-bearing, unlifted and lifted (Tests D'/E); and it **refuted**
the unlifted silent-wrong hypothesis (Test C). *Remaining for M1/M2 (not blocking):*
re-run the same two assertions on the **actual** `bubblepm` `μ=3i/2` spec (lives in
read-only `OLD_CODE/.../CC`; not yet wired into a TEST) to reproduce the literal
611/611 and the `≈ 1.27e-14 + 1.40e-15 i` reference — the mechanism is already
proven, this is fidelity to the physics point.

**Phase M1 — Detection + decomposition side (clears `liftcomplex`).** §4.1–4.2.
*Gate:* (a) real-`A` and prior complex-`B` goldens byte-identical (#25); (b) the
`bubblepm` `μ=3i/2` lifted decomposition produces **0** `liftcomplex` across all
611 cones; (c) #43 exactness guard passes on a complex-`A` lifted spec.

**Phase M2 — Monomial phase, convergent route (`bubblepm`).** §4.3–4.4, §4.6 for
the convergent emitter + `ValidateLiftedDecomposition`. *Gates:*
- **byte-identity:** every existing convergent `SplitRealImag` golden (incl.
  `cc_21b`) unchanged — the new term is absent when `Im(A)=0`.
- **#A-conv (vs NIntegrate):** lifted `SplitRealImag` reproduces NIntegrate of the
  original `bubblepm` integrand to < 1–2% (mirrors `cc_21b` validation a).
- **#A-load (load-bearing):** dropping the monomial phase term gives a *wrong*
  answer (mirrors `cc_21b` validation b).
- **#A-unlifted (regression — stays correct):** an **unlifted** complex-`A`
  `SplitRealImag` call still agrees with `Direct`/NIntegrate (it already did — M0
  Test C) *and* now decomposes on a real measure (variance, not correctness). This
  guards against §4.2-without-§4.4 silently corrupting the previously-correct
  unlifted result.

**Phase M3 — Divergent route (`bubblepp`), IBP + pinned-Subtraction.** §4.3 IBP
per-piece, §4.5, §4.6 for `ValidateIBP`/`ValidateSubtraction`. *Gates:*
- **#A-IBPvsSub:** complex `(pole, finite)` of `bubblepp` agree between the IBP
  route and the pinned-ε Subtraction route < 0.5%.
- **#A-vsNI:** reconstructed `I(ε*)` (complex) vs direct NIntegrate at
  `ε* ∈ {0.02, 0.01}` within tolerance.
- **#A-real-reduction:** with `Im(A)=0` the divergent route reproduces the
  planAXpDIV / planCXLIFTDIV real-`A` numbers exactly (phase term absent).

**Phase M4 — docs + cross-checks + message fix.** §7–8. Update the `liftcomplex`
message, SUMMARY, plan.md, and the manual; land `cc_47`.

---

## 7. Cross-checks (≥2 oracles each)

- **`cc_47` (new): complex-monomial-exponent bubble.** Both presectors of
  `bubblepm` (convergent) and `bubblepp` (divergent, IBP) at `μ = 3i/2`. Oracles:
  lifted `SplitRealImag` (this plan), pinned-ε Subtraction route, and NIntegrate
  of the original integrand. PASS = convergent value within 2% of NIntegrate
  (`≈ 1.27e-14 + 1.40e-15 i`); divergent complex pole agreement < 0.05 across
  routes. Include the **611/611 cone-success** assertion as a structural check.
- **`cc_21b` / `cc_45` regression:** unchanged (complex-`B`, real-`A`) — confirms
  the monomial term is inert when `Im(A)=0`.
- **Unlifted complex-`A` `SplitRealImag` vs `Direct`:** a small toy, agreement to
  MC tolerance (regression that the already-correct unlifted result is preserved —
  M0 Test C; *not* a bug fix).
- **#43 corpus:** complex-`A` lifted+divergent exact spec; `FreeQ[_, _Real]` on
  `Im(A).M`, `Const_A`, `ã^R` fields before `MmaToC`.
- **CUBA (#5) where available:** Cuhre on the direct complex integrand at fixed ε.

---

## 8. New messages / limitation updates

The current guard text is now **misleading** (it blames only `B`); fix it:
```
TropicalEval::liftcomplex =                                     (* tropical_eval.wl:209 *)
  "ProcessSectorLifted: cone `1` — all candidate pivots produce complex atilde;
   cannot emit a real-valued domain indicator.  In Direct mode, complex monomial
   (A) OR polynomial (B) exponents are unsupported with lifting — use
   ComplexExponentMode -> \"SplitRealImag\" (splits BOTH A and B: real parts build
   the domain map, imaginary parts ride the oscillatory phase)."
```
In `SplitRealImag` mode this message becomes **unreachable** (Re of both families
⇒ real `ã`); it remains the correct refusal for `Direct` + lift.

Update **SUMMARY.txt §2e / §6 L2(iii)**, **plan.md N3 / §9 D3**, the **manual**
(complex-exponent + lifting limitation), and **planCXLIFTDIV §9** from "complex
*polynomial* exponents with lifting supported via `SplitRealImag`" to: *"Complex
**monomial and polynomial** exponents with lifting (and with a `1/ε` pole) are
supported via `ComplexExponentMode -> "SplitRealImag"`, which splits both `A` and
`B` (cc_47); `Direct` mode with lifting remains refused for any complex exponent."*

---

## 9. Scope boundaries

**Delivers:** lifting with **complex monomial exponents `A`** (and complex `B`,
and a single `1/ε` pole), via `SplitRealImag`, on the convergent route
(`bubblepm`) and the divergent IBP + pinned-ε Subtraction routes (`bubblepp`);
the principal-series bubble; full exactness discipline; byte-identical real-`A`
and prior complex-`B` goldens; and, for the unlifted path, a real importance
measure for complex `A` (a **variance** improvement — the value was already
correct, M0 Test C).

**Does not deliver (explicit future work, refused with a clean `$Failed`):**
`Direct` mode with lifting for complex exponents (use `SplitRealImag` — inherits
planCXLIFTDIV C3 scope); the symbolic-ε **inline** Subtraction path with
lift × divergence × phase (`splitliftdiv`, `tropical_eval.wl:215`); Case B
(domain–pole coupling, planAXpDIV §3); `1/ε^{d≥2}` and multiple independent
auxiliary variables. **Unlifted-divergent** `SplitRealImag` stays gated off by the
`isLifted &&` conjunct at `:5543` (it falls back to the always-correct `Direct`
product — variance-suboptimal but never wrong), unchanged pre-existing scope. The
single-pole-per-sector limit and the `HasConstantTerm=False` route guidance (use
the Subtraction route) carry over unchanged.
