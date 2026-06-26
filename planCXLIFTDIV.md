# planCXLIFTDIV.md — Complex exponents × auxiliary lifting × divergence

**Goal.** Let an integral be **lifted**, **divergent** (a simple `1/ε` pole), **and**
carry **complex polynomial exponents** `B_j` — all in one call. This is the last
combination the v3 + planAXpDIV engine refuses. It closes the gap flagged as
future work in planAXpDIV.md §4.5 (`splitliftdiv`) and §4 (`liftcomplex` with
lifting), and removes limitation **L2(iii)** ("complex-B lifting inherits
`liftcomplex`").

**Thesis (no fundamental obstruction; mostly an implementation gap).** Lifting
acts on a coefficient, divergence on a monomial exponent, and a complex `B`
contributes an imaginary part to the post-δ effective exponents `ãtilde`. The
imaginary part is an *oscillatory phase* on a real measure — exactly what the
existing `SplitRealImag` mode already factors out (decompose on `Re(B)`,
reintroduce `Im(B)` as `exp(i·Im(B)·(log P + MonoFactorLog))`). cc_21b already
verifies this for **lifted, convergent** complex integrals. Two facts then make
the combination tractable:

1. **The Subtraction route already delivers it (verified).** `LaurentFromSubtraction`
   pins `ε` to several small values, which makes *every* lifted sector convergent
   (the divergent slot `ãtilde_k = c_k·ε > 0` at fixed `ε`), and the lifted
   `SplitRealImag` path then runs unchanged. Empirically, for
   `∫ x1^{-1+ε}(1 + x1 x2 + 10^6 x1^2 + x2^2)^{-2+0.3i} dx` with lift `{2,0},k=3`:

   ```
   ComplexExponentMode -> "SplitRealImag", route = LaurentFromSubtraction
     pole   = 0.7699 + 0.0882 i
     recon F(0.01) = 71.476 + 8.293 i   vs  NIntegrate = 71.612 + 8.301 i
     rel-err = 0.19%
   ```

   So **complex-B lift+divergence is, today, a one-liner via the Subtraction
   route** (this plan only formalizes/pins it with a cross-check).

2. **Only two paths still refuse it**, and both are implementation choices:
   - the **IBP route** (`Method -> "IBP"`) refuses via `TropicalEval::splitliftdiv`
     because the oscillatory-phase `MonoFactorLog` is sector-level, but IBP shifts
     exponents and re-flattens *per term/boundary*, so the phase coefficients must
     be **re-derived per IBP term/boundary**; and
   - the **Direct mode** (`ComplexExponentMode -> "Direct"`) refuses *any* complex-`B`
     lifting via `TropicalEval::liftcomplex`, because a complex `ãtilde` yields a
     complex domain-indicator coefficient `ic_j = m_other,j / ãtilde_j` and the
     half-space `Boole[log y_p* ≤ 0]` needs a **real** linear inequality.

This plan delivers the IBP route (Phase C1 — the principled second oracle),
formalizes the already-working Subtraction route (Phase C0), and optionally
unblocks Direct mode (Phase C3 — the genuine geometric subtlety, scoped honestly
and refused if not resolved). Written to slot beside `plan.md` / `planAXpDIV.md`:
same exactness discipline (§1.4 / §6.7), same "≥2 independent oracles per numeric
PASS" rule (§8.5), same "a known limitation produces a clean `$Failed`, never a
wrong number".

---

## 1. Where we are today (the three modes × the obstruction)

| route | mode | complex-B lift+div | why |
|---|---|---|---|
| `LaurentFromSubtraction` | `SplitRealImag` | **works** (verified, §thesis) | pinned ε ⇒ convergent ⇒ lifted SplitRealImag (cc_21b) |
| `LaurentFromSubtraction` | `Direct` | refused | pinned-ε sectors still have complex `ãtilde` ⇒ `liftcomplex` |
| `Method->"IBP"` (`evaluateTropicalIBPDriver`, lifted) | `SplitRealImag` | refused | `splitliftdiv`; per-term MonoFactorLog not re-derived |
| `Method->"IBP"` | `Direct` | refused | `liftcomplex` (complex `ãtilde`) |

`ProcessSectorLifted` classifies a sector with complex `ãtilde` as `"complex"` and
returns `$Failed` (`liftcomplex`); the **driver** in `SplitRealImag` mode decomposes
on `Re(B)` *before* calling it (`specForProc = MapAt[Re, liftedSpec,
{Key["PolynomialExponents"]}]`), so `ãtilde` is **real** and `liftcomplex` never
fires — but the divergent IBP codegen never receives the `Im(B)` phase, and the
driver guards the combination with `splitliftdiv`.

---

## 2. Why there is no fundamental obstruction (the math)

Write `B_j = B_j^R + i B_j^I`. The post-δ effective exponent of a lifted sector,
`ãtilde_j = a_other,j − a_p m_other,j/m_p + Σ_k B_k d̃_{k,j}`, splits as
`ãtilde_j = Re(ãtilde_j) + i·Im(ãtilde_j)` with
`Re(ãtilde_j) = a_other,j − a_p m_other,j/m_p + Σ_k B_k^R d̃_{k,j}` (the value the
driver gets from `Re(B)`) and `Im(ãtilde_j) = Σ_k B_k^I d̃_{k,j}`.

The sector integrand `∏_i y_i^{ãtilde_i − 1} ∏_j Q_j^{B_j}` factors as
**(real measure)·(oscillatory phase)**:

```
∏_i y_i^{Re(ãtilde_i) − 1} ∏_j Q_j^{B_j^R}     (flattened by Re(ãtilde) — the existing real path)
  ×  exp( i Σ_i Im(ãtilde_i) log y_i  +  i Σ_j B_j^I log Q_j ).
```

The phase's `log y_i` are re-expressed in flattened coordinates plus the *dropped
tropical monomial factor* — this is precisely the `MonoFactorLog` that
`SplitRealImag` already compiles (cc_21b sub-checks D/E confirm the lifted
`MonoFactorLog` is present and necessary). Because the **measure** is flattened by
`Re(ãtilde)` (real), the **domain indicator** `ic_j = m_other,j / Re(ãtilde_j)` is
real and the half-space `Boole` is well-defined. So in `SplitRealImag` the only
thing missing for divergence is wiring the phase through the divergence codegen.
The divergence machinery itself (IdentifyDivergences / IBPReduceSector / the pole
assembly) operates on `Re(ãtilde)` and is unchanged — the pole lives in the
`Re(ãtilde_k) = c_k ε → 0` direction exactly as in planAXpDIV (real-B) case.

---

## 3. The IBP subtlety — per-term `MonoFactorLog` (the only new math)

IBP reduces a divergent sector into a **boundary** integral (at `y_k = 1`, `n−1`
dims) and a set of **IBP terms** (each `n`-dim, with shifted monomial exponents
`α^{(t)}` and its own flattening `α0^{(t)} = α^{(t)}|_{ε→0}`). The cleared
polynomials `Q_j` — and therefore the tropical monomial factor `D_{k,j}` that
clearing removed from `Q_k` — are **shared** across boundary and all terms; only
the flattening differs. Hence the oscillatory phase of each piece is

```
exp( i Σ_j B_j^I ( log Q_j  +  MonoFactorLog_j^{(piece)} ) ),
   MonoFactorLog_j^{(piece)} = Const_j + Σ_i ( D_{j,i} / α0^{(piece)}_i ) log y'_i,
```

with `Const_j = (d^aug_{j,p}/m_p) log z0` (the same `z0` factor as the sector) and
`D_{j,i}` the un-flattened monomial-factor exponents (the numerators in
`ProcessSectorLifted`'s sector-level `MonoFactorLog`). So **re-deriving the phase
per piece is just re-dividing the fixed `D_{j,i}` by that piece's `α0`** (the
boundary drops the divergent slot `k`; at `y_k = 1` its `log y_k = 0`).

This is the entire content of `splitliftdiv` "per-term MonoFactorLog
re-derivation is future work". Nothing else in the IBP Laurent assembly changes:
`Re(B)` drives `c_k, r_k`, `PrefactorBase`, `DLogPrefactor` (all real); `Im(B)`
rides on the phase.

---

## 4. Design — data-flow changes

### 4.1 Carry the un-flattened monomial-factor data on a divergent lifted sector
`ProcessSectorLifted`, divergent branch: add a field
`"LiftedMonoFactor" -> <|"Const" -> {Const_j}, "DExp" -> {D_{j,i}}|>` (the
**un-flattened** numerators, length `n` per polynomial), so each IBP piece can
re-flatten. For unlifted/real sectors this is absent (codegen treats absent as
"no phase"), preserving byte-identity.

### 4.2 Re-derive per-piece `MonoFactorLog` in `IBPProcessSector`
For the boundary and each IBP term, build
`MonoFactorLog_j^{(piece)} = <|"Const"->Const_j, "Coeffs"->D_{j,i}/α0^{(piece)}_i|>`
(boundary: `i` over the `n−1` non-divergent coords, slot `k` dropped; terms: `i`
over all `n`). Attach to `BoundaryData` and each `IBPTerms[[t]]` as a
`"MonoFactorLog"` key, alongside the already-present per-piece `FlatPolys`/`Alpha0`.

### 4.3 Thread `Im(B)` + per-piece `MonoFactorLog` into the divergent codegen
`emitBaseFuncBody` already accepts `(domainConstraint, monoFactorLogs, imagExps)`
and emits the `SplitRealImag` phase block; the convergent emitter already passes
them. Change the **divergent** emitter calls in `emitIntegrandDefinitions`
(G0/G1/boundary-base/boundary-log/term-base/term-log; and the hand-written
remainder) to pass `Lookup[piece,"MonoFactorLog",None]` and the sector's
`Im(B)` (`Lookup[ibpSD/dd,"ImagPolyExponents",None]`). With `Im(B)=None`
(real path) the phase block is skipped ⇒ **byte-identical** emitted C++ (#25).

### 4.4 Attach `ImagPolyExponents` to divergent lifted sectors + route them
- The lifted driver already attaches `ImagPolyExponents = Im(B)` to convergent
  sectors in `cxSplit` mode; do the same for the **divergent** lifted sectors and
  carry it into `IBPProcessSector`/`ProcessDivergentSector` output.
- **Remove the `splitliftdiv` refusal** once C1 is green; route `cxSplit` lifted
  divergent calls to the (now phase-aware) IBP path.

### 4.5 Subtraction route (Phase C0 — already works)
No code change beyond a regression test: `LaurentFromSubtraction[liftedSpec, fan,
kp, "LiftData"->ld, "ComplexExponentMode"->"SplitRealImag", "EpsilonValues"->…]`
already returns the complex Laurent (§thesis). Add it as a pinned cross-check so
the capability cannot silently regress.

### 4.6 Direct-mode complex lifting (Phase C3 — optional, genuine subtlety)
Direct mode keeps the **complex** flattening (the spiral contour the unlifted
Direct path already uses; see `MANUAL/complex_flatten_spiral.pdf`). The only
blocker is the domain indicator. The domain constraint is geometrically a cut on
the **modulus** `|y_p*| ≤ 1`; with `y_j = (y'_j)^{1/ãtilde_j}`,
`log|y_j| = Re(1/ãtilde_j) log y'_j = (Re(ãtilde_j)/|ãtilde_j|^2) log y'_j`, so a
**real** indicator `ic_j = m_other,j · Re(ãtilde_j)/|ãtilde_j|^2` is available.
**Subtlety (must be sandboxed before trusting):** does the modulus cut on the
deformed (spiral) contour equal the true restriction of the original real
δ-resolution? This is the analog of the Case-B coupling question and is *not*
obvious. Phase C3 sandboxes it against NIntegrate; if it fails to reproduce the
reference, **keep `liftcomplex` for the Direct+lift case** and document
`SplitRealImag` as the canonical complex-lifting mode (C0/C1 already cover it).

---

## 5. Exactness discipline (plan.md §1.4 / §6.7)

- `Re(ãtilde)`, `Im(ãtilde)`, `D_{j,i}`, `Const_j` are built with `Re[…]`/`Im[…]`
  and exact arithmetic — **no `N[…]`** before the `MmaToC` boundary. The
  `SplitRealImag` codegen already numericizes `Im(B)`/`MonoFactorLog` only inside
  `mmaToCInternal`.
- Guard #43: extend the corpus with a **complex-B lifted+divergent** exact spec
  and assert `FreeQ[symbolicOutput, _Real]` for the boundary/term `MonoFactorLog`
  and `Re(ãtilde)`-derived fields before `MmaToC`.

---

## 6. Phased plan (with gates; ≥2 oracles each, no gate judged by its producer)

**Phase C0 — Pin the (already-working) Subtraction route.** Add a cross-check that
`LaurentFromSubtraction + SplitRealImag` reproduces the complex Laurent of a
lifted+divergent complex-B toy. *Gate:* recon `F(ε*)` within 1% of NIntegrate
(complex), pole non-zero and complex; matches the unlifted Direct reference of the
same integral with the extreme coefficient un-lifted.

**Phase C1 — IBP route via per-piece `MonoFactorLog` (the core).** Implement
§4.1–4.4. *Gates:*
- **byte-identity:** real-B and unlifted goldens unchanged (#25 byte-identical);
  every existing convergent `SplitRealImag` golden unchanged (the divergent
  emitter only adds a phase when `Im(B) ≠ None`).
- **#cx-A (IBP vs Subtraction, complex):** the (pole, finite) of the toy agree
  between the IBP route and the C0 Subtraction route (complex), `< 0.5%`.
- **#cx-B (vs NIntegrate at ε\*):** reconstructed `I(ε*)` (complex) vs direct
  NIntegrate of the original at `ε* ∈ {0.02, 0.01}` within tolerance.
- **#cx-C (real-B reduction):** with `Im(B) = 0` the IBP route reproduces the
  planAXpDIV real-B numbers exactly (the phase block is absent).

**Phase C2 — docs + cross-checks + limitation update.** Remove `splitliftdiv`;
update SUMMARY L2 / plan.md N3 / the manual (limitation item 2(iii)) to "complex-B
lift+divergence supported via `SplitRealImag` (both routes)". Extend `cc_45` with a
complex-B row at one dimension and a new `cc_46` (complex lift+div, IBP vs Sub vs
NIntegrate) + the #43 corpus entry.

**Phase C3 (optional, stretch) — Direct-mode complex lifting.** Sandbox the
modulus-based real domain indicator (§4.6) against NIntegrate on a lifted complex
toy with a *non-trivial* domain constraint. *Gate:* the Direct-mode lifted
decomposition reproduces NIntegrate (convergent) and the C0/C1 complex Laurent
(divergent). If the contour/modulus subtlety does **not** reproduce the reference,
abandon C3 and keep `liftcomplex` for Direct+lift (clean `$Failed`), with
`SplitRealImag` documented as the required complex-lifting mode.

---

## 7. Cross-checks (≥2 oracles each)

- **`cc_46` (new): complex-B lift × divergence.** Toy
  `∫ x1^{-1+ε}(1 + x1 x2 + 10^6 x1^2 + x2^2)^{B}`, `B = −2 + 0.3 i`, lift `{2,0},k=3`.
  Oracles: IBP route (C1), Subtraction route (C0), NIntegrate-at-ε\* of the
  original (complex). PASS = complex pole agreement `< 0.05` across routes and
  recon within 2% of NIntegrate. Add a 3D variant (e.g. the `cc_45` 3D-c
  numerator×denominator integrand with a complex denominator exponent) to exercise
  variable dimension.
- **Extend `cc_45`** with one complex-B row (Subtraction route) so the
  multi-dimensional battery covers complex exponents.
- **`cc_21b` regression:** unchanged (lifted convergent complex SplitRealImag) —
  confirms C1 did not perturb the convergent phase path.
- **#43 corpus:** complex-B lifted+divergent exact spec, `FreeQ[_, _Real]` on the
  boundary/term `MonoFactorLog` + `Re(ãtilde)` fields.
- **CUBA (#5) where available:** Cuhre on the direct complex integrand at fixed ε.

---

## 8. New messages / limitation updates

```
(* C1 removes: *)
TropicalEval::splitliftdiv   (* SplitRealImag x lift x divergence now supported *)

(* C3, if Direct-mode complex lifting is NOT delivered, keep and refine: *)
TropicalEval::liftcomplex =
  "ProcessSectorLifted: cone `1` — complex atilde has no real domain indicator in
   Direct mode.  Use ComplexExponentMode -> \"SplitRealImag\" for complex-exponent
   lifting (supported, incl. with a 1/eps pole)."
```

Update **SUMMARY.txt L2(iii)**, **plan.md N3 / §9 D3**, and **MANUAL** (limitation
item 2(iii), §5.7.4) from "complex-B lifting inherits `liftcomplex`" to:
*"Complex polynomial exponents with lifting (and with a `1/ε` pole) are supported
via `ComplexExponentMode -> "SplitRealImag"` on both the Subtraction and IBP
routes (cc_46); `Direct` mode with complex-B lifting remains refused
[unless C3 lands]."*

---

## 9. Scope boundaries

**Delivers:** lifting + a single `1/ε` pole + **complex** polynomial exponents,
via `SplitRealImag`, on **both** the Subtraction route (C0, already working) and
the IBP route (C1); full exactness discipline; byte-identical real-B and unlifted
goldens; variable dimension and non-trivial polynomials (cc_46 + cc_45 complex
row).

**Does not deliver (explicit future work, refused with a clean `$Failed`):**
Direct-mode complex-B lifting unless Phase C3's modulus-domain subtlety checks out;
Case B (domain–pole coupling, planAXpDIV §3); `1/ε^{d≥2}` in lifted sectors
(inherits L1/L8); multiple independent auxiliary variables. The single-pole-per-
sector limit and the `HasConstantTerm=False` route guidance (use the Subtraction
route) from planAXpDIV continue to apply.
