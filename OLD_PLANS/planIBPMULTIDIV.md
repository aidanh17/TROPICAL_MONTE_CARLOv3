# planIBPMULTIDIV — IBP on a multi-soft sector: **one real pole + N off-axis directions**

> **STATUS (implemented 2026-06-30).** P1–P4 done for the **UNLIFTED** case:
> `ibpDivClass` (§3.1), off-axis-aware `nestedIBP`/backstop (§3.2), the full
> `2^(N+1)`-corner iterated IBP (`IBPBuildCorners`/`ibpProcessLeaf`/
> `IBPProcessSectorMultiDiv`, §3.3) with the `[1/(c_kε)]·∏_j[1/(iθ_j)]` driver
> assembly, the off-axis-aware post-IBP `Catch`/`Throw` verify guard (§3.4), the
> `emitDomainIndicatorCpp` `.real()` complex-typing fix (§3.5), and the
> `AnalyticPole` `ck=0` guard (§3.6). Gates: **cc_51** (N=1 + an N=1 `c_m≠0`
> Part C exercising the K₁ correction) and **cc_52** (§3.5) pass vs the closed-form
> Dirichlet oracle; an N=2 reproduction passes; cc_45/46/47/48 + phase2/phase3
> goldens unchanged (byte-identity #25 preserved). **NOT delivered:** the
> *lifted* co-located one-pole+off-axis case — `ProcessSectorLifted`'s pivot
> admittance counts realified `Re≤0` directions (not poles) and refuses it
> (`liftnopivot`), a clean `$Failed` (cc_52 Part B). The off-axis-aware lifted
> pivot admittance (needed by the fourpt consumer) is the remaining future work;
> it was beyond this plan's §3 design, which assumed the lifted sector reaches
> `IBPProcessSector`.
>
> Spec for the next feature increment. **Nothing here is implemented yet.** Written
> 2026-06-30 in response to `CROSSCHECK/IBPCX_LIFTED_DIVERGENT_FEATURE.md` (the 4-point
> cross-check's request), after reading that document, `planIBPCX.md` (its §3.5/§8 name
> this exact item as future work), the M2b engine code, and **running a minimal
> oracle-backed reproduction** (§2 below) that pins down the actual failing guard.
>
> Engine references are line numbers in `tropical_eval.wl` (v3-build @ `3955950`,
> md5 `fa397a9e`, 8154 lines). They were re-verified against that commit; several line
> numbers in the incoming feature request are off (see §1.3) — use the ones here.

---

## 0. TL;DR

`planIBPCX` made the **single** divergent direction off-axis-aware: a lone direction with
`Re(α₀)=0, Im(α₀)=θ≠0` is no longer refused — it is finite (`1/(iθ)`), and IBP integrates
it (cc_47 Part C). What it did **not** do is make the **divergence *count*** off-axis-aware.
So a sector with **one genuine pole** (`Re=0 ∧ θ=0`) **plus one or more off-axis directions**
(`Re=0 ∧ θ≠0`) in the *same* cone is counted as `nDiv ≥ 2` and refused
(`TropicalEval::nestedIBP`) — even though it has exactly **one** true `1/ε` pole and the
off-axis directions are finite. This blocks lifting the *divergent* presectors of the 4-point
`J⁰⁰` figure.

**Three things to do, one of which is the bulk of the work:**

1. **Classify directions by `(Re, θ)`, not `Re` alone, everywhere directions are counted.**
   `genuine pole ⟺ Re(α₀)=0 ∧ θ=0`; `off-axis ⟺ Re(α₀)=0 ∧ θ≠0` (finite). Off-axis
   directions are excluded from the pole count. `ibpImagPole` (`:4950`) already computes
   `θ` for any direction — reuse it. (Small, mechanical.)

2. **Generalize the IBP boundary construction + driver Laurent assembly from a *single*
   divergent direction to *one pole + N off-axis*.** This is the real work: the off-axis
   directions are not absolutely integrable (`|y^{iθ-1}| = 1/y`, planIBPCX §1.2) so they
   each need an IBP step too, contributing a finite `1/(iθ_j)` boundary+bulk; the result is
   an **iterated IBP** over the divergent set with prefactor `[1/(c_k ε)]·∏_j[1/(iθ_j)]` and
   a multi-corner boundary/bulk sum. The current `IBPProcessSector` Step 2 builds **one**
   boundary at `y_k=1` for a single `k`; the driver uses a single `ck`/`ImagPole`.

3. **Fix the `emitDomainIndicatorCpp` complex-typing bug (`:2364–2379`).** It emits
   `double log_ypstar = (… cx(…) …)` when the lifted domain indicator is complex-typed,
   which g++ rejects. **Contrary to the feature request, this *is* on the IBP critical
   path** — the IBP codegen calls `emitDomainIndicatorCpp` at `:2786/:2811/:2845/:2862` — so
   it must be fixed for lifted-`SplitRealImag` divergent sectors. The indicator is provably
   real (the `SplitRealImag` domain map is built from `Re` parts); emit `.real()`.

Scope is the engine (`tropical_eval.wl`) and its `TEST/cc_*`/`phase2_selfgate.wl` gates.
This **delivers** "one pole + N off-axis" on the IBP route. It does **not** deliver genuine
nested poles (`>1` real pole, `1/ε^{d≥2}`), which stay refused (§7).

---

## 1. The misdiagnosis, corrected by reproduction

### 1.1 What the feature request got right

- It is a documented **feature gap**, not a defect: the engine cleanly `$Failed`s on a case
  `planIBPCX §8` explicitly lists as future work ("one real pole + one off-axis direction",
  `nestedIBP` counts both). Invariant #3 (known-limitation ⇒ clean `$Failed`) holds today.
- The physics: `J⁰⁰` has **one** simple UV pole, `c_k = 2/7`, **one divergent variable per
  divergent cone per the realified bookkeeping** — and the principal series `μ = iν̃` makes
  exponents complex, so off-axis directions appear (`Im(B)·D` / `Im(A)·M` routed onto a soft
  ray). IBP is the unique route (`planIBPCX §1.2`): Subtraction / pinned-ε share the
  vanishing-exponent blind spot.
- The lift is needed for the `k₁₂ ~ 10⁴` hierarchy at deep squeeze.

### 1.2 What the reproduction shows (the corrected diagnosis)

I built the minimal one-pole-plus-one-off-axis sector with a **closed-form oracle**
(`scratchpad/repro_offaxis.wl`):

```
∫₀^∞∫₀^∞ x₁^{2ε−1} x₂^{iθ−1} (1+x₁+x₂)^{−B} dx₁ dx₂
   = Γ(2ε) Γ(iθ) Γ(B−2ε−iθ) / Γ(B)            (Dirichlet)
```
`x₁` is a genuine real-soft pole (`α₁ = 2ε`, `Re=0 ∧ Im=0`); `x₂` is off-axis
(`α₂ = iθ`, `Re=0 ∧ Im=θ≠0`). `Γ(2ε)` gives the **single** `1/ε` pole; `Γ(iθ)` is a finite
complex constant — exactly **one pole + one finite off-axis direction**.

Running `EvaluateTropicalMC[…, "Method"->"IBP", …]` (both `Direct` and `SplitRealImag`),
`θ=1/2, B=3`:

```
1 convergent, 3 divergent sectors
TropicalEval::nestedIBP: Sector 1: 2 divergent variables. … Refusing this sector ($Failed)
IBPReduceSector: sector 2, divergent variables: y_{2}, effective exponents at eps=0: {3, I/2}    <- single off-axis, IBP runs
IBPReduceSector: sector 3, divergent variables: y_{2}, effective exponents at eps=0: {3-I/2, 0}   <- single real pole, IBP runs
ERROR: IBP processing failed for 1 of 3 divergent sector(s) …
```

**Sector 1 is the corner cone where both rays are soft**: effective exponents `{0, I/2}` →
both have `Re ≤ 0` → `nDiv = 2` → **refused by the `nestedIBP` front guard (`:5150–5157`),
before any IBP step runs.** Sectors 2 and 3 are single-divergence and proceed.

**Conclusions (these redirect the plan):**

- The failing guard is the **pre-IBP `nestedIBP` count** (`:5150–5157`), *not* a post-IBP
  "verify all `α₀>0`" trip. It is **case A** (off-axis miscount), **not case B** (no genuine
  higher-order residual exists — `Γ(iθ)` is finite, the single pole is `Γ(2ε)`).
- The off-axis direction is miscounted as a pole purely because the count uses `Re(α₀)≤0`
  and ignores `θ`. This holds in **both `Direct` and `SplitRealImag`** modes (the count is
  `Re`-only either way), so the fix is mode-independent.
- The discriminating diagnostic the feature request §4 asked for is hereby **run**: it is
  case A. (The feature request authored §4 honestly — "I have not confirmed it"; this is the
  confirmation.)

### 1.3 Line-number / mechanism corrections to the feature request

The feature request pins the same commit but its §3 attributes the failure to a post-IBP
"verify all `α₀>0`" guard "at `tropical_eval.wl:5148–5157`, error at `:5152`." That is wrong
on both counts, and it matters because it pointed at the wrong fix:

- `:5148–5157` **is the `nestedIBP` `nDiv>1` guard** (the one that actually fires), not a
  verify-`α₀` guard.
- The real post-IBP "verify all `α₀>0`" guard (the `ERROR: … alpha0[i] = … <= 0` string) is
  at **`:5310–5319`** (ERROR print `:5314`, `Return[$Failed]` `:5316`). It is reachable only
  *after* `nestedIBP` passes; in this feature it **never fires** (after full IBP reduction
  every divergent direction is raised, so all `α₀>0`). The feature request's "runs its single
  IBP step (1→31 terms) then `$Failed`" conflates the `IBPReduceSector` Prints of a
  *single-divergence* sector (like sectors 2/3 above) with the `nestedIBP` refusal of a
  *different* (corner) sector.

So we implement the **case-A classification fix at the count/assembly layer**, and there is
**no iterated-IBP-for-a-genuine-residual** (case B) requirement.

---

## 2. What is already correct (so the change targets the gap, not a rewrite)

1. **Per-direction `θ` already exists.** `ibpImagPole[sectorData, k, eps]` (`:4950–4963`)
   computes `θ` for **any** direction index `k`, unifying all three sources (Direct complex
   `Im[a0[[k]]]`; lifted `Im(B)·D` via `LiftedMonoFactor`; lifted `Im(A)` via
   `LiftedMonoPhase`). It is built from exact `Im[…]` (no `N[…]`), invariant #1 safe.

2. **The single-off-axis driver branch exists.** `evaluateTropicalIBPDriver` (`:6221–6242`)
   already branches: `PossibleZeroQ[ImagPole]` → real `1/(c_k ε)` pole (byte/numerically
   identical, #25); else off-axis finite `(bndBase−S0)/(iθ₀)`. The N-off-axis assembly
   **generalizes** this branch (multiply in `∏_j 1/(iθ_j)` and sum the off-axis corners).

3. **The bounded off-axis phases are already emitted.** Per-piece `MonoFactorLog` /
   `MonomialPhaseLog` divide the fixed `Im` numerators by the **raised** `α₀` (`:5341–5354`,
   boundary `:5266–5282`), so the C++ integrands are already numerically bounded for any
   off-axis direction. No integrand/codegen change for the *phases*.

4. **`IBPReduceSector` already iterates over all divergent directions.** Its `allDivVars`
   (`:4988–4996`) collects every `Re(α₀)≤0` direction and the outer `Do` (`:5019–5065`) IBP-
   reduces **each** — so the symbolic reduction for multiple soft directions is *already
   present*; it is the **count guards** (`:5150`, `:5174`) and the **single-direction
   boundary/driver** that block and mis-assemble it.

---

## 3. Design

### 3.1 Direction classification (the small, shared primitive)

Add one helper next to `ibpImagPole` (`:4963`):

```
ibpDivClass[sectorData, i, eps] :=
  Module[{a0i = sectorData["NewExponents"][[i]] /. eps -> 0,
          th  = ibpImagPole[sectorData, i, eps]},
    Which[
      ! (TrueQ[Re[a0i] <= 0] || (NumericQ[a0i] && Re[a0i] <= 0)), "Convergent",
      TrueQ[PossibleZeroQ[th]],                                    "Pole",
      True,                                                        "OffAxis"]]
```

`"Pole"` ⟺ `Re(α₀)=0 ∧ θ=0` (a genuine `1/ε`). `"OffAxis"` ⟺ `Re(α₀)≤0 ∧ θ≠0` (finite).
Decided **symbolically** (`PossibleZeroQ` before any `kinRules`), so it is exact and makes no
float decision in the decomposition (invariant #1; mirrors `planIBPCX §3.4`).

> **Edge case — power-divergent off-axis (`Re(α₀)<0 ∧ θ≠0`).** A single IBP step raises the
> exponent by `e≥1`; if `Re(α₀)<0` strictly it may not reach `Re>0` in one step. `J⁰⁰`'s
> off-axis directions are logarithmic (`Re=0`). Guard: classify `Re(α₀)<0 ∧ θ≠0` directions,
> and if a *raised* off-axis term still has `Re(α₀)≤0`, refuse cleanly (a new message, or
> reuse the post-IBP guard at `:5310` made off-axis-aware — see §3.4). Do not silently emit a
> divergent integrand.

### 3.2 Off-axis-aware divergence count (the unblock)

- **Front `nestedIBP` guard (`:5150–5157`).** Replace the `Re≤0` count with a **genuine-pole**
  count: `nPole = Count[Range[n], i /; ibpDivClass[sectorData,i,eps]==="Pole"]`. Refuse
  (`nestedIBP`) only when `nPole > 1`. Off-axis directions no longer inflate the count.
- **Backstop (`:5174–5177`).** Today it refuses when `Length[divVars] != 1`. With one pole +
  N off-axis, `IBPReduceSector` resolves `1+N` directions, so this must check the **pole**
  count (`==1`), not the total resolved count. Track which resolved direction is the pole.
- **`IBPReduceSector` (`:4988–4996`, `:5019–5065`).** Keep reducing **all** `Re≤0` directions
  (poles *and* off-axis both need raising), but tag each: return `"PoleVariable" -> k_pole`
  and `"OffAxisVariables" -> {m_1,…}` (classified via §3.1) alongside the existing
  `DivergentVariables`. The `badck` guard (`:5029–5033`) is already correct (fires only when
  `ck=0 ∧ θ=0`).

### 3.3 Boundary construction + driver assembly (one pole + N off-axis) — **the bulk**

The sector integral with divergent set `{k (pole), m₁…m_N (off-axis)}` is, by iterated IBP:

```
I(ε) = [1/(c_k ε)] · ∏_j [1/(iθ_j)] · Σ_{corners} (±1) · L_corner(ε)
```
where each *corner* chooses, per divergent direction, **boundary** (`y_d = 1`, drop the
coordinate, sign +) or **bulk** (raise `α_d`, integrate, sign −). The pole `1/(c_k ε)` is the
**only** ε-singular factor; `∏_j 1/(iθ_j)` is a finite constant. So the Laurent is still just
orders `−1` (pole) and `0` (finite) — **no `1/ε²`** (confirmed in §1.2). Two implementation
routes; **Route B recommended**:

- **Route A (full multi-corner):** generalize `IBPProcessSector` Step 2 to emit the `2^{N+1}`
  corner integrands and tag each with its prefactor sign and the boundary-fixed set; the
  driver sums them with `[1/(c_k ε)]·∏_j 1/(iθ_j)`. Most faithful, most code.

- **Route B (peel off-axis, then reuse the single-pole machinery) — recommended.** IBP-reduce
  the **off-axis** directions in a pre-pass that expands the sector into a list of
  sub-sectors, each tagged with a finite scalar `∏_j (±1/(iθ_j))` and either dropping `y_{m_j}`
  (boundary, dimension `−1`) or raising it (bulk). Each resulting sub-sector then has **exactly
  one** divergent direction (the pole `k`) and flows through the **existing**
  `IBPProcessSector`/driver single-pole path unchanged, with the off-axis scalar carried as a
  multiplicative `Coeff`-style constant into the driver assembly (`:6206–6219`, `:6221–6233`).
  This reuses the verified single-pole boundary/Laurent code and confines new logic to the
  pre-pass + a scalar multiply. The per-sub-sector dimensions differ (boundary corners drop
  coordinates), which the per-function dim handling (`integrandDims`) already supports.

Either route: the **pole** sub-pieces use the existing real-pole formulas (`:6225/:6228`),
multiplied by the carried off-axis constant; the driver's existing off-axis branch
(`:6240–6241`) is subsumed (single off-axis = the `N=1`, no-pole special case).

> The off-axis prefactor is `1/(iθ_j)` with `θ_j` the **exact symbolic** `ibpImagPole` value
> (kinematic-dependent allowed). Keep it symbolic through codegen; substitute `kinRules` and
> apply the mixed-grid `θ₀=0` fallback (`:6221–6222`) per kp in the driver, exactly as the
> single-off-axis branch already does.

### 3.4 Post-IBP verify guard, off-axis-aware (defensive)

The `:5310–5319` "verify all `α₀>0`" guard fires per IBP term. After §3.2/§3.3 fully raise
every divergent direction it should pass. Make it **off-axis-aware** as a backstop: pass a
`Re(α₀)=0` direction iff `ibpImagPole` for that term/direction is non-zero (finite off-axis);
fail only on `Re(α₀)≤0 ∧ θ=0` (a genuine unreduced residual, which would be a real bug or the
§3.1 power-divergent edge case). Replace the `Print["ERROR…"]; Return[$Failed]` with a proper
`Message[TropicalEval::…]` so the refusal is a clean `$Failed` (the current `Return` inside the
`Table`/`Module` does not even propagate cleanly out of `IBPProcessSector` — it drops a
`$Failed` into one `IBPTerms` cell — so this also fixes a latent failure-propagation bug).

### 3.5 The `emitDomainIndicatorCpp` complex-typing bug (`:2364–2379`)

`double log_ypstar = (logZ0str − (sumStr)) * (1.0/mpStr);` where `logZ0str`, `mpStr`, and the
`icList` entries come from `mmaToCInternal[N[…]]`. If any is `Complex[x, 0.]` (the machine-zero
imaginary that `realifyMonoA` warns about, `:2341–2347`), `mmaToCInternal` emits `cx(x,0.0)`
and the `double =` assignment fails to compile (`cannot convert std::complex<double> to
double`). The feature request observed this live on the lifted `SplitRealImag` fourpt indicator.

- **It is on the IBP critical path**, not "moot": the IBP boundary/term emitters call
  `emitDomainIndicatorCpp` at `:2786, :2811, :2845, :2862`. (cc_47 Part C happens to land a
  real-typed indicator and compiles — the bug is geometry-dependent, which is why it slipped.)
- **Fix:** the domain indicator is provably real in `SplitRealImag` mode (the map is built
  from `Re(A)/Re(B)`; `LogZ0`, `MP`, `IndicatorCoeffs` are real). Emit a real value: either
  wrap the assembled RHS with C++ `.real()` (so `double log_ypstar = (…).real();`), **or**
  `Re[…]` the three WL quantities before `mmaToCInternal` (`Re[N[dc["LogZ0"]]]`, etc.).
  Prefer the WL-side `Re[]` — it keeps the emitted token a plain `double` literal, smaller
  byte delta, and provably correct. Either way it changes bytes **only** when the indicator is
  complex-typed → real/`Direct` goldens byte-identical (#25).

### 3.6 Secondary robustness: `AnalyticPole -> 1/ck` (`:5433`)

When the lone (or pole) divergent direction is pure-imaginary off-axis with `ck=0` (e.g.
`α = iθ`, no `ε`), `1/ck` evaluates to `ComplexInfinity` and prints `Power::infy` (seen in the
§1.2 run, sectors 2/3). It is only a diagnostic field, but guard it:
`"AnalyticPole" -> If[TrueQ[ck==0], Indeterminate, 1/ck]`. Cosmetic, prevents a confusing
warning and a `ComplexInfinity` leaking into a returned association.

---

## 4. Exactness / byte-identity checklist

1. **Exactness (#1).** Classification and `θ_j` are symbolic (`PossibleZeroQ`, exact `Im[…]`);
   no `N[…]` upstream of `MmaToC`; the decomposition makes no float decision. The pole/finite
   branch lives in numeric driver assembly (as `planIBPCX §3.4` established).
2. **Byte-identity (#25).** Every new branch activates only when a sector actually has an
   off-axis divergent direction (`∃ i: class=="OffAxis"`), or when the domain indicator is
   complex-typed. A sector with `≤1` `Re≤0` direction and `θ=0` follows the **current** code
   path exactly → IBP codegen goldens (`TEST/baselines/codegen_goldens/ibpdiv_*.cpp`) and
   cc_45/46/47/48 numbers bit-identical.
3. **Known-limitation ⇒ clean `$Failed` (#3).** `>1` genuine pole, geometric Case B
   (`liftdivdomain`), eps-dependent `θ` (`splitdivmono`), and power-divergent off-axis that
   one IBP step cannot tame all stay refused with a `Message` + `$Failed`.
4. **≥2 oracles per PASS (#4).** §6.
5. **Never edit `OLD_CODE/` (#5).**

---

## 5. Where this stays refused (unchanged `$Failed`)

- **More than one genuine pole** (`>1` direction with `Re=0 ∧ θ=0`): real `1/ε^{d≥2}`,
  `nestedIBP`. Genuinely not implemented.
- **Geometric Case B** (`liftdivdomain`, `:215`): divergent variable couples to the lifted
  domain face. Unchanged.
- **eps-dependent `θ`** (`splitdivmono`, `:216`): `Im` that scales like `1/ε` reintroduces
  fast oscillation. Unchanged.
- **Inline symbolic-ε Subtraction** for lifted+complex+divergent (`splitliftdiv`, `:217`):
  IBP / `LaurentFromSubtraction` only; unchanged.

---

## 6. Validation gates

1. **Regression / byte-identity (θ=0).** cc_45/46/48 and cc_47 Parts A/B unchanged; IBP
   codegen goldens byte-identical; a θ=0 spot sector reproduces current `(pole, finite)`
   bit-for-bit.
2. **`cc_47` Part C unchanged.** Single off-axis (no co-located pole) still passes — it is the
   `N=1, no-pole` special case of the new assembly.
3. **New `cc_51` — unlifted one-pole + one-off-axis vs the Dirichlet oracle (§1.2).**
   `∫ x₁^{2ε−1} x₂^{iθ−1}(1+x₁+x₂)^{−B}`: assert IBP returns a result (no `$Failed`), with
   `pole = ½·Γ(iθ)Γ(B−iθ)/Γ(B)` (genuinely complex) and `finite` matching the `Series` of
   `Γ(2ε)Γ(iθ)Γ(B−2ε−iθ)/Γ(B)` (closed form — shares no code with the engine). Oracle #2:
   `Direct == SplitRealImag` on the IBP path. Confirm the sector realizes the structure by
   asserting `nPole==1 ∧ ∃ OffAxis` on the corner sector (a fixture that silently avoids the
   co-located case proves nothing — mirror cc_47 Part C's `ImagPole` self-check).
4. **New `cc_52` — lifted one-pole + one-off-axis.** Same integrand with a `LiftCoefficients`
   on an extreme coefficient (deep-squeeze analog of the fourpt hierarchy); IBP vs NIntegrate
   of the complex original (+ reseed). Exercises §3.5 (the domain-indicator `.real()` fix) and
   §3.3 with a lifted `PrefactorBase`/`DLogPrefactor`. Confirm the lifted indicator path
   compiles (the bug-fix gate).
5. **Fourpt acceptance (consumer-side, `CROSSCHECK/fourpt`).** Lifted divergent presector 1,
   `μ=I`, deep squeeze: IBP no longer `$Failed`; `(pole, finite)` **matches the unlifted IBP**
   at the same point (`probe6`: pole `≈ −5.4e-22+2.5e-21i`, finite `≈ −8.0e-19` at
   `r1=10^-3.5, r2=0.9`) — the lift is an exact identity so the pole is unchanged. Sum over the
   8 `++` presectors (lifted) reproduces the unlifted total within MC error at `r1=r2=0.9`.
6. Add cc_51/cc_52 (and the unchanged-goldens assertion) as hard gates where appropriate; keep
   every numeric PASS on ≥2 independent oracles (#4).

---

## 7. Phases

- **P0 — Diagnostic confirmation (mostly done).** §1.2 reproduction confirms case A
  (nestedIBP miscount), one pole + one off-axis, oracle-backed. Remaining: get the
  consumer's full-complex dump of fourpt `Fourptpp` presector-1 / sector-22 to confirm the
  fourpt cone matches the reproduction's structure (`nPole=1`, `∃ OffAxis`) and is not a
  separate power-divergent case (§3.1). Build cc_51 as the in-repo P0 artifact.
- **P1 — Classification + count unblock (§3.1, §3.2, §3.6).** `ibpDivClass`; off-axis-aware
  `nestedIBP`/backstop; tag pole vs off-axis in `IBPReduceSector`; guard `AnalyticPole`.
  **Gate:** the corner sector is no longer refused up front; reduction tags are correct.
- **P2 — One pole + N off-axis assembly (§3.3, §3.4) — the bulk.** Route B (peel off-axis →
  reuse single-pole machinery), or Route A. Off-axis-aware post-IBP guard with a clean
  `Message`. **Gate:** cc_51 numeric PASS vs the Dirichlet oracle (pole + finite, Re and Im);
  cc_47 Part C unchanged; θ=0 spot sector bit-identical; goldens byte-identical.
- **P3 — Lift + domain-indicator fix (§3.5).** `.real()`/`Re[]` the domain indicator; verify
  the lifted IBP path compiles. **Gate:** cc_52 numeric PASS; lifted-`Direct` goldens
  byte-identical.
- **P4 — Docs.** `SUMMARY.txt` "Known limitations" + `planIBPCX §8` (promote the "one pole +
  N off-axis" item to implemented; keep `>1` genuine pole refused); MANUAL §"Error Messages"
  and the IBP section; `why_not_IBP.md`. Note the fourpt acceptance (§6.5) as the motivating
  case.

---

## 8. Open questions / risks

1. **Fourpt cone structure (P0).** The reproduction is `nPole=1 ∧ one OffAxis`. The feature
   request's §3 quoted a single-zero realified vector (`{10,4/13,3/4,0,1/2}`) with an
   `nDiv=1`/post-IBP narrative that §1.3 argues is a misread of the debug log. Confirm the
   *actual* failing fourpt cone is the co-located pole+off-axis (`nDiv≥2 → nestedIBP`) and not
   a power-divergent off-axis (`Re<0`) before assuming P1/P2 alone suffice. The fix shape is
   robust either way (P1 unblocks; §3.1 guard catches power-divergent), but the gate fixture
   should match the real geometry.
2. **Route A vs B (P2).** Route B reuses verified code but assumes the off-axis pre-pass and
   the single-pole machinery compose cleanly across differing sub-sector dimensions. Validate
   on cc_51 (`N=1`) before trusting `N≥2`; if composition is awkward, fall back to Route A.
3. **`emitDomainIndicatorCpp` real-typing (P3).** Confirm `LogZ0/MP/IndicatorCoeffs` are
   provably real in lifted `SplitRealImag` (they should be — domain map from `Re` parts). If a
   genuinely non-real indicator can arise, `.real()` would be *wrong* and the case must be
   refused instead; check before applying the fix.
4. **NIntegrate oracle for cc_52.** The lifted off-axis original is oscillatory; reuse cc_47
   Part C's fit-`I(ε)=f+cε` approach (it converges there) and add a reseed envelope if shaky.
5. **`N` large.** Route A is `2^{N+1}` corners; Route B is `2^N` off-axis sub-sectors. `J⁰⁰`
   has small `N` (1, maybe 2). `log()` any cap; do not silently truncate corners.

## 9. Repro / reference pointers
- `SANDBOX/repro_offaxis_multidiv.wl` — the minimal one-pole+one-off-axis Dirichlet
  reproduction run in §1.2; basis for cc_51.
- `CROSSCHECK/IBPCX_LIFTED_DIVERGENT_FEATURE.md` — the incoming request (line numbers per §1.3
  superseded by this doc).
- `planIBPCX.md §1.2/§3.3/§3.4/§8`; `tropical_eval.wl` guards: `:5150` (nestedIBP), `:5174`
  (backstop), `:5310` (verify-α₀), `:4950` (`ibpImagPole`), `:6221–6242` (driver branch),
  `:2364–2379` (`emitDomainIndicatorCpp`), `:5433` (`AnalyticPole`).
- `TEST/cc_47.wl` Part C (single off-axis template), `TEST/cc_48.wl` (Dirichlet-oracle style).
