# planLIFTDIVZ0.md — Bug plan: dropped z0 factor in the lifted × divergent × IBP × SplitRealImag value assembly

**Author:** Claude Code (from `bug_summary` + `new_request`) · **Date:** 2026-07-01
**Engine under investigation:** branch `ibpcx-lifted-divergent-regulator-real` @ `6020ef9`
("IBP on lifted divergent sectors with the regulator inside the exponent"),
`md5(tropical_eval.wl) = 89ef1d2` — **this is exactly the version `bug_summary` instrumented.**
**Suspected defect class:** silent **wrong number** (violates invariants #3 and #4), *not* a `$Failed`.

---

## RESOLUTION (2026-07-01) — NO z0 BUG; `bug_summary` MISDIAGNOSED MC VARIANCE

Phases 1–2 were executed. **The lifted × divergent × IBP × `SplitRealImag` path is
CORRECT; there is nothing to fix.** `bug_summary`'s "dropped `z0^{~1.5}` factor" is the
**§4.1 trap it warned about**: it compared the (correct, low-variance) *lifted* result
against an *unlifted* reference that is **Monte-Carlo-variance-dominated noise** on the
`~z0` coefficient hierarchy — and never cross-checked that reference against an
independent oracle.

Evidence (all reproduced with repo primitives + the repo engine `6020ef9`, `md5 89ef1d26`):

1. **`TEST/cc_54.wl` (the Phase-1 decider, PASSES).** A repo-native fixture driving the
   *exact* implicated combination — lifted, divergent single-pole, `Im(B)≠0`, and (the
   axis `cc_53` lacks, §1.3) **`Im(A)≠0` routed onto the pivot so
   `LiftedMonoPhase["Const"] = −3·Log[z0]/4 ≠ 0`**. Lifted IBP `(pole, finite)` matches
   **both** the NIntegrate complex Laurent **and** the pinned-ε Subtraction route to
   ~0.3–0.6 % **in magnitude and phase**. A dropped `z0^{1.5}` would be a ~1500 % error.
   → The Im(A) z0-phase (the leading suspect) is threaded correctly.
2. **Unlifted path is sound *without* a hierarchy.** Same complex-A × complex-B divergent
   integrand with O(1) coefficients: unlifted IBP (Direct = SplitRealImag) matches
   NIntegrate to 0.3 %. So the unlifted algorithm is correct; only its *MC variance*
   degrades when a large coefficient is present.
3. **The actual fourpt presector-1 (`mu=I`, mild `k12=10`) confirms it.** Loaded the real
   spec from `CROSSCHECK/fourpt/.../mu_1.0000i/run_1/presector_integrands_fourpt.m`
   (5 vars, 6 polys, `A_i=-1/2+I`, complex `B`); `DetectExtremeCoefficients` at threshold
   35 gives **`z0 = 6.6667` — exactly `bug_summary`'s value**. Two *no-lift* MC draws at
   200k samples, seeds 42 vs 777, gave finite parts `+1.32e-7+7.83e-8i` vs
   `−1.56e-8−1.42e-7i` — **uncorrelated, opposite signs, spread larger than the values
   themselves.** The no-lift number is noise; `bug_summary`'s single draw was one sample
   from that distribution and its `|ratio|≈17.6≈z0^1.5` is (correct value)/(noise) that
   coincidentally landed near `z0^1.5`.
4. **Structural confirmation.** `probe10_samefn.wl` and `probe6_ibp_complex.wl`
   (the basis of the report) compare *only* lift-vs-no-lift; **neither uses an independent
   oracle.** The only independent-oracle probe, `probe9b_exact_mild.wl`, is **real-μ**, so
   it never exercised the complex/`SplitRealImag` path where the confusion arose.

**Deliverables:** `TEST/cc_54.wl` (new; passes; the permanent regression guard for this
path). **No change to `tropical_eval.wl`** (§5's fix is moot; touching the correct lifted
path would risk #25 / correctness). §5.3's `cc_53`/`cc_45–52`/`phase*_selfgate` gates are
untouched by construction.

**Follow-ups (not done here; separately tracked):** (a) the no-lift fallback of
`EvaluateTropicalMCLifted` silently returns variance-limited numbers on a hierarchy — a
reliability *warning* would stop a future lift-vs-no-lift comparison being misread as a
correctness bug; (b) the §6 Subtraction-route `.real()` typing item remains open.

The remainder of this document is the original (pre-resolution) investigation plan.

---

## 0. TL;DR

After the regulator-real fix (`6020ef9`) the lifted × divergent × IBP × `SplitRealImag`
path **runs** (it used to `$Failed` on a spurious `splitdivmono`). `bug_summary` reports that it
now returns a **wrong value**: with everything else identical and only the lift toggled on, the
finite part is off by a factor whose magnitude ≈ `z0^1.5` plus a residual phase — i.e. a complex
`z0` power is dropped somewhere in the value assembly, so `PrefactorBase`'s large negative `z0`
power (`z0^{-10}` in the reported case) is under-compensated. Because the lift is an **exact change
of variables**, lifted and unlifted must agree to Monte-Carlo precision; they don't.

This plan:

1. **§3 — builds a self-contained, repo-native example** (`TEST/cc_54.wl`) that decides the question
   with an *oracle-free identity* (lift-vs-no-lift, same driver, same point) plus a `NIntegrate`
   cross-check. All of `bug_summary`'s evidence lives in an external `CROSSCHECK/fourpt/` tree and a
   `scratchpad/bughunt/` sandbox that are **not in this repo**, so the bug must be reproduced here
   from repo primitives before anything is changed.
2. **§4 — determines whether the bug is real** in the clean repo and, if so, **localizes** the exact
   dropped `z0` factor (a `z0^{p+qi}` residual → the specific line).
3. **§5 — fixes it**, preserving byte-identity (#25), exactness (#1), and clean-`$Failed`-on-unsupported (#3).

**Do not start by "fixing" anything.** The first deliverable is the *example*; the bug is not yet
confirmed in the clean repo (see §1.3 and the honest-uncertainty note in §4.4).

---

## 1. Background and reconciliation of the two source documents

### 1.1 What path is implicated
Reproducing the 4-point loop seed `J^00` of arXiv:2211.03810, principal series `μ = i ν̃`, which makes
**both** the monomial exponents `A` **and** the polynomial exponents `B` complex, on the **divergent**
`++` presectors 1–2 (single simple UV pole, one divergent variable per cone), *with* the
auxiliary-variable uplift engaged (needed to tame the deep-squeeze `~10^4` coefficient hierarchy).
This is the `SplitRealImag × lifted × divergent × IBP` combination.

### 1.2 The two documents describe two *different* bugs, in sequence
- **`new_request` (2026-06-30):** the path used to `$Failed` on `TropicalEval::splitdivmono`
  ("theta_k is eps-DEPENDENT") because `Im[eps]` was left un-simplified (`theta_k = ν̃ − Im[eps]`
  *looked* regulator-dependent). Fix = assume the regulator is real where `theta_k`/`Im(B)`/`Im(A)`
  are formed. **This fix is already committed** (`6020ef9`: `imagPolyInfo`, `realifyPolyB`,
  `ibpImagPole` Refine; validated by `TEST/cc_53.wl`).
- **`bug_summary` (against `6020ef9`, the post-fix engine):** now that the path *runs*, it returns a
  **wrong number** — the `z0` bookkeeping in the value assembly is incomplete. **This is the bug this
  plan targets.** It is strictly downstream of `new_request`'s fix and more serious (a wrong number,
  not a refusal).

`bug_summary` also *disproves its own earlier* off-axis / MultiDiv hypothesis (point 2): instrumented
`ibpDivClass` on all 10 lifted divergent cones shows **exactly one `Pole` direction (theta = 0) and no
off-axis directions**. So the implicated code is the **standard single-pole lifted IBP path**
(`IBPProcessSector`, `tropical_eval.wl:5586`+), *not* `IBPProcessSectorMultiDiv`.

### 1.3 Why this is the *first* case to hit it, and why cc_53 did not catch it
`cc_53` Part B **is** a lifted × divergent × `SplitRealImag` IBP integral and it **passes at 8 %**
(it compares lifted-IBP to `NIntegrate` of the original and to the Subtraction route). The single
structural difference from the failing fourpt case:

| | monomial `A` | polynomial `B` | result |
|---|---|---|---|
| `cc_53` Part B | **real** (`{-1+eps, 0}`) | complex + regulator (`-2 + i/2 − eps`) | **PASS** |
| fourpt presector 1 (`μ=i`) | **complex** (`Im(A) ≠ 0`) | complex + regulator | **wrong number** |

So `cc_53` exercises the `Im(B)` phase (`MonoFactorLog` / `LiftedMonoFactor`) but leaves the
**`Im(A)` phase (`MonomialPhaseLog` / `LiftedMonoPhase`) completely dormant** (`imAexps === None`
⇒ both `None`, `tropical_eval.wl:1370`). The fourpt case is the first to drive a nonzero `Im(A)`
through the lifted **divergent** assembly. **This is the leading hypothesis and it dictates the
example design in §3: the reproducer must have `Im(A) ≠ 0`.** (It is a *hypothesis*, not a proven
fact — §4 tests it.)

### 1.4 The reported quantitative signature (the thing the example must reproduce)
Same driver `EvaluateTropicalMCLifted`, only the lift toggled (Threshold `1e9` = no-lift fallback vs
`35` = lift), divergent presector 1, `μ = i`, mild point `r1 = 0.1`:

```
no-lift: pole = -1.28e-8,   finite = -3.69e-8 + 1.40e-7 i
lift:    pole = -5.45e-13,  finite = -4.90e-9 + 6.59e-9 i
finite ratio no-lift/lift = 16.3 - 6.5 i,  |ratio| = 17.6 ≈ z0^1.5 (z0 = 6.667, z0^1.5 = 17.2)
```

A complex `z0` power (magnitude ≈ `z0^1.5`, plus a phase) is dropped. The pole is below MC noise at
this point; **the finite part is the clean signal.**

---

## 2. Map of the z0 machinery (so §4/§5 have precise coordinates)

`z0 = |C_primary|^{1/k}` is a **concrete algebraic number** (e.g. `10^{4/3}` in `cc_53`), *not* a
symbol. It is baked into coefficients/prefactors at `ProcessSectorLifted` time. (⚠️ **Red herring to
avoid:** an exploratory pass suggested "z0 is an unsubstituted C++ symbol at `:6588`". That is wrong
— if `z0` reached codegen as a free symbol, `cc_53` would fail to *compile*, but it passes. Do not
chase this.)

Where `z0` enters (all `tropical_eval.wl`):

| site | line | what carries `z0` |
|---|---|---|
| anchor computed | `941–942` | `z0 = Abs[Cprimary]^(1/kprimary)` (exact) |
| stored in `LiftData` | `1014` | `"z0" -> z0` |
| **cleared-poly coeffs** | **1093** | `cmono[[1]] * z0^(ep/mpLocal)` (per monomial; `ep` = pivot-column exponent) |
| **`PrefactorBase`** | **1298** | `(Abs[detM]/Abs[mp]) * z0^(ap/mp - 1)` (`ap = a[[pivotP]]`, the pivot effective exponent) |
| `MonoFactorLog` `Const` (Im(B) phase) | `1327` | `(dAug[[k,pivotP]]/mp) * logZ0` |
| `LiftedMonoFactor` (un-flattened Im(B)) | `1344–1355` | `Const = (dAug/mp) logZ0`, `DExp` |
| **`LiftedMonoPhase`** (un-flattened **Im(A)**) | **1370–1380** | `constA = (imEaug[[pivotP]]/mp) * logZ0`, `Num` — **`None` when `imAexps===None`** |
| divergent lifted `SectorData` return | `1389–1424` | carries all of the above |

How the **single-pole** path (`IBPProcessSector`, `5586`+) consumes them:
- `pfBase = Lookup[sectorData,"PrefactorBase",Abs[detM]]` (`5658`); boundary `bndPrefactor =
  (pfBase /. eps->0)/Times@@bndA0` (`5744`); per-term `flatPrefactor` (`5842`).
- boundary phase: `bndMonoFactorLog` threads `Const` (`5763`); `bndMonoPhaseLog` threads `Const` (`5772`).
- per-term phase: `termMonoFactorLog` threads `Const` (`5849`); `termMonoPhaseLog` threads `Const` (`5856`).
- `DLogPrefactor = D[Log[pfBase],eps]/.eps->0` (`5924`) — the finite-part correction for the
  eps-dependence of a lifted `PrefactorBase`.

**Observation that makes the bug non-obvious:** the phase `Const` (the `z0` log) *is* threaded in
both the boundary and the terms. So a naïve "the `Const` was forgotten" is not the whole story — the
dropped factor is more likely (a) an **`Im(A)`-specific** magnitude/phase piece that only appears when
`LiftedMonoPhase =!= None`, or (b) a mis-cancellation between the realified magnitude `z0^{Re(...)}`
and the phase `z0^{i·Im(...)}` when **both** `A` and `B` are complex, or (c) something specific to
the **`y_k = 1` boundary** where the divergent slot is dropped. §4 pins which.

The `SplitRealImag` realification that feeds this path (driver `EvaluateTropicalMC`, `4099`+):
`realifyPolyB[liftedSpec, imBList]` and `realifyMonoA[..., imAList]` subtract the eps-free imaginary
parts **before** `ProcessSectorLifted` (`4207–4229`, `4305–4311`, `4373–4392`), so the exponents that
build the `z0` magnitude powers are real and the imaginary parts must be re-supplied by
`MonoFactorLog`/`MonomialPhaseLog`. **The split of `z0^{complex exponent}` into a real-magnitude
power and an imaginary-phase constant is the exact place a factor can go missing.**

---

## 3. PHASE 1 — Build the example that decides it

**Deliverable:** `TEST/cc_54.wl`, a self-contained cross-check in the repo's `cc_*.wl` style
(mirror `cc_53.wl`'s scaffolding: `qr` quiet wrapper, `cc54Assert`, `INTERFILES/cc54/…` workdirs,
`CC54 PASS/FAIL` lines, `Quit[1]` on any FAIL so it can gate CI).

### 3.1 The decisive oracle — the lift-is-an-exact-identity test (probe10, repo-native)
The strongest discriminator needs **no external oracle**: run the *same* divergent integrand through
the *same* driver twice, toggling only the lift, and require agreement to MC precision.

```
ref  = EvaluateTropicalMCLifted[spec, pts, "Threshold" -> 10^12, (* no extreme coeff -> unlifted fallback *)
         "Method" -> Automatic, "ComplexExponentMode" -> "SplitRealImag",
         "Integrator" -> "MC", "NSamples" -> nS, "RunChecks" -> False, ...]
test = EvaluateTropicalMCLifted[spec, pts, "Threshold" -> tLift,  (* small enough to engage the lift *)
         (* all other options identical *) ...]
```

- `EvaluateTropicalMCLifted` (`tropical_eval.wl:4809`) with a huge `Threshold` finds no extreme
  coefficient and **falls back to plain `EvaluateTropicalMC` on the original spec** (`4827–4842`) —
  this is the trusted unlifted reference. With a small `Threshold` it lifts and routes through the
  suspect path. This is *exactly* `bug_summary`'s toggle (`1e9` vs `35`), reproduced from repo
  primitives.
- `Method -> Automatic` routes the divergent integral to IBP.
- **Assertion:** `|pole_test − pole_ref|` and `|(finite_test − finite_ref)/finite_ref|` are within
  MC tolerance (target the same order as `cc_53`, ~3–5 %). The bug predicts a **~1500 % finite-part
  error** (factor ≈ `z0^1.5`), which is unmissable.

### 3.2 Second, independent oracle (invariant #4)
Add a `NIntegrate` complex-Laurent fit of the **original** integrand (copy `cc_53`'s `niC`/`Fit`
block, `cc_53.wl:157–170`). Prediction under the bug: `ref ≈ NIntegrate` (both correct), `test`
off by `z0^{~1.5}·phase`. This distinguishes "the lift is wrong" from "the whole path is wrong."

### 3.3 Fixture requirements (what makes it different from cc_53)
The fixture must simultaneously:
1. **Be divergent** — a genuine single `1/eps` pole (mirror `cc_53`: `A_1 = -1 + eps` gives the
   `x_1 → 0` pole).
2. **Have an extreme coefficient so `z0 ≠ 1`** — e.g. `10^4 x_1^2` lifted with `k = 3`
   ⇒ `z0 = 10^{4/3} ≈ 21.5` (as in `cc_53`), or tune to `z0 ≈ 6.667` to match `bug_summary` exactly.
3. **Have `Im(B) ≠ 0`** (engages `MonoFactorLog`, as `cc_53` does) — regulator inside `B`
   (`B = b + i g − eps`) to also keep exercising the `6020ef9` regulator-real fix.
4. **★ Have `Im(A) ≠ 0`** — this is the new axis `cc_53` lacks (§1.3). Put the imaginary part on a
   **convergent** slot (e.g. `A_2 = i q`, `q ≠ 0`) so it activates `MonomialPhaseLog` /
   `LiftedMonoPhase` **without** turning the pole off-axis.
5. **Classify as one `Pole` (theta = 0), zero `OffAxis`** — matching `bug_summary` point 2. The test
   must *verify* this by calling `TropicalEval`Private`ibpDivClass` / `ibpImagPole` on the lifted
   divergent sector and printing `{poleDirs, offDirs}` (as `cc_53` Part A does at `cc_53.wl:101`).
   If the chosen `Im(A)` accidentally routes an off-axis direction (→ MultiDiv path), retune `q`
   until it is a clean single on-axis pole.

Concrete starting point (tune numerically in-session):
```
polys  = {1 + x[1] x[2] + 10^4 x[1]^2 + x[2]^2};
Avals  = {-1 + eps, I q};          (* q = 1/2, say; Im on the convergent slot only *)
pe     = {-2 + I g - eps};         (* g = 1/2 *)
lift {PolyIndex->1, ExponentVector->{2,0}, k->3}   (* z0 = 10^(4/3) *)
```
Then evaluate the two-threshold identity (§3.1) at a "mild" point (`cc_53` runs at a single default
point; a squeezed point is not required to expose the bug — `bug_summary`'s `r1 = 0.1` is mild).

### 3.4 PASS/FAIL semantics of the example itself
`cc_54.wl` is written to **PASS when the engine is correct** (lifted == unlifted == NIntegrate). So:
- On the **current** engine it is expected to **FAIL** (that failure *is* the bug confirmation, §4).
- After the §5 fix it must **PASS**. It then joins the permanent gate set.

This double-duty (repro now, regression later) is the standard `cc_*` pattern.

### 3.5 Escalation ladder (if the minimal fixture does NOT reproduce)
It is possible the minimal 2-variable, single-polynomial fixture does not trip the bug and the defect
needs fourpt's richer geometry. If §4 shows lifted == unlifted on the minimal fixture, escalate in
order, re-running §3.1 at each step:
1. Move `Im(A)` onto a slot that routes through the pivot in the flattened coords (inspect `mMatrix`).
2. Two polynomials / 3 variables (closer to a real presector).
3. The **actual fourpt `J^00` presector-1 integrand.** Its spec is not in this repo (it lived in the
   external `CROSSCHECK/fourpt/Fourptpp/`). **Action:** ask the user for the `IntegrandSpec`
   association (polynomials, `A`, `B`, kinematic symbols, the divergence map / lift rule) or for
   `probe10_samefn.wl` + `probe6_ibp_complex.wl`, and port them into `TEST/INTERFILES/cc54/` as a
   fixture. This is the guaranteed reproducer since it is the case `bug_summary` measured.

---

## 4. PHASE 2 — Determine whether the bug is real, and localize it

### 4.1 Confirm in the clean repo
Run `wolframscript -file TEST/cc_54.wl` on `6020ef9` (unmodified). Outcomes:
- **Lifted ≠ unlifted, and unlifted ≈ NIntegrate** ⇒ **bug CONFIRMED** in the clean repo (not a
  sandbox artifact). Proceed to §4.2.
- **All three agree** ⇒ bug not reproduced at this fixture size; climb the §3.5 ladder. Do **not**
  conclude "no bug" until the fourpt integrand itself has been run.
- **Unlifted ≠ NIntegrate** ⇒ the *unlifted* reference is itself suspect (a different, deeper
  problem); stop and reassess — the lift-identity test is only meaningful when the reference is sound.

### 4.2 Extract the residual factor
For each Laurent order (pole, finite) form `residual = value_unlifted / value_lifted`. `bug_summary`
gives `≈ 16.3 − 6.5 i` (finite), `|·| ≈ z0^1.5`. Fit `residual = z0^{p + q i}` (with the fixture's
known `z0`): solve `p = log|residual|/log z0`, `q = arg(residual)/log z0`. The pair `(p, q)` is the
**signature of the dropped term** and points directly at which exponent was dropped.

### 4.3 Localize with controlled toggles (cheap, read-only / minimal instrumentation)
Bisect the hypothesis space by toggling one axis at a time and re-running §3.1:

1. **`Im(A)` on/off** (`q → 0`): if the residual vanishes when `Im(A) = 0`, the dropped factor lives
   in the **`Im(A)` path** — `LiftedMonoPhase`/`MonomialPhaseLog` (`1370–1380`, `5771–5773`,
   `5855–5857`) or its interaction with `PrefactorBase`'s `ap` (`1298`). This is the §1.3 prediction.
2. **`Im(B)` on/off** (`g → 0`, keep `Im(A)`): isolates the `MonoFactorLog` path.
3. **boundary vs terms:** instrument `IBPProcessSector` (a scratch copy under
   `scratchpad/`, **never the repo file** per invariant #5-adjacent hygiene — mirror `bug_summary`'s
   sandbox discipline) to print, for the lifted sector, the numeric `z0` content of `bndPrefactor`
   (`5744`), each `flatPrefactor` (`5842`), and the `Const`/`Coeffs` of `bnd/termMonoFactorLog` and
   `bnd/termMonoPhaseLog`. Compare the *product* of `z0` powers actually emitted against the exact
   requirement that `PrefactorBase`'s `z0^{ap/mp−1}` be fully cancelled. The missing power is
   `(p, q)` from §4.2.
4. **Direct-mode control:** run the same lifted divergent integrand in `ComplexExponentMode ->
   "Direct"` (no split; complex `z0` powers stay attached to coefficients). If Direct is correct and
   `SplitRealImag` is wrong, the bug is unambiguously in the **real/imag split of the `z0` power**,
   confirming the §2 mechanism. (Direct may itself be limited for lifted divergent — if it refuses
   cleanly, note it and rely on the toggles above.)

### 4.4 Honest uncertainty
The exact dropped term is **not yet pinned** — the phase `Const` is demonstrably threaded (§2), so
the culprit is subtler than "forgot the constant." The most probable, evidence-backed localizations,
in order: **(i)** the `Im(A)` `z0` power (`LiftedMonoPhase`, only nonzero in the fourpt case);
**(ii)** a magnitude/phase mismatch in how `PrefactorBase`'s `z0^{ap/mp−1}` (realified `ap`) pairs
with the phase constants when `A` and `B` are *both* complex; **(iii)** a boundary-specific drop at
`y_k = 1`. §4.2–4.3 convert this ordering into a definite answer before any fix is written.

---

## 5. PHASE 3 — Fix

### 5.1 Shape of the fix (finalized only after §4 pins the site)
Thread the missing complex `z0` factor through the identified site so the exact-identity holds. The
fix will be one of:
- **(i) `Im(A)` path:** ensure `LiftedMonoPhase`/`MonomialPhaseLog` (and/or the realified `ap` in
  `PrefactorBase`) contribute the full `z0^{(imEaug/mp)·(…)}` — both the magnitude piece that pairs
  with the realified exponent and the phase constant — in the **boundary** (`5771–5773`) **and** the
  **terms** (`5855–5857`), matching what the convergent lifted path already does.
- **(ii) split mismatch:** correct the decomposition of `z0^{complex exponent}` into
  `z0^{Re}` (magnitude, into `PrefactorBase`/cleared coeffs) × `z0^{i·Im}` (phase, into the
  `…MonoFactorLog`/`…MonomialPhaseLog` `Const`) so the two halves reconstruct the Direct-mode product
  exactly.
- **(iii) boundary drop:** restore the `z0` factor lost when slot `k` is set to 1 / dropped in
  `bndPrefactor`/`bndMono*` (`5744`, `5761–5773`).

### 5.2 Invariants that gate the fix (from `CLAUDE.md` §"Invariants")
1. **#25 byte-identity / real-path untouched.** The new `z0` factor must be *identically 1* on the
   real-exponent / `Direct` path and on any sector where the relevant imaginary part is zero (guard on
   `lmp =!= None` / `imB` present, exactly as the existing phase code does). Re-run
   `TEST/phase2_selfgate.wl` and confirm `TEST/baselines/codegen_goldens/` are unchanged.
2. **#1 exactness.** `z0` powers stay exact rationals/radicals until codegen's `N[…]`; **no machine
   float upstream of `MmaToC`.** Assert `FreeQ[<payload>, _Real]` survives (phase2 already guards this).
3. **#3 clean `$Failed`.** The genuinely out-of-scope cases must still refuse: eps-dependent `theta`
   (`splitdivmono` negative control, `cc_53` Part A), nested/higher-order pole (`nestedIBP`), Case-B
   domain coupling (`liftdivdomain`). The fix must not turn any of these into a number.
4. **#4 two oracles.** `cc_54`'s PASS rests on lift==no-lift **and** IBP≈NIntegrate (and, if added,
   IBP≈Subtraction as in `cc_53` Part B).

### 5.3 Regression gate (definition of done)
- `wolframscript -file TEST/cc_54.wl` ⇒ `CC54 PASS` (was FAIL pre-fix).
- `TEST/cc_53.wl` still PASS (regulator-real fix intact; complex-B/real-A still correct).
- `TEST/cc_45.wl … cc_52.wl` still PASS (existing IBP / lift / complex cross-checks).
- `TEST/phase2_selfgate.wl` and `TEST/phase3_selfgate.wl` ⇒ exit 0 (goldens + exactness + complex-A/B).
- If the fourpt integrand was ported in (§3.5 step 3): its lifted divergent presector 1 IBP
  `(pole, finite)` matches the **unlifted** IBP reference `bug_summary` already has in hand
  (`probe6`: pole `-5.4e-22 + 2.5e-21 i`, finite `-8.0e-19` at `r1 = 10^{-3.5}, r2 = 0.9`).

### 5.4 Documentation
Update `SUMMARY.txt` (limitation/known-good matrix) and, if the fix changes the supported surface,
`planIBPCX.md` §8 / `planCXLIFTDIV.md`. Add a one-line pointer to this plan. Commit message should
credit that the fix completes the `6020ef9` feature (path now *correct*, not merely *running*).

---

## 6. Secondary, separately-tracked item (do NOT bundle into this fix)
`new_request` §0 flags a distinct real bug on the **Subtraction/Laurent** route (moot for IBP):
`emitDomainIndicatorCpp` emits `double log_ypstar = (… cx(…) …)` for a lifted `SplitRealImag` domain
indicator, which g++ rejects (`cannot convert std::complex<double> to double`); the indicator is
provably real, so a `.real()` fixes it (`tropical_eval.wl:2362–2377`, the `emitDomainIndicatorCpp`
region near `2406–2452`). **Out of scope here** — note it, file it, keep this PR to the z0 bug so the
byte-identity/regression story stays clean.

---

## 7. Key file/line index (all `tropical_eval.wl` unless noted)
- Lift anchor + `LiftData`: `941–942`, `1014`.
- Cleared-poly `z0` coeffs: `1093`.  `PrefactorBase`: `1298`.
- `MonoFactorLog`/`LiftedMonoFactor` (Im B): `1327`, `1344–1355`.
- `MonomialPhaseLog`/`LiftedMonoPhase` (Im A, `None` unless complex A): `1370–1380`.
- Divergent lifted `SectorData`: `1389–1424`.
- Realification in driver (`SplitRealImag`): `4207–4229`, `4305–4311`, `4373–4392`;
  `imagPolyInfo`/`realifyPolyB`: `~2366–2390`; `realifyMonoA`/`imagMonoInfo`: `~2336–2360`.
- Single-pole IBP path: `IBPProcessSector` `5586`+; classification/guards `5622–5648`, `5665–5695`;
  boundary `5711–5785` (`bndPrefactor` `5744`, phase `Const` `5763`/`5772`); terms `5787–5886`
  (`flatPrefactor` `5842`, phase `Const` `5849`/`5856`); `DLogPrefactor` `5924`.
- MultiDiv path (not implicated, for contrast): `ibpProcessLeaf` `5297`, `IBPBuildCorners` `5369`,
  `IBPProcessSectorMultiDiv` `5450`.
- `ibpImagPole` (regulator-real Refine): `5069`+.
- Driver `EvaluateTropicalMC` `4099`; `EvaluateTropicalMCLifted` (the two-threshold toggle) `4809–4893`.
- Existing near-miss test: `TEST/cc_53.wl` (Part B `146–219`).
