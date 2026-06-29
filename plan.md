# Implementation Plan: TROPICAL_MONTE_CARLOv3

**Status:** draft plan, ready for review.
**Working directory / target codebase:** `/Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3/`
**Reference codebases (READ-ONLY — never edit):** `OLD_CODE/TROPICAL_MONTE_CARLO/` (Tree A), `OLD_CODE/TROPICAL_MONTE_CARLO2/TROPICAL_MONTE_CARLOv2/` (Tree B).

This document is written for an implementing agent who has **not** read the old codebases. Every design decision, file, function, formula, fix, verification command, and PASS criterion needed to build v3 is included or pinned to a precise location in `OLD_CODE/`. It is modelled on the (excellent) v1→v2 migration plan at `OLD_CODE/TROPICAL_MONTE_CARLO2/PLAN_SUM/plan.md`, which used strict phase-gates; v3 keeps that discipline. **Phases must be executed in order, and each phase has a gate that must pass before the next begins.**

> One-sentence summary: **v3 = the divergence/IBP engine (Tree A) ∪ the lifting/complex engine (Tree B), refactored into a single non-duplicated package, with the duplicated batching *implementation* folded into the unified codegen (batched VEGAS kept as an option), the complex-lifting fixes from the bug logs baked in, and an extensive cross-check suite.**

---

## §0 How to read this document

- §1 — what v3 is (feature set) and is not (non-goals).
- §2 — source-of-truth analysis: exactly what code exists in `OLD_CODE`, which copy is canonical for each feature, and what is redundant.
- §3 — target architecture: module layout, files, public API, data structures.
- §4 — key design decisions, each with a recommendation and rationale (the few places a reviewer may want to overrule).
- §5 — the refactor / simplification / deletion catalog (the "remove redundant code" requirement, made concrete with line numbers).
- §6 — mathematical & codegen specs that MUST be transcribed correctly (the subtle, bug-prone parts — especially the complex-lifting fixes that exist only in the bug logs).
- §7 — the phased implementation plan with gates.
- §8 — the comprehensive cross-check suite (the "many cross checks" requirement), as concrete tests with PASS criteria.
- §9 — edge cases, risks, known limitations.
- §10 — execution conventions.
- Appendix A — file inventory of `OLD_CODE`. Appendix B — old→v3 function migration table.

---

## §1 What v3 is

### 1.1 The integral and the pipeline

v3 numerically evaluates generalized Euler integrals

```
I = ∫_{[0,∞)^n} dx_1 … dx_n  ∏_i x_i^{A_i}  ∏_j P_j(x)^{B_j}
```

with real or **complex** monomial exponents `A_i` and polynomial exponents `B_j` (possibly ε-dependent), polynomials `P_j` with real, complex, or kinematic-dependent coefficients. The strategy is unchanged from the old code:

```
(a) tropical fan  →  (b) symbolic sector processing  →  (c) flattening  →  (d) [optional] numerical integration
```

**The emphasis of v3 is the tropical decomposition (a–c), not the numerical back-end (d).** The decomposition is the product; numerical integration is one consumer of it. Concretely this means:
- The symbolic pipeline (fan, sector processing, divergence handling, lifting, flattening) is the core, is exact, and is validated symbolically (against `NIntegrate`, exact closed forms, FIESTA) independently of any sampler.
- The C++ Monte-Carlo / VEGAS code generator is a back-end that consumes the decomposition. It is kept, but simplified, and is **not** where new complexity is added.

### 1.2 Required features (all five must be present)

| # | Feature | Comes from | Notes |
|---|---|---|---|
| F1 | **Tropical decomposition** (fan + `ProcessSector` + flattening) | shared core (both trees, byte-identical helpers) | The heart of the package. |
| F2 | **1/ε divergences** | Tree A (Module 2 subtraction + Module 2b IBP) | Two methods that cross-check each other; both currently handle a *single* divergent variable per sector. |
| F3 | **Auxiliary-variable uplift** for extreme coefficients | Tree B (Module 1b/1c) + bug-log fixes | Real coefficients work in Tree B; complex coefficients & complex exponents need the bug-log fixes (§6). |
| F4 | **Complex exponents** (`A_i`, `B_j` ∈ ℂ) | Tree B (default path) + bug-log fixes (SplitRealImag) | The default `exp(B·log P)` path is correct; the VEGAS-only SplitRealImag optimization needs the §6.3 fixes. |
| F5 | **Optional VEGAS** numerical integration (incl. batched mode) | both trees (CUBA-Vegas) | MC is the dependency-free default; VEGAS opt-in via `Integrator -> "VEGAS"`, with an opt-in kinematic-scan `Batch` mode (§5.1). |

### 1.3 Explicit non-goals (removed or out of scope)

- **N1 — *Redundant* batching code (NOT the feature).** Batched VEGAS is **kept** as a driver option (§5.1, confirmed by the user). What is dropped is the *redundancy*: (i) the *duplicated* batch-`main()` implementation (`emitVegasBatchMain`), which is folded into the unified VEGAS emitter rather than kept as a separate tower; and (ii) the standalone batch *benchmark* harness (`TEST/bench_batch.cpp`, the 7200-set study — see N2). The batch *capability* (one low-discrepancy sample set shared across a chunk of kinematic points) survives, VEGAS-only.
- **N2 — The standalone numerical-integration benchmark study.** The `TEST/` sampler study (bench.cpp / bench_batch.cpp / 769-row sweeps) was a one-time investigation that already concluded "VEGAS wins, MC is the weak baseline." Its *conclusion* is baked into v3 (VEGAS is the recommended sampler); the heavy benchmark harness is **not** carried forward as a maintained part of v3 (it can stay archived in `OLD_CODE`).
- **N3 — Lifting a divergent integral simultaneously.** *Supported as of planAXpDIV.md (Case A).* A lifted sector may also carry a single 1/ε pole: after the lifting δ-resolution collapses the aux dimension, a lifted sector is structurally an ordinary `SectorData` of dimension `n` whose pole lives in one of the `n` post-δ effective exponents `ãtilde` — exactly the object the IBP/subtraction machinery consumes. `ProcessSectorLifted` is ε-aware; it carries a `PrefactorBase` and a `DomainConstraint`, and both the IBP and `LaurentFromSubtraction` routes resolve the Laurent (cross-check 41-C). **Complex polynomial exponents with lifting + a 1/ε pole are now supported** via `ComplexExponentMode -> "SplitRealImag"` on both routes (planCXLIFTDIV.md, cc_46): decompose on `Re(B)`, reintroduce `Im(B)` as an oscillatory phase (per-IBP-piece `MonoFactorLog` + the brought-down complex IBP coefficient — the `i·Im(B_j)` phase-derivative term). **Complex *monomial* exponents `A` with lifting + a 1/ε pole are also now supported** via the same mode (planAXpDIVv2.md, cc_47) — the principal-series bubble at imaginary `μ`, which previously aborted on every cone with `liftcomplex`. `Re(A)` is realified on the same footing as `Re(B)` (subtract the ε-free `Im(A)`, preserving the real regulator ε that lives in the monomial exponent), and `Im(A)` rides a new bare-monomial phase `exp(i (Const_A + Σ_a (Im(A)·M)_a/α0_a log y'_a))` — the monomial twin of the `MonoFactorLog` phase, with the lifted constant `Const_A = (Im(eaug)_p/m_p)·logZ0`. Convergent (cc_47 A) matches NIntegrate; divergent Case A (cc_47 B) agrees across IBP and pinned-ε Subtraction, including a genuinely complex pole. **Future work (refused with a clean `$Failed`, never a wrong number):** Case B — the divergent variable couples to the lifted domain face for *every* admissible pivot (`TropicalEval::liftdivdomain`; the caseA-first ranking makes this rare, cc_44); the inline symbolic-ε Subtraction path for SplitRealImag × lift × divergence (`TropicalEval::splitliftdiv`; use Automatic/IBP or pinned-ε Subtraction instead); its complex Case-B analogue, now covering a nonzero imaginary exponent in the divergent direction from **either** `Im(B)·D` **or** the bare monomial `Im(A)·M` (`TropicalEval::splitdivmono`, cc_47 C); and `Direct`-mode complex-`A`/`B` lifting (the pre-existing `liftcomplex` limitation — use SplitRealImag). The single-pole-per-sector limit (L1/L8) still applies.
- **N4 — Nested divergences (1/ε^{d≥2}).** Both old divergence methods refuse sectors with >1 divergent variable; v3 keeps that refusal (clean `$Failed`). Multi-order pole assembly is future work. (Full integral at a fixed ε is still reachable via `EpsilonValue`.)
- **N5 — Multiple independent auxiliary variables.** A single `z` with per-monomial powers and residual coefficients handles multiple extreme coefficients (§6.2); multiple independent `z`'s are future work.

### 1.4 The exactness invariant (core design principle — confirmed by the user; non-negotiable)

**All tropical-decomposition code is EXACT / symbolic. Floating-point arithmetic appears in exactly two places, and neither is the decomposition itself:**
1. the **numerical integration** of the flattened sector integrands (MC / VEGAS) — the deliberate, opt-in numeric step; and
2. **independent numerical cross-checks** (NIntegrate, CUBA) — which only *validate* a result and never feed back into it.

Every quantity the decomposition produces or *decides on* is kept exact (integers, exact rationals, exact radicals / `Power[C,1/k]`, exact complex numbers; ε and kinematics stay symbolic): ray matrices and `det M`, transformed/cleared exponent vectors, tropical minima, effective exponents, the flattening exponent-division, the lift anchor `z0` and residuals `c = C/z0^k`, the delta-elimination substitution, the convergence/divergence decision, the pivot-admissibility decision, the IBP reduction, the G0/G1/remainder construction, and the assembled Laurent (pole, finite) coefficients. **No decision inside the decomposition may be made by evaluating a float** — e.g. "is `Re[a] ≤ 0`?" or "is this `atilde` real?" must be decided *exactly*, never by a `10^-12` tolerance. Numericization happens only at the `MmaToC` codegen/integration boundary and in the validation layer.

**Why this is the linchpin for unattended running (the user's concern).** A subtly-wrong *numeric* value with a tight-but-lying error bar is the one thing that can slip a false "PASS" past an automated gate (cf. bug-log L1/L5) — and it cannot be eyeballed. The exactness invariant removes numerics from the place where correctness is *decided*, so the decomposition is verifiable by **exact symbolic identities** (which cannot be "subtly wrong"), and the only genuinely numeric output — the integrated value — is gated against **multiple independent references** by a **separate adversarial verifier agent** (§8.5). The lone heuristic permitted to touch a float is `DetectExtremeCoefficients` *choosing which* coefficient to lift; that choice **cannot affect correctness**, because the lift is an exact identity for any choice (it changes only variance). Implementation discipline in §6.7; verification mechanism in §8.5.

---

## §2 Source-of-truth analysis

### 2.1 The three engine copies in `OLD_CODE` (and which is canonical)

| Engine | Path | `tropical_eval.wl` size | Has | Lacks | Verdict |
|---|---|---|---|---|---|
| **Tree A** | `OLD_CODE/TROPICAL_MONTE_CARLO/` | 4584 lines | divergence (subtraction + IBP), VEGAS, **batching**, factored codegen (`emit*`) | lifting | **Canonical base for F2, F5, codegen.** Newest (Jun 14). |
| **Tree B** | `OLD_CODE/TROPICAL_MONTE_CARLO2/TROPICAL_MONTE_CARLOv2/` | 2666 lines | lifting (real-exp), complex exponents, VEGAS (no batch), `FlattenSector` refactor, `SeedBase` runtime override | divergence | **Canonical source for F3, F4, `FlattenSector`.** |
| **Tree B-inner** | `OLD_CODE/TROPICAL_MONTE_CARLO2/TROPICAL_MONTE_CARLO/` | 3834 lines | older pre-VEGAS divergence snapshot | VEGAS, lifting, batching | **REDUNDANT.** Superseded by Tree A. Do not carry anything forward from it. |

`tropical_fan.wl` is **byte-identical (md5 `256a8da…`) across all three** — the Polymake fan layer is stable. v3 carries it forward essentially unchanged (one additive hardening: K-scaled fans, §6.4).

### 2.2 The advanced complex-lifting work is NOT in the source — only in the bug logs

The two bug logs in Tree A —
`OLD_CODE/TROPICAL_MONTE_CARLO/TMCv2_BUG_LOG.md` and `OLD_CODE/TROPICAL_MONTE_CARLO/lift_error_log.md` —
document a *separate* branch (`complex-lifting-hd-vegas`, repo CALC2) that combined lifting + complex exponents + VEGAS in high dimension. **None of its symbols exist as source in `OLD_CODE`** (verified by grep: `ComplexExponentMode`, `SplitRealImag`, `computeFanScaled`, `resolveVegasSizing`, `MonoFactorLog`, `vegasbudget`, `test_complex_lifted` appear *only* inside those two `.md` files). They are therefore **specifications to (re)implement in v3**, not code to copy. §6.3–6.5 transcribe the fixes; §8 turns their tests into v3 cross-checks. This is the single most important scoping fact in this plan: **treat the bug logs as the design spec for correct complex lifting.**

### 2.3 Shared helpers that are already identical (low-risk to unify)

Verified byte-identical between Tree A and Tree B (diff of the corresponding blocks):
`ParsePolynomial`, `TransformExponents`, `MmaToC` / `mmaToCInternal` (all 18 down-values, incl. `Complex→cx`, `Re/Im`), `GenerateMonomialSumCpp`.
→ In v3 these become a single shared module (§3.1, Module 0) with high confidence.

### 2.4 Shared functions that DIVERGED (must be reconciled deliberately)

| Function | Tree A behavior | Tree B behavior | v3 reconciliation |
|---|---|---|---|
| `ProcessSector` | threads `eps`; divergence check at `eps→0`; flattening inline | no `eps` (hard `noregulator` guard); flattening extracted to `FlattenSector` | **Keep eps threading (A) AND the `FlattenSector` factorization (B).** ProcessSector calls `FlattenSector` with `eps` passed through. (§6.1) |
| `GenerateCppMonteCarlo` | factored `emit*` + `Which`-dispatch incl. batch; hard `cppBadTokens` gate; `"Integrator"->"MonteCarlo"` | 430-line monolith; warning-only bad-token scan; `"Integrator"->"MC"` | **Adopt A's factored skeleton, B's integrator vocabulary, A's hard `badcpp` gate** (§4.4, §5.2). |
| `CompileCpp` | memoized `detectCuba`, `TROPICAL_REQUIRES_CUBA` sentinel | `findCubaPrefix` (duplicated with `cuba_common.wl`) | **Keep A's memoized `detectCuba`; delete `findCubaPrefix`; single CUBA-probe source** (§5.3). |
| `EvaluateTropicalMC` | divergence routing, `EpsilonValue`, batch options | lifting routing (`LiftData`), `SeedBase` argv | **Unify into one driver with Method/Lift/Integrator routing + thin wrappers** (§3.4, §5.4). |
| Integrator vocabulary | `"MonteCarlo"` / `"Vegas"` | `"MC"` / `"VEGAS"` | **Pick `"MC"`/`"VEGAS"`** (§4.3). |

---

## §3 Target architecture

### 3.1 Module layout of `tropical_eval.wl` (one package, `TropicalEval``)

The biggest simplification over the old code: the parallel **convergent / IBP** code-generation and driver "towers" (which together account for ~1900 lines of near-duplicated codegen in Tree A) collapse into one parametrized path. Proposed modules:

| Module | Contents | Sourced from |
|---|---|---|
| **M0 — Primitives** | `ParsePolynomial`, `TransformExponents`, `MmaToC`/`mmaToCInternal`, `GenerateMonomialSumCpp` | shared (identical) |
| **M1 — Sector processing** | `ProcessSector`, `FlattenSector` (eps-aware), `CheckFlatteningMagnitude`, `evalFlattenedIntegrand` (new, extracted from triplicated log-exp eval) | A (eps) + B (`FlattenSector`) |
| **M1b — Lifting** | `DetectExtremeCoefficients`, `LiftCoefficients`, `ProcessSectorLifted`, anchor-selection (`kStar`) | B + bug-log fixes |
| **M2 — Divergence: tropical subtraction** | `IdentifyDivergences`, `ProcessDivergentSector` | A |
| **M2b — Divergence: IBP** | `IBPReduceSector`, `IBPProcessSector`, `IBPCheckBoundary` | A |
| **M3 — Unified C++ codegen** | one `emitIntegrandDefinitions` (parametrized by integrand-type tag), one `emitMain` (parametrized by `MC`/`VEGAS`), `GenerateCpp`, `CompileCpp`, `detectCuba`, `resolveVegasSizing` | A (skeleton) + bug-log (sizing) |
| **M4 — Driver** | one `EvaluateTropicalMC` (routes on `Method ∈ {Auto, None, Subtraction, IBP}`, `Lift`, `Integrator`); thin wrappers `EvaluateTropicalMCLifted`, `LaurentFromSubtraction` | unify A + B |
| **M5 — Validation** | `ValidateDecomposition` (parametrized to also serve the lifted case), `ValidateSubtraction`, `ValidateIBP`, `ValidateLiftedDecomposition` | A + B (de-duplicated) |

`tropical_fan.wl` stays a separate package (`TropicalFan``), with one additive private helper `computeFanScaled` (§6.4).

### 3.2 Files (created at the v3 repo root, alongside `OLD_CODE/`)

```
TROPICAL_MONTE_CARLO3/
  tropical_fan.wl            (from the identical OLD copy + computeFanScaled)
  tropical_eval.wl           (the merged, de-duplicated engine)
  SUMMARY.txt                (rewritten for v3)
  MANUAL/manual.tex          (merged from both manuals; convergent + divergence + lifting + complex)
  EXAMPLES/                  (worked examples — see below)
  TEST/                      (the cross-check suite — §8; lean, not the old benchmark harness)
  SANDBOX/                   (lifting/complex sandbox scripts, used during Phase 3–4 then retained)
  CC/                        (FIESTA cross-checks, ported with path fixes — optional, §8 item 8)
  INTERFILES/                (gitignored; generated C++, binaries, kinematic/result files)
```

> Naming note: keep the package contexts `TropicalEval`` and `TropicalFan`` (not "V3") so that scripts and any downstream callers (e.g. CALC2) are drop-in compatible. The directory is the version marker.

### 3.3 Public API (one consolidated surface)

Carry forward, after de-duplication:

```
(* fan — tropical_fan.wl *)
PolytopeVertices, ComputeDecomposition, ComputeDecompositiony, convexHullVertices, translateToOriginInteger

(* core decomposition *)
ParsePolynomial, ProcessSector, FlattenSector, CheckFlatteningMagnitude, ValidateDecomposition

(* lifting *)
DetectExtremeCoefficients, LiftCoefficients, ProcessSectorLifted, ValidateLiftedDecomposition

(* divergence *)
IdentifyDivergences, ProcessDivergentSector, ValidateSubtraction,
IBPReduceSector, IBPProcessSector, IBPCheckBoundary, ValidateIBP

(* codegen + back-end *)
MmaToC, GenerateCpp (renamed from GenerateCppMonteCarlo; the IBP variant is absorbed), CompileCpp, detectCuba

(* drivers *)
EvaluateTropicalMC (unified), EvaluateTropicalMCLifted (wrapper), LaurentFromSubtraction (wrapper)

(* tests *)
RunAllTests
```

**Removed from the public surface** (subsumed or deleted): `GenerateCppMonteCarloIBP`, `EvaluateTropicalMCIBP` (absorbed into `GenerateCpp` / `EvaluateTropicalMC` via a `Method` option); `findCubaPrefix` (private; replaced by `detectCuba`); `RunBenchmark2D` (dead). **Retained:** the `Batch` / `CubaMaxComp` options (batched VEGAS, §5.1).

### 3.4 Driver options (unified `EvaluateTropicalMC`)

```
"Method"        -> Automatic   (* Automatic: "None" if no divergent sectors, else "IBP" (the DEFAULT
                                    divergence method). Force with "None"|"IBP"|"Subtraction".
                                    Both IBP and Subtraction are kept and cross-check each other (#4). *)
"Lift"          -> None        (* None | Automatic | <lift-rules>; composes with a simple 1/eps pole (Case A; planAXpDIV.md / N3) and with complex B via ComplexExponentMode->"SplitRealImag" (planCXLIFTDIV.md) *)
"Integrator"    -> "MC"        (* "MC" (default, no deps) | "VEGAS" (needs CUBA) *)
"EpsilonValue"  -> None        (* substitute eps->value to get the full integral at fixed eps *)
"NSamples"      -> 10^6
"NThreads"      -> Automatic
"SeedBase"      -> 42          (* also passed as runtime argv → re-run same binary at new seed *)
"RunChecks"     -> True        (* symbolic Validate* gating before the expensive run *)
"PrecisionGoal" -> 3
"Verbose"       -> True
"WorkingDirectory" -> Automatic
(* VEGAS tuning, all Automatic-resolved by resolveVegasSizing (§6.5): *)
"VegasEpsRel" -> Automatic, "VegasEpsAbs" -> Automatic, "VegasSeed" -> 0,
"VegasNStart" -> Automatic, "VegasNIncrease" -> Automatic, "VegasNBatch" -> Automatic, "VegasMinEval" -> 0,
(* VEGAS kinematic-scan batching — VEGAS only; convergent / lifted-convergent path (divergent/IBP falls back to per-kp).
   NB: "Batch" (share a sample set across kinematic points) is distinct from "VegasNBatch" (a CUBA per-call tuning param). *)
"Batch" -> False, "CubaMaxComp" -> 512
(* complex exponents under VEGAS: *)
"ComplexExponentMode" -> Automatic   (* Automatic: "Direct" exp(B·logP); opt-in "SplitRealImag" for VEGAS variance, §6.3 *)
```

### 3.5 Data structures (unchanged Association schemas — see SUMMARY for full keys)

`IntegrandSpec`, `SectorData` (+lifted keys `DomainConstraint`, `PivotIndex`, `ZRow`, `HasConstantTerm`, `AugmentedA`, and new `MonoFactorLog`, §6.3), `DivergentSectorData`, `IBPSectorData`. Keep them stable for drop-in compatibility.

---

## §4 Key design decisions (recommendations a reviewer may want to confirm)

**D1 — Base engine = Tree A, graft lifting from Tree B.**
Rationale: Tree A is the superset (divergence + IBP + VEGAS + already-factored codegen) and is the newest. Lifting (Tree B) is an additive module. Going A→(remove batch *redundancy*)→(unify core)→(add lifting) is less destructive than B→(add divergence back). *Recommended.*

**D2 — Keep BOTH divergence methods; IBP is the DEFAULT, subtraction is the cross-check. (confirmed)**
Both are kept: they are independent algorithms, and their agreement on (pole, finite) is one of the strongest correctness signals (cross-check #4; old `test_divergent_crosscheck.wl` shows ≲1e-3 agreement). **`Method -> Automatic` resolves to `"IBP"`** for divergent integrals; tropical subtraction runs as the validating cross-check and is reachable on demand via `Method -> "Subtraction"` (and the `LaurentFromSubtraction` wrapper). IBP is the default because each IBP term is a *convergent* integral (effective exponents ≥1 at ε=0), avoiding the explicit pole/G0/G1/remainder bookkeeping of the subtraction scheme. The redundancy to remove is the *duplicated code generation* behind the two methods (§5.2), **not** either method.

**D3 — Lifting + a simple 1/ε pole is supported (Case A; N3, planAXpDIV.md), including complex monomial AND polynomial exponents (planCXLIFTDIV.md / planAXpDIVv2.md).** The (n+1)-dim lifted fan and the post-δ pole machinery compose; the divergent lifted sector reconstructs prefactors from `PrefactorBase` and threads a `DomainConstraint`. Complex-`B` (cc_46) AND complex-monomial-`A` (cc_47) lifting + divergence are supported via `ComplexExponentMode -> "SplitRealImag"` on both the Subtraction and IBP routes: `SplitRealImag` splits **both** families — `Re(A)`, `Re(B)` build the real measure/domain map, `Im(A)`, `Im(B)` ride the oscillatory phase. Case B (domain–pole coupling), the inline symbolic-ε Subtraction path for SplitRealImag × lift × divergence, its complex Case-B analogue (nonzero `Im(A)·M` or `Im(B)·D` in the divergent direction), and `Direct`-mode complex-exponent lifting are refused with clean `$Failed` messages. *Recommended.*

**D4 — Reconstruct complex lifting from the bug logs.** F3+F4 jointly require it. Implement complex *coefficients* (already mostly in Tree B's residual mechanism) and complex *exponents* in lifting (the SplitRealImag mode with the §6.3 fixes). *Recommended; this is the main genuinely-new engineering in v3.*

**D5 — Numerical integration is opt-in and minimal.** Default `Integrator -> "MC"` keeps zero new dependencies and byte-identical MC codegen; `"VEGAS"` needs CUBA and hard-fails (no silent fallback) if absent. VEGAS additionally offers an opt-in **kinematic-scan batch mode** (`Batch -> True`) that shares one low-discrepancy sample set across a chunk of kinematic points — a large speedup for big scans (old code: ~30–60× wall-clock on a 7200-point scan); implemented cleanly inside the unified codegen, not as a duplicate tower (§5.1). *Recommended (matches the user's "focus on decomposition, less on numerical integration").*

**D6 — Single integrator vocabulary `"MC"`/`"VEGAS"`** (§4.3 below).

**D7 — Exactness invariant (§1.4).** The tropical decomposition is exact/symbolic; floating-point lives only in MC/VEGAS integration and in independent cross-checks — never in a decomposition *decision*. This is what makes unattended gates trustworthy, via exact-identity verification by an independent agent (§8.5). *Confirmed by the user; non-negotiable.*

### 4.3 Vocabulary / convention unifications
- Integrator strings: **`"MC"` / `"VEGAS"`** (Tree B convention — newer, closer to v3's feature set).
- CUBA detection: **memoized `detectCuba[]`** (Tree A) is the single source of truth; `cuba_common.wl`'s probe (used by the CUBA cross-checks) calls into it rather than duplicating.
- `nocuba` behavior: **hard `$Failed`** when VEGAS requested but CUBA absent (both trees already do this — preserve, no silent MC fallback).
- Keep `SeedBase` **both** as a compile-time option and a runtime `argv[5]` override (Tree B) — lets a caller re-run the same binary under new seeds without recompiling (used by the determinism cross-check #22).

### 4.4 Codegen architecture
Adopt Tree A's **"shared Defs prefix + selectable `main()`"** factoring, but parametrize the integrand emitter by an **integrand-type tag** (`conv | g0 | g1 | rem | ibp-boundary | ibp-term`) and a **lifted-domain-indicator** flag, so there is exactly **one** `emitIntegrandDefinitions` and **one** `emitMain[kind]` with `kind ∈ {MC, VEGAS}`. This collapses the four Tree-A emit functions + their IBP twins into two parametrized functions (§5.2). Keep the hard `cppBadTokens`/`badcpp` gate (Tree A) — it prevents emitting non-compilable C++ when a nested divergence leaks `ComplexInfinity` or an unconverted symbolic head.

---

## §5 Refactor / simplification / deletion catalog

This section is the concrete realization of "refactor/simplify the old code base as much as possible … and remove redundant code, such as the old batching." Every item cites the Tree-A or Tree-B location. **Each item is verified non-breaking by the gate listed.**

### 5.1 REFACTOR (not delete): batched VEGAS — keep the capability, remove the redundancy

**Decision (confirmed by the user): v3 RETAINS batched VEGAS as a driver option** (`Integrator -> "VEGAS"`, `Batch -> True`) — it shares one low-discrepancy sample set across a chunk of kinematic points, the big speedup for large scans. The original "remove the old batching" instruction targets the *redundancy*, and that is what v3 removes — not the capability. All line numbers are in **Tree A** `tropical_eval.wl`.

- **FOLD, don't duplicate.** Tree A implements batching as a standalone `emitVegasBatchMain` (**lines 1926–2033**, ~107 lines) that re-derives the entire per-kp VEGAS `main()`. In v3 this becomes a **mode of the single VEGAS emitter** — `emitMain["VEGAS", "Batch" -> True]` — sharing the one `emitIntegrandDefinitions` Defs prefix (§5.2). The batch-specific pieces (chunked ncomp over the CUBA component axis, the `gb_*` globals, the `cubaBatch` callback, the `CubaMaxComp` ncomp ceiling) live only inside that mode. No duplicated tower; that elimination *is* the "remove redundant code."
- **KEEP the options** `"Batch" -> False` (default) and `"CubaMaxComp" -> 512` on `GenerateCpp` and `EvaluateTropicalMC`. Keep the `batch` local and its `vegasOpts` plumbing in the unified driver.
- **SCOPE: VEGAS + convergent / lifted-convergent only.** The divergent/IBP path falls back to per-kp VEGAS — Tree A's IBP "batch" was only a no-op note at **4205–4207** (delete the note); document the fallback rather than implementing a divergent batch main.
- **DROP the standalone batch *benchmark* harness** (`TEST/bench_batch.cpp` + the 7200-set sweep) — a one-time study, archived in `OLD_CODE`, not maintained in v3 (non-goal N2). This is the bulk of the genuinely redundant batching code.

The structural map confirmed batching only changes the emitted `main()`; the result-file format (4 floats/kp) is identical for both modes, so **Mathematica-side readback needs no change**.

**Gate:** `Integrator->"MC"` byte-identical to the MC golden; per-kp `"VEGAS"` matches the per-kp golden; **batched `"VEGAS"` reproduces the per-kp VEGAS value within combined error (cross-check #23)** and matches its own captured golden; full divergence suite green. (Cross-checks #23, #25.)

### 5.2 UNIFY: the duplicated convergent/IBP codegen towers

Tree A duplicates the entire codegen pipeline for convergent vs IBP (the single biggest redundancy in the codebase). Collapse:

| Tree A pair | → v3 |
|---|---|
| `emitIntegrandDefinitions` (1276–1640) + `emitIntegrandDefinitionsIBP` (3523–3940) | one `emitIntegrandDefinitions[sectors, spec, opts]` parametrized by per-integrand type tag |
| `emitMonteCarloMain` (1647–1801) + `emitMonteCarloMainIBP` (3947–4057) | one `emitMain["MC"]` (result-assembly driven by tags) |
| `emitVegasMain` (1815–1918) + `emitVegasMainIBP` (4067–4174) + `emitVegasBatchMain` (1926–2033) | one `emitMain["VEGAS"]` with a `Batch -> True` mode (§5.1) |
| `GenerateCppMonteCarlo` (2065) + `GenerateCppMonteCarloIBP` (4192) | one `GenerateCpp` |

**Gate:** golden-master codegen regression (capture v3's emitted C++ for a fixed set of specs once, then diff on every later change) — see cross-check #25; plus the full divergence suite must reproduce the old (pole, finite) values.

### 5.3 UNIFY / DELETE: CUBA detection
- Keep Tree A's memoized `detectCuba[]` (2134–2153). Delete Tree B's `findCubaPrefix` (2005). Make `EXAMPLES/cuba_common.wl`'s probe call `detectCuba[]` so there is **one** CUBA-locating code path (Tree B's comment at its 2001–2002 explicitly warns the two must not diverge — v3 removes the divergence by construction).

### 5.4 UNIFY: drivers
Collapse `EvaluateTropicalMC` + `EvaluateTropicalMCIBP` + `EvaluateTropicalMCLifted` into one `EvaluateTropicalMC` (the three share a ~10-step skeleton: workdir → fan unpack → per-sector processing → optional Validate → codegen → compile → run → readback → assemble). Routing on `Method`/`Lift`/`Integrator`. Keep `EvaluateTropicalMCLifted` and `LaurentFromSubtraction` as **thin wrappers** for API continuity.

### 5.5 DE-DUPLICATE: validation + log-exp evaluation
- Extract the triplicated log-exp integrand evaluation (Tree B: `CheckFlatteningMagnitude`, `ValidateDecomposition`, `ValidateLiftedDecomposition`) into one private `evalFlattenedIntegrand[sectorData, yVars, kinRules]`.
- Parametrize `ValidateDecomposition` so `ValidateLiftedDecomposition` reuses it with a sector-processor + domain-indicator argument (Tree B near-duplicates them).
- In `ProcessSectorLifted`, call `tryPivot[p]` **once** per pivot and cache (Tree B calls it up to 3×: candidate loop, complex-check loop, summary builder — Tree B 957–1006).

### 5.6 DELETE: dead code & vestigial declarations
- `RunBenchmark2D` — dead in both trees (Tree A 560–602; Tree B 1353). Delete.
- Tree B's vestigial divergence error messages (`divergent`, `nested`, `badck` at B 174–176) — but v3 *re-introduces* divergence, so these come back live from Tree A instead. Just don't carry Tree B's dead copies.
- Tree B-inner (the 3834-line copy): ignore entirely (§2.1).

### 5.7 SHARE: option-default boilerplate
The VEGAS option defaults are repeated across ≥4 option blocks (Tree A 1813/1926/2040/4065/4180/4256). Define one `$vegasOptionDefaults` list and `Join` it in. (After §5.1 this also stops `CubaMaxComp` from reappearing.)

### 5.8 HARDEN: the `MmaToC` CForm fallback
The recursive converter's CForm fallback (Tree A 1205–1213) is exactly the path that can leak `Power[`/`Times[` tokens that the `badcpp` gate later catches. Document the coupling; add a targeted test (cross-check #27) that feeds an exotic head and asserts a clean `$Failed` rather than uncompilable C++.

---

## §6 Specs that MUST be transcribed correctly

These are the subtle, bug-prone parts. The first is a refactor; the rest are the complex-lifting fixes that **exist only in the bug logs** (§2.2) and must be reconstructed.

### 6.1 ProcessSector ⊕ FlattenSector with eps threading (reconcile A & B)

Tree B factored the divergence-check + flatten + prefactor tail of `ProcessSector` into `FlattenSector[clearedPolys, effectiveAVals, prefactorBase]` (B 232–277) — but removed `eps`. Tree A keeps `eps` inline. v3 keeps **both**: `FlattenSector[clearedPolys, effectiveAVals, prefactorBase, eps:None]`, where the divergence predicate is evaluated at `eps→0` when `eps =!= None`:
```
a0 = If[eps =!= None, effectiveAVals /. eps -> 0, effectiveAVals];
isDivergent = Or @@ (TrueQ[Re[#] <= 0] || (NumericQ[#] && Re[#] <= 0) & /@ a0);
```
Preserve Tree A's "last violating index wins" loop semantics. `ProcessSector` calls it with `prefactorBase = Abs[detM]`; `ProcessSectorLifted` calls it with the lifted prefactor (§6.2). **Critical (from the v1→v2 plan §1.3 rule 4):** the `IsDivergent -> True` pre-flattening early-return branch of `ProcessSector` must survive — the lifted path consumes `ClearedPolys`/`NewExponents` from it, because (n+1)-dim pre-delta sectors are *routinely* flagged divergent even when the delta-constrained integral is finite.

### 6.2 Auxiliary-variable lifting (real and complex coefficients)

Lift an extreme coefficient `C` in monomial `C·x^α` to `z^k·x^α` with integer `k≥1`, anchor `z0 = |C|^(1/k)` (kept **exact** in WL), and **residual** `c = C/z0^k`. For real `C>0`, `c=1`; for real `C<0` or **complex `C`**, `c` is the (O(1)) phase/sign residual. Insert `∫dz δ(z−z0)`. Per-sector: run the (n+1)-dim sector machinery up to (not including) flattening, then resolve the delta by a pivot substitution (full formulas in `OLD_CODE/TROPICAL_MONTE_CARLO2/PLAN_SUM/plan.md` §3.2–3.4, and implemented in Tree B `ProcessSectorLifted` 849–1110). Carry forward verbatim:
- **Pivot ranking (amended, owner-approved):** `(1) HasConstantTerm=True, (2) |m_p|=1, (3) max min_j Re[ãtilde_j]`. Constant-term preservation is first because `HasConstantTerm=False` sectors can have **infinite MC variance** (§9 risk 1; bug-log L5).
- **Anchor selection (`kStar`):** `k* = max(1, ⌈|log|C|| / log τ⌉)` with `τ=1000` (the detector band). Decision order: **geometry gate first** (require a surviving sector with `HasConstantTerm=True` and all `Re[ãtilde]>0`; if none for any `k`, **do not lift**), then smallest `k` with `z0∈[1/τ,τ]`, then simpler geometry. Ship `DetectExtremeCoefficients`'s `SuggestedK -> kStar` (Tree B currently ships the placeholder `1`). Full procedure: `OLD_CODE/TROPICAL_MONTE_CARLO2/TROPICAL_MONTE_CARLOv2/AUXT/anchor_selection_procedure.md`.
- **BUG 2 fix (aux variable robustness):** build the auxiliary variable so it works even when base variables are *not* indexed (`x[i]`). Tree B / CALC2 used `auxVar = Head[vars[[1]]][n+1]`, which yields the invalid `Symbol[4]` (and a `Symbol::string` error) for plain-symbol specs. Use a guaranteed-fresh indexed symbol or normalize the spec's variables to indexed form, and emit a clear message otherwise. (`TMCv2_BUG_LOG.md` BUG 2.)
- **`liftidentity` check** must be a relative-tolerance comparison (not `Simplify[...]===0`) so it does not false-positive on floating `z0`. (`lift_error_log.md` pre-session fix.)

### 6.3 Complex EXPONENTS in codegen — the two fixes that made the bug-log branch correct

The **default** complex-exponent path is simply `result *= std::exp(B_j * std::log(P_j))` with complex `B_j` and complex `P_j` (already in Tree B). This is correct and is the v3 default (`ComplexExponentMode -> "Direct"`). **The opt-in `"SplitRealImag"` VEGAS-variance mode is where two real bugs lived** — implement it only with both fixes, and gate it with a complex-*base* cross-check:

- **Fix A (TMCv2_BUG_LOG BUG 1 — use the COMPLEX log).** The split identity is
  ```
  ∏_j P_j^{B_j} = ∏_j P_j^{Re(B_j)} · exp( i · Σ_j Im(B_j) · log P_j )     (log is the COMPLEX log)
  ```
  and `exp(i·Im(B)·log P) = exp(−Im(B)·arg P) · exp(i·Im(B)·log|P|)`. Using `std::log(std::abs(P_j))` (modulus) **silently drops the magnitude factor `exp(−Im(B)·arg P)`**, which is 1 only for real-positive `P_j`. Fix: emit `std::log(P_j)` (complex) and fold `cx(0, Im(B_j))` into the now-complex phase. *No-op for real-positive `P_j`* — which is exactly why every old lifting test missed it. **Therefore a v3 cross-check MUST use a complex-coefficient base (so `arg P ≠ 0`)** — see #21b.
- **Fix B (lift_error_log L3 — restore the tropical MONOMIAL factor).** `log P_k = log(monomial factor) + log Q_k`; tropical factoring removed `y^{d_k}` (and the lifted pivot substitution removed `z0^{…}·∏y^{…}`). The oscillatory phase must add back `Σ_k Im(B_k)·(Const_k + Σ_i Coeffs_{k,i}·log y'[i])`, where per sector & polynomial:
  - unlifted (`ProcessSector`): `Const=0`, `Coeffs_i = d_{k,i}/a_eff,i`;
  - lifted (`ProcessSectorLifted`): `Const = (d^aug_{k,p}/m_p)·log z0`, `Coeffs_j = [(d^aug_{k,rIdx[j]} − d^aug_{k,p}·m_{rIdx[j]}/m_p) + rcMin_{k,j}] / ãtilde_j`.
  Store this as a new SectorData key `"MonoFactorLog" -> {Const, Coeffs}`; codegen consumes it (with a safe fallback to the old form if absent). No effect on the real-exponent path.

### 6.4 High-D fan robustness — K-scaling (lift_error_log L2; reference impl already in TEST)

Automatic fan computation fails for ambient dimension ≥4 (lifted or not): thin lattice simplices (e.g. `conv{0, e_i}`) have no interior lattice point, so `ComputeDecomposition` leaks `$Failed` into the Polymake input. Fix: a private `computeFanScaled` that retries the normal-fan computation on a scaled copy `K·verts` with `K ∈ {1, n+2, 2n+4, 6n+6}` (the normal fan is scale-invariant). A reference implementation already exists in `OLD_CODE/TROPICAL_MONTE_CARLO/TEST/gen_sectors.wl` (`computeFanRobust`) — port it into `tropical_fan.wl` and wire it into the automatic-fan path of `EvaluateTropicalMCLifted`. After this, `liftdegenerate` fires only for a *genuinely* lower-dimensional polytope. (The 9D robust fan can take ~minutes — cache it and pass via `FanData` for repeated runs.)

### 6.5 High-D VEGAS sizing — `resolveVegasSizing` + `vegasbudget` (lift_error_log L1)

VEGAS defaults (`NStart=1000, NIncrease=500, NBatch=1000`) are tuned for n≲4 and return **confidently wrong** values (tight but lying error bars) for sector dimension ≳6 (measured: 8D unlifted 1.2% off at 7σ; lifted ~45% off). Fix: make those three options default to `Automatic`, resolved by `resolveVegasSizing[opt, maxSectorDim, base, kind]` — historical values for d≤4 (low-dim unchanged), `≈ base·3^(d−4)` for d≥5 (NStart ≈ 8.1e4 at d=8). Add `TropicalEval::vegasbudget`, fired when `NSamples < 20·NStart`. Residual guidance to document: in high-D, raise `NSamples` and **cross-check against a reference rather than trusting the VEGAS error bar**.

### 6.6 Harness hygiene (TMCv2_BUG_LOG harness note)
- Do **not** use over-broad `Check[expr, $Failed]` in test/driver harnesses — it caught benign WL messages (`General::munfl` on ~1e-204 sector magnitudes, then `Divide::infy`/`Greater::nord`) and reported a spurious `$Failed`. Capture the raw `Association` return; `Quiet`/guard tiny-magnitude messages in the package's WL post-processing.

### 6.7 Exactness discipline — purge decision-making numerics from the decomposition (enforces §1.4)

Audit every ported function and replace each float-based *decision* with an exact one. Known sites in the old engines:
- **Lifting pivot realness** (Tree B `ProcessSectorLifted` ~L963): the old filter `NumericQ[av] && Abs[Im[av]] < 10^-12 && Re[av] > 0` uses a `10^-12` tolerance to decide whether `atilde` is real-positive. Replace with an **exact** test — `Im[av]` is exactly zero (`PossibleZeroQ` / exact simplification, ε & kinematics symbolic) and `Re[av] > 0` decided exactly. A complex `atilde` then routes to `liftcomplex` deterministically, not by tolerance.
- **Divergence detection** (`IdentifyDivergences`, `FlattenSector`, IBP "convergent at ε=0" test): keep effective exponents exact, ε-expand symbolically, decide `Re[a] ≤ 0` from the exact leading coefficient. No `N[...]` on exponents.
- **`z0`, residuals, all lift data**: exact throughout (reinforce the existing discipline). `z0^k === C` and the round-trip `z→z0 ≡ original` must hold under exact simplification (verified symbolically, #40).
- **Prefactors / Jacobians**: `Abs[detM]` (exact integer), `Times@@effectiveAVals`, exponent division — all exact.
- **`DetectExtremeCoefficients`** is the *only* sanctioned numeric, and only as a heuristic *selector* of which coefficient to lift; compute magnitudes from exact `Abs[C]` when `C` is exact, and document that correctness is independent of the selection.
- **`MmaToC`** is the numericization boundary — exact rationals/radicals/complex become C++ literals here, and only here.
- **Ingest:** rationalize any inexact *input* coefficients at the boundary (e.g. `0.0001 -> 1/10000`, or carry as exact) so the decomposition runs on exact data.

Guard test (#43): for an exact-input spec, assert the symbolic decomposition output (sector exponents, prefactors, Laurent coefficients, lift data) is **free of inexact numbers** before the `MmaToC` boundary — i.e. `FreeQ[symbolicOutput, _Real]` — so any stray `N[...]` that crept into a decision path is caught mechanically.

---

## §7 Phased implementation plan (with gates)

> Discipline (from the proven v1→v2 plan): **order is law**; re-load the package after every edit (`wolframscript -code 'Get["tropical_eval.wl"]; Print["LOAD OK"]'`); keep `z0` and all lift data **exact** until the codegen/NIntegrate boundary (§1.4); every test prints `<NAME> PASS …` / `<NAME> FAIL expected=<e> got=<g>`; never edit anything under `OLD_CODE/`. **Every gate is evaluated by an independent verifier agent (§8.5) — never by the agent that produced the result — and numeric PASSes are accepted only against independent references.**

### Phase 0 — Setup & golden baselines
- **0.1** Create the v3 skeleton (§3.2). Copy `tropical_fan.wl` from the (identical) old copy.
- **0.2** Capture **golden baselines from BOTH old engines** (run read-only, in `OLD_CODE`, save transcripts under `TROPICAL_MONTE_CARLO3/TEST/baselines/`):
  - Tree A: `RunAllTests[]` (Tests 1–18) + IBP Tests 19–22 + `test_divergent_crosscheck.wl` → the divergence/IBP reference numbers (e.g. B=4 pole/finite = (1/6,−1/6); B=5 = (1/12,−1/8)).
  - Tree B: `run_validation_suite.wl` (Tests 1,2,3v2,5,6,7) + `test_lifted.wl` (Test 23) + `test_vegas.wl` (V1–V4, L1–L3) → the lifting/complex/VEGAS reference numbers (e.g. Test 6A = 7.850075e-4; lifted Case A VarRed ≈ 17.95×).
- **0.3** Smoke-test the Polymake plumbing in the new location (`PolytopeVertices` + `ComputeDecomposition` on `(1+x1^2+x2^2)^-1`).
- **Gate 0:** both baseline transcripts saved; fan smoke test passes.

### Phase 1 — Seed from Tree A; capture goldens
- **1.1** Copy Tree A `tropical_eval.wl` → v3. Confirm it loads and reproduces the Tree A golden (Tests 1–22), MC and VEGAS (per-kp **and** batched).
- **1.2** Capture v3 codegen goldens for every kind (conv / g0 / g1 / rem / ibp; MC, per-kp VEGAS, batched VEGAS) — the regression anchor for all later refactors.
- **Gate 1:** package loads; Tests 1–22 + `test_divergent_crosscheck` reproduce the Phase-0 Tree-A numbers; all goldens captured.

### Phase 2 — Refactor the shared core; unify codegen & drivers
- **2.1** Introduce Module 0 primitives as the single shared copy (they are already identical — pure reorganization).
- **2.2** Reconcile `ProcessSector` + `FlattenSector` with eps threading (§6.1). Verify the `IsDivergent` pre-flattening branch is preserved.
- **2.3** Unify the convergent/IBP codegen towers into one parametrized emitter (§5.2), **folding `emitVegasBatchMain` into the VEGAS emitter's `Batch` mode** (§5.1); unify the drivers (§5.4); de-duplicate validation + log-exp eval (§5.5); CUBA detection (§5.3); option defaults (§5.7).
- **Gate 2:** golden-master codegen regression green for **all** integrand kinds (conv, g0, g1, rem, ibp-boundary, ibp-term) under MC, per-kp VEGAS, **and batched VEGAS**; Tests 1–22 + crosscheck reproduce Phase-0 numbers **exactly**; `RunBenchmark2D`/`findCubaPrefix`/`emitVegasBatchMain`/IBP-codegen-twin symbols gone (folded or deleted; grep). This is the "simplify aggressively but prove nothing broke" gate.

### Phase 3 — Graft lifting (real coefficients), sandbox-first
- **3.1** Port Tree B's `SANDBOX/` toys (toy0 1D, toy1 eps, toy2 large-coeff) + `sandbox_lift_common.wl`. Run them read-only against the v3 core (this doubles as an integration proof of the §6.1 seam). **Hard gate 3a:** all toys PASS exactness (relerr < 1e-3 vs exact/NIntegrate); variance comparison recorded.
- **3.2** Port `DetectExtremeCoefficients`, `LiftCoefficients`, `ProcessSectorLifted`, `ValidateLiftedDecomposition` (Tree B M1b/1c). Apply the BUG 2 (aux-var) and `liftidentity` fixes (§6.2). Implement `SuggestedK -> kStar` + the geometry-gate anchor selection.
- **3.3** Wire `Lift` into the unified driver + the `EvaluateTropicalMCLifted` wrapper; codegen domain-indicator emission.
- **Gate 3:** Test 23 (A exactness <0.5% + DroppedSectors count; B end-to-end C++ <1%; C error paths fire cleanly; D EmptyDomain) reproduces Tree-B numbers; `bench_lift_variance` reproduces Case A VarRed ≈ 18× (and Case C correctly *fails* / does not lift). `HasConstantTerm` reported per sector.

### Phase 4 — Complex exponents end-to-end + complex lifting
- **4.1** Confirm the **default** complex path (`exp(B·log P)`, complex `A` and `B`) works through the unified codegen for Tree-B Examples 2,15,16,19 and exact-Gamma Example 21. (No SplitRealImag yet.)
- **4.2** Implement complex-*coefficient* lifting (anchor `z0=|C|^(1/k)`, residual `c=C/z0^k`) — mostly Tree B's residual mechanism; verify against a complex-coefficient toy.
- **4.3** Implement the opt-in `"SplitRealImag"` VEGAS mode with **both** §6.3 fixes (complex `log P` + `MonoFactorLog`). Add the `MonoFactorLog` SectorData key in both `ProcessSector` and `ProcessSectorLifted`.
- **Gate 4:** (a) complex exponents match exact Gamma/Beta (#2) and NIntegrate (#1) to old tolerances; (b) **the complex-BASE cross-check #21b PASSES** — a polynomial with complex coefficients (`arg P ≠ 0`) and complex `B`, lifted, where `SplitRealImag` and `Direct` agree and match a reference (this is the case that exposed BUG 1; without it the fix is unverified); (c) real-exponent codegen byte-identical (phase block only emitted when `Im(B)≠0`).

### Phase 5 — High-D hardening + anchor finalization
- **5.1** Port `computeFanRobust`→`computeFanScaled` into `tropical_fan.wl`; wire into the automatic lifted-fan path (§6.4).
- **5.2** Implement `resolveVegasSizing` + `vegasbudget` (§6.5).
- **Gate 5:** an 8D lifted-with-complex-coeff + VEGAS case (analogue of bug-log Test 26) converges with **default (Automatic)** sizing (≈0.1% vs reference) and `vegasbudget` fires when under-budgeted; low-dim (n=2) VEGAS unchanged (relerr ~1e-6).

### Phase 6 — Cross-check suite, benchmarks, docs
- **6.1** Assemble the full v3 cross-check suite (§8.2 **and** §8.3) under `TEST/`, with a single CI entry point (`run_validation_suite.wl` analogue) that exits non-zero on any FAIL. **Verify the §8.1 coverage matrix: every feature row must map to ≥1 implemented, passing check** — an empty row is a release blocker.
- **6.1b** Ship the new worked examples (§8.4) under `EXAMPLES/` (they double as cross-checks #8/#21b/#23/#33/#34/#35/#36/#41).
- **6.2** Port the FIESTA cross-checks (`CC/`) with path fixes (cross-check #8) — optional, gated on a FIESTA install; degrade gracefully if absent.
- **6.3** Rewrite `SUMMARY.txt` and merge the two `MANUAL/manual.tex` files (convergent + divergence + lifting + complex; carry the `complex_flatten_spiral` figure and the anchor-selection section).
- **Gate 6 (final):** the full suite is green on this machine; **every §8.1 feature row is covered by a passing check, and all new combination checks #31–#42 pass**; CUBA-dependent checks degrade gracefully when CUBA is absent; FIESTA checks skip cleanly when FIESTA is absent; docs build.

---

## §8 The comprehensive cross-check suite

The "many cross checks" requirement, in four parts: **§8.1** a feature-coverage matrix proving *every* feature is cross-checked (nothing slips through the merge); **§8.2** the cross-checks ported/adapted from the two old trees; **§8.3** *new* cross-checks v3 introduces — most testing feature **combinations** that were impossible before the merge (the features lived in separate codebases) plus the merge itself; **§8.4** new worked examples. Tier 1 = no external deps (WL + g++); Tier 2 = needs CUBA; Tier 3 = needs FIESTA. The CI gate requires all Tier-1 green; Tier-2/3 skip cleanly when their dependency is absent.

### §8.1 Feature-coverage matrix (no feature uncovered — an empty row is a release blocker)

| Feature / sub-feature | Cross-checks covering it |
|---|---|
| F1 tropical decomposition (fan, ProcessSector, flatten) | #1, #2, #25, #26, #29, #30, #37, #38 |
| F2a divergence — **IBP (default)** | #3, #4, #13, #14, #16, #27, #33, #35 |
| F2b divergence — tropical subtraction | #3, #4, #15, #16, #35 |
| F2c nested-divergence refusal (1/ε^{d≥2}) | #27 |
| F3a uplift — real coefficients | #17, #18, #19, #20, #28, #39, #40 |
| F3b uplift — complex coefficients | #34, #40 |
| F3c uplift — anchor / `k` selection | #20, #39 |
| F4a complex polynomial exponents B | #2, #21, #21b, #33, #34 |
| F4b complex monomial exponents A | #2, #21 |
| F4c SplitRealImag VEGAS mode (the BUG-1 regime) | #21b, #34 |
| F5a VEGAS (per-kp) | #5, #6, #7, #9, #42 |
| F5b VEGAS **batched** | #23, #36 |
| F5c high-D VEGAS sizing | #42, #34 |
| MC back-end | #1, #9, #22, #36, #39 |
| **the merge itself** (regression vs both old engines) | #31, #32 |
| integrand form: numerator factors (B>0) | #1, #2, #41 |
| **exactness invariant** (no decomposition numerics) | #43, §6.7; verified *exactly* by §8.5 layer 1 |
| **unattended numeric trustworthiness** | §8.5 independent adversarial verifier over #1, #2, #4, #5, #6, #9, #21b, #22 |
| external validation | #8 (FIESTA), #2 (exact Γ/Β), #5 (CUBA Cuhre) |

### §8.2 Cross-checks ported / adapted from the old trees

| # | Cross-check | Tier | PASS criterion (representative) |
|---|---|---|---|
| 1 | Tropical sector sum vs direct `NIntegrate` | 1 | relerr < 1e-3 (tighter for n≤4); the workhorse, every example |
| 2 | vs exact Gamma/Beta/Dirichlet closed form | 1 | e.g. Beta-integral Example 21 to ~1e-6 |
| 3 | Divergent (pole, finite) vs exact Laurent (Γ products) | 1 | B=4→(1/6,−1/6), B=5→(1/12,−1/8); tol ~5e-3 |
| 4 | **IBP (default) vs subtraction** agree on (pole, finite) | 1 | ≲1e-3 — subtraction validates the default IBP result (the reason both are kept, D2) |
| 5 | vs CUBA Cuhre on the direct (un-decomposed) integrand | 2 | sector sum within 1–5%; MC within 5σ |
| 6 | vs CUBA VEGAS (independent sampler) | 2 | agreement within combined error |
| 7 | CUBA on raw vs on decomposed+flattened integrand | 2 | decomposition decisively better (report ratio) |
| 8 | vs **FIESTA5** (external sector-decomposition tool) via GL(1)/Cheng–Wu | 3 | bubblepp pole exact, finite within ~2σ |
| 9 | **MC vs VEGAS parity** | 2 | `|MC−VEGAS| ≤ 5·√(errMC²+errVEGAS²)` |
| 13 | Pole coefficient vs numerical ε-slope estimate | 1 | relerr < 5% |
| 14 | G0 = B⁰ − Σ coeff·I (IBP ↔ subtraction pole identity) | 1 | the two pole-extraction routes share this analytic identity |
| 15 | ε/ε cross-term resolution (finite = TOTAL − G0 + γ·G0) | 1 | reproduces the `CC/diagnose_eps_eps` correction (<0.01%) |
| 16 | Subtracted-remainder convergence as ε→0 | 1 | remainder → 0 with documented rate |
| 17 | **Lifted vs unlifted exactness** | 1 | relerr < 0.5% vs exact; DroppedSectors count asserted |
| 18 | **Lifted vs unlifted variance** (VarRed = (σ_U/σ_L)²) | 1/2 | Case A ≥ 10× (≈18×); Case C correctly *not* lifted |
| 19 | Sample σ vs **true σ** (exact `NIntegrate` σ²=I2−I1²) | 1 | flags `HasConstantTerm=False` infinite-variance sectors |
| 20 | Anchor (integer `k`) sweep vs the `k*` rule | 1 | argmin relErr ≥ k* in sensitive cases (H1/H2/H3) |
| 28 | Lifted-sector domain / feasibility checks | 1 | EmptyDomain drop count, FeasibleFraction, HasConstantTerm tally all as asserted |
| 21 | **Complex-flattening** variance reduction (contour rotation) | 1 | sector value identical on both contours; ~20× variance drop |
| 21b | **Complex-BASE + complex-exponent lifting** (`arg P ≠ 0`) — *new, gates the BUG 1 fix* | 1/2 | `SplitRealImag` ≡ `Direct` ≡ reference; the case real-positive tests cannot catch (§6.3) |
| 22 | Determinism / seed reproducibility (`SeedBase` argv) | 1 | same seed → byte-identical; distinct seeds → distinct streams, each within 5σ |
| 23 | **Batched VEGAS == per-kp VEGAS** (batching invariance) | 2 | batched value within combined error of per-kp; both match the closed form (gates the §5.1 fold) |
| 25 | **Golden-master codegen regression** (byte-identical emitted C++) | 1 | conv / g0 / g1 / rem / ibp-boundary / ibp-term, MC & VEGAS — diff clean |
| 26 | GL(1)/Cheng–Wu affine↔simplex identity at finite ε | 1 | numeric equality of the chart change (validates #8) |
| 27 | **Known-limitation enforcement → graceful `$Failed`** | 1 | nested 1/ε^{d≥2} (`nested`/`nestedIBP`), missing CUBA (`nocuba`), bad integrand head (`badcpp`), no pivot (`liftnopivot`), complex-`atilde` (`liftcomplex`), degenerate polytope (`liftdegenerate`), under-budget VEGAS (`vegasbudget`) all fire cleanly, never silently wrong |
| 29 | Flattened-integrand magnitude spot-check | 1 | warns when max>1e3 or min<1e-6 (heavy-tail indicator) |
| 30 | Validation utilities exercised as assertions | 1 | `ValidateDecomposition`/`Subtraction`/`IBP`/`Lifted` themselves pass on known inputs |

**Highest-value, lowest-friction (port first):** 1, 2, 3, 4, 17, 19, 21b, 22, 25, 27. **Strongest independent validation:** 8 (FIESTA), 2 (exact Gamma), 5 (CUBA Cuhre).

### §8.3 NEW cross-checks introduced by v3 (combinations the split codebases could not test)

The genuinely new validation ideas (numbers ≥31 mark them v3-original). The theme: each old engine tested its own features in isolation; v3 can — and must — test the feature **crossings** and the **merge itself**. These are the cross-checks the user asked v3 to invent.

| # | New cross-check | Tier | Why it is new / what it catches | PASS criterion |
|---|---|---|---|---|
| 31 | **Merge fidelity vs Tree A** | 1 | the merge must not perturb the divergence engine | every Tree-A golden (Tests 1–22, divergent crosscheck) reproduced within MC noise (bitwise for symbolic) |
| 32 | **Merge fidelity vs Tree B** | 1 | …nor the lifting/complex engine | every Tree-B golden (Tests 1,2,3v2,5,6,7; Test 23; V/L VEGAS) reproduced |
| 33 | **Divergence × complex exponents, full matrix** | 1/2 | F2×F4×F5 — Tree A's crosscheck was real-exponent, Tree B had no divergence | one ε-integral with complex B through {IBP, subtraction}×{MC, VEGAS}: pairwise (pole,finite) ≲1e-2, vs exact Γ-product Laurent ≲5e-3 |
| 34 | **Lifting × complex coeff × VEGAS, high-D (≈8D)** | 2 | F3×F4×F5 in the exact regime that produced BUG-1 / L1 / L2 / L3 | `Direct` ≡ `SplitRealImag` ≡ reference ≲0.5%; K-scaled fan builds; `vegasbudget` fires if under-budgeted |
| 35 | **Laurent via three independent routes** | 1 | triangulates F2 | pole/finite from (a) IBP@ε=0, (b) subtraction@ε=0, (c) finite-ε `LaurentFromSubtraction` fit all agree ≲1e-2 |
| 36 | **Sampler × method × batch consistency** | 2 | exercises the §5.1 batch fold + kinematic plumbing together | one kinematic-scan integral: per-kp MC, per-kp VEGAS, batched VEGAS, and standalone one-point runs all agree within combined error; one-point == in-scan value |
| 37 | **Fan scale-invariance** | 1 | promotes a one-off "verified-correct" note to a gate; catches ray/K-scaling bugs (§6.4) | primitivized vs raw rays and K∈{1,n+2,…} fans → byte-identical generated C++ and bitwise-equal sector sum |
| 38 | **Variable-relabeling / regulator-position invariance** | 1 | catches index-ordering bugs across the merged code (generalizes Tree A Test 10) | permuting variables / moving the ε-regulator leaves the answer invariant to MC noise |
| 39 | **"The error bar lies" guardrail** | 1/2 | encodes the L1 / L5 / Example-20 heavy-tail lesson as an automated regression | a small-coeff case where unlifted plain-MC is confidently wrong: lifting/VEGAS/reference recover the truth AND a warning fires (HasConstantTerm / magnitude / vegasbudget) |
| 40 | **Lift identity — complex & negative & multi-rule (symbolic)** | 1 | hardens the residual mechanism for complex/negative C (BUG-2 regime) | z→z0 substitution ≡ original polynomials, exact (rel-tol for floats), for +, −, complex, and 2-rule incommensurate lifts |
| 41 | **Numerator (B>0) × {divergence, lifting}** | 1 | the old numerator examples were convergent-only | numerator factor + ε-pole, and numerator + lifting, vs NIntegrate/exact ≲1% |
| 42 | **High-D VEGAS sizing (`resolveVegasSizing`)** | 2 | makes the Phase-5 sizing fix a named maintained check | 8D convergent integrand: Automatic sizing reaches ~0.1%; low-D (n=2) unchanged (~1e-6); `vegasbudget` fires at `NSamples<20·NStart` |
| 43 | **Exactness guard** — no decision-making numerics in the decomposition | 1 | enforces §1.4 / §6.7 | symbolic decomposition output is `FreeQ[_, _Real]` before the `MmaToC` boundary; any stray `N[...]` in a decision path is flagged |

**Build-first among the new checks:** #31, #32 (catch any merge regression immediately), then #33, #34, #40 (the feature crossings that are the whole point of v3).

### §8.4 NEW worked examples (EXAMPLES/, doubling as the manual's worked-example chapter)

A small "capstone" set demonstrating the combinations only v3 supports; each also *is* one of the cross-checks above:
- **Ex-V3a — divergent + complex, full Laurent:** a 1-loop massive bubble (or the `CC/` bubblepp presector) with an ε-pole and complex exponents; (pole, finite) via IBP, cross-checked vs subtraction, exact Γ, and FIESTA (#8, #33, #35).
- **Ex-V3b — small complex coefficient, lifted + VEGAS:** an integrand with a tiny *complex* coefficient and complex B; lift via `SplitRealImag`, integrate with VEGAS, cross-check vs CUBA Cuhre on the direct integrand (#34, #21b).
- **Ex-V3c — kinematic scan, batched VEGAS:** one compile, many kinematic points, `Batch -> True`; cross-check per-point vs NIntegrate and vs a per-kp VEGAS run (#23, #36).
- **Ex-V3d — numerator + pole:** a genuine numerator `(…)^{+3/2}` times a divergent denominator, exercising B>0 alongside F2 (#41).

### §8.5 Independent adversarial verification — how every numeric PASS is gated (the user's explicit requirement)

The user cannot eyeball a subtly-wrong numeric value, and will not: **no gate is judged by the agent that produced the result, and no numeric PASS is ever accepted on a single value or a self-reported error bar.** For each result, the unattended per-phase workflow (§7) spawns a **separate verifier agent** — fresh context, did not produce the result — that is prompted to **refute** it. The verification is two-layer, exploiting the exactness invariant (§1.4):

- **Layer 1 — exact symbolic verification (the bulk).** Because the decomposition is exact, the verifier checks it with identities that *cannot be subtly wrong*, with **zero** floating point: the lift round-trip `z→z0 ≡ original` (#40), the tropical-factoring algebraic identity, `z0^k === C`, `det M` and exponent bookkeeping, the IBP↔subtraction pole identity (#14), and the GL(1) chart identity (#26). Any symbolic mismatch is a hard FAIL — and it is *decisive*, not statistical.
- **Layer 2 — numeric verification, only at the integration step, against ≥2 INDEPENDENT oracles.** The single genuinely-numeric output (the integrated sector value) is compared against references that share no code or algorithm with the producer: exact Γ/Β closed forms (#2), CUBA Cuhre on the *direct* integrand (#5), CUBA Vegas (#6), a second method/sampler (IBP vs subtraction #4, SplitRealImag vs Direct #21b, MC vs VEGAS #9), reseeded reproducibility (#22), and `NIntegrate` (#1). The verifier is instructed to *assume the value is wrong* and hunt for a reference that disagrees beyond combined error.
- **Fail-closed gate (no voting).** A result is accepted only if **every** verifier independently confirms it (`refuted=false`); **any** verifier that refutes — finds it wrong *or* cannot corroborate a numeric value — halts the run (no majority, no "most passed"). The complex-lifting phases (4–5) — the L1/L5 danger zone — run a **3-lens panel** (A: exact symbolic; B: an independent integrator; C: reseed/budget stability) and require **all three** to confirm; every disagreeing reference is surfaced in the phase report. (Orchestration realization + the orchestrator-implements-nothing rule: `ORCHESTRATION.md` §1.7.)

This is precisely what makes "leave each phase unattended" safe **without eyeballing**: correctness is decided by exact identities plus independent numeric oracles plus an adversarial skeptic — never by a human reading numbers, and never by the producer grading its own work. In the per-phase workflow scripts this is the canonical `producer agent() → independent verifier agent()/panel` pattern; a phase gate is GREEN only when **every** verifier signs off, and **any** refutation or non-corroboration halts the phase with diagnostics (fail-closed; no majority vote). (Cross-references: §1.4 exactness invariant, §6.7 audit, §7 gates, §10 execution.)

---

## §9 Edge cases, risks, known limitations

| # | Item | Handling |
|---|---|---|
| 1 | **`HasConstantTerm=False` → infinite MC variance** (lift_error_log L5) | The single biggest risk to lifting's *benefit* (not correctness). Pivot ranking prefers constant-term preservation; the correctness gate is the **exact** `ValidateLiftedDecomposition`, NOT the sampled 5σ; report `HasConstantTerm` and cross-check against a reference; never trust the sampled error bar for such sectors. |
| 1b | **Subtly-wrong numeric PASS slips an unattended gate** (the user's core concern — cannot be eyeballed) | Structurally closed: §1.4 keeps the decomposition exact (verified by exact identities, not numbers), and §8.5 gates the one numeric output with an **independent adversarial verifier** against ≥2 independent oracles, plus a 3-lens panel on the complex-lifting phases 4–5. No human eyeballing required; the producer never grades its own work. |
| 2 | **Over-zealous deletion** breaks lifting | `ProcessSector`'s `IsDivergent` pre-flattening branch is on the KEEP list (§6.1); Gate 2 + the Phase-3 toy rewiring catch a violation. |
| 3 | **Complex-base SplitRealImag** (the BUG 1 regime) | Implement with complex `log P` + `MonoFactorLog` (§6.3); gate with #21b (complex coefficients). Default `"Direct"` mode is always safe. |
| 4 | **Nested divergence (1/ε^{d≥2})** | Refuse cleanly (`nested`/`nestedIBP`); guarded in 3 places in Tree A (`IdentifyDivergences` 635–639; `IBPProcessSector` 2999–3006 & 3019–3024) — preserve all three. Full integral at fixed ε via `EpsilonValue`. |
| 5 | **High-D fan failure** (dim≥4) | `computeFanScaled` K-scaling (§6.4); reference impl in `TEST/gen_sectors.wl`. |
| 6 | **High-D VEGAS silent-wrong** | `resolveVegasSizing` + `vegasbudget` (§6.5). |
| 7 | **Kinematic-symbol coefficients can't auto-lift** | `DetectExtremeCoefficients` skips them (magnitude unknown); manual lift rules still allowed. `ãtilde` is kinematics-independent, so one compile serves all points. |
| 8 | **Lift + 1/ε pole** (N3) | Supported (Case A, planAXpDIV.md): lifted sector + simple pole via PrefactorBase/DomainConstraint, both routes. **Complex B supported via SplitRealImag** on both routes (planCXLIFTDIV.md, cc_46). Case B / inline-symbolic-ε SplitRealImag×lift×div / complex Case-B / Direct-mode complex-B-lift refused ($Failed). |
| 9 | **Harness false-`$Failed`** | No over-broad `Check`; `Quiet` tiny-magnitude messages (§6.6). |
| 10 | **CForm fallback leaks bad tokens** | `badcpp` gate + targeted test #27 (§5.8). |

**Stated limitations to document in SUMMARY/MANUAL:** single divergent variable per sector; lifting requires the lifted sectors to keep an O(1) constant term to *benefit*; complex-exponent lifting uses `SplitRealImag` (variance built on `Re(B)`, so very large `Im(B)` can inflate variance — cross-check against a reference); FIESTA/CC scripts need a FIESTA install + path edits.

---

## §10 Execution conventions

- **Run scripts:** `cd /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3 && wolframscript -file <script>`; scripts `SetDirectory`/use absolute paths when `Get`-ing packages.
- **Shared `INTERFILES/` hazard (from Tree B SUMMARY §5):** the C++ MC writes fixed filenames; run MC-executing scripts **one at a time** in a given `WorkingDirectory`, or give each its own. The symbolic, `NIntegrate`-gated tests are unaffected.
- **Toolchain (this machine):** `wolframscript` (Mathematica.app), `polymake` (`/opt/homebrew/bin`, auto-detected by `tropical_fan.wl`), `g++` (CompileCpp retries without `-fopenmp` automatically), CUBA (`brew install cuba`; `/opt/homebrew`). The brew `libcuba.a` macOS-version link warnings are harmless.
- **Numerics discipline:** keep `z0` and lift data exact through the symbolic pipeline; numericize only at the codegen/NIntegrate boundary — keeps the `z→z0` round-trip identity and the exactness gates sharp.
- **Order is law:** Phase 0 → Gate 0 → Phase 1 → … → Phase 6, each gate green before the next phase. Never edit `OLD_CODE/`. Never edit `tropical_fan.wl` except the additive `computeFanScaled`.

---

## Appendix A — `OLD_CODE` file inventory (what to mine vs ignore)

**Tree A (`OLD_CODE/TROPICAL_MONTE_CARLO/`) — divergence/IBP/VEGAS/batch base**
- `tropical_eval.wl` (4584) — engine base for v3. `tropical_fan.wl` (468) — identical everywhere.
- `SUMMARY.txt` (819), `MANUAL/manual.tex` (1192) — v3 docs source (divergence half).
- `TMCv2_BUG_LOG.md`, `lift_error_log.md` — **the complex-lifting spec** (§2.2, §6).
- `EXAMPLES/` — `tropical_eval_examples.wl` (Tests 1–22), `test_divergent_crosscheck.wl` (#4), `test_coverage.wl`, `test_nested_divergence.wl` (#27), `test_vegas.wl`, `test_mc_drivers.wl`, `vegas_*` — mine for cross-checks 1–16,22,25,27.
- `TEST/` — benchmark study (N2: archive, don't maintain) BUT `gen_sectors.wl` has `computeFanRobust` (§6.4) and `golden_capture.wl` is the codegen-golden pattern (#25). `CC/` — FIESTA (#8).
- `VEGAS_IMPLEMENTATION_PLAN.md`, `DIVERGENT_VALIDATION.md`, `plan_test.md` — background.

**Tree B (`OLD_CODE/TROPICAL_MONTE_CARLO2/TROPICAL_MONTE_CARLOv2/`) — lifting/complex/VEGAS source**
- `tropical_eval.wl` (2666) — source for `FlattenSector`, lifting M1b/1c, complex path, `SeedBase`.
- `SUMMARY.txt` (1117) — definitive lifting (§8) + complex docs. `MANUAL/manual.tex` (1992) + `complex_flatten_spiral.pdf` — complex-flattening §6.20.
- `EXAMPLES/` — `tropical_eval_examples2/3.wl` (complex + small-coeff Examples 15–23), `test_lifted.wl` (Test 23, #17), `test_cuba.wl` (#5), `test_vegas.wl` (V/L cases, #9,#18), `test_seedbase.wl` (#22), `cuba_common.wl` (#5/#6/#7 generator).
- `SANDBOX/` — lifting dev: `sandbox_lift_common.wl`, `sandbox_toy0/1/2`, `sandbox_true_sigma.wl` (#19), `bench_lift_variance.wl` (#18), `z0_sweep_*` (#20), `verify_complex_flatten.wl` (#21).
- `AUXT/` — `anchor_selection_procedure.md` (the `kStar` spec, §6.2), `z0_*` plans/results.

**Tree B-inner (`OLD_CODE/TROPICAL_MONTE_CARLO2/TROPICAL_MONTE_CARLO/`)** — REDUNDANT (§2.1); ignore.

## Appendix B — old → v3 function migration

| v3 function | Source | Action |
|---|---|---|
| `ParsePolynomial`, `TransformExponents`, `MmaToC`, `GenerateMonomialSumCpp` | A≡B (identical) | shared M0, as-is |
| `ProcessSector` | A (eps) + B (`FlattenSector` split) | merge (§6.1) |
| `FlattenSector` | B | + optional `eps` arg |
| `CheckFlatteningMagnitude`, `ValidateDecomposition` | A/B | de-dup, extract `evalFlattenedIntegrand` (§5.5) |
| `DetectExtremeCoefficients`, `LiftCoefficients`, `ProcessSectorLifted`, `ValidateLiftedDecomposition` | B | port + BUG2/identity/kStar fixes (§6.2) |
| `IdentifyDivergences`, `ProcessDivergentSector`, `ValidateSubtraction` | A | as-is |
| `IBPReduceSector`, `IBPProcessSector`, `IBPCheckBoundary`, `ValidateIBP` | A | as-is |
| `GenerateCpp` | A `GenerateCppMonteCarlo` + `…IBP` | **unify** (§5.2); batched VEGAS retained as a mode (§5.1) |
| `emitIntegrandDefinitions`, `emitMain[MC|VEGAS]` (VEGAS has a `Batch` mode) | A's 6 emit fns incl. `emitVegasBatchMain` | **collapse to 2** (§5.2) |
| `CompileCpp`, `detectCuba` | A | keep; delete B's `findCubaPrefix` (§5.3) |
| `resolveVegasSizing`, `computeFanScaled`, `MonoFactorLog` | bug logs + `gen_sectors.wl` | **reconstruct** (§6.3–6.5) |
| `EvaluateTropicalMC` | A + …IBP + B …Lifted | **unify** + wrappers (§5.4) |
| `RunBenchmark2D`, `findCubaPrefix`, `GenerateCppMonteCarloIBP`, `EvaluateTropicalMCIBP` | — | **delete / absorb** (`emitVegasBatchMain` is *folded*, not deleted — §5.1) |
