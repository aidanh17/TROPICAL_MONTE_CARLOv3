# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A Wolfram Language pipeline for numerically evaluating generalized Euler / Feynman integrals
`I = ∫_{[0,∞)^n} ∏_i x_i^{A_i} ∏_j P_j(x)^{B_j} dx` (exponents and polynomial coefficients may be
real or **complex**, and may depend on a regulator `eps`). The strategy: a **tropical fan
decomposition** partitions the domain into simplicial sectors, a symbolic **flattening** change of
variables renders each sector integrand O(1) on `[0,1]^n`, and the result is emitted as **C++** and
integrated by Monte Carlo / VEGAS. A single C++ compilation serves thousands of kinematic points.

The whole engine is **two files**:
- `tropical_fan.wl` — Newton-polytope → normal fan → simplicial triangulation, via **Polymake**.
  `computeFanScaled` adds K-scaling retries for high-dimensional lifted fans.
- `tropical_eval.wl` — everything else (~7970 lines, one `TropicalEval`` ` package). Internal layout
  (search these comment banners): **M0** primitives, **M1** `ProcessSector`/`FlattenSector`,
  **M1b** lifting, **M2** divergence subtraction, **M2b** IBP pole extraction, **M3** C++ codegen,
  **M4** unified driver, **M5** validation.

`SUMMARY.txt` (concise API + limitations) and `MANUAL/manual.tex` (full reference, 32pp) are the
authoritative prose; read them before changing decomposition/codegen behavior.

## Commands

`wolframscript` is on PATH (Mathematica 12+). All scripts `Get` both packages relative to their own
location, so run them from the repo root.

```bash
# A named cross-check (prints "CCnn PASS"/"CCnn FAIL" lines; these do NOT exit nonzero — grep output)
wolframscript -file TEST/cc_47.wl 2>&1 | grep -E "CC47 (PASS|FAIL)"

# Self-gate suites (these DO Exit[1] on any FAIL — usable as CI gates)
wolframscript -file TEST/phase2_selfgate.wl   # codegen goldens + exactness guards + complex-A/B
wolframscript -file TEST/phase3_selfgate.wl

# Rebuild the manual (run pdflatex twice for cross-refs)
cd MANUAL && pdflatex -interaction=nonstopmode manual.tex && pdflatex -interaction=nonstopmode manual.tex
```

There is no aggregate test runner. Tests are the `TEST/cc_*.wl` cross-checks (numbered ~1–48) and
the `TEST/phase*_selfgate.wl` gates. `EXAMPLES/` holds four narrated worked examples (`Ex-V3a`–`Ex-V3d.wl`,
each keyed to `OLD_PLANS/plan.md` §8.4 and specific cross-checks). **Note:** `MANUAL/manual.tex`'s "Command-Line
Execution" section instead lists `EXAMPLES/run_validation_suite.wl` / `EXAMPLES/test_*.wl`, which **do
not exist** — the real gate tests live in `TEST/`.

**Dependencies:** Polymake (required, on PATH or `/opt/homebrew`,`/usr/local`,`/usr`); g++ with
C++17 (required for any codegen test); **CUBA** (optional — only for `"Integrator" -> "VEGAS"`;
absent ⇒ plain MC still works, VEGAS fails fast with `TropicalEval::nocuba`); FIESTA5 (optional, only
for `CC/` scripts, which skip gracefully if absent).

## The unified driver and how features compose

`EvaluateTropicalMC[spec, fanData, kinematicPoints, opts]` is the entry point. `spec` is an
`IntegrandSpec` association: `<|"Polynomials"->{...}, "MonomialExponents"->{A_i},
"PolynomialExponents"->{B_j}, "Variables"->{x[1],...}, "KinematicSymbols"->{...},
"RegulatorSymbol"->eps|>`. Three orthogonal feature axes compose through it:

- **Divergence** (`1/eps` poles). `Method -> Automatic` (default) routes divergent integrals to the
  **IBP** path (`EvaluateTropicalMCIBP`); `"Subtraction"`/`"None"` use the inline tropical-subtraction
  path; `LaurentFromSubtraction` is the pinned-`eps` route. IBP and Subtraction are independent
  algorithms whose agreement is a primary correctness signal. Only **one** divergent variable per
  sector is supported (nested ⇒ clean `$Failed`).
- **Lifting** (`LiftData` from `LiftCoefficients`). An extreme polynomial coefficient is replaced by an
  auxiliary variable `z` with anchor `z_0=|C|^{1/k}` and a δ-constraint; the `(n+1)`-dim lifted fan
  adapts to it. Lifting composes with a `1/eps` pole (Case A; `ProcessSectorLifted` is eps-aware).
- **Complex exponents** (`ComplexExponentMode`). `"Direct"` (default, always correct) uses the
  complex `std::log`. `"SplitRealImag"` (opt-in, for VEGAS variance) decomposes on `Re(B)`/`Re(A)` and
  re-attaches the imaginary parts as oscillatory phases (`MonoFactorLog` for `Im(B)`,
  `MonomialPhaseLog` for `Im(A)`). See `manual.tex` §"Auxiliary Lifting with Complex Exponents".

## Invariants you must preserve (these gate every change)

These are project-wide contracts, enforced by the test suite and stated across `OLD_PLANS/plan.md` §1.4 /
`SUMMARY.txt` / the `OLD_PLANS/plan*.md` feature notes. Violating one is a regression even if numbers look right.

1. **Exactness.** The decomposition is symbolic-exact: **no decision inside it is made by evaluating
   a float.** Floats appear only in (a) the C++ integration and (b) numerical cross-checks. Tests
   assert `FreeQ[<symbolic payload>, _Real]` at the `MmaToC` boundary — keep exact rationals/radicals
   exact (e.g. `-1/2`, `Log[100]`), never introduce a machine float upstream of codegen.
2. **Byte-identity (#25).** The **real-exponent / `Direct` path must emit byte-identical C++** to the
   committed goldens in `TEST/baselines/codegen_goldens/`. Complex/`SplitRealImag` code only activates
   when an imaginary part is actually present, so it must not perturb the real path. After touching
   `MmaToC`/codegen, re-run `phase2_selfgate.wl` and confirm the goldens still match.
3. **Known limitation ⇒ clean `$Failed`, never a wrong number.** Unsupported cases must `Message`
   a `TropicalEval::*` and return `$Failed`. The message table is in `manual.tex` §"Error Messages";
   examples added recently: `splitliftdiv`, `splitdivmono`, `liftcomplex`, `liftdivdomain`.
4. **Every numeric PASS needs ≥2 independent oracles** (NIntegrate of the original, CUBA, IBP↔
   Subtraction agreement, `SplitRealImag`↔`Direct`, reseed). A single self-consistent number is not
   a pass.
5. **NEVER edit `OLD_CODE/`.** It holds the two parent trees (`TROPICAL_MONTE_CARLO`,
   `TROPICAL_MONTE_CARLO2`) that v3 merges; it is reference-only.

## Where the design is written down

- `OLD_PLANS/plan.md` — what v3 is, the merge of the two parent trees, per-phase scope/gates, design
  decisions (Dn) and non-goals (Nn).
- `OLD_PLANS/planAXpDIV.md` (lift + divergence), `OLD_PLANS/planCXLIFTDIV.md` (complex-`B` × lift ×
  divergence), `OLD_PLANS/planAXpDIVv2.md` (complex-**monomial** `A` × lift), `OLD_PLANS/planR.md`
  (review-remediation) — each is the spec for one feature; the matching `TEST/cc_4x.wl` is its
  cross-check. `OLD_PLANS/` also holds `planIBPCX.md`, `planIBPMULTIDIV.md`, and `planLIFTDIVZ0.md`
  (later feature specs, same convention).
- `ORCHESTRATION.md` + `ORCH/*.js` — how v3 was *built* (per-phase unattended multi-agent workflows
  with adversarial verifiers). Not needed for ordinary development.
- `fable_inventory.md` / `fable_refactor_plan.md` — the in-progress internal-cleanup refactor of
  `tropical_eval.wl` (dead-code removal, dedup, decomposition); behavior-preserving, gated by the
  same self-gate suites. Move to `OLD_PLANS/` once complete.
