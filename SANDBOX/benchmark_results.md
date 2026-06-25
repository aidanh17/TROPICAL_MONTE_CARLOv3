# TROPICAL_MONTE_CARLO v3: Lift-Variance Benchmark (HONEST exact-trueSigma)

**Date:** 2026-06-25  
**Method:** EXACT per-sample sigma per scheme, trueSigma^2 = I2 - I1^2 with
I1 = Int g and I2 = Int g^2 by converged NIntegrate (PrecisionGoal 4).
Total estimator variance = Sum_sectors trueSigma_s^2; INFINITE the moment ANY
sector has a divergent I2.  HasConstantTerm=False is the usual structural ROOT
CAUSE of a divergent I2 (re-clearing loses the constant term), but it is NOT a
perfect proxy: a HasConstantTerm=True sector can STILL have a divergent I2 (see
Case B sector 7).  The gate therefore keys on the EXACT I2 convergence and
raises a warning on every divergent scheme.  Sampled sigma (sigma = ReErr*Sqrt[N])
is shown for transparency only and is NEVER used as a pass criterion -- for a
sector with divergent I2 it is an optimistic (often wildly underestimated) proxy.

The CORRECTNESS gate is the EXACT `ValidateLiftedDecomposition` (Test 23 /
phase3_selfgate), NOT this variance number (plan.md Â§9 risk #1).

---

## Case A: `P = 1 + 10^6 x1^2 + x2^2 + x1 x2^2`, B={-2}  (lift {2,0}, k=2)

| Scheme | I1 sum | converged ref | I1 relErr | true sigma | classification |
|--------|--------|---------------|-----------|------------|----------------|
| Unlifted | 0.0007849842577801977 | 0.0007850074736127726`5. | 0.000029574027452369783 | 0.01658965697013982 | finite |
| Lifted k=2 | 0.0007850027299395547 | 0.0007850074736127726`5. | 6.042838287989023*^-6 | **INFINITE** | **class (ii)** |

Lifted HasConstantTerm per surviving sector: {True, True, True, True, True, False}  (z0 = 1000).
The HasConstantTerm=False sector (fan index {8})
has DIVERGENT I2, so the lifted estimator has **INFINITE true variance**.
**EXACT VarRed = INFINITE-VARIANCE / NOT-BENEFICIAL** (plan Â§9 risk #1, bug-log L5).

> The earlier "VarRed ~ 18x PASS" was computed from SAMPLED sigma (unlifted 0.016677064019884672 / lifted 0.003986604341149474); it is an
> artifact of the finite MC missing the heavy tail of the HasConstantTerm=False
> sector (feasFrac -> 0).  The HONEST result is infinite variance.

---

## Case B (report-only): mixed-coefficient 5-monomial, B={-2}  (1-rule k=2 {2,0})

| Scheme | I1 sum | converged ref | true sigma | EXACT VarRed |
|--------|--------|---------------|------------|--------------|
| Unlifted | 144.05793609582355 | 144.0578434006272533623`5. | INFINITE | - |
| Lifted k=2 | 143.38043296009027 | 144.0578434006272533623`5. | INFINITE | INFINITE |

Lifted HasConstantTerm per sector: {True, True, True, True, True, True}.
Despite ALL lifted sectors being HasConstantTerm=True, sector {7} has a DIVERGENT I2, so the lifted estimator has INFINITE true variance and the VarRed is reported INFINITE -- a sharper, more honest result than HasConstantTerm alone would suggest.
Case B is report-only (no pass gate).

---

## Case C: `P = 1 + 10^8 x1^3 x2 + x2^3`, B={-3}  (degenerate lifted polytope)

Automatic anchor selection (geometry gate) **REFUSES to lift** (no admissible kStar with a surviving HasConstantTerm=True sector).
Explicit-fan lifts (k=2 z0=10^4, k=4 z0=10^2) were run to expose the structure:

| Scheme | I1 sum | true sigma | HasConstantTerm/sector |
|--------|--------|------------|------------------------|
| Unlifted | 0.001684683411198275 | 0.039640142688450335 | (n/a) |
| Lifted k=2 | 0.0016841402792674741 | INFINITE | {False, False, False} |
| Lifted k=4 | 0.0016841402792674741 | INFINITE | {False, False, False} |

ALL lifted sectors (both k) are HasConstantTerm=False -> **INFINITE true variance** ->
**NOT BENEFICIAL** (structural, documented).  Correctly NOT auto-lifted.

---

## Case D (FINITE-VARIANCE demonstration): `P = 1 + 10^6 x1`, B={-2}, 1D  (lift {1}, k=2)

Hand-checkable Toy-0 with B=-2 (exact integral = 1e-6); f^2 IS integrable.  After
lifting (z0=1000), ALL surviving sectors have HasConstantTerm=True AND convergent I2.

| Scheme | I1 sum | exact | true sigma | EXACT VarRed |
|--------|--------|-------|------------|--------------|
| Unlifted | 1.0000000002120523*^-6 | 1.0e-6 | 0.0005773494031899867 | - |
| Lifted k=2 | 1.0000000000000076*^-6 | 1.0e-6 | 0.00010266413978186456 | **31.625681822076615** |

HasConstantTerm per sector: {True, True, True}.  This is the class-(i) branch:
a GENUINE finite variance reduction (VarRed = 31.625681822076615 >= 1), corroborated by EXACT
trueSigma -- lifting genuinely helps when the lifted sectors keep an O(1) constant term
and f^2 stays integrable.

---

## Summary (HONEST)

| Case | true sigma_unlifted | true sigma_lifted | EXACT VarRed | verdict |
|------|---------------------|-------------------|--------------|---------|
| A | 0.01658965697013982 | INFINITE | INFINITE | NOT-beneficial (class ii; divergent I2 in HasConstantTerm=False sector; warning fires) |
| B | INFINITE | INFINITE | INFINITE | report-only (INFINITE: divergent I2 even with HasConstantTerm=True) |
| C | 0.039640142688450335 | INFINITE | INFINITE | NOT-beneficial (class ii; all HasConstantTerm=False; not auto-lifted) |
| D | 0.0005773494031899867 | 0.00010266413978186456 | **31.625681822076615** | class (i): GENUINE finite reduction (all HasConstantTerm=True, all I2 converge) |

**Key correction vs the prior RED run:** Case A is reported as INFINITE true
variance (HasConstantTerm=False sector with divergent I2), NOT an 18x sampled-sigma
pass.  The ONLY genuine variance reduction (Case D, 31.625681822076615x) is corroborated by
EXACT trueSigma with all sectors finite-I2.  No finite VarRed is claimed from sampled sigma.