# Why the FINAL4 batched lift driver does not use the IBP route for divergent cones

> Written 2026-06-29, after the planRERUN.md aux-variable-lift implementation.
> Honest assessment: **IBP is not "impossible"** — for the bubble's *real-μ*
> divergent lifted cones it would even be *correct*. The reason the driver falls
> back to the standard path instead of using IBP is a combination of one decisive
> architectural incompatibility, one genuine physics refusal, and a couple of
> minor constraints. This file lays all of them out precisely, with engine
> references, so the decision is auditable and the path to enabling IBP later is
> clear.

Engine references below are line numbers in `Bubblepp/tropical_eval.wl`
(v3-build@9593338, md5 67ca04ce).

---

> **UPDATE (planIBPCX.md, implemented).** Two of the three obstacles below are now
> resolved at the engine level:
> - **§2 "IBP cannot batch over the grid" — RESOLVED.** `Integrator->"VEGAS"`,
>   `"Batch"->True` now runs chunked-ncomp Vegas over the IBP *function* table,
>   sharing one sample set across each kp-chunk, preserving the `(1+n_ibp)`-row-per-kp
>   output contract (cc_50 (B): batched IBP == per-kp IBP). One binary, all points.
> - **§3 "IBP $Faileds on imaginary-μ divergent cones (`splitdivmono`)" — RESOLVED.**
>   An off-axis divergent direction (`θ_k≠0`) is *not* a 1/ε divergence: the cone is
>   finite (pole at `ε=-iθ_k/c_k`, off the real axis). IBP is the **unique** route
>   that resolves it (it raises the divergent exponent before flattening). The
>   driver carries `θ_k` as `ImagPole` and branches the Laurent (cc_47 (C), cc_50).
>   `splitdivmono` now fires only for the narrow residue of an *ε-dependent* imaginary
>   exponent.
>
> What remains accurate: **§1** (the bespoke FINAL4 *Subtraction-combine* pipeline
> still consumes `TOTAL`/`G0`, not an IBP Laurent — wiring batched-IBP into it is the
> downstream P5 task, planIBPCX.md §7.5) and **§4** (single-pole-per-sector;
> `HasConstantTerm` clean pole). The §2/§3 reasons no longer block IBP as such.

---

## 0. TL;DR

The bubble's UV-divergent cones live only in **bubblepp presectors 1–2**. To lift
them (variance reduction at the squeezed corners, FLAG A) the engine offers three
routes for "lift × divergence":

| route | handles lifted ε-prefactor? | batched over 7381 kinematics? | used by FINAL4 runner? |
|---|---|---|---|
| inline symbolic-ε subtraction (`Method->"Subtraction"`) | **NO** (ε leaks; refused for complex via `splitliftdiv`) | yes | **yes** (the runner's whole design) |
| pinned-ε `LaurentFromSubtraction` | yes (ε pinned numerically) | yes (per pinned-ε) | no |
| **IBP** (`Method->"IBP"`, `IBPProcessSector`) | yes | **NO — falls back to per-kp Vegas** | no |

The runner is built on the **inline symbolic-ε subtraction** path. That path
cannot carry the lifted-divergent prefactor `z0^(a_p(ε)/m_p − 1)` (ε stays live →
the generated C++ references an undeclared `dϵ`). So lifted-divergent cones must
use one of the other two routes. **IBP is rejected primarily because it abandons
the batched-Vegas-over-kinematics architecture that makes the FINAL4 grid scan
feasible** (it reverts to one Vegas run *per kinematic point*). The driver
therefore falls back lifted-divergent presector-batches to the validated standard
(unlifted) path — correct pole and finite, no variance gain on those cones — and
keeps the full lift on every convergent presector (the FLAG-B headline, which is
entirely convergent and unaffected).

---

## 1. What the FINAL4 runner actually computes (the architecture IBP would have to replace)

`run_bispectrum_{pp,pm}.wl` → `bispectrum_run_core.wl` does NOT call the high-level
`EvaluateTropicalMC` driver. It calls the low-level primitives and assembles the
tropical **subtraction** by hand, per presector:

1. `ProcessSector` / `ProcessDivergentSector` per cone.
2. Compile **two** binaries per (presector, batch): a "main" (the TOTAL integrand)
   and a "G0" (the subtracted pole piece).
3. Combine (`CombineBubbleSectorResults`, `bispectrum_utils.wl`):
   - `finite = pref · (TOTAL + (γ−1)·G0)`
   - `pole   = pref · G0`     (the `PoleCoefficient` CSV column)
4. Across diagrams, the physics observable is `2 Re B₊₊ − 2 Re B₊₋`, and the
   pole-column gate asserts bubblepp poles match and bubblepm poles ≡ 0.

**IBP does not produce a `TOTAL` and a `G0`.** `IBPProcessSector`
(`tropical_eval.wl:108`) returns an `IBPSectorData` — integration-by-parts
*boundary* + *bulk* terms, each made convergent at ε→0 — that the IBP codegen
(`GenerateCppMonteCarloIBP`) integrates into a **Laurent pair (pole, finite)
directly**. There is no `(γ−1)·G0` channel to feed into the runner's combine.

Consequently, using IBP for presectors 1–2 means **replacing both the divergent
codegen path and the combine** for those presectors, then reconciling the IBP
Laurent output with (a) the per-presector prefactor handling, (b) the
`PoleCoefficient` column the CSV/`combine_multiruns.wl` expects, and (c) the
cross-diagram `2 Re B₊₊ − 2 Re B₊₋` assembly. This is a real refactor of the
divergent half of the pipeline, not a drop-in swap. It was scoped out so the
validated **convergent** lift (FLAG B) could ship first.

---

## 2. The decisive reason: IBP cannot batch over the 7381-point grid

This is the load-bearing efficiency property of the entire FINAL4 design
(`planAU.md §3.1`): the tropical decomposition is done **once per presector** and
serves **all 7381 kinematic points**, because `KinematicSymbols → params[]` makes
the kinematics *runtime data to one compiled batched-Vegas binary*. One binary,
one shared low-discrepancy sample set, all points.

`GenerateCppMonteCarloIBP::usage` (`tropical_eval.wl:120`) states verbatim:

> "Supports `Integrator->Vegas` (**per-kp CUBA Vegas**; same output layout);
> **`Batch` falls back to per-kp Vegas on the IBP path.**"

So on the IBP path, `"Batch"->True` is silently demoted to **one Vegas run per
kinematic point**. For a divergent presector that means ~7381 separate Vegas
initialisations/integrations (× the IBP boundary+bulk integrands) instead of one
batched pass. At grid scale this is orders of magnitude more expensive — it throws
away the exact property the FINAL4 pipeline was built around. The inline
subtraction path (and `LaurentFromSubtraction`) keep the batched-Vegas binary.

**This alone makes IBP the wrong tool for the full-grid scan**, independent of any
correctness question.

---

## 3. A genuine IBP refusal for imaginary-μ divergent cones (`splitdivmono`)

bubblepp at imaginary μ is **divergent AND complex** (the principal-series pp
diagram). For such a sector the divergent direction can carry a nonzero imaginary
exponent, and the IBP assembly refuses it. `tropical_eval.wl:214`:

> `TropicalEval::splitdivmono` — "IBPProcessSector: cone `1` — the divergent
> direction carries a nonzero imaginary exponent (Σⱼ Im(Bⱼ)·D_{j,k} … PLUS
> Im(A)·M …), so the 1/ε pole acquires an imaginary shift (the pole moves off
> ε=0) that the real-pole IBP assembly does not resolve. … Aborting ($Failed)."

So even if the architecture/batching issues were solved, **IBP would still
$Failed on the imaginary-μ divergent cones whose pole sits off ε=0** (complex
Case-B). Those would need the pinned-ε `LaurentFromSubtraction` route as oracle
anyway. (Real-μ divergent cones — FLAG A, `mu_1.2`, `mu_1.5` — are real-pole and
not subject to `splitdivmono`.)

---

## 4. Minor constraints (satisfied by the bubble, listed for completeness)

- **Single pole per sector.** `IBPProcessSector` supports one divergent variable;
  nested/higher-order poles are refused (`nestedIBP`, `tropical_eval.wl:201`). The
  bubble's UV divergence is a simple 1/ε pole, so this is **satisfied** — not a
  blocker, but must be asserted.

- **The "prefer Subtraction over IBP" guidance does NOT apply here.** Earlier notes
  (planAU §3.2, the v3 package memory) recommend Subtraction because *IBP
  undercounts the pole on `HasConstantTerm=False` sectors* (the L3 limitation).
  But Phase-0 Probe A established the bubble's divergent lifted sectors are
  **`HasConstantTerm=True`** (clean pole). So IBP would NOT undercount on the
  bubble — this particular caveat is **not** a reason against IBP here. (Stated so
  the record is accurate: the case against IBP rests on §2 and §3, not on L3.)

---

## 5. Why the current fallback is the right interim choice

For a lifted presector-batch that yields divergent sectors, `ProcessLiftedConeSet`
(`bispectrum_lift.wl`, `$LiftDivergentMode="fallback"`) returns
`Status->"FallbackStandard"`, and the runner processes that presector-batch via
the **standard unlifted 5-D path**. That path is byte-identical to the validated
base scan, so:

- the **pole and finite parts are correct** (no approximation),
- the **batched-Vegas architecture is preserved**,
- only the **variance gain on those specific cones is forgone** (FLAG A, pp 1–2).

Everything convergent — all of bubblepm (FLAG B, the headline), bubblepp
presectors 3–8, and pm/real-μ — gets the full lift. The squeezed-corner variance
fix that motivated the whole exercise is delivered where it matters most.

---

## 6. What it would take to actually use IBP (future work, P5)

If FLAG-A divergent-cone variance reduction is wanted later, the realistic path is
**not** the per-kp IBP route but the **pinned-ε `LaurentFromSubtraction`** route,
because it *keeps the batched-Vegas binary*:

1. For a lifted-divergent presector-batch, build the lifted sectors as today, but
   instead of the inline ε→0 G0/ck assembly, **pin ε to a few small values**,
   compile the batched-Vegas binary at each pinned ε (the prefactor
   `z0^(a_p(ε)/m_p−1)` is then numeric — no `dϵ` leak), run over the batch's
   points, and **fit `g(ε)=ε·F(ε)` → (pole, finite)** per point
   (`LaurentFromSubtraction`, `tropical_eval.wl:~4650`).
2. Feed the fitted (pole, finite) into the bubble combine in place of the
   `TOTAL`/`G0` pieces for those presectors.
3. For imaginary-μ divergent cones that hit `splitdivmono`, keep the standard
   fallback (or use NIntegrate as oracle).

This is the route `planAU.md §3.2/§11` and `planRERUN.md` P5 already designate as
the production path for lifted divergence. IBP proper (`Method->"IBP"`) is best
reserved as an *exactness oracle at spot points* (where its per-kp cost is
irrelevant), cross-checked against `LaurentFromSubtraction` — which is exactly how
`cc_45`/`cc_46` validated the two routes against each other.

---

## 7. One-line summary

IBP is correct in principle for the bubble's real-μ divergent cones, but it (a)
**defeats the batched-Vegas-over-7381-kinematics architecture** (`GenerateCpp
MonteCarloIBP` reverts to per-kp Vegas), (b) **returns a Laurent rather than the
TOTAL/G0 pieces** the bubble's subtraction combine consumes, and (c) **still
$Faileds on the imaginary-μ divergent cones** (`splitdivmono`). The pinned-ε
`LaurentFromSubtraction` route — not IBP — is the architecture-compatible way to
lift those cones, and is the documented P5 follow-up. Until then the fallback to
the validated standard path keeps them correct.
