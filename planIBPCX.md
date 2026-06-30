# planIBPCX — IBP with complex divergent directions, batched over kinematics, cross-validated against Subtraction

> Spec for the next feature increment. **Nothing here is implemented yet.** Written
> 2026-06-29 after reading `why_not_IBP.md`, the M2/M2b engine code, and the
> Case-A/B taxonomy in `planAXpDIV.md §3` / `planCXLIFTDIV.md §9` / `planAXpDIVv2.md §4.5`.
>
> Engine references are line numbers in `tropical_eval.wl` (v3-build @ 9593338).
> Read this together with `why_not_IBP.md` (the honest assessment it supersedes on
> two points) and `SUMMARY.txt` L2 "Known limitations" (lines 640–652, the claim §1 revises).

---

## 0. TL;DR

Two defects, one validation gap.

1. **IBP × complex exponents.** When the divergent direction carries a nonzero
   imaginary exponent `θ` (from `Im(B)·D` or `Im(A)·M`), `IBPProcessSector` refuses
   (`TropicalEval::splitdivmono`, `tropical_eval.wl:5032–5042`). The stated reason
   (`SUMMARY.txt:644–652`) is that "the real 1/ε pole moves off ε=0" and that "at
   pinned ε the divergent-slot phase coefficient ~1/ε makes the integrand oscillate
   too fast for **either route** to resolve." **The "either route" is wrong for IBP.**
   The cone is *not divergent* — `1/(c_k ε + iθ)` is regular at ε=0 (finite `1/(iθ)`,
   no pole) — and IBP computes that finite value through a **numerically tame**
   integrand, precisely because IBP *raises* the divergent exponent before flattening
   (bulk phase coefficient `θ/α₀ᵏ(raised)` = O(θ), not `θ/(c_kε)`). Verified
   symbolically and numerically (§1.3). The fix carries the full complex `s_k = c_k ε + iθ`
   through the prefactor and the Laurent assembly; `splitdivmono` is downgraded from a
   refusal to a branch.

2. **IBP × batching.** `emitMain` ignores `Batch` on the IBP path and falls back to
   per-kinematic-point Vegas (`tropical_eval.wl:5635–5640`, `GenerateCppMonteCarloIBP::usage`
   `:123`). The batched chunked-ncomp Vegas (`:3271–3359`) already exists for the
   convergent path; it **generalizes directly** to the IBP integrand table because the
   IBP per-kp output contract (one `conv-sum` row + `n_ibp` function rows per kp,
   `:5966–5999`) maps cleanly onto "one Vegas call per (function, kp-chunk) sharing
   samples across the chunk." No new combine; the driver's per-kp assembly (`:6012–6072`)
   is unchanged.

3. **Cross-validation (the confidence gap).** IBP and Subtraction must both run on
   integrals of the problematic type and agree, so that a single self-consistent number
   is never trusted (invariant #4). §5 defines the comparison matrix, the fixture
   geometry, and the new/updated cross-checks. Crucially it is **honest about which
   oracle is valid where**: the genuine-complex-pole family (θ=0) is a 4-way check
   (IBP / inline Subtraction / `LaurentFromSubtraction` / NIntegrate); the off-axis
   family (θ≠0) is IBP vs **NIntegrate** (+ reseed + batched-vs-per-kp), because the
   Subtraction and pinned-ε routes genuinely cannot resolve it (§1.2).

Scope is the engine (`tropical_eval.wl`) plus its `TEST/cc_*` cross-checks. Wiring the
IBP Laurent into a bespoke Subtraction-combine production pipeline (the bubble's
`TOTAL`/`G0` channel, `why_not_IBP.md §1`) is a separate downstream task (§7.5).

---

## 1. The misdiagnosis, precisely

### 1.1 Two different "Case B"s — do not conflate them

The `splitdivmono` text and `SUMMARY.txt:644` call this "the complex analogue of
Case B." That phrasing has caused the over-refusal. There are **two unrelated
couplings**, and only one is genuinely hard:

- **Geometric Case B** (`planAXpDIV.md §3`, `TropicalEval::liftdivdomain`, `:213`,
  `:1284`): the divergent variable enters the *lifted domain indicator* (`ic_k ≠ 0`).
  The `y_k→0` endpoint and the half-space face couple; IBP in `y_k` picks up an extra
  face boundary term. This is real and stays refused. **Not this plan.**

- **Complex pole location** (`splitdivmono`, `:214`, `:5032–5042`): the divergent
  *exponent* acquires a constant imaginary part `iθ`. Nothing geometric couples; the
  flattening, the boundary, and the IBP terms are structurally identical to the real
  case. The *only* thing that changes is the scalar prefactor `1/s_k`. **This is the
  defect.** Calling it "Case B" borrowed the difficulty of geometric Case B and applied
  it where it does not exist.

### 1.2 Why the divergent cone with θ≠0 is *not* divergent, and why only IBP can integrate it

The divergence criterion is geometric: `Re(α₀ᵏ) ≤ 0` (`:4832`, `:4988`). A cone with
imaginary exponent has `α₀ᵏ = iθ` → `Re = 0 ≤ 0` → flagged divergent. But the actual
ε-pole structure is governed by the *full* complex exponent `s_k = c_k ε + iθ`:

```
∫₀¹ y^{s_k − 1} g(y) dy  =  g(1)/s_k − (1/s_k) ∫₀¹ y^{s_k} g'(y) dy        (one IBP step)
```

- **θ = 0:** `1/s_k = 1/(c_k ε)` → genuine `1/ε` pole, residue `1/c_k`. (Current path.)
- **θ ≠ 0:** `1/s_k = 1/(c_k ε + iθ) = (1/iθ)(1 − (c_k/iθ)ε + …)` → **regular at ε=0**,
  no pole, finite value `g(1)/(iθ) − (1/iθ)∫₀¹ y^{iθ} g'(y) dy`.

So the user is right: such cones produce `1/(ε+iθ)` and are **finite** (the "pole" sits
at `ε = −iθ/c_k`, off the real axis). But finite ≠ trivial to integrate numerically:
`y^{iθ−1}g(y)` has `|·| = y^{−1}|g|`, **not absolutely integrable** — a plain
MC/Vegas of it has unbounded variance. Three routes, three fates:

| route | what it integrates near `y_k→0` | numerically? |
|---|---|---|
| convergent emitter (no subtraction) | `y_k^{iθ−1} g` directly | **diverges** (∥·∥ ~ 1/y) |
| inline Subtraction / pinned-ε `LaurentFromSubtraction` | flattens `y_k` by **vanishing** `α₀ᵏ = c_kε`: phase `y'^{iθ/(c_kε)}`, coefficient `θ/(c_kε)` | **fast oscillation, fails as ε→0** ✗ (this is what `SUMMARY.txt:650` correctly describes) |
| **IBP** | **raises** `α₀ᵏ` by `e_{m,k} ≥ 1` *before* flattening: bulk `y'^{iθ/(α₀ᵏ+e)}`, coefficient `θ/(α₀ᵏ+e)` = O(θ); boundary at `y_k=1` has no phase | **bounded, resolves cleanly** ✓ |

IBP is *the* route built for this: it trades the non-integrable `1/y` endpoint for a
finite boundary term plus a bounded bulk integrand. The refusal threw away the one
method that works.

### 1.3 Verified premise (sanity check already run, `scratchpad/verify_premise.wl`)

With `c_k=2`, `θ=7/10`, `g(y)=1+3y+2y²`:

- `∫₀¹ y^{c_kε+iθ−1} g dy = 10[(7i+20ε)⁻¹ + 3((10+7i)+20ε)⁻¹ + 2((20+7i)+20ε)⁻¹]`;
  `Series` at ε=0 has **pole coefficient 0**; value at ε=0 = `194300/66901 − 1475060i/468307` (finite).
- IBP form `g(1)/(iθ) − (1/iθ)∫₀¹ y^{iθ} g' dy` equals that value **exactly** (`MATCH → True`).
- `|y^{iθ−1}g|` at `y=10⁻³,10⁻⁶,10⁻⁹` = `{10³, 10⁶, 10⁹}` (blows up); `|y^{iθ}g'|` = `{3.0, 3.0, 3.0}` (bounded).
- Vanishing-exponent phase coefficient `θ/(c_kε)` at `ε=10⁻²,10⁻³` = `{35, 350}`; IBP raised-exponent coefficient `θ/e` for `e=1,2,3` = `{0.70, 0.35, 0.23}`.

---

## 2. What is already correct in the IBP path (so the change stays small)

The imaginary part already rides through the IBP integrands; only the scalar pole
prefactor is wrong. Concretely:

1. **Bulk coefficients are already complex.** `IBPExpandOneVariable` brings down the
   *full* complex `Bj = termPolyExps[[j]] + I·imB[[j]]` (`:4767–4768`), so the
   brought-down `S0/S1` carry `Im(B)`.
2. **Bulk & boundary phases already use the raised exponent.** `termMonoFactorLog` /
   `termMonoPhaseLog` divide the fixed `Im` numerators by **this term's** `α₀`
   (`:5182–5192`), which post-IBP has `Re ≥ 1` (the "verify all α₀ > 0" guard at
   `:5148–5157`). The boundary's `bndMonoFactorLog`/`bndMonoPhaseLog` drop slot `k`
   and divide by `bndA0 = a0[[ndVars]] > 0` (`:5108–5120`). → the C++ integrands
   (`GenerateCppMonteCarloIBP`) are **already emitted with bounded `θ`-phases**. No
   integrand/codegen change is needed for the integrands themselves.
3. **The `1/a_k` prefactor is applied in the driver, not baked into terms.**
   `ibpPrefactors = −1/ak` (`:4871`) is informational; the actual pole division is
   `poleCont = (bndBase − S0)/ckS`, `finiteCont = ((bndLog − S1) − rkS·(bndBase − S0))/ckS`
   in `evaluateTropicalIBPDriver` (`:6044–6048`), using the **real** `ck`. **This is the
   single locus of the bug.**

So the imaginary exponent reaches the integrands correctly; it just never reaches the
*prefactor*. The fix injects `iθ` into the prefactor and the Laurent assembly.

---

## 3. Design — IBP × complex (Deliverable 1)

### 3.1 Carry `θ_k` out of `IBPProcessSector`

Replace the refusal block (`:5032–5042`) with a **computation**: keep the existing
`ckImag` accumulation (it is already exactly the `θ_k` we need — `Σ_j imB_j·D_{j,k}` +
`lmp["Num"][[k]]`), and instead of `Message[splitdivmono]; Return[$Failed]`, store it:

```
"ImagPole" -> ckImag        (* θ_k; 0 for real/unlifted sectors -> byte-identical *)
```

added to the returned `IBPSectorData` (`:5223–5265`), alongside the existing
`"ck" -> ck`, `"rk" -> rk`. `ckImag` is built from exact `Im[…]` data (exactness
preserved, invariant #1; no `N[…]`). Keep `"AnalyticPole" -> 1/ck` as-is for the real
case but stop treating it as authoritative (see §3.3).

`ckImag` may be **kinematic-dependent and symbolic** (e.g. `Im(B)=μ` for imaginary μ).
Keep it symbolic here; the branch decision happens at assembly time (§3.4).

### 3.2 Relax the `badck` guard (`:4866–4869`)

`IBPReduceSector` aborts (`TropicalEval::badck`) when `c_k = d(a_k)/dε|₀ = 0`. That is
correct only when `θ_k = 0` (an unregulated genuine log divergence). When `θ_k ≠ 0` the
endpoint is regulated by the oscillation (prefactor `1/(iθ_k)`, no ε needed). Change the
guard to fire only when **both** `c_k = 0` **and** `θ_k = 0`. Pass `imB` (already in
scope, `:4825`) into the same `Σ imB_j·D_{j,k}` + `Im(A)` test used in §3.1, or defer the
combined guard to `IBPProcessSector` where `lmf`/`lmp` are available. (Edge case; the
bubble's regulated cones have `c_k ≠ 0`, so this is a completeness fix, not a blocker.)

### 3.3 The complex-pole Laurent in the driver (`:6021–6066`)

This is the heart of the change. Today:

```
ckS = ...["ck"] /. kinRules;   rkS = ...["rk"] /. kinRules;
poleCont   = (bndBase − S0)/ckS;
finiteCont = ((bndLog − S1) − rkS·(bndBase − S0))/ckS;
... finiteTotal += DLogPrefactor·poleCont ...
```

Generalize the prefactor from `1/(c_kε)` to `1/(c_kε + iθ_k)`. Let
`θ0 = (...["ImagPole"] /. kinRules)`:

- **θ0 == 0 (real pole):** unchanged. `poleCont`, `finiteCont`, and the
  `DLogPrefactor·poleCont` finite correction exactly as now → **byte-for-byte identical
  numbers** on every currently-supported sector (this is the #25 / numeric-regression
  guarantee).
- **θ0 ≠ 0 (off-axis, finite):**
  ```
  poleCont   = 0;
  finiteCont = (bndBase − S0) / (I·θ0);
  ```
  No `bndLog/S1` term (it is O(ε), dropped); no `DLogPrefactor` term (also O(ε)). This is
  the eps⁰ coefficient of `(pfBase(0)/(iθ0))·[(bndBase−S0) + ε(…)]`, derived in §1.2 and
  verified in §1.3.

`bndBase`, `S0` are exactly the MC outputs the driver already reads (`:6025–6041`); they
already include the bounded `θ`-phases (§2.2) and the complex brought-down `S0` (§2.1).
So **nothing upstream of this branch changes** for θ≠0 — the integrands are already right.

### 3.4 Where the pole/no-pole decision lives (exactness boundary)

The presence of the `1/ε` pole depends on `θ0`, which can vary across the kinematic grid
(real-μ slice: `θ0=0`, pole; imaginary-μ slice: `θ0≠0`, no pole). This is physics, not a
bug, and the branch is genuinely discontinuous. Resolve it without violating invariant #1
("no float decides the *decomposition*"):

- The **decomposition** (sectors, flattening, IBP terms, boundary, phases) is built once,
  symbolically, with `θ_k` carried as a symbol — *identical regardless of θ_k*. Invariant
  #1 is about this stage; it is untouched.
- The **branch** is in `evaluateTropicalIBPDriver`'s post-MC assembly (`:6012–6072`),
  which already substitutes numeric `kinRules` and multiplies MC floats. Deciding
  `θ0 == 0` there is a *numeric assembly* decision, not a decomposition decision.
- **Prefer the exact test:** decide per sector from the **symbolic** `ImagPole` via
  `PossibleZeroQ` *before* `kinRules`. If `PossibleZeroQ[ImagPole]` → real-pole branch for
  the whole call (provably real; preserves byte-identity). Else → complex branch.
- **Mixed grid guard:** if a single call's `ImagPole` is symbolically nonzero but a
  *specific* kp makes it vanish (e.g. μ ranges across 0), `(bndBase−S0)/(I·θ0)` would
  divide by zero. Detect `θ0==0` per kp inside the complex branch and fall back to the
  real-pole formula for that kp (removable-singularity-safe). Document that splitting the
  grid at the real slice is the clean alternative.

### 3.5 What stays refused (unchanged $Failed)

- **Nested / higher-order poles** (`>1` divergent variable per sector): still
  `TropicalEval::nestedIBP` (`:201`, `:4986–4992`, `:5010–5013`). A direction with
  `θ≠0` still counts as geometrically divergent (`Re(α₀)=0`), so a sector with one real
  pole *and* one off-axis direction is `nDiv=2` → refused today even though it is only one
  true pole. Note as future work (§8); out of scope here.
- **Geometric Case B** (`liftdivdomain`, `:213`): unchanged.
- **Inline symbolic-ε Subtraction with lift+complex** (`splitliftdiv`, `:215`):
  unchanged — that path's `G0/G1` emitters are not phase-aware and (per §1.2) flatten by
  the vanishing exponent, so it genuinely cannot resolve θ≠0. Leave the clean refusal.

### 3.6 Exactness / byte-identity checklist (Deliverable 1)

- `ImagPole` built from exact `Im[…]`; `FreeQ[_Real]` at the `MmaToC` boundary preserved.
- Complex branch only activates when `ImagPole` is not provably zero → real/Direct path
  unchanged → IBP goldens (`TEST/baselines/codegen_goldens/ibpdiv_*.cpp`) byte-identical.
- The driver assembly is Mathematica post-processing (not emitted C++); the θ0==0 branch
  reproduces current numbers bit-for-bit.

---

## 4. Design — IBP × batching (Deliverable 2)

### 4.1 What "batched" means and what already works

The FINAL4 architecture (`why_not_IBP.md §2`) does the tropical decomposition **once per
sector** and serves all kinematic points by making kinematics runtime `params[]` to one
batched-Vegas binary that shares one adaptive grid / sample set across points. The
convergent path implements this as chunked-ncomp Vegas (`:3271–3359`): outer loop over
sectors `for s in N_INTEGRANDS`, inner chunk over kp `for k0 in 0..n_kp step kpPerChunk`,
one `Vegas(dim_s, ncomp=2·cs, cubaBatch, …)` per `(sector, chunk)`; `cubaBatch`
(`:3272–3283`) evaluates `integrand_table[s]` at each kp in the chunk via global
`gb_kp0`/`gb_chunk`, packing `ff[2c]=re, ff[2c+1]=im`; `kpPerChunk = CubaMaxComp/2`
(`:3316`, default 256).

Note: the **plain-MC IBP** main already loops all kp in one binary
(`#pragma omp parallel for (kp)`, ~`:3155–3261`), each kp an independent RNG stream — so
MC is already "one binary, all points." The missing piece is **batched Vegas** (shared
grid), which is what gives the variance win at scale.

### 4.2 Why it falls back today, and why the fix is plumbing not a refactor

`emitMain` routes off `info["IsIBP"]` and **ignores `batch`** on that path (`:5639–5640`);
`GenerateCppMonteCarloIBP` prints "using per-kp VEGAS" (`:5635–5637`). The asserted
blocker (`why_not_IBP.md §1`, Explore survey) is that the IBP path emits **multiple
functions per sector** (boundary-base, boundary-log, per-term base/log) plus the
convergent sum, whereas the convergent batch sums to one integrand per kp.

But the convergent batch **already loops over a function index** (`s in N_INTEGRANDS`).
The IBP function table (`IBPFuncMap`, `NIBPFuncs`, built in `emitIntegrandDefinitions`
~`:2561`, `:2954`) is just a *longer* list of functions, each with its own integration
dim. So the generalization is: iterate the **IBP function table** = `[conv-sum]` ++
`[boundary-base, boundary-log, term-base, term-log …]` instead of `[sectors]`, one
`Vegas(dim_f, ncomp=2·cs, …)` per `(function, kp-chunk)`. Samples are shared across the
chunk's kp (the variance win); each function keeps its own dim (just like sectors do now).

### 4.3 Output contract — preserved exactly

The per-kp result block is `[conv-sum, IBP-func-0, IBP-func-1, …]`, `nLinesPerKP =
1 + nIBPFuncs` (`:5967`), read back by `mcCx[fid] := mcRawResults[[kpOffset + fid −
NConvergent + 2]]` (`:5997`). The batched emitter must write the **same** layout (row 0 =
convergent sum for that kp, rows 1..n_ibp = the IBP functions). Then the entire driver
combine (`:6012–6072`, including the §3.3 complex branch) is **unchanged** — it does not
care whether the rows came from per-kp or batched Vegas. The symbolic kinematic-dependent
`Coeff0/Coeff1` are applied in the driver after integration (`:6033–6038`), so they never
touch the batched integration.

### 4.4 Concrete changes

1. `emitMain` (`:2978`): add an `integrator === "VEGAS" && batch && isIBP` branch that
   reuses the `cubaBatch`/`gb_*` machinery (`:3271–3359`) but iterates the IBP function
   table and writes the `(1 + nIBPFuncs)`-row-per-kp layout. The convergent sectors are
   summed per kp into row 0 (mirror the per-kp IBP main's row-0 convention).
2. `GenerateCppMonteCarloIBP` (`:5621`): drop the "batch ignored" note (`:5635–5640`);
   pass `batch` through to the new branch.
3. `evaluateTropicalIBPDriver` (`:5783–6102`): allow `"Batch" -> True` to reach codegen
   (today `:5705` forces per-kp); read back the batched output (same parser, the layout is
   identical).

### 4.5 Exactness / regression checklist (Deliverable 2)

- Non-batched IBP (MC + per-kp Vegas) emission **byte-identical** — the new code path only
  activates under `Batch && VEGAS && IsIBP`.
- New gate "**batched IBP == per-kp IBP**" within MC tolerance + reseed (the §6.3 #23
  analog: `SUMMARY.txt:795` gates batched VEGAS == per-kp VEGAS for the convergent path —
  mirror it for IBP).
- Add the missing IBP codegen golden gate: `phase2_selfgate.wl` currently emits
  `ibpdiv_*.cpp` goldens but does **not** verify them (Explore survey C). Add an
  `ibpdiv_vegas_batch_*.cpp` golden and assert byte-identity, closing the gap.

---

## 5. Cross-validation — Subtraction vs IBP on the problematic type (Deliverable 3)

The user's requirement: never trust one engine. The honest design splits the problematic
type into two families because **the valid oracles differ**.

### 5.1 The comparison matrix

| family | divergent direction | true structure | IBP gives | independent oracle(s) | "Subtraction vs IBP"? |
|---|---|---|---|---|---|
| **F1: genuine complex pole** (θ=0, but `Im(B)`/`Im(A)` on *other* slots, complex finite part) | `θ_k = 0` | real `1/ε` pole, complex residue/finite | (pole, finite), both complex | inline Subtraction **and** `LaurentFromSubtraction` **and** NIntegrate(complex original) | **Yes — 4-way agreement.** This is the strong "subtraction vs IBP" check. |
| **F2: off-axis** (θ≠0, the imaginary-μ bubble cone) | `θ_k ≠ 0` | **no pole**, finite `~1/(iθ)` | pole≈0, finite | **NIntegrate(complex original)** + reseed + batched-vs-per-kp IBP | **No** — inline Subtraction (`splitliftdiv`/vanishing-exponent) and pinned-ε both fail (§1.2). IBP is the *unique* engine route; NIntegrate is the truth. State this plainly. |

F1 is what the user means by "compare the two methods" and it is fully supportable today
(the machinery exists — cc_46/cc_48 already do a 3-way version for θ=0; §5.4 strengthens
it). F2 is the genuinely-new capability; being explicit that Subtraction is *not* a valid
F2 oracle is itself a correctness statement — it stops anyone from "cross-checking" two
routes that share the same blind spot and getting false reassurance.

### 5.2 Fixture geometry (how to actually build θ≠0 cones)

From the `const-a-nonzero-fixture-geometry` memory and `planAXpDIVv2 §4.5`: the divergent
direction picks up `θ_k = Σ_j Im(B_j)·D_{j,k} + Im(A)·M_k`. To force **θ_k ≠ 0** (F2), put
the **divergence on the complex-weight variable** — i.e. the variable whose monomial/poly
exponent carries the imaginary part is the *same* variable that goes soft (`Re α₀=0`). To
force **θ_k = 0 with complex elsewhere** (F1), put the imaginary weight on a *non*-divergent
slot (the cc_48 pattern). Require `HasConstantTerm = True` so the pole is clean
(`why_not_IBP.md §4`, avoids the L3 undercount). Single divergent variable (avoid
`nestedIBP`).

### 5.3 Concrete fixtures (templates; tune coefficients so NIntegrate converges)

**F1 (θ=0, genuine complex pole)** — extend the cc_48 family
(`{1+x₁+x₂+x₁x₂}`, `A={2ε−1,0}`, `B={−(2+i/2)}`): divergence on `x₁`, imaginary `B` does
not net onto the `x₁` soft direction → `θ=0`. Run `Method->"Subtraction"`,
`Method->"IBP"`, `LaurentFromSubtraction`, and NIntegrate of the complex original;
assert all four (pole, finite) agree.

**F2 (θ≠0, off-axis)** — divergence carried on the imaginary-weight variable, e.g.
```
spec = <| "Polynomials"        -> { 1 + x[1] + x[2] + x[1] x[2] },
          "MonomialExponents"  -> { 2 eps - 1, 0 },         (* x1 soft -> divergent *)
          "PolynomialExponents"-> { -2 + (1/3) I },          (* Im(B) routes onto x1 *)
          "Variables"          -> { x[1], x[2] },
          "KinematicSymbols"   -> {}, "RegulatorSymbol" -> eps |>;
```
chosen so the fan sends the `Im(B)·D` weight onto the `x₁` soft slot → `θ_{x₁} ≠ 0`. Run
`Method->"IBP"` (expect pole ≈ 0, finite complex); oracle = NIntegrate of
`∫₀^∞∫₀^∞ x₁^{2ε−1} (1+x₁+x₂+x₁x₂)^{−2+i/3} dx` analytically continued, at small ε, plus a
reseed and a batched-vs-per-kp IBP comparison. **Confirm the cone actually realizes θ≠0**
by inspecting `ImagPole` (a fixture that silently lands θ=0 proves nothing). Build a
2–3 point fixture sweep (and at least one kinematic-`μ` variant so `ImagPole` is symbolic,
exercising §3.4).

### 5.4 Cross-checks to add / update

- **`cc_47` case (C)** (`TEST/cc_47.wl:24–29`): today asserts the `splitdivmono` **refusal**
  is the deliverable. After the fix, flip its expectation from "$Failed" to "IBP (pole≈0,
  finite) agrees with NIntegrate of the complex original" (F2). Keep a separate
  still-refused case for genuine geometric Case B / nested, so the refusal path stays tested.
- **New `cc_49`** (or extend `cc_48`): the F1 4-way agreement (IBP / inline Subtraction /
  `LaurentFromSubtraction` / NIntegrate), pole and finite, tol matching cc_46/cc_48
  (`|Δpole| < 5×10⁻³…0.05`, recon `< 2%`).
- **New `cc_50`**: the F2 off-axis fixture(s) — IBP vs NIntegrate (+ reseed + batched-vs-
  per-kp), explicitly documenting "Subtraction not applicable here (§1.2)."
- **Batching gate**: "batched IBP == per-kp IBP" (§4.5), and the `ibpdiv_vegas_batch`
  golden, added to `phase2_selfgate.wl` so it is a hard CI gate (`Exit[1]` on FAIL).
- Keep every numeric PASS on **≥2 independent oracles** (invariant #4): F1 has 4, F2 has
  NIntegrate + reseed + batched-vs-per-kp (≥2 truly independent of the IBP assembly).

---

## 6. Invariants this work must preserve

1. **Exactness (#1).** `ImagPole`/`θ_k` from exact `Im[…]`; no `N[…]` upstream of `MmaToC`;
   the decomposition makes no float decision (§3.4 keeps the branch in numeric assembly).
2. **Byte-identity (#25).** Real/Direct/θ=0 path emits byte-identical C++ and reproduces
   current numbers exactly; complex branch and batched-IBP branch activate only when their
   trigger (θ≠0 / `Batch&&VEGAS&&IsIBP`) is present.
3. **Known-limitation ⇒ clean $Failed (#3).** Nested poles, geometric Case B, inline-Sub
   lift+complex stay refused with their existing messages. Add no silently-wrong numbers.
4. **≥2 independent oracles per PASS (#4).** Enforced per §5.4.
5. **Never edit `OLD_CODE/` (#5).**

---

## 7. Phases & gates

- **P1 — IBP × complex (§3).** `ImagPole` plumbed; `splitdivmono` → branch; `badck`
  relaxed; driver complex Laurent. **Gate:** cc_47(C) flips to numeric PASS vs NIntegrate;
  cc_45/cc_46/cc_48 unchanged (regression); IBP goldens byte-identical; θ=0 numbers
  bit-identical on a spot sector.
- **P2 — IBP × batching (§4).** Batched-ncomp Vegas for the IBP function table; output
  contract preserved. **Gate:** "batched IBP == per-kp IBP" within MC tol + reseed;
  `ibpdiv_vegas_batch` golden added & byte-stable; non-batched emission byte-identical.
- **P3 — Cross-validation suite (§5).** F1 cc_49 (4-way), F2 cc_50 (IBP vs NIntegrate),
  batching gate into `phase2_selfgate.wl`. **Gate:** all PASS with ≥2 oracles each.
- **P4 — Docs.** Rewrite `SUMMARY.txt:640–652` from "out of scope — clean $Failed is the
  deliverable" to "off-axis (θ≠0) divergent directions are **supported on the IBP route**
  (finite, pole≈0); the Subtraction and pinned-ε routes remain inapplicable (fast
  oscillation); IBP batches over kinematics." Update the `splitdivmono` message to fire
  only for the still-refused residue (nested/geometric Case B), and the MANUAL §"Error
  Messages" + the IBP section. Update `why_not_IBP.md` (its §2 "cannot batch" and §3
  "still $Failed" claims are both resolved by this plan).
- **P5 — (optional, downstream) production-pipeline wiring (§7.5).**

### 7.5 Note on Subtraction-combine production pipelines

`why_not_IBP.md §1` observes that a bespoke pipeline built on the inline-Subtraction
`TOTAL`/`(γ−1)·G0` channel cannot directly consume IBP's `(pole, finite)` Laurent. That
impedance mismatch is real but **orthogonal** to this plan: the engine-level
cross-validation in §5 runs through `EvaluateTropicalMC`/`EvaluateTropicalMCIBP` and does
not require the production combine. Feeding batched-IBP Laurents into such a combine (or,
equivalently, adopting batched `LaurentFromSubtraction` for the divergent cones) is the
production-integration step and is scoped as P5, after the engine capability and its
cross-checks land.

---

## 8. Scope boundaries

**Delivers:** IBP on a **single** divergent direction carrying an imaginary exponent
(`θ≠0`, off-axis, finite) — no longer refused, computed via the bounded IBP bulk/boundary;
the genuine-complex-pole (`θ=0`) case strengthened to a 4-way Subtraction↔IBP↔
`LaurentFromSubtraction`↔NIntegrate check; IBP integration **batched over kinematics**
(shared-grid Vegas), preserving the per-kp output contract and combine.

**Now also implemented (planIBPMULTIDIV.md):**
- The mixed case "one real pole + N off-axis directions" in the same cone is **no longer
  refused** (UNLIFTED): `ibpDivClass` counts only genuine poles toward `nestedIBP`, and the
  full `2^(N+1)`-corner iterated IBP (`IBPBuildCorners`) assembles it. Verified vs the
  closed-form Dirichlet oracle for `N=1` (cc_51) and `N=2`. The driver multiplies the
  per-corner sum by `[1/(c_k ε)]·∏_j[1/(iθ_j)]` (still orders −1, 0 only — no `1/ε²`).

**Still refused / out of scope (clean `$Failed`, future work):**
- `>1` GENUINE pole (each `Re=0 ∧ θ=0`): a real `1/ε^{d≥2}` (`nestedIBP`).
- The LIFTED co-located "one pole + off-axis" case: `ProcessSectorLifted`'s pivot admittance
  counts realified `Re≤0` directions (not poles), so it refuses such a pivot
  (`liftnopivot`); the off-axis-aware lifted pivot admittance is future work (cc_52 (B)).
- Geometric Case B — divergent variable couples to the lifted domain face (`liftdivdomain`).
- Inline symbolic-ε Subtraction with lift+complex (`splitliftdiv`) and the pinned-ε route
  for θ≠0 — genuinely fail (fast oscillation); IBP is the supported route.
- An imaginary exponent that *itself* scales like `1/ε` (eps-dependent `Im`) — would
  reintroduce fast oscillation; guard and refuse if it ever arises (standard real-ε
  regulator with eps-free `Im(B)/Im(A)` is fine and is the case here).

---

## 9. Open questions / risks

1. **Mixed grid (§3.4).** A single call whose `μ` crosses the real axis mixes θ=0 (pole)
   and θ≠0 (no pole) kp. The per-kp fallback handles it, but the cleaner story may be to
   require the caller to split at the real slice. Decide during P1 from a μ-sweep fixture.
2. **F2 oracle quality.** NIntegrate of the analytically-continued oscillatory original may
   need method hints (`"DoubleExponentialOscillatory"` / large `MaxRecursion`) and a small
   pinned ε. Validate the oracle itself converges before trusting agreement; if NIntegrate
   is shaky, add a high-statistics IBP reseed envelope as the second oracle.
3. **Fixture realizes θ≠0.** The fan geometry decides whether `Im(B)·D` lands on the soft
   slot. The fixtures in §5.3 are templates — confirm `ImagPole ≠ 0` by inspection before
   asserting; otherwise adjust the polynomial so the soft ray couples to the imaginary
   weight.
4. **Batched-IBP CUBA limits.** `(1+n_ibp)` functions × kp can be many components; the
   existing `CubaMaxComp=512` chunking (`:3316`) already bounds ncomp per call, but verify
   per-function dim handling and that `n_ibp` (boundary+terms) does not explode for
   multi-term sectors. Log any silent truncation.
