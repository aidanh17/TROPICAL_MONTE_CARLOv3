# planAXpDIV.md — Lifting an integral that also has a 1/ε divergence

**Goal.** Remove known-limitation **L2 / plan.md N3** ("lifting and divergence are
mutually exclusive per call") for the case that actually matters, and replace it
with a real capability: a sector may be **lifted *and* carry a simple 1/ε pole**.

**Thesis (the user is right).** There is no fundamental obstruction. The two
features act on different objects — lifting acts on a *coefficient* (variance
reduction, +1 dimension, a δ-constraint), divergence acts on an *exponent*
`a_k(ε)→0` (a 1/ε endpoint pole). After the lifting δ-resolution collapses the
extra dimension, a lifted sector is **structurally an ordinary `SectorData`** of
dimension `n`, and the 1/ε pole — if any — necessarily lives in one of the `n`
post-δ effective exponents `ãtilde`, exactly the object the IBP/subtraction
machinery already consumes. The "incompatibility" is three implementation
choices, not a theorem. One genuine geometric subtlety (the lifted *domain
indicator* coupling to the divergent variable) is real and is scoped honestly
below: we implement the clean case in full and **refuse the hard case with a
clean `$Failed`**, exactly as v3 already refuses nested divergences (L1/L8) and
complex `ãtilde` (`liftcomplex`).

This document is written to slot in beside `plan.md` (same voice, same
exactness discipline §1.4 / §6.7, same "every numeric PASS gated by ≥2
independent oracles" rule from §8.5).

---

## 1. Why they are mutually exclusive *today* (the three real barriers)

All three are in `tropical_eval.wl`.

### Barrier A — `ProcessSectorLifted` refuses any divergent direction
`ProcessSectorLifted` (M1c, lines ~997–1253) does the (n+1)-dim sector +
δ-resolution and produces the post-δ effective exponents `ãtilde` (length `n`).
It then **demands every `ãtilde_j` be real and strictly positive**:

- `isRealPos[av] := PossibleZeroQ[Im[av]] && Re[av] > 0` (line 1080); a pivot is
  a *candidate* only if `And @@ (isRealPos /@ res["atilde"])` (line 1087).
- With a regulator present a divergent direction is `ãtilde_k = c_k·ε + O(ε²)`,
  for which `TrueQ[Re[c_k ε] > 0]` is **False** (undecidable symbolically). So
  the divergent pivot is silently dropped; the sector either picks an unrelated
  pivot or falls through to `liftnopivot` (line 1124).
- Even if a divergent `ãtilde` survives pivot selection, `FlattenSector` flags
  `IsDivergent` and the code aborts with `TropicalEval::liftdivergent`
  (lines 1195–1198): *"the lifted sector is divergent. Lifting supports
  convergent integrals only (plan.md N3)."*
- The emitted `SectorData` hardcodes `"IsDivergent" -> False`,
  `"DivergentVariable" -> 0` (lines 1241–1242).
- The candidate **ranking** `SortBy[..., -Min[Re[N[#["atilde"]]]] &]` (lines
  1132–1136) calls `N[ãtilde]` — with `ε` symbolic this is not a real number and
  the ordering is ill-defined.

### Barrier B — the divergence machinery cannot represent a lifted sector
Even if barrier A handed it a divergent lifted `SectorData`, the divergence
routines silently drop the two pieces of data that make a sector "lifted":

1. **Prefactor base is hardcoded to `Abs[detM]`.** A lifted sector's true
   prefactor base is `(Abs[detM]/Abs[mp]) · z0^(ap/mp − 1)` (line 1192), *not*
   `Abs[detM]`. But every divergence routine rebuilds the prefactor from
   `Abs[detM]` directly:
   - `ProcessDivergentSector`: `g0Prefactor = Abs[detM]/(Times@@g0aVals)` (1510),
     `remPrefactor = Abs[detM]/(Times@@remAvals)` (1571).
   - `IBPProcessSector`: `bndPrefactor = Abs[detM]/(Times@@bndA0)` (4590),
     `flatPrefactor = Abs[detM]/(Times@@alpha0)` (4656).
   - `ValidateSubtraction` / `ValidateIBP`: sector reconstruction uses
     `Abs[detM]` directly (e.g. 1677).
   For a lifted sector these all silently drop the `1/Abs[mp]` and
   `z0^(ap/mp−1)` factors → wrong pole/finite normalization.

2. **The lifted domain indicator is ignored.** A lifted convergent sector
   carries a `"DomainConstraint"` (the `Boole[log_ypstar ≤ 0]` half-space from
   the δ-resolution, lines 1152–1189) that the convergent codegen honors
   (`emitBaseFuncBody`, lines 2122–2143) and `ValidateLiftedDecomposition`
   honors (lines 1328–1334). **None** of the divergence integrands carry it:
   `emitBaseFuncBody` is called for `integrand_g0/g1/ibp_bnd/ibp_term` with the
   `domainConstraint` argument **defaulted to `None`** (lines 2306–2307,
   2325–2326, 2449–2450, 2493–2494, 2506–2507); the hand-written remainder
   integrand (2353–2424) emits no indicator at all.

### Barrier C — the driver forbids the combination
`EvaluateTropicalMC` (M4):
- `routeToIBP` short-circuits `isLifted, False` (lines 3581–3583): *"never route
  a lifted call to the IBP/divergence path."*
- The lifted branch hardcodes `divergentSectors = {}` (line 3686) and treats a
  `liftdivergent`/`liftnopivot` `$Failed` from any sector as a hard abort of the
  whole call (lines 3666–3684).
- `EvaluateTropicalMCLifted` (4118–4202) and `LiftCoefficients` (note line 949)
  already **preserve** `RegulatorSymbol` and ε-carrying exponents into the lifted
  spec — so nothing upstream strips ε. The driver is the only place the combo is
  forbidden.

---

## 2. Why there is no fundamental obstruction (the math)

Original integral, one extreme coefficient `C` in a monomial of some `P_j`, a
regulator `ε` producing a simple pole:

```
I(ε) = ∫_{[0,∞)^n} ∏_i x_i^{A_i(ε)} ∏_j P_j(x)^{B_j(ε)} dx
```

**Lift:** `C·x^α → z^k·x^α`, insert `∫dz δ(z−z0)`, `z0 = |C|^(1/k)`. The aux
variable is given **monomial exponent 0** (`LiftCoefficients`, line 941:
`newMonoExps = Append[monoExps, 0]`). Decompose the (n+1)-dim integral on the
(n+1)-dim fan (the fan is coefficient- and ε-independent — it depends only on the
support). Per sector: change of variables `M`, tropical factor →
augmented effective exponents `a` (length `n+1`); `ε` flows into `a` through
`rawA = (A+1)·M` and `Σ_k B_k(ε)·d_k`.

**δ-resolution** picks a pivot `p`, solves `δ(z−z0)` for `y_p`, leaving `n`
coordinates with effective exponents `ãtilde` (length `n`), prefactor
`(|detM|/|mp|)·z0^(ap/mp−1)`, re-cleared polynomials, and a half-space domain
indicator.

**Two structural facts that make the combination clean:**

1. **The pole cannot live in the pivot/aux direction.** The aux monomial
   exponent is 0 and `z` is fixed by the δ-function — the pivot coordinate is
   *integrated out*, not integrated over `[0,∞)`, so it produces no `∫_0 y^{a−1}`
   endpoint pole. **Every genuine 1/ε pole therefore lives in one of the `n`
   surviving `ãtilde` directions** — precisely the `NewExponents` that
   `IdentifyDivergences` / `IBPReduceSector` already ε-expand.

2. **A post-δ lifted sector is an ordinary `SectorData` plus two scalars/vectors
   of "lifted" metadata.** It has `Dimension = n`, `NewExponents = ãtilde`,
   `ClearedPolys = reclearedPolys`, `PolynomialExponents`, and additionally a
   non-`Abs[detM]` **prefactor base** and a **domain indicator**. The divergence
   routines consume exactly `{Dimension, NewExponents, ClearedPolys,
   PolynomialExponents}` for the *symbolic* reduction; they only mishandle the
   two extra pieces (Barrier B). The δ-resolution itself
   (`subPolys`/`atildeRaw`/`rcMin`, lines 1038–1062, and the `z0^(ep/mp)`
   coefficients) is purely symbolic in the exponents and already ε-agnostic — it
   works unchanged with ε-carrying `B_j`. **Only the convergence *decision*
   (candidate filter, ranking, FlattenSector flag) is not ε-aware.**

So: make the lifted path ε-aware enough to *emit* a divergent lifted sector
(Barrier A), teach the divergence routines the two extra pieces of metadata
(Barrier B), and let the driver route it (Barrier C). The pole machinery itself
is reused verbatim.

---

## 3. The one genuine subtlety — domain indicator × divergent variable

The lifted domain indicator is `Boole[log_ypstar ≤ 0]` with
```
log_ypstar = (logZ0 − Σ_i ic_i · log y'_i) / mp,   ic_i = mOther_i / ãtilde_i   (lines 1169–1185)
```
This is a half-space cut on the y′ coordinates. The 1/ε pole is the `y_k → 0`
endpoint of the divergent variable. There are two cases:

- **Case A — clean (implement in full).** Either `DomainConstraint === None`, **or**
  the divergent variable does not enter the indicator (`ic_k = mOther_k/ãtilde_k
  = 0`, i.e. `mOther_k = 0`). Then the indicator is constant in `y_k`; it
  commutes with the `y_k` IBP / subtraction and simply rides along on the
  remaining variables: G0 (at `y_k=0`) / boundary (at `y_k=1`) / IBP terms /
  remainder all keep the same `Boole[...]` evaluated on `y_{≠k}` (and at the
  fixed `y_k` value for boundary/G0). Mathematically identical to the unlifted
  pole extraction, just multiplied by a `y_k`-independent indicator.

- **Case B — hard (detect and refuse, `$Failed`).** `ic_k ≠ 0`: the divergent
  endpoint and the domain face couple. Either (i) the `y_k→0` region is cut off
  by the indicator (one sign of `ic_k·mp`) so there is **no real pole** — the
  "divergence" is spurious on the true domain — or (ii) the `y_k` integral has a
  `y_k`-dependent face and IBP in `y_k` picks up an extra boundary term, so the
  naive `(1/a_k)[B − Σ coeff·I_t]` formula is wrong. The lifting code already
  anticipates this: the `liftnopivot` message text says *"the domain constraint
  may cut off all divergent regions (log-space remap, future work)."* The honest
  v3 move is a new message `TropicalEval::liftdivdomain` and a clean `$Failed`,
  with the log-space remap left as documented future work.

**Pivot preference helps.** Because the candidate ranking is free to choose the
pivot, we add a tie-break that *prefers* a pivot giving `mOther_k = 0` on the
divergent slot (Case A) over one that couples (Case B). Many sectors that look
like Case B under one pivot are Case A under another. Only sectors that are Case
B for *every* admissible pivot get refused.

---

## 4. Design — minimal, exactness-preserving data-flow changes

The change set is deliberately small and **byte-compatible for every existing
(unlifted, real) path** so codegen golden #25 stays intact.

### 4.1 Introduce a `"PrefactorBase"` `SectorData` field (fixes Barrier B.1)
- `ProcessSector` (M1) and the unlifted convergent path: set
  `"PrefactorBase" -> Abs[detM]`.
- `ProcessSectorLifted`: set `"PrefactorBase" -> (Abs[detM]/Abs[mp])·z0^(ap/mp−1)`
  (the value already computed as `prefactorBase` at line 1192).
- Replace `Abs[detM]` with `Lookup[sectorData,"PrefactorBase",Abs[detM]]` at the
  ~6 reconstruction sites in §1/Barrier B.1 (`ProcessDivergentSector` 1510/1571;
  `IBPProcessSector` 4590/4656; `ValidateSubtraction`/`ValidateIBP`). The
  `Times@@(non-div aVals)` flatten-Jacobian is unchanged (it is the Jacobian of
  flattening the surviving variables and is independent of lifting).
- **Invariance:** for unlifted sectors `PrefactorBase = Abs[detM]`, so the value
  is identical and the emitted C++ is byte-for-byte unchanged. The codegen reads
  the prefactor from the already-computed `SectorData` keys (`G0Prefactor`,
  `Remainder["Prefactor"]`, `BoundaryData["Prefactor"]`, term `Prefactor`), so
  **no codegen change is needed for the prefactor** — only the WL processing
  functions change.

### 4.2 Thread `DomainConstraint` through the divergence integrands (Barrier B.2)
Only relevant when `DomainConstraint =!= None` and (Case A) `ic_k = 0`.
- `emitBaseFuncBody` already accepts `domainConstraint` and emits the indicator
  block (2122–2143). Change the *calls* for G0/G1/IBP-boundary/IBP-term to pass
  `Lookup[sd,"DomainConstraint",None]` instead of letting it default to `None`
  (sites: 2306–2307, 2325–2326, 2449–2450, 2493–2494, 2506–2507). For G0 and the
  IBP boundary the indicator is evaluated at the fixed `y_k` (0 and 1
  respectively); since `ic_k=0` in Case A the indicator does not reference
  `y_k`, so the existing `log_y[]`-over-`(n−1)`-dims emission is already correct
  — **but the indicator's `IndicatorCoeffs` must be re-indexed to the
  non-divergent coordinate ordering** (drop slot `k`). Provide a
  `dropDivVarFromDomain[dc, k]` helper that removes the `k`-th `IndicatorCoeffs`
  entry and (in Case A) asserts it was 0.
- The hand-written **remainder** integrand (2353–2424) needs the indicator block
  inserted manually (it is full-`n`-dim, so the indicator is emitted over all `n`
  `log_y[]` unchanged).
- `ValidateSubtraction` / `ValidateIBP` NIntegrate reconstructions must multiply
  by the `Boole[log_ypstar ≤ 0]` factor, copying the pattern from
  `ValidateLiftedDecomposition` (1328–1334).

### 4.3 Make `ProcessSectorLifted` ε-aware and able to emit a divergent sector (Barrier A)
Add an optional `eps` argument: `ProcessSectorLifted[..., liftData, eps:None]`
(mirrors `FlattenSector[..., eps:None]`, plan.md §6.1).

- **ε-aware classification of each `ãtilde_j`** (replaces `isRealPos`), decided
  *exactly* per §6.7 — ε-expand symbolically, never an `N`/tolerance test:
  ```
  a0_j = ãtilde_j /. eps -> 0
  class_j = "complex"     if  ¬PossibleZeroQ[Im[a0_j]]            (→ liftcomplex)
            "divergent"   if  PossibleZeroQ[Im[a0_j]] ∧ TrueQ[Re[a0_j] <= 0]
            "convergent"  if  PossibleZeroQ[Im[a0_j]] ∧ TrueQ[Re[a0_j] >  0]
  ```
- A pivot is **admissible** iff no direction is `"complex"` and **at most one**
  is `"divergent"` (mirrors the single-pole limit L1/L8). For a divergent
  direction additionally require `c_k = (D[ãtilde_k,eps] /. eps→0) ≠ 0`
  (`badck`).
- **Candidate ranking** becomes ε-safe: rank on `a0 = ãtilde /. eps→0` (use the
  leading real parts), keep the existing keys (HasConstantTerm, `|mp|=1`,
  `max min Re`), and add the **new tie-break: prefer pivots with `mOther_k = 0`
  on the divergent slot** (Case A over Case B, §3).
- **Domain-coupling gate (Case B refusal):** if the chosen divergent direction
  has `ic_k ≠ 0` (`mOther_k ≠ 0`) for *all* admissible pivots, emit
  `TropicalEval::liftdivdomain[coneIndex]` and `Return[$Failed]`.
- **Do not call `FlattenSector` blindly on the divergent direction.** When the
  sector is convergent (no divergent direction), behave exactly as today
  (byte-identical). When it has one divergent direction, **skip the
  `IsDivergent`→`liftdivergent` abort** and instead emit a *divergent lifted*
  `SectorData`:
  ```
  "IsDivergent"       -> True,
  "DivergentVariable" -> k,            (* index into the n surviving coords *)
  "NewExponents"      -> ãtilde,        (* ε-carrying; pole machinery ε-expands *)
  "ClearedPolys"      -> reclearedPolys,
  "PrefactorBase"     -> (Abs[detM]/Abs[mp])·z0^(ap/mp−1),
  "DomainConstraint"  -> domainClass,   (* Case A: ic_k = 0 *)
  ...all existing lifted keys (PivotIndex, ZRow, AugmentedA, HasConstantTerm,
     MonoFactorLog, LiftData)...
  ```
  The flattening of the **non-divergent** coordinates is then done downstream by
  the divergence routines (they flatten the surviving variables themselves), so
  `ProcessSectorLifted` does not pre-flatten in the divergent case.

### 4.4 Driver routing (Barrier C)
`EvaluateTropicalMC` (M4):
- Drop the `isLifted, False` short-circuit in `routeToIBP` (3581–3583).
- In lifted mode, process each sector with `ProcessSectorLifted[..., eps]`,
  **partition** the results into `convergentSectors` and `divergentSectors`
  (instead of hardcoding `divergentSectors = {}`), keeping `EmptyDomain`
  handling.
- If any divergent lifted sectors exist and `Method` resolves to a divergence
  method, route the **divergent** lifted sectors through `IBPProcessSector` /
  `ProcessDivergentSector` (now `PrefactorBase`/`DomainConstraint`-aware) and the
  **convergent** lifted sectors through the existing convergent codegen — the
  same conv+div coexistence the unlifted IBP path already supports.
- The Laurent assembly (pole, finite) is the existing IBP/subtraction assembly,
  unchanged — it sums sector contributions and is agnostic to whether the sector
  was lifted.
- `EvaluateTropicalMCLifted`: stop being convergent-only; pass `Method` through
  (`Automatic` resolves to IBP when divergent lifted sectors are present, exactly
  like the unlifted driver).
- Keep the single-pole / nested guard (`nestedIBP`) and add the new
  `liftdivdomain` guard. `>1` divergent `ãtilde` directions in a lifted sector →
  refuse (consistent with L1/L8).

### 4.5 Complex exponents (`MonoFactorLog`) — deferred sub-phase
The **Direct** complex-exponent path (`exp(B·log P)`, plan.md §6.3) is always
correct and needs no `MonoFactorLog`, so lift+divergence works for complex `B`
in Direct mode immediately. The **SplitRealImag** mode multiplies three features
(lift × divergence × oscillatory phase): the lifted `MonoFactorLog` (lines
1203–1225) is the *sector* monomial factor, but IBP shifts exponents per term and
re-flattens with each term's `alpha0`, so the phase's `MonoFactorLog` must be
re-derived per IBP term/boundary. **Scope:** implement and gate Direct mode in
the main phases; defer SplitRealImag×lift×divergence to a follow-up sub-phase
and, until then, refuse it with a clear message
(`TropicalEval::splitliftdiv`) rather than emitting a silently-wrong phase.

---

## 5. Exactness discipline (plan.md §1.4 / §6.7) — must hold

- All new classification decisions are **exact**: `PossibleZeroQ[Im[...]]` and
  `TrueQ[Re[a0]<=0]` on the ε→0 leading term; ε-expansion is symbolic
  (`/. eps->0`, `D[...,eps]`). **No `N[...]` in any decision** — the only `N` is
  the existing ranking heuristic on `a0` (numeric ordering does not affect
  correctness, only which valid pivot is chosen) and `MmaToC` numericization.
- `PrefactorBase = (Abs[detM]/Abs[mp])·z0^(ap/mp−1)` stays **exact** (`z0`,
  radicals, `Abs[detM]` integer) until the `MmaToC` boundary.
- Guard test #43 (`FreeQ[symbolicOutput, _Real]` before `MmaToC`) must pass for a
  lifted+divergent exact-input spec — add such a spec to the #43 corpus.

---

## 6. Phased implementation plan (with gates)

Sandbox-first, mirroring plan.md §7. Each gate is judged by ≥2 independent
oracles (plan.md §8.5); no gate is judged by the routine that produced it.

**Phase L0 — Sandbox the math.** In `SANDBOX/`, take one explicit
lifted+divergent toy (e.g. `∫_0^∞ dx1 dx2  x1^{-1+ε} (1 + 10^6·x1·x2 + x2^2)^{-2}`):
by hand / NIntegrate, confirm (a) which sectors are divergent post-δ, (b) the
pole lives in a surviving `ãtilde` direction, (c) the assembled
`pole/ε* + finite` matches NIntegrate of the original at small `ε*`. Identify a
Case-A sector and a Case-B sector explicitly. *Gate:* the toy's Laurent
reproduced two ways (subtraction-by-hand vs NIntegrate).

**Phase L1 — `PrefactorBase` plumbing (Barrier B.1).** Add the field; substitute
at the ~6 sites. *Gate:* all existing tests still green; codegen golden #25
**byte-identical** (unlifted `PrefactorBase = Abs[detM]`).

**Phase L2 — ε-aware `ProcessSectorLifted` (Barrier A).** Add `eps`,
classification, ranking tie-break, Case-B (`liftdivdomain`) refusal, divergent
lifted `SectorData` emission. *Gate:* (a) every existing lifted/convergent test
byte-identical (convergent path untouched); (b) on the L0 toy, the divergent
lifted `SectorData` has the expected `DivergentVariable`, `c_k`, `PrefactorBase`,
`DomainConstraint`; (c) #43 exactness guard passes on a lifted+divergent exact
spec.

**Phase L3 — divergence routines consume lifted metadata (Barrier B.2) + driver
(Barrier C), Direct mode.** Thread `DomainConstraint` into G0/G1/boundary/term
+ remainder integrands and the WL validators; wire driver routing. *Gate:*
- **#new-A (IBP vs subtraction in lifted mode):** the (pole, finite) of the L0
  toy agree between IBP and Subtraction routes (the same cross-check that
  justifies keeping both methods, plan.md #4) — now lifted.
- **#new-B (lifted+divergent vs exact Laurent):** where a Γ-product closed form
  exists, pole and finite match to < 0.5%.
- **#new-C (vs NIntegrate of the original at fixed ε):** reconstructed `I(ε*)`
  vs direct NIntegrate of the *unlifted original* at `ε* ∈ {0.02, 0.01}` within
  the established tolerance.
- **#new-D (variance):** lifted vs unlifted divergent sector variance — lifting
  reduces the MC variance of the finite part for the extreme-coefficient toy
  (with the L3/HasConstantTerm caveat honored: trust the exact validator, not
  the MC σ, for `HasConstantTerm=False` sectors).
- **Case-B refusal:** a sector that couples for all pivots returns clean
  `$Failed` with `liftdivdomain` (pinned like #27).

**Phase L4 — complex exponents (Direct) + docs.** Confirm lift+divergence with
complex `B_j` in Direct mode (extend the L0 toy to complex `B`); refuse
SplitRealImag×lift×divergence with `splitliftdiv`. Update SUMMARY L2 and plan.md
N3 to "supported (Case A); Case B and SplitRealImag are future work."

**Phase L5 (follow-up, optional) — SplitRealImag × lift × divergence.** Re-derive
the per-IBP-term `MonoFactorLog`; gate with the complex-*base* cross-check (#21b
analog) so Fix A / Fix B (plan.md §6.3) are exercised.

---

## 7. Cross-checks (extend the existing suite; ≥2 oracles each)

- **Extend `TEST/cc_41.wl`.** It already tests numerator×divergence (41-A) and
  numerator×lifting (41-B) *separately, precisely because they could not be
  combined.* Add **41-C: lifting × divergence** —
  `∫_0^∞ dx1 dx2  x1^{-1+ε} (1 + 10^6 x1 x2 + x2^2)^{-2}`. Oracles: IBP route,
  Subtraction route, and NIntegrate-at-ε* of the original; PASS = pole agreement
  < 0.05 across routes and finite within 2% of NIntegrate.
- **New `TEST/cc_44.wl` (Case-B refusal):** a lifted+divergent spec whose
  divergent variable couples to the domain face for every pivot → must return
  `$Failed` with `liftdivdomain` (limitation pinned, like #27/#34).
- **#19/#18 (variance) lifted-divergent variants:** exact `σ² = I2 − I1²` per
  surviving sector; `HasConstantTerm=False` sectors flagged infinite-variance
  (the existing honest-bench discipline).
- **#43 corpus:** add a lifted+divergent exact-input spec; assert
  `FreeQ[symbolicOutput,_Real]` before `MmaToC`.
- **CUBA (#5) where available:** Cuhre on the direct integrand of 41-C as an
  independent oracle at fixed ε.

---

## 8. New messages / limitation updates

```
TropicalEval::liftdivdomain =
  "ProcessSectorLifted: cone `1` — the divergent variable couples to the lifted
   domain constraint (ic_k != 0) for every admissible pivot.  The 1/eps pole and
   the domain face interact (Case B); a log-space remap is required (future
   work).  Aborting ($Failed).";

TropicalEval::splitliftdiv =
  "EvaluateTropicalMC: SplitRealImag combined with lifting AND divergence is not
   yet supported (per-term MonoFactorLog re-derivation is future work).  Use
   ComplexExponentMode -> \"Direct\" (always correct) for lifted+divergent
   complex-exponent integrals.";
```

Update **SUMMARY.txt §2d / §6 L2** and **plan.md N3 / §9 D3** from "mutually
exclusive" to: *"Supported when the divergent variable does not couple to the
lifted domain constraint (Case A); the coupled case (Case B) and
SplitRealImag×lift×divergence remain documented future work."* Keep the
single-pole-per-sector limit (L1/L8) in force for lifted sectors too.

---

## 9. Scope boundaries (what this plan does and does not deliver)

**Delivers:** lifting + a single 1/ε pole per sector, IBP **and** subtraction
routes, Direct complex exponents, full exactness discipline, byte-identical
unlifted goldens, refusal of the geometrically-coupled case and of nested
(`>1`-divergent) lifted sectors.

**Does not deliver (explicit future work):** Case B (domain-face coupling /
log-space remap); 1/ε^{d≥2} in lifted sectors (inherits L1/L8); SplitRealImag ×
lift × divergence (Phase L5). Each is refused with a clear message, never
silently mis-evaluated — consistent with the v3 rule that a known limitation must
produce a clean `$Failed`, not a wrong number.
