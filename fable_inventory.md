# fable_plan Step 1 — INVENTORY of tropical_eval.wl

Date: 2026-07-02. Branch `ibpcx-lifted-divergent-regulator-real` @ 1a84335.
Target file: `tropical_eval.wl` (8,766 lines). Note: fable_plan names the target
`trop_eval.wl` and guesses it is a min-plus/max-plus semiring package — no such
file exists; `tropical_eval.wl` (tropical-fan Feynman/Euler integral engine) is
the only tropical evaluation package in the repo and is taken as the target.
The min-plus semantics warning translates here to the real invariants:
symbolic exactness (no floats upstream of codegen) and byte-identical C++
emission on the real/Direct path.

## 1. Public API (26 symbols, all `::usage`-exported before Begin["`Private`"])

ProcessSector, CheckFlatteningMagnitude, ValidateDecomposition,
IdentifyDivergences, ProcessDivergentSector, ValidateSubtraction, MmaToC,
GenerateCppMonteCarlo, EvaluateTropicalMC, LaurentFromSubtraction, RunAllTests,
ParsePolynomial, CompileCpp, detectCuba, IBPReduceSector, IBPProcessSector,
IBPCheckBoundary, GenerateCppMonteCarloIBP, ValidateIBP, EvaluateTropicalMCIBP,
DetectExtremeCoefficients, LiftCoefficients, ProcessSectorLifted,
ValidateLiftedDecomposition, EvaluateTropicalMCLifted (+ $vegasOptionDefaults
et al. are private).

## 2. Module map (line ranges @ 1a84335)

| Region | Lines | Key symbols (public in caps context) |
|---|---|---|
| Package header + usage exports + messages | 1–232 | 21 messages TropicalEval::* |
| M0 primitives | 233–320 | ParsePolynomial, TransformExponents (priv), FlattenSector (priv) |
| M1 ProcessSector | 322–752 | ProcessSector (~200 ln), CheckFlatteningMagnitude, ValidateDecomposition |
| M1b lifting helpers | 758–1023 | makeAuxVar (priv), DetectExtremeCoefficients, LiftCoefficients (~152 ln) |
| M1c lifted sectors | 1026–1565 | ProcessSectorLifted (~420 ln, 6 local fns), ValidateLiftedDecomposition |
| M2 subtraction | 1568–2087 | IdentifyDivergences, ProcessDivergentSector (~190 ln), ValidateSubtraction (~250 ln) |
| M3 codegen | 2090–4055 | MmaToC/mmaToCInternal (15 clauses), GenerateMonomialSumCpp, $vegasOptionDefaults, resolveVegasSizing, normalizeIntegrator, imag*/realify* (4), emitDomainIndicatorCpp, dropDivVarFromDomain, liftedDomainBooleWL, emitBaseFuncBody, emitLogTail, emitIntegrandDefinitions (~465 ln), emitMain (~705 ln, 6 arms), cppBadHeads/cppBadTokens, GenerateCppMonteCarlo, detectCuba, CompileCpp |
| M4 driver | 4058–4960 | EvaluateTropicalMC (~694 ln, 12 steps), EvaluateTropicalMCLifted, LaurentFromSubtraction |
| M2b IBP | 4963–6366 | IBPExpandOneVariable, ibpImagPole, ibpDivClass, dropVarsFromDomain, IBPReduceSector, ibpProcessLeaf, IBPBuildCorners, IBPProcessSectorMultiDiv, IBPProcessSector (~350 ln), IBPCheckBoundary, ValidateIBP (~210 ln), GenerateCppMonteCarloIBP |
| M4b IBP driver | 6370–6909 | evaluateTropicalIBPDriver (priv, ~490 ln), EvaluateTropicalMCIBP (thin public wrapper) |
| M5 validation suite | 6913–8758 | RunAllTests + RunTest1..RunTest18 (priv) |
| End | 8760–8766 | End[], EndPackage[] |

Dependency: `tropical_fan.wl` (502 ln, TropicalFan` context, Polymake bridge) —
read in full; NOT a refactor target.

## 3. Test infrastructure (how to run)

- **Aggregate CI**: `wolframscript -file TEST/run_all.wl` — Phase A = in-package
  RunAllTests[] (18 tests); Phase B = every TEST/cc_*.wl (54 files, cc_1..cc_55
  with gaps 10,11,12,24) each in an isolated subprocess/INTERFILES dir.
  Exit 1 on any FAIL; SKIP allowed (Tier-2 CUBA / Tier-3 FIESTA).
  CLAUDE.md's "no aggregate runner" note is stale.
- **Fast gates** (~1 min total, both Exit[1] on FAIL):
  `TEST/phase2_selfgate.wl` — 9 codegen goldens byte-identity
  (TEST/baselines/codegen_goldens/{conv_nokat,subdiv,ibpdiv}×{mc,vegas_perkp,vegas_batch}.cpp),
  exactness guard FreeQ[...,_Real], FlattenSector eps-threading.
  `TEST/phase3_selfgate.wl` — lifting exactness, z→z0 round-trip, error routing
  (liftbadvar/liftcomplex/liftnopivot), exact ValidateLiftedDecomposition gate,
  domain-indicator C++ presence.
- Also byte-identity: cc_25 (goldens), cc_31 (treeA.txt structural), cc_5 (treeB.txt);
  cc_37/42/44/52 self-consistency byte checks.
- Tier-2 (CUBA, SKIP if absent): cc_5,6,7,9,18,21b,23,27,33,34,36,39,42,50.
  Tier-3 (FIESTA): cc_8 only.
- NOT in run_all: phase2/phase3 selfgates, phase5_*, test_lifted, test_complex_phase4
  — run separately.
- **Baseline status 2026-07-02**: phase2 ALL PASS (43s), phase3 ALL PASS (19s),
  cc_55 all PASS; full run_all.wl baseline in progress (log:
  scratchpad/run_all_baseline.log).

## 4. Coverage verdicts (public API)

Well-covered: ProcessSector, ValidateDecomposition, ProcessDivergentSector,
MmaToC, GenerateCppMonteCarlo, EvaluateTropicalMC(+IBP/+Lifted),
LaurentFromSubtraction, DetectExtremeCoefficients, LiftCoefficients,
ProcessSectorLifted, detectCuba, RunAllTests.
Covered: ValidateSubtraction, ParsePolynomial, CompileCpp,
GenerateCppMonteCarloIBP, ValidateLiftedDecomposition.
Thin: IdentifyDivergences (cc_55 only), IBPCheckBoundary (cc_55 only),
ValidateIBP (cc_30 only), IBPReduceSector (cc_14, cc_31),
CheckFlatteningMagnitude (cc_29 only).
No public symbol fully uncovered. cc_55 (untracked, passes) is a prior
fable_plan session's characterization test for the two zero-coverage symbols.

## 5. Refactor candidates (found during inventory)

Dead code / stale docs:
- D1 ValidateIBP 6186–6195: dead scaffolding Module referencing nonexistent
  ``ibpData`Private`dummyNotUsed`` — delete.
- D2 makeAuxVar 785: no-op `h = First[heads]` reassignment.
- D3 Unused locals: ProcessDivergentSector `aVals`; ValidateSubtraction
  `rawAVals`, `minExps`; RunTest3 `eps` binding.
- D4 Stale comments: 2242–2266 "GenerateCpp ONE public entry" (no such symbol);
  RunAllTests usage says "17 tests ... defined in tropical_eval_examples.wl"
  (actually 18, in-package).
- D5 resolveVegasSizing unused `kind_String` arg (doc-only; keep or prune with
  all 2 call sites).

Duplication (WL-side, no codegen bytes involved):
- W1 RunAllTests: 18 copy-pasted Module blocks → data-driven loop.
- W2 ValidateSubtraction: flat-poly evaluation block copy-pasted 5×; G1 block ≈
  G0 block × logInsertionSum; extract shared integrand builder.
- W3 ValidateIBP: 3 near-identical NIntegrate builders (+ D1 dead block).
- W4 ProcessDivergentSector: remFlatFullPolys/remFlatSimpPolys/G0-flat loops
  share one flatten shape; ndVars recomputed 3–7×.
- W5 Same flat-poly evaluator reappears in CheckFlatteningMagnitude (2×),
  ValidateDecomposition, RunTest12/14 — one private helper.
- W6 Divergent-test boilerplate shared by RunTest8/10/11/12/14/15/16/17/18.
- W7 Divergence predicate `TrueQ[Re[a0]<=0] || (NumericQ && Re<=0)` idiom
  repeated across FlattenSector/IdentifyDivergences/IBP paths (site-by-site
  equivalence check needed).
- W8 rk/ck ratio idiom 4× in M2b; SectorData assoc keys duplicated between
  divergent/convergent branches (M1, M1c) — preserve key order if merged.
- W9 EvaluateTropicalMC: SplitRealImag/realify/eps-pin block duplicated between
  lifted/unlifted branches; splitliftdiv guard copy-pasted 2×; validation-loop
  scaffolding 3×; kinematic Export 2×. Automatic-routing detector runs
  ProcessSector[Lifted] over all sectors then Step 2 re-processes (2× cost).

Codegen-region duplication (BYTE-IDENTITY-CRITICAL — every step gated on
goldens):
- C1 emitMain 6 arms: argv/usage parse, kin-data read, ofstream precision(17)
  tail near-verbatim ×6; Welford loop ×3 (+ CUBA accumulate ×2).
- C2 cubaBatch byte-identical ×2; cubaWrap+vegasOne near-identical ×2 (comment
  differs — extraction must keep emitted bytes identical, comment included).
- C3 remainder full_prod/simp_prod Do-blocks (2772–2800) differ only by name.
- C4 IBP base+log pair emission shape ×3 (MultiDiv corners, boundary, terms).
- C5 M2b: per-piece flatten+LogInsertions construction ×3 (ibpProcessLeaf was
  written as the generalized version; single-pole boundary/terms still
  hand-roll it).

Oversized functions (split into private stage helpers):
EvaluateTropicalMC 694 ln; emitMain 705 ln; evaluateTropicalIBPDriver 490 ln;
emitIntegrandDefinitions 465 ln; ProcessSectorLifted 420 ln; IBPProcessSector
350 ln; ValidateSubtraction 250 ln; ValidateIBP 210 ln; ProcessSector 200 ln.

## 6. Risk register (things a refactor must NOT change)

- R1 mmaToCInternal number/spacing formats (`n.0`, `(a.0/b.0)`, `" + "`,
  `" * "`, `std::pow(b, n.0)`) — byte-load-bearing for all 9 goldens.
- R2 Emission ordering conv→g0→g1→remainder / IBP boundary→terms→corners
  defines integrand_table[]; allFuncNames is regex-rebuilt from emitted text.
- R3 ProcessSectorLifted pivot ranking SortBy keys (caseA first) — correctness
  + byte-identity; do not reorder.
- R4 Exactness: all convergence/divergence decisions symbolic
  (TrueQ/PossibleZeroQ); no floats upstream of MmaToC. Float use exists BY
  DESIGN in: LiftCoefficients primary-rule ranking, classifyDomainFor
  EmptyDomain gate, diagnostics.
- R5 Return[] inside Module-inside-Table does NOT abort the outer function —
  the code works around this twice in IBPProcessSector (SelectFirst hoist;
  Catch/Throw "ibpverify"). Extracting blocks that contain Return[...] changes
  control flow — every extraction must convert to explicit value-check
  plumbing. (Known latent instance: ValidateLiftedDecomposition 1510 can leak
  $Failed into sectorResults silently — characterize, don't "fix" silently.)
- R6 Association key insertion order is observable (SameQ); preserve when
  merging assoc builders.
- R7 ibpImagPole's Refine[..., Reals] is load-bearing for the splitdivmono
  guard; do not simplify.
- R8 Options[...] defaults and $vegasOptionDefaults splice order feed emitted
  numbers (epsrel/epsabs/seed) — formatting via ToString[CForm[N[...]]] fixed.
- R9 Known limitation ⇒ Message + $Failed (21 message table entries) — every
  refusal path must keep firing the same message.
- R10 evaluateTropicalIBPDriver mcCx row indexing (kpOffset + fid − NConvergent
  + 2) is brittle to codegen layout: keep codegen ordering and driver parsing
  in lockstep.

## 7. Gating protocol for every refactor step

1. `wolframscript -file TEST/phase2_selfgate.wl` → ALL PASS (goldens byte-identical).
2. `wolframscript -file TEST/phase3_selfgate.wl` → ALL PASS.
3. Targeted cc tests for the touched region (e.g. cc_25/cc_31 for codegen,
   cc_30/cc_33 for subtraction/IBP numerics, cc_55 for IdentifyDivergences).
4. New characterization tests (cc_55, cc_56+) → PASS.
5. Milestone (end of each cluster): full `TEST/run_all.wl` → exit 0.
No step is accepted with any FAIL, skipped-that-should-run, or golden drift.
