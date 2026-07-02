# fable_plan Step 3 — ordered refactor steps for tropical_eval.wl

Rules binding every step (from fable_plan + CLAUDE.md):
- Behavior preserved exactly; public API signatures untouched; no floats
  introduced upstream of codegen; real/Direct C++ emission byte-identical to
  the 12 goldens (9 committed + cc_57's new lifted/split goldens).
- One logical change per step; suite green after every step; a step is
  rejected on ANY fail or golden drift.
- Refusal paths keep firing the same TropicalEval::* message + $Failed.
- WL trap: `Return[expr]` inside Module-inside-Table/Do exits only the inner
  Module. Extracting any block containing Return converts control flow to
  explicit value-check plumbing — reviewed line-by-line.
- Association key insertion order preserved (SameQ-observable).
- Engine edits only when no run_all.wl is in flight.
- Each accepted step = one commit on the working branch (small reviewable
  diffs; no pushes).

Gate tiers:
- T-fast (~1 min): phase2_selfgate + phase3_selfgate (+ cc_55, cc_56 when the
  touched region warrants; both are seconds).
- T-region: the cc tests named in the step.
- T-milestone (end of each phase): full TEST/run_all.wl exit 0 + phase5/test_
  lifted/test_complex_phase4 if the phase touched their regions.

## Phase R1 — dead code and stale docs (no semantic change)
- S1. Delete ValidateIBP dead scaffolding block (≈6186–6195; references
  nonexistent ``ibpData`Private`dummyNotUsed``). Gate: T-fast + cc_30.
- S2. Remove no-op `h = First[heads]` branch in makeAuxVar (≈785); delete
  unused locals (ProcessDivergentSector `aVals`; ValidateSubtraction
  `rawAVals`, `minExps`; RunTest3 `eps`); fix stale comments (2242–2266
  "GenerateCpp" banner → describe reality; RunAllTests::usage "17 tests ...
  tropical_eval_examples.wl" → 18 tests, in-package). Gate: T-fast + cc_15 +
  cc_30 + quick RunAllTests smoke (Tests 3, 8).

## Phase R2 — WL-side dedup in validators / test suite (no codegen text)
- S3. Extract private `flatPolyValueWL[flatPolys, logY]` (the
  `Total[Table[coeff*Exp[Total[exps*logY]]]]` shape) and use it at the 5
  ValidateSubtraction sites; keep produced expressions structurally identical.
  Gate: T-fast + cc_15, cc_30, cc_33, cc_38.
- S4. Reuse the same helper in ValidateDecomposition, CheckFlatteningMagnitude
  (2 sites), ValidateIBP (3 builders), collapsing ValidateIBP's near-identical
  NIntegrate builders. Gate: T-fast + cc_29, cc_30, cc_1, cc_2.
- S5. RunAllTests: replace 18 copy-pasted blocks with a data-driven loop over
  {RunTest1→"Test 1", ...}; printed output identical. Gate: full RunAllTests
  run diffed against pre-step output (module Test-internal randomness:
  compare the verdict/summary lines).
- S6. (optional, decide after S5) Shared harness for RunTest8/10/11/12/14/
  15/16/17/18 divergent-pipeline boilerplate. Gate: full RunAllTests diff.

## Phase R3 — M2/M2b shared-shape extraction (engine, non-codegen)
- S7. ProcessDivergentSector: hoist `ndVars` (3×), extract the shared
  "flatten by a0 over live vars" loop used by G0/remFull/remSimp; drop the
  `fullPolys = clearedPolys` alias. Structural outputs must be SameQ-identical
  (cc_31 treeA pins NewExponents/Prefactor structures). Gate: T-fast + cc_15,
  cc_30, cc_31, cc_33.
- S8. Extract divergence predicate `reNonPosQ` and substitute ONLY at sites
  verified textually identical (FlattenSector ≈292; IdentifyDivergences ≈1591,
  ≈1620; M2b sites if identical). Gate: T-fast + cc_55 + cc_31.
- S9. M2b small dedups: rk/ck ratio idiom (4×) → one helper. Do NOT remove the
  "redundant" second splitdivmono check (kept; message-path preserving).
  Gate: T-fast + cc_43, cc_48, cc_49 (subset), cc_55.

## Phase R4 — codegen C++-text dedup (byte-identity critical; one chunk/step)
- S10. emitMain: extract shared C++ text builders — argv/usage parse block,
  kinematic-file read block, ofstream/precision(17) result-write tail — as
  private WL constants/functions parameterized only where arms differ.
  Emitted bytes for all goldens unchanged. Gate: T-fast (phase2 goldens) +
  cc_25 + cc_56 + cc_57.
- S11. emitMain: extract cubaBatch (×2 identical), cubaWrap/vegasOne (×2,
  comment text parameterized to preserve bytes). Gate: same as S10.
- S12. emitIntegrandDefinitions: collapse remainder full_prod/simp_prod twin
  Do-blocks (name-parameterized); extract the IBP "base func + log func" pair
  emission (3 shapes) if byte-provable; extract Welford/accumulate C++ text
  (3×/2×). Gate: same as S10 + cc_31 + milestone run_all.

## Phase R5 — oversized function decomposition (semantic-neutral splits)
- S13. Single-pole IBP boundary/terms → route through ibpProcessLeaf (the
  written-for-this generalized helper) ONLY if the produced associations are
  SameQ-identical on the cc_30/cc_31/cc_43 corpus; else keep hand-rolled and
  just extract the common sub-expressions. Gate: T-fast + cc_25, cc_30,
  cc_31, cc_33, cc_43, cc_48–cc_52.
- S14. EvaluateTropicalMC (694 ln): extract per-Step private helpers
  (validateInputs, processAllSectors, pinEpsAcrossSectors, emitAndGuard,
  compileAndRun, parseResults, crossCheck) one step at a time; Return[$Failed]
  converted to checked returns. Gate per sub-step: T-fast + cc_31, cc_33,
  cc_41, cc_47; milestone run_all after the last.
- S15. evaluateTropicalIBPDriver (490 ln): same treatment; MultiDiv assembly
  block extracted whole (no internal reordering). Gate: cc_43, cc_48–cc_52 +
  T-fast.
- S16. (assess after S14/S15) ProcessSectorLifted / IBPProcessSector /
  ProcessSector internal stage extraction; SectorData assoc base+overlay
  merge preserving key order. Gate: T-fast + cc_18/20/27/28/32/39/44/52/54
  subset + milestone run_all.

## Phase R6 — coverage pass + closeout (fable_plan step 6)
- S17. Enumerate untested branches (error-message table vs tests; empty
  inputs; degenerate cones notsimplicial/degenerate; badcpp guard; vegasbudget
  guard; liftfandim/liftdegenerate; nKP==0 and bad-shape driver refusals) and
  delegate characterization tests (cc_58+) for gaps.
- S18. Final full gate: run_all.wl + phase2/3 selfgates + phase5_* +
  test_lifted + test_complex_phase4. Update stale docs (CLAUDE.md runner
  note; SUMMARY.txt if API prose drifted). Final report.

Explicit non-goals (fable_plan goal 2 interpreted for this repo):
- NO multi-file package split: plan.md/CLAUDE.md fix "the whole engine is two
  files" as a design decision; modularity is delivered as clean in-file
  modules + named private helpers with stable seams (cc_56).
- NO renames of public symbols; no option-default changes; no reordering of
  pivot-ranking keys, emission order, or assoc key order.
- NO behavior fixes discovered en route (e.g. ValidateLiftedDecomposition's
  silent $Failed leak at ≈1510, driver's zero-fill on short result rows):
  logged for the user, characterized if cheap, NOT changed.
