export const meta = {
  name: 'tmcv3-phase2',
  description: 'Build TROPICAL_MONTE_CARLOv3 Phase 2: refactor core, unify codegen+drivers, fold batching, begin exactness audit, with independent adversarial gate (golden-diff + exactness lenses)',
  phases: [ { title: 'Produce', model: 'opus' }, { title: 'Verify', model: 'opus' }, { title: 'Commit' } ],
}

// ============================ shared harness (ORCHESTRATION.md §2) ============================
const WORK = { type:'object', required:['summary','selfGatePass','commands','numerics','artifacts'],
  properties:{ summary:{type:'string'}, filesChanged:{type:'array',items:{type:'string'}},
    commands:{type:'array',items:{type:'object',properties:{cmd:{type:'string'},exitCode:{type:'number'},keyOutput:{type:'string'}}}},
    numerics:{type:'array',items:{type:'object',properties:{name:{type:'string'},value:{type:'string'},claimedRef:{type:'string'},deltaSigma:{type:'string'}}}},
    selfGatePass:{type:'boolean'}, artifacts:{type:'array',items:{type:'string'}} } }
const VERDICT = { type:'object', required:['lens','refuted','evidence'],
  properties:{ lens:{type:'string'}, refuted:{type:'boolean'}, confidence:{type:'string'}, evidence:{type:'string'},
    independentChecks:{type:'array',items:{type:'object',properties:{check:{type:'string'},result:{type:'string'}}}},
    disagreements:{type:'array',items:{type:'string'}}, notes:{type:'string'} } }
function verifyPrompt(lens, work) { return `
You are an INDEPENDENT VERIFIER. You did NOT produce this result. Your job is to REFUTE it; default to refuted=true unless your own independent checks PROVE it correct. READ-ONLY: do not edit any package or build file; never edit OLD_CODE/.
Lens: ${lens.id} — ${lens.instruction}
Gate spec for this phase: ${lens.gate}
Producer's claimed artifacts/results (do NOT trust these — re-derive): ${JSON.stringify(work)}
Procedure:
 1. Re-run the gate commands yourself (cd /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3; wolframscript / g++ / diff / grep / git). Do not trust the producer's exit codes.
 2. Run the EXACT symbolic identity checks (these cannot be subtly wrong): ${lens.exactChecks}
 3. For any NUMERIC value, compare against the independent oracle(s): ${lens.oracles}. Assume the value is wrong; hunt for a reference that disagrees beyond combined error.
 4. Return VERDICT: refuted=true with concrete disagreements if ANY check fails or any numeric lacks an independent corroborator; refuted=false ONLY if every check independently confirms.` }
async function runPhase({ producePrompt, produceModel = 'opus', produceEffort = 'max', verifierLenses, commitMsg }) {
  phase('Produce')
  const work = await agent(producePrompt, { label:'produce', phase:'Produce', model:produceModel, effort:produceEffort, schema:WORK })
  phase('Verify')
  const verdicts = (await parallel(verifierLenses.map(lens => () =>
    agent(verifyPrompt(lens, work), { label:`verify:${lens.id}`, phase:'Verify', model:'opus', effort:lens.effort||'high', schema:VERDICT }))
  )).filter(Boolean)
  const refuted = verdicts.filter(v => v.refuted).length
  const allRan  = verdicts.length === verifierLenses.length
  const green   = !!(work && work.selfGatePass) && allRan && refuted === 0
  log(`PHASE GATE: ${green ? 'GREEN' : 'RED'} (selfGate=${work && work.selfGatePass}, refuted=${refuted}/${verifierLenses.length}, allRan=${allRan})`)
  if (green && commitMsg) {
    phase('Commit')
    await agent(`The phase gate is GREEN. Commit the v3-build worktree. cd /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3, branch v3-build: git add -A (OLD_CODE/ + INTERFILES/ gitignored) then git commit -m ${JSON.stringify(commitMsg)}. Return the commit hash. Mechanical only.`,
      { label:'commit', phase:'Commit', model:'sonnet', effort:'low' })
  }
  return { green, work, verdicts }
}
// ========================================================================================

const r = await runPhase({
  produceModel: 'opus', produceEffort: 'max',
  producePrompt: `TROPICAL_MONTE_CARLOv3 Phase 2 — the HARD core-refactor surgery (plan.md §7 Phase 2 + §5 + §6.1 + §6.7; ORCHESTRATION.md §3 Phase 2). Root: /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3, branch v3-build. NEVER edit OLD_CODE/. tropical_fan.wl gets NO change this phase. Re-load the package after EVERY edit (wolframscript -code 'Get["tropical_eval.wl"]; Print["LOAD OK"]'); order is law.

Read plan.md §3.1 (module layout M0–M5), §5 (refactor/deletion catalog), §6.1, §6.7 FIRST. Then, editing ./tropical_eval.wl in place:
 2.1 Introduce Module-0 primitives (ParsePolynomial, TransformExponents, MmaToC/mmaToCInternal, GenerateMonomialSumCpp) as the single shared copy (already byte-identical — pure reorganization, §2.3).
 2.2 Reconcile ProcessSector ⊕ FlattenSector with EXACT eps threading (§6.1): keep Tree A's eps threading AND Tree B's FlattenSector factorization; ProcessSector calls FlattenSector passing eps through. CRITICAL: preserve the IsDivergent pre-flattening branch (plan.md §9 risk #2 KEEP list).
 2.3 Collapse the duplicated convergent/IBP codegen towers into ONE parametrized emitter (§5.2): one emitIntegrandDefinitions(type-tag) + one emitMain(MC|VEGAS). FOLD emitVegasBatchMain into the VEGAS emitter's Batch mode (§5.1). Unify drivers into one EvaluateTropicalMC routing on Method/Lift/Integrator + thin wrappers (§5.4). De-duplicate validation + log-exp eval (§5.5, extract evalFlattenedIntegrand). Single memoized detectCuba; DELETE findCubaPrefix (§5.3). Share option defaults (§5.7). Adopt integrator vocabulary "MC"/"VEGAS" (§4.3) and Tree A's hard badcpp gate.
 2.4 Begin the exactness audit (§6.7): ensure symbolic decomposition output is FreeQ[_,_Real] before the MmaToC boundary; remove decision-making numerics from the decomposition path.

SELF-GATE — the "simplify aggressively but PROVE nothing broke" gate (selfGatePass true ONLY if all hold):
 - golden-master codegen regression GREEN for ALL kinds {conv,g0,g1,rem,ibp-boundary,ibp-term} under {MC, per-kp VEGAS, batched VEGAS}: re-emit and diff against the Phase-1 TEST/baselines/codegen_goldens/ — clean (cross-check #25).
 - Tests 1–22 + test_divergent_crosscheck reproduce the Phase-0 numbers EXACTLY (Laurent bit-exact; numeric within MC noise).
 - Symbols GONE (grep across tropical_eval.wl, expect zero hits): RunBenchmark2D, findCubaPrefix, emitVegasBatchMain, and the IBP codegen-twin emitters folded into the unified emitter.
 - Exactness guard #43 passes on the touched paths: FreeQ[symbolic decomposition output, _Real] holds before MmaToC; no stray 10^-… tolerance in a decision path.
Long runs: run_in_background + poll; halt loudly on timeout; isolate INTERFILES per runner. Return WORK with filesChanged, all commands+exitCodes, the grep results (showing 0 hits), the reproduced numbers, and artifact paths (diff logs).`,
  verifierLenses: [
    {
      id: 'golden-diff-grep', effort: 'max',
      instruction: 'Re-emit ALL codegen kinds yourself and diff against the Phase-1 goldens; re-run Tests 1–22 + crosscheck; independently grep for the symbols that must be gone. Refute on ANY non-clean diff, ANY changed test number, or ANY surviving symbol.',
      gate: 'plan.md Gate 2: golden codegen regression green for all kinds×samplers; Tests 1–22+crosscheck reproduce Phase-0 numbers exactly; RunBenchmark2D/findCubaPrefix/emitVegasBatchMain/IBP-twin symbols gone.',
      exactChecks: 'diff(re-emitted codegen, TEST/baselines/codegen_goldens/*) is empty for every kind×{MC,per-kp VEGAS,batched VEGAS}; grep -n for RunBenchmark2D|findCubaPrefix|emitVegasBatchMain returns ZERO hits; B=4/B=5 Laurent reproduce bit-exact.',
      oracles: 'Phase-1 goldens in TEST/baselines/; Phase-0 treeA.txt; OLD_CODE Tree A SUMMARY.txt. Refute on any divergence.'
    },
    {
      id: 'exactness', effort: 'max',
      instruction: 'Independently assert the exactness invariant (§1.4/§6.7) on the touched paths. Construct your own decomposition inputs; verify FreeQ[_,_Real] holds on the symbolic output before MmaToC, and that no floating tolerance lives in a decision branch (IsDivergent, flattening, pivot, anchor).',
      gate: 'plan.md Gate 2 exactness component + cross-check #43.',
      exactChecks: 'FreeQ[symbolic decomposition output, _Real] === True at the MmaToC boundary for conv/g0/g1/rem/ibp inputs; grep the decision paths for N[…] / _Real / 10^- tolerances and confirm none gate a branch; the IsDivergent pre-flatten branch still exists and fires.',
      oracles: 'n/a — this lens is pure exact symbolic (Layer 1, §8.5). Refute if any decision numeric survives.'
    }
  ],
  commitMsg: 'Phase 2: unified core/codegen/drivers, folded batching, exactness audit begun (green gate)',
})
return r
