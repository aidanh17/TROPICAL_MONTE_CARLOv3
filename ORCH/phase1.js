export const meta = {
  name: 'tmcv3-phase1',
  description: 'Build TROPICAL_MONTE_CARLOv3 Phase 1: seed tropical_eval.wl from Tree A, reproduce Tree-A goldens, capture v3 codegen goldens, with independent adversarial gate',
  phases: [ { title: 'Produce', model: 'sonnet' }, { title: 'Verify', model: 'opus' }, { title: 'Commit' } ],
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
 1. Re-run the gate commands yourself (cd /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3; wolframscript / g++ / diff / git). Do not trust the producer's exit codes.
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
  produceModel: 'sonnet', produceEffort: 'medium',
  producePrompt: `TROPICAL_MONTE_CARLOv3 Phase 1 (plan.md §7 Phase 1; ORCHESTRATION.md §3 Phase 1). Root: /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3, branch v3-build. NEVER edit OLD_CODE/.

1.1 Copy Tree A's engine: cp OLD_CODE/TROPICAL_MONTE_CARLO/tropical_eval.wl -> ./tropical_eval.wl (verbatim — this is the canonical base, plan.md §2.1). Re-load to confirm: wolframscript -code 'Get["tropical_eval.wl"]; Print["LOAD OK"]'. Then confirm it reproduces the Phase-0 Tree-A goldens (TEST/baselines/treeA.txt): RunAllTests Tests 1–22 + test_divergent_crosscheck, in BOTH MC and VEGAS, per-kp AND batched. These are long; run_in_background and poll, halt loudly on timeout. Give each runner its own INTERFILES working dir (§10 clobber hazard).

1.2 Capture v3 codegen goldens — the byte-level regression anchor for ALL later refactors (cross-check #25). For EVERY integrand kind {conv, g0, g1, rem, ibp-boundary, ibp-term} emit the C++ under {MC, per-kp VEGAS, batched VEGAS} and save each generated source under TEST/baselines/codegen_goldens/. Use OLD_CODE/TROPICAL_MONTE_CARLO/TEST/golden_capture.wl as the pattern (do not edit it; copy/adapt into TEST/).

SELF-GATE (selfGatePass true ONLY if all hold): package loads ("LOAD OK"); Tests 1–22 + crosscheck reproduce the Phase-0 Tree-A numbers exactly (symbolic/Laurent) and within MC noise (numeric); all codegen goldens captured for every kind×sampler. Return WORK with commands+exitCodes, the reproduced reference numbers vs Phase-0, and artifact paths.`,
  verifierLenses: [{
    id: 'tree-a-regression', effort: 'high',
    instruction: 'Re-run Tests 1–22 + test_divergent_crosscheck yourself against the freshly-seeded v3 tropical_eval.wl; diff a SAMPLED set of the captured codegen goldens against a fresh re-emit; refute on ANY mismatch.',
    gate: 'plan.md Gate 1: package loads; Tests 1–22 + crosscheck reproduce Phase-0 Tree-A numbers; all goldens captured.',
    exactChecks: 'Laurent/symbolic outputs (B=4=(1/6,-1/6), B=5=(1/12,-1/8)) reproduce bit-exact; a re-emitted codegen golden is byte-identical (diff clean) to the saved one; v3 tropical_eval.wl is byte-identical to OLD_CODE/TROPICAL_MONTE_CARLO/tropical_eval.wl (it must be a verbatim copy at this phase).',
    oracles: 'TEST/baselines/treeA.txt (Phase-0 golden) and OLD_CODE/TROPICAL_MONTE_CARLO/SUMMARY.txt. Refute if any test number diverges beyond its documented tol or any golden differs.'
  }],
  commitMsg: 'Phase 1: seed tropical_eval.wl from Tree A; v3 codegen goldens captured (green gate)',
})
return r
