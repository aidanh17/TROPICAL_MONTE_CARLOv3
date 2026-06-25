export const meta = {
  name: 'tmcv3-phase0',
  description: 'Build TROPICAL_MONTE_CARLOv3 Phase 0: v3-build worktree, skeleton, capture both old goldens, fan smoke, with independent adversarial gate',
  phases: [ { title: 'Produce', model: 'sonnet' }, { title: 'Verify', model: 'opus' }, { title: 'Commit' } ],
}

// ============================ shared harness (ORCHESTRATION.md §2) ============================
const WORK = { type:'object', required:['summary','selfGatePass','commands','numerics','artifacts'],
  properties:{
    summary:{type:'string'},
    filesChanged:{type:'array',items:{type:'string'}},
    commands:{type:'array',items:{type:'object',properties:{cmd:{type:'string'},exitCode:{type:'number'},keyOutput:{type:'string'}}}},
    numerics:{type:'array',items:{type:'object',properties:{name:{type:'string'},value:{type:'string'},claimedRef:{type:'string'},deltaSigma:{type:'string'}}}},
    selfGatePass:{type:'boolean'},
    artifacts:{type:'array',items:{type:'string'}} } }

const VERDICT = { type:'object', required:['lens','refuted','evidence'],
  properties:{
    lens:{type:'string'},
    refuted:{type:'boolean'},                                // true = I found it wrong / unproven
    confidence:{type:'string'},
    evidence:{type:'string'},
    independentChecks:{type:'array',items:{type:'object',properties:{check:{type:'string'},result:{type:'string'}}}},
    disagreements:{type:'array',items:{type:'string'}},
    notes:{type:'string'} } }

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
  const work = await agent(producePrompt, {
    label: 'produce', phase: 'Produce', model: produceModel, effort: produceEffort, schema: WORK })

  phase('Verify')
  const verdicts = (await parallel(
    verifierLenses.map(lens => () =>
      agent(verifyPrompt(lens, work), {
        label: `verify:${lens.id}`, phase: 'Verify',
        model: 'opus', effort: lens.effort || 'high', schema: VERDICT }))
  )).filter(Boolean)

  const refuted = verdicts.filter(v => v.refuted).length
  const allRan  = verdicts.length === verifierLenses.length          // dead/absent verifier = NOT confirmed
  const green   = !!(work && work.selfGatePass) && allRan && refuted === 0   // FAIL-CLOSED (§1.7)
  log(`PHASE GATE: ${green ? 'GREEN' : 'RED'} (selfGate=${work && work.selfGatePass}, refuted=${refuted}/${verifierLenses.length}, allRan=${allRan})`)

  if (green && commitMsg) {
    phase('Commit')
    await agent(`The phase gate is GREEN. Commit the v3-build worktree state. cd /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3 and on branch v3-build run: git add -A (OLD_CODE/ and INTERFILES/ are gitignored and must stay untracked) then git commit -m ${JSON.stringify(commitMsg)}. Return the commit hash. Mechanical only — do not edit any source.`,
      { label: 'commit', phase: 'Commit', model: 'sonnet', effort: 'low' })
  }
  return { green, work, verdicts }
}
// ========================================================================================

const r = await runPhase({
  produceModel: 'sonnet', produceEffort: 'medium',
  producePrompt: `Set up the TROPICAL_MONTE_CARLOv3 build (plan.md Phase 0; ORCHESTRATION.md §3 Phase 0). Working root: /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3. NEVER edit anything under OLD_CODE/.

STEP A — git worktree (ORCHESTRATION.md §1.4). The project root is NOT yet a git repo (only the OLD_CODE subtrees have their own nested .git). So:
  - cd to the root; if no .git exists, 'git init'.
  - Create .gitignore that ignores 'OLD_CODE/' and 'INTERFILES/' (OLD_CODE must NEVER be tracked or edited; INTERFILES is generated).
  - Make an initial commit on the default branch containing only plan.md, ORCHESTRATION.md, ORCH/, .gitignore.
  - Create and check out branch 'v3-build'. ALL v3 build files accumulate here at the root, alongside OLD_CODE/.

STEP B — skeleton (plan.md §3.2): create the v3 directory skeleton at the root: tropical_fan.wl (copied verbatim from OLD_CODE/TROPICAL_MONTE_CARLO/tropical_fan.wl — it is byte-identical across all trees, md5 256a8da…; verify the md5 matches), and empty/placeholder dirs TEST/, TEST/baselines/, EXAMPLES/, SANDBOX/, CC/, INTERFILES/, MANUAL/. Do NOT yet copy tropical_eval.wl (that is Phase 1).

STEP C — capture goldens (plan.md §0.2), READ-ONLY against the OLD trees, IN PARALLEL, each in its OWN WorkingDirectory, run_in_background (they EXCEED the 10-min single-Bash cap; poll to completion — do NOT silently give up on timeout, halt loudly):
  - Tree A (OLD_CODE/TROPICAL_MONTE_CARLO/): RunAllTests[] (Tests 1–18) + IBP Tests 19–22 + EXAMPLES/test_divergent_crosscheck.wl. Save full transcript -> TEST/baselines/treeA.txt. Expected reference numbers to record: B=4 pole/finite=(1/6,-1/6); B=5=(1/12,-1/8).
  - Tree B (OLD_CODE/TROPICAL_MONTE_CARLO2/TROPICAL_MONTE_CARLOv2/): run_validation_suite.wl (Tests 1,2,3v2,5,6,7) + EXAMPLES/test_lifted.wl (Test 23) + EXAMPLES/test_vegas.wl (V1–V4,L1–L3). Save full transcript -> TEST/baselines/treeB.txt. Expected: Test 6A=7.850075e-4; lifted Case A VarRed ≈ 17.95×.
  Run scripts per plan.md §10: cd root && wolframscript -file <script>; use absolute paths when Get-ing packages; give each parallel runner its own INTERFILES working dir to avoid the fixed-filename clobber hazard (§10).

STEP D — fan smoke (plan.md §0.3): in the v3 location, Get the v3 tropical_fan.wl and run PolytopeVertices + ComputeDecomposition on (1+x1^2+x2^2)^-1. Record ray/sector counts.

SELF-GATE (set selfGatePass true ONLY if all hold): both baseline transcripts saved and contain the documented reference numbers (Tree A B=4/B=5 Laurent; Tree B Test 6A); fan smoke passes with sensible ray/sector counts; git branch v3-build exists with OLD_CODE/ untracked. Return WORK: filesChanged, every command with exitCode+keyOutput, numerics (each captured reference number with value+claimedRef), and artifacts (paths to treeA.txt, treeB.txt, .gitignore).`,
  verifierLenses: [{
    id: 'baseline', effort: 'high',
    instruction: 'Confirm both golden transcripts exist and the key reference numbers match the published OLD_CODE SUMMARYs; confirm the v3-build branch exists with OLD_CODE/ untracked and tropical_fan.wl byte-identical to the OLD copy. Re-run the fan smoke yourself.',
    gate: 'plan.md Phase 0 / Gate 0: both baseline transcripts saved with documented reference numbers; fan smoke passes.',
    exactChecks: 'md5 of v3 tropical_fan.wl === OLD_CODE/TROPICAL_MONTE_CARLO/tropical_fan.wl; fan ray/sector counts for (1+x1^2+x2^2)^-1 reproduce on a fresh run; git ls-files shows nothing under OLD_CODE/.',
    oracles: 'OLD_CODE/TROPICAL_MONTE_CARLO/SUMMARY.txt (Tree A B=4=(1/6,-1/6), B=5=(1/12,-1/8)) and OLD_CODE/TROPICAL_MONTE_CARLO2/TROPICAL_MONTE_CARLOv2/SUMMARY.txt (Tree B Test 6A=7.850075e-4, Case A VarRed≈17.95×). Refute if any transcript is missing or any number is off beyond stated tol.'
  }],
  commitMsg: 'Phase 0: v3-build skeleton + tropical_fan.wl + captured Tree A/B goldens (green gate)',
})
return r
