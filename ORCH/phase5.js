export const meta = {
  name: 'tmcv3-phase5',
  description: 'Build TROPICAL_MONTE_CARLOv3 Phase 5: high-D fan K-scaling (computeFanScaled) + resolveVegasSizing/vegasbudget, gated by a UNANIMOUS 3-lens panel (the L1 silent-wrong-answer trap)',
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
You are an INDEPENDENT VERIFIER on a 3-LENS PANEL guarding the L1 silent-wrong-answer trap (a "converged" high-D value that is quietly wrong). You did NOT produce this result. REFUTE it; default to refuted=true unless YOUR OWN independent checks PROVE it correct. The phase is GREEN only if ALL THREE lenses confirm — your single refutation halts the run. READ-ONLY; never edit OLD_CODE/.
Lens: ${lens.id} — ${lens.instruction}
Gate spec for this phase: ${lens.gate}
Producer's claimed artifacts/results (do NOT trust these — re-derive): ${JSON.stringify(work)}
Procedure:
 1. Re-run the gate commands yourself (cd /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3; wolframscript / g++ / diff / git). Do not trust the producer's exit codes.
 2. Run the EXACT symbolic identity checks (cannot be subtly wrong): ${lens.exactChecks}
 3. For any NUMERIC value, compare against the independent oracle(s): ${lens.oracles}. A "converged" 8D value with NO independent corroborator is an automatic REFUTE (the L1 trap).
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
  produceModel: 'opus', produceEffort: 'high',
  producePrompt: `TROPICAL_MONTE_CARLOv3 Phase 5 — high-D hardening + anchor finalization (plan.md §7 Phase 5 + §6.4 + §6.5; ORCHESTRATION.md §3 Phase 5). Root: /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3, branch v3-build. NEVER edit OLD_CODE/. tropical_fan.wl may ONLY receive the additive computeFanScaled — no other change (§10). The sizing fix is specified in the bug log lift_error_log.md (L1/L2). Re-load after every edit.

 5.1 Port computeFanRobust -> computeFanScaled (K-scaling) into ./tropical_fan.wl (reference impl in OLD_CODE/TROPICAL_MONTE_CARLO/TEST/gen_sectors.wl — read, do not edit). Wire it into the AUTOMATIC lifted-fan path (§6.4). Additive private helper only.
 5.2 Implement resolveVegasSizing + the vegasbudget guard (§6.5) in tropical_eval.wl.

SELF-GATE (selfGatePass true ONLY if all hold):
 - an 8D lifted-with-complex-coeff + VEGAS case (analogue of bug-log Test 26) converges with DEFAULT (Automatic) sizing to ≈0.1% vs reference.
 - vegasbudget FIRES when under-budgeted (NSamples < 20·NStart) — never silently wrong.
 - low-D (n=2) VEGAS is UNCHANGED (relerr ~1e-6) — no regression.
 - K-scaled vs raw-ray fans give IDENTICAL codegen + bitwise-equal sector sum (cross-check #37).
Long runs: the 9D robust fan (~260 s) and 8D VEGAS exceed 10 min — run_in_background + poll, halt loudly on timeout; VEGAS serial; isolate INTERFILES per runner. Return WORK with filesChanged, commands+exitCodes, the 8D value + reference + deltaSigma, the n=2 relerr, and #37 diff result; artifacts.`,
  verifierLenses: [
    {
      id: 'exact-fan-scale-invariance', effort: 'high',
      instruction: 'Independently verify fan scale-invariance (#37): primitivized-vs-raw rays and K∈{1,n+2,…} must give BYTE-IDENTICAL generated C++ and bitwise-equal sector sum. Pure exact / byte-level.',
      gate: 'plan.md Gate 5 (#37 component) + §8.5 Layer 1.',
      exactChecks: 'diff(codegen at K-scaled fan, codegen at raw-ray fan) is empty; sector sums are bitwise equal; computeFanScaled changes to tropical_fan.wl are purely additive (diff vs Phase-0 tropical_fan.wl touches only the new helper).',
      oracles: 'Phase-0 tropical_fan.wl; the two producer fan codegens. Refute on any byte difference or non-additive edit.'
    },
    {
      id: 'independent-integrator-8D', effort: 'high',
      instruction: 'Re-derive the 8D reference INDEPENDENTLY via CUBA Vegas at a LARGE budget (sharing no code with the producer). The L1 trap: a "converged" 8D value with no independent corroborator is automatically refuted. Confirm vegasbudget actually fires at NSamples<20·NStart.',
      gate: 'plan.md Gate 5 (numeric component) + §8.5 Layer 2.',
      exactChecks: 'n/a (numeric lens) — but confirm vegasbudget fires (a discrete, checkable event).',
      oracles: 'CUBA Vegas at large budget; CUBA Cuhre; any closed form if available. Refute if the 8D value lacks an independent corroborator within combined error.'
    },
    {
      id: 'reproducibility-sizing', effort: 'high',
      instruction: 'Verify sizing-sweep monotonicity and low-D non-regression: as the sizing budget increases the 8D estimate must converge monotonically toward the reference (no budget-dependent bias), and n=2 VEGAS must be unchanged (~1e-6). Refute on non-monotone drift or a low-D regression.',
      gate: 'plan.md Gate 5 (stability + low-D) + §8.5 Layer 2.',
      exactChecks: 'n/a (numeric lens).',
      oracles: 'a sizing sweep you run yourself; the Phase-1/Phase-4 n=2 VEGAS golden. Refute on non-monotonicity or low-D change.'
    }
  ],
  commitMsg: 'Phase 5: computeFanScaled K-scaling + resolveVegasSizing/vegasbudget (unanimous 3-lens green gate)',
})
return r
