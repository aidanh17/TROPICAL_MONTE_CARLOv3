export const meta = {
  name: 'tmcv3-phase3',
  description: 'Build TROPICAL_MONTE_CARLOv3 Phase 3: graft real-coefficient lifting sandbox-first, wire into driver, with independent adversarial gate (exact lift round-trip + independent integrator lenses)',
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
 1. Re-run the gate commands yourself (cd /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3; wolframscript / g++ / git). Do not trust the producer's exit codes.
 2. Run the EXACT symbolic identity checks (these cannot be subtly wrong): ${lens.exactChecks}
 3. For any NUMERIC value, compare against the independent oracle(s): ${lens.oracles}. Assume the value is wrong; hunt for a reference that disagrees beyond combined error. For variance claims use the EXACT trueSigma (NIntegrate sigma^2 = I2 - I1^2), NEVER the sampled sigma alone.
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
  producePrompt: `TROPICAL_MONTE_CARLOv3 Phase 3 — graft REAL-coefficient lifting, sandbox-first (plan.md §7 Phase 3 + §6.2; ORCHESTRATION.md §3 Phase 3). Root: /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3, branch v3-build. NEVER edit OLD_CODE/. Keep z0 and all lift data EXACT through the symbolic pipeline; numericize only at the codegen/NIntegrate boundary (§10). Re-load after every edit; order is law.

Read plan.md §6.2 and Tree B's lifting source FIRST: OLD_CODE/TROPICAL_MONTE_CARLO2/TROPICAL_MONTE_CARLOv2/ (tropical_eval.wl M1b/1c, SANDBOX/, AUXT/anchor_selection_procedure.md = the kStar spec).
 3.1 Port Tree B SANDBOX toys (sandbox_lift_common.wl, sandbox_toy0 1D, toy1 eps, toy2 large-coeff) into ./SANDBOX/. Run them READ-ONLY against the v3 core — this doubles as the integration proof of the §6.1 seam. HARD gate 3a: every toy PASSES exactness (relerr < 1e-3 vs exact/NIntegrate); record the variance comparison.
 3.2 Port DetectExtremeCoefficients, LiftCoefficients, ProcessSectorLifted, ValidateLiftedDecomposition (Tree B M1b/1c) into ./tropical_eval.wl. Apply: BUG-2 aux-var robustness fix; relative-tolerance liftidentity fix; SuggestedK -> kStar rename + the geometry-gate anchor selection (§6.2).
 3.3 Wire Lift into the unified EvaluateTropicalMC driver + the EvaluateTropicalMCLifted wrapper; emit the domain-indicator in codegen. Report HasConstantTerm per sector.

SELF-GATE (selfGatePass true ONLY if all hold):
 - Test 23 reproduces Tree-B numbers: (A) exactness <0.5% vs exact + the asserted DroppedSectors count; (B) end-to-end C++ <1%; (C) error paths fire cleanly (liftnopivot, etc.); (D) EmptyDomain handled.
 - bench_lift_variance reproduces Case A VarRed ≈ 18× and Case C is correctly NOT lifted.
 - HasConstantTerm reported per sector; the exact ValidateLiftedDecomposition (not the sampled 5σ) is the correctness check (plan.md §9 risk #1).
Long runs / NIntegrate: run_in_background + poll, halt loudly on timeout; isolate INTERFILES per runner. Return WORK with filesChanged, commands+exitCodes, all Test-23 + VarRed numbers (value+claimedRef+deltaSigma), artifacts.`,
  verifierLenses: [
    {
      id: 'exact-lift-roundtrip', effort: 'max',
      instruction: 'Independently verify the lift round-trip is exact (cross-check #40): z->z0 substitution reproduces the ORIGINAL polynomials symbolically, for +, -, and 2-rule multi-lift cases. Pure symbolic; zero floating point.',
      gate: 'plan.md Gate 3 (exactness component) + §8.5 Layer 1.',
      exactChecks: 'For each lift rule emitted by the producer, Simplify[ lifted /. z->z0 - original ] === 0 exactly (rel-tol only for any unavoidable float); z0^k === C exactly; DroppedSectors count and HasConstantTerm tallies match the asserted values.',
      oracles: 'n/a — pure exact symbolic. Refute if ANY round-trip is non-zero or any tally is off.'
    },
    {
      id: 'independent-integrator', effort: 'max',
      instruction: 'Re-derive Test-23A and Test-23B values with TWO independent integrators that share no code with the producer: NIntegrate AND CUBA Cuhre on the direct integrand. Refute if the lifted sector sum disagrees beyond combined error. Corroborate variance claims via EXACT trueSigma (sigma^2=I2-I1^2 by NIntegrate), never sampled sigma.',
      gate: 'plan.md Gate 3 (numeric component, §8.5 Layer 2, ≥2 independent oracles).',
      exactChecks: 'n/a (this is the numeric lens) — but confirm the producer used the exact validator, not the 5σ, as its correctness gate.',
      oracles: 'NIntegrate on the un-decomposed integrand; CUBA Cuhre (brew /opt/homebrew); Tree-B SUMMARY.txt Test-23 published values; exact trueSigma for VarRed. Refute on any disagreement beyond combined error.'
    }
  ],
  commitMsg: 'Phase 3: real-coefficient lifting grafted, sandbox-verified, wired into driver (green gate)',
})
return r
