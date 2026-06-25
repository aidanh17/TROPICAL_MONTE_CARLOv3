export const meta = {
  name: 'tmcv3-phase4',
  description: 'Build TROPICAL_MONTE_CARLOv3 Phase 4: complex exponents end-to-end + complex lifting + SplitRealImag (BUG-1/L3 zone), gated by a UNANIMOUS 3-lens adversarial panel',
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
You are an INDEPENDENT VERIFIER on a 3-LENS PANEL guarding the phase the user CANNOT eyeball (the BUG-1/L1/L3 silent-wrong-answer zone). You did NOT produce this result. REFUTE it; default to refuted=true unless YOUR OWN independent checks PROVE it correct. The phase is GREEN only if ALL THREE lenses confirm — your single refutation halts the run. READ-ONLY; never edit OLD_CODE/.
Lens: ${lens.id} — ${lens.instruction}
Gate spec for this phase: ${lens.gate}
Producer's claimed artifacts/results (do NOT trust these — re-derive): ${JSON.stringify(work)}
Procedure:
 1. Re-run the gate commands yourself (cd /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3; wolframscript / g++ / diff / git). Do not trust the producer's exit codes.
 2. Run the EXACT symbolic identity checks (cannot be subtly wrong): ${lens.exactChecks}
 3. For any NUMERIC value, compare against the independent oracle(s): ${lens.oracles}. Assume the value is wrong.
 4. Return VERDICT: refuted=true with concrete disagreements if ANY check fails or any numeric lacks an independent corroborator; refuted=false ONLY if every check independently confirms.` }
async function runPhase({ producePrompt, produceModel = 'opus', produceEffort = 'max', verifierLenses, commitMsg }) {
  phase('Produce')
  const work = await agent(producePrompt, { label:'produce', phase:'Produce', model:produceModel, effort:produceEffort, schema:WORK })
  phase('Verify')
  const verdicts = (await parallel(verifierLenses.map(lens => () =>
    agent(verifyPrompt(lens, work), { label:`verify:${lens.id}`, phase:'Verify', model:'opus', effort:lens.effort||'high', schema:VERDICT }))
  )).filter(Boolean)
  const refuted = verdicts.filter(v => v.refuted).length
  const allRan  = verdicts.length === verifierLenses.length   // UNANIMOUS panel: every lens must run AND confirm
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
  producePrompt: `TROPICAL_MONTE_CARLOv3 Phase 4 — the HARDEST phase: complex exponents end-to-end + complex lifting + the opt-in SplitRealImag VEGAS mode (plan.md §7 Phase 4 + §6.3; ORCHESTRATION.md §3 Phase 4). Root: /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3, branch v3-build. NEVER edit OLD_CODE/. The advanced complex-lifting work is NOT source anywhere — it is SPECIFIED ONLY in the bug logs OLD_CODE/TROPICAL_MONTE_CARLO/TMCv2_BUG_LOG.md and lift_error_log.md (plan.md §2.2). Read §6.3 and BOTH bug logs FIRST; treat them as the design spec. Keep exponents/coeffs exact through the symbolic pipeline; re-load after every edit.

 4.1 Confirm the DEFAULT complex path exp(B·log P) with complex A and B works through the unified codegen for Tree-B Examples 2,15,16,19 and exact-Gamma Example 21 (no SplitRealImag yet).
 4.2 Implement complex-COEFFICIENT lifting: anchor z0=|C|^(1/k), residual c=C/z0^k (Tree B's residual mechanism extended); verify on a complex-coefficient toy.
 4.3 Implement the opt-in "SplitRealImag" VEGAS mode with BOTH §6.3 fixes: (i) complex log P (BUG-1) and (ii) MonoFactorLog (L3). Add the MonoFactorLog SectorData key in BOTH ProcessSector AND ProcessSectorLifted. The phase block in codegen is emitted ONLY when Im(B)≠0.

SELF-GATE (selfGatePass true ONLY if all hold):
 - complex exponents match exact Gamma/Beta (#2) and NIntegrate (#1) to old tolerances.
 - the complex-BASE cross-check #21b PASSES: a polynomial with complex coefficients (arg P≠0) and complex B, lifted, where SplitRealImag ≡ Direct ≡ reference. This is the case that exposed BUG-1; without it the fix is unverified.
 - real-exponent codegen is BYTE-IDENTICAL to the Phase-2 goldens (phase block only when Im(B)≠0) — re-emit and diff.
Long runs: run_in_background + poll, halt loudly on timeout; isolate INTERFILES per runner. Return WORK with filesChanged, commands+exitCodes, every numeric (complex vs Gamma/Beta, #21b SplitRealImag vs Direct vs reference) with value+claimedRef+deltaSigma, artifacts.`,
  verifierLenses: [
    {
      id: 'A-exact-symbolic', effort: 'max',
      instruction: 'Independently verify the complex-split identity symbolically. Confirm prod P^(Re B) · exp(i·Σ Im(B)·log P) reproduces prod P^B symbolically at sampled points, and that the MonoFactorLog linear form is correct. Pure exact; zero stats.',
      gate: 'plan.md Gate 4 + §8.5 Layer 1; ORCHESTRATION.md §3 Phase 4 lens A.',
      exactChecks: 'Simplify[ ∏P^(Re B)·Exp[I·Σ Im(B)·Log P] - ∏P^B ] === 0 at several complex sample points; the MonoFactorLog key value equals the correct linear combination of monomial logs; real-exponent codegen diff vs Phase-2 goldens is empty.',
      oracles: 'n/a — pure exact symbolic. Refute if the identity fails or the linear form is wrong.'
    },
    {
      id: 'B-independent-integrator', effort: 'max',
      instruction: 'Re-integrate the #21b direct integrand (complex coefficients, arg P≠0, complex B, lifted) with CUBA Cuhre — an integrator sharing no code with the producer. Refute if SplitRealImag / Direct / your CUBA value disagree beyond combined error.',
      gate: 'plan.md Gate 4(b) + §8.5 Layer 2; ORCHESTRATION.md §3 Phase 4 lens B.',
      exactChecks: 'n/a (numeric lens).',
      oracles: 'CUBA Cuhre on the direct integrand (brew /opt/homebrew); exact Gamma/Beta for Example 21; NIntegrate. Refute if any pair disagrees beyond combined error.'
    },
    {
      id: 'C-reproducibility', effort: 'high',
      instruction: 'Verify reproducibility/stability: reseed determinism (#22 — same SeedBase → byte-identical stream; distinct seeds → distinct streams each within 5σ) and budget stability of the #21b VEGAS result. Refute on any seed non-determinism or budget-dependent drift beyond error.',
      gate: 'plan.md Gate 4 + §8.5 Layer 2 (reseed); ORCHESTRATION.md §3 Phase 4 lens C.',
      exactChecks: 'same SeedBase → byte-identical output stream (exact); distinct seeds → distinct streams.',
      oracles: 'reseeded runs; budget sweep of the #21b VEGAS value. Refute if the value drifts with budget beyond combined error.'
    }
  ],
  commitMsg: 'Phase 4: complex exponents + complex lifting + SplitRealImag (unanimous 3-lens green gate)',
})
return r
