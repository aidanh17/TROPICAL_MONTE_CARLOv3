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

RE-RUN NOTE (a prior Phase-4 attempt was RED). The complex-exponent / complex-coefficient-lifting / SplitRealImag PHYSICS was INDEPENDENTLY VERIFIED CORRECT (exact split identity, CUBA Cuhre on arg P≠0, #21b Direct≡Split≡reference, MonoFactorLog proven required). That work is ALREADY PRESENT in the working-tree tropical_eval.wl (5282 lines, uncommitted on top of Phase-3 commit 14b1499; MonoFactorLog + ComplexExponentMode present) and in TEST/test_complex_phase4.wl. PRESERVE it — do NOT reset, do NOT re-derive from scratch. Build the seed fix ON TOP, then re-run 4.1/4.2/4.3 to confirm they STILL pass unchanged. The refutation was a Phase-2 driver regression: the unified EvaluateTropicalMC DROPPED Tree B's SeedBase / argv[5] runtime seed override, so cross-check #22 is unsatisfiable (SeedBase->99 emits the same 'seed = 42ULL + kp' as ->42; grep argv[5] tropical_eval.wl = 0). plan §2.4 requires the unified driver to carry Tree B's SeedBase argv, and §8.6 mandates the argv[5] runtime override for #22.

 4.0 RESTORE SEEDING (the fix). (i) Add the SeedBase option to EvaluateTropicalMC (read it via OptionValue, forward it to GenerateCppMonteCarlo and into vegasOpts). (ii) In the emitted C++ main, restore Tree B's runtime override (see OLD_CODE/TROPICAL_MONTE_CARLO2/TROPICAL_MONTE_CARLOv2/tropical_eval.wl ~L1784): 'uint64_t seed_base = (argc>5) ? strtoull(argv[5],nullptr,10) : <SeedBaseDefault>ULL;' (default = the SeedBase option value, 42 if unset), and use seed_base+kp where the seed is currently hardcoded 42ULL+kp. One compile, many seeds. (iii) REGENERATE the 9 codegen goldens under TEST/baselines/codegen_goldens/ to this seed-enabled baseline, and PROVE the change is purely additive: diff each new golden vs its Phase-3-committed version and show the ONLY differing lines are the documented seed-override (argv[5] parse + seed_base) lines — all other bytes identical. Commit the regenerated goldens. (iv) Add cross-check #22 coverage to TEST/test_complex_phase4.wl: same SeedBase → byte-identical output; two distinct argv[5] seeds on ONE compiled binary → distinct streams, each within 5σ of the reference. (v) Fix the prior summary's mislabel: #21b is C++ MC (NSamples), not VEGAS.

 4.1 Confirm the DEFAULT complex path exp(B·log P) with complex A and B works through the unified codegen for Tree-B Examples 2,15,16,19 and exact-Gamma Example 21 (no SplitRealImag yet).
 4.2 Implement complex-COEFFICIENT lifting: anchor z0=|C|^(1/k), residual c=C/z0^k (Tree B's residual mechanism extended); verify on a complex-coefficient toy.
 4.3 Implement the opt-in "SplitRealImag" VEGAS mode with BOTH §6.3 fixes: (i) complex log P (BUG-1) and (ii) MonoFactorLog (L3). Add the MonoFactorLog SectorData key in BOTH ProcessSector AND ProcessSectorLifted. The phase block in codegen is emitted ONLY when Im(B)≠0.

SELF-GATE (selfGatePass true ONLY if all hold):
 - SEEDING RESTORED: SeedBase->99 produces a DIFFERENT emitted seed / output than ->42; one compiled binary run with two distinct argv[5] seeds gives distinct streams (each within 5σ of reference); same seed → byte-identical. grep argv[5] tropical_eval.wl > 0. #22 coverage present and passing in TEST/test_complex_phase4.wl.
 - complex exponents match exact Gamma/Beta (#2) and NIntegrate (#1) to old tolerances.
 - the complex-BASE cross-check #21b PASSES: a polynomial with complex coefficients (arg P≠0) and complex B, lifted, where SplitRealImag ≡ Direct ≡ reference. This is the case that exposed BUG-1; without it the fix is unverified.
 - real-exponent codegen equals the regenerated Phase-4 goldens, and those goldens differ from the Phase-3-committed goldens ONLY in the documented seed-override lines (all other bytes identical) — diff and show it. Complex phase block emitted only when Im(B)≠0.
Long runs: run_in_background + poll, halt loudly on timeout; isolate INTERFILES per runner. KEEP WORK COMPACT (a prior run risked the StructuredOutput cap): summary ≤ ~1200 chars; commands ≤ 8 with keyOutput ≤ 200 chars; numerics ≤ 14 gate-critical entries; write the full transcript (all output, the golden-diff showing seed-only changes, the #22 streams) to TEST/phase4_report.txt and reference it in artifacts.`,
  verifierLenses: [
    {
      id: 'A-exact-symbolic', effort: 'max',
      instruction: 'Independently verify the complex-split identity symbolically. Confirm prod P^(Re B) · exp(i·Σ Im(B)·log P) reproduces prod P^B symbolically at sampled points, and that the MonoFactorLog linear form is correct. Pure exact; zero stats.',
      gate: 'plan.md Gate 4 + §8.5 Layer 1; ORCHESTRATION.md §3 Phase 4 lens A.',
      exactChecks: 'Simplify[ ∏P^(Re B)·Exp[I·Σ Im(B)·Log P] - ∏P^B ] === 0 at several complex sample points; the MonoFactorLog key value equals the correct linear combination of monomial logs; real-exponent codegen equals the regenerated Phase-4 goldens.',
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
      instruction: 'Verify the RESTORED seeding and reproducibility. (1) Through the actual package API: SeedBase->99 must emit different C++ / output than ->42 (the prior bug was SeedBase being silently ignored). (2) On ONE compiled #21b binary, two distinct argv[5] seeds must give DISTINCT output streams, each within 5σ of the reference; the SAME seed must give byte-identical output. grep argv[5] in tropical_eval.wl must be > 0. (3) Confirm the regenerated codegen goldens differ from the Phase-3-committed goldens (git show <prev commit>:TEST/baselines/codegen_goldens/*) ONLY in the seed-override lines — refute if ANY non-seed byte changed (a smuggled codegen change). (4) Budget stability of the #21b MC result (it is C++ MC / NSamples, NOT VEGAS). Refute on: SeedBase ignored, identical streams for distinct seeds, any non-seed golden byte change, or budget-dependent drift beyond error.',
      gate: 'plan.md Gate 4 + §8.5 Layer 2 (reseed #22) + §8.6 argv[5]; ORCHESTRATION.md §3 Phase 4 lens C.',
      exactChecks: 'same seed → byte-identical output (exact); distinct argv[5] seeds on one binary → distinct streams; regenerated goldens vs Phase-3-committed goldens differ ONLY in seed-override lines.',
      oracles: 'reseeded runs on one compiled binary; the converged #21b reference; budget sweep of the #21b MC value. Refute if any check fails.'
    }
  ],
  commitMsg: 'Phase 4: complex exponents + complex lifting + SplitRealImag; restored SeedBase/argv[5] runtime seed override (#22) + regenerated seed-enabled codegen goldens (unanimous 3-lens green gate)',
})
return r
