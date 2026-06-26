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
  producePrompt: `CRITICAL OUTPUT FORMAT (a prior run FAILED 5× here): you MUST return the result by calling the StructuredOutput tool with a JSON OBJECT whose TOP-LEVEL keys are exactly: summary (a SHORT prose string), selfGatePass (boolean true/false), commands (array of objects), numerics (array of objects), filesChanged (array of strings), artifacts (array of strings). Do NOT emit XML/HTML tags such as <selfGatePass> or <summary>. Do NOT put the whole result inside the summary string — summary is a brief prose field ONLY; selfGatePass/commands/numerics/artifacts are SEPARATE top-level JSON keys. Keep it compact (summary ≤1200 chars; ≤8 commands; ≤14 numerics) and spool full detail to TEST/phase5_report.txt.

TROPICAL_MONTE_CARLOv3 Phase 5 — high-D hardening + anchor finalization (plan.md §7 Phase 5 + §6.4 + §6.5; ORCHESTRATION.md §3 Phase 5). Root: /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3, branch v3-build. NEVER edit OLD_CODE/. tropical_fan.wl may ONLY receive the additive computeFanScaled — no other change (§10). The sizing fix is specified in the bug log lift_error_log.md (L1/L2). Re-load after every edit.

RE-RUN NOTE (a prior Phase-5 attempt was RED). Your computeFanScaled (tropical_fan.wl, additive) and resolveVegasSizing + vegasbudget (tropical_eval.wl) are ALREADY PRESENT in the working tree (uncommitted on top of Phase-4 commit cb8da7d) and the sizing MACHINERY is verified sound (monotone scaling, low-D byte-identical, n=2 unchanged, and a plain convergent 8D case converges to 1e-6). PRESERVE that work; build the two fixes below on top. Two issues were found:

 ISSUE 1 (THE REAL ONE — fix it). The 8D lifted+complex+VEGAS test value is SYSTEMATICALLY WRONG: Re=0.00316713 vs independent truth 0.0031709 (CUBA Cuhre @1e9 and a Schwinger semi-analytic reduction agree to 7.8e-6) — 0.118% LOW, biased below ALL independent oracles, ~3437σ outside Cuhre's error bar, and vegasbudget did NOT fire (NSamples 4e6 > 20·81000) so it passed silently. The prior gate measured 0.126% against its OWN NIntegrate MC (same-method self-corroboration — PROHIBITED by §8.5 Layer 2). You MUST: (a) diagnose the root cause of the systematic low bias (candidates: an under-sampled heavy-tail / HasConstantTerm=False sector at high-D; SplitRealImag variance built on Re(B); NStart still too small for this integrand) — report which; (b) make the 8D DEFAULT-sizing value agree with an INDEPENDENT oracle (CUBA Cuhre on the direct integrand — NEVER your own NIntegrate MC) to ≈0.1%; (c) strengthen the under-resolution guard so this case CANNOT pass silently — vegasbudget (or an added convergence/cross-check guard) MUST fire whenever the VEGAS error bar understates the true error at high-D (the prior 20·NStart threshold was too loose: 4e6 samples passed yet the answer was 0.118% wrong with a tight bar). If the regime is genuinely unresolvable by SplitRealImag, the guard must say so and refuse to report a confident value.

 ISSUE 2 (#37 mis-specified — make the decomposition EXACT, per the user). The literal old #37 ("K=1 vs K=n+2 rays → byte-identical C++ + bitwise-equal sector sum") is FALSE and was the wrong thing to test: K-scaling the vertices legitimately changes the triangulation (translateToOriginInteger picks a scale-dependent lattice point) → different but equally-valid decompositions. What MUST hold (plan §1.4 exactness invariant) is that the sector decomposition the code produces is EXACT: for the fan actually used, the symbolic sum of sector contributions (with Jacobians/prefactors, under the change of variables) reproduces the ORIGINAL integrand as an EXACT symbolic identity, with the decomposition output FreeQ[_,_Real] (no floating-point decisions). Add/keep a check proving this EXACT sector-decomposition identity. The only byte-level guarantee that matters is the no-regression one: the wired computeFanScaled path (which tries K=1 first and returns it on success) is byte-identical to the prior raw ComputeDecomposition path — already true (#25, 9 goldens byte-identical). Do NOT chase byte-identity across different K.

SELF-GATE (selfGatePass true ONLY if all hold):
 - 8D lifted+complex+VEGAS, DEFAULT sizing, agrees with an INDEPENDENT oracle (CUBA Cuhre on the direct integrand) to ≈0.1% — corroborated independently, NOT against your own NIntegrate MC.
 - The under-resolution guard FIRES on the previously-silent 8D case (and on NSamples<20·NStart); it never reports a confident value when the VEGAS error bar understates the true high-D error. Demonstrate it firing on the exact case that previously passed silently.
 - low-D (n=2) VEGAS UNCHANGED (relerr ~1e-6); all 9 Phase-1 codegen goldens byte-identical.
 - EXACT sector decomposition (#37 redefined): the symbolic sector sum reproduces the original integrand exactly (FreeQ[_,_Real]); the wired computeFanScaled path is byte-identical to the raw ComputeDecomposition path.
Long runs: the 9D robust fan (~260 s) and 8D VEGAS/CUBA exceed 10 min — run_in_background + poll, halt loudly on timeout; VEGAS serial; isolate INTERFILES per runner. KEEP WORK COMPACT: summary ≤ ~1200 chars; commands ≤ 8 (keyOutput ≤ 200 chars); numerics ≤ 14; write the full transcript (root-cause diagnosis, the CUBA-corroborated 8D value, the guard firing, the exact-decomposition identity) to TEST/phase5_report.txt and reference it in artifacts.`,
  verifierLenses: [
    {
      id: 'exact-sector-decomposition', effort: 'high',
      instruction: 'Verify the sector decomposition is EXACT (#37 redefined, per §1.4). Do NOT test byte-identity across different K — that is legitimately false (K-scaling changes the triangulation). Instead: (1) For the fan the code actually uses, independently confirm the symbolic sum of sector contributions (with Jacobians/prefactors, under the change of variables) reproduces the ORIGINAL integrand as an EXACT symbolic identity, and the decomposition output is FreeQ[_,_Real] (no floating-point decisions in the decomposition). (2) Confirm the no-regression guarantee: the wired computeFanScaled path (K=1 returned on success) is byte-identical to the prior raw ComputeDecomposition path, and computeFanScaled edits to tropical_fan.wl are purely additive (diff vs Phase-0 tropical_fan.wl touches only the new helper). Refute if the decomposition is not exact, carries a stray _Real in a decision path, or the wired path diverges byte-wise from raw.',
      gate: 'plan.md §1.4 exactness invariant + Gate 5 #37 (redefined) + §8.5 Layer 1.',
      exactChecks: 'symbolic Σ(sector integrands) ≡ original integrand exactly (PossibleZeroQ / Expand) for the fan used; decomposition output FreeQ[_,_Real]; md5(wired computeFanScaled codegen) === md5(raw ComputeDecomposition codegen); tropical_fan.wl diff vs Phase-0 is additive-only.',
      oracles: 'n/a — pure exact symbolic. Refute if the exact identity fails, a decision numeric survives, or the wired path is not byte-identical to raw.'
    },
    {
      id: 'independent-integrator-8D', effort: 'max',
      instruction: 'THE L1 TRAP LENS — a prior run shipped an 8D value 0.118% low (3437σ) with a tight lying error bar; you previously caught it. Re-derive the 8D truth INDEPENDENTLY with ≥2 oracles sharing no code with the producer (build your OWN CUBA Cuhre at large budget; corroborate with CUBA Vegas and/or a Schwinger semi-analytic reduction and the c=0 closed form). Then: (1) the producer DEFAULT-sizing 8D value must agree with your independent truth to ≈0.1% and must NOT be systematically biased below all oracles. Reject any value corroborated only against the producer own NIntegrate MC (prohibited same-method self-corroboration). (2) Confirm the strengthened guard ACTUALLY FIRES on the previously-silent case (NSamples 4e6 that passed before): a tight VEGAS error bar that understates the true high-D error must trigger a warning / refusal. Refute if the 8D value still lacks an independent corroborator within combined error, is systematically biased, or the guard fails to fire when the error bar lies.',
      gate: 'plan.md Gate 5 (numeric, L1 silent-wrong trap) + §8.5 Layer 2 (≥2 independent oracles).',
      exactChecks: 'n/a (numeric lens) — but confirm the under-resolution guard fires on the previously-silent 8D case (a discrete, checkable event).',
      oracles: 'your OWN CUBA Cuhre at large budget; CUBA Vegas; Schwinger semi-analytic; c=0 closed form. Refute if the 8D value disagrees with the independent truth beyond combined error, is biased low vs all oracles, or relies on same-method self-corroboration.'
    },
    {
      id: 'reproducibility-sizing', effort: 'high',
      instruction: 'Verify sizing-sweep monotonicity and low-D non-regression: as the sizing budget increases the 8D estimate must converge monotonically toward the reference (no budget-dependent bias), and n=2 VEGAS must be unchanged (~1e-6). Refute on non-monotone drift or a low-D regression.',
      gate: 'plan.md Gate 5 (stability + low-D) + §8.5 Layer 2.',
      exactChecks: 'n/a (numeric lens).',
      oracles: 'a sizing sweep you run yourself; the Phase-1/Phase-4 n=2 VEGAS golden. Refute on non-monotonicity or low-D change.'
    }
  ],
  commitMsg: 'Phase 5: computeFanScaled K-scaling + resolveVegasSizing/vegasbudget; fixed 8D lifted+complex VEGAS bias (CUBA-corroborated ~0.1%) + strengthened under-resolution guard; #37 redefined to exact sector-decomposition identity (unanimous 3-lens green gate)',
})
return r
