export const meta = {
  name: 'tmcv3-phase6-patch',
  description: 'Phase 6 close-out: fix the 3 high-D/variance CHECK-HARNESS bugs (#18,#34,#42) + give 8D VEGAS adequate budget; re-verify only those; re-synthesize coverage to honest green. Engine is FINAL (commit 2961475) — do NOT edit tropical_eval.wl/tropical_fan.wl.',
  phases: [ { title: 'FixChecks', model: 'opus' }, { title: 'Verify', model: 'opus' }, { title: 'Synthesize', model: 'opus' }, { title: 'Commit' } ],
}

const VERDICT = { type:'object', required:['lens','refuted','evidence'],
  properties:{ lens:{type:'string'}, refuted:{type:'boolean'}, confidence:{type:'string'}, evidence:{type:'string'},
    checkId:{type:'string'}, skipped:{type:'boolean'},
    independentChecks:{type:'array',items:{type:'object',properties:{check:{type:'string'},result:{type:'string'}}}},
    disagreements:{type:'array',items:{type:'string'}}, notes:{type:'string'} } }
const WORK = { type:'object', required:['summary','selfGatePass','filesChanged'],
  properties:{ summary:{type:'string'}, selfGatePass:{type:'boolean'}, filesChanged:{type:'array',items:{type:'string'}},
    deliverablesOk:{type:'boolean'}, notes:{type:'string'} } }

const ROOT = '/Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3'

// The 3 remaining CHECK-HARNESS defects (engine numerics independently confirmed correct in the prior run).
const TARGETS = [
  { id:'18', desc:'lifted-vs-unlifted variance reduction (Tier-2)',
    bug:'weak NIntegrate (PrecisionGoal=4/MaxRecursion=20 + a TimeConstrained that maps a spurious timeout to $Failed→Infinity) made Case D trueSigma read INFINITE, though Case D variance is genuinely FINITE.',
    fix:'Compute trueSigma with a CONVERGED NIntegrate (high WorkingPrecision/MaxRecursion; no TimeConstrained→Infinity coercion). Apply the Phase-3 HONEST gate: a finite-variance case (all HasConstantTerm=True, e.g. Case D) reports a real exact VarRed (≈31×); an infinite-variance case (Case A/C, HasConstantTerm=False) is reported as the documented L5 limitation with the warning firing — never a sampled-σ N× pass.',
    oracle:'exact symbolic trueSigma σ²=I2−I1² (NIntegrate at high precision) and the Phase-3 result Case D VarRed≈31.6×' },
  { id:'34', desc:'lifting×complex-coeff×VEGAS 8D (Tier-2, v3-original #31–43)',
    bug:'sub-A (the BUG-34 lifted-complex fan) PASSES (10 sectors, no Missing[]); sub-B/C aborted ENVIRONMENTALLY because the 8D VEGAS ran at NSamples=4e6, UNDER the 8.1e6 vegasbudget threshold, so it never converged ($resB unassigned, rel-err 999) — and the check then mis-extracted the empty/under-budget result.',
    fix:'Fix the result extraction so an under-budget/aborted run is handled cleanly. Run the 8D lifted-complex VEGAS at an ADEQUATE budget (NSamples ≥ ~4e7, comfortably above 20·NStart so it converges) with run_in_background + poll and a GENEROUS timeout (it is slow). Demonstrate Direct ≡ SplitRealImag ≡ an independent CUBA/Schwinger oracle to ≲0.5%. Separately, a small 4e6 run may be used ONLY to demonstrate vegasbudget fires (expected). Do NOT change the engine.',
    oracle:'independent CUBA Cuhre / Schwinger reduction for the 8D integrand (the Phase-5 truth ≈0.0031709 for that polynomial family); Direct≡Split exact match' },
  { id:'42', desc:'resolveVegasSizing high-D (Tier-2, v3-original #31–43)',
    bug:'(a) the check used a hardcoded oracle for the WRONG polynomial; (b) it built the 8D fan via raw ComputeDecomposition (crashes polymake) instead of the Automatic computeFanScaled path; (c) its vegasbudget detector (a HandlerBlock) could not catch the fired message.',
    fix:'Compute the oracle for the ACTUAL check integrand via an independent Schwinger/NIntegrate reduction. Build the 8D fan through the Automatic computeFanScaled path (public driver) so K-scaling applies. Detect the vegasbudget warning robustly (Quiet+`MessageList`/`Check`, or capture General::message — NOT a broken HandlerBlock). Run 8D at adequate budget (≥~4e7) with run_in_background+poll+generous timeout for the convergence sub-check, and a 4e6 run to confirm vegasbudget fires. Do NOT change the engine.',
    oracle:'independent Schwinger/NIntegrate for the ACTUAL integrand; vegasbudget message fires at NSamples<threshold, silent above; n=2 unchanged' },
]

// ---- Stage 1: fix the 3 check scripts (engine FINAL, do not touch) ----
phase('FixChecks')
const fixWork = await agent(
  `Phase-6 close-out — repair THREE cross-check SCRIPTS only (the v3 ENGINE is FINAL at commit 2961475: do NOT edit tropical_eval.wl or tropical_fan.wl; if you change them you have failed). Root: ${ROOT}, branch v3-build. NEVER edit OLD_CODE/. The prior run independently confirmed the engine numerics are CORRECT in every case below — these are CHECK-HARNESS / budget / detection bugs.

Fix each of these TEST/cc_<id>.wl scripts and confirm it then prints "CC<id> PASS" and exits 0:
${TARGETS.map(t => ` • #${t.id} (${t.desc})\n     BUG: ${t.bug}\n     FIX: ${t.fix}`).join('\n')}

8D VEGAS is genuinely slow: launch those runs with run_in_background and poll to completion with a generous timeout; halt loudly (do not silently pass) if a run truly times out. Use an isolated INTERFILES/cc_<id>/ working dir per check (§10 clobber hazard). Keep z0/lift data exact; numericize only at the integration boundary.

SELF-GATE selfGatePass=true ONLY if: cc_18, cc_34, cc_42 all now print "CC<id> PASS" and exit 0 with their numerics corroborated against the stated independent oracle; tropical_eval.wl and tropical_fan.wl are UNCHANGED from commit 2961475 (verify with git diff --stat — must be empty for those two files); and the 9 codegen goldens remain byte-identical. Also set deliverablesOk=true iff SUMMARY.txt exists at the REPO ROOT (plan §3.2 — NOT under MANUAL/), and MANUAL/manual.tex+manual.pdf, EXAMPLES/Ex-V3a..d, TEST/run_all.wl all exist non-empty.
OUTPUT FORMAT: call StructuredOutput with a JSON OBJECT, top-level keys summary (≤1200 chars), selfGatePass (boolean), filesChanged (array), deliverablesOk (boolean), notes. Do NOT emit XML tags; do NOT nest fields inside summary. Spool detail to TEST/phase6_patch_report.txt.`,
  { label:'fix-checks', phase:'FixChecks', model:'opus', effort:'high', schema:WORK })
log(`FIXCHECKS: selfGatePass=${fixWork && fixWork.selfGatePass}; deliverablesOk=${fixWork && fixWork.deliverablesOk}`)

// ---- Stage 2: independent adversarial re-verification of ONLY the 3 fixed checks ----
phase('Verify')
const verdicts = (await parallel(TARGETS.map(t => () =>
  agent(`You are an INDEPENDENT VERIFIER (you did NOT write or fix this check). Re-run TEST/cc_${t.id}.wl in its OWN isolated dir ${ROOT}/INTERFILES/cc_${t.id}_v/ (cd ${ROOT}; §10 clobber hazard). Default refuted=true; refute unless YOUR OWN checks prove it correct. The ENGINE is final (commit 2961475) — confirm tropical_eval.wl/tropical_fan.wl are unchanged (git diff --stat empty for them); if the "fix" edited the engine, REFUTE. Independently corroborate the check's numeric result against: ${t.oracle}. Refute if the engine numerics are wrong, if the check still fails/aborts, if it passes VACUOUSLY, or if it relies on sampled-σ where exact trueSigma is mandated. 8D VEGAS is slow — run_in_background + poll + generous timeout.
OUTPUT FORMAT: StructuredOutput JSON object, top-level keys lens (string), refuted (boolean), evidence (string), checkId="${t.id}", skipped (boolean). No XML tags; do not nest fields in evidence.`,
    { label:`verify:${t.id}`, phase:'Verify', model:'opus', effort:'high', schema:VERDICT }))
)).filter(Boolean)

// ---- Stage 3: re-synthesize the §8.1 coverage gate (honest green) ----
phase('Synthesize')
const cov = await agent(
  `Final §8.1 coverage-gate re-synthesis for TROPICAL_MONTE_CARLOv3. CONTEXT: in the prior full Phase-6 run, ALL feature rows were covered by ≥1 passing check EXCEPT F5c (high-D VEGAS sizing), whose only checks #34 and #42 were refuted; #18 was also refuted (non-blocking — its rows F3a are covered by passing #17/#20/#28/#39/#40); #8 skipped cleanly (FIESTA absent); all deliverables exist (SUMMARY.txt at REPO ROOT per plan §3.2, MANUAL/manual.tex+pdf, EXAMPLES/Ex-V3a..d, TEST/run_all.wl). The three outstanding checks were just repaired and INDEPENDENTLY re-verified — here are the fresh verdicts: ${JSON.stringify(verdicts)}. The fix-stage deliverablesOk=${fixWork ? fixWork.deliverablesOk : 'unknown'}.
Independently CONFIRM on the filesystem: SUMMARY.txt at ${ROOT}/SUMMARY.txt (root, NOT MANUAL/), MANUAL/manual.tex + manual.pdf, EXAMPLES/Ex-V3a..d, TEST/run_all.wl all exist non-empty.
Decide the FINAL release gate: GREEN iff (a) #34 and #42 are now refuted=false (so F5c is covered and no #31–#43 combination check is refuted), (b) #18 is refuted=false (or its rows remain covered by passing checks — state which), (c) all deliverables exist with SUMMARY.txt at ROOT. Set selfGatePass accordingly; list any remaining emptyFeatureRows[] / refutedChecks[].
OUTPUT FORMAT: StructuredOutput JSON object, top-level keys summary (string), selfGatePass (boolean), emptyFeatureRows (array), refutedChecks (array). No XML tags; do not nest fields in summary.`,
  { label:'coverage', phase:'Synthesize', model:'opus', effort:'high',
    schema:{ type:'object', required:['summary','selfGatePass'], properties:{ summary:{type:'string'}, selfGatePass:{type:'boolean'}, emptyFeatureRows:{type:'array',items:{type:'string'}}, refutedChecks:{type:'array',items:{type:'string'}} } } })

const refutedTargets = verdicts.filter(v => v.refuted && !v.skipped).length
const engineUntouched = !!(fixWork && fixWork.selfGatePass)   // fix-stage asserts engine unchanged + checks pass
const emptyRows = (cov && cov.emptyFeatureRows) || []
const green = engineUntouched && verdicts.length === TARGETS.length && refutedTargets === 0
            && !!(cov && cov.selfGatePass) && emptyRows.length === 0
log(`PATCH GATE: ${green ? 'GREEN' : 'RED'} (engineUntouched=${engineUntouched}, refutedTargets=${refutedTargets}/${TARGETS.length}, covSelfGate=${cov && cov.selfGatePass}, emptyRows=${emptyRows.length})`)

if (green) {
  phase('Commit')
  await agent(`The final release gate is GREEN. Commit the v3-build worktree. cd ${ROOT}, branch v3-build: git add -A (OLD_CODE/ + INTERFILES/ gitignored) then git commit -m "Phase 6 FINAL (green): repaired high-D/variance check harness (#18/#34/#42); §8.1 coverage matrix complete — all feature rows covered, #31-#43 confirmed, deliverables present. v3 build complete." Return the commit hash. Mechanical only; do NOT edit any source.`,
    { label:'commit', phase:'Commit', model:'sonnet', effort:'low' })
}
return { green, fixWork, verdicts, coverage: cov }
