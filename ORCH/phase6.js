export const meta = {
  name: 'tmcv3-phase6',
  description: 'Build TROPICAL_MONTE_CARLOv3 Phase 6: full §8.2+§8.3 cross-check suite + §8.4 examples, fanned out (impl → run+independently-verify → coverage-matrix synthesis gate)',
  phases: [ { title: 'Implement', model: 'sonnet' }, { title: 'Run+Verify', model: 'opus' }, { title: 'Synthesize', model: 'opus' }, { title: 'Commit' } ],
}

// ============================ shared schemas (ORCHESTRATION.md §2) ============================
const WORK = { type:'object', required:['summary','selfGatePass','commands','numerics','artifacts'],
  properties:{ summary:{type:'string'}, filesChanged:{type:'array',items:{type:'string'}},
    commands:{type:'array',items:{type:'object',properties:{cmd:{type:'string'},exitCode:{type:'number'},keyOutput:{type:'string'}}}},
    numerics:{type:'array',items:{type:'object',properties:{name:{type:'string'},value:{type:'string'},claimedRef:{type:'string'},deltaSigma:{type:'string'}}}},
    selfGatePass:{type:'boolean'}, artifacts:{type:'array',items:{type:'string'}},
    emptyFeatureRows:{type:'array',items:{type:'string'}}, refutedChecks:{type:'array',items:{type:'string'}} } }
const VERDICT = { type:'object', required:['lens','refuted','evidence'],
  properties:{ lens:{type:'string'}, refuted:{type:'boolean'}, confidence:{type:'string'}, evidence:{type:'string'},
    checkId:{type:'string'}, tier:{type:'string'}, skipped:{type:'boolean'},
    independentChecks:{type:'array',items:{type:'object',properties:{check:{type:'string'},result:{type:'string'}}}},
    disagreements:{type:'array',items:{type:'string'}}, notes:{type:'string'} } }
// ========================================================================================

// The cross-check catalog (plan.md §8.2 #1-30 + §8.3 #31-43). tier: 1=WL+g++, 2=needs CUBA, 3=needs FIESTA.
// `new` marks the v3-original combination checks (#31-43) whose refutation hard-blocks release.
const CHECKS = [
  { id:'1',   tier:1, oracles:'NIntegrate on the direct integrand', exact:'n/a' },
  { id:'2',   tier:1, oracles:'exact Gamma/Beta/Dirichlet closed form (Example 21 ~1e-6)', exact:'closed-form Gamma identity' },
  { id:'3',   tier:1, oracles:'exact Laurent via Gamma products (B=4→(1/6,-1/6), B=5→(1/12,-1/8))', exact:'Gamma-product Laurent' },
  { id:'4',   tier:1, oracles:'IBP vs subtraction agree on (pole,finite)', exact:'IBP↔subtraction pole identity #14' },
  { id:'5',   tier:2, oracles:'CUBA Cuhre on the direct integrand (1-5%; MC within 5σ)', exact:'n/a' },
  { id:'6',   tier:2, oracles:'CUBA VEGAS (independent sampler), combined error', exact:'n/a' },
  { id:'7',   tier:2, oracles:'CUBA on raw vs decomposed+flattened (report ratio)', exact:'n/a' },
  { id:'8',   tier:3, oracles:'FIESTA5 via GL(1)/Cheng–Wu (skip cleanly if FIESTA absent)', exact:'GL(1) chart identity #26' },
  { id:'9',   tier:2, oracles:'MC vs VEGAS parity |MC−VEGAS|≤5·√(errMC²+errVEGAS²)', exact:'n/a' },
  { id:'13',  tier:1, oracles:'numerical ε-slope estimate (relerr<5%)', exact:'n/a' },
  { id:'14',  tier:1, oracles:'G0=B⁰−Σ coeff·I pole identity', exact:'IBP↔subtraction pole identity' },
  { id:'15',  tier:1, oracles:'CC/diagnose_eps_eps correction (<0.01%)', exact:'finite=TOTAL−G0+γ·G0 identity' },
  { id:'16',  tier:1, oracles:'subtracted-remainder → 0 as ε→0 at documented rate', exact:'n/a' },
  { id:'17',  tier:1, oracles:'exact lifted-vs-unlifted (relerr<0.5%; DroppedSectors asserted)', exact:'lift round-trip #40' },
  { id:'18',  tier:2, oracles:'lifted-vs-unlifted VarRed via EXACT trueSigma (Case A≈18×; Case C not lifted)', exact:'trueSigma σ²=I2−I1²' },
  { id:'19',  tier:1, oracles:'exact NIntegrate trueSigma; flags HasConstantTerm=False', exact:'σ²=I2−I1²' },
  { id:'20',  tier:1, oracles:'anchor k-sweep vs k* rule (argmin relErr≥k* in H1/H2/H3)', exact:'kStar rule' },
  { id:'21',  tier:1, oracles:'complex-flattening contour rotation (identical value, ~20× variance drop)', exact:'sector value contour-invariance' },
  { id:'21b', tier:2, oracles:'SplitRealImag ≡ Direct ≡ reference (the BUG-1 case)', exact:'complex split identity §6.3' },
  { id:'22',  tier:1, oracles:'SeedBase determinism (same seed→byte-identical; distinct→distinct, each 5σ)', exact:'byte-identical stream' },
  { id:'23',  tier:2, oracles:'batched VEGAS == per-kp VEGAS within combined error; both match closed form', exact:'n/a (gates §5.1 fold)' },
  { id:'25',  tier:1, oracles:'Phase-1 codegen goldens (byte-identical emitted C++)', exact:'diff-clean for all kinds×samplers' },
  { id:'26',  tier:1, oracles:'GL(1)/Cheng–Wu affine↔simplex numeric equality', exact:'chart-change identity (validates #8)' },
  { id:'27',  tier:1, oracles:'graceful $Failed: nested/nestedIBP, nocuba, badcpp, liftnopivot, liftcomplex, liftdegenerate, vegasbudget', exact:'each guard fires cleanly' },
  { id:'28',  tier:1, oracles:'EmptyDomain drop count, FeasibleFraction, HasConstantTerm tallies as asserted', exact:'feasibility bookkeeping' },
  { id:'29',  tier:1, oracles:'flattened magnitude warns when max>1e3 or min<1e-6', exact:'n/a' },
  { id:'30',  tier:1, oracles:'ValidateDecomposition/Subtraction/IBP/Lifted pass on known inputs', exact:'validator assertions' },
  { id:'31',  tier:1, new:true, oracles:'every Tree-A golden (Tests 1–22, divergent crosscheck) reproduced within MC noise (bitwise for symbolic)', exact:'merge fidelity vs Tree A' },
  { id:'32',  tier:1, new:true, oracles:'every Tree-B golden (Tests 1,2,3v2,5,6,7; Test 23; V/L VEGAS) reproduced', exact:'merge fidelity vs Tree B' },
  { id:'33',  tier:2, new:true, oracles:'complex-B ε-integral through {IBP,subtraction}×{MC,VEGAS}: pairwise (pole,finite)≲1e-2 vs exact Γ-product ≲5e-3', exact:'Γ-product Laurent' },
  { id:'34',  tier:2, new:true, oracles:'≈8D lifting×complex-coeff×VEGAS: Direct≡SplitRealImag≡reference ≲0.5%; K-scaled fan builds; vegasbudget fires if under-budgeted', exact:'complex split identity' },
  { id:'35',  tier:1, new:true, oracles:'Laurent via 3 routes (IBP@0, subtraction@0, finite-ε LaurentFromSubtraction) agree ≲1e-2', exact:'three-route pole/finite identity' },
  { id:'36',  tier:2, new:true, oracles:'per-kp MC, per-kp VEGAS, batched VEGAS, standalone one-point all agree; one-point==in-scan', exact:'n/a (batch fold + kinematics)' },
  { id:'37',  tier:1, new:true, oracles:'primitivized vs raw rays, K∈{1,n+2,…} → byte-identical C++ + bitwise-equal sector sum', exact:'fan scale-invariance' },
  { id:'38',  tier:1, new:true, oracles:'permuting variables / moving ε-regulator leaves answer invariant to MC noise', exact:'relabeling invariance' },
  { id:'39',  tier:2, new:true, oracles:'small-coeff case: unlifted plain-MC confidently wrong; lifting/VEGAS/reference recover truth AND warning fires', exact:'HasConstantTerm/magnitude/vegasbudget guard' },
  { id:'40',  tier:1, new:true, oracles:'z→z0 ≡ original exactly for +, −, complex, 2-rule incommensurate lifts (rel-tol for floats)', exact:'lift identity (symbolic)' },
  { id:'41',  tier:1, new:true, oracles:'numerator factor+ε-pole and numerator+lifting vs NIntegrate/exact ≲1%', exact:'n/a' },
  { id:'42',  tier:2, new:true, oracles:'8D Automatic sizing reaches ~0.1%; n=2 unchanged ~1e-6; vegasbudget fires at NSamples<20·NStart', exact:'n/a' },
  { id:'43',  tier:1, new:true, oracles:'symbolic decomposition output FreeQ[_,_Real] before MmaToC; stray N[…] in decision path flagged', exact:'exactness invariant §1.4/§6.7' },
]

// ---- Implement (cheap tier) → Run+independently-Verify (opus/high) : pipelined, no barrier ----
phase('Implement')
const results = await pipeline(CHECKS,
  c => agent(
    `Implement v3 cross-check #${c.id} (plan.md §8) as TEST/cc_${c.id}.wl in /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3 (branch v3-build; NEVER edit OLD_CODE/). It must Get the v3 ./tropical_eval.wl (and ./tropical_fan.wl as needed via absolute path), exercise the feature per plan.md §8.2/§8.3, and print "CC${c.id} PASS …" / "CC${c.id} FAIL expected=<e> got=<g>". Tier ${c.tier}: ${c.tier===1?'WL+g++ only.':c.tier===2?'needs CUBA — if CUBA absent, print "CC'+c.id+' SKIP (no CUBA)" and exit 0.':'needs FIESTA — if FIESTA absent, print "CC'+c.id+' SKIP (no FIESTA)" and exit 0.'} Mechanical port from the OLD trees where one exists; mirror the OLD test's logic. Return the file path.`,
    { label:`impl:${c.id}`, phase:'Implement', model:'sonnet', effort:'low' }),
  (_, c) => agent(
    `Run TEST/cc_${c.id}.wl in its OWN isolated working dir INTERFILES/cc_${c.id}/ (cd /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3; §10 fixed-filename clobber hazard — never share INTERFILES across parallel runners). Then, as an INDEPENDENT SKEPTIC who did NOT write the check, corroborate every numeric result against ${c.oracles}; default to refuted=true and refute if no independent oracle agrees within combined error. Exact checks that MUST hold symbolically (zero floating point): ${c.exact}. If the check legitimately SKIPs (dependency absent), set skipped=true, refuted=false. Long runs: run_in_background + poll, halt loudly on timeout. Return VERDICT with checkId="${c.id}", tier="${c.tier}".`,
    { label:`verify:${c.id}`, phase:'Run+Verify', model:'opus', effort:'high', schema:VERDICT })
)
const verdicts = results.filter(Boolean)

// ---- Synthesize the §8.1 coverage-matrix gate (the final release gate) ----
phase('Synthesize')
const coverage = await agent(
  `You are assembling the FINAL release gate for TROPICAL_MONTE_CARLOv3 (plan.md §8.1 feature-coverage matrix; ORCHESTRATION.md §3 Phase 6 Gate). Here are the independent verifier verdicts for every cross-check: ${JSON.stringify(verdicts)}.
Build the §8.1 feature-coverage report: for EACH feature row (F1, F2a, F2b, F2c, F3a, F3b, F3c, F4a, F4b, F4c, F5a, F5b, F5c, MC back-end, the-merge-itself, numerator B>0, exactness-invariant, external-validation) list the covering checks and whether ≥1 of them PASSED (refuted=false AND not skipped).
Also confirm a single CI entry point (TEST/run_all.wl analogue) exists that exits non-zero on any FAIL, and that the §8.4 worked examples (Ex-V3a/b/c/d under EXAMPLES/) exist and double as #8/#21b/#23/#33/#34/#35/#36/#41.
SELF-GATE selfGatePass=true ONLY if: NO feature row has zero passing checks (an empty row is a release blocker); NONE of the v3-original combination checks #31–#43 are refuted; Tier-2 checks that ran were corroborated and Tier-2/3 that lacked CUBA/FIESTA skipped CLEANLY (not silently wrong); docs build. Populate emptyFeatureRows[] and refutedChecks[] with any offenders. List every SKIP explicitly (no silent truncation, §5.4 safety rail).`,
  { label:'coverage', phase:'Synthesize', model:'opus', effort:'high', schema:WORK })

const refutedNew = verdicts.filter(v => v.refuted && !v.skipped && CHECKS.find(c => c.id===v.checkId && c.new))
const emptyRows  = coverage.emptyFeatureRows || []
const green = !!coverage.selfGatePass && refutedNew.length === 0 && emptyRows.length === 0
log(`PHASE 6 GATE: ${green ? 'GREEN' : 'RED'} (emptyRows=${emptyRows.length}, refuted #31-43=${refutedNew.length}, refuted total=${verdicts.filter(v=>v.refuted&&!v.skipped).length})`)

if (green) {
  phase('Commit')
  await agent(`The final release gate is GREEN. Commit the v3-build worktree. cd /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3, branch v3-build: git add -A (OLD_CODE/ + INTERFILES/ gitignored) then git commit -m "Phase 6: full cross-check suite + examples; §8.1 coverage matrix complete (final green gate)". Return the commit hash. Mechanical only.`,
    { label:'commit', phase:'Commit', model:'sonnet', effort:'low' })
}
return { green, coverage, verdicts }
