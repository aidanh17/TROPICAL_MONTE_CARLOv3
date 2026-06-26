export const meta = {
  name: 'tmcv3-phase6',
  description: 'Build TROPICAL_MONTE_CARLOv3 Phase 6 (re-run): fix engine bugs (#34/#43-IBP/#31-RunAllTests) + resolve #32, build deliverables (CI runner, EXAMPLES, MANUAL), re-fan corrected §8.2+§8.3 checks, coverage-matrix gate',
  phases: [ { title: 'Fix', model: 'opus' }, { title: 'Deliverables', model: 'sonnet' }, { title: 'Implement', model: 'sonnet' }, { title: 'Run+Verify', model: 'opus' }, { title: 'Synthesize', model: 'opus' }, { title: 'Commit' } ],
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
  { id:'9',   tier:2, oracles:'MC vs VEGAS parity |MC−VEGAS|≤5·√(errMC²+errVEGAS²)', exact:'n/a', fix:'CHECK BUG: the V4 reference list evaluated to a non-numeric symbolic (Im[s^-1]); supply numeric kinematics so the reference is numeric. Engine parity is sound.' },
  { id:'13',  tier:1, oracles:'numerical ε-slope estimate (relerr<5%)', exact:'n/a' },
  { id:'14',  tier:1, oracles:'G0=B⁰−Σ coeff·I pole identity', exact:'IBP↔subtraction pole identity' },
  { id:'15',  tier:1, oracles:'CC/diagnose_eps_eps correction (<0.01%)', exact:'finite=TOTAL−G0+γ·G0 identity', fix:'CHECK BUG: fatal ToExpression::sntx parse error (stray * in a line-154 comment) → zero assertions execute. Fix the syntax so the identity is actually tested.' },
  { id:'16',  tier:1, oracles:'subtracted-remainder → 0 as ε→0 at documented rate', exact:'n/a' },
  { id:'17',  tier:1, oracles:'exact lifted-vs-unlifted (relerr<0.5%; DroppedSectors asserted)', exact:'lift round-trip #40' },
  { id:'18',  tier:2, oracles:'lifted-vs-unlifted VarRed via EXACT trueSigma — use the HONEST gate: finite-variance case (all HasConstantTerm=True, e.g. Case D ≈31×) passes; infinite-variance case (Case A/C, HasConstantTerm=False) reported as L5 limitation + warning fires, NOT an N× pass', exact:'trueSigma σ²=I2−I1²', fix:'CHECK BUG: trueSigmaOne Set::shape destructuring bug made all trueSigma=INFINITE. Fix the destructuring; apply the Phase-3 honest variance gate (exact trueSigma, never sampled).' },
  { id:'19',  tier:1, oracles:'exact NIntegrate trueSigma; flags HasConstantTerm=False', exact:'σ²=I2−I1²', fix:'CHECK BUG: NIntegrate iterator-binding bug leaves I1 unevaluated + lifted fan non-simplicial. Fix the iterator binding (Block/With on the integration var); engine is sound.' },
  { id:'20',  tier:1, oracles:'anchor k-sweep vs k* rule (argmin relErr≥k* in H1/H2/H3)', exact:'kStar rule' },
  { id:'21',  tier:1, oracles:'complex-flattening contour rotation (identical value, ~20× variance drop)', exact:'sector value contour-invariance' },
  { id:'21b', tier:2, oracles:'SplitRealImag ≡ Direct ≡ reference (the BUG-1 case)', exact:'complex split identity §6.3' },
  { id:'22',  tier:1, oracles:'SeedBase determinism (same seed→byte-identical; distinct→distinct, each 5σ)', exact:'byte-identical stream' },
  { id:'23',  tier:2, oracles:'batched VEGAS == per-kp VEGAS within combined error; both match closed form', exact:'n/a (gates §5.1 fold)' },
  { id:'25',  tier:1, oracles:'Phase-1 codegen goldens (byte-identical emitted C++)', exact:'diff-clean for all kinds×samplers' },
  { id:'26',  tier:1, oracles:'GL(1)/Cheng–Wu affine↔simplex numeric equality', exact:'chart-change identity (validates #8)' },
  { id:'27',  tier:1, rev:'3', oracles:'graceful $Failed: nested/nestedIBP, nocuba, badcpp, liftnopivot, liftcomplex, liftdegenerate, vegasbudget', exact:'each guard fires cleanly', fix:'REAL ENGINE BUG (now fixed in the Fix stage, §5.8): cppBadTokens (tropical_eval.wl ~3237) matched BRACKET-form tokens (e.g. "Sin[") but CForm/MmaToC emits PAREN-form ("Sin("), so the badcpp guard never fired and uncompilable C++ was silently emitted. The Fix stage changes cppBadTokens to match the actually-emitted (paren/head) form so the guard fires on an unsupported head. ALSO (check-side): use a genuinely un-emittable head that reaches codegen (NOT Zeta[2], which auto-evaluates to a valid π²/6), and fix the nocuba sub-check TrueQ[detectCuba[]]["Found"] → TrueQ[detectCuba[]["Found"]]. After the engine fix, the badcpp guard must return $Failed on a real bad head on this CUBA-present machine.' },
  { id:'28',  tier:1, oracles:'EmptyDomain drop count, FeasibleFraction, HasConstantTerm tallies as asserted', exact:'feasibility bookkeeping', fix:'CHECK BUG: wrong expected Case-A HasConstantTerm value, and the FeasibleFraction 1%-MC gate rejects genuinely-feasible 1e-4-measure sectors. Correct the expected tallies against the actual engine output and relax/replace the MC feasibility gate with the exact feasibility test.' },
  { id:'29',  tier:1, oracles:'flattened magnitude warns when max>1e3 or min<1e-6', exact:'n/a' },
  { id:'30',  tier:1, oracles:'ValidateDecomposition/Subtraction/IBP/Lifted pass on known inputs', exact:'validator assertions' },
  { id:'31',  tier:1, new:true, oracles:'every Tree-A golden (Tests 1–22, divergent crosscheck) reproduced within MC noise (bitwise for symbolic)', exact:'merge fidelity vs Tree A', fix:'MIXED: (1,2) CHECK BUGS — reads key IBPTerms (engine emits Terms/NTerms via IBPReduceSector) and calls ProcessDivergentSector with 4 args (def takes 2); fix to the correct API. (D) REAL GAP — v3 has NO RunAllTests[] body; the Fix stage ports RunAllTests into v3, then this check calls it.' },
  { id:'32',  tier:1, new:true, oracles:'every Tree-B golden (Tests 1,2,3v2,5,6,7; Test 23; V/L VEGAS) reproduced', exact:'merge fidelity vs Tree B', fix:'T3v2-A: the Fix stage investigated whether v3 returning a numeric (vs Tree-B $Failed) is correct. Per that resolution, this check either accepts the numeric (cross-checked vs an independent oracle, with a documented capability-gain note) or asserts the restored $Failed refusal. Use the Fix stage outcome.' },
  { id:'33',  tier:2, new:true, rev:'3', oracles:'complex-B ε-integral through {IBP,subtraction}×{MC,VEGAS}: pairwise (pole,finite)≲1e-2 vs exact Γ-product ≲5e-3', exact:'Γ-product Laurent', fix:'CHECK BUG: cc_33.wl:97 printed finRef uses EulerGamma/2 instead of EulerGamma in the exact Γ-Laurent FINITE coefficient. Engine + all A1–A6 numerics are correct; fix the reference expression to EulerGamma so the mandated exact-symbolic finite coefficient matches.' },
  { id:'34',  tier:2, new:true, oracles:'≈8D lifting×complex-coeff×VEGAS: Direct≡SplitRealImag≡reference ≲0.5%; K-scaled fan builds; vegasbudget fires if under-budgeted', exact:'complex split identity', fix:'REAL ENGINE BUG (fixed in the Fix stage): computeFanScaled returned $Failed because the lifted complex monomial k-power coordinate resolved to Missing[KeyAbsent,k] instead of the integer exponent. After the engine fix, all 5 sub-checks must produce numbers and the K-scaled fan must build; corroborate vs an independent CUBA oracle.' },
  { id:'35',  tier:1, new:true, oracles:'Laurent via 3 routes (IBP@0, subtraction@0, finite-ε LaurentFromSubtraction) agree ≲1e-2', exact:'three-route pole/finite identity' },
  { id:'36',  tier:2, new:true, oracles:'per-kp MC, per-kp VEGAS, batched VEGAS, standalone one-point all agree; one-point==in-scan', exact:'n/a (batch fold + kinematics)' },
  { id:'37',  tier:1, new:true, oracles:'primitivized vs raw rays, K∈{1,n+2,…} → byte-identical C++ + bitwise-equal sector sum', exact:'fan scale-invariance' },
  { id:'38',  tier:1, new:true, oracles:'permuting variables / moving ε-regulator leaves answer invariant to MC noise', exact:'relabeling invariance' },
  { id:'39',  tier:2, new:true, oracles:'small-coeff case: unlifted plain-MC confidently wrong; lifting/VEGAS/reference recover truth AND warning fires', exact:'HasConstantTerm/magnitude/vegasbudget guard', fix:'CHECK BUGS: sub-check A reads a List as an Association; sub-check E asserts ALL sectors are variance-safe but 1 of 8 legitimately has HasConstantTerm=False (that warning firing IS the correct behavior per the Phase-3 L5 finding). Fix the data access and flip sub-E to EXPECT the warning to fire.' },
  { id:'40',  tier:1, new:true, oracles:'z→z0 ≡ original exactly for +, −, complex, 2-rule incommensurate lifts (rel-tol for floats)', exact:'lift identity (symbolic)' },
  { id:'41',  tier:1, new:true, oracles:'numerator factor+ε-pole and numerator+lifting vs NIntegrate/exact ≲1%', exact:'n/a' },
  { id:'42',  tier:2, new:true, oracles:'8D Automatic sizing reaches ~0.1% vs an oracle MATCHING the actual integrand; n=2 unchanged; vegasbudget fires at NSamples<20·NStart (note: code uses 100·NStart for d≥5)', exact:'n/a', fix:'CHECK BUGS: (a) the hardcoded oracle 0.00317086 is for a DIFFERENT polynomial than the check integrand (which includes x1^2; true value ≈0.00176159) — recompute the oracle for the ACTUAL integrand via independent Schwinger/NIntegrate; (b) the 8D fan must be built via the Automatic computeFanScaled path (raw ComputeDecomposition crashes polymake at 8D) — route through the public driver so K-scaling applies; (c) fix the vegasbudget interception.' },
  { id:'43',  tier:1, new:true, oracles:'symbolic decomposition output FreeQ[_,_Real] before MmaToC on BOTH convergent AND divergent/IBP paths; stray N[…] in decision path flagged', exact:'exactness invariant §1.4/§6.7', fix:'PARTLY REAL: convergent half genuinely passes. The IBP/divergent half is VACUOUS — the spec produces 0 divergent sectors and checkIBPResult keys on LaurentCoefficients which IBPProcessSector never returns. Fix the Fix stage to (1) use a spec that ACTUALLY produces divergent sectors, (2) key the exactness check on the field IBPProcessSector really returns, and (3) genuinely assert FreeQ[_,_Real] on the IBP decomposition output.' },
]

// ============================ STAGE 1 — FIX (opus): real engine bugs + #32 investigation ============================
phase('Fix')
const fix = await agent(
  `TROPICAL_MONTE_CARLOv3 Phase 6 engine-fix stage. Root: /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3, branch v3-build. NEVER edit OLD_CODE/; tropical_fan.wl only additive; re-load (Get) after every edit; keep the §1.4 exactness invariant.

RESUME SCOPE — IMPORTANT: items (A) #34 lifted-complex fan, (B) #43 IBP field/spec, (C) #31 RunAllTests, and (D) #32 T3v2-A are ALREADY FIXED AND VERIFIED in the working-tree tropical_eval.wl (uncommitted) from the prior run — do NOT redo them and do NOT re-run their long verifications. Your ONLY new code task this run is (E) below. Then do a QUICK regression check: package loads; the 9 codegen goldens are byte-identical; phase2/3/5 self-gate scripts still pass; cc_31/cc_32/cc_34/cc_43 still exit 0.

 (E) REAL ENGINE BUG #27 (§5.8 — the badcpp guard is INEFFECTIVE): cppBadTokens (tropical_eval.wl ~line 3237) matches BRACKET-form tokens like "Sin[", but CForm / MmaToC emits PAREN-form like "Sin(" — so the badcpp scan never matches and uncompilable C++ (unsupported heads) is SILENTLY emitted instead of returning $Failed. Fix cppBadTokens to detect the actually-emitted form (match the paren-form "name(" for the unsupported heads, or scan the post-CForm string for the bad head names as emitted), so the badcpp guard fires and GenerateCppMonteCarlo returns $Failed when an un-emittable head (e.g. a genuinely unsupported special function) reaches codegen. Verify: a genuinely un-emittable head now triggers $Failed (the badcpp guard fires); and a NORMAL integrand still compiles unchanged (the 9 codegen goldens remain byte-identical — the guard must not alter good output).

For REFERENCE ONLY (already done — do NOT redo) — the prior run fixed all of the following; they are context, not new tasks:

 (A) REAL BUG #34 — computeFanScaled returns $Failed on lifted COMPLEX 8D vertices: the lifted complex monomial's k-power (aux/z) coordinate resolves to Missing[KeyAbsent, k] instead of its integer exponent. Find where the lifted Newton-polytope vertices are built for complex-coefficient lifting and ensure the k exponent is correctly populated (no Missing[]). Verify: an 8D lifting×complex-coeff×VEGAS case builds a K-scaled fan AND Direct ≡ SplitRealImag ≡ an INDEPENDENT CUBA oracle to ≲0.5%.

 (B) REAL GAP #43/IBP-exactness — the exactness invariant is currently UNTESTED on the divergent/IBP path. (1) Determine the field IBPProcessSector actually returns its Laurent/pole data under (the check wrongly keyed on 'LaurentCoefficients', which IBPProcessSector never returns). (2) Ensure a divergent spec that ACTUALLY produces divergent sectors exists, and that IBPProcessSector's symbolic output is genuinely FreeQ[_,_Real] before the MmaToC boundary. Report the correct field name so cc_43 can be fixed.

 (C) REAL GAP #31D — v3 has NO RunAllTests[] body (0 DownValues). Port RunAllTests[] (Tree-A Tests 1–18 in-package suite) into v3 tropical_eval.wl so merge-fidelity #31 can call it; it MUST reproduce the Tree-A goldens (TEST/baselines/treeA.txt). Mirror OLD_CODE/TROPICAL_MONTE_CARLO/EXAMPLES test logic; do not weaken any test.

 (D) #32 T3v2-A INVESTIGATION (user: investigate, then decide). Tree B returned $Failed for negative test T3v2-A; v3's unified IBP path returns a NUMERIC result. Determine: what is the T3v2-A spec, WHY did Tree B refuse (read OLD_CODE/TROPICAL_MONTE_CARLO2/.../EXAMPLES), and is v3's numeric result actually CORRECT (cross-check vs independent NIntegrate/exact)? THEN ACT: if v3's number is correct → it is a legitimate capability gain: leave the behavior, and document it (one line in SUMMARY/notes) so cc_32 can be updated to ACCEPT the numeric (cross-checked vs the oracle); if v3's number is WRONG → it is a regression: restore the refusal so v3 returns $Failed like Tree B. Report your decision + the evidence (the spec, the oracle value, v3's value, your verdict).

SELF-GATE selfGatePass=true ONLY if: (E) the badcpp guard now FIRES (GenerateCppMonteCarlo returns $Failed) on a genuinely un-emittable head, while a normal integrand still compiles and the 9 codegen goldens remain BYTE-IDENTICAL; AND the regression check holds: package loads, phase2/3/5 self-gate scripts pass, cc_31/cc_32/cc_34/cc_43 still exit 0 (the already-done A/B/C/D from the prior run remain intact). Do not re-run the long A/B/C/D verifications.
OUTPUT FORMAT: call StructuredOutput with a JSON OBJECT, top-level keys exactly summary (short prose, ≤1500 chars; INCLUDE your #32 verdict + the correct IBP Laurent field name), selfGatePass (boolean), commands (≤8, keyOutput ≤200 chars), numerics (≤14), filesChanged (array), artifacts (array). Do NOT emit XML tags; do NOT nest the result inside summary. Spool full detail to TEST/phase6_fix_report.txt.`,
  { label:'fix-engine', phase:'Fix', model:'opus', effort:'max', schema:WORK })
log(`FIX STAGE: selfGatePass=${fix && fix.selfGatePass}; #32 + IBP-field decision in summary: ${fix ? fix.summary.slice(0,400) : 'FIX AGENT DIED'}`)

// ============================ STAGE 2 — DELIVERABLES (parallel): CI runner, EXAMPLES, MANUAL ============================
phase('Deliverables')
const deliverables = (await parallel([
  () => agent(`Build the single CI entry point TEST/run_all.wl for TROPICAL_MONTE_CARLOv3 (plan.md §6.1; root /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3, branch v3-build, NEVER edit OLD_CODE/). It must: run RunAllTests[] (Tree-A in-package suite) + every TEST/cc_*.wl, each in its own INTERFILES/cc_<id>/ working dir (§10 clobber hazard), collect PASS/FAIL/SKIP, print a summary, and EXIT NON-ZERO if any check FAILs. Tier-2 (CUBA) / Tier-3 (FIESTA) checks must SKIP cleanly (exit 0 contribution) when their dependency is absent — never silently wrong. Return the file path + a smoke-run result.`,
    { label:'deliv:ci-runner', phase:'Deliverables', model:'sonnet', effort:'medium' }),
  () => agent(`Build the §8.4 worked examples under EXAMPLES/ for TROPICAL_MONTE_CARLOv3 (root /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3, branch v3-build, NEVER edit OLD_CODE/). Each is a runnable .wl that Gets ./tropical_eval.wl and doubles as the named cross-checks: Ex-V3a (divergent + complex, full Laurent via IBP, cross-checked vs subtraction + exact Γ — #8/#33/#35); Ex-V3b (small complex coefficient, lifted + VEGAS, cross-check vs CUBA Cuhre — #34/#21b); Ex-V3c (kinematic scan, batched VEGAS, per-point vs NIntegrate + per-kp VEGAS — #23/#36); Ex-V3d (numerator (…)^{+3/2} × divergent denominator — #41). Each prints PASS/FAIL and exits non-zero on FAIL. Return the file paths + run results.`,
    { label:'deliv:examples', phase:'Deliverables', model:'sonnet', effort:'medium' }),
  () => agent(`Build the docs under MANUAL/ for TROPICAL_MONTE_CARLOv3 (plan.md §6.3; root /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3, branch v3-build, NEVER edit OLD_CODE/). Rewrite SUMMARY.txt for v3 (consolidated API, the convergent + divergence + lifting + complex feature set, the documented limitations incl. the L5 infinite-variance/HasConstantTerm note and the SplitRealImag Re(B) caveat) and merge the two OLD MANUAL/manual.tex files into MANUAL/manual.tex (convergent + divergence + lifting + complex; reference the complex_flatten_spiral figure and the anchor-selection section). Confirm manual.tex compiles (pdflatex or at least latex parse) OR report cleanly if no TeX toolchain; SUMMARY.txt must exist and be coherent. Return file paths + build result.`,
    { label:'deliv:manual', phase:'Deliverables', model:'sonnet', effort:'medium' }),
])).filter(Boolean)
log(`DELIVERABLES: ${deliverables.length}/3 stages returned`)

// ============================ STAGE 3 — RE-FAN corrected checks (impl sonnet/med → verify opus/high) ============================
phase('Implement')
const results = await pipeline(CHECKS,
  c => agent(
    `Implement/REPAIR v3 cross-check #${c.id} (plan.md §8) as TEST/cc_${c.id}.wl in /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3 (branch v3-build; NEVER edit OLD_CODE/). It must Get the v3 ./tropical_eval.wl (and ./tropical_fan.wl as needed via absolute path), exercise the feature per plan.md §8.2/§8.3, and print "CC${c.id} PASS …" / "CC${c.id} FAIL expected=<e> got=<g>". Tier ${c.tier}: ${c.tier===1?'WL+g++ only.':c.tier===2?'needs CUBA — if CUBA absent, print "CC'+c.id+' SKIP (no CUBA)" and exit 0.':'needs FIESTA — if FIESTA absent, print "CC'+c.id+' SKIP (no FIESTA)" and exit 0.'} Mirror the OLD test logic where one exists; call the engine with the CORRECT public API (verify key names/arity against tropical_eval.wl, do not invent keys).${c.fix ? ' KNOWN ISSUE TO FIX (a prior run refuted this check): ' + c.fix + ' The engine-fix stage has already landed the relevant engine changes; make the check correct and non-vacuous.' : ' This check passed previously — keep it correct; re-confirm it still passes after the engine fixes.'} Return the file path.`,
    { label:`impl:${c.id}`, phase:'Implement', model:'sonnet', effort:c.fix ? 'medium' : 'low' }),
  (_, c) => agent(
    `${c.rev?`(revision ${c.rev}: this check was REPAIRED after a prior refutation — re-verify FRESHLY from scratch, do not reuse any earlier verdict.) `:''}Run TEST/cc_${c.id}.wl in its OWN isolated working dir INTERFILES/cc_${c.id}/ (cd /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3; §10 fixed-filename clobber hazard — never share INTERFILES across parallel runners). Then, as an INDEPENDENT SKEPTIC who did NOT write the check, corroborate every numeric result against ${c.oracles}; default to refuted=true and refute if no independent oracle agrees within combined error. Exact checks that MUST hold symbolically (zero floating point): ${c.exact}. If the check legitimately SKIPs (dependency absent), set skipped=true, refuted=false. Long runs: run_in_background + poll, halt loudly on timeout.
OUTPUT FORMAT: return the result by calling StructuredOutput with a JSON OBJECT whose TOP-LEVEL keys include lens (string), refuted (boolean), evidence (string), checkId="${c.id}", tier="${c.tier}", skipped (boolean). Do NOT emit XML/HTML tags like <refuted>; do NOT nest the result inside evidence — refuted/skipped/checkId are SEPARATE top-level JSON keys.`,
    { label:`verify:${c.id}`, phase:'Run+Verify', model:'opus', effort:'high', schema:VERDICT })
)
const verdicts = results.filter(Boolean)

// ---- Synthesize the §8.1 coverage-matrix gate (the final release gate) ----
phase('Synthesize')
const coverage = await agent(
  `You are assembling the FINAL release gate for TROPICAL_MONTE_CARLOv3 (plan.md §8.1 feature-coverage matrix; ORCHESTRATION.md §3 Phase 6 Gate). Here are the independent verifier verdicts for every cross-check: ${JSON.stringify(verdicts)}.
The engine-fix stage reported: ${fix ? JSON.stringify({summary:fix.summary, selfGatePass:fix.selfGatePass}) : 'FIX STAGE DIED'}.
Build the §8.1 feature-coverage report: for EACH feature row (F1, F2a, F2b, F2c, F3a, F3b, F3c, F4a, F4b, F4c, F5a, F5b, F5c, MC back-end, the-merge-itself, numerator B>0, exactness-invariant, external-validation) list the covering checks and whether ≥1 of them PASSED (refuted=false AND not skipped).
Independently CONFIRM the deliverables EXIST and are non-empty by checking the filesystem yourself: TEST/run_all.wl (single CI entry point, exits non-zero on any FAIL); EXAMPLES/ Ex-V3a/b/c/d (§8.4, doubling as #8/#21b/#23/#33/#34/#35/#36/#41); MANUAL/SUMMARY.txt + MANUAL/manual.tex (docs build). An empty deliverable is a release blocker.
SELF-GATE selfGatePass=true ONLY if: the fix stage selfGatePass was true; NO feature row has zero passing checks (an empty row is a release blocker); NONE of the v3-original combination checks #31–#43 are refuted; Tier-2 checks that ran were corroborated and Tier-2/3 that lacked CUBA/FIESTA skipped CLEANLY (not silently wrong); TEST/run_all.wl + EXAMPLES/Ex-V3a..d + MANUAL/ all exist and are non-empty (docs build). Populate emptyFeatureRows[] (include any missing deliverable) and refutedChecks[] with any offenders. List every SKIP explicitly (no silent truncation, §5.4 safety rail).
OUTPUT FORMAT: return by calling StructuredOutput with a JSON OBJECT whose TOP-LEVEL keys are exactly summary (short prose string), selfGatePass (boolean), commands (array), numerics (array), filesChanged (array), artifacts (array), emptyFeatureRows (array), refutedChecks (array). Do NOT emit XML/HTML tags like <selfGatePass>; do NOT nest the whole result inside summary — these are SEPARATE top-level JSON keys.`,
  { label:'coverage', phase:'Synthesize', model:'opus', effort:'high', schema:WORK })

const fixOk     = !!(fix && fix.selfGatePass)              // fail-closed: dead/failed fix stage → RED
const refutedNew = verdicts.filter(v => v.refuted && !v.skipped && CHECKS.find(c => c.id===v.checkId && c.new))
const emptyRows  = (coverage && coverage.emptyFeatureRows) || []
const green = fixOk && !!(coverage && coverage.selfGatePass) && refutedNew.length === 0 && emptyRows.length === 0
log(`PHASE 6 GATE: ${green ? 'GREEN' : 'RED'} (fixOk=${fixOk}, emptyRows=${emptyRows.length}, refuted #31-43=${refutedNew.length}, refuted total=${verdicts.filter(v=>v.refuted&&!v.skipped).length})`)

if (green) {
  phase('Commit')
  await agent(`The final release gate is GREEN. Commit the v3-build worktree. cd /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3, branch v3-build: git add -A (OLD_CODE/ + INTERFILES/ gitignored) then git commit -m "Phase 6: engine fixes (#34/#43-IBP/#31-RunAllTests) + #32 resolved; full cross-check suite + EXAMPLES + MANUAL + CI runner; §8.1 coverage matrix complete (final green gate)". Return the commit hash. Mechanical only.`,
    { label:'commit', phase:'Commit', model:'sonnet', effort:'low' })
}
return { green, fix, coverage, verdicts }
