# Orchestrator Plan — building TROPICAL_MONTE_CARLOv3 via per-phase unattended ultracode workflows

**Companion to `plan.md`.** `plan.md` says *what* v3 is and *what* each phase must achieve + its gate. This document says *how* to execute it: one **unattended ultracode workflow per phase**, each with an **independent adversarial verifier** (plan.md §8.5) and **tiered effort/model**, with **checkpoints between phases**. Nothing here is run yet — these are the workflow scripts to launch when you say go.

Design constraints carried from the conversation:
- **Unattended within a phase, checkpoint between phases** — a phase runs start→gate with nobody watching; it halts at its gate (green or red) and reports.
- **No numeric PASS is self-graded.** Every gate is decided by a *separate* verifier agent that did not produce the result (plan.md §8.5). Most of the decomposition is verified by **exact symbolic identities** (plan.md §1.4) that cannot be subtly wrong; only the integrated value is checked numerically, against ≥2 independent oracles.
- **Tiered effort.** Opus+max only on the hard surgery and the verifier; cheap tiers for mechanical/runtime-bound legs.
- **Never edit `OLD_CODE/`.** The build accumulates in a dedicated worktree (§1.4).

---

## §1 Execution model

### 1.1 One workflow per phase, the main loop drives the gaps
Each plan.md phase (0–6) is **one `Workflow` invocation**. The main Claude loop (the orchestrator) fires phase *N*, waits for its `task-notification`, reads the structured `PHASE_REPORT`, then:
- **checkpoint mode (default, recommended through Phase 5):** post the report to you and wait for `next`;
- **auto-chain-on-green mode:** immediately fire phase *N+1* iff the gate was green; **halt** on any red.

You pick the mode once (and may differ it per phase — e.g. auto-chain 0→1, checkpoint 2→3→4→5, auto-chain 6). The phases are a strict dependency chain, so they are never run in parallel with each other; parallelism lives *inside* a phase. **The orchestrator itself implements nothing and the gate is fail-closed — see §1.7.**

### 1.2 The producer → verifier → gate contract (every phase)
Inside a phase workflow:
1. **Producer** (Opus, max/xhigh) performs that phase's edits/work and runs its own gate commands, returning a structured `WORK_RESULT` (files changed, commands+exit codes, self-reported numbers, artifacts). The producer may *not* declare the phase green.
2. **Independent verifier(s)** (Opus, high–max; fresh context; **read-only — verify, do not fix**) receive the gate spec + the producer's artifacts and are instructed to **refute**. They *re-run* the gate commands themselves and run the exact-symbolic identity checks. For the integrated value they consult ≥2 independent oracles (exact Γ/Β, CUBA Cuhre, CUBA Vegas, IBP↔subtraction, SplitRealImag↔Direct, reseed, NIntegrate). Each returns a `VERDICT`.
3. **Gate decision** is **plain code, no model** and **fail-closed (§1.7)**: GREEN iff the producer's self-gate passed *and* **every** verifier returns `refuted=false` (independently confirmed). **Any** verifier that refutes — finds something wrong *or* cannot independently corroborate a numeric value — makes the phase RED and **halts the run**. No majority vote. Complex-lifting phases (4–5) run a 3-lens panel and require **all** lenses to confirm.

### 1.3 Effort / model tiering (set per `agent()` call)
| Leg | model | effort |
|---|---|---|
| Core engine surgery, complex-lifting, exactness audit (plan §6.1/6.3/6.7) | `opus` | `max`/`xhigh` |
| Independent adversarial verifier + judge panels (plan §8.5) | `opus` | `high`–`max` |
| Cross-check authoring, synthesis, report assembly | `opus`/`sonnet` | `medium`–`high` |
| Mechanical ports, deletions, running suites, scaffolding, golden capture | `sonnet`/`haiku` | `low`–`medium` |

Rationale (from the effort discussion): much of the build is **runtime-bound, not reasoning-bound**; and on verbatim-port/deletion legs, high effort risks *over-thinking* and breaking the byte-identical golden invariant — so cheap is both faster and safer there. The verifier is the one place never to economize.

### 1.4 Isolation, working directories, toolchain
- Run the **entire build in a dedicated git worktree/branch** of this repo (e.g. `v3-build`), so the main tree stays clean until you review the diff. The producer edits there; you diff between phases.
- **Do not** use `isolation:'worktree'` per agent for the *sequential* producer (the build must accumulate in one tree). Use it only where parallel agents would mutate the same files — which in this plan is rare (most parallel legs are read-only verification / cross-check *runs*).
- **Parallel cross-check runners that compile+run C++ must not share `INTERFILES/`** (plan.md §10 hazard: fixed filenames clobber). Give each parallel runner its own `WorkingDirectory` (e.g. `INTERFILES/cc_<id>/`).
- Toolchain confirmed present on this machine: `wolframscript`, `polymake`, `g++`, CUBA. FIESTA absent → Tier-3 checks skip.

### 1.5 Long runs, timeouts, resume
- A single `Bash` call caps at **10 min**; `RunAllTests` (tens of min), the 9D robust fan (~260 s), and large VEGAS runs can exceed it. Long WL/compile/run steps must be launched with **`run_in_background`** and polled, or split into per-test invocations. Each phase script sets explicit per-step time budgets and **halts the phase** (not silently continues) on timeout.
- Every `Workflow` invocation returns a `runId`. If a phase is interrupted or you edit its script, **resume** with `Workflow({scriptPath, resumeFromRunId})` — the unchanged prefix of `agent()` calls returns cached results; only new/edited calls re-run.
- Scripts avoid `Date.now()`/`Math.random()` (unavailable in workflow scripts); any timestamp/seed is passed via `args`.

### 1.6 Running the orchestrator itself
The orchestrator is the main Claude loop — its model is the session model, its effort the session `/effort`. It is **decoupled** from the sub-agents: producer/verifier model+effort are set per `agent()` inside the scripts (§1.3), independent of the orchestrator.
- **Opus is a fine default** for the orchestrator (it is your interface and relays verifier findings clearly), but it does **not** need `max` effort. Under §1.7 it does no implementation, and the gate decision is **deterministic code in the workflow script**, not a model judgment — so "halt if any sub-agent finds something wrong" is guaranteed regardless of the orchestrator's model. Run the orchestrator at **medium** effort.
- Cheaper does not weaken the build: quality lives in the sub-agents (producers Opus/max, verifiers Opus/high–max), which keep their tiers no matter the session setting. (Sonnet would also suffice for pure coordination; Opus@medium is the safe default.)
- **Do not run the whole session at `max`** — that only inflates the cheap coordination turns. Net: **orchestrator = Opus @ medium; hard surgery & verifier = Opus @ max/high inside the scripts.**

### 1.7 The orchestrator implements NOTHING; the gate is FAIL-CLOSED (your requirements)
Two hard rules, both user-confirmed:
1. **Zero implementation by the orchestrator.** The orchestrator only launches each phase workflow, reads its `PHASE_REPORT`, enforces the gate, and relays to you. It **never** edits a package file, writes v3 code, or "fixes" a failing check. *All* producing/editing happens inside workflow sub-agents (the producers). If a phase needs a fix, that is a new producer task in a re-run of that phase — never the orchestrator acting directly.
2. **Fail-closed gate — ANY sub-agent finding a problem stops the run.** A phase is GREEN only if the producer's self-gate passed **and every** verifier returns `refuted=false` (independently *confirmed* correct). If **any** verifier refutes — finds something wrong **or** cannot independently corroborate a numeric value — or fails to run at all, the phase is RED and the run **halts immediately**. No majority vote, no "most checks passed." Stricter than the earlier majority rule, and the default. Enforced by deterministic code (`green = selfGatePass && allRan && refuted===0`), so it does not rely on the orchestrator model's judgment.

---

## §2 The reusable phase harness

Every phase script is an instance of this skeleton. (Workflow scripts are plain JS; `meta` must be a pure literal.)

```js
export const meta = {
  name: 'tmcv3-phaseN',
  description: 'Build TROPICAL_MONTE_CARLOv3 phase N: <goal> with independent adversarial gate',
  phases: [ { title: 'Produce' }, { title: 'Verify' } ],
}

// ---- structured-output schemas (validated at the tool boundary) ----
const WORK = { type:'object', required:['summary','selfGatePass','commands','numerics','artifacts'],
  properties:{
    summary:{type:'string'},
    filesChanged:{type:'array',items:{type:'string'}},
    commands:{type:'array',items:{type:'object'}},      // {cmd, exitCode, keyOutput}
    numerics:{type:'array',items:{type:'object'}},       // {name, value, claimedRef, deltaSigma}
    selfGatePass:{type:'boolean'},
    artifacts:{type:'array',items:{type:'string'}} } }

const VERDICT = { type:'object', required:['lens','refuted','evidence'],
  properties:{
    lens:{type:'string'},
    refuted:{type:'boolean'},                            // true = I found it wrong / unproven
    confidence:{type:'string'},
    independentChecks:{type:'array',items:{type:'object'}}, // {cmd|identity, result}
    disagreements:{type:'array',items:{type:'string'}},
    notes:{type:'string'} } }

// ---- the harness ----
async function runPhase({ producePrompt, verifierLenses, panel = 1 }) {
  phase('Produce')
  const work = await agent(producePrompt, {
    label: 'produce', phase: 'Produce', model: 'opus', effort: 'max', schema: WORK })

  phase('Verify')
  // independent verifiers: fresh context, READ-ONLY, instructed to refute; they re-run the gate themselves
  const verdicts = (await parallel(
    verifierLenses.map(lens => () =>
      agent(verifyPrompt(lens, work), {
        label: `verify:${lens.id}`, phase: 'Verify',
        model: 'opus', effort: lens.effort || 'high', schema: VERDICT }))
  )).filter(Boolean)

  const refuted = verdicts.filter(v => v.refuted).length
  const allRan  = verdicts.length === verifierLenses.length              // a dead/absent verifier counts as NOT confirmed
  const green   = work.selfGatePass && allRan && refuted === 0           // FAIL-CLOSED: any refute (or missing verifier) → RED
  log(`PHASE GATE: ${green ? 'GREEN' : 'RED'} (selfGate=${work.selfGatePass}, refuted=${refuted}/${verifierLenses.length}, allRan=${allRan})`)
  return { green, work, verdicts }
}

function verifyPrompt(lens, work) { return `
You are an INDEPENDENT VERIFIER. You did NOT produce this result. Your job is to REFUTE it; default to refuted=true unless your own independent checks prove it correct. READ-ONLY: do not edit package files.
Lens: ${lens.id} — ${lens.instruction}
Gate spec for this phase: ${lens.gate}
Producer's claimed artifacts/results: ${JSON.stringify(work)}
Procedure:
 1. Re-run the gate commands yourself (wolframscript / g++ / diff); do not trust the producer's exit codes.
 2. Run the EXACT symbolic identity checks (these cannot be subtly wrong): ${lens.exactChecks}
 3. For any NUMERIC value, compare against the independent oracle(s): ${lens.oracles}. Assume the value is wrong; hunt for a reference that disagrees beyond combined error.
 4. Return VERDICT: refuted=true with disagreements if ANY check fails or any numeric lacks an independent corroborator.` }
```

The phase report returned to the orchestrator is `{green, work, verdicts}`; the main loop renders it and checkpoints or auto-chains (§4).

---

## §3 Per-phase workflow specs

Each phase supplies a `producePrompt`, the gate checks, and the verifier lenses. Mechanical legs inside a producer prompt should themselves be delegated to cheap-tier sub-work where the script fans out; the tables below note effort.

### Phase 0 — Setup & golden baselines  *(mostly mechanical + runtime-bound)*
- **Produce (sonnet/medium):** create the `v3-build` worktree + v3 skeleton (`tropical_fan.wl` from the identical OLD copy); run **both** old suites to capture goldens (Tree A Tests 1–22 + divergent crosscheck; Tree B suite + Test 23 + V/L VEGAS) **in parallel, each in its own WorkingDirectory**, backgrounded (they exceed 10 min); fan smoke test.
- **Gate:** both baseline transcripts saved with the documented reference numbers (e.g. Tree A B=4→(1/6,−1/6); Tree B Test 6A=7.850075e-4); fan smoke passes.
- **Verifier lenses (1, opus/high):** re-run the fan smoke + spot-check 3 reference numbers against the published values in the OLD SUMMARYs; refute if any transcript is missing or a number is off.

```js
// phase0.js
export const meta = { name:'tmcv3-phase0', description:'v3 Phase 0: worktree, skeleton, capture both old goldens, fan smoke', phases:[{title:'Produce'},{title:'Verify'}] }
/* ...WORK/VERDICT/runPhase from §2... */
const r = await runPhase({
  producePrompt: `Set up the v3 build. In a NEW git worktree 'v3-build' of this repo:
   1. Create the v3 skeleton per plan.md §3.2; copy tropical_fan.wl from OLD_CODE/TROPICAL_MONTE_CARLO (identical everywhere).
   2. Capture goldens IN PARALLEL, each in its own WorkingDirectory, run_in_background (they exceed 10 min): Tree A RunAllTests+IBP+test_divergent_crosscheck → TEST/baselines/treeA.txt ; Tree B run_validation_suite+test_lifted+test_vegas → TEST/baselines/treeB.txt.
   3. Fan smoke: PolytopeVertices+ComputeDecomposition on (1+x1^2+x2^2)^-1.
   NEVER edit OLD_CODE/. Return WORK with the captured reference numbers.`,
  verifierLenses: [{ id:'baseline', effort:'high',
    instruction:'confirm both golden transcripts exist and key reference numbers match the OLD SUMMARYs',
    gate:'plan.md Phase 0 / Gate 0', exactChecks:'fan ray/sector counts for the smoke poly',
    oracles:'OLD_CODE SUMMARY.txt published values (Tree A B=4/B=5 Laurent; Tree B Test 6A)' }],
})
return r
```

### Phase 1 — Seed from Tree A; capture v3 goldens  *(mechanical)*
- **Produce (sonnet/medium):** copy Tree A `tropical_eval.wl` → v3; confirm it loads and reproduces Tree A goldens (MC + per-kp **and batched** VEGAS); capture v3 codegen goldens for every kind (conv/g0/g1/rem/ibp; MC, per-kp VEGAS, batched VEGAS).
- **Gate:** loads; Tests 1–22 + crosscheck reproduce Phase-0 Tree-A numbers; all goldens captured.
- **Verifier (1, opus/high):** re-run Tests 1–22; diff a sampled golden; refute on any mismatch.

### Phase 2 — Refactor core; unify codegen & drivers; fold batching  *(HARD surgery)*
- **Produce (opus/max):** Module-0 primitives single-source; reconcile `ProcessSector ⊕ FlattenSector` with **exact** ε threading (plan §6.1) — preserve the `IsDivergent` pre-flatten branch; collapse the convergent/IBP codegen towers into one parametrized emitter, **folding `emitVegasBatchMain` into the VEGAS `Batch` mode** (plan §5.1/§5.2); unify drivers (§5.4); de-dup validation+log-exp (§5.5); single `detectCuba` (§5.3); begin the **exactness audit** (§6.7).
- **Gate (the "simplify aggressively, prove nothing broke" gate):** golden-master codegen regression GREEN for all kinds under MC / per-kp VEGAS / **batched** VEGAS; Tests 1–22 + crosscheck reproduce Phase-0 numbers **exactly**; `RunBenchmark2D`/`findCubaPrefix`/`emitVegasBatchMain`/IBP-twin symbols gone (grep); exactness guard #43 passes on the touched paths.
- **Verifier lenses (2, opus/high):** (a) **golden diff + grep** — re-emit all codegen kinds and `diff` against Phase-1 goldens; re-run Tests 1–22; (b) **exactness** — assert `FreeQ[symbolic decomposition output, _Real]` before `MmaToC` and that no `10^-…` tolerance remains in a decision path. Either refutes → RED.

### Phase 3 — Graft lifting (real coefficients), sandbox-first  *(HARD + sandbox gate)*
- **Produce (opus/max):** port Tree B `SANDBOX/` toys + `sandbox_lift_common.wl`; run read-only vs the v3 core (proves the §6.1 seam); port `DetectExtremeCoefficients`/`LiftCoefficients`/`ProcessSectorLifted`/`ValidateLiftedDecomposition`; apply BUG-2 (aux-var robustness) + relative-tolerance `liftidentity` + `SuggestedK→kStar` + geometry-gate anchor (§6.2); wire `Lift` into the driver; domain-indicator codegen.
- **Gate:** Test 23 (A exactness <0.5% + DroppedSectors count; B end-to-end C++ <1%; C error paths fire; D EmptyDomain) reproduces Tree-B numbers; `bench_lift_variance` reproduces Case A VarRed ≈ 18× and Case C correctly **not** lifted; `HasConstantTerm` reported per sector.
- **Verifier lenses (2):** (a) **exact** — symbolic lift round-trip `z→z0 ≡ original` (#40) for +, −, multi-rule; (b) **independent integrator** — re-derive Test-23A/B values via NIntegrate **and** CUBA Cuhre, refute if the lifted sum disagrees. Variance claims are corroborated by exact `trueSigma` (NIntegrate σ²=I2−I1²), never by sample σ alone.

### Phase 4 — Complex exponents end-to-end + complex lifting  *(HARDEST — the L1/L3/BUG-1 zone; 3-lens panel)*
- **Produce (opus/max):** confirm the default `exp(B·log P)` complex path (Examples 2,15,16,19,21); complex-coefficient lifting (anchor `|C|^{1/k}`, residual `C/z0^k`); the opt-in `SplitRealImag` VEGAS mode with **both** fixes — complex `log P` (BUG-1) + `MonoFactorLog` (L3) — adding the `MonoFactorLog` SectorData key in `ProcessSector` & `ProcessSectorLifted` (plan §6.3).
- **Gate:** complex vs exact Γ/Β (#2) and NIntegrate (#1) to old tolerances; **the complex-BASE check #21b passes** (a polynomial with complex coefficients, `arg P≠0`, complex B, lifted, where `SplitRealImag ≡ Direct ≡ reference`); real-exponent codegen byte-identical (phase block only when `Im(B)≠0`).
- **Verifier panel (3, opus/max):** lens A **exact symbolic** (the split identity `∏P^{Re B}·exp(iΣ Im(B)·logP)` reproduces `∏P^B` symbolically at sample points; `MonoFactorLog` linear form correct); lens B **independent integrator** (CUBA Cuhre on the direct integrand with `arg P≠0`); lens C **reproducibility** (reseed #22 + budget stability). Accept only on **unanimous** non-refute — any lens that refutes or cannot corroborate halts the run. *This is the phase the user cannot eyeball — the panel is the safeguard.*

### Phase 5 — High-D hardening + anchor finalization  *(mixed; verifier-sensitive)*
- **Produce (opus/high for fan+sizing, sonnet for wiring):** port `computeFanRobust`→`computeFanScaled` (K-scaling) into `tropical_fan.wl` + wire the automatic lifted-fan path (§6.4); implement `resolveVegasSizing` + `vegasbudget` (§6.5).
- **Gate:** an 8D lifted-complex-coeff + VEGAS case converges with **Automatic** sizing (~0.1% vs reference); `vegasbudget` fires when under-budgeted; low-D (n=2) VEGAS unchanged (~1e-6); K-scaled vs raw-ray fans give identical codegen+result (#37).
- **Verifier panel (3, opus/high):** exact (fan scale-invariance #37 symbolic); independent integrator (8D reference via CUBA Vegas at large budget); reproducibility (sizing sweep monotonicity). Refute if the "converged" 8D value lacks an independent corroborator (this is the L1 silent-wrong-answer trap).

### Phase 6 — Cross-check suite, benchmarks, docs  *(FAN-OUT — the big parallel win)*
- **Produce:** implement the full §8.2+§8.3 suite under `TEST/` + the §8.4 examples; a single CI entry point exiting non-zero on any FAIL.
- **Fan-out (pipeline, mixed tiers):** the ~43 cross-checks are largely independent → run them concurrently, **each in its own `INTERFILES/` working dir**, cheap tier for the mechanical run-and-collect, opus/high for the adversarial verification leg of each numeric result. Then a synthesis agent (opus/high) builds the coverage report.
- **Gate (final):** every §8.1 feature-coverage-matrix row maps to ≥1 passing check (an empty row blocks release); all new combination checks #31–#42 + #43 pass; CUBA checks degrade gracefully, FIESTA checks skip cleanly; docs build.

```js
// phase6.js (fan-out shape)
export const meta = { name:'tmcv3-phase6', description:'v3 Phase 6: full cross-check suite + examples + coverage gate', phases:[{title:'Implement'},{title:'Run+Verify'},{title:'Synthesize'}] }
const CHECKS = [ /* {id, tier, script, oracles, exact?} for §8.2 + §8.3 (#1..#43) */ ]
const results = await pipeline(CHECKS,
  // stage 1: implement/port the check (cheap tier; mechanical)
  c => agent(`Implement cross-check ${c.id} per plan.md §8 as TEST/${c.id}.wl.`,
        {label:`impl:${c.id}`, phase:'Implement', model:'sonnet', effort:'low'}),
  // stage 2: RUN it in an isolated INTERFILES dir, then INDEPENDENTLY VERIFY the numeric result
  (_, c) => agent(`Run TEST/${c.id}.wl in INTERFILES/cc_${c.id}/. Then, as an independent skeptic, corroborate any numeric result against ${c.oracles}; refute if no independent oracle agrees within combined error. Exact checks (${c.exact||'n/a'}) must hold symbolically.`,
        {label:`verify:${c.id}`, phase:'Run+Verify', model:'opus', effort:'high', schema: VERDICT }))
const verdicts = results.filter(Boolean)
const coverage = await agent(`Build the §8.1 coverage report from these verdicts: ${JSON.stringify(verdicts)}. RED if any feature row has zero passing checks or any #31–#43 refuted.`,
  {label:'coverage', phase:'Synthesize', model:'opus', effort:'high', schema: WORK})
return { green: coverage.selfGatePass, coverage, verdicts }
```

---

## §4 The inter-phase orchestration loop (run by the main Claude loop)

```
for N in 0..6:
   fire Workflow({ scriptPath: "ORCH/phaseN.js" })          # one invocation per phase, runs unattended
   await task-notification                                  # phase ran start→gate with nobody watching
   report = read PHASE_REPORT (green/red, verdicts, artifacts, disagreements)
   render report to the user
   if report.red:           HALT (red = ANY verifier refuted, could not corroborate, or did not run).
                            Relay the disagreements + artifacts. The orchestrator does NOT fix and does NOT proceed.
   elif mode == checkpoint:  wait for the user's `next` (recommended through Phase 5)
   elif mode == auto-chain:  continue to N+1 only because green (every verifier confirmed)
```
- The orchestrator (main loop) is the *only* thing that crosses phase boundaries; it never overrides a RED gate.
- On RED: the report names which verifier lens refuted and why, with the disagreeing oracle values and the artifact paths — enough to fix without rerunning blind. Fix in the worktree, then **resume** that phase via `resumeFromRunId` (cached prefix, only the fixed leg re-runs).
- Between phases you can `git diff` the `v3-build` worktree to review exactly what changed.

---

## §5 Safety rails (why unattended is safe here)

1. **Orchestrator implements nothing; gate is fail-closed (§1.7).** The orchestrator only launches/reads/relays — never edits code, never fixes a failing check. A phase is GREEN only if **every** verifier confirms (`refuted=false`); **any** refutation, non-corroboration, or missing verifier halts the run. The producer never self-grades; the gate is independent (§1.2, plan §8.5).
2. **Exactness invariant** — the decomposition is symbolic, so most verification is exact-identity (cannot be subtly wrong); only the integrated value is numeric, and it needs ≥2 independent oracles (plan §1.4 / §6.7). The exactness guard (#43) is itself a gate.
3. **Complex-lifting phases get a 3-lens panel** (4–5) — the documented silent-wrong-answer zone (L1/L3/L5).
4. **No silent truncation** — if a phase caps coverage or skips a dependency-gated check, it `log()`s it; the §8.1 matrix gate blocks an empty feature row.
5. **Worktree isolation** — main repo untouched until you merge; `OLD_CODE/` never edited.
6. **Timeouts halt, never skip** — a step that exceeds budget fails the phase loudly (§1.5).
7. **Deterministic resume** — `resumeFromRunId` re-runs only changed legs.
8. **Budget guard** — each phase script may take `budget`/`+Nk` ceilings; the verifier legs are exempt from aggressive trimming.

---

## §6 Rough envelope (set expectations, not promises)

| Phase | Dominant cost | Wall-clock driver | Suggested mode |
|---|---|---|---|
| 0 | running both old suites | toolchain (tens of min, backgrounded) | auto-chain |
| 1 | re-running Tree-A suite | toolchain | auto-chain |
| 2 | engine surgery + golden diffs | model (opus/max) + compile | **checkpoint** |
| 3 | lifting port + sandbox NIntegrate | model + NIntegrate | **checkpoint** |
| 4 | complex-lifting + 3-lens panel | model (opus/max) heavy | **checkpoint** |
| 5 | 9D fan (~260 s) + 8D VEGAS | toolchain (serial VEGAS) | **checkpoint** |
| 6 | ~43 cross-checks fanned out | toolchain, parallelized (~10 lanes) | auto-chain |

Token cost concentrates in Phases 2–4 (opus/max surgery + verifier panels) and Phase 6's verification legs; Phases 0/1/5-mechanical are cheap-tier. This is exactly the tiering in §1.3 — max only where it changes the answer.

---

## §7 What I need from you to start
1. **Mode:** checkpoint between every phase, or auto-chain-on-green (my rec: auto-chain 0→1 and 6; checkpoint 2→3→4→5)?
2. **Worktree/branch name** (default `v3-build`) and whether to commit per phase.
3. Any **per-phase token/wall-clock ceiling**.
4. Confirm I should first **write the `ORCH/phaseN.js` scripts** (still not running them) for your review, or generate-and-launch Phase 0 directly when you say go.
