# ORCH/ — per-phase unattended ultracode workflows

Implementation of `ORCHESTRATION.md`. One self-contained Workflow script per plan.md phase.

## How to run (the inter-phase loop, §4)

The **orchestrator is the main Claude loop** (Opus @ medium). It implements NOTHING (§1.7):
it only fires a phase, reads the `PHASE_REPORT`, enforces the fail-closed gate, and relays.

```
for N in 0..6:
   Workflow({ scriptPath: "ORCH/phaseN.js" })   # runs unattended in background
   await task-notification
   read return value { green, work, verdicts }
   if !green:  HALT. Relay which verifier lens refuted + disagreements + artifacts.
   elif checkpoint phase:  wait for user `next`
   elif auto-chain phase:  fire phase N+1
```

**Mode chosen for this build:** auto-chain `0→1` and `6`; **checkpoint** `2→3→4→5`.

**Resume** a edited/interrupted phase: `Workflow({ scriptPath:"ORCH/phaseN.js", resumeFromRunId:"wf_…" })`
— the unchanged prefix of `agent()` calls returns cached; only new/edited legs re-run.

## The gate (identical in every script, §2 / §1.7)

```
green = work.selfGatePass && allRan && refuted === 0
```
FAIL-CLOSED: a producer cannot self-declare green; ANY verifier that refutes — finds something
wrong OR cannot independently corroborate a numeric value — OR fails to run at all → RED → halt.
No majority vote. Complex-lifting phases (4,5) run a 3-lens panel and require ALL lenses to confirm.

## Shared harness (inlined verbatim into each phaseN.js)

- `WORK` schema — producer's structured result (summary, filesChanged, commands, numerics, selfGatePass, artifacts).
- `VERDICT` schema — each verifier's result (lens, refuted, confidence, independentChecks, disagreements, notes).
- `verifyPrompt(lens, work)` — "you are an INDEPENDENT VERIFIER … default refuted=true … re-run the gate
  yourself … run the exact symbolic identities … corroborate every numeric vs ≥2 independent oracles."
- `runPhase({ producePrompt, produceModel, produceEffort, verifierLenses, commitMsg })` — Produce → Verify →
  deterministic gate. On green, commits the `v3-build` worktree (mechanical, cheap-tier; part of the build,
  not orchestrator action) when `commitMsg` is set.

## Tiering (§1.3)

| Leg | model | effort |
|---|---|---|
| Core surgery, complex-lifting, exactness audit (phases 2/3/4) | opus | max/xhigh |
| Independent adversarial verifier + judge panels | opus | high–max |
| Cross-check authoring, synthesis | opus/sonnet | medium–high |
| Mechanical ports, deletions, running suites, golden capture | sonnet/haiku | low–medium |

## Worktree (§1.4, decided this session)

Root is NOT a git repo (only the OLD_CODE subtrees are). Phase 0 therefore `git init`s at the
project root, gitignores `OLD_CODE/` + `INTERFILES/`, and builds on branch **`v3-build`**, committing
per green phase. `OLD_CODE/` is never tracked and never edited; `tropical_fan.wl` only ever gets the
additive `computeFanScaled`.
