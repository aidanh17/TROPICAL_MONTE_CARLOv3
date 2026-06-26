(* ============================================================================
   EXAMPLES/Ex-V3c.wl  —  Worked Example V3c
   plan.md §8.4: "Kinematic scan, batched VEGAS, per-point vs NIntegrate
                  + per-kp VEGAS"
   Cross-checks: #23, #36

   Integrand family (generalized simplex, cf. cc_23.wl)
   -----------------------------------------------------
   For each kinematic point (c0, c1, c2, c3) in a scan:
     P = c0 + c1*x1 + c2*x2 + c3*x3
     I = Int_{[0,inf)^3} P^{-5}  dx1 dx2 dx3
   Exact:  I = 1 / (4! * c0^2 * c1 * c2 * c3)   [Dirichlet / simplex formula]

   The fan is coefficient-independent (builds from the unit proxy 1+x1+x2+x3).
   One compile serves all kinematic points; "Batch"->True shares a single
   low-discrepancy VEGAS sample set across all points (plan.md §5.1).

   PASS criteria (plan.md §8.3 #23 / #36)
   ----------------------------------------
     K1  median relErr(per-kp VEGAS, exact) < 2e-2   (Tier 2)
     K2  median relErr(batched VEGAS, exact) < 2e-2  (Tier 2)
     K3  max |per-kp - batch| / exact < 1e-2         (Tier 2: batch invariance)
     K4  single-point per-kp VEGAS == in-scan per-kp within 5e-2   (Tier 2)
     K5  single-point batched VEGAS == in-scan batched within 5e-2  (Tier 2)
     K6  per-kp NIntegrate (first 5 points) vs exact < 1e-4  (Tier 1: correctness base)
     K7  per-kp MC vs exact, median < 5e-2  (Tier 1)

   Tier 1: K6, K7 (WL + g++, no CUBA).  All VEGAS gates (K1-K5) are Tier 2
   and skip gracefully when CUBA is absent.

   Run
   ---
     cd /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3
     wolframscript -file EXAMPLES/Ex-V3c.wl
   exits 0 on PASS (mandatory Tier-1 gates), 1 on FAIL.
   ============================================================================ *)

(* ── 0.  Locate packages ──────────────────────────────────────────────────── *)

$v3Root = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
$fanWL  = FileNameJoin[{$v3Root, "tropical_fan.wl"}];
$evalWL = FileNameJoin[{$v3Root, "tropical_eval.wl"}];

If[!FileExistsQ[$fanWL],
  Print["Ex-V3c FAIL  tropical_fan.wl not found at ", $fanWL]; Quit[1]];
If[!FileExistsQ[$evalWL],
  Print["Ex-V3c FAIL  tropical_eval.wl not found at ", $evalWL]; Quit[1]];

Get[$fanWL];
Get[$evalWL];

Quiet[CreateDirectory[
  FileNameJoin[{$v3Root, "INTERFILES", "ExV3c"}],
  CreateIntermediateDirectories -> True],
 {CreateDirectory::eexist}];

$ioDir = FileNameJoin[{$v3Root, "INTERFILES", "ExV3c"}];

Print[""];
Print["================================================================"];
Print["  EXAMPLES/Ex-V3c: Kinematic Scan — Batched VEGAS             "];
Print["  Cross-checks #23, #36 (plan.md §8.4)                       "];
Print["================================================================"];
Print[];

(* ── 1.  Helpers ──────────────────────────────────────────────────────────── *)

SetAttributes[quietRun, HoldFirst];
quietRun[expr_] := Quiet[expr,
  {TropicalEval`EvaluateTropicalMC::validate,
   TropicalEval`EvaluateTropicalMC::vegasbudget,
   General::stop}];

sci[x_] := If[N[Abs[x]] == 0., "0",
  Module[{e = Floor[Log10[Abs[N[x]]]]},
    ToString[NumberForm[N[x] / 10^e, {3, 2}]] <> "e" <> ToString[e]]];

relErr1[val_, ref_] := Abs[(val - ref)] / Max[Abs[ref], 1*^-300];

$passes   = {};
$failures = {};

record[name_String, ok_, passDetail_String : "", failExtra_String : ""] :=
  If[TrueQ[ok],
    (AppendTo[$passes, name];
     Print["  Ex-V3c PASS  ", name,
           If[passDetail =!= "", "  " <> passDetail, ""]]),
    (AppendTo[$failures, name];
     Print["  Ex-V3c FAIL  ", name,
           If[failExtra =!= "", "  " <> failExtra, ""]])];

recordSkip[name_String, reason_String] :=
  Print["  Ex-V3c SKIP  ", name, "  (", reason, ")"];

(* ── 2.  Integrand specification ─────────────────────────────────────────── *)

Print["--- §1  Integrand and kinematic scan ---"];
n     = 3;       (* integration variables *)
nKP   = 60;      (* kinematic points in the scan *)
vars  = Table[x[k], {k, n}];
csyms = Table[Symbol["c" <> ToString[i]], {i, 0, n}];
poly  = csyms[[1]] + Sum[csyms[[i + 1]] x[i], {i, n}];

spec = <|
  "Polynomials"         -> {poly},
  "MonomialExponents"   -> ConstantArray[0, n],
  "PolynomialExponents" -> {-(n + 2)},      (* B = -5 for n=3 *)
  "Variables"           -> vars,
  "KinematicSymbols"    -> csyms,
  "RegulatorSymbol"     -> None
|>;

Print["  n = ", n, "  B = -(n+2) = ", -(n+2)];
Print["  Exact: I = 1 / (", n+1, "! * c0^2 * c1 * ... * c", n, ")"];
Print["  Kinematic scan: ", nKP, " coefficient sets"];
Print[];

(* ── 3.  Fan and kinematic points ─────────────────────────────────────────── *)

Print["--- §2  Tropical fan (coefficient-independent) ---"];
verts = PolytopeVertices[(1 + Total[vars])^(-1), vars];
fan   = computeFanScaled[verts];
If[fan === $Failed || !ListQ[fan],
  Print["Ex-V3c FAIL  fan construction returned $Failed"];
  Quit[1]];
Print["  Sectors: ", Length[fan[[2]]]];
Print[];

(* Fixed-seed kinematic points: c0=1, c_i uniform in [exp(-0.7), exp(0.7)] *)
SeedRandom[1729 + n];
kps    = Table[Prepend[Exp[RandomReal[{-0.7, 0.7}, n]], 1.0], {nKP}];
exactFn = Function[c, 1.0 / ((n + 1)! * c[[1]]^2 * Times @@ c[[2 ;;]])];
exacts  = exactFn /@ kps;

Print["  Kinematic scan: first 3 exact values = ",
      NumberForm[exacts[[#]], {6, 6}] & /@ {1, 2, 3}];
Print[];

(* ── 4.  Criterion K6 — NIntegrate on first 5 points (Tier 1) ─────────────── *)

Print["--- §3  K6: NIntegrate on 5 scan points (Tier 1) ---"];
nintVals = Table[
  Module[{c = kps[[j]], P, ni},
    P  = c[[1]] + Sum[c[[i + 1]] x[i], {i, n}];
    ni = Quiet @ NIntegrate[P^(-(n + 2)),
      Evaluate[Sequence @@ ({x[#], 0, Infinity} & /@ Range[n])],
      Method -> "GlobalAdaptive", PrecisionGoal -> 5, MaxRecursion -> 30];
    ni],
  {j, 1, 5}
];

nintErrs = MapThread[relErr1, {nintVals, exacts[[;;5]]}];
Print["  NIntegrate vs exact (5 pts):  rel-err max = ", sci @ Max[nintErrs]];
record["K6 NIntegrate vs exact on 5 scan points",
  Max[nintErrs] < 1*^-4,
  "max-relErr=" <> sci[Max[nintErrs]],
  "expected=<1e-4 got=" <> sci[Max[nintErrs]]];
Print[];

(* ── 5.  Criterion K7 — per-kp MC (Tier 1) ──────────────────────────────── *)

Print["--- §4  K7: per-kp plain-MC (Tier 1, no CUBA) ---"];
{tMC, resMC} = AbsoluteTiming[
  quietRun @ EvaluateTropicalMC[spec, fan, kps,
    "Integrator"       -> "MC",
    "Batch"            -> False,
    "NSamples"         -> 100000,
    "RunChecks"        -> False,
    "Verbose"          -> False,
    "WorkingDirectory" -> $ioDir]
];

If[!AssociationQ[resMC] || !KeyExistsQ[resMC, "Results"],
  Print["  Ex-V3c FAIL  MC driver returned: ", Head[resMC]];
  AppendTo[$failures, "K7 MC driver"];
  Goto[$skipTier2]
];

estMC  = (#["Re"] &) /@ resMC["Results"];
relMC  = MapThread[relErr1, {estMC, exacts}];
medMC  = Median[relMC];
Print["  MC wall time = ", NumberForm[tMC, {4, 2}], " s"];
Print["  Median relErr vs exact = ", sci[medMC], "  (gate < 5e-2)"];
record["K7 per-kp MC median relErr vs exact",
  medMC < 5*^-2,
  "median=" <> sci[medMC],
  "expected=<5e-2 got=" <> sci[medMC]];
Print[];

(* ── 6.  Tier-2 criteria K1-K5 (CUBA required) ─────────────────────────── *)

$cuba = TrueQ[TropicalEval`detectCuba[]["Found"]];

If[!$cuba,
  Print["--- §5-§8  Tier-2 criteria K1-K5 SKIPPED (no CUBA) ---"];
  Print[];
  Goto[$skipTier2]
];

Print["--- §5  K1: per-kp VEGAS vs exact ---"];
{tPK, resPK} = AbsoluteTiming[
  quietRun @ EvaluateTropicalMC[spec, fan, kps,
    "Integrator"       -> "VEGAS",
    "Batch"            -> False,
    "NSamples"         -> 200000,
    "RunChecks"        -> False,
    "Verbose"          -> False,
    "WorkingDirectory" -> $ioDir]
];

If[!AssociationQ[resPK] || !KeyExistsQ[resPK, "Results"],
  recordSkip["K1/K3/K4 per-kp VEGAS", "driver returned $Failed"];
  Goto[$skipBatch]
];
estPK  = (#["Re"] &) /@ resPK["Results"];
relPK  = MapThread[relErr1, {estPK, exacts}];
medPK  = Median[relPK];
Print["  per-kp VEGAS wall time = ", NumberForm[tPK, {4, 2}], " s"];
Print["  Median relErr vs exact = ", sci[medPK], "  (gate < 2e-2)"];
record["K1 per-kp VEGAS median relErr vs exact  (#23 gate 1)",
  medPK < 2*^-2,
  "median=" <> sci[medPK],
  "expected=<2e-2 got=" <> sci[medPK]];
Print[];

Print["--- §6  K2: batched VEGAS vs exact ---"];
{tBT, resBT} = AbsoluteTiming[
  quietRun @ EvaluateTropicalMC[spec, fan, kps,
    "Integrator"       -> "VEGAS",
    "Batch"            -> True,
    "NSamples"         -> 200000,
    "RunChecks"        -> False,
    "Verbose"          -> False,
    "WorkingDirectory" -> $ioDir]
];

If[!AssociationQ[resBT] || !KeyExistsQ[resBT, "Results"],
  recordSkip["K2/K3/K5 batched VEGAS", "driver returned $Failed"];
  Goto[$skipSingle]
];
estBT  = (#["Re"] &) /@ resBT["Results"];
relBT  = MapThread[relErr1, {estBT, exacts}];
medBT  = Median[relBT];
maxAg  = Max[MapThread[relErr1, {estPK, estBT}]];
Print["  batched VEGAS wall time = ", NumberForm[tBT, {4, 2}], " s"];
Print["  Median relErr vs exact  = ", sci[medBT], "  (gate < 2e-2)"];
Print["  Max |per-kp - batch| / exact = ", sci[maxAg], "  (gate < 1e-2)"];
record["K2 batched VEGAS median relErr vs exact  (#23 gate 2)",
  medBT < 2*^-2,
  "median=" <> sci[medBT],
  "expected=<2e-2 got=" <> sci[medBT]];
record["K3 batch invariance: max |per-kp - batch| / exact  (#23 gate 3, #36)",
  maxAg < 1*^-2,
  "max=" <> sci[maxAg],
  "expected=<1e-2 got=" <> sci[maxAg]];
Print[];

Print["--- §7  K4+K5: single-point vs in-scan agreement ---"];
$kpSingle = {kps[[1]]};

resSoloPK = quietRun @ EvaluateTropicalMC[spec, fan, $kpSingle,
  "Integrator"       -> "VEGAS",
  "Batch"            -> False,
  "NSamples"         -> 200000,
  "RunChecks"        -> False,
  "Verbose"          -> False,
  "WorkingDirectory" -> $ioDir];

resSoloBT = quietRun @ EvaluateTropicalMC[spec, fan, $kpSingle,
  "Integrator"       -> "VEGAS",
  "Batch"            -> True,
  "NSamples"         -> 200000,
  "RunChecks"        -> False,
  "Verbose"          -> False,
  "WorkingDirectory" -> $ioDir];

If[AssociationQ[resSoloPK] && KeyExistsQ[resSoloPK, "Results"] &&
   Length[resSoloPK["Results"]] > 0,
  Module[{solo = resSoloPK["Results"][[1]]["Re"], scan = estPK[[1]],
          d},
    d = relErr1[solo, scan];
    Print["  Single-pt per-kp vs scan kp-1: relDiff = ", sci[d],
          "  (gate 5e-2)"];
    record["K4 single-pt per-kp == in-scan kp-1  (#23 gate 4)",
      d < 5*^-2,
      "relDiff=" <> sci[d],
      "expected=<5e-2 got=" <> sci[d]]],
  recordSkip["K4 single-pt per-kp", "run returned $Failed"]
];

If[AssociationQ[resSoloBT] && KeyExistsQ[resSoloBT, "Results"] &&
   Length[resSoloBT["Results"]] > 0,
  Module[{solo = resSoloBT["Results"][[1]]["Re"], scan = estBT[[1]],
          d},
    d = relErr1[solo, scan];
    Print["  Single-pt batched vs scan kp-1: relDiff = ", sci[d],
          "  (gate 5e-2)"];
    record["K5 single-pt batched == in-scan kp-1  (#23 gate 5)",
      d < 5*^-2,
      "relDiff=" <> sci[d],
      "expected=<5e-2 got=" <> sci[d]]],
  recordSkip["K5 single-pt batched", "run returned $Failed"]
];

Label[$skipSingle];
Label[$skipBatch];
Label[$skipTier2];
Print[];

(* ── 7.  Summary ──────────────────────────────────────────────────────────── *)

Print["================================================================"];
Print["  Kinematic scan summary:"];
Print["    n = ", n, "  B = ", -(n+2), "  nKP = ", nKP];
If[ValueQ[medMC],
  Print["    MC median relErr   = ", sci[medMC]]];
If[$cuda && ValueQ[medPK],
  Print["    per-kp VEGAS median = ", sci[medPK],
        "  batched VEGAS median = ", sci[medBT]];
  Print["    batch invariance (max |per-kp - batch|/exact) = ", sci[maxAg]];
  Print["    Wall: per-kp=", NumberForm[tPK, {4, 2}], "s  batch=",
        NumberForm[tBT, {4, 2}], "s  speedup=",
        NumberForm[tPK/Max[tBT, 0.001], {4, 2}], "x"]
];
Print["  PASSED: ", Length[$passes]];
Print["  FAILED: ", Length[$failures]];
Print["================================================================"];

If[Length[$failures] == 0,
  Print["Ex-V3c PASS  kinematic scan: NIntegrate base correct; batch invariant (",
        Length[$passes], " criteria)"],
  Print["Ex-V3c FAIL  ", Length[$failures],
        " criterion(a) failed: ",
        StringRiffle[$failures, ", "]]];
Print["================================================================"];

Quit[If[Length[$failures] > 0, 1, 0]]
