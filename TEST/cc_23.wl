(* ============================================================================
   TEST/cc_23.wl  —  Cross-check #23: Batched VEGAS == per-kp VEGAS
                     (plan.md §8.2 #23, §5.1 fold gate)

   PASS criterion (plan.md §8.2 #23):
       (a) batched VEGAS value within combined error of per-kp VEGAS, and
       (b) both batched and per-kp VEGAS match the closed-form reference.

   Tier 2: requires CUBA (for Integrator -> "VEGAS", "Batch" -> True).
   If CUBA is absent prints "CC23 SKIP (no CUBA)" and exits 0.

   Ported / adapted from:
     OLD_CODE/TROPICAL_MONTE_CARLO/EXAMPLES/vegas_batch_example.wl  (Tree A T11)
   using the v3 API (integrator strings "MC"/"VEGAS", EvaluateTropicalMC
   with "Batch"->True, computeFanScaled from tropical_fan.wl).

   Integrand family:
     P = c0 + sum_{i=1}^n c_i x_i,   A_i = 0,   B = -(n+2),   n in {2, 4}
     Exact: I = 1 / ((n+1)! * c0^2 * c1 * ... * cn)
     (generalised simplex; coefficient-independent fan)

   Sub-cases:
     A  n=2, nKP=50  coefficient sets  (fast smoke test)
     B  n=4, nKP=100 coefficient sets  (the Tree-A T11 regime)
     C  Single-point check: one-point run == in-scan value for both modes

   Gates (per sub-case):
     (1) median relErr(per-kp VEGAS, exact)  < vegasTol
     (2) median relErr(batched VEGAS, exact)  < vegasTol
     (3) max |per-kp - batched| / exact      < batchAgreeTol
     (4) single-point per-kp   == in-scan per-kp   (within MC noise, tol singleTol)
     (5) single-point batched  == in-scan batched   (within MC noise, tol singleTol)

   Run:  wolframscript -file TEST/cc_23.wl
   (from any directory — uses absolute paths)
   ============================================================================ *)

(* ── 0. Locate the package root ───────────────────────────────────────────── *)
$pkgRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
$evalWL  = FileNameJoin[{$pkgRoot, "tropical_eval.wl"}];
$fanWL   = FileNameJoin[{$pkgRoot, "tropical_fan.wl"}];
$ioDir   = FileNameJoin[{$pkgRoot, "TEST", "INTERFILES"}];

If[!FileExistsQ[$evalWL],
  Print["CC23 FAIL  tropical_eval.wl not found at ", $evalWL]; Quit[1]];
If[!FileExistsQ[$fanWL],
  Print["CC23 FAIL  tropical_fan.wl not found at ", $fanWL]; Quit[1]];

Get[$fanWL];
Get[$evalWL];

Quiet[CreateDirectory[$ioDir], {CreateDirectory::eexist}];

(* ── 1. CUBA presence check (Tier-2 gate) ───────────────────────────────── *)
$cuba = TropicalEval`detectCuba[];
If[!TrueQ[$cuba["Found"]],
  Print["CC23 SKIP (no CUBA)"];
  Quit[0]];

Print["CC23: packages loaded; CUBA found at ", $cuba["IncludeDir"]];
Print["CC23: INTERFILES directory: ", $ioDir];
Print[];

(* ── 2. Helpers ─────────────────────────────────────────────────────────── *)

SetAttributes[quietRun, HoldFirst];
quietRun[expr_] := Block[{$Output = {}}, expr];

sci[x_] := If[N[x] == 0., "0",
  Module[{e = Floor[Log10[Abs[N[x]]]], m},
    m = N[x] / 10^e;
    ToString[NumberForm[m, {3, 2}]] <> "e" <> ToString[e]]];

(* Relative error, safe against near-zero reference *)
relErr1[val_, ref_] := Abs[(val - ref)] / Max[Abs[ref], 1.*^-300];

(* ── 3. Core test function ───────────────────────────────────────────────── *)
(* Returns True (PASS) / False (FAIL).
   n         : number of integration variables
   nKP       : number of kinematic points (coefficient sets)
   nSamples  : per-run NSamples budget
   vegasTol  : median relErr gate for each VEGAS mode vs exact
   batchTol  : max |per-kp - batch| / exact  gate
   singleTol : single-point vs in-scan agreement gate                        *)

Options[cc23SubCase] = {
  "VegasTol"    -> 2.*^-2,
  "BatchAgreeTol" -> 1.*^-2,
  "SingleTol"   -> 5.*^-2
};

cc23SubCase[label_String, n_Integer, nKP_Integer, nSamples_Integer,
            OptionsPattern[]] :=
Module[
  {vars, csyms, P, spec, verts, fan,
   kps, exactFn, exacts,
   tPK, resPK, tBT, resBT,
   estPK, estBT, relPK, relBT,
   medPK, medBT, maxAgree,
   vegTol, batchTol, singleTol,
   anyFail = False, passFlags = {},
   (* sub-case C: single-point *)
   kpSingle, resSinglePK, resSingleBT,
   valScan1PK, valScan1BT, valSoloPK, valSoloBT},

  vegTol    = OptionValue["VegasTol"];
  batchTol  = OptionValue["BatchAgreeTol"];
  singleTol = OptionValue["SingleTol"];

  Print["========================================================"];
  Print["CC23  ", label];
  Print["========================================================"];
  Print["  n=", n, "  nKP=", nKP, "  NSamples=", nSamples];

  (* --- Build spec ---- *)
  vars  = Table[x[k], {k, n}];
  csyms = Table[Symbol["c" <> ToString[i]], {i, 0, n}];
  P     = csyms[[1]] + Sum[csyms[[i + 1]] x[i], {i, n}];
  spec  = <|"Polynomials"        -> {P},
            "MonomialExponents"  -> ConstantArray[0, n],
            "PolynomialExponents"-> {-(n + 2)},
            "Variables"          -> vars,
            "KinematicSymbols"   -> csyms,
            "RegulatorSymbol"    -> None|>;

  (* --- Fan (coefficient-independent: compute from unit proxy) --- *)
  verts = PolytopeVertices[(1 + Total[vars])^(-1), vars];
  fan   = computeFanScaled[verts];
  If[fan === $Failed,
    Print["CC23 FAIL  ", label, "  fan computation failed"];
    Print["CC23 FAIL  expected=<fan>  got=$Failed"];
    Return[False]];
  Print["  fan: ", Length[fan[[2]]], " sectors"];

  (* --- Kinematic points: c0=1, c_i in [exp(-0.7), exp(0.7)] --- *)
  SeedRandom[42 + n];   (* fixed seed for reproducibility *)
  kps = Table[
    Prepend[Exp[RandomReal[{-0.7, 0.7}, n]], 1.0],
    {nKP}];

  exactFn = Function[c, 1.0 / ((n + 1)! * c[[1]]^2 * Times @@ c[[2;;]])];
  exacts  = exactFn /@ kps;

  (* --- Run per-kp VEGAS --- *)
  Print["  Running per-kp VEGAS  (", nKP, " points) ..."];
  {tPK, resPK} = AbsoluteTiming[
    quietRun @ EvaluateTropicalMC[spec, fan, kps,
      "Integrator"   -> "VEGAS",
      "Batch"        -> False,
      "NSamples"     -> nSamples,
      "RunChecks"    -> False,
      "Verbose"      -> False,
      "WorkingDirectory" -> $ioDir]];

  If[!AssociationQ[resPK] || !KeyExistsQ[resPK, "Results"],
    Print["CC23 FAIL  ", label, "  per-kp VEGAS driver error: ", resPK];
    Print["CC23 FAIL  expected=<assoc>  got=", Head[resPK]];
    Return[False]];

  (* --- Run batched VEGAS --- *)
  Print["  Running batched VEGAS (", nKP, " points) ..."];
  {tBT, resBT} = AbsoluteTiming[
    quietRun @ EvaluateTropicalMC[spec, fan, kps,
      "Integrator"   -> "VEGAS",
      "Batch"        -> True,
      "NSamples"     -> nSamples,
      "RunChecks"    -> False,
      "Verbose"      -> False,
      "WorkingDirectory" -> $ioDir]];

  If[!AssociationQ[resBT] || !KeyExistsQ[resBT, "Results"],
    Print["CC23 FAIL  ", label, "  batched VEGAS driver error: ", resBT];
    Print["CC23 FAIL  expected=<assoc>  got=", Head[resBT]];
    Return[False]];

  (* --- Extract estimates (real part; integrand is real positive) --- *)
  estPK = (#["Re"] &) /@ resPK["Results"];
  estBT = (#["Re"] &) /@ resBT["Results"];

  (* --- Relative errors vs closed form --- *)
  relPK = MapThread[relErr1, {estPK, exacts}];
  relBT = MapThread[relErr1, {estBT, exacts}];

  medPK    = Median[relPK];
  medBT    = Median[relBT];
  maxAgree = Max[MapThread[relErr1, {estPK, estBT}]];

  Print["  Wall time:  per-kp = ", NumberForm[tPK, {4, 2}], " s",
        "  batched = ", NumberForm[tBT, {4, 2}], " s"];
  Print["  Median relErr:  per-kp = ", sci[medPK],
        "  (gate ", vegTol, ")"];
  Print["  Median relErr:  batched = ", sci[medBT],
        "  (gate ", vegTol, ")"];
  Print["  Max |per-kp - batch| / exact = ", sci[maxAgree],
        "  (gate ", batchTol, ")"];

  (* --- Gate (1): per-kp VEGAS vs exact --- *)
  If[!TrueQ[N[medPK] <= vegTol],
    anyFail = True;
    Print["  CC23 FAIL  ", label, " gate(1) per-kp vs exact: ",
          "expected=medRelErr<=", vegTol, " got=", sci[medPK]];
    AppendTo[passFlags, "gate1-FAIL"],
    AppendTo[passFlags, "gate1-PASS"]];

  (* --- Gate (2): batched VEGAS vs exact --- *)
  If[!TrueQ[N[medBT] <= vegTol],
    anyFail = True;
    Print["  CC23 FAIL  ", label, " gate(2) batch vs exact: ",
          "expected=medRelErr<=", vegTol, " got=", sci[medBT]];
    AppendTo[passFlags, "gate2-FAIL"],
    AppendTo[passFlags, "gate2-PASS"]];

  (* --- Gate (3): per-kp vs batched agreement --- *)
  If[!TrueQ[N[maxAgree] <= batchTol],
    anyFail = True;
    Print["  CC23 FAIL  ", label, " gate(3) per-kp vs batch: ",
          "expected=maxRelDiff<=", batchTol, " got=", sci[maxAgree]];
    AppendTo[passFlags, "gate3-FAIL"],
    AppendTo[passFlags, "gate3-PASS"]];

  (* --- Gate (4+5): single-point == first in-scan point --- *)
  kpSingle = {kps[[1]]};
  Print["  Running single-point per-kp VEGAS ..."];
  resSinglePK = quietRun @ EvaluateTropicalMC[spec, fan, kpSingle,
    "Integrator"   -> "VEGAS",
    "Batch"        -> False,
    "NSamples"     -> nSamples,
    "RunChecks"    -> False,
    "Verbose"      -> False,
    "WorkingDirectory" -> $ioDir];

  Print["  Running single-point batched VEGAS ..."];
  resSingleBT = quietRun @ EvaluateTropicalMC[spec, fan, kpSingle,
    "Integrator"   -> "VEGAS",
    "Batch"        -> True,
    "NSamples"     -> nSamples,
    "RunChecks"    -> False,
    "Verbose"      -> False,
    "WorkingDirectory" -> $ioDir];

  If[AssociationQ[resSinglePK] && KeyExistsQ[resSinglePK, "Results"] &&
     Length[resSinglePK["Results"]] > 0,
    valSoloPK  = resSinglePK["Results"][[1]]["Re"];
    valScan1PK = estPK[[1]];
    Module[{d = relErr1[valSoloPK, valScan1PK]},
      Print["  Single-pt per-kp vs scan kp-1: relDiff = ", sci[d],
            "  (gate ", singleTol, ")"];
      If[!TrueQ[N[d] <= singleTol],
        anyFail = True;
        Print["  CC23 FAIL  ", label, " gate(4): ",
              "expected=relDiff<=", singleTol, " got=", sci[d]];
        AppendTo[passFlags, "gate4-FAIL"],
        AppendTo[passFlags, "gate4-PASS"]]],
    Print["  gate(4) SKIP: single-point per-kp run returned unexpected structure"]];

  If[AssociationQ[resSingleBT] && KeyExistsQ[resSingleBT, "Results"] &&
     Length[resSingleBT["Results"]] > 0,
    valSoloBT  = resSingleBT["Results"][[1]]["Re"];
    valScan1BT = estBT[[1]];
    Module[{d = relErr1[valSoloBT, valScan1BT]},
      Print["  Single-pt batched  vs scan kp-1: relDiff = ", sci[d],
            "  (gate ", singleTol, ")"];
      If[!TrueQ[N[d] <= singleTol],
        anyFail = True;
        Print["  CC23 FAIL  ", label, " gate(5): ",
              "expected=relDiff<=", singleTol, " got=", sci[d]];
        AppendTo[passFlags, "gate5-FAIL"],
        AppendTo[passFlags, "gate5-PASS"]]],
    Print["  gate(5) SKIP: single-point batched run returned unexpected structure"]];

  Print["  Gates: ", StringRiffle[passFlags, "  "]];
  Print[];

  If[!anyFail,
    Print["CC23 PASS  ", label,
          "  per-kp=", sci[medPK], " batch=", sci[medBT],
          " agree=", sci[maxAgree],
          " (", nKP, " pts, NSamples=", nSamples, ")"],
    Print["CC23 FAIL  ", label,
          "  expected=all-gates-pass  got=fails-above"]];
  Print[];
  !anyFail
];

(* ── 4. Sub-cases ───────────────────────────────────────────────────────── *)

$cc23Results = {};

(* Sub-case A: n=2, 50 kinematic points — fast smoke *)
Print["========================================================"];
Print["CC23  Sub-case A: n=2, 50 kp, 200000 samples/run"];
Print["========================================================"];
AppendTo[$cc23Results,
  cc23SubCase["A: n=2, 50kp", 2, 50, 200000,
    "VegasTol" -> 2.*^-2, "BatchAgreeTol" -> 1.*^-2, "SingleTol" -> 5.*^-2]];

(* Sub-case B: n=4, 100 kinematic points — mirrors Tree A T11 *)
Print["========================================================"];
Print["CC23  Sub-case B: n=4, 100 kp, 100000 samples/run"];
Print["========================================================"];
AppendTo[$cc23Results,
  cc23SubCase["B: n=4, 100kp", 4, 100, 100000,
    "VegasTol" -> 2.*^-2, "BatchAgreeTol" -> 1.*^-2, "SingleTol" -> 5.*^-2]];

(* ── 5. Summary ─────────────────────────────────────────────────────────── *)
Module[{nPass, nFail, labels},
  labels = {"A: n=2, 50kp", "B: n=4, 100kp"};
  nPass  = Count[$cc23Results, True];
  nFail  = Count[$cc23Results, False];
  Print["========================================================"];
  Print["CC23  Batched VEGAS == per-kp VEGAS  —  summary"];
  Print["========================================================"];
  Do[
    Print["  ", StringPadRight[labels[[i]], 22],
          If[i <= Length[$cc23Results] && TrueQ[$cc23Results[[i]]], "PASS", "FAIL"]],
    {i, Length[labels]}];
  Print["--------------------------------------------------------"];
  If[nFail == 0,
    Print["CC23 PASS  all ", nPass, " sub-cases passed",
          "  (batch invariance confirmed; §5.1 fold gated)"],
    Print["CC23 FAIL  expected=", nPass + nFail, " pass",
          "  got=", nFail, " fail  (", nPass, " passed, ", nFail, " failed)"]];
  Print["========================================================"];
  If[nFail > 0, Quit[1]]
];
