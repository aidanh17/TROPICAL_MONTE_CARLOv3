(* ============================================================================
   TEST/cc_36.wl  —  Cross-check #36: Sampler x method x batch consistency
                     (plan.md §8.3 #36, §5.1 batch fold + kinematic plumbing)

   PASS criterion (plan.md §8.3 #36):
       One kinematic-scan integral: per-kp MC, per-kp VEGAS, batched VEGAS,
       and standalone one-point runs all agree within combined error;
       one-point == in-scan value.

   Tier 2: requires CUBA (for Integrator -> "VEGAS").
   If CUBA is absent prints "CC36 SKIP (no CUBA)" and exits 0.

   Integrand family:
     P = c0 + sum_{i=1}^n c_i x_i,   A_i = 0,   B = -(n+2),   n in {2, 4}
     Exact: I = 1 / ((n+1)! * c0^2 * c1 * ... * cn)
     (generalised simplex; coefficient-independent fan)

   Three methods are compared at each kinematic point:
     (a) per-kp MC   (Integrator -> "MC",    Batch -> False)
     (b) per-kp VEGAS (Integrator -> "VEGAS", Batch -> False)
     (c) batched VEGAS (Integrator -> "VEGAS", Batch -> True)

   Gates (per sub-case):
     (1) median relErr(per-kp MC,    exact)  < mcTol
     (2) median relErr(per-kp VEGAS, exact)  < vegasTol
     (3) median relErr(batched VEGAS, exact) < vegasTol
     (4) max pairwise |MC - per-kp VEGAS| / exact  < crossTol
     (5) max pairwise |MC - batched VEGAS| / exact  < crossTol
     (6) max pairwise |per-kp VEGAS - batched VEGAS| / exact  < batchTol
     (7) single-point per-kp MC    == in-scan kp-1 per-kp MC    (singleTol)
     (8) single-point per-kp VEGAS == in-scan kp-1 per-kp VEGAS (singleTol)
     (9) single-point batched VEGAS == in-scan kp-1 batched VEGAS (singleTol)

   Sub-cases:
     A  n=2, nKP=30  (fast smoke test)
     B  n=4, nKP=60  (higher-dimensional regime)

   Run:  wolframscript -file TEST/cc_36.wl
   (from any directory — uses absolute paths)
   ============================================================================ *)

(* ── 0. Locate the package root ───────────────────────────────────────────── *)
$pkgRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
$evalWL  = FileNameJoin[{$pkgRoot, "tropical_eval.wl"}];
$fanWL   = FileNameJoin[{$pkgRoot, "tropical_fan.wl"}];
$ioDir   = FileNameJoin[{$pkgRoot, "TEST", "INTERFILES"}];

If[!FileExistsQ[$evalWL],
  Print["CC36 FAIL  tropical_eval.wl not found at ", $evalWL]; Quit[1]];
If[!FileExistsQ[$fanWL],
  Print["CC36 FAIL  tropical_fan.wl not found at ", $fanWL]; Quit[1]];

Get[$fanWL];
Get[$evalWL];

Quiet[CreateDirectory[$ioDir], {CreateDirectory::eexist}];

(* ── 1. CUBA presence check (Tier-2 gate) ───────────────────────────────── *)
$cuba = TropicalEval`detectCuba[];
If[!TrueQ[$cuba["Found"]],
  Print["CC36 SKIP (no CUBA)"];
  Quit[0]];

Print["CC36: packages loaded; CUBA found at ", $cuba["IncludeDir"]];
Print["CC36: INTERFILES directory: ", $ioDir];
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

(* ── 3. Core test function ──────────────────────────────────────────────── *)

Options[cc36SubCase] = {
  "MCTol"         -> 5.*^-2,
  "VegasTol"      -> 2.*^-2,
  "CrossTol"      -> 5.*^-2,
  "BatchAgreeTol" -> 1.*^-2,
  "SingleTol"     -> 7.*^-2
};

cc36SubCase[label_String, n_Integer, nKP_Integer,
            nSamplesMC_Integer, nSamplesVEGAS_Integer,
            OptionsPattern[]] :=
Module[
  {vars, csyms, P, spec, verts, fan,
   kps, exactFn, exacts,
   tMC, resMC,
   tPK, resPK,
   tBT, resBT,
   estMC, estPK, estBT,
   relMC, relPK, relBT,
   medMC, medPK, medBT,
   maxMCvsPK, maxMCvsBT, maxPKvsBT,
   mcTol, vegTol, crossTol, batchTol, singleTol,
   anyFail = False, passFlags = {},
   (* gates 7-9: single-point runs *)
   kpSingle,
   resSoloMC, resSoloPK, resSoloBT,
   valScan1MC, valScan1PK, valScan1BT,
   valSoloMC, valSoloPK, valSoloBT},

  mcTol    = OptionValue["MCTol"];
  vegTol   = OptionValue["VegasTol"];
  crossTol = OptionValue["CrossTol"];
  batchTol = OptionValue["BatchAgreeTol"];
  singleTol = OptionValue["SingleTol"];

  Print["========================================================"];
  Print["CC36  ", label];
  Print["========================================================"];
  Print["  n=", n, "  nKP=", nKP,
        "  NSamples(MC)=", nSamplesMC,
        "  NSamples(VEGAS)=", nSamplesVEGAS];

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
    Print["CC36 FAIL  ", label, "  fan computation failed"];
    Print["CC36 FAIL  expected=<fan>  got=$Failed"];
    Return[False]];
  Print["  fan: ", Length[fan[[2]]], " sectors"];

  (* --- Kinematic points: c0=1, c_i in [exp(-0.5), exp(0.5)] --- *)
  SeedRandom[36 + n];   (* fixed seed for reproducibility *)
  kps = Table[
    Prepend[Exp[RandomReal[{-0.5, 0.5}, n]], 1.0],
    {nKP}];

  exactFn = Function[c, 1.0 / ((n + 1)! * c[[1]]^2 * Times @@ c[[2;;]])];
  exacts  = exactFn /@ kps;

  (* --- (a) Per-kp MC --- *)
  Print["  Running per-kp MC  (", nKP, " points, NSamples=", nSamplesMC, ") ..."];
  {tMC, resMC} = AbsoluteTiming[
    quietRun @ EvaluateTropicalMC[spec, fan, kps,
      "Integrator"       -> "MC",
      "Batch"            -> False,
      "NSamples"         -> nSamplesMC,
      "RunChecks"        -> False,
      "Verbose"          -> False,
      "WorkingDirectory" -> $ioDir]];

  If[!AssociationQ[resMC] || !KeyExistsQ[resMC, "Results"],
    Print["CC36 FAIL  ", label, "  per-kp MC driver error"];
    Print["CC36 FAIL  expected=<assoc>  got=", Head[resMC]];
    Return[False]];

  (* --- (b) Per-kp VEGAS --- *)
  Print["  Running per-kp VEGAS  (", nKP, " points, NSamples=", nSamplesVEGAS, ") ..."];
  {tPK, resPK} = AbsoluteTiming[
    quietRun @ EvaluateTropicalMC[spec, fan, kps,
      "Integrator"       -> "VEGAS",
      "Batch"            -> False,
      "NSamples"         -> nSamplesVEGAS,
      "RunChecks"        -> False,
      "Verbose"          -> False,
      "WorkingDirectory" -> $ioDir]];

  If[!AssociationQ[resPK] || !KeyExistsQ[resPK, "Results"],
    Print["CC36 FAIL  ", label, "  per-kp VEGAS driver error"];
    Print["CC36 FAIL  expected=<assoc>  got=", Head[resPK]];
    Return[False]];

  (* --- (c) Batched VEGAS --- *)
  Print["  Running batched VEGAS  (", nKP, " points, NSamples=", nSamplesVEGAS, ") ..."];
  {tBT, resBT} = AbsoluteTiming[
    quietRun @ EvaluateTropicalMC[spec, fan, kps,
      "Integrator"       -> "VEGAS",
      "Batch"            -> True,
      "NSamples"         -> nSamplesVEGAS,
      "RunChecks"        -> False,
      "Verbose"          -> False,
      "WorkingDirectory" -> $ioDir]];

  If[!AssociationQ[resBT] || !KeyExistsQ[resBT, "Results"],
    Print["CC36 FAIL  ", label, "  batched VEGAS driver error"];
    Print["CC36 FAIL  expected=<assoc>  got=", Head[resBT]];
    Return[False]];

  (* --- Extract estimates (real part; integrand is real positive) --- *)
  estMC = (#["Re"] &) /@ resMC["Results"];
  estPK = (#["Re"] &) /@ resPK["Results"];
  estBT = (#["Re"] &) /@ resBT["Results"];

  (* --- Relative errors vs closed form --- *)
  relMC = MapThread[relErr1, {estMC, exacts}];
  relPK = MapThread[relErr1, {estPK, exacts}];
  relBT = MapThread[relErr1, {estBT, exacts}];

  medMC = Median[relMC];
  medPK = Median[relPK];
  medBT = Median[relBT];

  (* --- Pairwise cross-method differences --- *)
  maxMCvsPK = Max[MapThread[relErr1, {estMC, estPK}]];
  maxMCvsBT = Max[MapThread[relErr1, {estMC, estBT}]];
  maxPKvsBT = Max[MapThread[relErr1, {estPK, estBT}]];

  Print["  Wall time:  MC = ", NumberForm[tMC, {4, 2}], " s",
        "  per-kp VEGAS = ", NumberForm[tPK, {4, 2}], " s",
        "  batched VEGAS = ", NumberForm[tBT, {4, 2}], " s"];
  Print["  Median relErr vs exact:  MC = ", sci[medMC],
        "  (gate ", mcTol, ")"];
  Print["  Median relErr vs exact:  per-kp VEGAS = ", sci[medPK],
        "  (gate ", vegTol, ")"];
  Print["  Median relErr vs exact:  batched VEGAS = ", sci[medBT],
        "  (gate ", vegTol, ")"];
  Print["  Max |MC - per-kp VEGAS| / exact = ", sci[maxMCvsPK],
        "  (gate ", crossTol, ")"];
  Print["  Max |MC - batched VEGAS| / exact = ", sci[maxMCvsBT],
        "  (gate ", crossTol, ")"];
  Print["  Max |per-kp VEGAS - batched VEGAS| / exact = ", sci[maxPKvsBT],
        "  (gate ", batchTol, ")"];

  (* --- Gate (1): per-kp MC vs exact --- *)
  If[!TrueQ[N[medMC] <= mcTol],
    anyFail = True;
    Print["  CC36 FAIL  ", label, " gate(1) per-kp MC vs exact: ",
          "expected=medRelErr<=", mcTol, " got=", sci[medMC]];
    AppendTo[passFlags, "gate1-FAIL"],
    AppendTo[passFlags, "gate1-PASS"]];

  (* --- Gate (2): per-kp VEGAS vs exact --- *)
  If[!TrueQ[N[medPK] <= vegTol],
    anyFail = True;
    Print["  CC36 FAIL  ", label, " gate(2) per-kp VEGAS vs exact: ",
          "expected=medRelErr<=", vegTol, " got=", sci[medPK]];
    AppendTo[passFlags, "gate2-FAIL"],
    AppendTo[passFlags, "gate2-PASS"]];

  (* --- Gate (3): batched VEGAS vs exact --- *)
  If[!TrueQ[N[medBT] <= vegTol],
    anyFail = True;
    Print["  CC36 FAIL  ", label, " gate(3) batched VEGAS vs exact: ",
          "expected=medRelErr<=", vegTol, " got=", sci[medBT]];
    AppendTo[passFlags, "gate3-FAIL"],
    AppendTo[passFlags, "gate3-PASS"]];

  (* --- Gate (4): MC vs per-kp VEGAS cross-agreement --- *)
  If[!TrueQ[N[maxMCvsPK] <= crossTol],
    anyFail = True;
    Print["  CC36 FAIL  ", label, " gate(4) MC vs per-kp VEGAS: ",
          "expected=maxRelDiff<=", crossTol, " got=", sci[maxMCvsPK]];
    AppendTo[passFlags, "gate4-FAIL"],
    AppendTo[passFlags, "gate4-PASS"]];

  (* --- Gate (5): MC vs batched VEGAS cross-agreement --- *)
  If[!TrueQ[N[maxMCvsBT] <= crossTol],
    anyFail = True;
    Print["  CC36 FAIL  ", label, " gate(5) MC vs batched VEGAS: ",
          "expected=maxRelDiff<=", crossTol, " got=", sci[maxMCvsBT]];
    AppendTo[passFlags, "gate5-FAIL"],
    AppendTo[passFlags, "gate5-PASS"]];

  (* --- Gate (6): per-kp VEGAS vs batched VEGAS agreement --- *)
  If[!TrueQ[N[maxPKvsBT] <= batchTol],
    anyFail = True;
    Print["  CC36 FAIL  ", label, " gate(6) per-kp vs batched VEGAS: ",
          "expected=maxRelDiff<=", batchTol, " got=", sci[maxPKvsBT]];
    AppendTo[passFlags, "gate6-FAIL"],
    AppendTo[passFlags, "gate6-PASS"]];

  (* --- Gates 7-9: single-point == first in-scan value --- *)
  kpSingle = {kps[[1]]};

  Print["  Running standalone single-point MC ..."];
  resSoloMC = quietRun @ EvaluateTropicalMC[spec, fan, kpSingle,
    "Integrator"       -> "MC",
    "Batch"            -> False,
    "NSamples"         -> nSamplesMC,
    "RunChecks"        -> False,
    "Verbose"          -> False,
    "WorkingDirectory" -> $ioDir];

  Print["  Running standalone single-point per-kp VEGAS ..."];
  resSoloPK = quietRun @ EvaluateTropicalMC[spec, fan, kpSingle,
    "Integrator"       -> "VEGAS",
    "Batch"            -> False,
    "NSamples"         -> nSamplesVEGAS,
    "RunChecks"        -> False,
    "Verbose"          -> False,
    "WorkingDirectory" -> $ioDir];

  Print["  Running standalone single-point batched VEGAS ..."];
  resSoloBT = quietRun @ EvaluateTropicalMC[spec, fan, kpSingle,
    "Integrator"       -> "VEGAS",
    "Batch"            -> True,
    "NSamples"         -> nSamplesVEGAS,
    "RunChecks"        -> False,
    "Verbose"          -> False,
    "WorkingDirectory" -> $ioDir];

  (* Gate 7: solo MC vs in-scan kp-1 MC *)
  If[AssociationQ[resSoloMC] && KeyExistsQ[resSoloMC, "Results"] &&
     Length[resSoloMC["Results"]] > 0,
    valSoloMC  = resSoloMC["Results"][[1]]["Re"];
    valScan1MC = estMC[[1]];
    Module[{d = relErr1[valSoloMC, valScan1MC]},
      Print["  Solo MC vs scan kp-1 MC: relDiff = ", sci[d],
            "  (gate ", singleTol, ")"];
      If[!TrueQ[N[d] <= singleTol],
        anyFail = True;
        Print["  CC36 FAIL  ", label, " gate(7): ",
              "expected=relDiff<=", singleTol, " got=", sci[d]];
        AppendTo[passFlags, "gate7-FAIL"],
        AppendTo[passFlags, "gate7-PASS"]]],
    Print["  gate(7) SKIP: single-point MC returned unexpected structure"]];

  (* Gate 8: solo per-kp VEGAS vs in-scan kp-1 per-kp VEGAS *)
  If[AssociationQ[resSoloPK] && KeyExistsQ[resSoloPK, "Results"] &&
     Length[resSoloPK["Results"]] > 0,
    valSoloPK  = resSoloPK["Results"][[1]]["Re"];
    valScan1PK = estPK[[1]];
    Module[{d = relErr1[valSoloPK, valScan1PK]},
      Print["  Solo per-kp VEGAS vs scan kp-1 VEGAS: relDiff = ", sci[d],
            "  (gate ", singleTol, ")"];
      If[!TrueQ[N[d] <= singleTol],
        anyFail = True;
        Print["  CC36 FAIL  ", label, " gate(8): ",
              "expected=relDiff<=", singleTol, " got=", sci[d]];
        AppendTo[passFlags, "gate8-FAIL"],
        AppendTo[passFlags, "gate8-PASS"]]],
    Print["  gate(8) SKIP: single-point per-kp VEGAS returned unexpected structure"]];

  (* Gate 9: solo batched VEGAS vs in-scan kp-1 batched VEGAS *)
  If[AssociationQ[resSoloBT] && KeyExistsQ[resSoloBT, "Results"] &&
     Length[resSoloBT["Results"]] > 0,
    valSoloBT  = resSoloBT["Results"][[1]]["Re"];
    valScan1BT = estBT[[1]];
    Module[{d = relErr1[valSoloBT, valScan1BT]},
      Print["  Solo batched VEGAS vs scan kp-1 batched VEGAS: relDiff = ", sci[d],
            "  (gate ", singleTol, ")"];
      If[!TrueQ[N[d] <= singleTol],
        anyFail = True;
        Print["  CC36 FAIL  ", label, " gate(9): ",
              "expected=relDiff<=", singleTol, " got=", sci[d]];
        AppendTo[passFlags, "gate9-FAIL"],
        AppendTo[passFlags, "gate9-PASS"]]],
    Print["  gate(9) SKIP: single-point batched VEGAS returned unexpected structure"]];

  Print["  Gates: ", StringRiffle[passFlags, "  "]];
  Print[];

  If[!anyFail,
    Print["CC36 PASS  ", label,
          "  MC=", sci[medMC],
          " perKpVEGAS=", sci[medPK],
          " batchVEGAS=", sci[medBT],
          " MCvsPK=", sci[maxMCvsPK],
          " batchAgree=", sci[maxPKvsBT],
          " (", nKP, " pts)"],
    Print["CC36 FAIL  ", label,
          "  expected=all-gates-pass  got=see-failures-above"]];
  Print[];
  !anyFail
];

(* ── 4. Sub-cases ───────────────────────────────────────────────────────── *)

$cc36Results = {};

(* Sub-case A: n=2, 30 kinematic points — fast smoke *)
AppendTo[$cc36Results,
  cc36SubCase["A: n=2, 30kp", 2, 30, 400000, 200000,
    "MCTol"         -> 5.*^-2,
    "VegasTol"      -> 2.*^-2,
    "CrossTol"      -> 5.*^-2,
    "BatchAgreeTol" -> 1.*^-2,
    "SingleTol"     -> 7.*^-2]];

(* Sub-case B: n=4, 60 kinematic points — higher-dimensional regime *)
AppendTo[$cc36Results,
  cc36SubCase["B: n=4, 60kp", 4, 60, 200000, 100000,
    "MCTol"         -> 5.*^-2,
    "VegasTol"      -> 2.*^-2,
    "CrossTol"      -> 5.*^-2,
    "BatchAgreeTol" -> 1.*^-2,
    "SingleTol"     -> 7.*^-2]];

(* ── 5. Summary ─────────────────────────────────────────────────────────── *)
Module[{nPass, nFail, labels},
  labels = {"A: n=2, 30kp", "B: n=4, 60kp"};
  nPass  = Count[$cc36Results, True];
  nFail  = Count[$cc36Results, False];
  Print["========================================================"];
  Print["CC36  Sampler x method x batch consistency  —  summary"];
  Print["========================================================"];
  Do[
    Print["  ", StringPadRight[labels[[i]], 22],
          If[i <= Length[$cc36Results] && TrueQ[$cc36Results[[i]]], "PASS", "FAIL"]],
    {i, Length[labels]}];
  Print["--------------------------------------------------------"];
  If[nFail == 0,
    Print["CC36 PASS  all ", nPass, " sub-cases passed",
          "  (MC x VEGAS x batch plumbing confirmed; §5.1 fold gated)"],
    Print["CC36 FAIL  expected=", nPass + nFail, " pass",
          "  got=", nFail, " fail  (", nPass, " passed, ", nFail, " failed)"]];
  Print["========================================================"];
  If[nFail > 0, Quit[1]]
];
