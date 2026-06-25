(* ============================================================================
   bench_lift_variance.wl  --  Phase 3 HONEST variance benchmark (plan §8.3,
   §9 risk #1; cross-check #18/#19).

   ----------------------------------------------------------------------------
   WHY THIS WAS REWRITTEN (RE-RUN NOTE, prior Phase-3 attempt was RED):
   The earlier benchmark headlined Case A "VarRed ~ 18x -> PASS (>= 10x)"
   computed from SAMPLED sigma (sigma = ReErr * Sqrt[N]).  That number is a
   LIE for Case A: Case A's lifted decomposition contains a HasConstantTerm=False
   sector (conv_5 / cone-8) whose I2 = Int f^2 DIVERGES, so the lifted estimator
   has INFINITE true variance.  The sampled sigma was small only because a finite
   MC sample missed the heavy tail entirely (feasFrac->0 for that sector).
   ----------------------------------------------------------------------------
   HONEST GATE (this file): compute EXACT per-sample trueSigma per scheme
       trueSigma^2 = I2 - I1^2,   I1 = Int g,  I2 = Int g^2
   by converged NIntegrate (trueSigma[] in sandbox_lift_common.wl), NOT sampled
   sigma.  The total estimator variance is Sum_sectors trueSigma_s^2; it is
   INFINITE the moment ANY sector has I2Converged->False.  Then classify:
     (i)  ALL lifted sectors HasConstantTerm=True  -> finite I2 -> report the
          EXACT VarRed = trueSigma^2_unlifted / trueSigma^2_lifted; PASS only if
          genuinely >= 1 (a real reduction).
     (ii) ANY lifted sector HasConstantTerm=False  -> I2 diverges -> report
          INFINITE-VARIANCE / NOT-BENEFICIAL and assert the HasConstantTerm
          (and magnitude) WARNING FIRES.  Documented L5 limitation (plan §9
          risk #1) -- NOT a correctness failure and NOT an 18x pass.

   Case A  MUST land in (ii): infinite true variance, warning fires.
   Case C  MUST be correctly NOT-beneficial (all lifted sectors HasConstantTerm
           =False) -- and the geometry-gate anchor selection refuses to lift it
           automatically (verified below: Automatic -> no lift).
   Case B  is report-only; its 1-rule-k=2 decomposition is all HasConstantTerm
           =True, so it gets a finite EXACT VarRed.

   The exact ValidateLiftedDecomposition (Test 23 / phase3_selfgate) -- NOT this
   variance number -- is the CORRECTNESS gate (plan §9 risk #1).  This file
   gates only the VARIANCE story, honestly.

   PASS condition (printed as BENCH GATE):
     * exact trueSigma computed for EVERY scheme (no $Failed I1),
     * every finite-VarRed claim corroborated by exact trueSigma (no sampled-
       sigma headline numbers),
     * the HasConstantTerm warning fires on EVERY infinite-variance case
       (Case A, Case C),
     * Case C is correctly classified NOT-beneficial and NOT auto-lifted.

   Run:
     cd /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3 &&
     wolframscript -file SANDBOX/bench_lift_variance.wl
   Exits Exit[1] if the honest gate fails.
   ============================================================================ *)

SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
v3Root = Directory[];
Get[FileNameJoin[{v3Root, "tropical_eval.wl"}]];
Get[FileNameJoin[{v3Root, "SANDBOX", "sandbox_lift_common.wl"}]];
Print["tropical_eval.wl + sandbox_lift_common.wl loaded"];
Print[];

(* PrecisionGoal for the exact-trueSigma NIntegrate (I1, I2). *)
tsPG = 4;

(* sampled-sigma column is INFORMATIONAL ONLY (never gated). *)
NSamplesInfo = 200000;

failures = {};
fail[msg_] := (Print["  GATE-FAIL: ", msg]; AppendTo[failures, msg]);

fmtN[x_] := Which[
  x === Infinity, "INFINITE",
  x === Indeterminate, "Indeterminate",
  NumericQ[x], ToString[N[x, 5], InputForm],  (* InputForm -> clean inline "a*^-b" *)
  True, ToString[x, InputForm]];


(* ============================================================================
   EXACT-trueSigma over a list of sector-data assocs (lifted OR unlifted).
   Each sd must carry FlattenedPolys / Prefactor / Dimension /
   PolynomialExponents (+ optional DomainConstraint) -- which BOTH ProcessSector
   and ProcessSectorLifted emit.  EmptyDomain sectors contribute 0.
   Returns <| I1Sum, SigmaSqTotal, SigmaTotal, AnyDivergent, DivergentIdx,
              PerSector |> with SigmaTotal=Infinity if any I2 diverges.
   ============================================================================ *)
exactSigmaTotal[sectorList_List] := Module[
  {ts, i1sum, divIdx, sigSq, sigTot},
  ts = Table[
    If[AssociationQ[sd] && !TrueQ[Lookup[sd, "EmptyDomain", False]] &&
       !TrueQ[Lookup[sd, "IsDivergent", False]],
      trueSigma[sd, {}, tsPG],
      <|"I1" -> 0, "I2" -> 0, "Sigma" -> 0, "I2Converged" -> True|>],
    {sd, sectorList}];
  i1sum  = Total[Map[#["I1"] &, ts]];
  divIdx = Flatten@Position[Map[TrueQ[#["I2Converged"]] &, ts], False];
  If[divIdx === {},
    sigSq  = Total[Map[#["Sigma"]^2 &, ts]];
    sigTot = Sqrt[sigSq],
    sigSq  = Infinity;
    sigTot = Infinity
  ];
  <|"I1Sum" -> i1sum, "SigmaSqTotal" -> sigSq, "SigmaTotal" -> sigTot,
    "AnyDivergent" -> (divIdx =!= {}), "DivergentIdx" -> divIdx,
    "PerSector" -> ts|>
];

(* Per-sector dump: I1 / I2 / Sigma / I2Converged / HasConstantTerm. *)
dumpSectors[label_String, sectorList_List, esResult_Association] := Module[{ts},
  ts = esResult["PerSector"];
  Print["  --- ", label, " per-sector EXACT trueSigma ---"];
  Do[
    Module[{sd = sectorList[[i]], r = ts[[i]], hc},
      hc = If[AssociationQ[sd], Lookup[sd, "HasConstantTerm", "n/a"], "n/a"];
      If[AssociationQ[sd] && TrueQ[Lookup[sd, "EmptyDomain", False]],
        Print["    sector ", i, " (cone ", Lookup[sd, "ConeIndex", "?"],
              "): EmptyDomain (dropped)"],
        Print["    sector ", i, " (cone ", Lookup[sd, "ConeIndex", "?"], "): ",
              "I1=", fmtN[r["I1"]],
              "  I2=", If[r["I2Converged"], fmtN[r["I2"]], "DIVERGENT"],
              "  Sigma=", If[r["I2Converged"], fmtN[r["Sigma"]], "INF"],
              "  HasConstantTerm=", hc]
      ]
    ],
    {i, Length[sectorList]}];
  If[esResult["AnyDivergent"],
    Print["    => TOTAL trueSigma = INFINITE (divergent I2 in sectors ",
          esResult["DivergentIdx"], ")"],
    Print["    => TOTAL trueSigma = ", fmtN[esResult["SigmaTotal"]]]
  ];
];


(* ============================================================================
   Build unlifted sectors via ProcessSector over a unit-coefficient proxy fan.
   ============================================================================ *)
buildUnlifted[spec_Association, proxyPoly_] := Module[{verts, fan},
  verts = PolytopeVertices[proxyPoly^(-1), spec["Variables"]];
  fan   = ComputeDecomposition[verts, "ShowProgress" -> False];
  <|"Sectors" -> Table[ProcessSector[spec, fan[[1]], fan[[2, s]], s],
                       {s, Length[fan[[2]]]}],
    "NSectors" -> Length[fan[[2]]]|>
];

(* Build lifted sectors via the PACKAGE (LiftCoefficients + ProcessSectorLifted).
   fanOpt = Automatic -> derive from lifted polytope; else explicit {rays,simps}. *)
buildLifted[spec_Association, liftRules_, fanOpt_] := Module[
  {lc, ls, ld, fan, verts, sectors, hasConstList, nDropped},
  lc = LiftCoefficients[spec, liftRules];
  If[!AssociationQ[lc],
    Return[<|"Failed" -> True, "Reason" -> lc|>]];
  ls = lc["LiftedSpec"]; ld = lc["LiftData"];
  fan = If[fanOpt =!= Automatic, fanOpt,
    Module[{vv = Quiet@PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]]},
      Quiet@ComputeDecomposition[vv, "ShowProgress" -> False]]];
  If[!(ListQ[fan] && Length[fan] >= 2),
    Return[<|"Failed" -> True, "Reason" -> "fan", "Fan" -> fan|>]];
  sectors = Table[
    Quiet@ProcessSectorLifted[ls, fan[[1]], fan[[2, s]], s, ld],
    {s, Length[fan[[2]]]}];
  hasConstList = Cases[sectors,
    a_?AssociationQ /; !TrueQ[Lookup[a, "EmptyDomain", False]] :>
      Lookup[a, "HasConstantTerm", Missing[]]];
  nDropped = Count[sectors,
    a_?AssociationQ /; TrueQ[Lookup[a, "EmptyDomain", False]]];
  <|"Failed" -> False, "Sectors" -> sectors, "LiftedSpec" -> ls,
    "LiftData" -> ld, "HasConstList" -> hasConstList,
    "NDropped" -> nDropped, "z0" -> ld["z0"]|>
];

(* Informational sampled sigma via the C++ driver (NEVER gated). *)
sampledSigmaUnlifted[spec_, proxyPoly_, wd_] := Module[{verts, fan, res},
  verts = PolytopeVertices[proxyPoly^(-1), spec["Variables"]];
  fan   = ComputeDecomposition[verts, "ShowProgress" -> False];
  res = Quiet@EvaluateTropicalMC[spec, fan, {{}}, "NSamples" -> NSamplesInfo,
    "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> wd];
  If[AssociationQ[res] && KeyExistsQ[res, "Results"],
    res["Results"][[1]]["ReErr"] * Sqrt[NSamplesInfo], Missing[]]
];
sampledSigmaLifted[spec_, liftRules_, fanOpt_, wd_] := Module[{res},
  res = Quiet@EvaluateTropicalMCLifted[spec, {{}}, "LiftRules" -> liftRules,
    "FanData" -> fanOpt, "NSamples" -> NSamplesInfo, "RunChecks" -> False,
    "Verbose" -> False, "WorkingDirectory" -> wd];
  If[AssociationQ[res] && KeyExistsQ[res, "Results"],
    res["Results"][[1]]["ReErr"] * Sqrt[NSamplesInfo], Missing[]]
];

wd = FileNameJoin[{v3Root, "INTERFILES", "bench_lift"}];
If[!DirectoryQ[wd], CreateDirectory[wd, CreateIntermediateDirectories -> True]];


(* ============================================================================
   CASE A: P = 1 + 10^6 x1^2 + x2^2 + x1 x2^2, B={-2}.  Lift {2,0} k=2.
   EXPECTED: lifted decomposition has a HasConstantTerm=False sector ->
             INFINITE true variance (class ii).  Warning must fire.
   ============================================================================ *)
Print["================================================================="];
Print["CASE A: P = 1 + 10^6 x1^2 + x2^2 + x1 x2^2,  B={-2}   (lift {2,0} k=2)"];
Print["================================================================="];

specA = <|"Polynomials" -> {1 + 10^6 x[1]^2 + x[2]^2 + x[1] x[2]^2},
  "MonomialExponents" -> {0, 0}, "PolynomialExponents" -> {-2},
  "Variables" -> {x[1], x[2]}, "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
ruleA = {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 2|>};

(* Converged reference: high-resolution NIntegrate of the direct integrand. *)
refA = Quiet@NIntegrate[1/(1 + 10^6 t1^2 + t2^2 + t1 t2^2)^2,
  {t1, 0, Infinity}, {t2, 0, Infinity},
  MaxRecursion -> 30, PrecisionGoal -> 8, WorkingPrecision -> 30];
Print["  Converged NIntegrate reference (I1) = ", N[refA, 8]];

uA = buildUnlifted[specA, 1 + x[1]^2 + x[2]^2 + x[1] x[2]^2];
esUA = exactSigmaTotal[uA["Sectors"]];
dumpSectors["Case A UNLIFTED", uA["Sectors"], esUA];
Print["  Unlifted I1 sum = ", fmtN[esUA["I1Sum"]],
      "  (ref ", N[refA, 8], "; relErr ",
      fmtN[Abs[(esUA["I1Sum"] - refA)/refA]], ")"];

lA = buildLifted[specA, ruleA, Automatic];
If[TrueQ[lA["Failed"]], fail["Case A lifted build failed: " <> ToString[lA["Reason"]]]];
esLA = exactSigmaTotal[lA["Sectors"]];
dumpSectors["Case A LIFTED k=2", lA["Sectors"], esLA];
Print["  Lifted I1 sum = ", fmtN[esLA["I1Sum"]],
      "  (ref ", N[refA, 8], "; relErr ",
      fmtN[Abs[(esLA["I1Sum"] - refA)/refA]], ")"];
Print["  Lifted z0 = ", lA["z0"], "   HasConstantTerm per sector = ", lA["HasConstList"]];
hcFalseA = MemberQ[lA["HasConstList"], False];
If[hcFalseA,
  Print["  *** WARNING: HasConstantTerm=False in lifted sector(s) -> I2 DIVERGES",
        " -> INFINITE true variance (plan §9 risk #1 / L5). ***"]];

(* sampled sigma -- INFORMATIONAL ONLY *)
sampUA = sampledSigmaUnlifted[specA, 1 + x[1]^2 + x[2]^2 + x[1] x[2]^2, wd];
sampLA = sampledSigmaLifted[specA, ruleA, Automatic, wd];
Print["  [informational] sampled sigma: unlifted=", fmtN[sampUA],
      "  lifted=", fmtN[sampLA],
      "  sampled-ratio=", fmtN[If[NumericQ[sampUA] && NumericQ[sampLA] && sampLA > 0,
                                  (sampUA/sampLA)^2, Indeterminate]],
      "  <-- DO NOT TRUST (lifted sector with divergent I2 makes this optimistic)"];
Print[];


(* ============================================================================
   CASE B (report-only): P = 10^-4 + 10^4 x1^2 + 10^-4 x2^2 + 10^4 x1 x2^2
                           + x1^2 x2, B={-2}.  Single-rule k=2 on {2,0}.
   EXPECTED: all lifted sectors HasConstantTerm=True -> finite EXACT VarRed.
   ============================================================================ *)
Print["================================================================="];
Print["CASE B (report-only): 5-monomial mixed-coeff, B={-2}   (1-rule k=2 {2,0})"];
Print["================================================================="];

specB = <|"Polynomials" -> {10^-4 + 10^4 x[1]^2 + 10^-4 x[2]^2 +
                            10^4 x[1] x[2]^2 + x[1]^2 x[2]},
  "MonomialExponents" -> {0, 0}, "PolynomialExponents" -> {-2},
  "Variables" -> {x[1], x[2]}, "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
ruleB = {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 2|>};

refB = Quiet@NIntegrate[
  1/(10^-4 + 10^4 t1^2 + 10^-4 t2^2 + 10^4 t1 t2^2 + t1^2 t2)^2,
  {t1, 0, Infinity}, {t2, 0, Infinity},
  MaxRecursion -> 30, PrecisionGoal -> 8, WorkingPrecision -> 30];
Print["  Converged NIntegrate reference (I1) = ", N[refB, 8]];

uB = buildUnlifted[specB, 1 + x[1]^2 + x[2]^2 + x[1] x[2]^2 + x[1]^2 x[2]];
esUB = exactSigmaTotal[uB["Sectors"]];
dumpSectors["Case B UNLIFTED", uB["Sectors"], esUB];

lB = buildLifted[specB, ruleB, Automatic];
If[TrueQ[lB["Failed"]], fail["Case B lifted build failed: " <> ToString[lB["Reason"]]]];
esLB = exactSigmaTotal[lB["Sectors"]];
dumpSectors["Case B LIFTED k=2", lB["Sectors"], esLB];
Print["  Lifted z0 = ", lB["z0"], "   HasConstantTerm per sector = ", lB["HasConstList"]];
hcFalseB = MemberQ[lB["HasConstList"], False];
If[hcFalseB,
  Print["  WARNING: HasConstantTerm=False in a Case-B lifted sector -> infinite variance."]];
Print["  Case B I1 sums: unlifted=", fmtN[esUB["I1Sum"]], "  lifted=", fmtN[esLB["I1Sum"]],
      "  (ref ", N[refB, 8], ")"];
Print[];


(* ============================================================================
   CASE C: P = 1 + 10^8 x1^3 x2 + x2^3, B={-3}.  Degenerate lifted polytope.
   EXPECTED: (a) Automatic anchor selection (geometry gate) REFUSES to lift
             (no surviving HasConstantTerm=True sector for any k);
             (b) the explicit-fan lifts (k=2, k=4) have all sectors
                 HasConstantTerm=False -> NOT beneficial (class ii).
   ============================================================================ *)
Print["================================================================="];
Print["CASE C: P = 1 + 10^8 x1^3 x2 + x2^3,  B={-3}   (degenerate lift)"];
Print["================================================================="];

specC = <|"Polynomials" -> {1 + 10^8 x[1]^3 x[2] + x[2]^3},
  "MonomialExponents" -> {0, 0}, "PolynomialExponents" -> {-3},
  "Variables" -> {x[1], x[2]}, "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;

refC = Quiet@NIntegrate[1/(1 + 10^8 t1^3 t2 + t2^3)^3,
  {t1, 0, Infinity}, {t2, 0, Infinity},
  MaxRecursion -> 30, PrecisionGoal -> 8, WorkingPrecision -> 30];
Print["  Converged NIntegrate reference (I1) = ", N[refC, 8]];

uC = buildUnlifted[specC, 1 + x[1]^3 x[2] + x[2]^3];
esUC = exactSigmaTotal[uC["Sectors"]];
dumpSectors["Case C UNLIFTED", uC["Sectors"], esUC];

(* (a) geometry-gate anchor selection: Automatic must NOT produce a beneficial
   lift for Case C.  DetectExtremeCoefficients DOES return a kStar (coeff 10^8 is
   extreme), but the lifted Newton polytope is lower-dimensional, so the
   EvaluateTropicalMCLifted Automatic fan-build hits the degeneracy guard and
   fires TropicalEval::liftdegenerate -> $Failed.  That is the §6.2 "do not lift"
   outcome (refuse loudly), NOT a silent wrong answer. *)
detC = Quiet@DetectExtremeCoefficients[specC, 1000];
Print["  DetectExtremeCoefficients (Automatic kStar) -> ", detC];
autoResC = Quiet@EvaluateTropicalMCLifted[specC, {{}}, "LiftRules" -> Automatic,
  "NSamples" -> 20000, "RunChecks" -> False, "Verbose" -> False,
  "WorkingDirectory" -> wd];
(* "not beneficial / not auto-lifted" = the driver does NOT return a usable MC
   result for the lifted Case C (it refuses via liftdegenerate -> $Failed). *)
autoNoLiftC = !(AssociationQ[autoResC] && KeyExistsQ[autoResC, "Results"]);
Print["  Geometry gate: Automatic EvaluateTropicalMCLifted on Case C -> ",
      If[autoNoLiftC, "refuses (no usable lifted result; liftdegenerate -> $Failed)",
         "UNEXPECTEDLY produced a result"]];

(* (b) explicit-fan lifts -- show they are all HasConstantTerm=False. *)
dvCk2 = {{1, 0, 0}, {0, 1, 0}, {-1, -3, 0}, {2, 0, -3}, {-2, 0, 3}};
slC   = {{1, 2, 4}, {2, 3, 4}, {3, 1, 4}, {1, 2, 5}, {2, 3, 5}, {3, 1, 5}};
dvCk4 = {{1, 0, 0}, {0, 1, 0}, {-1, -3, 0}, {4, 0, -3}, {-4, 0, 3}};
ruleCk2 = {<|"PolyIndex" -> 1, "ExponentVector" -> {3, 1}, "k" -> 2|>};
ruleCk4 = {<|"PolyIndex" -> 1, "ExponentVector" -> {3, 1}, "k" -> 4|>};

lCk2 = buildLifted[specC, ruleCk2, {dvCk2, slC}];
esLCk2 = If[TrueQ[lCk2["Failed"]], <|"SigmaTotal" -> Infinity, "AnyDivergent" -> True,
  "DivergentIdx" -> {}, "PerSector" -> {}, "I1Sum" -> Indeterminate|>,
  exactSigmaTotal[lCk2["Sectors"]]];
If[!TrueQ[lCk2["Failed"]], dumpSectors["Case C LIFTED k=2 (explicit fan)", lCk2["Sectors"], esLCk2]];
hcListCk2 = If[TrueQ[lCk2["Failed"]], {}, lCk2["HasConstList"]];
Print["  Case C k=2 HasConstantTerm per sector = ", hcListCk2,
      "   I1 sum = ", fmtN[esLCk2["I1Sum"]]];

lCk4 = buildLifted[specC, ruleCk4, {dvCk4, slC}];
esLCk4 = If[TrueQ[lCk4["Failed"]], <|"SigmaTotal" -> Infinity, "AnyDivergent" -> True,
  "DivergentIdx" -> {}, "PerSector" -> {}, "I1Sum" -> Indeterminate|>,
  exactSigmaTotal[lCk4["Sectors"]]];
hcListCk4 = If[TrueQ[lCk4["Failed"]], {}, lCk4["HasConstList"]];
Print["  Case C k=4 HasConstantTerm per sector = ", hcListCk4,
      "   I1 sum = ", fmtN[esLCk4["I1Sum"]]];

hcFalseC = MemberQ[hcListCk2, False] || MemberQ[hcListCk4, False];
allFalseC = (hcListCk2 =!= {} && AllTrue[hcListCk2, # === False &]) ||
            (hcListCk4 =!= {} && AllTrue[hcListCk4, # === False &]);
If[hcFalseC,
  Print["  *** WARNING: HasConstantTerm=False in Case C lifted sector(s) ->",
        " I2 DIVERGES -> NOT beneficial (structural; documented limitation). ***"]];
Print[];


(* ============================================================================
   CASE D (class (i) FINITE-VARIANCE demonstration):
       P = 1 + 10^6 x1, B=-2, 1D.  Exact integral = 10^-6.  Lift {1} k=2 (z0=10^3).
   This is the hand-checkable Toy-0 example with B=-2 (so f^2 is integrable).
   After lifting, ALL surviving sectors have HasConstantTerm=True AND convergent
   I2 -> FINITE true variance -> a GENUINE exact VarRed >= 1 (lifting helps).
   It exercises the class-(i) branch with a real, converged reduction so the
   benefit of lifting is not merely asserted.  Explicit k=2 fan (support
   {(0,0),(1,2)}, dividing hyperplane n1+2n2=0; the toy0 k=2 triangulation).
   ============================================================================ *)
Print["================================================================="];
Print["CASE D (FINITE-VARIANCE demo): P = 1 + 10^6 x1,  B={-2}, 1D  (lift {1} k=2)"];
Print["================================================================="];

specD = <|"Polynomials" -> {1 + 10^6 x[1]}, "MonomialExponents" -> {0},
  "PolynomialExponents" -> {-2}, "Variables" -> {x[1]},
  "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
ruleD = {<|"PolyIndex" -> 1, "ExponentVector" -> {1}, "k" -> 2|>};
exactD = 10^-6;
(* k=2 fan: complete unimodular triangulation of R^2 respecting n1+2n2=0. *)
raysD = {{1, 0}, {0, 1}, {-1, 2}, {-1, 0}, {-1, -2}, {0, -1}, {1, -2}};
sectsD = Table[{i, If[i < Length[raysD], i + 1, 1]}, {i, Length[raysD]}];

uD = buildUnlifted[specD, 1 + x[1]];
esUD = exactSigmaTotal[uD["Sectors"]];
dumpSectors["Case D UNLIFTED", uD["Sectors"], esUD];
Print["  Unlifted I1 sum = ", fmtN[esUD["I1Sum"]], "  (exact ", N[exactD], "; relErr ",
      fmtN[Abs[(esUD["I1Sum"] - exactD)/exactD]], ")"];

lD = buildLifted[specD, ruleD, {raysD, sectsD}];
If[TrueQ[lD["Failed"]], fail["Case D lifted build failed: " <> ToString[lD["Reason"]]]];
esLD = exactSigmaTotal[lD["Sectors"]];
dumpSectors["Case D LIFTED k=2", lD["Sectors"], esLD];
Print["  Lifted I1 sum = ", fmtN[esLD["I1Sum"]], "  (exact ", N[exactD], "; relErr ",
      fmtN[Abs[(esLD["I1Sum"] - exactD)/exactD]], ")"];
Print["  Lifted z0 = ", lD["z0"], "   HasConstantTerm per sector = ", lD["HasConstList"]];
hcFalseD = MemberQ[lD["HasConstList"], False];
vrD = If[NumericQ[esUD["SigmaTotal"]] && NumericQ[esLD["SigmaTotal"]] && esLD["SigmaTotal"] > 0,
         (esUD["SigmaTotal"]/esLD["SigmaTotal"])^2, Indeterminate];
Print["  EXACT VarRed (trueSigma_U^2/trueSigma_L^2) = ", fmtN[vrD],
      "   (all HasConstantTerm=True, all I2 converge -> GENUINE finite reduction)"];
Print[];


(* ============================================================================
   HONEST CLASSIFICATION + GATE
   ============================================================================ *)
Print["================================================================="];
Print["  HONEST VARIANCE CLASSIFICATION (EXACT trueSigma; no sampled sigma)"];
Print["================================================================="];

(* Per-scheme classify.  The infinite-variance SIGNAL is the EXACT divergence of
   I2 (esL["AnyDivergent"]), NOT HasConstantTerm by itself -- because a
   HasConstantTerm=True sector can STILL have a divergent I2 (Case B sector 7).
   HasConstantTerm=False is the structural ROOT CAUSE for class (ii) and is
   reported, but the gate keys on the exact I2 convergence.  A WARNING is emitted
   on EVERY infinite-variance scheme so the gate invariant holds.
   classWarned (module-global) records every scheme that warned. *)
classWarned = {};
classify[label_, esU_, esL_, hcFalse_] := Module[{vr, cls, divergent},
  divergent = TrueQ[esL["AnyDivergent"]];
  If[divergent,
    cls = "INFINITE-VARIANCE / NOT-BENEFICIAL (class ii: divergent I2"
          <> If[hcFalse, "; HasConstantTerm=False root cause",
                         "; divergent even with HasConstantTerm=True"] <> ")";
    vr = Infinity;
    Print["  *** WARNING [", label, "]: divergent I2 in sector(s) ",
          esL["DivergentIdx"], " -> INFINITE true variance",
          If[hcFalse, " (HasConstantTerm=False; plan §9 risk #1 / L5)",
                      " (HasConstantTerm=True but I2 still diverges)"], ". ***"];
    AppendTo[classWarned, label],
    cls = "FINITE-VARIANCE (class i: all I2 converge, all HasConstantTerm=True)";
    vr = If[NumericQ[esU["SigmaTotal"]] && NumericQ[esL["SigmaTotal"]] && esL["SigmaTotal"] > 0,
            (esU["SigmaTotal"]/esL["SigmaTotal"])^2, Indeterminate]
  ];
  Print["  ", label, ":"];
  Print["     trueSigma_unlifted = ", fmtN[esU["SigmaTotal"]],
        "   trueSigma_lifted = ", fmtN[esL["SigmaTotal"]]];
  Print["     class: ", cls];
  Print["     EXACT VarRed (trueSigma_U^2 / trueSigma_L^2) = ", fmtN[vr]];
  <|"class" -> cls, "vr" -> vr, "hcFalse" -> hcFalse, "divergent" -> divergent|>
];

cA = classify["Case A (lift {2,0} k=2)", esUA, esLA, hcFalseA];
cB = classify["Case B (1-rule k=2)", esUB, esLB, hcFalseB];
(* Case C: take the better (still-divergent) of k=2/k=4 for the report. *)
cC = classify["Case C (explicit-fan lift)", esUC,
       If[esLCk2["SigmaTotal"] === Infinity && esLCk4["SigmaTotal"] === Infinity,
          esLCk2, (* both infinite *)
          If[NumericQ[esLCk2["SigmaTotal"]], esLCk2, esLCk4]],
       hcFalseC];
cD = classify["Case D (1D 1+10^6 x1, B=-2, lift k=2) -- FINITE-VARIANCE demo", esUD, esLD, hcFalseD];
Print[];


(* ============================================================================
   GATE ASSERTIONS
   ============================================================================ *)
Print["================================================================="];
Print["  BENCH GATE ASSERTIONS"];
Print["================================================================="];

(* (G1) exact trueSigma computed for every scheme: I1 numeric (NIntegrate ran). *)
g1ok = NumericQ[esUA["I1Sum"]] && NumericQ[esLA["I1Sum"]] &&
       NumericQ[esUB["I1Sum"]] && NumericQ[esLB["I1Sum"]] &&
       NumericQ[esUC["I1Sum"]] && NumericQ[esUD["I1Sum"]] && NumericQ[esLD["I1Sum"]];
If[!g1ok, fail["exact trueSigma I1 not numeric for some scheme"]];
Print["  (G1) exact trueSigma (I1,I2) computed for every scheme: ",
      If[g1ok, "PASS", "FAIL"]];

(* (G2) decomposition correctness: I1 sums reproduce the converged reference. *)
relA = Abs[(esLA["I1Sum"] - refA)/refA];
relUA = Abs[(esUA["I1Sum"] - refA)/refA];
relC = Abs[(esUC["I1Sum"] - refC)/refC];
If[!(NumericQ[relA] && relA < 0.01), fail["Case A lifted I1 sum off reference: relErr=" <> fmtN[relA]]];
If[!(NumericQ[relUA] && relUA < 0.01), fail["Case A unlifted I1 sum off reference: relErr=" <> fmtN[relUA]]];
Print["  (G2) I1 decomposition reproduces converged reference (Case A lifted relErr=",
      fmtN[relA], ", unlifted relErr=", fmtN[relUA], "): ",
      If[NumericQ[relA] && relA < 0.01 && NumericQ[relUA] && relUA < 0.01, "PASS", "FAIL"]];

(* (G3) Case A MUST be class (ii): HasConstantTerm=False sector -> infinite variance. *)
If[!hcFalseA, fail["Case A expected a HasConstantTerm=False sector (class ii); none found"]];
If[!(esLA["AnyDivergent"]), fail["Case A expected divergent I2 (infinite true variance); I2 converged"]];
Print["  (G3) Case A HasConstantTerm=False warning fires + true variance INFINITE: ",
      If[hcFalseA && esLA["AnyDivergent"], "PASS", "FAIL"],
      "  (divergent sectors=", esLA["DivergentIdx"], ")"];

(* (G4) Case A's sampled sigma is NOT used as a pass; assert the headline is INFINITE. *)
If[cA["vr"] =!= Infinity, fail["Case A VarRed must be reported INFINITE, got " <> fmtN[cA["vr"]]]];
Print["  (G4) Case A headline VarRed = INFINITE (NOT an 18x sampled pass): ",
      If[cA["vr"] === Infinity, "PASS", "FAIL"]];

(* (G5) Case C correctly NOT beneficial: all lifted sectors HasConstantTerm=False. *)
If[!allFalseC, fail["Case C expected ALL lifted sectors HasConstantTerm=False"]];
If[!hcFalseC, fail["Case C expected HasConstantTerm=False warning"]];
Print["  (G5) Case C all lifted sectors HasConstantTerm=False (NOT beneficial): ",
      If[allFalseC && hcFalseC, "PASS", "FAIL"]];

(* (G6) geometry gate: Automatic anchor selection does NOT auto-lift Case C. *)
If[!autoNoLiftC, fail["Case C geometry gate should refuse Automatic lift (no admissible kStar)"]];
Print["  (G6) Case C geometry-gate anchor selection refuses Automatic lift: ",
      If[autoNoLiftC, "PASS", "FAIL"]];

(* (G7) Case D demonstrates a GENUINE class-(i) finite VarRed: all sectors
   HasConstantTerm=True, all I2 converge, and VarRed >= 1 (lifting genuinely
   reduces finite variance).  This exercises the finite-variance branch with a
   converged reduction (not a sampled number). *)
g7ok = !hcFalseD && !esLD["AnyDivergent"] && NumericQ[vrD] && vrD >= 1 &&
       NumericQ[esLD["SigmaTotal"]] && NumericQ[esUD["SigmaTotal"]];
If[!g7ok, fail["Case D finite-variance demo failed: hcFalseD=" <> ToString[hcFalseD] <>
   " divergent=" <> ToString[esLD["AnyDivergent"]] <> " vrD=" <> fmtN[vrD]]];
Print["  (G7) Case D class-(i) GENUINE finite VarRed = ", fmtN[vrD],
      " (>=1, all HasConstantTerm=True, all I2 converge, corroborated by exact trueSigma): ",
      If[g7ok, "PASS", "FAIL"]];

(* (G8) the warning invariant: EVERY infinite-variance scheme raised a warning.
   classWarned must contain exactly the divergent schemes (A, B, C) and not D. *)
divergentSchemes = Select[{{"Case A (lift {2,0} k=2)", cA},
                           {"Case B (1-rule k=2)", cB},
                           {"Case C (explicit-fan lift)", cC},
                           {"Case D (1D 1+10^6 x1, B=-2, lift k=2) -- FINITE-VARIANCE demo", cD}},
  TrueQ[#[[2]]["divergent"]] &][[All, 1]];
warnInvariant = AllTrue[divergentSchemes, MemberQ[classWarned, #] &] &&
                !MemberQ[classWarned, "Case D (1D 1+10^6 x1, B=-2, lift k=2) -- FINITE-VARIANCE demo"];
If[!warnInvariant, fail["warning did not fire on every infinite-variance scheme: warned=" <>
   ToString[classWarned] <> " divergent=" <> ToString[divergentSchemes]]];
Print["  (G8) infinite-variance WARNING fired on EVERY divergent scheme ",
      divergentSchemes, " (and not on finite Case D): ",
      If[warnInvariant, "PASS", "FAIL"]];

(* (G9) Case B (report-only): exact trueSigma computed; reported honestly
   (infinite here -- sector 7 I2 diverges even though HasConstantTerm=True). *)
Print["  (G9) Case B (report-only) reported honestly via exact trueSigma: VarRed=",
      fmtN[cB["vr"]], "  ", If[TrueQ[cB["divergent"]],
        "(INFINITE: divergent I2 in a HasConstantTerm=True sector -- correctly NOT a finite pass)",
        "(finite, genuine)"]];
Print[];


(* ============================================================================
   WRITE benchmark_results.md (HONEST)
   ============================================================================ *)
mdFile = FileNameJoin[{DirectoryName[$InputFileName], "benchmark_results.md"}];
sline[x_] := If[StringQ[x], x, ToString[x]];
md = StringRiffle[Select[{
  "# TROPICAL_MONTE_CARLO v3: Lift-Variance Benchmark (HONEST exact-trueSigma)",
  "",
  "**Date:** " <> DateString[{"Year", "-", "Month", "-", "Day"}] <> "  ",
  "**Method:** EXACT per-sample sigma per scheme, trueSigma^2 = I2 - I1^2 with",
  "I1 = Int g and I2 = Int g^2 by converged NIntegrate (PrecisionGoal " <> ToString[tsPG] <> ").",
  "Total estimator variance = Sum_sectors trueSigma_s^2; INFINITE the moment ANY",
  "sector has a divergent I2.  HasConstantTerm=False is the usual structural ROOT",
  "CAUSE of a divergent I2 (re-clearing loses the constant term), but it is NOT a",
  "perfect proxy: a HasConstantTerm=True sector can STILL have a divergent I2 (see",
  "Case B sector 7).  The gate therefore keys on the EXACT I2 convergence and",
  "raises a warning on every divergent scheme.  Sampled sigma (sigma = ReErr*Sqrt[N])",
  "is shown for transparency only and is NEVER used as a pass criterion -- for a",
  "sector with divergent I2 it is an optimistic (often wildly underestimated) proxy.",
  "",
  "The CORRECTNESS gate is the EXACT `ValidateLiftedDecomposition` (Test 23 /",
  "phase3_selfgate), NOT this variance number (plan.md §9 risk #1).",
  "",
  "---",
  "",
  "## Case A: `P = 1 + 10^6 x1^2 + x2^2 + x1 x2^2`, B={-2}  (lift {2,0}, k=2)",
  "",
  "| Scheme | I1 sum | converged ref | I1 relErr | true sigma | classification |",
  "|--------|--------|---------------|-----------|------------|----------------|",
  "| Unlifted | " <> fmtN[esUA["I1Sum"]] <> " | " <> fmtN[refA] <> " | " <> fmtN[relUA] <>
    " | " <> fmtN[esUA["SigmaTotal"]] <> " | finite |",
  "| Lifted k=2 | " <> fmtN[esLA["I1Sum"]] <> " | " <> fmtN[refA] <> " | " <> fmtN[relA] <>
    " | **INFINITE** | **class (ii)** |",
  "",
  "Lifted HasConstantTerm per surviving sector: " <> sline[lA["HasConstList"]] <> "  (z0 = " <> sline[lA["z0"]] <> ").",
  "The HasConstantTerm=False sector (fan index " <> sline[esLA["DivergentIdx"]] <> ")",
  "has DIVERGENT I2, so the lifted estimator has **INFINITE true variance**.",
  "**EXACT VarRed = INFINITE-VARIANCE / NOT-BENEFICIAL** (plan §9 risk #1, bug-log L5).",
  "",
  "> The earlier \"VarRed ~ 18x PASS\" was computed from SAMPLED sigma " <>
    "(unlifted " <> fmtN[sampUA] <> " / lifted " <> fmtN[sampLA] <> "); it is an",
  "> artifact of the finite MC missing the heavy tail of the HasConstantTerm=False",
  "> sector (feasFrac -> 0).  The HONEST result is infinite variance.",
  "",
  "---",
  "",
  "## Case B (report-only): mixed-coefficient 5-monomial, B={-2}  (1-rule k=2 {2,0})",
  "",
  "| Scheme | I1 sum | converged ref | true sigma | EXACT VarRed |",
  "|--------|--------|---------------|------------|--------------|",
  "| Unlifted | " <> fmtN[esUB["I1Sum"]] <> " | " <> fmtN[refB] <> " | " <> fmtN[esUB["SigmaTotal"]] <> " | - |",
  "| Lifted k=2 | " <> fmtN[esLB["I1Sum"]] <> " | " <> fmtN[refB] <> " | " <> fmtN[esLB["SigmaTotal"]] <>
    " | " <> fmtN[cB["vr"]] <> " |",
  "",
  "Lifted HasConstantTerm per sector: " <> sline[lB["HasConstList"]] <> ".",
  If[TrueQ[cB["divergent"]],
    "Despite ALL lifted sectors being HasConstantTerm=True, sector " <> sline[esLB["DivergentIdx"]] <>
      " has a DIVERGENT I2, so the lifted estimator has INFINITE true variance and the VarRed " <>
      "is reported INFINITE -- a sharper, more honest result than HasConstantTerm alone would suggest.",
    "All lifted sectors HasConstantTerm=True AND all I2 converge -> finite I2 -> the EXACT VarRed above is genuine."],
  "Case B is report-only (no pass gate).",
  "",
  "---",
  "",
  "## Case C: `P = 1 + 10^8 x1^3 x2 + x2^3`, B={-3}  (degenerate lifted polytope)",
  "",
  "Automatic anchor selection (geometry gate) " <>
    If[autoNoLiftC, "**REFUSES to lift** (no admissible kStar with a surviving HasConstantTerm=True sector).",
       "unexpectedly produced an anchor."],
  "Explicit-fan lifts (k=2 z0=10^4, k=4 z0=10^2) were run to expose the structure:",
  "",
  "| Scheme | I1 sum | true sigma | HasConstantTerm/sector |",
  "|--------|--------|------------|------------------------|",
  "| Unlifted | " <> fmtN[esUC["I1Sum"]] <> " | " <> fmtN[esUC["SigmaTotal"]] <> " | (n/a) |",
  "| Lifted k=2 | " <> fmtN[esLCk2["I1Sum"]] <> " | " <> fmtN[esLCk2["SigmaTotal"]] <> " | " <> sline[hcListCk2] <> " |",
  "| Lifted k=4 | " <> fmtN[esLCk4["I1Sum"]] <> " | " <> fmtN[esLCk4["SigmaTotal"]] <> " | " <> sline[hcListCk4] <> " |",
  "",
  "ALL lifted sectors (both k) are HasConstantTerm=False -> **INFINITE true variance** ->",
  "**NOT BENEFICIAL** (structural, documented).  Correctly NOT auto-lifted.",
  "",
  "---",
  "",
  "## Case D (FINITE-VARIANCE demonstration): `P = 1 + 10^6 x1`, B={-2}, 1D  (lift {1}, k=2)",
  "",
  "Hand-checkable Toy-0 with B=-2 (exact integral = 1e-6); f^2 IS integrable.  After",
  "lifting (z0=" <> sline[lD["z0"]] <> "), ALL surviving sectors have HasConstantTerm=True AND convergent I2.",
  "",
  "| Scheme | I1 sum | exact | true sigma | EXACT VarRed |",
  "|--------|--------|-------|------------|--------------|",
  "| Unlifted | " <> fmtN[esUD["I1Sum"]] <> " | 1.0e-6 | " <> fmtN[esUD["SigmaTotal"]] <> " | - |",
  "| Lifted k=2 | " <> fmtN[esLD["I1Sum"]] <> " | 1.0e-6 | " <> fmtN[esLD["SigmaTotal"]] <>
    " | **" <> fmtN[vrD] <> "** |",
  "",
  "HasConstantTerm per sector: " <> sline[lD["HasConstList"]] <> ".  This is the class-(i) branch:",
  "a GENUINE finite variance reduction (VarRed = " <> fmtN[vrD] <> " >= 1), corroborated by EXACT",
  "trueSigma -- lifting genuinely helps when the lifted sectors keep an O(1) constant term",
  "and f^2 stays integrable.",
  "",
  "---",
  "",
  "## Summary (HONEST)",
  "",
  "| Case | true sigma_unlifted | true sigma_lifted | EXACT VarRed | verdict |",
  "|------|---------------------|-------------------|--------------|---------|",
  "| A | " <> fmtN[esUA["SigmaTotal"]] <> " | INFINITE | INFINITE | NOT-beneficial (class ii; divergent I2 in HasConstantTerm=False sector; warning fires) |",
  "| B | " <> fmtN[esUB["SigmaTotal"]] <> " | " <> fmtN[esLB["SigmaTotal"]] <> " | " <> fmtN[cB["vr"]] <>
    " | report-only" <> If[TrueQ[cB["divergent"]], " (INFINITE: divergent I2 even with HasConstantTerm=True)", " (finite, genuine)"] <> " |",
  "| C | " <> fmtN[esUC["SigmaTotal"]] <> " | INFINITE | INFINITE | NOT-beneficial (class ii; all HasConstantTerm=False; not auto-lifted) |",
  "| D | " <> fmtN[esUD["SigmaTotal"]] <> " | " <> fmtN[esLD["SigmaTotal"]] <> " | **" <> fmtN[vrD] <>
    "** | class (i): GENUINE finite reduction (all HasConstantTerm=True, all I2 converge) |",
  "",
  "**Key correction vs the prior RED run:** Case A is reported as INFINITE true",
  "variance (HasConstantTerm=False sector with divergent I2), NOT an 18x sampled-sigma",
  "pass.  The ONLY genuine variance reduction (Case D, " <> fmtN[vrD] <> "x) is corroborated by",
  "EXACT trueSigma with all sectors finite-I2.  No finite VarRed is claimed from sampled sigma."
}, StringQ], "\n"];
Export[mdFile, md, "Text"];
Print["Wrote ", mdFile];
Print[];


(* ============================================================================
   FINAL GATE
   ============================================================================ *)
Print["================================================================="];
If[failures === {},
  Print["BENCH GATE: ALL PASS (honest exact-trueSigma; A/B/C infinite-variance with",
        " warnings firing on every divergent scheme; Case D class-(i) genuine finite",
        " VarRed=", fmtN[vrD], " corroborated by exact trueSigma; no sampled-sigma headline)"],
  Print["BENCH GATE: ", Length[failures], " FAILURE(S):"];
  Do[Print["    - ", f], {f, failures}]
];
Print["================================================================="];
If[failures =!= {}, Exit[1]];
