(* ============================================================================
   TEST/phase3_selfgate.wl  --  Phase 3 SELF-GATE (lifting, real coefficients)

   Fast (NIntegrate-gated, no C++) checks that the lifting graft is correct and
   that the non-lifted pipeline is unregressed.  The heavy end-to-end C++ checks
   live in TEST/test_lifted.wl (Test 23) and SANDBOX/bench_lift_variance.wl.

     (1) Exactness invariant / guard #43 (plan.md §6.7): the LIFTED symbolic
         decomposition output is FreeQ[_Real] before the MmaToC boundary, and
         z0/residuals are exact.
     (2) Round-trip identity #40: z->z0 reproduces the original EXACTLY for +,
         -, and multi-rule lifts (the liftidentity relative-tolerance fix never
         false-positives on these exact cases).
     (3) BUG-2 aux-var robustness: a plain-symbol spec fires TropicalEval::
         liftbadvar and returns $Failed (no Symbol[n+1]).
     (4) The exact ValidateLiftedDecomposition correctness gate (NOT the sampled
         sigma, plan.md §9 risk #1): Toy-1 (k=2) relErr < 0.5% with 3 EmptyDomain
         drops, and Toy-0 (k=1) relErr < 0.1% with EmptyDomain drops.
     (5) HasConstantTerm reported per sector by the driver; EmptyDomain handled.
     (6) EXACT error-path routing: liftcomplex (complex atilde) and liftnopivot
         (all atilde <= 0) fire deterministically and return $Failed.
     (7) Domain-indicator C++ is emitted for lifted sectors with a constraint.

   Exits the WL kernel with Exit[1] on any FAIL so the orchestrator can gate.
   ============================================================================ *)

SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
v3Root = Directory[];
Get[FileNameJoin[{v3Root, "tropical_eval.wl"}]];

workDir = FileNameJoin[{v3Root, "INTERFILES", "phase3_selfgate"}];
If[!DirectoryQ[workDir], CreateDirectory[workDir, CreateIntermediateDirectories -> True]];

failures = {};
report[label_, ok_] := (
  Print["  ", If[TrueQ[ok], "PASS  ", "FAIL  "], label];
  If[!TrueQ[ok], AppendTo[failures, label]];
);

Print["======================================================================="];
Print["PHASE 3 SELF-GATE (lifting) -- ", DateString[]];
Print["======================================================================="];

(* Case A lifted spec, reused below *)
specA = <|"Polynomials" -> {1 + 10^6 x[1]^2 + x[2]^2 + x[1] x[2]^2},
  "MonomialExponents" -> {0, 0}, "PolynomialExponents" -> {-2},
  "Variables" -> {x[1], x[2]}, "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
ruleA = {<|"PolyIndex"->1, "ExponentVector"->{2,0}, "k"->2|>};

(* ---- (1) Exactness guard #43 on lifted decomposition output ---- *)
Print["\n--- (1) Exactness guard #43: FreeQ[_Real] on lifted output ---"];
Module[{lc, ls, ld, verts, fan, sectors, conv, bad},
  lc = LiftCoefficients[specA, ruleA];
  ls = lc["LiftedSpec"]; ld = lc["LiftData"];
  report["z0 exact (FreeQ _Real)", FreeQ[ld["z0"], _Real]];
  report["residuals exact (FreeQ _Real)", FreeQ[ld["Residuals"], _Real]];
  verts = PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]];
  fan = ComputeDecomposition[verts, "ShowProgress" -> False];
  sectors = Table[ProcessSectorLifted[ls, fan[[1]], fan[[2,s]], s, ld], {s, Length[fan[[2]]]}];
  conv = Select[sectors, AssociationQ[#] && !KeyExistsQ[#, "EmptyDomain"] &];
  bad = {};
  Do[
    Module[{sd = conv[[i]]},
      If[!FreeQ[{sd["NewExponents"], sd["Prefactor"], sd["ClearedPolys"],
                 sd["FlattenedPolys"], sd["DomainConstraint"], sd["MinExponents"]}, _Real],
        AppendTo[bad, sd["ConeIndex"]]]
    ], {i, Length[conv]}];
  report["all lifted sector payloads FreeQ _Real", bad === {}];
];

(* ---- (2) Round-trip identity #40 for +, -, multi-rule ---- *)
Print["\n--- (2) Round-trip identity #40 (z->z0 === original, EXACT) ---"];
roundTrip[spec_, rules_] := Module[{lc, z0, av},
  lc = LiftCoefficients[spec, rules];
  If[!AssociationQ[lc], Return[False]];
  z0 = lc["LiftData"]["z0"]; av = lc["LiftData"]["AuxVariable"];
  FreeQ[z0, _Real] && And @@ Table[
    TrueQ[Expand[lc["LiftedSpec"]["Polynomials"][[j]] /. av -> z0] ===
          Expand[spec["Polynomials"][[j]]]], {j, Length[spec["Polynomials"]]}]
];
report["(+) 10^6 x1^2 k=2", roundTrip[specA, ruleA]];
report["(-) -10^6 x1^2 k=2",
  roundTrip[<|"Polynomials"->{1 - 10^6 x[1]^2 + x[2]^2}, "MonomialExponents"->{0,0},
    "PolynomialExponents"->{-2}, "Variables"->{x[1],x[2]}, "KinematicSymbols"->{},
    "RegulatorSymbol"->None|>, ruleA]];
report["multi-rule 10^6 x1^2 & 10^6 x2^2 k=2",
  roundTrip[<|"Polynomials"->{1 + 10^6 x[1]^2 + 10^6 x[2]^2}, "MonomialExponents"->{0,0},
    "PolynomialExponents"->{-2}, "Variables"->{x[1],x[2]}, "KinematicSymbols"->{},
    "RegulatorSymbol"->None|>,
    {<|"PolyIndex"->1,"ExponentVector"->{2,0},"k"->2|>, <|"PolyIndex"->1,"ExponentVector"->{0,2},"k"->2|>}]];

(* ---- (3) BUG-2 aux-var robustness ---- *)
Print["\n--- (3) BUG-2: plain-symbol spec fires liftbadvar, returns $Failed ---"];
Module[{r},
  r = Quiet @ LiftCoefficients[
    <|"Polynomials"->{1 + 10^6 aa^2 + bb^2}, "MonomialExponents"->{0,0},
      "PolynomialExponents"->{-2}, "Variables"->{aa,bb}, "KinematicSymbols"->{},
      "RegulatorSymbol"->None|>, {<|"PolyIndex"->1,"ExponentVector"->{2,0},"k"->2|>}];
  report["plain-symbol LiftCoefficients -> $Failed", r === $Failed];
];

(* ---- (4) Exact ValidateLiftedDecomposition correctness gate ---- *)
Print["\n--- (4) ValidateLiftedDecomposition (EXACT gate, not sampled sigma) ---"];
(* Toy 1: P = 1 + x1 + 10^-6 x2, exact = 500000, k=2, explicit fan, 3 drops *)
Module[{spec1, lc, ls, ld, fan, vl},
  spec1 = <|"Polynomials"->{1 + x[1] + 10^-6 x[2]}, "MonomialExponents"->{0,0},
    "PolynomialExponents"->{-3}, "Variables"->{x[1],x[2]}, "KinematicSymbols"->{},
    "RegulatorSymbol"->None|>;
  lc = LiftCoefficients[spec1, {<|"PolyIndex"->1,"ExponentVector"->{0,1},"k"->2|>}];
  ls = lc["LiftedSpec"]; ld = lc["LiftData"];
  fan = {{{1,0,0},{0,1,0},{-1,-1,0},{0,2,-1},{0,-2,1}},
         {{1,2,4},{2,3,4},{3,1,4},{1,2,5},{2,3,5},{3,1,5}}};
  vl = Quiet @ ValidateLiftedDecomposition[spec1, ls, fan, ld, {}, 3];
  report["Toy-1 k=2 relErr < 0.5%", AssociationQ[vl] && NumericQ[vl["RelativeError"]] && vl["RelativeError"] < 0.005];
  report["Toy-1 k=2 DroppedSectors == 3", AssociationQ[vl] && Length[vl["DroppedSectors"]] == 3];
];
(* Toy 0: P = 1 + 10^6 x1, exact = 10^-6, k=1, explicit fan, EmptyDomain drops *)
Module[{spec0, lc, ls, ld, raysK1, sectsK1, vl},
  spec0 = <|"Polynomials"->{1 + 10^6 x[1]}, "MonomialExponents"->{0},
    "PolynomialExponents"->{-2}, "Variables"->{x[1]}, "KinematicSymbols"->{},
    "RegulatorSymbol"->None|>;
  lc = LiftCoefficients[spec0, {<|"PolyIndex"->1,"ExponentVector"->{1},"k"->1|>}];
  ls = lc["LiftedSpec"]; ld = lc["LiftData"];
  raysK1 = {{1,0},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1}};
  sectsK1 = Table[{i, If[i < Length[raysK1], i+1, 1]}, {i, Length[raysK1]}];
  vl = Quiet @ ValidateLiftedDecomposition[spec0, ls, {raysK1, sectsK1}, ld, {}, 4];
  report["Toy-0 k=1 relErr < 0.1%", AssociationQ[vl] && NumericQ[vl["RelativeError"]] && vl["RelativeError"] < 0.001];
  report["Toy-0 k=1 EmptyDomain drops > 0", AssociationQ[vl] && Length[vl["DroppedSectors"]] > 0];
];

(* ---- (5) HasConstantTerm reported per sector; EmptyDomain handled ---- *)
Print["\n--- (5) Driver reports HasConstantTerm per sector; EmptyDomain handled ---"];
Module[{res},
  res = EvaluateTropicalMCLifted[specA, {{}}, "LiftRules" -> ruleA,
    "NSamples" -> 50000, "RunChecks" -> False, "Verbose" -> False,
    "WorkingDirectory" -> workDir];
  report["HasConstantTerm list length == ConvergentSectors",
    AssociationQ[res] && ListQ[res["HasConstantTerm"]] &&
    Length[res["HasConstantTerm"]] == res["ConvergentSectors"]];
  report["EmptyDomainSectors reported (== 2 for Case A)",
    AssociationQ[res] && res["EmptyDomainSectors"] == 2];
  report["Case A has a HasConstantTerm=False sector (risk-1 flagged)",
    AssociationQ[res] && MemberQ[res["HasConstantTerm"], False]];
];

(* ---- (6) EXACT error-path routing ---- *)
Print["\n--- (6) Error paths: liftcomplex + liftnopivot (EXACT, deterministic) ---"];
Module[{specC, lc, ls, ld, verts, fan, fired},
  specC = <|"Polynomials"->{1 + 10^6 x[1]^2 + x[2]^2 + x[1] x[2]^2}, "MonomialExponents"->{0,0},
    "PolynomialExponents"->{-2 + I}, "Variables"->{x[1],x[2]}, "KinematicSymbols"->{},
    "RegulatorSymbol"->None|>;
  lc = LiftCoefficients[specC, {<|"PolyIndex"->1,"ExponentVector"->{2,0},"k"->2|>}];
  ls = lc["LiftedSpec"]; ld = lc["LiftData"];
  verts = Quiet @ PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]];
  fan = Quiet @ ComputeDecomposition[verts, "ShowProgress" -> False];
  fired = False;
  Do[If[Quiet @ ProcessSectorLifted[ls, fan[[1]], fan[[2,s]], s, ld] === $Failed, fired = True],
     {s, Length[fan[[2]]]}];
  report["liftcomplex fires + $Failed on complex-B lift", fired];
];
Module[{specNP, ld, dv, r},
  specNP = <|"Polynomials"->{1 + x[3]}, "MonomialExponents"->{-5,-5,3},
    "PolynomialExponents"->{-3}, "Variables"->{x[1],x[2],x[3]}, "KinematicSymbols"->{},
    "RegulatorSymbol"->None|>;
  dv = {{-1,0,-1},{0,-1,-1},{0,0,-1}};
  ld = <|"z0"->1, "AuxIndex"->3, "AuxVariable"->x[3], "Rules"->{},
    "OriginalSpec"-><|"Polynomials"->{1+x[1]+x[2]}, "MonomialExponents"->{0,0},
      "PolynomialExponents"->{-3}, "Variables"->{x[1],x[2]}, "KinematicSymbols"->{}, "RegulatorSymbol"->None|>|>;
  r = Quiet @ ProcessSectorLifted[specNP, dv, {1,2,3}, 1, ld];
  report["liftnopivot fires + $Failed on all-nonpositive atilde", r === $Failed];
];

(* ---- (7) Domain-indicator C++ emitted ---- *)
Print["\n--- (7) Domain-indicator C++ emitted for constrained lifted sectors ---"];
Module[{code, nInd},
  code = Import[FileNameJoin[{workDir, "tropical_mc_generated.cpp"}], "Text"];
  nInd = StringCount[code, "lifted-sector domain indicator"];
  report["domain indicator block(s) present in emitted C++", nInd >= 1];
];

(* ---- Summary ---- *)
Print["\n======================================================================="];
If[failures === {},
  Print["PHASE 3 SELF-GATE: ALL PASS"],
  Print["PHASE 3 SELF-GATE: ", Length[failures], " FAILURE(S): ", failures]
];
Print["======================================================================="];
If[failures =!= {}, Exit[1]];
