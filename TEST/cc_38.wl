(* ============================================================================
   TEST/cc_38.wl  —  v3 Cross-Check #38
   plan.md §8.3 row #38:  Variable-relabeling / regulator-position invariance.

   Tier 1: WL + g++ only (no CUBA, no FIESTA).

   Motivation (plan.md §8.1 / §8.3 row 38):
     Index-ordering bugs in the merged codebase could silently produce sector
     sums that depend on which variable carries the regulator.  This is a
     non-trivial regression path because the two OLD trees used different
     internal conventions.  The check generalises Tree A Test 10, which
     verified divergence detection for the specific relabelling
       x[1] <-> x[2]  with regulator on x[2].

   Strategy:
     Consider the 2D integral
       I(eps) = Int_0^inf Int_0^inf  x1^{2eps-1} (1+x1+x2)^{-2} dx1 dx2
     with regulator carried by x1 (standard orientation).  The exact result is
       I(eps) = Gamma(2 eps) Gamma(1-2 eps) / Gamma(1) * 1/(2-1)
     i.e. a simple pole at eps=0 with residue 1/2.

     Now form three equivalent specifications:
       A  standard: MonomialExponents = {2eps-1, 0}, Variables = {x[1],x[2]}
       B  swapped variables: MonomialExponents = {0, 2eps-1}, Variables = {x[2],x[1]}
          (same integrand, just the coordinate labels x[1]<->x[2] interchanged)
       C  regulator on x[2]: MonomialExponents = {0, 2eps-1}, Variables = {x[1],x[2]}
          (i.e. x2 plays the role of the regulated variable)

     In each case:
       1. Compute the tropical fan via PolytopeVertices + ComputeDecomposition.
       2. Run ProcessSector over all sectors; collect divergent sectors.
       3. Check that the number of divergent sectors is >= 1.
       4. Run NIntegrate-based ValidateSubtraction on the first divergent sector
          at eps = 0.1 (well inside the convergence radius).
       5. Record relErr.

     PASS sub-criteria (all must hold for CC38 PASS):
       P1  divCount(A) >= 1  (at least one divergent sector found)
       P2  divCount(B) >= 1  (same after relabelling)
       P3  divCount(C) >= 1  (same with regulator on the other variable)
       P4  relErr(A) < 0.05  (subtraction matches NIntegrate)
       P5  relErr(B) < 0.05  (relabelled version also passes)
       P6  relErr(C) < 0.05  (regulator-moved version also passes)
       P7  divCount(A) == divCount(B) == divCount(C)
           (identical divergent-sector count is the strongest cross-check)

   Ported / adapted from:
     OLD_CODE/.../EXAMPLES/tropical_eval_examples.wl  Test 10
     (Divergent variable index permutation, lines 1780-1858)

   Run:
     wolframscript -file TEST/cc_38.wl
   from the TROPICAL_MONTE_CARLOv3 root.  No external deps required.
   ============================================================================ *)

(* --------------------------------------------------------------------------
   Load v3 packages by absolute path so this file runs from any working dir.
   -------------------------------------------------------------------------- *)
$v3Root = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{$v3Root, "tropical_fan.wl"}]];
Get[FileNameJoin[{$v3Root, "tropical_eval.wl"}]];

Print[""];
Print["================================================================"];
Print["  v3 Cross-Check #38: Variable-relabeling / regulator-position"];
Print["                      invariance (plan.md §8.3 row #38)"];
Print["================================================================"];
Print[""];

(* --------------------------------------------------------------------------
   Helpers
   -------------------------------------------------------------------------- *)

(* K-scaled fan robust helper — mirrors gen_sectors.wl computeFanRobust (§6.4) *)
computeFanRobust38[verts_] := Module[{n, fd},
  n = Length[First[verts]];
  Do[fd = Quiet @ ComputeDecomposition[K * verts, "ShowProgress" -> False];
     If[ListQ[fd] && Length[fd] == 2 && FreeQ[fd, $Failed] && Length[fd[[1]]] > 0,
       Return[fd, Module]],
     {K, {1, n + 2, 2 n + 4, 6 n + 6}}];
  $Failed];

(* Run the full ProcessSector + ProcessDivergentSector + ValidateSubtraction
   pipeline for one spec.  Returns an Association with keys:
     "DivCount"  -> number of divergent sectors found
     "RelErr"    -> RelativeError from ValidateSubtraction on first div sector
     "Pass"      -> True/False overall sub-pass
     "Note"      -> human-readable string  *)
runVariant38[label_String, spec_, poly_, vars_, testEps_] :=
  Module[
    {verts, fanData, dualVerts, simplexList, allSectors, divSectors,
     sd, divData, vsResult, relErr, divCount, note, pass},

    verts = PolytopeVertices[poly^(-2), vars];
    fanData = computeFanRobust38[verts];

    If[fanData === $Failed,
      note = label <> ": fan $Failed";
      Print["  [", label, "] Fan computation failed"];
      Return[<|"DivCount" -> 0, "RelErr" -> Missing[], "Pass" -> False,
               "Note" -> note|>]
    ];

    dualVerts   = fanData[[1]];
    simplexList = fanData[[2]];

    Print["  [", label, "] Fan: ", Length[dualVerts], " rays, ",
          Length[simplexList], " sectors"];

    allSectors = Table[
      Quiet @ ProcessSector[spec, dualVerts, simplexList[[i]], i],
      {i, Length[simplexList]}
    ];

    divSectors = Select[allSectors, AssociationQ[#] && TrueQ[#["IsDivergent"]] &];
    divCount   = Length[divSectors];

    Print["  [", label, "] Divergent sectors: ", divCount];

    If[divCount == 0,
      note = label <> ": divCount=0";
      Print["  [", label, "] No divergent sectors found — FAIL"];
      Return[<|"DivCount" -> 0, "RelErr" -> Missing[], "Pass" -> False,
               "Note" -> note|>]
    ];

    (* Print which variables carry the divergence *)
    Print["  [", label, "] DivergentVariables: ",
          Union[#["DivergentVariable"] & /@ divSectors]];

    (* ValidateSubtraction on the first divergent sector at testEps *)
    sd      = First[divSectors];
    divData = Quiet @ ProcessDivergentSector[sd, spec];

    If[!AssociationQ[divData],
      note = label <> ": ProcessDivergentSector failed";
      Print["  [", label, "] ProcessDivergentSector failed — FAIL"];
      Return[<|"DivCount" -> divCount, "RelErr" -> Missing[], "Pass" -> False,
               "Note" -> note|>]
    ];

    vsResult = Quiet @ ValidateSubtraction[divData, sd, spec, {}, testEps];

    If[!AssociationQ[vsResult],
      note = label <> ": ValidateSubtraction failed";
      Print["  [", label, "] ValidateSubtraction failed — FAIL"];
      Return[<|"DivCount" -> divCount, "RelErr" -> Missing[], "Pass" -> False,
               "Note" -> note|>]
    ];

    relErr = Lookup[vsResult, "RelativeError", Missing[]];
    Print["  [", label, "] ValidateSubtraction relErr = ", relErr];

    pass = NumericQ[relErr] && relErr < 0.05;
    note = label <> ": divCount=" <> ToString[divCount] <>
           " relErr=" <> ToString[relErr];

    Print["  [", label, "] ", If[pass, "PASS", "FAIL"]];
    <|"DivCount" -> divCount, "RelErr" -> relErr, "Pass" -> pass,
      "Note" -> note|>
  ];


(* ============================================================================
   The shared polynomial.  Both variables enter symmetrically:
     (1 + x1 + x2)  is symmetric under x1 <-> x2.
   Exact:  Int_0^inf x_reg^{2eps-1} (1+x1+x2)^{-2} dx1 dx2
           = Gamma(2eps) Gamma(1-2eps)  [has 1/eps pole, residue = 1/2]
   ============================================================================ *)

eps38 = Symbol["eps38"];
poly38 = 1 + x[1] + x[2];
testEps38 = 0.1;    (* evaluate at eps = 0.1 for subtraction check *)

(* --------------------------------------------------------------------------
   Variant A — regulator on x[1], Variables = {x[1], x[2]}  (standard)
   -------------------------------------------------------------------------- *)
Print["--- Variant A: regulator on x[1], vars = {x[1], x[2]} ---"];
specA = <|
  "Polynomials"        -> {poly38},
  "MonomialExponents"  -> {2 eps38 - 1, 0},
  "PolynomialExponents"-> {-2},
  "Variables"          -> {x[1], x[2]},
  "KinematicSymbols"   -> {},
  "RegulatorSymbol"    -> eps38
|>;
rA = runVariant38["A", specA, poly38, {x[1], x[2]}, testEps38];
Print[];

(* --------------------------------------------------------------------------
   Variant B — swap variable labels: x[1] <-> x[2].
   Polynomial is symmetric, so the integrand is identical.
   MonomialExponents are permuted to keep the regulator on the *same*
   integration variable (physical x1) but now called x[2] in the spec.
   This exercises the index-relabelling path directly (cf. Test 10 in OLD_CODE).
   -------------------------------------------------------------------------- *)
Print["--- Variant B: vars = {x[2], x[1]}  (relabelled), regulator on x[2] ---"];
(* After relabelling: x[1]->x[2], x[2]->x[1], poly = 1+x[2]+x[1] = poly38.
   The regulated variable is now x[2] in the spec list. *)
specB = <|
  "Polynomials"        -> {1 + x[2] + x[1]},
  "MonomialExponents"  -> {0, 2 eps38 - 1},   (* second slot = x[2] = regulated *)
  "PolynomialExponents"-> {-2},
  "Variables"          -> {x[1], x[2]},        (* x[2] carries the regulator now *)
  "KinematicSymbols"   -> {},
  "RegulatorSymbol"    -> eps38
|>;
(* Fan must be computed over the same polynomial (symmetric) *)
rB = runVariant38["B", specB, 1 + x[2] + x[1], {x[1], x[2]}, testEps38];
Print[];

(* --------------------------------------------------------------------------
   Variant C — same as A but regulator explicitly on x[2], vars canonical.
   MonomialExponents = {0, 2eps-1}.  Same integrand by symmetry of poly38.
   -------------------------------------------------------------------------- *)
Print["--- Variant C: MonomialExponents={0,2eps-1}, regulator on x[2] ---"];
specC = <|
  "Polynomials"        -> {poly38},
  "MonomialExponents"  -> {0, 2 eps38 - 1},
  "PolynomialExponents"-> {-2},
  "Variables"          -> {x[1], x[2]},
  "KinematicSymbols"   -> {},
  "RegulatorSymbol"    -> eps38
|>;
rC = runVariant38["C", specC, poly38, {x[1], x[2]}, testEps38];
Print[];


(* ============================================================================
   Evaluate all PASS criteria
   ============================================================================ *)
Print["================================================================"];
Print["  Cross-check #38 summary"];
Print["================================================================"];

passP1 = rA["DivCount"] >= 1;
passP2 = rB["DivCount"] >= 1;
passP3 = rC["DivCount"] >= 1;
passP4 = TrueQ[rA["Pass"]];
passP5 = TrueQ[rB["Pass"]];
passP6 = TrueQ[rC["Pass"]];
passP7 = (rA["DivCount"] == rB["DivCount"] == rC["DivCount"]);

printCriterion[label_String, p_] :=
  Print["  ", label, ": ", If[p, "PASS", "FAIL"]];

printCriterion["P1 divCount(A)>=1", passP1];
printCriterion["P2 divCount(B)>=1", passP2];
printCriterion["P3 divCount(C)>=1", passP3];
printCriterion["P4 relErr(A)<0.05", passP4];
printCriterion["P5 relErr(B)<0.05", passP5];
printCriterion["P6 relErr(C)<0.05", passP6];
printCriterion[
  "P7 divCount(A)==divCount(B)==divCount(C)  [" <>
  ToString[rA["DivCount"]] <> "," <>
  ToString[rB["DivCount"]] <> "," <>
  ToString[rC["DivCount"]] <> "]",
  passP7];
Print[];

allPass = And[passP1, passP2, passP3, passP4, passP5, passP6, passP7];

If[allPass,
  Print["CC38 PASS  divCount=", rA["DivCount"],
        "  relErr(A)=", rA["RelErr"],
        "  relErr(B)=", rB["RelErr"],
        "  relErr(C)=", rC["RelErr"]],
  (* Collect failures *)
  With[{fails = Select[
          {If[!passP1, "P1:divCount(A)<1",          Nothing],
           If[!passP2, "P2:divCount(B)<1",          Nothing],
           If[!passP3, "P3:divCount(C)<1",          Nothing],
           If[!passP4, "P4:relErr(A)>=0.05 got="<>ToString[rA["RelErr"]], Nothing],
           If[!passP5, "P5:relErr(B)>=0.05 got="<>ToString[rB["RelErr"]], Nothing],
           If[!passP6, "P6:relErr(C)>=0.05 got="<>ToString[rC["RelErr"]], Nothing],
           If[!passP7,
              "P7:divCounts differ A="<>ToString[rA["DivCount"]]<>
              " B="<>ToString[rB["DivCount"]]<>
              " C="<>ToString[rC["DivCount"]],
              Nothing]},
          StringQ]},
    Print["CC38 FAIL  expected=all-criteria-pass  got=",
          StringRiffle[fails, "; "]]
  ];
  Quit[1]
];
