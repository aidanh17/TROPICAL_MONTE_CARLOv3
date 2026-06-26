(* ============================================================================
   TEST/phase2_selfgate.wl  --  Phase 2 SELF-GATE checks (non-runtime parts)

   (1) NEW-vocabulary codegen: re-emit every kind under the canonical
       "MC" / "VEGAS" Integrator strings (the Phase-2 vocabulary) and assert it
       equals the Phase-1 golden produced with the old "MonteCarlo"/"Vegas"
       strings -- i.e. the vocabulary migration changed NO emitted byte and the
       single emitIntegrandDefinitions + emitMain reproduce all 6 kinds.
   (2) Exactness guard (#43, plan.md §6.7) on the touched decomposition paths:
       FreeQ[symbolic ProcessSector / FlattenSector output, _Real] for an
       exact-input spec, before the MmaToC boundary.
   (3) FlattenSector eps-threading sanity: convergent vs divergent classification
       matches the inline pre-refactor semantics on a known divergent spec.

   Exits the WL kernel with a nonzero code (via Exit[1]) on any FAIL so the
   orchestrator can gate on it.
   ============================================================================ *)

SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
v3Root = Directory[];
Get[FileNameJoin[{v3Root, "tropical_eval.wl"}]];

workDir = FileNameJoin[{v3Root, "INTERFILES", "phase2_selfgate"}];
If[!DirectoryQ[workDir], CreateDirectory[workDir, CreateIntermediateDirectories -> True]];
goldenDir = FileNameJoin[{v3Root, "TEST", "baselines", "codegen_goldens"}];

eps = Symbol["eps"];
failures = {};
report[label_, ok_] := (
  Print["  ", If[TrueQ[ok], "PASS  ", "FAIL  "], label];
  If[!TrueQ[ok], AppendTo[failures, label]]);

(* diff a freshly-emitted file against its committed golden *)
diffGolden[label_, srcFile_] := Module[{a, b},
  If[!FileExistsQ[srcFile], report[label <> " (emitted)", False]; Return[]];
  a = Import[srcFile, "Text"];
  b = Import[FileNameJoin[{goldenDir, label}], "Text"];
  report[label <> " == golden (NEW vocab)", a === b]
];

Print["======================================================================="];
Print["PHASE 2 SELF-GATE -- ", DateString[]];
Print["======================================================================="];

(* ---------------------------------------------------------------------------
   (1) NEW-vocabulary codegen golden regression
   --------------------------------------------------------------------------- *)
Print["\n--- (1) NEW-vocabulary codegen (\"MC\"/\"VEGAS\") vs Phase-1 goldens ---"];

(* Spec 1: convergent, no kinematics *)
Module[{poly, vars, spec, verts, fan, sd, conv, f},
  poly = 1 + x[1]^2 + x[2]^2; vars = {x[1], x[2]};
  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {0, 0},
    "PolynomialExponents" -> {-3}, "Variables" -> vars,
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
  verts = PolytopeVertices[poly^(-3), vars];
  fan = ComputeDecomposition[verts, "ShowProgress" -> False];
  sd = Table[ProcessSector[spec, fan[[1]], fan[[2, s]], s], {s, Length[fan[[2]]]}];
  conv = Select[sd, AssociationQ[#] && !#["IsDivergent"] &];

  f = FileNameJoin[{workDir, "conv_nokat_mc.cpp"}];
  GenerateCppMonteCarlo[conv, {}, spec, f, "Integrator" -> "MC"];
  diffGolden["conv_nokat_mc.cpp", f];

  f = FileNameJoin[{workDir, "conv_nokat_vegas_perkp.cpp"}];
  GenerateCppMonteCarlo[conv, {}, spec, f, "Integrator" -> "VEGAS", "Batch" -> False];
  diffGolden["conv_nokat_vegas_perkp.cpp", f];

  f = FileNameJoin[{workDir, "conv_nokat_vegas_batch.cpp"}];
  GenerateCppMonteCarlo[conv, {}, spec, f, "Integrator" -> "VEGAS", "Batch" -> True];
  diffGolden["conv_nokat_vegas_batch.cpp", f];
];

(* Spec 2: subtraction divergent 1D (g0/g1/rem) *)
Module[{poly, vars, spec, verts, fan, sd, conv, divSD, procDiv, fMC, fV, fB},
  poly = 1 + x[1]; vars = {x[1]};
  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {-1 + eps},
    "PolynomialExponents" -> {-2}, "Variables" -> vars,
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> eps|>;
  verts = PolytopeVertices[poly^(-1), vars];
  fan = ComputeDecomposition[verts, "ShowProgress" -> False];
  sd = Table[ProcessSector[spec, fan[[1]], fan[[2, s]], s], {s, Length[fan[[2]]]}];
  conv   = Select[sd, AssociationQ[#] && !#["IsDivergent"] &];
  divSD  = Select[sd, AssociationQ[#] && #["IsDivergent"] &];
  procDiv = Select[Quiet[ProcessDivergentSector[#, spec] & /@ divSD], AssociationQ];

  fMC = FileNameJoin[{workDir, "subdiv_mc.cpp"}];
  GenerateCppMonteCarlo[conv /. eps -> 0, procDiv, spec, fMC, "Integrator" -> "MC"];
  diffGolden["subdiv_mc.cpp", fMC];

  fV = FileNameJoin[{workDir, "subdiv_vegas_perkp.cpp"}];
  GenerateCppMonteCarlo[conv /. eps -> 0, procDiv, spec, fV, "Integrator" -> "VEGAS", "Batch" -> False];
  diffGolden["subdiv_vegas_perkp.cpp", fV];

  fB = FileNameJoin[{workDir, "subdiv_vegas_batch.cpp"}];
  GenerateCppMonteCarlo[conv /. eps -> 0, procDiv, spec, fB, "Integrator" -> "VEGAS", "Batch" -> True];
  diffGolden["subdiv_vegas_batch.cpp", fB];
];

(* Spec 3: IBP divergent 1D (ibp-boundary, ibp-term) *)
Module[{poly, vars, spec, verts, fan, sd, conv, divSD, ibpSD, fMC, fV, fB},
  poly = 1 + x[1]; vars = {x[1]};
  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {-1 + eps},
    "PolynomialExponents" -> {-2}, "Variables" -> vars,
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> eps|>;
  verts = PolytopeVertices[poly^(-1), vars];
  fan = ComputeDecomposition[verts, "ShowProgress" -> False];
  sd = Table[ProcessSector[spec, fan[[1]], fan[[2, s]], s], {s, Length[fan[[2]]]}];
  conv  = Select[sd, AssociationQ[#] && !#["IsDivergent"] &];
  divSD = Select[sd, AssociationQ[#] && #["IsDivergent"] &];
  ibpSD = Select[Quiet[IBPProcessSector[#, spec] & /@ divSD], AssociationQ];

  fMC = FileNameJoin[{workDir, "ibpdiv_mc.cpp"}];
  GenerateCppMonteCarloIBP[conv /. eps -> 0, ibpSD, spec, fMC, "Integrator" -> "MC"];
  diffGolden["ibpdiv_mc.cpp", fMC];

  fV = FileNameJoin[{workDir, "ibpdiv_vegas_perkp.cpp"}];
  GenerateCppMonteCarloIBP[conv /. eps -> 0, ibpSD, spec, fV, "Integrator" -> "VEGAS", "Batch" -> False];
  diffGolden["ibpdiv_vegas_perkp.cpp", fV];

  fB = FileNameJoin[{workDir, "ibpdiv_vegas_batch.cpp"}];
  GenerateCppMonteCarloIBP[conv /. eps -> 0, ibpSD, spec, fB, "Integrator" -> "VEGAS", "Batch" -> True];
  diffGolden["ibpdiv_vegas_batch.cpp", fB];
];

(* ---------------------------------------------------------------------------
   (2) Exactness guard #43 (plan.md §6.7) on the touched decomposition paths.
   For an exact-input spec, the symbolic ProcessSector output (exponents,
   prefactor, flattened polynomials) must be free of inexact reals before the
   MmaToC boundary.  eps and kinematics stay symbolic, so we keep them symbolic.
   --------------------------------------------------------------------------- *)
Print["\n--- (2) Exactness guard #43: FreeQ[ProcessSector output, _Real] ---"];

Module[{poly, vars, spec, verts, fan, sd, conv, keysToCheck, payload},
  poly = 1 + 2 x[1]^2 + x[2]^2 + x[1] x[2]^2 + 3 x[1]^2 x[2];
  vars = {x[1], x[2]};
  (* exact rational + symbolic-eps polynomial exponent; nothing inexact in *)
  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {0, 0},
    "PolynomialExponents" -> {-3 + eps}, "Variables" -> vars,
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> eps|>;
  verts = PolytopeVertices[poly^(-3), vars];
  fan = ComputeDecomposition[verts, "ShowProgress" -> False];
  sd = Table[ProcessSector[spec, fan[[1]], fan[[2, s]], s], {s, Length[fan[[2]]]}];
  conv = Select[sd, AssociationQ[#] && !#["IsDivergent"] &];

  (* the codegen-bound payload of every convergent sector *)
  payload = Table[{c["FlattenedPolys"], c["Prefactor"], c["NewExponents"],
                   c["PolynomialExponents"], c["MinExponents"]}, {c, conv}];
  report["ProcessSector convergent payload FreeQ _Real",
         FreeQ[payload, _Real] && Length[conv] > 0];

  (* FlattenSector direct: exact in -> exact out *)
  Module[{cleared, eff, fl},
    cleared = conv[[1]]["ClearedPolys"];
    eff     = conv[[1]]["NewExponents"];
    (* re-flatten via the public FlattenSector and check exactness *)
    fl = FlattenSector[cleared, eff, Abs[conv[[1]]["DetM"]], eps];
    report["FlattenSector output FreeQ _Real",
           FreeQ[{fl["FlattenedPolys"], fl["Prefactor"]}, _Real]]
  ];
];

(* #43 corpus, complex-B lifted+divergent (planCXLIFTDIV.md §5): an exact spec
   (B = -2 + 3/10 i, eps symbolic) must carry EXACT per-IBP-piece MonoFactorLog
   and Re(atilde)-derived fields all the way to the MmaToC boundary — the
   SplitRealImag codegen numericizes only inside mmaToCInternal. *)
Module[{poly, vars, eps2, spec, lift, ls, ld, specRe, imB, fan, verts,
        sd, divSecs, ibp, payload},
  eps2 = Symbol["epsCx43"];
  poly = 1 + x[1] x[2] + 10^6 x[1]^2 + x[2]^2;
  vars = {x[1], x[2]};
  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {-1 + eps2, 0},
    "PolynomialExponents" -> {-2 + 3/10 I}, "Variables" -> vars,
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> eps2|>;
  lift = LiftCoefficients[spec,
    {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 3|>}];
  If[!AssociationQ[lift],
    report["complex-B lift+div: lift built", False],
    ls  = lift["LiftedSpec"];  ld = lift["LiftData"];
    imB = Im[ls["PolynomialExponents"]];
    specRe = MapAt[Re, ls, {Key["PolynomialExponents"]}];  (* decompose on Re(B) *)
    verts = PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]];
    fan = computeFanScaled[verts];
    sd = Table[
      Module[{s0},
        s0 = ProcessSectorLifted[specRe, fan[[1]], fan[[2, s]], s, ld,
               "Eps" -> eps2, "Verbose" -> False];
        If[AssociationQ[s0] && !KeyExistsQ[s0, "EmptyDomain"],
          s0["ImagPolyExponents"] = imB]; s0],
      {s, Length[fan[[2]]]}];
    divSecs = Select[sd, AssociationQ[#] && TrueQ[#["IsDivergent"]] &];
    report["complex-B lift+div: a divergent sector exists",
           Length[divSecs] > 0];
    If[Length[divSecs] > 0,
      ibp = IBPProcessSector[divSecs[[1]], specRe];
      report["complex-B lift+div: IBPProcessSector succeeds",
             AssociationQ[ibp]];
      If[AssociationQ[ibp],
        (* the codegen-bound phase + magnitude payload: per-piece MonoFactorLog
           (boundary + terms), Re(atilde)-derived B0/a0, the un-flattened
           LiftedMonoFactor, and the IBP prefactors — all EXACT (no inexact real)
           on an exact spec.  Im(B) (3/10) is exact too. *)
        payload = {
          divSecs[[1]]["LiftedMonoFactor"],
          ibp["BoundaryData"]["MonoFactorLog"],
          #["MonoFactorLog"] & /@ ibp["IBPTerms"],
          (* the brought-down complex IBP coefficient (Re(B)+i Im(B)) must be
             exact too — Im(B)=3/10 is an exact rational, not an inexact real. *)
          #["Coeff0"] & /@ ibp["IBPTerms"],
          #["Coeff1"] & /@ ibp["IBPTerms"],
          ibp["B0"], ibp["a0"], ibp["BoundaryData"]["Avals"]};
        report["complex-B lift+div MonoFactorLog/Re(atilde)/coeff payload FreeQ _Real",
               FreeQ[payload, _Real]]
      ]
    ]
  ]
];

(* ---------------------------------------------------------------------------
   (3) FlattenSector eps-threading classification sanity.
   1+x1 with A=-1+eps, B=-2: the sector containing the x1->0 boundary is
   divergent at eps->0 (Re[a_eff|eps=0] <= 0); FlattenSector must flag it.
   --------------------------------------------------------------------------- *)
Print["\n--- (3) FlattenSector eps-threading classification ---"];
Module[{poly, vars, spec, verts, fan, sd, nDiv, nConv},
  poly = 1 + x[1]; vars = {x[1]};
  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {-1 + eps},
    "PolynomialExponents" -> {-2}, "Variables" -> vars,
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> eps|>;
  verts = PolytopeVertices[poly^(-1), vars];
  fan = ComputeDecomposition[verts, "ShowProgress" -> False];
  sd = Table[ProcessSector[spec, fan[[1]], fan[[2, s]], s], {s, Length[fan[[2]]]}];
  nDiv  = Count[sd, _?(AssociationQ[#] && #["IsDivergent"] &)];
  nConv = Count[sd, _?(AssociationQ[#] && !#["IsDivergent"] &)];
  Print["    sectors: ", Length[sd], "  divergent: ", nDiv, "  convergent: ", nConv];
  (* the 1D 1+x1 fan has exactly one divergent sector (the x1->0 cone) *)
  report["1+x1 has >=1 divergent sector at eps->0", nDiv >= 1];
  (* divergent branch preserves ClearedPolys/NewExponents (KEEP list) *)
  report["divergent sector keeps ClearedPolys + NewExponents",
    AllTrue[Select[sd, AssociationQ[#] && #["IsDivergent"] &],
            (KeyExistsQ[#, "ClearedPolys"] && KeyExistsQ[#, "NewExponents"]) &]];
];

Print["\n======================================================================="];
If[failures === {},
  Print["PHASE 2 SELF-GATE: ALL PASS"],
  Print["PHASE 2 SELF-GATE: FAIL  -> ", failures]];
Print["======================================================================="];
If[failures =!= {}, Exit[1]];
