(* ============================================================================
   TEST/phase1_codegen_goldens.wl    --  Phase 1.2 codegen golden capture
   Captures byte-level C++ golden files for every integrand kind x sampler.
   Kinds: conv, g0, g1, rem, ibp-boundary, ibp-term
   Samplers: MC, per-kp VEGAS, batched VEGAS

   Saves to TEST/baselines/codegen_goldens/
   Uses INTERFILES/phase1_codegen/ as working dir.
   ============================================================================ *)

SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
v3Root = Directory[];
Get[FileNameJoin[{v3Root, "tropical_eval.wl"}]];

workDir = FileNameJoin[{v3Root, "INTERFILES", "phase1_codegen"}];
If[!DirectoryQ[workDir], CreateDirectory[workDir, CreateIntermediateDirectories -> True]];

goldenDir = FileNameJoin[{v3Root, "TEST", "baselines", "codegen_goldens"}];
If[!DirectoryQ[goldenDir], CreateDirectory[goldenDir, CreateIntermediateDirectories -> True]];

Print["======================================================================="];
Print["PHASE 1 CODEGEN GOLDENS -- started: ", DateString[]];
Print["v3 root: ", v3Root];
Print["work dir: ", workDir];
Print["golden dir: ", goldenDir];
Print["======================================================================="];

eps = Symbol["eps"];

(* -----------------------------------------------------------------------
   Helper: save and report a generated file
   ----------------------------------------------------------------------- *)
saveGolden[label_, srcFile_] := Module[{content},
  If[FileExistsQ[srcFile],
    content = Import[srcFile, "Text"];
    Export[FileNameJoin[{goldenDir, label}], content, "Text"];
    Print["  SAVED ", label, "  (", StringLength[content], " chars)"];
    {label, StringLength[content]},
    Print["  MISSING source: ", srcFile, " -- SKIPPED"];
    {label, 0}
  ]
];

captures = {};

(* =======================================================================
   SPEC 1: Convergent, no kinematics  (Pi/8)
   P = 1+x1^2+x2^2, B=-3
   ======================================================================= *)
Print["\n--- Spec 1: conv, no kinematics (Pi/8) ---"];
Module[{poly, vars, spec, verts, fan, sd, conv, f},
  poly = 1 + x[1]^2 + x[2]^2; vars = {x[1], x[2]};
  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {0, 0},
    "PolynomialExponents" -> {-3}, "Variables" -> vars,
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
  verts = PolytopeVertices[poly^(-3), vars];
  fan = ComputeDecomposition[verts, "ShowProgress" -> False];
  sd = Table[ProcessSector[spec, fan[[1]], fan[[2, s]], s, "Verbose" -> False],
    {s, Length[fan[[2]]]}];
  conv = Select[sd, AssociationQ[#] && !#["IsDivergent"] &];

  (* MC *)
  f = FileNameJoin[{workDir, "conv_nokat_mc.cpp"}];
  GenerateCppMonteCarlo[conv, {}, spec, f, "Integrator" -> "MonteCarlo"];
  AppendTo[captures, saveGolden["conv_nokat_mc.cpp", f]];

  (* Per-kp VEGAS *)
  f = FileNameJoin[{workDir, "conv_nokat_vegas_perkp.cpp"}];
  GenerateCppMonteCarlo[conv, {}, spec, f, "Integrator" -> "Vegas", "Batch" -> False];
  AppendTo[captures, saveGolden["conv_nokat_vegas_perkp.cpp", f]];

  (* Batched VEGAS *)
  f = FileNameJoin[{workDir, "conv_nokat_vegas_batch.cpp"}];
  GenerateCppMonteCarlo[conv, {}, spec, f, "Integrator" -> "Vegas", "Batch" -> True];
  AppendTo[captures, saveGolden["conv_nokat_vegas_batch.cpp", f]];
];

(* =======================================================================
   SPEC 2: Divergent 1D  (for g0, g1, rem)
   P = 1+x1, A = {-1+eps}, B = -2
   ======================================================================= *)
Print["\n--- Spec 2: divergent 1D (g0, g1, rem) ---"];
Module[{poly, vars, spec, verts, fan, sd, conv, divSD, procDiv, ibpSD,
        fMC, fVegas, fBatch},
  poly = 1 + x[1]; vars = {x[1]};
  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {-1 + eps},
    "PolynomialExponents" -> {-2}, "Variables" -> vars,
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> eps|>;
  verts = PolytopeVertices[(poly)^(-1), vars];
  fan = ComputeDecomposition[verts, "ShowProgress" -> False];
  sd = Table[ProcessSector[spec, fan[[1]], fan[[2, s]], s, "Verbose" -> False],
    {s, Length[fan[[2]]]}];
  conv   = Select[sd, AssociationQ[#] && !#["IsDivergent"] &];
  divSD  = Select[sd, AssociationQ[#] && #["IsDivergent"] &];
  procDiv = Quiet[ProcessDivergentSector[#, spec] & /@ divSD];
  procDiv = Select[procDiv, AssociationQ];

  (* MC: includes conv + g0 + g1 + rem *)
  fMC = FileNameJoin[{workDir, "subdiv_mc.cpp"}];
  (* eps->0 on convergent sectors for subtraction method *)
  GenerateCppMonteCarlo[conv /. eps -> 0, procDiv, spec, fMC,
                        "Integrator" -> "MonteCarlo"];
  AppendTo[captures, saveGolden["subdiv_mc.cpp", fMC]];

  (* Per-kp VEGAS *)
  fVegas = FileNameJoin[{workDir, "subdiv_vegas_perkp.cpp"}];
  GenerateCppMonteCarlo[conv /. eps -> 0, procDiv, spec, fVegas,
                        "Integrator" -> "Vegas", "Batch" -> False];
  AppendTo[captures, saveGolden["subdiv_vegas_perkp.cpp", fVegas]];

  (* Batched VEGAS (falls back to per-kp on divergent path, but still exercises the option) *)
  fBatch = FileNameJoin[{workDir, "subdiv_vegas_batch.cpp"}];
  GenerateCppMonteCarlo[conv /. eps -> 0, procDiv, spec, fBatch,
                        "Integrator" -> "Vegas", "Batch" -> True];
  AppendTo[captures, saveGolden["subdiv_vegas_batch.cpp", fBatch]];
];

(* =======================================================================
   SPEC 3: IBP divergent 1D  (ibp-boundary, ibp-term)
   Same spec as above, but IBP path
   ======================================================================= *)
Print["\n--- Spec 3: IBP divergent 1D (ibp-boundary, ibp-term) ---"];
Module[{poly, vars, spec, verts, fan, sd, conv, divSD, ibpSD,
        fMC, fVegas, fBatch},
  poly = 1 + x[1]; vars = {x[1]};
  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {-1 + eps},
    "PolynomialExponents" -> {-2}, "Variables" -> vars,
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> eps|>;
  verts = PolytopeVertices[(poly)^(-1), vars];
  fan = ComputeDecomposition[verts, "ShowProgress" -> False];
  sd = Table[ProcessSector[spec, fan[[1]], fan[[2, s]], s, "Verbose" -> False],
    {s, Length[fan[[2]]]}];
  conv  = Select[sd, AssociationQ[#] && !#["IsDivergent"] &];
  divSD = Select[sd, AssociationQ[#] && #["IsDivergent"] &];
  ibpSD = Quiet[IBPProcessSector[#, spec] & /@ divSD];
  ibpSD = Select[ibpSD, AssociationQ];

  (* MC *)
  fMC = FileNameJoin[{workDir, "ibpdiv_mc.cpp"}];
  GenerateCppMonteCarloIBP[conv /. eps -> 0, ibpSD, spec, fMC,
                            "Integrator" -> "MonteCarlo"];
  AppendTo[captures, saveGolden["ibpdiv_mc.cpp", fMC]];

  (* Per-kp VEGAS *)
  fVegas = FileNameJoin[{workDir, "ibpdiv_vegas_perkp.cpp"}];
  GenerateCppMonteCarloIBP[conv /. eps -> 0, ibpSD, spec, fVegas,
                            "Integrator" -> "Vegas", "Batch" -> False];
  AppendTo[captures, saveGolden["ibpdiv_vegas_perkp.cpp", fVegas]];

  (* Batched VEGAS *)
  fBatch = FileNameJoin[{workDir, "ibpdiv_vegas_batch.cpp"}];
  GenerateCppMonteCarloIBP[conv /. eps -> 0, ibpSD, spec, fBatch,
                            "Integrator" -> "Vegas", "Batch" -> True];
  AppendTo[captures, saveGolden["ibpdiv_vegas_batch.cpp", fBatch]];
];

(* Summary *)
Print[""];
Print["======================================================================="];
Print["CODEGEN GOLDENS SUMMARY"];
Print["======================================================================="];
Print["Total files attempted: ", Length[captures]];
nOK = Length[Select[captures, #[[2]] > 0 &]];
nFail = Length[captures] - nOK;
Do[Print["  ", c[[1]], "  ", If[c[[2]] > 0, ToString[c[[2]]] <> " chars", "MISSING"]],
  {c, captures}];
Print[""];
Print["OK: ", nOK, "  MISSING: ", nFail];
Print["Golden dir: ", goldenDir];
Print[""];
Print["======================================================================="];
Print["PHASE 1 CODEGEN GOLDENS COMPLETE -- ended: ", DateString[]];
Print["======================================================================="];
