(* ============================================================================
   TEST/cc_25.wl  —  Cross-check #25 (plan.md §8.2)
   Golden-master codegen regression: byte-identical emitted C++.

   For every integrand kind × sampler combination captured during Phase 1
   (conv/subdiv/ibp × MC/per-kp-VEGAS/batched-VEGAS = 9 golden files),
   regenerate the C++ from the v3 engine and diff against the stored golden.
   A single differing character is a FAIL.

   PASS criterion (plan.md §8.2 #25):
     conv / g0 / g1 / rem / ibp-boundary / ibp-term, MC & VEGAS — diff clean.

   Tier 1: WL + g++ only (goldens stored; g++ not actually invoked here).

   Ported from:
     OLD_CODE/TROPICAL_MONTE_CARLO/TEST/golden_capture.wl  (check mode)
     TEST/phase1_codegen_goldens.wl                         (Phase-1 golden capture)

   Run:
     wolframscript -file TEST/cc_25.wl
   from the TROPICAL_MONTE_CARLO3 directory.
   ============================================================================ *)

(* --- locate root via absolute path; load v3 packages --- *)
With[{root = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]]},
  Get[FileNameJoin[{root, "tropical_fan.wl"}]];
  Get[FileNameJoin[{root, "tropical_eval.wl"}]];
];

Print["CC25: packages loaded."];
Print[];

(* ============================================================================
   Locate golden directory (written by TEST/phase1_codegen_goldens.wl)
   ============================================================================ *)
cc25Root    = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
goldenDir   = FileNameJoin[{cc25Root, "TEST", "baselines", "codegen_goldens"}];
workDir     = FileNameJoin[{cc25Root, "TEST", "INTERFILES", "cc25"}];

If[!DirectoryQ[goldenDir],
  Print["CC25 FAIL  golden directory not found: ", goldenDir];
  Print["           Run TEST/phase1_codegen_goldens.wl first to capture goldens."];
  Exit[1];
];
If[!DirectoryQ[workDir],
  CreateDirectory[workDir, CreateIntermediateDirectories -> True]];

eps = Symbol["eps"];

(* ============================================================================
   Result accumulation
   nFail counts differing or missing files; nPass counts byte-identical files.
   ============================================================================ *)
nPass = 0; nFail = 0;

(* Helper: compare generated file against stored golden; print result. *)
checkGolden[label_String, srcFile_String] := Module[{golden, generated},
  If[!FileExistsQ[FileNameJoin[{goldenDir, label}]],
    Print["CC25 FAIL  MISSING golden: ", label,
          "  (run TEST/phase1_codegen_goldens.wl to capture)"];
    nFail++; Return[];
  ];
  If[!FileExistsQ[srcFile],
    Print["CC25 FAIL  generated file missing: ", srcFile];
    nFail++; Return[];
  ];
  golden    = Import[FileNameJoin[{goldenDir, label}], "Text"];
  generated = Import[srcFile, "Text"];
  If[golden === generated,
    Print["CC25 PASS  [byte-identical] ", label,
          "  (", StringLength[generated], " chars)"];
    nPass++,
    Print["CC25 FAIL  [diff] ", label,
          "  expected=", StringLength[golden], " chars",
          "  got=",      StringLength[generated], " chars"];
    nFail++;
  ];
];

(* ============================================================================
   SPEC 1: Convergent, no kinematics  (Pi/8)
   P = 1+x1^2+x2^2,  A = {0,0},  B = -3
   Goldens: conv_nokat_mc.cpp / conv_nokat_vegas_perkp.cpp / conv_nokat_vegas_batch.cpp
   ============================================================================ *)
Print["--- Spec 1: convergent, no kinematics (Pi/8) ---"];
Module[{poly, vars, spec, verts, fan, sd, conv, f},
  poly = 1 + x[1]^2 + x[2]^2;
  vars = {x[1], x[2]};
  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents"-> {-3},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;
  verts = PolytopeVertices[poly^(-3), vars];
  fan   = ComputeDecomposition[verts, "ShowProgress" -> False];
  sd    = Table[
    ProcessSector[spec, fan[[1]], fan[[2, s]], s, "Verbose" -> False],
    {s, Length[fan[[2]]]}
  ];
  conv = Select[sd, AssociationQ[#] && !#["IsDivergent"] &];

  (* MC *)
  f = FileNameJoin[{workDir, "conv_nokat_mc.cpp"}];
  GenerateCppMonteCarlo[conv, {}, spec, f, "Integrator" -> "MonteCarlo"];
  checkGolden["conv_nokat_mc.cpp", f];

  (* Per-kp VEGAS *)
  f = FileNameJoin[{workDir, "conv_nokat_vegas_perkp.cpp"}];
  GenerateCppMonteCarlo[conv, {}, spec, f, "Integrator" -> "Vegas", "Batch" -> False];
  checkGolden["conv_nokat_vegas_perkp.cpp", f];

  (* Batched VEGAS *)
  f = FileNameJoin[{workDir, "conv_nokat_vegas_batch.cpp"}];
  GenerateCppMonteCarlo[conv, {}, spec, f, "Integrator" -> "Vegas", "Batch" -> True];
  checkGolden["conv_nokat_vegas_batch.cpp", f];
];

(* ============================================================================
   SPEC 2: Subtraction divergent 1D  (g0 / g1 / rem)
   P = 1+x1,  A = {-1+eps},  B = -2
   Goldens: subdiv_mc.cpp / subdiv_vegas_perkp.cpp / subdiv_vegas_batch.cpp
   ============================================================================ *)
Print[];
Print["--- Spec 2: subtraction divergent 1D (g0/g1/rem) ---"];
Module[{poly, vars, spec, verts, fan, sd, conv, divSD, procDiv,
        fMC, fVegas, fBatch},
  poly = 1 + x[1];
  vars = {x[1]};
  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {-1 + eps},
    "PolynomialExponents"-> {-2},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;
  verts   = PolytopeVertices[poly^(-1), vars];
  fan     = ComputeDecomposition[verts, "ShowProgress" -> False];
  sd      = Table[
    ProcessSector[spec, fan[[1]], fan[[2, s]], s, "Verbose" -> False],
    {s, Length[fan[[2]]]}
  ];
  conv    = Select[sd, AssociationQ[#] && !#["IsDivergent"] &];
  divSD   = Select[sd, AssociationQ[#] && #["IsDivergent"] &];
  procDiv = Quiet[ProcessDivergentSector[#, spec] & /@ divSD];
  procDiv = Select[procDiv, AssociationQ];

  (* MC (eps->0 on convergent sectors for subtraction method) *)
  fMC = FileNameJoin[{workDir, "subdiv_mc.cpp"}];
  GenerateCppMonteCarlo[conv /. eps -> 0, procDiv, spec, fMC,
                        "Integrator" -> "MonteCarlo"];
  checkGolden["subdiv_mc.cpp", fMC];

  (* Per-kp VEGAS *)
  fVegas = FileNameJoin[{workDir, "subdiv_vegas_perkp.cpp"}];
  GenerateCppMonteCarlo[conv /. eps -> 0, procDiv, spec, fVegas,
                        "Integrator" -> "Vegas", "Batch" -> False];
  checkGolden["subdiv_vegas_perkp.cpp", fVegas];

  (* Batched VEGAS (falls back to per-kp on divergent path; still exercises the option) *)
  fBatch = FileNameJoin[{workDir, "subdiv_vegas_batch.cpp"}];
  GenerateCppMonteCarlo[conv /. eps -> 0, procDiv, spec, fBatch,
                        "Integrator" -> "Vegas", "Batch" -> True];
  checkGolden["subdiv_vegas_batch.cpp", fBatch];
];

(* ============================================================================
   SPEC 3: IBP divergent 1D  (ibp-boundary / ibp-term)
   Same polynomial/exponents as Spec 2, IBP path.
   Goldens: ibpdiv_mc.cpp / ibpdiv_vegas_perkp.cpp / ibpdiv_vegas_batch.cpp
   ============================================================================ *)
Print[];
Print["--- Spec 3: IBP divergent 1D (ibp-boundary/ibp-term) ---"];
Module[{poly, vars, spec, verts, fan, sd, conv, divSD, ibpSD,
        fMC, fVegas, fBatch},
  poly = 1 + x[1];
  vars = {x[1]};
  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {-1 + eps},
    "PolynomialExponents"-> {-2},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;
  verts = PolytopeVertices[poly^(-1), vars];
  fan   = ComputeDecomposition[verts, "ShowProgress" -> False];
  sd    = Table[
    ProcessSector[spec, fan[[1]], fan[[2, s]], s, "Verbose" -> False],
    {s, Length[fan[[2]]]}
  ];
  conv  = Select[sd, AssociationQ[#] && !#["IsDivergent"] &];
  divSD = Select[sd, AssociationQ[#] && #["IsDivergent"] &];
  ibpSD = Quiet[IBPProcessSector[#, spec] & /@ divSD];
  ibpSD = Select[ibpSD, AssociationQ];

  (* MC *)
  fMC = FileNameJoin[{workDir, "ibpdiv_mc.cpp"}];
  GenerateCppMonteCarloIBP[conv /. eps -> 0, ibpSD, spec, fMC,
                            "Integrator" -> "MonteCarlo"];
  checkGolden["ibpdiv_mc.cpp", fMC];

  (* Per-kp VEGAS *)
  fVegas = FileNameJoin[{workDir, "ibpdiv_vegas_perkp.cpp"}];
  GenerateCppMonteCarloIBP[conv /. eps -> 0, ibpSD, spec, fVegas,
                            "Integrator" -> "Vegas", "Batch" -> False];
  checkGolden["ibpdiv_vegas_perkp.cpp", fVegas];

  (* Batched VEGAS *)
  fBatch = FileNameJoin[{workDir, "ibpdiv_vegas_batch.cpp"}];
  GenerateCppMonteCarloIBP[conv /. eps -> 0, ibpSD, spec, fBatch,
                            "Integrator" -> "Vegas", "Batch" -> True];
  checkGolden["ibpdiv_vegas_batch.cpp", fBatch];
];

(* ============================================================================
   Summary
   ============================================================================ *)
Print[];
Print["======================================================================="];
Print["CC25 SUMMARY: ", nPass, " PASS  /  ", nFail, " FAIL  (", nPass + nFail, " files checked)"];
Print["Golden dir: ", goldenDir];
Print["======================================================================="];

If[nFail == 0,
  Print["CC25 PASS  all ", nPass, " golden files byte-identical"],
  Print["CC25 FAIL  expected=0 diffs  got=", nFail, " file(s) differ"]
];
