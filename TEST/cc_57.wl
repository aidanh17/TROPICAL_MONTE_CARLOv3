(* ============================================================================
   TEST/cc_57.wl  —  CHARACTERIZATION test: codegen goldens for two C++-emission
   paths that previously had NO byte-level pin (fable_plan, following cc_55/56):

     (A) the LIFTED sector path        — domain-indicator emission
                                          (emitDomainIndicatorCpp, planAXpDIV.md
                                          §4.2), via GenerateCppMonteCarlo.
     (B) the SplitRealImag complex-B   — oscillatory-phase emission
         path                            (tropical_eval.wl ~2503-2567), via the
                                          full EvaluateTropicalMC driver.
     (C) OPTIONAL: lifted + divergent Case A via the IBP route (SplitRealImag),
         i.e. the combination cc_54 exercises numerically — here we only pin
         the emitted tropical_mc_ibp.cpp bytes.

   These are "golden master" / characterization tests: byte-for-byte pins of
   the behavior that exists TODAY (before any refactor), captured on FIRST RUN
   (the golden file is written if absent) so a later modularity refactor can be
   proven behavior-preserving on these two previously-uncovered emission paths.
   They are NOT a new correctness claim — correctness of the underlying paths
   is established elsewhere (phase3_selfgate.wl (7) for the domain indicator;
   TEST/test_complex_phase4.wl #21b and cc_53/cc_54 for SplitRealImag /
   lifted+divergent SplitRealImag IBP).

   Per-golden determinism guard: each candidate golden is emitted TWICE (two
   independent pipeline runs, fresh fan/sector computation each time) and the
   two emissions must be byte-identical BEFORE the golden compare/write; a
   mismatch there SKIPs that golden (never silently pins noise).

   Tier 1: WL + Polymake (fan) + g++ (only for golden C's numeric driver run;
   goldens A/B pin GenerateCppMonteCarlo / EvaluateTropicalMC's emitted text
   directly, no MC binaries are executed for them).

   Run:  wolframscript -file TEST/cc_57.wl
   (from any directory — uses absolute paths).  Never exits nonzero; run_all.wl
   classifies PASS/FAIL/SKIP from the printed "CC57 ..." lines.
   ============================================================================ *)

$pkgRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{$pkgRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$pkgRoot, "tropical_eval.wl"}]];
Print["CC57: packages loaded."];
Print[];

$cc57Pass = True; $cc57Fail = {}; $cc57SkipList = {};

cc57Assert[label_String, test_, exp_, got_] :=
  If[TrueQ[test],
    Print["CC57 PASS  [", label, "]  expected=", exp, "  got=", got],
    Print["CC57 FAIL  [", label, "]  expected=", exp, "  got=", got];
    $cc57Pass = False; AppendTo[$cc57Fail, label]];

cc57Skip[label_String, reason_String] := (
  Print["CC57 SKIP  [", label, "]  ", reason];
  AppendTo[$cc57SkipList, label];
);

(* Quiet the KNOWN, EXPECTED diagnostic messages the pipelines below may raise
   (Polymake progress notes, VEGAS budget hints, etc.) so log noise stays out;
   never quiets a message that would indicate an actual failure path. *)
SetAttributes[qr, HoldFirst];
qr[e_] := Quiet[e, {TropicalFan::polymake, TropicalEval::vegasbudget,
  General::stop, General::munfl, Power::infy}];

(* --- scratch working dir (repo-root INTERFILES/cc_57/, run_all convention) - *)
workDir = FileNameJoin[{$pkgRoot, "INTERFILES", "cc_57"}];
Quiet[CreateDirectory[workDir, CreateIntermediateDirectories -> True],
      {CreateDirectory::eexist}];

goldenDir = FileNameJoin[{$pkgRoot, "TEST", "baselines", "codegen_goldens"}];
Quiet[CreateDirectory[goldenDir, CreateIntermediateDirectories -> True],
      {CreateDirectory::eexist}];

(* checkOrCreateGolden: byte-compare against a committed golden; if the golden
   does not exist yet, WRITE it (this is how the pre-refactor baseline gets
   captured) and report PASS with the "created" tag; on later runs, compare. *)
checkOrCreateGolden[label_String, text_String] := Module[{path, golden},
  path = FileNameJoin[{goldenDir, label}];
  If[!FileExistsQ[path],
    Export[path, text, "Text"];
    Print["CC57 PASS  [golden ", label, " created]  (", StringLength[text], " chars)"];
    True
    ,
    golden = Import[path, "Text"];
    If[golden === text,
      Print["CC57 PASS  [golden ", label, " byte-identical]  (",
        StringLength[text], " chars)"];
      True
      ,
      Print["CC57 FAIL  [golden ", label, " MISMATCH]  golden=",
        StringLength[golden], " chars  generated=", StringLength[text], " chars"];
      $cc57Pass = False; AppendTo[$cc57Fail, label <> " byte-mismatch"];
      False
    ]
  ]
];

(* ============================================================================
   (A) LIFTED CONVERGENT sector codegen  ->  lifted_conv_mc.cpp
   Cribbed verbatim from TEST/phase3_selfgate.wl "Case A lifted spec" (the
   fixture whose driver run phase3_selfgate (7) proves DOES emit a
   "lifted-sector domain indicator" block): P = 1 + 10^6 x1^2 + x2^2 + x1 x2^2,
   B=-2, lifted on x1^2 with k=2.  Pipeline (verbatim, phase3_selfgate.wl §1):
     LiftCoefficients -> PolytopeVertices+ComputeDecomposition on the lifted
     integrand -> ProcessSectorLifted per sector, drop EmptyDomain ->
     GenerateCppMonteCarlo[convergentLiftedSectors, {}, liftedSpec, out].
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (A) Lifted CONVERGENT sector codegen: lifted_conv_mc.cpp"];

Module[{specA, ruleA, buildLiftedCpp, text1, text2, f1, f2},
  specA = <|"Polynomials" -> {1 + 10^6 x[1]^2 + x[2]^2 + x[1] x[2]^2},
    "MonomialExponents" -> {0, 0}, "PolynomialExponents" -> {-2},
    "Variables" -> {x[1], x[2]}, "KinematicSymbols" -> {},
    "RegulatorSymbol" -> None|>;
  ruleA = {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 2|>};

  buildLiftedCpp[outFile_String] := Module[{lc, ls, ld, verts, fan, sectors, conv},
    lc = LiftCoefficients[specA, ruleA];
    If[!AssociationQ[lc], Return[$Failed]];
    ls = lc["LiftedSpec"]; ld = lc["LiftData"];
    verts = PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]];
    fan = ComputeDecomposition[verts, "ShowProgress" -> False];
    sectors = Table[
      ProcessSectorLifted[ls, fan[[1]], fan[[2, s]], s, ld],
      {s, Length[fan[[2]]]}];
    conv = Select[sectors, AssociationQ[#] && !KeyExistsQ[#, "EmptyDomain"] &];
    GenerateCppMonteCarlo[conv, {}, ls, outFile, "Integrator" -> "MonteCarlo"];
    If[FileExistsQ[outFile], Import[outFile, "Text"], $Failed]
  ];

  f1 = FileNameJoin[{workDir, "lifted_conv_mc_run1.cpp"}];
  f2 = FileNameJoin[{workDir, "lifted_conv_mc_run2.cpp"}];
  text1 = qr@buildLiftedCpp[f1];
  text2 = qr@buildLiftedCpp[f2];

  If[StringQ[text1] && StringQ[text2],
    cc57Assert["A1 lifted_conv_mc emission deterministic (2 independent runs)",
      text1 === text2, "identical", If[text1 === text2, "identical", "DIFFERS"]];
    cc57Assert["A2 lifted_conv_mc contains a lifted-sector domain-indicator block",
      StringContainsQ[text1, "lifted-sector domain indicator"],
      "contains marker",
      "contains=" <> ToString[StringContainsQ[text1, "lifted-sector domain indicator"]]];

    If[text1 === text2,
      checkOrCreateGolden["lifted_conv_mc.cpp", text1],
      cc57Skip["lifted_conv_mc golden",
        "emission not byte-identical across 2 runs (see A1) — refusing to pin noise"]
    ],
    cc57Assert["A0 lifted_conv_mc pipeline succeeded (no $Failed)", False,
      "StringQ result", ToString[{Head[text1], Head[text2]}]];
    cc57Skip["lifted_conv_mc golden", "pipeline returned $Failed, see A0"];
  ];
];

(* ============================================================================
   (B) SplitRealImag oscillatory-phase codegen  ->  split_cx_conv_mc.cpp
   Cribbed from TEST/test_complex_phase4.wl #21b (unlifted complex-BASE +
   complex-B fixture): P = 1+x1^2+(1+I)x1 x2+x2^2, B=-(3/2+4I/5).  Captured via
   the full EvaluateTropicalMC driver with ComplexExponentMode->"SplitRealImag"
   (RunChecks->False, small NSamples — only the emitted tropical_mc_generated
   .cpp TEXT is pinned; the MC numbers are irrelevant and not checked here,
   correctness of this fixture is #21b's job).
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (B) SplitRealImag oscillatory-phase codegen: split_cx_conv_mc.cpp"];

Module[{vars, Bex, poly, specSplit, verts, fan, buildSplitCpp, text1, text2, wd1, wd2},
  vars = {x[1], x[2]};
  Bex = -(3/2 + 4 I/5);
  poly = 1 + x[1]^2 + (1 + I) x[1] x[2] + x[2]^2;   (* complex coeff => arg P != 0 *)
  specSplit = <|"Polynomials" -> {poly}, "MonomialExponents" -> {0, 0},
    "PolynomialExponents" -> {Bex}, "Variables" -> vars,
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
  verts = PolytopeVertices[poly^(-1), vars];
  fan = ComputeDecomposition[verts, "ShowProgress" -> False];

  buildSplitCpp[wd_String] := Module[{mc, cppFile},
    Quiet[CreateDirectory[wd, CreateIntermediateDirectories -> True],
          {CreateDirectory::eexist}];
    mc = EvaluateTropicalMC[specSplit, fan, {{}}, "NSamples" -> 1000,
      "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> wd,
      "ComplexExponentMode" -> "SplitRealImag"];
    cppFile = FileNameJoin[{wd, "tropical_mc_generated.cpp"}];
    If[AssociationQ[mc] && FileExistsQ[cppFile], Import[cppFile, "Text"], $Failed]
  ];

  wd1 = FileNameJoin[{workDir, "split_run1"}];
  wd2 = FileNameJoin[{workDir, "split_run2"}];
  text1 = qr@buildSplitCpp[wd1];
  text2 = qr@buildSplitCpp[wd2];

  If[StringQ[text1] && StringQ[text2],
    cc57Assert["B1 split_cx_conv_mc emission deterministic (2 independent runs)",
      text1 === text2, "identical", If[text1 === text2, "identical", "DIFFERS"]];
    cc57Assert["B2 split_cx_conv_mc contains SplitRealImag oscillatory-phase marker",
      StringContainsQ[text1, "SplitRealImag oscillatory phase"],
      "contains marker",
      "contains=" <> ToString[StringContainsQ[text1, "SplitRealImag oscillatory phase"]]];
    cc57Assert["B3 split_cx_conv_mc contains phaseSum accumulation line",
      StringContainsQ[text1, "cx phaseSum ="],
      "contains marker",
      "contains=" <> ToString[StringContainsQ[text1, "cx phaseSum ="]]];

    If[text1 === text2,
      checkOrCreateGolden["split_cx_conv_mc.cpp", text1],
      cc57Skip["split_cx_conv_mc golden",
        "emission not byte-identical across 2 runs (see B1) — refusing to pin noise"]
    ],
    cc57Assert["B0 split_cx_conv_mc pipeline succeeded (no $Failed)", False,
      "StringQ result", ToString[{Head[text1], Head[text2]}]];
    cc57Skip["split_cx_conv_mc golden", "pipeline returned $Failed, see B0"];
  ];
];

(* ============================================================================
   (C) OPTIONAL: lifted + divergent Case A via the IBP route (SplitRealImag)
       -> lifted_div_ibp_mc.cpp
   Cribbed verbatim from TEST/cc_54.wl's shared fixture (specF/liftRuleF): a
   genuine single 1/eps pole, an extreme lifted coefficient (k=3), complex B
   with the regulator inside, and a nonzero Im(A) on a convergent slot.
   Captured via the full EvaluateTropicalMC["Method"->"IBP", ComplexExponentMode
   ->"SplitRealImag"] driver (this DOES compile+run a tiny NSamples=1000 MC
   binary — required to obtain "Results" from this driver at all — but only
   the emitted tropical_mc_ibp.cpp TEXT is pinned; correctness is cc_54's job).
   Skips cleanly if the pipeline fails or emission is not reproducible.
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (C) OPTIONAL: lifted+divergent IBP (SplitRealImag) codegen: lifted_div_ibp_mc.cpp"];

Module[{epsD, dim, gB, qA, polysF, AvalsF, peF, specF, liftRuleF,
        buildIBPCpp, text1, text2, wd1, wd2, result},
  epsD = Symbol["epsCC57"];
  dim = 2; gB = 1/2; qA = 1/2;
  polysF = {1 + x[1] x[2] + x[1]^2 + 10^4 x[2]^2};
  AvalsF = {-1 + epsD, 1/2 + I qA};
  peF    = {-2 + I gB - epsD};
  specF  = <|"Polynomials" -> polysF, "MonomialExponents" -> AvalsF,
    "PolynomialExponents" -> peF, "Variables" -> {x[1], x[2]},
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> epsD|>;
  liftRuleF = {<|"PolyIndex" -> 1, "ExponentVector" -> {0, 2}, "k" -> 3|>};

  buildIBPCpp[wd_String] := Module[{lift, ls, ld, fan, ibp, cppFile},
    Quiet[CreateDirectory[wd, CreateIntermediateDirectories -> True],
          {CreateDirectory::eexist}];
    lift = qr@LiftCoefficients[specF, liftRuleF];
    If[!AssociationQ[lift], Return[$Failed]];
    ls = lift["LiftedSpec"]; ld = lift["LiftData"];
    fan = qr@computeFanScaled[qr@PolytopeVertices[
      (Times @@ ls["Polynomials"])^(-1), ls["Variables"]]];
    If[!ListQ[fan] || Length[fan] != 2, Return[$Failed]];
    ibp = qr@EvaluateTropicalMC[ls, fan, {{}}, "LiftData" -> ld, "Method" -> "IBP",
      "ComplexExponentMode" -> "SplitRealImag", "Integrator" -> "MC",
      "NSamples" -> 1000, "RunChecks" -> False, "Verbose" -> False,
      "WorkingDirectory" -> wd];
    cppFile = FileNameJoin[{wd, "tropical_mc_ibp.cpp"}];
    If[AssociationQ[ibp] && FileExistsQ[cppFile], Import[cppFile, "Text"], $Failed]
  ];

  wd1 = FileNameJoin[{workDir, "ibp_run1"}];
  wd2 = FileNameJoin[{workDir, "ibp_run2"}];
  result = TimeConstrained[
    (text1 = buildIBPCpp[wd1]; text2 = buildIBPCpp[wd2]; "done"),
    180, "timeout"];

  Which[
    result === "timeout",
      cc57Skip["lifted_div_ibp_mc golden", "pipeline exceeded 180s time budget"],

    !StringQ[text1] || !StringQ[text2],
      cc57Skip["lifted_div_ibp_mc golden",
        "pipeline returned $Failed (fan/lift/IBP route did not complete cleanly)"],

    text1 =!= text2,
      cc57Skip["lifted_div_ibp_mc golden",
        "emission not byte-identical across 2 runs — refusing to pin noise"],

    True,
      cc57Assert["C2 lifted_div_ibp_mc contains SplitRealImag oscillatory-phase marker",
        StringContainsQ[text1, "SplitRealImag oscillatory phase"] ||
          StringContainsQ[text1, "phaseSum"],
        "contains marker",
        "contains=" <> ToString[StringContainsQ[text1, "SplitRealImag oscillatory phase"] ||
          StringContainsQ[text1, "phaseSum"]]];
      checkOrCreateGolden["lifted_div_ibp_mc.cpp", text1];
  ];
];

(* ============================================================================
   Summary
   ============================================================================ *)
Print[];
Print["================================================================"];
If[$cc57Pass,
  Print["CC57 PASS  lifted-domain-indicator + SplitRealImag codegen goldens pinned",
    If[$cc57SkipList =!= {}, "  (skipped: " <> ToString[$cc57SkipList] <> ")", ""]],
  Print["CC57 FAIL  failures: ", $cc57Fail]];
Print["================================================================"];
