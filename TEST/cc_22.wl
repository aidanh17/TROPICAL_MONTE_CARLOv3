(* ============================================================================
   TEST/cc_22.wl  —  Cross-check #22: Determinism / seed reproducibility
                      (plan.md §8.2 #22, §4.3 SeedBase note)

   PASS criteria (plan.md §8.2 #22):
     (A) Determinism   : same seed_base -> byte-identical result (run twice)
     (B) Override wired: omitting argv[5] == passing the compile-time SeedBase=42
     (C) Distinct stream: seeds 42, 99, 12345 give different MC estimates
     (D) Correctness   : every stream within 5 sigma of Pi/8

   Integrand:  ∫_{[0,∞)²} dx1 dx2 / (1 + x1² + x2²)³ = Pi/8
   (simple, two-dimensional, no kinematics — isolates the seed plumbing)

   Ported from OLD_CODE/.../EXAMPLES/test_seedbase.wl (Tree B, canonical
   source for SeedBase / argv[5] override).

   Tier 1: WL + g++ only (no CUBA required; uses the MC back-end).

   Run:  wolframscript -file TEST/cc_22.wl
         (from the TROPICAL_MONTE_CARLO3 root, or anywhere — uses absolute paths)
   ============================================================================ *)

(* ── 0. Locate the package root via this file's own path ─────────────────── *)
$pkgRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
$evalWL  = FileNameJoin[{$pkgRoot, "tropical_eval.wl"}];
$fanWL   = FileNameJoin[{$pkgRoot, "tropical_fan.wl"}];
$ioDir   = FileNameJoin[{$pkgRoot, "TEST", "INTERFILES"}];

If[!FileExistsQ[$evalWL],
  Print["CC22 FAIL  tropical_eval.wl not found at ", $evalWL]; Quit[1]];

Get[$fanWL];
Get[$evalWL];

Quiet[CreateDirectory[$ioDir], {CreateDirectory::eexist}];

Print["CC22: tropical_eval.wl loaded"];
Print["CC22: working directory for generated files: ", $ioDir];
Print[];

(* ── 1. Build the integrand spec and fan (once) ─────────────────────────── *)
Module[{poly, vars, spec, verts, fanData, dualVerts, simplices, sectors,
        cppFile, binary, kinFile,
        nSamples = 400000,
        exact = N[Pi/8],
        runSeed,
        rNo5, r42, r42b, r99, r12345,
        within,
        passAll = True,
        nPass, nFail},

  poly = 1 + x[1]^2 + x[2]^2;
  vars = {x[1], x[2]};
  spec = <|"Polynomials"        -> {poly},
           "MonomialExponents"  -> {0, 0},
           "PolynomialExponents"-> {-3},
           "Variables"          -> vars,
           "KinematicSymbols"   -> {},
           "RegulatorSymbol"    -> None|>;

  Print["CC22: computing tropical fan for (1 + x1^2 + x2^2)^(-3) ..."];
  verts   = PolytopeVertices[poly^(-3), vars];
  fanData = computeFanScaled[verts];
  If[fanData === $Failed,
    Print["CC22 FAIL  expected=fan  got=$Failed (fan computation failed)"]; Quit[1]];
  {dualVerts, simplices} = fanData;

  sectors = Select[
    Table[ProcessSector[spec, dualVerts, simplices[[s]], s], {s, Length[simplices]}],
    AssociationQ];
  Print["CC22: fan ", Length[dualVerts], " rays, ", Length[simplices],
        " simplices, ", Length[sectors], " valid sectors"];
  Print["CC22: exact = Pi/8 = ", exact];
  Print[];

  (* ── 2. Compile ONCE with compile-time SeedBase = 42 ────────────────── *)
  cppFile = FileNameJoin[{$ioDir, "cc22_seedtest_generated.cpp"}];
  binary  = FileNameJoin[{$ioDir, "cc22_seedtest"}];
  kinFile = FileNameJoin[{$ioDir, "cc22_kin.txt"}];

  GenerateCppMonteCarlo[sectors, {}, spec, cppFile,
    "NSamples" -> nSamples, "SeedBase" -> 42];
  If[CompileCpp[cppFile, binary, False] === $Failed,
    Print["CC22 FAIL  expected=compiled  got=$Failed (compilation failed)"]; Quit[1]];
  Print["CC22: compiled binary: ", FileNameTake[binary]];

  (* One (trivial) kinematic point — no kinematics in this integral *)
  Export[kinFile, "1\n", "Text"];

  (* ── 3. Helper: run the binary with an optional seed argv ────────────── *)
  (* runSeed[seedArgOrNone] -> {re, im, reErr, imErr} from the results file *)
  runSeed[seedArg_] := Module[{outFile, proc, tbl},
    outFile = FileNameJoin[{$ioDir,
               "cc22_res_" <> ToString[seedArg] <> ".txt"}];
    proc = RunProcess[Join[
             {binary, kinFile, outFile,
              ToString[nSamples], "4"},
             If[seedArg === None, {}, {ToString[seedArg]}]]];
    If[proc["ExitCode"] =!= 0,
      Print["CC22: binary run FAILED (seed ", seedArg, "): ",
            proc["StandardError"]];
      Return[$Failed]];
    tbl = Import[outFile, "Table"];
    If[!ListQ[tbl] || Length[tbl] < 1,
      Print["CC22: result file unreadable"]; Return[$Failed]];
    tbl[[1]]
  ];

  (* ── 4. Run five times reusing the compiled binary ───────────────────── *)
  Print["CC22: running five seed trials ..."];
  rNo5   = runSeed[None];    (* no argv[5] -> compile-time 42 *)
  r42    = runSeed[42];
  r42b   = runSeed[42];      (* second run of same seed *)
  r99    = runSeed[99];
  r12345 = runSeed[12345];

  If[MemberQ[{rNo5, r42, r42b, r99, r12345}, $Failed],
    Print["CC22 FAIL  expected=all-runs  got=$Failed (one or more seed trials failed)"];
    Quit[1]];

  Print[];
  Print["  no-argv5 : ", rNo5];
  Print["  seed  42 : ", r42];
  Print["  seed  42': ", r42b];
  Print["  seed  99 : ", r99];
  Print["  seed 12345: ", r12345];
  Print[];

  nPass = 0; nFail = 0;

  (* ── (A) Determinism: seed 42 run twice -> byte-identical results ──── *)
  If[r42 === r42b,
    nPass++;
    Print["  (A) determinism (seed 42 == seed 42'): PASS"],
    nFail++;
    passAll = False;
    Print["  (A) determinism: FAIL  expected=", r42, "  got=", r42b]];

  (* ── (B) Override wired: no argv[5] == argv[5]=42 ──────────────────── *)
  If[rNo5 === r42,
    nPass++;
    Print["  (B) argv[5] override (no-arg == 42): PASS"],
    nFail++;
    passAll = False;
    Print["  (B) argv[5] override: FAIL  expected=", r42, "  got=", rNo5]];

  (* ── (C) Distinct streams: 42, 99, 12345 give different Re ─────────── *)
  If[r42[[1]] =!= r99[[1]] && r42[[1]] =!= r12345[[1]] && r99[[1]] =!= r12345[[1]],
    nPass++;
    Print["  (C) distinct streams (42 / 99 / 12345 differ): PASS"],
    nFail++;
    passAll = False;
    Print["  (C) distinct streams: FAIL  expected=distinct  got=",
          {r42[[1]], r99[[1]], r12345[[1]]}]];

  (* ── (D) Correctness: every stream within 5 sigma of Pi/8 ─────────── *)
  within[r_] := ListQ[r] && Length[r] >= 3 && r[[3]] > 0 &&
                Abs[r[[1]] - exact] <= 5 r[[3]];
  If[AllTrue[{rNo5, r42, r99, r12345}, within],
    nPass++;
    Print["  (D) all streams within 5 sigma of Pi/8: PASS"],
    nFail++;
    passAll = False;
    Print["  (D) correctness: FAIL  expected=within-5-sigma  got=devs-in-sigma=",
          Map[If[#[[3]] > 0, Abs[#[[1]] - exact] / #[[3]], Infinity] &,
              {rNo5, r42, r99, r12345}]]];

  Print[];
  Print["================================================================"];
  Print["CC22  Determinism / seed reproducibility  —  summary"];
  Print["================================================================"];
  Print["  (A) determinism         ", If[r42 === r42b,        "PASS", "FAIL"]];
  Print["  (B) argv[5] override    ", If[rNo5 === r42,        "PASS", "FAIL"]];
  Print["  (C) distinct streams    ",
        If[r42[[1]] =!= r99[[1]] && r42[[1]] =!= r12345[[1]] &&
           r99[[1]] =!= r12345[[1]],                          "PASS", "FAIL"]];
  Print["  (D) correctness (5 sig) ", If[AllTrue[{rNo5, r42, r99, r12345}, within],
                                        "PASS", "FAIL"]];
  Print["----------------------------------------------------------------"];
  If[passAll,
    Print["CC22 PASS  all 4 sub-checks passed  (determinism + seed override + ",
          "distinct streams + correctness within 5 sigma)"],
    Print["CC22 FAIL  expected=", nPass + nFail, " pass  got=",
          nFail, " fail  (", nPass, " passed, ", nFail, " failed)"]];
  Print["================================================================"];

  If[!passAll, Quit[1]]
]
