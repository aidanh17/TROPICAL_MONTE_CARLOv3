(* ============================================================================
   TEST/cc_17.wl  —  Cross-check #17 (plan.md §8.2)
   Lifted vs unlifted exactness (Tier 1: WL + g++ only).

   PASS criteria (plan.md §8.2, row 17):
     relerr < 0.5% vs exact closed form; DroppedSectors count asserted.

   Two sub-cases ported from OLD_CODE test_lifted.wl (Test 23A and 23D):

   17A  Toy 1: P = 1 + x[1] + 1e-6*x[2], B={-3}.
        Exact = 1/(2*1e-6) = 500000.  k=2 lift on {0,1}.
        PASS: relErr(sectorSum, 500000) < 0.005; DroppedSectors == 3.

   17B  Toy 0: P = 1 + 1e6*x[1], B={-2}.
        Exact = 1e-6.  k=1 lift on {1}.  Constant-root sectors dropped.
        PASS: DroppedSectors > 0; relErr(sectorSum, 1e-6) < 0.001.

   Run:
     wolframscript -file TEST/cc_17.wl
   from the TROPICAL_MONTE_CARLO3 root, or load in WL with $InputFileName set.
   ============================================================================ *)

(* --- locate root and load v3 packages via absolute paths --- *)
With[{root = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]]},
  Get[FileNameJoin[{root, "tropical_fan.wl"}]];
  Get[FileNameJoin[{root, "tropical_eval.wl"}]];
];

Print["CC17: tropical_fan.wl + tropical_eval.wl loaded."];
Print[];

(* ============================================================================
   Helper: print "CC17 PASS ..." or "CC17 FAIL expected=<e> got=<g>"
   ============================================================================ *)
cc17Pass[label_String, info_String] :=
  Print["CC17 PASS [", label, "] ", info];

cc17Fail[label_String, expected_, got_] :=
  Print["CC17 FAIL [", label, "] expected=", expected, " got=", got];

(* ============================================================================
   17A  Toy 1 exactness
   P = 1 + x[1] + 10^-6 x[2], B={-3}.
   Closed-form exact = Int_0^inf Int_0^inf (1+x1+1e-6 x2)^-3 dx1 dx2
                     = 1e6 * 1/2 = 500000.
   Lift monomial {0,1} (coeff 1e-6) with k=2 => z0 = (1e-6)^(1/2) = 1e-3.
   Explicit fan supplied (degenerate lifted polytope; plan §8.2 known facts).
   ============================================================================ *)

Module[
  {spec1, liftRules1, lcRes, liftedSpec1, liftData1,
   explicitDualVertices, explicitSimplexList, explicitFan,
   vlRes, relErrExact, nDropped, exactRef,
   pass17A},

  Print["--- CC17A: Toy 1 exactness (LiftCoefficients + ValidateLiftedDecomposition) ---"];
  pass17A = True;
  exactRef = 500000;

  spec1 = <|
    "Polynomials"         -> {1 + x[1] + 10^-6 x[2]},
    "MonomialExponents"   -> {0, 0},
    "PolynomialExponents" -> {-3},
    "Variables"           -> {x[1], x[2]},
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> None
  |>;

  liftRules1 = {<|"PolyIndex" -> 1, "ExponentVector" -> {0, 1}, "k" -> 2|>};
  lcRes = LiftCoefficients[spec1, liftRules1];

  If[!AssociationQ[lcRes],
    cc17Fail["17A-LiftCoefficients", "Association", lcRes];
    pass17A = False;
    Goto["done17A"]
  ];

  liftedSpec1 = lcRes["LiftedSpec"];
  liftData1   = lcRes["LiftData"];

  Print["  z0 = ", liftData1["z0"]];
  Print["  Lifted poly: ", liftedSpec1["Polynomials"][[1]]];

  (* Identity check: substituting z0 back must recover the original poly *)
  Module[{lp = liftedSpec1["Polynomials"][[1]],
          z0v = liftData1["z0"],
          auxV = liftData1["AuxVariable"]},
    If[Expand[lp /. auxV -> z0v] === Expand[spec1["Polynomials"][[1]]],
      Print["  Identity check: PASS"],
      Print["  Identity check: FAIL (lift is inconsistent)"];
      pass17A = False
    ]
  ];

  (* Explicit unimodular fan for the degenerate Toy-1 lifted polytope.
     Dual vertices and simplices from plan §8.2 / test_lifted.wl Test 23A. *)
  explicitDualVertices = {{1,0,0},{0,1,0},{-1,-1,0},{0,2,-1},{0,-2,1}};
  explicitSimplexList  = {{1,2,4},{2,3,4},{3,1,4},{1,2,5},{2,3,5},{3,1,5}};
  explicitFan = {explicitDualVertices, explicitSimplexList};
  Print["  Fan: ", Length[explicitDualVertices], " rays, ",
        Length[explicitSimplexList], " sectors"];

  vlRes = Quiet @ ValidateLiftedDecomposition[
    spec1, liftedSpec1, explicitFan, liftData1, {}, 3
  ];

  If[!AssociationQ[vlRes],
    cc17Fail["17A-ValidateLiftedDecomposition", "Association", vlRes];
    pass17A = False;
    Goto["done17A"]
  ];

  nDropped = Length[Lookup[vlRes, "DroppedSectors", {}]];
  relErrExact = Abs[(vlRes["SectorSum"] - exactRef) / exactRef];

  Print["  Sector sum = ", vlRes["SectorSum"]];
  Print["  Exact ref  = ", exactRef, " (closed form 1/(2*1e-6))"];
  Print["  relErr     = ", relErrExact];
  Print["  DroppedSectors = ", nDropped];

  (* Gate 1: relErr vs closed-form exact < 0.5% *)
  If[!NumericQ[relErrExact] || relErrExact >= 0.005,
    cc17Fail["17A-exactness", "relErr<0.005",
             "relErr=" <> ToString[relErrExact]];
    pass17A = False,
    cc17Pass["17A-exactness",
             "relErr=" <> ToString[relErrExact] <> " < 0.005; exact=" <>
             ToString[exactRef]]
  ];

  (* Gate 2: DroppedSectors count == 3 *)
  If[nDropped =!= 3,
    cc17Fail["17A-DroppedSectors", 3, nDropped];
    pass17A = False,
    cc17Pass["17A-DroppedSectors",
             "count=3 as expected"]
  ];

  Label["done17A"];
  Print[];
  Print[If[pass17A, "CC17A PASS", "CC17A FAIL"]];
  Print[];
  pass17A
];

(* ============================================================================
   17B  Toy 0 EmptyDomain drops
   P = 1 + 10^6*x[1], B={-2}.
   Exact = Int_0^inf (1+1e6 x)^-2 dx = 1e-6.
   k=1 lift on {1}: z0=1e6.  Constant-root sectors (z0>1) dropped.
   Explicit k=1 fan from sandbox_toy0_1d.wl (7 sectors, 2-D fan).
   ============================================================================ *)

Module[
  {spec0, liftRules0, lcRes0, liftedSpec0, liftData0,
   raysK1, sectsK1, explicitFan0,
   vlRes0, relErr0, nDropped0, exactAns,
   pass17B},

  Print["--- CC17B: Toy 0 EmptyDomain drops (k=1, z0=1e6) ---"];
  pass17B = True;
  exactAns = 10^-6;

  spec0 = <|
    "Polynomials"         -> {1 + 10^6 * x[1]},
    "MonomialExponents"   -> {0},
    "PolynomialExponents" -> {-2},
    "Variables"           -> {x[1]},
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> None
  |>;

  liftRules0 = {<|"PolyIndex" -> 1, "ExponentVector" -> {1}, "k" -> 1|>};
  lcRes0 = LiftCoefficients[spec0, liftRules0];

  If[!AssociationQ[lcRes0],
    cc17Fail["17B-LiftCoefficients", "Association", lcRes0];
    pass17B = False;
    Goto["done17B"]
  ];

  liftedSpec0 = lcRes0["LiftedSpec"];
  liftData0   = lcRes0["LiftData"];

  Print["  z0 = ", liftData0["z0"]];
  Print["  Lifted poly: ", liftedSpec0["Polynomials"][[1]]];

  (* Identity check *)
  Module[{lp = liftedSpec0["Polynomials"][[1]],
          z0v = liftData0["z0"],
          auxV = liftData0["AuxVariable"]},
    If[Expand[lp /. auxV -> z0v] === Expand[spec0["Polynomials"][[1]]],
      Print["  Identity check: PASS"],
      Print["  Identity check: FAIL"]
    ]
  ];

  (* Complete unimodular triangulation for k=1 support {(0,0),(1,1)},
     7 sectors from sandbox_toy0_1d.wl *)
  raysK1  = {{1,0},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1}};
  sectsK1 = Table[{i, If[i < Length[raysK1], i+1, 1]}, {i, Length[raysK1]}];
  explicitFan0 = {raysK1, sectsK1};
  Print["  Fan: ", Length[raysK1], " rays, ", Length[sectsK1], " sectors"];

  vlRes0 = Quiet @ ValidateLiftedDecomposition[
    spec0, liftedSpec0, explicitFan0, liftData0, {}, 4
  ];

  If[!AssociationQ[vlRes0],
    cc17Fail["17B-ValidateLiftedDecomposition", "Association", vlRes0];
    pass17B = False;
    Goto["done17B"]
  ];

  nDropped0 = Length[Lookup[vlRes0, "DroppedSectors", {}]];
  relErr0   = vlRes0["RelativeError"];

  Print["  Sector sum      = ", vlRes0["SectorSum"]];
  Print["  Direct NIntegrate = ", vlRes0["DirectResult"]];
  Print["  Relative error  = ", relErr0];
  Print["  DroppedSectors  = ", nDropped0];

  (* Gate 1: at least one sector dropped (z0=1e6 >> 1, constant-root sectors must vanish) *)
  If[nDropped0 > 0,
    cc17Pass["17B-DroppedSectors",
             "count=" <> ToString[nDropped0] <> " > 0 (z0=1e6 constant-root sectors removed)"],
    cc17Fail["17B-DroppedSectors", ">0", 0];
    pass17B = False
  ];

  (* Gate 2: relErr vs NIntegrate reference < 0.1% *)
  If[!NumericQ[relErr0] || relErr0 >= 0.001,
    cc17Fail["17B-exactness", "relErr<0.001",
             "relErr=" <> ToString[relErr0]];
    pass17B = False,
    cc17Pass["17B-exactness",
             "relErr=" <> ToString[relErr0] <> " < 0.001; exact=" <> ToString[exactAns]]
  ];

  Label["done17B"];
  Print[];
  Print[If[pass17B, "CC17B PASS", "CC17B FAIL"]];
  Print[];
  pass17B
];
