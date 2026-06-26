(* ============================================================================
   TEST/cc_30.wl  --  Cross-check #30 (plan.md §8.2)
   Validation utilities exercised as assertions (Tier 1: WL + g++ only).

   PASS criterion (plan.md §8.2, row 30):
     ValidateDecomposition / ValidateSubtraction / ValidateIBP /
     ValidateLiftedDecomposition all PASS on known inputs.

   Four sub-cases, one per utility:

   30A  ValidateDecomposition:
        P = (1+x[1]^2+x[2]^2)^{-3}, no regulator.
        Exact closed form: pi^2/24 (Dirichlet / Gamma products).
        Ported from cc_1.wl case 1a (which passes in the suite).
        PASS: relErr(SectorSum, NIntegrate) < 1e-3.

   30B  ValidateSubtraction:
        x[1]^{2eps-1} * (1+x[1]+x[2]+x[1]*x[2])^{-2}, 2D, eps regulator.
        Ported from OLD_CODE Test 8 (tropical_eval_examples.wl).
        Run on all divergent sectors; each must give relErr < 0.02.
        PASS: all divergent sectors give AssociationQ result with
              RelativeError < 0.02.

   30C  ValidateIBP:
        Same 2D integrand as 30B.
        Ported from OLD_CODE Test 20 (tropical_eval_examples.wl).
        PASS: all divergent sectors give AssociationQ result with
              RelativeError < 0.05.

   30D  ValidateLiftedDecomposition:
        P = (1+10^6*x[1])^{-2}, n=1, k=1 lift.  Exact = 10^{-6}.
        Ported from cc_17B / OLD_CODE test_lifted.wl Test 23B.
        PASS: RelativeError < 1e-3.

   Run:
     wolframscript -file TEST/cc_30.wl
   from the TROPICAL_MONTE_CARLO3 root.
   ============================================================================ *)

(* --- locate root and load v3 packages via absolute paths --- *)
With[{root = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]]},
  Get[FileNameJoin[{root, "tropical_fan.wl"}]];
  Get[FileNameJoin[{root, "tropical_eval.wl"}]];
];

Print["CC30: tropical_fan.wl + tropical_eval.wl loaded."];
Print[];

(* ============================================================================
   Helpers
   ============================================================================ *)
cc30Pass[label_String, info_String] :=
  Print["CC30 PASS [", label, "] ", info];

cc30Fail[label_String, expected_, got_] :=
  Print["CC30 FAIL [", label, "] expected=", expected, " got=", got];


(* ============================================================================
   30A -- ValidateDecomposition on a known convergent integral
   Int_0^inf Int_0^inf (1+x1^2+x2^2)^{-3} dx1 dx2
   Ported from cc_1.wl case 1a; ValidateDecomposition must return
   relErr < 1e-3 (same tolerance as #1).
   ============================================================================ *)
Module[
  {spec, verts, fanData, vr, relErr, pass30A},

  Print["--- CC30A: ValidateDecomposition (2D, (1+x1^2+x2^2)^{-3}) ---"];
  pass30A = True;

  spec = <|
    "Polynomials"         -> {1 + x[1]^2 + x[2]^2},
    "MonomialExponents"   -> {0, 0},
    "PolynomialExponents" -> {-3},
    "Variables"           -> {x[1], x[2]},
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> None
  |>;

  verts = PolytopeVertices[
    (1 + x[1]^2 + x[2]^2)^(-3), {x[1], x[2]}
  ];
  fanData = Quiet @ ComputeDecomposition[verts, "ShowProgress" -> False];

  If[fanData === $Failed || !ListQ[fanData] || Length[fanData] < 2,
    cc30Fail["30A-fan", "valid fanData", fanData];
    pass30A = False;
    Goto["done30A"]
  ];

  Print["  Fan: ", Length[fanData[[1]]], " rays, ",
        Length[fanData[[2]]], " sectors"];

  vr = Quiet @ ValidateDecomposition[spec, fanData, {}, 4];

  If[!AssociationQ[vr],
    cc30Fail["30A-ValidateDecomposition", "Association", vr];
    pass30A = False;
    Goto["done30A"]
  ];

  relErr = vr["RelativeError"];

  If[!NumericQ[relErr],
    cc30Fail["30A-relErr-numeric", "NumericQ relErr", relErr];
    pass30A = False;
    Goto["done30A"]
  ];

  Print["  SectorSum     = ", vr["SectorSum"]];
  Print["  DirectResult  = ", vr["DirectResult"]];
  Print["  RelativeError = ", relErr];

  If[Abs[relErr] >= 1*^-3,
    cc30Fail["30A-relErr", "relErr<1e-3",
             "relErr=" <> ToString[CForm[relErr]]];
    pass30A = False,
    cc30Pass["30A-relErr",
             "relErr=" <> ToString[CForm[relErr]] <> " < 1e-3"]
  ];

  Label["done30A"];
  Print[];
  Print[If[pass30A, "CC30A PASS", "CC30A FAIL"]];
  Print[];
  pass30A
];


(* ============================================================================
   30B -- ValidateSubtraction on 2D eps-divergent integrand
   x[1]^{2eps-1} (1+x1+x2+x1*x2)^{-2}
   Ported from OLD_CODE Test 8 (tropical_eval_examples.wl).
   For each divergent sector: call ProcessDivergentSector then
   ValidateSubtraction at testEps=0.05.
   PASS: all sectors give AssociationQ with RelativeError < 0.02.
   ============================================================================ *)
Module[
  {eps, poly, vars, spec, verts, fanData, dualVerts, simplexList,
   allSectors, divSectors, allPass30B, nDiv},

  Print["--- CC30B: ValidateSubtraction (2D, x1^{2eps-1}(1+x1+x2+x1x2)^{-2}) ---"];
  allPass30B = True;

  eps  = Symbol["eps30B"];
  poly = 1 + x[1] + x[2] + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"         -> {poly},
    "MonomialExponents"   -> {2 eps - 1, 0},
    "PolynomialExponents" -> {-2},
    "Variables"           -> vars,
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> eps
  |>;

  verts   = PolytopeVertices[poly^(-2), vars];
  fanData = Quiet @ ComputeDecomposition[verts, "ShowProgress" -> False];

  If[fanData === $Failed || !ListQ[fanData] || Length[fanData] < 2,
    cc30Fail["30B-fan", "valid fanData", fanData];
    allPass30B = False;
    Goto["done30B"]
  ];

  dualVerts   = fanData[[1]];
  simplexList = fanData[[2]];

  allSectors = Table[
    Quiet @ ProcessSector[spec, dualVerts, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  divSectors = Select[allSectors, AssociationQ[#] && #["IsDivergent"] &];
  nDiv = Length[divSectors];

  If[nDiv == 0,
    cc30Fail["30B-divCount", ">0", 0];
    allPass30B = False;
    Goto["done30B"]
  ];

  Print["  Found ", nDiv, " divergent sector(s)"];

  Do[
    Module[{sd, divData, vsResult, relErr, label},
      sd    = divSectors[[idx]];
      label = "30B-sector" <> ToString[sd["ConeIndex"]];

      divData = Quiet @ ProcessDivergentSector[sd, spec];
      If[!AssociationQ[divData],
        cc30Fail[label <> "-ProcessDivergentSector", "Association", divData];
        allPass30B = False;
        Return[]
      ];

      vsResult = Quiet @ ValidateSubtraction[divData, sd, spec, {}, 0.05];
      If[!AssociationQ[vsResult],
        cc30Fail[label <> "-ValidateSubtraction", "Association", vsResult];
        allPass30B = False;
        Return[]
      ];

      relErr = vsResult["RelativeError"];
      Print["  Sector ", sd["ConeIndex"], ": RelativeError = ", relErr];

      If[!NumericQ[relErr],
        cc30Fail[label <> "-numeric", "NumericQ", relErr];
        allPass30B = False;
        Return[]
      ];

      If[relErr >= 0.02,
        cc30Fail[label <> "-relErr", "relErr<0.02",
                 "relErr=" <> ToString[CForm[relErr]]];
        allPass30B = False,
        cc30Pass[label <> "-relErr",
                 "relErr=" <> ToString[CForm[relErr]] <> " < 0.02"]
      ]
    ],
    {idx, nDiv}
  ];

  Label["done30B"];
  Print[];
  Print[If[allPass30B, "CC30B PASS", "CC30B FAIL"]];
  Print[];
  allPass30B
];


(* ============================================================================
   30C -- ValidateIBP on the same 2D eps-divergent integrand
   x[1]^{2eps-1} (1+x1+x2+x1*x2)^{-2}
   Ported from OLD_CODE Test 20 (tropical_eval_examples.wl).
   PASS: all divergent sectors give AssociationQ with RelativeError < 0.05.
   ============================================================================ *)
Module[
  {eps, poly, vars, spec, verts, fanData, dualVerts, simplexList,
   allSectors, divSectors, allPass30C, nDiv},

  Print["--- CC30C: ValidateIBP (2D, x1^{2eps-1}(1+x1+x2+x1x2)^{-2}) ---"];
  allPass30C = True;

  eps  = Symbol["eps30C"];
  poly = 1 + x[1] + x[2] + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"         -> {poly},
    "MonomialExponents"   -> {2 eps - 1, 0},
    "PolynomialExponents" -> {-2},
    "Variables"           -> vars,
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> eps
  |>;

  verts   = PolytopeVertices[poly^(-2), vars];
  fanData = Quiet @ ComputeDecomposition[verts, "ShowProgress" -> False];

  If[fanData === $Failed || !ListQ[fanData] || Length[fanData] < 2,
    cc30Fail["30C-fan", "valid fanData", fanData];
    allPass30C = False;
    Goto["done30C"]
  ];

  dualVerts   = fanData[[1]];
  simplexList = fanData[[2]];

  allSectors = Table[
    Quiet @ ProcessSector[spec, dualVerts, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  divSectors = Select[allSectors, AssociationQ[#] && #["IsDivergent"] &];
  nDiv = Length[divSectors];

  If[nDiv == 0,
    cc30Fail["30C-divCount", ">0", 0];
    allPass30C = False;
    Goto["done30C"]
  ];

  Print["  Found ", nDiv, " divergent sector(s)"];

  Do[
    Module[{sd, ibpSD, vResult, relErr, label},
      sd    = divSectors[[idx]];
      label = "30C-sector" <> ToString[sd["ConeIndex"]];

      ibpSD = Quiet @ IBPProcessSector[sd, spec];
      If[!AssociationQ[ibpSD],
        cc30Fail[label <> "-IBPProcessSector", "Association", ibpSD];
        allPass30C = False;
        Return[]
      ];

      vResult = Quiet @ ValidateIBP[ibpSD, sd, spec, {}, 0.05];
      If[!AssociationQ[vResult],
        cc30Fail[label <> "-ValidateIBP", "Association", vResult];
        allPass30C = False;
        Return[]
      ];

      relErr = vResult["RelativeError"];
      Print["  Sector ", sd["ConeIndex"], ": RelativeError = ", relErr];

      If[!NumericQ[relErr],
        cc30Fail[label <> "-numeric", "NumericQ", relErr];
        allPass30C = False;
        Return[]
      ];

      If[relErr >= 0.05,
        cc30Fail[label <> "-relErr", "relErr<0.05",
                 "relErr=" <> ToString[CForm[relErr]]];
        allPass30C = False,
        cc30Pass[label <> "-relErr",
                 "relErr=" <> ToString[CForm[relErr]] <> " < 0.05"]
      ]
    ],
    {idx, nDiv}
  ];

  Label["done30C"];
  Print[];
  Print[If[allPass30C, "CC30C PASS", "CC30C FAIL"]];
  Print[];
  allPass30C
];


(* ============================================================================
   30D -- ValidateLiftedDecomposition on 1D large-coefficient integrand
   P = (1 + 10^6 * x[1])^{-2}, n=1, k=1 lift.  Exact = 10^{-6}.
   Ported from cc_17B / OLD_CODE test_lifted.wl / sandbox_toy0_1d.wl.
   Explicit 2-D fan for the 1D+1 lifted support {(0,0),(1,1)}.
   PASS: AssociationQ; RelativeError < 1e-3.
   ============================================================================ *)
Module[
  {spec0, liftRules0, lcRes, liftedSpec, liftData,
   raysK1, sectsK1, explicitFan,
   vlRes, relErr, nDropped, pass30D},

  Print["--- CC30D: ValidateLiftedDecomposition (1D, k=1, exact=10^{-6}) ---"];
  pass30D = True;

  spec0 = <|
    "Polynomials"         -> {1 + 10^6 * x[1]},
    "MonomialExponents"   -> {0},
    "PolynomialExponents" -> {-2},
    "Variables"           -> {x[1]},
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> None
  |>;

  liftRules0 = {<|"PolyIndex" -> 1, "ExponentVector" -> {1}, "k" -> 1|>};

  lcRes = LiftCoefficients[spec0, liftRules0];

  If[!AssociationQ[lcRes],
    cc30Fail["30D-LiftCoefficients", "Association", lcRes];
    pass30D = False;
    Goto["done30D"]
  ];

  liftedSpec = lcRes["LiftedSpec"];
  liftData   = lcRes["LiftData"];

  Print["  z0 = ", liftData["z0"]];

  (* Lift identity check: z -> z0 must recover the original polynomial *)
  Module[{lp = liftedSpec["Polynomials"][[1]],
          z0v = liftData["z0"],
          auxV = liftData["AuxVariable"]},
    If[Expand[lp /. auxV -> z0v] === Expand[spec0["Polynomials"][[1]]],
      Print["  Lift identity: PASS"],
      Print["  Lift identity: FAIL -- z->z0 substitution does not recover original"];
      pass30D = False
    ]
  ];

  (* Explicit complete unimodular triangulation of the 2-D Newton polytope
     for support {(0,0),(1,1)} (k=1 lifted 1D polytope).
     7 maximal cones; from sandbox_toy0_1d.wl / cc_17B. *)
  raysK1  = {{1,0},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1}};
  sectsK1 = Table[{i, If[i < Length[raysK1], i+1, 1]}, {i, Length[raysK1]}];
  explicitFan = {raysK1, sectsK1};

  Print["  Fan: ", Length[raysK1], " rays, ",
        Length[sectsK1], " sectors (explicit, from sandbox_toy0_1d)"];

  vlRes = Quiet @ ValidateLiftedDecomposition[
    spec0, liftedSpec, explicitFan, liftData, {}, 4
  ];

  If[!AssociationQ[vlRes],
    cc30Fail["30D-ValidateLiftedDecomposition", "Association", vlRes];
    pass30D = False;
    Goto["done30D"]
  ];

  nDropped = Length[Lookup[vlRes, "DroppedSectors", {}]];
  relErr   = Lookup[vlRes, "RelativeError", Missing[]];

  Print["  SectorSum      = ", Lookup[vlRes, "SectorSum", Missing[]]];
  Print["  DirectResult   = ", Lookup[vlRes, "DirectResult", Missing[]]];
  Print["  RelativeError  = ", relErr];
  Print["  DroppedSectors = ", nDropped, " (z0=10^6; constant-root sectors removed)"];

  If[!NumericQ[relErr],
    cc30Fail["30D-numeric", "NumericQ relErr", relErr];
    pass30D = False;
    Goto["done30D"]
  ];

  If[relErr >= 1*^-3,
    cc30Fail["30D-relErr", "relErr<1e-3",
             "relErr=" <> ToString[CForm[relErr]]];
    pass30D = False,
    cc30Pass["30D-relErr",
             "relErr=" <> ToString[CForm[relErr]] <> " < 1e-3"]
  ];

  Label["done30D"];
  Print[];
  Print[If[pass30D, "CC30D PASS", "CC30D FAIL"]];
  Print[];
  pass30D
];
