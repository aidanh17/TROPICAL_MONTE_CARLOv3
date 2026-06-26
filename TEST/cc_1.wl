(* ============================================================================
   TEST/cc_1.wl  —  Cross-check #1 (plan.md §8.2)
   Tropical sector sum vs direct NIntegrate (Tier 1: WL + g++ only).

   For each test case:
     - Build the tropical fan via tropical_fan.wl / tropical_eval.wl
     - Call ValidateDecomposition to get:
         sectorSum   = Σ_s ∫_{[0,1]^n} sector integrand
         directResult = NIntegrate of original integrand over [0,∞)^n
     - PASS if relErr = |sectorSum - directResult| / |directResult| < tol
     - Print  "CC1 PASS …" / "CC1 FAIL expected=<e> got=<g>"

   Ported from:
     OLD_CODE/TROPICAL_MONTE_CARLO/EXAMPLES/tropical_eval_examples.wl
     Example 1 (2D, real exponent), Example 2 (2D, complex exponent),
     Example 3 (2D, kinematic parameter), plus a 3D case analogous to
     old test_coverage.wl.

   Run:
     wolframscript -file TEST/cc_1.wl
   from the TROPICAL_MONTE_CARLO3 directory, or load directly in WL with
   $InputFileName set appropriately.
   ============================================================================ *)

(* --- locate root and load v3 package via absolute paths --- *)
With[{root = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]]},
  Get[FileNameJoin[{root, "tropical_fan.wl"}]];
  Get[FileNameJoin[{root, "tropical_eval.wl"}]];
];

Print["CC1: tropical_eval.wl loaded."];
Print[];

(* ============================================================================
   Helper: run one sub-case and print PASS / FAIL.
   relTol  — maximum acceptable relative error
   ============================================================================ *)
cc1Check[label_String, spec_Association, kinRules_List, relTol_Real] :=
Module[{verts, fanData, vr, relErr, pass},

  (* Build fan using the polynomial with real exponent for geometry *)
  verts   = PolytopeVertices[
    (Times @@ MapThread[Power, {spec["Polynomials"], Re /@ spec["PolynomialExponents"]}]),
    spec["Variables"]
  ];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  If[fanData === $Failed || !ListQ[fanData] || Length[fanData] < 2,
    Print["CC1 FAIL  [", label, "]  fan computation returned $Failed"];
    Return[$Failed];
  ];

  (* ValidateDecomposition computes sectorSum + directResult *)
  vr = Quiet @ ValidateDecomposition[spec, fanData, kinRules, 4];

  If[!AssociationQ[vr],
    Print["CC1 FAIL  [", label, "]  ValidateDecomposition returned non-Association: ", vr];
    Return[$Failed];
  ];

  relErr = vr["RelativeError"];

  If[!NumericQ[relErr],
    Print["CC1 FAIL  [", label, "]  RelativeError is not numeric: ", relErr];
    Return[$Failed];
  ];

  pass = TrueQ[Abs[relErr] < relTol];

  If[pass,
    Print["CC1 PASS  [", label, "]  relErr=", ScientificForm[Abs[relErr], 3],
          "  directResult=", vr["DirectResult"],
          "  sectorSum=", vr["SectorSum"],
          "  sectors=", Length[fanData[[2]]]],
    Print["CC1 FAIL  [", label, "]  expected=relErr<", relTol,
          "  got=", ScientificForm[Abs[relErr], 3],
          "  directResult=", vr["DirectResult"],
          "  sectorSum=", vr["SectorSum"]]
  ];

  pass
];


(* ============================================================================
   Case 1a — 2D real-exponent integral (ported from OLD Example 1)
   ∫_{[0,∞)^2} dx1 dx2 / (1 + x1^2 + x2^2)^3
   ============================================================================ *)
Print["--- Case 1a: 2D basic convergent, (1+x1^2+x2^2)^{-3} ---"];
cc1Check[
  "1a_2D_real",
  <|
    "Polynomials"         -> {1 + x[1]^2 + x[2]^2},
    "MonomialExponents"   -> {0, 0},
    "PolynomialExponents" -> {-3},
    "Variables"           -> {x[1], x[2]},
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> None
  |>,
  {},   (* no kinematic substitution *)
  0.001 (* relTol < 1e-3 per plan §8.2 #1 *)
];
Print[];


(* ============================================================================
   Case 1b — 2D complex polynomial exponent (ported from OLD Example 2)
   ∫_{[0,∞)^2} dx1 dx2 / (1 + 2x1^2 + x2^2 + x1*x2^2 + 3x1^2*x2)^{2+I}
   Demonstrates F4 (complex polynomial exponent) going through the fan.
   ============================================================================ *)
Print["--- Case 1b: 2D complex exponent, P^{-(2+I)} ---"];
cc1Check[
  "1b_2D_complex_exp",
  <|
    "Polynomials"         -> {1 + 2 x[1]^2 + x[2]^2 + x[1] x[2]^2 + 3 x[1]^2 x[2]},
    "MonomialExponents"   -> {0, 0},
    "PolynomialExponents" -> {-(2 + I)},
    "Variables"           -> {x[1], x[2]},
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> None
  |>,
  {},
  0.001
];
Print[];


(* ============================================================================
   Case 1c — 2D kinematic-dependent integral at a fixed kinematic value
              (ported from OLD Example 3, evaluated at lam=2)
   ∫_{[0,∞)^2} dx1 dx2 / (1 + lam*x1^2 + x2^2 + x1*x2^2)^2  at lam=2
   Demonstrates F1 with a KinematicSymbols entry, cross-checked at lam=2.
   ============================================================================ *)
Print["--- Case 1c: 2D kinematic param lam, evaluated at lam=2 ---"];
Block[{lam},
  cc1Check[
    "1c_2D_kinematic_lam2",
    <|
      "Polynomials"         -> {1 + lam x[1]^2 + x[2]^2 + x[1] x[2]^2},
      "MonomialExponents"   -> {0, 0},
      "PolynomialExponents" -> {-2},
      "Variables"           -> {x[1], x[2]},
      "KinematicSymbols"    -> {lam},
      "RegulatorSymbol"     -> None
    |>,
    {lam -> 2},  (* substitute lam=2 for the NIntegrate cross-check *)
    0.001
  ];
];
Print[];


(* ============================================================================
   Case 1d — 3D real-exponent integral (higher-dimensional smoke test)
   ∫_{[0,∞)^3} dx1 dx2 dx3 / (1 + x1^2 + x2^2 + x3^2)^4
   Analogue of old test_coverage.wl dimension extension test.
   Tolerance relaxed to 1e-2 for n=3 per plan §8.2 #1 note.
   ============================================================================ *)
Print["--- Case 1d: 3D basic convergent, (1+x1^2+x2^2+x3^2)^{-4} ---"];
cc1Check[
  "1d_3D_real",
  <|
    "Polynomials"         -> {1 + x[1]^2 + x[2]^2 + x[3]^2},
    "MonomialExponents"   -> {0, 0, 0},
    "PolynomialExponents" -> {-4},
    "Variables"           -> {x[1], x[2], x[3]},
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> None
  |>,
  {},
  0.01   (* plan: "tighter for n≤4"; 1e-2 for n=3 *)
];
Print[];


(* ============================================================================
   Case 1e — 2D two-polynomial product (numerator + denominator)
   ∫_{[0,∞)^2} dx1 dx2  (1 + x1 x2)^{1/2} / (1 + x1^2 + x2^2)^2
   Exercises the B>0 (numerator factor) code path (plan §8.1 "numerator factors").
   ============================================================================ *)
Print["--- Case 1e: 2D two polynomials, numerator B=+1/2 ---"];
cc1Check[
  "1e_2D_numerator",
  <|
    "Polynomials"         -> {1 + x[1] x[2], 1 + x[1]^2 + x[2]^2},
    "MonomialExponents"   -> {0, 0},
    "PolynomialExponents" -> {1/2, -2},
    "Variables"           -> {x[1], x[2]},
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> None
  |>,
  {},
  0.001
];
Print[];


Print["CC1: all cases complete."];
