(* ============================================================================
   EXAMPLES/Ex-V3d.wl  —  Worked Example V3d
   plan.md §8.4: "Numerator (...)^{+3/2} x divergent denominator"
   Cross-check: #41

   Integrand
   ---------
   I(eps) = Int_{[0,inf)^2} x1^{-1+eps} (1+x1)^{3/2} (1+x1+x2)^{-4}  dx1 dx2

   The factor (1+x1)^{3/2} has a POSITIVE polynomial exponent (B>0: numerator),
   while x1^{-1+eps} and (1+x1+x2)^{-4} make the integral divergent as eps->0.
   This is the first demonstration that IBP handles B>0 alongside a 1/eps pole.

   References
   ----------
   (i)  NIntegrate at fixed eps* = 0.02 (no analytic closed form for the full
        Laurent; the convergent sub-check below provides an analytic anchor).
   (ii) Exact analytic for the CONVERGENT sub-check (eps=0, convergent spec):
        I_conv = Int_{[0,inf)^2} (1+x1)^{3/2} (1+x1+x2)^{-4}  dx1 dx2  =  2/3
        (established in OLD_CODE/TEST/test_coverage.wl entry B5 and
         DIVERGENT_VALIDATION.md B5).

   Two routes for (pole, finite)
   ------------------------------
   Route A — EvaluateTropicalMCIBP (IBP @ eps=0, MC sampler)
   Route B — LaurentFromSubtraction (finite-eps fit, MC sampler)

   PASS criteria: all D1-D5 must hold.  Tier 1: WL + g++ only.
     D1  Convergent sub-check rel-err vs 2/3 < 1e-2
     D2  IBP pole is non-zero
     D3  relErr F_IBP at eps-star vs NIntRef < 3e-2
     D4  relErr F_Sub at eps-star vs NIntRef < 3e-2
     D5  pole-A minus pole-B < 0.1

   Run
   ---
     cd /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3
     wolframscript -file EXAMPLES/Ex-V3d.wl
   exits 0 on PASS, 1 on FAIL.
   ============================================================================ *)

(* ── 0.  Locate packages ──────────────────────────────────────────────────── *)

$v3Root = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
$fanWL  = FileNameJoin[{$v3Root, "tropical_fan.wl"}];
$evalWL = FileNameJoin[{$v3Root, "tropical_eval.wl"}];

If[!FileExistsQ[$fanWL],
  Print["Ex-V3d FAIL  tropical_fan.wl not found at ", $fanWL]; Quit[1]];
If[!FileExistsQ[$evalWL],
  Print["Ex-V3d FAIL  tropical_eval.wl not found at ", $evalWL]; Quit[1]];

Get[$fanWL];
Get[$evalWL];

Quiet[CreateDirectory[
  FileNameJoin[{$v3Root, "INTERFILES", "ExV3d"}],
  CreateIntermediateDirectories -> True],
 {CreateDirectory::eexist}];

$ioDir = FileNameJoin[{$v3Root, "INTERFILES", "ExV3d"}];

Print[""];
Print["================================================================"];
Print["  EXAMPLES/Ex-V3d: Numerator (...)^{+3/2} x Divergent Denom   "];
Print["  Cross-check #41 (plan.md §8.4)                              "];
Print["================================================================"];
Print[];

(* ── 1.  Helpers ──────────────────────────────────────────────────────────── *)

SetAttributes[quietRun, HoldFirst];
quietRun[expr_] := Quiet[expr,
  {TropicalEval`EvaluateTropicalMCIBP::validate,
   TropicalEval`LaurentFromSubtraction::validate,
   TropicalEval`ValidateDecomposition::validate,
   General::stop}];

sci[x_] := ToString[ScientificForm[N[x], 3]];

$passes   = {};
$failures = {};

record[name_String, ok_, passDetail_String : "", failExtra_String : ""] :=
  If[TrueQ[ok],
    (AppendTo[$passes, name];
     Print["  Ex-V3d PASS  ", name,
           If[passDetail =!= "", "  " <> passDetail, ""]]),
    (AppendTo[$failures, name];
     Print["  Ex-V3d FAIL  ", name,
           If[failExtra =!= "", "  " <> failExtra, ""]])];

(* ── 2.  Integrand specification ─────────────────────────────────────────── *)

Print["--- §1  Integrand specification ---"];
eps   = Symbol["epsV3d"];
vars  = {x[1], x[2]};
pNum  = 1 + x[1];          (* NUMERATOR factor: B>0 = 3/2 *)
pDen  = 1 + x[1] + x[2];   (* DENOMINATOR factor: B<0 = -4 *)

spec = <|
  "Polynomials"         -> {pNum, pDen},
  "MonomialExponents"   -> {-1 + eps, 0},   (* x1^{-1+eps}: divergence source *)
  "PolynomialExponents" -> {3/2, -4},       (* numerator 3/2; denominator -4 *)
  "Variables"           -> vars,
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> eps
|>;

Print["  Polynomials: pNum = ", pNum, "^{3/2}   pDen = ", pDen, "^{-4}"];
Print["  MonomialExponents: {-1+eps, 0}   (divergent as eps->0)"];
Print[];

(* ── 3.  Fan ──────────────────────────────────────────────────────────────── *)

Print["--- §2  Tropical fan ---"];
verts = PolytopeVertices[(pNum * pDen)^(-1), vars];
fan   = ComputeDecomposition[verts, "ShowProgress" -> False];

If[!ListQ[fan] || Length[fan] < 2,
  Print["Ex-V3d FAIL  fan computation returned $Failed"];
  Quit[1]];
Print["  Sectors: ", Length[fan[[2]]]];
Print[];

(* ── 4.  Criterion D1 — convergent sub-check (exact B5 = 2/3) ─────────────── *)

Print["--- §3  D1: Convergent sub-check: (1+x1)^{3/2} (1+x1+x2)^{-4} = 2/3 ---"];
specConv = <|
  "Polynomials"         -> {pNum, pDen},
  "MonomialExponents"   -> {0, 0},    (* no eps: convergent *)
  "PolynomialExponents" -> {3/2, -4},
  "Variables"           -> vars,
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> None
|>;

fanConv = ComputeDecomposition[
  PolytopeVertices[(pNum * pDen)^(-1), vars],
  "ShowProgress" -> False];

vrConv = quietRun @ ValidateDecomposition[specConv, fanConv, {}, 4];

If[AssociationQ[vrConv],
  Print["  Sector sum = ", vrConv["SectorSum"],
        "  exact = 2/3 = ", N[2/3],
        "  rel-err = ", sci @ vrConv["RelativeError"]];
  record["D1 convergent B5 numerator vs exact 2/3",
    vrConv["RelativeError"] < 1*^-2,
    "rel-err=" <> sci[vrConv["RelativeError"]],
    "expected=<1e-2 got=" <> sci[vrConv["RelativeError"]]],
  Print["  ValidateDecomposition returned: ", vrConv];
  AppendTo[$failures, "D1 ValidateDecomposition"];
];
Print[];

(* ── 5.  NIntegrate reference at fixed eps* ──────────────────────────────── *)

Print["--- §4  NIntegrate reference at eps* = 0.02 ---"];
epsStar  = 0.02;

(* Compactify: xi = ti/(1-ti) -> dxi = 1/(1-ti)^2 dti, ti in [0,1] *)
nintRef = Quiet @ NIntegrate[
  (u1/(1-u1))^(-1 + epsStar)
    * (1 + u1/(1-u1))^(3/2)
    * (1 + u1/(1-u1) + u2/(1-u2))^(-4)
    / ((1-u1)^2 * (1-u2)^2),
  {u1, 0, 1}, {u2, 0, 1},
  Method -> "GlobalAdaptive", PrecisionGoal -> 5, MaxRecursion -> 40];
Print["  NIntegrate(eps*=", epsStar, ") = ", N[nintRef]];
Print[];

(* ── 6.  Route A — EvaluateTropicalMCIBP (MC) ────────────────────────────── *)

Print["--- §5  Route A: IBP @ eps=0 (EvaluateTropicalMCIBP, MC) ---"];
resA = quietRun @ EvaluateTropicalMCIBP[spec, fan, {{}},
  "Integrator"       -> "MC",
  "NSamples"         -> 600000,
  "RunChecks"        -> False,
  "Verbose"          -> False,
  "WorkingDirectory" -> FileNameJoin[{$ioDir, "routeA"}]];

If[!AssociationQ[resA] || !KeyExistsQ[resA, "Results"] ||
   Length[resA["Results"]] == 0,
  Print["  EvaluateTropicalMCIBP returned: ", resA];
  AppendTo[$failures, "D2/D3 IBP driver"];
  AppendTo[$failures, "D5 IBP driver"];
  Goto[$skipRouteA]
];

poleA   = Re @ resA["Results"][[1]]["PoleCoefficient"];
finA    = Re @ resA["Results"][[1]]["FinitePart"];
fIBPeps = poleA / epsStar + finA;    (* reconstructed F at epsStar from Laurent *)

Print["  IBP-MC:  pole = ", poleA, "  finite = ", finA];
Print["  F_IBP(eps*=", epsStar, ") = pole/eps* + finite = ", fIBPeps,
      "  NInt = ", N[nintRef]];

(* D2: pole is non-zero *)
record["D2 IBP pole is non-zero (genuine 1/eps)",
  NumericQ[poleA] && Abs[poleA] > 0.01,
  "pole=" <> sci[poleA],
  "pole=" <> sci[poleA] <> " (expected |pole|>0.01)"];

(* D3: F_IBP at epsStar vs NIntegrate *)
If[NumericQ[nintRef] && NumericQ[fIBPeps] && Abs[nintRef] > 0,
  Module[{relE = Abs[(fIBPeps - nintRef) / nintRef]},
    Print["  rel-err(F_IBP vs NInt) = ", sci[relE], "  (gate < 3e-2)"];
    record["D3 F_IBP(eps*) vs NIntegrate",
      relE < 3*^-2,
      "rel-err=" <> sci[relE],
      "expected=<3e-2 got=" <> sci[relE]]],
  record["D3 F_IBP(eps*) vs NIntegrate",
    False, "",
    "nintRef or fIBPeps non-numeric: " <> ToString[{N[nintRef], fIBPeps}]]
];

Label[$skipRouteA];
Print[];

(* ── 7.  Route B — LaurentFromSubtraction ────────────────────────────────── *)

Print["--- §6  Route B: LaurentFromSubtraction (finite-eps fit, MC) ---"];
epsVals = {0.01, 0.02, 0.04};
resB = quietRun @ LaurentFromSubtraction[spec, fan, {{}},
  "Integrator"       -> "MC",
  "NSamples"         -> 400000,
  "EpsilonValues"    -> epsVals,
  "WorkingDirectory" -> FileNameJoin[{$ioDir, "routeB"}]];

If[!AssociationQ[resB] || !KeyExistsQ[resB, "Results"] ||
   Length[resB["Results"]] == 0,
  Print["  LaurentFromSubtraction returned: ", resB];
  AppendTo[$failures, "D4/D5 Sub driver"];
  Goto[$skipRouteB]
];

poleB = Re @ resB["Results"][[1]]["Pole"];
finB  = Re @ resB["Results"][[1]]["Finite"];
(* eps* = 0.02 is the second element of epsVals;
   LaurentFromSubtraction stores the per-eps integrals in "FullIntegrals" *)
$fullInts = Quiet[resB["Results"][[1]]["FullIntegrals"]];
fSubeps   = If[ListQ[$fullInts] && Length[$fullInts] >= 2,
  Re[$fullInts[[2]]],      (* index 2 = eps=0.02 *)
  poleB / epsStar + finB   (* fallback: reconstruct from Laurent fit *)
];

Print["  Fit-MC: pole = ", poleB, "  finite = ", finB];
Print["  F_Sub(eps*=", epsStar, ") = ", fSubeps,
      "  NInt = ", N[nintRef]];

(* D4: F_Sub at epsStar vs NIntegrate *)
If[NumericQ[nintRef] && NumericQ[fSubeps] && Abs[nintRef] > 0,
  Module[{relE = Abs[(fSubeps - nintRef) / nintRef]},
    Print["  rel-err(F_Sub vs NInt) = ", sci[relE], "  (gate < 3e-2)"];
    record["D4 F_Sub(eps*) vs NIntegrate",
      relE < 3*^-2,
      "rel-err=" <> sci[relE],
      "expected=<3e-2 got=" <> sci[relE]]],
  record["D4 F_Sub(eps*) vs NIntegrate",
    False, "",
    "nintRef or fSubeps non-numeric"]
];

(* D5: IBP vs subtraction pole agreement (only if both routes succeeded) *)
If[ValueQ[poleA] && ValueQ[poleB],
  Module[{dpole = Abs[poleA - poleB]},
    Print["  |pole_A - pole_B| = ", sci[dpole], "  (gate < 0.1)"];
    record["D5 IBP vs Sub pole agreement  (#41 cross-check)",
      dpole < 0.1,
      "|dpole|=" <> sci[dpole],
      "expected=<0.1 got=" <> sci[dpole]]]
];

Label[$skipRouteB];
Print[];

(* ── 8.  Summary ──────────────────────────────────────────────────────────── *)

Print["================================================================"];
Print["  Results summary:"];
Print["    I(eps) = Int x1^{-1+eps} (1+x1)^{3/2} (1+x1+x2)^{-4} dx1 dx2"];
If[ValueQ[vrConv], Print["    convergent (eps=0) sub-check:  sector sum = ",
  vrConv["SectorSum"], "  exact = 2/3"]];
Print["    NIntegrate(eps*=", epsStar, ") = ", N[nintRef]];
If[ValueQ[poleA], Print["    IBP:  pole = ", poleA, "  finite = ", finA,
  "  F(eps*)~", fIBPeps]];
If[ValueQ[poleB], Print["    Fit:  pole = ", poleB, "  finite = ", finB,
  "  F(eps*)~", fSubeps]];
Print["  PASSED: ", Length[$passes]];
Print["  FAILED: ", Length[$failures]];
Print["================================================================"];

If[Length[$failures] == 0,
  Print["Ex-V3d PASS  numerator^{+3/2} x divergent denom:",
        " IBP and Fit both reproduce NIntegrate at eps*  (",
        Length[$passes], " criteria)"],
  Print["Ex-V3d FAIL  ", Length[$failures],
        " criterion(a) failed: ",
        StringRiffle[$failures, ", "]]];
Print["================================================================"];

Quit[If[Length[$failures] > 0, 1, 0]]
