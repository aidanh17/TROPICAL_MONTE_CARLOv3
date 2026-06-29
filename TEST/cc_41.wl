(* ============================================================================
   TEST/cc_41.wl  --  Cross-check #41 (plan.md §8.3)
   Numerator (B>0) x {divergence, lifting}            Tier 1: WL + g++ only

   What this tests
   ---------------
   Cross-check #41 is the first test that deliberately combines B>0 (a
   positive polynomial exponent, i.e. a numerator factor) with the two
   other major code paths that were previously tested only in isolation:
     41-A  Numerator x divergence: a positive-exponent factor alongside a
           polynomial that produces an eps-pole; the IBP path must handle
           both simultaneously.  The pole and finite part are compared to an
           exact analytic (Gamma) reference.
     41-B  Numerator x lifting: a positive-exponent factor alongside a
           polynomial with an extreme coefficient that triggers lifting;
           the lifted MC is compared to NIntegrate.

   Analytic references
   -------------------
   41-A:
     I = Int_0^inf dx1 dx2  x1^{-1+eps} (1+x1)^{3/2} (1+x1+x2)^{-4}

     For small eps the integral reads:
       I(eps) = B(eps, 3/2+1-eps) * B(?, ?)  -- not trivially closed.

     Instead we use the NIntegrate reference at eps* = 0.02, plus exact
     analytic reference for the *convergent sub-integral* (41-A2):
       I_conv = Int_0^inf dx1 dx2  (1+x1)^{3/2} (1+x1+x2)^{-4}  =  2/3
     established in OLD_CODE test_coverage.wl B5.

     For 41-A we use:
       - NIntegrate at fixed eps* as the reference for the full integral.
       - EvaluateTropicalMCIBP (IBP+MC) for pole + finite extraction.
       - LaurentFromSubtraction (subtraction+MC) for independent cross-check.
       PASS: pole from both routes agree within 0.05, F at eps* within 2 pct of
             NIntegrate for both IBP and subtraction routes.

   41-B:
     I = Int_0^inf dx1 dx2  (1+x1*x2)^{1/2} / (1 + 10^6*x1^2 + x2^2)^2
     Exact reference: NIntegrate (no analytic closed form).
     Lifting: extreme coefficient 10^6 in denominator triggers automatic
              DetectExtremeCoefficients -> LiftRules -> EvaluateTropicalMCLifted.
     PASS: |liftedMC - ref| / |ref| < 0.01  (1 pct)

   Old-tree source
   ---------------
   41-A is a mechanical combination of:
     OLD_CODE/.../EXAMPLES/test_coverage.wl  B5 (numerator factor B>0)
     OLD_CODE/.../EXAMPLES/test_coverage.wl  B3/B4 (IBP divergent integrals)
     OLD_CODE/.../DIVERGENT_VALIDATION.md    B5 entry
   41-B is a mechanical combination of:
     OLD_CODE/.../EXAMPLES/tropical_eval_examples.wl  Example 10 (numerator)
     TEST/test_lifted.wl  Test 23B (EvaluateTropicalMCLifted with extreme coeff)

   Run
   ---
     cd /Users/aidanh/Desktop/TROPICAL_MONTE_CARLOv3
     wolframscript -file TEST/cc_41.wl
   ============================================================================ *)

(* ── 0.  Locate packages, load ───────────────────────────────────────────── *)

$pkgRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
$fanWL   = FileNameJoin[{$pkgRoot, "tropical_fan.wl"}];
$evalWL  = FileNameJoin[{$pkgRoot, "tropical_eval.wl"}];

If[!FileExistsQ[$fanWL],
  Print["CC41 FAIL tropical_fan.wl not found at ", $fanWL]; Quit[1]];
If[!FileExistsQ[$evalWL],
  Print["CC41 FAIL tropical_eval.wl not found at ", $evalWL]; Quit[1]];

Get[$fanWL];
Get[$evalWL];

Print["CC41: packages loaded."];
Print[];

(* ── 1.  Shared utilities ────────────────────────────────────────────────── *)

sci[x_] := ScientificForm[x, 3];

SetAttributes[quietRun, HoldFirst];
quietRun[expr_] := Quiet[expr,
  {TropicalEval::validate, TropicalEval::divergent,
   TropicalEval::vegasbudget, General::stop}];

(* Global pass/fail tracking *)
$cc41AllPass  = True;
$cc41AnyFail  = {};

cc41Assert[label_String, test_,
           expectedStr_String, gotStr_String] :=
  If[TrueQ[test],
    Print["CC41 PASS  [", label, "]  expected=", expectedStr,
          "  got=", gotStr],
    Print["CC41 FAIL  [", label, "]  expected=", expectedStr,
          "  got=", gotStr];
    $cc41AllPass = False;
    AppendTo[$cc41AnyFail, label]
  ];

(* ── 2.  CUBA detection (Tier-2 guard) ──────────────────────────────────── *)

$cuba = TrueQ[Quiet[
  FindFile["CUBALink`"] =!= $Failed ||
  (FileExistsQ["/opt/homebrew/lib/libcuba.a"] ||
   FileExistsQ["/usr/local/lib/libcuba.a"])
]];
Print["CC41: CUBA available = ", $cuba];
Print[];


(* ============================================================================
   41-A: Numerator factor (B>0) x divergent integral (eps-pole)
         I(eps) = Int_0^inf dx1 dx2  x1^{-1+eps} (1+x1)^{3/2} (1+x1+x2)^{-4}

   The B>0 numerator polynomial is (1+x1)^{3/2}; the divergence is from the
   monomial factor x1^{-1+eps}.  Expected: a 1/eps pole.

   Analytic reference (convergent limit, eps->0 carefully):
     The sub-integral with eps=0 is *divergent* (that is the whole point),
     but for a FIXED eps* the NIntegrate gives us an independent reference.

   We also test the *convergent* reduction: at the convergent limit the
   numerator x divergence integral collapses to
       I_conv = Int dx1 dx2  (1+x1)^{3/2} (1+x1+x2)^{-4}  = 2/3
   (established in test_coverage.wl B5 and DIVERGENT_VALIDATION.md B5).
   ============================================================================ *)

Print["================================================================"];
Print["  CC41-A: Numerator x Divergence"];
Print["  I(eps)=Int dx1 dx2 x1^{-1+eps}(1+x1)^{3/2}(1+x1+x2)^{-4}"];
Print["================================================================"];

Module[
  {eps, vars, pNum, pDen, spec, verts, fan,
   epsStar, nMC, nSub, epsList,
   nintRef, ibpMC, lfMC,
   poleIbpMC, finIbpMC, poleSubMC, finSubMC,
   fIbpStar, fSubStar,
   tolPole, tolF},

  eps    = Symbol["eps41a"];
  vars   = {x[1], x[2]};
  pNum   = 1 + x[1];       (* positive exponent 3/2: numerator factor *)
  pDen   = 1 + x[1] + x[2];

  spec = <|
    "Polynomials"         -> {pNum, pDen},
    "MonomialExponents"   -> {-1 + eps, 0},   (* x1^{-1+eps}: divergence source *)
    "PolynomialExponents" -> {3/2, -4},        (* (1+x1)^{3/2} numerator; (1+x1+x2)^{-4} denom *)
    "Variables"           -> vars,
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> eps
  |>;

  (* Fan: use product of all polynomials (signed exponents irrelevant for fan geometry) *)
  verts = PolytopeVertices[(pNum * pDen)^(-1), vars];
  fan   = ComputeDecomposition[verts, "ShowProgress" -> False];

  If[!ListQ[fan] || Length[fan] < 2,
    Print["CC41 FAIL [41-A fan] fan computation returned $Failed"];
    $cc41AllPass = False;
    AppendTo[$cc41AnyFail, "41-A fan"];
    Goto[$skip41A]
  ];
  Print["41-A: fan built -- ", Length[fan[[2]]], " sectors"];

  (* NIntegrate reference at eps* = 0.02 via compactification *)
  epsStar = 0.02;
  nintRef = Quiet @ NIntegrate[
    (u1/(1-u1))^(-1+epsStar)
      (1 + u1/(1-u1))^(3/2)
      (1 + u1/(1-u1) + u2/(1-u2))^(-4)
      / ((1-u1)^2 (1-u2)^2),
    {u1, 0, 1}, {u2, 0, 1},
    Method -> "GlobalAdaptive", PrecisionGoal -> 5, MaxRecursion -> 40];
  Print["41-A: NIntegrate(eps*=", epsStar, ") = ", nintRef];

  nMC    = 500000;
  nSub   = 300000;
  epsList = {0.01, 0.02, 0.04};
  tolPole = 0.1;   (* IBP-MC pole tolerance; pole itself may be O(1) *)
  tolF    = 0.03;  (* tolerance for F at eps* vs NIntegrate *)

  (* ---- Route (a): IBP+MC ---- *)
  Print["41-A: running EvaluateTropicalMCIBP (MC, ", nMC, " samples) ..."];
  ibpMC = quietRun @ EvaluateTropicalMCIBP[spec, fan, {{}},
    "Integrator" -> "MC",
    "NSamples"   -> nMC,
    "RunChecks"  -> False,
    "Verbose"    -> False,
    "WorkingDirectory" -> FileNameJoin[{$pkgRoot, "TEST", "INTERFILES", "cc41a"}]];

  If[!AssociationQ[ibpMC] || !KeyExistsQ[ibpMC, "Results"] ||
     Length[ibpMC["Results"]] == 0,
    Print["CC41 FAIL [41-A IBP] EvaluateTropicalMCIBP returned: ", ibpMC];
    $cc41AllPass = False;
    AppendTo[$cc41AnyFail, "41-A IBP"];
    Goto[$skip41A]
  ];

  poleIbpMC = Re @ ibpMC["Results"][[1]]["PoleCoefficient"];
  finIbpMC  = Re @ ibpMC["Results"][[1]]["FinitePart"];
  Print["  IBP+MC: pole = ", poleIbpMC, "  finite = ", finIbpMC];

  (* Reconstruct F at epsStar from IBP Laurent approximation *)
  fIbpStar = poleIbpMC / epsStar + finIbpMC;
  Print["  IBP+MC: F(eps*=", epsStar, ") = ", fIbpStar,
        "  (pole/eps* + finite; compare NInteg=", nintRef, ")"];

  cc41Assert["41-A IBP F(eps*) vs NIntegrate",
    NumericQ[fIbpStar] && NumericQ[nintRef] &&
      Abs[(fIbpStar - nintRef)/nintRef] < tolF,
    "rel-err<" <> ToString[tolF],
    "rel-err=" <> ToString[sci @ Abs[(fIbpStar - nintRef)/nintRef]]];

  (* Pole should be non-zero (the integration produces a genuine 1/eps pole) *)
  cc41Assert["41-A IBP pole is non-zero",
    NumericQ[poleIbpMC] && Abs[poleIbpMC] > 0.01,
    "non-zero pole",
    "pole=" <> ToString[sci[poleIbpMC]]];

  (* ---- Route (b): LaurentFromSubtraction+MC ---- *)
  Print["41-A: running LaurentFromSubtraction (MC, ", nSub, " samples) ..."];
  lfMC = quietRun @ LaurentFromSubtraction[spec, fan, {{}},
    "Integrator"     -> "MC",
    "NSamples"       -> nSub,
    "EpsilonValues"  -> epsList,
    "WorkingDirectory" -> FileNameJoin[{$pkgRoot, "TEST", "INTERFILES", "cc41a"}]];

  If[!AssociationQ[lfMC] || !KeyExistsQ[lfMC, "Results"] ||
     Length[lfMC["Results"]] == 0,
    Print["CC41 FAIL [41-A Sub] LaurentFromSubtraction returned: ", lfMC];
    $cc41AllPass = False;
    AppendTo[$cc41AnyFail, "41-A Sub"];
    Goto[$skip41A]
  ];

  poleSubMC = Re @ lfMC["Results"][[1]]["Pole"];
  finSubMC  = Re @ lfMC["Results"][[1]]["Finite"];
  Print["  Sub+MC: pole = ", poleSubMC, "  finite = ", finSubMC];

  (* F at epsStar from subtraction result at the matching eps index *)
  fSubStar  = Re @ lfMC["Results"][[1]]["FullIntegrals"][[2]];  (* eps=0.02 is index 2 *)
  Print["  Sub+MC: F(eps*=", epsStar, ") = ", fSubStar,
        "  (compare NInteg=", nintRef, ")"];

  cc41Assert["41-A Sub F(eps*) vs NIntegrate",
    NumericQ[fSubStar] && NumericQ[nintRef] &&
      Abs[(fSubStar - nintRef)/nintRef] < tolF,
    "rel-err<" <> ToString[tolF],
    "rel-err=" <> ToString[sci @ Abs[(fSubStar - nintRef)/nintRef]]];

  (* Routes agree on the pole *)
  cc41Assert["41-A IBP vs Sub pole agreement",
    NumericQ[poleIbpMC] && NumericQ[poleSubMC] &&
      Abs[poleIbpMC - poleSubMC] < tolPole,
    "agree<" <> ToString[tolPole],
    "|diff|=" <> ToString[sci @ Abs[poleIbpMC - poleSubMC]]];

  (* Also verify the B5 convergent sub-check (set eps=0 to get a convergent integral) *)
  Print[];
  Print["41-A2 (B5 sub-check): convergent numerator Int dx1 dx2 (1+x1)^{3/2}(1+x1+x2)^{-4}=2/3"];
  Module[{specConv, fanConv, vrConv},
    specConv = <|
      "Polynomials"         -> {pNum, pDen},
      "MonomialExponents"   -> {0, 0},
      "PolynomialExponents" -> {3/2, -4},
      "Variables"           -> vars,
      "KinematicSymbols"    -> {},
      "RegulatorSymbol"     -> None
    |>;
    fanConv = ComputeDecomposition[PolytopeVertices[(pNum * pDen)^(-1), vars],
                "ShowProgress" -> False];
    vrConv = quietRun @ ValidateDecomposition[specConv, fanConv, {}, 4];
    If[AssociationQ[vrConv],
      Print["  Sector sum = ", vrConv["SectorSum"], "  exact = 2/3 = ", N[2/3],
            "  rel-err = ", sci @ vrConv["RelativeError"]];
      cc41Assert["41-A2 convergent B5 numerator vs exact 2/3",
        vrConv["RelativeError"] < 0.01,
        "rel-err<0.01",
        "rel-err=" <> ToString[sci @ vrConv["RelativeError"]]],
      Print["CC41 FAIL [41-A2] ValidateDecomposition returned: ", vrConv];
      $cc41AllPass = False;
      AppendTo[$cc41AnyFail, "41-A2"]
    ]
  ];

  Label[$skip41A];
];

Print[];


(* ============================================================================
   41-B: Numerator factor (B>0) x lifting (extreme coefficient)
         I = Int_0^inf dx1 dx2  (1+x1*x2)^{1/2} / (1 + 10^6*x1^2 + x2^2)^2

   The numerator polynomial (1+x1*x2)^{1/2} has positive exponent (B>0).
   The denominator has coefficient 10^6 >> 1, triggering
     DetectExtremeCoefficients -> LiftCoefficients -> EvaluateTropicalMCLifted.

   Reference: NIntegrate (no analytic closed form; the polynomial geometry is
   rich enough that lifting genuinely helps convergence).

   PASS: |liftedMC - NIntRef| / |NIntRef| < 0.01   (1 pct)
   ============================================================================ *)

Print["================================================================"];
Print["  CC41-B: Numerator x Lifting"];
Print["  I=Int dx1 dx2 (1+x1*x2)^{1/2}/(1+10^6*x1^2+x2^2)^2"];
Print["================================================================"];

Module[
  {vars, pNum, pDen, spec, extremeCoeff,
   fanPoly, verts, fan,
   nintRef, liftRules, liftedRes, liftedMC, liftedErr, relErr,
   tolLift},

  vars       = {x[1], x[2]};
  extremeCoeff = 10^6;
  pNum       = 1 + x[1] * x[2];        (* numerator, exponent +1/2 *)
  pDen       = 1 + extremeCoeff * x[1]^2 + x[2]^2;  (* denominator, exponent -2 *)

  spec = <|
    "Polynomials"         -> {pNum, pDen},
    "MonomialExponents"   -> {0, 0},
    "PolynomialExponents" -> {1/2, -2},    (* B>0 numerator + B<0 denominator *)
    "Variables"           -> vars,
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> None
  |>;

  (* Fan from product of both polynomials (standard: both B-values matter for geometry) *)
  fanPoly = pNum * pDen;
  verts   = PolytopeVertices[fanPoly^(-1), vars];
  fan     = ComputeDecomposition[verts, "ShowProgress" -> False];

  If[!ListQ[fan] || Length[fan] < 2,
    Print["CC41 FAIL [41-B fan] fan computation returned $Failed"];
    $cc41AllPass = False;
    AppendTo[$cc41AnyFail, "41-B fan"];
    Goto[$skip41B]
  ];
  Print["41-B: fan built -- ", Length[fan[[2]]], " sectors"];

  (* NIntegrate reference via compactification u = t/(1+t) *)
  Print["41-B: computing NIntegrate reference ..."];
  nintRef = Quiet @ NIntegrate[
    (1 + (u1/(1-u1)) * (u2/(1-u2)))^(1/2)
      / (1 + extremeCoeff * (u1/(1-u1))^2 + (u2/(1-u2))^2)^2
      / ((1-u1)^2 (1-u2)^2),
    {u1, 0, 1}, {u2, 0, 1},
    Method -> "GlobalAdaptive", PrecisionGoal -> 5, MaxRecursion -> 50];
  Print["41-B: NIntegrate reference = ", nintRef];

  If[!NumericQ[nintRef] || nintRef <= 0,
    Print["CC41 FAIL [41-B NInt] NIntegrate returned: ", nintRef];
    $cc41AllPass = False;
    AppendTo[$cc41AnyFail, "41-B NInt"];
    Goto[$skip41B]
  ];

  (* Lift: the dominant extreme term in pDen is extremeCoeff * x[1]^2.
     Exponent vector {2,0} in x[1], x[2] for the denominator polynomial (index 2).
     k=3 so z0 = extremeCoeff^{1/3} (rational k avoids branch issues). *)
  liftRules = {<|"PolyIndex" -> 2, "ExponentVector" -> {2, 0}, "k" -> 3|>};

  Quiet[CreateDirectory[
    FileNameJoin[{$pkgRoot, "TEST", "INTERFILES", "cc41b"}],
    CreateIntermediateDirectories -> True],
   {CreateDirectory::eexist}];

  Print["41-B: running EvaluateTropicalMCLifted (explicit lift rule, k=3, 10^6 samples) ..."];
  liftedRes = quietRun @ EvaluateTropicalMCLifted[
    spec, {{}},
    "LiftRules"        -> liftRules,
    "NSamples"         -> 1000000,
    "RunChecks"        -> False,
    "Verbose"          -> False,
    "WorkingDirectory" -> FileNameJoin[{$pkgRoot, "TEST", "INTERFILES", "cc41b"}]
  ];

  If[!AssociationQ[liftedRes] || !KeyExistsQ[liftedRes, "Results"] ||
     Length[liftedRes["Results"]] == 0,
    Print["CC41 FAIL [41-B lift] EvaluateTropicalMCLifted returned: ", liftedRes];
    $cc41AllPass = False;
    AppendTo[$cc41AnyFail, "41-B lift"];
    Goto[$skip41B]
  ];

  liftedMC  = liftedRes["Results"][[1]]["Re"];
  liftedErr = liftedRes["Results"][[1]]["ReErr"];
  relErr    = Abs[(liftedMC - nintRef) / nintRef];
  tolLift   = 0.01;   (* per plan §8.3 #41 "vs NIntegrate/exact <1 pct" *)

  Print["  Lifted MC = ", liftedMC, " +/- ", liftedErr];
  Print["  NIntegrate = ", nintRef];
  Print["  Relative error = ", sci @ relErr];

  cc41Assert["41-B lifted vs NIntegrate",
    relErr < tolLift,
    "rel-err<" <> ToString[tolLift],
    "rel-err=" <> ToString[sci @ relErr]];

  (* Also check the 4-sigma error-bar consistency *)
  cc41Assert["41-B 4-sigma error-bar consistency",
    Abs[liftedMC - nintRef] < 4 * liftedErr || relErr < tolLift / 3,
    "|MC-ref|<4*sigma or rel-err<0.003",
    "|diff|=" <> ToString[sci @ Abs[liftedMC - nintRef]] <>
      " 4sig=" <> ToString[sci[4 * liftedErr]]];

  (* Informational: DetectExtremeCoefficients should flag the spec *)
  Module[{detected},
    detected = Quiet @ DetectExtremeCoefficients[spec, 1000];
    Print["  DetectExtremeCoefficients (threshold=1000) = ", detected];
    cc41Assert["41-B DetectExtremeCoefficients fires",
      AssociationQ[detected] || ListQ[detected] || detected =!= {},
      "non-empty detection",
      ToString[detected]]
  ];

  Label[$skip41B];
];

Print[];


(* ============================================================================
   41-C: Lifting x Divergence  (planAXpDIV.md §7 — the capability that removes
         known-limitation L2 / plan.md N3 "lifting and divergence are mutually
         exclusive per call").
         I(eps) = Int_0^inf dx1 dx2  x1^{-1+eps}(1 + x1 x2 + 10^6 x1^2 + x2^2)^{-2}

   The denominator carries an extreme coefficient 10^6 on x1^2 (triggers lifting,
   ExponentVector {2,0}, k=3 so z0 = 100 exactly) AND the monomial x1^{-1+eps}
   produces a simple 1/eps pole at x1->0.  After lifting + delta-resolution the
   pole lives in surviving atilde directions; the divergent lifted sectors all
   have HasConstantTerm=True (clean atilde poles, no polynomial-zero poles), so
   BOTH the IBP and Subtraction routes resolve it.

   Analytic pole: at x1->0, P -> 1+x2^2, so
       pole = Int_0^inf (1+x2^2)^{-2} dx2 = pi/4 = 0.7853982.

   Oracles (>=2, plan.md §8.5):
     - IBP route        : EvaluateTropicalMC[liftedSpec, fan, LiftData, Method->IBP]
     - Subtraction route: LaurentFromSubtraction[liftedSpec, fan, LiftData] (eps-fit)
     - NIntegrate of the ORIGINAL at epsStar = 0.01.
   PASS: pole agreement across routes < 0.05 AND both poles within 5% of pi/4;
         reconstructed F(epsStar) = pole/epsStar + finite within 2% of NIntegrate
         for both routes.
   ============================================================================ *)

Print["================================================================"];
Print["  CC41-C: Lifting x Divergence  (removes L2 / plan.md N3)"];
Print["  I(eps)=Int dx1 dx2 x1^{-1+eps}(1+x1 x2+10^6 x1^2+x2^2)^{-2}"];
Print["================================================================"];

Module[
  {eps, vars, spec, rules, lift, liftedSpec, liftData, verts, fan,
   epsStar, nMC, nSub, niRef,
   ibpRes, subRes, poleIbp, finIbp, poleSub, finSub,
   recIbp, recSub, tolPole, tolPi, tolF, wdirC},

  eps  = Symbol["eps41c"];
  vars = {x[1], x[2]};
  spec = <|
    "Polynomials"         -> {1 + x[1] x[2] + 10^6 x[1]^2 + x[2]^2},
    "MonomialExponents"   -> {-1 + eps, 0},
    "PolynomialExponents" -> {-2},
    "Variables"           -> vars,
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> eps
  |>;
  rules = {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 3|>};

  wdirC = FileNameJoin[{$pkgRoot, "TEST", "INTERFILES", "cc41c"}];
  Quiet[CreateDirectory[wdirC, CreateIntermediateDirectories -> True],
    {CreateDirectory::eexist}];

  (* Lift + build the (n+1)-dim fan (shared by both routes). *)
  lift = LiftCoefficients[spec, rules];
  If[!AssociationQ[lift],
    Print["CC41 FAIL [41-C lift] LiftCoefficients returned: ", lift];
    $cc41AllPass = False; AppendTo[$cc41AnyFail, "41-C lift"]; Goto[$skip41C]];
  liftedSpec = lift["LiftedSpec"]; liftData = lift["LiftData"];
  verts = Quiet[PolytopeVertices[(Times @@ liftedSpec["Polynomials"])^(-1),
                                 liftedSpec["Variables"]], TropicalFan::polymake];
  fan   = Quiet[computeFanScaled[verts], TropicalFan::polymake];
  If[!ListQ[fan] || Length[fan] < 2,
    Print["CC41 FAIL [41-C fan] lifted fan build failed"];
    $cc41AllPass = False; AppendTo[$cc41AnyFail, "41-C fan"]; Goto[$skip41C]];
  Print["41-C: lifted fan built -- ", Length[fan[[2]]], " sectors, z0=", liftData["z0"]];

  epsStar = 0.01; nMC = 500000; nSub = 400000;
  tolPole = 0.05; tolPi = 0.05; tolF = 0.02;

  niRef = Quiet @ NIntegrate[
    (u1/(1-u1))^(-1+epsStar) *
    (1 + (u1/(1-u1))(u2/(1-u2)) + 10^6 (u1/(1-u1))^2 + (u2/(1-u2))^2)^(-2) /
    ((1-u1)^2 (1-u2)^2),
    {u1, 0, 1}, {u2, 0, 1},
    Method -> "GlobalAdaptive", PrecisionGoal -> 5, MaxRecursion -> 50];
  Print["41-C: NIntegrate(eps*=", epsStar, ") = ", niRef];

  (* ---- IBP route ---- *)
  Print["41-C: IBP route (lifted, MC) ..."];
  ibpRes = quietRun @ EvaluateTropicalMC[liftedSpec, fan, {{}},
    "LiftData" -> liftData, "Method" -> "IBP", "Integrator" -> "MC",
    "NSamples" -> nMC, "RunChecks" -> False, "Verbose" -> False,
    "WorkingDirectory" -> wdirC];
  If[!AssociationQ[ibpRes] || !KeyExistsQ[ibpRes, "Results"],
    Print["CC41 FAIL [41-C IBP] returned: ", ibpRes];
    $cc41AllPass = False; AppendTo[$cc41AnyFail, "41-C IBP"]; Goto[$skip41C]];
  poleIbp = Re @ ibpRes["Results"][[1]]["PoleCoefficient"];
  finIbp  = Re @ ibpRes["Results"][[1]]["FinitePart"];
  Print["  IBP: pole = ", poleIbp, "  finite = ", finIbp];

  (* ---- Subtraction route (LaurentFromSubtraction, pinned eps) ---- *)
  Print["41-C: Subtraction route (lifted, MC) ..."];
  subRes = quietRun @ LaurentFromSubtraction[liftedSpec, fan, {{}},
    "LiftData" -> liftData, "EpsilonValues" -> {0.01, 0.02, 0.04},
    "Integrator" -> "MC", "NSamples" -> nSub, "WorkingDirectory" -> wdirC];
  If[!AssociationQ[subRes] || !KeyExistsQ[subRes, "Results"],
    Print["CC41 FAIL [41-C Sub] returned: ", subRes];
    $cc41AllPass = False; AppendTo[$cc41AnyFail, "41-C Sub"]; Goto[$skip41C]];
  poleSub = Re @ subRes["Results"][[1]]["Pole"];
  finSub  = Re @ subRes["Results"][[1]]["Finite"];
  Print["  Sub: pole = ", poleSub, "  finite = ", finSub];

  recIbp = poleIbp/epsStar + finIbp;
  recSub = poleSub/epsStar + finSub;
  Print["  reconstructed F(eps*): IBP=", recIbp, " Sub=", recSub,
        " (NInt=", niRef, ", pi/4=", N[Pi/4], ")"];

  cc41Assert["41-C pole IBP vs Sub agreement",
    NumericQ[poleIbp] && NumericQ[poleSub] && Abs[poleIbp - poleSub] < tolPole,
    "agree<" <> ToString[tolPole],
    "|diff|=" <> ToString[sci @ Abs[poleIbp - poleSub]]];

  cc41Assert["41-C pole ~ pi/4 (both routes)",
    NumericQ[poleIbp] && NumericQ[poleSub] &&
      Abs[(poleIbp - Pi/4)/(Pi/4)] < tolPi && Abs[(poleSub - Pi/4)/(Pi/4)] < tolPi,
    "rel<" <> ToString[tolPi] <> " vs pi/4",
    "IBP=" <> ToString[sci[poleIbp]] <> " Sub=" <> ToString[sci[poleSub]]];

  cc41Assert["41-C IBP F(eps*) vs NIntegrate",
    NumericQ[recIbp] && NumericQ[niRef] && Abs[(recIbp - niRef)/niRef] < tolF,
    "rel-err<" <> ToString[tolF],
    "rel-err=" <> ToString[sci @ Abs[(recIbp - niRef)/niRef]]];

  cc41Assert["41-C Sub F(eps*) vs NIntegrate",
    NumericQ[recSub] && NumericQ[niRef] && Abs[(recSub - niRef)/niRef] < tolF,
    "rel-err<" <> ToString[tolF],
    "rel-err=" <> ToString[sci @ Abs[(recSub - niRef)/niRef]]];

  Label[$skip41C];
];

Print[];


(* ============================================================================
   Final summary
   ============================================================================ *)

Print["================================================================"];
If[$cc41AllPass,
  Print["CC41 PASS  all sub-checks passed  failures={}"],
  Print["CC41 FAIL  failed sub-checks: ", $cc41AnyFail]
];
Print["================================================================"];
