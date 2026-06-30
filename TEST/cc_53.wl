(* ============================================================================
   TEST/cc_53.wl  —  Cross-check #53  (new_request: IBPCX_LIFTED_DIVERGENT_FEATURE
   UPDATE 2026-06-30)

   The regulator eps is physically REAL, but when it is an undeclared symbol that
   sits inside a polynomial / monomial exponent (the principal-series propagator
   power delta = 3/2 + i nu - eps), Im[.] and Re[.] of that exponent leave stray
   Im[eps] / Re[eps] terms:
       Im[3/2 + I nu - eps] = nu - Im[eps]      (instead of nu)
       Re[3/2 + I nu - eps] = 3/2 - Re[eps]     (instead of 3/2 - eps)
   Consequences before the fix (all observed on the 4-point lifted divergent
   presectors, and reproduced here):
     1. The off-axis divergence guard sees theta_k = nu - Im[eps] as eps-DEPENDENT
        and refuses a CONSTANT off-axis direction (TropicalEval::splitdivmono),
        even though it is the supported off-axis case (s_k = c_k eps + i theta_k,
        regular at eps=0).
     2. Im(B) = nu - Im[eps] leaks into C++ as "(eps).imag()" (uncompilable).
     3. Re(B) = 3/2 - Re[eps] leaks "Derivative(Re)" into the log-insertion C++.

   The fix assumes the regulator is REAL wherever B/theta is formed:
     - ibpImagPole  : Refine[theta, eps in Reals]  (the splitdivmono guard reads it)
     - imagPolyInfo : eps-free Im(B) for the oscillatory phase (ImagPolyExponents)
     - realifyPolyB : Re(B) by subtracting I*Im(B) (keeps eps a TRUE real), not
                      MapAt[Re] (which wraps it in Re[eps]).
   A GENUINELY eps-dependent imaginary part (theta ~ 1/eps, the true out-of-scope
   case) keeps a bare eps and is still refused (splitdivmono) — preserved here as a
   negative control.

   Parts:
     (A) UNIT — ibpImagPole/ibpDivClass collapse the spurious Im[eps] on a lifted
         off-axis slot (nu - Im[eps] -> nu), classify it OffAxis (guard passes),
         keep a true pole at 0, and STILL flag a genuine 1/eps Im as eps-dependent.
     (B) END-TO-END — a lifted divergent integral whose polynomial exponent B
         carries the regulator (B = b + i g - eps): IBP (pole, finite) vs the
         NIntegrate complex Laurent of the original.  Pre-fix this $Failed at
         codegen; post-fix it matches the independent oracle.

   Oracles per numeric PASS (invariant #4): NIntegrate of the complex original
   (Laurent fit), plus the closed structural collapse in (A).

   Tier 1: WL + g++ (MC).  Run:  wolframscript -file TEST/cc_53.wl
   ============================================================================ *)

$pkgRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{$pkgRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$pkgRoot, "tropical_eval.wl"}]];
Print["CC53: packages loaded."];
Print[];

sci[v_] := If[NumericQ[v], ScientificForm[N[v], 3], v];
$cc53Pass = True; $cc53Fail = {};
cc53Assert[label_String, test_, exp_String, got_String] :=
  If[TrueQ[test],
    Print["CC53 PASS  [", label, "]  expected=", exp, "  got=", got],
    Print["CC53 FAIL  [", label, "]  expected=", exp, "  got=", got];
    $cc53Pass = False; AppendTo[$cc53Fail, label]];

SetAttributes[qr, HoldFirst];
qr[e_] := Quiet[e, {TropicalEval::validate, TropicalEval::divergent, TropicalEval::vegasbudget,
   TropicalEval::liftcomplex, TropicalEval::splitdivmono, General::stop, NIntegrate::slwcon,
   NIntegrate::ncvb, NIntegrate::inumr, NIntegrate::eincr, NIntegrate::izero, NIntegrate::precw,
   TropicalFan::polymake, General::munfl, Power::infy, Infinity::indet}];

wdir[tag_] := Module[{d = FileNameJoin[{$pkgRoot, "TEST", "INTERFILES", "cc53", tag}]},
  Quiet[CreateDirectory[d, CreateIntermediateDirectories -> True], {CreateDirectory::eexist}]; d];

(* ============================================================================
   PART A — UNIT: ibpImagPole collapses the spurious Im[eps]
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (A) UNIT: ibpImagPole / ibpDivClass regulator-real collapse"];

Module[{ibpImagPole, ibpDivClass, imagPolyInfo, realifyPolyB, epsR, sdOff, sdBad,
        sdA, th1, th2, thBad, thA, nu, specK, imBK, reBK},
  ibpImagPole  = TropicalEval`Private`ibpImagPole;
  ibpDivClass  = TropicalEval`Private`ibpDivClass;
  imagPolyInfo = TropicalEval`Private`imagPolyInfo;
  realifyPolyB = TropicalEval`Private`realifyPolyB;
  cc53Assert["A private hooks resolve",
    AllTrue[{ibpImagPole, ibpDivClass, imagPolyInfo, realifyPolyB},
      Head[#] === Symbol &], "Symbol", ToString[Head[ibpImagPole]]];

  epsR = Symbol["epsCC53A"];
  (* Lifted SplitRealImag sector: NewExponents realified (Re only).  slot 1 =
     off-axis soft (Re=0), slot 2 = genuine pole (Re=0).  Im(B_1) = nu - Im[eps]
     (nu=1, the delta_+ = 3/2 + I nu - eps shape) routed onto slot 1 via the
     dropped tropical monomial factor D_{1,1}=1. *)
  sdOff = <|"ConeIndex" -> 1, "Dimension" -> 2, "NewExponents" -> {0, 2 epsR},
    "LiftedMonoFactor" -> <|"DExp" -> {{1, 0}}|>,
    "ImagPolyExponents" -> {1 - Im[epsR]}, "LiftedMonoPhase" -> None|>;
  th1 = ibpImagPole[sdOff, 1, epsR]; th2 = ibpImagPole[sdOff, 2, epsR];
  Print["   theta(off-axis slot1) = ", th1, "   theta(pole slot2) = ", th2];

  cc53Assert["A off-axis theta loses Im[eps] (-> constant nu)",
    FreeQ[th1, epsR] && TrueQ[Simplify[th1 == 1]], "1 (FreeQ eps)",
    ToString[th1, InputForm]];
  cc53Assert["A off-axis theta is nonzero (genuinely off-axis, not a pole)",
    ! TrueQ[PossibleZeroQ[th1]], "PossibleZeroQ=False", ToString[PossibleZeroQ[th1]]];
  cc53Assert["A pole slot theta == 0 (genuine 1/eps pole preserved)",
    TrueQ[PossibleZeroQ[th2]], "0", ToString[th2, InputForm]];
  cc53Assert["A ibpDivClass: slot1 OffAxis, slot2 Pole",
    ibpDivClass[sdOff, 1, epsR] === "OffAxis" && ibpDivClass[sdOff, 2, epsR] === "Pole",
    "{OffAxis,Pole}",
    ToString[{ibpDivClass[sdOff, 1, epsR], ibpDivClass[sdOff, 2, epsR]}]];

  (* Im(A) twin of the same contamination (LiftedMonoPhase numerator). *)
  sdA = <|"ConeIndex" -> 1, "Dimension" -> 2, "NewExponents" -> {0, 2 epsR},
    "LiftedMonoFactor" -> None, "ImagPolyExponents" -> None,
    "LiftedMonoPhase" -> <|"Num" -> {2 - 3 Im[epsR], 0}|>|>;
  thA = ibpImagPole[sdA, 1, epsR];
  cc53Assert["A Im(A) twin also collapses (2 - 3 Im[eps] -> 2)",
    FreeQ[thA, epsR] && TrueQ[Simplify[thA == 2]], "2 (FreeQ eps)",
    ToString[thA, InputForm]];

  (* NEGATIVE CONTROL: a genuinely eps-dependent Im(B) ~ 1/eps must NOT collapse
     (theta ~ 1/eps reintroduces the fast oscillation -> still out of scope). *)
  sdBad = sdOff; sdBad["ImagPolyExponents"] = {1/epsR};
  thBad = ibpImagPole[sdBad, 1, epsR];
  cc53Assert["A negative control: genuine eps-dependent theta stays eps-dependent",
    ! FreeQ[thBad, epsR], "!FreeQ eps (still refused)", ToString[thBad, InputForm]];

  (* KINEMATIC SYMBOL in B's imaginary part (planR.md R2): the principal-series
     delta = 3/2 + I nu - eps with nu a SYMBOL.  imagPolyInfo must declare nu real
     so Im(B) = nu (NOT Re[nu]); otherwise realifyPolyB leaves B complex and the
     phase coefficient leaks "(nu).real()" into C++.  This guards the kinematic gap
     a numeric-only fixture would miss. *)
  nu = Symbol["nuCC53A"];
  specK = <|"PolynomialExponents" -> {3/2 + I nu - epsR},
    "KinematicSymbols" -> {nu}, "RegulatorSymbol" -> epsR|>;
  imBK = imagPolyInfo[specK];
  reBK = realifyPolyB[specK, imBK]["PolynomialExponents"];
  cc53Assert["A kinematic nu in B: imagPolyInfo gives {nu}, not {Re[nu]}",
    FreeQ[imBK, Re] && TrueQ[Simplify[imBK == {nu}, nu \[Element] Reals]],
    "{nu}", ToString[imBK, InputForm]];
  cc53Assert["A kinematic nu in B: realifyPolyB output is genuinely real (no Complex)",
    FreeQ[Simplify[reBK, nu \[Element] Reals], Complex],
    "FreeQ Complex", ToString[reBK, InputForm]];
];

(* ============================================================================
   PART B — END-TO-END: lifted divergent with the regulator inside B
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (B) END-TO-END: lifted divergent, B = b + i g - eps  IBP vs NIntegrate"];

Module[{epsD, dim, g, es, polys, Avals, pe, niC, gv, fitRe, fitIm, poleNI, finNI,
        spec, lift, ls, ld, fan, ibp, poleI, finI, relP, relF, ee,
        sub, poleS, finS},
  epsD = Symbol["epsCC53B"]; dim = 2; g = 1/2; es = {0.01, 0.02, 0.04, 0.08};
  polys = {1 + x[1] x[2] + 10^4 x[1]^2 + x[2]^2};
  Avals = {-1 + epsD, 0};
  pe    = {-2 + I g - epsD};   (* regulator AND imaginary part inside B *)

  cc53Assert["B fixture: Im[B] raw carries Im[eps] (exercises the contamination)",
    ! FreeQ[Im[pe[[1]]], epsD], "Im[B] contains Im[eps]", ToString[Im[pe[[1]]], InputForm]];

  niC[esv_] := Module[{us, subr, jac, igc, Aev, Bev},
    Aev = Avals /. epsD -> esv; Bev = pe /. epsD -> esv;
    us = Table[Unique["u"], {dim}]; subr = Table[x[i] -> us[[i]]/(1 - us[[i]]), {i, dim}];
    jac = Times @@ Table[1/(1 - us[[i]])^2, {i, dim}];
    igc = (Times @@ MapThread[#1^#2 &, {Table[x[i], {i, dim}], Aev}] *
           Times @@ MapThread[#1^#2 &, {polys, Bev}] /. subr) jac;
    qr@NIntegrate[igc, Evaluate[Sequence @@ ({#, 0, 1} & /@ us)], Method -> "GlobalAdaptive",
      PrecisionGoal -> 8, WorkingPrecision -> 40, MaxRecursion -> 80]];
  gv = niC /@ N[es, 40];
  fitRe = Fit[MapThread[{#1, Re[#2]} &, {N[es, 40], gv}], {1/ee, 1, ee}, ee];
  fitIm = Fit[MapThread[{#1, Im[#2]} &, {N[es, 40], gv}], {1/ee, 1, ee}, ee];
  poleNI = Coefficient[fitRe, ee, -1] + I Coefficient[fitIm, ee, -1];
  finNI  = Coefficient[fitRe, ee,  0] + I Coefficient[fitIm, ee,  0];
  Print["   NIntegrate Laurent: pole=", N@poleNI, "  finite=", N@finNI];

  spec = <|"Polynomials" -> polys, "MonomialExponents" -> Avals, "PolynomialExponents" -> pe,
    "Variables" -> {x[1], x[2]}, "KinematicSymbols" -> {}, "RegulatorSymbol" -> epsD|>;
  lift = qr@LiftCoefficients[spec, {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 3|>}];
  ls = lift["LiftedSpec"]; ld = lift["LiftData"];
  fan = qr@computeFanScaled[qr@PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]]];

  ibp = qr@EvaluateTropicalMC[ls, fan, {{}}, "LiftData" -> ld, "Method" -> "IBP",
    "ComplexExponentMode" -> "SplitRealImag", "Integrator" -> "MC", "NSamples" -> 4000000,
    "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> wdir["B"]];
  (* Second INDEPENDENT engine oracle (invariant #4): the pinned-eps Subtraction
     route.  This is a genuine pole (theta=0), so Subtraction is a valid oracle
     (cf. cc_46); it exercises the same imagPolyInfo/realifyPolyB fix via a
     different algorithm, and IBP<->Sub agreement does not depend on the
     finite-difference Laurent fit. *)
  sub = qr@LaurentFromSubtraction[ls, fan, {{}}, "LiftData" -> ld,
    "ComplexExponentMode" -> "SplitRealImag", "EpsilonValues" -> es[[;; 3]],
    "Integrator" -> "MC", "NSamples" -> 2000000, "Verbose" -> False,
    "WorkingDirectory" -> wdir["Bsub"]];

  If[!(AssociationQ[ibp] && KeyExistsQ[ibp, "Results"]),
    cc53Assert["B IBP no longer $Failed (compiles + runs)", False,
      "Association[Results]", ToString[Head[ibp]]],
    poleI = ibp["Results"][[1]]["PoleCoefficient"]; finI = ibp["Results"][[1]]["FinitePart"];
    relP = Abs[(poleI - poleNI)/poleNI]; relF = Abs[(finI - finNI)/finNI];
    Print["   IBP: pole=", N@poleI, "  finite=", N@finI];
    cc53Assert["B IBP no longer $Failed (compiles + runs)", True,
      "Association[Results]", "Association"];
    (* Oracle 1 (independent engine): IBP vs Subtraction. *)
    If[AssociationQ[sub] && KeyExistsQ[sub, "Results"],
      poleS = sub["Results"][[1]]["Pole"]; finS = sub["Results"][[1]]["Finite"];
      Print["   Sub: pole=", N@poleS, "  finite=", N@finS];
      cc53Assert["B IBP pole == Subtraction pole (2nd independent engine, #4)",
        NumericQ[Abs[poleI - poleS]] && Abs[poleI - poleS] < 0.05,
        "|d|<0.05", "|d|=" <> ToString[sci@Abs[poleI - poleS]]];
      cc53Assert["B IBP finite == Subtraction finite (2nd independent engine, #4)",
        NumericQ[Abs[finI - finS]] && Abs[(finI - finS)/finS] < 0.05,
        "rel<5e-2", "rel=" <> ToString[sci@Abs[(finI - finS)/finS]]],
      cc53Assert["B Subtraction route returns a result (2nd oracle)", False,
        "Association[Results]", ToString[Head[sub]]]];
    (* Oracle 2 (external): NIntegrate Laurent.  Pole is well-determined; the
       finite part carries a ~1% truncation bias from the finite-eps Laurent fit,
       so it is a looser external sanity check than the IBP<->Sub agreement above. *)
    cc53Assert["B IBP pole vs NIntegrate Laurent",
      NumericQ[relP] && relP < 0.03, "rel<3e-2", "rel=" <> ToString[sci@relP]];
    cc53Assert["B IBP finite vs NIntegrate Laurent (looser: fit bias)",
      NumericQ[relF] && relF < 0.08, "rel<8e-2", "rel=" <> ToString[sci@relF]];
  ];
];

Print[];
Print["================================================================"];
If[$cc53Pass,
  Print["CC53 PASS  regulator-real off-axis / lifted-divergent-with-regulator-in-B ",
        "verified: ibpImagPole collapses the spurious Im[eps] (off-axis no longer ",
        "refused, genuine eps-dependence still refused); the regulator-in-B lifted ",
        "divergent integral compiles and its IBP Laurent matches NIntegrate ",
        "(new_request UPDATE 2026-06-30).  failures={}"],
  Print["CC53 FAIL  failed: ", $cc53Fail]];
Print["================================================================"];
If[!$cc53Pass, Quit[1]];
