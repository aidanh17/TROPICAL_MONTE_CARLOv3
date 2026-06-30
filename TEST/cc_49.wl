(* ============================================================================
   TEST/cc_49.wl  --  v3 Cross-Check #49  (planIBPCX.md §5, family F1)
   GENUINE COMPLEX POLE (theta_k = 0):  the divergent direction is ON-AXIS (real
   1/eps pole), but the imaginary exponent on a NON-divergent slot makes the pole
   RESIDUE and finite part genuinely complex.  This is the "compare the two
   engines" family the user asks for: IBP and Subtraction are independent
   algorithms that must AGREE, and both must match the analytic oracle.

   The strong 4-way agreement (planIBPCX.md §5.1 F1), pole AND finite:
     (1) IBP            (Method->"IBP";  EvaluateTropicalMCIBP)
     (2) Subtraction    (LaurentFromSubtraction; pinned-eps, an independent algo)
     (3) NIntegrate     (eps-fit of the COMPLEX original; shares no engine code)
     (4) exact Gamma-Laurent (closed form)

   Fixture (Dirichlet, a2 = 1 + I c on a NON-divergent slot -> theta_{x1} = 0):
       I(eps) = Int x1^{2eps-1} x2^{I c} (1+x1+x2)^{-B} dx
              = Gamma(2eps) Gamma(1+I c) Gamma(B-2eps-1-I c) / Gamma(B)
       pole   = (1/2) Gamma(1+I c) Gamma(B-1-I c) / Gamma(B)     (genuinely complex)
       finite = Gamma(1+I c)Gamma(B-1-I c)/Gamma(B) (-EulerGamma - psi(B-1-I c))
     B = 3 (real), c = 1/2.  Single 1/eps pole; imaginary weight off the pole var.

   (Contrast cc_50: there the imaginary weight is ON the divergent slot -> theta!=0
    -> OFF-AXIS, no pole; Subtraction is NOT a valid oracle there.  Here theta=0,
    so all four routes are valid and must agree -- the honest distinction of
    planIBPCX.md §5.1.)

   Tier 1: WL + g++ (MC).  Run:  wolframscript -file TEST/cc_49.wl
   ============================================================================ *)

$pkgRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{$pkgRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$pkgRoot, "tropical_eval.wl"}]];
Print["CC49: packages loaded."];

SetAttributes[qr, HoldFirst];
qr[e_] := Quiet[e, {TropicalEval::validate, TropicalEval::divergent, TropicalEval::vegasbudget,
   General::stop, NIntegrate::slwcon, NIntegrate::ncvb, NIntegrate::inumr, NIntegrate::eincr,
   NIntegrate::izero, NIntegrate::precw, TropicalFan::polymake}];

wdir[tag_] := Module[{d = FileNameJoin[{$pkgRoot, "TEST", "INTERFILES", "cc49", tag}]},
  Quiet[CreateDirectory[d, CreateIntermediateDirectories -> True], {CreateDirectory::eexist}]; d];

sci[x_] := Module[{a = N[Abs[x]], e}, If[a == 0., "0",
   e = Floor[Log10[a]]; ToString[NumberForm[a/10^e, {3, 2}]] <> "e" <> ToString[e]]];
fmtC[z_] := ToString[NumberForm[N[Re[z]], {6, 5}]] <> If[Im[z] >= 0, "+", "-"] <>
   ToString[NumberForm[N[Abs[Im[z]]], {6, 5}]] <> "I";

$pass = True; $fail = {};
assert[label_, ok_, detail_] := (
  If[TrueQ[ok], Print["CC49 PASS ", label, "  (", detail, ")"],
     Print["CC49 FAIL ", label, "  (", detail, ")"]; AppendTo[$fail, label]];
  $pass = $pass && TrueQ[ok]);

cplxClose[a_, b_, tol_] := NumericQ[Abs[a - b]] && Abs[Re[a] - Re[b]] < tol && Abs[Im[a] - Im[b]] < tol;

NS = 2000000; epsVals = {0.01, 0.02, 0.04};

Print[""];
Print["----------------------------------------------------------------"];
Print["  F1 genuine COMPLEX pole (theta=0):  IBP / Subtraction / NIntegrate / exact"];
Module[{cc = 1/2, Bv = 3, eps, vars, poly, spec, fan,
        poleRef, finRef, niC, gv, fd, poleNI, finNI, esStar = 0.01, dim = 2,
        rIBP, rSub, pI, fI, pS, fS},
  eps = Symbol["epsCC49"]; vars = {x[1], x[2]}; poly = 1 + x[1] + x[2];
  (* (4) exact Gamma-Laurent *)
  poleRef = N[(1/2) Gamma[1 + I cc] Gamma[Bv - 1 - I cc]/Gamma[Bv], 16];
  finRef  = N[Gamma[1 + I cc] Gamma[Bv - 1 - I cc]/Gamma[Bv]*(-EulerGamma - PolyGamma[0, Bv - 1 - I cc]), 16];
  Print["   exact pole   = ", fmtC[poleRef]];
  Print["   exact finite = ", fmtC[finRef]];

  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {2 eps - 1, I cc},
    "PolynomialExponents" -> {-Bv}, "Variables" -> vars,
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> eps|>;
  fan = qr@ComputeDecomposition[qr@PolytopeVertices[poly^(-1), vars]];
  assert["F1-0 fan built", ListQ[fan] && Length[fan] == 2 && Length[fan[[2]]] > 0,
     "sectors=" <> ToString[If[ListQ[fan], Length[fan[[2]]], fan]]];

  (* (3) NIntegrate of the complex original -> eps-fit Laurent (pole + finite) *)
  niC[es_] := Module[{us, subr, jac, igc},
    us = Table[Unique["u"], {dim}]; subr = Table[x[i] -> us[[i]]/(1 - us[[i]]), {i, dim}];
    jac = Times @@ Table[1/(1 - us[[i]])^2, {i, dim}];
    igc = (x[1]^(2 es - 1) x[2]^(I cc) (poly)^(-Bv) /. subr) jac;
    qr@NIntegrate[igc, Evaluate[Sequence @@ ({#, 0, 1} & /@ us)], Method -> "GlobalAdaptive",
      PrecisionGoal -> 6, WorkingPrecision -> 30, MaxRecursion -> 60]];
  gv = MapThread[#1 #2 &, {N[epsVals, 30], niC /@ N[epsVals, 30]}];  (* eps*I(eps) -> pole=intercept, finite=slope *)
  fd = MapThread[{#1, #2} &, {N[epsVals, 30], gv}];
  poleNI = (Fit[{#1, Re[#2]} & @@@ fd, {1, ee}, ee] /. ee -> 0) +
           I (Fit[{#1, Im[#2]} & @@@ fd, {1, ee}, ee] /. ee -> 0);
  finNI  = Coefficient[Fit[{#1, Re[#2]} & @@@ fd, {1, ee}, ee], ee] +
           I Coefficient[Fit[{#1, Im[#2]} & @@@ fd, {1, ee}, ee], ee];
  Print["   NIntegrate pole = ", fmtC[poleNI], "  finite = ", fmtC[finNI]];

  (* (1) IBP *)
  rIBP = qr@EvaluateTropicalMCIBP[spec, fan, {{}}, "Integrator" -> "MC", "NSamples" -> NS,
    "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> wdir["ibp"]];
  (* (2) Subtraction (pinned-eps, independent algorithm) *)
  rSub = qr@LaurentFromSubtraction[spec, fan, {{}}, "Integrator" -> "MC", "NSamples" -> NS,
    "EpsilonValues" -> epsVals, "WorkingDirectory" -> wdir["sub"], "Verbose" -> False];

  If[!(AssociationQ[rIBP] && KeyExistsQ[rIBP, "Results"]),
    assert["F1 IBP returns a result", False, ToString[Head[rIBP]]],
    (
    If[!(AssociationQ[rSub] && KeyExistsQ[rSub, "Results"]),
      assert["F1 Subtraction returns a result", False, ToString[Head[rSub]]],
      (
      pI = rIBP["Results"][[1]]["PoleCoefficient"]; fI = rIBP["Results"][[1]]["FinitePart"];
      pS = rSub["Results"][[1]]["Pole"];            fS = rSub["Results"][[1]]["Finite"];
      Print["   IBP: pole=", fmtC[pI], " finite=", fmtC[fI]];
      Print["   Sub: pole=", fmtC[pS], " finite=", fmtC[fS]];

      assert["F1-1 IBP pole vs exact", cplxClose[pI, poleRef, 5*^-3],
         "dRe=" <> sci[Re[pI - poleRef]] <> " dIm=" <> sci[Im[pI - poleRef]]];
      assert["F1-2 Subtraction pole vs exact", cplxClose[pS, poleRef, 5*^-3],
         "dRe=" <> sci[Re[pS - poleRef]] <> " dIm=" <> sci[Im[pS - poleRef]]];
      assert["F1-3 NIntegrate pole vs exact", cplxClose[poleNI, poleRef, 1*^-2],
         "dRe=" <> sci[Re[poleNI - poleRef]] <> " dIm=" <> sci[Im[poleNI - poleRef]]];
      assert["F1-4 IBP <-> Subtraction pole agree (two engines)", Abs[pI - pS] < 1*^-2,
         "|dpole|=" <> sci[Abs[pI - pS]]];
      assert["F1-5 IBP finite vs exact", cplxClose[fI, finRef, 2.5*^-2],
         "dRe=" <> sci[Re[fI - finRef]] <> " dIm=" <> sci[Im[fI - finRef]]];
      assert["F1-6 Subtraction finite vs exact", cplxClose[fS, finRef, 2.5*^-2],
         "dRe=" <> sci[Re[fS - finRef]] <> " dIm=" <> sci[Im[fS - finRef]]];
      assert["F1-7 pole genuinely complex (Im != 0)", Abs[Im[poleRef]] > 0.01,
         "Im(pole)=" <> sci[Im[poleRef]]];
      )
    ];
    )
  ];
];

Print[""];
Print["================================================================"];
If[$pass,
  Print["CC49 PASS  genuine complex pole (theta=0): IBP, Subtraction, NIntegrate ",
        "and exact Gamma-Laurent agree on pole & finite (two independent engines + ",
        "two analytic oracles).  failures={}"],
  Print["CC49 FAIL  failed: ", $fail]];
Print["================================================================"];
If[!$pass, Quit[1]];
