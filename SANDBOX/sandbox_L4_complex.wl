(* L4: lift + divergence with a COMPLEX polynomial exponent in Direct mode
   (planAXpDIV.md §4.5).  Toy B with B = -2 + (3/10) I:
     I(eps) = Int x1^{-1+eps}(1 + x1 x2 + 10^6 x1^2 + x2^2)^{-2 + (3/10)I} dx
   Direct mode (exp(B log P)) is always correct; cxSplit is OFF so the
   splitliftdiv guard does not fire.  Check IBP vs Subtraction vs NIntegrate
   (complex). *)
$pkgRoot = "/Users/aidanh/Desktop/TROPICAL_MONTE_CARLOv3";
Get[FileNameJoin[{$pkgRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$pkgRoot, "tropical_eval.wl"}]];
wdir = FileNameJoin[{$pkgRoot, "INTERFILES", "sandbox_L4_complex"}];
Quiet[CreateDirectory[wdir, CreateIntermediateDirectories -> True], CreateDirectory::eexist];
SetAttributes[qr, HoldFirst];
qr[e_] := Quiet[e, {TropicalEval::validate, TropicalEval::vegasbudget, General::stop,
   NIntegrate::slwcon, NIntegrate::ncvb, TropicalFan::polymake}];

eps  = Symbol["epsX"];
vars = {x[1], x[2]};
Bc   = -2 + (3/10) I;
spec = <|"Polynomials" -> {1 + x[1] x[2] + 10^6 x[1]^2 + x[2]^2},
  "MonomialExponents" -> {-1 + eps, 0}, "PolynomialExponents" -> {Bc},
  "Variables" -> vars, "KinematicSymbols" -> {}, "RegulatorSymbol" -> eps|>;
rules = {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 3|>};
lift = LiftCoefficients[spec, rules];
liftedSpec = lift["LiftedSpec"]; liftData = lift["LiftData"];
verts = qr@PolytopeVertices[(Times @@ liftedSpec["Polynomials"])^(-1), liftedSpec["Variables"]];
fan   = qr@computeFanScaled[verts];

es = 0.01;
ni = qr@NIntegrate[
  (u1/(1-u1))^(-1+es) *
  (1 + (u1/(1-u1))(u2/(1-u2)) + 10^6 (u1/(1-u1))^2 + (u2/(1-u2))^2)^Bc /
  ((1-u1)^2 (1-u2)^2),
  {u1,0,1},{u2,0,1}, Method->"GlobalAdaptive", PrecisionGoal->5, MaxRecursion->50];
Print["L4cx: NIntegrate(", es, ") = ", ni];

ibp = qr@EvaluateTropicalMC[liftedSpec, fan, {{}},
  "LiftData"->liftData, "Method"->"IBP", "ComplexExponentMode"->"Direct",
  "Integrator"->"MC", "NSamples"->500000, "RunChecks"->False, "Verbose"->False,
  "WorkingDirectory"->wdir];
poleI = ibp["Results"][[1]]["PoleCoefficient"];
finI  = ibp["Results"][[1]]["FinitePart"];
Print["L4cx IBP: pole=", poleI, "  finite=", finI];

sub = qr@LaurentFromSubtraction[liftedSpec, fan, {{}},
  "LiftData"->liftData, "EpsilonValues"->{0.01,0.02,0.04},
  "ComplexExponentMode"->"Direct", "Integrator"->"MC", "NSamples"->400000,
  "WorkingDirectory"->wdir];
poleS = sub["Results"][[1]]["Pole"];
finS  = sub["Results"][[1]]["Finite"];
Print["L4cx Sub: pole=", poleS, "  finite=", finS];

recI = poleI/es + finI; recS = poleS/es + finS;
Print["L4cx recon(", es, "): IBP=", recI, " Sub=", recS, " NInt=", ni];
gP = Abs[poleI - poleS] < 0.05;
gI = Abs[(recI - ni)/ni] < 0.03;
gS = Abs[(recS - ni)/ni] < 0.03;
Print["L4cx GATES: poleIBP~Sub=", gP, " (|d|=", Abs[poleI-poleS], ")",
      "  recIBP/NInt=", gI, " (", Abs[(recI-ni)/ni], ")",
      "  recSub/NInt=", gS, " (", Abs[(recS-ni)/ni], ")"];
Print["L4cx OVERALL: ", If[gP && gI && gS, "PASS", "FAIL"]];
