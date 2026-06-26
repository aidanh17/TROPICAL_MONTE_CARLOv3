(* L3 MC gate on toy T (lifted + divergent), Case A:
     I(eps) = Int_0^inf dx1 dx2 x1^{-1+eps}(1 + 10^6 x1 x2 + x1^2 + x2^2)^{-2}
   Gates:
     #new-A : IBP route pole/finite  ==  Subtraction route pole/finite
     #new-C : reconstructed I(es) = pole/es + finite  ==  NIntegrate(es)
     pole sanity: pole ~ pi/4 = 0.7853982
*)
$pkgRoot = "/Users/aidanh/Desktop/TROPICAL_MONTE_CARLOv3";
Get[FileNameJoin[{$pkgRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$pkgRoot, "tropical_eval.wl"}]];

wdir = FileNameJoin[{$pkgRoot, "INTERFILES", "sandbox_L3_mc"}];
Quiet[CreateDirectory[wdir, CreateIntermediateDirectories -> True], CreateDirectory::eexist];

eps  = Symbol["epsL3m"];
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
lift = LiftCoefficients[spec, rules];
liftedSpec = lift["LiftedSpec"]; liftData = lift["LiftData"];
verts = Quiet[PolytopeVertices[(Times @@ liftedSpec["Polynomials"])^(-1),
                               liftedSpec["Variables"]], TropicalFan::polymake];
fan   = Quiet[computeFanScaled[verts], TropicalFan::polymake];

(* NIntegrate references of the ORIGINAL at fixed eps *)
NIorig[es_] := Quiet@NIntegrate[
  (u1/(1-u1))^(-1+es) *
  (1 + (u1/(1-u1))(u2/(1-u2)) + 10^6 (u1/(1-u1))^2 + (u2/(1-u2))^2)^(-2) /
  ((1-u1)^2 (1-u2)^2),
  {u1,0,1},{u2,0,1}, Method->"GlobalAdaptive", PrecisionGoal->5, MaxRecursion->50];
ni01 = NIorig[0.01]; ni02 = NIorig[0.02];
Print["L3 MC: NIntegrate(0.01)=", ni01, "  NIntegrate(0.02)=", ni02];

SetAttributes[qr, HoldFirst];
qr[e_] := Quiet[e, {TropicalEval::validate, TropicalEval::divergent,
   TropicalEval::vegasbudget, General::stop, NIntegrate::slwcon, NIntegrate::ncvb}];

(* ---- IBP route ---- *)
Print["L3 MC: IBP route ..."];
ibp = qr@EvaluateTropicalMC[liftedSpec, fan, {{}},
  "LiftData" -> liftData, "Method" -> "IBP", "Integrator" -> "MC",
  "NSamples" -> 400000, "RunChecks" -> False, "Verbose" -> False,
  "WorkingDirectory" -> wdir];
poleIBP = Re@ibp["Results"][[1]]["PoleCoefficient"];
finIBP  = Re@ibp["Results"][[1]]["FinitePart"];
Print["  IBP: pole=", poleIBP, "  finite=", finIBP];

(* ---- Subtraction route (LaurentFromSubtraction, pinned eps) ---- *)
Print["L3 MC: Subtraction route ..."];
sub = qr@LaurentFromSubtraction[liftedSpec, fan, {{}},
  "LiftData" -> liftData, "EpsilonValues" -> {0.01, 0.02, 0.04},
  "Integrator" -> "MC", "NSamples" -> 400000,
  "WorkingDirectory" -> wdir];
poleSub = Re@sub["Results"][[1]]["Pole"];
finSub  = Re@sub["Results"][[1]]["Finite"];
Print["  Sub: pole=", poleSub, "  finite=", finSub];

(* ---- gates ---- *)
recIBP01 = poleIBP/0.01 + finIBP;
recSub01 = poleSub/0.01 + finSub;
Print["L3 MC: pi/4 = ", N[Pi/4]];
Print["L3 MC: reconstructed I(0.01): IBP=", recIBP01, " Sub=", recSub01, " NInt=", ni01];

gPole   = Abs[poleIBP - poleSub] < 0.05;
gPolePi = Abs[poleIBP - Pi/4]/(Pi/4) < 0.05 && Abs[poleSub - Pi/4]/(Pi/4) < 0.05;
gFin    = Abs[finIBP - finSub] < 0.5;
gRecIBP = Abs[(recIBP01 - ni01)/ni01] < 0.03;
gRecSub = Abs[(recSub01 - ni01)/ni01] < 0.03;
Print["L3 MC GATES:"];
Print["  pole IBP~Sub (<0.05):        ", gPole, "  |diff|=", Abs[poleIBP-poleSub]];
Print["  pole ~ pi/4 (<5%):           ", gPolePi];
Print["  finite IBP~Sub (<0.5):       ", gFin, "  |diff|=", Abs[finIBP-finSub]];
Print["  recon I(0.01) IBP vs NInt:   ", gRecIBP, "  rel=", Abs[(recIBP01-ni01)/ni01]];
Print["  recon I(0.01) Sub vs NInt:   ", gRecSub, "  rel=", Abs[(recSub01-ni01)/ni01]];
Print["L3 MC OVERALL: ", If[gPole&&gPolePi&&gFin&&gRecIBP&&gRecSub, "PASS", "FAIL"]];
