(* SANDBOX/sandbox_CXLIFTDIV_C1.wl  --  planCXLIFTDIV.md Phase C1
   Complex-B x lifting x divergence via the IBP route (the principled second
   oracle).  Compares IBP vs the C0 Subtraction route vs NIntegrate (complex).

   Toy:  I(eps) = Int_0^inf dx1 dx2  x1^{-1+eps} (1 + x1 x2 + 10^6 x1^2 + x2^2)^{B}
         B = -2 + 0.3 i,  lift {2,0}, k=3.

   Gates (planCXLIFTDIV.md §6/§7):
     #cx-A  |pole_IBP - pole_Sub| < 0.05
     #cx-B  recon I at eps-star (IBP) within 2% of NIntegrate (complex)
     #cx-C  is exercised separately (real-B reduction) by the byte-identity gate.
*)

$pkgRoot = "/Users/aidanh/Desktop/TROPICAL_MONTE_CARLOv3";
Get[FileNameJoin[{$pkgRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$pkgRoot, "tropical_eval.wl"}]];
Print["C1: packages loaded."];

SetAttributes[qr, HoldFirst];
qr[e_] := Quiet[e, {TropicalEval::validate, TropicalEval::divergent,
   TropicalEval::vegasbudget, General::stop, NIntegrate::slwcon, NIntegrate::ncvb,
   NIntegrate::inumr, NIntegrate::eincr, NIntegrate::izero, TropicalFan::polymake}];

dim  = 2;
poly = 1 + x[1] x[2] + 10^6 x[1]^2 + x[2]^2;
B    = -2 + 3/10 I;
vars = {x[1], x[2]};

(* ---- complex NIntegrate Laurent reference ---- *)
niOriginalCx[es_] := Module[{us, sub, jac, ig},
  us  = Table[Unique["u"], {dim}];
  sub = Table[x[i] -> us[[i]]/(1 - us[[i]]), {i, dim}];
  jac = Times @@ Table[1/(1 - us[[i]])^2, {i, dim}];
  ig  = (us[[1]]/(1 - us[[1]]))^(-1 + es) * (poly^B /. sub) * jac;
  qr@NIntegrate[ig, Evaluate[Sequence @@ ({#, 0, 1} & /@ us)],
    Method -> "GlobalAdaptive", PrecisionGoal -> 6, WorkingPrecision -> 30,
    MaxRecursion -> 60]];
epsFit = {1/100, 2/100, 4/100};
niVals = niOriginalCx /@ N[epsFit, 30];
gVals  = MapThread[#1 #2 &, {N[epsFit, 30], niVals}];
fitData = MapThread[{#1, #2} &, {N[epsFit, 30], gVals}];
poleNI = (Fit[{#1, Re[#2]} & @@@ fitData, {1, ee}, ee] /. ee -> 0) +
       I (Fit[{#1, Im[#2]} & @@@ fitData, {1, ee}, ee] /. ee -> 0);
finNI  = Coefficient[Fit[{#1, Re[#2]} & @@@ fitData, {1, ee}, ee], ee] +
       I Coefficient[Fit[{#1, Im[#2]} & @@@ fitData, {1, ee}, ee], ee];
Print["C1 NIntegrate: pole ~ ", poleNI, "  finite ~ ", finNI];

(* ---- build lifted spec / fan ---- *)
eps  = Symbol["eeC1"];
spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {-1 + eps, 0},
  "PolynomialExponents" -> {B}, "Variables" -> vars,
  "KinematicSymbols" -> {}, "RegulatorSymbol" -> eps|>;
lift = qr@LiftCoefficients[spec, {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 3|>}];
ls = lift["LiftedSpec"]; ld = lift["LiftData"];
verts = qr@PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]];
fan   = qr@computeFanScaled[verts];
Print["C1 lifted fan: ", Length[fan[[2]]], " sectors, z0=", ld["z0"]];
wdir = FileNameJoin[{$pkgRoot, "INTERFILES", "sandbox_cxliftdiv_c1"}];
Quiet[CreateDirectory[wdir, CreateIntermediateDirectories -> True], CreateDirectory::eexist];

(* ---- Subtraction route (C0) ---- *)
sub = qr@LaurentFromSubtraction[ls, fan, {{}}, "LiftData" -> ld,
  "ComplexExponentMode" -> "SplitRealImag",
  "EpsilonValues" -> {0.01, 0.02, 0.04}, "Integrator" -> "MC",
  "NSamples" -> 800000, "WorkingDirectory" -> wdir, "Verbose" -> False];
poleSub = sub["Results"][[1]]["Pole"];  finSub = sub["Results"][[1]]["Finite"];
Print["C1 Subtraction: pole = ", poleSub, "  finite = ", finSub];

(* ---- IBP route (C1, the new capability) ---- *)
ibp = qr@EvaluateTropicalMC[ls, fan, {{}}, "LiftData" -> ld, "Method" -> "IBP",
  "ComplexExponentMode" -> "SplitRealImag", "Integrator" -> "MC",
  "NSamples" -> 800000, "RunChecks" -> False, "Verbose" -> False,
  "WorkingDirectory" -> wdir];
If[!AssociationQ[ibp] || !KeyExistsQ[ibp, "Results"],
  Print["C1 FAIL: IBP route returned ", Head[ibp], "  ", ibp]; Quit[1]];
poleIBP = ibp["Results"][[1]]["PoleCoefficient"];
finIBP  = ibp["Results"][[1]]["FinitePart"];
Print["C1 IBP: pole = ", poleIBP, "  finite = ", finIBP];

(* ---- gates ---- *)
esStar   = 0.01;
reconIBP = poleIBP/esStar + finIBP;
reconNI  = poleNI/esStar + finNI;
gA = Abs[poleIBP - poleSub] < 0.05;
gB = Abs[(reconIBP - reconNI)/reconNI] < 0.02;
gPoleNI = Abs[poleIBP - poleNI] < 0.05;
Print["----"];
Print["C1 #cx-A |pole_IBP - pole_Sub| = ", N@Abs[poleIBP - poleSub], "  (<0.05: ", gA, ")"];
Print["C1 #cx-B recon_IBP = ", reconIBP];
Print["         recon_NI  = ", reconNI];
Print["         rel = ", N@Abs[(reconIBP - reconNI)/reconNI], "  (<0.02: ", gB, ")"];
Print["C1 (aux) |pole_IBP - pole_NI| = ", N@Abs[poleIBP - poleNI], "  (<0.05: ", gPoleNI, ")"];
Print["C1 OVERALL: ", If[gA && gB, "PASS", "FAIL"]];
If[!(gA && gB), Quit[1]];
