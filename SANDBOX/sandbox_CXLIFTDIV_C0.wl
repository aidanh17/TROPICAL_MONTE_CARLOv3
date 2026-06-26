(* SANDBOX/sandbox_CXLIFTDIV_C0.wl  --  planCXLIFTDIV.md Phase C0
   Complex-B  x  lifting  x  divergence, via the (already-working) Subtraction
   route, cross-checked against a complex NIntegrate Laurent fit.

   Toy (planCXLIFTDIV.md §7, cc_46):
       I(eps) = Int_0^inf dx1 dx2  x1^{-1+eps} (1 + x1 x2 + 10^6 x1^2 + x2^2)^{B}
       B = -2 + 0.3 i,   lift {2,0} on poly 1 with new var index k=3.

   Goals:
     (a) Establish the EXACT-ish complex Laurent of the toy from NIntegrate of
         the complex original at several small eps*, fit P_{-1}/eps + P_0.
     (b) Confirm LaurentFromSubtraction + SplitRealImag + LiftData reproduces the
         complex pole + finite (the route the plan claims works today, §thesis).
     (c) Record reference numbers for the C1 (IBP) gate.
*)

$pkgRoot = "/Users/aidanh/Desktop/TROPICAL_MONTE_CARLOv3";
Get[FileNameJoin[{$pkgRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$pkgRoot, "tropical_eval.wl"}]];
Print["C0: packages loaded."];

SetAttributes[qr, HoldFirst];
qr[e_] := Quiet[e, {TropicalEval::validate, TropicalEval::divergent,
   TropicalEval::vegasbudget, General::stop, NIntegrate::slwcon, NIntegrate::ncvb,
   NIntegrate::inumr, NIntegrate::eincr, NIntegrate::izero, TropicalFan::polymake}];

dim    = 2;
poly   = 1 + x[1] x[2] + 10^6 x[1]^2 + x[2]^2;
B      = -2 + 3/10 I;
vars   = {x[1], x[2]};

(* ---- (a) complex NIntegrate of the original at several eps*, fit Laurent ---- *)
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
Print["C0 NIntegrate(eps*) (complex):"];
Do[Print["   eps=", N[epsFit[[i]]], "  I = ", niVals[[i]]], {i, Length[epsFit]}];

(* g(eps) = eps*I(eps) = P_{-1} + P_0 eps + ...  Linear fit -> P_{-1}, P_0. *)
gVals = MapThread[#1 #2 &, {N[epsFit, 30], niVals}];
fitData = MapThread[{#1, #2} &, {N[epsFit, 30], gVals}];
{poleRe, p0Re} = With[{f = Fit[{#1, Re[#2]} & @@@ fitData, {1, ee}, ee]},
  {f /. ee -> 0, Coefficient[f, ee]}];
{poleIm, p0Im} = With[{f = Fit[{#1, Im[#2]} & @@@ fitData, {1, ee}, ee]},
  {f /. ee -> 0, Coefficient[f, ee]}];
poleNI = poleRe + I poleIm;
finNI  = p0Re + I p0Im;
Print["C0 NIntegrate Laurent fit:"];
Print["   pole   ~ ", poleNI];
Print["   finite ~ ", finNI];

(* ---- (b) Subtraction route (lifted, SplitRealImag) ---- *)
eps  = Symbol["eeC0"];
spec = <|"Polynomials" -> {poly},
  "MonomialExponents" -> {-1 + eps, 0},
  "PolynomialExponents" -> {B}, "Variables" -> vars,
  "KinematicSymbols" -> {}, "RegulatorSymbol" -> eps|>;

lift = qr@LiftCoefficients[spec, {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 3|>}];
If[!AssociationQ[lift], Print["C0 FAIL: lift failed"]; Quit[1]];
ls = lift["LiftedSpec"]; ld = lift["LiftData"];
Print["C0 lifted: z0=", ld["z0"]];

verts = qr@PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]];
fan   = qr@computeFanScaled[verts];
Print["C0 lifted fan: ", Length[fan[[2]]], " sectors"];

wdir = FileNameJoin[{$pkgRoot, "INTERFILES", "sandbox_cxliftdiv_c0"}];
Quiet[CreateDirectory[wdir, CreateIntermediateDirectories -> True], CreateDirectory::eexist];

sub = qr@LaurentFromSubtraction[ls, fan, {{}}, "LiftData" -> ld,
  "ComplexExponentMode" -> "SplitRealImag",
  "EpsilonValues" -> {0.01, 0.02, 0.04}, "Integrator" -> "MC",
  "NSamples" -> 800000, "WorkingDirectory" -> wdir, "Verbose" -> False];

If[!AssociationQ[sub] || !KeyExistsQ[sub, "Results"],
  Print["C0 FAIL: Subtraction route returned ", Head[sub]]; Quit[1]];
poleSub = sub["Results"][[1]]["Pole"];
finSub  = sub["Results"][[1]]["Finite"];
Print["C0 Subtraction (SplitRealImag) route:"];
Print["   pole   = ", poleSub];
Print["   finite = ", finSub];

(* ---- (c) compare ---- *)
esStar = 0.01;
reconSub = poleSub/esStar + finSub;
reconNI  = poleNI/esStar + finNI;
relPole  = Abs[poleSub - poleNI]/Abs[poleNI];
relRecon = Abs[reconSub - reconNI]/Abs[reconNI];
Print["C0 comparison @ eps*=", esStar, ":"];
Print["   |pole_Sub - pole_NI|        = ", N@Abs[poleSub - poleNI]];
Print["   rel pole                    = ", N@relPole];
Print["   recon_Sub  = ", reconSub];
Print["   recon_NI   = ", reconNI];
Print["   rel recon                   = ", N@relRecon];

gPole  = Abs[poleSub - poleNI] < 0.05;
gRecon = relRecon < 0.02;
Print["C0 GATES: pole agree(<0.05)=", gPole, "  recon rel(<0.02)=", gRecon];
Print["C0 OVERALL: ", If[gPole && gRecon, "PASS", "FAIL"]];
