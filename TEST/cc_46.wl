(* ============================================================================
   TEST/cc_46.wl  —  Cross-check #46 (planCXLIFTDIV.md §6/§7)
   Complex polynomial exponents  x  lifting  x  divergence (a simple 1/eps pole),
   via ComplexExponentMode -> "SplitRealImag" on BOTH routes:
     - Subtraction (C0): LaurentFromSubtraction on the lifted spec at pinned eps
       (every sector convergent => the lifted SplitRealImag path runs unchanged).
     - IBP (C1): EvaluateTropicalMC[.., Method->"IBP", "ComplexExponentMode"->
       "SplitRealImag"] — per-IBP-piece MonoFactorLog + the brought-down complex
       IBP coefficient reinstate Im(B) as the oscillatory phase.
   Independent oracle: NIntegrate of the COMPLEX original at several small eps,
   fit P_{-1}/eps + P_0 (complex Laurent).

   Each integrand:  I(eps) = Int_[0,inf)^n x1^{-1+eps} prod_j P_j(x)^{B_j} dx,
   with at least one COMPLEX B_j (nonzero Im) and an extreme coefficient (>=10^4)
   that triggers lifting.

   PASS per example (planCXLIFTDIV.md §7):
     (a) complex pole agreement across routes: |pole_IBP - pole_Sub| < 0.05;
     (b) recon F(epsStar)=pole/epsStar+finite within 2% of NIntegrate (complex)
         for BOTH routes;
     (c) the pole is genuinely complex (nonzero Im) and nonzero.

   Tier 1: WL + g++ (MC).  Run:  wolframscript -file TEST/cc_46.wl
   ============================================================================ *)

$pkgRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{$pkgRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$pkgRoot, "tropical_eval.wl"}]];
Print["CC46: packages loaded."];
Print[];

sci[v_] := If[NumericQ[v], ScientificForm[N[v], 3], v];
$cc46Pass = True; $cc46Fail = {};
cc46Assert[label_String, test_, exp_String, got_String] :=
  If[TrueQ[test],
    Print["CC46 PASS  [", label, "]  expected=", exp, "  got=", got],
    Print["CC46 FAIL  [", label, "]  expected=", exp, "  got=", got];
    $cc46Pass = False; AppendTo[$cc46Fail, label]];

SetAttributes[qr, HoldFirst];
qr[e_] := Quiet[e, {TropicalEval::validate, TropicalEval::divergent,
   TropicalEval::vegasbudget, General::stop, NIntegrate::slwcon, NIntegrate::ncvb,
   NIntegrate::inumr, NIntegrate::eincr, NIntegrate::izero, TropicalFan::polymake}];

(* Complex NIntegrate Laurent reference: NIntegrate the COMPLEX original at a few
   small eps via xi=ui/(1-ui), fit g(eps)=eps*I(eps)=P_{-1}+P_0 eps + ... .  *)
niCxLaurent[polys_, polyExps_, dim_, epsList_] := Module[
  {niCx, gv, fd, poleR, poleI, finR, finI},
  niCx[es_] := Module[{us, sub, jac, ig},
    us  = Table[Unique["u"], {dim}];
    sub = Table[x[i] -> us[[i]]/(1 - us[[i]]), {i, dim}];
    jac = Times @@ Table[1/(1 - us[[i]])^2, {i, dim}];
    ig  = (us[[1]]/(1 - us[[1]]))^(-1 + es) *
          (Times @@ MapThread[#1^#2 &, {polys, polyExps}] /. sub) * jac;
    qr@NIntegrate[ig, Evaluate[Sequence @@ ({#, 0, 1} & /@ us)],
      Method -> "GlobalAdaptive", PrecisionGoal -> 6, WorkingPrecision -> 30,
      MaxRecursion -> 60]];
  gv = MapThread[#1 #2 &, {N[epsList, 30], niCx /@ N[epsList, 30]}];
  fd = MapThread[{#1, #2} &, {N[epsList, 30], gv}];
  poleR = Fit[{#1, Re[#2]} & @@@ fd, {1, ee}, ee];
  poleI = Fit[{#1, Im[#2]} & @@@ fd, {1, ee}, ee];
  {(poleR /. ee -> 0) + I (poleI /. ee -> 0),
   Coefficient[poleR, ee] + I Coefficient[poleI, ee]}];

(* Run one complex-B lift x divergence example through all three oracles. *)
runCxExample[label_String, polys_List, polyExps_List, lr_Association,
             esStar_:0.01, nSamp_:800000] := Module[
  {dim, eps, vars, spec, lift, ls, ld, verts, fan, dimOK,
   poleNI, finNI, sub, poleS, finS, recS,
   ibp, poleI, finI, recI, reconNI, wdir},
  dim = Length[lr["ExponentVector"]];
  eps = Symbol["ee46x" <> ToString[Hash[label]]];
  vars = Table[x[i], {i, dim}];
  wdir = FileNameJoin[{$pkgRoot, "TEST", "INTERFILES", "cc46",
           StringReplace[label, {" " -> "_", ":" -> "", "(" -> "", ")" -> "",
                                 "*" -> "", "+" -> "p"}]}];
  Quiet[CreateDirectory[wdir, CreateIntermediateDirectories -> True],
        {CreateDirectory::eexist}];

  Print["----------------------------------------------------------------"];
  Print["  ", label, "   (n=", dim, ")"];

  spec = <|"Polynomials" -> polys,
    "MonomialExponents" -> Join[{-1 + eps}, ConstantArray[0, dim - 1]],
    "PolynomialExponents" -> polyExps, "Variables" -> vars,
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> eps|>;

  lift = qr@LiftCoefficients[spec, {lr}];
  If[!AssociationQ[lift],
    cc46Assert[label <> " lift", False, "LiftedSpec", "lift failed"]; Return[]];
  ls = lift["LiftedSpec"]; ld = lift["LiftData"];
  verts = qr@PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]];
  fan = qr@computeFanScaled[verts];
  dimOK = ListQ[fan] && Length[fan] >= 2 && AllTrue[fan[[2]], Length[#] == dim + 1 &];
  If[!dimOK,
    cc46Assert[label <> " fan", False, "full-dim (n+1) lifted fan", "degenerate"];
    Return[]];
  Print["   lifted fan: ", Length[fan[[2]]], " sectors, z0=", ld["z0"]];

  (* Oracle 1: complex NIntegrate Laurent. *)
  {poleNI, finNI} = niCxLaurent[polys, polyExps, dim, {esStar, 2 esStar, 4 esStar}];
  reconNI = poleNI/esStar + finNI;
  Print["   NIntegrate: pole=", N@poleNI, "  recon=", N@reconNI];

  (* Oracle 2: Subtraction route (C0). *)
  sub = qr@LaurentFromSubtraction[ls, fan, {{}}, "LiftData" -> ld,
    "ComplexExponentMode" -> "SplitRealImag",
    "EpsilonValues" -> {esStar, 2 esStar, 4 esStar}, "Integrator" -> "MC",
    "NSamples" -> nSamp, "WorkingDirectory" -> wdir, "Verbose" -> False];
  If[!AssociationQ[sub] || !KeyExistsQ[sub, "Results"],
    cc46Assert[label <> " Sub", False, "Association[Results]", ToString[Head[sub]]];
    Return[]];
  poleS = sub["Results"][[1]]["Pole"];  finS = sub["Results"][[1]]["Finite"];
  recS  = poleS/esStar + finS;
  Print["   Sub: pole=", poleS, "  recon=", recS];

  (* Oracle 3: IBP route (C1). *)
  ibp = qr@EvaluateTropicalMC[ls, fan, {{}}, "LiftData" -> ld, "Method" -> "IBP",
    "ComplexExponentMode" -> "SplitRealImag", "Integrator" -> "MC",
    "NSamples" -> nSamp, "RunChecks" -> False, "Verbose" -> False,
    "WorkingDirectory" -> wdir];
  If[!AssociationQ[ibp] || !KeyExistsQ[ibp, "Results"],
    cc46Assert[label <> " IBP", False, "Association[Results]", ToString[Head[ibp]]];
    Return[]];
  poleI = ibp["Results"][[1]]["PoleCoefficient"];  finI = ibp["Results"][[1]]["FinitePart"];
  recI  = poleI/esStar + finI;
  Print["   IBP: pole=", poleI, "  recon=", recI];

  (* Gates. *)
  cc46Assert[label <> " IBP pole == Sub pole",
    NumericQ[Abs[poleI - poleS]] && Abs[poleI - poleS] < 0.05,
    "|pole_IBP-pole_Sub|<0.05", "|diff|=" <> ToString[sci@Abs[poleI - poleS]]];
  cc46Assert[label <> " Sub recon vs NIntegrate",
    NumericQ[Abs[recS - reconNI]] && Abs[(recS - reconNI)/reconNI] < 0.02,
    "rel<0.02", "rel=" <> ToString[sci@Abs[(recS - reconNI)/reconNI]]];
  cc46Assert[label <> " IBP recon vs NIntegrate",
    NumericQ[Abs[recI - reconNI]] && Abs[(recI - reconNI)/reconNI] < 0.02,
    "rel<0.02", "rel=" <> ToString[sci@Abs[(recI - reconNI)/reconNI]]];
  cc46Assert[label <> " genuinely complex pole",
    NumericQ[Im[poleS]] && Abs[Im[poleS]] > 0.01 && Abs[poleS] > 0.01,
    "|Im(pole)|>0.01", "Im(pole)=" <> ToString[sci@Im[poleS]]];
];

(* ── Battery: complex-B x lift x divergence ─────────────────────────────── *)

(* 2D canonical toy (planCXLIFTDIV.md §thesis/§7): B = -2 + 0.3 i. *)
runCxExample["2D: (1+x1 x2+10^6 x1^2+x2^2)^(-2+0.3i)",
  {1 + x[1] x[2] + 10^6 x[1]^2 + x[2]^2}, {-2 + 3/10 I},
  <|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 3|>];

(* 3D, two polynomials: real numerator x COMPLEX extreme-coeff denominator. *)
runCxExample["3D: (1+x1+x2+x3)(1+10^5 x1^2+x2^2+x3^2)^(-3+0.2i)",
  {1 + x[1] + x[2] + x[3], 1 + 10^5 x[1]^2 + x[2]^2 + x[3]^2}, {1, -3 + 2/10 I},
  <|"PolyIndex" -> 2, "ExponentVector" -> {2, 0, 0}, "k" -> 1|>];

Print[];
Print["================================================================"];
If[$cc46Pass,
  Print["CC46 PASS  complex-B x lift x divergence verified: IBP route agrees ",
        "with the Subtraction route and NIntegrate (complex Laurent).  failures={}"],
  Print["CC46 FAIL  failed: ", $cc46Fail]];
Print["================================================================"];
If[!$cc46Pass, Quit[1]];
