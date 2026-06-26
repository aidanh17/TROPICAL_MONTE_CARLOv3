(* ============================================================================
   TEST/cc_45.wl  —  Cross-check #45 (planAXpDIV.md §7, §9)
   Lifting x divergence on a BATTERY of variable-dimension, non-trivial-polynomial
   integrands.  Companion to cc_41-C (the single canonical toy); cc_45 stresses
   the capability across n = 2, 3, 4 with cross terms, linear terms, numerators
   (B>0) and multiple polynomials, and an extreme coefficient that triggers lifting.

   Each integrand is   I(eps) = Int_[0,inf)^n  x1^{-1+eps} prod_j P_j(x)^{B_j} dx,
   i.e. a simple 1/eps pole at x1->0 (regulated by eps) AND an extreme coefficient
   (>=10^4) on one monomial that triggers DetectExtremeCoefficients-style lifting.

   Routes / oracles (planAXpDIV.md §8.5, >=2 per PASS):
     - SUBTRACTION (universal): LaurentFromSubtraction on the lifted spec at pinned
       eps.  This sums the FULL lifted decomposition (exact), so it is robust to
       the HasConstantTerm=False sectors that non-trivial polynomials produce.
     - IBP (clean only): EvaluateTropicalMC[..,Method->IBP].  The per-sector pole
       extraction is exact only when the pole localizes to an effective-exponent
       endpoint (HasConstantTerm-clean).  Examples flagged "clean" additionally
       assert IBP pole == Subtraction pole.
     - NIntegrate of the ORIGINAL at epsStar = 0.01 (independent oracle).

   PASS per example:
     (a) reconstructed F(epsStar) = pole/epsStar + finite (Subtraction route)
         within 2% of NIntegrate;
     (b) pole > 0 (a genuine 1/eps pole is present);
     (c) for "clean" examples, |pole_IBP - pole_Sub| < 0.05 (both routes agree).

   Tier 1: WL + g++ (MC).  Run:  wolframscript -file TEST/cc_45.wl
   ============================================================================ *)

$pkgRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{$pkgRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$pkgRoot, "tropical_eval.wl"}]];
Print["CC45: packages loaded."];
Print[];

sci[v_] := If[NumericQ[v], ScientificForm[N[v], 3], v];
$cc45Pass = True; $cc45Fail = {};
cc45Assert[label_String, test_, exp_String, got_String] :=
  If[TrueQ[test],
    Print["CC45 PASS  [", label, "]  expected=", exp, "  got=", got],
    Print["CC45 FAIL  [", label, "]  expected=", exp, "  got=", got];
    $cc45Pass = False; AppendTo[$cc45Fail, label]];

SetAttributes[qr, HoldFirst];
qr[e_] := Quiet[e, {TropicalEval::validate, TropicalEval::divergent,
   TropicalEval::vegasbudget, General::stop, NIntegrate::slwcon, NIntegrate::ncvb,
   NIntegrate::inumr, NIntegrate::eincr, TropicalFan::polymake}];

(* NIntegrate of the original (monomial x1^{-1+eps}) at fixed eps, via xi=ui/(1-ui) *)
niOriginal[polys_, polyExps_, dim_, es_] := Module[{us, sub, jac, ig},
  us  = Table[Unique["u"], {dim}];
  sub = Table[x[i] -> us[[i]]/(1 - us[[i]]), {i, dim}];
  jac = Times @@ Table[1/(1 - us[[i]])^2, {i, dim}];
  ig  = (us[[1]]/(1 - us[[1]]))^(-1 + es) *
        (Times @@ MapThread[#1^#2 &, {polys, polyExps}] /. sub) * jac;
  qr@NIntegrate[ig, Evaluate[Sequence @@ ({#, 0, 1} & /@ us)],
    Method -> "GlobalAdaptive", PrecisionGoal -> 5, MaxRecursion -> 50]];

(* Run one battery example. *)
runExample[label_String, polys_List, polyExps_List, lr_Association,
           clean_, esStar_:0.01, nSamp_:400000] := Module[
  {dim, eps, vars, spec, lift, ls, ld, verts, fan, ni,
   sub, poleS, finS, recS, ibp, poleI, relS, wdir, dimOK},
  dim = Length[lr["ExponentVector"]];
  eps = Symbol["ee45x" <> ToString[Hash[label]]];
  vars = Table[x[i], {i, dim}];
  wdir = FileNameJoin[{$pkgRoot, "TEST", "INTERFILES", "cc45",
                       StringReplace[label, {" " -> "_", ":" -> "", "(" -> "", ")" -> ""}]}];
  Quiet[CreateDirectory[wdir, CreateIntermediateDirectories -> True], {CreateDirectory::eexist}];

  Print["----------------------------------------------------------------"];
  Print["  ", label, "   (n=", dim, ", clean=", clean, ")"];

  spec = <|"Polynomials" -> polys,
    "MonomialExponents" -> Join[{-1 + eps}, ConstantArray[0, dim - 1]],
    "PolynomialExponents" -> polyExps, "Variables" -> vars,
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> eps|>;

  lift = qr@LiftCoefficients[spec, {lr}];
  If[!AssociationQ[lift],
    cc45Assert[label <> " lift", False, "LiftedSpec", "lift failed"]; Return[]];
  ls = lift["LiftedSpec"]; ld = lift["LiftData"];
  verts = qr@PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]];
  fan = qr@computeFanScaled[verts];
  dimOK = ListQ[fan] && Length[fan] >= 2 && AllTrue[fan[[2]], Length[#] == dim + 1 &];
  If[!dimOK,
    cc45Assert[label <> " fan", False, "full-dim (n+1) lifted fan", "degenerate/failed"];
    Return[]];
  Print["   lifted fan: ", Length[fan[[2]]], " sectors (full-dim), z0=", ld["z0"]];

  ni = niOriginal[polys, polyExps, dim, esStar];
  Print["   NIntegrate(", esStar, ") = ", ni];

  (* Subtraction route (universal). *)
  sub = qr@LaurentFromSubtraction[ls, fan, {{}}, "LiftData" -> ld,
    "EpsilonValues" -> {0.01, 0.02, 0.04}, "Integrator" -> "MC",
    "NSamples" -> nSamp, "WorkingDirectory" -> wdir];
  If[!AssociationQ[sub] || !KeyExistsQ[sub, "Results"],
    cc45Assert[label <> " Sub", False, "Association[Results]", ToString[Head[sub]]];
    Return[]];
  poleS = Re@sub["Results"][[1]]["Pole"]; finS = Re@sub["Results"][[1]]["Finite"];
  recS  = poleS/esStar + finS;
  relS  = If[NumericQ[ni] && ni != 0, Abs[(recS - ni)/ni], Infinity];
  Print["   Sub: pole=", poleS, " finite=", finS, "  recon=", recS,
        "  rel-err=", sci[relS]];

  cc45Assert[label <> " Sub recon vs NIntegrate",
    NumericQ[relS] && relS < 0.02, "rel-err<0.02", "rel-err=" <> ToString[sci[relS]]];
  cc45Assert[label <> " genuine 1/eps pole",
    NumericQ[poleS] && poleS > 0.01, "pole>0.01", "pole=" <> ToString[sci[poleS]]];

  (* IBP route (assert agreement only for clean examples). *)
  If[clean,
    ibp = qr@EvaluateTropicalMC[ls, fan, {{}}, "LiftData" -> ld, "Method" -> "IBP",
      "Integrator" -> "MC", "NSamples" -> nSamp, "RunChecks" -> False,
      "Verbose" -> False, "WorkingDirectory" -> wdir];
    If[AssociationQ[ibp] && KeyExistsQ[ibp, "Results"],
      poleI = Re@ibp["Results"][[1]]["PoleCoefficient"];
      Print["   IBP: pole=", poleI, "  |pole_IBP - pole_Sub|=", Abs[poleI - poleS]];
      cc45Assert[label <> " IBP pole == Sub pole (clean)",
        NumericQ[poleI] && Abs[poleI - poleS] < 0.05,
        "|diff|<0.05", "|diff|=" <> ToString[sci @ Abs[poleI - poleS]]],
      cc45Assert[label <> " IBP route", False, "Association[Results]", ToString[Head[ibp]]]
    ]
  ];
];

(* ── Battery: variable dimension, non-trivial polynomials ──────────────── *)

(* 2D, denominator with a cross term + an extreme x1^2 (clean: pole localizes). *)
runExample["2D-B: (1+x1 x2+10^6 x1^2+x2^2)^-2",
  {1 + x[1] x[2] + 10^6 x[1]^2 + x[2]^2}, {-2},
  <|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 3|>, True];

(* 2D, 6-monomial denominator with linear + cross + extreme x1 x2 (messy). *)
runExample["2D-a: (1+x1+x2+10^4 x1 x2+x1^2+x2^2)^-2",
  {1 + x[1] + x[2] + 10^4 x[1] x[2] + x[1]^2 + x[2]^2}, {-2},
  <|"PolyIndex" -> 1, "ExponentVector" -> {1, 1}, "k" -> 1|>, False];

(* 3D, two cross terms + extreme x1^2 (messy). *)
runExample["3D-a: (1+x1 x2+x2 x3+10^5 x1^2+x2^2+x3^2)^-2",
  {1 + x[1] x[2] + x[2] x[3] + 10^5 x[1]^2 + x[2]^2 + x[3]^2}, {-2},
  <|"PolyIndex" -> 1, "ExponentVector" -> {2, 0, 0}, "k" -> 1|>, False];

(* 3D, TWO polynomials: numerator (B=+1) x extreme-coeff denominator (B=-3) (clean). *)
runExample["3D-c: (1+x1+x2+x3)(1+10^5 x1^2+x2^2+x3^2)^-3",
  {1 + x[1] + x[2] + x[3], 1 + 10^5 x[1]^2 + x[2]^2 + x[3]^2}, {1, -3},
  <|"PolyIndex" -> 2, "ExponentVector" -> {2, 0, 0}, "k" -> 1|>, True];

(* 4D, two cross terms + extreme x1^2 (messy). *)
runExample["4D-a: (1+x1 x2+x3 x4+10^4 x1^2+x2^2+x3^2+x4^2)^-3",
  {1 + x[1] x[2] + x[3] x[4] + 10^4 x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2}, {-3},
  <|"PolyIndex" -> 1, "ExponentVector" -> {2, 0, 0, 0}, "k" -> 1|>, False, 0.01, 500000];

Print[];
Print["================================================================"];
If[$cc45Pass,
  Print["CC45 PASS  lift x divergence verified across n=2,3,4 (non-trivial ",
        "polynomials); Subtraction route robust, IBP agrees where clean.  failures={}"],
  Print["CC45 FAIL  failed: ", $cc45Fail]];
Print["================================================================"];
If[!$cc45Pass, Quit[1]];
