(* ============================================================================
   SANDBOX/sandbox_liftdiv_L0.wl  --  planAXpDIV.md Phase L0
   Lifted + divergent toy:
       I(eps) = Int_0^inf dx1 dx2  x1^{-1+eps} (1 + 10^6 x1 x2 + x2^2)^{-2}

   Goal of L0 (math sandbox, no engine changes):
     (a) Establish the EXACT Laurent of the toy in two independent ways:
           - closed form  I(eps) = (1/2) (10^6)^{-eps} Gamma(eps)
                                     Gamma((1-eps)/2) Gamma((3-eps)/2)
           - NIntegrate of the original at small eps*, fit pole/eps + finite.
     (b) Confirm the CURRENT engine refuses the combination (liftdivergent),
         i.e. document the barrier this plan removes.
   ============================================================================ *)

$pkgRoot = "/Users/aidanh/Desktop/TROPICAL_MONTE_CARLOv3";
Get[FileNameJoin[{$pkgRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$pkgRoot, "tropical_eval.wl"}]];
Print["L0: packages loaded."];

(* ---- (a1) closed form -------------------------------------------------- *)
Iexact[eps_] := (1/2) (10^6)^(-eps) Gamma[eps] Gamma[(1-eps)/2] Gamma[(3-eps)/2];

(* Laurent about eps=0: g(eps) = eps*I(eps) is analytic, g(0)=pole, g'(0)=finite *)
gser = Series[eps Iexact[eps], {eps, 0, 1}];
poleExact   = SeriesCoefficient[gser, 0];   (* coefficient of eps^0 in eps*I = pole *)
finiteExact = SeriesCoefficient[gser, 1];   (* coefficient of eps^1 in eps*I = finite *)
Print["L0 closed form:   pole   = ", poleExact, "  = ", N[poleExact]];
Print["L0 closed form:   finite = ", N[finiteExact]];
Print["L0 (pi/4 check):  ", N[Pi/4], "  match=", PossibleZeroQ[poleExact - Pi/4]];

(* ---- (a2) NIntegrate of the ORIGINAL at small eps* --------------------- *)
(* compactify x = u/(1-u), dx = du/(1-u)^2 *)
Iorig[epsStar_?NumericQ] := NIntegrate[
   (u1/(1-u1))^(-1+epsStar) *
   (1 + 10^6 (u1/(1-u1)) (u2/(1-u2)) + (u2/(1-u2))^2)^(-2) /
   ((1-u1)^2 (1-u2)^2),
   {u1, 0, 1}, {u2, 0, 1},
   Method -> "GlobalAdaptive", PrecisionGoal -> 6, MaxRecursion -> 60];

epsList = {0.005, 0.01, 0.02};
Fvals = Iorig /@ epsList;
Print["L0 NIntegrate(original) F(eps*) for eps*=", epsList, ":"];
Print["   ", Fvals];

(* fit g=eps*F as linear+quadratic *)
gFit = MapThread[Times, {epsList, Fvals}];
design = Table[ev^p, {ev, epsList}, {p, 0, 2}];
coef = LeastSquares[design, gFit];
poleNI   = coef[[1]];
finiteNI = coef[[2]];
Print["L0 NIntegrate fit: pole   = ", poleNI, " (exact ", N[poleExact], ")"];
Print["L0 NIntegrate fit: finite = ", finiteNI, " (exact ", N[finiteExact], ")"];
Print["L0 GATE pole rel-err   = ", Abs[(poleNI - poleExact)/poleExact]];
Print["L0 GATE finite rel-err = ", Abs[(finiteNI - finiteExact)/finiteExact]];

(* Also: direct closed-form vs NIntegrate at eps*=0.01 (no fit) *)
Print["L0 cross: Iexact(0.01)=", N[Iexact[0.01]], "  NIntegrate(0.01)=", Fvals[[2]],
      "  rel-err=", Abs[(N[Iexact[0.01]] - Fvals[[2]])/Fvals[[2]]]];

(* ---- (b) confirm the CURRENT engine refuses lifting+divergence --------- *)
Print[];
Print["L0 (b): does the CURRENT engine refuse lift+divergence?"];
Module[{eps, vars, spec, res},
  eps  = Symbol["epsL0"];
  vars = {x[1], x[2]};
  spec = <|
    "Polynomials"         -> {1 + 10^6 x[1] x[2] + x[2]^2},
    "MonomialExponents"   -> {-1 + eps, 0},
    "PolynomialExponents" -> {-2},
    "Variables"           -> vars,
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> eps
  |>;
  res = Quiet @ EvaluateTropicalMCLifted[spec, {{}},
    "LiftRules" -> {<|"PolyIndex" -> 1, "ExponentVector" -> {1, 1}, "k" -> 1|>},
    "NSamples" -> 10000, "RunChecks" -> False, "Verbose" -> False,
    "WorkingDirectory" -> FileNameJoin[{$pkgRoot, "INTERFILES", "sandbox_liftdiv_L0"}]];
  Print["L0 (b): EvaluateTropicalMCLifted returned: ",
        If[res === $Failed, "$Failed (barrier confirmed: refuses lift+divergence)",
           "an Association (already supported?!)"]];
];

Print[];
Print["L0 DONE."];
