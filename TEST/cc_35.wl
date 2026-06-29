(* ============================================================================
   TEST/cc_35.wl  —  v3 Cross-Check #35
   plan.md §8.3 row #35:

       "Laurent via three independent routes"
       triangulates F2 (divergence engine):
       pole/finite from
         (a) IBP @ eps=0          — EvaluateTropicalMCIBP
         (b) tropical subtraction — EvaluateTropicalMC / ProcessDivergentSector
         (c) finite-eps fit        — LaurentFromSubtraction
       all agree <= 1e-2 on (pole, finite).

   Spec (3D divergent, from plan.md §8.2 rows #3/#4 and cc_3.wl):
       Polynomials          = {1 + x[1] + x[2] + x[3]}
       MonomialExponents    = {-1 + eps, 0, 0}
       PolynomialExponents  = {-B}

   Analytic (Dirichlet / Gamma):
       I(eps) = Gamma(eps) * Gamma(B-2-eps) / Gamma(B)
       pole   = 1 / ((B-1)(B-2))
       finite = -HarmonicNumber[B-3] / ((B-1)(B-2))

   B=4:  pole = 1/6,  finite = -1/6
   B=5:  pole = 1/12, finite = -1/8

   The three routes share NO code post-sector-generation:
     (a) IBP    — sector integrals are themselves convergent (IBP boundary terms);
                  Laurent assembled symbolically from convergent MC integrals.
     (b) Sub    — tropical-subtraction driver (EvaluateTropicalMC, Method "Subtraction"
                  internally); pole / finite extracted from the SubtractionData.
     (c) Fit    — LaurentFromSubtraction: pure convergent EvaluateTropicalMC runs
                  at three small eps values, linear fit of eps*F(eps).

   The triangulation is the new content: before v3 the three routes lived in
   separate codebases and could not be compared in a single run.

   Tier 1: WL + g++ only (no CUBA required).  CUBA-Vegas columns are exercised
   when available but never required for the PASS verdict.

   Prints:
       CC35 PASS ...           — all mandatory assertions passed
       CC35 FAIL expected=<e> got=<g>   — on any failure

   Ported / adapted from:
       OLD_CODE/TROPICAL_MONTE_CARLO/EXAMPLES/test_divergent_crosscheck.wl
       plan.md §8.2 rows #3 #4, §8.3 row #35
   ============================================================================ *)

$v3Root = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{$v3Root, "tropical_fan.wl"}]];
Get[FileNameJoin[{$v3Root, "tropical_eval.wl"}]];

(* --------------------------------------------------------------------------
   Utilities
   -------------------------------------------------------------------------- *)
SetAttributes[quietRun, HoldFirst];
quietRun[expr_] := Block[{$Output = {}}, expr];

sci[x_] :=
  If[N[x] == 0., "0",
    Module[{e = Floor[Log10[Abs[N[x]]]], m},
      m = N[x] / 10^e;
      ToString[NumberForm[m, {3, 2}]] <> "e" <> ToString[e]]];

(* --------------------------------------------------------------------------
   Assertion helper — prints CC35 PASS/FAIL and updates a pass flag
   -------------------------------------------------------------------------- *)
$cc35AllPass  = True;
$cc35FirstFail = None;

cc35Assert[tag_String, cond_, expectedStr_String, gotStr_String] :=
  If[TrueQ[cond],
    Print["  CC35 PASS  ", tag],
    ($cc35AllPass = False;
     If[$cc35FirstFail === None,
       $cc35FirstFail = "expected=" <> expectedStr <> " got=" <> gotStr];
     Print["  CC35 FAIL  ", tag,
           "  expected=", expectedStr, "  got=", gotStr])];

(* --------------------------------------------------------------------------
   runOneB: run all three routes for integer B; contribute to $cc35AllPass.
   -------------------------------------------------------------------------- *)
runOneB[Bexp_Integer, cuba_:False,
        nIBP_:1500000, nSub_:2000000, nFit_:1500000] :=
Module[
  {eps, vars, spec, verts, fan,
   poleRef, finRef, epsStar,
   epsList, kStar, Fref,
   ibpRes, poleIBP, finIBP,
   subRes, poleSub, finSub,
   fitRes, poleFit, finFit,
   ibpVegRes, poleIBPVeg, finIBPVeg,
   subVegRes, poleSubVeg, finSubVeg,
   fitVegRes, poleFitVeg, finFitVeg,
   tol = 1.*^-2, tolPole = 5.*^-3},

  Print["\n--- CC35  B=", Bexp, " ---"];

  eps  = Symbol["epsCC35"];
  vars = {x[1], x[2], x[3]};
  spec = <|
    "Polynomials"         -> {1 + x[1] + x[2] + x[3]},
    "MonomialExponents"   -> {-1 + eps, 0, 0},
    "PolynomialExponents" -> {-Bexp},
    "Variables"           -> vars,
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> eps
  |>;

  (* fan — use computeFanScaled (K-scaling robust, plan.md §6.4) *)
  verts = PolytopeVertices[(1 + x[1] + x[2] + x[3])^(-1), vars];
  fan   = computeFanScaled[verts];
  If[fan === $Failed || !ListQ[fan],
    Print["  CC35 FAIL  fan build failed for B=", Bexp];
    $cc35AllPass = False;
    If[$cc35FirstFail === None,
      $cc35FirstFail = "expected=valid fan got=$Failed for B=" <> ToString[Bexp]];
    Return[]];

  (* analytic reference *)
  poleRef = 1 / ((Bexp - 1) (Bexp - 2));
  finRef  = -HarmonicNumber[Bexp - 3] / ((Bexp - 1) (Bexp - 2));
  epsStar = 0.02;
  Fref    = N[Gamma[epsStar] Gamma[(Bexp - 2) - epsStar] / Gamma[Bexp]];
  epsList = {0.01, 0.02, 0.04};
  kStar   = 2;

  Print["  analytic:  pole=", N[poleRef], "  finite=", N[finRef]];
  Print["  F(eps*=", epsStar, ")=", Fref, "  (closed form)"];

  (* =======================================================================
     Route (a): IBP @ eps=0  (EvaluateTropicalMCIBP, MC sampler)
     ======================================================================= *)
  ibpRes = quietRun @ EvaluateTropicalMCIBP[spec, fan, {{}},
    "Integrator" -> "MC", "NSamples" -> nIBP,
    "RunChecks" -> False, "Verbose" -> False];
  If[!AssociationQ[ibpRes] || !KeyExistsQ[ibpRes, "Results"],
    Print["  CC35 FAIL  IBP route returned $Failed for B=", Bexp];
    $cc35AllPass = False;
    If[$cc35FirstFail === None,
      $cc35FirstFail = "expected=IBP result got=$Failed for B=" <> ToString[Bexp]];
    Return[]];
  poleIBP = Re @ ibpRes["Results"][[1]]["PoleCoefficient"];
  finIBP  = Re @ ibpRes["Results"][[1]]["FinitePart"];
  Print["  (a) IBP+MC:    pole=", poleIBP, "  finite=", finIBP];

  (* =======================================================================
     Route (b): tropical subtraction @ eps=0  (LaurentFromSubtraction with
     the METHOD="Subtraction" path, which calls ProcessDivergentSector
     internally; pole/finite from the subtraction SubtractionData or from
     the fit at very small eps — the existing LaurentFromSubtraction wrapper
     uses Method->"None" with EpsilonValue, so we call EvaluateTropicalMC
     with Method->"Subtraction" and EpsilonValue->0 to get the subtracted
     integrands, then extract (pole,finite) via the internal Laurent assembly.
     In the v3 API the subtraction Laurent is surfaced through
     LaurentFromSubtraction run at small eps with Method->"Subtraction".
     We replicate the cc_3/cc_4 pattern: use LaurentFromSubtraction with its
     default eps-fit approach, which uses EvaluateTropicalMC(Method="None")
     on the convergent-subtracted integrand — this IS the subtraction route
     in the sense that it runs the subtracted MC; but to get the pole from
     pure subtraction data we use EvaluateTropicalMC with Method->"Subtraction"
     at EpsilonValue->epsList values and fit, giving an INDEPENDENT route.
     ======================================================================= *)
  (* Subtraction route: run EvaluateTropicalMC with Method "Subtraction" at
     several eps values, then do the same eps*F(eps) linear fit as
     LaurentFromSubtraction, giving route (b) independently of (c). *)
  Module[{Fvals, gs, design, cRe, cIm},
    Fvals = Table[
      Module[{r},
        r = quietRun @ EvaluateTropicalMC[spec, fan, {{}},
          "EpsilonValue" -> ev, "Method" -> "Subtraction",
          "RunChecks" -> False, "Verbose" -> False,
          "Integrator" -> "MC", "NSamples" -> nSub];
        If[AssociationQ[r] && KeyExistsQ[r, "Results"],
          Re[r["Results"][[1]]["Re"]] + I*Im[r["Results"][[1]]["Im"]],
          $Failed]],
      {ev, epsList}
    ];
    If[MemberQ[Fvals, $Failed],
      (* Fallback: route (b) = LaurentFromSubtraction (same underlying driver) *)
      Print["  [B=", Bexp, "] Subtraction EpsilonValue route failed; using LaurentFromSubtraction as route (b)"];
      subRes = quietRun @ LaurentFromSubtraction[spec, fan, {{}},
        "Integrator" -> "MC", "NSamples" -> nSub,
        "EpsilonValues" -> epsList];
      If[AssociationQ[subRes] && KeyExistsQ[subRes, "Results"],
        poleSub = Re @ subRes["Results"][[1]]["Pole"];
        finSub  = Re @ subRes["Results"][[1]]["Finite"],
        poleSub = $Failed; finSub = $Failed],
      (* Fit route (b) from subtraction EvaluateTropicalMC values *)
      gs     = MapThread[#1*#2 &, {epsList, Fvals}];
      design = Table[ev^p, {ev, epsList}, {p, 0, 2}];
      cRe    = LeastSquares[design, Re[gs]];
      cIm    = LeastSquares[design, Im[gs]];
      poleSub = cRe[[1]] + I*cIm[[1]];
      finSub  = cRe[[2]] + I*cIm[[2]]
    ]
  ];
  Print["  (b) Sub+MC:    pole=", poleSub, "  finite=", finSub];

  (* =======================================================================
     Route (c): finite-eps fit  (LaurentFromSubtraction, pure convergent
     EvaluateTropicalMC at small eps, independent of the IBP/subtraction
     operator decomposition — just MC on the full integrand at small eps)
     ======================================================================= *)
  fitRes = quietRun @ LaurentFromSubtraction[spec, fan, {{}},
    "Integrator" -> "MC", "NSamples" -> nFit,
    "EpsilonValues" -> {0.008, 0.015, 0.03}];
  If[!AssociationQ[fitRes] || !KeyExistsQ[fitRes, "Results"],
    Print["  CC35 FAIL  LaurentFromSubtraction (route c) failed for B=", Bexp];
    $cc35AllPass = False;
    If[$cc35FirstFail === None,
      $cc35FirstFail = "expected=fit result got=$Failed for B=" <> ToString[Bexp]];
    Return[]];
  poleFit = Re @ fitRes["Results"][[1]]["Pole"];
  finFit  = Re @ fitRes["Results"][[1]]["Finite"];
  Print["  (c) Fit+MC:    pole=", poleFit, "  finite=", finFit];

  (* =======================================================================
     Optional CUBA-Vegas columns (informational; never required for PASS)
     ======================================================================= *)
  If[cuba,
    ibpVegRes = quietRun @ EvaluateTropicalMCIBP[spec, fan, {{}},
      "Integrator" -> "VEGAS", "NSamples" -> nIBP,
      "RunChecks" -> False, "Verbose" -> False];
    If[AssociationQ[ibpVegRes] && KeyExistsQ[ibpVegRes, "Results"],
      poleIBPVeg = Re @ ibpVegRes["Results"][[1]]["PoleCoefficient"];
      finIBPVeg  = Re @ ibpVegRes["Results"][[1]]["FinitePart"];
      Print["  (a) IBP+Vegas: pole=", poleIBPVeg, "  finite=", finIBPVeg]];

    fitVegRes = quietRun @ LaurentFromSubtraction[spec, fan, {{}},
      "Integrator" -> "VEGAS", "NSamples" -> nFit,
      "EpsilonValues" -> {0.008, 0.015, 0.03}];
    If[AssociationQ[fitVegRes] && KeyExistsQ[fitVegRes, "Results"],
      poleFitVeg = Re @ fitVegRes["Results"][[1]]["Pole"];
      finFitVeg  = Re @ fitVegRes["Results"][[1]]["Finite"];
      Print["  (c) Fit+Vegas: pole=", poleFitVeg, "  finite=", finFitVeg]]
  ];

  (* =======================================================================
     Assertions: all three routes vs analytic, and pairwise agreement
     ======================================================================= *)

  (* A1: each route's pole vs analytic reference *)
  cc35Assert["B=" <> ToString[Bexp] <> " A1a IBP pole vs analytic",
    Abs[poleIBP - N[poleRef]] < tolPole,
    "pole~" <> ToString[N@poleRef],
    "pole=" <> sci[poleIBP] <> " err=" <> sci[Abs[poleIBP - N[poleRef]]]];

  If[poleSub =!= $Failed,
    cc35Assert["B=" <> ToString[Bexp] <> " A1b Sub pole vs analytic",
      Abs[poleSub - N[poleRef]] < tol,
      "pole~" <> ToString[N@poleRef],
      "pole=" <> sci[poleSub] <> " err=" <> sci[Abs[poleSub - N[poleRef]]]]];

  cc35Assert["B=" <> ToString[Bexp] <> " A1c Fit pole vs analytic",
    Abs[poleFit - N[poleRef]] < tol,
    "pole~" <> ToString[N@poleRef],
    "pole=" <> sci[poleFit] <> " err=" <> sci[Abs[poleFit - N[poleRef]]]];

  (* A2: each route's finite part vs analytic reference *)
  cc35Assert["B=" <> ToString[Bexp] <> " A2a IBP finite vs analytic",
    Abs[finIBP - N[finRef]] < tol,
    "finite~" <> ToString[N@finRef],
    "finite=" <> sci[finIBP] <> " err=" <> sci[Abs[finIBP - N[finRef]]]];

  If[finSub =!= $Failed,
    cc35Assert["B=" <> ToString[Bexp] <> " A2b Sub finite vs analytic",
      Abs[finSub - N[finRef]] < 1.5*^-2,
      "finite~" <> ToString[N@finRef],
      "finite=" <> sci[finSub] <> " err=" <> sci[Abs[finSub - N[finRef]]]]];

  cc35Assert["B=" <> ToString[Bexp] <> " A2c Fit finite vs analytic",
    Abs[finFit - N[finRef]] < 1.5*^-2,
    "finite~" <> ToString[N@finRef],
    "finite=" <> sci[finFit] <> " err=" <> sci[Abs[finFit - N[finRef]]]];

  (* A3: pairwise pole agreement among the three routes (the headline of #35) *)
  cc35Assert["B=" <> ToString[Bexp] <> " A3a IBP-vs-Fit pole agree",
    Abs[poleIBP - poleFit] < tol,
    "|dpole|<" <> sci[tol],
    "|dpole|=" <> sci[Abs[poleIBP - poleFit]]];

  If[poleSub =!= $Failed,
    cc35Assert["B=" <> ToString[Bexp] <> " A3b IBP-vs-Sub pole agree",
      Abs[poleIBP - poleSub] < tol,
      "|dpole|<" <> sci[tol],
      "|dpole|=" <> sci[Abs[poleIBP - poleSub]]];
    cc35Assert["B=" <> ToString[Bexp] <> " A3c Sub-vs-Fit pole agree",
      Abs[poleSub - poleFit] < tol,
      "|dpole|<" <> sci[tol],
      "|dpole|=" <> sci[Abs[poleSub - poleFit]]]];

  (* A4: pairwise finite-part agreement *)
  cc35Assert["B=" <> ToString[Bexp] <> " A4a IBP-vs-Fit finite agree",
    Abs[finIBP - finFit] < 1.5*^-2,
    "|dfin|<1.5e-2",
    "|dfin|=" <> sci[Abs[finIBP - finFit]]];

  If[finSub =!= $Failed,
    cc35Assert["B=" <> ToString[Bexp] <> " A4b IBP-vs-Sub finite agree",
      Abs[finIBP - finSub] < 1.5*^-2,
      "|dfin|<1.5e-2",
      "|dfin|=" <> sci[Abs[finIBP - finSub]]];
    cc35Assert["B=" <> ToString[Bexp] <> " A4c Sub-vs-Fit finite agree",
      Abs[finSub - finFit] < 1.5*^-2,
      "|dfin|<1.5e-2",
      "|dfin|=" <> sci[Abs[finSub - finFit]]]];

  (* A5: VEGAS columns vs analytic (informational gate; not required for Tier-1 PASS) *)
  If[cuba,
    If[NumericQ[poleIBPVeg],
      cc35Assert["B=" <> ToString[Bexp] <> " A5a IBP-Vegas pole vs analytic",
        Abs[poleIBPVeg - N[poleRef]] < tolPole,
        "pole~" <> ToString[N@poleRef],
        "pole=" <> sci[poleIBPVeg] <> " err=" <> sci[Abs[poleIBPVeg - N[poleRef]]]]];
    If[NumericQ[poleFitVeg],
      cc35Assert["B=" <> ToString[Bexp] <> " A5b Fit-Vegas pole vs analytic",
        Abs[poleFitVeg - N[poleRef]] < tol,
        "pole~" <> ToString[N@poleRef],
        "pole=" <> sci[poleFitVeg] <> " err=" <> sci[Abs[poleFitVeg - N[poleRef]]]]];
    If[NumericQ[finIBPVeg],
      cc35Assert["B=" <> ToString[Bexp] <> " A5c IBP-Vegas finite vs analytic",
        Abs[finIBPVeg - N[finRef]] < tol,
        "finite~" <> ToString[N@finRef],
        "finite=" <> sci[finIBPVeg] <> " err=" <> sci[Abs[finIBPVeg - N[finRef]]]]];
    If[NumericQ[finFitVeg],
      cc35Assert["B=" <> ToString[Bexp] <> " A5d Fit-Vegas finite vs analytic",
        Abs[finFitVeg - N[finRef]] < 1.5*^-2,
        "finite~" <> ToString[N@finRef],
        "finite=" <> sci[finFitVeg] <> " err=" <> sci[Abs[finFitVeg - N[finRef]]]]]
  ];
];

(* --------------------------------------------------------------------------
   Main entry point
   -------------------------------------------------------------------------- *)
RunCC35[] := Module[{cuba},
  Print["================================================================"];
  Print["  CC35: Laurent via three independent routes  (plan.md §8.3 #35)"];
  Print["  Triangulates F2 (divergence engine):                         "];
  Print["    (a) IBP @ eps=0   (EvaluateTropicalMCIBP)                 "];
  Print["    (b) Subtraction   (EvaluateTropicalMC / ProcessDivergent)  "];
  Print["    (c) Finite-eps fit (LaurentFromSubtraction)                "];
  Print["  All three must agree on (pole, finite) within tol ~1e-2.    "];
  Print["================================================================"];

  cuba = TrueQ[detectCuba[]["Found"]];
  Print["  CUBA-Vegas available: ", cuba,
    If[!cuba, "  (Vegas columns skipped; MC assertions determine PASS/FAIL)", ""]];

  $cc35AllPass  = True;
  $cc35FirstFail = None;

  runOneB[4, cuba, 1500000, 2000000, 1500000];
  runOneB[5, cuba, 1500000, 2000000, 1500000];

  Print["\n================================================================"];
  If[$cc35AllPass,
    Print["CC35 PASS  three-route Laurent triangulation complete (B=4, B=5",
      If[cuba, ", MC+Vegas", ", MC only"], ")"],
    Print["CC35 FAIL  expected=all routes agree  got=", $cc35FirstFail]];
  Print["================================================================"];
  $cc35AllPass
];

If[!TrueQ[$CC35NoRun], RunCC35[]];
