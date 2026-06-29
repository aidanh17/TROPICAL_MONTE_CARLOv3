(* ============================================================================
   TEST/cc_3.wl  —  v3 Cross-check #3  (plan.md §8.2, §8.3)

   "Divergent (pole, finite) vs exact Laurent (Γ products)"

   Spec (3D-div):
       Polynomials          = {1 + x1 + x2 + x3}
       MonomialExponents    = {-1 + eps, 0, 0}
       PolynomialExponents  = {-B}

   Analytic (Dirichlet identity):
       I(eps) = Gamma(eps) * Gamma(B-2-eps) / Gamma(B)
       pole  = 1 / ((B-1)(B-2))
       finite = -HarmonicNumber[B-3] / ((B-1)(B-2))

   For B=4:  (pole, finite) = (1/6,  -1/6 )
   For B=5:  (pole, finite) = (1/12, -1/8 )

   Method: IBP (default) and tropical subtraction, both via MC sampler.
   Tier 1: WL + g++ only (no CUBA required).  CUBA-Vegas columns are run
   when CUBA is available and reported, but are never required for PASS/FAIL.

   Prints:
       CC3 PASS <summary>         — all mandatory assertions passed
       CC3 FAIL expected=<e> got=<g>  — first failed assertion

   Ported from:
       OLD_CODE/TROPICAL_MONTE_CARLO/EXAMPLES/test_divergent_crosscheck.wl
   with path/API updates for v3 and the CC3 PASS/FAIL output contract.
   ============================================================================ *)

(* --- locate & load v3 packages (absolute paths, branch-safe) --- *)
$v3Root = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{$v3Root, "tropical_fan.wl"}]];
Get[FileNameJoin[{$v3Root, "tropical_eval.wl"}]];

(* --- helpers --- *)
SetAttributes[quietRun, HoldFirst];
quietRun[expr_] := Block[{$Output = {}}, expr];

sci[x_] := If[N[x] == 0., "0",
  Module[{e = Floor[Log10[Abs[N[x]]]], m},
    m = N[x] / 10^e;
    ToString[NumberForm[m, {3, 2}]] <> "e" <> ToString[e]]];

(* --- overall PASS/FAIL bookkeeping --- *)
$cc3AllPass  = True;
$cc3FirstFail = None;

cc3Assert[tag_String, cond_, expected_, got_] :=
  If[TrueQ[cond],
    Print["  [", tag, "] CC3 PASS"],
    ($cc3AllPass = False;
     If[$cc3FirstFail === None,
       $cc3FirstFail = "expected=" <> ToString[expected] <> " got=" <> ToString[got]];
     Print["  [", tag, "] CC3 FAIL expected=", expected, " got=", got])];

(* ============================================================
   Per-B cross-check: IBP+MC, subtraction+MC (± CUBA Vegas)
   ============================================================ *)
runOneB[Bexp_Integer, cuba_, nMC_, nSub_] := Module[
  {eps, vars, spec, verts, fan,
   poleRef, finRef, epsStar, Fref, epsList, kStar,
   ibpMC, poleIbpMC, finIbpMC,
   ibpVeg, poleIbpVeg, finIbpVeg,
   lfMC, poleSubMC, finSubMC, FsubMC,
   lfVeg, poleSubVeg, finSubVeg, FsubVeg,
   nint, laurFromIBP,
   tolPole = 5.*^-3, tolFin = 1.5*^-2, tolAgree = 1.5*^-2, tolFull = 8.*^-3},

  Print["\n--- CC3 B=", Bexp, " ---"];

  eps  = Symbol["epsCC3"];
  vars = {x[1], x[2], x[3]};
  spec = <|
    "Polynomials"        -> {1 + x[1] + x[2] + x[3]},
    "MonomialExponents"  -> {-1 + eps, 0, 0},
    "PolynomialExponents"-> {-Bexp},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;

  (* fan — use computeFanScaled (K-scaling robust, §6.4) *)
  verts = PolytopeVertices[(1 + x[1] + x[2] + x[3])^(-1), vars];
  fan   = computeFanScaled[verts];
  If[fan === $Failed || !ListQ[fan],
    Print["  CC3 FAIL fan build failed for B=", Bexp];
    $cc3AllPass = False;
    If[$cc3FirstFail === None,
      $cc3FirstFail = "expected=valid fan got=$Failed for B=" <> ToString[Bexp]];
    Return[]];

  (* analytic reference *)
  poleRef = 1 / ((Bexp - 1) (Bexp - 2));
  finRef  = -HarmonicNumber[Bexp - 3] / ((Bexp - 1) (Bexp - 2));
  epsStar = 0.02;
  Fref    = N[Gamma[epsStar] Gamma[(Bexp - 2) - epsStar] / Gamma[Bexp]];
  epsList = {0.01, 0.02, 0.04}; kStar = 2;
  Print["  analytic: pole=", N[poleRef], "  finite=", N[finRef],
        "  F(eps*=0.02)=", Fref];

  (* direct NIntegrate for full-integral consistency *)
  nint = Quiet @ NIntegrate[
    (u1/(1-u1))^(-1+epsStar)
      (1 + u1/(1-u1) + u2/(1-u2) + u3/(1-u3))^(-Bexp)
      / ((1-u1)^2 (1-u2)^2 (1-u3)^2),
    {u1, 0, 1}, {u2, 0, 1}, {u3, 0, 1},
    Method -> "GlobalAdaptive", PrecisionGoal -> 5, MaxRecursion -> 40];
  Print["  NIntegrate(original,eps*)=", nint, "  closed-form=", Fref];

  (* ---- (1) IBP + MC ---- *)
  ibpMC = quietRun @ EvaluateTropicalMCIBP[spec, fan, {{}},
    "Integrator" -> "MC", "NSamples" -> nMC,
    "RunChecks" -> False, "Verbose" -> False];
  poleIbpMC = Re @ ibpMC["Results"][[1]]["PoleCoefficient"];
  finIbpMC  = Re @ ibpMC["Results"][[1]]["FinitePart"];
  Print["  IBP+MC:    pole=", poleIbpMC, "  finite=", finIbpMC];

  (* ---- (2) Subtraction + MC (Laurent fit) ---- *)
  lfMC = quietRun @ LaurentFromSubtraction[spec, fan, {{}},
    "Integrator" -> "MC", "NSamples" -> nSub, "EpsilonValues" -> epsList];
  poleSubMC = Re @ lfMC["Results"][[1]]["Pole"];
  finSubMC  = Re @ lfMC["Results"][[1]]["Finite"];
  FsubMC    = Re @ lfMC["Results"][[1]]["FullIntegrals"][[kStar]];
  Print["  Sub+MC:    pole=", poleSubMC, "  finite=", finSubMC];

  (* ---- optional CUBA-Vegas columns ---- *)
  If[cuba,
    ibpVeg = quietRun @ EvaluateTropicalMCIBP[spec, fan, {{}},
      "Integrator" -> "VEGAS", "NSamples" -> nMC,
      "RunChecks" -> False, "Verbose" -> False];
    poleIbpVeg = Re @ ibpVeg["Results"][[1]]["PoleCoefficient"];
    finIbpVeg  = Re @ ibpVeg["Results"][[1]]["FinitePart"];
    Print["  IBP+Vegas: pole=", poleIbpVeg, "  finite=", finIbpVeg];

    lfVeg = quietRun @ LaurentFromSubtraction[spec, fan, {{}},
      "Integrator" -> "VEGAS", "NSamples" -> nSub, "EpsilonValues" -> epsList];
    poleSubVeg = Re @ lfVeg["Results"][[1]]["Pole"];
    finSubVeg  = Re @ lfVeg["Results"][[1]]["Finite"];
    FsubVeg    = Re @ lfVeg["Results"][[1]]["FullIntegrals"][[kStar]];
    Print["  Sub+Vegas: pole=", poleSubVeg, "  finite=", finSubVeg]];

  (* ===================== PASS/FAIL assertions ===================== *)

  (* A1: pole vs analytic — all MC runs *)
  cc3Assert["B=" <> ToString[Bexp] <> " A1 IBP-MC pole vs analytic",
    Abs[poleIbpMC - poleRef] < tolPole,
    "pole~" <> ToString[N@poleRef], "pole=" <> sci[poleIbpMC] <> " err=" <> sci[Abs[poleIbpMC - poleRef]]];

  cc3Assert["B=" <> ToString[Bexp] <> " A1 Sub-MC pole vs analytic",
    Abs[poleSubMC - poleRef] < tolPole,
    "pole~" <> ToString[N@poleRef], "pole=" <> sci[poleSubMC] <> " err=" <> sci[Abs[poleSubMC - poleRef]]];

  If[cuba,
    cc3Assert["B=" <> ToString[Bexp] <> " A1 IBP-Vegas pole vs analytic",
      Abs[poleIbpVeg - poleRef] < tolPole,
      "pole~" <> ToString[N@poleRef], "pole=" <> sci[poleIbpVeg] <> " err=" <> sci[Abs[poleIbpVeg - poleRef]]];
    cc3Assert["B=" <> ToString[Bexp] <> " A1 Sub-Vegas pole vs analytic",
      Abs[poleSubVeg - poleRef] < tolPole,
      "pole~" <> ToString[N@poleRef], "pole=" <> sci[poleSubVeg] <> " err=" <> sci[Abs[poleSubVeg - poleRef]]]];

  (* A2: finite vs analytic *)
  cc3Assert["B=" <> ToString[Bexp] <> " A2 IBP-MC finite vs analytic",
    Abs[finIbpMC - finRef] < tolFin,
    "finite~" <> ToString[N@finRef], "finite=" <> sci[finIbpMC] <> " err=" <> sci[Abs[finIbpMC - finRef]]];

  cc3Assert["B=" <> ToString[Bexp] <> " A2 Sub-MC finite vs analytic",
    Abs[finSubMC - finRef] < tolFin,
    "finite~" <> ToString[N@finRef], "finite=" <> sci[finSubMC] <> " err=" <> sci[Abs[finSubMC - finRef]]];

  If[cuba,
    cc3Assert["B=" <> ToString[Bexp] <> " A2 IBP-Vegas finite vs analytic",
      Abs[finIbpVeg - finRef] < tolFin,
      "finite~" <> ToString[N@finRef], "finite=" <> sci[finIbpVeg] <> " err=" <> sci[Abs[finIbpVeg - finRef]]];
    cc3Assert["B=" <> ToString[Bexp] <> " A2 Sub-Vegas finite vs analytic",
      Abs[finSubVeg - finRef] < tolFin,
      "finite~" <> ToString[N@finRef], "finite=" <> sci[finSubVeg] <> " err=" <> sci[Abs[finSubVeg - finRef]]]];

  (* A3: IBP vs subtraction agree (the plan.md §8.2 #3/#4 headline) *)
  cc3Assert["B=" <> ToString[Bexp] <> " A3 IBP-MC vs Sub-MC pole agree",
    Abs[poleIbpMC - poleSubMC] < tolPole,
    "|dpole|<" <> sci[tolPole], "|dpole|=" <> sci[Abs[poleIbpMC - poleSubMC]]];

  cc3Assert["B=" <> ToString[Bexp] <> " A3 IBP-MC vs Sub-MC finite agree",
    Abs[finIbpMC - finSubMC] < tolAgree,
    "|dfin|<" <> sci[tolAgree], "|dfin|=" <> sci[Abs[finIbpMC - finSubMC]]];

  If[cuba,
    cc3Assert["B=" <> ToString[Bexp] <> " A3 IBP-Vegas vs Sub-Vegas pole agree",
      Abs[poleIbpVeg - poleSubVeg] < tolPole,
      "|dpole|<" <> sci[tolPole], "|dpole|=" <> sci[Abs[poleIbpVeg - poleSubVeg]]];
    cc3Assert["B=" <> ToString[Bexp] <> " A3 IBP-Vegas vs Sub-Vegas finite agree",
      Abs[finIbpVeg - finSubVeg] < tolAgree,
      "|dfin|<" <> sci[tolAgree], "|dfin|=" <> sci[Abs[finIbpVeg - finSubVeg]]]];

  (* A4: MC vs Vegas agree within each method — informational when CUBA present *)
  If[cuba,
    Print["  [B=", Bexp, " A4 MC vs Vegas agree] ",
      If[Abs[poleIbpMC - poleIbpVeg] < tolPole && Abs[finIbpMC - finIbpVeg] < tolAgree &&
         Abs[poleSubMC - poleSubVeg] < tolPole && Abs[finSubMC - finSubVeg] < tolAgree,
        "PASS (informational)", "NOTE: disagreement (informational, not required for CC3)"]]];

  (* A5: full-integral consistency at eps*=0.02 *)
  laurFromIBP = poleIbpMC / epsStar + finIbpMC;
  cc3Assert["B=" <> ToString[Bexp] <> " A5 Sub-MC full-integral vs closed-form",
    Abs[FsubMC - Fref] / Abs[Fref] < tolFull,
    "F~" <> ToString[Fref], "F=" <> ToString[FsubMC] <> " relerr=" <> sci[Abs[FsubMC - Fref]/Abs[Fref]]];

  cc3Assert["B=" <> ToString[Bexp] <> " A5 Sub-MC vs IBP-Laurent consistency",
    Abs[FsubMC - laurFromIBP] / Abs[Fref] < 1.5*^-2,
    "F~Laurent", "relerr=" <> sci[Abs[FsubMC - laurFromIBP]/Abs[Fref]]];

  If[NumericQ[nint],
    cc3Assert["B=" <> ToString[Bexp] <> " A5 NIntegrate vs closed-form",
      Abs[nint - Fref] / Abs[Fref] < 5.*^-3,
      "NInt~Fref", "relerr=" <> sci[Abs[nint - Fref]/Abs[Fref]]]];
];

(* ============================================================
   Main entry point
   ============================================================ *)
RunCC3[] := Module[{cuba},
  Print["================================================================"];
  Print["  CC3: 3D divergent integral — (pole,finite) vs exact Gamma Laurent"];
  Print["  plan.md §8.2 checks #3 (analytic) and #4 (IBP vs subtraction)"];
  Print["================================================================"];

  cuba = TrueQ[detectCuba[]["Found"]];
  Print["  CUBA-Vegas available: ", cuba,
    If[!cuba, "  (Vegas columns skipped; MC assertions determine PASS/FAIL)", ""]];

  runOneB[4, cuba, 1500000, 2000000];
  runOneB[5, cuba, 1500000, 2000000];

  Print["\n================================================================"];
  If[$cc3AllPass,
    Print["CC3 PASS all assertions passed (B=4 and B=5, IBP+MC and Sub+MC",
      If[cuba, " and Vegas", ""], ")"],
    Print["CC3 FAIL ", $cc3FirstFail]];
  Print["================================================================"];
  $cc3AllPass
];

If[!TrueQ[$CC3NoRun], RunCC3[]];
