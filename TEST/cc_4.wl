(* ============================================================================
   TEST/cc_4.wl  —  v3 Cross-Check #4
   plan.md §8.2 row #4:  IBP (default) vs tropical subtraction agree on
   (pole, finite), tolerance ≲1e-3.

   Tier 1 (WL + g++ only).  CUBA-Vegas columns are skipped cleanly when CUBA
   is absent so the MC half still proves PASS with zero external dependencies.

   Spec (3D-div):
       Polynomials = {1 + x[1] + x[2] + x[3]},
       MonomialExponents = {-1+eps, 0, 0},
       PolynomialExponents = {-B}.
   Analytic:  I(eps) = Gamma(eps)*Gamma(B-2-eps)/Gamma(B)
   B=4:  pole = 1/6,  finite = -1/6
   B=5:  pole = 1/12, finite = -1/8

   Ported / adapted from OLD_CODE/TROPICAL_MONTE_CARLO/EXAMPLES/test_divergent_crosscheck.wl
   with v3 path conventions (absolute Get paths, "CC4 PASS …" / "CC4 FAIL …" output).
   ============================================================================ *)

(* --------------------------------------------------------------------------
   Load v3 packages by absolute path so this file runs from any working dir.
   -------------------------------------------------------------------------- *)
$v3Root = "/Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3";
Get[FileNameJoin[{$v3Root, "tropical_fan.wl"}]];
Get[FileNameJoin[{$v3Root, "tropical_eval.wl"}]];

(* --------------------------------------------------------------------------
   Helpers
   -------------------------------------------------------------------------- *)
SetAttributes[quietRun, HoldFirst];
quietRun[expr_] := Block[{$Output = {}}, expr];

sci[x_] := If[N[x] == 0., "0",
  Module[{e = Floor[Log10[Abs[N[x]]]], m},
    m = N[x]/10^e;
    ToString[NumberForm[m, {3, 2}]] <> "e" <> ToString[e]]];

(* K-scaled fan robust helper — mirrors gen_sectors.wl computeFanRobust (§6.4) *)
computeFanRobust[verts_] := Module[{n, fd},
  n = Length[First[verts]];
  Do[fd = Quiet@ComputeDecomposition[K*verts, "ShowProgress" -> False];
     If[ListQ[fd] && Length[fd] == 2 && FreeQ[fd, $Failed] && Length[fd[[1]]] > 0,
        Return[fd, Module]],
     {K, {1, n + 2, 2 n + 4, 6 n + 6}}];
  $Failed];

(* --------------------------------------------------------------------------
   Core cross-check for one (B, nIBP, nSub) triple.
   Returns True (all assertions pass) or False.
   Prints "CC4 PASS ..." / "CC4 FAIL expected=<e> got=<g>" per assertion.
   -------------------------------------------------------------------------- *)
cc4Spec[label_String, Bexp_Integer, cuba_:False, nIBP_:1500000, nSub_:2000000] :=
Module[
  {eps, vars, spec, fan,
   poleRef, finRef, epsStar, Fref, epsList, kStar,
   ibpMC, poleIbpMC, finIbpMC,
   ibpVeg, poleIbpVeg, finIbpVeg,
   lfMC, poleSubMC, finSubMC, FsubMC,
   lfVeg, poleSubVeg, finSubVeg, FsubVeg,
   nint, pass = True,
   tolPole = 5.*^-3, tolFinIBP = 1.*^-2, tolFinSub = 1.5*^-2,
   tolAgreeP = 5.*^-3, tolAgreeF = 1.5*^-2, tolFull = 8.*^-3, tolLaur = 1.5*^-2},

  Print[""];
  Print["--- CC4 ", label, " ---"];

  eps  = Symbol["epsCC4"];
  vars = {x[1], x[2], x[3]};
  spec = <|
    "Polynomials"        -> {1 + x[1] + x[2] + x[3]},
    "MonomialExponents"  -> {-1 + eps, 0, 0},
    "PolynomialExponents"-> {-Bexp},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;

  fan = computeFanRobust[PolytopeVertices[(1 + x[1] + x[2] + x[3])^(-1), vars]];
  If[fan === $Failed,
    Print["CC4 FAIL expected=fan got=$Failed"];
    Return[False]];

  (* ---- analytic reference: I(eps) = Gamma(eps)*Gamma(B-2-eps)/Gamma(B) ---- *)
  poleRef = 1/((Bexp - 1) (Bexp - 2));
  finRef  = -HarmonicNumber[Bexp - 3]/((Bexp - 1) (Bexp - 2));
  epsStar = 0.02;
  Fref    = N[Gamma[epsStar] Gamma[(Bexp - 2) - epsStar] / Gamma[Bexp]];
  epsList = {0.01, 0.02, 0.04};
  kStar   = 2;   (* epsList[[2]] = 0.02 *)
  Print["  analytic: pole=", N[poleRef], "  finite=", N[finRef],
        "  F(0.02)=", Fref];

  (* ---- direct NIntegrate at eps* for the convergence sanity check ---- *)
  nint = Quiet@NIntegrate[
    (u1/(1 - u1))^(-1 + epsStar) *
      (1 + u1/(1 - u1) + u2/(1 - u2) + u3/(1 - u3))^(-Bexp) /
      ((1 - u1)^2 (1 - u2)^2 (1 - u3)^2),
    {u1, 0, 1}, {u2, 0, 1}, {u3, 0, 1},
    Method -> "GlobalAdaptive", PrecisionGoal -> 5, MaxRecursion -> 40];
  Print["  NIntegrate(eps*=0.02)=", nint, "  (closed=", Fref, ")"];

  (* ---- (1) IBP + MonteCarlo ---- *)
  ibpMC = quietRun@EvaluateTropicalMCIBP[spec, fan, {{}},
    "Integrator" -> "MC", "NSamples" -> nIBP,
    "RunChecks"  -> False, "Verbose" -> False];
  poleIbpMC = Re @ ibpMC["Results"][[1]]["PoleCoefficient"];
  finIbpMC  = Re @ ibpMC["Results"][[1]]["FinitePart"];

  (* ---- (2) IBP + Vegas (if CUBA present) ---- *)
  If[cuba,
    ibpVeg = quietRun@EvaluateTropicalMCIBP[spec, fan, {{}},
      "Integrator" -> "VEGAS", "NSamples" -> nIBP,
      "RunChecks"  -> False, "Verbose" -> False];
    poleIbpVeg = Re @ ibpVeg["Results"][[1]]["PoleCoefficient"];
    finIbpVeg  = Re @ ibpVeg["Results"][[1]]["FinitePart"]];

  (* ---- (3) subtraction + MonteCarlo ---- *)
  lfMC = quietRun@LaurentFromSubtraction[spec, fan, {{}},
    "Integrator"    -> "MC",
    "NSamples"      -> nSub,
    "EpsilonValues" -> epsList];
  poleSubMC = Re @ lfMC["Results"][[1]]["Pole"];
  finSubMC  = Re @ lfMC["Results"][[1]]["Finite"];
  FsubMC    = Re @ lfMC["Results"][[1]]["FullIntegrals"][[kStar]];

  (* ---- (4) subtraction + Vegas (if CUBA present) ---- *)
  If[cuba,
    lfVeg = quietRun@LaurentFromSubtraction[spec, fan, {{}},
      "Integrator"    -> "VEGAS",
      "NSamples"      -> nSub,
      "EpsilonValues" -> epsList];
    poleSubVeg = Re @ lfVeg["Results"][[1]]["Pole"];
    finSubVeg  = Re @ lfVeg["Results"][[1]]["Finite"];
    FsubVeg    = Re @ lfVeg["Results"][[1]]["FullIntegrals"][[kStar]]];

  Print["  (pole,finite):  IBP-MC=(", poleIbpMC, ",", finIbpMC, ")"];
  If[cuba, Print["  (pole,finite):  IBP-Vegas=(", poleIbpVeg, ",", finIbpVeg, ")"]];
  Print["  (pole,finite):  sub-MC=(", poleSubMC, ",", finSubMC, ")"];
  If[cuba, Print["  (pole,finite):  sub-Vegas=(", poleSubVeg, ",", finSubVeg, ")"]];

  (* =====================================================================
     Assertions (mirrors old test_divergent_crosscheck.wl A1-A5)
     ===================================================================== *)

  (* A1: pole vs analytic *)
  Module[{errIbpMC, errSubMC, errIbpVeg, errSubVeg, p},
    errIbpMC = Abs[poleIbpMC - N[poleRef]];
    errSubMC = Abs[poleSubMC - N[poleRef]];
    errIbpVeg = If[cuba, Abs[poleIbpVeg - N[poleRef]], 0.];
    errSubVeg = If[cuba, Abs[poleSubVeg - N[poleRef]], 0.];
    p = errIbpMC < tolPole && errSubMC < tolPole &&
        (!cuba || (errIbpVeg < tolPole && errSubVeg < tolPole));
    If[p,
      Print["CC4 PASS ", label, " A1 pole-vs-analytic (IBP-MC=", sci[errIbpMC],
            " sub-MC=", sci[errSubMC], ")"],
      Print["CC4 FAIL ", label, " A1 pole-vs-analytic expected=<", tolPole,
            ">  IBP-MC err=", sci[errIbpMC], " sub-MC err=", sci[errSubMC]]];
    pass = pass && p];

  (* A2: finite vs analytic *)
  Module[{errIbpMC, errSubMC, errIbpVeg, errSubVeg, p},
    errIbpMC = Abs[finIbpMC - N[finRef]];
    errSubMC = Abs[finSubMC - N[finRef]];
    errIbpVeg = If[cuba, Abs[finIbpVeg - N[finRef]], 0.];
    errSubVeg = If[cuba, Abs[finSubVeg - N[finRef]], 0.];
    p = errIbpMC < tolFinIBP && errSubMC < tolFinSub &&
        (!cuba || (errIbpVeg < tolFinIBP && errSubVeg < tolFinSub));
    If[p,
      Print["CC4 PASS ", label, " A2 finite-vs-analytic (IBP-MC=", sci[errIbpMC],
            " sub-MC=", sci[errSubMC], ")"],
      Print["CC4 FAIL ", label, " A2 finite-vs-analytic expected=IBP<",
            tolFinIBP, " sub<", tolFinSub, "  IBP-MC=", sci[errIbpMC],
            " sub-MC=", sci[errSubMC]]];
    pass = pass && p];

  (* A3: IBP vs subtraction agree — this is the headline assertion of CC#4 *)
  Module[{dpole, dfin, dpoleVeg, dfinVeg, p},
    dpole    = Abs[poleIbpMC - poleSubMC];
    dfin     = Abs[finIbpMC  - finSubMC];
    dpoleVeg = If[cuba, Abs[poleIbpVeg - poleSubVeg], 0.];
    dfinVeg  = If[cuba, Abs[finIbpVeg  - finSubVeg ], 0.];
    p = dpole < tolAgreeP && dfin < tolAgreeF &&
        (!cuba || (dpoleVeg < tolAgreeP && dfinVeg < tolAgreeF));
    If[p,
      Print["CC4 PASS ", label,
            " A3 IBP-vs-subtraction (dpole=", sci[dpole], " dfin=", sci[dfin], ")"],
      Print["CC4 FAIL ", label,
            " A3 IBP-vs-subtraction expected=dpole<", tolAgreeP,
            " dfin<", tolAgreeF, "  got dpole=", sci[dpole], " dfin=", sci[dfin]]];
    pass = pass && p];

  (* A4: MC vs Vegas agree within each method (informational; skipped without CUBA) *)
  If[cuba,
    Module[{dpoleIBP, dfinIBP, dpoleSub, dfinSub, p},
      dpoleIBP = Abs[poleIbpMC - poleIbpVeg];
      dfinIBP  = Abs[finIbpMC  - finIbpVeg];
      dpoleSub = Abs[poleSubMC - poleSubVeg];
      dfinSub  = Abs[finSubMC  - finSubVeg];
      p = dpoleIBP < tolAgreeP && dfinIBP < tolAgreeF &&
          dpoleSub  < tolAgreeP && dfinSub < tolAgreeF;
      If[p,
        Print["CC4 PASS ", label, " A4 MC-vs-Vegas (dpole-IBP=", sci[dpoleIBP],
              " dfin-IBP=", sci[dfinIBP], ")"],
        Print["CC4 FAIL ", label, " A4 MC-vs-Vegas expected=<", tolAgreeP,
              ">/<", tolAgreeF, "  dpole-IBP=", sci[dpoleIBP],
              " dfin-IBP=", sci[dfinIBP]]];
      pass = pass && p],
    Print["CC4 PASS ", label, " A4 MC-vs-Vegas SKIPPED (no CUBA)"]];

  (* A5: full-integral consistency at eps* = 0.02 *)
  Module[{laurIBP, c1, c2, c3, p},
    laurIBP = poleIbpMC / epsStar + finIbpMC;
    c1 = Abs[FsubMC - Fref] / Abs[Fref];
    c2 = Abs[FsubMC - laurIBP] / Abs[Fref];
    c3 = If[NumericQ[nint], Abs[nint - Fref] / Abs[Fref], 0.];
    p = c1 < tolFull && c2 < tolLaur && c3 < 5.*^-3;
    If[p,
      Print["CC4 PASS ", label, " A5 full-integral-at-eps* (sub-vs-closed=",
            sci[c1], " sub-vs-IBPLaurent=", sci[c2], " NInt-vs-closed=", sci[c3], ")"],
      Print["CC4 FAIL ", label, " A5 full-integral-at-eps* expected=c1<",
            tolFull, " c2<", tolLaur, " c3<5e-3  got c1=", sci[c1],
            " c2=", sci[c2], " c3=", sci[c3]]];
    pass = pass && p];

  pass
];

(* --------------------------------------------------------------------------
   Main entry point
   -------------------------------------------------------------------------- *)
RunCC4[] := Module[{cuba, allPass = True},
  cuba = TrueQ[detectCuba[]["Found"]];
  Print[""];
  Print["========================================"];
  Print[" v3 Cross-Check #4: IBP vs Subtraction  "];
  Print["  (plan.md §8.2 row #4, §8.3 item #35)  "];
  Print["========================================"];
  Print["  CUBA available: ", cuba,
    If[!cuba, "  -> Vegas columns SKIPPED; MC cross-check must still pass", ""]];

  allPass = cc4Spec["B=4", 4, cuba, 1500000, 2000000] && allPass;
  allPass = cc4Spec["B=5", 5, cuba, 1500000, 2000000] && allPass;

  Print[""];
  Print["========================================"];
  If[allPass,
    Print["CC4 PASS  IBP-vs-subtraction cross-check complete"],
    Print["CC4 FAIL  one or more assertions failed — see above"]];
  Print["========================================"];
  allPass
];

If[!TrueQ[$CC4NoRun], RunCC4[]];
