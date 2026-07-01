(* ============================================================================
   TEST/cc_54.wl  —  Cross-check #54  (planLIFTDIVZ0.md)

   PURPOSE.  planLIFTDIVZ0.md was written from an external `bug_summary` that
   claimed the lifted x divergent x IBP x SplitRealImag path drops a complex
   `z0^{~1.5}` factor in the value assembly (a silent WRONG NUMBER).  This
   cross-check is the repo-native decider the plan asked for (§3): it drives the
   *exact* implicated combination and pins the result against TWO independent
   oracles.

   THE AXIS cc_53 LACKS (planLIFTDIVZ0.md §1.3).  cc_53 Part B already exercises
   a lifted x divergent x SplitRealImag IBP integral, but with a REAL monomial
   exponent A (Im(A)=0), so its Im(A) z0-phase machinery (LiftedMonoPhase /
   MonomialPhaseLog) is dormant.  This fixture puts a NONZERO Im(A) on a
   convergent slot AND arranges the geometry so Im(A) projects onto the lift
   pivot => LiftedMonoPhase["Const"] = (imEaug_pivot/m_p) log z0 is NONZERO.  So
   the divergent-lifted Im(A) z0-phase term — the leading suspect in the plan —
   is genuinely driven here for the first time.

   OUTCOME (the honest reassessment, planLIFTDIVZ0.md §4.1/§4.4).  The bug is
   NOT reproduced at repo-native fixture sizes: the lifted IBP (pole, finite)
   agrees with BOTH the NIntegrate complex Laurent of the original AND the
   pinned-eps Subtraction route — in magnitude AND phase — to MC precision.  The
   lifted path is CORRECT, including the Im(A) z0-phase.  What `bug_summary`
   measured (lift != no-lift, "|ratio| ~ z0^1.5") is the §4.1 trap: its no-lift
   REFERENCE is the unsound one.  The unlifted decomposition on an integrand with
   an extreme (~1e4) coefficient has blown-up MC variance — precisely the
   hierarchy lifting exists to tame — so a naive lift-vs-no-lift MC identity
   false-fails on unlifted noise, not a lifted z0 drop.  (Verified in-session:
   with the hierarchy removed the unlifted complex-A/B divergent IBP matches
   NIntegrate to <0.3%; the lifted result matches NIntegrate regardless of the
   hierarchy.)  Hence this test pins the LIFTED result against oracles that do
   NOT depend on the variance-limited unlifted path, and it PASSES on the current
   engine (6020ef9).  It stands as the regression guard for the correctness of
   the lifted x divergent x complex-A/B x SplitRealImag IBP assembly.

   Parts:
     (A) UNIT — build the lifted divergent sectors and assert the fixture really
         exercises the implicated code: >=1 divergent cone that is a single
         on-axis Pole (theta=0) with ZERO off-axis directions (matching
         bug_summary point 2), carrying a NONZERO LiftedMonoPhase["Const"]
         (Im(A) z0-phase live), with z0 != 1 (a genuine lift).
     (B) END-TO-END — lifted IBP (pole, finite) vs NIntegrate complex Laurent
         (external oracle) AND vs the pinned-eps Subtraction route (2nd
         independent engine, invariant #4).  Agreement in magnitude AND phase.

   Oracles per numeric PASS (invariant #4): NIntegrate of the complex original
   (Laurent fit) + the pinned-eps Subtraction route (shares no MC assembly with
   IBP).  Two independent engines + one external analytic reference.

   Tier 1: WL + g++ (MC).  Run:  wolframscript -file TEST/cc_54.wl
   ============================================================================ *)

$pkgRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{$pkgRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$pkgRoot, "tropical_eval.wl"}]];
Print["CC54: packages loaded."];
Print[];

sci[v_] := If[NumericQ[v], ScientificForm[N[v], 3], v];
$cc54Pass = True; $cc54Fail = {};
cc54Assert[label_String, test_, exp_String, got_String] :=
  If[TrueQ[test],
    Print["CC54 PASS  [", label, "]  expected=", exp, "  got=", got],
    Print["CC54 FAIL  [", label, "]  expected=", exp, "  got=", got];
    $cc54Pass = False; AppendTo[$cc54Fail, label]];

SetAttributes[qr, HoldFirst];
qr[e_] := Quiet[e, {TropicalEval::validate, TropicalEval::divergent, TropicalEval::vegasbudget,
   TropicalEval::liftcomplex, TropicalEval::splitdivmono, TropicalEval::liftnopivot,
   TropicalEval::liftdivdomain, General::stop, NIntegrate::slwcon,
   NIntegrate::ncvb, NIntegrate::inumr, NIntegrate::eincr, NIntegrate::izero, NIntegrate::precw,
   TropicalFan::polymake, General::munfl, Power::infy, Infinity::indet}];

wdir[tag_] := Module[{d = FileNameJoin[{$pkgRoot, "TEST", "INTERFILES", "cc54", tag}]},
  Quiet[CreateDirectory[d, CreateIntermediateDirectories -> True], {CreateDirectory::eexist}]; d];

(* ----------------------------------------------------------------------------
   Shared fixture.  A genuine single 1/eps pole (A_1 = -1 + eps), an EXTREME
   coefficient (1e4 on x2^2) lifted with k=3 => z0 = 10^(4/3) != 1, a complex
   polynomial exponent with the regulator inside (B = -2 + i/2 - eps, exercises
   the 6020ef9 regulator-real fix), AND — the new axis — a nonzero Im(A) = 1/2 on
   a genuinely CONVERGENT slot (Re(A_2)=1/2>0, so no oscillatory endpoint and no
   off-axis pole).  The {0,2} lift makes Im(A) project onto the pivot.
   ---------------------------------------------------------------------------- *)
epsD  = Symbol["epsCC54"];
dim   = 2;
gB    = 1/2;                                  (* Im(B) *)
qA    = 1/2;                                  (* Im(A) on the convergent slot *)
esV   = {0.01, 0.02, 0.04, 0.08};
polysF = {1 + x[1] x[2] + x[1]^2 + 10^4 x[2]^2};
AvalsF = {-1 + epsD, 1/2 + I qA};
peF    = {-2 + I gB - epsD};
specF  = <|"Polynomials" -> polysF, "MonomialExponents" -> AvalsF,
           "PolynomialExponents" -> peF, "Variables" -> {x[1], x[2]},
           "KinematicSymbols" -> {}, "RegulatorSymbol" -> epsD|>;
liftRuleF = {<|"PolyIndex" -> 1, "ExponentVector" -> {0, 2}, "k" -> 3|>};

(* ============================================================================
   PART A — UNIT: the fixture really drives the implicated code
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (A) UNIT: lifted divergent single-pole cone with LIVE Im(A) z0-phase"];

Module[{lift, ls, ld, imB, imA, hasImA, lsReal, verts, fan, dv, sl, divInfo, z0},
  lift = qr@LiftCoefficients[specF, liftRuleF];
  cc54Assert["A LiftCoefficients succeeds", AssociationQ[lift],
    "Association", ToString[Head[lift]]];
  ls = lift["LiftedSpec"]; ld = lift["LiftData"]; z0 = ld["z0"];
  cc54Assert["A z0 != 1 (genuine lift of an extreme coefficient)",
    FreeQ[z0, _Real] && !TrueQ[z0 == 1] && TrueQ[Simplify[z0 == 10^(4/3)]],
    "10^(4/3)", ToString[z0, InputForm]];

  imB = TropicalEval`Private`imagPolyInfo[ls];
  {imA, hasImA} = TropicalEval`Private`imagMonoInfo[ls, epsD];
  cc54Assert["A Im(A) is present (the axis cc_53 lacks)",
    TrueQ[hasImA] && AnyTrue[imA, (!TrueQ[PossibleZeroQ[#]]) &],
    "hasImagA=True", ToString[{hasImA, imA}, InputForm]];
  lsReal = TropicalEval`Private`realifyMonoA[
             TropicalEval`Private`realifyPolyB[ls, imB], imA];

  verts = qr@PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]];
  fan   = qr@computeFanScaled[verts];
  {dv, sl} = fan;

  divInfo = {};
  Do[
    Module[{sd, cls, pd, od},
      sd = qr@ProcessSectorLifted[lsReal, dv, sl[[c]], c, ld,
             "Eps" -> epsD, "ImagMonoExps" -> imA];
      If[AssociationQ[sd] && TrueQ[Lookup[sd, "IsDivergent", False]],
        cls = Table[TropicalEval`Private`ibpDivClass[sd, i, epsD], {i, sd["Dimension"]}];
        pd  = Select[Range[sd["Dimension"]], cls[[#]] === "Pole" &];
        od  = Select[Range[sd["Dimension"]], cls[[#]] === "OffAxis" &];
        AppendTo[divInfo, <|"cone" -> c, "pole" -> pd, "off" -> od,
          "lmpConst" -> Lookup[sd["LiftedMonoPhase"], "Const", 0]|>]
      ]
    ],
    {c, Length[sl]}
  ];

  cc54Assert["A >=1 divergent lifted cone exists",
    Length[divInfo] >= 1, ">=1", ToString[Length[divInfo]]];
  (* bug_summary point 2: each divergent cone is exactly ONE on-axis Pole,
     ZERO off-axis directions (the standard single-pole lifted IBP path). *)
  cc54Assert["A every divergent cone is single-Pole / zero-OffAxis",
    Length[divInfo] >= 1 &&
      AllTrue[divInfo, (Length[#["pole"]] == 1 && Length[#["off"]] == 0) &],
    "all {1 pole, 0 off}",
    ToString[{#["pole"], #["off"]} & /@ divInfo, InputForm]];
  (* The leading-suspect term (planLIFTDIVZ0.md §1.3) is actually LIVE here. *)
  cc54Assert["A some divergent cone has a NONZERO LiftedMonoPhase Const (Im(A) z0-phase)",
    AnyTrue[divInfo, (!TrueQ[PossibleZeroQ[#["lmpConst"]]]) &],
    "exists nonzero (imEaug_pivot/m_p) log z0",
    ToString[N[#["lmpConst"]] & /@ divInfo]];
];

(* ============================================================================
   PART B — END-TO-END: lifted IBP vs NIntegrate vs Subtraction
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (B) END-TO-END: lifted IBP (pole,finite) vs 2 independent oracles"];

Module[{niC, gv, fitRe, fitIm, poleNI, finNI, ee, lift, ls, ld, fan,
        ibp, poleI, finI, relPn, relFn, sub, poleS, finS, relPs, relFs},
  (* Oracle 1 (external): NIntegrate complex Laurent of the ORIGINAL integrand. *)
  niC[esv_] := Module[{us, subr, jac, igc, Aev, Bev},
    Aev = AvalsF /. epsD -> esv; Bev = peF /. epsD -> esv;
    us = Table[Unique["u"], {dim}]; subr = Table[x[i] -> us[[i]]/(1 - us[[i]]), {i, dim}];
    jac = Times @@ Table[1/(1 - us[[i]])^2, {i, dim}];
    igc = (Times @@ MapThread[#1^#2 &, {Table[x[i], {i, dim}], Aev}] *
           Times @@ MapThread[#1^#2 &, {polysF, Bev}] /. subr) jac;
    qr@NIntegrate[igc, Evaluate[Sequence @@ ({#, 0, 1} & /@ us)], Method -> "GlobalAdaptive",
      PrecisionGoal -> 8, WorkingPrecision -> 40, MaxRecursion -> 80]];
  gv = niC /@ N[esV, 40];
  fitRe = Fit[MapThread[{#1, Re[#2]} &, {N[esV, 40], gv}], {1/ee, 1, ee}, ee];
  fitIm = Fit[MapThread[{#1, Im[#2]} &, {N[esV, 40], gv}], {1/ee, 1, ee}, ee];
  poleNI = Coefficient[fitRe, ee, -1] + I Coefficient[fitIm, ee, -1];
  finNI  = Coefficient[fitRe, ee,  0] + I Coefficient[fitIm, ee,  0];
  Print["   NIntegrate Laurent: pole=", N@poleNI, "  finite=", N@finNI];

  lift = qr@LiftCoefficients[specF, liftRuleF];
  ls = lift["LiftedSpec"]; ld = lift["LiftData"];
  fan = qr@computeFanScaled[qr@PolytopeVertices[
    (Times @@ ls["Polynomials"])^(-1), ls["Variables"]]];

  ibp = qr@EvaluateTropicalMC[ls, fan, {{}}, "LiftData" -> ld, "Method" -> "IBP",
    "ComplexExponentMode" -> "SplitRealImag", "Integrator" -> "MC", "NSamples" -> 4000000,
    "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> wdir["B"]];

  sub = qr@LaurentFromSubtraction[ls, fan, {{}}, "LiftData" -> ld,
    "ComplexExponentMode" -> "SplitRealImag", "EpsilonValues" -> esV[[;; 3]],
    "Integrator" -> "MC", "NSamples" -> 2000000, "Verbose" -> False,
    "WorkingDirectory" -> wdir["Bsub"]];

  If[!(AssociationQ[ibp] && KeyExistsQ[ibp, "Results"]),
    cc54Assert["B IBP compiles + runs", False, "Association[Results]", ToString[Head[ibp]]],
    poleI = ibp["Results"][[1]]["PoleCoefficient"]; finI = ibp["Results"][[1]]["FinitePart"];
    Print["   IBP: pole=", N@poleI, "  finite=", N@finI];
    cc54Assert["B IBP compiles + runs (no $Failed)", True, "Association[Results]", "Association"];

    relPn = Abs[(poleI - poleNI)/poleNI]; relFn = Abs[(finI - finNI)/finNI];
    (* Oracle 2 (independent engine): pinned-eps Subtraction. *)
    If[AssociationQ[sub] && KeyExistsQ[sub, "Results"],
      poleS = sub["Results"][[1]]["Pole"]; finS = sub["Results"][[1]]["Finite"];
      Print["   Sub: pole=", N@poleS, "  finite=", N@finS];
      relPs = Abs[(poleI - poleS)/poleS]; relFs = Abs[(finI - finS)/finS];
      cc54Assert["B IBP pole == Subtraction pole (2nd engine, #4)",
        NumericQ[relPs] && relPs < 0.05, "rel<5e-2", "rel=" <> ToString[sci@relPs]];
      cc54Assert["B IBP finite == Subtraction finite (2nd engine, #4; magnitude+phase)",
        NumericQ[relFs] && relFs < 0.06, "rel<6e-2", "rel=" <> ToString[sci@relFs]],
      cc54Assert["B Subtraction route returns a result (2nd oracle)", False,
        "Association[Results]", ToString[Head[sub]]]];

    (* Oracle 1: NIntegrate Laurent.  Pole well-determined; finite carries a mild
       finite-eps fit bias -> looser tolerance (cf. cc_53). *)
    cc54Assert["B IBP pole vs NIntegrate Laurent",
      NumericQ[relPn] && relPn < 0.05, "rel<5e-2", "rel=" <> ToString[sci@relPn]];
    cc54Assert["B IBP finite vs NIntegrate Laurent (looser: fit bias)",
      NumericQ[relFn] && relFn < 0.08, "rel<8e-2", "rel=" <> ToString[sci@relFn]];
    (* If bug_summary's z0^~1.5 drop were real, relFn would be ~1500%, not <8%. *)
    cc54Assert["B NO dropped z0 factor (finite agrees, not off by z0^~1.5)",
      NumericQ[relFn] && relFn < 0.5, "rel<<z0^1.5-1", "rel=" <> ToString[sci@relFn]];
  ];
];

Print[];
Print["================================================================"];
If[$cc54Pass,
  Print["CC54 PASS  lifted x divergent x complex-A/B x SplitRealImag IBP is ",
        "CORRECT: the Im(A) z0-phase (LiftedMonoPhase, dormant in cc_53) is live ",
        "on a single-pole cone, and the lifted IBP (pole,finite) matches BOTH the ",
        "NIntegrate complex Laurent and the Subtraction route in magnitude AND ",
        "phase.  planLIFTDIVZ0's bug_summary z0^~1.5 drop is NOT reproduced at ",
        "repo-native fixture sizes (§4.1: the no-lift reference it trusted is the ",
        "variance-unsound one).  failures={}"],
  Print["CC54 FAIL  failed: ", $cc54Fail]];
Print["================================================================"];
If[!$cc54Pass, Quit[1]];
