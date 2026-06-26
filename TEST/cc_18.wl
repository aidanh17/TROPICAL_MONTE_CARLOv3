(* ============================================================================
   TEST/cc_18.wl  —  v3 Cross-Check #18
   plan.md §8.2 row #18:  Lifted vs unlifted variance
   VarRed = (sigma_U / sigma_L)^2.

   §8.2 PASS criterion:
     Case A  (1 + 10^6 x1^2 + x2^2 + x1 x2^2)^-2 :
       Lifted sectors must be correctly classified.  The v3 HONEST bench
       (SANDBOX/bench_lift_variance.wl) showed that Case A has a
       HasConstantTerm=False sector -> divergent I2 -> INFINITE true variance,
       so the class-(ii) classification and the WARNING fire are the PASS gates.
     Case C  (1 + 10^8 x1^3 x2 + x2^3)^-3  (degenerate lifted polytope):
       ALL lifted sectors HasConstantTerm=False -> correctly NOT beneficial.
       Automatic fan detection fires liftdegenerate (geometry gate refuses).
     Case D  (1 + 10^6 x1)^-2  (1D, B=-2):
       Class-(i) finite-variance demo: all sectors HasConstantTerm=True, I2
       converges, exact VarRed >= 1 (genuine improvement).

   Tier 2  (CUBA required for the VEGAS-lifted vs MC-lifted parity sub-gate).
   If CUBA is absent this script prints "CC18 SKIP (no CUBA)" and exits 0.

   KNOWN BUG FIXED: The earlier version had a Set::shape destructuring bug in
   trueSigmaOne — it tried to use FlattenedPolys as {expr, exp} pairs, but
   FlattenedPolys is a list of polynomials each being a list of {coeff, expVec}
   monomials.  The fix reconstructs the polynomial value from the monomial list
   and reads PolynomialExponents from the sector Association (Phase-3 honest
   variance gate: exact trueSigma, never sampled).

   Ported / adapted from
     SANDBOX/bench_lift_variance.wl       (v3 honest trueSigma bench)
     OLD_CODE/.../SANDBOX/bench_lift_variance.wl (Old v2 bench — used for
       the spec definitions and fan data; logic replaced by honest trueSigma)
     OLD_CODE/.../EXAMPLES/test_vegas.wl  (L1/L2 lifted+VEGAS vs lifted+MC)

   Run from any directory:
     wolframscript -file /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3/TEST/cc_18.wl
   ============================================================================ *)

(* --------------------------------------------------------------------------
   Load v3 packages by absolute path.
   -------------------------------------------------------------------------- *)
$v3Root = "/Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3";
Get[FileNameJoin[{$v3Root, "tropical_fan.wl"}]];
Get[FileNameJoin[{$v3Root, "tropical_eval.wl"}]];

(* --------------------------------------------------------------------------
   CUBA presence check (same dirs as the package's detectCuba).
   Tier 2: if CUBA absent, skip cleanly.
   -------------------------------------------------------------------------- *)
$cubaPrefix = SelectFirst[
  {"/opt/homebrew", "/usr/local", "/usr"},
  FileExistsQ[FileNameJoin[{#, "include", "cuba.h"}]] &&
  (FileExistsQ[FileNameJoin[{#, "lib", "libcuba.a"}]] ||
   FileExistsQ[FileNameJoin[{#, "lib", "libcuba.dylib"}]] ||
   FileExistsQ[FileNameJoin[{#, "lib", "libcuba.so"}]]) &, $Failed];

If[$cubaPrefix === $Failed,
  Print["CC18 SKIP (no CUBA)"];
  Quit[0]
];

Print["CC18: CUBA found at ", $cubaPrefix];
Print["CC18: tropical_eval.wl loaded"];
Print[];

(* --------------------------------------------------------------------------
   Working directory for C++ codegen artifacts.
   -------------------------------------------------------------------------- *)
$wd = FileNameJoin[{$v3Root, "INTERFILES", "cc18"}];
If[!DirectoryQ[$wd], CreateDirectory[$wd, CreateIntermediateDirectories -> True]];

(* --------------------------------------------------------------------------
   Failure accumulator.
   -------------------------------------------------------------------------- *)
$cc18Failures = {};
cc18Fail[msg_String] := (
  Print["  CC18 GATE-FAIL: ", msg];
  AppendTo[$cc18Failures, msg]
);

(* --------------------------------------------------------------------------
   Numeric formatter (no stray backticks from InputForm on reals).
   -------------------------------------------------------------------------- *)
fmtN[x_] := Which[
  x === Infinity,      "INFINITE",
  x === Indeterminate, "Indeterminate",
  NumericQ[x],         ToString[N[x, 5], InputForm],
  True,                ToString[x, InputForm]
];


(* ==========================================================================
   HONEST trueSigma infrastructure
   Mirrors SANDBOX/bench_lift_variance.wl exactSigmaTotal / buildUnlifted /
   buildLifted, self-contained here so cc_18.wl runs independently.

   CRITICAL FIX (Set::shape bug):
   FlattenedPolys from ProcessSector / ProcessSectorLifted is a list of
   polynomials, each polynomial being a list of {coeff, {e1,...,en}} monomials:
     flatPolys[[j]] = { {c1, {e1,...}}, {c2, {e1,...}}, ... }
   The polynomial VALUE at point yv is computed by summing over monomials:
     Sum[ mono[[1]] * Exp[Total[mono[[2]] * Log[yv]]] , {mono, flatPolys[[j]]} ]
   The exponent for polynomial j is PolynomialExponents[[j]] from the sector.
   The old (buggy) code tried Power[#[[1]], #[[2]]] treating flatPolys as
   {expr, exp} pairs, causing Set::shape errors and making all I2 diverge.
   ========================================================================== *)

(* PrecisionGoal for the exact-sigma NIntegrate (I1, I2).
   HONEST gate (Phase-3): the reference sigma must be CONVERGED, never sampled
   and never a spurious-timeout->Infinity.  A genuinely-divergent I2 (the
   HasConstantTerm=False sectors of Case A / Case C) is detected by NIntegrate's
   OWN non-convergence messages (ncvb/slwcon/eincr/inumr/...), passed to Check as
   SEPARATE message-name arguments (NOT an Alternatives expression, which left the
   old Check unevaluated -> non-numeric -> every sector read INFINITE, even the
   trivially-finite Case D).  Convergent sectors finish in << 1 s; the long safety
   TimeConstrained below therefore never fires on a finite-variance case. *)
$tsPG = 6;

(* -----------------------------------------------------------------------
   trueSigmaOne[sd, pg]:
   Compute exact per-sample sigma for one sector sd via NIntegrate:
     sigma^2 = I2 - I1^2,  I1 = Int g,  I2 = Int g^2
   Works for both lifted and unlifted sector associations.
   EmptyDomain / IsDivergent sectors return zero contribution.
   ----------------------------------------------------------------------- *)
trueSigmaOne[sd_Association, pg_Integer] :=
Module[
  {flatPolys, polyExps, pre, dim, domConstr, yVars,
   polyVals, integrand, integrand2,
   logZ0, mp, indCoeffs, logYpStar,
   i1, i2result, i2, i2cvgd},

  (* Skip empty / divergent sectors *)
  If[TrueQ[Lookup[sd, "EmptyDomain", False]] ||
     TrueQ[Lookup[sd, "IsDivergent",  False]],
    Return[<|"I1" -> 0, "I2" -> 0, "Sigma" -> 0, "I2Converged" -> True|>]
  ];

  flatPolys = Lookup[sd, "FlattenedPolys", {}];
  polyExps  = Lookup[sd, "PolynomialExponents", {}];
  pre       = Lookup[sd, "Prefactor", 1];
  dim       = Lookup[sd, "Dimension", Length[flatPolys]];
  domConstr = Lookup[sd, "DomainConstraint", None];

  If[!ListQ[flatPolys] || Length[flatPolys] == 0,
    Return[<|"I1" -> 0, "I2" -> 0, "Sigma" -> 0, "I2Converged" -> False|>]
  ];

  (* Fresh integration variables — one per dimension (same pattern as engine's ValidateDecomposition) *)
  yVars = Table[Unique["cc18yv"], {dim}];

  (* Build polynomial values from the monomial list representation.
     FlattenedPolys[[j]] = { {coeff, {e1,...,en}}, ... } — each entry is a
     {coefficient, exponent-vector} pair.  Evaluate the polynomial sum at yVars
     using the log-exp pattern (same as trueSigma in sandbox_lift_common.wl):
       P_j(y) = Sum_m  coeff_m * Exp[ Sum_i alpha_{m,i} * Log[y_i] ]  *)
  polyVals = Table[
    Total[Table[
      mono[[1]] * Exp[Total[mono[[2]] * (Log /@ yVars)]],
      {mono, flatPolys[[j]]}
    ]],
    {j, Length[flatPolys]}
  ];

  (* integrand = pre * Prod_j P_j(y)^{polyExps[[j]]} *)
  integrand = pre * Times @@ MapThread[
    Function[{pv, be}, Exp[be * Log[pv]]],
    {polyVals, polyExps}
  ];

  (* Domain indicator for lifted sectors *)
  If[domConstr =!= None,
    logZ0     = domConstr["LogZ0"];
    mp        = domConstr["MP"];
    indCoeffs = domConstr["IndicatorCoeffs"];
    logYpStar = (logZ0 - Total[indCoeffs * (Log /@ yVars)]) / mp;
    integrand = integrand * Boole[logYpStar <= 0]
  ];

  integrand2 = integrand^2;

  (* I1 = Int g  (converged adaptive NIntegrate). *)
  i1 = Quiet @ NIntegrate[
    integrand,
    Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
    MaxRecursion -> 30, PrecisionGoal -> pg,
    Method -> {"GlobalAdaptive", "SingularityHandler" -> "IMT"}
  ];

  (* I2 = Int g^2.  HONEST divergence detection: a genuinely-divergent I2 (the
     HasConstantTerm=False sectors of Case A / Case C) makes NIntegrate emit its
     OWN non-convergence messages (ncvb/slwcon/eincr/inumr), which we hand to
     Check as SEPARATE message-name arguments.  This is the bug fix: the prior
     version OR-ed the tags into a single Alternatives ("a|b|c"), which is invalid
     Check syntax — Check stayed UNEVALUATED, so i2result was never numeric and
     EVERY sector (including the trivially-finite Case A unlifted and Case D) was
     mis-classified divergent -> INFINITE.  No TimeConstrained->Infinity coercion:
     a finite-variance sector converges in << 1 s and is detected as finite; a
     genuinely-divergent sector is detected by the messages above, not by a clock.
     A long safety TimeConstrained guards only a true hang. *)
  i2result = Quiet @ Check[
    TimeConstrained[
      NIntegrate[
        integrand2,
        Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
        MaxRecursion -> 30, PrecisionGoal -> pg,
        Method -> {"GlobalAdaptive", "SingularityHandler" -> "IMT"}
      ],
      300,        (* generous safety net; genuine divergence is caught first *)
      $TimedOut
    ],
    $Failed,      (* failure value when a non-convergence message fires *)
    NIntegrate::ncvb,  NIntegrate::slwcon, NIntegrate::eincr,
    NIntegrate::inumr, NIntegrate::izero,  NIntegrate::deltam
  ];

  Which[
    i2result === $Failed || i2result === $TimedOut,
      i2cvgd = False; i2 = Infinity,
    !NumericQ[i2result],
      i2cvgd = False; i2 = Infinity,
    Abs[i2result] > 10^30,
      i2cvgd = False; i2 = Infinity,
    True,
      i2cvgd = True; i2 = i2result
  ];

  <|"I1"          -> i1,
    "I2"          -> i2,
    "Sigma"       -> If[i2cvgd, Sqrt[Max[0, i2 - i1^2]], Infinity],
    "I2Converged" -> i2cvgd
  |>
];


(* Total exact sigma over a list of sector associations (lifted OR unlifted). *)
exactSigmaTotal[sectorList_List] :=
Module[{ts, i1sum, divIdx, sigSq, sigTot},
  ts = Table[
    If[AssociationQ[sd] &&
       !TrueQ[Lookup[sd, "EmptyDomain", False]] &&
       !TrueQ[Lookup[sd, "IsDivergent",  False]],
      trueSigmaOne[sd, $tsPG],
      <|"I1" -> 0, "I2" -> 0, "Sigma" -> 0, "I2Converged" -> True|>
    ],
    {sd, sectorList}
  ];
  i1sum  = Total[Map[#["I1"] &, ts]];
  divIdx = Flatten @ Position[Map[TrueQ[#["I2Converged"]] &, ts], False];
  If[divIdx === {},
    sigSq  = Total[Map[#["Sigma"]^2 &, ts]];
    sigTot = Sqrt[sigSq],
    sigSq  = Infinity;
    sigTot = Infinity
  ];
  <|"I1Sum"        -> i1sum,
    "SigmaSqTotal" -> sigSq,
    "SigmaTotal"   -> sigTot,
    "AnyDivergent" -> (divIdx =!= {}),
    "DivergentIdx" -> divIdx,
    "PerSector"    -> ts|>
];


(* Build unlifted sector list via unit-coefficient proxy fan. *)
buildUnlifted[spec_Association, proxyPoly_] :=
Module[{verts, fan},
  verts = PolytopeVertices[proxyPoly^(-1), spec["Variables"]];
  fan   = ComputeDecomposition[verts, "ShowProgress" -> False];
  <|"Sectors"  -> Table[ProcessSector[spec, fan[[1]], fan[[2, s]], s],
                        {s, Length[fan[[2]]]}],
    "NSectors" -> Length[fan[[2]]]|>
];


(* Build lifted sector list via LiftCoefficients + ProcessSectorLifted.
   fanOpt = Automatic derives the fan from the lifted Newton polytope. *)
buildLifted[spec_Association, liftRules_, fanOpt_] :=
Module[{lc, ls, ld, fan, sectors, hasConstList, nDropped},
  lc = LiftCoefficients[spec, liftRules];
  If[!AssociationQ[lc],
    Return[<|"Failed" -> True, "Reason" -> lc|>]];
  ls = lc["LiftedSpec"]; ld = lc["LiftData"];
  fan = If[fanOpt =!= Automatic, fanOpt,
    Module[{vv = Quiet @ PolytopeVertices[
      (Times @@ ls["Polynomials"])^(-1), ls["Variables"]]},
      Quiet @ ComputeDecomposition[vv, "ShowProgress" -> False]]];
  If[!(ListQ[fan] && Length[fan] >= 2),
    Return[<|"Failed" -> True, "Reason" -> "fan-build", "Fan" -> fan|>]];
  sectors = Table[
    Quiet @ ProcessSectorLifted[ls, fan[[1]], fan[[2, s]], s, ld],
    {s, Length[fan[[2]]]}];
  hasConstList = Cases[sectors,
    a_?AssociationQ /; !TrueQ[Lookup[a, "EmptyDomain", False]] :>
      Lookup[a, "HasConstantTerm", Missing[]]];
  nDropped = Count[sectors,
    a_?AssociationQ /; TrueQ[Lookup[a, "EmptyDomain", False]]];
  <|"Failed"      -> False,
    "Sectors"     -> sectors,
    "LiftedSpec"  -> ls,
    "LiftData"    -> ld,
    "HasConstList"-> hasConstList,
    "NDropped"    -> nDropped,
    "z0"          -> ld["z0"]|>
];


(* Warning tracker (mirrors SANDBOX/bench_lift_variance.wl classify[]). *)
$classWarned = {};
classify[label_String, esU_Association, esL_Association, hcFalse_] :=
Module[{vr, cls, divergent},
  divergent = TrueQ[esL["AnyDivergent"]];
  If[divergent,
    cls = "INFINITE-VARIANCE / NOT-BENEFICIAL (class ii)";
    vr  = Infinity;
    Print["  *** WARNING [", label, "]: divergent I2 in sector(s) ",
          esL["DivergentIdx"], " -> INFINITE true variance",
          If[hcFalse, " (HasConstantTerm=False; plan §9 risk #1 / L5)", ""], ". ***"];
    AppendTo[$classWarned, label],
    (* class (i): finite I2, genuine VarRed *)
    cls = "FINITE-VARIANCE (class i: all I2 converge)";
    vr  = If[NumericQ[esU["SigmaTotal"]] && NumericQ[esL["SigmaTotal"]] &&
             esL["SigmaTotal"] > 0,
             (esU["SigmaTotal"]/esL["SigmaTotal"])^2,
             Indeterminate]
  ];
  Print["  ", label, ":"];
  Print["     trueSigma_unlifted = ", fmtN[esU["SigmaTotal"]],
        "   trueSigma_lifted = ", fmtN[esL["SigmaTotal"]]];
  Print["     class: ", cls];
  Print["     EXACT VarRed (trueSigma_U^2 / trueSigma_L^2) = ", fmtN[vr]];
  <|"class" -> cls, "vr" -> vr, "hcFalse" -> hcFalse, "divergent" -> divergent|>
];


(* ==========================================================================
   CASE A:  P = 1 + 10^6 x1^2 + x2^2 + x1 x2^2,  B={-2}.
   Lift rule: explicit k=2 on {2,0}.
   Expected: HasConstantTerm=False sector -> INFINITE true variance (class ii).
             HasConstantTerm=False warning must fire.
   ========================================================================== *)
Print["================================================================="];
Print["CC18 CASE A:  P = 1 + 10^6 x1^2 + x2^2 + x1 x2^2,  B={-2}"];
Print["             lift {2,0} k=2"];
Print["================================================================="];

$specA = <|"Polynomials"        -> {1 + 10^6 x[1]^2 + x[2]^2 + x[1] x[2]^2},
            "MonomialExponents"  -> {0, 0},
            "PolynomialExponents"-> {-2},
            "Variables"          -> {x[1], x[2]},
            "KinematicSymbols"   -> {},
            "RegulatorSymbol"    -> None|>;
$ruleA = {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 2|>};

(* Converged NIntegrate reference for the direct integrand. *)
$refA = Quiet @ NIntegrate[
  1/(1 + 10^6 t1^2 + t2^2 + t1 t2^2)^2,
  {t1, 0, Infinity}, {t2, 0, Infinity},
  MaxRecursion -> 30, PrecisionGoal -> 8, WorkingPrecision -> 30];
Print["  Converged NIntegrate reference = ", N[$refA, 8]];

$uA   = buildUnlifted[$specA, 1 + x[1]^2 + x[2]^2 + x[1] x[2]^2];
$esUA = exactSigmaTotal[$uA["Sectors"]];
Print["  Unlifted: ", $uA["NSectors"], " sectors, I1=", fmtN[$esUA["I1Sum"]],
      "  trueSigma=", fmtN[$esUA["SigmaTotal"]]];

$lA = buildLifted[$specA, $ruleA, Automatic];
If[TrueQ[$lA["Failed"]],
  cc18Fail["Case A lifted build failed: " <> ToString[$lA["Reason"]]];
  $esLA = <|"I1Sum" -> 0, "SigmaTotal" -> Infinity, "AnyDivergent" -> True,
             "DivergentIdx" -> {}, "PerSector" -> {}|>,
  $esLA = exactSigmaTotal[$lA["Sectors"]];
  Print["  Lifted:   z0=", $lA["z0"], "  sectors=", Length[$lA["Sectors"]],
        "  dropped=", $lA["NDropped"]];
  Print["  HasConstantTerm per non-dropped sector: ", $lA["HasConstList"]];
  Print["  I1 sum=", fmtN[$esLA["I1Sum"]], "  trueSigma=", fmtN[$esLA["SigmaTotal"]]];
];

$hcFalseA = MemberQ[$lA["HasConstList"], False];
$cA = classify["Case A (lift {2,0} k=2)", $esUA, $esLA, $hcFalseA];
Print[];


(* ==========================================================================
   CASE C:  P = 1 + 10^8 x1^3 x2 + x2^3,  B={-3}.  Degenerate lifted polytope.
   Expected:
     (a) ALL lifted sectors HasConstantTerm=False -> NOT beneficial.
     (b) Automatic anchor selection refuses (geometry gate -> liftdegenerate).
   ========================================================================== *)
Print["================================================================="];
Print["CC18 CASE C:  P = 1 + 10^8 x1^3 x2 + x2^3,  B={-3}"];
Print["             degenerate lifted polytope (3 monomials in 3D)"];
Print["================================================================="];

$specC = <|"Polynomials"        -> {1 + 10^8 x[1]^3 x[2] + x[2]^3},
            "MonomialExponents"  -> {0, 0},
            "PolynomialExponents"-> {-3},
            "Variables"          -> {x[1], x[2]},
            "KinematicSymbols"   -> {},
            "RegulatorSymbol"    -> None|>;

$refC = Quiet @ NIntegrate[1/(1 + 10^8 t1^3 t2 + t2^3)^3,
  {t1, 0, Infinity}, {t2, 0, Infinity},
  MaxRecursion -> 30, PrecisionGoal -> 8, WorkingPrecision -> 30];
Print["  Converged NIntegrate reference = ", N[$refC, 8]];

$uC   = buildUnlifted[$specC, 1 + x[1]^3 x[2] + x[2]^3];
$esUC = exactSigmaTotal[$uC["Sectors"]];
Print["  Unlifted: ", $uC["NSectors"], " sectors, I1=", fmtN[$esUC["I1Sum"]],
      "  trueSigma=", fmtN[$esUC["SigmaTotal"]]];

(* (a) geometry-gate: Automatic should refuse (liftdegenerate -> $Failed). *)
$autoResC = Quiet @ EvaluateTropicalMCLifted[$specC, {{}},
  "LiftRules" -> Automatic, "NSamples" -> 20000, "RunChecks" -> False,
  "Verbose" -> False, "WorkingDirectory" -> $wd];
$autoNoLiftC = !(AssociationQ[$autoResC] && KeyExistsQ[$autoResC, "Results"]);
Print["  Geometry gate (Automatic): ",
      If[$autoNoLiftC,
         "refuses to lift (liftdegenerate -> $Failed) [EXPECTED]",
         "UNEXPECTEDLY produced a lifted result"]];

(* (b) explicit-fan lifts — both k=2 and k=4 should have all HasConstantTerm=False. *)
$dvCk2 = {{1, 0, 0}, {0, 1, 0}, {-1, -3, 0}, {2, 0, -3}, {-2, 0, 3}};
$slC   = {{1, 2, 4}, {2, 3, 4}, {3, 1, 4}, {1, 2, 5}, {2, 3, 5}, {3, 1, 5}};
$dvCk4 = {{1, 0, 0}, {0, 1, 0}, {-1, -3, 0}, {4, 0, -3}, {-4, 0, 3}};

$lCk2 = buildLifted[$specC, {<|"PolyIndex"->1,"ExponentVector"->{3,1},"k"->2|>},
           {$dvCk2, $slC}];
$lCk4 = buildLifted[$specC, {<|"PolyIndex"->1,"ExponentVector"->{3,1},"k"->4|>},
           {$dvCk4, $slC}];

$hcListCk2 = If[TrueQ[$lCk2["Failed"]], {}, $lCk2["HasConstList"]];
$hcListCk4 = If[TrueQ[$lCk4["Failed"]], {}, $lCk4["HasConstList"]];

$esLCk2 = If[TrueQ[$lCk2["Failed"]],
  <|"SigmaTotal"->Infinity,"AnyDivergent"->True,"DivergentIdx"->{},"PerSector"->{},"I1Sum"->0|>,
  exactSigmaTotal[$lCk2["Sectors"]]];
$esLCk4 = If[TrueQ[$lCk4["Failed"]],
  <|"SigmaTotal"->Infinity,"AnyDivergent"->True,"DivergentIdx"->{},"PerSector"->{},"I1Sum"->0|>,
  exactSigmaTotal[$lCk4["Sectors"]]];

Print["  k=2 HasConstantTerm per sector = ", $hcListCk2];
Print["  k=4 HasConstantTerm per sector = ", $hcListCk4];

$hcFalseC  = MemberQ[$hcListCk2, False] || MemberQ[$hcListCk4, False];
$allFalseC = ($hcListCk2 =!= {} && AllTrue[$hcListCk2, # === False &]) ||
             ($hcListCk4 =!= {} && AllTrue[$hcListCk4, # === False &]);
$esLCbest  = If[$esLCk2["SigmaTotal"] === Infinity, $esLCk2, $esLCk4];
$cC = classify["Case C (explicit-fan lift)", $esUC, $esLCbest, $hcFalseC];
Print[];


(* ==========================================================================
   CASE D  (class-(i) FINITE-VARIANCE demonstration):
   P = 1 + 10^6 x1,  B={-2},  1D.  Lift {1} k=2 (z0 = 10^3).
   Expected: all HasConstantTerm=True, all I2 converge, exact VarRed >= 1.
   This is the hand-checkable Toy-0 example from SANDBOX/bench_lift_variance.wl.
   ========================================================================== *)
Print["================================================================="];
Print["CC18 CASE D:  P = 1 + 10^6 x1,  B={-2}  (1D, lift {1} k=2)"];
Print["             Finite-variance demo — class-(i) gate"];
Print["================================================================="];

$specD = <|"Polynomials"        -> {1 + 10^6 x[1]},
            "MonomialExponents"  -> {0},
            "PolynomialExponents"-> {-2},
            "Variables"          -> {x[1]},
            "KinematicSymbols"   -> {},
            "RegulatorSymbol"    -> None|>;
$ruleD  = {<|"PolyIndex" -> 1, "ExponentVector" -> {1}, "k" -> 2|>};
$exactD = 10^-6;

(* k=2 explicit fan: complete unimodular triangulation of R^2 for support
   {(0,0),(1,2)} — the 1D-lift fan from SANDBOX/bench_lift_variance.wl. *)
$raysD  = {{1, 0}, {0, 1}, {-1, 2}, {-1, 0}, {-1, -2}, {0, -1}, {1, -2}};
$sectsD = Table[{i, If[i < Length[$raysD], i+1, 1]}, {i, Length[$raysD]}];

$uD   = buildUnlifted[$specD, 1 + x[1]];
$esUD = exactSigmaTotal[$uD["Sectors"]];
Print["  Unlifted: ", $uD["NSectors"], " sectors, I1=", fmtN[$esUD["I1Sum"]],
      "  (exact ", N[$exactD], ", relErr=",
      fmtN[Abs[($esUD["I1Sum"] - $exactD)/$exactD]], ")",
      "  trueSigma=", fmtN[$esUD["SigmaTotal"]]];

$lD = buildLifted[$specD, $ruleD, {$raysD, $sectsD}];
If[TrueQ[$lD["Failed"]],
  cc18Fail["Case D lifted build failed: " <> ToString[$lD["Reason"]]];
  $esLD = <|"I1Sum"->0,"SigmaTotal"->Infinity,"AnyDivergent"->True,
             "DivergentIdx"->{},"PerSector"->{}|>,
  $esLD = exactSigmaTotal[$lD["Sectors"]];
  Print["  Lifted:   z0=", $lD["z0"], "  sectors=", Length[$lD["Sectors"]],
        "  dropped=", $lD["NDropped"]];
  Print["  HasConstantTerm per non-dropped sector: ", $lD["HasConstList"]];
  Print["  I1 sum=", fmtN[$esLD["I1Sum"]],
        "  (exact ", N[$exactD], ", relErr=",
        fmtN[Abs[($esLD["I1Sum"] - $exactD)/$exactD]], ")",
        "  trueSigma=", fmtN[$esLD["SigmaTotal"]]];
];

$hcFalseD = MemberQ[$lD["HasConstList"], False];
$vrD = If[NumericQ[$esUD["SigmaTotal"]] && NumericQ[$esLD["SigmaTotal"]] &&
          $esLD["SigmaTotal"] > 0,
         ($esUD["SigmaTotal"]/$esLD["SigmaTotal"])^2, Indeterminate];
Print["  EXACT VarRed (trueSigma_U^2/trueSigma_L^2) = ", fmtN[$vrD]];
$cD = classify["Case D (1D 1+10^6 x1, B=-2, lift k=2)", $esUD, $esLD, $hcFalseD];
Print[];


(* ==========================================================================
   VEGAS-LIFTED sub-gate (Tier-2, ported from OLD .../EXAMPLES/test_vegas.wl
   cases L1 and L2):
   Run EvaluateTropicalMCLifted with Integrator->"VEGAS" and compare against
   the NIntegrate reference.  L1 (large coeff 10^6): VEGAS codegen must be
   reached and VEGAS result within 5% of reference.  L2 (small coeff 10^-4):
   VEGAS must also be within 2% and be more accurate than lifted MC.
   ========================================================================== *)
Print["================================================================="];
Print["CC18 VEGAS sub-gate  (cases L1 and L2 from test_vegas.wl)"];
Print["================================================================="];

(* Robust NIntegrate reference (avoids HoldAll issues). *)
niRef2D[t1_, t2_, poly_] := Module[{r},
  r = Quiet @ NIntegrate[poly,
    {t1, 0, Infinity}, {t2, 0, Infinity},
    MaxRecursion -> 40, PrecisionGoal -> 8, WorkingPrecision -> 30];
  If[NumericQ[r], N[r], $Failed]];

runLiftedIntegrator[spec_, liftRules_, integrator_, neval_, label_] :=
Module[{res, val, err},
  Print["  ", label, " (", integrator, ", neval=", neval, ") ..."];
  res = Quiet @ EvaluateTropicalMCLifted[spec, {{}},
    "LiftRules"        -> liftRules,
    "Integrator"       -> integrator,
    "NSamples"         -> neval,
    "VegasEpsRel"      -> 1.*^-9,
    "RunChecks"        -> False,
    "Verbose"          -> False,
    "WorkingDirectory" -> $wd];
  If[AssociationQ[res] && KeyExistsQ[res, "Results"],
    val = res["Results"][[1]]["Re"];
    err = res["Results"][[1]]["ReErr"];
    Print["    ", integrator, " result = ", val, " +/- ", err];
    <|"value" -> val, "err" -> err, "ok" -> True|>,
    Print["    ", integrator, " run FAILED: ", res];
    <|"ok" -> False|>
  ]
];

(* L1: P = 1 + 10^6 x1^2 + x2^2 + x1 x2^2,  B=-2,  lift {2,0} k=2. *)
Print["--- L1: large coeff 10^6 ---"];
$specL1 = $specA;  (* same spec as Case A *)
$ruleL1 = $ruleA;
$refL1  = $refA;
Print["  NIntegrate reference = ", N[$refL1, 7]];

$rVEL1 = runLiftedIntegrator[$specL1, $ruleL1, "VEGAS", 500000, "L1"];
$rMCL1 = runLiftedIntegrator[$specL1, $ruleL1, "MC",    500000, "L1"];

$vegasOkL1 = TrueQ[$rVEL1["ok"]];
If[$vegasOkL1,
  $devVEL1 = Abs[$rVEL1["value"]/$refL1 - 1];
  $devMCL1 = If[TrueQ[$rMCL1["ok"]], Abs[$rMCL1["value"]/$refL1 - 1], Indeterminate];
  Print["  L1 relErr: VEGAS=", fmtN[$devVEL1], "  MC=", fmtN[$devMCL1],
        "   (plan: VEGAS codegen reached + within 5%)"];
  If[!($devVEL1 < 0.05),
    cc18Fail["L1 VEGAS relErr = " <> fmtN[$devVEL1] <> " >= 0.05"]],
  cc18Fail["L1 VEGAS run failed"]
];
Print[];

(* L2: P = 1 + 10^-4 x1^2 + x2^2 + x1 x2^2,  B=-2,  lift {2,0} k=2. *)
Print["--- L2: small coeff 10^-4 ---"];
$specL2 = <|"Polynomials"        -> {1 + 10^-4 x[1]^2 + x[2]^2 + x[1] x[2]^2},
             "MonomialExponents"  -> {0, 0},
             "PolynomialExponents"-> {-2},
             "Variables"          -> {x[1], x[2]},
             "KinematicSymbols"   -> {},
             "RegulatorSymbol"    -> None|>;
$ruleL2 = {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 2|>};
$refL2  = niRef2D[t1L2, t2L2, (1 + 10^-4 t1L2^2 + t2L2^2 + t1L2 t2L2^2)^(-2)];
If[$refL2 === $Failed,
  cc18Fail["L2 NIntegrate reference failed"],
  Print["  NIntegrate reference = ", N[$refL2, 7]];
  $rVEL2 = runLiftedIntegrator[$specL2, $ruleL2, "VEGAS", 500000, "L2"];
  $rMCL2 = runLiftedIntegrator[$specL2, $ruleL2, "MC",    500000, "L2"];
  If[TrueQ[$rVEL2["ok"]],
    $devVEL2 = Abs[$rVEL2["value"]/$refL2 - 1];
    $devMCL2 = If[TrueQ[$rMCL2["ok"]], Abs[$rMCL2["value"]/$refL2 - 1], Indeterminate];
    Print["  L2 relErr: VEGAS=", fmtN[$devVEL2], "  MC=", fmtN[$devMCL2],
          "   (plan: VEGAS within 2% AND more accurate than MC)"];
    If[!($devVEL2 < 0.02),
      cc18Fail["L2 VEGAS relErr = " <> fmtN[$devVEL2] <> " >= 0.02"]];
    If[TrueQ[$rMCL2["ok"]] && !($devVEL2 < $devMCL2),
      cc18Fail["L2 VEGAS not more accurate than lifted MC (devVEGAS=" <>
               fmtN[$devVEL2] <> " devMC=" <> fmtN[$devMCL2] <> ")"]],
    cc18Fail["L2 VEGAS run failed"]
  ]
];
Print[];


(* ==========================================================================
   GATE ASSERTIONS
   ========================================================================== *)
Print["================================================================="];
Print["  CC18 GATE ASSERTIONS"];
Print["================================================================="];

(* G1: I1 sums reproduce converged reference. *)
$relAlifted   = If[NumericQ[$esLA["I1Sum"]] && NumericQ[$refA],
                   Abs[($esLA["I1Sum"] - $refA)/$refA], Indeterminate];
$relAunlifted = If[NumericQ[$esUA["I1Sum"]] && NumericQ[$refA],
                   Abs[($esUA["I1Sum"] - $refA)/$refA], Indeterminate];
$g1ok = (NumericQ[$relAlifted]   && $relAlifted   < 0.01) &&
        (NumericQ[$relAunlifted] && $relAunlifted < 0.01);
If[!$g1ok,
  cc18Fail["Case A I1 off reference: lifted relErr=" <> fmtN[$relAlifted] <>
           " unlifted relErr=" <> fmtN[$relAunlifted]]];
Print["  (G1) Case A I1 reproduces NIntegrate ref (lifted relErr=",
      fmtN[$relAlifted], ", unlifted relErr=", fmtN[$relAunlifted], "): ",
      If[$g1ok, "PASS", "FAIL"]];

(* G2: Case A must be class-(ii): HasConstantTerm=False + divergent I2. *)
$g2ok = $hcFalseA && TrueQ[$esLA["AnyDivergent"]];
If[!$g2ok,
  cc18Fail["Case A expected class-(ii) (HasConstantTerm=False + divergent I2)"]];
Print["  (G2) Case A class-(ii): HasConstantTerm=False warning + INFINITE true variance: ",
      If[$g2ok, "PASS", "FAIL"],
      "  (divergent sectors=", $esLA["DivergentIdx"], ")"];

(* G3: Case A headline VarRed = INFINITE (not a sampled 18x pass). *)
$g3ok = ($cA["vr"] === Infinity);
If[!$g3ok, cc18Fail["Case A VarRed must be INFINITE, got " <> fmtN[$cA["vr"]]]];
Print["  (G3) Case A headline VarRed = INFINITE (not a sampled-sigma pass): ",
      If[$g3ok, "PASS", "FAIL"]];

(* G4: Case C all lifted sectors HasConstantTerm=False (NOT beneficial). *)
$g4ok = $allFalseC && $hcFalseC;
If[!$g4ok, cc18Fail["Case C expected all lifted sectors HasConstantTerm=False"]];
Print["  (G4) Case C all lifted sectors HasConstantTerm=False (NOT beneficial): ",
      If[$g4ok, "PASS", "FAIL"]];

(* G5: Geometry gate refuses Automatic lift for Case C. *)
$g5ok = $autoNoLiftC;
If[!$g5ok, cc18Fail["Case C geometry gate should refuse Automatic lift"]];
Print["  (G5) Case C geometry-gate refuses Automatic lift (liftdegenerate): ",
      If[$g5ok, "PASS", "FAIL"]];

(* G6: Case D class-(i) finite VarRed: all HasConstantTerm=True, VarRed >= 1. *)
$g6ok = !$hcFalseD && !TrueQ[$esLD["AnyDivergent"]] &&
        NumericQ[$vrD] && $vrD >= 1 &&
        NumericQ[$esLD["SigmaTotal"]] && NumericQ[$esUD["SigmaTotal"]];
If[!$g6ok,
  cc18Fail["Case D finite-variance demo failed: hcFalseD=" <> ToString[$hcFalseD] <>
           " divergent=" <> ToString[TrueQ[$esLD["AnyDivergent"]]] <>
           " vrD=" <> fmtN[$vrD]]];
Print["  (G6) Case D class-(i) genuine finite VarRed=", fmtN[$vrD],
      " (all HasConstantTerm=True, all I2 converge, VarRed>=1): ",
      If[$g6ok, "PASS", "FAIL"]];

(* G7: Warning invariant — every divergent scheme warned; Case D did not. *)
$divergentSchemes = {"Case A (lift {2,0} k=2)", "Case C (explicit-fan lift)"};
$warnOk = AllTrue[$divergentSchemes, MemberQ[$classWarned, #] &] &&
          !MemberQ[$classWarned, "Case D (1D 1+10^6 x1, B=-2, lift k=2)"];
If[!$warnOk,
  cc18Fail["Warning invariant: warned=" <> ToString[$classWarned] <>
           " divergent=" <> ToString[$divergentSchemes]]];
Print["  (G7) WARNING fired on every divergent scheme (A, C) and NOT on finite D: ",
      If[$warnOk, "PASS", "FAIL"]];

(* G8: Case D I1 sum close to exact 10^-6. *)
$relDlifted   = If[NumericQ[$esLD["I1Sum"]],
                   Abs[($esLD["I1Sum"] - $exactD)/$exactD], Indeterminate];
$relDunlifted = If[NumericQ[$esUD["I1Sum"]],
                   Abs[($esUD["I1Sum"] - $exactD)/$exactD], Indeterminate];
$g8ok = NumericQ[$relDlifted]   && $relDlifted   < 0.01 &&
        NumericQ[$relDunlifted] && $relDunlifted < 0.01;
If[!$g8ok,
  cc18Fail["Case D I1 off exact: lifted relErr=" <> fmtN[$relDlifted] <>
           " unlifted relErr=" <> fmtN[$relDunlifted]]];
Print["  (G8) Case D I1 reproduces exact 10^-6 (lifted relErr=",
      fmtN[$relDlifted], ", unlifted relErr=", fmtN[$relDunlifted], "): ",
      If[$g8ok, "PASS", "FAIL"]];

Print[];

(* ==========================================================================
   FINAL VERDICT
   ========================================================================== *)
$cc18AllPass = Length[$cc18Failures] == 0;

If[$cc18AllPass,
  Print["CC18 PASS  {caseA=class-ii-infinite, caseC=not-beneficial, caseD=vrD=",
        fmtN[$vrD], ", L1-VEGAS-relErr=",
        If[TrueQ[$rVEL1["ok"]], fmtN[$devVEL1], "N/A"],
        ", L2-VEGAS-relErr=",
        If[TrueQ[$rVEL2["ok"]], fmtN[$devVEL2], "N/A"], "}"],
  Print["CC18 FAIL  expected=<all gates PASS>  got=<",
        Length[$cc18Failures], " failure(s): ",
        StringRiffle[$cc18Failures, "; "], ">"];
  Quit[1]
];
