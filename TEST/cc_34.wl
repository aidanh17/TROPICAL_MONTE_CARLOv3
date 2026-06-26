(* ============================================================================
   TEST/cc_34.wl  —  Cross-check #34 (plan.md §8.3)
                     Lifting × complex coeff × VEGAS, high-D (≈8D)

   What this tests
   ---------------
   The "F3×F4×F5" feature crossing: tropical lifting of a tiny complex
   coefficient (F3b, plan.md §8.1) combined with complex polynomial exponents
   (F4a) and high-dimensional VEGAS integration (F5c), in exactly the regime
   that produced BUG-1 / L1 / L2 / L3 in the bug logs.

   Integrand (8D, from phase5_8d.wl oracle):
       P = 1 + c x[1]^2 + x[1] + Sum_{i=2}^{8} x[i]^2,  c = (1+I) 1e-6
       B = -6  (real exponent, complex coefficient in P)

   Independent oracle (Schwinger WP=40 + CUBA 1e9, confirmed in phase5_8d.wl):
       Re = 0.00317086,  Im ≈ 0  (negligible for real B)

   Five sub-checks are run in order:

     A  K-scaled fan builds without error (computeFanScaled, plan.md §6.4).
        DetectExtremeCoefficients returns rules with "SuggestedK" key.
        LiftCoefficients (post-fix) accepts "SuggestedK" via ruleK[], so the
        lifted Newton polytope vertex carries an integer exponent rather than
        Missing[KeyAbsent,k].  computeFanScaled then builds the (n+1)-dim fan.
        Asserts: AssociationQ return, LiftedSpec key present, polytope vertices
        list non-empty, fan returns {dualVerts, simplexList} with sector count > 0,
        and no Missing value in any lifted polynomial exponent vector.

     B  SplitRealImag VEGAS vs oracle: rel-err < 0.5% (plan.md §8.3 #34).
        EvaluateTropicalMCLifted with ComplexExponentMode -> "SplitRealImag".
        FanData -> Automatic lets the wrapper build the K-scaled fan.

     C  Direct VEGAS vs oracle: rel-err < 0.5%.
        EvaluateTropicalMCLifted with ComplexExponentMode -> "Direct".
        Must independently agree with both oracle and B.

     D  Direct ≡ SplitRealImag within combined VEGAS error + 0.5% oracle.
        The two modes must be mutually consistent (not just individually vs oracle).
        This gates plan.md §6.3 Fix A / BUG-1: without the complex-log fix,
        SplitRealImag silently dropped exp(-Im(B)*arg P) and would disagree
        with Direct for complex-coefficient polynomials where arg P ≠ 0.

     E  vegasbudget fires when NSamples < guardFactor*NStart (plan.md §6.5).
        For d=9 (the lifted n+1 dimension) the guard factor is 100; with
        NSamples=1000 and VegasNStart=100 the resolved NStart=100 and
        threshold=10000, so NSamples=1000 < 10000 fires TropicalEval::vegasbudget.

   PASS criterion summary
   ----------------------
     A  fan built; Length[sectors] > 0; no Missing in lifted exponents
     B  rel-err(SplitRealImag VEGAS vs oracle) < 5e-3
     C  rel-err(Direct VEGAS vs oracle) < 5e-3
     D  |SplitRealImag - Direct| < combined_sigma + 5e-3 * oracle
     E  vegasbudget message fires (handler or StringContains fallback)

   Tier 2: requires CUBA (for the VEGAS integration).
   Sub-checks B, C, D, E require CUBA; sub-check A is Tier 1.
   If CUBA is absent, prints "CC34 SKIP (no CUBA)" and exits 0.

   Engine fix (BUG #34): prior to the Phase-6 engine fix, LiftCoefficients
   read r["k"] on rules produced by DetectExtremeCoefficients (which carry
   "SuggestedK", not "k"), getting Missing[KeyAbsent,k].  The lifted Newton
   polytope vertex then carried Missing[KeyAbsent,k] as an exponent, making
   computeFanScaled return $Failed.  The fix in LiftCoefficients accepts both
   "k" and "SuggestedK" via ruleK[].  After the fix all 5 sub-checks produce
   numbers and the K-scaled fan builds correctly.

   Source / port notes
   -------------------
   #34 is a genuinely new v3 cross-check (plan.md §8.3 "Build-first" set).
   The integrand and oracle are ported from TEST/phase5_8d.wl.
   The SplitRealImag vs Direct comparison mirrors cc_21b.wl sub-checks A+B.
   The vegasbudget sub-check mirrors the budgetFired probe in phase5_8d.wl,
   using the same message-match pattern as that file.
   The K-scaled fan check mirrors phase5_lowd.wl and plan.md §6.4.

   Run
   ---
     wolframscript -file TEST/cc_34.wl
   (from any directory; absolute paths are derived from $InputFileName)
   ============================================================================ *)


(* ── 0.  Locate package root from this file's own path ───────────────────── *)

$pkgRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
$evalWL  = FileNameJoin[{$pkgRoot, "tropical_eval.wl"}];
$fanWL   = FileNameJoin[{$pkgRoot, "tropical_fan.wl"}];
$ioDir   = FileNameJoin[{$pkgRoot, "TEST", "INTERFILES", "cc34"}];

If[!FileExistsQ[$evalWL],
  Print["CC34 FAIL  tropical_eval.wl not found at ", $evalWL];
  Quit[1]];

Get[$fanWL];
Get[$evalWL];

Quiet[CreateDirectory[$ioDir, CreateIntermediateDirectories -> True],
      {CreateDirectory::eexist}];


(* ── 1.  CUBA presence check (Tier-2 gate) ───────────────────────────────── *)

$cubaInfo = TropicalEval`detectCuba[];
If[!TrueQ[$cubaInfo["Found"]],
  Print["CC34 SKIP (no CUBA)"];
  Quit[0]];

Print["CC34: tropical_eval.wl loaded; CUBA found at ", $cubaInfo["IncludeDir"]];
Print["CC34: INTERFILES directory: ", $ioDir];
Print[];


(* ── 2.  Shared helpers ──────────────────────────────────────────────────── *)

(* Relative error, real-valued (|got - ref| / |ref|; fallback to abs if tiny). *)
relErr[got_, ref_] :=
  If[Abs[N[ref]] < 1*^-300,
     Abs[N[got - ref]],
     Abs[N[got - ref]] / Abs[N[ref]]];

$passes   = {};
$failures = {};

recordPass[name_String, detail_String : ""] :=
  (AppendTo[$passes, name];
   Print["CC34 PASS  ", name, If[detail =!= "", "  " <> detail, ""]]);

recordFail[name_String, extra_String : ""] :=
  (AppendTo[$failures, name];
   Print["CC34 FAIL  ", name, If[extra =!= "", "  " <> extra, ""]]);

record[name_String, ok_, passDetail_String : "", failExtra_String : ""] :=
  If[TrueQ[ok], recordPass[name, passDetail], recordFail[name, failExtra]];


(* ── 3.  Integrand and oracle ─────────────────────────────────────────────── *)
(* 8D integrand from phase5_8d.wl.  c = (1+I)*1e-6 is the tiny complex
   coefficient that drives lifting (F3b).  B = -6 is real, so the integrand's
   imaginary part is negligible; we gate only on Re.
   The complex coefficient means arg(P) != 0 in general, which is the condition
   that exposed BUG-1 in the bug logs (plan.md §6.3 Fix A). *)

$n    = 8;
$vars = Table[x[i], {i, $n}];
$cc   = (1 + I) 10^-6;
$poly = 1 + $cc x[1]^2 + x[1] + Sum[x[i]^2, {i, 2, $n}];
$Bex  = -6;

$spec = <|
  "Polynomials"         -> {$poly},
  "MonomialExponents"   -> ConstantArray[0, $n],
  "PolynomialExponents" -> {$Bex},
  "Variables"           -> $vars,
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> None
|>;

(* Oracle from phase5_8d.wl (Schwinger WP=40 + CUBA 1e9 independent reference).
   Im is numerically 0 for real B; we check only Re. *)
$oracleRe = 0.00317086`;

Print["CC34: 8D integrand  c = ", $cc, "  B = ", $Bex];
Print["CC34: oracle  Re = ", $oracleRe, "  (Schwinger WP=40, phase5_8d.wl)"];
Print["CC34: This exercises BUG-1 / BUG-34 fix regime: arg(P) != 0 for complex c"];
Print[];


(* ── 4A.  K-scaled fan builds (computeFanScaled, plan.md §6.4) ───────────── *)
(* Standalone sub-check: detect extreme coefficients, lift the spec, build the
   (n+1)-dim Newton polytope vertices, and call computeFanScaled.

   This is the non-vacuous exercise of the BUG #34 fix: DetectExtremeCoefficients
   emits rules with "SuggestedK"; LiftCoefficients post-fix accepts "SuggestedK"
   via ruleK[], so the lifted spec carries integer exponents (not Missing values)
   and computeFanScaled can build the fan.  Before the fix, this path returned
   $Failed because Missing[KeyAbsent,k] poisoned the aux-variable exponent. *)

Print["=== Sub-check A: K-scaled fan construction (BUG-34 fix path) ==="];

(* Step 1: detect extreme coefficients — returns rules with "SuggestedK" key. *)
$detectRules = TropicalEval`DetectExtremeCoefficients[$spec, 1000];
Print["  DetectExtremeCoefficients: ", Length[$detectRules], " rules found"];
If[Length[$detectRules] > 0,
  Print["  First rule keys: ", Keys[$detectRules[[1]]]]];

(* Step 2: lift the spec directly using DetectExtremeCoefficients rules as-is
   (with "SuggestedK" key, not "k").  This is the exact path that BUG #34
   broke: LiftCoefficients must accept "SuggestedK" via ruleK[]. *)
$liftResult = TropicalEval`LiftCoefficients[$spec, $detectRules];
$liftOK     = AssociationQ[$liftResult] && KeyExistsQ[$liftResult, "LiftedSpec"];

If[$liftOK,
  $liftedSpec = $liftResult["LiftedSpec"];
  $liftData   = $liftResult["LiftData"];
  $auxVar     = $liftData["AuxVariable"];
  Print["  LiftCoefficients OK; aux variable = ", $auxVar],
  Print["  LiftCoefficients returned: ", $liftResult];
  Print["  (FAIL: before fix this returned $Failed due to Missing[KeyAbsent,k])"]
];

(* Step 3: check that no Missing value appears in the lifted polynomial
   exponent structure — the specific signature of BUG #34. *)
$noMissing = If[$liftOK,
  FreeQ[$liftedSpec["Polynomials"], _Missing],
  False];
If[$liftOK,
  Print["  Lifted polynomials free of Missing values: ", $noMissing]];

(* Step 4: Newton polytope vertices of the lifted integrand (n+1 = 9 dims). *)
$liftedVerts = If[$liftOK,
  Quiet[
    TropicalFan`PolytopeVertices[
      (Times @@ $liftedSpec["Polynomials"])^(-1),
      $liftedSpec["Variables"]],
    {TropicalFan::polymake}],
  $Failed];

$vertsOK = ListQ[$liftedVerts] && Length[$liftedVerts] > 0;
Print["  PolytopeVertices: ", If[$vertsOK, Length[$liftedVerts], 0], " vertices; ",
      If[$vertsOK,
         "ambient dim = " <> ToString[Length[$liftedVerts[[1]]]] <>
           " (expected " <> ToString[$n + 1] <> ")",
         "FAILED"]];

(* Step 5: K-scaled fan (plan.md §6.4 / lift_error_log L2).
   Uses TropicalFan`computeFanScaled which retries ComputeDecomposition at
   K in {1, n+2, 2n+4, 6n+6} to handle thin lattice simplices in dim >= 4. *)
$liftedFan = If[$vertsOK,
  Quiet[TropicalFan`computeFanScaled[$liftedVerts], {TropicalFan::polymake}],
  $Failed];

$fanOK    = ListQ[$liftedFan] && Length[$liftedFan] >= 2;
$nSectors = If[$fanOK, Length[$liftedFan[[2]]], 0];

If[$fanOK,
  Print["  computeFanScaled OK: ", $nSectors, " sectors in 9D fan"],
  Print["  computeFanScaled returned: ", $liftedFan]];

record["A K-scaled fan builds (post BUG-34 fix); no Missing in lifted exponents",
  $liftOK && $noMissing && $fanOK && $nSectors > 0,
  "sectors=" <> ToString[$nSectors] <> " noMissing=" <> ToString[$noMissing],
  "liftOK=" <> ToString[$liftOK] <>
    " noMissing=" <> ToString[$noMissing] <>
    " fanOK=" <> ToString[$fanOK] <>
    " sectors=" <> ToString[$nSectors] <> " (expected > 0)"];
Print[];


(* ── 4B.  SplitRealImag VEGAS vs oracle ─────────────────────────────────── *)
(* EvaluateTropicalMCLifted with FanData -> Automatic so the wrapper performs
   the full lifting + K-scaled fan construction internally, exercising the
   end-to-end lifted-VEGAS path.  ComplexExponentMode -> "SplitRealImag"
   decomposes on Re(B) and reintroduces Im(B) as an oscillatory phase
   (plan.md §6.3).  Here B=-6 is real, so Im(B)=0 and SplitRealImag ≡ Direct
   at the codegen level — but the polynomial P has complex coefficients (arg P != 0),
   which is what exposed BUG-1 (the complex-log fix). *)

Print["=== Sub-check B: SplitRealImag VEGAS vs oracle ==="];
$wdB = FileNameJoin[{$ioDir, "B_split_vegas"}];
Quiet[CreateDirectory[$wdB, CreateIntermediateDirectories -> True],
      {CreateDirectory::eexist}];

$resB = TropicalEval`EvaluateTropicalMCLifted[$spec, {{}},
  "Integrator"          -> "VEGAS",
  "NSamples"            -> 4000000,
  "RunChecks"           -> False,
  "Verbose"             -> False,
  "WorkingDirectory"    -> $wdB,
  "FanData"             -> Automatic,
  "ComplexExponentMode" -> "SplitRealImag"];

$resB_ok = AssociationQ[$resB] && KeyExistsQ[$resB, "Results"] &&
           Length[$resB["Results"]] >= 1;
$vBRe  = If[$resB_ok, $resB["Results"][[1]]["Re"],  Missing["noResult"]];
$vBIm  = If[$resB_ok, $resB["Results"][[1]]["Im"],  Missing["noResult"]];
$eBRe  = If[$resB_ok, $resB["Results"][[1]]["ReErr"], 0.];
$rBoracle = If[$resB_ok && NumericQ[$vBRe], relErr[$vBRe, $oracleRe], 999.];

Print["  SplitRealImag VEGAS: Re = ", $vBRe, " +/- ", $eBRe,
      "  Im = ", $vBIm];
Print["  rel-err vs oracle:   ", ScientificForm[$rBoracle, 3],
      "  (gate: < 5e-3)"];
record["B SplitRealImag VEGAS vs oracle",
  NumericQ[$vBRe] && $rBoracle < 5*^-3,
  "expected=<5e-3 got=" <> ToString[ScientificForm[$rBoracle, 3]],
  "expected=<5e-3 got=" <> ToString[ScientificForm[$rBoracle, 3]]];
Print[];


(* ── 4C.  Direct VEGAS vs oracle ─────────────────────────────────────────── *)
(* ComplexExponentMode -> "Direct" uses exp(B*log P) with complex log.
   For real B=-6 this is numerically equivalent to SplitRealImag (since Im(B)=0
   makes the oscillatory-phase term vanish).  Cross-checking both modes against
   the oracle independently validates the codegen paths for this lifted case. *)

Print["=== Sub-check C: Direct VEGAS vs oracle ==="];
$wdC = FileNameJoin[{$ioDir, "C_direct_vegas"}];
Quiet[CreateDirectory[$wdC, CreateIntermediateDirectories -> True],
      {CreateDirectory::eexist}];

$resC = TropicalEval`EvaluateTropicalMCLifted[$spec, {{}},
  "Integrator"          -> "VEGAS",
  "NSamples"            -> 4000000,
  "RunChecks"           -> False,
  "Verbose"             -> False,
  "WorkingDirectory"    -> $wdC,
  "FanData"             -> Automatic,
  "ComplexExponentMode" -> "Direct"];

$resC_ok = AssociationQ[$resC] && KeyExistsQ[$resC, "Results"] &&
           Length[$resC["Results"]] >= 1;
$vCRe  = If[$resC_ok, $resC["Results"][[1]]["Re"],  Missing["noResult"]];
$vCIm  = If[$resC_ok, $resC["Results"][[1]]["Im"],  Missing["noResult"]];
$eCRe  = If[$resC_ok, $resC["Results"][[1]]["ReErr"], 0.];
$rCoracle = If[$resC_ok && NumericQ[$vCRe], relErr[$vCRe, $oracleRe], 999.];

Print["  Direct VEGAS: Re = ", $vCRe, " +/- ", $eCRe,
      "  Im = ", $vCIm];
Print["  rel-err vs oracle:  ", ScientificForm[$rCoracle, 3],
      "  (gate: < 5e-3)"];
record["C Direct VEGAS vs oracle",
  NumericQ[$vCRe] && $rCoracle < 5*^-3,
  "expected=<5e-3 got=" <> ToString[ScientificForm[$rCoracle, 3]],
  "expected=<5e-3 got=" <> ToString[ScientificForm[$rCoracle, 3]]];
Print[];


(* ── 4D.  Direct ≡ SplitRealImag within combined error ───────────────────── *)
(* Mutual consistency check (plan.md §8.3 #34: "Direct ≡ SplitRealImag ≡ ref").
   For real B the two modes should agree to VEGAS statistical precision.
   This sub-check would catch a systematic bias in either codegen path. *)

Print["=== Sub-check D: Direct ≡ SplitRealImag (mutual consistency) ==="];
$modeDiff      = If[NumericQ[$vBRe] && NumericQ[$vCRe],
                    Abs[$vBRe - $vCRe], 999.];
$sigmaCombined = Abs[$eBRe] + Abs[$eCRe];
$tol           = $sigmaCombined + 5*^-3 * Abs[$oracleRe];
$modeOK        = $modeDiff < $tol;

Print["  |SplitRealImag - Direct| = ", ScientificForm[$modeDiff, 3]];
Print["  tolerance (|err_B| + |err_C| + 0.5% oracle) = ",
      ScientificForm[$tol, 3]];
record["D Direct == SplitRealImag within combined error",
  NumericQ[$vBRe] && NumericQ[$vCRe] && $modeOK,
  "|diff|=" <> ToString[ScientificForm[$modeDiff, 3]] <>
    " < tol=" <> ToString[ScientificForm[$tol, 3]],
  "|diff|=" <> ToString[ScientificForm[$modeDiff, 3]] <>
    " >= tol=" <> ToString[ScientificForm[$tol, 3]]];
Print[];


(* ── 4E.  vegasbudget fires when under-budgeted ──────────────────────────── *)
(* plan.md §6.5 / lift_error_log L1: TropicalEval::vegasbudget fires when
   NSamples < guardFactor * resolvedNStart.  For lifted 9D sectors (d >= 5),
   guardFactor = 100.  With NSamples=1000 and VegasNStart=100 (explicit, so
   resolveVegasSizing returns 100 verbatim), threshold = 100*100 = 10000.
   NSamples=1000 < 10000 -> message fires.

   Message-match pattern mirrors phase5_8d.wl (the working reference):
     Hold[Message[MessageName[TropicalEval, "vegasbudget"], ___], _]
   where TropicalEval resolves to TropicalEval`TropicalEval after Get.
   A StringContainsQ fallback catches any context-path variation.             *)

Print["=== Sub-check E: vegasbudget fires when NSamples=1000 < 10000 (9D) ==="];
$wdE = FileNameJoin[{$ioDir, "E_budget_probe"}];
Quiet[CreateDirectory[$wdE, CreateIntermediateDirectories -> True],
      {CreateDirectory::eexist}];

$budgetFired = False;

Internal`HandlerBlock[
  {"Message",
   Function[{msgHeld},
     If[
       (* Primary pattern: matches MessageName[TropicalEval, "vegasbudget"]
          (same as phase5_8d.wl — works because TropicalEval is in context
          after Get["tropical_eval.wl"]). *)
       MatchQ[msgHeld,
         Hold[Message[MessageName[TropicalEval`TropicalEval, "vegasbudget"], ___],
              _]] ||
       (* Secondary: fully-qualified context form. *)
       MatchQ[msgHeld,
         Hold[Message[
           HoldPattern[MessageName[TropicalEval`TropicalEval, "vegasbudget"]],
           ___], _]] ||
       (* Fallback: string-match on the tag (catches re-exported symbols). *)
       StringContainsQ[ToString[msgHeld], "vegasbudget"],
       $budgetFired = True
     ]
   ]},
  Quiet[
    TropicalEval`EvaluateTropicalMCLifted[$spec, {{}},
      "Integrator"          -> "VEGAS",
      "NSamples"            -> 1000,
      "VegasNStart"         -> 100,
      "RunChecks"           -> False,
      "Verbose"             -> False,
      "WorkingDirectory"    -> $wdE,
      "FanData"             -> Automatic,
      "ComplexExponentMode" -> "Direct"],
    (* Quiet all messages EXCEPT vegasbudget — but the handler fires first,
       so Quiet only affects message printing, not handler interception. *)
    All]
];

(* Fallback probe: check "VegasBudgetWarning" or "Warnings" key in return value
   in case the message symbol lives in an unexpected context. *)
If[!$budgetFired,
  Module[{resE},
    resE = Quiet[
      TropicalEval`EvaluateTropicalMCLifted[$spec, {{}},
        "Integrator"          -> "VEGAS",
        "NSamples"            -> 1000,
        "VegasNStart"         -> 100,
        "RunChecks"           -> False,
        "Verbose"             -> False,
        "WorkingDirectory"    -> $wdE,
        "FanData"             -> Automatic,
        "ComplexExponentMode" -> "Direct"],
      All];
    If[AssociationQ[resE] &&
       (TrueQ[Lookup[resE, "VegasBudgetWarning", False]] ||
        (KeyExistsQ[resE, "Warnings"] &&
         AnyTrue[resE["Warnings"], StringContainsQ[#, "budget"] &])),
      $budgetFired = True]
  ]
];

Print["  NSamples=1000, VegasNStart=100 -> resolvedNStart=100, ",
      "threshold=10000 (guardFactor=100 for d>=5)"];
Print["  vegasbudget message detected: ", $budgetFired];
record["E vegasbudget fires for NSamples=1000 < 10000 (lifted 9D, guardFactor=100)",
  $budgetFired,
  "fired=True",
  "expected=True got=False  " <>
    "(check resolveVegasSizing for 9D and guardFactor logic in EvaluateTropicalMC)"];
Print[];


(* ── 5.  Final summary ───────────────────────────────────────────────────── *)

$nPass = Length[$passes];
$nFail = Length[$failures];
$nTot  = $nPass + $nFail;

Print["================================================================"];
If[$nFail == 0,
  Print["CC34 PASS  all ", $nTot, " sub-checks passed",
        "  (lifting x complex-coeff x VEGAS 8D: Direct/Split agree; ",
        "K-fan OK; vegasbudget fires)"],
  Print["CC34 FAIL  ", $nFail, " of ", $nTot, " sub-checks failed:  ",
        StringRiffle[$failures, ", "]]];
Print["================================================================"];

Quit[If[$nFail > 0, 1, 0]]
