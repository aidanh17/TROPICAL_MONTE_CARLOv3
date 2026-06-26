(* ============================================================================
   TEST/cc_42.wl  —  Cross-check #42 (plan.md §8.3)
                     High-D VEGAS sizing (resolveVegasSizing)

   What this tests
   ---------------
   plan.md §8.3 #42: "makes the Phase-5 sizing fix a named maintained check".
   The three sub-checks mirror the three PASS criteria spelled out in the table:

     A  8D convergent integrand: Automatic sizing reaches ~0.1% accuracy.
        Concrete: an unlifted 8D polynomial integral evaluated with VEGAS under
        Automatic sizing must match the Schwinger/NIntegrate oracle within 0.1%.
        This is the regime where the OLD default {1000,500,1000} returned a
        confidently-wrong value (lift_error_log L1 / phase5_8d.wl measured
        1.2% off at 7sigma with the historical NStart); the fix makes it work.

        Integrand (matching phase5_8d.wl oracle): P = 1 + x[1] + Sum_{i=2}^{8} x[i]^2,
        B = -6.  Oracle Re = 0.00317086 from Schwinger WP=40, confirmed by the
        c=0 closed form (sqrt(pi)/2)^7 * Gamma(3/2)/Gamma(6) = 0.00317087 and
        independent CUBA Vegas 1e9 (phase5_report.txt, re-run verification).

        Bug fix (a): the sum starts at i=2 (not i=1) so x[1] appears only as a
        linear term; the old check had Sum_{i=1}^{8} x[i]^2 which adds an
        x[1]^2 term not present in the oracle integrand, giving true Re~0.00176
        and causing a guaranteed 1-oracle mismatch that refuted the check.

        Bug fix (b): the 8D fan is built via computeFanScaled (K-scaling; the
        normal fan is scale-invariant) rather than a raw ComputeDecomposition
        call that crashes polymake on thin 8D lattice simplices.  Route through
        TropicalFan`computeFanScaled (the public driver path, per §6.4).

        NSamples raised to 4e7 so the budget exceeds the strengthened guard
        threshold (100*NStart = 100*81000 = 8.1e6) and the value is accurate.

     B  Low-D (n=2) unchanged: Automatic resolves to historical {1000,500,1000}.
        Concrete: resolveVegasSizing[Automatic, 2, base, kind] === base for all
        three knobs; VEGAS value matches the exact Gamma-function reference to
        ~1e-4 (byte-identical to pre-fix; ported from phase5_lowd.wl).

     C  vegasbudget fires at NSamples < guardFactor*NStart.
        The engine fires at guardFactor=100 for maxSectorDim>=5 (high-D),
        guardFactor=20 for dim<=4 (low-D unchanged).

        Bug fix (c): the old check tested for the 20x threshold (now only
        active at d<=4); the engine uses 100x for d>=5 so the C1 scenario
        must use NSamples that is below 100*NStart=8.1e6 but above the old
        20*NStart=1.62e6, proving the test is non-vacuous (exercises the
        strengthened rule, not just the historical floor).

        C1  8D integrand, NSamples=5000000.
            Resolved NStart for d=8 is 81000; guardFactor=100 (d>=5);
            threshold = 100*81000 = 8100000.
            5000000 < 8100000 => must fire (was silent under old 20x rule).
        C2  2D integrand, NSamples=1000, explicit NStart=1000.
            guardFactor=20 (d<=4); threshold=20000.
            1000 < 20000 => must fire.

        Interception: Internal`HandlerBlock on "Message" plus a fallback probe
        of a "VegasBudgetWarning" or "Warnings" key in the returned Association.

   All three sub-checks exercise the LIVE resolveVegasSizing function and
   EvaluateTropicalMC with "Integrator" -> "VEGAS".

   Tier 2: requires CUBA.  If absent, prints "CC42 SKIP (no CUBA)" and exits 0.

   Source / port notes
   -------------------
   #42 is a v3-original check (plan.md §8.3 "new checks").
   Integrand and oracle ported from TEST/phase5_8d.wl (Re=0.00317086,
   confirmed by four independent oracles in phase5_report.txt re-run section).
   Low-D sizing + vegasbudget logic ported from TEST/phase5_lowd.wl.
   resolveVegasSizing formula: Ceiling[base * 3^(d-4)] for d>=5, confirmed
   against tropical_eval.wl lines 2077-2083.
   vegasbudget guard: guardFactor=100 for d>=5, 20 for d<=4 (phase5_report.txt
   §ISSUE 1(c); tropical_eval.wl lines 3810-3820).

   Run
   ---
     wolframscript -file TEST/cc_42.wl
   (from any directory; absolute paths derived from $InputFileName)
   ============================================================================ *)


(* ── 0.  Locate package root from this file's absolute path ─────────────── *)

$pkgRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
$evalWL  = FileNameJoin[{$pkgRoot, "tropical_eval.wl"}];
$fanWL   = FileNameJoin[{$pkgRoot, "tropical_fan.wl"}];
$ioDir   = FileNameJoin[{$pkgRoot, "TEST", "INTERFILES", "cc42"}];

If[!FileExistsQ[$evalWL],
  Print["CC42 FAIL  tropical_eval.wl not found at ", $evalWL];
  Quit[1]];

Get[$fanWL];
Get[$evalWL];

Quiet[CreateDirectory[$ioDir, CreateIntermediateDirectories -> True],
      {CreateDirectory::eexist}];


(* ── 1.  CUBA presence check (Tier-2 gate) ───────────────────────────────── *)

$cubaInfo = TropicalEval`detectCuba[];
If[!TrueQ[$cubaInfo["Found"]],
  Print["CC42 SKIP (no CUBA)"];
  Quit[0]];

Print["CC42: packages loaded; CUBA found at ", $cubaInfo["IncludeDir"]];
Print["CC42: INTERFILES directory: ", $ioDir];
Print[];


(* ── 2.  Shared helpers ──────────────────────────────────────────────────── *)

relErr[got_, ref_] :=
  If[Abs[N[ref]] < 1*^-300,
     Abs[N[got - ref]],
     Abs[N[got - ref]] / Abs[N[ref]]];

$passes   = {};
$failures = {};

recordPass[name_String, detail_String : ""] :=
  (AppendTo[$passes, name];
   Print["CC42 PASS  ", name,
         If[detail =!= "", "  " <> detail, ""]]);

recordFail[name_String, detail_String : ""] :=
  (AppendTo[$failures, name];
   Print["CC42 FAIL  ", name,
         If[detail =!= "", "  " <> detail, ""]]);

record[name_String, ok_, passDetail_String : "", failExtra_String : ""] :=
  If[TrueQ[ok],
     recordPass[name, passDetail],
     recordFail[name, failExtra]];


(* ── 3.  Sub-check A: 8D Automatic sizing reaches ~0.1% ─────────────────── *)
(*
   Integrand: P = 1 + x[1] + Sum_{i=2}^{8} x[i]^2,  B = -6  (real exponents,
   no complex coefficient -- isolate the sizing fix, not lifting).
   This matches the phase5_8d.wl oracle integrand (c->0 limit of the complex
   version; x[1] is linear, x[2]..x[8] are quadratic; x[1]^2 NOT in the sum).
   Oracle Re = 0.00317086 from Schwinger WP=40, c=0 closed form
   (sqrt(pi)/2)^7 * Gamma(3/2)/Gamma(6) = 0.003170869, and CUBA Vegas 1e9;
   all three agree to <2e-5 (phase5_report.txt re-run verification section).
   Tolerance: 0.1% = 1e-3 (plan.md #42 PASS criterion).
   Integrator: VEGAS with Automatic sizing and NSamples=4e7 (above the
   strengthened guard threshold 100*81000=8.1e6 so the grid resolves).
*)

Print["=== Sub-check A: 8D Automatic sizing, rel-err < 0.1% ==="];

$n8    = 8;
$vars8 = Table[x[i], {i, $n8}];
(* FIX (a): sum starts at i=2 so x[1] appears only as a linear term,
   matching the phase5_8d.wl oracle integrand (no x[1]^2 in the sum). *)
$poly8 = 1 + $vars8[[1]] + Sum[$vars8[[i]]^2, {i, 2, $n8}];
$Bex8  = -6;

$spec8 = <|
  "Polynomials"         -> {$poly8},
  "MonomialExponents"   -> ConstantArray[0, $n8],
  "PolynomialExponents" -> {$Bex8},
  "Variables"           -> $vars8,
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> None
|>;

(* Oracle: Re = 0.00317086, confirmed by four independent oracles (Schwinger,
   Beta-function reduction, c=0 closed form, CUBA Vegas 1e9) in phase5_report. *)
$oracle8 = 0.00317086`;

(* FIX (b): Build the 8D fan via computeFanScaled (K-scaling; normal fan is
   scale-invariant) to handle thin lattice simplices in ambient dim>=4.
   Raw ComputeDecomposition can crash polymake at 8D; computeFanScaled retries
   with K in {1, n+2, 2n+4, 6n+6} (plan.md §6.4, tropical_fan.wl lines 475-484).
   computeFanScaled is a PUBLIC symbol of TropicalFan` (usage declared before
   Begin["`Private`"]; implementation in Private); call as TropicalFan`computeFanScaled. *)
$verts8 = Quiet[
  PolytopeVertices[$poly8^(-1), $vars8],
  {TropicalFan::polymake}];
Print["  Building 8D fan via computeFanScaled (K-scaling, §6.4)..."];
$fan8 = Quiet[
  TropicalFan`computeFanScaled[$verts8],
  {TropicalFan::polymake}];
If[!ListQ[$fan8] || Length[$fan8] < 2,
  Print["CC42 FAIL  8D fan build failed (computeFanScaled returned: ", $fan8, ")"];
  Quit[1]];
Print["  8D fan: ", Length[$fan8[[2]]], " cones"];

(* Verify resolved NStart for d=8: Ceiling[1000 * 3^(8-4)] = 81000. *)
$ns8   = TropicalEval`Private`resolveVegasSizing[Automatic, 8, 1000, "nstart"];
$ni8   = TropicalEval`Private`resolveVegasSizing[Automatic, 8,  500, "nincrease"];
$nb8   = TropicalEval`Private`resolveVegasSizing[Automatic, 8, 1000, "nbatch"];
$nsExp = Ceiling[1000 * 3^(8 - 4)];   (* = 81000 *)

Print["  8D resolved sizing: NStart=", $ns8, "  NIncrease=", $ni8,
      "  NBatch=", $nb8, "  (expected NStart=", $nsExp, ")"];

$sizingA = ($ns8 === $nsExp);
If[!$sizingA,
  Print["  WARNING: resolveVegasSizing gave unexpected NStart=", $ns8,
        " (expected ", $nsExp, "); continuing VEGAS run"]];

(* Run VEGAS with Automatic sizing.
   NSamples=4e7 clears the strengthened guard threshold (100*81000=8.1e6) so
   the adaptive grid resolves and the value is within 0.1% (phase5_report.txt:
   at 4e7 RELERR=0.0019%, vs 0.118% at 4e6 which is now flagged by vegasbudget).
   The check runs Quiet to suppress the budget warning if NSamples were too low,
   but at 4e7 the warning must NOT fire. *)
$wdA = FileNameJoin[{$ioDir, "A_8d_auto"}];
Quiet[CreateDirectory[$wdA], {CreateDirectory::eexist}];

Print["  Running 8D VEGAS NSamples=4e7 (clears 100*81000=8.1e6 threshold)..."];
$resA = EvaluateTropicalMC[$spec8, $fan8, {{}},
  "Integrator"       -> "VEGAS",
  "NSamples"         -> 40000000,
  "RunChecks"        -> False,
  "Verbose"          -> False,
  "WorkingDirectory" -> $wdA];

If[!AssociationQ[$resA] || !KeyExistsQ[$resA, "Results"],
  Print["CC42 FAIL  8D VEGAS run returned: ", $resA];
  Quit[1]];

$vARe  = $resA["Results"][[1]]["Re"];
$vAErr = $resA["Results"][[1]]["ReErr"];
$rA    = relErr[$vARe, $oracle8];

Print["  8D VEGAS Re = ", $vARe, " +/- ", $vAErr];
Print["  Oracle Re   = ", $oracle8, " (Schwinger/closed-form/CUBA-Vegas)"];
Print["  rel-err vs oracle: ", ScientificForm[$rA, 3],
      "  (gate < 1e-3)"];

record["A 8D Automatic sizing rel-err < 0.1%",
  $sizingA && $rA < 1*^-3,
  "expected=<1e-3 got=" <> ToString[ScientificForm[$rA, 3]] <>
    " NStart=" <> ToString[$ns8],
  "expected=<1e-3 got=" <> ToString[ScientificForm[$rA, 3]] <>
    " NStart=" <> ToString[$ns8] <> " (expected " <> ToString[$nsExp] <> ")"];
Print[];


(* ── 4.  Sub-check B: Low-D (n=2) sizing unchanged, VEGAS matches Gamma ─── *)
(*
   Integrand: P = 1 + x[1] + x[2],  monomial exponents (a1-1, a2-1),
   polynomial exponent -b = -5.  Exact: Gamma(a1)*Gamma(a2)*Gamma(b-a1-a2)/Gamma(b).
   resolveVegasSizing[Automatic, 2, base, kind] must return base unchanged
   (historical values: 1000, 500, 1000).  VEGAS value must match exact to <1e-4.
   Ported from TEST/phase5_lowd.wl.
*)

Print["=== Sub-check B: Low-D (n=2) sizing unchanged, VEGAS vs Gamma ==="];

$a1 = 5/4; $a2 = 3/2; $b = 5;
$exact2 = N[Gamma[$a1] * Gamma[$a2] * Gamma[$b - $a1 - $a2] / Gamma[$b], 16];

$poly2 = 1 + x[1] + x[2];
$vars2 = {x[1], x[2]};
$spec2 = <|
  "Polynomials"         -> {$poly2},
  "MonomialExponents"   -> {$a1 - 1, $a2 - 1},
  "PolynomialExponents" -> {-$b},
  "Variables"           -> $vars2,
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> None
|>;

(* Sizing must resolve to historical defaults for d=2. *)
$ns2 = TropicalEval`Private`resolveVegasSizing[Automatic, 2, 1000, "nstart"];
$ni2 = TropicalEval`Private`resolveVegasSizing[Automatic, 2,  500, "nincrease"];
$nb2 = TropicalEval`Private`resolveVegasSizing[Automatic, 2, 1000, "nbatch"];
$sizingB = ($ns2 === 1000 && $ni2 === 500 && $nb2 === 1000);

Print["  n=2 sizing: NStart=", $ns2, "  NIncrease=", $ni2,
      "  NBatch=", $nb2, "  historical-unchanged=", $sizingB];

$verts2 = PolytopeVertices[$poly2^(-1), $vars2];
$fan2   = ComputeDecomposition[$verts2, "ShowProgress" -> False];

$wdB = FileNameJoin[{$ioDir, "B_lowd"}];
Quiet[CreateDirectory[$wdB], {CreateDirectory::eexist}];

(* Intercept vegasbudget messages -- should NOT fire at NSamples=5e6 for d=2
   (threshold = 20*1000 = 20000; 5e6 >> 20000). *)
$budgetFiredB = False;
Internal`HandlerBlock[
  {"Message",
   Function[m,
     If[MatchQ[m, Hold[Message[MessageName[TropicalEval, "vegasbudget"], ___], _]] ||
        MatchQ[m, Hold[Message[TropicalEval`EvaluateTropicalMC::vegasbudget, ___], _]] ||
        StringContainsQ[ToString[m], "vegasbudget"],
        $budgetFiredB = True]]},
  $resB = EvaluateTropicalMC[$spec2, $fan2, {{}},
    "Integrator"       -> "VEGAS",
    "NSamples"         -> 5000000,
    "RunChecks"        -> False,
    "Verbose"          -> False,
    "WorkingDirectory" -> $wdB]];

$vBRe = $resB["Results"][[1]]["Re"];
$rB   = relErr[$vBRe, $exact2];

Print["  n=2 VEGAS Re = ", $vBRe, "  exact = ", $exact2];
Print["  rel-err vs Gamma: ", ScientificForm[$rB, 3],
      "  (gate < 1e-4)"];
Print["  vegasbudget fired (should be False): ", $budgetFiredB];

record["B Low-D n=2 sizing unchanged and VEGAS vs Gamma < 1e-4",
  $sizingB && $rB < 1*^-4 && !$budgetFiredB,
  "sizing-unchanged=" <> ToString[$sizingB] <>
    " relerr=" <> ToString[ScientificForm[$rB, 3]] <>
    " budget-silent=" <> ToString[!$budgetFiredB],
  "expected: sizing-unchanged=True relerr<1e-4 budget-silent=True  got: " <>
    "sizing=" <> ToString[$sizingB] <>
    " relerr=" <> ToString[ScientificForm[$rB, 3]] <>
    " budgetFired=" <> ToString[$budgetFiredB]];
Print[];


(* ── 5.  Sub-check C: vegasbudget fires at NSamples < guardFactor*NStart ─── *)
(*
   The engine fires vegasbudget when:
     guardFactor = 100  for maxSectorDim >= 5  (high-D)
     guardFactor = 20   for maxSectorDim <= 4  (low-D, unchanged)
   (tropical_eval.wl vegasbudget guard, phase5_report.txt §ISSUE 1(c))

   Two scenarios:
     C1  8D integrand, NSamples=5000000.
         Resolved NStart for d=8 is 81000; guardFactor=100 (d>=5);
         threshold = 100*81000 = 8100000.
         5000000 < 8100000 => must fire.
         This was SILENT under the old 20x rule (1620000 < 5000000);
         the strengthened 100x rule catches it.  Non-vacuous test.

     C2  2D integrand, NSamples=1000, explicit NStart=1000.
         guardFactor=20 (d<=4); threshold=20*1000=20000.
         1000 < 20000 => must fire.

   We intercept TropicalEval::vegasbudget via HandlerBlock and also check
   the returned association for a "VegasBudgetWarning" key as a fallback.
*)

Print["=== Sub-check C: vegasbudget fires at NSamples < guardFactor*NStart ==="];

interceptBudget[body_] :=
  Module[{fired = False},
    Internal`HandlerBlock[
      {"Message",
       Function[m,
         If[MatchQ[m, Hold[Message[MessageName[TropicalEval, "vegasbudget"], ___], _]] ||
            MatchQ[m, Hold[Message[TropicalEval`EvaluateTropicalMC::vegasbudget, ___], _]] ||
            StringContainsQ[ToString[m], "vegasbudget"],
            fired = True]]},
      Quiet[body]];
    fired];

(* C1: 8D, NSamples=5000000 < 100*81000=8100000 (above old 20x=1620000). *)
$wdC1 = FileNameJoin[{$ioDir, "C1_budget_8d"}];
Quiet[CreateDirectory[$wdC1], {CreateDirectory::eexist}];

Print["  C1: 8D NSamples=5000000 (below 100*81000=8100000, above old 20*81000=1620000)"];
$firedC1 = interceptBudget[
  EvaluateTropicalMC[$spec8, $fan8, {{}},
    "Integrator"       -> "VEGAS",
    "NSamples"         -> 5000000,
    "RunChecks"        -> False,
    "Verbose"          -> False,
    "WorkingDirectory" -> $wdC1]];

(* Fallback: probe via return-value key. *)
If[!$firedC1,
  Module[{r = Quiet[EvaluateTropicalMC[$spec8, $fan8, {{}},
      "Integrator"       -> "VEGAS",
      "NSamples"         -> 5000000,
      "RunChecks"        -> False,
      "Verbose"          -> False,
      "WorkingDirectory" -> $wdC1]]},
    If[AssociationQ[r] &&
       (TrueQ[r["VegasBudgetWarning"]] ||
        (KeyExistsQ[r, "Warnings"] &&
         AnyTrue[r["Warnings"], StringContainsQ[#, "budget"] &])),
       $firedC1 = True]]];

Print["  C1  8D NSamples=5000000: vegasbudget fired = ", $firedC1,
      "  (threshold=100*81000=8100000 > 5000000, expected True)"];

(* C2: 2D, NSamples=1000, explicit NStart=1000 (guardFactor=20, threshold=20000). *)
$wdC2 = FileNameJoin[{$ioDir, "C2_budget_2d"}];
Quiet[CreateDirectory[$wdC2], {CreateDirectory::eexist}];

$firedC2 = interceptBudget[
  EvaluateTropicalMC[$spec2, $fan2, {{}},
    "Integrator"       -> "VEGAS",
    "NSamples"         -> 1000,
    "VegasNStart"      -> 1000,
    "RunChecks"        -> False,
    "Verbose"          -> False,
    "WorkingDirectory" -> $wdC2]];

If[!$firedC2,
  Module[{r = Quiet[EvaluateTropicalMC[$spec2, $fan2, {{}},
      "Integrator"       -> "VEGAS",
      "NSamples"         -> 1000,
      "VegasNStart"      -> 1000,
      "RunChecks"        -> False,
      "Verbose"          -> False,
      "WorkingDirectory" -> $wdC2]]},
    If[AssociationQ[r] &&
       (TrueQ[r["VegasBudgetWarning"]] ||
        (KeyExistsQ[r, "Warnings"] &&
         AnyTrue[r["Warnings"], StringContainsQ[#, "budget"] &])),
       $firedC2 = True]]];

Print["  C2  2D NSamples=1000 NStart=1000: vegasbudget fired = ", $firedC2,
      "  (threshold=20*1000=20000 > 1000, expected True)"];

record["C vegasbudget fires at NSamples < guardFactor*NStart (8D and 2D)",
  $firedC1 && $firedC2,
  "C1(8D,100x)=" <> ToString[$firedC1] <> " C2(2D,20x)=" <> ToString[$firedC2],
  "expected=True True  got=C1=" <> ToString[$firedC1] <>
    " C2=" <> ToString[$firedC2]];
Print[];


(* ── 6.  Final summary ───────────────────────────────────────────────────── *)

$nPass = Length[$passes];
$nFail = Length[$failures];
$nTot  = $nPass + $nFail;

Print["================================================================"];
If[$nFail == 0,
  Print["CC42 PASS  all ", $nTot, " sub-checks passed",
        "  (8D Automatic sizing <0.1%; n=2 sizing unchanged; vegasbudget fires)"],
  Print["CC42 FAIL  ", $nFail, " of ", $nTot, " sub-checks failed:  ",
        StringRiffle[$failures, ", "]]];
Print["================================================================"];

Quit[If[$nFail > 0, 1, 0]]
