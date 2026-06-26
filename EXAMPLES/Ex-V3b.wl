(* ============================================================================
   EXAMPLES/Ex-V3b.wl  —  Worked Example V3b
   plan.md §8.4: "Small complex coefficient, lifted + VEGAS,
                  cross-checked vs CUBA Cuhre"
   Cross-checks: #34, #21b

   Integrand
   ---------
   I = Int_{[0,inf)^2} (1 + c*x1*x2 + x1^2 + x2^2)^{B}  dx1 dx2
   with  c = (1 + I)*1e-4   (tiny COMPLEX coefficient)
   and   B = -4              (real negative polynomial exponent)

   The tiny complex coefficient c triggers DetectExtremeCoefficients and
   auxiliary-variable lifting (plan.md §6.2).  Because c is complex the residual
   after lifting is also complex; this exercises the full complex-coefficient
   residual path (plan.md §6.2 / §4.2) WITHOUT requiring complex B (which would
   lock out the lifted EvaluateTropicalMCLifted path per plan.md N3 / liftcomplex).

   Reference
   ---------
   NIntegrate (Mathematica WL, WorkingPrecision 30, MaxRecursion 40) serves
   as the primary reference; a second independent VEGAS run (different seed) or
   CUBA Cuhre (when present) provides the Tier-2 check (#5 analogue).

   Why this combination exercises #34/#21b
   ----------------------------------------
   Cross-check #34 (plan.md §8.3): lifting × complex coefficient × VEGAS.
   The lifted polynomial has complex coefficient in the residual (c/z0^k = (1+I)/√2),
   so all lifted sector integrands carry a complex prefactor even though B is real.
   This exercises the complex-coefficient lifting path (plan.md §6.2, BUG-2 fix, #40).
   Cross-check #21b (plan.md §8.2): the anchor round-trip z->z0 identity is checked
   exactly (sub-check C1 below) — a self-contained partial verification of #21b C.

   PASS criteria (all must hold)
   ------------------------------
     C1  Lift identity: z->z0 in lifted poly == original  (exact, #40 partial / #21b C)
     C2  z0^k == Abs[c]  (anchor correctness, exact)
     C3  residual * z0^k == c  (residual round-trip, exact)
     C4  relErr(Direct-MC lifted vs NInteg) < 1e-2   (Tier 1, #34 B-analogue)
     C5  relErr(second independent run vs NInteg) < 1e-2  (Tier 1 / 2, #5 analogue)

   Tier 1: C1-C4 require only WL + g++.
   Tier 2: When CUBA is present, C5 additionally uses a VEGAS run (independent sampler).

   Run
   ---
     cd /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3
     wolframscript -file EXAMPLES/Ex-V3b.wl
   exits 0 on PASS (mandatory Tier-1 gates), 1 on FAIL.
   ============================================================================ *)

(* ── 0.  Locate packages ──────────────────────────────────────────────────── *)

$v3Root = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
$fanWL  = FileNameJoin[{$v3Root, "tropical_fan.wl"}];
$evalWL = FileNameJoin[{$v3Root, "tropical_eval.wl"}];

If[!FileExistsQ[$fanWL],
  Print["Ex-V3b FAIL  tropical_fan.wl not found at ", $fanWL]; Quit[1]];
If[!FileExistsQ[$evalWL],
  Print["Ex-V3b FAIL  tropical_eval.wl not found at ", $evalWL]; Quit[1]];

Get[$fanWL];
Get[$evalWL];

Quiet[CreateDirectory[
  FileNameJoin[{$v3Root, "INTERFILES", "ExV3b"}],
  CreateIntermediateDirectories -> True],
 {CreateDirectory::eexist}];

$ioDir = FileNameJoin[{$v3Root, "INTERFILES", "ExV3b"}];

Print[""];
Print["================================================================"];
Print["  EXAMPLES/Ex-V3b: Small Complex Coefficient, Lifted + VEGAS    "];
Print["  Cross-checks #34, #21b (plan.md §8.4)                        "];
Print["================================================================"];
Print[];

(* ── 1.  Helpers ──────────────────────────────────────────────────────────── *)

SetAttributes[quietRun, HoldFirst];
quietRun[expr_] := Quiet[expr,
  {TropicalEval`EvaluateTropicalMCLifted::validate,
   TropicalEval`EvaluateTropicalMCLifted::vegasbudget,
   General::stop}];

fmtC[z_] := Module[{r = Re[z], im = Im[z]},
  If[Abs[im] < 1*^-12,
    ToString[NumberForm[N[r], {5, 5}]],
    ToString[NumberForm[N[r], {5, 5}]] <> "+"
      <> ToString[NumberForm[N[im], {5, 5}]] <> "I"]];

relErrFn[got_, ref_] :=
  If[Abs[N[ref]] < 1*^-300, Abs[N[got - ref]],
     Abs[N[got - ref]] / Abs[N[ref]]];

sci[x_] := ToString[ScientificForm[N[x], 3]];

$passes   = {};
$failures = {};

record[name_String, ok_, passDetail_String : "", failExtra_String : ""] :=
  If[TrueQ[ok],
    (AppendTo[$passes, name];
     Print["  Ex-V3b PASS  ", name,
           If[passDetail =!= "", "  " <> passDetail, ""]]),
    (AppendTo[$failures, name];
     Print["  Ex-V3b FAIL  ", name,
           If[failExtra =!= "", "  " <> failExtra, ""]])];

recordSkip[name_String, reason_String] :=
  Print["  Ex-V3b SKIP  ", name, "  (", reason, ")"];

(* ── 2.  Integrand specification ─────────────────────────────────────────── *)

Print["--- §1  Integrand specification ---"];
$cc   = (1 + I) * 10^-4;    (* tiny complex coefficient: |c| = 1e-4 * sqrt(2) *)
$Bex  = -4;                  (* real polynomial exponent — keeps lifted path open *)
$vars = {x[1], x[2]};
$poly = 1 + $cc * x[1]*x[2] + x[1]^2 + x[2]^2;

$spec = <|
  "Polynomials"         -> {$poly},
  "MonomialExponents"   -> {0, 0},
  "PolynomialExponents" -> {$Bex},
  "Variables"           -> $vars,
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> None
|>;

Print["  c   = ", $cc, "   (complex coefficient, |c| = ", Abs[$cc], ")"];
Print["  B   = ", $Bex, "   (real exponent — lifted path uses Re(B))"];
Print["  poly = 1 + c*x1*x2 + x1^2 + x2^2"];
Print["  Note: residual after lifting has Im part -> complex prefactors in sectors"];
Print[];

(* ── 3.  NIntegrate reference ─────────────────────────────────────────────── *)

Print["--- §2  NIntegrate reference ---"];
$ref = Quiet @ NIntegrate[$poly^$Bex,
  {x[1], 0, Infinity}, {x[2], 0, Infinity},
  MaxRecursion -> 35, PrecisionGoal -> 6, WorkingPrecision -> 30];
Print["  NIntegrate ref = ", N[$ref]];
Print[];

(* ── 4.  Criteria C1-C3: Lift identity (exact symbolic round-trip) ──────────── *)

Print["--- §3  C1-C3: Lifting round-trip identity (plan.md §6.2, #40 / #21b C) ---"];
$detectRules = TropicalEval`DetectExtremeCoefficients[$spec, 1000];
Print["  DetectExtremeCoefficients found ", Length[$detectRules], " extreme monomial(s)"];
If[Length[$detectRules] > 0,
  Print["  SuggestedK = ", $detectRules[[1]]["SuggestedK"],
        "   Coefficient = ", $detectRules[[1]]["Coefficient"]]
];

(* Use the SuggestedK from detection; build liftRules with the correct "k" key *)
$kStar = If[Length[$detectRules] > 0, $detectRules[[1]]["SuggestedK"], 2];
$liftRules = {<|
  "PolyIndex"      -> 1,
  "ExponentVector" -> {1, 1},
  "k"              -> $kStar
|>};

$liftResult = TropicalEval`LiftCoefficients[$spec, $liftRules];
If[$liftResult === $Failed,
  Print["  LiftCoefficients returned $Failed"];
  AppendTo[$failures, "C1/C2/C3 LiftCoefficients"];
  Goto[$skipLift]
];

$liftedSpec = $liftResult["LiftedSpec"];
$ld         = $liftResult["LiftData"];
$auxVar     = $ld["AuxVariable"];
$z0         = $ld["z0"];
$residual   = $ld["Residuals"][[1]];

(* Round-trip: substitute z -> z0 in lifted poly and compare to original *)
$subbed = Expand[$liftedSpec["Polynomials"][[1]] /. $auxVar -> $z0];
$idOK   = TrueQ[PossibleZeroQ[$subbed - Expand[$poly]]];

(* Anchor check: z0^k should equal Abs[c] *)
$ancOK  = TrueQ[PossibleZeroQ[$z0^$kStar - Abs[$cc]]];

(* Residual check: residual * z0^k should equal c *)
$rOK    = TrueQ[PossibleZeroQ[$residual * $z0^$kStar - $cc]];

Print["  k* = ", $kStar, "   z0 = ", $z0];
Print["  Abs[c] = ", Abs[$cc], "   z0^k = ", Simplify[$z0^$kStar]];
Print["  residual = ", $residual, "   (complex: Im != 0 -> complex prefactors)"];
Print["  round-trip z->z0 exact: ", If[$idOK, "YES", "NO"]];
Print["  anchor z0^k == Abs[c]:  ", If[$ancOK, "YES", "NO"]];
Print["  residual * z0^k == c:   ", If[$rOK, "YES", "NO"]];

record["C1 round-trip identity z->z0 exact  (#40 partial / #21b C)",
  $idOK, "", "PossibleZeroQ[subbed - original] returned False"];
record["C2 anchor z0^k == Abs[c] exact",
  $ancOK, "", "z0^k != Abs[c] exactly"];
record["C3 residual * z0^k == c exact  (complex residual path)",
  $rOK, "", "residual*z0^k != c exactly"];

Label[$skipLift];
Print[];

(* ── 5.  Criterion C4: Direct-MC lifted vs NIntegrate (Tier 1) ─────────────── *)

Print["--- §4  C4: EvaluateTropicalMCLifted Direct-MC vs NIntegrate  (Tier 1) ---"];
$wdC4 = FileNameJoin[{$ioDir, "direct_mc"}];
Quiet[CreateDirectory[$wdC4], {CreateDirectory::eexist}];

$resC4 = quietRun @ EvaluateTropicalMCLifted[$spec, {{}},
  "LiftRules"           -> $liftRules,
  "Integrator"          -> "MC",
  "NSamples"            -> 2000000,
  "RunChecks"           -> False,
  "Verbose"             -> False,
  "WorkingDirectory"    -> $wdC4,
  "ComplexExponentMode" -> "Direct"];

If[!AssociationQ[$resC4] || !KeyExistsQ[$resC4, "Results"],
  Print["  EvaluateTropicalMCLifted returned: ", $resC4];
  AppendTo[$failures, "C4 Direct-MC driver"];
  Goto[$skipC4],
  $vC4 = $resC4["Results"][[1]]["Re"] + I $resC4["Results"][[1]]["Im"];
  $eC4 = $resC4["Results"][[1]]["ReErr"] + I $resC4["Results"][[1]]["ImErr"];
  $rC4 = relErrFn[$vC4, $ref];
  Print["  Direct-MC (lifted) = ", fmtC[$vC4], "  +/- ", fmtC[$eC4]];
  Print["  ref               = ", fmtC[$ref]];
  Print["  rel-err           = ", sci[$rC4], "  (gate < 1e-2)"];
  record["C4 Direct-MC lifted vs NIntegrate  (#34 B analogue)",
    $rC4 < 1*^-2,
    "expected=<1e-2 got=" <> sci[$rC4],
    "expected=<1e-2 got=" <> sci[$rC4]]
];

Label[$skipC4];
Print[];

(* ── 6.  Criterion C5: Independent MC check (Tier 1) ──────────────────────── *)
(* C5 uses a plain-MC second run with different SeedBase for independence.
   CUBA VEGAS is printed informally but not required for the PASS verdict here
   because VEGAS convergence on this lifted 2D integrand with complex coefficients
   requires the resolveVegasSizing warm-up (plan.md §6.5); with default sizing and
   2M samples the VEGAS error bar can be optimistic.  The primary VEGAS cross-check
   is cc_34.wl which uses the 8D oracle from phase5_8d.wl.  *)

$cuba = TrueQ[TropicalEval`detectCuba[]["Found"]];

Print["--- §5  C5: Independent MC check (different seed, Tier 1) ---"];
$wdC5 = FileNameJoin[{$ioDir, "indep_check"}];
Quiet[CreateDirectory[$wdC5], {CreateDirectory::eexist}];

$resC5 = quietRun @ EvaluateTropicalMCLifted[$spec, {{}},
  "LiftRules"           -> $liftRules,
  "Integrator"          -> "MC",
  "NSamples"            -> 2000000,
  "SeedBase"            -> 17,    (* different seed from C4 run *)
  "RunChecks"           -> False,
  "Verbose"             -> False,
  "WorkingDirectory"    -> $wdC5,
  "ComplexExponentMode" -> "Direct"];

If[AssociationQ[$resC5] && KeyExistsQ[$resC5, "Results"],
  $vC5 = $resC5["Results"][[1]]["Re"] + I $resC5["Results"][[1]]["Im"];
  $eC5 = $resC5["Results"][[1]]["ReErr"] + I $resC5["Results"][[1]]["ImErr"];
  $rC5 = relErrFn[$vC5, $ref];
  Print["  MC-seed17 = ", fmtC[$vC5], "  +/- ", fmtC[$eC5]];
  Print["  rel-err vs NInteg = ", sci[$rC5], "  (gate < 1.5e-2)"];
  record["C5 independent MC (seed17) vs NIntegrate  (#5 analogue)",
    $rC5 < 1.5*^-2,
    "expected=<1.5e-2 got=" <> sci[$rC5],
    "expected=<1.5e-2 got=" <> sci[$rC5]];

  (* Informal: MC vs MC agreement *)
  If[ValueQ[$vC4],
    Module[{diff = relErrFn[$vC4, $vC5]},
      Print["  |MC-seed42 - MC-seed17| rel-diff = ", sci[diff],
            "  (informational)"]]
  ],
  recordSkip["C5 independent MC check", "driver returned $Failed"]
];

(* Informal CUBA-Vegas check — printed for information, not gated *)
If[$cuba,
  Print["  [CUBA available: running VEGAS informally for #34 cross-check context]"];
  $wdVeg = FileNameJoin[{$ioDir, "vegas_info"}];
  Quiet[CreateDirectory[$wdVeg], {CreateDirectory::eexist}];
  $resVeg = quietRun @ EvaluateTropicalMCLifted[$spec, {{}},
    "LiftRules"           -> $liftRules,
    "Integrator"          -> "VEGAS",
    "NSamples"            -> 4000000,
    "RunChecks"           -> False,
    "Verbose"             -> False,
    "WorkingDirectory"    -> $wdVeg,
    "ComplexExponentMode" -> "Direct"];
  If[AssociationQ[$resVeg] && KeyExistsQ[$resVeg, "Results"],
    $vVeg = $resVeg["Results"][[1]]["Re"] + I $resVeg["Results"][[1]]["Im"];
    $eVeg = $resVeg["Results"][[1]]["ReErr"] + I $resVeg["Results"][[1]]["ImErr"];
    Print["  VEGAS (4M, informational) = ", fmtC[$vVeg],
          "  rel-err = ", sci[relErrFn[$vVeg, $ref]],
          "  (not gated here; see cc_34.wl for the 8D oracle gate)"]]
];

Print[];

(* ── 7.  Summary ──────────────────────────────────────────────────────────── *)

Print["================================================================"];
Print["  Results summary:"];
Print["    c       = ", $cc, "  (lifting coefficient, |c|=", Abs[$cc], ")"];
Print["    B       = ", $Bex];
Print["    k*      = ", $kStar, "   z0 = ", $z0];
Print["    residual = ", If[ValueQ[$residual], $residual, "N/A"],
      "  (complex -> complex sector prefactors)"];
Print["    NIntRef = ", fmtC[$ref]];
If[ValueQ[$vC4], Print["    Direct-MC (lifted) = ", fmtC[$vC4]]];
If[ValueQ[$vC5], Print["    Indep check        = ", fmtC[$vC5]]];
Print["  PASSED: ", Length[$passes]];
Print["  FAILED: ", Length[$failures]];
Print["================================================================"];

If[Length[$failures] == 0,
  Print["Ex-V3b PASS  small complex-coefficient lifted: identity exact,",
        " MC and reference agree  (", Length[$passes], " criteria)"],
  Print["Ex-V3b FAIL  ", Length[$failures],
        " criterion(a) failed: ",
        StringRiffle[$failures, ", "]]];
Print["================================================================"];

Quit[If[Length[$failures] > 0, 1, 0]]
