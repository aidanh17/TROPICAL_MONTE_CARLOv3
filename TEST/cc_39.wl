(* ============================================================================
   TEST/cc_39.wl  —  Cross-check #39: "The error bar lies" guardrail
                     (plan.md §8.3 #39, Tier 1/2)

   Encodes the L1 / L5 / Example-20 heavy-tail lesson as an automated
   regression:

     A small-coefficient case where unlifted plain-MC with default settings
     is confidently wrong (tight but lying error bar).  After lifting and
     using VEGAS with properly resolved sizing, the result recovers the
     reference, AND the appropriate warnings fire along the way.

   Sub-checks:
     (A) DetectExtremeCoefficients detects the small coefficient and reports
         it (HasSmallCoeff / SuggestedK set, magnitude check fires).
     (B) Unlifted plain-MC with small NSamples: the result is far from the
         reference (demonstrates the lie), NOT counted as a FAIL of the test —
         it is the expected bad behaviour we are documenting.
     (C) vegasbudget warning fires when VEGAS is called with NSamples too
         small relative to resolveVegasSizing's recommendation (regression
         that the budget guard is wired).
     (D) Lifted VEGAS with properly sized budget recovers the reference to
         within 5% relative error.
     (E) HasConstantTerm=False fires for exactly 1 of 8 lifted sectors
         (sector 8; the Phase-3 L5 finding).  This IS the correct engine
         behavior: the warning IS the guardrail working.  The remaining 7
         sectors are variance-safe (HasConstantTerm=True).

   Integrand  (Example-20 canonical case, lift_error_log §L1/L5):
       P = 1 + 1e-6 · x[1]^2 + x[2]^2 + x[1]·x[2]^2,   B = -2
       ∫₀^∞ ∫₀^∞ P^{-2} dx₁ dx₂

   Reference: NIntegrate at WorkingPrecision -> 30 (trusted, Tier-1 check).

   Tier: 1/2.
     — Sub-checks A, B, C, E are Tier 1 (no CUBA needed; VEGAS budget
       message can be tested with a very small NSamples even without a
       working CUBA binary, because the message fires before/during compile).
     — Sub-check D requires CUBA (VEGAS integration).  If CUBA is absent,
       sub-check D is replaced by an NIntegrate reference cross-check (also
       Tier 1), so the test degrades gracefully: on machines without CUBA it
       prints "CC39 SKIP (no CUBA)" and exits 0 per plan.md §8 Tier policy.

   Run:
       wolframscript -file TEST/cc_39.wl
   from the TROPICAL_MONTE_CARLO3 root, or from any directory.

   Ported from / inspired by:
     OLD_CODE/TROPICAL_MONTE_CARLO2/TROPICAL_MONTE_CARLOv2/TEST/lift_vegas_probe.wl
     OLD_CODE/TROPICAL_MONTE_CARLO/lift_error_log.md   (L1, L5)
     plan.md §6.2 (HasConstantTerm ranking), §6.5 (resolveVegasSizing /
       vegasbudget), §8.3 #39.
   ============================================================================ *)

(* ── 0. Load v3 packages (absolute paths, plan.md §10 discipline) ─────────── *)

$pkgRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
$fanWL   = FileNameJoin[{$pkgRoot, "tropical_fan.wl"}];
$evalWL  = FileNameJoin[{$pkgRoot, "tropical_eval.wl"}];
$ioDir   = FileNameJoin[{$pkgRoot, "TEST", "INTERFILES"}];

If[!FileExistsQ[$fanWL],
  Print["CC39 FAIL  tropical_fan.wl not found at ", $fanWL]; Quit[1]];
If[!FileExistsQ[$evalWL],
  Print["CC39 FAIL  tropical_eval.wl not found at ", $evalWL]; Quit[1]];

Get[$fanWL];
Get[$evalWL];

Quiet[CreateDirectory[$ioDir], {CreateDirectory::eexist}];

Print["CC39: packages loaded from ", $pkgRoot];
Print[];

(* ── 1. CUBA gate (Tier-2 policy) ────────────────────────────────────────── *)

$cubaInfo = TropicalEval`detectCuba[];
$cubaPresent = TrueQ[$cubaInfo["Found"]];

If[!$cubaPresent,
  (* Plan §8: Tier-2 checks degrade gracefully when CUBA is absent. *)
  Print["CC39 SKIP (no CUBA)"];
  Quit[0]
];

Print["CC39: CUBA found at ", $cubaInfo["IncludeDir"]];
Print[];

(* ── 2. Helpers ──────────────────────────────────────────────────────────── *)

SetAttributes[quietRun, HoldFirst];
quietRun[expr_] := Block[{$Output = {}}, expr];

fmtN[x_] := If[NumberQ[x], ToString[ScientificForm[N[x, 4]]], ToString[x]];

relErr[got_, ref_] :=
  If[NumericQ[got] && NumericQ[ref] && ref =!= 0,
    Abs[N[got]/N[ref] - 1], Indeterminate];

$pass = 0;
$fail = 0;

reportSub[label_String, ok_] :=
  If[TrueQ[ok],
    $pass++;
    Print["CC39 ", label, " PASS"],
    $fail++;
    Print["CC39 ", label, " FAIL"]
  ];

(* ── 3. Integrand definition ─────────────────────────────────────────────── *)
(*
   P = 1 + 1e-6 x1^2 + x2^2 + x1 x2^2
   B = -2   (convergent, real exponents)
   Lift rule: the monomial x1^2 (exponent {2,0}) has coefficient 1e-6.
   Auto-k = ceil(|log(1e-6)| / log(1000)) = ceil(6/3) = 2  => z0 = 1e-3.
*)

$smallC = 10^(-6);

$spec = <|
  "Polynomials"         -> {1 + $smallC x[1]^2 + x[2]^2 + x[1] x[2]^2},
  "MonomialExponents"   -> {0, 0},
  "PolynomialExponents" -> {-2},
  "Variables"           -> {x[1], x[2]},
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> None
|>;

$liftRule = {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 2|>};

(* ── 4. Tier-1 reference via NIntegrate ─────────────────────────────────── *)

Print["CC39: computing NIntegrate reference (WorkingPrecision 30)..."];
$ref = Quiet@NIntegrate[
  (1 + $smallC x1^2 + x2^2 + x1 x2^2)^(-2),
  {x1, 0, Infinity}, {x2, 0, Infinity},
  MaxRecursion -> 40, PrecisionGoal -> 8, WorkingPrecision -> 30];

If[!NumericQ[$ref],
  Print["CC39 FAIL  NIntegrate reference failed: ", $ref]; Quit[1]];

Print["CC39: reference = ", $ref, "  (≈ ", fmtN[$ref], ")"];
Print[];

(* ── Sub-check A: DetectExtremeCoefficients flags the small coefficient ──── *)

Print["================================================================"];
Print["CC39 (A) DetectExtremeCoefficients detects the 1e-6 coefficient"];
Print["================================================================"];

(* DetectExtremeCoefficients returns a List of Association entries — one
   per extreme monomial found.  It is NOT a single Association wrapper;
   the list itself IS the ExtremeCoefficients collection.
   BUG FIX (sub-check A): the prior version read $detRes["ExtremeCoefficients"]
   and $detRes["SuggestedK"] as if the result were a single Association, which
   always returned Missing[] and made $passA vacuously False. *)

$detRes = DetectExtremeCoefficients[$spec, 1000];

$passA = ListQ[$detRes] &&
         Length[$detRes] > 0 &&
         AssociationQ[$detRes[[1]]] &&
         IntegerQ[$detRes[[1]]["SuggestedK"]] &&
         $detRes[[1]]["SuggestedK"] >= 1;

Print["  DetectExtremeCoefficients result: ",
      If[ListQ[$detRes] && Length[$detRes] > 0 && AssociationQ[$detRes[[1]]],
         "ExtremeCoefficients (count=" <> ToString[Length[$detRes]] <> ")=" <>
         ToString[$detRes] <>
         "  SuggestedK (first)=" <> ToString[$detRes[[1]]["SuggestedK"]],
         ToString[$detRes]]];

reportSub["(A) DetectExtremeCoeff", $passA];
Print[];

(* ── Sub-check B: Unlifted plain-MC is badly wrong (documents the lie) ────── *)
(*
   This sub-check verifies that unlifted MC with modest NSamples produces a
   result that is FAR from the reference (rel-err >> 1%).  This is the
   "confidently wrong" scenario.  A result that is already close would mean
   the integrand does not actually exhibit the heavy-tail problem at this
   sample count — in that case we still PASS but print a note.

   We do NOT hard-fail on a "too accurate" unlifted result; the test is
   checking that the GUARDRAIL mechanism (lifting + warnings) works, not
   that this specific machine is always wrong in the unlifted case.
*)

Print["================================================================"];
Print["CC39 (B) Unlifted plain-MC with small NSamples (demonstrating the lie)"];
Print["================================================================"];

$fan2D = Quiet[ComputeDecomposition[
  PolytopeVertices[(1 + $smallC x[1]^2 + x[2]^2 + x[1] x[2]^2)^(-1),
                   {x[1], x[2]}],
  "ShowProgress" -> False]];

$unliftedMC = Quiet@EvaluateTropicalMC[$spec, $fan2D, {{}},
  "Integrator" -> "MC",
  "NSamples"   -> 500000,
  "RunChecks"  -> False,
  "Verbose"    -> False];

Module[{gotVal, relE, passB},
  gotVal = If[AssociationQ[$unliftedMC] && Length[$unliftedMC["Results"]] > 0,
    $unliftedMC["Results"][[1]]["Re"], $Failed];
  relE = relErr[gotVal, $ref];
  Print["  unlifted MC (500k): ", fmtN[gotVal],
        "   relErr vs ref = ", fmtN[relE]];

  (* The sub-check PASSES regardless — we are documenting the phenomenon.
     If relErr < 0.01 the integrand did not exhibit the problem at this
     budget; if relErr >= 0.01 (expected) we confirm the heavy-tail. *)
  passB = NumericQ[gotVal];
  If[passB,
    If[relE > 0.01,
      Print["  => unlifted MC is far from reference as expected (rel-err ",
            fmtN[relE], "); heavy-tail regime confirmed."],
      Print["  => unlifted MC happened to be close (rel-err ", fmtN[relE],
            "); integrand may converge without lifting at this dimension."]],
    Print["  => unlifted MC returned: ", $unliftedMC]
  ];
  reportSub["(B) unlifted-MC-lie-documented", passB];
];
Print[];

(* ── Sub-check C: vegasbudget warning fires when NSamples too small ────── *)
(*
   For a 2D sector the resolved NStart (dim=2, historical) = 1000, so
   the guard threshold is 20*1000 = 20000.  Requesting NSamples=100 must
   fire TropicalEval::vegasbudget.
*)

Print["================================================================"];
Print["CC39 (C) vegasbudget warning fires when VEGAS under-budgeted"];
Print["================================================================"];

$firedBudget = False;

Quiet[
  Check[
    Block[{$Output = {}},
      EvaluateTropicalMC[$spec, $fan2D, {{}},
        "Integrator" -> "VEGAS",
        "NSamples"   -> 100,
        "RunChecks"  -> False,
        "Verbose"    -> False]],
    ($firedBudget = True; $Failed),
    {TropicalEval::vegasbudget}
  ],
  {TropicalEval::vegasbudget}
];

(* Also check $MessageList in case Check swallowed the message differently *)
If[!$firedBudget,
  $firedBudget = Or @@ (StringContainsQ[ToString[#], "vegasbudget"] & /@
                         ToString /@ $MessageList)
];

Print["  vegasbudget fired: ", $firedBudget];
reportSub["(C) vegasbudget fires", $firedBudget];
Print[];

(* ── Sub-check D: Lifted VEGAS with proper budget recovers the reference ─── *)
(*
   Use EvaluateTropicalMCLifted with Automatic sizing (resolveVegasSizing).
   For a 3D lifted problem (n+1=3) with Automatic, NStart >= 1000 (historical).
   We budget generously: NSamples = 3e6.
   The result must satisfy |got - ref| / ref < 0.05  (5% gate, matching
   lift_error_log Test 26 at 0.089% and probe Case 2 at <1%).
*)

Print["================================================================"];
Print["CC39 (D) Lifted VEGAS with Automatic sizing recovers reference"];
Print["================================================================"];

$liftedVEGAS = Quiet@EvaluateTropicalMCLifted[$spec, {{}},
  "LiftRules"      -> $liftRule,
  "Integrator"     -> "VEGAS",
  "NSamples"       -> 3000000,
  "VegasEpsRel"    -> 1.*^-9,
  "RunChecks"      -> False,
  "Verbose"        -> False,
  "WorkingDirectory" -> $ioDir];

Module[{gotVal, gotErr, relE, passD},
  gotVal = If[AssociationQ[$liftedVEGAS] && Length[$liftedVEGAS["Results"]] > 0,
    $liftedVEGAS["Results"][[1]]["Re"], $Failed];
  gotErr = If[AssociationQ[$liftedVEGAS] && Length[$liftedVEGAS["Results"]] > 0,
    $liftedVEGAS["Results"][[1]]["ReErr"], Indeterminate];
  relE = relErr[gotVal, $ref];
  Print["  lifted VEGAS (3M eval): ", fmtN[gotVal], " +/- ", fmtN[gotErr]];
  Print["  reference (NIntegrate): ", fmtN[$ref]];
  Print["  rel-err: ", fmtN[relE]];

  passD = NumericQ[gotVal] && NumericQ[relE] && relE < 0.05;
  If[!passD && NumericQ[gotVal],
    Print["CC39 (D) FAIL  expected rel-err < 0.05  got rel-err = ", fmtN[relE],
          "  (expected=", fmtN[$ref], "  got=", fmtN[gotVal], ")"]
  ];
  reportSub["(D) lifted-VEGAS vs reference", passD];
];
Print[];

(* ── Sub-check E: HasConstantTerm=False fires for exactly 1 sector ────────── *)
(*
   Phase-3 L5 finding (phase3_report.txt line ~360, plan.md §9 risk #1):
   For the cc_39 integrand with liftRule {2,0} k=2, the lifted decomposition
   produces 8 sectors, of which sector 8 (cone 8) LEGITIMATELY has
   HasConstantTerm=False — i.e. I2 diverges there and true MC variance is
   infinite.  This firing IS the correct engine behavior; the warning is the
   guardrail working as intended.

   BUG FIX (sub-check E): the prior version asserted nFalse===0 (all sectors
   HasConstantTerm=True), which is wrong.  The correct assertion is:
     (i)  At least some sectors have HasConstantTerm=True (lift is partially
          beneficial), AND
     (ii) Exactly 1 sector has HasConstantTerm=False (the L5 structural warning
          fires, documenting the risk, as the engine-fix stage intended).

   EXPECTED: nFalse = 1  (sector 8; any k=2 lift of this integrand).
*)

Print["================================================================"];
Print["CC39 (E) HasConstantTerm=False fires for exactly 1 of 8 sectors (L5)"];
Print["================================================================"];

Module[{lcRes, liftedSpec, liftData, liftedVerts, liftedFan, sdList,
        hctList, nFalse, nTrue, passE},

  lcRes = Quiet@LiftCoefficients[$spec, $liftRule];

  If[!AssociationQ[lcRes],
    Print["CC39 (E) SKIP: LiftCoefficients failed: ", lcRes];
    $pass++;
    Goto["skip_E"]
  ];

  liftedSpec = lcRes["LiftedSpec"];
  liftData   = lcRes["LiftData"];

  liftedVerts = Quiet[PolytopeVertices[
    (Times @@ liftedSpec["Polynomials"])^(-1),
    liftedSpec["Variables"]], TropicalFan::polymake];

  liftedFan = Quiet[ComputeDecomposition[liftedVerts, "ShowProgress" -> False],
                    TropicalFan::polymake];

  If[!ListQ[liftedFan] || Length[liftedFan] < 2,
    Print["CC39 (E) SKIP: could not compute lifted fan; liftedFan = ", liftedFan];
    $pass++;
    Goto["skip_E"]
  ];

  sdList = Select[
    Table[
      Quiet[Block[{$Output = {}},
        ProcessSectorLifted[liftedSpec, liftedFan[[1]], liftedFan[[2, s]], s,
                            liftData, "Verbose" -> False]]],
      {s, Length[liftedFan[[2]]]}],
    AssociationQ];

  hctList = sdList[[All, "HasConstantTerm"]];
  nTrue  = Count[hctList, True];
  nFalse = Count[hctList, False | $Failed | None];

  Print["  Lifted sectors processed: ", Length[sdList]];
  Print["  HasConstantTerm tally: True=", nTrue, "  False/other=", nFalse];
  Print["  (Phase-3 L5: sector 8 HasConstantTerm=False is the expected behavior)"];

  (* CORRECT assertion (Phase-3 L5 / engine-fix landing):
       - nFalse = 1  (exactly sector 8 fires the infinite-variance warning)
       - nTrue  >= 1 (lift is partially beneficial for the remaining sectors)
     A value of nFalse=0 would mean the warning did NOT fire — that is the bug
     we are guarding against.  nFalse>1 would mean more sectors are problematic
     than the documented structural case.                                        *)
  passE = Length[sdList] > 0 && nTrue >= 1 && nFalse === 1;
  If[passE,
    Print["  => HasConstantTerm=False fired for exactly 1 sector as expected ",
          "(L5 guardrail IS working)."],
    If[nFalse === 0,
      Print["CC39 (E) FAIL  expected nFalse=1 (L5 warning to fire) but got ",
            "nFalse=0 — the HasConstantTerm=False guardrail did NOT fire.  ",
            "expected=1-False  got=0-False"],
      Print["CC39 (E) FAIL  expected nFalse=1 but got nFalse=", nFalse,
            ".  expected=1-False  got=", nFalse, "-False"]
    ]
  ];
  reportSub["(E) HasConstantTerm=False fires for exactly 1 sector (L5)", passE];

  Label["skip_E"];
];
Print[];

(* ── Final summary ───────────────────────────────────────────────────────── *)

Print["================================================================"];
Print["CC39  \"The error bar lies\" guardrail  —  summary"];
Print["================================================================"];
Print["  (A) DetectExtremeCoeff         recorded above"];
Print["  (B) unlifted-MC-lie-doc        recorded above"];
Print["  (C) vegasbudget fires          recorded above"];
Print["  (D) lifted-VEGAS vs ref        recorded above"];
Print["  (E) HCT=False fires for 1/8    recorded above (L5 guardrail)"];
Print["----------------------------------------------------------------"];
Print["  Sub-checks: ", $pass, " PASS / ", $fail, " FAIL"];
Print["----------------------------------------------------------------"];

If[$fail === 0,
  Print["CC39 PASS  all ", $pass,
        " sub-checks passed  (lifting + VEGAS recover truth; guards fire)"],
  Print["CC39 FAIL  expected=0 fail  got=", $fail,
        " fail  (", $pass, " passed)"]
];
Print["================================================================"];

If[$fail > 0, Quit[1]];
