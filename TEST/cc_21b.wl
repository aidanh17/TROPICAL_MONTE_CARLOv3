(* ============================================================================
   TEST/cc_21b.wl  —  Cross-check #21b: complex-BASE + complex-exponent lifting
                       SplitRealImag == Direct == reference  (plan.md §8.2 #21b)

   What this tests
   ---------------
   A polynomial P with complex coefficients (so arg P != 0 on the domain) raised
   to a complex power B.  This is the exact regime that exposed BUG 1 in the bug
   log (TMCv2_BUG_LOG.md): the SplitRealImag VEGAS-variance mode silently dropped
   the magnitude factor exp(-Im(B)*arg P) because it used log|P| (modulus) rather
   than log P (complex log).  For real-positive P the two agree; for complex P the
   error can reach dozens of orders of magnitude.

   Five sub-checks are run in order:

     A  Unlifted: Direct vs NIntegrate (Tier-1 sanity, no CUBA needed, but we
        gate the whole check on CUBA so the VEGAS sub-check always runs together).

     B  Unlifted: SplitRealImag VEGAS vs Direct MC vs NIntegrate ref.
        SplitRealImag and Direct must agree within combined MC+VEGAS error AND
        both must match NIntegrate to < 1%.  This is the primary BUG-1 regression.

     C  Lifted: round-trip identity z -> z0 is EXACT (plan.md §6.2, #40 partial).
        Complex coefficient Cc = (1+I)*1e-4 triggers lifting; anchor z0 = |Cc|^(1/k)
        must be exact; residual c = Cc/z0^k must satisfy PossibleZeroQ.

     D  Lifted: SplitRealImag sector sum WITH MonoFactorLog matches ref via
        WL NIntegrate (no MC noise in the gate).  The MonoFactorLog term (L3 fix,
        lift_error_log.md) accounts for the log of the cleared real-positive
        monomial factor that tropical factoring removed.

     E  Lifted: dropping MonoFactorLog gives a WRONG answer (rel-err > 5%).
        Confirms the L3 term is both present and necessary, not a no-op.

   PASS criterion summary
   ----------------------
     A  rel-err(Direct vs NInteg) < 1e-3
     B  rel-err(Split vs NInteg) < 1e-2  AND  |Split-Direct| < sigma_combined + 1% ref
     C  PossibleZeroQ identity holds, z0^k === Abs[Cc], residual round-trip exact
     D  rel-err(lifted SplitRealImag with MFL vs ref) < 5e-3
     E  rel-err(lifted SplitRealImag without MFL vs ref) > 5e-2

   Tier 2: requires CUBA (for the VEGAS lane in sub-check B).
   If CUBA is absent, prints "CC21b SKIP (no CUBA)" and exits 0.

   Source / port notes
   -------------------
   Ported and unified from:
     OLD_CODE/TROPICAL_MONTE_CARLO/TMCv2_BUG_LOG.md   (BUG 1 root cause + fix)
     OLD_CODE/TROPICAL_MONTE_CARLO/lift_error_log.md   (L3 MonoFactorLog fix)
     TEST/test_complex_phase4.wl  #21b blocks          (Phase 4 implementation)
   Sub-checks D and E are a mechanical port of the liftedSplit[useMFL_] pattern
   in test_complex_phase4.wl, isolated as an independent named check.

   Run
   ---
     wolframscript -file TEST/cc_21b.wl
   (from any directory; absolute paths are derived from $InputFileName)
   ============================================================================ *)


(* ── 0.  Locate package root from this file's own path ───────────────────── *)

$pkgRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
$evalWL  = FileNameJoin[{$pkgRoot, "tropical_eval.wl"}];
$fanWL   = FileNameJoin[{$pkgRoot, "tropical_fan.wl"}];
$ioDir   = FileNameJoin[{$pkgRoot, "TEST", "INTERFILES", "cc21b"}];

If[!FileExistsQ[$evalWL],
  Print["CC21b FAIL  tropical_eval.wl not found at ", $evalWL];
  Quit[1]];

Get[$fanWL];
Get[$evalWL];

Quiet[CreateDirectory[$ioDir, CreateIntermediateDirectories -> True],
      {CreateDirectory::eexist}];


(* ── 1.  CUBA presence check (Tier-2 gate) ───────────────────────────────── *)

$cuba = TropicalEval`detectCuba[];
If[!TrueQ[$cuba["Found"]],
  Print["CC21b SKIP (no CUBA)"];
  Quit[0]];

Print["CC21b: tropical_eval.wl loaded; CUBA found at ", $cuba["IncludeDir"]];
Print["CC21b: INTERFILES directory: ", $ioDir];
Print[];


(* ── 2.  Shared helpers ──────────────────────────────────────────────────── *)

(* Relative error for complex values — max of |re-err| and |im-err|, normalized
   by |ref|.  Falls back to absolute error when ref is numerically tiny. *)
relErr[got_, ref_] :=
  If[Abs[N[ref]] < 1*^-300,
     Abs[N[got - ref]],
     Abs[N[got - ref]] / Abs[N[ref]]];

(* Accumulate pass/fail outcomes; print the final summary line. *)
$passes   = {};
$failures = {};

recordPass[name_String, detail_String : ""] :=
  (AppendTo[$passes, name];
   Print["CC21b PASS  ", name, If[detail =!= "", "  " <> detail, ""]]);

recordFail[name_String, extra_String : ""] :=
  (AppendTo[$failures, name];
   Print["CC21b FAIL  ", name, If[extra =!= "", "  " <> extra, ""]]);

record[name_String, ok_, detail_String : "", failExtra_String : ""] :=
  If[TrueQ[ok], recordPass[name, detail], recordFail[name, failExtra]];


(* ── 3.  Build the shared fan (unlifted) ─────────────────────────────────── *)
(* Integrand: P = 1 + x1^2 + (1+I)*x1*x2 + x2^2,  B = -(3/2 + 4I/5).
   The mixed term (1+I)*x1*x2 gives arg P != 0 on the positive orthant.
   The fan is built from the real part (Re B = -3/2) for tropical geometry.   *)

$vars  = {x[1], x[2]};
$poly  = 1 + x[1]^2 + (1 + I) x[1] x[2] + x[2]^2;
$Bex   = -(3/2 + 4 I / 5);

$spec = <|
  "Polynomials"         -> {$poly},
  "MonomialExponents"   -> {0, 0},
  "PolynomialExponents" -> {$Bex},
  "Variables"           -> $vars,
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> None
|>;

(* Build fan from |P|^(-Re B) — the shadow polynomial whose tropical geometry
   matches the convergence structure of Re B. *)
$shadowPoly = 1 + x[1]^2 + x[1] x[2] + x[2]^2;   (* |coeffs| only — real *)
$verts = PolytopeVertices[$shadowPoly^(-Re[$Bex]), $vars];
$fan   = ComputeDecomposition[$verts, "ShowProgress" -> False];

Print["Unlifted fan: ", Length[$fan[[2]]], " sectors"];
Print[];


(* ── 4.  NIntegrate reference (high-precision WL, no C++) ───────────────── *)

Print["Computing NIntegrate reference (WorkingPrecision 30) ..."];
$ref = NIntegrate[$poly^$Bex,
  {x[1], 0, Infinity}, {x[2], 0, Infinity},
  MaxRecursion -> 35, PrecisionGoal -> 8, WorkingPrecision -> 30];
Print["  ref = ", N[$ref]];
Print[];


(* ── 5A.  Unlifted: Direct MC vs NIntegrate ─────────────────────────────── *)
(* Gate: rel-err < 1e-3.  Verifies the basic complex-B pipeline. *)

Print["=== Sub-check A: unlifted Direct MC vs NIntegrate ==="];
$wdA = FileNameJoin[{$ioDir, "A_direct_mc"}];
Quiet[CreateDirectory[$wdA], {CreateDirectory::eexist}];

$mcA = TropicalEval`EvaluateTropicalMC[$spec, $fan, {{}},
  "Integrator"          -> "MC",
  "NSamples"            -> 4000000,
  "RunChecks"           -> False,
  "Verbose"             -> False,
  "WorkingDirectory"    -> $wdA,
  "ComplexExponentMode" -> "Direct"];

$vA = $mcA["Results"][[1]]["Re"] + I $mcA["Results"][[1]]["Im"];
$eA = $mcA["Results"][[1]]["ReErr"] + I $mcA["Results"][[1]]["ImErr"];
$rA = relErr[$vA, $ref];

Print["  Direct MC = ", N[$vA], " +/- ", N[$eA]];
Print["  rel-err vs NInteg = ", ScientificForm[$rA, 3], "  (gate < 1e-3)"];
record["A unlifted Direct MC vs NInteg",
  $rA < 1*^-3,
  "expected=<1e-3 got=" <> ToString[ScientificForm[$rA, 3]],
  "expected=<1e-3 got=" <> ToString[ScientificForm[$rA, 3]]];
Print[];


(* ── 5B.  Unlifted: SplitRealImag VEGAS vs Direct MC vs NIntegrate ───────── *)
(* Primary BUG-1 regression.  SplitRealImag uses the complex log P (the fix);
   before the fix this would be off by exp(Im(B)*arg P) ~ O(1) here but wildly
   wrong on the CALC2 bubble.  Gates: Split vs NInteg < 1%, Split vs Direct
   within combined sigma + 1% of |ref|.                                        *)

Print["=== Sub-check B: unlifted SplitRealImag VEGAS vs Direct MC vs NIntegrate ==="];
$wdB = FileNameJoin[{$ioDir, "B_split_vegas"}];
Quiet[CreateDirectory[$wdB], {CreateDirectory::eexist}];

$mcB = TropicalEval`EvaluateTropicalMC[$spec, $fan, {{}},
  "Integrator"          -> "VEGAS",
  "NSamples"            -> 2000000,
  "VegasEpsRel"         -> 1*^-9,   (* exhaust budget; don't stop early *)
  "RunChecks"           -> False,
  "Verbose"             -> False,
  "WorkingDirectory"    -> $wdB,
  "ComplexExponentMode" -> "SplitRealImag"];

$vB = $mcB["Results"][[1]]["Re"] + I $mcB["Results"][[1]]["Im"];
$eB = $mcB["Results"][[1]]["ReErr"] + I $mcB["Results"][[1]]["ImErr"];
$rBref   = relErr[$vB, $ref];
$rBdirect = relErr[$vB, $vA];
$sigmaCombined = Abs[$eA] + Abs[$eB];

Print["  SplitRealImag VEGAS = ", N[$vB], " +/- ", N[$eB]];
Print["  rel-err vs NInteg   = ", ScientificForm[$rBref,    3], "  (gate < 1e-2)"];
Print["  rel-err vs Direct   = ", ScientificForm[$rBdirect,  3]];
Print["  |Split-Direct|      = ", ScientificForm[Abs[N[$vB - $vA]], 3],
      "   sigma_combined = ", ScientificForm[$sigmaCombined, 3],
      "   1% |ref| = ", ScientificForm[0.01 Abs[N[$ref]], 3]];

$okBref   = $rBref   < 1*^-2;
$okBagree = Abs[N[$vB - $vA]] < $sigmaCombined + 0.01 Abs[N[$ref]];
record["B SplitRealImag VEGAS vs NInteg (BUG-1 regression)",
  $okBref,
  "expected=<1e-2 got=" <> ToString[ScientificForm[$rBref, 3]],
  "expected=<1e-2 got=" <> ToString[ScientificForm[$rBref, 3]]];
record["B SplitRealImag VEGAS vs Direct MC within sigma",
  $okBagree,
  "|got-direct|=" <> ToString[ScientificForm[Abs[N[$vB - $vA]], 3]] <>
    " sigma_tol=" <> ToString[ScientificForm[$sigmaCombined + 0.01 Abs[N[$ref]], 3]],
  "|got-direct|=" <> ToString[ScientificForm[Abs[N[$vB - $vA]], 3]] <>
    " > sigma_tol=" <> ToString[ScientificForm[$sigmaCombined + 0.01 Abs[N[$ref]], 3]]];
Print[];


(* ── 5C.  Lifted: round-trip identity (exact symbolic) ───────────────────── *)
(* Cc = (1+I)*1e-4 is a tiny complex coefficient; lifting replaces Cc*x1^2 by
   z^k * x1^2 with z0 = |Cc|^(1/k) and residual c = Cc/z0^k.  The round-trip
   z -> z0 in the lifted polynomial must recover the original exactly.         *)

Print["=== Sub-check C: lifted round-trip identity (exact, §6.2 / #40 partial) ==="];
$Cc    = (1 + I) 1*^-4;
$polyL = 1 + $Cc x[1]^2 + x[1] + x[2]^2;    (* poly to lift *)
$specL = <|
  "Polynomials"         -> {$polyL},
  "MonomialExponents"   -> {0, 0},
  "PolynomialExponents" -> {-(2 + I/2)},
  "Variables"           -> $vars,
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> None
|>;

$det    = TropicalEval`DetectExtremeCoefficients[$specL, 1000];
$kStar  = $det[[1]]["SuggestedK"];
$liftRules = {<|
  "PolyIndex"      -> 1,
  "ExponentVector" -> $det[[1]]["ExponentVector"],
  "k"              -> $kStar
|>};

$lifted = TropicalEval`LiftCoefficients[$specL, $liftRules];
$lspec  = $lifted["LiftedSpec"];
$ld     = $lifted["LiftData"];
$auxVar = $ld["AuxVariable"];
$z0     = $ld["z0"];

(* Check 1: z -> z0 in lifted poly == original poly (exact) *)
$subbed  = Expand[$lspec["Polynomials"][[1]] /. $auxVar -> $z0];
$idOK    = TrueQ[PossibleZeroQ[$subbed - Expand[$polyL]]];

(* Check 2: z0^k === |Cc| (anchor selection, exact) *)
$anchorOK = TrueQ[PossibleZeroQ[$z0^$kStar - Abs[$Cc]]];

(* Check 3: residual * z0^k === Cc (round-trip on the coefficient, exact) *)
$res        = $ld["Residuals"][[1]];
$residualOK = TrueQ[PossibleZeroQ[$res $z0^$kStar - $Cc]];

Print["  k* = ", $kStar, "   z0 = ", $z0];
Print["  |Cc| = ", Abs[$Cc], "   z0^k* = ", Simplify[$z0^$kStar]];
Print["  round-trip identity (z->z0): ", If[$idOK, "EXACT", "FAILED"]];
Print["  anchor   z0^k == |Cc|:       ", If[$anchorOK, "EXACT", "FAILED"]];
Print["  residual c*z0^k == Cc:       ", If[$residualOK, "EXACT", "FAILED"]];

record["C lifted round-trip identity z->z0 exact",    $idOK,
  "", "PossibleZeroQ[subbed - original] returned False"];
record["C lifted anchor z0^k == |Cc| exact",          $anchorOK,
  "", "z0^k != Abs[Cc] exactly"];
record["C lifted residual c*z0^k == Cc exact",        $residualOK,
  "", "residual*z0^k != Cc exactly"];
Print[];


(* ── 5D & 5E.  Lifted: SplitRealImag with / without MonoFactorLog ─────────
   The lifted polynomial has complex B, so we cannot use EvaluateTropicalMC
   directly (the Direct path refuses complex B on the lifted domain; SplitRealImag
   is the only admitted path).  Instead we reconstruct the SplitRealImag sector
   sum in WL via NIntegrate — this removes all MC noise from the gate.

   The lifted spec uses Re(B) for the sector geometry and carries Im(B) separately.
   Each sector contributes:

     mag_s  = prefactor * prod_j  Q_j^{Re(B_j)}
     phase_s = sum_j  I*Im(B_j) * [ log(Q_j) + MFL_j["Const"]
                                   + sum_i MFL_j["Coeffs"][[i]] * log(y_i) ]

   where Q_j is the flattened polynomial for the j-th factor and MFL is the
   MonoFactorLog association (plan.md §6.3 / lift_error_log L3).

   Sub-check D: with MFL   => rel-err < 5e-3 (the L3 fix must be correct).
   Sub-check E: without MFL => rel-err > 5e-2 (confirms MFL is load-bearing).
   --------------------------------------------------------------------------- *)

Print["=== Sub-checks D+E: lifted SplitRealImag WITH/WITHOUT MonoFactorLog ==="];

(* NIntegrate reference for the lifted spec. *)
$refL = NIntegrate[$polyL^(-(2 + I/2)),
  {x[1], 0, Infinity}, {x[2], 0, Infinity},
  MaxRecursion -> 40, PrecisionGoal -> 8];
Print["  NIntegrate reference (lifted poly) = ", N[$refL]];

(* Build fan for the lifted polynomial using the real shadow (no complex coeffs). *)
$shadowL = 1 + x[1]^2 + x[1] + x[2]^2;    (* |coeff| -> 1 for the extreme monomial *)
$vertsL  = PolytopeVertices[($shadowL + $auxVar^$kStar x[1]^2)^(-1), $lspec["Variables"]];
$fanL    = ComputeDecomposition[$vertsL, "ShowProgress" -> False];
Print["  Lifted fan: ", Length[$fanL[[2]]], " sectors"];

(* specRe: like lspec but with Re(B) so ProcessSectorLifted uses the right fan *)
$specRe  = $lspec;
$specRe["PolynomialExponents"] = Re[$lspec["PolynomialExponents"]];
$imB     = Im[$lspec["PolynomialExponents"]];
$nVars   = Length[$vars];   (* 2 — the y-integration variables after delta removal *)

(* Sector-by-sector WL integral, optionally including the MonoFactorLog term. *)
computeLiftedSplit[useMFL_] :=
  Sum[
    Module[{sd, flatRe, pfRe, BRe, mfl, dc, yv, lg, Q, mag, mterm, integ, lyps},
      sd = Quiet @ TropicalEval`ProcessSectorLifted[
        $specRe, $fanL[[1]], $fanL[[2, s]], s, $ld];
      Which[
        sd === $Failed,                     0,
        KeyExistsQ[sd, "EmptyDomain"],      0,
        True,
          flatRe = sd["FlattenedPolys"];
          pfRe   = sd["Prefactor"];
          BRe    = sd["PolynomialExponents"];
          mfl    = sd["MonoFactorLog"];
          dc     = sd["DomainConstraint"];
          (* fresh unique symbols for the y-integration variables *)
          yv     = Table[Unique["y"], {$nVars}];
          lg     = Log /@ yv;
          (* reconstruct the cleared (flattened) polynomial Q_j from its
             monomial-sum representation stored in FlattenedPolys *)
          Q      = Table[
            Total[(#[[1]] Exp[Total[#[[2]] lg]]) & /@ flatRe[[j]]],
            {j, Length[flatRe]}];
          (* real-exponent magnitude *)
          mag    = pfRe Product[Exp[BRe[[j]] Log[Q[[j]]]], {j, Length[flatRe]}];
          (* oscillatory phase: with or without MonoFactorLog correction *)
          mterm  = If[useMFL,
            Sum[I $imB[[j]] (Log[Q[[j]]]
                              + mfl[[j]]["Const"]
                              + Sum[mfl[[j]]["Coeffs"][[i]] lg[[i]], {i, $nVars}]),
                {j, Length[flatRe]}],
            Sum[I $imB[[j]] Log[Q[[j]]],
                {j, Length[flatRe]}]];
          integ  = mag Exp[mterm];
          (* apply delta-constraint indicator if present *)
          If[dc =!= None,
            lyps  = (N[dc["LogZ0"]] - Total[N[dc["IndicatorCoeffs"]] lg]) / N[dc["MP"]];
            integ = integ Boole[lyps <= 0]];
          Quiet @ NIntegrate[Evaluate[integ],
            Evaluate[Sequence @@ ({#, 0, 1} & /@ yv)],
            MaxRecursion -> 18, PrecisionGoal -> 4, Method -> "GlobalAdaptive"]
      ]],
    {s, Length[$fanL[[2]]]}];

Print["  Computing lifted SplitRealImag WITH MonoFactorLog (D) ..."];
$withMFL = computeLiftedSplit[True];
$rD      = relErr[$withMFL, $refL];
Print["  WITH  MFL = ", N[$withMFL], "  rel-err = ", ScientificForm[$rD, 3],
      "  (gate < 5e-3)"];

Print["  Computing lifted SplitRealImag WITHOUT MonoFactorLog (E) ..."];
$noMFL   = computeLiftedSplit[False];
$rE      = relErr[$noMFL, $refL];
Print["  WITHOUT MFL = ", N[$noMFL], "  rel-err = ", ScientificForm[$rE, 3],
      "  (must be > 5e-2 to confirm MFL is load-bearing)"];

record["D lifted SplitRealImag WITH MonoFactorLog vs ref",
  $rD < 5*^-3,
  "expected=<5e-3 got=" <> ToString[ScientificForm[$rD, 3]],
  "expected=<5e-3 got=" <> ToString[ScientificForm[$rD, 3]]];
record["E MonoFactorLog is load-bearing (dropping it is wrong)",
  $rE > 5*^-2,
  "expected=>5e-2 got=" <> ToString[ScientificForm[$rE, 3]],
  "expected=>5e-2 got=" <> ToString[ScientificForm[$rE, 3]]];
Print[];


(* ── 6.  Final summary ───────────────────────────────────────────────────── *)

$nPass = Length[$passes];
$nFail = Length[$failures];
$nTot  = $nPass + $nFail;

Print["================================================================"];
If[$nFail == 0,
  Print["CC21b PASS  all ", $nTot, " sub-checks passed",
        "  (Direct/Split/ref agree; MonoFactorLog present and correct)"],
  Print["CC21b FAIL  ", $nFail, " of ", $nTot, " sub-checks failed:",
        "  ", StringRiffle[$failures, ", "]]];
Print["================================================================"];

Quit[If[$nFail > 0, 1, 0]]
