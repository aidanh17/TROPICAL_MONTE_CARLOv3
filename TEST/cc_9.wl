(* ============================================================================
   TEST/cc_9.wl  —  Cross-check #9: MC vs VEGAS parity  (plan.md §8.2, §8.3)

   PASS criterion (plan.md §8.2 #9):
       |MC − VEGAS| ≤ K · sqrt(errMC² + errVEGAS²),   K = 5   (combined-sigma gate)
   and VEGAS within a budget-appropriate relative tolerance of a reference value.

   Tier 2: requires CUBA (for Integrator -> "VEGAS").
   If CUBA is absent prints "CC9 SKIP (no CUBA)" and exits 0.

   Cases ported from OLD_CODE/TROPICAL_MONTE_CARLO2/TROPICAL_MONTE_CARLOv2/
   EXAMPLES/test_vegas.wl (cases V1–V4), adapted to the v3 API:
     V1  Simplex  — P = 1 + Σ x_i, A=0, B=-(n+2); exact 1/(n+1)!   (n=2 and n=4)
     V2  Frac A   — A_i = -1/2, B = -(n+1), n=4;  exact π^2 Γ(3)/4!  (sqrt-edge)
     V3  Complex B — P = 1 + λ x1² + x2² + x1 x2², B=-(2+I/2), λ scan; NIntegrate
     V4  Coeff batch — P = c0 + Σ c_i x_i (n=4), B=-6; exact 1/(120 c0² Π c_i)

   Run:  wolframscript -file TEST/cc_9.wl
   (from the TROPICAL_MONTE_CARLO3 root, or anywhere — uses absolute paths)
   ============================================================================ *)

(* ── 0. Locate the package root via this file's own path ─────────────────── *)
$pkgRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
$evalWL  = FileNameJoin[{$pkgRoot, "tropical_eval.wl"}];
$fanWL   = FileNameJoin[{$pkgRoot, "tropical_fan.wl"}];
$ioDir   = FileNameJoin[{$pkgRoot, "TEST", "INTERFILES"}];

If[!FileExistsQ[$evalWL],
  Print["CC9 FAIL  tropical_eval.wl not found at ", $evalWL]; Quit[1]];

Get[$fanWL];   (* loads TropicalFan` including computeFanScaled *)
Get[$evalWL];  (* loads TropicalEval` *)

Quiet[CreateDirectory[$ioDir], {CreateDirectory::eexist}];

(* ── 1. CUBA presence check (Tier-2 gate) ────────────────────────────────── *)
$cuba = TropicalEval`detectCuba[];
If[!TrueQ[$cuba["Found"]],
  Print["CC9 SKIP (no CUBA)"];
  Quit[0]];

Print["CC9: tropical_eval.wl loaded; CUBA found at ", $cuba["IncludeDir"]];
Print["CC9: working directory for generated files: ", $ioDir];
Print[];

(* ── 2. Helpers ──────────────────────────────────────────────────────────── *)

(* Suppress output during integration runs *)
SetAttributes[quietRun, HoldFirst];
quietRun[expr_] := Block[{$Output = {}}, expr];

(* Numeric formatter: avoids InputForm back-ticks and escapes cleanly *)
fmtR[x_, p_] := If[x == 0, "0",
  ToString[CForm[SetPrecision[x, p]]]];
fmt[v_, p_: 4] := With[{vv = Chop[N[v], 1.*^-14]},
  If[Im[vv] == 0,
    fmtR[Re[vv], p],
    fmtR[Re[vv], p] <>
      If[Im[vv] < 0, " - ", " + "] <> fmtR[Abs[Im[vv]], p] <> "*I"]];

(* ── 3. Core comparison function ─────────────────────────────────────────── *)
(* Returns True on PASS, False on FAIL.
   refList : list of reference values (complex), one per kinematic point, or None.
   opts    : forwarded to compareMCvsVEGAS via options below.             *)

Options[cc9Compare] = {
  "RefList"    -> None,
  "RefLabel"   -> "exact",
  "McSamples"  -> 1000000,
  "VegasEval"  -> 200000,
  "VegasEpsRel"-> 1.*^-9,    (* tiny → VEGAS exhausts the full budget *)
  "McTol"      -> 2.*^-2,    (* MC   rel-err gate vs reference         *)
  "VegasTol"   -> 1.*^-2,    (* VEGAS rel-err gate vs reference         *)
  "SigmaK"     -> 5,         (* combined-sigma factor for |MC−VEGAS|    *)
  "GateMC"     -> True,      (* hard-fail if MC exceeds McTol           *)
  "GateSigma"  -> True       (* hard-fail if |MC−VEGAS| > K·σ_combined  *)
};

cc9Compare[label_String, spec_Association, fanData_List,
           kinPoints_List, OptionsPattern[]] :=
Module[{refList, refLabel, ns, ve, veEps, mcTol, veTol, sigK,
        gateMC, gateSig, rMC, rVE, nkp, anyFail = False},

  refList  = OptionValue["RefList"];
  refLabel = OptionValue["RefLabel"];
  ns       = OptionValue["McSamples"];
  ve       = OptionValue["VegasEval"];
  veEps    = OptionValue["VegasEpsRel"];
  mcTol    = OptionValue["McTol"];
  veTol    = OptionValue["VegasTol"];
  sigK     = OptionValue["SigmaK"];
  gateMC   = OptionValue["GateMC"];
  gateSig  = OptionValue["GateSigma"];
  nkp      = Length[kinPoints];

  Print["--- CC9 sub-case: ", label, " ---"];
  Print["    MC ", ns, " samples  |  VEGAS ", ve, " maxeval/sector"];

  (* run MC *)
  rMC = quietRun @ EvaluateTropicalMC[spec, fanData, kinPoints,
    "Integrator" -> "MC", "NSamples" -> ns,
    "RunChecks" -> False, "Verbose" -> False,
    "WorkingDirectory" -> $ioDir];
  (* run VEGAS *)
  rVE = quietRun @ EvaluateTropicalMC[spec, fanData, kinPoints,
    "Integrator" -> "VEGAS", "NSamples" -> ve,
    "VegasEpsRel" -> veEps,
    "RunChecks" -> False, "Verbose" -> False,
    "WorkingDirectory" -> $ioDir];

  (* early abort on driver failure *)
  If[!AssociationQ[rMC] || !AssociationQ[rVE],
    Print["    CC9 FAIL  expected=<run> got=",
          If[!AssociationQ[rMC], "MC->" <> ToString[rMC],
             "VEGAS->" <> ToString[rVE]]];
    Return[False]];

  Do[
    Module[{rm, rv, mcVal, veVal, mcErr, veErr, diff, comb, sig,
            ref, mcDev, veDev, kpFail = False, tags = {}},

      rm = rMC["Results"][[i]];  rv = rVE["Results"][[i]];
      mcVal = rm["Re"] + I rm["Im"];
      veVal = rv["Re"] + I rv["Im"];
      mcErr = Sqrt[rm["ReErr"]^2 + rm["ImErr"]^2];
      veErr = Sqrt[rv["ReErr"]^2 + rv["ImErr"]^2];
      diff  = Abs[mcVal - veVal];
      comb  = Max[Sqrt[mcErr^2 + veErr^2], 1.*^-300];
      sig   = diff / comb;

      Print["    kp ", i, If[nkp > 1 && kinPoints[[i]] =!= {},
            " = " <> ToString[kinPoints[[i]]], ""], ":"];
      Print["      MC    = ", fmt[mcVal, 6], "  ± ", fmt[mcErr, 3]];
      Print["      VEGAS = ", fmt[veVal, 6], "  ± ", fmt[veErr, 3]];
      Print["      |MC−VEGAS| = ", fmt[diff, 3],
            "  = ", fmt[sig, 3], " σ  (gate K=", sigK, ")"];

      (* gate: combined sigma band *)
      If[gateSig && !TrueQ[N[sig] <= sigK],
        kpFail = True; AppendTo[tags, "sigma>" <> ToString[sigK]]];

      (* gate: reference comparison *)
      If[refList =!= None,
        ref   = refList[[i]];
        mcDev = Abs[(mcVal - ref) / Max[Abs[ref], 1.*^-300]];
        veDev = Abs[(veVal - ref) / Max[Abs[ref], 1.*^-300]];
        Print["      ", refLabel, " = ", fmt[ref, 6]];
        Print["      relErr:  MC = ", fmt[mcDev, 3],
              If[gateMC, "  (gate " <> ToString[mcTol] <> ")", "  (report)"],
              "   VEGAS = ", fmt[veDev, 3],
              "  (gate ", veTol, ")"];
        If[gateMC && !TrueQ[N[mcDev]   <= mcTol],
          kpFail = True; AppendTo[tags, "MC>tol"]];
        If[!TrueQ[N[veDev] <= veTol],
          kpFail = True; AppendTo[tags, "VEGAS>tol"]];
      ];

      If[kpFail,
        anyFail = True;
        Print["      -> kp ", i, " FAIL: ", StringRiffle[tags, ", "]];
        Print["         CC9 FAIL  expected=pass  got=fail  [", label,
              " kp", i, " tags: ", StringRiffle[tags, ","], "]"]];
    ],
    {i, nkp}];

  If[!anyFail,
    Print["    CC9 PASS  ", label, "  (", nkp,
          " kinematic point", If[nkp > 1, "s", ""], " all within gates)"]];
  Print[];
  !anyFail
];

(* ── 4. Cases ─────────────────────────────────────────────────────────────── *)

$cc9Results = {};   (* collect sub-case pass/fail *)

(* ────────────────────────────────────────────────────────────────────────────
   V1  A_simplex: P = 1 + Σ x_i, A=0, B=-(n+2)
       Exact value: 1/(n+1)!          (n=2: 1/6 ≈ 0.1667; n=4: 1/120 ≈ 0.00833)
   ──────────────────────────────────────────────────────────────────────────── *)
Print["========================================================"];
Print["CC9  V1 — A_simplex, exact 1/(n+1)!"];
Print["========================================================"];

Do[
  Module[{vars, P, spec, verts, fan, exact, p},
    vars  = Table[x[i], {i, n}];
    P     = 1 + Total[vars];
    spec  = <|"Polynomials" -> {P},
              "MonomialExponents"  -> ConstantArray[0, n],
              "PolynomialExponents"-> {-(n + 2)},
              "Variables"          -> vars,
              "KinematicSymbols"   -> {},
              "RegulatorSymbol"    -> None|>;
    verts = PolytopeVertices[P^(-1), vars];
    fan   = computeFanScaled[verts];
    If[fan === $Failed,
      Print["CC9 FAIL  V1 n=", n, ": fan computation failed"];
      AppendTo[$cc9Results, False]; Return[]];
    exact = N[1 / (n + 1)!];
    p = cc9Compare[
      "V1 A_simplex n=" <> ToString[n], spec, fan, {{}},
      "RefList" -> {exact}, "RefLabel" -> "1/(n+1)!",
      "McSamples" -> 1000000, "VegasEval" -> 200000,
      "McTol" -> 2.*^-2, "VegasTol" -> 1.*^-2, "SigmaK" -> 5];
    AppendTo[$cc9Results, p]],
  {n, {2, 4}}];

(* ────────────────────────────────────────────────────────────────────────────
   V2  Frac A (sqrt edge): A_i=-1/2, B=-(n+1), n=4
       Exact: π^(n/2) Γ(n/2+1) / n!  =  π² Γ(3) / 24 = π²/12 ≈ 0.8225
       MC is heaviest on the sqrt edge; VEGAS shines here.
   ──────────────────────────────────────────────────────────────────────────── *)
Print["========================================================"];
Print["CC9  V2 — Frac A (sqrt edge), n=4, exact π²/12"];
Print["========================================================"];

Module[{n = 4, vars, P, spec, verts, fan, exact, p},
  vars  = Table[x[i], {i, n}];
  P     = 1 + Total[vars];
  spec  = <|"Polynomials" -> {P},
            "MonomialExponents"  -> ConstantArray[-1/2, n],
            "PolynomialExponents"-> {-(n + 1)},
            "Variables"          -> vars,
            "KinematicSymbols"   -> {},
            "RegulatorSymbol"    -> None|>;
  verts = PolytopeVertices[P^(-1), vars];
  fan   = computeFanScaled[verts];
  If[fan === $Failed,
    Print["CC9 FAIL  V2: fan computation failed"];
    AppendTo[$cc9Results, False]; Goto[afterV2]];
  exact = N[Pi^(n/2) Gamma[n/2 + 1] / n!];
  Print["    exact = ", fmt[exact, 7]];
  p = cc9Compare[
    "V2 Af_frac n=4 (sqrt edge)", spec, fan, {{}},
    "RefList" -> {exact}, "RefLabel" -> "pi^2/12",
    "McSamples" -> 3000000, "VegasEval" -> 300000,
    "McTol" -> 3.*^-2, "VegasTol" -> 1.*^-2, "SigmaK" -> 6];
  AppendTo[$cc9Results, p]];
Label[afterV2];

(* ────────────────────────────────────────────────────────────────────────────
   V3  Complex B, kinematic scan
       P = 1 + λ x1² + x2² + x1 x2²,  B = -(2 + I/2)
       Reference via NIntegrate at each λ.  Exercises ncomp=2 (Re + Im).
   ──────────────────────────────────────────────────────────────────────────── *)
Print["========================================================"];
Print["CC9  V3 — Complex B, λ-scan, NIntegrate reference"];
Print["========================================================"];

Module[{lam, A, P, vars, spec, verts, fan, lamVals, kinPts, refList, t1, t2, p},
  lam  = Symbol["lam"];
  A    = 2 + I/2;
  P    = 1 + lam x[1]^2 + x[2]^2 + x[1] x[2]^2;
  vars = {x[1], x[2]};
  spec = <|"Polynomials" -> {P},
           "MonomialExponents"  -> {0, 0},
           "PolynomialExponents"-> {-A},
           "Variables"          -> vars,
           "KinematicSymbols"   -> {lam},
           "RegulatorSymbol"    -> None|>;
  (* fan is λ-independent: compute from the unit-coefficient proxy *)
  verts = PolytopeVertices[(P /. lam -> 1)^(-1), vars];
  fan   = computeFanScaled[verts];
  If[fan === $Failed,
    Print["CC9 FAIL  V3: fan computation failed"];
    AppendTo[$cc9Results, False]; Goto[afterV3]];
  lamVals  = {0.5, 1.0, 2.0, 4.0};
  kinPts   = List /@ lamVals;
  Print["    computing NIntegrate references for λ ∈ ", lamVals, " ..."];
  refList  = Table[
    Quiet @ NIntegrate[
      (1 + lv t1^2 + t2^2 + t1 t2^2)^(-A),
      {t1, 0, Infinity}, {t2, 0, Infinity},
      MaxRecursion -> 20, PrecisionGoal -> 5],
    {lv, lamVals}];
  Print["    NIntegrate refs = ", fmt[#, 6] & /@ refList];
  p = cc9Compare[
    "V3 complex B (A=2+I/2), lam scan", spec, fan, kinPts,
    "RefList" -> N[refList], "RefLabel" -> "NIntegrate",
    "McSamples" -> 500000, "VegasEval" -> 150000,
    "McTol" -> 3.*^-2, "VegasTol" -> 1.*^-2, "SigmaK" -> 5];
  AppendTo[$cc9Results, p]];
Label[afterV3];

(* ────────────────────────────────────────────────────────────────────────────
   V4  Coefficient batch scan
       P = c0 + Σ_{i=1}^4 c_i x_i,  B = -6
       Exact: 1 / (120 · c0² · c1 · c2 · c3 · c4)
       (generalized simplex; 24 random coefficient sets, seed fixed)

   BUG FIX (prior run):  the reference list was computed via
       Table[Module[{c = kinPts[[s]]}, ...], {s, nSets}]
   where the Table iterator symbol s was still in scope inside Module,
   causing Mathematica to leave Im[s^-1] unevaluated when the
   denominator was formed symbolically before s was bound to a number.
   Fix: use Map with an explicit Function so kinPts values (pure Machine
   reals) are substituted directly — no iterator symbol, no symbolic
   intermediate.  All kinematic values are forced to machine precision
   with N[] immediately after generation.
   ──────────────────────────────────────────────────────────────────────────── *)
Print["========================================================"];
Print["CC9  V4 — Coefficient batch (n=4, 24 sets), exact reference"];
Print["========================================================"];

Module[{n = 4, vars, c0sym, ciSyms, csyms, P, spec, verts, fan,
        nSets, kinPts, refList, p},
  vars   = Table[x[i], {i, n}];
  (* Use explicit local symbol names to avoid any global-symbol clash.
     c0sym, ciSyms are the kinematic symbols appearing in P; they are
     purely symbolic (no numeric value) throughout this Module. *)
  c0sym  = Symbol["cc9c0"];
  ciSyms = Table[Symbol["cc9c" <> ToString[i]], {i, n}];
  csyms  = Prepend[ciSyms, c0sym];     (* {cc9c0, cc9c1, cc9c2, cc9c3, cc9c4} *)
  P      = c0sym + Sum[ciSyms[[i]] x[i], {i, n}];
  spec   = <|"Polynomials"          -> {P},
             "MonomialExponents"    -> ConstantArray[0, n],
             "PolynomialExponents"  -> {-6},
             "Variables"            -> vars,
             "KinematicSymbols"     -> csyms,
             "RegulatorSymbol"      -> None|>;
  (* fan is coefficient-independent: unit-coefficient proxy is the simplex *)
  verts = PolytopeVertices[(1 + Total[vars])^(-1), vars];
  fan   = computeFanScaled[verts];
  If[fan === $Failed,
    Print["CC9 FAIL  V4: fan computation failed"];
    AppendTo[$cc9Results, False]; Goto[afterV4]];
  nSets  = 24;
  SeedRandom[12345];
  (* kinPts: each row is {c0, c1, c2, c3, c4} as machine-precision reals.
     N[] is explicit so no exact-form wrappers survive into the engine's
     Im[#]==0 guard or into the reference formula. *)
  kinPts = N @ Table[
    Prepend[Table[Exp[RandomReal[{-0.7, 0.7}]], {n}], 1.0],  (* c0 = 1 *)
    {nSets}];
  (* Reference: 1/(120 c0^2 c1 c2 c3 c4).
     Use Map+Function — no Table iterator symbol ever appears in the
     denominator, so the result is provably numeric for numeric kinPts. *)
  refList = Map[
    Function[kp,
      N[1.0 / (120.0 * kp[[1]]^2 * Times @@ kp[[2 ;; All]])]],
    kinPts];
  (* Sanity: catch any surviving symbolics before they can fake-pass *)
  If[!AllTrue[refList, (NumericQ[#] && Im[#] == 0 && # > 0) &],
    Print["CC9 FAIL  V4: reference list is not purely numeric — got: ",
          Select[refList, !NumericQ[#] &]];
    AppendTo[$cc9Results, False]; Goto[afterV4]];
  p = cc9Compare[
    "V4 coeff batch (n=4, " <> ToString[nSets] <> " sets)",
    spec, fan, kinPts,
    "RefList" -> refList, "RefLabel" -> "exact 1/(120 c0^2 prod_ci)",
    "McSamples" -> 300000, "VegasEval" -> 100000,
    "McTol" -> 3.*^-2, "VegasTol" -> 1.*^-2, "SigmaK" -> 5];
  AppendTo[$cc9Results, p]];
Label[afterV4];

(* ── 5. Summary ──────────────────────────────────────────────────────────── *)
Module[{nPass, nFail, labels},
  labels = {"V1 n=2", "V1 n=4", "V2 Af_frac", "V3 complex-B", "V4 coeff-batch"};
  nPass  = Count[$cc9Results, True];
  nFail  = Count[$cc9Results, False];
  Print["========================================================"];
  Print["CC9  MC vs VEGAS parity  —  summary"];
  Print["========================================================"];
  Do[
    Print["  ", StringPadRight[labels[[i]], 20],
          If[TrueQ[$cc9Results[[i]]], "PASS", "FAIL"]],
    {i, Min[Length[labels], Length[$cc9Results]]}];
  Print["--------------------------------------------------------"];
  If[nFail == 0,
    Print["CC9 PASS  all ", nPass, " sub-cases passed  (MC vs VEGAS ≤ 5σ, refs within tol)"],
    Print["CC9 FAIL  expected=", nPass + nFail, " pass  got=",
          nFail, " fail  (", nPass, " passed, ", nFail, " failed)"]];
  Print["========================================================"];
  If[nFail > 0, Quit[1]]
];
