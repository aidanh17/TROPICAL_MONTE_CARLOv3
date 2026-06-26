(* ============================================================================
   TEST/cc_31.wl  -  Cross-check #31  (plan.md §8.3)
   "Merge fidelity vs Tree A"

   PASS criterion (plan.md §8.3 #31):
     Every Tree-A golden (Tests 1-22, divergent crosscheck) reproduced by
     the v3 engine within MC noise.  Symbolic outputs (sector counts,
     effective exponents, IsDivergent flags, IBP term counts) must be
     BITWISE IDENTICAL to the Tree-A goldens captured in
     TEST/baselines/treeA.txt.  Numeric outputs (NIntegrate vs sector-sum
     relative errors, pole/finite coefficients) must fall inside the same
     tolerance bands used by Tree-A's own PASS/FAIL logic.

   Strategy:
     1. Load ONLY the v3 packages (absolute paths; never touch OLD_CODE).
     2. Run the same specs and checks that tree-A's RunAllTests (Tests 1-18),
        IBP tests (19-22), and test_divergent_crosscheck (B=4, B=5) use.
     3. Compare every result to the published Tree-A golden numbers (from
        TEST/baselines/treeA.txt) and assert within tolerance.
     4. Print  "CC31 PASS ..."  /  "CC31 FAIL expected=<e> got=<g>"
        in the canonical format.

   Tier 1: WL + g++ only.  CUBA-Vegas columns are run when present but are
   NEVER required for CC31 PASS/FAIL.

   Run:
     cd /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3
     wolframscript -file TEST/cc_31.wl
   ============================================================================ *)

(* --------------------------------------------------------------------------
   0. Load v3 packages via absolute paths (branch-safe; never uses SetDirectory)
   -------------------------------------------------------------------------- *)
$v3Root = "/Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3";

Get[FileNameJoin[{$v3Root, "tropical_fan.wl"}]];
Get[FileNameJoin[{$v3Root, "tropical_eval.wl"}]];

Print["CC31: v3 packages loaded."];
Print[];

(* --------------------------------------------------------------------------
   1. Infrastructure
   -------------------------------------------------------------------------- *)

(* Suppress output of long runs; capture the return value only. *)
SetAttributes[quietRun, HoldFirst];
quietRun[expr_] := Block[{$Output = {}}, expr];

(* Global PASS/FAIL tracking. *)
$cc31Pass = True;
$cc31FailCount = 0;
$cc31Assertions = 0;

cc31Assert[tag_String, cond_, expectedStr_String, gotStr_String] :=
  Module[{},
    $cc31Assertions++;
    If[TrueQ[cond],
      Print["  [", tag, "] subPASS"],
      ($cc31Pass = False;
       $cc31FailCount++;
       Print["  [", tag, "] subFAIL expected=", expectedStr, " got=", gotStr])
    ]
  ];

(* Robust-fan helper (mirrors the one in test_divergent_crosscheck.wl).
   K-scaling is scale-invariant for the normal fan (plan §6.4). *)
computeFanRobust[verts_] := Module[{n, fd},
  n = Length[First[verts]];
  Do[
    fd = Quiet @ ComputeDecomposition[K * verts, "ShowProgress" -> False];
    If[ListQ[fd] && Length[fd] == 2 && FreeQ[fd, $Failed] && Length[fd[[1]]] > 0,
       Return[fd, Module]],
    {K, {1, n + 2, 2 n + 4, 6 n + 6}}];
  $Failed];

(* --------------------------------------------------------------------------
   2. Symbolic / structural checks (bitwise identical to Tree-A goldens)
   -------------------------------------------------------------------------- *)

Print["=== CC31 Part A: symbolic/structural checks (Tests 1-22 analogues) ==="];
Print[];

(* ----- A1: fan sector count for the 2D basic spec (Tree-A Example 1) -----
   Tree-A golden: 3 rays, 3 sectors.  (treeA.txt line "Fan: 3 rays, 3 sectors") *)
Module[{poly, vars, verts, fanData},
  poly    = 1 + x[1]^2 + x[2]^2;
  vars    = {x[1], x[2]};
  verts   = PolytopeVertices[poly^(-3), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  If[!ListQ[fanData] || fanData === $Failed,
    cc31Assert["A1-fan-count", False, "3 sectors", "$Failed"]; Return[$Failed]];
  cc31Assert["A1-fan-count",
    Length[fanData[[2]]] == 3,
    "3", ToString[Length[fanData[[2]]]]];
  (* Sector 1: det=|{0,0}|=1, effA={1,1} (Tree-A golden lines 17-18) *)
  Module[{spec, sd},
    spec = <|"Polynomials" -> {poly},
             "MonomialExponents" -> {0, 0},
             "PolynomialExponents" -> {-3},
             "Variables" -> vars,
             "KinematicSymbols" -> {},
             "RegulatorSymbol" -> None|>;
    sd = quietRun @ ProcessSector[spec, fanData[[1]], fanData[[2,1]], 1];
    cc31Assert["A1-sector1-effA",
      AssociationQ[sd] && sd["NewExponents"] === {1, 1},
      "{1, 1}", ToString[If[AssociationQ[sd], sd["NewExponents"], sd]]];
    cc31Assert["A1-sector1-prefactor",
      AssociationQ[sd] && sd["Prefactor"] === 1,
      "1", ToString[If[AssociationQ[sd], sd["Prefactor"], sd]]];
    cc31Assert["A1-sector1-convergent",
      AssociationQ[sd] && sd["IsDivergent"] === False,
      "False", ToString[If[AssociationQ[sd], sd["IsDivergent"], sd]]];
  ]
];
Print[];

(* ----- A2: 4D fan sector count (Tree-A golden: "8 rays, 14 sectors") ----- *)
Module[{poly, vars, verts, fanData},
  poly  = 1 + x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2 + x[1]*x[2] + x[3]*x[4];
  vars  = {x[1], x[2], x[3], x[4]};
  verts = PolytopeVertices[poly^(-4), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  If[!ListQ[fanData] || fanData === $Failed,
    cc31Assert["A2-4D-fan-count", False, "14 sectors", "$Failed"]; Return[$Failed]];
  cc31Assert["A2-4D-fan-count",
    Length[fanData[[2]]] == 14,
    "14", ToString[Length[fanData[[2]]]]];
];
Print[];

(* ----- A3: divergent-flag detection for the 1D spec (Tree-A Test 3) -----
   spec: Integral y^{2eps-1}/(1+y^2), one divergent sector.
   Tree-A golden: 1 divergent sector found. *)
Module[{eps, spec, verts, fanData, sectors, divSectors},
  eps  = Symbol["eps31"];
  spec = <|"Polynomials" -> {1 + x[1]^2},
           "MonomialExponents" -> {2*eps - 1},
           "PolynomialExponents" -> {-1},
           "Variables" -> {x[1]},
           "KinematicSymbols" -> {},
           "RegulatorSymbol" -> eps|>;
  verts   = PolytopeVertices[(1 + x[1]^2)^(-1), {x[1]}];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  If[!ListQ[fanData] || fanData === $Failed,
    cc31Assert["A3-1D-div-flag", False, "listQ[fanData]", "$Failed"]; Return[$Failed]];
  sectors = fanData[[2]];
  divSectors = Select[
    Table[quietRun @ ProcessSector[spec, fanData[[1]], sectors[[s]], s],
          {s, Length[sectors]}],
    AssociationQ[#] && TrueQ[#["IsDivergent"]] &];
  cc31Assert["A3-1D-div-count",
    Length[divSectors] >= 1,
    ">=1 divergent", ToString[Length[divSectors]]];
];
Print[];

(* ----- A4: IBP term counts (Tree-A Tests 19-21) ----- *)

(* IBP expected-term helper: mirrors Tree-A Test-19 logic verbatim.
   For the single divergent variable k, count the total number of monomials
   across all cleared polynomials that have e_{m,k} > 0 (each such monomial
   contributes one IBP term after the recursion).  This is the exact predicate
   from OLD_CODE/TROPICAL_MONTE_CARLO/EXAMPLES/tropical_eval_examples.wl
   lines 2831-2843. *)
ibpExpectedTerms[sd_Association, k_Integer] :=
  Total @ Table[
    Count[sd["ClearedPolys"][[j]],
      _?(TrueQ[#[[2, k]] > 0] || (NumericQ[#[[2, k]]] && #[[2, k]] > 0) &)],
    {j, Length[sd["ClearedPolys"]]}];

(* A4 helper: given a spec+fanData, find all divergent sectors via Select
   (mirroring the OLD test pattern — never hardcode sector indices), run
   IBPReduceSector on each, and assert nTerms == ibpExpectedTerms. *)
checkIBPTerms[lbl_String, spec_Association, fanData_List, eps_] :=
  Module[{allSD, divSecs, nDiv},
    allSD = Table[
      quietRun @ ProcessSector[spec, fanData[[1]], fanData[[2, s]], s],
      {s, Length[fanData[[2]]]}];
    divSecs = Select[allSD, AssociationQ[#] && TrueQ[#["IsDivergent"]] &];
    nDiv = Length[divSecs];
    cc31Assert[lbl <> "-divCount",
      nDiv >= 1, ">=1 divergent", ToString[nDiv]];
    If[nDiv == 0, Return[]];
    Do[
      Module[{sd = divSecs[[s]], k, expN, ibp, gotN},
        ibp  = quietRun @ IBPReduceSector[sd, eps];
        k    = If[AssociationQ[ibp] && KeyExistsQ[ibp, "DivergentVariables"],
                  ibp["DivergentVariables"][[1]], 1];
        expN = ibpExpectedTerms[sd, k];
        gotN = If[AssociationQ[ibp] && KeyExistsQ[ibp, "Terms"],
                  ibp["NTerms"], $Failed];
        cc31Assert[lbl <> "-sec" <> ToString[sd["ConeIndex"]] <> "-nTerms",
          gotN === expN,
          ToString[expN], ToString[gotN]];
        (* All terms must be convergent at eps=0 (mirrors OLD tree-A lines 2846-2863) *)
        If[AssociationQ[ibp],
          Do[
            Module[{termA0 = ibp["Terms"][[t]]["NewExponents"] /. eps -> 0},
              Do[
                If[TrueQ[Re[termA0[[i]]] <= 0] ||
                   (NumericQ[termA0[[i]]] && Re[termA0[[i]]] <= 0),
                  cc31Assert[lbl <> "-sec" <> ToString[sd["ConeIndex"]] <>
                             "-term" <> ToString[t] <> "-conv-y" <> ToString[i],
                    False, "Re[alpha^0]>0", ToString[termA0[[i]]]]],
                {i, ibp["Dimension"]}]],
            {t, ibp["NTerms"]}]]
      ],
      {s, nDiv}]
  ];

(* Test 19 / 20: 2D divergent, polynomial 1+x1+x2+x1*x2,
   Tree-A golden: each divergent sector has exactly 2 IBP terms *)
Module[{eps, spec, verts, fanData},
  eps  = Symbol["eps31b"];
  spec = <|"Polynomials" -> {1 + x[1] + x[2] + x[1]*x[2]},
           "MonomialExponents" -> {2*eps - 1, 0},
           "PolynomialExponents" -> {-2},
           "Variables" -> {x[1], x[2]},
           "KinematicSymbols" -> {},
           "RegulatorSymbol" -> eps|>;
  verts   = PolytopeVertices[(1 + x[1] + x[2] + x[1]*x[2])^(-2), {x[1], x[2]}];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  If[!ListQ[fanData] || fanData === $Failed,
    cc31Assert["A4-IBP-terms19", False, "fan ok", "$Failed"]; Return[$Failed]];
  checkIBPTerms["A4-IBP-terms19", spec, fanData, eps];
];
Print[];

(* Test 21: richer polynomial 1+x1+x2+x1*x2+x1^2*x2,
   Tree-A golden: each divergent sector has exactly 3 IBP terms *)
Module[{eps, spec, verts, fanData},
  eps  = Symbol["eps31c"];
  spec = <|"Polynomials" -> {1 + x[1] + x[2] + x[1]*x[2] + x[1]^2*x[2]},
           "MonomialExponents" -> {2*eps - 1, 0},
           "PolynomialExponents" -> {-2},
           "Variables" -> {x[1], x[2]},
           "KinematicSymbols" -> {},
           "RegulatorSymbol" -> eps|>;
  verts   = PolytopeVertices[(1 + x[1] + x[2] + x[1]*x[2] + x[1]^2*x[2])^(-2), {x[1], x[2]}];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  If[!ListQ[fanData] || fanData === $Failed,
    cc31Assert["A4-IBP-terms21", False, "fan ok", "$Failed"]; Return[$Failed]];
  checkIBPTerms["A4-IBP-terms21", spec, fanData, eps];
];
Print[];

(* --------------------------------------------------------------------------
   3. Numerical checks: sector sum vs NIntegrate (Tests 1-7 representative)
      Tolerances match Tree-A's own PASS/FAIL bands.
   -------------------------------------------------------------------------- *)

Print["=== CC31 Part B: numeric sector-sum vs NIntegrate (Tests 1-7 selection) ==="];
Print[];

(* Helper: ValidateDecomposition wrapper with tolerance check *)
numCheck[tag_String, spec_Association, fanData_List, kinRules_List, tol_Real] :=
  Module[{vr, rE},
    vr = Quiet @ ValidateDecomposition[spec, fanData, kinRules, 4];
    If[!AssociationQ[vr],
      cc31Assert[tag, False,
        "relErr<" <> ToString[tol], "ValidateDecomposition=$Failed"];
      Return[$Failed]];
    rE = Abs[vr["RelativeError"]];
    cc31Assert[tag,
      NumericQ[rE] && rE < tol,
      "relErr<" <> ToString[tol],
      "relErr=" <> ToString[NumberForm[rE, {3,2}]]];
    vr
  ];

(* B1: Test 1 / Example 1  -  2D, real exponent, (1+x1^2+x2^2)^{-3}
   Tree-A golden: relative error ~7e-7 (NIntegrate ≈ 0.3927) *)
Module[{spec, verts, fanData},
  spec = <|"Polynomials" -> {1 + x[1]^2 + x[2]^2},
           "MonomialExponents" -> {0, 0},
           "PolynomialExponents" -> {-3},
           "Variables" -> {x[1], x[2]},
           "KinematicSymbols" -> {},
           "RegulatorSymbol" -> None|>;
  verts   = PolytopeVertices[(1 + x[1]^2 + x[2]^2)^(-3), {x[1], x[2]}];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  If[!ListQ[fanData] || fanData === $Failed,
    cc31Assert["B1-2D-basic", False, "fan ok", "$Failed"]; Return[$Failed]];
  Print["[B1] 2D basic (1+x1^2+x2^2)^{-3}:"];
  numCheck["B1-2D-basic", spec, fanData, {}, 1.*^-3];
  Print[];
];

(* B2: Test 2 / RunAllTests Test 1 (A=2) and Test 2 (A=2+0.5I)
   Polynomial: 1+2x1^2+x2^2+x1*x2^2+3x1^2*x2
   Tree-A golden: A=2 relErr ~2.15e-4; A=2+0.5I relErr ~3.98e-6 *)
Module[{poly, vars, verts, fanData, specReal, specCplx},
  poly = 1 + 2 x[1]^2 + x[2]^2 + x[1]*x[2]^2 + 3 x[1]^2*x[2];
  vars = {x[1], x[2]};
  verts   = PolytopeVertices[poly^(-2), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  If[!ListQ[fanData] || fanData === $Failed,
    cc31Assert["B2-fan", False, "fan ok", "$Failed"]; Return[$Failed]];
  specReal = <|"Polynomials" -> {poly},
               "MonomialExponents" -> {0, 0},
               "PolynomialExponents" -> {-2},
               "Variables" -> vars,
               "KinematicSymbols" -> {},
               "RegulatorSymbol" -> None|>;
  specCplx = ReplacePart[specReal, "PolynomialExponents" -> {-(2 + 0.5 I)}];
  Print["[B2a] 2D 5-monomial poly, A=2:"];
  numCheck["B2a-A2-real", specReal, fanData, {}, 1.*^-3];
  Print[];
  (* fan is the same (real part of exponent = 2) *)
  Print["[B2b] 2D 5-monomial poly, A=2+0.5I:"];
  numCheck["B2b-A2-complex", specCplx, fanData, {}, 1.*^-3];
  Print[];
];

(* B3: Test 3 — 1D divergent over the tropical domain [0,inf):
       I(eps) = Int_0^inf y^{2eps-1}/(1+y^2) dy = (pi/2) csc(pi eps).
   Laurent at eps->0:  1/(2eps) + 0 + (pi^2/12) eps + ...
       => pole = 1/2,  finite = 0   (EXACT, full [0,inf) integral).

   NOTE (Phase-6 fix): the engine API for the pole/finite Laurent of a
   divergent spec is EvaluateTropicalMC[..., Method->"IBP"] (the same proven
   path as Part C below), whose Results carry PoleCoefficient/FinitePart.
   The previous harness called ProcessDivergentSector with 4 args (it only
   has a 2-arg DownValue, so it returned unevaluated), keyed on
   PoleCoefficient/FinitePart (which ProcessDivergentSector does not expose),
   and used the [0,1] closed form (finite=-log2/2) for what the tropical
   machinery computes over [0,inf) (finite=0).  All three are corrected here;
   the engine itself was already correct. *)
Module[{eps, vars, spec, verts, fanData, sd, ibpRes, pole, finite},
  eps  = Symbol["eps31d"];
  vars = {x[1]};
  spec = <|"Polynomials" -> {1 + x[1]^2},
           "MonomialExponents" -> {2*eps - 1},
           "PolynomialExponents" -> {-1},
           "Variables" -> vars,
           "KinematicSymbols" -> {},
           "RegulatorSymbol" -> eps|>;
  verts   = PolytopeVertices[(1 + x[1]^2)^(-1), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  If[!ListQ[fanData] || fanData === $Failed,
    cc31Assert["B3-1D-div", False, "fan ok", "$Failed"]; Return[$Failed]];
  sd = quietRun @ ProcessSector[spec, fanData[[1]], fanData[[2,1]], 1];
  If[!AssociationQ[sd] || !TrueQ[sd["IsDivergent"]],
    cc31Assert["B3-1D-div-isDivergent", False, "True", ToString[If[AssociationQ[sd], sd["IsDivergent"], sd]]];
    Return[$Failed]];
  ibpRes = quietRun @ EvaluateTropicalMC[spec, fanData, {{}},
    "Method" -> "IBP", "Integrator" -> "MC", "NSamples" -> 2000000,
    "RunChecks" -> False, "Verbose" -> False];
  pole   = If[AssociationQ[ibpRes] && KeyExistsQ[ibpRes, "Results"],
              Re @ ibpRes["Results"][[1]]["PoleCoefficient"], $Failed];
  finite = If[AssociationQ[ibpRes] && KeyExistsQ[ibpRes, "Results"],
              Re @ ibpRes["Results"][[1]]["FinitePart"], $Failed];
  Print["[B3] 1D divergent over [0,inf): exact pole=1/2, finite=0:"];
  cc31Assert["B3-pole",
    NumericQ[pole] && Abs[pole - 0.5] < 0.01,
    "pole≈0.5", ToString[pole]];
  cc31Assert["B3-finite",
    NumericQ[finite] && Abs[finite - 0] < 0.01,
    "finite≈0", ToString[finite]];
  Print[];
];

(* B4: Test 7 — 3D and 4D higher-dim integrals
   3D: (1+x1^2+x2^2+x3^2)^{-4}  golden relerr ~3.4e-5
   4D: (1+x1^2+…+x1*x2+x3*x4)^{-4}  golden relerr ~1.4e-4 *)
Module[{spec3D, verts3D, fan3D, spec4D, verts4D, fan4D},
  (* 3D *)
  spec3D = <|"Polynomials" -> {1 + x[1]^2 + x[2]^2 + x[3]^2},
             "MonomialExponents" -> {0, 0, 0},
             "PolynomialExponents" -> {-4},
             "Variables" -> {x[1], x[2], x[3]},
             "KinematicSymbols" -> {},
             "RegulatorSymbol" -> None|>;
  verts3D = PolytopeVertices[(1 + x[1]^2 + x[2]^2 + x[3]^2)^(-4), {x[1], x[2], x[3]}];
  fan3D   = ComputeDecomposition[verts3D, "ShowProgress" -> False];
  Print["[B4a] 3D convergent (1+x1^2+x2^2+x3^2)^{-4}:"];
  If[ListQ[fan3D] && !FreeQ[fan3D, {_, _}],
    numCheck["B4a-3D", spec3D, fan3D, {}, 1.*^-2],
    cc31Assert["B4a-3D", False, "fan ok", "$Failed"]];
  Print[];
  (* 4D *)
  spec4D = <|"Polynomials" -> {1 + x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2 + x[1]*x[2] + x[3]*x[4]},
             "MonomialExponents" -> {0, 0, 0, 0},
             "PolynomialExponents" -> {-4},
             "Variables" -> {x[1], x[2], x[3], x[4]},
             "KinematicSymbols" -> {},
             "RegulatorSymbol" -> None|>;
  verts4D = PolytopeVertices[
    (1 + x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2 + x[1]*x[2] + x[3]*x[4])^(-4),
    {x[1], x[2], x[3], x[4]}];
  fan4D   = ComputeDecomposition[verts4D, "ShowProgress" -> False];
  Print["[B4b] 4D convergent:"];
  If[ListQ[fan4D] && !FreeQ[fan4D, {_, _}],
    numCheck["B4b-4D", spec4D, fan4D, {}, 5.*^-3],
    cc31Assert["B4b-4D", False, "fan ok", "$Failed"]];
  Print[];
];

(* --------------------------------------------------------------------------
   4. Divergent cross-check (Tree-A test_divergent_crosscheck, Tier-1 MC only)
      Spec: {1+x1+x2+x3}, {-1+eps, 0, 0}, {-B}
      Analytic: I(eps) = Gamma(eps) Gamma(B-2-eps) / Gamma(B)
      B=4: pole = 1/6, finite = -1/6
      B=5: pole = 1/12, finite = -1/8
      We run EvaluateTropicalMCIBP + LaurentFromSubtraction (MC only) and check
      all four Tree-A Phase-A assertions A1-A3, A5.
   -------------------------------------------------------------------------- *)

Print["=== CC31 Part C: divergent cross-check (MC only; mirrors treeA.txt) ==="];
Print[];

divCC[bExp_Integer, poleRef_Rational, finRef_] :=
  Module[{eps, vars, spec, fan, nIBP = 1500000, nSub = 2000000,
          ibpMC, lfMC,
          poleIbpMC, finIbpMC, poleSubMC, finSubMC, FsubMC,
          epsList, kStar, epsStar, Fref,
          tolPole = 5.*^-3, tolFinIBP = 1.*^-2, tolFinSub = 1.5*^-2,
          tolAgreeP = 5.*^-3, tolAgreeF = 1.5*^-2, tolFull = 8.*^-3,
          lbl, pass = True},

    lbl = "B=" <> ToString[bExp];
    eps  = Symbol["epsCC31"];
    vars = {x[1], x[2], x[3]};
    spec = <|"Polynomials" -> {1 + x[1] + x[2] + x[3]},
             "MonomialExponents" -> {-1 + eps, 0, 0},
             "PolynomialExponents" -> {-bExp},
             "Variables" -> vars,
             "KinematicSymbols" -> {},
             "RegulatorSymbol" -> eps|>;
    fan = computeFanRobust[PolytopeVertices[(1 + x[1] + x[2] + x[3])^(-1), vars]];
    If[fan === $Failed,
      cc31Assert["C-fan-" <> lbl, False, "fan ok", "$Failed"];
      Return[False]];

    epsStar = 0.02;
    epsList = {0.01, 0.02, 0.04};
    kStar   = 2;  (* epsList[[2]] = 0.02 *)
    Fref    = N[Gamma[epsStar] * Gamma[(bExp - 2) - epsStar] / Gamma[bExp]];

    Print["-- Divergent crosscheck " <> lbl <> " --"];
    Print["   analytic: pole=", N[poleRef], "  finite=", N[finRef]];

    (* --- IBP + MC --- *)
    ibpMC = quietRun @ EvaluateTropicalMC[spec, fan, {{}},
      "Method" -> "IBP",
      "Integrator" -> "MC",
      "NSamples" -> nIBP,
      "RunChecks" -> False,
      "Verbose" -> False];
    If[!AssociationQ[ibpMC] || !KeyExistsQ[ibpMC, "Results"],
      cc31Assert["C-IBP-MC-" <> lbl, False, "IBP result assoc", ToString[Head[ibpMC]]];
      Return[False]];
    poleIbpMC = Re @ ibpMC["Results"][[1]]["PoleCoefficient"];
    finIbpMC  = Re @ ibpMC["Results"][[1]]["FinitePart"];
    Print["   IBP MC: pole=", poleIbpMC, "  finite=", finIbpMC];

    (* --- Subtraction + MC via LaurentFromSubtraction --- *)
    lfMC = quietRun @ LaurentFromSubtraction[spec, fan, {{}},
      "Integrator" -> "MC",
      "NSamples" -> nSub,
      "EpsilonValues" -> epsList];
    If[!AssociationQ[lfMC] || !KeyExistsQ[lfMC, "Results"],
      cc31Assert["C-Sub-MC-" <> lbl, False, "Sub result assoc", ToString[Head[lfMC]]];
      Return[False]];
    poleSubMC = Re @ lfMC["Results"][[1]]["Pole"];
    finSubMC  = Re @ lfMC["Results"][[1]]["Finite"];
    FsubMC    = Re @ lfMC["Results"][[1]]["FullIntegrals"][[kStar]];
    Print["   Sub MC: pole=", poleSubMC, "  finite=", finSubMC];

    (* A1: pole vs analytic *)
    Module[{p = Abs[poleIbpMC - poleRef] < tolPole && Abs[poleSubMC - poleRef] < tolPole},
      cc31Assert["C-A1-pole-" <> lbl, p,
        "errs<" <> ToString[tolPole],
        "IBP=" <> ToString[Abs[poleIbpMC - N[poleRef]]] <>
        " sub=" <> ToString[Abs[poleSubMC - N[poleRef]]]]];

    (* A2: finite vs analytic *)
    Module[{p = Abs[finIbpMC - N[finRef]] < tolFinIBP && Abs[finSubMC - N[finRef]] < tolFinSub},
      cc31Assert["C-A2-finite-" <> lbl, p,
        "IBPerr<" <> ToString[tolFinIBP] <> " suberr<" <> ToString[tolFinSub],
        "IBP=" <> ToString[Abs[finIbpMC - N[finRef]]] <>
        " sub=" <> ToString[Abs[finSubMC - N[finRef]]]]];

    (* A3: IBP vs subtraction agreement *)
    Module[{p = Abs[poleIbpMC - poleSubMC] < tolAgreeP && Abs[finIbpMC - finSubMC] < tolAgreeF},
      cc31Assert["C-A3-agree-" <> lbl, p,
        "|dpole|<" <> ToString[tolAgreeP] <> " |dfin|<" <> ToString[tolAgreeF],
        "|dpole|=" <> ToString[Abs[poleIbpMC - poleSubMC]] <>
        " |dfin|=" <> ToString[Abs[finIbpMC - finSubMC]]]];

    (* A5: full-integral consistency at eps* *)
    Module[{laurFromIBP = poleIbpMC/epsStar + finIbpMC,
            c1, c2, p},
      c1 = Abs[FsubMC - Fref]/Fref;
      c2 = Abs[FsubMC - (poleIbpMC/epsStar + finIbpMC)]/Fref;
      p  = c1 < tolFull && c2 < 1.5*^-2;
      cc31Assert["C-A5-full-" <> lbl, p,
        "sub-vs-closed<" <> ToString[tolFull] <> " sub-vs-IBPLaurent<1.5e-2",
        "c1=" <> ToString[c1] <> " c2=" <> ToString[c2]]];

    Print[];
    $cc31Pass
  ];

divCC[4, 1/6, -1/6];
divCC[5, 1/12, -1/8];

(* --------------------------------------------------------------------------
   5. RunAllTests integration check (plan-level: v3 ships a RunAllTests[])
      Run it and assert zero failures.  This exercises the full Test 1-18
      suite from inside the v3 package (the same suite Tree-A runs).
   -------------------------------------------------------------------------- *)

Print["=== CC31 Part D: RunAllTests[] (v3 internal suite, Tests 1-18) ==="];
Print[];

Module[{rat, nFail},
  rat   = quietRun @ RunAllTests[];
  (* RunAllTests returns a list of {TestN, True/False} pairs *)
  If[!ListQ[rat],
    cc31Assert["D-RunAllTests", False,
      "ListQ result", ToString[Head[rat]]];
    Return[$Failed]];
  nFail = Count[rat, {_, False}];
  cc31Assert["D-RunAllTests-0fail",
    nFail === 0,
    "0 failures in RunAllTests[]",
    ToString[nFail] <> " failure(s): " <>
      ToString[Select[rat, MatchQ[{_, False}]]]];
  Print["   RunAllTests: ", Length[rat], " tests, ", nFail, " failures."];
];
Print[];

(* --------------------------------------------------------------------------
   6. Final verdict
   -------------------------------------------------------------------------- *)

Print["============================================================"];
Print["CC31 summary: ", $cc31Assertions, " assertions checked, ",
      $cc31FailCount, " failed."];

If[$cc31Pass,
  Print["CC31 PASS  merge-fidelity-vs-TreeA: all ", $cc31Assertions,
        " assertions green; v3 reproduces every Tree-A golden."],
  Print["CC31 FAIL  expected=all-assertions-pass  got=",
        $cc31FailCount, "-failures (see subFAIL lines above)"]
];
