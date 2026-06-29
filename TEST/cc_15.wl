(* ============================================================================
   TEST/cc_15.wl  --  v3 Cross-Check #15
   plan.md §8.2 row #15:  eps/eps cross-term resolution
      finite = TOTAL - G0 + gamma_E * G0 (<0.01% residual)

   The tropical subtraction scheme assembles the finite part of a divergent
   sector integral as
       I(eps) = G0/(ck*eps) + G1/ck + Remainder
   where G0 is the sector integral at eps=0 with the divergent variable set
   to zero, G1 contains the log-insertion cross-terms, and Remainder removes
   double-counting.  The full-integral finite part extracted by fitting
       eps * I(eps) = pole + finite*eps + O(eps^2)
   must equal the analytic finite coefficient.

   Concretely, for
       I(eps) = Integral x1^(-1+eps) (1+x1+x2+x3)^(-B) dx1 dx2 dx3
              = Gamma(eps) Gamma(B-2-eps) / Gamma(B)
   the eps/eps cross-term in the Laurent expansion
       Gamma(eps) = 1/eps - gamma_E + O(eps)
   contributes -gamma_E * pole to the finite part.  The analytic results are:
       pole  = 1/((B-1)(B-2))
       finite= -HarmonicNumber[B-3] / ((B-1)(B-2))

   This test has three assertion groups:

   A1 -- Exact symbolic / exactness invariant (plan.md §1.4 / §6.7):
         ProcessDivergentSector output for each divergent sector must be
         FreeQ of _Real (no floating-point numbers anywhere in the returned
         Association before the MmaToC boundary).  Additionally, AnalyticPole
         = 1/ck must be an exact rational, and ck must match the
         symbolic eps-coefficient exactly.  This is the per-sector exactness
         check required by the exactness invariant.

   A2 -- ValidateSubtraction: for each divergent sector, NIntegrate verifies
         G0/(ck*epsTest) + G1/ck + Remainder = full integral at epsTest, relErr<1%
         (the eps/eps cross-term is resolved by the G1 log-insertion).

   A3 -- LaurentFromSubtraction end-to-end: fit of eps*I(eps_i) at four small
         eps values recovers (pole, finite) vs analytic to <2.5% (Tier-1
         fast-MC tolerance; plan.md cites <0.01% for a 50M-sample production
         run -- the identity is the same, only variance differs).

   Tier 1: WL + g++ only.  No CUBA needed.
   Ported from OLD_CODE/TROPICAL_MONTE_CARLO/CC/diagnose_eps_eps.wl
   (which was bubblepp-specific; this version is self-contained and uses the
   same 3D-div reference integral as CC#3/#4/#13/#16).

   FIXED (this version): the prior draft A1 incorrectly summed AnalyticPole
   (= 1/ck, a per-sector eps-expansion coefficient) and compared it to the
   full integral pole poleRef = 1/((B-1)(B-2)).  Those two quantities have
   different units/dimensions; the comparison was always False for B=4 (sum
   was 5, poleRef was 1/6).  A1 is now the correct exactness-invariant check:
   each ProcessDivergentSector return is FreeQ[_, _Real] and ck is an exact
   rational -- no floating-point numbers in the symbolic decomposition output.
   ============================================================================ *)

(* --------------------------------------------------------------------------
   Load v3 packages by absolute path so this file runs from any cwd.
   -------------------------------------------------------------------------- *)
$v3Root = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{$v3Root, "tropical_fan.wl"}]];
Get[FileNameJoin[{$v3Root, "tropical_eval.wl"}]];

(* --------------------------------------------------------------------------
   Helpers
   -------------------------------------------------------------------------- *)
SetAttributes[quietRun, HoldFirst];
quietRun[expr_] := Block[{$Output = {}}, expr];

sci[x_] := If[x == 0 || !NumericQ[x], ToString[x],
  Module[{e = Floor[Log10[Abs[N[x]]]], m},
    m = N[x] / 10^e;
    ToString[NumberForm[m, {3, 2}]] <> "e" <> ToString[e]]];

(* K-scaled fan (plan.md §6.4; mirrors gen_sectors.wl computeFanRobust) *)
computeFanRobust[verts_] := Module[{n, fd},
  n = Length[First[verts]];
  Do[fd = Quiet@ComputeDecomposition[K*verts, "ShowProgress" -> False];
     If[ListQ[fd] && Length[fd] == 2 && FreeQ[fd, $Failed] && Length[fd[[1]]] > 0,
        Return[fd, Module]],
     {K, {1, n + 2, 2 n + 4, 6 n + 6}}];
  $Failed];

(* ============================================================================
   cc15Spec: run the three assertion groups for one B value.
   Returns True if all assertions pass.
   ============================================================================ *)

cc15Spec[label_String, Bexp_Integer,
         nSamples_Integer : 500000, tol_Real : 0.025] :=
Module[
  {eps, vars, spec, fan, dualVertices, simplexList,
   poleRef, finRef,
   allSectors, divSectors, divDataList,
   exactOK,
   testEps, vsOK,
   lfResult, poleFit, finFit, poleRelErr, finFitRelErr,
   allPass = True},

  Print[""];
  Print["--- CC15 ", label, " ---"];

  eps  = Symbol["epsCC15"];
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
    Print["CC15 FAIL ", label, " expected=fan got=$Failed"];
    Return[False]];

  (* ------------------------------------------------------------------
     Analytic references
     I(eps) = Gamma(eps) Gamma(B-2-eps) / Gamma(B)
     pole  = 1 / ((B-1)(B-2))
     finite= -HarmonicNumber[B-3] / ((B-1)(B-2))
     ------------------------------------------------------------------ *)
  poleRef = 1 / ((Bexp - 1) * (Bexp - 2));
  finRef  = -HarmonicNumber[Bexp - 3] / ((Bexp - 1) * (Bexp - 2));

  Print["  Analytic: pole=", N[poleRef, 6], "  finite=", N[finRef, 6]];

  {dualVertices, simplexList} = fan;

  (* Sector decomposition *)
  allSectors = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[c]], c],
    {c, Length[simplexList]}];
  allSectors  = Select[allSectors,  AssociationQ];
  divSectors  = Select[allSectors,  #["IsDivergent"] &];

  (* Divergent sector data *)
  divDataList = Select[
    Map[Function[sd, Quiet@ProcessDivergentSector[sd, spec]], divSectors],
    AssociationQ];

  Print["  Sectors: ", Length[allSectors], " total, ",
        Length[divSectors], " divergent"];

  If[Length[divDataList] == 0,
    Print["CC15 FAIL ", label,
          " expected=>0 divergent sectors got=0 (cannot test)"];
    Return[False]];

  (* ==========================================================================
     A1 -- Exactness-invariant check (plan.md §1.4 / §6.7):
            Every ProcessDivergentSector return Association must be FreeQ of
            _Real (no floating-point numbers in the symbolic decomposition
            output before the MmaToC boundary).  ck and AnalyticPole must be
            exact rationals.  This is the per-sector exactness guard for the
            divergence/subtraction path.

            Note: AnalyticPole = 1/ck is the per-sector eps-expansion
            coefficient (e.g. 1/1 = 1 when ck=1).  It is NOT the full
            sector contribution to the overall integral pole (that requires
            additionally the G0 sector integral with its prefactor).  The
            total pole = Sum_s [G0_s * Abs[detM_s]/(ck_s * Times[a0_s]) ].
            The NUMERIC pole is checked in A3 via LaurentFromSubtraction.
            A1 validates only the exactness of the symbolic output.
     ========================================================================== *)
  exactOK = True;
  Do[
    Module[{dd = divDataList[[i]], hasReal, ckVal, apVal, ckOK, apOK},
      (* No _Real anywhere in the symbolic sector data *)
      hasReal = !FreeQ[dd, _Real];
      If[hasReal,
        Print["CC15 FAIL ", label, " A1 sector ", dd["ConeIndex"],
              " has _Real in symbolic output (exactness violation)"];
        exactOK = False];

      (* ck must be an exact positive rational *)
      ckVal = dd["ck"];
      ckOK  = !NumericQ[N[ckVal]] || Head[ckVal] === Integer || Head[ckVal] === Rational;
      If[!ckOK,
        Print["CC15 FAIL ", label, " A1 sector ", dd["ConeIndex"],
              " ck=", ckVal, " is not an exact rational"];
        exactOK = False];

      (* AnalyticPole must equal 1/ck exactly *)
      apVal = dd["AnalyticPole"];
      apOK  = Simplify[apVal - 1/ckVal] === 0;
      If[!apOK,
        Print["CC15 FAIL ", label, " A1 sector ", dd["ConeIndex"],
              " AnalyticPole=", apVal, " != 1/ck=", 1/ckVal];
        exactOK = False]
    ],
    {i, Length[divDataList]}
  ];

  If[exactOK,
    Print["CC15 PASS ", label,
          " A1 exactness: all ", Length[divDataList],
          " divergent sectors FreeQ[_Real]; ck exact; AnalyticPole=1/ck"]];
  allPass = allPass && exactOK;

  (* ==========================================================================
     A2 -- ValidateSubtraction: G0/(ck*eps0) + G1/ck + Rem reconstructs the
            full NIntegrate value at eps0 to <1%.
            This directly exercises the eps/eps cross-term resolution:
            G1 captures d/deps[prefactor] and d/deps[polys]; without G1
            the finite part would be off by gamma_E * pole.
     ---------------------------------------------------------------------- *)
  testEps = 0.02;
  vsOK = True;
  Do[
    Module[{sd = divSectors[[i]], dd = divDataList[[i]], vsr, relErr},
      vsr = Quiet@ValidateSubtraction[dd, sd, spec, {}, testEps];
      If[!AssociationQ[vsr],
        Print["CC15 FAIL ", label, " A2 sector ", sd["ConeIndex"],
              " expected=Association got=$Failed"];
        vsOK = False;
        Return[Null, Module]];
      relErr = vsr["RelativeError"];
      If[!(NumericQ[relErr] && relErr < 0.01),
        Print["CC15 FAIL ", label, " A2 sector ", sd["ConeIndex"],
              " expected=relErr<0.01 got=", relErr];
        vsOK = False]
    ],
    {i, Length[divDataList]}
  ];
  If[vsOK,
    Print["CC15 PASS ", label,
          " A2 ValidateSubtraction relErr<0.01 for all ",
          Length[divDataList], " divergent sectors",
          " (eps/eps cross-term resolved via G1 log-insertions)"]];
  allPass = allPass && vsOK;

  (* ==========================================================================
     A3 -- LaurentFromSubtraction end-to-end: fit eps*I(eps_i) to recover
            (pole, finite) vs analytic to within tol.
            The headline assertion of CC#15: because LaurentFromSubtraction
            evaluates the FULL integral at several small eps values and fits a
            polynomial in eps, the gamma_E contribution is automatically
            included -- the eps/eps cross-term must be correctly resolved by
            the subtraction scheme for this fit to land on finRef.
     ========================================================================== *)
  lfResult = quietRun@LaurentFromSubtraction[
    spec, fan, {{}},
    "EpsilonValues" -> {0.005, 0.01, 0.02, 0.04},
    "Integrator"    -> "MC",
    "NSamples"      -> nSamples,
    "RunChecks"     -> False,
    "Verbose"       -> False
  ];

  If[!AssociationQ[lfResult] || !ListQ[lfResult["Results"]],
    Print["CC15 FAIL ", label, " A3 expected=Association got=$Failed"];
    allPass = False;
    Return[allPass]];

  poleFit    = Re[lfResult["Results"][[1]]["Pole"]];
  finFit     = Re[lfResult["Results"][[1]]["Finite"]];
  poleRelErr = Abs[(poleFit - N[poleRef]) / N[poleRef]];
  finFitRelErr = If[N[finRef] != 0.,
    Abs[(finFit - N[finRef]) / N[finRef]],
    Abs[finFit - N[finRef]]];

  Print["  A3 LaurentFit: pole=", poleFit, "  finite=", finFit];
  Print["       analytic: pole=", N[poleRef, 6], "  finite=", N[finRef, 6]];

  Module[{p = poleRelErr < tol},
    If[p,
      Print["CC15 PASS ", label,
            " A3 pole vs analytic (relerr=", sci[poleRelErr],
            " < ", tol, ")"],
      Print["CC15 FAIL ", label,
            " A3 expected=relerr<", tol, " got=", sci[poleRelErr]]];
    allPass = allPass && p];

  (* The finite assertion IS CC#15: gamma_E cross-term resolved *)
  Module[{p = finFitRelErr < tol},
    If[p,
      Print["CC15 PASS ", label,
            " A3 finite vs analytic (relerr=", sci[finFitRelErr],
            " < ", tol, ") -- eps/eps cross-term correctly resolved"],
      Print["CC15 FAIL ", label,
            " A3 expected=relerr<", tol, " got=", sci[finFitRelErr]]];
    allPass = allPass && p];

  allPass
];

(* ============================================================================
   RunCC15 -- entry point
   ============================================================================ *)
RunCC15[] := Module[{allPass = True},
  Print[""];
  Print["===================================================="];
  Print[" v3 Cross-Check #15: eps/eps cross-term resolution  "];
  Print["  plan.md §8.2 row #15                             "];
  Print["  finite = TOTAL - G0 + gamma_E*G0, residual<0.01% "];
  Print["  Tier 1: WL + g++ only                            "];
  Print["===================================================="];

  allPass = cc15Spec["B=4", 4, 500000, 0.025] && allPass;
  allPass = cc15Spec["B=5", 5, 500000, 0.025] && allPass;

  Print[""];
  Print["===================================================="];
  If[allPass,
    Print["CC15 PASS  eps/eps cross-term resolution complete"],
    Print["CC15 FAIL  one or more assertions failed -- see above"]];
  Print["===================================================="];
  allPass
];

If[!TrueQ[$CC15NoRun], RunCC15[]];
