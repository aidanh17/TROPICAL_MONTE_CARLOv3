(* ============================================================================
   TEST/cc_16.wl  ---  v3 Cross-Check #16
   plan.md §8.2 row #16:  Subtracted-remainder convergence as eps->0.
   PASS criterion: remainder R(eps) -> 0 with documented O(eps) rate.

   Tier 1 (WL + g++ only; no CUBA required).

   Mathematical basis
   ------------------
   For a single-divergence spec the tropical subtraction writes

       F(eps) = pole/eps + finite + c_1*eps + c_2*eps^2 + ...

   Define the "subtracted remainder"

       R(eps) := F(eps) - pole/eps - finite

   Then R(eps) = O(eps) as eps -> 0 (linear rate; slope c_1 is finite).

   For the canonical 3-D spec used here

       Polynomials = {1 + x[1] + x[2] + x[3]},
       MonomialExponents = {-1+eps, 0, 0},
       PolynomialExponents = {-4},

   the EXACT closed form is

       F(eps) = Gamma(eps) * Gamma(2-eps) / Gamma(4)
              = 1/(6 eps) - 1/6 + (pi^2/18 - EulerGamma/3 - HarmonicNumber[2]/3)*eps + O(eps^2)

   so  pole_exact = 1/6,   finite_exact = -1/6.

   Test protocol
   -------------
   Layer A (exact -- documents the O(eps) rate analytically):
     A1. Compute R_exact(eps) = F_exact(eps) - pole_exact/eps - finite_exact
         at a grid of eps values.  Assert |R_exact(eps)| = O(eps):
           |R_exact(eps)| < tol * |pole_exact| * eps  for each eps.
     A2. Ratio test on R_exact: R_exact(2*eps)/R_exact(eps) ~ 2 (linear rate).

   Layer B (numeric -- confirms the MC driver respects the subtraction identity):
     B1. Run EvaluateTropicalMC (subtraction mode) at several eps values.
         Fit the Laurent via eps*F_MC(eps) = pole_fit + finite_fit*eps + ...
         Assert |pole_fit - pole_exact| < 1% and |finite_fit - finite_exact| < 3%.
     B2. Assert F_MC(eps_i) agrees with F_exact(eps_i) to < 2% (relerr).

   The O(eps) convergence *rate* is documented by Layer A (exact, no MC noise);
   Layer B validates the driver is computing the same function.

   Ported / adapted from OLD_CODE/TROPICAL_MONTE_CARLO/EXAMPLES/
   test_divergent_crosscheck.wl and the LaurentFromSubtraction logic.
   ============================================================================ *)

(* --------------------------------------------------------------------------
   Load v3 packages by absolute path so this file runs from any working dir.
   -------------------------------------------------------------------------- *)
$v3Root = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{$v3Root, "tropical_fan.wl"}]];
Get[FileNameJoin[{$v3Root, "tropical_eval.wl"}]];

(* --------------------------------------------------------------------------
   Helpers
   -------------------------------------------------------------------------- *)
SetAttributes[quietRun, HoldFirst];
quietRun[expr_] := Block[{$Output = {}}, expr];

sci[x_] := If[!NumericQ[x] || x === 0, "0",
  Module[{e = Floor[Log10[Abs[N[x]]]], m},
    m = N[x] / 10^e;
    ToString[NumberForm[m, {3, 2}]] <> "e" <> ToString[e]]];

(* K-scaled fan robust helper -- mirrors gen_sectors.wl computeFanRobust (plan §6.4) *)
computeFanRobust[verts_] := Module[{n, fd},
  n = Length[First[verts]];
  Do[fd = Quiet@ComputeDecomposition[K * verts, "ShowProgress" -> False];
     If[ListQ[fd] && Length[fd] == 2 && FreeQ[fd, $Failed] && Length[fd[[1]]] > 0,
        Return[fd, Module]], {K, {n + 2, 2 n + 4, 6 n + 6}}];
  $Failed];

(* --------------------------------------------------------------------------
   Exact closed form for the 3-D spec
       I(eps) = Gamma(eps)*Gamma(2-eps)/Gamma(4)
   -------------------------------------------------------------------------- *)
Fexact[e_] := N[Gamma[e] * Gamma[2 - e] / Gamma[4]];

poleExact = N[1/6];
finExact  = N[-1/6];

(* ============================================================================
   LAYER A:  Exact O(eps) rate documentation (no MC, no fan needed)
   ========================================================================== *)
Print["================================================================"];
Print["CC16 Layer A: Exact remainder R(eps) = F_exact(eps) - pole/eps - finite"];
Print["================================================================"];

epsA = {0.01, 0.02, 0.04, 0.08, 0.16, 0.32};

RexactVals = Table[Fexact[ev] - poleExact/ev - finExact, {ev, epsA}];

Print["CC16-A: eps grid  = ", epsA];
Print["CC16-A: R_exact   = ", Map[sci, RexactVals]];
Print["CC16-A: R_exact/eps = ", MapThread[sci[#1/#2] &, {RexactVals, epsA}]];

(* A1: magnitude bound  |R_exact(eps)| < 5 * |poleExact| * eps *)
tolA1 = 5.0;
a1Res = Table[
  Module[{ok, r, b},
    r  = Abs[RexactVals[[i]]];
    b  = tolA1 * Abs[poleExact] * epsA[[i]];
    ok = r < b;
    Print["  A1 eps=", epsA[[i]], "  |R|=", sci[r], "  bound=", sci[b], "  ok=", ok];
    ok],
  {i, Length[epsA]}];
a1Pass = And @@ a1Res;

(* A2: ratio test -- R_exact(eps_{i+1})/R_exact(eps_i) ~ eps_{i+1}/eps_i *)
tolA2 = 0.15;
Print["\nCC16-A: Ratio test (R_exact(2 eps)/R_exact(eps) ~ 2), tol=", tolA2];
a2Res = Table[
  Module[{rRatio, eRatio, reldev, ok},
    rRatio  = If[RexactVals[[i]] != 0, RexactVals[[i+1]] / RexactVals[[i]], Indeterminate];
    eRatio  = epsA[[i+1]] / epsA[[i]];
    reldev  = If[NumericQ[rRatio], Abs[rRatio - eRatio] / Abs[eRatio], 999.];
    ok      = reldev < tolA2;
    Print["  A2 R(", epsA[[i+1]], ")/R(", epsA[[i]], ") = ", sci[rRatio],
          "  eps-ratio=", sci[eRatio], "  reldev=", sci[reldev], "  ok=", ok];
    ok],
  {i, Length[epsA] - 1}];
a2Pass = And @@ a2Res;

Print["\nCC16 Layer A summary: A1=", If[a1Pass, "PASS", "FAIL"],
      "  A2=", If[a2Pass, "PASS", "FAIL"]];

(* ============================================================================
   LAYER B:  Numeric MC driver validation
   ========================================================================== *)
Print["\n================================================================"];
Print["CC16 Layer B: MC driver F_MC(eps) vs F_exact(eps)"];
Print["================================================================"];

(* Build fan *)
Print["CC16-B: building fan..."];
vars  = {x[1], x[2], x[3]};
eps   = Symbol["epsCC16"];
spec  = <|"Polynomials"         -> {1 + x[1] + x[2] + x[3]},
           "MonomialExponents"   -> {-1 + eps, 0, 0},
           "PolynomialExponents" -> {-4},
           "Variables"           -> vars,
           "KinematicSymbols"    -> {},
           "RegulatorSymbol"     -> eps|>;

verts  = PolytopeVertices[(1 + x[1] + x[2] + x[3])^(-1), vars];
fanRaw = quietRun@computeFanRobust[verts];
If[fanRaw === $Failed,
  Print["CC16 FAIL expected=fan_ok got=$Failed"]; Exit[1]];

(* Evaluate F_MC at three eps values -- keep them moderate so MC noise < |R| *)
epsB  = {0.02, 0.05, 0.10};
nSamp = 800000;

Print["CC16-B: evaluating F_MC at eps = ", epsB, "  NSamples=", nSamp, "..."];
FmcVals = Table[
  Module[{r},
    r = quietRun@EvaluateTropicalMC[spec, fanRaw, {{}},
          "EpsilonValue"  -> ev,
          "Method"        -> "None",
          "Integrator"    -> "MonteCarlo",
          "NSamples"      -> nSamp,
          "RunChecks"     -> False,
          "Verbose"       -> False];
    If[AssociationQ[r] && KeyExistsQ[r, "Results"],
      Re[r["Results"][[1]]["Re"]],
      Missing["EvalFailed"]]],
  {ev, epsB}];

If[MemberQ[FmcVals, _Missing],
  Print["CC16 FAIL expected=all_F_MC_numeric got=", FmcVals]; Exit[1]];

Print["CC16-B: F_exact(eps) = ", Map[Fexact, epsB]];
Print["CC16-B: F_MC(eps)    = ", FmcVals];

(* B1: relative error F_MC vs F_exact *)
tolB1 = 0.02;
b1Res = Table[
  Module[{fe, fm, relerr, ok},
    fe     = Fexact[epsB[[i]]];
    fm     = FmcVals[[i]];
    relerr = Abs[(fm - fe) / fe];
    ok     = relerr < tolB1;
    Print["  B1 eps=", epsB[[i]], "  F_exact=", sci[fe], "  F_MC=", sci[fm],
          "  relerr=", sci[relerr], "  ok=", ok];
    ok],
  {i, Length[epsB]}];
b1Pass = And @@ b1Res;

(* B2: Laurent fit from MC data: pole_fit and finite_fit vs analytic *)
gVals  = epsB * FmcVals;
deg    = Min[Length[epsB] - 1, 2];
design = Table[ev^p, {ev, epsB}, {p, 0, deg}];
coeffs = LeastSquares[design, gVals];
poleFit = coeffs[[1]];
finFit  = If[Length[coeffs] >= 2, coeffs[[2]], Missing[]];

tolB2pole = 0.01;
tolB2fin  = 0.03;
b2pPass = Abs[poleFit - poleExact] < tolB2pole * Abs[poleExact];
b2fPass = If[NumericQ[finFit],
  Abs[finFit - finExact] < tolB2fin * Abs[finExact],
  False];
b2Pass  = b2pPass && b2fPass;

Print["\n  B2 Laurent fit from MC:  pole_fit=", sci[poleFit],
      "  (exact=", sci[poleExact], ")  relerr=", sci[Abs[poleFit - poleExact]/Abs[poleExact]],
      "  ok=", b2pPass];
Print["  B2 Laurent fit from MC:  fin_fit=", sci[finFit],
      "  (exact=", sci[finExact], ")  relerr=", sci[If[NumericQ[finFit],
        Abs[finFit - finExact]/Abs[finExact], Indeterminate]],
      "  ok=", b2fPass];

Print["\nCC16 Layer B summary: B1=", If[b1Pass, "PASS", "FAIL"],
      "  B2=", If[b2Pass, "PASS", "FAIL"]];

(* ============================================================================
   Final verdict
   ========================================================================== *)
allPass = a1Pass && a2Pass && b1Pass && b2Pass;

Print["\n================================================================"];
Print["CC16 A1 (exact |R| magnitude bound)  : ", If[a1Pass, "PASS", "FAIL"]];
Print["CC16 A2 (exact linear-rate ratio)    : ", If[a2Pass, "PASS", "FAIL"]];
Print["CC16 B1 (F_MC vs F_exact relerr<2%)  : ", If[b1Pass, "PASS", "FAIL"]];
Print["CC16 B2 (Laurent fit pole/finite)     : ", If[b2Pass, "PASS", "FAIL"]];
Print["================================================================"];

If[allPass,
  Print["CC16 PASS subtracted-remainder converges O(eps); ",
        "rate documented analytically; MC driver consistent with closed form"],
  Module[{failing},
    failing = StringRiffle[
      Select[{"A1", "A2", "B1", "B2"},
             Function[tag, Not[Switch[tag,
               "A1", a1Pass, "A2", a2Pass, "B1", b1Pass, "B2", b2Pass]]]],
      ","];
    Print["CC16 FAIL expected=O(eps)_convergence got=failures_in_", failing]]];
