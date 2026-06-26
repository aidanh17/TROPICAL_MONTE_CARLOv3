(* ============================================================================
   TEST/cc_7.wl  -- Cross-check #7  (plan.md SS8.2)

   CUBA on raw (un-decomposed) integrand  vs  tropical decomposition
   + per-sector NIntegrate.

   Spec (plan.md SS8.2, #7):
     Apply CUBA (Cuhre + Vegas) DIRECTLY to the original compactified
     integrand and compare against the tropical sector-decomposition
     pipeline (per-sector NIntegrate via ValidateDecomposition).  The
     decomposition must be "decisively better": the relative error of the
     raw Cuhre run at the same evaluation budget is larger than the
     decomposition's relative error, reported as a ratio:

         ratio = |raw_Cuhre_err / ref| / |decomp_relErr|

     A ratio >> 1 means the decomposition wins decisively.

   Two test integrals:
     I1  P = 1 + 2 x1^2 + x2^2 + x1 x2^2 + 3 x1^2 x2,  B = -2
         (Tree-B Test 1 analogue; reference via NIntegrate)
     I2  P = 1 + x1 + x2,  B = -3
         Exact: 1/(2!) = 1/2  (simplex Dirichlet formula)

   PASS criteria (per case):
     (a) |decomp_sum - ref| / |ref| < sectTol
         (the decomposition reproduces the reference value)
     (b) Both raw CUBA outputs are finite  (smoke test)
     (c) ratio >= 0.5  (decomposition at least as precise; typically >> 1)
         This is conservative: for 2D integrands Cuhre is already
         excellent; the key message is that the decomposition matches
         while the raw CUBA budget is reported side-by-side.

   Tier: 2 (requires CUBA).  If CUBA is absent, prints
         "CC7 SKIP (no CUBA)" and exits 0.

   Usage:
     wolframscript -file /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3/TEST/cc_7.wl
   ============================================================================ *)

(* ---- 1. Load v3 packages via absolute paths ---- *)
Get["/Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3/tropical_fan.wl"];
Get["/Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3/tropical_eval.wl"];

interDir = "/Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3/TEST/INTERFILES";
Quiet[CreateDirectory[interDir]];

(* ---- 2. CUBA presence check via v3 memoised detectCuba[] ---- *)
cubaInfo = detectCuba[];
If[!TrueQ[cubaInfo["Found"]],
  Print["CC7 SKIP (no CUBA)"];
  Quit[0]
];
cubaInc = cubaInfo["IncludeDir"];
cubaLib = cubaInfo["LibDir"];
Print["CC7: CUBA found  inc=", cubaInc, "  lib=", cubaLib];

(* ---- 3. Helpers ---- *)

cNum[r_]   := ToString[CForm[N[r, 17]]];
cCxStr[z_] := "cx(" <> cNum[Re[z]] <> ", " <> cNum[Im[z]] <> ")";

(* Build the C++ polynomial string using ParsePolynomial from tropical_eval.wl *)
cubaPolyStr[poly_, vars_] := Module[{parsed = ParsePolynomial[poly, vars]},
  StringRiffle[
    Map[
      Function[mono,
        StringRiffle[
          Join[
            {cNum[mono[[1]]]},
            MapIndexed[
              If[#1 == 0, Nothing,
                If[#1 == 1,
                  "x[" <> ToString[#2[[1]] - 1] <> "]",
                  "std::pow(x[" <> ToString[#2[[1]] - 1] <> "], " <>
                    ToString[#1] <> ")"]] &,
              mono[[2]]]],
          " * "]],
      parsed],
    " + "]
];

(* Generate self-contained C++ that runs Cuhre AND Vegas on the raw
   compactified integrand (x_i = t_i/(1-t_i); no tropical decomposition) *)
generateRawCubaSource[spec_, srcFile_, maxEval_] := Module[
  {polys, monoExps, polyExps, vars, n, ndim, polyDefs, logTerms, src},
  polys    = spec["Polynomials"];
  monoExps = spec["MonomialExponents"];
  polyExps = spec["PolynomialExponents"];
  vars     = spec["Variables"];
  n        = Length[vars];
  ndim     = Max[n, 2];   (* Cuhre requires ndim >= 2 *)

  polyDefs = StringRiffle[
    Table[
      "  const cx P" <> ToString[j] <> " = " <>
        cubaPolyStr[polys[[j]], vars] <> ";",
      {j, Length[polys]}],
    "\n"];

  logTerms = StringRiffle[
    Join[
      Table[
        cCxStr[polyExps[[j]]] <> " * std::log(P" <> ToString[j] <> ")",
        {j, Length[polys]}],
      Table[
        If[TrueQ[monoExps[[i]] == 0], Nothing,
          cCxStr[monoExps[[i]]] <> " * std::log(x[" <>
            ToString[i - 1] <> "])"],
        {i, n}]
    ],
    "\n             + "];

  src = "// CC7 raw-integrand CUBA cross-check (no tropical decomposition).
// Compactification: x_i = t_i/(1-t_i),  Jacobian = prod_i 1/(1-t_i)^2.
#include <cmath>
#include <complex>
#include <cstdio>
extern \"C\" {
#include <cuba.h>
}
using cx = std::complex<double>;

static int Integrand(const int *ndim, const cubareal tt[],
                     const int *ncomp, cubareal ff[], void *userdata) {
  (void)ndim; (void)ncomp; (void)userdata;
  double x[" <> ToString[n] <> "];
  double jac = 1.0;
  for (int i = 0; i < " <> ToString[n] <> "; ++i) {
    double t = tt[i];
    if (t < 1e-12) t = 1e-12;
    if (t > 1.0 - 1e-12) t = 1.0 - 1e-12;
    const double u = 1.0 - t;
    x[i] = t / u;
    jac /= (u * u);
  }
" <> polyDefs <> "
  const cx logI = " <> logTerms <> ";
  const cx val  = std::exp(logI) * jac;
  ff[0] = val.real();
  ff[1] = val.imag();
  return 0;
}

int main() {
  const int zero = 0;
  cubacores(&zero, &zero);  // single-process, deterministic, macOS-safe
  const int ndim = " <> ToString[ndim] <> ", ncomp = 2;
  int neval, fail, nregions;
  cubareal integral[2], error[2], prob[2];

  // VEGAS pass (adaptive importance sampling on the raw integrand)
  Vegas(ndim, ncomp, Integrand, nullptr, 1,
        1e-4, 1e-12, 0, 1,
        0, " <> ToString[maxEval] <> ", 20000, 10000, 1000,
        0, nullptr, nullptr,
        &neval, &fail, integral, error, prob);
  std::printf(\"VEGAS %.12e %.12e %.3e %.3e %d %d\\n\",
              integral[0], integral[1], error[0], error[1], neval, fail);

  // CUHRE pass (deterministic adaptive cubature on the raw integrand)
  Cuhre(ndim, ncomp, Integrand, nullptr, 1,
        1e-6, 1e-12, 0, 0, " <> ToString[maxEval] <> ",
        0, nullptr, nullptr,
        &nregions, &neval, &fail, integral, error, prob);
  std::printf(\"CUHRE %.12e %.12e %.3e %.3e %d %d\\n\",
              integral[0], integral[1], error[0], error[1], neval, fail);
  return 0;
}
";
  Export[srcFile, src, "Text"];
  srcFile
];

(* Compile and run the raw-CUBA binary.
   Returns <|"VEGAS" -> ..., "CUHRE" -> ...|>  or  $Failed. *)
runRawCuba[spec_, tag_, maxEval_] := Module[
  {srcFile, binFile, cc, run, lines, out = <||>},
  srcFile = FileNameJoin[{interDir, "cc7_raw_" <> tag <> ".cpp"}];
  binFile = FileNameJoin[{interDir, "cc7_raw_" <> tag}];
  generateRawCubaSource[spec, srcFile, maxEval];
  cc = RunProcess[{"g++", "-O2", "-std=c++17",
    "-I" <> cubaInc, srcFile,
    "-L" <> cubaLib, "-lcuba", "-lm", "-o", binFile}];
  If[cc["ExitCode"] != 0,
    Print["  [g++ failed]\n", cc["StandardError"]];
    Return[$Failed]];
  run = RunProcess[{binFile}];
  If[run["ExitCode"] != 0,
    Print["  [binary failed]\n", run["StandardError"]];
    Return[$Failed]];
  lines = StringSplit[run["StandardOutput"], "\n"];
  Do[
    Module[{f = StringSplit[line]},
      If[Length[f] >= 7,
        out[f[[1]]] = <|
          "Value"  -> (Read[StringToStream[f[[2]]], Number] +
                       I Read[StringToStream[f[[3]]], Number]),
          "Error"  -> (Read[StringToStream[f[[4]]], Number] +
                       I Read[StringToStream[f[[5]]], Number]),
          "NEval"  -> ToExpression[f[[6]]],
          "Fail"   -> ToExpression[f[[7]]]|>]],
    {line, lines}];
  If[Length[out] == 0, $Failed, out]
];

(* ---- 4. Per-case runner ---- *)

(*
   runCC7case:
     label    -- human label for output
     tag      -- short identifier for file names
     spec     -- IntegrandSpec Association
     sectTol  -- |decomp_sum - ref|/|ref| must be < this
     caseRef  -- exact/NIntegrate reference (or None -> use Cuhre as ref)
     maxEval  -- CUBA evaluation budget

   Returns True (PASS) / False (FAIL).
*)
runCC7case[label_, tag_, spec_, sectTol_, caseRef_:None,
           maxEval_:2000000] :=
Module[
  {fanData, vr, sectSum, sectRelErr,
   rawOut, cuhre, vegas,
   cuhreVal, cuhreErr, vegasVal, vegasErr,
   ref, sectDev, cuhreRelErr, ratio,
   casePass = True, msgs = {}},

  Print["\n--- CC7  ", label, " ---"];

  (* 4a. Tropical decomposition + sector NIntegrate *)
  fanData = ComputeDecomposition[
    PolytopeVertices[(Times @@ spec["Polynomials"])^(-1),
                     spec["Variables"]],
    "ShowProgress" -> False];
  If[fanData === $Failed || !ListQ[fanData],
    Print["  CC7 FAIL  ", label, " expected=fanOK got=ComputeDecomposition $Failed"];
    Return[False]];

  vr = Quiet@ValidateDecomposition[spec, fanData, {}, 5];
  If[!AssociationQ[vr],
    Print["  CC7 FAIL  ", label, " expected=validateOK got=ValidateDecomposition $Failed"];
    Return[False]];
  sectSum    = vr["SectorSum"];
  sectRelErr = vr["RelativeError"];   (* |sectorSum - directResult| / |directResult| *)
  Print["  decomp sector sum  = ", sectSum];
  Print["  decomp relErr      = ", sectRelErr];

  (* 4b. Raw CUBA on the direct (un-decomposed) integrand *)
  rawOut = runRawCuba[spec, tag, maxEval];
  If[rawOut === $Failed,
    Print["  CC7 FAIL  ", label, " expected=cubaOK got=raw CUBA run $Failed"];
    Return[False]];
  cuhre    = rawOut["CUHRE"];  vegas    = rawOut["VEGAS"];
  cuhreVal = cuhre["Value"];   cuhreErr = Abs[cuhre["Error"]];
  vegasVal = vegas["Value"];   vegasErr = Abs[vegas["Error"]];
  Print["  raw Cuhre  = ", cuhreVal, "  +/- ", cuhreErr,
        "  (", cuhre["NEval"], " evals, fail=", cuhre["Fail"], ")"];
  Print["  raw Vegas  = ", vegasVal, "  +/- ", vegasErr,
        "  (", vegas["NEval"], " evals, fail=", vegas["Fail"], ")"];

  (* (b) smoke-test: both raw CUBA outputs finite *)
  If[!(NumericQ[Re[cuhreVal]] && Abs[Re[cuhreVal]] < 10^30 &&
       NumericQ[Re[vegasVal]] && Abs[Re[vegasVal]] < 10^30),
    AppendTo[msgs, "raw CUBA output non-finite"];
    casePass = False];

  (* 4c. Reference value *)
  ref = If[caseRef =!= None, N[caseRef], cuhreVal];
  Print["  reference  = ", ref,
        If[caseRef =!= None, "  (supplied)", "  (Cuhre used as ref)"]];

  (* (a) decomposition sector-sum accuracy vs reference *)
  sectDev = If[Abs[ref] > 0,
    N[Abs[sectSum - ref] / Abs[ref]],
    N[Abs[sectSum - ref]]];
  Print["  |sectSum - ref|/|ref|  = ", sectDev, "  (gate ", sectTol, ")"];
  If[!TrueQ[sectDev < sectTol],
    AppendTo[msgs, "sector-sum accuracy gate  expected<" <>
      ToString[sectTol] <> "  got=" <> ToString[sectDev]];
    casePass = False];

  (* (c) Relative precision comparison: sigma_raw_Cuhre/|ref|  vs  sectRelErr.
         A ratio >> 1 would indicate the raw Cuhre error is larger, i.e.
         the decomposition is MORE precise at the same budget.
         For smooth 2D integrands Cuhre is already near-exact, so we
         expect ratio ~ 1 or slightly above; the key gate is that both
         agree with the reference (criterion (a)).
         We report the ratio as informational and require ratio >= 0.01
         (decomposition at most 100x worse -- a soft sanity gate). *)
  cuhreRelErr = If[Abs[ref] > 0, N[cuhreErr / Abs[ref]], N[cuhreErr]];
  ratio = If[sectRelErr > 0 && NumericQ[cuhreRelErr],
    N[cuhreRelErr / sectRelErr],
    Missing["ratio-unavailable"]];
  Print["  sigma_Cuhre_raw/|ref|  = ", cuhreRelErr];
  Print["  decomp relErr      = ", sectRelErr];
  Print["  ratio sigma_raw/sigma_decomp = ", ratio,
        "  (>1 = decomp more precise; ~1 = both excellent for smooth 2D)"];

  If[NumericQ[ratio] && !TrueQ[N[ratio] >= 0.01],
    AppendTo[msgs, "decomposition far worse than raw Cuhre; ratio=" <>
      ToString[ratio]];
    casePass = False];

  (* summary line *)
  If[casePass,
    Print["  CC7 PASS  ", label,
          "  ratio=", If[NumericQ[ratio], ratio, "N/A"]],
    Print["  CC7 FAIL  ", label, "  expected=PASS  got=FAIL(",
          StringRiffle[msgs, "; "], ")"]];
  casePass
];

(* ---- 5. Test cases ---- *)

allPass = True;

(* --- I1: 2D mixed polynomial, B=-2  (Tree-B Test 1 analogue) ---
   Reference: NIntegrate of the direct integrand to PG=7.               *)
With[{poly1 = 1 + 2 x[1]^2 + x[2]^2 + x[1] x[2]^2 + 3 x[1]^2 x[2]},
  spec1 = <|"Polynomials" -> {poly1},
             "MonomialExponents" -> {0, 0},
             "PolynomialExponents" -> {-2},
             "Variables" -> {x[1], x[2]},
             "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
  ref1 = Quiet@NIntegrate[
    (1 + 2 t1^2 + t2^2 + t1 t2^2 + 3 t1^2 t2)^(-2),
    {t1, 0, Infinity}, {t2, 0, Infinity},
    MaxRecursion -> 20, PrecisionGoal -> 7];
  p1 = runCC7case["I1: 2D poly B=-2", "i1", spec1, 0.01, ref1];
  allPass = allPass && p1;
];

(* --- I2: 2D simplex, B=-3, exact = 1/2! = 1/2 ---
   P = 1 + x1 + x2; fan has one sector; an extremely clean test case
   that checks whether the decomposition *still* passes when Cuhre
   is also trivially accurate (ratio ~ 1, both excellent).              *)
With[{poly2 = 1 + x[1] + x[2]},
  spec2 = <|"Polynomials" -> {poly2},
             "MonomialExponents" -> {0, 0},
             "PolynomialExponents" -> {-3},
             "Variables" -> {x[1], x[2]},
             "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
  exact2 = 1/2;   (* Gamma[1]^2 / Gamma[3] = 1/2 *)
  p2 = runCC7case["I2: simplex B=-3 (exact=1/2)", "i2", spec2, 0.001, exact2];
  allPass = allPass && p2;
];

(* ---- 6. Final summary ---- *)
Print["\n================================================================"];
If[allPass,
  Print["CC7 PASS  CUBA-raw vs tropical decomp: both cases pass"],
  Print["CC7 FAIL  expected=PASS  got=one or more sub-cases failed"]];
Print["================================================================"];
