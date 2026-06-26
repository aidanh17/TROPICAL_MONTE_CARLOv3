(* ============================================================================
   TEST/cc_5.wl  —  v3 Cross-check #5 (plan.md §8.2)

   CUBA Cuhre vs tropical sector sum (direct integrand, no decomposition).

   For one representative convergent integral the test:
     (a) runs CUBA Cuhre deterministic cubature on the ORIGINAL integrand
         directly on [0,inf)^n (compactified via t = x/(1+x));
     (b) computes the tropical-decomposition sector sum via
         ValidateDecomposition (NIntegrate over each flattened sector);
     (c) asserts  |sectorSum - cuhreVal| / |cuhreVal| < 0.01  (1 %).

   Integral used:
       P = 1 + 2 x[1]^2 + x[2]^2 + x[1] x[2]^2 + 3 x[1]^2 x[2],  B = -2
   This is Test 1 / Case A=2 of Tree-B's validation suite (treeB.txt
   baseline) — a clean 2D convergent case with no degeneracies.

   Tier: 2 (needs CUBA).
   If CUBA is absent the script prints "CC5 SKIP (no CUBA)" and exits 0.

   Run:
       cd /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3
       wolframscript -file TEST/cc_5.wl

   Reference (plan.md §8, #5):
       "vs CUBA Cuhre on the direct (un-decomposed) integrand —
        sector sum within 1–5%; MC within 5 sigma"
   ============================================================================ *)

(* --------------------------------------------------------------------------
   0.  Load v3 packages using absolute paths (plan.md §10: use absolute paths).
   -------------------------------------------------------------------------- *)

$v3Root = "/Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3";
Get[FileNameJoin[{$v3Root, "tropical_fan.wl"}]];
Get[FileNameJoin[{$v3Root, "tropical_eval.wl"}]];

Print["CC5: tropical_fan.wl + tropical_eval.wl loaded."];

(* --------------------------------------------------------------------------
   1.  CUBA availability gate (Tier-2: skip cleanly if absent).
   -------------------------------------------------------------------------- *)

$cubaInfo = detectCuba[];

If[!TrueQ[$cubaInfo["Found"]],
  Print["CC5 SKIP (no CUBA)"];
  Exit[0]
];

$cubaInc = $cubaInfo["IncludeDir"];
$cubaLib = $cubaInfo["LibDir"];

Print["CC5: CUBA found — include: ", $cubaInc, "  lib: ", $cubaLib];

(* --------------------------------------------------------------------------
   2.  Integral specification.
   -------------------------------------------------------------------------- *)

$spec = <|
  "Polynomials"         -> {1 + 2 x[1]^2 + x[2]^2 + x[1] x[2]^2 + 3 x[1]^2 x[2]},
  "MonomialExponents"   -> {0, 0},
  "PolynomialExponents" -> {-2},
  "Variables"           -> {x[1], x[2]},
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> None
|>;

(* --------------------------------------------------------------------------
   3.  Tropical fan + sector decomposition (symbolic — no numerics here).
   -------------------------------------------------------------------------- *)

$verts   = PolytopeVertices[(Times @@ $spec["Polynomials"])^(-1),
                             $spec["Variables"]];
$fanData = ComputeDecomposition[$verts, "ShowProgress" -> False];
Print["CC5: fan computed — ", Length[$fanData[[2]]], " sectors."];

(* --------------------------------------------------------------------------
   4.  Tropical sector sum via ValidateDecomposition + NIntegrate.
   -------------------------------------------------------------------------- *)

$vr = Quiet @ ValidateDecomposition[$spec, $fanData, {}, 3];

If[!AssociationQ[$vr] || !NumericQ[$vr["SectorSum"]],
  Print["CC5 FAIL expected=<cuhre> got=$Failed (ValidateDecomposition failed)"];
  Exit[1]
];

$sectorSum = $vr["SectorSum"];
Print["CC5: sector sum = ", $sectorSum];

(* --------------------------------------------------------------------------
   5.  Generate + compile + run a CUBA Cuhre program for the DIRECT integrand.

   Compactification: x_i = t_i / (1 - t_i), dx_i = dt_i / (1-t_i)^2.
   Integrand on [0,1]^n (clamped away from boundaries):
       prod_i (1/(1-t_i))^2  *  prod_j P_j(x)^{B_j}  *  prod_i x_i^{A_i}
   Implemented in std::complex arithmetic so the framework generalises;
   two Cuba components (Re, Im) with only the real part used for real B.

   No kinematic symbols; numeric literals only — no MmaToC dependency.
   -------------------------------------------------------------------------- *)

(* Helpers matching cuba_common.wl pattern from Tree B. *)
cNum[r_] := ToString[CForm[N[r, 17]]];
cubaComplexStr[z_] := "cx(" <> cNum[Re[z]] <> ", " <> cNum[Im[z]] <> ")";

(* Emit the C++ literal for one monomial coeff * prod x[i]^exp  (C-0-indexed) *)
cubaMonoStr[coeff_, expVec_] :=
  StringRiffle[
    Join[
      {cubaComplexStr[coeff]},
      MapIndexed[
        Function[{e, idx},
          Which[
            e == 0, Nothing,
            e == 1, "x[" <> ToString[idx[[1]] - 1] <> "]",
            True,   "std::pow(x[" <> ToString[idx[[1]] - 1] <> "], " <>
                    ToString[e] <> ")"
          ]
        ],
        expVec
      ]
    ],
    " * "
  ];

(* Emit the C++ expression for a polynomial (sum of monomials). *)
cubaPolyStr[poly_, vars_] := Module[{parsed},
  parsed = ParsePolynomial[poly, vars];
  StringRiffle[
    Table[cubaMonoStr[mono[[1]], mono[[2]]], {mono, parsed}],
    " + "
  ]
];

(* Build the full C++ source for Cuhre on the direct integrand. *)
cc5GenerateCubaSource[spec_, srcFile_,
    cubaIncDir_, cubaLibDir_,
    maxEval_: 2000000, cuhreEpsRel_: 10^-6] :=
Module[{polys, monoExps, polyExps, vars, n, ndim,
        polyDefs, logTerms, cuhreCall, src},

  polys    = spec["Polynomials"];
  monoExps = spec["MonomialExponents"];
  polyExps = spec["PolynomialExponents"];
  vars     = spec["Variables"];
  n        = Length[vars];
  ndim     = Max[n, 2];    (* Cuhre requires ndim >= 2 *)

  (* P_j definitions *)
  polyDefs = StringRiffle[
    Table[
      "    const cx P" <> ToString[j] <> " = " <>
      cubaPolyStr[polys[[j]], vars] <> ";",
      {j, Length[polys]}
    ], "\n"
  ];

  (* log-integrand:  sum_j B_j log P_j  +  sum_i A_i log x_i  (skip A=0) *)
  logTerms = StringRiffle[
    Join[
      Table[
        cubaComplexStr[polyExps[[j]]] <> " * std::log(P" <>
        ToString[j] <> ")",
        {j, Length[polys]}
      ],
      Table[
        If[TrueQ[monoExps[[i]] == 0], Nothing,
          cubaComplexStr[monoExps[[i]]] <> " * std::log(x[" <>
          ToString[i - 1] <> "])"],
        {i, n}
      ]
    ],
    "\n             + "
  ];

  cuhreCall = StringJoin[
    "  Cuhre(ndim, ncomp, Integrand, nullptr, 1,\n",
    "         ", cNum[cuhreEpsRel], ", 0.0, 0,\n",
    "         0, ", ToString[maxEval], ", 0,\n",
    "         nullptr, nullptr,\n",
    "         &nregions, &neval, &fail, val, err, prob);\n"
  ];

  src = "// cc_5.wl — auto-generated CUBA Cuhre cross-check (direct integrand).
// Compactification x_i = t_i/(1-t_i); no tropical decomposition.
#include <cmath>
#include <complex>
#include <cstdio>
#include <algorithm>
#include <cuba.h>

typedef std::complex<double> cx;

static const int ndim  = " <> ToString[ndim] <> ";
static const int ncomp = 2;   // Re, Im

static int Integrand(const int *ndim_, const double t[],
                     const int *ncomp_, double out[], void *) {
  // compactify: x_i = t_i/(1-t_i), clamp t away from 0 and 1
  double x_raw[" <> ToString[ndim] <> "];
  double jacobian = 1.0;
  for (int i = 0; i < " <> ToString[n] <> "; ++i) {
    double ti = std::max(1e-12, std::min(1.0-1e-12, t[i]));
    double xi = ti / (1.0 - ti);
    x_raw[i]  = xi;
    jacobian  *= 1.0 / ((1.0 - ti) * (1.0 - ti));
  }
  const double *x = x_raw;  (void)x;
" <> polyDefs <> "
  cx logval = " <> logTerms <> ";
  cx val = cx(jacobian, 0.0) * std::exp(logval);
  out[0] = val.real();
  out[1] = val.imag();
  return 0;
}

int main() {
  int    nregions, neval, fail;
  double val[2], err[2], prob[2];
" <> cuhreCall <>
"  printf(\"CUHRE  Re=%.15e  Im=%.15e  errRe=%.6e  errIm=%.6e  neval=%d  fail=%d\\n\",
         val[0], val[1], err[0], err[1], neval, fail);
  return fail;
}
";
  Export[srcFile, src, "Text"];
  srcFile
];

(* Write source, compile, run, parse. *)
$cc5SrcDir = FileNameJoin[{$v3Root, "TEST", "INTERFILES"}];
If[!DirectoryQ[$cc5SrcDir], CreateDirectory[$cc5SrcDir]];

$cc5Src = FileNameJoin[{$cc5SrcDir, "cc5_direct.cpp"}];
$cc5Bin = FileNameJoin[{$cc5SrcDir, "cc5_direct"}];

cc5GenerateCubaSource[$spec, $cc5Src, $cubaInc, $cubaLib,
  2000000, 10^-6];

(* Compile. Try with -fopenmp first, fall back without it. *)
$compileCmd = {
  "g++", "-std=c++17", "-O2", "-fopenmp",
  "-I", $cubaInc, "-L", $cubaLib,
  $cc5Src, "-o", $cc5Bin, "-lcuba", "-lm"
};

$compileResult = RunProcess[$compileCmd];

If[$compileResult["ExitCode"] =!= 0,
  (* retry without -fopenmp *)
  $compileCmd2 = Select[$compileCmd, (# =!= "-fopenmp") &];
  $compileResult = RunProcess[$compileCmd2];
];

If[$compileResult["ExitCode"] =!= 0,
  Print["CC5 FAIL expected=<cuhre> got=$Failed (compile error: ",
        $compileResult["StandardError"], ")"];
  Exit[1]
];

(* Run — set DYLD_LIBRARY_PATH so the dynamic lib is found on macOS. *)
$runResult = RunProcess[
  {$cc5Bin},
  ProcessEnvironment -> Association[
    Normal[GetEnvironment[]],
    "DYLD_LIBRARY_PATH" -> $cubaLib,
    "LD_LIBRARY_PATH"   -> $cubaLib
  ]
];

(* CUBA Cuhre exits non-zero when fail=1 (max evals reached but result is
   still returned).  We tolerate fail=1 — what matters is whether the
   printed value agrees with the sector sum.  A genuine crash produces no
   output, which the parse step below catches. *)
If[$runResult["ExitCode"] > 1,
  Print["CC5 FAIL expected=<cuhre> got=$Failed (runtime error: ",
        $runResult["StandardError"], ")"];
  Exit[1]
];

$cuhreOut = $runResult["StandardOutput"];
Print["CC5: Cuhre output: ", StringTrim[$cuhreOut]];

(* Parse "Re=<v>  Im=<v>  ..." — values are in scientific notation e.g. 3.61e-01 *)
$sciNumPat = RegularExpression["[+-]?[0-9]+(?:\\.[0-9]+)?(?:[eE][+-]?[0-9]+)?"];
$reMatch  = StringCases[$cuhreOut, "Re="  ~~ n:$sciNumPat :> n];
$errMatch = StringCases[$cuhreOut, "errRe=" ~~ n:$sciNumPat :> n];

If[$reMatch === {} || $errMatch === {},
  Print["CC5 FAIL expected=<cuhre> got=$Failed (could not parse Cuhre output)"];
  Exit[1]
];

(* ToExpression mis-parses "3.61e-01" — use ImportString["CSV"] instead *)
parseSciNum[s_String] := First @ Flatten @ ImportString[s, "CSV"];
$cuhreVal = parseSciNum[$reMatch[[1]]];
$cuhreErr = parseSciNum[$errMatch[[1]]];

Print["CC5: Cuhre value = ", $cuhreVal, "  +/- ", $cuhreErr];

(* --------------------------------------------------------------------------
   6.  Gate: |sectorSum - cuhreVal| / |cuhreVal| < 0.01  (plan §8 #5: 1–5%)
   -------------------------------------------------------------------------- *)

$relDev = Abs[($sectorSum - $cuhreVal) / $cuhreVal];

Print["CC5: sector sum = ", $sectorSum];
Print["CC5: |rel dev|  = ", $relDev, "  (gate 0.01)"];

If[TrueQ[$relDev < 0.01],
  Print["CC5 PASS  sectorSum=", $sectorSum,
        "  cuhre=", $cuhreVal, "  reldev=", $relDev],
  Print["CC5 FAIL expected=", $cuhreVal, " got=", $sectorSum,
        "  reldev=", $relDev, "  (gate 0.01 exceeded)"];
  Exit[1]
];
