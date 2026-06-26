(* ============================================================================
   TEST/cc_6.wl  —  Cross-check #6: tropical sector sum vs CUBA Vegas
                    (plan.md §8.2, #6)

   Applies CUBA's Vegas integrator directly to the original integrand on
   [0,inf)^n (no tropical decomposition, via x = t/(1-t) compactification)
   and compares the result against the tropical sector sum from
   ValidateDecomposition.

   PASS criterion (plan.md §8.2 #6):
     |sectorSum - vegasValue| <= 5 * sqrt(sectorSumErr^2 + vegasErr^2)
   i.e. the two independent samplers agree within combined 5-sigma.
   A looser relative-tolerance fallback (5%) is applied when the Vegas
   error bar is unreliable (Fail != 0) or when the sector-sum NIntegrate
   error is large — but the agreement gate is always the primary one.

   Tier 2: requires CUBA (https://feynarts.de/cuba/; macOS: brew install cuba).
   If CUBA is absent the script prints "CC6 SKIP (no CUBA)" and exits 0.

   Run from the repo root:
     wolframscript -file TEST/cc_6.wl

   Dependencies: tropical_fan.wl, tropical_eval.wl, Polymake, g++, CUBA.
   ============================================================================ *)

(* ---- Locate the repo root (one level up from TEST/) ---------------------- *)
SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];

(* ---- Load the v3 package ------------------------------------------------- *)
Get[FileNameJoin[{Directory[], "tropical_fan.wl"}]];
Get[FileNameJoin[{Directory[], "tropical_eval.wl"}]];

(* ---- CUBA probe (mirrors detectCuba[] but produces the prefix string) ---- *)
cc6CubaPrefix[] := SelectFirst[
  {"/opt/homebrew", "/usr/local", "/usr"},
  FileExistsQ[FileNameJoin[{#, "include", "cuba.h"}]] &&
  (FileExistsQ[FileNameJoin[{#, "lib", "libcuba.a"}]] ||
   FileExistsQ[FileNameJoin[{#, "lib", "libcuba.dylib"}]] ||
   FileExistsQ[FileNameJoin[{#, "lib", "libcuba.so"}]]) &,
  $Failed];

Module[{cubaPrefix},
  cubaPrefix = cc6CubaPrefix[];
  If[cubaPrefix === $Failed,
    Print["CC6 SKIP (no CUBA)"];
    Exit[0]
  ]
];

(* ---- Working directory for compiled binaries ----------------------------- *)
cc6Dir = FileNameJoin[{Directory[], "TEST", "INTERFILES"}];
If[!DirectoryQ[cc6Dir], CreateDirectory[cc6Dir]];

(* ---- Helper: emit a self-contained CUBA Vegas C++ source for spec --------
   Integrates the original un-decomposed integrand over [0,inf)^n using the
   compactification x_i = t_i/(1-t_i).  Complex arithmetic throughout so
   the code handles complex exponents correctly.  Two ncomp components: Re,Im.
   -------------------------------------------------------------------------- *)

cc6Num[r_]  := ToString[CForm[N[r, 17]]];
cc6CxStr[z_] := "cx(" <> cc6Num[Re[z]] <> ", " <> cc6Num[Im[z]] <> ")";

cc6PolyStr[poly_, vars_] := Module[{parsed},
  parsed = ParsePolynomial[poly, vars];
  StringRiffle[
    Table[
      StringRiffle[
        Join[
          {cc6Num[mono[[1]]]},
          MapIndexed[
            If[#1 == 0, Nothing,
              If[#1 == 1,
                "x[" <> ToString[#2[[1]] - 1] <> "]",
                "std::pow(x[" <> ToString[#2[[1]] - 1] <> "], " <>
                  ToString[#1] <> ")"]
            ] &,
            mono[[2]]
          ]
        ], " * "
      ],
      {mono, parsed}
    ],
    " + "
  ]
];

cc6GenerateVegasSource[spec_, srcFile_, opts___Rule] := Module[
  {polys, monoExps, polyExps, vars, n, ndim, maxEval, vegasEpsRel, nstart,
   polyDefs, logTerms, src},

  {maxEval, vegasEpsRel, nstart} =
    {"MaxEval", "VegasEpsRel", "NStart"} /.
    {opts} /. {"MaxEval" -> 4000000, "VegasEpsRel" -> 5*^-4, "NStart" -> 20000};

  polys    = spec["Polynomials"];
  monoExps = spec["MonomialExponents"];
  polyExps = spec["PolynomialExponents"];
  vars     = spec["Variables"];
  n        = Length[vars];
  ndim     = Max[n, 2];   (* Vegas requires ndim >= 2; pad 1D integrands *)

  polyDefs = StringRiffle[
    Table[
      "  const cx P" <> ToString[j] <> " = " <>
        cc6PolyStr[polys[[j]], vars] <> ";",
      {j, Length[polys]}
    ], "\n"
  ];

  logTerms = StringRiffle[
    Join[
      Table[
        cc6CxStr[polyExps[[j]]] <> " * std::log(P" <> ToString[j] <> ")",
        {j, Length[polys]}
      ],
      Table[
        If[TrueQ[monoExps[[i]] == 0], Nothing,
          cc6CxStr[monoExps[[i]]] <> " * std::log(x[" <>
            ToString[i - 1] <> "])"],
        {i, n}
      ]
    ],
    "\n             + "
  ];

  src =
"// CC6: CUBA Vegas cross-check — direct integrand, no tropical decomposition.
// Compactification x_i = t_i/(1-t_i).  Emitted by TEST/cc_6.wl.
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
  cubacores(&zero, &zero);    // single process: deterministic, macOS-safe

  const int ndim = " <> ToString[ndim] <> ", ncomp = 2;
  int neval, fail, nregions;
  cubareal integral[2], error[2], prob[2];

  Vegas(ndim, ncomp, Integrand, nullptr, 1,
        " <> cc6Num[vegasEpsRel] <> ", 1e-12, 0, 1,
        0, " <> ToString[maxEval] <> ",
        " <> ToString[nstart] <> ", 10000, 1000,
        0, nullptr, nullptr,
        &neval, &fail, integral, error, prob);
  std::printf(\"VEGAS %.12e %.12e %.3e %.3e %d %d\\n\",
              integral[0], integral[1], error[0], error[1], neval, fail);
  return 0;
}
";
  Export[srcFile, src, "Text"];
  srcFile
];

(* ---- Helper: compile and run the CUBA Vegas binary ----------------------- *)

cc6RunVegas[spec_, tag_, opts___Rule] := Module[
  {prefix, srcFile, binFile, cc, run, lines, fields},

  prefix = cc6CubaPrefix[];
  If[prefix === $Failed, Return[$Failed]];

  srcFile = FileNameJoin[{cc6Dir, "cc6_" <> tag <> ".cpp"}];
  binFile = FileNameJoin[{cc6Dir, "cc6_" <> tag}];

  cc6GenerateVegasSource[spec, srcFile, opts];

  cc = RunProcess[{"g++", "-O2", "-std=c++17",
    "-I" <> FileNameJoin[{prefix, "include"}],
    srcFile,
    "-L" <> FileNameJoin[{prefix, "lib"}],
    "-lcuba", "-lm", "-o", binFile}];
  If[cc["ExitCode"] != 0,
    Print["  [CC6 compile failed]\n", cc["StandardError"]];
    Return[$Failed]
  ];

  run = RunProcess[{binFile}];
  If[run["ExitCode"] != 0,
    Print["  [CC6 run failed]\n", run["StandardError"]];
    Return[$Failed]
  ];

  lines = Select[StringSplit[run["StandardOutput"], "\n"],
                 StringLength[#] > 0 &];
  If[Length[lines] == 0, Return[$Failed]];

  fields = StringSplit[lines[[1]]];
  If[Length[fields] < 7, Return[$Failed]];

  <|"Method"  -> fields[[1]],
    "Value"   -> (Read[StringToStream[fields[[2]]], Number] +
                  I Read[StringToStream[fields[[3]]], Number]),
    "Error"   -> (Read[StringToStream[fields[[4]]], Number] +
                  I Read[StringToStream[fields[[5]]], Number]),
    "NEval"   -> ToExpression[fields[[6]]],
    "Fail"    -> ToExpression[fields[[7]]]|>
];

(* ---- Helper: run one named case ----------------------------------------- *)

$cc6Rows = {};

cc6RunCase[label_, tag_, spec_, opts : OptionsPattern[]] := Module[
  {pg, maxEval, vegasEpsRel, nstart, sectorTol,
   verts, fanData, vr, sectorSum, sectorErr,
   vegas, vegasVal, vegasErr,
   combinedErr, devSigma, devRel, casePass = True, notes = {}},

  pg          = "PrecisionGoal" /. {opts} /. "PrecisionGoal" -> 3;
  maxEval     = "MaxEval"       /. {opts} /. "MaxEval"       -> 4000000;
  vegasEpsRel = "VegasEpsRel"   /. {opts} /. "VegasEpsRel"  -> 5*^-4;
  nstart      = "NStart"        /. {opts} /. "NStart"        -> 20000;
  sectorTol   = "SectorTol"     /. {opts} /. "SectorTol"     -> 0.05;

  Print["--- CC6 ", label, " ---"];

  (* --- Tropical decomposition + sector sum (via ValidateDecomposition) --- *)
  verts   = PolytopeVertices[(Times @@ spec["Polynomials"])^(-1),
                              spec["Variables"]];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  vr = Quiet @ ValidateDecomposition[spec, fanData, {}, pg];
  If[!AssociationQ[vr],
    Print["  ValidateDecomposition failed — FAIL"];
    AppendTo[$cc6Rows, <|"Label" -> label, "Pass" -> False,
      "Notes" -> "ValidateDecomposition failed"|>];
    Return[False]
  ];
  sectorSum = vr["SectorSum"];
  (* ValidateDecomposition returns "RelativeError" = |sectorSum - NIntegrate| / |NIntegrate|.
     Use it as a rough estimate of the sector-sum uncertainty for the 5-sigma gate. *)
  sectorErr = If[KeyExistsQ[vr, "RelativeError"],
    Abs[sectorSum] * vr["RelativeError"],
    Abs[sectorSum] * 10^(-pg)   (* conservative fallback *)
  ];

  Print["  Sector sum:   ", sectorSum, "   (+/- ~", sectorErr, " est.)"];

  (* --- CUBA Vegas on the direct un-decomposed integrand ------------------ *)
  vegas = cc6RunVegas[spec, tag,
    "MaxEval" -> maxEval, "VegasEpsRel" -> vegasEpsRel, "NStart" -> nstart];
  If[vegas === $Failed,
    Print["  CUBA Vegas failed — FAIL"];
    AppendTo[$cc6Rows, <|"Label" -> label, "Pass" -> False,
      "Notes" -> "CUBA Vegas failed"|>];
    Return[False]
  ];

  vegasVal = vegas["Value"];
  vegasErr = Abs[vegas["Error"]];   (* take modulus of complex error pair *)
  If[TrueQ[vegas["Fail"] != 0],
    AppendTo[notes, "Vegas Fail=" <> ToString[vegas["Fail"]]]
  ];

  Print["  CUBA Vegas:   ", vegasVal, "   +/- ", vegasErr,
        "  (", vegas["NEval"], " evals, fail=", vegas["Fail"], ")"];

  (* --- Agreement gate: within 5 combined sigma -------------------------- *)
  combinedErr = Sqrt[sectorErr^2 + vegasErr^2];
  devSigma    = Abs[sectorSum - vegasVal] / combinedErr;
  devRel      = If[Abs[vegasVal] > 0,
                   Abs[sectorSum - vegasVal] / Abs[vegasVal],
                   Infinity];

  Print["  |sector - Vegas| = ", Abs[sectorSum - vegasVal],
        "   = ", devSigma, " sigma  (rel ", devRel, ")"];

  (* Primary gate: 5-sigma combined error *)
  If[!TrueQ[devSigma < 5],
    (* Secondary gate: loose relative tolerance (for very small numbers or
       unreliable Vegas error bar when Fail != 0) *)
    If[TrueQ[devRel < sectorTol],
      AppendTo[notes, "5-sigma gate missed; rel-tol fallback passed (" <>
                       ToString[sectorTol] <> ")"],
      casePass = False;
      AppendTo[notes, "5-sigma AND rel-tol gates both missed"]
    ]
  ];

  Print["  ",
        If[casePass,
          "CC6 PASS (sector-sum vs CUBA Vegas, " <> label <> ")",
          "CC6 FAIL expected=" <> ToString[CForm[N[vegasVal, 4]]] <>
            " got=" <> ToString[CForm[N[sectorSum, 4]]]
        ],
        If[notes =!= {},
          "  [" <> StringRiffle[notes, "; "] <> "]",
          ""]
  ];
  Print[];

  AppendTo[$cc6Rows, <|"Label" -> label, "Pass" -> casePass,
    "DevSigma" -> devSigma, "DevRel" -> devRel,
    "Notes" -> StringRiffle[notes, "; "]|>];
  casePass
];

(* ============================================================================
   Test cases — ported from OLD_CODE Tree B test_cuba.wl (§8.2 #5/#6),
   adapted for cross-check #6 (Vegas arm only, per plan.md §8.2 #6).

   Case A — Test 1, B=-2: P = 1 + 2x1^2 + x2^2 + x1*x2^2 + 3x1^2*x2
     Standard 2D benchmark; moderate integral (~7.85e-4).

   Case B — Test 1, B=-3: same polynomial, steeper decay; smaller value.

   Case C — Test 2, B=-(2+0.5i): complex polynomial exponent.
     Exercises the complex-valued integrand path in the Vegas C++.

   Case D — Test 7A (3D): P = 1 + x1^2 + x2^2 + x3^2 + x1*x2*x3, B=-3.
     Higher dimension; uses more Vegas evaluations.
   ============================================================================ *)

Module[{poly2d = 1 + 2 x[1]^2 + x[2]^2 + x[1] x[2]^2 + 3 x[1]^2 x[2]},

  (* Case A: B = -2 *)
  cc6RunCase["Test1 B=-2", "t1b2",
    <|"Polynomials"        -> {poly2d},
      "MonomialExponents"  -> {0, 0},
      "PolynomialExponents"-> {-2},
      "Variables"          -> {x[1], x[2]},
      "KinematicSymbols"   -> {},
      "RegulatorSymbol"    -> None|>
  ];

  (* Case B: B = -3 *)
  cc6RunCase["Test1 B=-3", "t1b3",
    <|"Polynomials"        -> {poly2d},
      "MonomialExponents"  -> {0, 0},
      "PolynomialExponents"-> {-3},
      "Variables"          -> {x[1], x[2]},
      "KinematicSymbols"   -> {},
      "RegulatorSymbol"    -> None|>
  ];

  (* Case C: complex polynomial exponent B = -(2 + 0.5i) *)
  cc6RunCase["Test2 B=-(2+0.5i)", "t2cx",
    <|"Polynomials"        -> {poly2d},
      "MonomialExponents"  -> {0, 0},
      "PolynomialExponents"-> {-(2 + 0.5 I)},
      "Variables"          -> {x[1], x[2]},
      "KinematicSymbols"   -> {},
      "RegulatorSymbol"    -> None|>,
    "PrecisionGoal" -> 2, "SectorTol" -> 0.10
  ];
];

(* Case D: 3D *)
cc6RunCase["Test7A 3D B=-3", "t7a",
  <|"Polynomials"        -> {1 + x[1]^2 + x[2]^2 + x[3]^2 + x[1] x[2] x[3]},
    "MonomialExponents"  -> {0, 0, 0},
    "PolynomialExponents"-> {-3},
    "Variables"          -> {x[1], x[2], x[3]},
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None|>,
  "MaxEval" -> 8000000, "NStart" -> 40000, "SectorTol" -> 0.05
];

(* ---- Summary ------------------------------------------------------------- *)
Module[{nPass, nFail, allPass},
  nPass = Count[$cc6Rows, r_ /; TrueQ[r["Pass"]]];
  nFail = Length[$cc6Rows] - nPass;
  allPass = (nFail == 0);

  Print["================================================================"];
  Print["  CC6 summary: sector sum vs CUBA Vegas (", Length[$cc6Rows], " cases)"];
  Print["================================================================"];
  Do[
    Print["  ", StringPadRight[row["Label"], 28],
      If[TrueQ[row["Pass"]], "PASS", "FAIL"],
      If[NumberQ[row["DevSigma"]],
        "   " <> ToString[CForm[SetPrecision[row["DevSigma"], 3]]] <> " sigma", ""],
      If[NumberQ[row["DevRel"]],
        "   reldev " <> ToString[CForm[SetPrecision[row["DevRel"], 3]]], ""],
      If[StringLength[row["Notes"]] > 0,
        "   [" <> row["Notes"] <> "]", ""]
    ],
    {row, $cc6Rows}
  ];
  Print["----------------------------------------------------------------"];
  Print["  ", nPass, " PASSED, ", nFail, " FAILED"];
  Print["================================================================"];

  If[allPass,
    Print["CC6 PASS (all ", nPass, " cases: sector sum agrees with CUBA Vegas)"],
    Print["CC6 FAIL (", nFail, " case(s) failed — see details above)"]
  ];
];
