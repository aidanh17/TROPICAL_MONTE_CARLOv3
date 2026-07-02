(* ============================================================================
   TEST/cc_56.wl  —  CHARACTERIZATION test (fable_plan step 2, following cc_55)

   Pins the CURRENT behavior of low-level public API + string-emission
   behavior of tropical_eval.wl's Module 3 (C++ codegen) primitives, ahead of
   a modularity refactor.  These are "golden master" tests: they assert the
   behavior that exists TODAY (before any refactor), so a later structural
   split can be proven behavior-preserving.  They are NOT a claim that the
   behavior is "correct" in an external sense (e.g. cppBadTokens returning a
   bare List, or normalizeIntegrator silently accepting non-strings, are
   pinned as-is) — they lock down the status quo byte-for-byte / value-for-
   value.  Every expected string/value below was captured by RUNNING the
   code (not guessed) via a throwaway probe script (scratch, not committed).

   Covers, one PART per seam:
     A. MmaToC[expr, paramMap]           (PUBLIC) — one assertion per clause
        of the internal mmaToCInternal dispatcher.
     B. ParsePolynomial[poly, vars]      (PUBLIC).
     C. TropicalEval`Private`normalizeIntegrator[s]        (PRIVATE seam).
     D. TropicalEval`Private`resolveVegasSizing[...]       (PRIVATE seam).
     E. TropicalEval`Private`emitDomainIndicatorCpp[...],
        TropicalEval`Private`dropDivVarFromDomain[...],
        TropicalEval`Private`liftedDomainBooleWL[...]      (PRIVATE seam).
     F. TropicalEval`Private`cppBadTokens[code]             (PRIVATE seam).
     G. TropicalEval`Private`imagMonoInfo / realifyMonoA,
        TropicalEval`Private`imagPolyInfo / realifyPolyB    (PRIVATE seam).
     H. TropicalEval`Private`emitLogTail[...]               (PRIVATE seam).

   These PRIVATE names are the seams the fable_plan refactor has committed
   to keep stable (same name, same signature, same return value) even though
   they live in `Private` — any later restructuring of tropical_eval.wl must
   keep them at TropicalEval`Private`<name> with byte-identical behavior, or
   this test must be consciously updated (not silently broken).

   Tier 0: pure WL, no Polymake fan, no g++, fully deterministic (no
   SeedRandom / no float-tolerance compares — every assertion below is exact
   string or exact-value equality, safe to run in the same process as other
   Tier-0 checks).
   ============================================================================ *)

$pkgRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{$pkgRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$pkgRoot, "tropical_eval.wl"}]];
Print["CC56: packages loaded."];
Print[];

$cc56Pass = True; $cc56Fail = {};
cc56Assert[label_String, test_, exp_, got_] :=
  If[TrueQ[test],
    Print["CC56 PASS  [", label, "]  expected=", exp, "  got=", got],
    Print["CC56 FAIL  [", label, "]  expected=", exp, "  got=", got];
    $cc56Pass = False; AppendTo[$cc56Fail, label]];

(* Local aliases for the private seams under test (style per cc_52/cc_53/cc_54). *)
normalizeIntegrator   = TropicalEval`Private`normalizeIntegrator;
resolveVegasSizing    = TropicalEval`Private`resolveVegasSizing;
emitDomainIndicatorCpp = TropicalEval`Private`emitDomainIndicatorCpp;
dropDivVarFromDomain  = TropicalEval`Private`dropDivVarFromDomain;
liftedDomainBooleWL   = TropicalEval`Private`liftedDomainBooleWL;
cppBadTokensFn        = TropicalEval`Private`cppBadTokens;
imagMonoInfo          = TropicalEval`Private`imagMonoInfo;
realifyMonoA          = TropicalEval`Private`realifyMonoA;
imagPolyInfo          = TropicalEval`Private`imagPolyInfo;
realifyPolyB          = TropicalEval`Private`realifyPolyB;
emitLogTail           = TropicalEval`Private`emitLogTail;

(* ============================================================================
   PART A — MmaToC[expr, paramMap]  (PUBLIC).  One assertion per clause of
   mmaToCInternal (tropical_eval.wl lines ~2099-2193), in source order.
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (A) MmaToC characterization"];

cc56Assert["A1 Integer -> \"N.0\"",
  MmaToC[3] === "3.0", "3.0", MmaToC[3]];

cc56Assert["A2 negative Integer -> \"-N.0\"",
  MmaToC[-3] === "-3.0", "-3.0", MmaToC[-3]];

cc56Assert["A3 Rational -1/2 -> \"(-1.0/2.0)\"",
  MmaToC[-1/2] === "(-1.0/2.0)", "(-1.0/2.0)", MmaToC[-1/2]];

cc56Assert["A4 Real 1.5 -> \"1.5\"",
  MmaToC[1.5] === "1.5", "1.5", MmaToC[1.5]];

cc56Assert["A5 Complex 2+3I -> \"cx(2.0, 3.0)\"",
  MmaToC[2 + 3 I] === "cx(2.0, 3.0)", "cx(2.0, 3.0)", MmaToC[2 + 3 I]];

cc56Assert["A6 pure-imaginary Complex 3I -> \"cx(0.0, 3.0)\"",
  MmaToC[3 I] === "cx(0.0, 3.0)", "cx(0.0, 3.0)", MmaToC[3 I]];

cc56Assert["A7 plain Symbol, no paramMap -> literal name",
  MmaToC[b] === "b", "b", MmaToC[b]];

cc56Assert["A8 plain Symbol, WITH paramMap -> mapped string verbatim",
  MmaToC[b, <|b -> "params[0]"|>] === "params[0]",
  "params[0]", MmaToC[b, <|b -> "params[0]"|>]];

cc56Assert["A9 indexed symbol x[2], no paramMap -> \"x[2]\"",
  MmaToC[x[2]] === "x[2]", "x[2]", MmaToC[x[2]]];

cc56Assert["A10 indexed symbol x[2], WITH paramMap -> mapped array form",
  MmaToC[x[2], <|x[2] -> "kin[2]"|>] === "kin[2]",
  "kin[2]", MmaToC[x[2], <|x[2] -> "kin[2]"|>]];

cc56Assert["A11 Power[b,-1] -> \"(1.0/b)\"",
  MmaToC[Power[b, -1]] === "(1.0/b)", "(1.0/b)", MmaToC[Power[b, -1]]];

cc56Assert["A12 Sqrt[b] -> \"std::sqrt(b)\"",
  MmaToC[Sqrt[b]] === "std::sqrt(b)", "std::sqrt(b)", MmaToC[Sqrt[b]]];

cc56Assert["A13 1/Sqrt[b] -> \"(1.0/std::sqrt(b))\"",
  MmaToC[1/Sqrt[b]] === "(1.0/std::sqrt(b))",
  "(1.0/std::sqrt(b))", MmaToC[1/Sqrt[b]]];

cc56Assert["A14 Power[E,x] -> \"std::exp(x[1])\"",
  MmaToC[Power[E, x[1]]] === "std::exp(x[1])",
  "std::exp(x[1])", MmaToC[Power[E, x[1]]]];

cc56Assert["A15 Power[b,4] (integer power) -> \"std::pow(b, 4.0)\"",
  MmaToC[Power[b, 4]] === "std::pow(b, 4.0)",
  "std::pow(b, 4.0)", MmaToC[Power[b, 4]]];

cc56Assert["A16 Power[b,c] (general power) -> \"std::pow(b, c)\"",
  MmaToC[Power[b, c]] === "std::pow(b, c)",
  "std::pow(b, c)", MmaToC[Power[b, c]]];

cc56Assert["A17 Log[b] -> \"std::log(b)\"",
  MmaToC[Log[b]] === "std::log(b)", "std::log(b)", MmaToC[Log[b]]];

cc56Assert["A18 Exp[b] -> \"std::exp(b)\"",
  MmaToC[Exp[b]] === "std::exp(b)", "std::exp(b)", MmaToC[Exp[b]]];

cc56Assert["A19 Abs[b] -> \"std::abs(b)\"",
  MmaToC[Abs[b]] === "std::abs(b)", "std::abs(b)", MmaToC[Abs[b]]];

cc56Assert["A20 Re[b] -> \"(b).real()\"",
  MmaToC[Re[b]] === "(b).real()", "(b).real()", MmaToC[Re[b]]];

cc56Assert["A21 Im[b] -> \"(b).imag()\"",
  MmaToC[Im[b]] === "(b).imag()", "(b).imag()", MmaToC[Im[b]]];

cc56Assert["A22 Plus (3 args) -> \"(a + b + c)\" (\" + \" spacing)",
  MmaToC[a + b + c] === "(a + b + c)", "(a + b + c)", MmaToC[a + b + c]];

cc56Assert["A23 Times (3 args) -> \"(a * b * c)\" (\" * \" spacing)",
  MmaToC[a*b*c] === "(a * b * c)", "(a * b * c)", MmaToC[a*b*c]];

cc56Assert["A24 Times[-1,x[1]] negation -> \"(-x[1])\"",
  MmaToC[Times[-1, x[1]]] === "(-x[1])", "(-x[1])", MmaToC[Times[-1, x[1]]]];

cc56Assert["A25 fallback clause: Sin[x[1]] -> CForm-derived \"Sin(x(1))\"",
  MmaToC[Sin[x[1]]] === "Sin(x(1))", "Sin(x(1))", MmaToC[Sin[x[1]]]];

(* ============================================================================
   PART B — ParsePolynomial[poly, vars]  (PUBLIC).
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (B) ParsePolynomial characterization"];

cc56Assert["B1 simple 2-var poly (incl. constant term)",
  ParsePolynomial[3 + 2*x[1] + 5*x[1]*x[2]^2, {x[1], x[2]}] ===
    {{3, {0, 0}}, {2, {1, 0}}, {5, {1, 2}}},
  {{3, {0, 0}}, {2, {1, 0}}, {5, {1, 2}}},
  ParsePolynomial[3 + 2*x[1] + 5*x[1]*x[2]^2, {x[1], x[2]}]];

cc56Assert["B2 poly with symbolic (kinematic) coefficients",
  ParsePolynomial[m^2*x[1] + s*x[2]^2, {x[1], x[2]}] ===
    {{m^2, {1, 0}}, {s, {0, 2}}},
  {{m^2, {1, 0}}, {s, {0, 2}}},
  ParsePolynomial[m^2*x[1] + s*x[2]^2, {x[1], x[2]}]];

cc56Assert["B3 single monomial",
  ParsePolynomial[7*x[1]^3*x[2], {x[1], x[2]}] === {{7, {3, 1}}},
  {{7, {3, 1}}}, ParsePolynomial[7*x[1]^3*x[2], {x[1], x[2]}]];

cc56Assert["B4 bare constant",
  ParsePolynomial[5, {x[1], x[2]}] === {{5, {0, 0}}},
  {{5, {0, 0}}}, ParsePolynomial[5, {x[1], x[2]}]];

(* ============================================================================
   PART C — TropicalEval`Private`normalizeIntegrator[s]  (PRIVATE seam).
   Always returns a STRING "MC" or "VEGAS" (Head === String pinned below,
   not just value equality: this distinguishes it from returning a bare
   Symbol, which would print identically under some Print idioms).
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (C) normalizeIntegrator characterization"];

cc56Assert["C1 \"MC\" -> \"MC\"",
  normalizeIntegrator["MC"] === "MC" && Head[normalizeIntegrator["MC"]] === String,
  "MC", normalizeIntegrator["MC"]];

cc56Assert["C2 \"MonteCarlo\" (alias) -> \"MC\"",
  normalizeIntegrator["MonteCarlo"] === "MC",
  "MC", normalizeIntegrator["MonteCarlo"]];

cc56Assert["C3 \"VEGAS\" -> \"VEGAS\"",
  normalizeIntegrator["VEGAS"] === "VEGAS" &&
    Head[normalizeIntegrator["VEGAS"]] === String,
  "VEGAS", normalizeIntegrator["VEGAS"]];

cc56Assert["C4 \"Vegas\" (alias) -> \"VEGAS\"",
  normalizeIntegrator["Vegas"] === "VEGAS",
  "VEGAS", normalizeIntegrator["Vegas"]];

cc56Assert["C5 unrecognized string \"anything-else\" -> \"MC\" (silent default)",
  normalizeIntegrator["anything-else"] === "MC",
  "MC", normalizeIntegrator["anything-else"]];

cc56Assert["C6 non-string Symbol Automatic -> \"MC\" (tolerated, no error)",
  normalizeIntegrator[Automatic] === "MC",
  "MC", normalizeIntegrator[Automatic]];

cc56Assert["C7 non-string Integer 42 -> \"MC\" (tolerated, no error)",
  normalizeIntegrator[42] === "MC",
  "MC", normalizeIntegrator[42]];

(* ============================================================================
   PART D — TropicalEval`Private`resolveVegasSizing[opt, maxSectorDim, base, kind]
   (PRIVATE seam).  Automatic: base unchanged for maxSectorDim<=4, else
   Ceiling[base * 3^(maxSectorDim-4)].  Non-Automatic: passthrough verbatim.
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (D) resolveVegasSizing characterization"];

cc56Assert["D1 Automatic, dim=2, base=1000 -> 1000",
  resolveVegasSizing[Automatic, 2, 1000, "kind"] === 1000,
  1000, resolveVegasSizing[Automatic, 2, 1000, "kind"]];

cc56Assert["D2 Automatic, dim=4, base=1000 -> 1000",
  resolveVegasSizing[Automatic, 4, 1000, "kind"] === 1000,
  1000, resolveVegasSizing[Automatic, 4, 1000, "kind"]];

cc56Assert["D3 Automatic, dim=5, base=1000 -> 3000",
  resolveVegasSizing[Automatic, 5, 1000, "kind"] === 3000,
  3000, resolveVegasSizing[Automatic, 5, 1000, "kind"]];

cc56Assert["D4 Automatic, dim=7, base=1000 -> 27000",
  resolveVegasSizing[Automatic, 7, 1000, "kind"] === 27000,
  27000, resolveVegasSizing[Automatic, 7, 1000, "kind"]];

cc56Assert["D5 non-Automatic passthrough (12345 unchanged)",
  resolveVegasSizing[12345, 7, 1000, "kind"] === 12345,
  12345, resolveVegasSizing[12345, 7, 1000, "kind"]];

(* ============================================================================
   PART E — emitDomainIndicatorCpp / dropDivVarFromDomain / liftedDomainBooleWL
   (PRIVATE seam).  DomainConstraint association keys per tropical_eval.wl
   ~line 1170 (classifyDomainFor): "LogZ0" -> Log[z0], "MP" -> mp,
   "IndicatorCoeffs" -> {...}.
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (E) emitDomainIndicatorCpp / dropDivVarFromDomain / liftedDomainBooleWL"];

cc56Assert["E1 emitDomainIndicatorCpp[None, <||>] -> \"\"",
  emitDomainIndicatorCpp[None, <||>] === "",
  "\"\"", emitDomainIndicatorCpp[None, <||>]];

$dcE = <|"LogZ0" -> Log[2], "MP" -> 3, "IndicatorCoeffs" -> {1, -1}|>;
$expectedE2 = "    // lifted-sector domain indicator\n" <>
  "    double log_ypstar = (0.6931471805599453 - (1. * log_y[0] + -1. * log_y[1])) * (1.0/3.);\n" <>
  "    if (log_ypstar > 0.0) return cx(0.0, 0.0);\n\n";
(* Round-trip sanity: rebuild the same string with explicit line pieces and
   confirm it matches the \n-escaped literal above (guards against a stray
   editor-introduced whitespace difference in this test file itself). *)
$expectedE2RoundTrip = StringRiffle[{
  "    // lifted-sector domain indicator",
  "    double log_ypstar = (0.6931471805599453 - (1. * log_y[0] + -1. * log_y[1])) * (1.0/3.);",
  "    if (log_ypstar > 0.0) return cx(0.0, 0.0);",
  "",
  ""}, "\n"];
cc56Assert["E1b expectedE2 literal round-trips via explicit \\n pieces",
  $expectedE2 === $expectedE2RoundTrip, $expectedE2, $expectedE2RoundTrip];

cc56Assert["E2 emitDomainIndicatorCpp[assoc, <||>] pinned C++ block",
  emitDomainIndicatorCpp[$dcE, <||>] === $expectedE2,
  $expectedE2, emitDomainIndicatorCpp[$dcE, <||>]];

cc56Assert["E3 dropDivVarFromDomain[None, 1] -> None",
  dropDivVarFromDomain[None, 1] === None, None, dropDivVarFromDomain[None, 1]];

cc56Assert["E4 dropDivVarFromDomain[assoc, 2] drops+reindexes IndicatorCoeffs",
  dropDivVarFromDomain[
      <|"LogZ0" -> Log[2], "MP" -> 3, "IndicatorCoeffs" -> {1, -1, 5}|>, 2] ===
    <|"LogZ0" -> Log[2], "MP" -> 3, "IndicatorCoeffs" -> {1, 5}|>,
  <|"LogZ0" -> Log[2], "MP" -> 3, "IndicatorCoeffs" -> {1, 5}|>,
  dropDivVarFromDomain[
      <|"LogZ0" -> Log[2], "MP" -> 3, "IndicatorCoeffs" -> {1, -1, 5}|>, 2]];

cc56Assert["E5 liftedDomainBooleWL[None, {...}] -> 1",
  liftedDomainBooleWL[None, {1, 2, 3}] === 1,
  1, liftedDomainBooleWL[None, {1, 2, 3}]];

cc56Assert["E6 liftedDomainBooleWL[assoc, {2,1}] -> 1 (on boundary)",
  liftedDomainBooleWL[
      <|"LogZ0" -> Log[2], "MP" -> 3, "IndicatorCoeffs" -> {1, -1}|>, {2, 1}] === 1,
  1, liftedDomainBooleWL[
      <|"LogZ0" -> Log[2], "MP" -> 3, "IndicatorCoeffs" -> {1, -1}|>, {2, 1}]];

(* ============================================================================
   PART F — TropicalEval`Private`cppBadTokens[code]  (PRIVATE seam).
   Returns a LIST of the offending token strings (NOT a string; "" would only
   ever appear as one hypothetical element, never as the clean-case return —
   the clean case returns {}).  Pinned as-observed.
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (F) cppBadTokens characterization"];

cc56Assert["F1 clean snippet -> {} (empty list, not \"\")",
  cppBadTokensFn["cx P0 = exp(log_y[0]);"] === {},
  {}, cppBadTokensFn["cx P0 = exp(log_y[0]);"]];

cc56Assert["F2 detects \"Gamma(\" -> {\"Gamma(\"}",
  cppBadTokensFn["cx foo = Gamma(x);"] === {"Gamma("},
  {"Gamma("}, cppBadTokensFn["cx foo = Gamma(x);"]];

cc56Assert["F3 detects \"DirectedInfinity\" -> {\"DirectedInfinity\"}",
  cppBadTokensFn["cx foo = DirectedInfinity;"] === {"DirectedInfinity"},
  {"DirectedInfinity"}, cppBadTokensFn["cx foo = DirectedInfinity;"]];

cc56Assert["F4 no false-positive on \"P0(\"",
  cppBadTokensFn["cx result = P0(1.0) + Pfull0(2.0);"] === {},
  {}, cppBadTokensFn["cx result = P0(1.0) + Pfull0(2.0);"]];

cc56Assert["F5 no false-positive on \"Vegas(\" (CUBA entry point)",
  cppBadTokensFn["int main() { Vegas(f, 3); }"] === {},
  {}, cppBadTokensFn["int main() { Vegas(f, 3); }"]];

(* ============================================================================
   PART G — imagMonoInfo/realifyMonoA (Im(A) helpers) and
            imagPolyInfo/realifyPolyB (Im(B) helpers)  (PRIVATE seam).
   Spec: MonomialExponents has a complex entry -1+eps+2I (regulator eps
   folded in), PolynomialExponents has a complex entry -1+I*nu with nu a
   declared-real KinematicSymbol.
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (G) imagMonoInfo/realifyMonoA + imagPolyInfo/realifyPolyB"];

$specG = <|
  "MonomialExponents"    -> {-1 + eps + 2 I, 3},
  "PolynomialExponents"  -> {-1 + I*nu},
  "KinematicSymbols"     -> {nu},
  "RegulatorSymbol"      -> eps
|>;

Module[{r},
  r = imagMonoInfo[$specG, eps];
  cc56Assert["G1 imagMonoInfo: {imA, hasImagA} = {{2,0}, True}",
    r === {{2, 0}, True}, {{2, 0}, True}, r];
];

Module[{imA, r},
  imA = imagMonoInfo[$specG, eps][[1]];
  r = realifyMonoA[$specG, imA];
  cc56Assert["G2 realifyMonoA: MonomialExponents -> {-1+eps, 3}",
    r["MonomialExponents"] === {-1 + eps, 3},
    {-1 + eps, 3}, r["MonomialExponents"]];
];

Module[{r},
  r = imagPolyInfo[$specG];
  cc56Assert["G3 imagPolyInfo: {nu} (eps-free Im(B), nu declared real)",
    r === {nu}, {nu}, r];
];

Module[{imB, r},
  imB = imagPolyInfo[$specG];
  r = realifyPolyB[$specG, imB];
  cc56Assert["G4 realifyPolyB: PolynomialExponents -> {-1}",
    r["PolynomialExponents"] === {-1}, {-1}, r["PolynomialExponents"]];
];

(* ============================================================================
   PART H — TropicalEval`Private`emitLogTail[varTerms, polyTerms, paramMap,
   withComment]  (PRIVATE seam).  withComment=True emits the "// Log
   insertion factors" header (G1 / IBP boundary-log function); False omits it
   (IBP term-log function).
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (H) emitLogTail characterization"];

$varTermsH = {{2, 1}, {-1, 3}};
$polyTermsH = {{1, 1}};

$expectedHWith = "\n    // Log insertion factors\n    cx log_sum(0.0, 0.0);\n" <>
  "    log_sum += 2.0 * log_y[0];\n" <>
  "    log_sum += -1.0 * log_y[2];\n" <>
  "    log_sum += 1.0 * std::log(P0);\n";

$expectedHNo = "\n    cx log_sum(0.0, 0.0);\n" <>
  "    log_sum += 2.0 * log_y[0];\n" <>
  "    log_sum += -1.0 * log_y[2];\n" <>
  "    log_sum += 1.0 * std::log(P0);\n";

cc56Assert["H1 emitLogTail withComment=True pinned block",
  emitLogTail[$varTermsH, $polyTermsH, <||>, True] === $expectedHWith,
  $expectedHWith, emitLogTail[$varTermsH, $polyTermsH, <||>, True]];

cc56Assert["H2 emitLogTail withComment=False pinned block (no header line)",
  emitLogTail[$varTermsH, $polyTermsH, <||>, False] === $expectedHNo,
  $expectedHNo, emitLogTail[$varTermsH, $polyTermsH, <||>, False]];

(* ============================================================================
   Summary
   ============================================================================ *)
Print[];
Print["================================================================"];
If[$cc56Pass,
  Print["CC56 PASS  MmaToC/ParsePolynomial + codegen-primitive private seams characterization pinned"],
  Print["CC56 FAIL  failures: ", $cc56Fail]];
Print["================================================================"];
