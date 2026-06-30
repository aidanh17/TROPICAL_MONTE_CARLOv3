(* ============================================================================
   TEST/cc_52.wl  —  Cross-check #52 (planIBPMULTIDIV.md §3.5 / §6.4)
   The lifted-sector domain-indicator complex-typing fix, plus the documented
   scope boundary of the lifted divergent multi-soft case.

   §3.5: emitDomainIndicatorCpp emitted
       double log_ypstar = (logZ0 - (sum)) * (1.0/mp);
   When any of logZ0/mp/IndicatorCoeffs is a machine-complex literal Complex[x,0.]
   (the noise realifyMonoA warns about), mmaToCInternal emits cx(x,0.0) and the
   `double = <complex>` assignment does NOT compile (g++).  The SplitRealImag
   domain map is built from Re(A)/Re(B), so the indicator is PROVABLY REAL; the
   fix takes .real() — but ONLY when a cx(...) literal is actually present, so the
   real-typed path stays byte-identical (#25).

   Part A (the bug-fix gate):
     - real-typed dc  -> NO ".real()" (byte-identical to the pre-fix form), and
     - complex-typed dc -> ".real()" present AND the emitted block COMPILES.
   Part B (scope boundary, planIBPMULTIDIV.md §5 / §8): a LIFTED, co-located
     one-pole + one-off-axis sector is refused cleanly at ProcessSectorLifted
     (TropicalEval::liftnopivot, ndiv>1) -> $Failed, never a wrong number
     (invariant #3).  (The UNLIFTED one-pole + N-off-axis case is delivered and
     gated by cc_51; lifted single-off-axis end-to-end is gated by cc_47 Part C.)

   Tier 1: WL + g++.  Run:  wolframscript -file TEST/cc_52.wl
   ============================================================================ *)

$pkgRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{$pkgRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$pkgRoot, "tropical_eval.wl"}]];
Print["CC52: packages loaded."];
Print[];

$cc52Pass = True; $cc52Fail = {};
cc52Assert[label_String, test_, exp_String, got_String] :=
  If[TrueQ[test],
    Print["CC52 PASS  [", label, "]  expected=", exp, "  got=", got],
    Print["CC52 FAIL  [", label, "]  expected=", exp, "  got=", got];
    $cc52Pass = False; AppendTo[$cc52Fail, label]];

SetAttributes[qr, HoldFirst];
qr[e_] := Quiet[e, {TropicalEval::validate, TropicalEval::liftnopivot,
   TropicalEval::splitdivmono, TropicalEval::nestedIBP, General::stop,
   TropicalFan::polymake, Power::infy}];

emitDC = TropicalEval`Private`emitDomainIndicatorCpp;

(* ============================================================================
   PART A — the §3.5 domain-indicator complex-typing fix
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (A) emitDomainIndicatorCpp: real byte-identity + complex .real() compiles"];

Module[{dcReal, dcCx, outReal, outCx, cppFile, binFile, harness, compileOK},
  dcReal = <|"LogZ0" -> Log[100], "MP" -> 2,
             "IndicatorCoeffs" -> {1/2, -1/3}|>;
  (* Machine-complex twin: Complex[x, 0.] in every numeric quantity (the noise
     realifyMonoA leaves after eps-pinning). *)
  dcCx = <|"LogZ0" -> Complex[N[Log[100]], 0.], "MP" -> Complex[2., 0.],
           "IndicatorCoeffs" -> {Complex[0.5, 0.], Complex[-1./3, 0.]}|>;

  outReal = emitDC[dcReal, <||>];
  outCx   = emitDC[dcCx,   <||>];
  Print["   real-dc block:\n", outReal];
  Print["   cx-dc block:\n",   outCx];

  (* (1) Real-typed path is byte-identical to the pre-fix form: no ".real()". *)
  cc52Assert["A real-typed indicator unchanged (no .real(), #25)",
    !StringContainsQ[outReal, ".real()"] &&
      StringContainsQ[outReal, "double log_ypstar = ("],
    "no .real()", "contains .real()=" <> ToString[StringContainsQ[outReal, ".real()"]]];

  (* (2) Complex-typed path takes .real() (so `double = <complex>` compiles). *)
  cc52Assert["A complex-typed indicator takes .real()",
    StringContainsQ[outCx, ".real()"] && StringContainsQ[outCx, "cx("],
    ".real() present", "contains .real()=" <> ToString[StringContainsQ[outCx, ".real()"]]];

  (* (3) The complex-typed block actually COMPILES (the live bug-fix gate). *)
  harness = "#include <complex>\n#include <cmath>\n" <>
    "using cx = std::complex<double>;\n" <>
    "cx test_indicator(const double* log_y) {\n" <> outCx <>
    "    return cx(1.0, 0.0);\n}\n" <>
    "int main(){ double y[2]={0.1,0.2}; volatile cx v = test_indicator(y); (void)v; return 0; }\n";
  cppFile = FileNameJoin[{$pkgRoot, "TEST", "INTERFILES", "cc52_indicator.cpp"}];
  binFile = FileNameJoin[{$pkgRoot, "TEST", "INTERFILES", "cc52_indicator"}];
  Quiet[CreateDirectory[DirectoryName[cppFile],
    CreateIntermediateDirectories -> True], {CreateDirectory::eexist}];
  Export[cppFile, harness, "Text"];
  compileOK = (Run["g++ -std=c++17 -O2 -o " <> binFile <> " " <> cppFile <>
                   " 2> " <> binFile <> ".err"] == 0);
  cc52Assert["A complex-typed indicator block compiles (g++)",
    compileOK, "g++ exit 0", "compileOK=" <> ToString[compileOK]];
];

(* ============================================================================
   PART B — scope boundary: lifted co-located one-pole + one-off-axis refuses
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (B) lifted co-located one-pole + one-off-axis -> clean $Failed (liftnopivot)"];

Module[{th, eps, polys, pe, Avals, spec, lift, ls, ld, fan, res, wd},
  th = 1/2; eps = Symbol["epsCC52B"];
  polys = {1 + x[1] + x[2] + 10^4 x[2]^2};   (* extreme coeff on x2^2 -> lift *)
  pe    = {-3};
  Avals = {2 eps - 1, I th - 1};             (* x1 real pole; x2 off-axis *)
  spec = <|"Polynomials" -> polys, "MonomialExponents" -> Avals,
    "PolynomialExponents" -> pe, "Variables" -> {x[1], x[2]},
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> eps|>;
  lift = qr@LiftCoefficients[spec,
    {<|"PolyIndex" -> 1, "ExponentVector" -> {0, 2}, "k" -> 3|>}];
  ls = lift["LiftedSpec"]; ld = lift["LiftData"];
  fan = qr@computeFanScaled[
    qr@PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]]];
  wd = FileNameJoin[{$pkgRoot, "TEST", "INTERFILES", "cc52", "liftB"}];
  Quiet[CreateDirectory[wd, CreateIntermediateDirectories -> True],
        {CreateDirectory::eexist}];

  res = qr@EvaluateTropicalMC[ls, fan, {{}}, "LiftData" -> ld, "Method" -> "IBP",
    "ComplexExponentMode" -> "SplitRealImag", "Integrator" -> "MC",
    "NSamples" -> 100000, "RunChecks" -> False, "Verbose" -> False,
    "WorkingDirectory" -> wd];
  cc52Assert["B lifted co-located pole+off-axis refuses cleanly ($Failed)",
    res === $Failed, "$Failed", ToString[Head[res]]];
];

(* ============================================================================
   PART C — scope boundary: a strictly power-divergent off-axis direction
   (Re(alpha0) < 0, theta != 0) is genuinely divergent (the oscillation does not
   regulate a magnitude divergence), so it must refuse cleanly rather than
   assemble the logarithmic 1/(i theta) prefactor (planIBPMULTIDIV.md §3.1).
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (C) power-divergent off-axis (Re<0) -> clean $Failed (offaxispow)"];

Module[{th, eps, polys, pe, Avals, spec, fan, sds, div, res},
  th = 1/2; eps = Symbol["epsCC52C"];
  polys = {1 + x[1] + x[2]}; pe = {-3};
  (* x2 off-axis with Re(alpha0) = -1/2 < 0 (power-divergent), theta = 1/2 *)
  Avals = {2 eps - 1, -1/2 + I th - 1};
  spec = <|"Polynomials" -> polys, "MonomialExponents" -> Avals,
    "PolynomialExponents" -> pe, "Variables" -> {x[1], x[2]},
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> eps|>;
  fan = qr@computeFanScaled[
    qr@PolytopeVertices[(Times @@ polys)^(-1), {x[1], x[2]}]];
  sds = Table[ProcessSector[spec, fan[[1]], fan[[2, s]], s], {s, Length[fan[[2]]]}];
  div = Select[sds, AssociationQ[#] && #["IsDivergent"] &];
  res = qr[IBPProcessSector[#, spec] & /@ div];
  cc52Assert["C power-divergent off-axis refuses cleanly ($Failed)",
    MemberQ[res, $Failed], "some $Failed", "any$Failed=" <> ToString[MemberQ[res, $Failed]]];
];

Print[];
Print["================================================================"];
If[$cc52Pass,
  Print["CC52 PASS  domain-indicator complex-typing fix verified (real byte-",
        "identical, complex .real() compiles); lifted co-located one-pole+",
        "off-axis refused cleanly (documented scope boundary).  failures={}"],
  Print["CC52 FAIL  failed: ", $cc52Fail]];
Print["================================================================"];
If[!$cc52Pass, Quit[1]];
