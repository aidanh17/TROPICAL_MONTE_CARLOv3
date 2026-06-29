(* ============================================================================
   SANDBOX/sandbox_M0_complexmono.wl  --  planAXpDIVv2.md Phase M0
   Verifies the complex-MONOMIAL-exponent (A) obstruction and the proposed fix.

   Test A (lifted, headline):  ProcessSectorLifted with Re(B)-only (A complex)
        => liftcomplex on every cone;  with Re(A)+Re(B) => every cone succeeds.
        (Reproduces LIFTCOMPLISSUE/liftcomplex_issue.tex table; the engine's
         611/611 on the real bubble, here on a small toy.)

   Test D' (unlifted, the new math):  reconstruct the SplitRealImag sector sum
        and confirm the monomial phase  exp(i * sum_a (Im(A).M)_a/aEff_a * log y'_a)
        (a) matches NIntegrate of the complex original  [formula is CORRECT], and
        (b) dropping it gives a wrong answer            [phase is LOAD-BEARING].
        B is kept REAL so the ONLY imaginary structure comes from Im(A).

   Test C (unlifted, latent bug):  the SplitRealImag DRIVER takes Re(B) only and
        leaves A complex (tropical_eval.wl:3968) -> does it silently mis-evaluate
        a complex-A integrand?  Compare Direct (always correct) vs SplitRealImag
        vs NIntegrate.
   ============================================================================ *)

$pkgRoot = "/Users/aidanh/Desktop/TROPICAL_MONTE_CARLOv3";
Get[FileNameJoin[{$pkgRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$pkgRoot, "tropical_eval.wl"}]];
Print["M0: packages loaded.\n"];

SetAttributes[qr, HoldFirst];
qr[e_] := Quiet[e, {TropicalEval::validate, TropicalEval::divergent,
   TropicalEval::liftcomplex, TropicalEval::liftnopivot, TropicalEval::liftdivergent,
   TropicalEval::vegasbudget, General::stop, TropicalFan::polymake,
   NIntegrate::slwcon, NIntegrate::ncvb, NIntegrate::inumr, NIntegrate::eincr,
   NIntegrate::izero, NIntegrate::deuc}];

relErr[got_, ref_] := If[Abs[N[ref]] < 1*^-300, Abs[N[got - ref]],
   Abs[N[got - ref]] / Abs[N[ref]]];

$pass = {}; $fail = {};
rec[name_, ok_, detail_] := (If[TrueQ[ok], AppendTo[$pass, name], AppendTo[$fail, name]];
   Print["  ", If[TrueQ[ok], "PASS", "FAIL"], "  ", name, "  ", detail]);

vars = {x[1], x[2]};
ImAvals = {3/2, -3/2};                       (* Im of monomial exponents (bubble-like) *)
ReAvals = {-1/2, -1/2};
Avals   = ReAvals + I ImAvals;               (* COMPLEX monomial exponents *)

(* high-precision NIntegrate of the ORIGINAL complex integral, via x=u/(1-u) *)
niOriginal[poly_, Aexp_, Bexp_] := Module[{us, sub, jac, ig},
  us  = Table[Unique["u"], {2}];
  sub = Table[x[i] -> us[[i]]/(1 - us[[i]]), {i, 2}];
  jac = Times @@ Table[1/(1 - us[[i]])^2, {i, 2}];
  ig  = (x[1]^Aexp[[1]] x[2]^Aexp[[2]] poly^Bexp /. sub) jac;
  qr@NIntegrate[ig, Evaluate[Sequence @@ ({#, 0, 1} & /@ us)],
    Method -> "GlobalAdaptive", PrecisionGoal -> 6, WorkingPrecision -> 30,
    MaxRecursion -> 80]];

(* ========================================================================
   TEST A — lifted: liftcomplex(Re B only)  vs  success(Re A + Re B)
   ======================================================================== *)
Print["=== TEST A: lifted, complex monomial exponents ==="];
polyA = 1 + x[1] x[2] + 10^6 x[1]^2 + x[2]^2;   (* real coeffs; lift the 10^6 *)
Bcx   = -2 + (3/10) I;                            (* complex B, bubble-like *)
specA = <|"Polynomials" -> {polyA}, "MonomialExponents" -> Avals,
  "PolynomialExponents" -> {Bcx}, "Variables" -> vars,
  "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
liftA = qr@LiftCoefficients[specA, {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 3|>}];
lsA = liftA["LiftedSpec"]; ldA = liftA["LiftData"];
fanA = qr@computeFanScaled[qr@PolytopeVertices[(Times @@ lsA["Polynomials"])^(-1), lsA["Variables"]]];
nA = Length[fanA[[2]]];
Print["  lifted fan: ", nA, " cones;  z0 = ", ldA["z0"]];

specBonly = lsA; specBonly["PolynomialExponents"] = Re[lsA["PolynomialExponents"]];
specBoth  = specBonly; specBoth["MonomialExponents"] = Re[lsA["MonomialExponents"]];
cls[r_]  := Which[r === $Failed, "FAIL", AssociationQ[r] && TrueQ[Lookup[r, "EmptyDomain", False]], "empty",
   AssociationQ[r], "ok", True, "?"];
rB = Table[qr@ProcessSectorLifted[specBonly, fanA[[1]], fanA[[2, s]], s, ldA], {s, nA}];
r2 = Table[qr@ProcessSectorLifted[specBoth,  fanA[[1]], fanA[[2, s]], s, ldA], {s, nA}];
Print["  Re(B)-only  (A complex): ", Tally[cls /@ rB]];
Print["  Re(A)+Re(B) (fix)      : ", Tally[cls /@ r2]];
rec["A1 Re(B)-only fails every cone (liftcomplex)", Count[cls /@ rB, "FAIL"] == nA,
  "fails=" <> ToString[Count[cls /@ rB, "FAIL"]] <> "/" <> ToString[nA]];
rec["A2 Re(A)+Re(B) succeeds every cone", Count[cls /@ r2, "FAIL"] == 0,
  "fails=" <> ToString[Count[cls /@ r2, "FAIL"]] <> "/" <> ToString[nA]];
Print[];

(* ========================================================================
   TEST D' — unlifted: monomial phase  exp(i sum_a (Im(A).M)_a/aEff_a log y'_a)
   B is REAL so the ONLY imaginary part comes from Im(A).
   ======================================================================== *)
Print["=== TEST D': unlifted, monomial phase form + load-bearing ==="];
polyD = 1 + x[1] x[2] + x[1]^2 + x[2]^2;   (* real coeffs, no lifting *)
Breal = -2;
specReD = <|"Polynomials" -> {polyD}, "MonomialExponents" -> ReAvals,
  "PolynomialExponents" -> {Breal}, "Variables" -> vars,
  "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
fanD = qr@computeFanScaled[qr@PolytopeVertices[polyD^(-1), vars]];
nD = Length[fanD[[2]]];
Print["  unlifted fan: ", nD, " cones"];

refD = niOriginal[polyD, Avals, Breal];
Print["  NIntegrate(original, complex A) = ", N[refD]];

computeSplit[useMonoPhase_] := Sum[
  Module[{sd, flat, pf, BRe, aEffR, M, ImAdotM, yv, lg, Q, mag, mph, integ},
    sd = qr@ProcessSector[specReD, fanD[[1]], fanD[[2, s]], s];
    Which[
      sd === $Failed, 0,
      TrueQ[sd["IsDivergent"]], (Print["  WARN: sector ", s, " divergent"]; 0),
      True,
        flat  = sd["FlattenedPolys"]; pf = sd["Prefactor"];
        BRe   = sd["PolynomialExponents"]; aEffR = sd["NewExponents"];
        M     = -Transpose[fanD[[1]][[#]] & /@ fanD[[2, s]]];
        ImAdotM = ImAvals . M;                       (* (Im(A).M)_a, length 2 *)
        yv  = Table[Unique["y"], {2}]; lg = Log /@ yv;
        Q   = Table[Total[(#[[1]] Exp[#[[2]] . lg]) & /@ flat[[j]]], {j, Length[flat]}];
        mag = pf Product[Exp[BRe[[j]] Log[Q[[j]]]], {j, Length[flat]}];
        mph = If[useMonoPhase, I Sum[(ImAdotM[[a]]/aEffR[[a]]) lg[[a]], {a, 2}], 0];
        integ = mag Exp[mph];
        qr@NIntegrate[Evaluate[integ], Evaluate[Sequence @@ ({#, 0, 1} & /@ yv)],
          Method -> "GlobalAdaptive", PrecisionGoal -> 5, MaxRecursion -> 40]
    ]], {s, nD}];

withP = computeSplit[True];  rW = relErr[withP, refD];
noP   = computeSplit[False]; rN = relErr[noP,   refD];
Print["  WITH    mono phase = ", N[withP], "   rel-err = ", ScientificForm[rW, 3]];
Print["  WITHOUT mono phase = ", N[noP],   "   rel-err = ", ScientificForm[rN, 3]];
rec["D'1 mono-phase formula matches NIntegrate (<1%)", rW < 1*^-2,
  "rel-err=" <> ToString[ScientificForm[rW, 3]]];
rec["D'2 mono phase is load-bearing (dropping it is wrong, >5%)", rN > 5*^-2,
  "rel-err=" <> ToString[ScientificForm[rN, 3]]];
Print[];

(* ========================================================================
   TEST C — unlifted DRIVER: Direct vs SplitRealImag (complex A AND B)
   ======================================================================== *)
Print["=== TEST C: unlifted driver, latent silent-wrong check ==="];
specC = <|"Polynomials" -> {polyD}, "MonomialExponents" -> Avals,
  "PolynomialExponents" -> {Bcx}, "Variables" -> vars,
  "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
fanC = fanD;
refC = niOriginal[polyD, Avals, Bcx];
Print["  NIntegrate(original) = ", N[refC]];
wd = FileNameJoin[{$pkgRoot, "INTERFILES", "sandbox_m0"}];
Quiet[CreateDirectory[wd, CreateIntermediateDirectories -> True], CreateDirectory::eexist];
runDriver[mode_] := Module[{r},
  r = qr@Check[EvaluateTropicalMC[specC, fanC, {{}}, "Integrator" -> "MC",
        "NSamples" -> 3000000, "RunChecks" -> False, "Verbose" -> False,
        "WorkingDirectory" -> wd, "ComplexExponentMode" -> mode], $Failed];
  If[AssociationQ[r] && KeyExistsQ[r, "Results"],
    r["Results"][[1]]["Re"] + I r["Results"][[1]]["Im"], $Failed]];
vDir = runDriver["Direct"];
vSpl = runDriver["SplitRealImag"];
Print["  Direct        = ", N[vDir], If[vDir =!= $Failed, "   rel-err vs NI = " <> ToString[ScientificForm[relErr[vDir, refC], 3]], ""]];
Print["  SplitRealImag = ", N[vSpl], If[vSpl =!= $Failed, "   rel-err vs NI = " <> ToString[ScientificForm[relErr[vSpl, refC], 3]], ""]];
If[vDir =!= $Failed && vSpl =!= $Failed,
  rec["C1 Direct matches NIntegrate (sanity, <3%)", relErr[vDir, refC] < 3*^-2,
    "rel-err=" <> ToString[ScientificForm[relErr[vDir, refC], 3]]];
  rec["C2 SplitRealImag DISAGREES with Direct/NI (>5%) => latent bug confirmed",
    relErr[vSpl, refC] > 5*^-2,
    "rel-err=" <> ToString[ScientificForm[relErr[vSpl, refC], 3]]],
  Print["  (driver returned $Failed for a mode; see above)"]];
Print[];

(* ========================================================================
   SUMMARY
   ======================================================================== *)
Print["================================================================"];
Print["M0 PASS: ", Length[$pass], "   FAIL: ", Length[$fail]];
If[Length[$fail] > 0, Print["  failed: ", $fail]];
Print["================================================================"];
