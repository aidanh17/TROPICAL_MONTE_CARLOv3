(* ============================================================================
   TEST/test_complex_phase4.wl  —  Phase 4 cross-checks (plan.md §7 Phase 4, §6.3)
   Complex exponents end-to-end (Direct), complex-coefficient lifting, and the
   opt-in SplitRealImag VEGAS mode with both bug-log fixes (complex log P =
   TMCv2_BUG_LOG BUG 1; MonoFactorLog = lift_error_log L3).

   Covers cross-checks #1 (NIntegrate), #2 (exact Gamma/Beta), #21b
   (complex-BASE + complex-exponent lifting: SplitRealImag == Direct == ref).

   Tier 1 (symbolic + C++ MC).  Each test prints "<NAME> PASS ..." / "<NAME> FAIL".
   Run:  wolframscript -file TEST/test_complex_phase4.wl
   ============================================================================ *)

SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
v3Root = Directory[];
Get[FileNameJoin[{v3Root, "tropical_eval.wl"}]];
workDir = FileNameJoin[{v3Root, "INTERFILES", "phase4_test"}];
If[!DirectoryQ[workDir], CreateDirectory[workDir, CreateIntermediateDirectories -> True]];

passes = {};
record[name_, ok_] := (AppendTo[passes, ok];
  Print[name, If[ok, " PASS", " FAIL"]]);
rel[a_, b_] := N[Abs[(a - b)/b]];
reportRel[name_, got_, ref_, tol_] := Module[{r = rel[got, ref]},
  Print[name, If[r < tol, " PASS", " FAIL"], "  got=", got, "  ref=", ref,
    "  relerr=", ScientificForm[r, 3], "  tol=", tol];
  AppendTo[passes, r < tol]; r < tol];

(* ---------------------------------------------------------------------------
   #1 / #2 — Direct complex path: Examples 2,15,16,19,21 (symbolic sector sum
   vs NIntegrate; Example 21 also vs exact Gamma).
   --------------------------------------------------------------------------- *)
Print["=== Direct complex path (#1 NIntegrate, #2 exact Gamma) ==="];
Module[{poly, vars, A, spec, verts, fan, vr},
  poly = 1 + 2 x[1]^2 + x[2]^2 + x[1] x[2]^2 + 3 x[1]^2 x[2]; vars = {x[1], x[2]}; A = 2 + I;
  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {0, 0},
    "PolynomialExponents" -> {-A}, "Variables" -> vars, "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
  verts = PolytopeVertices[poly^(-Re[A]), vars]; fan = ComputeDecomposition[verts, "ShowProgress" -> False];
  vr = Quiet@ValidateDecomposition[spec, fan, {}, 3];
  reportRel["EX2  Direct #1", vr["SectorSum"], vr["DirectResult"], 1*^-3]];

Module[{p1, p2, vars, spec, verts, fan, vr},
  p1 = 1 + x[1] + x[2]; p2 = 1 + x[1] x[2]; vars = {x[1], x[2]};
  spec = <|"Polynomials" -> {p1, p2}, "MonomialExponents" -> {0, 0},
    "PolynomialExponents" -> {-(3/2 + I/2), -(1 - I/3)}, "Variables" -> vars,
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
  verts = PolytopeVertices[(p1 p2)^(-1), vars]; fan = ComputeDecomposition[verts, "ShowProgress" -> False];
  vr = Quiet@ValidateDecomposition[spec, fan, {}, 3];
  reportRel["EX15 Direct #1", vr["SectorSum"], vr["DirectResult"], 1*^-3]];

Module[{poly, vars, spec, verts, fan, vr},
  poly = 1 + x[1]^2 + x[2]^2 + x[1] x[2]; vars = {x[1], x[2]};
  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {1/2 + I/3, -1/4 - I/5},
    "PolynomialExponents" -> {-(2 + I)}, "Variables" -> vars, "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
  verts = PolytopeVertices[poly^(-1), vars]; fan = ComputeDecomposition[verts, "ShowProgress" -> False];
  vr = Quiet@ValidateDecomposition[spec, fan, {}, 3];
  reportRel["EX16 Direct #1", vr["SectorSum"], vr["DirectResult"], 1*^-3]];

Module[{poly, vars, spec, verts, fan, vr},
  poly = 1 + 10^-3 x[1]^2 + x[2]^2 + 10^-2 x[1] x[2]^2; vars = {x[1], x[2]};
  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {0, 0},
    "PolynomialExponents" -> {-(2 + I/2)}, "Variables" -> vars, "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
  verts = PolytopeVertices[(1 + x[1]^2 + x[2]^2 + x[1] x[2]^2)^(-1), vars]; fan = ComputeDecomposition[verts, "ShowProgress" -> False];
  vr = Quiet@ValidateDecomposition[spec, fan, {}, 4];
  reportRel["EX19 Direct #1", vr["SectorSum"], vr["DirectResult"], 5*^-3]];

Module[{a1, a2, b, exact, poly, vars, spec, verts, fan, vr},
  a1 = 5/4 + I/3; a2 = 3/2 - I/5; b = 5 + I/2;
  exact = N[Gamma[a1] Gamma[a2] Gamma[b - a1 - a2]/Gamma[b]];
  poly = 1 + x[1] + x[2]; vars = {x[1], x[2]};
  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {a1 - 1, a2 - 1},
    "PolynomialExponents" -> {-b}, "Variables" -> vars, "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
  verts = PolytopeVertices[poly^(-1), vars]; fan = ComputeDecomposition[verts, "ShowProgress" -> False];
  vr = Quiet@ValidateDecomposition[spec, fan, {}, 3];
  reportRel["EX21 Direct #2 exactGamma", vr["SectorSum"], exact, 1*^-4]];

(* Example 21 also through the C++ MC pipeline vs exact Gamma (#2 end-to-end). *)
Module[{a1, a2, b, exact, poly, vars, spec, verts, fan, mc, r, got},
  a1 = 5/4 + I/3; a2 = 3/2 - I/5; b = 5 + I/2;
  exact = N[Gamma[a1] Gamma[a2] Gamma[b - a1 - a2]/Gamma[b]];
  poly = 1 + x[1] + x[2]; vars = {x[1], x[2]};
  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {a1 - 1, a2 - 1},
    "PolynomialExponents" -> {-b}, "Variables" -> vars, "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
  verts = PolytopeVertices[poly^(-1), vars]; fan = ComputeDecomposition[verts, "ShowProgress" -> False];
  mc = TropicalEval`EvaluateTropicalMC[spec, fan, {{}}, "NSamples" -> 2000000,
    "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> workDir, "ComplexExponentMode" -> "Direct"];
  got = If[AssociationQ[mc], mc["Results"][[1]]["Re"] + I mc["Results"][[1]]["Im"], $Failed];
  reportRel["EX21 Direct #2 C++ MC", got, exact, 1*^-2]];

(* ---------------------------------------------------------------------------
   #40 — complex-coefficient lift round-trip identity (exact) + lifted
   decomposition exactness (Direct, real B).
   --------------------------------------------------------------------------- *)
Print["=== Complex-coefficient lifting (#40 identity, exactness) ==="];
Module[{vars, Cc, poly, spec, det, k, rules, lifted, lspec, ld, auxVar, z0, subbed, orig,
        id, magOK, reconOK, verts, fan, vr},
  vars = {x[1], x[2]}; Cc = (1 + I) 1*^-4;
  poly = 1 + Cc x[1]^2 + x[2]^2 + x[1] x[2]^2;
  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {0, 0},
    "PolynomialExponents" -> {-2}, "Variables" -> vars, "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
  det = DetectExtremeCoefficients[spec, 1000]; k = det[[1]]["SuggestedK"];
  rules = {<|"PolyIndex" -> 1, "ExponentVector" -> det[[1]]["ExponentVector"], "k" -> k|>};
  lifted = LiftCoefficients[spec, rules]; lspec = lifted["LiftedSpec"]; ld = lifted["LiftData"];
  auxVar = ld["AuxVariable"]; z0 = ld["z0"];
  subbed = Expand[lspec["Polynomials"][[1]] /. auxVar -> z0]; orig = Expand[poly];
  id = TrueQ[PossibleZeroQ[subbed - orig]];
  magOK = TrueQ[PossibleZeroQ[z0^k - Abs[Cc]]];
  reconOK = TrueQ[PossibleZeroQ[ld["Residuals"][[1]] z0^k - Cc]];
  record["#40 cx-lift round-trip identity exact", id && magOK && reconOK];
  verts = PolytopeVertices[(1 + x[1]^2 + x[2]^2 + x[1] x[2]^2 + auxVar^k x[1]^2)^(-1), lspec["Variables"]];
  fan = ComputeDecomposition[verts, "ShowProgress" -> False];
  vr = Quiet@ValidateLiftedDecomposition[spec, lspec, fan, ld, {}, 6];
  reportRel["#40 lifted decomposition exact (pg6)", vr["SectorSum"], vr["DirectResult"], 1*^-5]];

(* ---------------------------------------------------------------------------
   #21b — complex-BASE (arg P != 0) + complex B, SplitRealImag == Direct == ref.
   UNLIFTED (isolates the codegen phase / BUG-1 + L3): the complex coefficient
   (1+I) on x1*x2 makes the cleared poly complex (arg P != 0).  Both modes run
   through the C++ MC pipeline; they must agree and match NIntegrate.
   --------------------------------------------------------------------------- *)
Print["=== #21b complex-BASE SplitRealImag vs Direct vs ref (unlifted C++ MC) ==="];
Module[{vars, Bex, poly, spec, verts, fan, ref, mcD, mcS, vD, vS, eD, eS, dirD, wd1, wd2},
  vars = {x[1], x[2]}; Bex = -(3/2 + 4 I/5);
  poly = 1 + x[1]^2 + (1 + I) x[1] x[2] + x[2]^2;   (* complex coeff => arg P != 0 *)
  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {0, 0},
    "PolynomialExponents" -> {Bex}, "Variables" -> vars, "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
  verts = PolytopeVertices[poly^(-1), vars]; fan = ComputeDecomposition[verts, "ShowProgress" -> False];
  ref = NIntegrate[poly^Bex, {x[1], 0, Infinity}, {x[2], 0, Infinity},
    MaxRecursion -> 30, PrecisionGoal -> 8, WorkingPrecision -> 30];
  wd1 = FileNameJoin[{workDir, "21b_direct"}]; wd2 = FileNameJoin[{workDir, "21b_split"}];
  Scan[If[!DirectoryQ[#], CreateDirectory[#]] &, {wd1, wd2}];
  mcD = TropicalEval`EvaluateTropicalMC[spec, fan, {{}}, "NSamples" -> 4000000, "RunChecks" -> False,
    "Verbose" -> False, "WorkingDirectory" -> wd1, "ComplexExponentMode" -> "Direct"];
  mcS = TropicalEval`EvaluateTropicalMC[spec, fan, {{}}, "NSamples" -> 4000000, "RunChecks" -> False,
    "Verbose" -> False, "WorkingDirectory" -> wd2, "ComplexExponentMode" -> "SplitRealImag"];
  vD = mcD["Results"][[1]]["Re"] + I mcD["Results"][[1]]["Im"];
  vS = mcS["Results"][[1]]["Re"] + I mcS["Results"][[1]]["Im"];
  eD = mcD["Results"][[1]]["ReErr"] + I mcD["Results"][[1]]["ImErr"];
  eS = mcS["Results"][[1]]["ReErr"] + I mcS["Results"][[1]]["ImErr"];
  Print["  ref=", N[ref]];
  Print["  Direct=", vD, " +/-", eD, "   Split=", vS, " +/-", eS];
  Print["  Direct vs ref=", ScientificForm[rel[vD, ref], 3],
        "   Split vs ref=", ScientificForm[rel[vS, ref], 3],
        "   Split vs Direct=", ScientificForm[rel[vS, vD], 3]];
  record["#21b Split==Direct==ref (unlifted)",
    rel[vD, ref] < 1*^-2 && rel[vS, ref] < 1*^-2 && Abs[vS - vD] < Abs[eS] + Abs[eD] + 0.01 Abs[ref]]];

(* ---------------------------------------------------------------------------
   #21b — LIFTED complex-base + complex B (the BUG-1 regime, lift_error_log L3).
   tiny COMPLEX coefficient forces lifting (residual (1+I)/sqrt2, arg P != 0);
   complex B.  The lifted-Direct path correctly REFUSES complex B (liftcomplex,
   needs a real domain indicator) — so SplitRealImag is the only way to lift a
   complex-exponent integrand.  We validate the lifted SplitRealImag sector sum
   (exact decomposition + phase, via NIntegrate — NO MC noise) against the
   unlifted reference, and confirm that DROPPING MonoFactorLog gives a wrong
   answer (so the L3 term is both present and correct).  This is the decisive
   lifted check; the sampled MC value carries the documented L5 heavy-tail
   variance (one HasConstantTerm=False sector) and is not the gate.
   --------------------------------------------------------------------------- *)
Print["=== #21b LIFTED complex-base SplitRealImag (MonoFactorLog / BUG-1, exact) ==="];
Module[{vars, Cc, Bex, poly, shadow, spec, ref, det, k, rules, lifted, lspec, ld,
        auxVar, id40, verts, fan, n, specRe, imB, liftedSplit, withMFL, noMFL},
  vars = {x[1], x[2]}; Cc = (1 + I) 1*^-4; Bex = -(2 + I/2);
  poly = 1 + Cc x[1]^2 + x[1] + x[2]^2; shadow = 1 + x[1]^2 + x[1] + x[2]^2; n = 2;
  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {0, 0},
    "PolynomialExponents" -> {Bex}, "Variables" -> vars, "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
  ref = NIntegrate[poly^Bex, {x[1], 0, Infinity}, {x[2], 0, Infinity}, MaxRecursion -> 40, PrecisionGoal -> 8];
  det = DetectExtremeCoefficients[spec, 1000]; k = det[[1]]["SuggestedK"];
  rules = {<|"PolyIndex" -> 1, "ExponentVector" -> det[[1]]["ExponentVector"], "k" -> k|>};
  lifted = LiftCoefficients[spec, rules]; lspec = lifted["LiftedSpec"]; ld = lifted["LiftData"]; auxVar = ld["AuxVariable"];
  id40 = TrueQ[PossibleZeroQ[Expand[lspec["Polynomials"][[1]] /. auxVar -> ld["z0"]] - Expand[poly]]];
  record["#21b lifted #40 round-trip identity exact", id40];
  verts = PolytopeVertices[(shadow + auxVar^k x[1]^2)^(-1), lspec["Variables"]];
  fan = ComputeDecomposition[verts, "ShowProgress" -> False];
  specRe = lspec; specRe["PolynomialExponents"] = Re[lspec["PolynomialExponents"]];
  imB = Im[lspec["PolynomialExponents"]];
  liftedSplit[useMFL_] := Sum[
    Module[{sd, flatRe, pfRe, BRe, mfl, dc, yv, lg, Q, mag, mterm, integ, lyps},
      sd = Quiet@ProcessSectorLifted[specRe, fan[[1]], fan[[2, s]], s, ld];
      Which[sd === $Failed, 0, KeyExistsQ[sd, "EmptyDomain"], 0, True,
        flatRe = sd["FlattenedPolys"]; pfRe = sd["Prefactor"]; BRe = sd["PolynomialExponents"];
        mfl = sd["MonoFactorLog"]; dc = sd["DomainConstraint"];
        yv = Table[Unique["y"], {n}]; lg = Log /@ yv;
        Q = Table[Total[(#[[1]] Exp[Total[#[[2]] lg]]) & /@ flatRe[[j]]], {j, Length[flatRe]}];
        mag = pfRe Product[Exp[BRe[[j]] Log[Q[[j]]]], {j, Length[flatRe]}];
        mterm = If[useMFL,
          Sum[I imB[[j]] (Log[Q[[j]]] + mfl[[j]]["Const"] + Sum[mfl[[j]]["Coeffs"][[i]] lg[[i]], {i, n}]), {j, Length[flatRe]}],
          Sum[I imB[[j]] Log[Q[[j]]], {j, Length[flatRe]}]];
        integ = mag Exp[mterm];
        If[dc =!= None, lyps = (N[dc["LogZ0"]] - Total[N[dc["IndicatorCoeffs"]] lg])/N[dc["MP"]];
          integ = integ Boole[lyps <= 0]];
        Quiet@NIntegrate[Evaluate[integ], Evaluate[Sequence @@ ({#, 0, 1} & /@ yv)],
          MaxRecursion -> 18, PrecisionGoal -> 4, Method -> "GlobalAdaptive"]
      ]], {s, Length[fan[[2]]]}];
  withMFL = liftedSplit[True]; noMFL = liftedSplit[False];
  Print["  ref=", N[ref]];
  Print["  lifted Split WITH  MonoFactorLog = ", N[withMFL], "  relErr=", ScientificForm[rel[withMFL, ref], 3]];
  Print["  lifted Split WITHOUT MonoFactorLog = ", N[noMFL], "  relErr=", ScientificForm[rel[noMFL, ref], 3], "  (must be wrong)"];
  record["#21b lifted SplitRealImag == ref (exact decomposition+phase)", rel[withMFL, ref] < 5*^-3];
  record["#21b MonoFactorLog is REQUIRED (dropping it is wrong)", rel[noMFL, ref] > 5*^-2]];

(* ---------------------------------------------------------------------------
   #22 — Determinism / seed reproducibility via the runtime argv[5] override
   (plan.md §3.4, §8.6, cross-check #22; Tree B SeedBase).  Restored in Phase 4.0.
   ONE compiled binary is exercised three ways:
     (a) same seed (no argv[5], twice)        -> BYTE-IDENTICAL output;
     (b) two distinct argv[5] seeds           -> DISTINCT streams;
     (c) each stream within 5 sigma of the exact reference.
   Testbed: P = 1+x1^2+x2^2, B=-3, no kinematics  => exact integral = Pi/8 (real).
   This needs the seed_base argv[5] parse + `seed = seed_base + kp` in the MC
   main; before the fix SeedBase->99 and ->42 emitted the identical seed and
   the override was a no-op (#22 unsatisfiable). grep argv[5] tropical_eval.wl>0.
   --------------------------------------------------------------------------- *)
Print["=== #22 seed reproducibility (runtime argv[5] override, ONE binary) ==="];
Module[{poly, vars, spec, verts, fan, sd, conv, wd, src, bin, kin, exact,
        runSeed, runDefault, oA, oB, dflt1, dflt2, s123, s777, p123, p777,
        re123, ree123, re777, ree777, dist, byteEq, within5,
        compiledOK, argv5InSrc},
  exact = N[Pi/8];
  poly = 1 + x[1]^2 + x[2]^2; vars = {x[1], x[2]};
  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {0, 0},
    "PolynomialExponents" -> {-3}, "Variables" -> vars,
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
  verts = PolytopeVertices[poly^(-3), vars];
  fan = ComputeDecomposition[verts, "ShowProgress" -> False];
  sd = Table[ProcessSector[spec, fan[[1]], fan[[2, s]], s, "Verbose" -> False],
    {s, Length[fan[[2]]]}];
  conv = Select[sd, AssociationQ[#] && !#["IsDivergent"] &];
  wd = FileNameJoin[{workDir, "seed22"}];
  If[!DirectoryQ[wd], CreateDirectory[wd]];
  src = FileNameJoin[{wd, "seed22.cpp"}];
  bin = FileNameJoin[{wd, "seed22_bin"}];
  kin = FileNameJoin[{wd, "kin.txt"}];
  Export[kin, "1\n", "Text"];   (* N_PARAMS=0 => count of kinematic points *)

  (* Compile ONE binary at the default SeedBase (42). *)
  GenerateCppMonteCarlo[conv, {}, spec, src, "Integrator" -> "MonteCarlo", "SeedBase" -> 42];
  argv5InSrc = StringContainsQ[Import[src, "Text"], "argv[5]"];
  record["#22 emitted MC main parses argv[5] (seed override present)", argv5InSrc];
  compiledOK = (CompileCpp[src, bin] =!= $Failed) && FileExistsQ[bin];
  record["#22 single binary compiled", compiledOK];

  If[compiledOK,
    (* Helper: run the ONE binary, optional seed arg, return {re, im, reErr, imErr}. *)
    runSeed[seedStr_] := Module[{out, r},
      out = FileNameJoin[{wd, "out_" <> seedStr <> ".txt"}];
      r = RunProcess[{bin, kin, out, "4000000", "1", seedStr}];
      If[r["ExitCode"] != 0, $Failed, ToExpression /@ StringSplit[StringTrim@Import[out, "Text"]]]];
    runDefault[tag_] := Module[{out, r},
      out = FileNameJoin[{wd, "outdef_" <> tag <> ".txt"}];
      r = RunProcess[{bin, kin, out, "4000000", "1"}];   (* NO argv[5]: compile-time seed *)
      If[r["ExitCode"] != 0, $Failed, {out, StringTrim@Import[out, "Text"]}]];

    (* (a) same seed (no override), twice => byte-identical output file. *)
    dflt1 = runDefault["1"]; dflt2 = runDefault["2"];
    byteEq = (dflt1 =!= $Failed && dflt2 =!= $Failed && dflt1[[2]] === dflt2[[2]]);
    Print["  default-seed run #1 line: ", dflt1[[2]]];
    Print["  default-seed run #2 line: ", dflt2[[2]]];
    record["#22 same seed => byte-identical output", byteEq];

    (* (b)+(c) two distinct argv[5] seeds on the SAME binary. *)
    s123 = runSeed["123"]; s777 = runSeed["777"];
    If[s123 =!= $Failed && s777 =!= $Failed,
      {re123, p123, ree123} = {s123[[1]], s123[[2]], s123[[3]]};
      {re777, p777, ree777} = {s777[[1]], s777[[2]], s777[[3]]};
      dist = (re123 =!= re777);   (* distinct random streams => different estimate *)
      within5 = (Abs[re123 - exact] <= 5 ree123) && (Abs[re777 - exact] <= 5 ree777);
      Print["  exact = Pi/8 = ", exact];
      Print["  seed=123: Re=", re123, " +/-", ree123, "  (", Abs[re123 - exact]/ree123, " sigma)"];
      Print["  seed=777: Re=", re777, " +/-", ree777, "  (", Abs[re777 - exact]/ree777, " sigma)"];
      Print["  streams distinct (Re differ): ", dist];
      record["#22 distinct argv[5] seeds => distinct streams", dist];
      record["#22 each seeded stream within 5 sigma of Pi/8", within5],
      record["#22 distinct argv[5] seeds => distinct streams", False];
      record["#22 each seeded stream within 5 sigma of Pi/8", False]],
    record["#22 same seed => byte-identical output", False];
    record["#22 distinct argv[5] seeds => distinct streams", False];
    record["#22 each seeded stream within 5 sigma of Pi/8", False]]];

(* --------------------------------------------------------------------------- *)
Print[];
Print["PHASE4-TEST ", If[And @@ passes, "PASS", "FAIL"], "  (",
  Count[passes, True], "/", Length[passes], " checks passed)"];
