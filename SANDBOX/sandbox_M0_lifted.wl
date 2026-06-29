(* SANDBOX/sandbox_M0_lifted.wl -- planAXpDIVv2.md Phase M0, Test E (lifted).
   Verifies the LIFTED monomial phase, incl. the new Const_A (the z0 term):
     Const_A   = (Im(eaug)_pivot / mp) * logZ0,        eaug = Im(A_aug) . M
     Coeffs_jj = (Im(eaug)_{rIdx[jj]} - Im(eaug)_pivot * mOther_jj/mp) / atilde_jj
   (the analog of the engine's lifted MonoFactorLog for B, with d^aug -> Im(A).M
    and NO rcMin term, since the bare monomial is not a polynomial).
   B is REAL so the ONLY imaginary structure is Im(A).  Reconstruct the lifted
   SplitRealImag sector sum WITH / WITHOUT the monomial phase vs NIntegrate. *)

$pkgRoot = "/Users/aidanh/Desktop/TROPICAL_MONTE_CARLOv3";
Get[FileNameJoin[{$pkgRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$pkgRoot, "tropical_eval.wl"}]];
Print["M0-lifted: packages loaded."];

SetAttributes[qr, HoldFirst];
qr[e_] := Quiet[e, {TropicalEval::validate, TropicalEval::divergent,
   TropicalEval::liftcomplex, TropicalEval::liftnopivot, TropicalEval::liftdivergent,
   General::stop, TropicalFan::polymake, NIntegrate::slwcon, NIntegrate::ncvb,
   NIntegrate::inumr, NIntegrate::eincr, NIntegrate::izero, NIntegrate::deuc}];
relErr[got_, ref_] := If[Abs[N[ref]] < 1*^-300, Abs[N[got - ref]], Abs[N[got - ref]]/Abs[N[ref]]];

vars   = {x[1], x[2]};
ImAvals = {3/2, -3/2}; ReAvals = {-1/2, -1/2};
Avals  = ReAvals + I ImAvals;
poly   = 1 + x[1] x[2] + 10^4 x[1]^2 + x[2]^2;   (* modest extreme coeff: NIntegrate-friendly *)
Breal  = -2;

(* NIntegrate reference of the ORIGINAL complex integral *)
niOriginal[Aexp_, Bexp_] := Module[{us, sub, jac, ig},
  us = Table[Unique["u"], {2}]; sub = Table[x[i] -> us[[i]]/(1 - us[[i]]), {i, 2}];
  jac = Times @@ Table[1/(1 - us[[i]])^2, {i, 2}];
  ig = (x[1]^Aexp[[1]] x[2]^Aexp[[2]] poly^Bexp /. sub) jac;
  qr@NIntegrate[ig, Evaluate[Sequence @@ ({#, 0, 1} & /@ us)],
    Method -> "GlobalAdaptive", PrecisionGoal -> 6, WorkingPrecision -> 30, MaxRecursion -> 80]];
ref = niOriginal[Avals, Breal];
Print["NIntegrate(original, complex A, B=-2) = ", N[ref]];

(* lifted spec + fan *)
spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> Avals,
  "PolynomialExponents" -> {Breal}, "Variables" -> vars,
  "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
lift = qr@LiftCoefficients[spec, {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 3|>}];
ls = lift["LiftedSpec"]; ld = lift["LiftData"]; logZ0 = Log[ld["z0"]];
fan = qr@computeFanScaled[qr@PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]]];
nC = Length[fan[[2]]];
Print["lifted fan: ", nC, " cones;  z0 = ", ld["z0"]];

(* Re(A) AND Re(B): the proposed real decomposition (no liftcomplex) *)
specRe = ls;
specRe["MonomialExponents"]   = Re[ls["MonomialExponents"]];
specRe["PolynomialExponents"] = Re[ls["PolynomialExponents"]];
n = 2; n1 = 3;
ImAaug = Append[ImAvals, 0];   (* aux variable has monomial exponent 0 *)

computeLiftedSplit[useMonoPhase_] := Sum[
  Module[{sd, flat, pf, BRe, atilde, M, pivotP, mVec, rIdx, mp, mOther,
          Imeaug, ConstA, CoeffsA, yv, lg, Q, mag, mph, integ, dc, lyps},
    sd = qr@ProcessSectorLifted[specRe, fan[[1]], fan[[2, s]], s, ld];
    Which[
      sd === $Failed, 0,
      AssociationQ[sd] && TrueQ[Lookup[sd, "EmptyDomain", False]], 0,
      True,
        flat = sd["FlattenedPolys"]; pf = sd["Prefactor"]; BRe = sd["PolynomialExponents"];
        atilde = sd["NewExponents"]; M = sd["RayMatrix"];
        pivotP = sd["PivotIndex"]; mVec = sd["ZRow"]; dc = sd["DomainConstraint"];
        rIdx = DeleteCases[Range[n1], pivotP]; mp = mVec[[pivotP]]; mOther = mVec[[rIdx]];
        Imeaug = ImAaug . M;                              (* length n+1 *)
        ConstA  = (Imeaug[[pivotP]]/mp) logZ0;
        CoeffsA = Table[(Imeaug[[rIdx[[jj]]]] - Imeaug[[pivotP]] mOther[[jj]]/mp)/atilde[[jj]], {jj, n}];
        yv = Table[Unique["y"], {n}]; lg = Log /@ yv;
        Q  = Table[Total[(#[[1]] Exp[#[[2]] . lg]) & /@ flat[[j]]], {j, Length[flat]}];
        mag = pf Product[Exp[BRe[[j]] Log[Q[[j]]]], {j, Length[flat]}];
        mph = If[useMonoPhase, I (ConstA + Sum[CoeffsA[[a]] lg[[a]], {a, n}]), 0];
        integ = mag Exp[mph];
        If[dc =!= None,
          lyps = (N[dc["LogZ0"]] - Total[N[dc["IndicatorCoeffs"]] lg])/N[dc["MP"]];
          integ = integ Boole[lyps <= 0]];
        qr@NIntegrate[Evaluate[integ], Evaluate[Sequence @@ ({#, 0, 1} & /@ yv)],
          Method -> "GlobalAdaptive", PrecisionGoal -> 5, MaxRecursion -> 50]
    ]], {s, nC}];

withP = computeLiftedSplit[True];  rW = relErr[withP, ref];
noP   = computeLiftedSplit[False]; rN = relErr[noP, ref];
Print["lifted WITH    mono phase (incl Const_A) = ", N[withP], "   rel-err = ", ScientificForm[rW, 3]];
Print["lifted WITHOUT mono phase                = ", N[noP],   "   rel-err = ", ScientificForm[rN, 3]];
Print["----"];
Print["E1 lifted mono phase matches NIntegrate (<2%): ", If[rW < 2*^-2, "PASS", "FAIL"], "  (", ScientificForm[rW, 3], ")"];
Print["E2 mono phase is load-bearing (>5% without):   ", If[rN > 5*^-2, "PASS", "FAIL"], "  (", ScientificForm[rN, 3], ")"];
