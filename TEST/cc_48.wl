(* ============================================================================
   TEST/cc_48.wl  --  v3 Cross-Check #48
   plan.md section 8.1 feature-coverage matrix -- fills the empty F2 x F4b cell:
   COMPLEX EXPONENTS x DIVERGENCE, *UNLIFTED* (Direct == SplitRealImag == exact Gamma).

   Why this is new (gap analysis, 2026-06):
     * F4b (complex *monomial* exponents A) was covered ONLY by convergent
       cross-checks (#2, #21).  The single divergent complex-A test (cc_47) is
       LIFTED.  There was NO unlifted divergent complex-monomial-A test at all.
     * The unlifted SplitRealImag x divergence path (the Im(A) MonomialPhaseLog /
       Im(B) MonoFactorLog phase machinery WITHOUT the lifting machinery) was
       only ever exercised bundled with lifting (cc_46/cc_47).  This isolates it.

   Both parts use a Dirichlet/Beta integrand with an EXACT closed-form Laurent
   (independent oracle -- shares no code with the engine), and additionally
   cross-check the two independent divergence algorithms (IBP vs Subtraction).
   Single 1/eps pole per sector (L1/L8); the imaginary exponent sits on a
   NON-divergent direction (Case A), so the pole stays at eps -> 0.

   PART A -- complex MONOMIAL A (real B): the actual missing F2 x F4b cell.
       I_A(eps) = Int x1^{2eps-1} x2^{I c} (1+x1+x2)^{-B} dx1 dx2
       Dirichlet:  = Gamma(2eps) Gamma(1+I c) Gamma(B-2eps-1-I c) / Gamma(B)
       pole   = (1/2) Gamma(1+I c) Gamma(B-1-I c) / Gamma(B)   (genuinely complex)
       finite = Gamma(1+I c) Gamma(B-1-I c)/Gamma(B) * (-EulerGamma - psi(B-1-I c))
     B = 3 (real), c = 1/2 -- the ONLY complex feature is the monomial exponent.

   PART B -- complex POLYNOMIAL B (real A), UNLIFTED, via SplitRealImag:
       I_B(eps) = Int x1^{2eps-1} ((1+x1)(1+x2))^{-B} dx1 dx2, with
       (1+x1)(1+x2) = 1+x1+x2+x1 x2 (a 4-term square Newton polytope), factorizes:
                   = Gamma(2eps) Gamma(B-2eps) Gamma(B-1) / Gamma(B)^2
       pole   = 1/(2(B-1));   finite = (-EulerGamma - psi(B)) / (B-1)
     B = 2 + I/2.  Isolates the unlifted complex-B SplitRealImag divergent path
     that cc_46 only ever tested WITH lifting.

   PASS criteria (Re and Im parts independently):
     pole vs exact Gamma   < 5e-3   (every route)
     finite vs exact Gamma < 2.5e-2 (IBP tight; Subtraction looser via eps-fit)
     IBP vs Subtraction pole agree < 1e-2
     Direct == SplitRealImag on the IBP path (option must not perturb the result)
     pole genuinely complex: |Im(pole)| > 0.01

   Tier 1 (no external deps): WL + g++.  CUBA/VEGAS not required.
   ============================================================================ *)

$v3Root = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{$v3Root, "tropical_fan.wl"}]];
Get[FileNameJoin[{$v3Root, "tropical_eval.wl"}]];

SetAttributes[qr, HoldFirst]; qr[e_] := Quiet[e];

wdir[tag_] := Module[{d = FileNameJoin[{$v3Root, "TEST", "INTERFILES", "cc48", tag}]},
  Quiet[CreateDirectory[d, CreateIntermediateDirectories -> True], {CreateDirectory::eexist}]; d];

sci[x_] := Module[{a = N[Abs[x]], e}, If[a == 0., "0",
   e = Floor[Log10[a]]; ToString[NumberForm[a/10^e, {3, 2}]] <> "e" <> ToString[e]]];
fmtC[z_] := ToString[NumberForm[N[Re[z]], {5, 4}]] <> If[Im[z] >= 0, "+", "-"] <>
   ToString[NumberForm[N[Abs[Im[z]]], {5, 4}]] <> "I";

allPass = True;
assert[label_, ok_, detail_] := (
  If[TrueQ[ok], Print["CC48 PASS ", label, "  (", detail, ")"],
     Print["CC48 FAIL ", label, "  (", detail, ")"]];
  allPass = allPass && TrueQ[ok]);

tolPole = 5.*^-3; tolFin = 2.5*^-2; tolAgree = 1.*^-2; tolSplit = 1.*^-6;
epsVals = {0.01, 0.02, 0.04};
NS = 1500000;

poleFinIBP[res_] := If[AssociationQ[res] && KeyExistsQ[res, "Results"],
   {res["Results"][[1]]["PoleCoefficient"], res["Results"][[1]]["FinitePart"]}, $Failed];
poleFinSub[res_] := If[AssociationQ[res] && KeyExistsQ[res, "Results"],
   {res["Results"][[1]]["Pole"], res["Results"][[1]]["Finite"]}, $Failed];

cplxClose[a_, b_, tol_] := NumericQ[Abs[a - b]] &&
   Abs[Re[a] - Re[b]] < tol && Abs[Im[a] - Im[b]] < tol;

Print[""];
Print["================================================================"];
Print["  v3 Cross-Check #48: complex exponents x divergence (UNLIFTED)"];
Print["  fills plan 8.1 F2 x F4b ; isolates SplitRealImag from lifting"];
Print["================================================================"];

(* ===========================================================================
   PART A -- complex MONOMIAL A, real B.  The missing F2 x F4b cell.
   =========================================================================== *)
Print[""];
Print["----------------------------------------------------------------"];
Print["  (A) complex MONOMIAL exponent A x divergence (real B=3, c=1/2)"];
Module[{cc = 1/2, Bv = 3, eps, vars, poly, spec, fan,
        poleRef, finRef, dIBP, dSub, sIBP, sSub, ok, pdI, pdS, psI, psS},
  poleRef = N[(1/2) Gamma[1 + I cc] Gamma[Bv - 1 - I cc]/Gamma[Bv], 20];
  finRef  = N[Gamma[1 + I cc] Gamma[Bv - 1 - I cc]/Gamma[Bv] *
              (-EulerGamma - PolyGamma[0, Bv - 1 - I cc]), 20];
  Print["   exact pole   = ", fmtC[poleRef]];
  Print["   exact finite = ", fmtC[finRef]];

  eps = Symbol["epsA48"]; vars = {x[1], x[2]}; poly = 1 + x[1] + x[2];
  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {2 eps - 1, I cc},
     "PolynomialExponents" -> {-Bv}, "Variables" -> vars,
     "KinematicSymbols" -> {}, "RegulatorSymbol" -> eps|>;
  fan = qr@ComputeDecomposition[qr@PolytopeVertices[poly^(-1), vars]];
  assert["A0 fan built", ListQ[fan] && Length[fan] == 2 && Length[fan[[2]]] > 0,
     "sectors=" <> ToString[If[ListQ[fan], Length[fan[[2]]], fan]]];

  dIBP = poleFinIBP@qr@EvaluateTropicalMCIBP[spec, fan, {{}}, "Integrator" -> "MC",
     "NSamples" -> NS, "RunChecks" -> False, "Verbose" -> False,
     "WorkingDirectory" -> wdir["A_dirIBP"]];
  dSub = poleFinSub@qr@LaurentFromSubtraction[spec, fan, {{}}, "Integrator" -> "MC",
     "NSamples" -> NS, "EpsilonValues" -> epsVals, "WorkingDirectory" -> wdir["A_dirSub"]];
  sIBP = poleFinIBP@qr@EvaluateTropicalMC[spec, fan, {{}}, "Method" -> "IBP",
     "ComplexExponentMode" -> "SplitRealImag", "Integrator" -> "MC", "NSamples" -> NS,
     "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> wdir["A_splIBP"]];
  sSub = poleFinSub@qr@LaurentFromSubtraction[spec, fan, {{}},
     "ComplexExponentMode" -> "SplitRealImag", "Integrator" -> "MC", "NSamples" -> NS,
     "EpsilonValues" -> epsVals, "WorkingDirectory" -> wdir["A_splSub"]];

  ok = AllTrue[{dIBP, dSub, sIBP, sSub}, ListQ];
  assert["A-routes returned results", ok,
     "DirIBP/DirSub/SplIBP/SplSub all returned Association[Results]"];
  If[ok, (
    {pdI, pdS, psI, psS} = {dIBP[[1]], dSub[[1]], sIBP[[1]], sSub[[1]]};
    Print["   DirIBP pole=", fmtC[pdI], "  DirSub pole=", fmtC[pdS]];
    Print["   SplIBP pole=", fmtC[psI], "  SplSub pole=", fmtC[psS]];
    assert["A1 Direct-IBP pole vs exact Gamma", cplxClose[pdI, poleRef, tolPole],
      "dRe=" <> sci[Re[pdI - poleRef]] <> " dIm=" <> sci[Im[pdI - poleRef]]];
    assert["A2 Direct-Sub pole vs exact Gamma", cplxClose[pdS, poleRef, tolPole],
      "dRe=" <> sci[Re[pdS - poleRef]] <> " dIm=" <> sci[Im[pdS - poleRef]]];
    assert["A3 Split-IBP pole vs exact Gamma (complex-A phase path runs)",
      cplxClose[psI, poleRef, tolPole],
      "dRe=" <> sci[Re[psI - poleRef]] <> " dIm=" <> sci[Im[psI - poleRef]]];
    assert["A4 Split-Sub pole vs exact Gamma", cplxClose[psS, poleRef, tolPole],
      "dRe=" <> sci[Re[psS - poleRef]] <> " dIm=" <> sci[Im[psS - poleRef]]];
    assert["A5 IBP vs Subtraction pole agree (Direct)", Abs[pdI - pdS] < tolAgree,
      "|dpole|=" <> sci[Abs[pdI - pdS]]];
    assert["A6 SplitRealImag == Direct on IBP path (option inert here)",
      Abs[psI - pdI] < tolSplit, "|dpole|=" <> sci[Abs[psI - pdI]]];
    assert["A7 Direct-IBP finite vs exact Gamma", cplxClose[dIBP[[2]], finRef, tolFin],
      "dRe=" <> sci[Re[dIBP[[2]] - finRef]] <> " dIm=" <> sci[Im[dIBP[[2]] - finRef]]];
    assert["A8 pole genuinely complex (Im != 0)", Abs[Im[poleRef]] > 0.01,
      "Im(pole)=" <> sci[Im[poleRef]]];
  )];
];

(* ===========================================================================
   PART B -- complex POLYNOMIAL B, real A, UNLIFTED, via SplitRealImag.
   =========================================================================== *)
Print[""];
Print["----------------------------------------------------------------"];
Print["  (B) complex POLYNOMIAL exponent B x divergence, UNLIFTED SplitRealImag"];
Module[{Bv = 2 + I/2, eps, vars, poly, spec, fan,
        poleRef, finRef, dIBP, sIBP, sSub, ok, pdI, psI, psS},
  poleRef = N[1/(2 (Bv - 1)), 20];
  finRef  = N[(-EulerGamma - PolyGamma[0, Bv])/(Bv - 1), 20];
  Print["   exact pole   = ", fmtC[poleRef]];
  Print["   exact finite = ", fmtC[finRef]];

  eps = Symbol["epsB48"]; vars = {x[1], x[2]}; poly = 1 + x[1] + x[2] + x[1] x[2];
  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {2 eps - 1, 0},
     "PolynomialExponents" -> {-Bv}, "Variables" -> vars,
     "KinematicSymbols" -> {}, "RegulatorSymbol" -> eps|>;
  fan = qr@ComputeDecomposition[qr@PolytopeVertices[poly^(-1), vars]];
  assert["B0 fan built (4-term square polytope)",
     ListQ[fan] && Length[fan] == 2 && Length[fan[[2]]] > 0,
     "sectors=" <> ToString[If[ListQ[fan], Length[fan[[2]]], fan]]];

  dIBP = poleFinIBP@qr@EvaluateTropicalMCIBP[spec, fan, {{}}, "Integrator" -> "MC",
     "NSamples" -> NS, "RunChecks" -> False, "Verbose" -> False,
     "WorkingDirectory" -> wdir["B_dirIBP"]];
  sIBP = poleFinIBP@qr@EvaluateTropicalMC[spec, fan, {{}}, "Method" -> "IBP",
     "ComplexExponentMode" -> "SplitRealImag", "Integrator" -> "MC", "NSamples" -> NS,
     "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> wdir["B_splIBP"]];
  sSub = poleFinSub@qr@LaurentFromSubtraction[spec, fan, {{}},
     "ComplexExponentMode" -> "SplitRealImag", "Integrator" -> "MC", "NSamples" -> NS,
     "EpsilonValues" -> epsVals, "WorkingDirectory" -> wdir["B_splSub"]];

  ok = AllTrue[{dIBP, sIBP, sSub}, ListQ];
  assert["B-routes returned results", ok,
     "DirIBP/SplIBP/SplSub all returned Association[Results]"];
  If[ok, (
    {pdI, psI, psS} = {dIBP[[1]], sIBP[[1]], sSub[[1]]};
    Print["   DirIBP pole=", fmtC[pdI], "  SplIBP pole=", fmtC[psI],
          "  SplSub pole=", fmtC[psS]];
    assert["B1 Direct-IBP pole vs exact Gamma", cplxClose[pdI, poleRef, tolPole],
      "dRe=" <> sci[Re[pdI - poleRef]] <> " dIm=" <> sci[Im[pdI - poleRef]]];
    assert["B2 Split-IBP pole vs exact Gamma", cplxClose[psI, poleRef, tolPole],
      "dRe=" <> sci[Re[psI - poleRef]] <> " dIm=" <> sci[Im[psI - poleRef]]];
    assert["B3 Split-Sub pole vs exact Gamma (pinned-eps, independent algo)",
      cplxClose[psS, poleRef, tolPole],
      "dRe=" <> sci[Re[psS - poleRef]] <> " dIm=" <> sci[Im[psS - poleRef]]];
    assert["B4 SplitRealImag == Direct on IBP path", Abs[psI - pdI] < tolSplit,
      "|dpole|=" <> sci[Abs[psI - pdI]]];
    assert["B5 Split-IBP finite vs exact Gamma", cplxClose[sIBP[[2]], finRef, tolFin],
      "dRe=" <> sci[Re[sIBP[[2]] - finRef]] <> " dIm=" <> sci[Im[sIBP[[2]] - finRef]]];
    assert["B6 pole genuinely complex (Im != 0)", Abs[Im[poleRef]] > 0.01,
      "Im(pole)=" <> sci[Im[poleRef]]];
  )];
];

Print[""];
Print["================================================================"];
If[allPass,
  Print["CC48 PASS  complex exponents x divergence (unlifted): monomial-A and ",
        "polynomial-B, Direct == SplitRealImag == exact Gamma-Laurent, IBP<->Sub agree"],
  Print["CC48 FAIL  one or more assertions failed -- see above"]];
Print["================================================================"];
If[!allPass, Quit[1]];
