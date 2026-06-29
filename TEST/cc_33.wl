(* ============================================================================
   TEST/cc_33.wl  --  v3 Cross-Check #33
   plan.md §8.3 row #33:  Divergence x complex exponents, full matrix.
   F2a x F4 x F5 -- IBP (default) and subtraction both handle a 2D integral
   with a single 1/eps pole and complex polynomial exponent B.

   Tier 1 (no external deps):  {IBP, Subtraction} x MC;
                                pole vs exact Gamma-Laurent; IBP vs Sub agree.
   Tier 2 (needs CUBA):        add VEGAS columns; MC vs VEGAS agree.
   Degrades cleanly: if CUBA absent, Tier-2 assertions are skipped (not failed).

   Integrand:
       I(eps) = Int_0^inf x1^{2*eps-1} (1+x1+x2)^{-B} dx1 dx2,
   with B = 2 + I/2  (complex; Im(B) != 0 so arg(P) != 0 in general sectors).

   Exact Laurent (Mellin + Beta function):
     Step 1:  Int_0^inf (1+x1+x2)^{-B} dx2 = (1+x1)^{1-B}/(B-1)  [Re(B)>1]
     Step 2:  I(eps) = 1/(B-1) * Int_0^inf x1^{2eps-1} (1+x1)^{1-B} dx1
                     = 1/(B-1) * Beta(2 eps, B-1-2 eps)
                     = Gamma(2 eps) Gamma(B-1-2 eps) / ((B-1) Gamma(B-1))
     Pole:    lim_{eps->0} eps * I(eps)
              = [lim eps*Gamma(2 eps)] * Gamma(B-1) / ((B-1) Gamma(B-1))
              = (1/2) / (B-1)    [since eps*Gamma(2eps) -> 1/2]
              = 1/(2(B-1))
     Finite:  subleading: d/deps [eps * I(eps)] |_{eps=0}
              ~ (-EulerGamma - PolyGamma[0, B-1]) / (B-1)
              Derivation: Gamma(2 eps) = 1/(2 eps) - EulerGamma + O(eps);
              Gamma(B-1-2 eps) = Gamma(B-1)[1 - 2 eps*PolyGamma[0,B-1] + O(eps^2)].
              Product / ((B-1)*Gamma(B-1)) gives finite part (-EulerGamma - PolyGamma[0,B-1])/(B-1).

   PASS criteria:
     A1  |pole_IBP_MC  - poleRef| < 5e-3  (real and imag parts)
     A2  |pole_Sub_MC  - poleRef| < 5e-3
     A3  |pole_IBP_MC  - pole_Sub_MC| < 1e-2   (IBP vs Sub pole agree)
     A4  |fin_IBP_MC   - fin_Sub_MC|  < 2e-2   (IBP vs Sub finite agree)
     A5  (CUBA) |pole_IBP_VEGAS - pole_IBP_MC| < 1e-2
     A6  (CUBA) |pole_Sub_VEGAS - pole_Sub_MC| < 1e-2

   Reference/mirror:
     OLD_CODE Tree-A EXAMPLES/tropical_eval_examples.wl Test 18 (divergent +
     complex exponents, ValidateSubtraction) and plan.md §8.4 Ex-V3a.
   ============================================================================ *)

(* --------------------------------------------------------------------------
   Load v3 packages by absolute path so this file runs from any cwd.
   -------------------------------------------------------------------------- *)
$v3Root = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{$v3Root, "tropical_fan.wl"}]];
Get[FileNameJoin[{$v3Root, "tropical_eval.wl"}]];

(* --------------------------------------------------------------------------
   Helpers
   -------------------------------------------------------------------------- *)
SetAttributes[quietRun, HoldFirst];
quietRun[expr_] := Block[{$Output = {}}, expr];

fmtC[z_] := Module[{r = Re[z], im = Im[z]},
  If[Abs[im] < 1.*^-12,
    ToString[NumberForm[N[r], {4, 4}]],
    ToString[NumberForm[N[r], {4, 4}]] <> "+" <>
    ToString[NumberForm[N[im], {4, 4}]] <> "I"]];

sci[x_] := Module[{a = N[Abs[x]], e},
  If[a == 0., "0",
    e = Floor[Log10[a]];
    ToString[NumberForm[a/10^e, {3, 2}]] <> "e" <> ToString[e]]];

(* K-scaled fan -- mirrors the computeFanRobust helper used elsewhere in v3 *)
computeFanRobust33[verts_] := Module[{n, fd},
  n = Length[First[verts]];
  Do[fd = Quiet@ComputeDecomposition[K*verts, "ShowProgress" -> False];
     If[ListQ[fd] && Length[fd] == 2 && FreeQ[fd, $Failed] && Length[fd[[1]]] > 0,
        Return[fd, Module]],
     {K, {1, n + 2, 2 n + 4, 6 n + 6}}];
  $Failed];

(* --------------------------------------------------------------------------
   Detect CUBA (memoised in detectCuba[])
   -------------------------------------------------------------------------- *)
hasCuba = TrueQ[detectCuba[]["Found"]];

Print[""];
Print["================================================================"];
Print["  v3 Cross-Check #33: Divergence x Complex Exponents (full matrix)"];
Print["  plan.md §8.3 row #33 -- F2a x F4 x F5"];
Print["  Tier 1 (MC), Tier 2 (VEGAS, CUBA=" <> ToString[hasCuba] <> ")"];
Print["================================================================"];

(* --------------------------------------------------------------------------
   Problem specification
   Integrand: x1^{2*eps-1} * (1 + x1 + x2)^{-B},  B = 2 + I/2
   -------------------------------------------------------------------------- *)
eps33  = Symbol["eps33"];
Bval   = 2 + 1/2 * I;
vars33 = {x[1], x[2]};
poly33 = 1 + x[1] + x[2];

(* Exact Laurent coefficients *)
poleRef = N[1 / (2 (Bval - 1))];
finRef  = N[(-EulerGamma - PolyGamma[0, Bval - 1]) / (Bval - 1)];

Print["  B       = ", Bval];
Print["  poleRef = 1/(2(B-1)) = ", poleRef];
Print["  finRef  = ", finRef, "  (subleading Laurent coeff)"];
Print[""];

spec33 = <|
  "Polynomials"        -> {poly33},
  "MonomialExponents"  -> {2 eps33 - 1, 0},
  "PolynomialExponents"-> {-Bval},
  "Variables"          -> vars33,
  "KinematicSymbols"   -> {},
  "RegulatorSymbol"    -> eps33
|>;

(* Build fan *)
verts33 = PolytopeVertices[poly33^(-2), vars33];
fan33   = computeFanRobust33[verts33];
If[fan33 === $Failed,
  Print["CC33 FAIL expected=fan got=$Failed"];
  Quit[1]];
Print["  Fan: ", Length[fan33[[2]]], " sectors"];
Print[""];

allPass = True;
tolPoleRef = 5.*^-3;
tolAgreeP  = 1.*^-2;
tolAgreeF  = 2.*^-2;
tolVegas   = 1.*^-2;

(* --------------------------------------------------------------------------
   Tier 1: IBP + MC
   -------------------------------------------------------------------------- *)
Print["  [Tier 1] Running IBP + MC ..."];
ibpMC = quietRun @ EvaluateTropicalMCIBP[spec33, fan33, {{}},
  "Integrator" -> "MC", "NSamples" -> 2000000,
  "RunChecks" -> False, "Verbose" -> False];
If[!AssociationQ[ibpMC] || !KeyExistsQ[ibpMC, "Results"],
  Print["CC33 FAIL expected=IBP-MC-result got=$Failed"];
  Quit[1]];
poleIbpMC = ibpMC["Results"][[1]]["PoleCoefficient"];
finIbpMC  = ibpMC["Results"][[1]]["FinitePart"];
Print["  IBP-MC:   pole=", fmtC[poleIbpMC], "  finite=", fmtC[finIbpMC]];

(* --------------------------------------------------------------------------
   Tier 1: Subtraction + MC
   -------------------------------------------------------------------------- *)
Print["  [Tier 1] Running Subtraction + MC ..."];
epsVals33 = {0.01, 0.02, 0.04};
subMC = quietRun @ LaurentFromSubtraction[spec33, fan33, {{}},
  "Integrator"    -> "MC",
  "NSamples"      -> 2500000,
  "EpsilonValues" -> epsVals33];
If[!AssociationQ[subMC] || !KeyExistsQ[subMC, "Results"],
  Print["CC33 FAIL expected=Sub-MC-result got=$Failed"];
  Quit[1]];
poleSubMC = subMC["Results"][[1]]["Pole"];
finSubMC  = subMC["Results"][[1]]["Finite"];
Print["  Sub-MC:   pole=", fmtC[poleSubMC], "  finite=", fmtC[finSubMC]];
Print[""];

(* --------------------------------------------------------------------------
   Tier-1 assertions (A1-A4, no CUBA needed)
   -------------------------------------------------------------------------- *)

(* A1: IBP-MC pole vs exact Gamma-Laurent *)
Module[{eRe = Abs[Re[poleIbpMC] - Re[poleRef]],
         eIm = Abs[Im[poleIbpMC] - Im[poleRef]], p},
  p = eRe < tolPoleRef && eIm < tolPoleRef;
  If[p,
    Print["CC33 PASS A1 IBP-MC pole vs exact (Re err=", sci[eRe],
          " Im err=", sci[eIm], ")"],
    Print["CC33 FAIL A1 IBP-MC pole vs exact expected=<", tolPoleRef,
          " got Re err=", sci[eRe], " Im err=", sci[eIm]]];
  allPass = allPass && p];

(* A2: Sub-MC pole vs exact *)
Module[{eRe = Abs[Re[poleSubMC] - Re[poleRef]],
         eIm = Abs[Im[poleSubMC] - Im[poleRef]], p},
  p = eRe < tolPoleRef && eIm < tolPoleRef;
  If[p,
    Print["CC33 PASS A2 Sub-MC pole vs exact (Re err=", sci[eRe],
          " Im err=", sci[eIm], ")"],
    Print["CC33 FAIL A2 Sub-MC pole vs exact expected=<", tolPoleRef,
          " got Re err=", sci[eRe], " Im err=", sci[eIm]]];
  allPass = allPass && p];

(* A3: IBP vs Sub pole agreement -- headline assertion of #33 *)
Module[{dp = Abs[poleIbpMC - poleSubMC], p},
  p = dp < tolAgreeP;
  If[p,
    Print["CC33 PASS A3 IBP vs Sub pole agree (|dpole|=", sci[dp], ")"],
    Print["CC33 FAIL A3 IBP vs Sub pole agree expected=<", tolAgreeP,
          " got=", sci[dp]]];
  allPass = allPass && p];

(* A4: IBP vs Sub finite agreement *)
Module[{df = Abs[finIbpMC - finSubMC], p},
  p = df < tolAgreeF;
  If[p,
    Print["CC33 PASS A4 IBP vs Sub finite agree (|dfin|=", sci[df], ")"],
    Print["CC33 FAIL A4 IBP vs Sub finite agree expected=<", tolAgreeF,
          " got=", sci[df]]];
  allPass = allPass && p];

(* --------------------------------------------------------------------------
   Tier 2: VEGAS (CUBA), skip gracefully if absent
   -------------------------------------------------------------------------- *)
If[!hasCuba,
  Print[""];
  Print["  [Tier 2] CC33 SKIP (no CUBA) -- A5/A6 not evaluated"],
  (* else: CUBA present -- run VEGAS columns *)
  Print[""];
  Print["  [Tier 2] Running IBP + VEGAS ..."];
  ibpVeg = quietRun @ EvaluateTropicalMCIBP[spec33, fan33, {{}},
    "Integrator" -> "VEGAS", "NSamples" -> 2000000,
    "RunChecks" -> False, "Verbose" -> False];
  If[!AssociationQ[ibpVeg] || !KeyExistsQ[ibpVeg, "Results"],
    Print["CC33 FAIL A5 expected=IBP-VEGAS-result got=$Failed"];
    allPass = False,
    (* else *)
    poleIbpVeg = ibpVeg["Results"][[1]]["PoleCoefficient"];
    finIbpVeg  = ibpVeg["Results"][[1]]["FinitePart"];
    Print["  IBP-VEGAS: pole=", fmtC[poleIbpVeg], "  finite=", fmtC[finIbpVeg]];
    Module[{dp = Abs[poleIbpMC - poleIbpVeg],
             df = Abs[finIbpMC  - finIbpVeg], p},
      p = dp < tolVegas && df < tolAgreeF;
      If[p,
        Print["CC33 PASS A5 IBP: MC vs VEGAS (|dpole|=", sci[dp],
              " |dfin|=", sci[df], ")"],
        Print["CC33 FAIL A5 IBP: MC vs VEGAS expected pole<", tolVegas,
              " fin<", tolAgreeF, " got |dpole|=", sci[dp],
              " |dfin|=", sci[df]]];
      allPass = allPass && p]];

  Print["  [Tier 2] Running Subtraction + VEGAS ..."];
  subVeg = quietRun @ LaurentFromSubtraction[spec33, fan33, {{}},
    "Integrator"    -> "VEGAS",
    "NSamples"      -> 2500000,
    "EpsilonValues" -> epsVals33];
  If[!AssociationQ[subVeg] || !KeyExistsQ[subVeg, "Results"],
    Print["CC33 FAIL A6 expected=Sub-VEGAS-result got=$Failed"];
    allPass = False,
    (* else *)
    poleSubVeg = subVeg["Results"][[1]]["Pole"];
    finSubVeg  = subVeg["Results"][[1]]["Finite"];
    Print["  Sub-VEGAS: pole=", fmtC[poleSubVeg], "  finite=", fmtC[finSubVeg]];
    Module[{dp = Abs[poleSubMC - poleSubVeg],
             df = Abs[finSubMC  - finSubVeg], p},
      p = dp < tolVegas && df < tolAgreeF;
      If[p,
        Print["CC33 PASS A6 Sub: MC vs VEGAS (|dpole|=", sci[dp],
              " |dfin|=", sci[df], ")"],
        Print["CC33 FAIL A6 Sub: MC vs VEGAS expected pole<", tolVegas,
              " fin<", tolAgreeF, " got |dpole|=", sci[dp],
              " |dfin|=", sci[df]]];
      allPass = allPass && p]]];

(* --------------------------------------------------------------------------
   Summary
   -------------------------------------------------------------------------- *)
Print[""];
Print["================================================================"];
If[allPass,
  Print["CC33 PASS  divergence x complex exponents full matrix (B=", Bval,
        ")  {IBP,Sub}x{MC" <> If[hasCuba, ",VEGAS", ""] <>
        "} pairwise agree and match exact Gamma-Laurent"],
  Print["CC33 FAIL  one or more assertions failed -- see above"]];
Print["================================================================"];

If[!allPass, Quit[1]];
