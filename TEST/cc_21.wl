(* ============================================================================
   TEST/cc_21.wl  —  Cross-check #21 (plan.md §8.2)
   Complex-flattening variance reduction (contour rotation).  Tier 1: WL only.

   Integral:  I = Int_0^inf (1+x)^B dx,   B = -2 + I*b,   b = 6
   Exact result:  1 / (1 - I*b)

   The large-x sector (ray +1) carries a complex effective exponent
       a = 1 - I*b
   whose imaginary part causes the real-axis integrand g(y) to oscillate
   (spiralling contour image).  Flattening y -> (y')^(1/a) rotates the
   contour back to a smooth path; the plan calls for a ~20x variance drop.

   PASS criteria (both must hold):
     P1  |ItotFlat - exact| / |exact| < 1e-4          (value agrees on both contours)
     P2  min(VarRed_Re, VarRed_Im) >= 10              (~20x variance reduction, >=10x gate)

   Ported from:
     OLD_CODE/TROPICAL_MONTE_CARLO2/TROPICAL_MONTE_CARLOv2/SANDBOX/verify_complex_flatten.wl
   with path fixes, absolute-path Get, and CC21 PASS/FAIL print convention.

   Run:
     wolframscript -file TEST/cc_21.wl
   from the TROPICAL_MONTE_CARLO3 directory, or with $InputFileName set.
   ============================================================================ *)

(* --- locate root and load v3 packages via absolute paths --- *)
With[{root = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]]},
  Get[FileNameJoin[{root, "tropical_fan.wl"}]];
  Get[FileNameJoin[{root, "tropical_eval.wl"}]];
];

(* ------------------------------------------------------------------ *)
(* Parameters                                                           *)
(* ------------------------------------------------------------------ *)
bval  = 6;
Bc    = -2 + I bval;
exact = 1 / (1 - I bval);

(* Tolerance gates *)
tol1       = 1*^-4;    (* relative error of the sector sum vs exact *)
varRedGate = 10;       (* minimum variance-reduction factor (both Re and Im) *)

(* MC sample count for the variance comparison *)
nMC = 200000;

Print["================================================================"];
Print["  CC21: complex-flattening contour-rotation check"];
Print["  I = Int_0^inf (1+x)^(", Bc, ") dx"];
Print["  Exact = 1/(1-I*", bval, ") = ", N[exact]];
Print["================================================================"];

(* ------------------------------------------------------------------ *)
(* 1D fan: rays +1 (large-x sector) and -1 (small-x sector)            *)
(*   fan = { rayMatrix, { {sectorIndex,...}, ... } }                    *)
(* ------------------------------------------------------------------ *)
fanRays    = {{1}, {-1}};
fanSectors = {{1}, {2}};
fan        = {fanRays, fanSectors};

spec = <|
  "Polynomials"         -> {1 + x[1]},
  "MonomialExponents"   -> {0},
  "PolynomialExponents" -> {Bc},
  "Variables"           -> {x[1]},
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> None
|>;

(* Process both sectors *)
sd = Table[
  ProcessSector[spec, fanRays, {fanSectors[[s, 1]]}, s],
  {s, 2}
];

Do[
  Print["--- Sector ", s, "  ray = ", fanRays[[fanSectors[[s, 1]]]], " ---"];
  Print["   NewExponents (a) = ", sd[[s]]["NewExponents"]];
  Print["   Prefactor        = ", N[sd[[s]]["Prefactor"]]];
  Print["   FlattenedPolys   = ", sd[[s]]["FlattenedPolys"]];
  Print["   PolyExps (B)     = ", sd[[s]]["PolynomialExponents"]];,
  {s, 2}
];

(* ------------------------------------------------------------------ *)
(* 2. Focus on the complex (large-x) sector: sector 1, ray +1          *)
(* ------------------------------------------------------------------ *)
sCx = 1;
aP  = sd[[sCx]]["NewExponents"][[1]];    (* complex effective exponent *)
pf  = sd[[sCx]]["Prefactor"];

Print[""];
Print["Complex sector: a = ", aP, "   Im(a) = ", Im[aP]];

Print[""];
Print["Contour image y = (y')^(1/a) for y' = 1, .5, .1, .01, .001:"];
Do[
  Print["   y' = ", yp, "  ->  y = ", N[yp^(1/aP)],
        "   (|y|=", N[Abs[yp^(1/aP)]], ", arg=", N[Arg[yp^(1/aP)]], ")"],
  {yp, {1, 0.5, 0.1, 0.01, 0.001}}
];

(* ------------------------------------------------------------------ *)
(* 3. Integrands                                                         *)
(*   g(y) — real-axis (un-rotated) integrand for sector +               *)
(*   f(y') — flattened (rotated contour) integrand for sector +         *)
(* ------------------------------------------------------------------ *)
g[y_]  := y^(aP - 1) * (1 + y)^Bc;
f[yp_] := pf * (1 + yp^(1/aP))^Bc;

(* Other sector (ray -1): flattened integrand *)
aOther    = sd[[2]]["NewExponents"][[1]];
pfOther   = sd[[2]]["Prefactor"];
gOther[y_] := pfOther * (1 + y^(1/aOther))^Bc;

(* ------------------------------------------------------------------ *)
(* 4. Numeric integrals (both contours must give the same answer)       *)
(* ------------------------------------------------------------------ *)
IgReal = NIntegrate[g[y], {y, 0, 1}, MaxRecursion -> 40, WorkingPrecision -> 20];
IfFlat = NIntegrate[f[yp], {yp, 0, 1}, MaxRecursion -> 40, WorkingPrecision -> 20];

ItotFlat = IfFlat + NIntegrate[
  gOther[y], {y, 0, 1}, MaxRecursion -> 40, WorkingPrecision -> 20
];

Print[""];
Print["Sector+ value, real contour     Int_0^1 g  = ", N[IgReal]];
Print["Sector+ value, rotated contour  Int_0^1 f  = ", N[IfFlat]];
Print["Total (both flattened sectors)             = ", N[ItotFlat]];
Print["Exact 1/(1 - I*", bval, ")                = ", N[exact]];

relErr = N[Abs[(ItotFlat - exact) / exact]];
Print["Relative error                             = ", relErr];

(* ------------------------------------------------------------------ *)
(* 5. Variance comparison: real contour g vs flattened f                *)
(*    (uniform MC over 2e5 points)                                      *)
(* ------------------------------------------------------------------ *)
SeedRandom[42];
ys = RandomReal[{0, 1}, nMC];
gv = g /@ ys;
fv = f /@ ys;

varReG = Variance[Re[gv]];
varImG = Variance[Im[gv]];
varReF = Variance[Re[fv]];
varImF = Variance[Im[fv]];

varRedRe = N[varReG / varReF];
varRedIm = N[varImG / varImF];

Print[""];
Print["Per-sample statistics over ", nMC, " uniform points (sector +):"];
Print["   real contour  g : mean=", N[Mean[gv]],
      "  Var[Re]=", N[varReG], "  Var[Im]=", N[varImG]];
Print["   rotated (flat) f: mean=", N[Mean[fv]],
      "  Var[Re]=", N[varReF], "  Var[Im]=", N[varImF]];
Print["   |g| range = ", {N@Min@Abs@gv, N@Max@Abs@gv},
      "    |f| range = ", {N@Min@Abs@fv, N@Max@Abs@fv}];
Print["   Var reduction Re: ", varRedRe, "   Im: ", varRedIm];

(* ------------------------------------------------------------------ *)
(* 6. PASS / FAIL                                                        *)
(* ------------------------------------------------------------------ *)
pass1 = (relErr < tol1);
pass2 = (Min[varRedRe, varRedIm] >= varRedGate);

If[pass1,
  Print["CC21 PASS P1: total vs exact  relErr=", relErr, " < ", tol1],
  Print["CC21 FAIL P1: expected relErr<", tol1, " got=", relErr,
        " expected=", N[exact], " got=", N[ItotFlat]]
];

If[pass2,
  Print["CC21 PASS P2: variance reduction Re=", varRedRe,
        " Im=", varRedIm, " both >= ", varRedGate],
  Print["CC21 FAIL P2: expected VarRed>=", varRedGate,
        " got Re=", varRedRe, " Im=", varRedIm]
];

If[pass1 && pass2,
  Print["CC21 PASS: complex-flattening contour-rotation check OK"],
  Print["CC21 FAIL: one or more criteria not met (see above)"]
];
