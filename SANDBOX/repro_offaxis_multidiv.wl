(* Minimal reproduction: one REAL pole (x1: 2eps-1) + one OFF-AXIS direction
   (x2: I*theta-1) in the SAME cone, UNLIFTED.  Closed-form Dirichlet oracle:
     I(eps) = Gamma(2eps) Gamma(I th) Gamma(B-2eps-I th)/Gamma(B)
   => simple 1/eps pole (from Gamma(2eps)); Gamma(I th) is a finite complex const.
   So there is exactly ONE true pole + one off-axis (finite) direction.
   Question: does the engine refuse today, and via WHICH guard? *)

$pkgRoot = "/Users/aidanh/Desktop/TROPICAL_MONTE_CARLOv3";
Get[FileNameJoin[{$pkgRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$pkgRoot, "tropical_eval.wl"}]];
Print["loaded."];

th = 1/2; B = 3;
eps = Symbol["epsR"];
polys = {1 + x[1] + x[2]};
pe    = {-B};
Avals = {2 eps - 1, I th - 1};   (* x1 real-soft pole; x2 off-axis soft *)

spec = <|"Polynomials" -> polys, "MonomialExponents" -> Avals,
  "PolynomialExponents" -> pe, "Variables" -> {x[1], x[2]},
  "KinematicSymbols" -> {}, "RegulatorSymbol" -> eps|>;

fan = computeFanScaled[PolytopeVertices[(Times @@ polys)^(-1), {x[1], x[2]}]];
Print["fan cones = ", Length[fan]];

wd = FileNameJoin[{$pkgRoot, "TEST", "INTERFILES", "repro_offaxis"}];
Quiet[CreateDirectory[wd, CreateIntermediateDirectories -> True], {CreateDirectory::eexist}];

Print["==== Method -> IBP, Direct ===="];
rD = EvaluateTropicalMC[spec, fan, {{}}, "Method" -> "IBP",
  "ComplexExponentMode" -> "Direct", "Integrator" -> "MC", "NSamples" -> 200000,
  "RunChecks" -> False, "Verbose" -> True, "WorkingDirectory" -> wd];
Print["Direct head: ", Head[rD]];
If[AssociationQ[rD] && KeyExistsQ[rD, "Results"],
  Print["Direct pole=", rD["Results"][[1]]["PoleCoefficient"],
        " finite=", rD["Results"][[1]]["FinitePart"]]];

Print["==== Method -> IBP, SplitRealImag ===="];
rS = EvaluateTropicalMC[spec, fan, {{}}, "Method" -> "IBP",
  "ComplexExponentMode" -> "SplitRealImag", "Integrator" -> "MC", "NSamples" -> 200000,
  "RunChecks" -> False, "Verbose" -> True, "WorkingDirectory" -> wd];
Print["Split head: ", Head[rS]];

(* oracle *)
orc[e_] := Gamma[2 e] Gamma[I th] Gamma[B - 2 e - I th]/Gamma[B];
ser = Series[orc[ee], {ee, 0, 0}];
Print["oracle pole = ", N[SeriesCoefficient[ser, -1]]];
Print["oracle finite = ", N[SeriesCoefficient[ser, 0]]];
