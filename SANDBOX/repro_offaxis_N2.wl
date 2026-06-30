(* N=2 validation: ONE real pole (x1: 2eps) + TWO off-axis (x2: i th2, x3: i th3)
   in the same corner cone.  Dirichlet oracle:
     I(eps) = Gamma(2eps) Gamma(i th2) Gamma(i th3) Gamma(B-2eps-i th2-i th3)/Gamma(B)
   Confirms the generalized corner assembly is correct for N=2 (planIBPMULTIDIV
   §8.2/§8.5: validate N>=2 before trusting it). *)

$pkgRoot = "/Users/aidanh/Desktop/TROPICAL_MONTE_CARLOv3";
Get[FileNameJoin[{$pkgRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$pkgRoot, "tropical_eval.wl"}]];
Print["loaded."];

th2 = 1/2; th3 = 1/3; B = 4;
eps = Symbol["epsN2"];
polys = {1 + x[1] + x[2] + x[3]};
pe    = {-B};
Avals = {2 eps - 1, I th2 - 1, I th3 - 1};

spec = <|"Polynomials" -> polys, "MonomialExponents" -> Avals,
  "PolynomialExponents" -> pe, "Variables" -> {x[1], x[2], x[3]},
  "KinematicSymbols" -> {}, "RegulatorSymbol" -> eps|>;
fan = computeFanScaled[PolytopeVertices[(Times @@ polys)^(-1), {x[1], x[2], x[3]}]];
wd = FileNameJoin[{$pkgRoot, "TEST", "INTERFILES", "repro_N2"}];
Quiet[CreateDirectory[wd, CreateIntermediateDirectories -> True], {CreateDirectory::eexist}];

orc[e_] := Gamma[2 e] Gamma[I th2] Gamma[I th3] Gamma[B - 2 e - I th2 - I th3]/Gamma[B];
ser = Series[orc[ee], {ee, 0, 0}];
Print["oracle pole   = ", N[SeriesCoefficient[ser, -1]]];
Print["oracle finite = ", N[SeriesCoefficient[ser, 0]]];

r = EvaluateTropicalMC[spec, fan, {{}}, "Method" -> "IBP",
  "ComplexExponentMode" -> "Direct", "Integrator" -> "MC", "NSamples" -> 4000000,
  "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> wd];
If[AssociationQ[r] && KeyExistsQ[r, "Results"],
  Print["engine pole   = ", N[r["Results"][[1]]["PoleCoefficient"]]];
  Print["engine finite = ", N[r["Results"][[1]]["FinitePart"]]];
  Print["MultiDiv sectors (HasPole, OffAxisDirs): ",
    ({#["HasPole"], #["OffAxisDirs"]} & /@
       Select[r["IBPProcessedSectors"], TrueQ[#["MultiDiv"]] &])],
  Print["engine head = ", Head[r]]];
