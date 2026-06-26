(* L2 inspection on the liftable+divergent toy T:
     I(eps) = Int_0^inf dx1 dx2  x1^{-1+eps} (1 + 10^6 x1 x2 + x1^2 + x2^2)^{-2}
   The added x1^2 lifts the Newton polytope off the (z x1 x2)-plane so the lifted
   3D polytope is full-dimensional.  Pole is unchanged at x1->0 (P->1+x2^2):
     pole = Int_0^inf (1+x2^2)^{-2} dx2 = pi/4. *)
$pkgRoot = "/Users/aidanh/Desktop/TROPICAL_MONTE_CARLOv3";
Get[FileNameJoin[{$pkgRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$pkgRoot, "tropical_eval.wl"}]];

eps  = Symbol["epsL2"];
vars = {x[1], x[2]};
spec = <|
  "Polynomials"         -> {1 + 10^6 x[1] x[2] + x[1]^2 + x[2]^2},
  "MonomialExponents"   -> {-1 + eps, 0},
  "PolynomialExponents" -> {-2},
  "Variables"           -> vars,
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> eps
|>;
rules = {<|"PolyIndex" -> 1, "ExponentVector" -> {1, 1}, "k" -> 1|>};

lift = LiftCoefficients[spec, rules];
liftedSpec = lift["LiftedSpec"]; liftData = lift["LiftData"];
Print["L2: z0 = ", liftData["z0"]];

verts = Quiet[PolytopeVertices[(Times @@ liftedSpec["Polynomials"])^(-1),
                               liftedSpec["Variables"]], TropicalFan::polymake];
fan   = Quiet[computeFanScaled[verts], TropicalFan::polymake];
{dv, sl} = fan;
Print["L2: lifted fan has ", Length[sl], " sectors; simplex sizes = ",
      Length /@ sl];

nDiv = 0; nConv = 0; nEmpty = 0; nFailed = 0; caseAcount = 0; exactOK = True;
Do[
  Module[{sd},
    sd = Quiet@ProcessSectorLifted[liftedSpec, dv, sl[[s]], s, liftData, "Eps" -> eps];
    Which[
      sd === $Failed, nFailed++,
      KeyExistsQ[sd, "EmptyDomain"] && sd["EmptyDomain"], nEmpty++,
      TrueQ[sd["IsDivergent"]],
        nDiv++;
        Module[{k = sd["DivergentVariable"], at = sd["NewExponents"], ck, dc},
          ck = (D[at[[k]], eps] /. eps -> 0);
          dc = sd["DomainConstraint"];
          Print["  DIV sector ", s, ": divVar=", k, "  atilde|eps0=", at /. eps -> 0,
                "  c_k=", ck, "  PrefBase=", sd["PrefactorBase"],
                "  Domain=", If[dc === None, "None(CaseA)",
                  "ic=" <> ToString[dc["IndicatorCoeffs"]] <>
                  " ic_k=" <> ToString[dc["IndicatorCoeffs"][[k]]]]];
          If[dc === None || PossibleZeroQ[dc["IndicatorCoeffs"][[k]]], caseAcount++];
          If[!FreeQ[{at, sd["PrefactorBase"]}, _Real], exactOK = False];
        ],
      True, nConv++;
        If[!FreeQ[{sd["NewExponents"], sd["Prefactor"], sd["PrefactorBase"]}, _Real],
           exactOK = False]
    ]
  ], {s, Length[sl]}];

Print["L2 summary: ", nConv, " conv, ", nDiv, " div, ",
      nEmpty, " empty, ", nFailed, " failed."];
Print["L2: divergent sectors Case A = ", caseAcount, " / ", nDiv];
Print["L2: SectorData exact (FreeQ _Real) = ", exactOK];
Print["L2 GATE: ", If[nDiv >= 1 && nFailed == 0 && exactOK, "PASS", "CHECK"]];
