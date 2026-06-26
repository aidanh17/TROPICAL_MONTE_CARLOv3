(* L3 fast WL-only checks: IBPProcessSector on the lifted divergent sector of
   toy T; verify DLogPrefactor = -ln(10^6), DomainConstraint propagation, and
   ValidateIBP/ValidateSubtraction reconstruct the lifted divergent sector. *)
$pkgRoot = "/Users/aidanh/Desktop/TROPICAL_MONTE_CARLOv3";
Get[FileNameJoin[{$pkgRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$pkgRoot, "tropical_eval.wl"}]];
Print["L3 WL: loaded."];

eps  = Symbol["epsL3"];
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
verts = Quiet[PolytopeVertices[(Times @@ liftedSpec["Polynomials"])^(-1),
                               liftedSpec["Variables"]], TropicalFan::polymake];
fan   = Quiet[computeFanScaled[verts], TropicalFan::polymake];
{dv, sl} = fan;

(* find the divergent lifted sector *)
divSD = {};
Do[Module[{sd = Quiet@ProcessSectorLifted[liftedSpec, dv, sl[[s]], s, liftData, "Eps"->eps]},
  If[AssociationQ[sd] && TrueQ[sd["IsDivergent"]], AppendTo[divSD, sd]]], {s, Length[sl]}];
Print["L3 WL: ", Length[divSD], " divergent lifted sector(s)."];

Do[
  Module[{sd = divSD[[s]], ibp, dlp},
    ibp = IBPProcessSector[sd, liftedSpec];
    If[ibp === $Failed, Print["  IBPProcessSector FAILED"]; Continue[]];
    dlp = ibp["DLogPrefactor"];
    Print["  sector ", sd["ConeIndex"], ": DLogPrefactor = ", dlp,
          "  (expect -Log[10^6] = ", -Log[10^6], " = ", N[-Log[10^6]], ")",
          "  match=", PossibleZeroQ[dlp - (-Log[10^6])]];
    Print["    DomainConstraint propagated = ", ibp["DomainConstraint"]];
    (* ValidateIBP reconstruction at fixed eps *)
    Module[{vr = Quiet@ValidateIBP[ibp, sd, liftedSpec, {}, 0.01]},
      If[AssociationQ[vr],
        Print["    ValidateIBP rel-err = ", vr["RelativeError"]]]];
    (* Subtraction route validator *)
    Module[{pd, vs},
      pd = ProcessDivergentSector[sd, liftedSpec];
      If[pd =!= $Failed,
        vs = Quiet@ValidateSubtraction[pd, sd, liftedSpec, {}, 0.01];
        If[AssociationQ[vs],
          Print["    ValidateSubtraction rel-err = ", vs["RelativeError"]]]]];
  ], {s, Length[divSD]}];
Print["L3 WL DONE."];
