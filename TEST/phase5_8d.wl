(* ============================================================================
   TEST/phase5_8d.wl — 8D lifted+complex+VEGAS, DEFAULT sizing (Issue 1).
   P = 1 + c x1^2 + x1 + Sum_{i>=2} xi^2,  c=(1+I)1e-6,  B=-6 (PolynomialExponents).
   Independent truth (Schwinger WP=40 + CUBA Cuhre/Vegas 1e9): Re=0.00317086.
   Env: NSAMPLES (per-sector maxeval); VNSTART (explicit VegasNStart, optional).
   ============================================================================ *)
SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
v3Root = Directory[];
Get[FileNameJoin[{v3Root, "tropical_eval.wl"}]];

nsamp = ToExpression[Environment["NSAMPLES"]];
If[!IntegerQ[nsamp], nsamp = 4000000];
vnstart = ToExpression[Environment["VNSTART"]];   (* Null or integer *)
wd = FileNameJoin[{v3Root, "INTERFILES",
   "phase5_8d_" <> ToString[nsamp] <> "_" <> ToString[vnstart]}];
If[!DirectoryQ[wd], CreateDirectory[wd, CreateIntermediateDirectories -> True]];

oracleRe = 0.00317086`;

n = 8; vars = Table[x[i], {i, n}]; cc = (1 + I) 10^-6;
poly = 1 + cc x[1]^2 + x[1] + Sum[x[i]^2, {i, 2, n}];
spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> ConstantArray[0, n],
   "PolynomialExponents" -> {-6}, "Variables" -> vars,
   "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;

Print["=== 8D lifted+complex+VEGAS NSamples=", nsamp, " VNStart=", vnstart, " ==="];

(* Optional cached lifted fan to avoid the ~260s rebuild across sizing points. *)
fanCacheFile = FileNameJoin[{v3Root, "INTERFILES", "phase5_lifted_fan.mx"}];
liftedFan = Automatic;
If[FileExistsQ[fanCacheFile],
  liftedFan = Import[fanCacheFile];
  Print["loaded cached lifted fan: ", Length[liftedFan[[2]]], " cones"],
  Module[{lr, ls, lf, lv},
    lr = DetectExtremeCoefficients[spec, 1000];
    lr = LiftCoefficients[spec, lr]; ls = lr["LiftedSpec"];
    lv = Quiet[PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]],
               TropicalFan::polymake];
    lf = Quiet[TropicalFan`computeFanScaled[lv], TropicalFan::polymake];
    If[ListQ[lf] && Length[lf] >= 2, liftedFan = lf; Export[fanCacheFile, lf];
       Print["built + cached lifted fan: ", Length[lf[[2]]], " cones"]]
  ]
];

(* Capture vegasbudget firing via a transparent message hook. *)
budgetFired = False;
Internal`HandlerBlock[{"Message", Function[m,
    If[MatchQ[m, Hold[Message[MessageName[TropicalEval, "vegasbudget"], ___], _]],
       budgetFired = True]]},
  opts = {"Integrator" -> "VEGAS", "ComplexExponentMode" -> "SplitRealImag",
     "NSamples" -> nsamp, "RunChecks" -> False, "Verbose" -> False,
     "WorkingDirectory" -> wd, "FanData" -> liftedFan};
  If[IntegerQ[vnstart], opts = Append[opts, "VegasNStart" -> vnstart]];
  t0 = AbsoluteTime[];
  res = EvaluateTropicalMCLifted[spec, {{}}, Sequence @@ opts];
  el = Round[AbsoluteTime[] - t0];
];

Print["elapsed=", el, "s"];
Print["VEGASBUDGET_FIRED=", budgetFired];
If[AssociationQ[res],
  re = res["Results"][[1]]["Re"]; im = res["Results"][[1]]["Im"];
  reErr = res["Results"][[1]]["ReErr"];
  Print["VALUE Re=", CForm[re], "  Im=", CForm[im]];
  Print["VEGAS_ReErr=", CForm[reErr]];
  Print["RELERR_vs_oracle=", CForm[Abs[(re - oracleRe)/oracleRe]]];
  Print["DELTA_Re=", CForm[re - oracleRe], " (negative => biased LOW)"];
  Print["SIGMA_vs_VEGASbar=", CForm[If[reErr > 0, Abs[(re - oracleRe)/reErr], -1.]]],
  Print["RESULT=$Failed"]
];
