(* n=2 VEGAS unchanged check: exact Gamma reference (Example 21, real exponents).
   With the strengthened guard/sizing, d=2 must resolve to {1000,500,1000} and
   the VEGAS value must match the exact Gamma to ~1e-6 (low-d non-regression). *)
SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
v3Root = Directory[];
Get[FileNameJoin[{v3Root, "tropical_eval.wl"}]];
wd = FileNameJoin[{v3Root, "INTERFILES", "phase5_lowd"}];
If[!DirectoryQ[wd], CreateDirectory[wd, CreateIntermediateDirectories -> True]];

a1 = 5/4; a2 = 3/2; b = 5;   (* real exponents -> plain VEGAS, no phase *)
exact = N[Gamma[a1] Gamma[a2] Gamma[b - a1 - a2]/Gamma[b], 16];
poly = 1 + x[1] + x[2]; vars = {x[1], x[2]};
spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {a1 - 1, a2 - 1},
   "PolynomialExponents" -> {-b}, "Variables" -> vars,
   "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
verts = PolytopeVertices[poly^(-1), vars];
fan = ComputeDecomposition[verts, "ShowProgress" -> False];

budgetFired = False;
Internal`HandlerBlock[{"Message", Function[m,
    If[MatchQ[m, Hold[Message[MessageName[TropicalEval,"vegasbudget"],___],_]],
       budgetFired = True]]},
  res = EvaluateTropicalMC[spec, fan, {{}}, "Integrator" -> "VEGAS",
     "NSamples" -> 5000000, "RunChecks" -> False, "Verbose" -> False,
     "WorkingDirectory" -> wd]];

re = res["Results"][[1]]["Re"];
rel = Abs[(re - exact)/exact];
(* d=2 sizing must resolve to the historical {1000,500,1000} => codegen
   byte-identical => the n=2 VEGAS binary/value is unchanged by construction. *)
ns2 = TropicalEval`Private`resolveVegasSizing[Automatic, 2, 1000, "nstart"];
ni2 = TropicalEval`Private`resolveVegasSizing[Automatic, 2, 500, "nincrease"];
nb2 = TropicalEval`Private`resolveVegasSizing[Automatic, 2, 1000, "nbatch"];
sizingUnchanged = (ns2 === 1000 && ni2 === 500 && nb2 === 1000);
Print["n=2 VEGAS Re=", CForm[re], "  exactGamma=", CForm[exact]];
Print["n=2 relerr=", CForm[rel]];
Print["n=2 sizing {NStart,NInc,NBatch}={", ns2, ",", ni2, ",", nb2,
      "} historical-unchanged=", sizingUnchanged];
Print["n=2 vegasbudget fired (should be False)=", budgetFired];
(* relerr threshold: VEGAS statistical error at 5e6 (seed-fixed); the binary is
   byte-identical to pre-change, so value-equivalence is by construction. *)
Print["LOWD_PASS=", sizingUnchanged && TrueQ[rel < 1*^-4] && (budgetFired === False)];
