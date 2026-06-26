(* Fast NIntegrate-only scan of lift+divergent candidates across dimensions.
   For each: lifted fan (full-dim?), divergent sector count + HasConstantTerm,
   IBP pole (sum of divergent-sector G0/ck via NIntegrate), and the TRUE pole
   (= Int over x2..xn of prod_j P_j(0,x2..xn)^{B_j}, since the measure is
   x1^{-1+eps}).  match => IBP route valid (HasConstantTerm-clean). *)
$pkgRoot = "/Users/aidanh/Desktop/TROPICAL_MONTE_CARLOv3";
Get[FileNameJoin[{$pkgRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$pkgRoot, "tropical_eval.wl"}]];
SetAttributes[qr, HoldFirst];
qr[e_] := Quiet[e, {General::stop, TropicalFan::polymake, NIntegrate::slwcon,
   NIntegrate::inumr, NIntegrate::ncvb, NIntegrate::eincr}];

(* compactify xi = ui/(1-ui) on [0,1], Jacobian prod 1/(1-ui)^2 *)
poleRef[polys_, polyExps_, dim_] := Module[{us, sub, jac, integrand, p0},
  us = Table[Unique["u"], {dim - 1}];   (* x2..xn *)
  (* set x1=0, map x_{i} (i>=2) -> us[[i-1]]/(1-us[[i-1]]) *)
  sub = Join[{x[1] -> 0},
    Table[x[i] -> us[[i-1]]/(1 - us[[i-1]]), {i, 2, dim}]];
  jac = Times @@ Table[1/(1 - us[[i]])^2, {i, dim - 1}];
  integrand = (Times @@ MapThread[#1^#2 &, {polys, polyExps}] /. sub) * jac;
  qr@NIntegrate[integrand, Evaluate[Sequence @@ ({#, 0, 1} & /@ us)],
    Method -> "GlobalAdaptive", PrecisionGoal -> 5, MaxRecursion -> 50]
];

scan[label_, polys_, polyExps_, lr_] := Module[
  {eps, dim, vars, spec, lift, ls, ld, verts, fan, dv, sl,
   ibpPole = 0, hc = {}, nd = 0, ne = 0, nf = 0, full, pr},
  eps = Symbol["es"]; dim = Length[lr["ExponentVector"]];
  vars = Table[x[i], {i, dim}];
  spec = <|"Polynomials" -> polys,
    "MonomialExponents" -> Join[{-1 + eps}, ConstantArray[0, dim - 1]],
    "PolynomialExponents" -> polyExps, "Variables" -> vars,
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> eps|>;
  lift = qr@LiftCoefficients[spec, {lr}];
  If[!AssociationQ[lift], Print[label, ": LIFT FAILED"]; Return[]];
  ls = lift["LiftedSpec"]; ld = lift["LiftData"];
  verts = qr@PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]];
  fan = qr@computeFanScaled[verts];
  If[!ListQ[fan] || Length[fan] < 2, Print[label, ": FAN FAILED (lower-dim?)"]; Return[]];
  {dv, sl} = fan;
  full = AllTrue[sl, Length[#] == dim + 1 &];
  Do[
    Module[{sd = qr@ProcessSectorLifted[ls, dv, sl[[s]], s, ld, "Eps" -> eps]},
      Which[
        sd === $Failed, nf++,
        KeyExistsQ[sd, "EmptyDomain"] && sd["EmptyDomain"], ne++,
        TrueQ[sd["IsDivergent"]],
          nd++; AppendTo[hc, sd["HasConstantTerm"]];
          Module[{pd = ProcessDivergentSector[sd, ls], g0fp, g0pf, b0, gd, yv, pv, ig, g0, ck},
            g0fp = pd["G0FlatPolys"]; g0pf = pd["G0Prefactor"]; b0 = pd["B0"]; ck = pd["ck"];
            gd = pd["G0Dimension"]; yv = Table[Unique["g"], {gd}];
            pv = Table[Total[Table[m[[1]] Exp[Total[m[[2]] Log /@ yv]], {m, g0fp[[j]]}]], {j, Length[g0fp]}];
            ig = (g0pf /. eps -> 0) Times @@ MapThread[Exp[#2 Log[#1]] &, {pv, b0}];
            g0 = qr@NIntegrate[ig, Evaluate[Sequence @@ ({#, 0, 1} & /@ yv)],
                   Method -> "GlobalAdaptive", PrecisionGoal -> 5, MaxRecursion -> 40];
            ibpPole += g0/ck]]],
    {s, Length[sl]}];
  pr = poleRef[polys, polyExps, dim];
  Print[label];
  Print["   dim=", dim, " fan=", Length[sl], " sec (full-dim=", full, "), ",
        nd, " div, ", ne, " empty, ", nf, " failed;  HasConst(div)=", hc];
  Print["   IBP pole(G0/ck)=", ibpPole, "   truePole(Int P(0)^B)=", pr,
        "   IBP-match=", If[NumericQ[ibpPole] && NumericQ[pr] && Abs[pr] > 10^-9,
          Abs[(ibpPole - pr)/pr] < 0.02, False]];
];

es = Symbol["es"];
(* 2D non-trivial *)
scan["2D-a: 1+x1+x2+1e4 x1x2+x1^2+x2^2, B=-2, lift{1,1}",
  {1 + x[1] + x[2] + 10^4 x[1] x[2] + x[1]^2 + x[2]^2}, {-2}, <|"PolyIndex"->1,"ExponentVector"->{1,1},"k"->1|>];
scan["2D-b: (1+x2)^1 (1+1e5 x1^2+x2^2+x1 x2)^-2, lift{2,0}",
  {1 + x[2], 1 + 10^5 x[1]^2 + x[2]^2 + x[1] x[2]}, {1, -2}, <|"PolyIndex"->2,"ExponentVector"->{2,0},"k"->1|>];
scan["2D-c: 1+1e6 x1^3+x1 x2+x2^2, B=-2, lift{3,0}",
  {1 + 10^6 x[1]^3 + x[1] x[2] + x[2]^2}, {-2}, <|"PolyIndex"->1,"ExponentVector"->{3,0},"k"->1|>];
(* 3D non-trivial *)
scan["3D-a: 1+x1x2+x2x3+1e5 x1^2+x2^2+x3^2, B=-2, lift{2,0,0}",
  {1 + x[1] x[2] + x[2] x[3] + 10^5 x[1]^2 + x[2]^2 + x[3]^2}, {-2}, <|"PolyIndex"->1,"ExponentVector"->{2,0,0},"k"->1|>];
scan["3D-b: 1+x1+x2+x3+1e4 x1x2+x3^2, B=-3, lift{1,1,0}",
  {1 + x[1] + x[2] + x[3] + 10^4 x[1] x[2] + x[3]^2}, {-3}, <|"PolyIndex"->1,"ExponentVector"->{1,1,0},"k"->1|>];
scan["3D-c: (1+x1+x2+x3)^1 (1+1e5 x1^2+x2^2+x3^2)^-3, lift{2,0,0}",
  {1 + x[1] + x[2] + x[3], 1 + 10^5 x[1]^2 + x[2]^2 + x[3]^2}, {1, -3}, <|"PolyIndex"->2,"ExponentVector"->{2,0,0},"k"->1|>];
(* 4D non-trivial *)
scan["4D-a: 1+x1x2+x3x4+1e4 x1^2+x2^2+x3^2+x4^2, B=-3, lift{2,0,0,0}",
  {1 + x[1] x[2] + x[3] x[4] + 10^4 x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2}, {-3}, <|"PolyIndex"->1,"ExponentVector"->{2,0,0,0},"k"->1|>];
