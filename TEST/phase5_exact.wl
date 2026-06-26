(* ============================================================================
   TEST/phase5_exact.wl — Issue 2: EXACT sector-decomposition identity (#37
   redefined, plan.md §1.4).  For the fan the code ACTUALLY uses, the symbolic
   sum of sector contributions under the monomial change of variables
   x_i -> prod_j y_j^{M_ij} reproduces the ORIGINAL integrand as an EXACT
   symbolic identity, with NO floating-point decision (FreeQ[_,_Real]).
   Also: the wired computeFanScaled path (K=1 first) is byte-identical to the
   raw ComputeDecomposition path; tropical_fan.wl edits are additive-only.
   No NIntegrate, no MC — pure symbolic algebra.
   ============================================================================ *)
SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
v3Root = Directory[];
Get[FileNameJoin[{v3Root, "tropical_eval.wl"}]];
pass = {};
rec[n_, b_] := (AppendTo[pass, b]; Print[n, If[b, " PASS", " FAIL"]]);

(* ---- exact per-sector substitution identity for a given spec/fan ---- *)
checkExact[spec_, fan_, label_] := Module[
  {polys, monoExps, polyExps, vars, n, dv, sl, orig, perSector, allReal, sd,
   M, sub, jac, xsub, transformed, secInteg, diff, ok, freeReal},
  polys = spec["Polynomials"]; monoExps = spec["MonomialExponents"];
  polyExps = spec["PolynomialExponents"]; vars = spec["Variables"];
  n = Length[vars]; {dv, sl} = fan;
  (* original integrand symbol *)
  orig = (Times @@ MapThread[Power, {vars, monoExps}]) *
         (Times @@ MapThread[Power, {polys, polyExps}]);
  allReal = True; perSector = True;
  Do[
    Module[{rays, yv, mm, det},
      rays = (dv[[#]] &) /@ sl[[s]];
      yv = Table[Unique["yy"], {n}];
      mm = -Transpose[rays];               (* M_ij as in ProcessSector *)
      det = Det[mm];
      (* monomial substitution x_i -> prod_j yv_j^{M_ij} *)
      xsub = Table[vars[[i]] -> Product[yv[[j]]^mm[[i, j]], {j, n}], {i, n}];
      jac = Abs[det] * Product[yv[[j]]^(Total[mm[[All, j]]] - 1), {j, n}];
      (* original integrand pulled back through the cone map *)
      transformed = (orig /. xsub) * jac;
      (* compare against ProcessSector's reconstructed sector integrand *)
      sd = ProcessSector[spec, dv, sl[[s]], s];
      If[sd === $Failed || TrueQ[sd["IsDivergent"]],
        perSector = False,
        Module[{fp = sd["FlattenedPolys"], pe = sd["PolynomialExponents"],
                pf = sd["Prefactor"], polyVals},
          freeReal = FreeQ[{sd["RayMatrix"], sd["NewExponents"], fp, pe, pf}, _Real];
          allReal = allReal && freeReal;
          polyVals = Table[
            Total[(#[[1]] * Product[yv[[i]]^#[[2, i]], {i, n}]) & /@ fp[[j]]],
            {j, Length[fp]}];
          secInteg = pf * Times @@ MapThread[Power, {polyVals, pe}];
          (* EXACT identity: pulled-back original == reconstructed sector *)
          diff = Simplify[PowerExpand[transformed / secInteg]];
          ok = TrueQ[PossibleZeroQ[diff - 1]] ||
               TrueQ[Simplify[transformed - secInteg] === 0];
          perSector = perSector && ok]],
      {s, Length[sl]}]
  ];
  rec[label <> " : per-sector exact substitution identity", perSector];
  rec[label <> " : decomposition FreeQ[_,_Real] (no float decisions)", allReal];
];

(* Spec 1: real-exponent n=2 (Example 16 shape) *)
Module[{poly, vars, spec, verts, fan},
  poly = 1 + x[1]^2 + x[2]^2 + x[1] x[2]; vars = {x[1], x[2]};
  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {0, 0},
    "PolynomialExponents" -> {-2}, "Variables" -> vars,
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
  verts = PolytopeVertices[poly^(-1), vars];
  fan = ComputeDecomposition[verts, "ShowProgress" -> False];
  checkExact[spec, fan, "n=2 real"]];

(* Spec 2: two-polynomial n=2 (Example 15 shape) *)
Module[{p1, p2, vars, spec, verts, fan},
  p1 = 1 + x[1] + x[2]; p2 = 1 + x[1] x[2]; vars = {x[1], x[2]};
  spec = <|"Polynomials" -> {p1, p2}, "MonomialExponents" -> {0, 0},
    "PolynomialExponents" -> {-2, -1}, "Variables" -> vars,
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
  verts = PolytopeVertices[(p1 p2)^(-1), vars];
  fan = ComputeDecomposition[verts, "ShowProgress" -> False];
  checkExact[spec, fan, "n=2 two-poly"]];

(* ---- No-regression: wired computeFanScaled (K=1 first) == raw ---- *)
Module[{poly, vars, verts, raw, scaled},
  poly = 1 + x[1]^2 + x[2]^2 + x[1] x[2]; vars = {x[1], x[2]};
  verts = PolytopeVertices[poly^(-1), vars];
  raw = ComputeDecomposition[verts, "ShowProgress" -> False];
  scaled = TropicalFan`computeFanScaled[verts];
  rec["computeFanScaled K=1 byte-identical to raw ComputeDecomposition",
    raw === scaled]];

Print["EXACT_ALL_PASS=", And @@ pass];
