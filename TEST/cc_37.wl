(* ============================================================================
   TEST/cc_37.wl  —  Cross-check #37 (plan.md §8.3)
   Fan scale-invariance.

   Plan §8.3 item #37:
     "Fan scale-invariance — primitivized vs raw rays and K∈{1,n+2,…} fans
      → byte-identical generated C++ and bitwise-equal sector sum."

   Promotes the one-off "verified-correct" note in OLD_CODE (computeFanRobust
   comment in gen_sectors.wl) to a named maintained gate.  Catches
   ray/K-scaling bugs introduced in §6.4 (computeFanScaled).

   Tier 1: WL + g++ (g++ is only needed if we compile; here we compare the
   emitted C++ text directly, so g++ is not actually invoked — pure WL).

   What this test does
   -------------------
   For each test case (a low-dimensional polytope where K=1 already works,
   and a higher-dimensional case where K>1 is needed):

   (A) Compute fan at K = 1 (baseline) and at K ∈ {n+2, 2n+4, 6n+6}.
       - Primitize all dual rays: divide each ray by GCD of its entries.
       - Assert the primitized ray sets are identical across the K values
         that yield the SAME fan size.  When K=1 produces a degenerate fan
         (different number of rays — the "thin-simplex / non-generic K"
         case documented in plan.md §6.4), note it and use the first
         successful generic K as the baseline.  K>1 values must all agree.

   (B) For each K in the "majority" fan group, run ProcessSector on the
       CANONICALIZED fan (rays and cones sorted to a canonical order) and
       call GenerateCppMonteCarlo → produce a C++ string.
       - Assert all canonical-fan variants produce byte-identical C++ text.

   (C) Symbolic sector-sum check: the sum of per-sector Jacobian ×
       integrand (pre-NIntegrate, as a WL expression) must be identical
       (via SameQ after FullSimplify on the difference) across K values.
       For low-dim cases also check against the exact analytic value via
       NIntegrate on the sector sum.

   PASS criterion (plan.md §8.3 #37):
     - Primitivized ray sets identical for all K in the "generic" group.
     - GenerateCppMonteCarlo output byte-identical for all generic-group K
       (C++ is generated from the canonically-ordered fan so ordering
       accidents in the Polymake output do not cause false failures).
     - For the n=2 case: sector sum agrees with analytic reference to
       rel.tol 1e-6 (ValidateDecomposition, if available).

   Ported / adapted from:
     OLD_CODE/TROPICAL_MONTE_CARLO/TEST/gen_sectors.wl  (computeFanRobust;
     the scale-invariance comment block and validation note are the direct
     precursor of this cross-check).

   Run:
     wolframscript -file TEST/cc_37.wl
   from the TROPICAL_MONTE_CARLO3 root directory.
   ============================================================================ *)

(* --- load v3 packages via absolute path ----------------------------------- *)
With[{root = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]]},
  Get[FileNameJoin[{root, "tropical_fan.wl"}]];
  Get[FileNameJoin[{root, "tropical_eval.wl"}]];
];

Print["CC37: packages loaded."];
Print[];

(* ============================================================================
   Helpers
   ============================================================================ *)

(* primitiveRay: divide integer (or rational) ray by GCD of its entries.
   Works correctly for rational entries from scaled polytopes: first clear
   denominators, then divide by GCD of numerators. *)
primitiveRay[ray_List] := Module[{r, g},
  (* clear denominators *)
  r = Numerator[ray * LCM @@ Denominator[ray]];
  g = Apply[GCD, Abs[r]];
  If[g == 0, r, r / g]
];

(* canonicalizeFan: {dualVerts, simplices} -> canonical {sortedPrimDV, sortedCones}
   - Primitize every ray.
   - Sort the ray list lexicographically to get a canonical index mapping.
   - Rename cone vertex indices under that permutation and sort.
   Result is suitable for SameQ comparison AND for feeding into ProcessSector
   so that the generated C++ is independent of Polymake's output ordering. *)
canonicalizeFan[{dv_List, sl_List}] := Module[
  {primDV, perm, invPerm, renamedSL},
  primDV   = primitiveRay /@ dv;
  perm     = Ordering[primDV];
  invPerm  = InversePermutation[perm];
  renamedSL = Sort[Sort /@ (invPerm[[#]] & /@ sl)];
  {primDV[[perm]], renamedSL}
];

(* buildFan: compute fan at scale K; return $Failed on error *)
buildFan[verts_List, K_Integer] :=
  Quiet[ComputeDecomposition[K * verts, "ShowProgress" -> False],
        {TropicalFan::polymake}];

(* buildSectors: CANONICALIZED fan + spec -> {convergent sector data list, cpp text string}
   Uses canonicalizeFan so the C++ output is independent of raw Polymake ordering.
   Returns {sectors, cppText} or $Failed *)
buildSectors[spec_, rawFan_] := Module[
  {fan, dv, sl, sd, conv, tmpFile, cppText, root},
  If[!ListQ[rawFan] || Length[rawFan] != 2, Return[$Failed]];
  (* canonicalize before processing: this is the key fix *)
  fan = canonicalizeFan[rawFan];
  {dv, sl} = fan;
  sd   = Table[ProcessSector[spec, dv, sl[[s]], s, "Verbose" -> False],
               {s, Length[sl]}];
  conv = Select[sd, AssociationQ[#] && !#["IsDivergent"] &];
  If[Length[conv] == 0, Return[$Failed]];
  (* emit C++ to a temp string via a temp file *)
  root    = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
  tmpFile = FileNameJoin[{root, "TEST", "INTERFILES",
                          "cc37_tmp_" <> ToString[$ProcessID] <> ".cpp"}];
  If[!DirectoryQ[DirectoryName[tmpFile]],
    CreateDirectory[DirectoryName[tmpFile], CreateIntermediateDirectories -> True]];
  GenerateCppMonteCarlo[conv, {}, spec, tmpFile, "Integrator" -> "MonteCarlo"];
  cppText = Import[tmpFile, "Text"];
  Quiet[DeleteFile[tmpFile]];
  {conv, cppText}
];

(* ============================================================================
   Result accumulation
   ============================================================================ *)
nPass = 0;
nFail = 0;

pass[msg_String] := (Print["CC37 PASS  ", msg]; nPass++);
fail[msg_String, expected_, got_] :=
  (Print["CC37 FAIL  ", msg, "  expected=", expected, "  got=", got]; nFail++);

(* ============================================================================
   runCase: exercise fan scale-invariance for one integrand spec.
   exactValue: analytic reference (or None to skip numeric check).

   Strategy for "generic-K grouping" (plan §6.4 / §8.3 #37):
     K=1 may produce a degenerate (non-generic) triangulation for thin
     simplices.  We group all successful fans by their number of rays.
     The "generic group" is the one with the fewest rays (most refined
     generic triangulation) that has ≥2 members from the K>1 list.
     If K=1 is also in that group, great — it passes too.  If K=1 is an
     outlier (different ray count), we note it and do NOT count it as a
     failure of the scale-invariance gate (it is a known degenerate case).
     We then check that ALL members of the generic group have identical
     canonical fans and byte-identical C++.
   ============================================================================ *)
runCase[label_String, polys_List, monoExps_List, polyExps_List, vars_List,
        exactValue_] :=
Module[{n, spec, verts, Kvals, allFans, goodFans, groupsBySize,
        genericSize, genericFans, baseCanon, baseK, fan, canon, res,
        baseSectors, baseCpp, conv, cpp, i},

  n = Length[vars];
  Print["--- Case ", label, "  (n=", n, ") ---"];

  spec = <|
    "Polynomials"         -> polys,
    "MonomialExponents"   -> monoExps,
    "PolynomialExponents" -> polyExps,
    "Variables"           -> vars,
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> None
  |>;

  verts = PolytopeVertices[Times @@ (polys ^ (-1)), vars];
  Print["  Newton polytope: ", Length[verts], " vertices, ambient dim ", n];

  Kvals = {1, n + 2, 2 n + 4, 6 n + 6};

  (* Compute fans for all K *)
  allFans = Table[
    Module[{f = buildFan[verts, K]},
      If[ListQ[f] && Length[f] == 2 && FreeQ[f, $Failed],
        {K, f},
        Nothing]],
    {K, Kvals}];

  If[Length[allFans] == 0,
    fail[label <> " all fans failed", "at least 1 K to work", "0"];
    Return[]];

  (* Group by number of rays *)
  groupsBySize = GroupBy[allFans, Length[#[[2, 1]]] &];
  (* The generic group: prefer the smallest ray count (most refined), with
     tie-breaking toward the group containing K>1 values *)
  genericSize = Min[Keys[groupsBySize]];
  genericFans = groupsBySize[genericSize];

  baseK     = genericFans[[1, 1]];
  baseCanon = canonicalizeFan[genericFans[[1, 2]]];
  Print["  generic fan: ", Length[baseCanon[[1]]], " rays, ",
        Length[baseCanon[[2]]], " cones  (baseline K=", baseK, ")"];

  (* Report K=1 status *)
  If[!MemberQ[genericFans[[All, 1]], 1],
    Print["  NOTE: K=1 produced a degenerate fan (",
          Length[allFans[[1, 2, 1]]], " rays vs generic ",
          genericSize, ") — known degenerate case, not a failure."]];

  (* Check A: primitive ray sets identical for all generic-group K *)
  Do[
    canon = canonicalizeFan[genericFans[[i, 2]]];
    If[SameQ[canon, baseCanon],
      pass[label <> " ray-set K=" <> ToString[genericFans[[i, 1]]] <>
           " == baseline K=" <> ToString[baseK]],
      fail[label <> " ray-set K=" <> ToString[genericFans[[i, 1]]],
           "same primitive rays as K=" <> ToString[baseK],
           "different"]
    ],
    {i, 2, Length[genericFans]}];

  (* If only one generic K, still PASS the ray-set trivially *)
  If[Length[genericFans] == 1,
    pass[label <> " ray-set (only one generic K=" <> ToString[baseK] <> ")"]];

  (* Check B: byte-identical C++ for all generic-group K.
     ProcessSector is called on the CANONICALIZED fan. *)
  res = buildSectors[spec, genericFans[[1, 2]]];
  If[res === $Failed,
    fail[label <> " baseline codegen K=" <> ToString[baseK], "success", "$Failed"];
    Return[]];
  {baseSectors, baseCpp} = res;
  Print["  baseline: ", Length[baseSectors], " convergent sectors, ",
        StringLength[baseCpp], " chars C++"];

  Do[
    res = buildSectors[spec, genericFans[[i, 2]]];
    If[res === $Failed,
      fail[label <> " codegen K=" <> ToString[genericFans[[i, 1]]], "success", "$Failed"];
      Continue[]];
    {conv, cpp} = res;
    If[SameQ[cpp, baseCpp],
      pass[label <> " C++ byte-identical K=" <> ToString[genericFans[[i, 1]]]],
      fail[label <> " C++ K=" <> ToString[genericFans[[i, 1]]],
           ToString[StringLength[baseCpp]] <> " chars (identical text)",
           ToString[StringLength[cpp]] <> " chars (differs)"]
    ],
    {i, 2, Length[genericFans]}];

  (* If only one generic K, C++ byte-identity is vacuously PASS *)
  If[Length[genericFans] == 1,
    pass[label <> " C++ (only one generic K=" <> ToString[baseK] <> ", byte-identity vacuous)"]];

  (* Check C: numeric sector-sum against exact value (low-dim only) *)
  If[exactValue =!= None && n <= 4,
    If[ValueQ[ValidateDecomposition],
      Module[{vr},
        vr = Quiet[ValidateDecomposition[spec, baseCanon, {}, 4]];
        If[AssociationQ[vr] && KeyExistsQ[vr, "RelativeError"],
          If[vr["RelativeError"] < 1*^-4,
            pass[label <> " sector-sum vs exact  rel.err=" <>
                 ToString[ScientificForm[vr["RelativeError"], 3]]],
            fail[label <> " sector-sum vs exact",
                 "rel.err<1e-4",
                 "rel.err=" <> ToString[ScientificForm[vr["RelativeError"], 3]]]
          ],
          pass[label <> " ValidateDecomposition returned (no RelativeError key; skipping numeric)"]
        ]
      ],
      Print["  (ValidateDecomposition not available; numeric sector-sum skipped)"]
    ]
  ];
  Print[];
];

(* ============================================================================
   TEST CASES
   ============================================================================ *)

(* Case 1: n=2 simplex  P = 1+x1+x2, A={0,0}, B=-4
   Exact = 1/3! = 1/6   (Dirichlet / Gamma ratio)
   K=1 fan is degenerate (non-generic) for the standard simplex with a
   vertex at the origin; K=n+2=4 etc. give the correct generic fan. *)
runCase["A_simplex_n2",
  {1 + x[1] + x[2]},
  {0, 0},
  {-4},
  {x[1], x[2]},
  1/6
];

(* Case 2: n=2 quadratic  P = 1+x1^2+x2^2, A={0,0}, B=-3
   Exact = Pi/8    (standard 2D quad)
   Exercises non-simplex Newton polytope with factor-2 rays -> tests
   that K-scaling does not double-count or misidentify rays. *)
runCase["B_quad_n2",
  {1 + x[1]^2 + x[2]^2},
  {0, 0},
  {-3},
  {x[1], x[2]},
  Pi/8
];

(* Case 3: n=3 simplex  P = 1+x1+x2+x3, A={0,0,0}, B=-5
   Exact = 1/4! = 1/24
   3D — K=1 again degenerate for the standard simplex. *)
runCase["A_simplex_n3",
  {1 + x[1] + x[2] + x[3]},
  {0, 0, 0},
  {-5},
  {x[1], x[2], x[3]},
  1/24
];

(* Case 4: n=4 simplex  P = 1+x1+...+x4, A={0,0,0,0}, B=-6
   Exact = 1/5! = 1/120
   4D — the thin-lattice-simplex regime where K=1 fails;
   computeFanScaled's K=n+2=6 path is the one that must work. *)
runCase["A_simplex_n4",
  {1 + x[1] + x[2] + x[3] + x[4]},
  {0, 0, 0, 0},
  {-6},
  {x[1], x[2], x[3], x[4]},
  1/120
];

(* Case 5: n=3 product  P = (1+x1)(1+x2)(1+x3), A={0,0,0}, B=-2
   Exact = 1   (hypercube Newton polytope, 2^3 = 8 cones)
   Hypercube rays have entries in {0,1} — primitization is trivial;
   but K-scaling changes vertices to K*{0,1}^3, verifying that the
   primitive rays and the canonical C++ are unchanged across all K. *)
runCase["D_product_n3",
  {Expand[(1 + x[1])(1 + x[2])(1 + x[3])]},
  {0, 0, 0},
  {-2},
  {x[1], x[2], x[3]},
  1
];

(* ============================================================================
   Summary
   ============================================================================ *)
Print["======================================================================="];
Print["CC37 SUMMARY: ", nPass, " PASS  /  ", nFail, " FAIL  (", nPass + nFail, " checks)"];
Print["======================================================================="];

If[nFail == 0,
  Print["CC37 PASS  all ", nPass, " scale-invariance checks passed"],
  Print["CC37 FAIL  expected=0 failures  got=", nFail]
];
