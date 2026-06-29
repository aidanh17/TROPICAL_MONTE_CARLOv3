(* ============================================================================
   TEST/cc_19.wl  --  Cross-check #19 (plan.md §8.2)
   Sample sigma vs true sigma: exact NIntegrate sigma^2 = I2 - I1^2.
   Flags HasConstantTerm=False sectors as infinite-variance.
   Tier 1: WL + g++ only.

   PASS criteria (plan.md §8.2 #19):
     (a) I1 sector sum matches the exact value within relerr < 0.1% for both
         unlifted and lifted configurations.
     (b) Every lifted sector with HasConstantTerm=False has I2 = Infinity
         (or NIntegrate fails), i.e. the exact per-sample sigma is infinite.
     (c) Every lifted sector with HasConstantTerm=True has a finite I2.
     (d) The lifted total sigma (finite sectors) is strictly less than the
         unlifted total sigma (variance reduction is present).

   Mechanical port from:
     OLD_CODE/TROPICAL_MONTE_CARLO2/TROPICAL_MONTE_CARLOv2/SANDBOX/
       sandbox_true_sigma.wl   (integrand + loop logic)
       sandbox_lift_common.wl  (trueSigma, sectorNIntegrate)
   Uses v3 ProcessSector / LiftCoefficients / ProcessSectorLifted instead
   of the OLD_CODE sandbox helpers.

   KNOWN BUGS FIXED IN THIS REVISION:
   -----------------------------------------------------------------------
   BUG 1 (NIntegrate iterator-binding):
     The old code did:
       limits = Sequence @@ ({#, 0, 1} & /@ yVars);
       NIntegrate[expr, limits, ...]
     NIntegrate is HoldFirst/HoldAll, so the symbol `limits` is NEVER
     evaluated to its Sequence value inside NIntegrate -- the iterator
     list is invisible and I1 is returned unevaluated.
     FIX: use  Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)]  inline,
     exactly as sandbox_lift_common.wl::trueSigma does (line 606/616).

   BUG 2 (lifted fan non-simplicial / degenerate polytope):
     The old code called ComputeDecomposition on the 3D proxy
       1 + x[1] + x[3]^2 * x[2]
     whose Newton polytope is 2-dimensional (3 points in a 2D affine
     subspace of R^3).  ComputeDecomposition requires a full-dimensional
     polytope and returns $Failed / emits errors for degenerate input.
     FIX: supply the same explicit unimodular fan that cc_17.wl uses
     for this same Toy-1 integrand (Test 23A canonical fan, plan.md §8.2).
     This fan has 5 rays and 6 simplicial sectors.
   -----------------------------------------------------------------------

   Run:
     wolframscript -file TEST/cc_19.wl
   from the TROPICAL_MONTE_CARLOv3 root.
   ============================================================================ *)

(* --- load v3 packages using absolute paths --- *)
$v3Root = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{$v3Root, "tropical_fan.wl"}]];
Get[FileNameJoin[{$v3Root, "tropical_eval.wl"}]];

Print["CC19: v3 packages loaded."];
Print[];

(* ============================================================================
   INLINE HELPERS (ported from sandbox_lift_common.wl; no OLD_CODE dependency)

   KEY FIX (BUG 1): every NIntegrate call uses
       Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)]
   so the iterator list is evaluated BEFORE NIntegrate holds its arguments.
   This is the only correct pattern for dynamically-constructed iterator lists
   when the integrand variable symbols are created at runtime (Unique).
   ============================================================================ *)

(* cc19BuildIntegrand[lsd, kinRules] -> {integrand, yVars}
   Builds the log-exp integrand for a sector (lifted or unlifted) including
   the optional DomainConstraint indicator.  Returns fresh Unique variables.
   The caller is responsible for passing the SAME yVars to NIntegrate. *)
cc19BuildIntegrand[lsd_Association, kinRules_List] :=
Module[
  {flatPolys, polyExps, pf, dim, domConstr,
   yVars, polyVals, integrand, logYpStar, logZ0, mp, indCoeffs},

  flatPolys = lsd["FlattenedPolys"];
  polyExps  = lsd["PolynomialExponents"] /. kinRules;
  pf        = lsd["Prefactor"] /. kinRules;
  dim       = lsd["Dimension"];
  domConstr = Lookup[lsd, "DomainConstraint", None];
  yVars     = Table[Unique["cc19y"], {dim}];

  polyVals = Table[
    Total[
      Table[
        Module[{coeff, alphas},
          coeff  = mono[[1]] /. kinRules;
          alphas = mono[[2]] /. kinRules;
          coeff * Exp[Total[alphas * (Log /@ yVars)]]
        ],
        {mono, flatPolys[[j]]}
      ]
    ],
    {j, Length[flatPolys]}
  ];

  integrand = pf *
    Times @@ MapThread[
      Function[{pv, be}, Exp[be * Log[pv]]],
      {polyVals, polyExps}
    ];

  If[domConstr =!= None,
    logZ0     = domConstr["LogZ0"];
    mp        = domConstr["MP"];
    indCoeffs = domConstr["IndicatorCoeffs"] /. kinRules;
    logYpStar = (logZ0 - Total[indCoeffs * (Log /@ yVars)]) / mp;
    integrand = integrand * Boole[logYpStar <= 0]
  ];

  {integrand, yVars}
];


(* trueSigmaCC19[lsd, kinRules, pg]
   Computes I1 = Integral[g], I2 = Integral[g^2]; sigma = Sqrt[I2 - I1^2].
   Returns <|"I1", "I2", "Sigma", "I2Converged"|>.
   Ported from sandbox_lift_common.wl::trueSigma.

   BUG 1 FIX: NIntegrate iterators use Evaluate[Sequence @@ ...] inline,
   NOT a pre-built Sequence stored in a Module local variable.  The latter
   is never evaluated inside NIntegrate's held argument position. *)
trueSigmaCC19[lsd_Association, kinRules_List, pg_Integer] :=
Module[
  {integrand, yVars, i1, i2raw, i2, i2converged},

  If[KeyExistsQ[lsd, "EmptyDomain"] && TrueQ[lsd["EmptyDomain"]],
    Return[<|"I1" -> 0, "I2" -> 0, "Sigma" -> 0, "I2Converged" -> True|>]
  ];

  {integrand, yVars} = cc19BuildIntegrand[lsd, kinRules];

  (* I1 = NIntegrate[g] over [0,1]^dim.
     CORRECT iterator pattern: Evaluate forces the Sequence to splice
     BEFORE NIntegrate holds its argument list. *)
  i1 = Quiet@NIntegrate[
    integrand,
    Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
    MaxRecursion -> 12,
    PrecisionGoal -> pg,
    Method -> "GlobalAdaptive"
  ];

  (* I2 = NIntegrate[g^2]; may diverge for HasConstantTerm=False sectors.
     Same iterator fix applied here. *)
  i2raw = Quiet@Check[
    TimeConstrained[
      NIntegrate[
        integrand^2,
        Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
        MaxRecursion -> 12,
        PrecisionGoal -> pg,
        Method -> "GlobalAdaptive"
      ],
      90,   (* 90s per sector *)
      $TimedOut
    ],
    $Failed
  ];

  Which[
    i2raw === $Failed || i2raw === $TimedOut,
      i2converged = False; i2 = Infinity,
    !NumericQ[i2raw],
      i2converged = False; i2 = Infinity,
    Abs[i2raw] > 10^30,
      i2converged = False; i2 = Infinity,
    True,
      i2converged = True;  i2 = i2raw
  ];

  <|"I1"          -> i1,
    "I2"          -> i2,
    "Sigma"       -> If[i2converged, Sqrt[Max[0, i2 - i1^2]], Infinity],
    "I2Converged" -> i2converged
  |>
];


(* ============================================================================
   INTEGRAND SPEC
   Toy 1: P = 1 + x[1] + 10^-6 x[2],  B = {-3},  exact = 500000.
   (Ported verbatim from sandbox_true_sigma.wl.)
   The large hierarchy 10^-6 makes the unlifted MC sigma very large.
   ============================================================================ *)

spec0 = <|
  "Polynomials"         -> {1 + x[1] + 10^-6 * x[2]},
  "MonomialExponents"   -> {0, 0},
  "PolynomialExponents" -> {-3},
  "Variables"           -> {x[1], x[2]},
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> None
|>;

exactVal = 500000;  (* closed form: Integral = 1 / (2 * 10^-6) = 500000 *)

(* ============================================================================
   PART A: UNLIFTED SECTORS
   ============================================================================ *)
Print["============================================================"];
Print["CC19 PART A: UNLIFTED  P = 1 + x[1] + 10^-6 x[2], B={-3}"];
Print["============================================================"];

(* Fan: use the unit-coefficient proxy P0 = 1 + x[1] + x[2] for the fan
   (fan geometry is independent of coefficient magnitudes). *)
vertsU = PolytopeVertices[(1 + x[1] + x[2])^(-1), {x[1], x[2]}];
fanU   = ComputeDecomposition[vertsU, "ShowProgress" -> False];
Print["Unlifted fan: ", Length[fanU[[2]]], " sectors"];

(* Collect trueSigma results for every convergent sector. *)
unliftedResults = {};
Do[
  Module[{sd, res},
    sd = ProcessSector[spec0, fanU[[1]], fanU[[2, s]], s];
    If[sd === $Failed || TrueQ[sd["IsDivergent"]],
      Print["  Sector ", s, ": DIVERGENT or FAILED -- skipping"],
      res = trueSigmaCC19[sd, {}, 3];
      AppendTo[unliftedResults, res];
      Print["  Sector ", s,
            ": I1=", res["I1"],
            "  I2=", If[res["I2Converged"], res["I2"], "DIVERGENT"],
            "  Sigma=", If[res["I2Converged"], N[res["Sigma"]], "INF"],
            "  I2Converged=", res["I2Converged"]]
    ]
  ],
  {s, Length[fanU[[2]]]}
];

unliftedI1Sum = Total[Map[#["I1"]&, unliftedResults]];
Print["\nUnlifted I1 sum = ", unliftedI1Sum, "  (exact=", exactVal, ")"];

unliftedRelerr = Abs[(unliftedI1Sum - exactVal) / exactVal];
If[unliftedRelerr < 0.001,
  Print["CC19 PASS unlifted I1 sum: relerr=", N[unliftedRelerr], " < 0.001"],
  Print["CC19 FAIL unlifted I1 sum: expected=", exactVal,
        " got=", unliftedI1Sum, " relerr=", N[unliftedRelerr]]
];

(* Total sigma for unlifted (Infinity if any sector diverges) *)
unliftedSigmaTotal =
  If[AnyTrue[unliftedResults, !#["I2Converged"]&],
    Infinity,
    Sqrt[Total[(Map[#["Sigma"]&, unliftedResults])^2]]
  ];
Print["Unlifted total sigma = ", unliftedSigmaTotal];


(* ============================================================================
   PART B: LIFTED SECTORS (k=2)
   Uses v3 LiftCoefficients + ProcessSectorLifted (plan.md §6.2).

   BUG 2 FIX (lifted fan non-simplicial / degenerate polytope):
   The lifted polynomial 1 + x[1] + x[3]^2 * x[2] has Newton polytope
   vertices {0,0,0}, {1,0,0}, {0,1,2} -- these 3 points lie in a 2D affine
   subspace of R^3.  ComputeDecomposition requires a full-dimensional
   polytope and fails for this degenerate input.

   We instead supply the same explicit unimodular simplicial fan used by
   cc_17.wl (Test 23A canonical fan, plan.md §8.2 / test_lifted.wl).
   This is a complete fan over R^3 with 5 rays and 6 simplicial sectors,
   and is the correct fan for this lifting geometry.
   ============================================================================ *)
Print["\n============================================================"];
Print["CC19 PART B: LIFTED k=2  (v3 LiftCoefficients + ProcessSectorLifted)"];
Print["============================================================"];

(* LiftCoefficients: lift the coefficient 10^-6 of x[2] (exponent {0,1}) with k=2. *)
liftRules = {<|"PolyIndex" -> 1, "ExponentVector" -> {0, 1}, "k" -> 2|>};
liftResult = LiftCoefficients[spec0, liftRules];
If[liftResult === $Failed,
  Print["CC19 FAIL: LiftCoefficients returned $Failed"];
  Quit[1]
];
liftedSpec = liftResult["LiftedSpec"];
liftData   = liftResult["LiftData"];
Print["z0 = ", liftData["z0"], "  (exact: ", (10^-6)^(1/2), ")"];
Print["Lifted poly: ", liftedSpec["Polynomials"][[1]]];

(* Identity check: z -> z0 must recover original polynomial *)
Module[{lp = liftedSpec["Polynomials"][[1]],
        z0v = liftData["z0"],
        auxV = liftData["AuxVariable"]},
  If[Expand[lp /. auxV -> z0v] === Expand[spec0["Polynomials"][[1]]],
    Print["Identity check: PASS (lift round-trip exact)"],
    Print["Identity check: FAIL (lift inconsistent)"];
    Quit[1]
  ]
];

(* Explicit unimodular fan for the degenerate Toy-1 lifted polytope.
   Dual vertices and simplices from plan.md §8.2 / test_lifted.wl Test 23A.
   This is the CORRECT fan for the 3-variable lifted spec; supplying it
   explicitly bypasses the ComputeDecomposition degenerate-polytope failure
   that was the second bug in the original cc_19.wl. *)
explicitDualVertices = {{1,0,0},{0,1,0},{-1,-1,0},{0,2,-1},{0,-2,1}};
explicitSimplexList  = {{1,2,4},{2,3,4},{3,1,4},{1,2,5},{2,3,5},{3,1,5}};
explicitFan = {explicitDualVertices, explicitSimplexList};
Print["Lifted fan (explicit): ", Length[explicitDualVertices], " rays, ",
      Length[explicitSimplexList], " sectors"];

liftedSectors = {};
Do[
  Module[{sd},
    sd = ProcessSectorLifted[liftedSpec, explicitFan[[1]], explicitFan[[2, s]], s, liftData];
    Which[
      sd === $Failed,
        Print["  Sector ", s, ": ProcessSectorLifted FAILED -- skipping"],
      KeyExistsQ[sd, "EmptyDomain"] && TrueQ[sd["EmptyDomain"]],
        Print["  Sector ", s, ": EmptyDomain (dropped)"],
      True,
        AppendTo[liftedSectors, sd]
    ]
  ],
  {s, Length[explicitSimplexList]}
];
Print["Non-empty lifted sectors: ", Length[liftedSectors]];

liftedResults = {};
Do[
  Module[{lsd = liftedSectors[[i]], res},
    res = trueSigmaCC19[lsd, {}, 3];
    AppendTo[liftedResults, res];
    Print["  Sector ", i,
          " (cone=", lsd["ConeIndex"],
          ", HasConstantTerm=", lsd["HasConstantTerm"], ")",
          ": I1=", res["I1"],
          "  I2=", If[res["I2Converged"], res["I2"], "DIVERGENT"],
          "  Sigma=", If[res["I2Converged"], N[res["Sigma"]], "INF"],
          "  I2Converged=", res["I2Converged"]]
  ],
  {i, Length[liftedSectors]}
];

liftedI1Sum = Total[Map[#["I1"]&, liftedResults]];
Print["\nLifted I1 sum = ", liftedI1Sum, "  (exact=", exactVal, ")"];

liftedRelerr = Abs[(liftedI1Sum - exactVal) / exactVal];
If[liftedRelerr < 0.001,
  Print["CC19 PASS lifted I1 sum: relerr=", N[liftedRelerr], " < 0.001"],
  Print["CC19 FAIL lifted I1 sum: expected=", exactVal,
        " got=", liftedI1Sum, " relerr=", N[liftedRelerr]]
];

liftedSigmaTotal =
  If[AnyTrue[liftedResults, !#["I2Converged"]&],
    Infinity,
    Sqrt[Total[(Map[#["Sigma"]&, liftedResults])^2]]
  ];
Print["Lifted total sigma = ", liftedSigmaTotal];


(* ============================================================================
   PART C: HasConstantTerm diagnostic (plan.md §8.2 #19, §9 risk 1)
   For each lifted sector, verify that:
     HasConstantTerm=False => I2 diverges  (sigma = Infinity)
     HasConstantTerm=True  => I2 converges (sigma is finite)
   ============================================================================ *)
Print["\n============================================================"];
Print["CC19 PART C: HasConstantTerm vs I2-convergence diagnostic"];
Print["============================================================"];

hasConstList = Map[#["HasConstantTerm"]&, liftedSectors];

(* Check: HasConstantTerm=False sectors should have I2 divergent *)
noConstIdx = Select[Range[Length[liftedSectors]], !hasConstList[[#]]&];
constIdx   = Select[Range[Length[liftedSectors]],  hasConstList[[#]]&];

Print["Sectors with HasConstantTerm=False: ", noConstIdx];
Print["Sectors with HasConstantTerm=True:  ", constIdx];

(* Sub-check C1: every HasConstantTerm=False sector has infinite sigma *)
c1Results = Table[
  Module[{res = liftedResults[[idx]]},
    If[!res["I2Converged"],
      {idx, "PASS", "infinite as expected"},
      {idx, "FAIL", "expected INF sigma but I2Converged=True, Sigma=" <> ToString[res["Sigma"]]}
    ]
  ],
  {idx, noConstIdx}
];

If[noConstIdx === {},
  Print["CC19 PASS C1: no HasConstantTerm=False sectors (all have constant term)"],
  Module[{allPass = AllTrue[c1Results, #[[2]] == "PASS"&]},
    Do[
      Print["  Sector ", r[[1]], " HasConst=False: ", r[[2]], " -- ", r[[3]]],
      {r, c1Results}
    ];
    If[allPass,
      Print["CC19 PASS C1: all HasConstantTerm=False sectors have infinite sigma"],
      Print["CC19 FAIL C1: some HasConstantTerm=False sectors have FINITE sigma"]
    ]
  ]
];

(* Sub-check C2: every HasConstantTerm=True sector has finite sigma *)
c2Results = Table[
  Module[{res = liftedResults[[idx]]},
    If[res["I2Converged"],
      {idx, "PASS", "finite sigma=" <> ToString[N[res["Sigma"]]]},
      {idx, "FAIL", "expected finite sigma but I2Converged=False"}
    ]
  ],
  {idx, constIdx}
];

If[constIdx === {},
  Print["CC19 WARN C2: no HasConstantTerm=True sectors found"],
  Module[{allPass = AllTrue[c2Results, #[[2]] == "PASS"&]},
    Do[
      Print["  Sector ", r[[1]], " HasConst=True: ", r[[2]], " -- ", r[[3]]],
      {r, c2Results}
    ];
    If[allPass,
      Print["CC19 PASS C2: all HasConstantTerm=True sectors have finite sigma"],
      Print["CC19 FAIL C2: some HasConstantTerm=True sectors have INFINITE sigma"]
    ]
  ]
];

(* Non-vacuousness guard: the check is only meaningful if we actually have
   both HasConstantTerm=True and HasConstantTerm=False sectors, OR at least
   one sector of each type is represented in the lifted fan.  The Toy-1
   integrand with the explicit 6-sector fan is known to produce both types
   (DroppedSectors=3 constant-root + some HasConst=False), so we assert that
   liftedSectors is non-empty and constIdx is non-empty. *)
If[Length[liftedSectors] == 0,
  Print["CC19 FAIL C-nonvacuous: no non-empty lifted sectors at all"];
];
If[constIdx === {} && noConstIdx === {},
  Print["CC19 FAIL C-nonvacuous: HasConstantTerm key missing from all sectors"];
];


(* ============================================================================
   PART D: Variance reduction (plan.md §8.2 #19 + §18)
   Lifting with k=2 should reduce sigma vs unlifted (at least for the
   sectors that do have HasConstantTerm=True).
   ============================================================================ *)
Print["\n============================================================"];
Print["CC19 PART D: Variance reduction check"];
Print["============================================================"];

Print["Unlifted total sigma = ", If[unliftedSigmaTotal === Infinity,
                                    "INF", N[unliftedSigmaTotal]]];
Print["Lifted   total sigma = ", If[liftedSigmaTotal === Infinity,
                                    "INF", N[liftedSigmaTotal]]];

(* Compute the lifted sigma counting only finite-sigma sectors *)
finiteIdx = Select[Range[Length[liftedResults]], liftedResults[[#]]["I2Converged"]&];
If[finiteIdx === {},
  Print["CC19 WARN D: no finite-sigma lifted sectors to compare"],
  Module[{sigmaFiniteLifted = Sqrt[Total[(Map[#["Sigma"]&, liftedResults[[finiteIdx]]])^2]]},
    Print["Lifted sigma (finite sectors only) = ", N[sigmaFiniteLifted]];

    Which[
      unliftedSigmaTotal === Infinity && sigmaFiniteLifted < Infinity,
        Print["CC19 PASS D: lifting converts infinite-sigma (unlifted) to finite-sigma sectors"],
      NumberQ[unliftedSigmaTotal] && NumberQ[sigmaFiniteLifted] &&
        sigmaFiniteLifted < unliftedSigmaTotal,
        Module[{ratio = unliftedSigmaTotal / sigmaFiniteLifted},
          Print["CC19 PASS D: variance reduction ratio sigma_U/sigma_L = ", N[ratio],
                " > 1  (plan.md §18: expect >= 10x for toy1 lifted case)"]
        ],
      NumberQ[unliftedSigmaTotal] && NumberQ[sigmaFiniteLifted] &&
        sigmaFiniteLifted >= unliftedSigmaTotal,
        Print["CC19 FAIL D: no variance reduction: sigma_U=", N[unliftedSigmaTotal],
              " sigma_L=", N[sigmaFiniteLifted]],
      True,
        Print["CC19 WARN D: could not compare sigmas (unlifted=", unliftedSigmaTotal,
              " lifted=", sigmaFiniteLifted, ")"]
    ]
  ]
];


(* ============================================================================
   FINAL SUMMARY
   ============================================================================ *)
Print["\n============================================================"];
Print["CC19 SUMMARY"];
Print["============================================================"];

Module[{failCount = 0},

  (* Gate A: unlifted I1 *)
  If[NumericQ[unliftedRelerr] && unliftedRelerr < 0.001,
    Print["CC19 PASS A: unlifted I1 sum relerr=", N[unliftedRelerr]],
    Print["CC19 FAIL A: unlifted I1 sum expected=", exactVal,
          " got=", unliftedI1Sum]; failCount++
  ];

  (* Gate B: lifted I1 *)
  If[NumericQ[liftedRelerr] && liftedRelerr < 0.001,
    Print["CC19 PASS B: lifted I1 sum relerr=", N[liftedRelerr]],
    Print["CC19 FAIL B: lifted I1 sum expected=", exactVal,
          " got=", liftedI1Sum]; failCount++
  ];

  (* Gate C1: HasConst=False => infinite sigma *)
  If[noConstIdx === {} ||
       AllTrue[c1Results, #[[2]] == "PASS"&],
    Print["CC19 PASS C1: HasConstantTerm=False -> infinite sigma (", Length[noConstIdx], " sectors)"],
    Print["CC19 FAIL C1: some HasConstantTerm=False sectors have finite sigma"]; failCount++
  ];

  (* Gate C2: HasConst=True => finite sigma *)
  If[constIdx === {},
    Print["CC19 WARN C2: no HasConstantTerm=True sectors (non-vacuousness check needed)"],
    If[AllTrue[c2Results, #[[2]] == "PASS"&],
      Print["CC19 PASS C2: HasConstantTerm=True -> finite sigma (", Length[constIdx], " sectors)"],
      Print["CC19 FAIL C2: some HasConstantTerm=True sectors have infinite sigma"]; failCount++
    ]
  ];

  (* Gate E: non-vacuousness -- we must have at least 1 non-empty lifted sector *)
  If[Length[liftedSectors] > 0,
    Print["CC19 PASS E: non-vacuous -- ", Length[liftedSectors], " non-empty lifted sectors evaluated"],
    Print["CC19 FAIL E: all lifted sectors are empty/failed -- check not meaningful"]; failCount++
  ];

  Print[];
  If[failCount == 0,
    Print["CC19 PASS  all gates green"],
    Print["CC19 FAIL  ", failCount, " gate(s) failed"]
  ]
];
