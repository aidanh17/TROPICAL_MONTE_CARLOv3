(* ============================================================
   CC13 — Pole coefficient vs numerical eps-slope estimate
   plan.md §8 cross-check #13 (§8.2): F2a (IBP default / subtraction)
   Tier 1: WL + g++ only (no CUBA required; uses NIntegrate).

   Mechanical port of OLD_CODE TROPICAL_MONTE_CARLO/EXAMPLES/
   tropical_eval_examples.wl RunTest14[] (lines 2206-2368),
   adapted to v3 paths and naming conventions.

   Integral: ∫₀^∞ x1^{2ε-1} (1+x1+x2)^{-2} dx1 dx2
   → has a 1/ε pole from the x1→0 sector.

   PASS criterion (plan.md §8.2 #13): relerr < 5% for every
   divergent sector.
   ============================================================ *)

(* --- load v3 packages using absolute paths ------------------- *)
$v3Root = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{$v3Root, "tropical_fan.wl"}]];
Get[FileNameJoin[{$v3Root, "tropical_eval.wl"}]];

(* --- integrand spec ------------------------------------------ *)
eps   = Symbol["cc13eps"];
poly  = 1 + x[1] + x[2];
vars  = {x[1], x[2]};

spec = <|
  "Polynomials"         -> {poly},
  "MonomialExponents"   -> {2 eps - 1, 0},
  "PolynomialExponents" -> {-2},
  "Variables"           -> vars,
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> eps
|>;

(* --- tropical fan -------------------------------------------- *)
verts       = PolytopeVertices[poly^(-2), vars];
fanData     = ComputeDecomposition[verts, "ShowProgress" -> False];
dualVertices = fanData[[1]];
simplexList  = fanData[[2]];

(* --- sector processing --------------------------------------- *)
allSectorData = Table[
  Quiet @ ProcessSector[spec, dualVertices, simplexList[[i]], i],
  {i, Length[simplexList]}
];

divSectors = Select[allSectorData,
  (AssociationQ[#] && TrueQ[#["IsDivergent"]]) &];

If[Length[divSectors] == 0,
  Print["CC13 FAIL expected=divergent-sectors got=none"];
  Quit[1]
];

(* --- per-sector check ---------------------------------------- *)
allPass   = True;
sectorLog = {};

Do[
  Module[{divData, ck, g0Val,
          g0FlatPolys, g0Pf, g0Dim, B0, g0yVars,
          g0PolyVals, g0Integrand,
          poleCoeff,
          eps1, eps2, Ieps1, Ieps2,
          yVars, clearedPolys, aNum1, aNum2,
          polyVals1, polyVals2, integ1, integ2, polyExps,
          slopeCoeff, relErr, sectorPass},

    (* ProcessDivergentSector: tropical subtraction scheme *)
    divData = Quiet @ ProcessDivergentSector[sd, spec];

    If[!AssociationQ[divData],
      Print["CC13 FAIL expected=Association got=$Failed for sector ",
            sd["ConeIndex"]];
      allPass = False;
      Return[]   (* next Do iteration *)
    ];

    ck = divData["ck"];

    (* --- G0 via NIntegrate ---------------------------------- *)
    g0FlatPolys = divData["G0FlatPolys"];
    g0Pf        = divData["G0Prefactor"];
    g0Dim       = divData["G0Dimension"];
    B0          = divData["G0PolyExponents"];
    g0yVars     = Table[Unique["g0y"], {g0Dim}];

    g0PolyVals = Table[
      Total[Table[
        mono[[1]] * Exp[Total[mono[[2]] * Log /@ g0yVars]],
        {mono, g0FlatPolys[[j]]}
      ]],
      {j, Length[g0FlatPolys]}
    ];

    g0Integrand = g0Pf *
      Times @@ MapThread[
        Function[{pv, be}, Exp[be * Log[pv]]],
        {g0PolyVals, B0}
      ];

    g0Val = Quiet @ NIntegrate[
      g0Integrand,
      Evaluate[Sequence @@ ({#, 0, 1} & /@ g0yVars)],
      MaxRecursion -> 20, PrecisionGoal -> 4,
      Method -> "GlobalAdaptive"
    ];

    poleCoeff = g0Val / ck;
    Print["  Sector ", sd["ConeIndex"], ": G0/ck = ", poleCoeff];

    (* --- Sector integrand at two small eps values ----------- *)
    eps1 = 0.01;  eps2 = 0.02;
    yVars       = Table[Unique["py"], {sd["Dimension"]}];
    clearedPolys = divData["ClearedPolys"];
    polyExps     = sd["PolynomialExponents"];

    (* helper: sector integrand at a numeric eps value *)
    sectorIntegrand[epsVal_] := Module[{aNum, pvs},
      aNum = sd["NewExponents"] /. eps -> epsVal;
      pvs  = Table[
        Total[Table[
          mono[[1]] * Exp[Total[mono[[2]] * Log /@ yVars]],
          {mono, clearedPolys[[j]]}
        ]],
        {j, Length[clearedPolys]}
      ];
      Abs[sd["DetM"]] *
        Exp[Total[(aNum - 1) * Log /@ yVars]] *
        Times @@ MapThread[
          Function[{pv, be}, Exp[be * Log[pv]]],
          {pvs, polyExps /. eps -> epsVal}
        ]
    ];

    Ieps1 = Quiet @ NIntegrate[
      sectorIntegrand[eps1],
      Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
      MaxRecursion -> 20, PrecisionGoal -> 3,
      Method -> "GlobalAdaptive"
    ];

    Ieps2 = Quiet @ NIntegrate[
      sectorIntegrand[eps2],
      Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
      MaxRecursion -> 20, PrecisionGoal -> 3,
      Method -> "GlobalAdaptive"
    ];

    (* Numerical slope: I(eps) ~ C/eps + D + ...
       d I / d(1/eps) = C
       slopeCoeff = (I(eps1) - I(eps2)) / (1/eps1 - 1/eps2)  *)
    slopeCoeff = (Ieps1 - Ieps2) / (1/eps1 - 1/eps2);
    Print["  Sector ", sd["ConeIndex"], ": numerical slope = ", slopeCoeff];

    relErr = If[NumericQ[poleCoeff] && Abs[poleCoeff] > 0,
      Abs[(poleCoeff - slopeCoeff) / poleCoeff],
      Infinity
    ];
    Print["  relErr = ", relErr];

    sectorPass = NumericQ[relErr] && relErr < 0.05;
    AppendTo[sectorLog, <|"Sector" -> sd["ConeIndex"],
                          "PoleCoeff" -> poleCoeff,
                          "Slope"     -> slopeCoeff,
                          "RelErr"    -> relErr,
                          "Pass"      -> sectorPass|>];

    If[sectorPass,
      Print["CC13 PASS sector=", sd["ConeIndex"],
            " expected=<5% got=", relErr],
      Print["CC13 FAIL expected=<0.05 got=", relErr,
            " sector=", sd["ConeIndex"]];
      allPass = False
    ]
  ],
  {sd, divSectors}
];

(* --- final verdict ------------------------------------------- *)
If[allPass,
  Print["CC13 PASS all=", Length[divSectors],
        " divergent sectors pole-vs-slope relerr<5%"],
  Print["CC13 FAIL see per-sector output above"]
];
