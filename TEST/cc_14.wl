(* ============================================================
   TEST/cc_14.wl  —  v3 Cross-Check #14
   plan.md §8.2 row #14:
       G0 = B^0 − Σ coeff·I   (IBP ↔ subtraction pole identity)

   The two pole-extraction routes share one analytic identity.
   For every divergent sector:
       subtraction pole coefficient  = G0(ck=1) / ck
       IBP pole coefficient          = (1/ck) * [B_boundary − Σ_t coeff_t · I_t]
   These must agree to relative tolerance < 1%.

   Tier 1: WL + g++ only (no CUBA).

   Integral:  ∫ x1^{2ε−1} (1+x1+x2)^{−2} dx1 dx2
   Exact pole coefficient: 1/2 (from Γ(2ε)Γ(0)/Γ(2) → ck=2, G0=1).

   Mechanical port from OLD_CODE RunTest14[] (lines 2206–2368) extended
   to exercise IBP via IBPReduceSector and verify the B^0−Σ identity
   against the subtraction-route G0.
   ============================================================ *)

(* --- load v3 packages using absolute paths ------------------- *)
$v3Root = "/Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3";
Get[FileNameJoin[{$v3Root, "tropical_fan.wl"}]];
Get[FileNameJoin[{$v3Root, "tropical_eval.wl"}]];

(* --- integrand spec ------------------------------------------ *)
eps  = Symbol["cc14eps"];
poly = 1 + x[1] + x[2];
vars = {x[1], x[2]};

spec = <|
  "Polynomials"         -> {poly},
  "MonomialExponents"   -> {2 eps - 1, 0},
  "PolynomialExponents" -> {-2},
  "Variables"           -> vars,
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> eps
|>;

(* --- fan decomposition --------------------------------------- *)
verts       = PolytopeVertices[poly^(-2), vars];
fanData     = Quiet @ ComputeDecomposition[verts, "ShowProgress" -> False];
dualVertices = fanData[[1]];
simplexList  = fanData[[2]];

allSectorData = Table[
  Quiet @ ProcessSector[spec, dualVertices, simplexList[[i]], i],
  {i, Length[simplexList]}
];

divSectors = Select[allSectorData,
  (AssociationQ[#] && TrueQ[#["IsDivergent"]]) &];

If[Length[divSectors] == 0,
  Print["CC14 FAIL expected=divergent-sectors got=none"];
  Quit[1]
];

(* --- helpers ------------------------------------------------- *)

(* Numerically integrate a flat-poly sector integrand over [0,1]^dim.
   flatPolys : list of monomial-list representations {{coeff, {exp1,...}}, ...}
   aVals     : Listable effective exponents (must be numeric)
   polyExps  : polynomial exponents (must be numeric)
   Returns a numeric value or $Failed. *)
integrateFlatSector[flatPolys_List, aVals_List, polyExps_List, detM_] :=
Module[{dim = Length[aVals], yVars, polyVals, integrand},
  yVars = Table[Unique["iy"], {dim}];
  polyVals = Table[
    Total[Table[
      mono[[1]] * Exp[Total[mono[[2]] * Log /@ yVars]],
      {mono, flatPolys[[j]]}
    ]],
    {j, Length[flatPolys]}
  ];
  integrand = Abs[detM] *
    Exp[Total[(aVals - 1) * Log /@ yVars]] *
    Times @@ MapThread[
      Function[{pv, be}, Exp[be * Log[pv]]],
      {polyVals, polyExps}
    ];
  Quiet @ NIntegrate[
    integrand,
    Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
    MaxRecursion -> 20, PrecisionGoal -> 4,
    Method -> "GlobalAdaptive"
  ]
];

(* Compute G0 (subtraction route) from ProcessDivergentSector output. *)
computeG0[divData_] :=
Module[{g0FlatPolys, g0Pf, g0Dim, B0, g0yVars, g0PolyVals, g0Integrand},
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

  Quiet @ NIntegrate[
    g0Integrand,
    Evaluate[Sequence @@ ({#, 0, 1} & /@ g0yVars)],
    MaxRecursion -> 20, PrecisionGoal -> 4,
    Method -> "GlobalAdaptive"
  ]
];

(* --------------------------------------------------------------------------
   Per-sector identity check.
   The B^0 − Σ coeff·I identity (plan.md §8.2 #14):

     Subtraction: pole = G0 / ck      (G0 is the (n-1)-dim integral at eps=0)
     IBP:         pole = B^0 / ck − Σ_t (termCoeff_t / ck) · I_t

   where:
     B^0    = boundary integral  ∫ f(y1,...,y_{k-1},1,y_{k+1},...) at eps=0
     I_t    = reduced-term integrals at eps=0
     ck     = d/dε a_k |_{ε=0}   (the IBP denominator coefficient)

   The identity asserts   G0 = B^0 − Σ_t termCoeff_t · I_t
   (up to numerical MC noise). PASS: |G0 − IBP_num| / |G0| < 0.01.
   -------------------------------------------------------------------------- *)

allPass   = True;
nChecked  = 0;

Do[
  Module[{divData, ibpData,
          ck, g0Val, poleSubtraction,
          k, a0, B0, clearedPolys, detM, n,
          ndVars, bndPolys, bndA0, bndVal,
          terms, termIntegrals, sumTerms, ibpNumerator, poleIBP,
          relErr, pass},

    (* ---- subtraction route ---------------------------------- *)
    divData = Quiet @ ProcessDivergentSector[sd, spec];
    If[!AssociationQ[divData],
      Print["CC14 FAIL expected=Association got=$Failed sector=",
            sd["ConeIndex"]];
      allPass = False;
      Return[]   (* next Do iteration *)
    ];

    ck    = divData["ck"];
    g0Val = computeG0[divData];
    If[!NumericQ[g0Val],
      Print["CC14 FAIL expected=numeric-G0 got=$Failed sector=",
            sd["ConeIndex"]];
      allPass = False;
      Return[]
    ];
    poleSubtraction = g0Val / ck;
    Print["Sector ", sd["ConeIndex"],
          ": subtraction pole = G0/ck = ", N[poleSubtraction, 6]];

    (* ---- IBP route ------------------------------------------ *)
    ibpData = Quiet @ IBPReduceSector[sd, eps];
    If[ibpData === $Failed,
      Print["CC14 FAIL expected=IBPData got=$Failed sector=",
            sd["ConeIndex"]];
      allPass = False;
      Return[]
    ];

    n          = ibpData["Dimension"];
    clearedPolys = ibpData["ClearedPolys"];
    detM       = ibpData["DetM"];
    terms      = ibpData["Terms"];
    (* divergent variable index *)
    k = ibpData["DivergentVariables"][[1]];

    (* B^0: boundary integral — evaluate clearedPolys at y_k=1, eps=0 *)
    (* construct flat boundary polys: set y_k -> 1 in each monomial   *)
    ndVars = DeleteCases[Range[n], k];
    bndPolys = Table[
      Table[
        Module[{mcoeff = mono[[1]], mexp = mono[[2]]},
          (* set y_k exponent to 0 (y_k=1); keep others *)
          {mcoeff, mexp[[ndVars]]}
        ],
        {mono, clearedPolys[[j]]}
      ],
      {j, Length[clearedPolys]}
    ];
    a0  = sd["NewExponents"] /. eps -> 0;
    bndA0 = a0[[ndVars]];   (* monomial exponents of non-divergent vars *)
    B0  = sd["PolynomialExponents"] /. eps -> 0;

    bndVal = integrateFlatSector[bndPolys, bndA0, B0, detM];
    If[!NumericQ[bndVal],
      Print["CC14 FAIL expected=numeric-B0 got=$Failed sector=",
            sd["ConeIndex"]];
      allPass = False;
      Return[]
    ];
    Print["Sector ", sd["ConeIndex"],
          ": IBP boundary B^0 = ", N[bndVal, 6]];

    (* Σ_t termCoeff_t · I_t: sum over IBP-reduced terms at eps=0 *)
    termIntegrals = Table[
      Module[{term, tc, taVals, tpExps, ta0, tp0, tVal},
        term   = terms[[t]];
        tc     = term["Coefficient"];
        taVals = term["NewExponents"] /. eps -> 0;
        tpExps = term["PolyExponents"] /. eps -> 0;
        tVal   = integrateFlatSector[clearedPolys, taVals, tpExps, detM];
        If[NumericQ[tVal], tc * tVal, $Failed]
      ],
      {t, Length[terms]}
    ];

    If[MemberQ[termIntegrals, $Failed],
      Print["CC14 FAIL expected=numeric-term-integrals got=$Failed sector=",
            sd["ConeIndex"]];
      allPass = False;
      Return[]
    ];

    sumTerms   = Total[termIntegrals];
    ibpNumerator = bndVal - sumTerms;   (* B^0 − Σ coeff·I *)
    poleIBP    = ibpNumerator / ck;
    Print["Sector ", sd["ConeIndex"],
          ": IBP pole = (B^0 − Σ coeff·I)/ck = ", N[poleIBP, 6]];

    (* ---- identity check ------------------------------------- *)
    relErr = If[Abs[N[poleSubtraction]] > 1*^-10,
      Abs[N[poleSubtraction - poleIBP]] / Abs[N[poleSubtraction]],
      Abs[N[poleSubtraction - poleIBP]]
    ];
    Print["Sector ", sd["ConeIndex"],
          ": relErr |sub−ibp|/|sub| = ", N[relErr, 4]];

    pass = NumericQ[relErr] && relErr < 0.01;
    If[pass,
      Print["CC14 PASS sector=", sd["ConeIndex"],
            " expected=<0.01 got=", N[relErr, 4]],
      Print["CC14 FAIL expected=<0.01 got=", N[relErr, 4],
            " sector=", sd["ConeIndex"]];
      allPass = False
    ];

    nChecked++;
  ],
  {sd, divSectors}
];

(* --- final verdict ------------------------------------------- *)
If[nChecked == 0,
  Print["CC14 FAIL expected=at-least-one-sector got=0"];
  Quit[1]
];

If[allPass,
  Print["CC14 PASS all=", nChecked,
        " divergent sectors IBP-pole=subtraction-pole within 1%"],
  Print["CC14 FAIL see per-sector output above"]
];
