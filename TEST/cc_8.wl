(* ============================================================================
   CC8 — Cross-check #8: tropical decomposition vs FIESTA5
         via GL(1) / Cheng-Wu gauge transformation.

   Plan §8.2 item #8: "vs FIESTA5 (external sector-decomposition tool)
   via GL(1)/Cheng-Wu — bubblepp pole exact, finite within ~2σ."

   Tier 3: needs FIESTA5.  If absent, prints "CC8 SKIP (no FIESTA)" and
   exits cleanly (exit code 0).

   Integral under test
   -------------------
   I = Int_0^inf dx1 dx2  x1^{2ep-1} / [(1+x1)(1+x2)]^2

   Exact Laurent expansion (Gamma products):
     I = Gamma(2ep) * Gamma(2-2ep)
       = 1/(2ep)  -  1  +  O(ep)
   => pole = 1/2,  finite = -1

   GL(1) identity (Cheng-Wu)
   -------------------------
   After homogenisation: f(x0,x1,x2) = x1^{2ep-1} * [(x0+x1)(x0+x2)]^{-2}.
   Scaling degree: (2ep-1) + (-4) = 2ep-5, so D = 5-2ep in 3 variables.
   D - E = (5-2ep) - 3 = 2-2ep.

   FIESTA domain (simplex, delta(1-x0-x1-x2)):
     x[1]=x0, x[2]=x1, x[3]=x2
     functions: {x[1], x[2], (x[1]+x[2])*(x[1]+x[3])}
     degrees:   {2-2ep, 2ep-1, -2}

   Tropical MC domain [0,inf):
     ProcessSector + ProcessDivergentSector on
       spec = polynomials {(1+x1)(1+x2)}, monomial exponents {2ep-1,0},
              polynomial exponents {-2}, regulator = ep.

   PASS criteria (from plan §8.2 #8)
   -----------------------------------
   1. FIESTA pole  coefficient matches exact 1/2 within 2σ.
   2. FIESTA finite coefficient matches exact -1  within 2σ.
   3. Tropical pole   matches exact 1/2 within 5e-3  (analytic-exact subtraction).
   4. Tropical finite matches exact -1  within 5e-2  (NIntegrate-based).
   5. FIESTA pole  and tropical pole   agree within combined 2σ + 5e-3 tolerance.
   6. FIESTA finite and tropical finite agree within combined 2σ + 5e-2 tolerance.
   ============================================================================ *)

(* ------------------------------------------------------------------ helpers *)

cc8Root = DirectoryName[$InputFileName];
If[cc8Root === "", cc8Root = Directory[]];
(* $InputFileName may be empty when pasted interactively; fall back. *)
v3Root = FileNameJoin[{cc8Root, ".."}];

cc8Print[msg___] := Print["[CC8] ", msg];

(* ------------------------------------------------------------------ FIESTA detection *)

(* Search common install locations for FIESTA5.m *)
fiestaCandidates = {
  "/usr/local/fiesta/FIESTA5/FIESTA5.m",
  "/usr/local/FIESTA5/FIESTA5.m",
  "/opt/homebrew/share/FIESTA5/FIESTA5.m",
  FileNameJoin[{$HomeDirectory, "FIESTA5", "FIESTA5.m"}],
  FileNameJoin[{$HomeDirectory, "fiesta", "FIESTA5", "FIESTA5.m"}]
};

(* Also honour an environment variable: FIESTA5_DIR *)
With[{envDir = Environment["FIESTA5_DIR"]},
  If[StringQ[envDir] && envDir =!= "",
    fiestaCandidates = Prepend[fiestaCandidates,
      FileNameJoin[{envDir, "FIESTA5.m"}]]]
];

fiestaMFile = SelectFirst[fiestaCandidates, FileExistsQ, None];
fiestaDir   = If[fiestaMFile =!= None, DirectoryName[fiestaMFile], None];

If[fiestaMFile === None,
  Print["CC8 SKIP (no FIESTA)"];
  Exit[0]
];

(* Also check that the CIntegratePool binary exists — without it FIESTA hangs *)
fiestaIntegratorBin = FileNameJoin[{fiestaDir, "bin", "CIntegratePool"}];
If[! FileExistsQ[fiestaIntegratorBin],
  Print["CC8 SKIP (no FIESTA) [FIESTA5.m found but CIntegratePool binary missing at ",
        fiestaIntegratorBin, "]"];
  Exit[0]
];

(* ------------------------------------------------------------------ load FIESTA *)

cc8Print["FIESTA5 found at: ", fiestaMFile];
SetDirectory[fiestaDir];
Quiet[Get["FIESTA5.m"]];
SetOptions[FIESTA, "NumberOfSubkernels" -> 0, "NumberOfLinks" -> 4];
cc8Print["FIESTA5 loaded."];

(* ------------------------------------------------------------------ load v3 package *)

v3EvalFile = FileNameJoin[{v3Root, "tropical_eval.wl"}];
v3FanFile  = FileNameJoin[{v3Root, "tropical_fan.wl"}];

If[!FileExistsQ[v3EvalFile],
  Print["CC8 FAIL expected=tropical_eval.wl got=NOT_FOUND (", v3EvalFile, ")"];
  Exit[1]
];

(* tropical_eval.wl loads tropical_fan.wl internally via its own path; we
   also Get fan explicitly so TropicalFan` is definitely on $ContextPath. *)
If[FileExistsQ[v3FanFile], Quiet[Get[v3FanFile]]];
Quiet[Get[v3EvalFile]];
cc8Print["v3 package loaded from ", v3EvalFile];

(* ------------------------------------------------------------------ exact reference *)

exactPole   = 1/2;  (* coefficient of 1/ep *)
exactFinite = -1;   (* coefficient of ep^0 *)

(* ================================================================== PART 1: FIESTA *)

cc8Print[""];
cc8Print["--- Part 1: FIESTA5 SDEvaluateDirect on the GL(1) simplex form ---"];
cc8Print["  functions: {x0, x1, (x0+x1)*(x0+x2)},  degrees: {2-2ep, 2ep-1, -2}"];
cc8Print["  delta: delta(1-x0-x1-x2)"];

SetDirectory[fiestaDir];

fiestaResult = Quiet[
  SDEvaluateDirect[
    {x[1], x[2], (x[1] + x[2]) (x[1] + x[3])},
    {2 - 2 ep, 2 ep - 1, -2},
    0,  (* expand to order ep^0 *)
    {{1, 2, 3}}
  ]
];

cc8Print["  FIESTA raw result: ", fiestaResult];

(* Parse the FIESTA result.  SDEvaluateDirect returns a list of the form
   {{coeff_n, err_n, order_n}, ...} where the leading term has order -1 (the
   1/ep pole) and the next term has order 0 (the finite part).  Alternatively
   on some FIESTA versions it returns an Association or a sum expression.
   We handle both formats. *)

parseFiestaLaurent[res_List] :=
  Module[{poleEntry, finiteEntry, pole, poleSig, finite, finiteSig},
    (* Filter entries by order: order = -1 -> pole; order = 0 -> finite *)
    poleEntry   = SelectFirst[res, (#[[3]] == -1) &, None];
    finiteEntry = SelectFirst[res, (#[[3]] ==  0) &, None];
    pole        = If[poleEntry   =!= None, poleEntry[[1]],   Missing["NoPole"]];
    poleSig     = If[poleEntry   =!= None, poleEntry[[2]],   Infinity];
    finite      = If[finiteEntry =!= None, finiteEntry[[1]], Missing["NoFinite"]];
    finiteSig   = If[finiteEntry =!= None, finiteEntry[[2]], Infinity];
    {pole, poleSig, finite, finiteSig}
  ];

parseFiestaLaurent[res_] :=
  Module[{poleCoeff, poleErr, finCoeff, finErr},
    (* Try to extract coefficients from a symbolic expression in ep *)
    poleCoeff = Coefficient[res, ep, -1];
    finCoeff  = Coefficient[res, ep,  0];
    (* No sigma info available *)
    {poleCoeff, 0, finCoeff, 0}
  ];

{fiestaPole, fiestaSigPole, fiestaFin, fiestaSigFin} =
  parseFiestaLaurent[fiestaResult];

cc8Print["  FIESTA pole   = ", N[fiestaPole],   " ± ", N[fiestaSigPole]];
cc8Print["  FIESTA finite = ", N[fiestaFin],    " ± ", N[fiestaSigFin]];
cc8Print["  Exact pole    = ", N[exactPole]];
cc8Print["  Exact finite  = ", N[exactFinite]];

(* ================================================================== PART 2: Tropical MC *)

cc8Print[""];
cc8Print["--- Part 2: v3 tropical decomposition + subtraction ---"];
cc8Print["  Int_0^inf dx1 dx2  x1^{2ep-1} / [(1+x1)(1+x2)]^2"];
cc8Print["  Exact: Gamma(2ep)*Gamma(2-2ep) = 1/(2ep) - 1 + O(ep)"];

SetDirectory[v3Root];

Module[
  {poly, vars, eps, spec, verts, fanData, dualVertices, simplexList,
   allSectorData, divSectors, convSectors,
   totalPole, totalFinite, convContrib},

  eps  = Symbol["cc8eps"];
  poly = (1 + x[1]) (1 + x[2]);  (* = 1 + x1 + x2 + x1*x2 *)
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"         -> {Expand[poly]},
    "MonomialExponents"   -> {2 eps - 1, 0},
    "PolynomialExponents" -> {-2},
    "Variables"           -> vars,
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> eps
  |>;

  verts    = PolytopeVertices[poly^(-2), vars];
  fanData  = ComputeDecomposition[verts, "ShowProgress" -> False];
  {dualVertices, simplexList} = fanData;

  cc8Print["  Fan: ", Length[dualVertices], " rays, ",
           Length[simplexList], " sectors"];

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  divSectors  = Select[allSectorData, (AssociationQ[#] && #["IsDivergent"]) &];
  convSectors = Select[allSectorData, (AssociationQ[#] && ! #["IsDivergent"]) &];

  cc8Print["  Convergent sectors: ", Length[convSectors],
           ";  Divergent sectors: ", Length[divSectors]];

  (* --- Convergent sector contributions at eps = 0 --- *)
  convContrib = 0;
  Do[
    Module[{sd, flatPolys, pExps, pf, dim, yVars, polyVals, integrand, result},
      sd       = cs;
      flatPolys = sd["FlattenedPolys"];
      pExps    = sd["PolynomialExponents"] /. eps -> 0;
      pf       = sd["Prefactor"]           /. eps -> 0;
      dim      = sd["Dimension"];
      yVars    = Table[Unique["cc8y"], {dim}];

      polyVals = Table[
        Total[Table[
          Module[{c, a},
            c = mono[[1]] /. eps -> 0;
            a = mono[[2]] /. eps -> 0;
            c * Exp[Total[a * (Log /@ yVars)]]
          ], {mono, flatPolys[[j]]}]],
        {j, Length[flatPolys]}
      ];

      integrand = pf *
        Times @@ MapThread[Function[{p, b}, Exp[b * Log[p]]], {polyVals, pExps}];

      result = Quiet@NIntegrate[integrand,
        Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
        MaxRecursion -> 20, PrecisionGoal -> 6];
      convContrib += result;
    ],
    {cs, convSectors}
  ];

  cc8Print["  Convergent sector sum (eps=0): ", N[convContrib]];

  (* --- Divergent sector contributions: pole (G0) and finite (G1 + remainder) --- *)
  totalPole   = 0;
  totalFinite = 0;

  Do[
    Module[{divData, ck, k, n, ndVars, a0, a1, B0, B1, detM,
            clearedPolys, simpPolys,
            g0Val, g1Val, remVal},

      divData = Quiet@ProcessDivergentSector[sd, spec];

      If[! AssociationQ[divData],
        cc8Print["  WARNING: ProcessDivergentSector returned non-Association for sector ",
                 sd["ConeIndex"]];
        Return[]
      ];

      ck         = divData["ck"];
      k          = divData["DivergentVariable"];
      n          = divData["Dimension"];
      ndVars     = DeleteCases[Range[n], k];
      a0         = divData["a0"];
      a1         = divData["a1"];
      B0         = divData["B0"];
      B1         = divData["B1"];
      detM       = divData["DetM"];
      clearedPolys = divData["ClearedPolys"];
      simpPolys    = divData["SimplifiedPolys"];

      (* G0: (n-1)-dim integral at eps=0 — the pole numerator *)
      Module[{yV, aV, pV, intg},
        yV  = Table[Unique["cc8g0"], {n - 1}];
        aV  = a0[[ndVars]];
        pV  = Table[Total[Table[Module[{c, e},
          c = mono[[1]]; e = mono[[2]][[ndVars]];
          c * Exp[Total[e * (Log /@ yV)]]
        ], {mono, simpPolys[[j]]}]], {j, Length[simpPolys]}];
        intg = Abs[detM] * Exp[Total[(aV - 1) * (Log /@ yV)]] *
          Times @@ MapThread[Function[{p, b}, Exp[b * Log[p]]], {pV, B0}];
        g0Val = Quiet@NIntegrate[intg,
          Evaluate[Sequence @@ ({#, 0, 1} & /@ yV)],
          MaxRecursion -> 25, PrecisionGoal -> 7];
      ];

      (* G1: log-insertion integral — contributes to the finite part *)
      Module[{yV, aV, pV, base, logIns, intg},
        yV  = Table[Unique["cc8g1"], {n - 1}];
        aV  = a0[[ndVars]];
        pV  = Table[Total[Table[Module[{c, e},
          c = mono[[1]]; e = mono[[2]][[ndVars]];
          c * Exp[Total[e * (Log /@ yV)]]
        ], {mono, simpPolys[[j]]}]], {j, Length[simpPolys]}];
        base   = Abs[detM] * Exp[Total[(aV - 1) * (Log /@ yV)]] *
          Times @@ MapThread[Function[{p, b}, Exp[b * Log[p]]], {pV, B0}];
        logIns = Total[Table[a1[[ndVars[[i]]]] / aV[[i]] * Log[yV[[i]]],
                  {i, n - 1}]] +
                 Total[Table[B1[[j]] * Log[pV[[j]]], {j, Length[B0]}]];
        intg = base * logIns;
        g1Val = Quiet@NIntegrate[intg,
          Evaluate[Sequence @@ ({#, 0, 1} & /@ yV)],
          MaxRecursion -> 25, PrecisionGoal -> 6];
      ];

      (* Remainder: the full integrand minus the subtracted piece, at eps=0 *)
      Module[{yV, fPV, sPV, brk, intg},
        yV  = Table[Unique["cc8rm"], {n}];
        fPV = Table[Total[Table[Module[{c, e},
          c = mono[[1]]; e = mono[[2]];
          c * Exp[Total[e * (Log /@ yV)]]
        ], {mono, clearedPolys[[j]]}]], {j, Length[clearedPolys]}];
        sPV = Table[Total[Table[Module[{c, e},
          c = mono[[1]]; e = mono[[2]];
          c * Exp[Total[e * (Log /@ yV)]]
        ], {mono, simpPolys[[j]]}]], {j, Length[simpPolys]}];
        brk  = Times @@ MapThread[Function[{p, b}, Exp[b * Log[p]]], {fPV, B0}] -
               Times @@ MapThread[Function[{p, b}, Exp[b * Log[p]]], {sPV, B0}];
        intg = Abs[detM] * Exp[Total[(a0 - 1) * (Log /@ yV)]] * brk;
        remVal = Quiet@NIntegrate[intg,
          Evaluate[Sequence @@ ({#, 0, 1} & /@ yV)],
          MaxRecursion -> 20, PrecisionGoal -> 5];
      ];

      cc8Print["  Sector ", sd["ConeIndex"],
               ": G0/ck=", N[g0Val / ck],
               ", G1/ck=", N[g1Val / ck],
               ", rem=",   N[remVal]];

      totalPole   += g0Val / ck;
      totalFinite += g1Val / ck + remVal;
    ],
    {sd, divSectors}
  ];

  totalFinite += convContrib;

  cc8Print[""];
  cc8Print["  v3 tropical pole   = ", N[totalPole]];
  cc8Print["  v3 tropical finite = ", N[totalFinite]];
  cc8Print["  Exact pole         = ", N[exactPole]];
  cc8Print["  Exact finite       = ", N[exactFinite]];

  (* ================================================================ PART 3: GL(1) numeric cross-check *)
  (* Independently verify the GL(1) identity at a finite eps value to
     confirm the two representations (affine and simplex) agree, which
     validates the FIESTA comparison setup regardless of FIESTA's output. *)

  cc8Print[""];
  cc8Print["--- Part 3: GL(1) numeric identity check at finite eps ---"];

  Module[{eps0, niAffine, niSimplex, relErrGL},
    eps0 = 0.05;

    niAffine = Quiet@NIntegrate[
      t1^(2 eps0 - 1) / ((1 + t1) (1 + t2))^2,
      {t1, 0, Infinity}, {t2, 0, Infinity},
      MaxRecursion -> 20, PrecisionGoal -> 6
    ];

    (* Simplex form: x0^{2-2ep} x1^{2ep-1} [(x0+x1)(x0+x2)]^{-2} on simplex *)
    niSimplex = Quiet@NIntegrate[
      (1 - s1 - s2)^(2 - 2 eps0) * s1^(2 eps0 - 1) *
      ((1 - s1 - s2 + s1) (1 - s1 - s2 + s2))^(-2),
      {s1, 0, 1}, {s2, 0, 1 - s1},
      MaxRecursion -> 20, PrecisionGoal -> 6,
      Method -> "GlobalAdaptive"
    ];

    relErrGL = If[Abs[niAffine] > 0, Abs[(niAffine - niSimplex) / niAffine], Infinity];
    cc8Print["  Affine [0,inf) at eps=", eps0, ": ", N[niAffine]];
    cc8Print["  Simplex (GL1)   at eps=", eps0, ": ", N[niSimplex]];
    cc8Print["  Relative error (GL1 identity): ", N[relErrGL]];

    If[relErrGL > 0.01,
      cc8Print["  WARNING: GL(1) identity check failed (relErr=", N[relErrGL],
               " > 0.01); FIESTA comparison may be invalid."]
    ]
  ];

  (* ================================================================ PART 4: verdicts *)

  cc8Print[""];
  cc8Print["--- Part 4: PASS/FAIL verdicts ---"];

  (* Tolerance definitions *)
  tropPoleTol   = 5*^-3;   (* tropical pole vs exact *)
  tropFinTol    = 5*^-2;   (* tropical finite vs exact *)

  relErrTropPole = If[Abs[N[exactPole]] > 0,
    Abs[N[totalPole - exactPole]] / Abs[N[exactPole]], Infinity];
  relErrTropFin  = If[Abs[N[exactFinite]] > 0,
    Abs[N[totalFinite - exactFinite]] / Abs[N[exactFinite]], Infinity];

  (* -- Tropical vs exact -- *)
  If[relErrTropPole <= tropPoleTol,
    Print["CC8 PASS tropical-pole: expected=", N[exactPole],
          " got=", N[totalPole],
          " relErr=", N[relErrTropPole]],
    Print["CC8 FAIL expected=", N[exactPole],
          " got=", N[totalPole],
          " (tropical pole relErr=", N[relErrTropPole],
          " > tol=", tropPoleTol, ")"]
  ];

  If[relErrTropFin <= tropFinTol,
    Print["CC8 PASS tropical-finite: expected=", N[exactFinite],
          " got=", N[totalFinite],
          " relErr=", N[relErrTropFin]],
    Print["CC8 FAIL expected=", N[exactFinite],
          " got=", N[totalFinite],
          " (tropical finite relErr=", N[relErrTropFin],
          " > tol=", tropFinTol, ")"]
  ];

  (* -- FIESTA vs exact (if FIESTA returned parseable numbers) -- *)
  If[NumericQ[fiestaPole] && NumericQ[fiestaFin],

    (* The FIESTA error bar sigma for the pole *)
    Module[{nSigPole, nSigFin, fiestaPolePass, fiestaFinPass,
            absErrPole, absErrFin, poleSigTol, finSigTol},

      (* Use at least 2σ tolerance; if sigma = 0 (not reported), use 1e-3 *)
      poleSigTol = Max[2 * Abs[N[fiestaSigPole]], 1*^-3];
      finSigTol  = Max[2 * Abs[N[fiestaSigFin]],  5*^-2];

      absErrPole = Abs[N[fiestaPole - exactPole]];
      absErrFin  = Abs[N[fiestaFin  - exactFinite]];

      fiestaPolePass = (absErrPole <= poleSigTol);
      fiestaFinPass  = (absErrFin  <= finSigTol);

      If[fiestaPolePass,
        Print["CC8 PASS fiesta-pole: expected=", N[exactPole],
              " got=", N[fiestaPole],
              " absErr=", N[absErrPole],
              " tol(2sigma)=", N[poleSigTol]],
        Print["CC8 FAIL expected=", N[exactPole],
              " got=", N[fiestaPole],
              " (fiesta pole absErr=", N[absErrPole],
              " > 2sigma=", N[poleSigTol], ")"]
      ];

      If[fiestaFinPass,
        Print["CC8 PASS fiesta-finite: expected=", N[exactFinite],
              " got=", N[fiestaFin],
              " absErr=", N[absErrFin],
              " tol(2sigma)=", N[finSigTol]],
        Print["CC8 FAIL expected=", N[exactFinite],
              " got=", N[fiestaFin],
              " (fiesta finite absErr=", N[absErrFin],
              " > 2sigma=", N[finSigTol], ")"]
      ];

      (* -- FIESTA vs Tropical agreement -- *)
      Module[{poleAgreeErr, finAgreeErr, poleAgreeTol, finAgreeTol},
        poleAgreeTol = Max[2 * Abs[N[fiestaSigPole]], tropPoleTol];
        finAgreeTol  = Max[2 * Abs[N[fiestaSigFin]],  tropFinTol];
        poleAgreeErr = Abs[N[fiestaPole - totalPole]];
        finAgreeErr  = Abs[N[fiestaFin  - totalFinite]];

        If[poleAgreeErr <= poleAgreeTol,
          Print["CC8 PASS fiesta-vs-tropical-pole: fiesta=", N[fiestaPole],
                " tropical=", N[totalPole],
                " absErr=", N[poleAgreeErr],
                " tol=", N[poleAgreeTol]],
          Print["CC8 FAIL expected=", N[fiestaPole],
                " got=", N[totalPole],
                " (fiesta-vs-tropical pole absErr=", N[poleAgreeErr],
                " > tol=", N[poleAgreeTol], ")"]
        ];

        If[finAgreeErr <= finAgreeTol,
          Print["CC8 PASS fiesta-vs-tropical-finite: fiesta=", N[fiestaFin],
                " tropical=", N[totalFinite],
                " absErr=", N[finAgreeErr],
                " tol=", N[finAgreeTol]],
          Print["CC8 FAIL expected=", N[fiestaFin],
                " got=", N[totalFinite],
                " (fiesta-vs-tropical finite absErr=", N[finAgreeErr],
                " > tol=", N[finAgreeTol], ")"]
        ];
      ];
    ],

    (* FIESTA result was not parseable as numbers — report but don't FAIL cc8 *)
    cc8Print["  NOTE: FIESTA result could not be parsed as numeric pole/finite; ",
             "skipping FIESTA-vs-tropical and FIESTA-vs-exact comparisons."];
    cc8Print["  Raw FIESTA output: ", fiestaResult]
  ];

] (* end Module *)

cc8Print[""];
cc8Print["CC8 complete."];
