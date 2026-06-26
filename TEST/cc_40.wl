(* ============================================================================
   TEST/cc_40.wl  --  Cross-check #40 (plan.md §8.3)
   "Lift identity -- complex & negative & multi-rule (symbolic)"

   Plan.md §8.3 #40 (F3a/F3b uplift / exactness invariant / §6.2):
       "z->z0 substitution ≡ original polynomials, exact (rel-tol for floats),
        for +, −, complex, and 2-rule incommensurate lifts."

   Strategy (Tier 1: WL + g++ only, no CUBA):
     For each of four cases, call LiftCoefficients, then:
       (1) substitute auxVar -> z0 into the lifted polynomial;
       (2) check Expand[subbed] === Expand[original] or PossibleZeroQ of diff;
       (3) for float z0 (complex/incommensurate), check numeric relative
           residual at two independent sample points < 1e-8;
       (4) check z0^k === C exactly (the anchor identity);
       (5) check residual c = C / z0^k exactly (the residual identity).

   Cases:
     C40-A  Positive real coefficient:  C = 10^6,  k=2,  z0 = 1000 (exact Int).
     C40-B  Negative real coefficient:  C = -10^6, k=2,  z0 = 1000; residual = -1.
     C40-C  Complex coefficient:        C = 3 + 4*I, k=1, z0 = 5; residual = (3+4I)/5.
     C40-D  Two-rule incommensurate:    C1 = 10^6 (exp {2,0}), C2 = 10^-3 (exp {0,1}),
                k1=2, k2=1, primary = C1; residuals c1=-1 (test with neg) and c2 exact.

   PASS criterion (per plan.md):
     Each sub-case prints "CC40 PASS [C40-X] ..." on success.
     Any failure prints "CC40 FAIL expected=<e> got=<g>".
     The final summary prints "CC40 PASS ..." iff ALL sub-cases pass.

   Ported from:
     OLD_CODE/TROPICAL_MONTE_CARLO/TMCv2_BUG_LOG.md   (BUG 2, liftidentity logic)
     OLD_CODE/TROPICAL_MONTE_CARLO/lift_error_log.md   (pre-session liftidentity fix)
     plan.md §6.2 (residual formula, anchor identity z0^k===C)

   Run:  wolframscript -file TEST/cc_40.wl
   (from TROPICAL_MONTE_CARLO3 root or anywhere; uses absolute paths)
   ============================================================================ *)

(* ---- 0. Locate root and load packages via absolute paths ---- *)
$cc40Root = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
$cc40Fan  = FileNameJoin[{$cc40Root, "tropical_fan.wl"}];
$cc40Eval = FileNameJoin[{$cc40Root, "tropical_eval.wl"}];

If[!FileExistsQ[$cc40Fan],
  Print["CC40 FAIL expected=tropical_fan.wl got=file-not-found at ", $cc40Fan];
  Quit[1]
];
If[!FileExistsQ[$cc40Eval],
  Print["CC40 FAIL expected=tropical_eval.wl got=file-not-found at ", $cc40Eval];
  Quit[1]
];

Get[$cc40Fan];
Get[$cc40Eval];

Print["CC40: packages loaded from ", $cc40Root];
Print[];

(* ---- 1. Shared helpers ---- *)

(* Numeric relative residual at two sample points (avoids single-point accidents).
   Returns the maximum over two probes. *)
numericRelResid[diff_, original_, liftedSubbed_, vars_] :=
Module[{subs, maxResid = 0},
  Do[
    Module[{sub, num, den, r},
      sub = Thread[vars -> Table[1.7 + seed*0.37 I, {i, Length[vars]}]];
      num = Abs[N[diff /. sub]];
      den = Abs[N[original /. sub]] + Abs[N[liftedSubbed /. sub]] + 1;
      r = If[den > 0, num / den, num];
      maxResid = Max[maxResid, r]
    ],
    {seed, {1, 2}}
  ];
  maxResid
];

(* Core identity verifier: z->z0 in liftedPoly must recover originalPoly.
   Returns True (pass) or False (fail), and prints diagnostics. *)
checkLiftIdentity[label_String, liftedPoly_, originalPoly_, auxVar_, z0_,
                  vars_List] :=
Module[{subbed, diff, relR},
  subbed = Expand[liftedPoly /. auxVar -> z0];
  diff   = Simplify[subbed - Expand[originalPoly]];
  Which[
    TrueQ[diff === 0],
      Print["  identity: exact symbolic match"];
      True,
    TrueQ[PossibleZeroQ[diff]],
      Print["  identity: PossibleZeroQ (structurally zero)"];
      True,
    True,
      relR = numericRelResid[diff, Expand[originalPoly], subbed, vars];
      Print["  identity numeric relResid = ", relR];
      If[NumericQ[relR] && relR < 10.^-8,
        Print["  identity: numeric residual < 1e-8 (PASS)"];
        True,
        Print["  identity: numeric residual = ", relR, " >= 1e-8 (FAIL)"];
        False
      ]
  ]
];

(* Check anchor identity z0^{k_primary} === |C_primary| for the PRIMARY rule only.
   z0 = |C_primary|^(1/k_primary) so z0^{k_primary} = |C_primary| exactly.
   For secondary rules z0 is shared and z0^{k_i} need not equal |C_i| in general
   -- the difference is absorbed into residual c_i = C_i / z0^{k_i} (plan.md §6.2).
   primaryIdx: 1-based index of the primary rule (default 1). *)
checkAnchorIdentity[label_String, z0_, liftRules_List, ruleCoeffs_List,
                    primaryIdx_Integer : 1] :=
Module[{allOK = True},
  Do[
    If[i =!= primaryIdx,
      Print["  anchor z0^k=|C| rule ", i,
            ": secondary rule -- anchor applies only to primary (rule ", primaryIdx, ")"];
      Continue[]
    ];
    Module[{k = liftRules[[i]]["k"], C = ruleCoeffs[[i]], lhs, rhs, diff, relR},
      lhs = Simplify[z0^k];
      rhs = Abs[C];
      diff = Simplify[lhs - rhs];
      Which[
        TrueQ[diff === 0] || TrueQ[PossibleZeroQ[diff]],
          Print["  anchor z0^k=|C| rule ", i, " (primary): exact"],
        True,
          Module[{num = Abs[N[diff]], den = Abs[N[lhs]] + Abs[N[rhs]] + 1},
            relR = num / den;
            If[NumericQ[relR] && relR < 10.^-10,
              Print["  anchor z0^k=|C| rule ", i, " (primary): numeric relResid=", relR, " OK"],
              Print["  anchor z0^k=|C| rule ", i, " (primary): FAIL relResid=", relR];
              allOK = False
            ]
          ]
      ]
    ],
    {i, Length[liftRules]}
  ];
  allOK
];

(* Check residual c_i = C_i / z0^{k_i} exactly. *)
checkResidualIdentity[label_String, z0_, liftRules_List, ruleCoeffs_List,
                      residuals_List] :=
Module[{allOK = True},
  Do[
    Module[{k = liftRules[[i]]["k"], C = ruleCoeffs[[i]],
            c = residuals[[i]], diff, relR},
      diff = Simplify[c - C / z0^k];
      Which[
        TrueQ[diff === 0] || TrueQ[PossibleZeroQ[diff]],
          Print["  residual c=C/z0^k rule ", i, ": exact"],
        True,
          Module[{num = Abs[N[diff]],
                  den = Abs[N[c]] + Abs[N[C / z0^k]] + 1},
            relR = num / den;
            If[NumericQ[relR] && relR < 10.^-10,
              Print["  residual c=C/z0^k rule ", i, ": numeric relResid=", relR, " OK"],
              Print["  residual c=C/z0^k rule ", i, ": FAIL relResid=", relR];
              allOK = False
            ]
          ]
      ]
    ],
    {i, Length[liftRules]}
  ];
  allOK
];

(* Run one sub-case end-to-end.  Returns True/False. *)
runSubCase[label_String, spec_Association, liftRules_List] :=
Module[{lcRes, liftedSpec, liftData, z0, auxVar, vars, liftedPoly, origPoly,
        polyIdx, ruleCoeffs, residuals, identOK, anchorOK, residOK,
        pass = True},
  Print["--------------------------------------------------"];
  Print["CC40 sub-case ", label];

  lcRes = Quiet @ LiftCoefficients[spec, liftRules];
  If[!AssociationQ[lcRes],
    Print["  LiftCoefficients returned: ", lcRes];
    Print["CC40 FAIL expected=Association got=", lcRes, " [", label, "]"];
    Return[False]
  ];

  liftedSpec = lcRes["LiftedSpec"];
  liftData   = lcRes["LiftData"];
  z0         = liftData["z0"];
  auxVar     = liftData["AuxVariable"];
  vars       = spec["Variables"];
  residuals  = liftData["Residuals"];

  Print["  z0 = ", z0, "  (exact: ", Head[z0] =!= Real, ")"];
  Print["  residuals = ", residuals];

  (* Collect ruleCoeffs from parsedPolys so the anchor/residual checks work *)
  Module[{parsedPolys = ParsePolynomial[#, vars] & /@ spec["Polynomials"]},
    ruleCoeffs = Table[
      Module[{r = liftRules[[i]], j0, alpha, matchPos},
        j0       = r["PolyIndex"];
        alpha    = r["ExponentVector"];
        matchPos = Position[parsedPolys[[j0, All, 2]], alpha];
        If[matchPos === {}, Missing["notfound", i],
           parsedPolys[[j0, matchPos[[1, 1]], 1]]]
      ],
      {i, Length[liftRules]}
    ]
  ];

  (* Check polynomial identity for each polynomial in the spec *)
  Do[
    liftedPoly = liftedSpec["Polynomials"][[j]];
    origPoly   = spec["Polynomials"][[j]];
    Print["  poly ", j, " identity:"];
    identOK = checkLiftIdentity[label, liftedPoly, origPoly, auxVar, z0, vars];
    If[!identOK,
      Print["CC40 FAIL expected=identity-holds got=residual-too-large [", label,
            " poly ", j, "]"];
      pass = False
    ],
    {j, Length[spec["Polynomials"]]}
  ];

  (* Determine primary rule index (argmax |Log[|C||) -- mirrors LiftCoefficients logic *)
  Module[{logMags = Abs[Log[Abs[N[#]]]] & /@ ruleCoeffs,
          primIdx},
    primIdx = First @ Ordering[logMags, -1];

    (* Anchor identity z0^{k_primary} = |C_primary| *)
    Print["  Anchor identities (primary rule=", primIdx, "):"];
    anchorOK = checkAnchorIdentity[label, z0, liftRules, ruleCoeffs, primIdx]
  ];
  If[!anchorOK,
    Print["CC40 FAIL expected=z0^k=|C| got=mismatch [", label, "]"];
    pass = False
  ];

  (* Residual identity c = C / z0^k *)
  Print["  Residual identities:"];
  residOK = checkResidualIdentity[label, z0, liftRules, ruleCoeffs, residuals];
  If[!residOK,
    Print["CC40 FAIL expected=c=C/z0^k got=mismatch [", label, "]"];
    pass = False
  ];

  If[pass,
    Print["CC40 PASS [", label, "] lift identity, anchor, and residual all exact"],
    Print["CC40 FAIL [", label, "] one or more identity checks failed"]
  ];
  Print[];
  pass
];

(* ---- 2. Sub-case C40-A: positive real coefficient ---- *)
(* C = 10^6, k=2 => z0 = (10^6)^(1/2) = 1000 exactly (IntegerQ).
   P = 1 + 10^6 * x[1]^2 + x[2]^2, B=-2.
   Lifted: P -> 1 + z^2 * x[1]^2 + x[2]^2 (residual c = 10^6 / 1000^2 = 1). *)
passA = runSubCase["C40-A (positive real C=10^6 k=2)",
  <|"Polynomials"        -> {1 + 10^6 * x[1]^2 + x[2]^2},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents"-> {-2},
    "Variables"          -> {x[1], x[2]},
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None|>,
  {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 2|>}
];

(* ---- 3. Sub-case C40-B: negative real coefficient ---- *)
(* C = -10^6, k=2 => z0 = |C|^(1/2) = 1000; residual c = -10^6 / 1000^2 = -1.
   P = 1 - 10^6 * x[1]^2 + x[2]^2, B=-2.
   The sign/phase goes into the residual (plan.md §6.2 explicitly). *)
passB = runSubCase["C40-B (negative real C=-10^6 k=2)",
  <|"Polynomials"        -> {1 - 10^6 * x[1]^2 + x[2]^2},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents"-> {-2},
    "Variables"          -> {x[1], x[2]},
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None|>,
  {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 2|>}
];

(* ---- 4. Sub-case C40-C: complex coefficient ---- *)
(* C = (3 + 4*I) (|C| = 5), k=1 => z0 = 5; residual c = (3+4I)/5.
   P = 1 + (3+4I) * x[1] + x[2], B=-1.
   This is the BUG-2-regime: complex C with exact Gaussian integer magnitude. *)
passC = runSubCase["C40-C (complex C=3+4I k=1)",
  <|"Polynomials"        -> {1 + (3 + 4*I) * x[1] + x[2]},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents"-> {-1},
    "Variables"          -> {x[1], x[2]},
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None|>,
  {<|"PolyIndex" -> 1, "ExponentVector" -> {1, 0}, "k" -> 1|>}
];

(* ---- 5. Sub-case C40-D: two-rule incommensurate lift ---- *)
(* Two monomials lifted simultaneously:
     rule 1: C1 = -10^6, alpha1 = {2,0}, k1 = 2  => primary (larger |log|C||)
     rule 2: C2 = 10^-3, alpha2 = {0,1}, k2 = 1
   primary z0 = |C1|^(1/2) = 1000.
   residuals:  c1 = -10^6 / 1000^2 = -1
               c2 = 10^-3 / 1000^1 = 10^-6
   P = 1 - 10^6 * x[1]^2 + 10^-3 * x[2] + x[1]*x[2], B=-2.
   This tests the incommensurate k case documented in plan.md §6.2 and the
   multi-rule residual identity.  "Incommensurate" means k1 != k2. *)
passD = runSubCase["C40-D (two-rule incommensurate: C1=-10^6 k=2, C2=10^-3 k=1)",
  <|"Polynomials"        -> {1 - 10^6 * x[1]^2 + 10^-3 * x[2] + x[1]*x[2]},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents"-> {-2},
    "Variables"          -> {x[1], x[2]},
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None|>,
  {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 2|>,
   <|"PolyIndex" -> 1, "ExponentVector" -> {0, 1}, "k" -> 1|>}
];

(* ---- 6. Summary ---- *)
Print["=================================================="];
If[TrueQ[passA] && TrueQ[passB] && TrueQ[passC] && TrueQ[passD],
  Print["CC40 PASS  Lift identity exact for: positive(A), negative(B), ",
        "complex(C), two-rule-incommensurate(D); z0^k=|C| and c=C/z0^k ",
        "all hold symbolically."],
  Print["CC40 FAIL  One or more sub-cases failed.  ",
        "A=", passA, " B=", passB, " C=", passC, " D=", passD]
];
Print["=================================================="];
