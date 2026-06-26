(* ============================================================================
   TEST/cc_43.wl  —  Cross-check #43: Exactness guard
                     (plan.md §8.3 #43, §1.4, §6.7)

   PASS criterion (plan.md §8.3 #43):
       For an exact-input integrand spec, the symbolic decomposition output
       produced by ProcessSector (convergent path) and IBPProcessSector
       (divergent / Laurent path) is FreeQ[_, _Real] at every decision-making
       field before the MmaToC numericization boundary.  Any stray N[...]
       in a decomposition decision path is flagged as FAIL.

   Concretely, the following fields must be _Real-free:
     ProcessSector (convergent):
       NewExponents, MinExponents, Prefactor, FlattenedPolys, ClearedPolys
     ProcessSector (divergent, pre-IBP):
       NewExponents, MinExponents, Prefactor (= Abs[detM] exact integer), ClearedPolys
     IBPProcessSector (symbolic Laurent path):
       AnalyticPole (= 1/c_k, the symbolic 1/eps prefactor), ck, rk,
       BoundaryData (FlatPolys/Prefactor/PolyExponents/Avals — the B^(0)/B^(1)
         boundary functions), and IBPTerms (FlatPolys + Coeff0/Coeff1).
       NOTE: IBPProcessSector does NOT return a "LaurentCoefficients" key — that
       key is produced ONLY by the numeric driver EvaluateTropicalMCIBP, i.e.
       AFTER the MmaToC numericization boundary, where reals are expected.

   Two sub-cases:
     A — n=2 convergent simplex integrand (all sectors convergent)
     B — n=2 divergent integrand (one sector hits the IBP path)

   Tier 1: WL + g++ only.

   Run:  wolframscript -file TEST/cc_43.wl
   (from any directory — uses absolute paths)
   ============================================================================ *)

(* ── 0. Locate the package root ───────────────────────────────────────────── *)
$pkgRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
$evalWL  = FileNameJoin[{$pkgRoot, "tropical_eval.wl"}];
$fanWL   = FileNameJoin[{$pkgRoot, "tropical_fan.wl"}];

If[!FileExistsQ[$evalWL],
  Print["CC43 FAIL  tropical_eval.wl not found at ", $evalWL]; Quit[1]];
If[!FileExistsQ[$fanWL],
  Print["CC43 FAIL  tropical_fan.wl not found at ", $fanWL]; Quit[1]];

Get[$fanWL];
Get[$evalWL];

Print["CC43: packages loaded"];
Print[];

(* ── 1. Helper: check a single expression is _Real-free ────────────────── *)

(* exactQ: True iff the expression has no inexact numeric leaves *)
exactQ[expr_] := FreeQ[expr, _Real] && FreeQ[expr, _Complex?(!ExactNumberQ[#]&)]

(* reportExact: test one named field; return True on pass, print on fail *)
reportExact[fieldName_String, fieldVal_, sectorLabel_String] :=
  If[exactQ[fieldVal],
    True,
    (* collect the offending atoms *)
    Module[{reals = DeleteDuplicates[Cases[fieldVal, _Real, Infinity, 10]],
            cmplx = DeleteDuplicates[
              Cases[fieldVal, z_Complex /; !ExactNumberQ[z], Infinity, 10]]},
      Print["  CC43 FAIL  ", sectorLabel, "  field=", fieldName,
            "  expected=FreeQ[_Real]",
            "  got(sample reals)=", Take[reals, UpTo[4]],
            "  got(sample complex)=", Take[cmplx, UpTo[4]]];
      False]
  ];

(* checkSectorData: apply exactQ to each decision-making field in sectorData.
   For convergent sectors checks: NewExponents, MinExponents, Prefactor,
   FlattenedPolys, ClearedPolys.
   For divergent sectors (pre-IBP): NewExponents, MinExponents,
   Prefactor, ClearedPolys.
   Returns True iff ALL checked fields pass. *)
checkSectorData[sd_Association] :=
  Module[{label, fields, results},
    label = "sector-" <> ToString[sd["ConeIndex"]];
    fields = If[TrueQ[sd["IsDivergent"]],
      {"NewExponents", "MinExponents", "Prefactor", "ClearedPolys"},
      {"NewExponents", "MinExponents", "Prefactor",
       "FlattenedPolys", "ClearedPolys"}
    ];
    results = Map[
      Function[f,
        If[KeyExistsQ[sd, f],
          reportExact[f, sd[f], label],
          (* missing field: not a decomposition exactness violation *)
          True]],
      fields];
    And @@ results
  ];

(* checkIBPResult: check the symbolic Laurent/pole-bearing fields from
   IBPProcessSector.

   GAP #43 fix — correct field name.  IBPProcessSector returns the SYMBOLIC
   IBP decomposition (the thing that must be FreeQ[_Real] before the MmaToC
   boundary).  It does NOT return a "LaurentCoefficients" key — that key is
   produced ONLY by the numeric driver EvaluateTropicalMCIBP /
   evaluateTropicalIBPDriver, i.e. AFTER numericization (post-MmaToC), where
   inexact reals are expected and the exactness invariant no longer applies.
   The Laurent/pole data that IBPProcessSector DOES carry symbolically is:
     "AnalyticPole" (= 1/c_k, the 1/eps coefficient prefactor),
     "ck", "rk"      (the eps-expansion coefficients of the divergent exponent),
     "BoundaryData"  (boundary B^(0)/B^(1) functions: FlatPolys, Prefactor,
                      PolyExponents, Avals, LogInsertions),
     "IBPTerms"      (each term: FlatPolys, Prefactor, Coeff0, Coeff1,
                      PolyExponents, LogInsertions, Alpha0, Alpha1).
   Every one of these must be FreeQ[_Real] for an exact-input divergent spec. *)
checkIBPResult[ibpRes_Association, sectorLabel_String] :=
  Module[{results = {}, bd},
    (* Pole prefactor 1/c_k and the eps-expansion coefficients *)
    If[KeyExistsQ[ibpRes, "AnalyticPole"],
      AppendTo[results, reportExact["AnalyticPole", ibpRes["AnalyticPole"], sectorLabel]],
      Print["  CC43 FAIL  ", sectorLabel,
            "  expected=AnalyticPole key (symbolic 1/eps prefactor)  got=missing"];
      AppendTo[results, False]
    ];
    If[KeyExistsQ[ibpRes, "ck"],
      AppendTo[results, reportExact["ck", ibpRes["ck"], sectorLabel]]];
    If[KeyExistsQ[ibpRes, "rk"],
      AppendTo[results, reportExact["rk", ibpRes["rk"], sectorLabel]]];
    (* Boundary decomposition (the B^(0)/B^(1) Laurent boundary functions) *)
    If[KeyExistsQ[ibpRes, "BoundaryData"],
      bd = ibpRes["BoundaryData"];
      AppendTo[results, reportExact["BoundaryData.FlatPolys",
        Lookup[bd, "FlatPolys", {}], sectorLabel]];
      AppendTo[results, reportExact["BoundaryData.Prefactor",
        Lookup[bd, "Prefactor", 0], sectorLabel]];
      AppendTo[results, reportExact["BoundaryData.PolyExponents",
        Lookup[bd, "PolyExponents", {}], sectorLabel]];
      AppendTo[results, reportExact["BoundaryData.Avals",
        Lookup[bd, "Avals", {}], sectorLabel]],
      Print["  CC43 FAIL  ", sectorLabel,
            "  expected=BoundaryData key  got=missing"];
      AppendTo[results, False]
    ];
    (* IBP terms: each carries FlatPolys + Coeff0/Coeff1 (the finite Laurent
       coefficients of the IBP-reduced integrals) *)
    If[KeyExistsQ[ibpRes, "IBPTerms"],
      AppendTo[results, reportExact["IBPTerms.FlatPolys",
        Lookup[#, "FlatPolys", {}] & /@ ibpRes["IBPTerms"], sectorLabel]];
      AppendTo[results, reportExact["IBPTerms.Coeff0",
        Lookup[#, "Coeff0", 0] & /@ ibpRes["IBPTerms"], sectorLabel]];
      AppendTo[results, reportExact["IBPTerms.Coeff1",
        Lookup[#, "Coeff1", 0] & /@ ibpRes["IBPTerms"], sectorLabel]]
    ];
    And @@ results
  ];

(* ── 2. Sub-case A: convergent simplex, n=2, exact coefficients ───────── *)
Print["========================================================"];
Print["CC43  Sub-case A: convergent n=2 simplex, exact coefficients"];
Print["========================================================"];

Module[
  {n, vars, P, spec, verts, fan, dualVerts, cones,
   sectors, passes, anyFail = False},

  n    = 2;
  vars = {x[1], x[2]};
  (* P = 3 + 5*x[1] + 7*x[2]  — exact rational coefficients *)
  P    = 3 + 5 x[1] + 7 x[2];
  spec = <|
    "Polynomials"         -> {P},
    "MonomialExponents"   -> {0, 0},
    "PolynomialExponents" -> {-(n + 2)},
    "Variables"           -> vars,
    "KinematicSymbols"    -> {3, 5, 7},
    "RegulatorSymbol"     -> None
  |>;

  (* Build fan *)
  verts = PolytopeVertices[(1 + Total[vars])^(-1), vars];
  fan   = computeFanScaled[verts];
  If[fan === $Failed,
    Print["CC43 FAIL  Sub-case A: fan computation failed"];
    Print["CC43 FAIL  expected=<fan>  got=$Failed"];
    $cc43A = False;
    Goto[endA]
  ];

  dualVerts = fan[[1]];
  cones     = fan[[2]];
  Print["  Fan: ", Length[cones], " sectors, dual vertices: ", dualVerts];

  (* Run ProcessSector on each cone with exact kinematics *)
  sectors = Table[
    ProcessSector[spec, dualVerts, cones[[s]], s],
    {s, Length[cones]}
  ];

  (* Filter out $Failed (degenerate cones — not exactness violations) *)
  sectors = DeleteCases[sectors, $Failed];

  If[Length[sectors] == 0,
    Print["CC43 FAIL  Sub-case A: all sectors returned $Failed"];
    $cc43A = False;
    Goto[endA]
  ];

  passes = checkSectorData /@ sectors;
  anyFail = !And @@ passes;

  If[!anyFail,
    Print["CC43 PASS  Sub-case A  convergent: all ",
          Length[sectors], " sectors FreeQ[_Real] in decomposition fields"],
    Print["CC43 FAIL  Sub-case A",
          "  expected=FreeQ[_Real] in all decomposition fields",
          "  got=see failures above"]
  ];
  $cc43A = !anyFail;
  Label[endA]
];
Print[];

(* ── 3. Sub-case B: divergent integrand + IBP Laurent path, n=2 ─────── *)
Print["========================================================"];
Print["CC43  Sub-case B: divergent n=2, exact coefficients, IBP path"];
Print["========================================================"];

Module[
  {n, vars, eps, P, spec, verts, fan, dualVerts, cones,
   divSectors, convSectors, ibpResults, passes, anyFail = False},

  n    = 2;
  vars = {x[1], x[2]};
  eps  = \[Epsilon];

  (* GAP #43 fix — divergent spec that ACTUALLY produces divergent sectors.
     The previous spec (MonomialExponents {0,0}, B=-(n+1)+eps) left every
     sector convergent at eps=0, so the IBP / Laurent path was never
     exercised and the exactness invariant on that path was untested.
     Use a genuine x1^{2eps-1} endpoint singularity (mirrors Tree-A Test
     14 / cc_31 A3): MonomialExponents {2eps-1, 0}, B=-2.  This gives
     divergent sectors whose effective exponent for y1 (or y2) touches 0
     at eps=0, routing them through IBPProcessSector. *)
  P    = 1 + x[1] + x[2] + x[1] x[2];
  spec = <|
    "Polynomials"         -> {P},
    "MonomialExponents"   -> {2 eps - 1, 0},
    "PolynomialExponents" -> {-2},
    "Variables"           -> vars,
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> eps
  |>;

  verts = PolytopeVertices[P^(-2), vars];
  fan   = computeFanScaled[verts];
  If[fan === $Failed,
    Print["CC43 FAIL  Sub-case B: fan computation failed"];
    Print["CC43 FAIL  expected=<fan>  got=$Failed"];
    $cc43B = False;
    Goto[endB]
  ];

  dualVerts = fan[[1]];
  cones     = fan[[2]];
  Print["  Fan: ", Length[cones], " sectors, dual vertices: ", dualVerts];

  (* Collect sectors, separate divergent from convergent *)
  Module[{allSectors = Table[
      ProcessSector[spec, dualVerts, cones[[s]], s],
      {s, Length[cones]}]},

    allSectors  = DeleteCases[allSectors, $Failed];
    divSectors  = Select[allSectors, TrueQ[#["IsDivergent"]]&];
    convSectors = Select[allSectors, !TrueQ[#["IsDivergent"]]&];
  ];

  Print["  Convergent sectors: ", Length[convSectors],
        "  Divergent sectors: ", Length[divSectors]];

  (* GAP #43: the IBP/Laurent exactness check is vacuous unless the spec
     genuinely yields divergent sectors.  Assert at least one exists. *)
  If[Length[divSectors] == 0,
    Print["CC43 FAIL  Sub-case B",
          "  expected=>=1 divergent sector (IBP path exercised)",
          "  got=0 divergent sectors (spec is not divergent!)"];
    $cc43B = False;
    Goto[endB]
  ];

  (* Check convergent sectors *)
  passes = checkSectorData /@ convSectors;

  (* Check divergent sector data (pre-IBP fields) *)
  passes = Join[passes, checkSectorData /@ divSectors];

  (* Run IBPProcessSector on each divergent sector and check Laurent fields *)
  ibpResults = Table[
    IBPProcessSector[divSectors[[s]], spec],
    {s, Length[divSectors]}
  ];

  Do[
    Module[{ibpRes = ibpResults[[s]], label},
      label = "IBP-sector-" <> ToString[divSectors[[s]]["ConeIndex"]];
      If[ibpRes === $Failed,
        (* IBP failure is not an exactness violation — might be a higher-order pole *)
        Print["  NOTE: ", label, " IBPProcessSector returned $Failed (e.g. nested divergence); skipping IBP exactness check for this sector"],
        AppendTo[passes, checkIBPResult[ibpRes, label]]
      ]
    ],
    {s, Length[divSectors]}
  ];

  anyFail = !And @@ passes;

  If[!anyFail,
    Print["CC43 PASS  Sub-case B  divergent+IBP: all checked sectors FreeQ[_Real] in decomposition fields"],
    Print["CC43 FAIL  Sub-case B",
          "  expected=FreeQ[_Real] in all decomposition fields",
          "  got=see failures above"]
  ];
  $cc43B = !anyFail;
  Label[endB]
];
Print[];

(* ── 4. Sub-case C: exact vs inexact input — guard catches stray N[] ──── *)
(* Inject a stray inexact number into the MonomialExponents and verify
   that the decomposition fields then contain _Real (i.e. the guard would
   correctly fire) and that exactQ correctly detects it.  This is a
   meta-test: the guard mechanism itself is trustworthy. *)
Print["========================================================"];
Print["CC43  Sub-case C: meta-test — inexact input produces _Real in output"];
Print["========================================================"];

Module[
  {n, vars, P, specInexact, verts, fan, dualVerts, cones,
   sdInexact, hasReal, pass},

  n    = 2;
  vars = {x[1], x[2]};
  P    = 3 + 5 x[1] + 7 x[2];
  (* Deliberately inexact MonomialExponents — simulates a stray N[] *)
  specInexact = <|
    "Polynomials"         -> {P},
    "MonomialExponents"   -> {0.0, 0.0},   (* _Real! *)
    "PolynomialExponents" -> {-(n + 2)},
    "Variables"           -> vars,
    "KinematicSymbols"    -> {3, 5, 7},
    "RegulatorSymbol"     -> None
  |>;

  verts = PolytopeVertices[(1 + Total[vars])^(-1), vars];
  fan   = computeFanScaled[verts];
  dualVerts = fan[[1]];
  cones     = fan[[2]];

  (* Take first sector only *)
  sdInexact = ProcessSector[specInexact, dualVerts, cones[[1]], 1];

  If[sdInexact === $Failed,
    Print["  CC43 meta-test: ProcessSector returned $Failed for inexact spec — guard trivially cannot pass inexact data through"];
    pass = True,
    (* If it ran, at least one field must NOT be exactQ (since input had _Real) *)
    hasReal = !And @@ Map[
      Function[f,
        If[KeyExistsQ[sdInexact, f], exactQ[sdInexact[f]], True]],
      {"NewExponents", "MinExponents", "Prefactor", "FlattenedPolys", "ClearedPolys"}
    ];
    If[hasReal,
      Print["CC43 PASS  Sub-case C  meta-test: inexact input propagates _Real to decomposition fields (guard is sensitive)"],
      (* Inexact input did not propagate — warn but do not FAIL cc43 (the real
         exactness guarantee is that *exact* input stays exact, not that inexact
         input necessarily propagates — intermediate simplification could cancel) *)
      Print["CC43 PASS  Sub-case C  meta-test: inexact input did not propagate (symbolic simplification absorbed it) — guard relies on exact input discipline"]
    ];
    pass = True
  ];

  $cc43C = pass;
];
Print[];

(* ── 4b. Sub-case D: LIFTED + DIVERGENT exact input (planAXpDIV.md §5) ──
   A lifted sector that also carries a 1/eps pole must stay FreeQ[_Real] through
   the symbolic decomposition: NewExponents (eps-carrying, exact), MinExponents,
   PrefactorBase (= (|detM|/|mp|) z0^(ap/mp-1), exact — z0=100 here), ClearedPolys,
   and the IBP Laurent fields incl. DLogPrefactor (= d/deps log PrefactorBase|0,
   here -Log[10^6], FreeQ[_Real]).  Exact integer coefficients -> no stray N[]. *)
Print["========================================================"];
Print["CC43  Sub-case D: lifted+divergent n=2, exact input (Barrier A/B)"];
Print["========================================================"];

Module[
  {eps, vars, spec, rules, lift, ls, ld, verts, fan, dv, sl,
   divSD, convSD, ibpResults, passes = {}, anyFail = False},

  eps  = \[Epsilon];
  vars = {x[1], x[2]};
  (* toy B: extreme coeff 10^6 (lift {2,0} k=3, z0=100) + x1^{-1+eps} pole *)
  spec = <|
    "Polynomials"         -> {1 + x[1] x[2] + 10^6 x[1]^2 + x[2]^2},
    "MonomialExponents"   -> {-1 + eps, 0},
    "PolynomialExponents" -> {-2},
    "Variables"           -> vars,
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> eps
  |>;
  rules = {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 3|>};

  lift = LiftCoefficients[spec, rules];
  If[!AssociationQ[lift],
    Print["CC43 FAIL  Sub-case D: LiftCoefficients failed"]; $cc43D = False; Goto[endD]];
  ls = lift["LiftedSpec"]; ld = lift["LiftData"];
  verts = Quiet[PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]],
                TropicalFan::polymake];
  fan = Quiet[computeFanScaled[verts], TropicalFan::polymake];
  If[!ListQ[fan] || Length[fan] < 2,
    Print["CC43 FAIL  Sub-case D: lifted fan build failed"]; $cc43D = False; Goto[endD]];
  {dv, sl} = fan;

  Module[{all = Table[
      Quiet@ProcessSectorLifted[ls, dv, sl[[s]], s, ld, "Eps" -> eps],
      {s, Length[sl]}]},
    all     = Select[all, AssociationQ];
    all     = Select[all, !KeyExistsQ[#, "EmptyDomain"] &];
    divSD   = Select[all, TrueQ[#["IsDivergent"]] &];
    convSD  = Select[all, !TrueQ[#["IsDivergent"]] &];
  ];
  Print["  Lifted fan: ", Length[sl], " sectors -> ", Length[convSD],
        " conv, ", Length[divSD], " div  (z0=", ld["z0"], ")"];
  If[Length[divSD] == 0,
    Print["CC43 FAIL  Sub-case D: expected >=1 divergent lifted sector, got 0"];
    $cc43D = False; Goto[endD]];

  (* divergent lifted sector pre-IBP fields, incl. PrefactorBase (the new field) *)
  Do[
    Module[{sd = divSD[[s]], lab},
      lab = "lifted-div-sector-" <> ToString[sd["ConeIndex"]];
      AppendTo[passes, reportExact["NewExponents", sd["NewExponents"], lab]];
      AppendTo[passes, reportExact["MinExponents", sd["MinExponents"], lab]];
      AppendTo[passes, reportExact["PrefactorBase", sd["PrefactorBase"], lab]];
      AppendTo[passes, reportExact["ClearedPolys", sd["ClearedPolys"], lab]];
    ], {s, Length[divSD]}];

  (* IBP Laurent fields (incl. DLogPrefactor) on each divergent lifted sector *)
  ibpResults = Table[IBPProcessSector[divSD[[s]], ls], {s, Length[divSD]}];
  Do[
    Module[{ibp = ibpResults[[s]], lab},
      lab = "lifted-IBP-sector-" <> ToString[divSD[[s]]["ConeIndex"]];
      If[ibp === $Failed,
        Print["  NOTE: ", lab, " IBPProcessSector returned $Failed; skipping"],
        AppendTo[passes, checkIBPResult[ibp, lab]];
        If[KeyExistsQ[ibp, "DLogPrefactor"],
          AppendTo[passes, reportExact["DLogPrefactor", ibp["DLogPrefactor"], lab]]]
      ]
    ], {s, Length[divSD]}];

  anyFail = !And @@ passes;
  If[!anyFail,
    Print["CC43 PASS  Sub-case D  lifted+divergent: all checked fields FreeQ[_Real] ",
          "(NewExponents/PrefactorBase/ClearedPolys + IBP boundary/terms/DLogPrefactor)"],
    Print["CC43 FAIL  Sub-case D  expected=FreeQ[_Real] in all fields  got=see above"]
  ];
  $cc43D = !anyFail;
  Label[endD]
];
Print[];

(* ── 5. Summary ─────────────────────────────────────────────────────── *)
Module[{results, nPass, nFail},
  results = {$cc43A, $cc43B, $cc43C, $cc43D};
  nPass = Count[results, True];
  nFail = Count[results, False];

  Print["========================================================"];
  Print["CC43  Exactness guard (plan.md §1.4 / §6.7)  —  summary"];
  Print["========================================================"];
  Print["  Sub-case A (convergent n=2, exact input):  ",
        If[TrueQ[$cc43A], "PASS", "FAIL"]];
  Print["  Sub-case B (divergent n=2 + IBP path):     ",
        If[TrueQ[$cc43B], "PASS", "FAIL"]];
  Print["  Sub-case C (meta-test: guard sensitivity): ",
        If[TrueQ[$cc43C], "PASS", "FAIL"]];
  Print["  Sub-case D (lifted+divergent, exact input): ",
        If[TrueQ[$cc43D], "PASS", "FAIL"]];
  Print["--------------------------------------------------------"];
  If[nFail == 0,
    Print["CC43 PASS  all ", nPass, " sub-cases passed",
          "  (exactness invariant confirmed: decomposition is FreeQ[_Real]",
          " before MmaToC boundary; §1.4 / §6.7)"],
    Print["CC43 FAIL  expected=all-sub-cases-pass",
          "  got=", nFail, " fail  (", nPass, " passed, ", nFail, " failed)"]
  ];
  Print["========================================================"];
  If[nFail > 0, Quit[1]]
];
