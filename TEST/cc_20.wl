(* ============================================================================
   TEST/cc_20.wl  --  Cross-check #20: Anchor (integer k) sweep vs the k* rule.

   Plan.md §8.2 #20 (F3a uplift / F3c anchor selection):
       "Anchor (integer k) sweep vs the k* rule  --  argmin relErr >= k* in
        sensitive cases (H1/H2/H3)"

   Strategy (Tier 1: WL + g++ only, no CUBA):
     For each of three sensitive cases (H1: large coeff, H2: degenerate-fan
     regime, H3: small coeff) compute
       * k* = max(1, Ceil[|log10|C|| / 3])   (tau = 1000, anchor_selection_procedure.md §3)
       * for each k in {1 .. k*+2}: lift, compute exact NIntegrate reference,
         run ProcessSectorLifted geometry scan, record relErr (MC via v3
         EvaluateTropicalMCLifted when the sector machinery is working, else
         geometry-only + NIntegrate for the reference)
       * PASS: the k that minimises relErr (or achieves the best geometry gate
         score) satisfies  k_best >= k*  for H1 and H3 (where mechanism (b)
         dominates), and the k* geometry gate passes for H1/H3 while for H2
         (Case C / degenerate) the gate correctly FAILS for all k, causing
         "do not lift" -- the gate fires cleanly.

   Ported from the logic of:
     OLD_CODE/TROPICAL_MONTE_CARLO2/TROPICAL_MONTE_CARLOv2/SANDBOX/z0_sweep_common.wl
     OLD_CODE/TROPICAL_MONTE_CARLO2/TROPICAL_MONTE_CARLOv2/SANDBOX/z0_sweep_battery.wl
     OLD_CODE/TROPICAL_MONTE_CARLO2/TROPICAL_MONTE_CARLOv2/AUXT/anchor_selection_procedure.md

   The old sweep was a large CSV study; here we distil it to the three
   decisive assertions the plan requires, each printing PASS/FAIL.

   Run:  wolframscript -file TEST/cc_20.wl
   (from the TROPICAL_MONTE_CARLO3 root, or anywhere -- uses absolute paths)
   ============================================================================ *)

(* ── 0. Locate package root and load ──────────────────────────────────────── *)
$pkgRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
$evalWL  = FileNameJoin[{$pkgRoot, "tropical_eval.wl"}];
$fanWL   = FileNameJoin[{$pkgRoot, "tropical_fan.wl"}];
$ioDir   = FileNameJoin[{$pkgRoot, "TEST", "INTERFILES"}];

If[!FileExistsQ[$evalWL],
  Print["CC20 FAIL  tropical_eval.wl not found at ", $evalWL]; Quit[1]];
If[!FileExistsQ[$fanWL],
  Print["CC20 FAIL  tropical_fan.wl not found at ", $fanWL]; Quit[1]];

Get[$fanWL];
Get[$evalWL];

Quiet[CreateDirectory[$ioDir], {CreateDirectory::eexist}];

Print["CC20: packages loaded from ", $pkgRoot];
Print["CC20: INTERFILES directory: ", $ioDir];
Print[];

(* ── 1. k* formula (anchor_selection_procedure.md §3, tau = 1000) ─────────── *)
(* k* = max(1, Ceil[ |log10|C|| / log10(tau) ])  with tau=1000 => log10(tau)=3 *)
kStarFor[magnitude_] :=
  Max[1, Ceiling[Abs[Log[10, N[Abs[magnitude]]]] / 3]];

(* ── 2. Geometry gate scan ──────────────────────────────────────────────────
   For a given lifted spec + fan, scan all sectors with ProcessSectorLifted
   and return: {gateOK, nSectors, nSurvivors, hasConstAny, minMargin}
   gateOK = exists survivor with HasConstantTerm=True and all Re[atilde]>0.
   --------------------------------------------------------------------------- *)
geometryGate[liftedSpec_Association, fan_List, liftData_Association] :=
Module[{dv, sl, sectors, survivors, hasConstList, margins, gateOK},
  {dv, sl} = fan;
  sectors = Table[
    Quiet @ ProcessSectorLifted[liftedSpec, dv, sl[[s]], s, liftData],
    {s, Length[sl]}
  ];
  survivors = Select[sectors,
    AssociationQ[#] && !TrueQ[Lookup[#, "EmptyDomain", False]] &];
  hasConstList = Lookup[#, "HasConstantTerm", False] & /@ survivors;
  margins = Table[
    Min[Re[N[Lookup[s, "NewExponents", {Infinity}]]]], {s, survivors}];
  gateOK = MemberQ[hasConstList, True] &&
            AnyTrue[survivors, Min[Re[N[Lookup[#, "NewExponents", {1}]]]] > 0 &];
  <|"gateOK" -> gateOK,
    "nSectors"   -> Length[sl],
    "nSurvivors" -> Length[survivors],
    "hasConstAny"-> MemberQ[hasConstList, True],
    "minMargin"  -> If[margins === {}, Indeterminate, Min[margins]]|>
];

(* ── 3. Lifted NIntegrate reference ─────────────────────────────────────────
   Direct NIntegrate of the ORIGINAL (unlifted) integrand.  This is the truth
   value independent of any k choice.
   --------------------------------------------------------------------------- *)
nIntegrateRef[spec_Association] :=
Module[{vars, integrand, res},
  vars = spec["Variables"];
  integrand = (Times @@ MapThread[Power,
    {spec["Polynomials"], spec["PolynomialExponents"]}]) *
    (Times @@ MapThread[Power, {vars, spec["MonomialExponents"]}]);
  res = Quiet @ NIntegrate[integrand,
    Evaluate[Sequence @@ ({#, 0, Infinity} & /@ vars)],
    MaxRecursion -> 30, PrecisionGoal -> 5, Method -> "GlobalAdaptive"];
  If[NumericQ[res], res, $Failed]
];

(* ── 4. Per-k sweep for one case ────────────────────────────────────────────
   For each k:
     (a) lift the spec
     (b) build the fan (auto, via PolytopeVertices + computeFanScaled)
     (c) run geometry gate
     (d) run MC via EvaluateTropicalMCLifted (using "MC" integrator, no CUBA)
     (e) compute relErr vs truth
   Returns a list of <|k, gateOK, hasConstAny, minMargin, mcVal, relErr|>
   --------------------------------------------------------------------------- *)
sweepOneCase[spec_Association, proxyPoly_, liftRule_, truth_, kList_List,
             nSamples_Integer : 500000] :=
Module[{vars, results = {}},
  vars = spec["Variables"];
  Do[
    Module[{rule, liftRes, liftedSpec, liftData, verts, fan, gate,
            mcRes, mcVal, relErr, z0},
      rule = Append[liftRule, "k" -> k];
      liftRes = Quiet @ LiftCoefficients[spec, {rule}];
      If[!AssociationQ[liftRes],
        AppendTo[results, <|"k" -> k, "gateOK" -> False,
          "hasConstAny" -> False, "minMargin" -> Indeterminate,
          "mcVal" -> Indeterminate, "relErr" -> Indeterminate,
          "skipReason" -> "LiftCoefficients failed"|>];
        Continue[]
      ];
      liftedSpec = liftRes["LiftedSpec"];
      liftData   = liftRes["LiftData"];
      z0         = liftData["z0"];

      (* Build fan from the lifted spec's Newton polytope *)
      verts = Quiet @ PolytopeVertices[
        (Times @@ liftedSpec["Polynomials"])^(-1),
        liftedSpec["Variables"]];
      If[!ListQ[verts],
        AppendTo[results, <|"k" -> k, "gateOK" -> False,
          "hasConstAny" -> False, "minMargin" -> Indeterminate,
          "mcVal" -> Indeterminate, "relErr" -> Indeterminate,
          "skipReason" -> "PolytopeVertices failed", "z0" -> N[z0]|>];
        Continue[]
      ];
      fan = Quiet @ computeFanScaled[verts];
      If[fan === $Failed || !ListQ[fan] || Length[fan] < 2,
        AppendTo[results, <|"k" -> k, "gateOK" -> False,
          "hasConstAny" -> False, "minMargin" -> Indeterminate,
          "mcVal" -> Indeterminate, "relErr" -> Indeterminate,
          "skipReason" -> "fan degenerate/failed", "z0" -> N[z0]|>];
        Continue[]
      ];
      (* Check simplex size consistency -- degeneracy detection *)
      With[{nv = Length[liftedSpec["Variables"]], sl = fan[[2]]},
        If[!AllTrue[sl, Length[#] == nv &],
          AppendTo[results, <|"k" -> k, "gateOK" -> False,
            "hasConstAny" -> False, "minMargin" -> Indeterminate,
            "mcVal" -> Indeterminate, "relErr" -> Indeterminate,
            "skipReason" -> "fan simplex dim mismatch (degenerate polytope)",
            "z0" -> N[z0]|>];
          Continue[]
        ]
      ];

      (* Geometry gate *)
      gate = geometryGate[liftedSpec, fan, liftData];

      (* MC via EvaluateTropicalMCLifted -- Integrator="MC" so no CUBA needed *)
      mcRes = Quiet @ EvaluateTropicalMCLifted[spec, {{}},
        "LiftRules" -> {rule},
        "FanData"   -> fan,
        "Integrator"-> "MC",
        "NSamples"  -> nSamples,
        "RunChecks" -> False,
        "Verbose"   -> False,
        "WorkingDirectory" -> $ioDir];

      If[AssociationQ[mcRes] && KeyExistsQ[mcRes, "Results"] &&
         Length[mcRes["Results"]] >= 1,
        mcVal = mcRes["Results"][[1]]["Re"];
        relErr = If[NumericQ[truth] && truth =!= 0 && NumericQ[mcVal],
          Abs[(mcVal - truth) / truth], Indeterminate],
        mcVal  = Indeterminate;
        relErr = Indeterminate
      ];

      AppendTo[results, <|
        "k"          -> k,
        "z0"         -> N[z0],
        "gateOK"     -> gate["gateOK"],
        "hasConstAny"-> gate["hasConstAny"],
        "minMargin"  -> gate["minMargin"],
        "nSectors"   -> gate["nSectors"],
        "nSurvivors" -> gate["nSurvivors"],
        "mcVal"      -> mcVal,
        "relErr"     -> relErr,
        "skipReason" -> ""|>]
    ],
    {k, kList}
  ];
  results
];

(* ── 5. Reporting helper ─────────────────────────────────────────────────── *)
fmtN[x_] := If[NumericQ[x], ToString[N[x, 4]], ToString[x]];

printSweepRow[row_Association] :=
  Print["    k=", row["k"], "  z0=", If[KeyExistsQ[row,"z0"],fmtN[row["z0"]],"?"],
        "  gateOK=", row["gateOK"],
        "  hasConst=", row["hasConstAny"],
        "  minMargin=", fmtN[Lookup[row,"minMargin",Indeterminate]],
        "  relErr=", fmtN[row["relErr"]],
        If[StringLength[Lookup[row,"skipReason",""]]>0,
           "  [" <> row["skipReason"] <> "]", ""]
  ];

(* ── 6. CASE H1 -- large coefficient (mechanism b dominates; k*=2 for |C|=10^6) ──
   P = 1 + 10^6 x[1]^2 + x[2]^2 + x[1] x[2]^2,  B=-2
   kStar = Ceil[6/3] = 2
   PASS: k_best (argmin relErr among gate-passing k) >= kStar=2.
   --------------------------------------------------------------------------- *)
Print["========================================================"];
Print["CC20  H1 -- large coefficient  |C|=10^6  (expect k*=2)"];
Print["========================================================"];

Module[{spec, mag, kStar, kList, truth, rows, gateRows, bestK, bestRelErr, pass},
  mag   = 10^6;
  kStar = kStarFor[mag];
  kList = Range[1, kStar + 2];   (* 1..4 *)
  Print["  k* = ", kStar, "  (formula: Ceil[|log10(", mag, ")|/3])"];

  spec = <|
    "Polynomials"        -> {1 + mag x[1]^2 + x[2]^2 + x[1] x[2]^2},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents"-> {-2},
    "Variables"          -> {x[1], x[2]},
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None|>;

  Print["  Computing NIntegrate reference ..."];
  truth = nIntegrateRef[spec];
  Print["  truth = ", fmtN[truth]];
  If[truth === $Failed,
    Print["CC20 FAIL  expected=<truth> got=$Failed (NIntegrate failed for H1)"];
    Quit[1]
  ];

  rows = sweepOneCase[spec, 1 + x[1]^2 + x[2]^2 + x[1] x[2]^2,
    <|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}|>, truth, kList, 600000];

  Print["  Sweep results:"];
  printSweepRow /@ rows;

  (* PASS criterion: argmin relErr among gate-passing k is >= kStar *)
  gateRows = Select[rows, TrueQ[#["gateOK"]] && NumericQ[#["relErr"]] &];
  If[gateRows === {},
    Print["CC20 FAIL  expected=gate-pass  got=no gate-passing k for H1"];
    Quit[1]
  ];
  {bestRelErr, bestK} = MinimalBy[gateRows, #["relErr"] &][[1]] /.
    r_Association :> {r["relErr"], r["k"]};

  pass = (bestK >= kStar);
  Print["  bestK=", bestK, "  bestRelErr=", fmtN[bestRelErr],
        "  kStar=", kStar];
  If[pass,
    Print["CC20 PASS  H1: argmin-relErr k=", bestK,
          " >= k*=", kStar, " (relErr=", fmtN[bestRelErr], ")"],
    Print["CC20 FAIL  H1: expected=bestK>=", kStar,
          " got=bestK=", bestK, " (relErr=", fmtN[bestRelErr], ")"]
  ];
  Print[];
  If[!pass, Quit[1]]
];

(* ── 7. CASE H2 -- degenerate geometry (Case C; gate must refuse all k) ─────
   P = 1 + 10^8 x[1]^3 x[2] + x[2]^3,  B=-3
   The lifted Newton polytope has a lineality pair for any k (benchmark Case C).
   PASS: gateOK=False for every k in {2,3,4}  (lifting is correctly refused).
   --------------------------------------------------------------------------- *)
Print["========================================================"];
Print["CC20  H2 -- degenerate (Case C)  |C|=10^8  (gate must refuse all k)"];
Print["========================================================"];

Module[{spec, mag, kStar, kList, rows, anyGateOK, pass},
  mag   = 10^8;
  kStar = kStarFor[mag];
  kList = {2, 3, 4};
  Print["  k* = ", kStar, "  sweep k in ", kList];

  spec = <|
    "Polynomials"        -> {1 + mag x[1]^3 x[2] + x[2]^3},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents"-> {-3},
    "Variables"          -> {x[1], x[2]},
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None|>;

  (* Use explicit fans provided by the old sweep for Case C:
     the lineality pair {k,0,-3}/{-k,0,3} means the auto-fan yields a
     lower-dimensional polytope (degenerateQ fires).  We test via the auto-fan
     path; all k should degenerate (return $Failed fan or gateOK=False). *)
  rows = Table[
    Module[{rule, liftRes, liftedSpec, liftData, verts, fan, gate, z0},
      rule = <|"PolyIndex" -> 1, "ExponentVector" -> {3, 1}, "k" -> k|>;
      liftRes = Quiet @ LiftCoefficients[spec, {rule}];
      If[!AssociationQ[liftRes],
        <|"k" -> k, "gateOK" -> False, "skipReason" -> "LiftCoefficients failed"|>,
        liftedSpec = liftRes["LiftedSpec"];
        liftData   = liftRes["LiftData"];
        z0         = liftData["z0"];
        verts = Quiet @ PolytopeVertices[
          (Times @@ liftedSpec["Polynomials"])^(-1),
          liftedSpec["Variables"]];
        fan = If[ListQ[verts],
          Quiet @ computeFanScaled[verts], $Failed];
        If[fan === $Failed || !ListQ[fan] || Length[fan] < 2,
          <|"k" -> k, "z0" -> N[z0], "gateOK" -> False,
            "skipReason" -> "fan degenerate (expected for Case C)"|>,
          (* check simplex dim *)
          With[{nv = Length[liftedSpec["Variables"]], sl = fan[[2]]},
            If[!AllTrue[sl, Length[#] == nv &],
              <|"k" -> k, "z0" -> N[z0], "gateOK" -> False,
                "skipReason" -> "simplex dim mismatch (degenerate, expected)"|>,
              gate = geometryGate[liftedSpec, fan, liftData];
              Join[gate, <|"k" -> k, "z0" -> N[z0], "skipReason" -> ""|>]
            ]
          ]
        ]
      ]
    ],
    {k, kList}
  ];

  Print["  Sweep results:"];
  Do[
    Print["    k=", r["k"],
          If[KeyExistsQ[r,"z0"], "  z0=" <> fmtN[r["z0"]], ""],
          "  gateOK=", r["gateOK"],
          If[StringLength[Lookup[r,"skipReason",""]]>0,
             "  [" <> r["skipReason"] <> "]", ""]],
    {r, rows}
  ];

  anyGateOK = AnyTrue[rows, TrueQ[Lookup[#, "gateOK", False]] &];
  pass = !anyGateOK;
  If[pass,
    Print["CC20 PASS  H2: all k in ", kList,
          " correctly gate-refused (degenerate polytope; do-not-lift is correct)"],
    Print["CC20 FAIL  H2: expected=gateOK=False for all k  got=gateOK=True for some k"]
  ];
  Print[];
  If[!pass, Quit[1]]
];

(* ── 8. CASE H3 -- small coefficient (|C|=10^-4; k*=2 same as H1 by symmetry) ──
   P = 1 + 10^-4 x[1]^2 + x[2]^2 + x[1] x[2]^2,  B=-2
   kStar = Ceil[|-4|/3] = Ceil[4/3] = 2
   PASS: k_best >= kStar = 2 among gate-passing k.
   --------------------------------------------------------------------------- *)
Print["========================================================"];
Print["CC20  H3 -- small coefficient  |C|=10^-4  (expect k*=2)"];
Print["========================================================"];

Module[{spec, mag, kStar, kList, truth, rows, gateRows, bestK, bestRelErr, pass},
  mag   = 10^-4;
  kStar = kStarFor[mag];
  kList = Range[1, kStar + 2];   (* 1..4 *)
  Print["  k* = ", kStar, "  (formula: Ceil[|log10(", N[mag], ")|/3])"];

  spec = <|
    "Polynomials"        -> {1 + mag x[1]^2 + x[2]^2 + x[1] x[2]^2},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents"-> {-2},
    "Variables"          -> {x[1], x[2]},
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None|>;

  Print["  Computing NIntegrate reference ..."];
  truth = nIntegrateRef[spec];
  Print["  truth = ", fmtN[truth]];
  If[truth === $Failed,
    Print["CC20 FAIL  expected=<truth> got=$Failed (NIntegrate failed for H3)"];
    Quit[1]
  ];

  rows = sweepOneCase[spec, 1 + x[1]^2 + x[2]^2 + x[1] x[2]^2,
    <|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}|>, truth, kList, 600000];

  Print["  Sweep results:"];
  printSweepRow /@ rows;

  gateRows = Select[rows, TrueQ[#["gateOK"]] && NumericQ[#["relErr"]] &];
  If[gateRows === {},
    Print["CC20 FAIL  expected=gate-pass  got=no gate-passing k for H3"];
    Quit[1]
  ];
  {bestRelErr, bestK} = MinimalBy[gateRows, #["relErr"] &][[1]] /.
    r_Association :> {r["relErr"], r["k"]};

  pass = (bestK >= kStar);
  Print["  bestK=", bestK, "  bestRelErr=", fmtN[bestRelErr],
        "  kStar=", kStar];
  If[pass,
    Print["CC20 PASS  H3: argmin-relErr k=", bestK,
          " >= k*=", kStar, " (relErr=", fmtN[bestRelErr], ")"],
    Print["CC20 FAIL  H3: expected=bestK>=", kStar,
          " got=bestK=", bestK, " (relErr=", fmtN[bestRelErr], ")"]
  ];
  Print[];
  If[!pass, Quit[1]]
];

(* ── 9. k* formula unit-tests (pure Tier-1 symbolic; no MC needed) ─────────
   These check the formula kStarFor directly against the manual table in
   anchor_selection_procedure.md §3 (Table: Case A, Case B, Example 20).
   --------------------------------------------------------------------------- *)
Print["========================================================"];
Print["CC20  k* formula unit-tests (anchor_selection_procedure.md §3 table)"];
Print["========================================================"];

Module[{tests, nPass = 0, nFail = 0},
  tests = {
    (*  |C|   expected k*   label *)
    {10^6,  2, "Case A (|C|=10^6)"},
    {10^4,  2, "Case B (|C|=10^4)"},
    {10^-4, 2, "Example 20 (|C|=10^-4)"},
    {10^3,  1, "on-band (|C|=10^3, k*=1)"},
    {10^9,  3, "very large (|C|=10^9)"},
    {10^-7, 3, "very small (|C|=10^-7)"},
    {1,     1, "unit coeff (k*=1)"}
  };
  Do[
    Module[{mag, expected, label, got, ok},
      {mag, expected, label} = t;
      got = kStarFor[mag];
      ok  = (got === expected);
      If[ok,
        Print["  PASS  ", label, "  k*=", got];
        nPass++,
        Print["  FAIL  ", label, "  expected=", expected, " got=", got];
        nFail++]
    ],
    {t, tests}
  ];
  Print[];
  If[nFail == 0,
    Print["CC20 PASS  k* formula: all ", nPass, " unit-tests correct"],
    Print["CC20 FAIL  k* formula: expected=", nPass + nFail, " pass  got=",
          nFail, " fail"]
  ];
  Print[];
  If[nFail > 0, Quit[1]]
];

(* ── 10. Final summary ───────────────────────────────────────────────────── *)
Print["========================================================"];
Print["CC20 PASS  Anchor k-sweep vs k* rule: H1 k_best>=k*, H2 degenerate refused,"];
Print["           H3 k_best>=k*, k* formula correct on all test cases."];
Print["========================================================"];
