(* ============================================================================
   TEST/cc_55.wl  —  CHARACTERIZATION test (fable_plan step 2 / step 6)

   Pins the CURRENT behavior of the two public symbols that had ZERO test
   coverage before this file:

     * IdentifyDivergences[sectorData, eps]   (Module 2)
     * IBPCheckBoundary[sectorData, spec, n]  (Module 2b)

   These are "golden master" / characterization tests: they assert the behavior
   that exists TODAY (before any refactor), so a later modularity refactor can be
   proven behavior-preserving.  They are NOT a claim that the behavior is
   "correct" in an external sense — they lock down the status quo.

   Both functions consume a sectorData association.  For a REAL-exponent,
   UNLIFTED sector the internal helpers collapse to trivial forms
   (ibpImagPole -> Im[a0]=0, ibpDivClass -> "Pole" on a Re(a0)<=0 slot), so a
   minimal hand-built sectorData is sufficient and deterministic — no fan /
   Polymake / g++ needed.  Tier 0 (pure WL); runs in isolation under run_all.wl.
   ============================================================================ *)

$pkgRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{$pkgRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$pkgRoot, "tropical_eval.wl"}]];
Print["CC55: packages loaded."];
Print[];

$cc55Pass = True; $cc55Fail = {};
cc55Assert[label_String, test_, exp_, got_] :=
  If[TrueQ[test],
    Print["CC55 PASS  [", label, "]  expected=", exp, "  got=", got],
    Print["CC55 FAIL  [", label, "]  expected=", exp, "  got=", got];
    $cc55Pass = False; AppendTo[$cc55Fail, label]];

(* Quiet the KNOWN, EXPECTED messages these functions raise on the $Failed
   branches (so General::stop noise stays out of the log). *)
SetAttributes[qr, HoldFirst];
qr[e_] := Quiet[e, {TropicalEval::nested, TropicalEval::badck, General::stop}];

(* Suppress the diagnostic Print[] side effects of the functions under test by
   pointing $Output at no channels for the duration of the call, then restore.
   Messages (routed via $Messages) are unaffected — handled by qr above. *)
SetAttributes[silent, HoldFirst];
silent[e_] := Module[{r}, Block[{$Output = {}}, r = e]; r];

eps55 = Symbol["epsCC55"];

(* ============================================================================
   PART A — IdentifyDivergences[sectorData, eps]

   Reads only "NewExponents", "Dimension", "ConeIndex".
     a0 = NewExponents /. eps->0 ;  a1 = D[NewExponents, eps] /. eps->0
   Divergent slot i  <=>  Re(a0[[i]]) <= 0.
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (A) IdentifyDivergences characterization"];

(* A1. Fully convergent sector -> <|"IsDivergent" -> False|>, nothing else. *)
Module[{sd, r},
  sd = <|"NewExponents" -> {1, 2}, "Dimension" -> 2, "ConeIndex" -> 1|>;
  r  = silent@qr@IdentifyDivergences[sd, eps55];
  cc55Assert["A1 convergent -> IsDivergent=False (sole key)",
    r === <|"IsDivergent" -> False|>,
    <|"IsDivergent" -> False|>, r];
];

(* A2. Single logarithmic pole: a0 = {0, 2}, a1 = {1, 0}.
       Divergent variable k=1, ck = a1[[1]] = 1, ak0 = 0. *)
Module[{sd, r},
  sd = <|"NewExponents" -> {eps55, 2}, "Dimension" -> 2, "ConeIndex" -> 7|>;
  r  = silent@qr@IdentifyDivergences[sd, eps55];
  cc55Assert["A2 single pole -> IsDivergent=True",
    AssociationQ[r] && TrueQ[r["IsDivergent"]], True, Head[r]];
  cc55Assert["A2 DivergentVariable = 1", TrueQ[r["DivergentVariable"] === 1],
    1, r["DivergentVariable"]];
  cc55Assert["A2 ck = 1 (exact, no float)",
    r["ck"] === 1 && FreeQ[r["ck"], _Real], 1, r["ck"]];
  cc55Assert["A2 ak0 = 0", r["ak0"] === 0, 0, r["ak0"]];
  cc55Assert["A2 a0 = {0,2}, a1 = {1,0}",
    r["a0"] === {0, 2} && r["a1"] === {1, 0},
    {{0, 2}, {1, 0}}, {r["a0"], r["a1"]}];
];

(* A3. Two divergent variables -> Message[TropicalEval::nested] + $Failed. *)
Module[{sd, r},
  sd = <|"NewExponents" -> {eps55, eps55}, "Dimension" -> 2, "ConeIndex" -> 3|>;
  r  = silent@qr@IdentifyDivergences[sd, eps55];
  cc55Assert["A3 two divergent vars (nested) -> $Failed",
    r === $Failed, $Failed, r];
];

(* A4. Divergent slot with ck = 0 (a0={0,2}, a1={0,0}) -> badck + $Failed. *)
Module[{sd, r},
  sd = <|"NewExponents" -> {0, 2}, "Dimension" -> 2, "ConeIndex" -> 4|>;
  r  = silent@qr@IdentifyDivergences[sd, eps55];
  cc55Assert["A4 ck=0 (higher-order pole) -> $Failed",
    r === $Failed, $Failed, r];
];

(* ============================================================================
   PART B — IBPCheckBoundary[sectorData, integrandSpec, nPoints]

   Diagnostic: for each genuine 1/eps POLE slot (Re(a0)<=0 AND theta=0) it
   samples nPoints random points and returns
     { <|"Variable"->k, "BoundaryMags"-> {nPoints reals}|>, ... }.
   Off-axis / convergent slots contribute nothing.  Stochastic (RandomReal):
   we pin STRUCTURE, not sampled values, and seed for reproducibility.
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (B) IBPCheckBoundary characterization"];

specB = <|"RegulatorSymbol" -> eps55|>;

(* B1. Convergent sector (no pole slot) -> {} . *)
Module[{sd, r},
  sd = <|"Dimension" -> 2, "NewExponents" -> {1, 2},
         "ClearedPolys" -> {{{1, {0, 0}}, {1, {1, 0}}}},
         "PolynomialExponents" -> {-1}, "DetM" -> 1, "ConeIndex" -> 1|>;
  r  = silent@qr@IBPCheckBoundary[sd, specB, 5];
  cc55Assert["B1 convergent sector -> {} (no boundary checks)",
    r === {}, {}, r];
];

(* B2. Single real 1/eps pole slot (a0={0,2}, theta=0 -> "Pole").
       Structure: one entry, keyed Variable/BoundaryMags, Variable=1,
       BoundaryMags has nPoints non-negative real magnitudes. *)
Module[{sd, r, e1},
  SeedRandom[55];
  sd = <|"Dimension" -> 2, "NewExponents" -> {eps55, 2},
         "ClearedPolys" -> {{{1, {0, 0}}, {1, {1, 0}}}},
         "PolynomialExponents" -> {-1}, "DetM" -> 1, "ConeIndex" -> 2|>;
  r  = silent@qr@IBPCheckBoundary[sd, specB, 5];
  cc55Assert["B2 one pole slot -> length-1 result list",
    ListQ[r] && Length[r] === 1, "list len 1",
    If[ListQ[r], Length[r], Head[r]]];
  e1 = If[ListQ[r] && Length[r] >= 1, r[[1]], <||>];
  cc55Assert["B2 entry keys = {Variable, BoundaryMags}",
    AssociationQ[e1] && Sort[Keys[e1]] === {"BoundaryMags", "Variable"},
    {"BoundaryMags", "Variable"}, If[AssociationQ[e1], Keys[e1], e1]];
  cc55Assert["B2 Variable = 1", e1["Variable"] === 1, 1, e1["Variable"]];
  cc55Assert["B2 BoundaryMags: 5 non-negative reals",
    ListQ[e1["BoundaryMags"]] && Length[e1["BoundaryMags"]] === 5 &&
      AllTrue[e1["BoundaryMags"], (NumericQ[#] && Re[#] >= 0 && Im[#] == 0) &],
    "5 reals >= 0",
    If[ListQ[e1["BoundaryMags"]], Length[e1["BoundaryMags"]], e1["BoundaryMags"]]];
];

(* ============================================================================
   Summary
   ============================================================================ *)
Print[];
Print["================================================================"];
If[$cc55Pass,
  Print["CC55 PASS  IdentifyDivergences + IBPCheckBoundary characterization pinned"],
  Print["CC55 FAIL  failures: ", $cc55Fail]];
Print["================================================================"];
