(* ============================================================================
   TEST/cc_44.wl  —  Cross-check #44 (planAXpDIV.md §3, §7, §8)
   Case-B (domain x divergent-variable coupling) — refusal + safety invariant.

   Background (planAXpDIV.md §3).  In a LIFTED + DIVERGENT sector the 1/eps pole
   lives in a surviving effective exponent atilde_k.  The lifted DOMAIN indicator
   is a half-space Boole[log_ypstar <= 0] with IndicatorCoeffs ic_j = mOther_j /
   atilde_j.  Two cases:

     Case A (supported): the divergent slot decouples from the domain face,
       i.e. DomainConstraint === None OR ic_k = 0 (mOther_k = 0).  The indicator
       commutes with the y_k IBP/subtraction and rides along on the other coords;
       the 1/eps extraction is identical to the unlifted pole, just multiplied by
       a y_k-independent indicator.  ProcessSectorLifted emits these.

     Case B (refused): ic_k != 0 for EVERY admissible pivot — the divergent
       endpoint and the domain face couple.  A log-space remap is required
       (future work).  ProcessSectorLifted issues TropicalEval::liftdivdomain
       and returns $Failed (a clean refusal, never a wrong number).

   The candidate ranking PREFERS a decoupled (Case A) pivot whenever one exists
   (the caseA-first key, planAXpDIV.md §4.3), so Case B only fires when NO
   admissible pivot decouples.  Empirically this is rare: across a diverse
   battery the pivot preference resolves every would-be-Case-B sector to Case A.

   This cross-check verifies the v3 honesty contract for the coupled case:
     INVARIANT (the safety property): every divergent lifted SectorData that
       ProcessSectorLifted EMITS is Case A — DomainConstraint === None, or its
       IndicatorCoeffs[[DivergentVariable]] is exactly 0.  Hence the engine NEVER
       silently emits a coupled (Case B) divergent sector that would yield a
       wrong pole; it either emits a decoupled sector or refuses.
     REFUSAL: TropicalEval::liftdivdomain is defined and is the wired refusal
       (a clean $Failed), so any genuine Case-B sector aborts the call rather
       than producing a finite-looking but wrong Laurent — consistent with the
       v3 rule pinned by #27 / #34 (a known limitation => clean $Failed).

   Tier 1: WL only (no g++, no MC).
   Run:  wolframscript -file TEST/cc_44.wl
   ============================================================================ *)

$pkgRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{$pkgRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$pkgRoot, "tropical_eval.wl"}]];
Print["CC44: packages loaded."];
Print[];

$cc44AllPass = True; $cc44Fail = {};
cc44Assert[label_String, test_, exp_String, got_String] :=
  If[TrueQ[test],
    Print["CC44 PASS  [", label, "]  expected=", exp, "  got=", got],
    Print["CC44 FAIL  [", label, "]  expected=", exp, "  got=", got];
    $cc44AllPass = False; AppendTo[$cc44Fail, label]];

SetAttributes[qr, HoldFirst];
qr[e_] := Quiet[e, {General::stop, TropicalFan::polymake, NIntegrate::slwcon,
   TropicalEval::liftdivdomain}];

(* caseAQ: True iff this emitted divergent lifted sector is decoupled. *)
caseAQ[sd_Association] := Module[{dc = sd["DomainConstraint"], k = sd["DivergentVariable"]},
  dc === None || (k >= 1 && k <= Length[dc["IndicatorCoeffs"]] &&
                  PossibleZeroQ[dc["IndicatorCoeffs"][[k]]])];

(* ── 1.  Refusal guard is defined and wired ─────────────────────────── *)
cc44Assert["liftdivdomain message defined",
  StringQ[TropicalEval::liftdivdomain] &&
    StringContainsQ[TropicalEval::liftdivdomain, "Case B"],
  "message string mentioning Case B",
  If[StringQ[TropicalEval::liftdivdomain], "defined", "undefined"]];

(* ── 2.  Safety invariant over a battery of lifted+divergent specs ────── *)
eps = Symbol["e44"];
battery = {
  {"2D x1x2 k1", {1 + 10^6 x[1] x[2] + x[1]^2 + x[2]^2}, {-1 + eps, 0}, {-2},
    <|"PolyIndex" -> 1, "ExponentVector" -> {1, 1}, "k" -> 1|>},
  {"2D x1^2 k3", {1 + x[1] x[2] + 10^6 x[1]^2 + x[2]^2}, {-1 + eps, 0}, {-2},
    <|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 3|>},
  {"2D x2 k1", {1 + 10^6 x[2] + x[1]^2 + x[2]^2}, {-1 + eps, 0}, {-2},
    <|"PolyIndex" -> 1, "ExponentVector" -> {0, 1}, "k" -> 1|>},
  {"2D x1 k1", {1 + 10^6 x[1] + x[1]^2 + x[2]^2}, {-1 + eps, 0}, {-2},
    <|"PolyIndex" -> 1, "ExponentVector" -> {1, 0}, "k" -> 1|>},
  {"2D x1x2 k2 1e9", {1 + 10^9 x[1] x[2] + x[1]^3 + x[2]}, {-1 + eps, 0}, {-2},
    <|"PolyIndex" -> 1, "ExponentVector" -> {1, 1}, "k" -> 2|>},
  {"3D x1x2 k1", {1 + 10^6 x[1] x[2] + x[2] x[3] + x[1]^2 + x[3]^2},
    {-1 + eps, 0, 0}, {-2}, <|"PolyIndex" -> 1, "ExponentVector" -> {1, 1, 0}, "k" -> 1|>},
  {"3D x1x3 k1", {1 + 10^6 x[1] x[3] + x[2] + x[1]^2 + x[3]^2},
    {-1 + eps, 0, 0}, {-2}, <|"PolyIndex" -> 1, "ExponentVector" -> {1, 0, 1}, "k" -> 1|>}
};

totDiv = 0; totRefuse = 0; allInvariant = True;
Do[
  Module[{label, polys, mexp, pexp, lr, vars, spec, lift, ls, ld, verts, fan, dv, sl,
          refused = False},
    {label, polys, mexp, pexp, lr} = cand;
    vars = Table[x[i], {i, Length[mexp]}];
    spec = <|"Polynomials" -> polys, "MonomialExponents" -> mexp,
      "PolynomialExponents" -> pexp, "Variables" -> vars,
      "KinematicSymbols" -> {}, "RegulatorSymbol" -> eps|>;
    lift = qr@LiftCoefficients[spec, {lr}];
    If[!AssociationQ[lift], Print["  ", label, ": lift failed, skipping"]; Continue[]];
    ls = lift["LiftedSpec"]; ld = lift["LiftData"];
    verts = qr@PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]];
    fan = qr@computeFanScaled[verts];
    If[!ListQ[fan] || Length[fan] < 2, Print["  ", label, ": fan failed, skipping"]; Continue[]];
    {dv, sl} = fan;
    Do[
      Check[
        Module[{sd = Quiet@ProcessSectorLifted[ls, dv, sl[[s]], s, ld, "Eps" -> eps]},
          Which[
            sd === $Failed, Null,  (* some other refusal (liftcomplex/liftnopivot) *)
            AssociationQ[sd] && TrueQ[sd["IsDivergent"]],
              totDiv++;
              If[!caseAQ[sd],
                allInvariant = False;
                Print["  INVARIANT VIOLATED: ", label, " sector ", s,
                      " emitted a COUPLED (Case B) divergent sector!"]]]],
        (* liftdivdomain fired -> clean $Failed refusal *)
        refused = True; totRefuse++,
        {TropicalEval::liftdivdomain}],
      {s, Length[sl]}];
    Print["  ", label, ": done", If[refused, " (liftdivdomain refusal seen)", ""]];
  ], {cand, battery}];

Print[];
Print["CC44: across battery — ", totDiv, " divergent lifted sectors emitted, ",
      totRefuse, " liftdivdomain refusals."];

cc44Assert["safety invariant: every emitted divergent lifted sector is Case A",
  allInvariant && totDiv > 0,
  "all emitted divergent sectors decoupled (ic_k=0 or domain None), >=1 emitted",
  "emitted=" <> ToString[totDiv] <> " invariantHeld=" <> ToString[allInvariant]];

(* ── 3.  The refusal is a CLEAN $Failed (no crash / no wrong number) ──── *)
(* Drive ProcessSectorLifted through a battery sector and confirm that whenever
   it would refuse Case B it returns $Failed (verified above via Check); here we
   confirm a refused call surfaces $Failed end-to-end through the driver for a
   spec, OR — when the pivot preference avoids Case B everywhere (the common
   outcome) — the call succeeds with only Case-A divergent sectors. *)
Module[{spec, res},
  spec = <|"Polynomials" -> {1 + x[1] x[2] + 10^6 x[1]^2 + x[2]^2},
    "MonomialExponents" -> {-1 + eps, 0}, "PolynomialExponents" -> {-2},
    "Variables" -> {x[1], x[2]}, "KinematicSymbols" -> {}, "RegulatorSymbol" -> eps|>;
  res = qr@EvaluateTropicalMCLifted[spec, {{}},
    "LiftRules" -> {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 3|>},
    "Method" -> "IBP", "Integrator" -> "MC", "NSamples" -> 20000,
    "RunChecks" -> False, "Verbose" -> False,
    "WorkingDirectory" -> FileNameJoin[{$pkgRoot, "TEST", "INTERFILES", "cc44"}]];
  cc44Assert["driver returns a clean result or clean $Failed (never a crash)",
    res === $Failed || (AssociationQ[res] && KeyExistsQ[res, "Results"]),
    "$Failed or Association[Results]",
    If[res === $Failed, "$Failed", If[AssociationQ[res], "Association", ToString[Head[res]]]]];
];

Print[];
Print["================================================================"];
If[$cc44AllPass,
  Print["CC44 PASS  Case-B refusal wired + safety invariant holds  failures={}"],
  Print["CC44 FAIL  failed: ", $cc44Fail]];
Print["================================================================"];
If[!$cc44AllPass, Quit[1]];
