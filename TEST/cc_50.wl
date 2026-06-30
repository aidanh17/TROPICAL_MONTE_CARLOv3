(* ============================================================================
   TEST/cc_50.wl  --  v3 Cross-Check #50  (planIBPCX.md §5, family F2)
   OFF-AXIS divergent directions:  IBP on a single divergent direction that
   carries a NONZERO imaginary exponent theta_k.

   Physics (planIBPCX.md §1.2):  the divergent endpoint exponent is
       s_k = c_k eps + i theta_k,   theta_k != 0
   so the "1/eps pole" sits at eps = -i theta_k / c_k, OFF the real axis: the
   cone is FINITE (pole = 0, value ~ 1/(i theta_k)).  IBP is the UNIQUE engine
   route that resolves it numerically -- it RAISES the divergent exponent before
   flattening, so the bulk theta-phase coefficient is theta/alpha0_raised = O(theta),
   bounded.  The inline-Subtraction and pinned-eps routes flatten the divergent
   direction by its VANISHING exponent (phase coeff ~ theta/(c_k eps)) and oscillate
   too fast to resolve -- they SHARE a blind spot, so neither is a valid F2 oracle
   (planIBPCX.md §5.1 F2).  The independent oracles here are the exact closed-form
   Gamma-Laurent and NIntegrate of the complex original (plus a reseed and a
   batched-vs-per-kp IBP comparison).

   PART A  off-axis Direct (unlifted), exact-Gamma oracle + NIntegrate + reseed
       I(eps) = Int x1^{2eps-1+I c} x2^0 (1+x1+x2)^{-B} dx
              = Gamma(2eps+I c) Gamma(B-2eps-I c-1) / Gamma(B)   (Dirichlet, a2=1)
       eps=0:  pole = 0 (Gamma(I c) regular),  finite = Gamma(I c)Gamma(B-1-I c)/Gamma(B).
     The imaginary weight is ON the divergent variable x1 -> theta_{x1} = c != 0.
   PART B  batching gate (Tier 2, CUBA): batched IBP == per-kp IBP, on the
     off-axis fixture and on a real-divergent multi-kp fixture (planIBPCX.md §4/§4.5).

   Independent oracles per numeric PASS (invariant #4): exact Gamma-Laurent AND
   NIntegrate of the complex original AND a reseed (Part A); per-kp IBP AND exact
   (Part B).  Subtraction is deliberately NOT used as an oracle (it cannot resolve
   off-axis -- stated plainly, planIBPCX.md §5.1).

   Tier 1 (Part A): WL + g++ (MC).  Tier 2 (Part B): + CUBA (skips if absent).
   Run:  wolframscript -file TEST/cc_50.wl
   ============================================================================ *)

$pkgRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{$pkgRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$pkgRoot, "tropical_eval.wl"}]];
Print["CC50: packages loaded."];

SetAttributes[qr, HoldFirst];
qr[e_] := Quiet[e, {TropicalEval::validate, TropicalEval::divergent, TropicalEval::vegasbudget,
   TropicalEval::splitdivmono, General::stop, NIntegrate::slwcon, NIntegrate::ncvb,
   NIntegrate::inumr, NIntegrate::eincr, NIntegrate::izero, NIntegrate::precw,
   TropicalFan::polymake}];

wdir[tag_] := Module[{d = FileNameJoin[{$pkgRoot, "TEST", "INTERFILES", "cc50", tag}]},
  Quiet[CreateDirectory[d, CreateIntermediateDirectories -> True], {CreateDirectory::eexist}]; d];

sci[x_] := Module[{a = N[Abs[x]], e}, If[a == 0., "0",
   e = Floor[Log10[a]]; ToString[NumberForm[a/10^e, {3, 2}]] <> "e" <> ToString[e]]];
fmtC[z_] := ToString[NumberForm[N[Re[z]], {6, 5}]] <> If[Im[z] >= 0, "+", "-"] <>
   ToString[NumberForm[N[Abs[Im[z]]], {6, 5}]] <> "I";

$pass = True; $fail = {};
assert[label_, ok_, detail_] := (
  If[TrueQ[ok], Print["CC50 PASS ", label, "  (", detail, ")"],
     Print["CC50 FAIL ", label, "  (", detail, ")"]; AppendTo[$fail, label]];
  $pass = $pass && TrueQ[ok];
  TrueQ[ok]);

NS = 4000000;

(* ============================================================================
   PART A -- off-axis Direct, exact Gamma + NIntegrate + reseed
   ============================================================================ *)
Print[""];
Print["----------------------------------------------------------------"];
Print["  (A) OFF-AXIS Direct  A={2eps-1+(1/3)I, 0}  poly=1+x1+x2  B=-3"];
Module[{cc = 1/3, Bv = 3, eps, vars, poly, spec, fan, ex, poleRef, finRef,
        rA, rA2, pA, fA, pA2, fA2, imagPoles, ip, niC, gv, fd, finNI, esStar = 0.01, dim = 2},
  eps = Symbol["epsA50"]; vars = {x[1], x[2]}; poly = 1 + x[1] + x[2];
  ex[es_] := Gamma[2 es + I cc] Gamma[Bv - 2 es - I cc - 1]/Gamma[Bv];
  poleRef = N[SeriesCoefficient[Series[ex[ee], {ee, 0, 0}], -1], 12];   (* = 0 *)
  finRef  = N[ex[0], 12];
  Print["   exact pole   = ", fmtC[poleRef], "  (off-axis: no 1/eps pole)"];
  Print["   exact finite = ", fmtC[finRef]];

  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {2 eps - 1 + I cc, 0},
    "PolynomialExponents" -> {-Bv}, "Variables" -> vars,
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> eps|>;
  fan = qr@ComputeDecomposition[qr@PolytopeVertices[poly^(-1), vars]];
  assert["A0 fan built", ListQ[fan] && Length[fan] == 2 && Length[fan[[2]]] > 0,
     "sectors=" <> ToString[If[ListQ[fan], Length[fan[[2]]], fan]]];

  (* Oracle 2: NIntegrate of the complex original; finite = I(eps->0) via eps-fit. *)
  niC[es_] := Module[{us, subr, jac, igc},
    us = Table[Unique["u"], {dim}]; subr = Table[x[i] -> us[[i]]/(1 - us[[i]]), {i, dim}];
    jac = Times @@ Table[1/(1 - us[[i]])^2, {i, dim}];
    igc = (x[1]^(2 es - 1 + I cc) (poly)^(-Bv) /. subr) jac;
    qr@NIntegrate[igc, Evaluate[Sequence @@ ({#, 0, 1} & /@ us)], Method -> "GlobalAdaptive",
      PrecisionGoal -> 6, WorkingPrecision -> 30, MaxRecursion -> 60]];
  gv = niC /@ N[{esStar, 2 esStar, 4 esStar}, 30];
  fd = MapThread[{#1, #2} &, {N[{esStar, 2 esStar, 4 esStar}, 30], gv}];
  finNI = (Fit[{#1, Re[#2]} & @@@ fd, {1, ee}, ee] /. ee -> 0) +
          I (Fit[{#1, Im[#2]} & @@@ fd, {1, ee}, ee] /. ee -> 0);
  Print["   NIntegrate finite (eps-fit) = ", fmtC[finNI]];

  rA = qr@EvaluateTropicalMCIBP[spec, fan, {{}}, "Integrator" -> "MC", "NSamples" -> NS,
    "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> wdir["A"]];
  rA2 = qr@EvaluateTropicalMCIBP[spec, fan, {{}}, "Integrator" -> "MC", "NSamples" -> NS,
    "SeedBase" -> 909, "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> wdir["A2"]];
  If[!(AssociationQ[rA] && KeyExistsQ[rA, "Results"]),
    assert["A IBP returns a result", False, ToString[Head[rA]]],
    (
    pA = rA["Results"][[1]]["PoleCoefficient"]; fA = rA["Results"][[1]]["FinitePart"];
    imagPoles = (#["ImagPole"] & /@ rA["IBPProcessedSectors"]);
    ip = AnyTrue[imagPoles, (NumericQ[#] && Abs[#] > 0.01) || (!NumericQ[#] && !TrueQ[PossibleZeroQ[#]]) &];
    Print["   IBP pole=", fmtC[pA], "  finite=", fmtC[fA], "   ImagPoles=", imagPoles];

    assert["A1 fixture realizes off-axis (theta_k != 0)", ip,
       "ImagPoles=" <> ToString[imagPoles, InputForm]];
    assert["A2 IBP pole ~ 0 (cone finite, no 1/eps)", NumericQ[Abs[pA]] && Abs[pA] < 5*^-3,
       "|pole|=" <> sci[Abs[pA]]];
    assert["A3 IBP finite vs EXACT Gamma-Laurent",
       NumericQ[Abs[fA - finRef]] && Abs[(fA - finRef)/finRef] < 2*^-2,
       "rel=" <> sci[Abs[(fA - finRef)/finRef]]];
    assert["A4 IBP finite vs NIntegrate (independent oracle)",
       NumericQ[Abs[fA - finNI]] && Abs[(fA - finNI)/finNI] < 3*^-2,
       "rel=" <> sci[Abs[(fA - finNI)/finNI]]];
    If[AssociationQ[rA2] && KeyExistsQ[rA2, "Results"],
      pA2 = rA2["Results"][[1]]["PoleCoefficient"]; fA2 = rA2["Results"][[1]]["FinitePart"];
      assert["A5 reseed consistency (finite stable under reseed)",
         Abs[fA2 - fA] < 2*^-2, "|dfin|=" <> sci[Abs[fA2 - fA]]]];
    )
  ];
  Print["   NOTE: Subtraction / pinned-eps are NOT valid oracles here -- they",
        " flatten the divergent slot by its vanishing exponent and oscillate",
        " too fast (planIBPCX.md §1.2/§5.1).  IBP is the supported route."];
];

(* ============================================================================
   PART B -- batching gate: batched IBP == per-kp IBP  (Tier 2, CUBA)
   ============================================================================ *)
Print[""];
Print["----------------------------------------------------------------"];
Print["  (B) batched IBP == per-kp IBP  (planIBPCX.md §4)"];
Module[{cuba},
  cuba = TropicalEval`detectCuba[];
  If[!TrueQ[cuba["Found"]],
    Print["   CC50 (B) SKIP (no CUBA) -- batching gate needs Integrator->VEGAS"],
    Print["   CUBA found: ", cuba["IncludeDir"]];

    (* B1: off-axis Direct, single kp -- batched vs per-kp (and pole~0 preserved) *)
    Module[{cc = 1/3, Bv = 3, eps, vars, poly, spec, fan, rPK, rB},
      eps = Symbol["epsB50a"]; vars = {x[1], x[2]}; poly = 1 + x[1] + x[2];
      spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {2 eps - 1 + I cc, 0},
        "PolynomialExponents" -> {-Bv}, "Variables" -> vars,
        "KinematicSymbols" -> {}, "RegulatorSymbol" -> eps|>;
      fan = qr@ComputeDecomposition[qr@PolytopeVertices[poly^(-1), vars]];
      rPK = qr@EvaluateTropicalMCIBP[spec, fan, {{}}, "Integrator" -> "VEGAS", "Batch" -> False,
        "NSamples" -> 300000, "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> wdir["Bpk"]];
      rB = qr@EvaluateTropicalMCIBP[spec, fan, {{}}, "Integrator" -> "VEGAS", "Batch" -> True,
        "NSamples" -> 300000, "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> wdir["Bb"]];
      If[AssociationQ[rPK] && AssociationQ[rB] && KeyExistsQ[rPK, "Results"] && KeyExistsQ[rB, "Results"],
        Module[{ppk = rPK["Results"][[1]]["PoleCoefficient"], fpk = rPK["Results"][[1]]["FinitePart"],
                pb = rB["Results"][[1]]["PoleCoefficient"], fb = rB["Results"][[1]]["FinitePart"]},
          Print["   off-axis: per-kp fin=", fmtC[fpk], "  batch fin=", fmtC[fb]];
          assert["B1 off-axis batched pole ~ per-kp pole (both ~0)", Abs[pb - ppk] < 5*^-3,
             "|dpole|=" <> sci[Abs[pb - ppk]]];
          assert["B2 off-axis batched finite == per-kp finite", Abs[fb - fpk] < 1*^-2,
             "|dfin|=" <> sci[Abs[fb - fpk]]]],
        assert["B1/B2 off-axis VEGAS runs returned results",
           AssociationQ[rPK] && AssociationQ[rB],
           "pk=" <> ToString[Head[rPK]] <> " b=" <> ToString[Head[rB]]]];
    ];

    (* B3: real-divergent, MULTI-KP -- batched grid shared across points (the
       headline batching property), vs per-kp. *)
    Module[{eps, vars, poly, spec, fan, kps, rPK, rB, maxd},
      eps = Symbol["epsB50b"]; vars = {x[1], x[2]}; poly = 1 + p1 x[1] + x[2];
      spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {2 eps - 1, 0},
        "PolynomialExponents" -> {-2}, "Variables" -> vars,
        "KinematicSymbols" -> {p1}, "RegulatorSymbol" -> eps|>;
      fan = qr@ComputeDecomposition[qr@PolytopeVertices[poly^(-1), vars]];
      kps = {{1.0}, {1.5}, {2.0}};
      rPK = qr@EvaluateTropicalMCIBP[spec, fan, kps, "Integrator" -> "VEGAS", "Batch" -> False,
        "NSamples" -> 200000, "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> wdir["B3pk"]];
      rB = qr@EvaluateTropicalMCIBP[spec, fan, kps, "Integrator" -> "VEGAS", "Batch" -> True,
        "NSamples" -> 200000, "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> wdir["B3b"]];
      If[AssociationQ[rPK] && AssociationQ[rB] && KeyExistsQ[rPK, "Results"] && KeyExistsQ[rB, "Results"],
        maxd = Max[Table[
          Abs[rB["Results"][[i]]["PoleCoefficient"] - rPK["Results"][[i]]["PoleCoefficient"]] +
          Abs[rB["Results"][[i]]["FinitePart"] - rPK["Results"][[i]]["FinitePart"]],
          {i, Length[kps]}]];
        Print["   real-div multi-kp: max |batched - per-kp| (pole+fin) = ", sci[maxd]];
        assert["B3 batched IBP == per-kp IBP across kp (real divergent)", maxd < 1*^-2,
           "max|d|=" <> sci[maxd]],
        assert["B3 real-div VEGAS runs returned results",
           AssociationQ[rPK] && AssociationQ[rB],
           "pk=" <> ToString[Head[rPK]] <> " b=" <> ToString[Head[rB]]]];
    ];
  ];
];

Print[""];
Print["================================================================"];
If[$pass,
  Print["CC50 PASS  off-axis (theta!=0) divergent directions resolved by IBP ",
        "(pole~0, finite) vs exact Gamma-Laurent & NIntegrate; batched IBP == ",
        "per-kp IBP.  failures={}"],
  Print["CC50 FAIL  failed: ", $fail]];
Print["================================================================"];
If[!$pass, Quit[1]];
