(* ============================================================================
   TEST/cc_47.wl  —  Cross-check #47 (planAXpDIVv2.md §6/§7)
   Complex *MONOMIAL* exponents A  x  lifting  (x divergence), via
   ComplexExponentMode -> "SplitRealImag", which now splits BOTH A and B:
     Re(A), Re(B) build the real importance measure / domain indicator;
     Im(A) rides a NEW bare-monomial oscillatory phase
        exp( i (Const_A + Sum_a (Im(A).M)_a / alpha0_a * log y'_a) ),
     the monomial twin of the Im(B) MonoFactorLog phase (planCXLIFTDIV).

   This is the principal-series structure (bubble at imaginary mu): the FIRST
   lifted integrand whose MONOMIAL exponents are complex.  Before this plan the
   lift aborted on every cone with TropicalEval::liftcomplex; SplitRealImag now
   realifies A on the same footing as B and clears it.

   Parts:
     (A) CONVERGENT complex-A lifted (bubblepm analog): EvaluateTropicalMC with
         SplitRealImag vs NIntegrate of the COMPLEX original.  Gate < 2%.
         + Direct-mode lifting still refuses (liftcomplex) — SplitRealImag is the
           only mode that lifts a complex integrand (the message fix).
     (B) DIVERGENT complex-A lifted (1/eps pole), Case A (the imaginary monomial
         phase decouples from the divergent direction), genuinely COMPLEX pole:
         IBP route vs pinned-eps Subtraction route vs NIntegrate complex Laurent.
         Pole agreement < 0.05; recon within ~3% of NIntegrate.
     (C) DIVERGENT complex-A, OFF-AXIS (imaginary phase couples to the divergent
         slot — imaginary part on the pole variable shifts the endpoint exponent to
         s_k = c_k eps + i theta_k, theta_k != 0): the "pole" sits off eps=0, the
         cone is FINITE (pole = 0, value ~ 1/(i theta_k)).  IBP is the UNIQUE route
         that resolves it (planIBPCX.md §1.2: it raises the divergent exponent
         before flattening -> bounded bulk theta-phase).  IBP (pole~0, finite) vs
         NIntegrate of the complex original.  (Subtraction / pinned-eps share the
         vanishing-exponent blind spot and are NOT valid oracles here — §5.1 F2.)

   Independent oracles per numeric PASS: NIntegrate of the complex original AND a
   second engine route (Subtraction vs IBP).  (planAXpDIVv2.md §8.5.)

   Tier 1: WL + g++ (MC).  Run:  wolframscript -file TEST/cc_47.wl
   ============================================================================ *)

$pkgRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{$pkgRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$pkgRoot, "tropical_eval.wl"}]];
Print["CC47: packages loaded."];
Print[];

sci[v_] := If[NumericQ[v], ScientificForm[N[v], 3], v];
$cc47Pass = True; $cc47Fail = {};
cc47Assert[label_String, test_, exp_String, got_String] :=
  If[TrueQ[test],
    Print["CC47 PASS  [", label, "]  expected=", exp, "  got=", got],
    Print["CC47 FAIL  [", label, "]  expected=", exp, "  got=", got];
    $cc47Pass = False; AppendTo[$cc47Fail, label]];

SetAttributes[qr, HoldFirst];
qr[e_] := Quiet[e, {TropicalEval::validate, TropicalEval::divergent, TropicalEval::vegasbudget,
   TropicalEval::liftcomplex, TropicalEval::splitdivmono, General::stop, NIntegrate::slwcon,
   NIntegrate::ncvb, NIntegrate::inumr, NIntegrate::eincr, NIntegrate::izero, NIntegrate::precw,
   TropicalFan::polymake}];

wdir[tag_] := Module[{d = FileNameJoin[{$pkgRoot, "TEST", "INTERFILES", "cc47", tag}]},
  Quiet[CreateDirectory[d, CreateIntermediateDirectories -> True], {CreateDirectory::eexist}]; d];

(* ============================================================================
   PART A — convergent complex-A lifted (bubblepm analog)
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (A) CONVERGENT complex-A lifted  vs NIntegrate"];

Module[{vars, Avals, poly, Breal, us, sub, jac, ig, ref, spec, lift, ls, ld, fan,
        res, val, relA, direct},
  vars  = {x[1], x[2]};
  Avals = {-1/2 + 3/2 I, -1/2 - 3/2 I};        (* complex monomial exponents *)
  poly  = 1 + x[1] x[2] + 10^4 x[1]^2 + x[2]^2;  (* extreme coeff -> triggers lift *)
  Breal = -2;

  (* Oracle 1: NIntegrate of the COMPLEX original. *)
  us = Table[Unique["u"], {2}]; sub = Table[x[i] -> us[[i]]/(1 - us[[i]]), {i, 2}];
  jac = Times @@ Table[1/(1 - us[[i]])^2, {i, 2}];
  ig  = (x[1]^Avals[[1]] x[2]^Avals[[2]] poly^Breal /. sub) jac;
  ref = qr@NIntegrate[ig, Evaluate[Sequence @@ ({#, 0, 1} & /@ us)],
    Method -> "GlobalAdaptive", PrecisionGoal -> 6, WorkingPrecision -> 30, MaxRecursion -> 80];
  Print["   NIntegrate(complex original) = ", N[ref]];

  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> Avals,
    "PolynomialExponents" -> {Breal}, "Variables" -> vars, "KinematicSymbols" -> {},
    "RegulatorSymbol" -> None|>;
  lift = qr@LiftCoefficients[spec, {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 3|>}];
  ls = lift["LiftedSpec"]; ld = lift["LiftData"];
  fan = qr@computeFanScaled[qr@PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]]];
  Print["   lifted fan: ", Length[fan[[2]]], " cones, z0=", ld["z0"]];

  (* Oracle 2: the engine (SplitRealImag) — full codegen + g++ + MC. *)
  res = qr@EvaluateTropicalMC[ls, fan, {{}}, "LiftData" -> ld,
    "ComplexExponentMode" -> "SplitRealImag", "Integrator" -> "MC", "NSamples" -> 3000000,
    "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> wdir["convA"]];
  If[AssociationQ[res] && KeyExistsQ[res, "Results"],
    val = res["Results"][[1]]["Re"] + I res["Results"][[1]]["Im"];
    relA = Abs[(val - ref)/ref];
    Print["   engine SplitRealImag = ", N[val], "   rel-err = ", sci@relA];
    cc47Assert["A conv complex-A lifted SplitRealImag vs NIntegrate",
      NumericQ[relA] && relA < 2*^-2, "rel<2e-2", "rel=" <> ToString[sci@relA]],
    cc47Assert["A conv complex-A lifted SplitRealImag vs NIntegrate", False,
      "Association[Results]", ToString[Head[res]]]];

  (* Direct + lift still refuses a complex integrand (liftcomplex): SplitRealImag
     is the ONLY mode that lifts complex exponents (the message-fix contract). *)
  direct = qr@EvaluateTropicalMC[ls, fan, {{}}, "LiftData" -> ld,
    "ComplexExponentMode" -> "Direct", "Integrator" -> "MC", "NSamples" -> 100000,
    "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> wdir["convA_direct"]];
  cc47Assert["A Direct+lift refuses complex-A (liftcomplex)",
    direct === $Failed, "$Failed", ToString[Head[direct]]];
];

(* ============================================================================
   PART B — divergent complex-A lifted, Case A (genuinely COMPLEX pole):
            IBP vs Subtraction vs NIntegrate
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (B) DIVERGENT complex-A lifted (Case A)  IBP vs Sub vs NIntegrate"];

Module[{epsD, esStar, dim, polys, pe, Avals, spec, lift, ls, ld, fan,
        niCx, gv, fd, poleNI, finNI, reconNI, sub, ibp, poleS, finS, poleI, finI, recS, recI},
  epsD = Symbol["epsCC47B"]; esStar = 0.01; dim = 2;
  (* 1/eps pole on x1; the COMPLEX monomial exponent on x2 decouples from the
     divergent direction (Case A) and produces a genuinely complex pole. *)
  polys = {1 + x[1] x[2] + 10^4 x[1]^2 + x[2]^2}; pe = {-2};
  Avals = {-1 + epsD, -1/2 - 3/2 I};

  (* Oracle 1: NIntegrate complex Laurent. *)
  niCx[es_] := Module[{us, subr, jac, igc, Aev}, Aev = Avals /. epsD -> es;
    us = Table[Unique["u"], {dim}]; subr = Table[x[i] -> us[[i]]/(1 - us[[i]]), {i, dim}];
    jac = Times @@ Table[1/(1 - us[[i]])^2, {i, dim}];
    igc = (Times @@ MapThread[#1^#2 &, {Table[x[i], {i, dim}], Aev}] *
           Times @@ MapThread[#1^#2 &, {polys, pe}] /. subr) jac;
    qr@NIntegrate[igc, Evaluate[Sequence @@ ({#, 0, 1} & /@ us)], Method -> "GlobalAdaptive",
      PrecisionGoal -> 5, WorkingPrecision -> 25, MaxRecursion -> 50]];
  gv = MapThread[#1 #2 &, {N[{esStar, 2 esStar, 4 esStar}, 25], niCx /@ N[{esStar, 2 esStar, 4 esStar}, 25]}];
  fd = MapThread[{#1, #2} &, {N[{esStar, 2 esStar, 4 esStar}, 25], gv}];
  poleNI = (Fit[{#1, Re[#2]} & @@@ fd, {1, ee}, ee] /. ee -> 0) +
           I (Fit[{#1, Im[#2]} & @@@ fd, {1, ee}, ee] /. ee -> 0);
  finNI  = Coefficient[Fit[{#1, Re[#2]} & @@@ fd, {1, ee}, ee], ee] +
           I Coefficient[Fit[{#1, Im[#2]} & @@@ fd, {1, ee}, ee], ee];
  reconNI = poleNI/esStar + finNI;
  Print["   NIntegrate: pole=", N@poleNI, "  recon=", N@reconNI];

  spec = <|"Polynomials" -> polys, "MonomialExponents" -> Avals, "PolynomialExponents" -> pe,
    "Variables" -> {x[1], x[2]}, "KinematicSymbols" -> {}, "RegulatorSymbol" -> epsD|>;
  lift = qr@LiftCoefficients[spec, {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 3|>}];
  ls = lift["LiftedSpec"]; ld = lift["LiftData"];
  fan = qr@computeFanScaled[qr@PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]]];

  (* Oracle 2: pinned-eps Subtraction route. *)
  sub = qr@LaurentFromSubtraction[ls, fan, {{}}, "LiftData" -> ld, "ComplexExponentMode" -> "SplitRealImag",
    "EpsilonValues" -> {esStar, 2 esStar, 4 esStar}, "Integrator" -> "MC", "NSamples" -> 800000,
    "WorkingDirectory" -> wdir["divB_sub"], "Verbose" -> False];
  (* Oracle 3: IBP route. *)
  ibp = qr@EvaluateTropicalMC[ls, fan, {{}}, "LiftData" -> ld, "Method" -> "IBP",
    "ComplexExponentMode" -> "SplitRealImag", "Integrator" -> "MC", "NSamples" -> 800000,
    "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> wdir["divB_ibp"]];

  If[!(AssociationQ[sub] && KeyExistsQ[sub, "Results"]),
    cc47Assert["B Sub route", False, "Association[Results]", ToString[Head[sub]]]; Return[]];
  If[!(AssociationQ[ibp] && KeyExistsQ[ibp, "Results"]),
    cc47Assert["B IBP route", False, "Association[Results]", ToString[Head[ibp]]]; Return[]];
  poleS = sub["Results"][[1]]["Pole"]; finS = sub["Results"][[1]]["Finite"]; recS = poleS/esStar + finS;
  poleI = ibp["Results"][[1]]["PoleCoefficient"]; finI = ibp["Results"][[1]]["FinitePart"]; recI = poleI/esStar + finI;
  Print["   Sub: pole=", N@poleS, "  recon=", N@recS];
  Print["   IBP: pole=", N@poleI, "  recon=", N@recI];

  cc47Assert["B IBP pole == Sub pole",
    NumericQ[Abs[poleI - poleS]] && Abs[poleI - poleS] < 0.05,
    "|pole_IBP-pole_Sub|<0.05", "|diff|=" <> ToString[sci@Abs[poleI - poleS]]];
  cc47Assert["B Sub recon vs NIntegrate",
    NumericQ[Abs[recS - reconNI]] && Abs[(recS - reconNI)/reconNI] < 0.03,
    "rel<0.03", "rel=" <> ToString[sci@Abs[(recS - reconNI)/reconNI]]];
  cc47Assert["B IBP recon vs NIntegrate",
    NumericQ[Abs[recI - reconNI]] && Abs[(recI - reconNI)/reconNI] < 0.03,
    "rel<0.03", "rel=" <> ToString[sci@Abs[(recI - reconNI)/reconNI]]];
  cc47Assert["B genuinely complex pole (Im(pole) != 0)",
    NumericQ[Im[poleNI]] && Abs[Im[poleNI]] > 0.01,
    "|Im(pole)|>0.01", "Im(pole_NI)=" <> ToString[sci@Im[poleNI]]];
];

(* ============================================================================
   PART C — divergent complex-A, OFF-AXIS (planIBPCX.md §1/§3): the imaginary
   phase couples to the divergent direction (imaginary part 3/2 ON the pole
   variable x1), so the divergent endpoint exponent is s_k = c_k eps + i theta_k
   with theta_k != 0.  The "1/eps pole" is then REGULAR at eps=0 (it sits at
   eps = -i theta_k/c_k, off the real axis): the cone is FINITE, pole = 0, value
   ~ 1/(i theta_k).  This was formerly refused as "complex Case B"
   (splitdivmono -> $Failed); planIBPCX.md shows IBP is the UNIQUE route that
   resolves it (it RAISES the divergent exponent before flattening, so the bulk
   theta-phase coefficient is O(theta), bounded — the inline-Subtraction and
   pinned-eps routes flatten by the vanishing exponent and oscillate too fast).
   Oracle: NIntegrate of the COMPLEX original (Subtraction is NOT a valid oracle
   here — same vanishing-exponent blind spot; planIBPCX.md §5.1 F2).
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (C) DIVERGENT complex-A OFF-AXIS  IBP (pole~0, finite) vs NIntegrate"];

Module[{epsD, esStar = 0.01, dim = 2, polys, pe, Avals, spec, lift, ls, ld, fan,
        niC, gv, fd, finNI, ibp, poleI, finI, imagPoles, ip},
  epsD = Symbol["epsCC47C"];
  polys = {1 + x[1] x[2] + 10^4 x[1]^2 + x[2]^2}; pe = {-2};
  Avals = {-1 + epsD + 3/2 I, 0};   (* imaginary part ON the pole variable x1 -> off-axis *)

  spec = <|"Polynomials" -> polys, "MonomialExponents" -> Avals, "PolynomialExponents" -> pe,
    "Variables" -> {x[1], x[2]}, "KinematicSymbols" -> {}, "RegulatorSymbol" -> epsD|>;
  lift = qr@LiftCoefficients[spec, {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 3|>}];
  ls = lift["LiftedSpec"]; ld = lift["LiftData"];
  fan = qr@computeFanScaled[qr@PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]]];

  (* Oracle: NIntegrate of the COMPLEX original.  Off-axis -> no pole; the value
     is regular at eps=0, so fit I(eps)=f+c*eps and take the eps^0 part f. *)
  niC[es_] := Module[{us, subr, jac, igc, Aev}, Aev = Avals /. epsD -> es;
    us = Table[Unique["u"], {dim}]; subr = Table[x[i] -> us[[i]]/(1 - us[[i]]), {i, dim}];
    jac = Times @@ Table[1/(1 - us[[i]])^2, {i, dim}];
    igc = (Times @@ MapThread[#1^#2 &, {Table[x[i], {i, dim}], Aev}] *
           Times @@ MapThread[#1^#2 &, {polys, pe}] /. subr) jac;
    qr@NIntegrate[igc, Evaluate[Sequence @@ ({#, 0, 1} & /@ us)], Method -> "GlobalAdaptive",
      PrecisionGoal -> 6, WorkingPrecision -> 30, MaxRecursion -> 60]];
  gv = niC /@ N[{esStar, 2 esStar, 4 esStar}, 30];
  fd = MapThread[{#1, #2} &, {N[{esStar, 2 esStar, 4 esStar}, 30], gv}];
  finNI = (Fit[{#1, Re[#2]} & @@@ fd, {1, ee}, ee] /. ee -> 0) +
          I (Fit[{#1, Im[#2]} & @@@ fd, {1, ee}, ee] /. ee -> 0);
  Print["   NIntegrate finite (I at eps->0) = ", N@finNI, "  (values bounded -> no pole)"];

  ibp = qr@EvaluateTropicalMC[ls, fan, {{}}, "LiftData" -> ld, "Method" -> "IBP",
    "ComplexExponentMode" -> "SplitRealImag", "Integrator" -> "MC", "NSamples" -> 2000000,
    "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> wdir["divC_ibp"]];
  If[!(AssociationQ[ibp] && KeyExistsQ[ibp, "Results"]),
    cc47Assert["C IBP off-axis returns a result", False, "Association[Results]",
      ToString[Head[ibp]]]; Return[]];
  poleI = ibp["Results"][[1]]["PoleCoefficient"]; finI = ibp["Results"][[1]]["FinitePart"];
  (* Confirm the fixture REALLY realizes the off-axis condition (planIBPCX.md §9
     risk #3): at least one divergent sector must carry an off-axis theta != 0.
     A single off-axis direction now routes through the generalized corner path
     (planIBPMULTIDIV.md), which records OffAxisThetas instead of a scalar
     ImagPole; accept either representation. *)
  imagPoles = Flatten[(If[TrueQ[#["MultiDiv"]], #["OffAxisThetas"],
                          {#["ImagPole"]}]) & /@ ibp["IBPProcessedSectors"]];
  ip = AnyTrue[imagPoles, (NumericQ[#] && Abs[#] > 0.01) || (!NumericQ[#] && !TrueQ[PossibleZeroQ[#]]) &];
  Print["   IBP: pole=", N@poleI, "  finite=", N@finI, "   ImagPoles=", imagPoles];

  cc47Assert["C fixture realizes off-axis (ImagPole != 0)", ip,
    "some ImagPole != 0", "ImagPoles=" <> ToString[imagPoles, InputForm]];
  cc47Assert["C IBP off-axis pole ~ 0 (cone is finite, no 1/eps)",
    NumericQ[Abs[poleI]] && Abs[poleI] < 0.02,
    "|pole_IBP|<0.02", "|pole|=" <> ToString[sci@Abs[poleI]]];
  cc47Assert["C IBP off-axis finite vs NIntegrate complex original",
    NumericQ[Abs[finI - finNI]] && Abs[(finI - finNI)/finNI] < 0.03,
    "rel<0.03", "rel=" <> ToString[sci@Abs[(finI - finNI)/finNI]]];
];

Print[];
Print["================================================================"];
If[$cc47Pass,
  Print["CC47 PASS  complex MONOMIAL exponents A x lifting verified: SplitRealImag ",
        "clears liftcomplex; convergent matches NIntegrate; divergent IBP & ",
        "Subtraction agree (Case A); off-axis IBP (pole~0,finite) matches NIntegrate ",
        "(planIBPCX.md).  failures={}"],
  Print["CC47 FAIL  failed: ", $cc47Fail]];
Print["================================================================"];
If[!$cc47Pass, Quit[1]];
