(* ============================================================================
   TEST/cc_51.wl  —  Cross-check #51 (planIBPMULTIDIV.md §6.3)
   IBP on a multi-soft sector: ONE real 1/eps pole + ONE off-axis direction in
   the SAME cone.  Before planIBPMULTIDIV the corner cone (both rays soft) was
   refused (TropicalEval::nestedIBP counted the off-axis direction as a second
   pole); it is now resolved by the generalized 2^(N+1)-corner iterated IBP.

     I(eps) = int_0^inf int_0^inf x1^{2eps-1} x2^{i th -1} (1+x1+x2)^{-B} dx
            = Gamma(2 eps) Gamma(i th) Gamma(B - 2 eps - i th) / Gamma(B)   (Dirichlet)

   x1 is a genuine real-soft pole (alpha1 = 2 eps : Re=0, Im=0); x2 is off-axis
   (alpha2 = i th : Re=0, Im=th != 0).  Gamma(2 eps) -> the SINGLE 1/eps pole;
   Gamma(i th) is a finite COMPLEX constant.  Exactly one pole + one finite
   off-axis direction (planIBPMULTIDIV.md §1.2).

   Oracles per numeric PASS (invariant #4):
     (1) the closed-form Dirichlet Series (shares NO code with the engine), and
     (2) Direct == SplitRealImag on the IBP path.
   Plus a STRUCTURAL self-check (planIBPMULTIDIV.md §6.3): the corner sector must
   realize nPole==1 AND an OffAxis direction (a fixture that silently avoids the
   co-located case would prove nothing).

   Tier 1: WL + g++ (MC).  Run:  wolframscript -file TEST/cc_51.wl
   ============================================================================ *)

$pkgRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{$pkgRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$pkgRoot, "tropical_eval.wl"}]];
Print["CC51: packages loaded."];
Print[];

sci[v_] := If[NumericQ[v], ScientificForm[N[v], 3], v];
$cc51Pass = True; $cc51Fail = {};
cc51Assert[label_String, test_, exp_String, got_String] :=
  If[TrueQ[test],
    Print["CC51 PASS  [", label, "]  expected=", exp, "  got=", got],
    Print["CC51 FAIL  [", label, "]  expected=", exp, "  got=", got];
    $cc51Pass = False; AppendTo[$cc51Fail, label]];

SetAttributes[qr, HoldFirst];
qr[e_] := Quiet[e, {TropicalEval::validate, TropicalEval::divergent,
   TropicalEval::vegasbudget, TropicalEval::splitdivmono, General::stop,
   TropicalFan::polymake, Power::infy}];

wdir[tag_] := Module[{d = FileNameJoin[{$pkgRoot, "TEST", "INTERFILES", "cc51", tag}]},
  Quiet[CreateDirectory[d, CreateIntermediateDirectories -> True],
        {CreateDirectory::eexist}]; d];

(* ---- The closed-form Dirichlet oracle (no engine code) ---- *)
th = 1/2; Bexp = 3;
orc[e_] := Gamma[2 e] Gamma[I th] Gamma[Bexp - 2 e - I th]/Gamma[Bexp];
ser       = Series[orc[ee], {ee, 0, 0}];
poleOrc   = N[SeriesCoefficient[ser, -1]];
finOrc    = N[SeriesCoefficient[ser, 0]];
Print["  Dirichlet oracle: pole = ", N@poleOrc, "  finite = ", N@finOrc];

(* ---- The integrand spec ---- *)
epsR  = Symbol["epsCC51"];
polys = {1 + x[1] + x[2]};
pe    = {-Bexp};
Avals = {2 epsR - 1, I th - 1};   (* x1 real-soft pole; x2 off-axis soft *)
spec  = <|"Polynomials" -> polys, "MonomialExponents" -> Avals,
  "PolynomialExponents" -> pe, "Variables" -> {x[1], x[2]},
  "KinematicSymbols" -> {}, "RegulatorSymbol" -> epsR|>;
fan = qr@computeFanScaled[
  qr@PolytopeVertices[(Times @@ polys)^(-1), {x[1], x[2]}]];

runIBP[mode_String, tag_String, ns_Integer] := qr@EvaluateTropicalMC[
  spec, fan, {{}}, "Method" -> "IBP", "ComplexExponentMode" -> mode,
  "Integrator" -> "MC", "NSamples" -> ns, "RunChecks" -> False,
  "Verbose" -> False, "WorkingDirectory" -> wdir[tag]];

(* ============================================================================
   PART A — Direct mode vs the Dirichlet oracle + structural self-check
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (A) Direct IBP  vs Dirichlet oracle  (+ nPole==1, exists OffAxis)"];

Module[{rD, poleD, finD, procs, mdSector, hasStruct},
  rD = runIBP["Direct", "direct", 3000000];
  If[!(AssociationQ[rD] && KeyExistsQ[rD, "Results"]),
    cc51Assert["A Direct returns a result (no nestedIBP refusal)", False,
      "Association[Results]", ToString[Head[rD]]]; Return[]];
  poleD = rD["Results"][[1]]["PoleCoefficient"];
  finD  = rD["Results"][[1]]["FinitePart"];
  Print["   Direct IBP: pole = ", N@poleD, "  finite = ", N@finD];

  (* Structural self-check: a divergent sector must be a one-pole + off-axis
     corner (MultiDiv, HasPole, OffAxisDirs != {}). *)
  procs = rD["IBPProcessedSectors"];
  mdSector = SelectFirst[procs, TrueQ[#["MultiDiv"]] &&
               TrueQ[#["HasPole"]] && Length[#["OffAxisDirs"]] >= 1 &, None];
  hasStruct = (mdSector =!= None);
  Print["   corner sector: ", If[hasStruct,
     "MultiDiv HasPole=True OffAxisDirs=" <> ToString[mdSector["OffAxisDirs"]],
     "NOT FOUND"]];
  cc51Assert["A structural: one-pole + off-axis corner present",
    hasStruct, "MultiDiv & HasPole & exists OffAxis", ToString[hasStruct]];

  cc51Assert["A Direct pole vs Dirichlet oracle",
    NumericQ[Abs[poleD - poleOrc]] && Abs[(poleD - poleOrc)/poleOrc] < 0.02,
    "rel<0.02", "rel=" <> ToString[sci@Abs[(poleD - poleOrc)/poleOrc]]];
  cc51Assert["A Direct finite vs Dirichlet oracle",
    NumericQ[Abs[finD - finOrc]] && Abs[(finD - finOrc)/finOrc] < 0.02,
    "rel<0.02", "rel=" <> ToString[sci@Abs[(finD - finOrc)/finOrc]]];
  cc51Assert["A pole genuinely complex (Im != 0)",
    Abs[Im[poleOrc]] > 0.01, "|Im(pole)|>0.01",
    "Im=" <> ToString[sci@Im[poleOrc]]];
  $ccDirect = {poleD, finD};
];

(* ============================================================================
   PART B — SplitRealImag mode vs oracle, and Direct == SplitRealImag
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (B) SplitRealImag IBP  vs oracle  and  Direct == SplitRealImag"];

Module[{rS, poleS, finS},
  rS = runIBP["SplitRealImag", "split", 3000000];
  If[!(AssociationQ[rS] && KeyExistsQ[rS, "Results"]),
    cc51Assert["B SplitRealImag returns a result", False,
      "Association[Results]", ToString[Head[rS]]]; Return[]];
  poleS = rS["Results"][[1]]["PoleCoefficient"];
  finS  = rS["Results"][[1]]["FinitePart"];
  Print["   SplitRealImag IBP: pole = ", N@poleS, "  finite = ", N@finS];

  cc51Assert["B Split pole vs Dirichlet oracle",
    NumericQ[Abs[poleS - poleOrc]] && Abs[(poleS - poleOrc)/poleOrc] < 0.02,
    "rel<0.02", "rel=" <> ToString[sci@Abs[(poleS - poleOrc)/poleOrc]]];
  cc51Assert["B Split finite vs Dirichlet oracle",
    NumericQ[Abs[finS - finOrc]] && Abs[(finS - finOrc)/finOrc] < 0.02,
    "rel<0.02", "rel=" <> ToString[sci@Abs[(finS - finOrc)/finOrc]]];

  If[ListQ[$ccDirect],
    cc51Assert["B Direct == SplitRealImag (pole)",
      Abs[poleS - $ccDirect[[1]]] < 0.02 * (Abs[poleOrc] + 0.01),
      "|dPole|<2% scale", "|d|=" <> ToString[sci@Abs[poleS - $ccDirect[[1]]]]];
    cc51Assert["B Direct == SplitRealImag (finite)",
      Abs[finS - $ccDirect[[2]]] < 0.02 * (Abs[finOrc] + 0.01),
      "|dFinite|<2% scale", "|d|=" <> ToString[sci@Abs[finS - $ccDirect[[2]]]]];
  ];
];

(* ============================================================================
   PART C — off-axis direction whose REAL part carries eps (c_m != 0): exercises
   the off-axis-prefactor O(eps) correction (K1 = K0 sum_j(-c_j/(i theta_j))),
   which Parts A/B (eps-free off-axis, c_m=0 -> K1=0) do NOT touch.
     I(eps) = Gamma(2eps) Gamma(2eps + i th) Gamma(B - 4eps - i th)/Gamma(B)
   x2 is off-axis (Re=2eps -> 0, Im=th) with c_m = 2 != 0.
   ============================================================================ *)
Print["----------------------------------------------------------------"];
Print["  (C) off-axis with c_m != 0 (exercises the K1 finite correction)"];

Module[{epsC, AvalsC, specC, fanC, orcC, serC, poleOrcC, finOrcC, rC, poleC, finC,
        ck},
  epsC = Symbol["epsCC51C"];
  AvalsC = {2 epsC - 1, 2 epsC + I th - 1};   (* x2 off-axis, c_m = 2 *)
  orcC[e_] := Gamma[2 e] Gamma[2 e + I th] Gamma[Bexp - 4 e - I th]/Gamma[Bexp];
  serC = Series[orcC[ee], {ee, 0, 0}];
  poleOrcC = N[SeriesCoefficient[serC, -1]];
  finOrcC  = N[SeriesCoefficient[serC, 0]];
  Print["   oracle: pole = ", N@poleOrcC, "  finite = ", N@finOrcC];

  specC = <|"Polynomials" -> polys, "MonomialExponents" -> AvalsC,
    "PolynomialExponents" -> pe, "Variables" -> {x[1], x[2]},
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> epsC|>;
  fanC = qr@computeFanScaled[
    qr@PolytopeVertices[(Times @@ polys)^(-1), {x[1], x[2]}]];
  rC = qr@EvaluateTropicalMC[specC, fanC, {{}}, "Method" -> "IBP",
    "ComplexExponentMode" -> "Direct", "Integrator" -> "MC", "NSamples" -> 4000000,
    "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> wdir["cm"]];
  If[!(AssociationQ[rC] && KeyExistsQ[rC, "Results"]),
    cc51Assert["C off-axis c_m!=0 returns a result", False,
      "Association[Results]", ToString[Head[rC]]]; Return[]];
  poleC = rC["Results"][[1]]["PoleCoefficient"];
  finC  = rC["Results"][[1]]["FinitePart"];
  (* confirm a nonzero off-axis c_m is actually present (K1 != 0 path) *)
  ck = Flatten[#["OffAxisCks"] & /@
        Select[rC["IBPProcessedSectors"], TrueQ[#["MultiDiv"]] &]];
  Print["   engine: pole = ", N@poleC, "  finite = ", N@finC,
        "   OffAxisCks=", ck];
  cc51Assert["C off-axis c_m != 0 present (K1 path exercised)",
    AnyTrue[ck, (NumericQ[#] && Abs[#] > 0.01) &],
    "some c_m != 0", "OffAxisCks=" <> ToString[ck, InputForm]];
  cc51Assert["C pole vs Dirichlet oracle (c_m!=0)",
    NumericQ[Abs[poleC - poleOrcC]] && Abs[(poleC - poleOrcC)/poleOrcC] < 0.02,
    "rel<0.02", "rel=" <> ToString[sci@Abs[(poleC - poleOrcC)/poleOrcC]]];
  cc51Assert["C finite vs Dirichlet oracle (c_m!=0, K1 term)",
    NumericQ[Abs[finC - finOrcC]] && Abs[(finC - finOrcC)/finOrcC] < 0.02,
    "rel<0.02", "rel=" <> ToString[sci@Abs[(finC - finOrcC)/finOrcC]]];
];

Print[];
Print["================================================================"];
If[$cc51Pass,
  Print["CC51 PASS  one real pole + one off-axis IBP verified: corner cone no ",
        "longer refused (nestedIBP); pole & finite (Re and Im) match the closed-",
        "form Dirichlet oracle; Direct == SplitRealImag.  failures={}"],
  Print["CC51 FAIL  failed: ", $cc51Fail]];
Print["================================================================"];
If[!$cc51Pass, Quit[1]];
