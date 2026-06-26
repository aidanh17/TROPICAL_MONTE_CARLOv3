(* ============================================================================
   EXAMPLES/Ex-V3a.wl  —  Worked Example V3a
   plan.md §8.4: "Divergent + complex, full Laurent via IBP,
                  cross-checked vs subtraction + exact Gamma"
   Cross-checks: #8 (FIESTA, optional/skipped if absent), #33, #35

   Integrand
   ---------
   I(eps) = Int_{[0,inf)^2} x1^{2*eps - 1} (1 + x1 + x2)^{-B} dx1 dx2
   with  B = 2 + I/2  (complex polynomial exponent, Im B != 0).

   Exact Laurent (derivation from §Mellin/Gamma, see cc_33.wl):
     Step 1: integrate x2:  J(x1) = (1+x1)^{1-B} / (B-1)    [Re B > 1]
     Step 2: 1D integral = Beta(2*eps, B-1-2*eps) / (B-1)
     so
       I(eps) = Gamma(2*eps) Gamma(B-1-2*eps) / ((B-1) Gamma(B-1))
       pole   = residue_at_eps=0  =  1 / (2*(B-1))
       finite = d/d(eps) [eps * I(eps)] |_{eps=0}
              = [-EulerGamma/2 - PolyGamma[0, B-1]] / (B-1)

   Three independent routes for (pole, finite)
   -------------------------------------------
     Route A — EvaluateTropicalMCIBP  (IBP @ eps=0, MC sampler)
     Route B — LaurentFromSubtraction (finite-eps fit, MC sampler)
     Route C — Exact Gamma Laurent    (closed-form reference)

   PASS criteria (all must hold)
   ------------------------------
     P1  |Re(pole_A) - Re(poleRef)| < 5e-3  AND  |Im(pole_A) - Im(poleRef)| < 5e-3
     P2  |Re(pole_B) - Re(poleRef)| < 1e-2  AND  |Im(pole_B) - Im(poleRef)| < 1e-2
     P3  |pole_A - pole_B|  < 1e-2   (IBP vs fit agree: cross-check #33 A3 analogue)
     P4  |finite_A - finite_B| < 2e-2

   Tier 1: WL + g++ only (no external deps).
   CUBA-Vegas columns printed informally when CUBA is present but never gate PASS.
   FIESTA cross-check (#8) is skipped when FIESTA is absent.

   Run
   ---
     cd /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3
     wolframscript -file EXAMPLES/Ex-V3a.wl
   exits 0 on PASS, 1 on FAIL.
   ============================================================================ *)

(* ── 0.  Locate packages ──────────────────────────────────────────────────── *)

$v3Root  = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
$fanWL   = FileNameJoin[{$v3Root, "tropical_fan.wl"}];
$evalWL  = FileNameJoin[{$v3Root, "tropical_eval.wl"}];

If[!FileExistsQ[$fanWL],
  Print["Ex-V3a FAIL  tropical_fan.wl not found at ", $fanWL]; Quit[1]];
If[!FileExistsQ[$evalWL],
  Print["Ex-V3a FAIL  tropical_eval.wl not found at ", $evalWL]; Quit[1]];

Get[$fanWL];
Get[$evalWL];

Quiet[CreateDirectory[
  FileNameJoin[{$v3Root, "INTERFILES", "ExV3a"}],
  CreateIntermediateDirectories -> True],
 {CreateDirectory::eexist}];

$ioDir = FileNameJoin[{$v3Root, "INTERFILES", "ExV3a"}];

Print[""];
Print["================================================================"];
Print["  EXAMPLES/Ex-V3a: Divergent + Complex — Full Laurent via IBP  "];
Print["  Cross-checks #33, #35 (plan.md §8.4)                        "];
Print["================================================================"];
Print[];

(* ── 1.  Helpers ──────────────────────────────────────────────────────────── *)

SetAttributes[quietRun, HoldFirst];
quietRun[expr_] := Quiet[expr,
  {TropicalEval`EvaluateTropicalMCIBP::validate,
   TropicalEval`LaurentFromSubtraction::validate,
   General::stop}];

fmtC[z_] := Module[{r = Re[z], im = Im[z]},
  If[Abs[im] < 1*^-12,
    ToString[NumberForm[N[r], {5, 5}]],
    ToString[NumberForm[N[r], {5, 5}]] <> "+"
      <> ToString[NumberForm[N[im], {5, 5}]] <> "I"]];

sci[x_] := If[N[Abs[x]] == 0., "0",
  Module[{e = Floor[Log10[Abs[N[x]]]]},
    ToString[NumberForm[N[x] / 10^e, {3, 2}]] <> "e" <> ToString[e]]];

(* Accumulate pass/fail *)
$passes   = {};
$failures = {};

record[name_String, ok_, passDetail_String : "", failExtra_String : ""] :=
  If[TrueQ[ok],
    (AppendTo[$passes, name];
     Print["  Ex-V3a PASS  ", name,
           If[passDetail =!= "", "  " <> passDetail, ""]]),
    (AppendTo[$failures, name];
     Print["  Ex-V3a FAIL  ", name,
           If[failExtra =!= "", "  " <> failExtra, ""]])];

(* ── 2.  Integrand specification ─────────────────────────────────────────── *)

Print["--- §1  Integrand and exact reference ---"];
eps    = Symbol["epsV3a"];
Bval   = 2 + 1/2 * I;     (* complex polynomial exponent *)
vars   = {x[1], x[2]};
poly   = 1 + x[1] + x[2];

spec = <|
  "Polynomials"         -> {poly},
  "MonomialExponents"   -> {2 eps - 1, 0},
  "PolynomialExponents" -> {-Bval},
  "Variables"           -> vars,
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> eps
|>;

(* Exact Laurent coefficients *)
poleRef = 1 / (2 (Bval - 1));
finRef  = N[(-EulerGamma/2 - PolyGamma[0, Bval - 1]) / (Bval - 1)];

Print["  B       = ", Bval];
Print["  poleRef = 1/(2(B-1)) = ", N[poleRef]];
Print["  finRef  = ", N[finRef]];
Print[];

(* ── 3.  Fan construction ─────────────────────────────────────────────────── *)

Print["--- §2  Tropical fan ---"];
verts = PolytopeVertices[poly^(-2), vars];
fan   = computeFanScaled[verts];
If[fan === $Failed || !ListQ[fan],
  Print["Ex-V3a FAIL  fan construction returned $Failed"];
  Quit[1]];
Print["  Sectors: ", Length[fan[[2]]]];
Print[];

(* ── 4.  Route A — IBP + MC ───────────────────────────────────────────────── *)

Print["--- §3  Route A: IBP @ eps=0 (EvaluateTropicalMCIBP, MC) ---"];
resA = quietRun @ EvaluateTropicalMCIBP[spec, fan, {{}},
  "Integrator"       -> "MC",
  "NSamples"         -> 2000000,
  "RunChecks"        -> False,
  "Verbose"          -> False,
  "WorkingDirectory" -> FileNameJoin[{$ioDir, "routeA"}]];

If[!AssociationQ[resA] || !KeyExistsQ[resA, "Results"],
  Print["Ex-V3a FAIL  EvaluateTropicalMCIBP returned $Failed: ", resA];
  Quit[1]];

poleA = resA["Results"][[1]]["PoleCoefficient"];
finA  = resA["Results"][[1]]["FinitePart"];
Print["  IBP-MC:  pole = ", fmtC[poleA], "  finite = ", fmtC[finA]];
Print[];

(* ── 5.  Route B — LaurentFromSubtraction (finite-eps fit) ────────────────── *)

Print["--- §4  Route B: Laurent from subtraction (finite-eps fit) ---"];
epsVals = {0.01, 0.02, 0.04};
resB = quietRun @ LaurentFromSubtraction[spec, fan, {{}},
  "Integrator"       -> "MC",
  "NSamples"         -> 2000000,
  "EpsilonValues"    -> epsVals,
  "WorkingDirectory" -> FileNameJoin[{$ioDir, "routeB"}]];

If[!AssociationQ[resB] || !KeyExistsQ[resB, "Results"],
  Print["Ex-V3a FAIL  LaurentFromSubtraction returned $Failed: ", resB];
  Quit[1]];

poleB = resB["Results"][[1]]["Pole"];
finB  = resB["Results"][[1]]["Finite"];
Print["  Fit-MC:  pole = ", fmtC[poleB], "  finite = ", fmtC[finB]];
Print[];

(* ── 6.  Optional CUBA-Vegas columns (informational) ─────────────────────── *)

$cuba = TrueQ[TropicalEval`detectCuba[]["Found"]];
If[$cuba,
  Print["--- §5  CUBA-Vegas columns (informational) ---"];
  resAVeg = quietRun @ EvaluateTropicalMCIBP[spec, fan, {{}},
    "Integrator"       -> "VEGAS",
    "NSamples"         -> 2000000,
    "RunChecks"        -> False,
    "Verbose"          -> False,
    "WorkingDirectory" -> FileNameJoin[{$ioDir, "routeA_veg"}]];
  If[AssociationQ[resAVeg] && KeyExistsQ[resAVeg, "Results"],
    poleAV = resAVeg["Results"][[1]]["PoleCoefficient"];
    finAV  = resAVeg["Results"][[1]]["FinitePart"];
    Print["  IBP-VEGAS: pole = ", fmtC[poleAV], "  finite = ", fmtC[finAV]];
    Print["  |IBP-MC - IBP-VEGAS| (pole) = ",
          sci @ Abs[poleA - poleAV]]];
  Print[];
];

(* ── 7.  FIESTA cross-check (#8) — optional, skip if absent ─────────────── *)

$fiestaPresent = FileExistsQ[FileNameJoin[{$v3Root, "CC", "fiesta_common.wl"}]];
If[$fiestaPresent,
  Print["--- §6  FIESTA cross-check (#8) ---"];
  Get[FileNameJoin[{$v3Root, "CC", "fiesta_common.wl"}]];
  (* Run the FIESTA check if the helper exports RunFIESTACheck *)
  If[MemberQ[Names["CC`*"], "CC`RunFIESTACheck"],
    CC`RunFIESTACheck[spec, fan, poleA, finA],
    Print["  FIESTA helper found but RunFIESTACheck not defined — skipping."]];
  Print[],
  Print["--- §6  FIESTA cross-check (#8) SKIPPED (CC/fiesta_common.wl absent) ---"];
  Print[];
];

(* ── 8.  Assertions ───────────────────────────────────────────────────────── *)

Print["--- §7  Assertions ---"];

(* P1: IBP pole vs exact Gamma Laurent *)
Module[{errRe = Abs[Re[poleA] - Re[poleRef]],
        errIm = Abs[Im[poleA] - Im[poleRef]]},
  record["P1 IBP pole vs exact  (#33 A1 analogue)",
    errRe < 5*^-3 && errIm < 5*^-3,
    "Re err=" <> sci[errRe] <> "  Im err=" <> sci[errIm],
    "expected Re<5e-3 Im<5e-3  got Re=" <> sci[errRe] <> " Im=" <> sci[errIm]]];

(* P2: Fit pole vs exact Gamma Laurent *)
Module[{errRe = Abs[Re[poleB] - Re[poleRef]],
        errIm = Abs[Im[poleB] - Im[poleRef]]},
  record["P2 Fit pole vs exact  (#35 A1c analogue)",
    errRe < 1*^-2 && errIm < 1*^-2,
    "Re err=" <> sci[errRe] <> "  Im err=" <> sci[errIm],
    "expected Re<1e-2 Im<1e-2  got Re=" <> sci[errRe] <> " Im=" <> sci[errIm]]];

(* P3: IBP vs fit agree on pole (#33 A3 analogue) *)
Module[{dpole = Abs[poleA - poleB]},
  record["P3 IBP vs Fit pole agreement  (#33 A3, #35 A3a)",
    dpole < 1*^-2,
    "|dpole|=" <> sci[dpole],
    "expected<1e-2  got=" <> sci[dpole]]];

(* P4: IBP vs fit agree on finite part (#33 A4 analogue) *)
Module[{dfin = Abs[finA - finB]},
  record["P4 IBP vs Fit finite agreement  (#33 A4)",
    dfin < 2*^-2,
    "|dfin|=" <> sci[dfin],
    "expected<2e-2  got=" <> sci[dfin]]];

Print[];

(* ── 9.  Summary ──────────────────────────────────────────────────────────── *)

Print["================================================================"];
Print["  Results summary:"];
Print["    B         = ", Bval];
Print["    pole_ref  = ", fmtC[N[poleRef]]];
Print["    pole_IBP  = ", fmtC[poleA]];
Print["    pole_Fit  = ", fmtC[poleB]];
Print["    fin_ref   = ", fmtC[finRef]];
Print["    fin_IBP   = ", fmtC[finA]];
Print["    fin_Fit   = ", fmtC[finB]];
Print["  PASSED: ", Length[$passes]];
Print["  FAILED: ", Length[$failures]];
Print["================================================================"];

If[Length[$failures] == 0,
  Print["Ex-V3a PASS  divergent + complex Laurent: IBP and Fit agree with",
        " exact Gamma  (", Length[$passes], " criteria)"],
  Print["Ex-V3a FAIL  ", Length[$failures],
        " criterion(a) failed: ",
        StringRiffle[$failures, ", "]]];
Print["================================================================"];

Quit[If[Length[$failures] > 0, 1, 0]]
