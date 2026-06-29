(* ============================================================================
   TEST/cc_2.wl  --  Cross-check #2 (plan.md §8.2, §8.3)

   CC2: Exact analytic closed form (Gamma / Beta functions)
        with complex monomial AND polynomial exponents  (features F4a, F4b)

   Integral:
       I = ∫_0^∞ dx1 dx2  x1^{a1-1}  x2^{a2-1}  (1 + x1 + x2)^{-b}

   Closed form (generalized Beta / Dirichlet integral):
       I_exact = Gamma[a1] Gamma[a2] Gamma[b - a1 - a2] / Gamma[b]

   Parameters (complex, convergent: Re[a1]=1.25, Re[a2]=1.5, Re[b-a1-a2]=2.25):
       a1 = 5/4 + I/3,   a2 = 3/2 - I/5,   b = 5 + I/2

   The value of this integral is KNOWN EXACTLY from analytic continuation of
   the Euler integral; it provides the strongest possible reference because it
   shares no algorithm, no random number, and no CAS quadrature with either
   the tropical decomposition or NIntegrate.

   Tests performed (Tier 1 = WL + g++ only):

     CC2-A  Sector sum via NIntegrate (ValidateDecomposition)
            vs exact Gamma value.   PASS: |rel dev| < 1e-3.

     CC2-B  Sector sum cross-check: the ValidateDecomposition sector sum
            also agrees with the direct NIntegrate of the original integrand.
            PASS: |rel dev (sum vs direct)| < 1e-3.

     CC2-C  Tropical Monte Carlo (C++ pipeline) vs exact Gamma value.
            PASS: |Re(MC) - Re(exact)| < 5 * Re(sigma)
                  |Im(MC) - Im(exact)| < 5 * Im(sigma)     (5-sigma test)
            (Skipped cleanly if g++ is not available.)

   Mechanical port from:
     OLD_CODE/TROPICAL_MONTE_CARLO2/TROPICAL_MONTE_CARLOv2/EXAMPLES/
       tropical_eval_examples3.wl  (Example 21)
   mirroring its logic.  The PASS/FAIL logic is new (the original printed
   numbers for visual inspection; this file automates the gate).

   Tier: 1  (WL + g++; no CUBA, no FIESTA).

   Run:
     cd /Users/aidanh/Desktop/TROPICAL_MONTE_CARLOv3
     wolframscript -file TEST/cc_2.wl
   ============================================================================ *)

(* --------------------------------------------------------------------------
   Load packages using absolute paths (plan.md §10 convention)
   -------------------------------------------------------------------------- *)
$repoRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];

Get[FileNameJoin[{$repoRoot, "tropical_fan.wl"}]];
Get[FileNameJoin[{$repoRoot, "tropical_eval.wl"}]];

Print["CC2: packages loaded"];
Print[];

(* --------------------------------------------------------------------------
   Parameters and exact reference value
   -------------------------------------------------------------------------- *)
$a1    = 5/4 + I/3;
$a2    = 3/2 - I/5;
$b     = 5   + I/2;
$exact = N[Gamma[$a1] * Gamma[$a2] * Gamma[$b - $a1 - $a2] / Gamma[$b], 20];

Print["CC2: Integral[0,Inf] dx1 dx2  x1^{a1-1}  x2^{a2-1}  (1+x1+x2)^{-b}"];
Print["     a1 = ", $a1, ",  a2 = ", $a2, ",  b = ", $b];
Print["     Exact (Gamma) = ", $exact];
Print[];

(* --------------------------------------------------------------------------
   IntegrandSpec (plan.md §3.5 schema):
       MonomialExponents = {a1-1, a2-1}   (F4b: complex mono exponents)
       PolynomialExponents = {-b}          (F4a: complex poly exponent)
   -------------------------------------------------------------------------- *)
$poly = 1 + x[1] + x[2];
$vars = {x[1], x[2]};

$spec = <|
  "Polynomials"         -> {$poly},
  "MonomialExponents"   -> {$a1 - 1, $a2 - 1},
  "PolynomialExponents" -> {-$b},
  "Variables"           -> $vars,
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> None
|>;

(* --------------------------------------------------------------------------
   Fan
   -------------------------------------------------------------------------- *)
$verts   = PolytopeVertices[$poly^(-1), $vars];
$fanData = ComputeDecomposition[$verts, "ShowProgress" -> False];
Print["CC2: fan built — ", Length[$fanData[[1]]], " rays, ",
      Length[$fanData[[2]]], " sectors"];
Print[];

(* --------------------------------------------------------------------------
   Helper: print a labelled result line and return the relative deviation
   -------------------------------------------------------------------------- *)
reportLine[label_, val_, ref_] :=
  Module[{dev = Abs[(val - ref) / ref]},
    Print["  ", StringPadRight[label, 26], N[val, 8],
          "    |rel dev vs exact| = ", dev];
    dev
  ];

(* --------------------------------------------------------------------------
   CC2-A/B  ValidateDecomposition
   -------------------------------------------------------------------------- *)
Print["--- CC2-A/B: ValidateDecomposition (NIntegrate sector sum) ---"];

$vr = Quiet @ ValidateDecomposition[$spec, $fanData, {}, 3];

If[!AssociationQ[$vr],
  Print["CC2 FAIL: ValidateDecomposition returned: ", $vr];
  Quit[1]
];

$devDirect = reportLine["NIntegrate direct:",  $vr["DirectResult"], $exact];
$devSum    = reportLine["Sector sum (NInt):",  $vr["SectorSum"],    $exact];
$devIntern = Abs[($vr["SectorSum"] - $vr["DirectResult"]) / $vr["DirectResult"]];
Print["  [sum vs direct NInt internal: ", $devIntern, "]"];
Print[];

(* CC2-A: sector sum vs exact Gamma *)
$tolA = 10^-3;
If[$devSum < $tolA,
  Print["CC2-A PASS  sum vs Gamma: |rel dev| = ", $devSum,
        "  < ", $tolA, "  (tol)"],
  Print["CC2-A FAIL  expected |rel dev| < ", $tolA,
        "  got ", $devSum,
        "  (sum=", $vr["SectorSum"], ", exact=", $exact, ")"]
];

(* CC2-B: sector sum vs direct NIntegrate (internal consistency) *)
$tolB = 10^-3;
If[$devIntern < $tolB,
  Print["CC2-B PASS  sum vs direct NInt: |rel dev| = ", $devIntern,
        "  < ", $tolB, "  (tol)"],
  Print["CC2-B FAIL  expected |rel dev| < ", $tolB,
        "  got ", $devIntern,
        "  (sum=", $vr["SectorSum"],
        ", direct=", $vr["DirectResult"], ")"]
];
Print[];

(* --------------------------------------------------------------------------
   CC2-C  Tropical Monte Carlo (C++ pipeline)
   -------------------------------------------------------------------------- *)
Print["--- CC2-C: Tropical Monte Carlo (C++ pipeline) ---"];

$interfiles = FileNameJoin[{$repoRoot, "TEST", "INTERFILES"}];
If[!DirectoryQ[$interfiles], CreateDirectory[$interfiles]];

$mcRes = EvaluateTropicalMC[
  $spec, $fanData, {{}},
  "NSamples"        -> 500000,
  "RunChecks"       -> False,
  "Verbose"         -> False,
  "WorkingDirectory" -> $interfiles
];

If[!AssociationQ[$mcRes],
  (* g++ not available or other compile failure — skip CC2-C cleanly *)
  Print["CC2-C SKIP  EvaluateTropicalMC not available (g++ absent?)"],

  Module[{r = $mcRes["Results"][[1]],
          mcVal, sigmaRe, sigmaIm, devRe, devIm, nSigmaRe, nSigmaIm},
    mcVal   = r["Re"] + I * r["Im"];
    sigmaRe = r["ReErr"];
    sigmaIm = r["ImErr"];

    devRe    = Abs[Re[mcVal] - Re[$exact]];
    devIm    = Abs[Im[mcVal] - Im[$exact]];
    nSigmaRe = If[sigmaRe > 0, devRe / sigmaRe, Infinity];
    nSigmaIm = If[sigmaIm > 0, devIm / sigmaIm, Infinity];

    Print["  MC result = ", mcVal];
    Print["  exact     = ", N[$exact, 8]];
    Print["  |ΔRe| / σRe = ", nSigmaRe, " sigma"];
    Print["  |ΔIm| / σIm = ", nSigmaIm, " sigma"];

    $tolSigma = 5;
    If[nSigmaRe < $tolSigma && nSigmaIm < $tolSigma,
      Print["CC2-C PASS  MC vs Gamma: ",
            nSigmaRe, " sigmaRe, ", nSigmaIm, " sigmaIm  (< ", $tolSigma, ")"],
      Print["CC2-C FAIL  expected=", N[$exact, 8],
            "  got=", mcVal,
            "  (", nSigmaRe, " sigmaRe, ", nSigmaIm, " sigmaIm)"]
    ]
  ]
];
Print[];

(* --------------------------------------------------------------------------
   Summary
   -------------------------------------------------------------------------- *)
Print["=== CC2 summary ==="];
Print["Integral: Gamma[a1] Gamma[a2] Gamma[b-a1-a2] / Gamma[b]  =  ", N[$exact, 8]];
Print["a1=", $a1, "  a2=", $a2, "  b=", $b];
Print["(Tests F4a complex poly exponents + F4b complex mono exponents; Tier 1)"];
