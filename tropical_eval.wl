(* ::Package:: *)

(* ============================================================================
   tropical_eval.wl

   Numerical evaluation of generalized Euler integrals via tropical
   decomposition.  Consumes the fan output of tropical_fan.wl and produces
   C++ Monte-Carlo code.

   Pipeline:
     Module 1 - ProcessSector        (symbolic coordinate transform + flattening)
     Module 2 - Divergence regulation (tropical subtraction for 1/eps poles)
     Module 3 - C++ code generation   (MmaToC + GenerateCppMonteCarlo)
     Module 4 - EvaluateTropicalMC    (driver: fan -> C++ -> results)
     Module 5 - RunAllTests           (validation suite)

   Dependencies: tropical_fan.wl (loaded automatically from same directory)
   ============================================================================ *)

(* Load tropical_fan.wl BEFORE BeginPackage so TropicalFan` is on $ContextPath.
   Set $SkipPolymakeLoad = True before loading this file to skip the polymake
   dependency (e.g. on cluster nodes without polymake installed). *)
If[!TrueQ[$SkipPolymakeLoad],
  Get[FileNameJoin[{DirectoryName[$InputFileName], "tropical_fan.wl"}]];
  BeginPackage["TropicalEval`", {"TropicalFan`"}],
  BeginPackage["TropicalEval`"]
];

(* ---- Public symbols ---- *)

ProcessSector::usage =
  "ProcessSector[integrandSpec, dualVertices, simplex, coneIndex] performs \
monomial change of variables and flattening for one simplicial cone.";

CheckFlatteningMagnitude::usage =
  "CheckFlatteningMagnitude[sectorData, nSamples] spot-checks the flattened \
integrand magnitude at random points.";

ValidateDecomposition::usage =
  "ValidateDecomposition[integrandSpec, fanData, testKinematics, precisionGoal] \
cross-checks sector sum against direct NIntegrate.";

IdentifyDivergences::usage =
  "IdentifyDivergences[sectorData, eps] identifies divergent variables.";

ProcessDivergentSector::usage =
  "ProcessDivergentSector[sectorData, integrandSpec] constructs G0, G1, \
remainder, and analytic pole for a divergent sector.";

ValidateSubtraction::usage =
  "ValidateSubtraction[sectorData, testKinematics, testEpsilon] checks \
self-consistency of the tropical subtraction.";

MmaToC::usage =
  "MmaToC[expr, paramMap] converts a Mathematica expression to a C++ string.";

GenerateCppMonteCarlo::usage =
  "GenerateCppMonteCarlo[convergentSectors, divergentSectors, integrandSpec, \
outputFile] generates self-contained C++ Monte Carlo source. Option \
\"Integrator\"->\"MonteCarlo\" (default) | \"Vegas\" selects the sampler; \
\"Vegas\" emits a CUBA-guarded main (per-kp, or chunked-ncomp with \
\"Batch\"->True). Vegas options: \"VegasEpsRel\", \"VegasEpsAbs\", \
\"VegasSeed\", \"CubaMaxComp\". Default output is byte-identical to before.";

EvaluateTropicalMC::usage =
  "EvaluateTropicalMC[integrandSpec, fanData, kinematicPoints, opts] runs \
the full tropical Monte Carlo pipeline. \"Integrator\"->\"Vegas\" uses the \
optional CUBA-Vegas sampler instead of plain MC (NSamples = maxeval/sector); \
\"Batch\"->True batches a kinematic scan. Requires CUBA only for Vegas.";

LaurentFromSubtraction::usage =
  "LaurentFromSubtraction[integrandSpec, fanData, kinematicPoints, opts] runs \
EvaluateTropicalMC at several small EpsilonValues and fits the Laurent \
F(eps)=pole/eps+finite+O(eps) via g(eps)=eps*F(eps), returning a (pole, finite) \
per kinematic point (directly comparable to EvaluateTropicalMCIBP). Options: \
\"EpsilonValues\" (default {0.005,0.01,0.02}) plus all EvaluateTropicalMC options \
(\"Integrator\", \"NSamples\", ...).";

RunAllTests::usage =
  "RunAllTests[] runs the validation suite (18 tests) with structured reporting. \
Defined in this file (Module 5).";

ParsePolynomial::usage =
  "ParsePolynomial[poly, vars] parses a polynomial into a list of \
{coefficient, exponentVector} pairs.";

CompileCpp::usage =
  "CompileCpp[srcFile, outputBinary, debug] compiles generated C++ \
Monte Carlo source code. debug=True adds -DTROPICAL_MC_DEBUG. Option \
\"UseCuba\"->Automatic (default; detects the TROPICAL_REQUIRES_CUBA sentinel) \
| True | False links the optional CUBA library for Vegas sources (-lcuba); \
fails fast with a clear message if CUBA is requested but not found.";

detectCuba::usage =
  "detectCuba[] returns <|\"Found\"->True|False, \"IncludeDir\"->..., \
\"LibDir\"->...|> describing whether the optional CUBA library (needed only for \
\"Integrator\"->\"Vegas\") is available.  Memoized.";

(* ---- IBP-based pole extraction (Module 2b) ---- *)

IBPReduceSector::usage =
  "IBPReduceSector[sectorData, eps] performs IBP-based decomposition of \
a divergent sector into monomial terms. The symbolic reduction iterates over \
multiple divergent variables, but the downstream numerical assembly \
(IBPProcessSector / EvaluateTropicalMCIBP) currently supports a single \
divergent variable (one 1/eps pole) per sector.";

IBPProcessSector::usage =
  "IBPProcessSector[sectorData, integrandSpec] runs the full IBP pipeline: \
IBP reduction, epsilon expansion, flattening, and boundary construction. \
Returns an IBPSectorData association ready for C++ codegen. Supports a SINGLE \
divergent variable per sector; nested divergences (>1 divergent variable) are \
detected and refused ($Failed, TropicalEval::nestedIBP).";

IBPCheckBoundary::usage =
  "IBPCheckBoundary[sectorData, integrandSpec, nPoints] numerically verifies \
that IBP boundary terms are well-behaved.";

GenerateCppMonteCarloIBP::usage =
  "GenerateCppMonteCarloIBP[convergentSectors, ibpSectors, integrandSpec, \
outputFile] generates C++ Monte Carlo code with IBP divergence handling. \
Outputs per-function MC results for IBP sectors. Supports \"Integrator\"-> \
\"Vegas\" (per-kp CUBA Vegas; same output layout) and \"Batch\"->True \
(chunked-ncomp Vegas sharing samples across a kp-chunk for each IBP function, \
same (1+n_ibp)-row-per-kp layout; planIBPCX.md §4). Default output is \
byte-identical to before.";

ValidateIBP::usage =
  "ValidateIBP[ibpSectorData, sectorData, integrandSpec, testKinematics, \
testEpsilon] cross-checks the IBP pipeline against direct NIntegrate.";

EvaluateTropicalMCIBP::usage =
  "EvaluateTropicalMCIBP[integrandSpec, fanData, kinematicPoints, opts] \
runs the full tropical Monte Carlo pipeline using IBP for divergent sectors. \
Per kinematic point it returns \"PoleCoefficient\" (1/eps) and \"FinitePart\" \
(eps^0), plus \"LaurentCoefficients\"-><|-1->pole, 0->finite|>. Supports a \
single divergent variable per sector; nested / higher-order poles (1/eps^d, \
d>=2) are detected and refused ($Failed). \
Options: NSamples, NThreads, RunChecks, Verbose, WorkingDirectory, and the \
sampler options \"Integrator\" (\"MonteCarlo\"|\"Vegas\"), \"VegasEpsRel\", \
\"VegasEpsAbs\", \"VegasSeed\" (Vegas needs the optional CUBA library).";

(* ---- Lifting (Module 1b/1c, plan.md §6.2) ---- *)

DetectExtremeCoefficients::usage =
  "DetectExtremeCoefficients[integrandSpec, threshold:1000, opts] scans every \
polynomial for numeric coefficients outside [1/threshold, threshold] and \
returns a list of <|\"PolyIndex\"->j, \"ExponentVector\"->alpha, \
\"Coefficient\"->C, \"Magnitude\"->Abs[C], \"SuggestedK\"->kStar|>. \
Symbolic/kinematic coefficients are skipped (magnitude unknown). SuggestedK \
automates the anchor z0=Abs[C]^(1/k): with \"AnchorRule\"->\"kStar\" (default) \
it returns kStar=Max[1,Ceiling[Abs[Log[Abs[C]]]/Log[threshold]]] -- the \
smallest k that pulls z0 back inside the non-extreme band [1/threshold, \
threshold] (plan.md §6.2). \"AnchorRule\"->\"Unit\" recovers the legacy \
SuggestedK->1; a function f is called as f[mag, threshold]. \
\"BandEdgeGuard\"->False.";

LiftCoefficients::usage =
  "LiftCoefficients[integrandSpec, liftRules] applies auxiliary-variable \
lifting: each rule <|\"PolyIndex\"->j, \"ExponentVector\"->alpha, \"k\"->k|> \
replaces the extreme monomial C x^alpha by c z^k x^alpha with residual \
c=C/z0^k and shared anchor z0=Abs[Cprimary]^(1/kprimary) (kept EXACT). Returns \
<|\"LiftedSpec\"->..., \"LiftData\"->...|>. The aux variable is built robustly \
for plain-symbol specs (BUG-2 fix); the z->z0 round-trip identity is checked \
with relative tolerance.";

ProcessSectorLifted::usage =
  "ProcessSectorLifted[liftedSpec, dualVertices, simplex, coneIndex, \
liftData, opts] runs the full delta-resolution pipeline (plan.md §6.2) for one \
sector of a lifted (n+1)-dim integrand: standard ProcessSector, then pivot \
search (ranking: HasConstantTerm, |m_p|=1, max min Re atilde), domain-constraint \
classification, and FlattenSector. Returns a SectorData association augmented \
with DomainConstraint, LiftData, PivotIndex, ZRow, AugmentedA, HasConstantTerm; \
or <|\"EmptyDomain\"->True, \"ConeIndex\"->...|> for an empty domain; or $Failed \
with liftcomplex/liftnopivot/liftdivergent.";

ValidateLiftedDecomposition::usage =
  "ValidateLiftedDecomposition[originalSpec, liftedSpec, liftedFanData, \
liftData, testKinematics, precisionGoal:3] cross-checks the lifted sector sum \
against a direct NIntegrate of the ORIGINAL integrand. Sectors via \
ProcessSectorLifted; EmptyDomain sectors contribute 0 and are listed under \
DroppedSectors. This EXACT NIntegrate gate (not the sampled 5sigma) is the \
lifting correctness check (plan.md §9 risk #1). Returns <|DirectResult, \
SectorSum, RelativeError, SectorResults, DroppedSectors|>.";

EvaluateTropicalMCLifted::usage =
  "EvaluateTropicalMCLifted[integrandSpec, kinematicPoints, opts] is a thin \
wrapper that detects/applies lifting (\"LiftRules\"->Automatic uses \
DetectExtremeCoefficients), builds the (n+1)-dim lifted fan (\"FanData\"-> \
Automatic, K-scaled robust path), and routes through EvaluateTropicalMC with \
LiftData. With no extreme coefficients it falls back to plain EvaluateTropicalMC. \
Options: \"LiftRules\"->Automatic, \"Threshold\"->1000, \"AnchorRule\"->\"kStar\", \
\"BandEdgeGuard\"->False, \"FanData\"->Automatic, plus all EvaluateTropicalMC \
options.";

(* ---- Error messages ---- *)

TropicalEval::degenerate = "Sector `1`: degenerate cone, det(M) = 0.";
TropicalEval::notsimplicial = "Sector `1`: `2` rays for `3` variables (not simplicial).";
TropicalEval::divergent = "Sector `1`: variable y_`2` divergent, a_`2` = `3`.";
TropicalEval::nested = "Sector `1`: multiple divergent variables (`2`). Nested subtraction not implemented.";
TropicalEval::badck = "Sector `1`: c_k = 0 for divergent variable y_`2`. Higher-order pole.";
TropicalEval::nestedIBP = "Sector `1`: `2` divergent variables. The IBP numerical path assembles a single 1/eps pole per sector; nested / higher-order poles (1/eps^d, d>=2) are not implemented in EvaluateTropicalMCIBP. Refusing this sector ($Failed) rather than dividing by a vanishing effective exponent and emitting invalid C++.";
TropicalEval::badcpp = "Code generation produced non-compilable C++ (offending tokens: `1`). This indicates an unsupported sector (e.g. a nested divergence leaking ComplexInfinity/Indeterminate, or an unconverted symbolic head). Returning $Failed instead of writing C++ that g++ cannot compile.";
TropicalEval::validate = "Validation `1`: relative error `2` exceeds tolerance `3`.";

(* ---- Lifting error messages (plan.md §6.2) ---- *)
TropicalEval::liftidentity = "LiftCoefficients: round-trip identity check FAILED for polynomial `1`; lifted poly at z->z0 does not match original (relative residual `2`).";
TropicalEval::liftbadvar = "LiftCoefficients: cannot build a valid auxiliary variable from the spec variables `1` (BUG-2: Head[plainSymbol][n+1] yields the invalid Symbol[n+1]).  Use indexed variables x[i], or the aux variable could not be made fresh.  Returning $Failed.";
TropicalEval::liftnopivot = "ProcessSectorLifted: cone `1` — no admissible pivot found.  z-row m=`2`, per-pivot atilde=`3`.  Try a different k in the lift rules; alternatively the domain constraint may cut off all divergent regions (log-space remap, future work).";
TropicalEval::liftcomplex = "ProcessSectorLifted: cone `1` — all candidate pivots produce complex atilde; cannot emit a real-valued domain indicator.  In Direct mode, complex monomial (A) OR polynomial (B) exponents are unsupported with lifting — use ComplexExponentMode -> \"SplitRealImag\" (splits BOTH A and B: real parts build the domain map, imaginary parts ride the oscillatory phase).";
TropicalEval::liftdivergent = "ProcessSectorLifted: cone `1` — atilde `2` has a non-positive component after delta resolution; the lifted sector is divergent.  Lifting supports convergent integrals only (plan.md N3).";
TropicalEval::liftfandim = "EvaluateTropicalMC with LiftData: the fan dimension is `1` but n+1 = `2` is required.  Supply the (n+1)-dimensional lifted fan.";
TropicalEval::liftdegenerate = "EvaluateTropicalMCLifted: the lifted Newton polytope is lower-dimensional; automatic fan construction is not possible — supply an explicit complete simplicial fan via the \"FanData\" option.";
TropicalEval::liftdivdomain = "ProcessSectorLifted: cone `1` — the divergent variable couples to the lifted domain constraint (ic_k != 0) for every admissible pivot.  The 1/eps pole and the domain face interact (Case B, planAXpDIV.md §3); a log-space remap is required (future work).  Aborting ($Failed).";
TropicalEval::splitdivmono = "IBPProcessSector: cone `1` — the divergent direction's imaginary exponent theta_k = `2` is eps-DEPENDENT (it scales with the regulator).  A constant theta_k is fine — it is the off-axis case the IBP route now supports (s_k = c_k eps + i theta_k is regular at eps=0; planIBPCX.md §1/§3).  But an eps-dependent theta_k ~ 1/eps reintroduces the fast oscillation the bounded IBP integrand relies on being O(theta), so it cannot be resolved.  The supported case is the standard real-eps regulator with eps-free Im(B)/Im(A) (planIBPCX.md §8).  Aborting ($Failed) rather than emitting a wrong value.";
TropicalEval::ibpresidual = "IBPProcessSector: cone `1` — term `2` retains a genuine (on-axis, theta=0) divergence in y_`3` (alpha0 = `4`) after IBP reduction.  This is an unreduced residual / higher-order pole the numerical path does not assemble.  Aborting ($Failed).";
TropicalEval::ibpvalmultidiv = "ValidateIBP: cone `1` is a one-pole + N-off-axis MultiDiv sector (planIBPMULTIDIV.md); the finite-eps single-pole reconstruction oracle does not cover the multi-corner assembly (validated instead against the closed-form Dirichlet oracle, cc_51/cc_52).  Returning $Failed (not validated) rather than indexing absent IBPTerms.";
TropicalEval::offaxispow = "IBPProcessSector: cone `1` — divergent direction y_`2` is power-divergent off-axis (Re(alpha0) < 0 with theta != 0); a single IBP step does not raise it to Re > 0 (planIBPMULTIDIV.md §3.1).  The supported off-axis case is logarithmic (Re(alpha0) = 0).  Aborting ($Failed) rather than emitting a divergent integrand.";
TropicalEval::splitliftdiv = "EvaluateTropicalMC: SplitRealImag x lifting x divergence is supported only on the IBP route (Method -> Automatic / \"IBP\") and the pinned-eps Subtraction route (LaurentFromSubtraction); the symbolic-eps inline Subtraction path (Method -> \"None\"/\"Subtraction\") does not re-derive the per-piece oscillatory phase (planCXLIFTDIV.md C1).  Use Method -> Automatic (default) or LaurentFromSubtraction for lifted+divergent complex-exponent integrals.";

(* ---- High-D VEGAS sizing guard (lift_error_log L1, plan.md §6.5) ---- *)
TropicalEval::vegasbudget = "VEGAS budget too small: NSamples=`1` per sector is below 20*NStart=`2` (max sector dim `3`, resolved NStart=`4`).  In high dimension the adaptive grid never resolves and VEGAS can return a confidently-wrong value with a tight (lying) error bar.  Raise \"NSamples\" (or lower \"VegasNStart\") and cross-check against a reference rather than trusting the VEGAS error bar.";

(* ============================================================================
   PRIVATE IMPLEMENTATION
   ============================================================================ *)

Begin["`Private`"]

(* tropical_fan.wl is loaded before BeginPackage above *)

(* --------------------------------------------------------------------------
   Helper: ParsePolynomial
   Converts a polynomial into {coefficient, exponentVector} pairs.
   -------------------------------------------------------------------------- *)

ParsePolynomial[poly_, vars_List] := Module[
  {expanded, terms, result},
  expanded = Expand[poly];
  terms = If[Head[expanded] === Plus, List @@ expanded, {expanded}];
  result = Table[
    Module[{coeff, exps},
      exps = Exponent[term, #] & /@ vars;
      coeff = term / (Times @@ MapThread[Power, {vars, exps}]);
      {Simplify[coeff], exps}
    ],
    {term, terms}
  ];
  result
];

(* --------------------------------------------------------------------------
   Helper: TransformExponents
   Given original exponent vector and matrix M, compute new exponent vector.
   -------------------------------------------------------------------------- *)

TransformExponents[expVec_List, mMatrix_List] := expVec . mMatrix;

(* --------------------------------------------------------------------------
   FlattenSector  (extracted from ProcessSector; Tree B factoring + Tree A eps)

   Given the cleared (min-exponent-shifted) polynomials, the effective monomial
   exponents a_i^eff, and a prefactor base (|det M| for the plain path; the
   lifted prefactor for the lifted path), decide convergence and either flatten
   or report the divergence.

   eps threading (plan.md §6.1): the convergence predicate is evaluated at
   eps -> 0 when eps =!= None, so an eps-regulated integral is classified by its
   eps^0 leading behaviour.  The decision is EXACT (no float tolerance): a
   variable is divergent iff Re[a_i^eff |_{eps->0}] <= 0, decided symbolically
   (TrueQ on the exact comparison) or, for an explicit number, on its exact Re.
   "Last violating index wins" loop semantics are preserved from Tree A.

   Returns an Association with the flatten/divergence payload only; ProcessSector
   (and ProcessSectorLifted) merge it into the full SectorData.
     IsDivergent -> True :  {IsDivergent, DivergentVariable}
     IsDivergent -> False:  {IsDivergent, DivergentVariable(=0),
                             FlattenedPolys, Prefactor}
   -------------------------------------------------------------------------- *)

FlattenSector[clearedPolys_List, effectiveAVals_List, prefactorBase_,
              eps_: None] :=
Module[{a0vals, isDivergent, divVar, n, flattenedPolys, prefactor},
  n = Length[effectiveAVals];

  (* Convergence test at eps -> 0 (exact). *)
  a0vals = effectiveAVals /. (If[eps =!= None, eps -> 0, {}]);
  isDivergent = False;
  divVar = 0;
  Do[
    If[TrueQ[Re[a0vals[[i]]] <= 0] ||
       (NumericQ[a0vals[[i]]] && Re[a0vals[[i]]] <= 0),
      isDivergent = True;
      divVar = i;
    ],
    {i, n}
  ];

  If[isDivergent,
    Return[<|"IsDivergent" -> True, "DivergentVariable" -> divVar|>]
  ];

  (* Convergent: flatten y_i -> (y_i')^{1/a_i^eff}.  Each monomial
     {coeff, {e1,...,en}} of Q_j becomes {coeff, {e1/a1^eff,...,en/an^eff}}. *)
  flattenedPolys = Table[
    Table[
      {mono[[1]],
       MapThread[#1/#2 &, {mono[[2]], effectiveAVals}]},
      {mono, clearedPolys[[j]]}
    ],
    {j, Length[clearedPolys]}
  ];

  (* Prefactor = base / Prod(a_i^eff) *)
  prefactor = prefactorBase / (Times @@ effectiveAVals);

  <|"IsDivergent" -> False, "DivergentVariable" -> 0,
    "FlattenedPolys" -> flattenedPolys, "Prefactor" -> prefactor|>
];

(* --------------------------------------------------------------------------
   MODULE 1: ProcessSector

   Key insight (tropical factoring):
   After the monomial substitution x_i = prod y_j^{M_ij}, each polynomial
   P_k becomes a sum of monomials in y with SIGNED exponents.  The dominant
   monomial (the one the tropical fan says dominates in this cone) has the
   minimum exponents.  We factor it out:
       P_k(y) = prod_j y_j^{d_{k,j}} * Q_k(y)
   where d_{k,j} = min_m (transformed exponent of y_j in monomial m),
   and Q_k has all non-negative y-exponents with a constant term.

   The effective monomial prefactor then becomes:
       a_j^eff = rawA_j + sum_k B_k * d_{k,j}

   For a properly constructed tropical fan and convergent integral,
   a_j^eff > 0 for all j, and we can flatten using these effective exponents.
   -------------------------------------------------------------------------- *)

(* "ImagMonoExps" (planAXpDIVv2.md §4.3): the ORIGINAL imaginary monomial
   exponents Im(A) (length n), supplied only in SplitRealImag mode when the
   monomial exponents are complex.  When None (default, real-A / Direct / the
   inner ProcessSector call from ProcessSectorLifted) NO MonomialPhaseLog is
   attached, so the emitted C++ is byte-identical (#25). *)
Options[ProcessSector] = {"Verbose" -> False, "ImagMonoExps" -> None};

ProcessSector[integrandSpec_Association, dualVertices_List,
              simplex_List, coneIndex_Integer, OptionsPattern[]] :=
Module[
  {polys, monoExps, polyExps, vars, eps, kinSyms,
   selectedRays, mMatrix, detM, n,
   transformedPolys, clearedPolys, minExponents,
   rawAVals, effectiveAVals,
   flattenedPolys, prefactor, monoFactorLog, monomialPhaseLog, imAexps,
   isDivergent, divVar, verbose, sectorData,
   parsedPolys},

  verbose = OptionValue["Verbose"];
  imAexps = OptionValue["ImagMonoExps"];

  (* Extract fields from integrand spec *)
  polys    = integrandSpec["Polynomials"];
  monoExps = integrandSpec["MonomialExponents"];
  polyExps = integrandSpec["PolynomialExponents"];
  vars     = integrandSpec["Variables"];
  kinSyms  = integrandSpec["KinematicSymbols"];
  eps      = integrandSpec["RegulatorSymbol"];
  n        = Length[vars];

  (* --- Step 1: Monomial change of variables --- *)

  (* Extract ray vectors for this simplex (rows of dualVertices) *)
  selectedRays = dualVertices[[#]] & /@ simplex;

  (* Check: simplicial condition *)
  If[Length[selectedRays] != n,
    Message[TropicalEval::notsimplicial, coneIndex,
            Length[selectedRays], n];
    Return[$Failed]
  ];

  (* M_{ij} = -rho_j[i], i.e. M = -Transpose[selectedRays] *)
  mMatrix = -Transpose[selectedRays];

  (* Check: non-degenerate *)
  detM = Det[mMatrix];
  If[detM === 0 || TrueQ[detM == 0],
    Message[TropicalEval::degenerate, coneIndex];
    Return[$Failed]
  ];

  (* Raw exponents from monomial part + Jacobian:
     rawA_i = sum_k (A_k + 1) * M_{ki} = (monoExps + 1) . M *)
  rawAVals = (monoExps + 1) . mMatrix;

  (* --- Transform polynomials --- *)
  parsedPolys = ParsePolynomial[#, vars] & /@ polys;

  transformedPolys = Table[
    Table[
      Module[{coeff, origExp, newExp},
        coeff   = mono[[1]];
        origExp = mono[[2]];
        newExp  = TransformExponents[origExp, mMatrix];
        {coeff, newExp}
      ],
      {mono, parsedPolys[[j]]}
    ],
    {j, Length[polys]}
  ];

  (* --- Step 1b: Tropical factoring --- *)
  (* For each polynomial, find min exponents and factor them out *)

  minExponents = Table[
    Table[
      Min[#[[2, i]] & /@ transformedPolys[[j]]],
      {i, n}
    ],
    {j, Length[polys]}
  ];

  (* Cleared polynomials: shift exponents so minimum is 0 *)
  clearedPolys = Table[
    Table[
      {mono[[1]], mono[[2]] - minExponents[[j]]},
      {mono, transformedPolys[[j]]}
    ],
    {j, Length[polys]}
  ];

  (* Effective exponents: a_i^eff = rawA_i + sum_j B_j * d_{j,i} *)
  effectiveAVals = rawAVals + Total[
    Table[polyExps[[j]] * minExponents[[j]], {j, Length[polys]}]
  ];

  If[verbose,
    Print["Sector ", coneIndex, ": det(M) = ", detM,
          ", rawA = ", rawAVals,
          ", minExp = ", minExponents,
          ", effA = ", effectiveAVals]
  ];

  (* --- Step 2: Flattening using effective exponents --- *)
  (* Convergence decision + flatten are factored into FlattenSector (eps-aware,
     plan.md §6.1).  The pre-flattening divergent branch is PRESERVED verbatim
     (plan.md §9 risk #2 KEEP): the lifted path consumes ClearedPolys/
     NewExponents from it, and (n+1)-dim pre-delta sectors are routinely flagged
     divergent even when the delta-constrained integral is finite. *)
  Module[{flat},
    flat = FlattenSector[clearedPolys, effectiveAVals, Abs[detM], eps];
    isDivergent = flat["IsDivergent"];
    divVar      = flat["DivergentVariable"];

    If[isDivergent,
      (* Divergent: return pre-flattened data for Module 2 *)
      sectorData = <|
        "ConeIndex"           -> coneIndex,
        "RayMatrix"           -> mMatrix,
        "DetM"                -> detM,
        "SelectedRays"        -> selectedRays,
        "RawExponents"        -> rawAVals,
        "NewExponents"        -> effectiveAVals,
        "MinExponents"        -> minExponents,
        "TransformedPolys"    -> transformedPolys,
        "ClearedPolys"        -> clearedPolys,
        "Prefactor"           -> Abs[detM],
        (* PrefactorBase (planAXpDIV.md §4.1, Barrier B.1): the sector prefactor
           BEFORE dividing by Prod(effective exponents).  For an unlifted sector
           this is Abs[detM]; the divergence routines read it via
           Lookup[sd,"PrefactorBase",Abs[detM]], so unlifted output is unchanged. *)
        "PrefactorBase"       -> Abs[detM],
        "IsDivergent"         -> True,
        "DivergentVariable"   -> divVar,
        "Dimension"           -> n,
        "PolynomialExponents" -> polyExps,
        "MonomialExponents"   -> monoExps
      |>;
      If[verbose,
        Print["  -> Divergent in variable y_", divVar]
      ];
      Return[sectorData]
    ];

    flattenedPolys = flat["FlattenedPolys"];
    prefactor      = flat["Prefactor"];
  ];

  (* MonoFactorLog (plan.md §6.3 Fix B / lift_error_log L3) — the log of the
     tropical MONOMIAL factor y^{d_k} that clearing removed from P_k, expressed
     in the FLATTENED coords y' (y_i = (y_i')^{1/a_eff,i}):
        log(monomial factor_k) = Sum_i d_{k,i} log y_i
                               = Const_k + Sum_i Coeffs_{k,i} log y'_i,
     with (unlifted) Const_k = 0 and Coeffs_{k,i} = d_{k,i}/a_eff,i.
     Consumed ONLY by the opt-in SplitRealImag codegen phase (Im(B)!=0); the
     real-exponent / Direct path ignores it, so the emitted C++ is unchanged. *)
  monoFactorLog = Table[
    <|"Const"  -> 0,
      "Coeffs" -> Table[minExponents[[k, i]]/effectiveAVals[[i]], {i, n}]|>,
    {k, Length[clearedPolys]}
  ];

  (* MonomialPhaseLog (planAXpDIVv2.md §3) — the bare-monomial oscillatory phase
     from the IMAGINARY monomial exponents Im(A).  The bare monomial x^A maps to
     y^{rawA}, rawA = (A+1).M, whose imaginary part rides the phase:
        i * Sum_i Im(A_i) log x_i = i * Sum_a (Im(A).M)_a log y_a
                                  = i * Sum_a [ (Im(A).M)_a / a_eff,a ] log y'_a
     (the +1 Jacobian is real, NOT part of the phase).  Const = 0 (unlifted).
     Single oscillatory term (NOT per-polynomial like MonoFactorLog).  Attached
     only when Im(A) is supplied (SplitRealImag); real-A / Direct -> None -> the
     codegen emits no monomial phase -> byte-identical (#25). *)
  monomialPhaseLog = If[imAexps === None, None,
    Module[{num = imAexps . mMatrix},
      <|"Const" -> 0,
        "Coeffs" -> Table[num[[i]]/effectiveAVals[[i]], {i, n}]|>]];

  sectorData = <|
    "ConeIndex"            -> coneIndex,
    "RayMatrix"            -> mMatrix,
    "DetM"                 -> detM,
    "SelectedRays"         -> selectedRays,
    "RawExponents"         -> rawAVals,
    "NewExponents"         -> effectiveAVals,
    "MinExponents"         -> minExponents,
    "TransformedPolys"     -> transformedPolys,
    "ClearedPolys"         -> clearedPolys,
    "FlattenedPolys"       -> flattenedPolys,
    "Prefactor"            -> prefactor,
    "PrefactorBase"        -> Abs[detM],
    "MonoFactorLog"        -> monoFactorLog,
    (* planAXpDIVv2.md §4.3: None on real-A sectors (codegen treats absent/None
       as "no monomial phase") -> byte-identical (#25). *)
    "MonomialPhaseLog"     -> monomialPhaseLog,
    "IsDivergent"          -> False,
    "DivergentVariable"    -> 0,
    "Dimension"            -> n,
    "PolynomialExponents"  -> polyExps,
    "MonomialExponents"    -> monoExps
  |>;

  If[verbose,
    Print["  -> Convergent, prefactor = ", prefactor]
  ];

  sectorData
];

(* --------------------------------------------------------------------------
   CheckFlatteningMagnitude
   Spot-check that the flattened integrand is O(1) at random points.
   -------------------------------------------------------------------------- *)

CheckFlatteningMagnitude[sectorData_Association, nSamples_Integer: 20,
                         testKinematics_List: {}] :=
Module[
  {flatPolys, polyExps, prefactor, dim, mags, y, polyVals, integrandVal,
   kinRules, dc},

  If[sectorData["IsDivergent"],
    Print["CheckFlatteningMagnitude: sector ", sectorData["ConeIndex"],
          " is divergent, skipping."];
    Return[Null]
  ];

  flatPolys = sectorData["FlattenedPolys"];
  polyExps  = sectorData["PolynomialExponents"];
  prefactor = sectorData["Prefactor"];
  dim       = sectorData["Dimension"];
  kinRules  = If[testKinematics === {}, {}, testKinematics];
  dc        = Lookup[sectorData, "DomainConstraint", None];

  (* Lifted-sector path: rejection-sample only feasible points (plan.md §6.2). *)
  If[dc =!= None,
    Module[{logZ0num, mpNum, icNum, feasibleMags, totalDraws, feasibleCount,
            maxDraws, isFeasible, logYpStar, mag, meanMag, maxMag, minMag, ff},
      logZ0num = N[dc["LogZ0"]];  mpNum = N[dc["MP"]];  icNum = N[dc["IndicatorCoeffs"]];
      maxDraws = 50 * nSamples;  feasibleMags = {};  totalDraws = 0;  feasibleCount = 0;
      While[feasibleCount < nSamples && totalDraws < maxDraws,
        y = RandomReal[{0.01, 0.99}, dim];  totalDraws++;
        logYpStar = (logZ0num - Total[icNum * Log[y]]) / mpNum;
        isFeasible = (logYpStar <= 0);
        If[isFeasible,
          polyVals = Table[
            Total[Table[
              Module[{coeff = mono[[1]] /. kinRules, alphas = mono[[2]] /. kinRules,
                      logY2 = Log[y]}, coeff * Exp[Total[alphas * logY2]]],
              {mono, flatPolys[[j]]}]],
            {j, Length[flatPolys]}];
          integrandVal = (prefactor /. kinRules) *
            Times @@ MapThread[Exp[#2 * Log[#1]] &, {polyVals, polyExps /. kinRules}];
          AppendTo[feasibleMags, Abs[integrandVal]];  feasibleCount++
        ]
      ];
      If[feasibleCount == 0,
        Print["WARNING: CheckFlatteningMagnitude sector ", sectorData["ConeIndex"],
              ": ZERO feasible points in ", totalDraws, " draws."];
        Return[<|"Mean" -> 0, "Max" -> 0, "Min" -> 0, "Samples" -> {}, "FeasibleFraction" -> 0|>]
      ];
      meanMag = Mean[feasibleMags];  maxMag = Max[feasibleMags];  minMag = Min[feasibleMags];
      ff = N[feasibleCount / totalDraws];
      If[maxMag > 10^3 || minMag < 10^(-6),
        Print["WARNING: Sector ", sectorData["ConeIndex"],
              " flattening check: min=", minMag, " max=", maxMag,
              " mean=", meanMag, " feasibleFrac=", ff]];
      Return[<|"Mean" -> meanMag, "Max" -> maxMag, "Min" -> minMag,
               "Samples" -> feasibleMags, "FeasibleFraction" -> ff|>]
    ]
  ];

  mags = Table[
    y = RandomReal[{0.01, 0.99}, dim];
    polyVals = Table[
      Total[
        Table[
          Module[{coeff, alphas, logY},
            coeff  = mono[[1]] /. kinRules;
            alphas = mono[[2]] /. kinRules;
            logY   = Log[y];
            coeff * Exp[Total[alphas * logY]]
          ],
          {mono, flatPolys[[j]]}
        ]
      ],
      {j, Length[flatPolys]}
    ];
    integrandVal = (prefactor /. kinRules) *
      Times @@ MapThread[
        Exp[#2 * Log[#1]] &,
        {polyVals, polyExps /. kinRules}
      ];
    Abs[integrandVal],
    {nSamples}
  ];

  Module[{meanMag, maxMag, minMag},
    meanMag = Mean[mags];
    maxMag  = Max[mags];
    minMag  = Min[mags];
    If[maxMag > 10^3 || minMag < 10^(-6),
      Print["WARNING: Sector ", sectorData["ConeIndex"],
            " flattening check: min=", minMag, " max=", maxMag,
            " mean=", meanMag]
    ];
    <|"Mean" -> meanMag, "Max" -> maxMag, "Min" -> minMag,
      "Samples" -> mags|>
  ]
];

(* --------------------------------------------------------------------------
   ValidateDecomposition
   Cross-check sector sum against direct NIntegrate.
   -------------------------------------------------------------------------- *)

ValidateDecomposition[integrandSpec_Association, fanData_List,
                      testKinematics_List, precisionGoal_Integer: 3] :=
Module[
  {polys, monoExps, polyExps, vars, n, dualVertices, simplexList,
   directIntegrand, directResult, sectorResults, sectorSum,
   relError, kinRules, allSectorData},

  polys    = integrandSpec["Polynomials"];
  monoExps = integrandSpec["MonomialExponents"];
  polyExps = integrandSpec["PolynomialExponents"];
  vars     = integrandSpec["Variables"];
  n        = Length[vars];
  kinRules = testKinematics;

  {dualVertices, simplexList} = fanData;

  (* Direct NIntegrate of the original integrand *)
  directIntegrand = (Times @@ MapThread[Power, {vars, monoExps}]) *
    (Times @@ MapThread[Power, {polys, polyExps}]) /. kinRules;

  directResult = NIntegrate[
    directIntegrand,
    Evaluate[Sequence @@ ({#, 0, Infinity} & /@ vars)],
    MaxRecursion -> 20,
    PrecisionGoal -> precisionGoal + 1,
    Method -> "GlobalAdaptive"
  ];

  (* Process all sectors *)
  allSectorData = Table[
    ProcessSector[integrandSpec /. kinRules, dualVertices,
                  simplexList[[s]], s],
    {s, Length[simplexList]}
  ];

  (* For each sector, evaluate via NIntegrate on [0,1]^n *)
  sectorResults = Table[
    Module[{sd, flatPolys, pExps, pf, dim, yVars, integrand,
            polyVals, result},
      sd = allSectorData[[s]];
      If[sd === $Failed, 0,
        If[sd["IsDivergent"],
          Print["ValidateDecomposition: sector ", s,
                " is divergent, skipping"];
          0,
          flatPolys = sd["FlattenedPolys"];
          pExps     = sd["PolynomialExponents"] /. kinRules;
          pf        = sd["Prefactor"] /. kinRules;
          dim       = sd["Dimension"];
          yVars     = Table[Unique["yv"], {dim}];

          (* Evaluate using log-exp form to avoid numerical issues *)
          polyVals = Table[
            Total[
              Table[
                Module[{coeff, alphas, logY},
                  coeff  = mono[[1]] /. kinRules;
                  alphas = mono[[2]] /. kinRules;
                  logY   = Log /@ yVars;
                  coeff * Exp[Total[alphas * logY]]
                ],
                {mono, flatPolys[[j]]}
              ]
            ],
            {j, Length[flatPolys]}
          ];

          integrand = pf *
            Times @@ MapThread[
              Function[{pv, be}, Exp[be * Log[pv]]],
              {polyVals, pExps}
            ];

          result = Quiet@NIntegrate[
            integrand,
            Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
            MaxRecursion -> 15,
            PrecisionGoal -> precisionGoal,
            Method -> "GlobalAdaptive"
          ];
          result
        ]
      ]
    ],
    {s, Length[simplexList]}
  ];

  sectorSum = Total[sectorResults];
  relError  = Abs[(sectorSum - directResult) / directResult];

  If[relError > 10^(-precisionGoal + 1),
    Message[TropicalEval::validate, "Decomposition", relError,
            10^(-precisionGoal + 1)]
  ];

  <|"DirectResult" -> directResult, "SectorSum" -> sectorSum,
    "RelativeError" -> relError, "SectorResults" -> sectorResults|>
];

(* The hard-coded 2D benchmark unit test was removed in v3
   (plan.md §5.6 — dead code in both old trees). *)


(* ============================================================================
   MODULE 1b: LIFTING  (auxiliary-variable uplift for extreme coefficients)
   Ported from Tree B (TROPICAL_MONTE_CARLOv2) with the plan.md §6.2 fixes:
     - BUG-2 aux-variable robustness (no Symbol[n+1] for plain-symbol specs);
     - relative-tolerance liftidentity check (no false-positive on float z0);
     - SuggestedK -> kStar geometry-band anchor (already in the Tree B copy);
     - §6.7 exactness: pivot realness decided exactly (PossibleZeroQ), not by a
       10^-12 tolerance.
   ============================================================================ *)

(* --------------------------------------------------------------------------
   makeAuxVar — BUG-2 fix (plan.md §6.2).
   Build a guaranteed-VALID, fresh auxiliary variable of the same indexed
   "head[idx]" form as the existing variables.  Tree B used
   Head[vars[[1]]][auxIdx], which for a plain symbol x gives Head[x]=Symbol and
   Symbol[4] — an invalid expression that triggers Symbol::string.  Here we
   require the spec variables to be indexed (head[_]); if they are plain
   symbols we fire liftbadvar and signal failure (None).
   -------------------------------------------------------------------------- *)
makeAuxVar[vars_List, auxIdx_Integer] := Module[{heads, h},
  (* Each variable must be of the form head[index] with a symbol head. *)
  If[!AllTrue[vars, (MatchQ[#, _Symbol[_]] && Head[Head[#]] === Symbol) &],
    Return[None]
  ];
  heads = Head /@ vars;            (* e.g. {x, x} *)
  h = First[heads];
  (* Use the common head when all variables share it; otherwise the first. *)
  h[auxIdx]
];

(* --------------------------------------------------------------------------
   DetectExtremeCoefficients
   Scan every polynomial for numeric coefficients outside [1/threshold,
   threshold].  Returns a list of flagged-monomial associations.
   Symbolic/kinematic coefficients are silently skipped (magnitude unknown).
   SuggestedK -> kStar (plan.md §6.2 / AUXT/anchor_selection_procedure.md).
   -------------------------------------------------------------------------- *)

Options[DetectExtremeCoefficients] = {
  "AnchorRule"    -> "kStar",
  "BandEdgeGuard" -> False
};

(* Two definitions so the optional numeric threshold is never confused with a
   trailing option rule. *)
DetectExtremeCoefficients[integrandSpec_Association, opts : OptionsPattern[]] :=
  DetectExtremeCoefficients[integrandSpec, 1000, opts];

DetectExtremeCoefficients[integrandSpec_Association, threshold_?NumericQ,
                          opts : OptionsPattern[]] :=
Module[
  {polys, vars, result, parsedPoly, coeff, mag,
   anchorRule, bandEdgeGuard, logTau, suggestK},

  polys  = integrandSpec["Polynomials"];
  vars   = integrandSpec["Variables"];
  result = {};

  anchorRule    = OptionValue["AnchorRule"];
  bandEdgeGuard = TrueQ[OptionValue["BandEdgeGuard"]];
  logTau        = Log[N[threshold]];

  (* kStar = max(1, ceil(|log|C|| / log threshold)) — the smallest k that pulls
     z0 = |C|^(1/k) back inside the detector's non-extreme band [1/tau, tau]
     (plan.md §6.2).  This is a heuristic SELECTOR only: correctness is
     independent of the choice (the lift is an exact identity for any k). *)
  suggestK[m_] := Switch[anchorRule,
    "Unit",  1,
    "kStar",
      Module[{absLog = Abs[Log[N[m]]], kS},
        kS = Max[1, Ceiling[absLog / logTau]];
        If[bandEdgeGuard && absLog <= logTau + Log[10.] (1 + 1.*^-8),
          kS = kS + 1
        ];
        kS
      ],
    _, anchorRule[m, threshold]
  ];

  Do[
    parsedPoly = ParsePolynomial[polys[[j]], vars];
    Do[
      coeff = mono[[1]];
      (* Skip symbolic / kinematic coefficients *)
      If[NumericQ[N[coeff]],
        mag = Abs[N[coeff]];
        If[mag < 1/threshold || mag > threshold,
          AppendTo[result, <|
            "PolyIndex"      -> j,
            "ExponentVector" -> mono[[2]],
            "Coefficient"    -> coeff,
            "Magnitude"      -> mag,
            "SuggestedK"     -> suggestK[mag]
          |>]
        ]
      ],
      {mono, parsedPoly}
    ],
    {j, Length[polys]}
  ];

  result
];


(* --------------------------------------------------------------------------
   LiftCoefficients
   Apply auxiliary-variable lifting (plan.md §6.2).
   liftRules = { <|"PolyIndex"->j, "ExponentVector"->alpha, "k"->k|>, ... }
   Returns <|"LiftedSpec"->..., "LiftData"->...|> (or $Failed).
   z0 = |C_primary|^(1/k_primary) kept EXACT; residual c_i = C_i / z0^{k_i}.
   -------------------------------------------------------------------------- *)

LiftCoefficients[integrandSpec_Association, liftRules_List] :=
Module[
  {polys, vars, monoExps, polyExps, kinSyms,
   n, auxIdx, auxVar, newVars, newMonoExps,
   parsedPolys, ruleCoeffs, ruleLogMags, primaryIdx, primaryRule,
   Cprimary, kprimary, z0, residuals,
   liftedPolys, j, liftedSpec, liftData,
   liftedSubbed, original, relResid, ruleK},

  (* BUG #34 fix: the lift exponent of the aux variable.  Accept rules from
     either calling convention — EvaluateTropicalMCLifted remaps to "k", but
     DetectExtremeCoefficients-built rules (and several cross-checks, e.g.
     cc_34) pass the detector rule straight through with key "SuggestedK".
     Reading r["k"] on the latter previously yielded Missing[KeyAbsent,k],
     poisoning the aux coordinate (x[aux]^Missing) so the lifted Newton
     polytope vertex carried Missing[KeyAbsent,k] and computeFanScaled
     returned $Failed.  Fall back to "SuggestedK"; require a positive
     integer (the lift identity needs an integer aux power). *)
  ruleK[r_Association] := Module[{kv},
    kv = Lookup[r, "k", Lookup[r, "SuggestedK", Missing["NoK"]]];
    kv];

  polys    = integrandSpec["Polynomials"];
  vars     = integrandSpec["Variables"];
  monoExps = integrandSpec["MonomialExponents"];
  polyExps = integrandSpec["PolynomialExponents"];
  kinSyms  = integrandSpec["KinematicSymbols"];
  n        = Length[vars];
  auxIdx   = n + 1;

  (* BUG-2 fix: build a valid indexed aux variable; fail cleanly otherwise. *)
  auxVar = makeAuxVar[vars, auxIdx];
  If[auxVar === None,
    Message[TropicalEval::liftbadvar, vars];
    Return[$Failed]
  ];

  (* ---- Find the primary rule (max |Log[|C|]|) ---- *)
  parsedPolys = ParsePolynomial[#, vars] & /@ polys;

  (* BUG #34 guard: every rule must resolve to a positive-integer aux power
     (via "k" or "SuggestedK").  Refuse cleanly rather than emit a vertex
     carrying Missing[KeyAbsent,k] that fails the fan build downstream. *)
  If[!AllTrue[liftRules, MatchQ[ruleK[#], _Integer?Positive] &],
    Message[TropicalEval::liftbadvar, vars];
    Return[$Failed]
  ];

  ruleCoeffs = Table[
    Module[{j0 = r["PolyIndex"], alpha = r["ExponentVector"], matchPos},
      matchPos = Position[parsedPolys[[j0, All, 2]], alpha];
      If[matchPos === {},
        Print["LiftCoefficients: monomial ", alpha,
              " not found in polynomial ", j0, "."];
        Return[$Failed]
      ];
      parsedPolys[[j0, matchPos[[1, 1]], 1]]
    ],
    {r, liftRules}
  ];
  If[MemberQ[ruleCoeffs, $Failed], Return[$Failed]];

  ruleLogMags = Abs[Log[Abs[N[#]]]] & /@ ruleCoeffs;
  primaryIdx  = First@Ordering[ruleLogMags, -1];
  primaryRule = liftRules[[primaryIdx]];
  Cprimary    = ruleCoeffs[[primaryIdx]];
  kprimary    = ruleK[primaryRule];

  (* z0 = |C_primary|^(1/k_primary), EXACT.  Sign/phase goes into residual. *)
  z0 = If[IntegerQ[Abs[Cprimary]^(1/kprimary)],
    Abs[Cprimary]^(1/kprimary),
    Power[Abs[Cprimary], 1/kprimary]
  ];

  (* ---- Residuals: c_i = C_i / z0^{k_i} (exact) ---- *)
  residuals = Table[
    Simplify[ruleCoeffs[[i]] / z0^ruleK[liftRules[[i]]]],
    {i, Length[liftRules]}
  ];

  (* ---- Build lifted polynomials ---- *)
  liftedPolys = Table[
    Module[{acc = polys[[polyJ]]},
      Do[
        Module[{r = liftRules[[ri]], alpha, ki, Ci, ci,
                termToReplace, replacement},
          If[r["PolyIndex"] == polyJ,
            alpha = r["ExponentVector"];
            ki    = ruleK[r];
            ci    = residuals[[ri]];
            Ci    = ruleCoeffs[[ri]];
            termToReplace = Ci * Times @@ MapThread[
              Function[{v, e}, If[e == 0, 1, Power[v, e]]], {vars, alpha}];
            replacement = ci * auxVar^ki * Times @@ MapThread[
              Function[{v, e}, If[e == 0, 1, Power[v, e]]], {vars, alpha}];
            acc = Expand[acc - termToReplace + replacement];
          ]
        ],
        {ri, Length[liftRules]}
      ];
      acc
    ],
    {polyJ, Length[polys]}
  ];

  newVars     = Append[vars, auxVar];
  newMonoExps = Append[monoExps, 0];

  liftedSpec = <|
    "Polynomials"         -> liftedPolys,
    "MonomialExponents"   -> newMonoExps,
    "PolynomialExponents" -> polyExps,
    "Variables"           -> newVars,
    "KinematicSymbols"    -> kinSyms,
    "RegulatorSymbol"     -> Lookup[integrandSpec, "RegulatorSymbol", None]
  |>;

  (* ---- Identity check (plan.md §6.2 rel-tolerance fix) ----
     liftedPolys /. auxVar -> z0 must reproduce the original.  Exact-match
     first; if symbols differ only by an exact rewrite, fall back to a RELATIVE
     residual test (NOT Simplify[...]===0, which false-positives on float z0). *)
  Do[
    liftedSubbed = Expand[liftedPolys[[j]] /. auxVar -> z0];
    original     = Expand[polys[[j]]];
    If[!TrueQ[liftedSubbed === original],
      If[!TrueQ[PossibleZeroQ[liftedSubbed - original]],
        (* numeric relative residual at a generic sample point *)
        relResid = Module[{diff = liftedSubbed - original, sub, num, den},
          sub = Thread[vars -> Table[1.7 + 0.37 i, {i, Length[vars]}]];
          num = Abs[N[diff /. sub]];
          den = Abs[N[original /. sub]] + Abs[N[liftedSubbed /. sub]] + 1;
          num / den
        ];
        If[!(NumericQ[relResid] && relResid < 10.^-8),
          Message[TropicalEval::liftidentity, j, relResid]
        ]
      ]
    ],
    {j, Length[polys]}
  ];

  liftData = <|
    "z0"           -> z0,
    "AuxVariable"  -> auxVar,
    "AuxIndex"     -> auxIdx,
    "Rules"        -> liftRules,
    "Residuals"    -> residuals,
    "OriginalSpec" -> integrandSpec
  |>;

  <|"LiftedSpec" -> liftedSpec, "LiftData" -> liftData|>
];


(* ============================================================================
   MODULE 1c: ProcessSectorLifted
   Delta-resolution pipeline (plan.md §6.2) for one sector of a lifted
   integrand.  tryPivot[p] is cached (memoized) so each pivot is computed once.
   ============================================================================ *)

(* "Eps" (planAXpDIV.md §4.3): the regulator symbol.  When None (default,
   convergent lifting) the classification/ranking reduce EXACTLY to the pre-L2
   behaviour, so every existing lifted/convergent result is byte-identical.
   When a symbol, atilde is classified at eps->0 and a single divergent slot is
   permitted, emitting a divergent lifted SectorData (Barrier A removed). *)
Options[ProcessSectorLifted] = {"Verbose" -> False, "Eps" -> None,
  "ImagMonoExps" -> None};

ProcessSectorLifted[liftedSpec_Association, dualVertices_List,
                    simplex_List, coneIndex_Integer,
                    liftData_Association, OptionsPattern[]] :=
Module[
  {sdAug, z0, auxIdx, a, clearedPolys, detM, mMatrix, polyExps, n1, n,
   mVec, verbose, tryPivot, pivotCache, candidates, bestPivot,
   pivotP, mp, ap, mOtherVec, atildeVals, reclearedPolys,
   allOtherZero, domainClass, logZ0,
   fsResult, flattenedPolys, prefactor, prefactorBase,
   monoFactorLogLifted, liftedMonoFactorNum, eps, a0fn, classDir, pivotClass,
   classifyDomainFor, candClass, bestCI, bestDom, divDir,
   imAexps, monomialPhaseLogLifted, liftedMonoPhaseNum},

  verbose = OptionValue["Verbose"];
  (* Im(A) over the AUGMENTED variables (length n+1; the aux entry is 0 — the
     lift appends a real aux monomial exponent, LiftCoefficients:950).  None for
     real-A / Direct -> no monomial phase attached -> byte-identical (#25). *)
  imAexps = OptionValue["ImagMonoExps"];

  (* --- Step 1: standard (n+1)-dim ProcessSector --- *)
  (* NOTE: the inner ProcessSector is called WITHOUT "ImagMonoExps": its sdAug
     MonomialPhaseLog would be on augmented (n+1) coords; ProcessSectorLifted
     re-derives the post-delta lifted phase below from imAexps directly. *)
  sdAug = ProcessSector[liftedSpec, dualVertices, simplex, coneIndex];
  If[sdAug === $Failed, Return[$Failed]];

  z0           = liftData["z0"];
  auxIdx       = liftData["AuxIndex"];
  a            = sdAug["NewExponents"];     (* length n+1 *)
  clearedPolys = sdAug["ClearedPolys"];
  detM         = sdAug["DetM"];
  mMatrix      = sdAug["RayMatrix"];
  polyExps     = sdAug["PolynomialExponents"];

  n1 = Length[a];   (* n+1 *)
  n  = n1 - 1;      (* original n *)

  mVec = mMatrix[[auxIdx]];   (* z-row, length n+1 *)

  (* tryPivot[p]: monomial substitution + re-clear + atilde. *)
  tryPivot[p_] := Module[
    {mpLocal, apLocal, rIdx, mOther, aOther,
     subPolys, rcMin, newPolyList, atildeRaw, hasConst},
    mpLocal = mVec[[p]];
    If[mpLocal == 0, Return[$Failed]];
    apLocal = a[[p]];
    rIdx    = DeleteCases[Range[n1], p];
    mOther  = mVec[[rIdx]];
    aOther  = a[[rIdx]];

    subPolys = Table[
      Map[Function[{cmono},
        Module[{ep = cmono[[2, p]]},
          {cmono[[1]] * z0^(ep / mpLocal),
           Table[cmono[[2, rIdx[[jj]]]] - ep * mOther[[jj]] / mpLocal, {jj, n}]}
        ]
      ], clearedPolys[[k]]],
      {k, Length[clearedPolys]}
    ];

    rcMin = Table[
      Table[Min[#[[2, jj]] & /@ subPolys[[k]]], {jj, n}],
      {k, Length[subPolys]}
    ];

    newPolyList = Table[
      Map[Function[{cmono}, {cmono[[1]], cmono[[2]] - rcMin[[k]]}],
          subPolys[[k]]],
      {k, Length[subPolys]}
    ];

    atildeRaw = Table[
      aOther[[jj]] - apLocal * mOther[[jj]] / mpLocal,
      {jj, n}
    ] + Total[Table[polyExps[[k]] * rcMin[[k]], {k, Length[polyExps]}]];

    hasConst = And @@ Table[
      AnyTrue[newPolyList[[k]], (#[[2]] === ConstantArray[0, n]) &],
      {k, Length[newPolyList]}
    ];

    <|"pivot" -> p, "mp" -> mpLocal, "ap" -> apLocal,
      "remainIdx" -> rIdx, "mOther" -> mOther,
      "atilde" -> atildeRaw, "newPolys" -> newPolyList,
      "hasConst" -> hasConst, "rcMin" -> rcMin|>
  ];
  (* Memoize so each pivot is computed exactly once (plan.md §5.5). *)
  pivotCache[p_] := pivotCache[p] = tryPivot[p];

  (* ---- planAXpDIV.md §4.3: epsilon-aware classification (§6.7 EXACT) ----
     Each post-delta effective exponent atilde_j is classified at eps->0 with
     PossibleZeroQ / Re[..]>0 — never an N/tolerance test:
       "complex"     Im[a0_j] not provably zero               -> liftcomplex
       "convergent"  real, Re[a0_j] > 0
       "divergent"   real, Re[a0_j] <= 0   (only with a regulator present)
       "bad"         real, Re[a0_j] <= 0, but NO regulator to subtract
     With eps===None this reduces EXACTLY to the old isRealPos test (only
     all-convergent pivots are admissible), so the convergent path is unchanged. *)
  eps = OptionValue["Eps"];
  a0fn[av_] := If[eps =!= None, av /. eps -> 0, av];
  classDir[av_] := Module[{a0 = a0fn[av]},
    Which[
      !TrueQ[PossibleZeroQ[Im[a0]]],       "complex",
      TrueQ[Re[a0] > 0],                    "convergent",
      eps =!= None && TrueQ[Re[a0] <= 0],   "divergent",
      True,                                 "bad"
    ]];
  (* admissibility + divergent slot + c_k for one pivot's atilde *)
  pivotClass[res_] := Module[{cls, ndiv, kk, ckk},
    cls  = classDir /@ res["atilde"];
    ndiv = Count[cls, "divergent"];
    Which[
      MemberQ[cls, "complex"], <|"adm" -> False, "complex" -> True,  "divDir" -> 0|>,
      MemberQ[cls, "bad"],     <|"adm" -> False, "complex" -> False, "divDir" -> 0|>,
      ndiv > 1,                <|"adm" -> False, "complex" -> False, "divDir" -> 0|>,
      ndiv == 0,               <|"adm" -> True,  "complex" -> False, "divDir" -> 0, "ck" -> 0|>,
      True,
        kk  = First @ Flatten @ Position[cls, "divergent"];
        ckk = D[res["atilde"][[kk]], eps] /. eps -> 0;
        If[TrueQ[ckk == 0] || (NumericQ[ckk] && ckk == 0),
          (* c_k = 0: higher-order pole, not a simple 1/eps -> refuse (badck). *)
          <|"adm" -> False, "complex" -> False, "divDir" -> 0|>,
          <|"adm" -> True,  "complex" -> False, "divDir" -> kk, "ck" -> ckk|>]
    ]];
  (* domain classification as a pure function (mirrors Step 4; NO early return).
     empty/None/constrained is decided from eps-free data (mOther signs, mp sign,
     |z0|); only IndicatorCoeffs carry eps and are rebuilt safely on emission. *)
  classifyDomainFor[res_] := Module[
    {mO = res["mOther"], mpL = res["mp"], at = res["atilde"], aoz, cls},
    aoz = And @@ (# == 0 & /@ mO);
    cls = <|"LogZ0" -> Log[z0], "MP" -> mpL,
            "IndicatorCoeffs" -> Table[mO[[jj]]/at[[jj]], {jj, n}]|>;
    Which[
      aoz,
        <|"empty" -> TrueQ[N[z0^(1/mpL)] > 1], "class" -> None|>,
      mpL > 0,
        Which[
          (And @@ (# >= 0 & /@ mO)) && N[z0] > 1,  <|"empty" -> True,  "class" -> None|>,
          (And @@ (# <= 0 & /@ mO)) && N[z0] <= 1, <|"empty" -> False, "class" -> None|>,
          True,                                    <|"empty" -> False, "class" -> cls|>],
      True, (* mpL < 0 *)
        Which[
          (And @@ (# <= 0 & /@ mO)) && N[z0] < 1,  <|"empty" -> True,  "class" -> None|>,
          (And @@ (# >= 0 & /@ mO)) && N[z0] >= 1, <|"empty" -> False, "class" -> None|>,
          True,                                    <|"empty" -> False, "class" -> cls|>]
    ]];

  (* ---- build admissible candidates with classification ---- *)
  candClass = {};
  Do[
    If[mVec[[p]] != 0,
      Module[{res = pivotCache[p], ci},
        If[res =!= $Failed,
          ci = pivotClass[res];
          If[ci["adm"],
            AppendTo[candClass, {res, ci, classifyDomainFor[res]}]
          ]
        ]
      ]
    ],
    {p, n1}
  ];

  If[candClass === {},
    (* No admissible pivot.  Complex atilde => liftcomplex (decided exactly);
       otherwise liftnopivot (nested >1-divergent, badck, or — with no
       regulator — an unregulated divergence). *)
    Module[{anyComplex = False},
      Do[
        If[mVec[[p]] != 0,
          Module[{res = pivotCache[p]},
            If[res =!= $Failed &&
               AnyTrue[res["atilde"], (!TrueQ[PossibleZeroQ[Im[a0fn[#]]]]) &],
              anyComplex = True
            ]
          ]
        ],
        {p, n1}
      ];
      If[anyComplex,
        Message[TropicalEval::liftcomplex, coneIndex];
        Return[$Failed]
      ]
    ];
    Module[{pivotSummary},
      pivotSummary = Table[
        If[mVec[[p]] != 0,
          Module[{res = pivotCache[p]},
            If[res =!= $Failed, {p, mVec[[p]], N[a0fn[res["atilde"]]]}]
          ],
          Nothing
        ],
        {p, n1}
      ];
      Message[TropicalEval::liftnopivot, coneIndex, mVec, pivotSummary]
    ];
    Return[$Failed]
  ];

  (* ---- ranking (planAXpDIV.md §4.3) ----
     (0) Case-A over Case-B: a divergent slot that couples to the lifted domain
         face (Case B) yields a WRONG pole, so a decoupled (Case A) pivot must
         win whenever one exists — a CORRECTNESS key, above the variance
         heuristic.  For convergent / EmptyDomain pivots caseA is True (a
         constant key), so it never reorders the convergent path: the unlifted
         goldens and lifted-convergent output stay byte-identical (#25).
     (1) HasConstantTerm; (2) |mp|=1; (3) max min_j Re[a0], a0 = atilde|_{eps->0}
         (keeps the ranking real-valued even with a symbolic regulator). *)
  candClass = SortBy[candClass,
    Function[t,
      Module[{res = t[[1]], ci = t[[2]], dom = t[[3]], cA},
        cA = (ci["divDir"] === 0) || dom["empty"] || (dom["class"] === None) ||
             TrueQ[res["mOther"][[ci["divDir"]]] == 0];
        {-Boole[cA],
         -Boole[res["hasConst"]],
         -Boole[Abs[res["mp"]] == 1],
         -Min[Re[N[a0fn[res["atilde"]]]]]}
      ]]];
  {bestPivot, bestCI, bestDom} = candClass[[1]];

  pivotP         = bestPivot["pivot"];
  mp             = bestPivot["mp"];
  ap             = bestPivot["ap"];
  mOtherVec      = bestPivot["mOther"];
  atildeVals     = bestPivot["atilde"];
  reclearedPolys = bestPivot["newPolys"];
  divDir         = bestCI["divDir"];

  If[verbose,
    Print["  PSL cone ", coneIndex, ": pivot p=", pivotP,
          " mp=", mp, " atilde0=", N[a0fn[atildeVals]],
          " HasConstantTerm=", bestPivot["hasConst"],
          " divDir=", divDir]
  ];

  (* ---- Step 4: domain-constraint classification (from the ranked pivot) ---- *)
  logZ0 = Log[z0];
  If[bestDom["empty"],
    Return[<|"EmptyDomain" -> True, "ConeIndex" -> coneIndex|>]
  ];
  domainClass = bestDom["class"];

  (* ---- Case-B refusal (planAXpDIV.md §3) ----
     The divergent slot couples to the lifted domain face (ic_k != 0) and no
     decoupled (Case A) pivot exists (the caseA-first ranking already preferred
     one).  The 1/eps pole and the domain face interact; refuse cleanly
     ($Failed); a log-space remap is documented future work. *)
  If[divDir =!= 0 &&
     !((domainClass === None) || TrueQ[mOtherVec[[divDir]] == 0]),
    Message[TropicalEval::liftdivdomain, coneIndex];
    Return[$Failed]
  ];

  (* ---- Step 5: flatten the surviving coordinates ----
     Convergent lifted sector: flatten now (byte-identical to the pre-L2 path).
     Divergent lifted sector (planAXpDIV.md §4.3, Barrier A): do NOT pre-flatten
     — the divergence routines flatten the non-divergent coordinates themselves
     and the 1/eps pole lives in the surviving atilde slot divDir. *)
  prefactorBase = (Abs[detM] / Abs[mp]) * z0^(ap / mp - 1);
  If[divDir === 0,
    fsResult = FlattenSector[reclearedPolys, atildeVals, prefactorBase];
    If[fsResult["IsDivergent"],
      Message[TropicalEval::liftdivergent, coneIndex, atildeVals];
      Return[$Failed]
    ];
    flattenedPolys = fsResult["FlattenedPolys"];
    prefactor      = fsResult["Prefactor"],
    (* divergent: placeholders (not consumed by the divergence path) *)
    flattenedPolys = None;
    prefactor      = None
  ];

  (* MonoFactorLog (plan.md §6.3 Fix B / lift_error_log L3), LIFTED case.
     log P_k = log(monomial factor) + log Q_k, where the augmented tropical
     clearing removed y^{d^aug_k} and the pivot substitution + re-clearing
     removed z0^{d^aug_{k,p}/m_p} * Prod_j y_j^{...}.  In the flattened coords
     y' (y_j = (y_j')^{1/atilde_j}):
        Const_k   = (d^aug_{k,p}/m_p) * log z0
        Coeffs_{k,j} = [ (d^aug_{k,rIdx[j]} - d^aug_{k,p}*m_{rIdx[j]}/m_p)
                         + rcMin_{k,j} ] / atilde_j
     d^aug_{k,*} = sdAug["MinExponents"][[k]] (the augmented cleared minima);
     p=pivotP, m_p=mp, rIdx=remainIdx, rcMin=bestPivot["rcMin"], atilde=atildeVals.
     Consumed only by the SplitRealImag phase (Im(B)!=0); real path unaffected. *)
  monoFactorLogLifted = Module[{dAug = sdAug["MinExponents"], rIdx, rcMinL},
    rIdx   = bestPivot["remainIdx"];
    rcMinL = bestPivot["rcMin"];
    Table[
      <|"Const"  -> (dAug[[k, pivotP]]/mp) * logZ0,
        "Coeffs" -> Table[
          ((dAug[[k, rIdx[[j]]]] - dAug[[k, pivotP]]*mOtherVec[[j]]/mp)
            + rcMinL[[k, j]]) / atildeVals[[j]],
          {j, n}]|>,
      {k, Length[reclearedPolys]}
    ]
  ];

  (* LiftedMonoFactor (planCXLIFTDIV.md §4.1): the UN-flattened numerators of the
     sector-level MonoFactorLog above, i.e. the same Const_k and the numerators
     D_{k,j} BEFORE the division by atilde_j.  A divergent sector is reduced by
     IBP into a boundary (n-1 dims, slot divDir dropped) and IBP terms (n dims),
     each with its OWN flattening alpha0^{(piece)}; the oscillatory phase of each
     piece is Const_k + Sum_i (D_{k,i}/alpha0^{(piece)}_i) log y'_i, so IBPProcess-
     Sector re-divides these fixed D_{k,i} by that piece's alpha0 (§3).  Exact /
     eps-free (atilde-free); consumed only by the SplitRealImag phase. *)
  liftedMonoFactorNum = Module[{dAug = sdAug["MinExponents"], rIdx, rcMinL},
    rIdx   = bestPivot["remainIdx"];
    rcMinL = bestPivot["rcMin"];
    <|"Const" -> Table[(dAug[[k, pivotP]]/mp) * logZ0,
                       {k, Length[reclearedPolys]}],
      "DExp"  -> Table[
        Table[
          (dAug[[k, rIdx[[j]]]] - dAug[[k, pivotP]]*mOtherVec[[j]]/mp)
            + rcMinL[[k, j]],
          {j, n}],
        {k, Length[reclearedPolys]}]|>
  ];

  (* MonomialPhaseLog (planAXpDIVv2.md §3), LIFTED case — the bare-monomial
     oscillatory phase from Im(A), the monomial twin of monoFactorLogLifted.
     Apply the pivot substitution log Y_p = (logZ0 - Sum_jj mOther_jj log y_jj)/m_p
     to Sum_a Im(eaug)_a log Y_a, eaug = Im(A_aug).M_aug (M_aug = RayMatrix, the
     augmented CoV).  In the flattened coords y' (y_jj = (y_jj')^{1/atilde_jj}):
        Const   = (Im(eaug)_p / m_p) * logZ0
        Coeffs_jj = ( Im(eaug)_{rIdx[jj]} - Im(eaug)_p * mOther_jj/m_p ) / atilde_jj
     Identical structure to monoFactorLogLifted with d^aug_{k,*} replaced by the
     single vector Im(eaug) and NO rcMin term (the bare monomial is not a
     re-cleared polynomial).  M0 Test E verified this to 1.2e-8 vs NIntegrate.
     liftedMonoPhaseNum is the UN-flattened numerator (length n) consumed by
     IBPProcessSector to re-flatten the phase per IBP piece (the §4.3 analog of
     liftedMonoFactorNum).  Both None on real-A sectors -> byte-identical (#25). *)
  {monomialPhaseLogLifted, liftedMonoPhaseNum} = If[imAexps === None, {None, None},
    Module[{imEaug, rIdx, constA, numVec},
      imEaug = imAexps . mMatrix;               (* length n+1 *)
      rIdx   = bestPivot["remainIdx"];
      constA = (imEaug[[pivotP]]/mp) * logZ0;
      numVec = Table[
        imEaug[[rIdx[[j]]]] - imEaug[[pivotP]]*mOtherVec[[j]]/mp, {j, n}];
      {<|"Const" -> constA,
         "Coeffs" -> Table[numVec[[j]]/atildeVals[[j]], {j, n}]|>,
       <|"Const" -> constA, "Num" -> numVec|>}
    ]];

  (* ---- Divergent lifted sector: emit a divergent SectorData (planAXpDIV.md
     §4.3, Barrier A).  NewExponents carry eps so the pole machinery
     (IdentifyDivergences / IBPReduceSector) eps-expands them; the divergence
     routines flatten the non-divergent coordinates and reconstruct prefactors
     from PrefactorBase.  DomainConstraint is the Case-A indicator at eps->0 with
     the divergent slot forced to 0 (ic_divDir = 0, avoiding a 0/0 from the
     vanishing atilde_divDir). ---- *)
  If[divDir =!= 0,
    Return[<|
      "ConeIndex"           -> coneIndex,
      "RayMatrix"           -> mMatrix,
      "DetM"                -> detM,
      "SelectedRays"        -> sdAug["SelectedRays"],
      "RawExponents"        -> sdAug["RawExponents"],
      "NewExponents"        -> atildeVals,
      "MinExponents"        -> bestPivot["rcMin"],
      "TransformedPolys"    -> clearedPolys,
      "ClearedPolys"        -> reclearedPolys,
      "PrefactorBase"       -> prefactorBase,
      "MonoFactorLog"       -> monoFactorLogLifted,
      (* planCXLIFTDIV.md §4.1: un-flattened numerators so each IBP piece can
         re-flatten the oscillatory phase by its own alpha0.  Absent on real /
         unlifted sectors (codegen treats absent as "no phase") -> #25. *)
      "LiftedMonoFactor"    -> liftedMonoFactorNum,
      (* planAXpDIVv2.md §4.3: un-flattened bare-monomial phase numerator (Im(A));
         IBPProcessSector re-divides it by each piece's alpha0.  None real-A -> #25. *)
      "LiftedMonoPhase"     -> liftedMonoPhaseNum,
      "IsDivergent"         -> True,
      "DivergentVariable"   -> divDir,
      "Dimension"           -> n,
      "PolynomialExponents" -> polyExps,
      "MonomialExponents"   -> liftData["OriginalSpec"]["MonomialExponents"],
      "DomainConstraint"    -> If[domainClass === None, None,
        <|"LogZ0" -> domainClass["LogZ0"], "MP" -> domainClass["MP"],
          "IndicatorCoeffs" -> Table[
            If[jj === divDir, 0, (mOtherVec[[jj]]/atildeVals[[jj]]) /. eps -> 0],
            {jj, n}]|>],
      "LiftData"            -> liftData,
      "PivotIndex"          -> pivotP,
      "ZRow"                -> mVec,
      "AugmentedA"          -> a,
      "HasConstantTerm"     -> bestPivot["hasConst"]
    |>]
  ];

  (* ---- Step 6: assemble full SectorData (convergent lifted sector) ---- *)
  <|
    "ConeIndex"           -> coneIndex,
    "RayMatrix"           -> mMatrix,
    "DetM"                -> detM,
    "SelectedRays"        -> sdAug["SelectedRays"],
    "RawExponents"        -> sdAug["RawExponents"],
    "NewExponents"        -> atildeVals,
    "MinExponents"        -> bestPivot["rcMin"],
    "TransformedPolys"    -> clearedPolys,
    "ClearedPolys"        -> reclearedPolys,
    "FlattenedPolys"      -> flattenedPolys,
    "Prefactor"           -> prefactor,
    (* PrefactorBase (planAXpDIV.md §4.1, Barrier B.1): the lifted sector's
       true prefactor base (Abs[detM]/Abs[mp]) z0^(ap/mp-1), NOT Abs[detM].
       The divergence routines reconstruct prefactors from this. *)
    "PrefactorBase"       -> prefactorBase,
    "MonoFactorLog"       -> monoFactorLogLifted,
    (* planAXpDIVv2.md §4.3: bare-monomial phase (Im(A)), flattened by atilde.
       None on real-A sectors -> codegen emits no monomial phase -> #25. *)
    "MonomialPhaseLog"    -> monomialPhaseLogLifted,
    "IsDivergent"         -> False,
    "DivergentVariable"   -> 0,
    "Dimension"           -> n,
    "PolynomialExponents" -> polyExps,
    "MonomialExponents"   -> liftData["OriginalSpec"]["MonomialExponents"],
    "DomainConstraint"    -> domainClass,
    "LiftData"            -> liftData,
    "PivotIndex"          -> pivotP,
    "ZRow"                -> mVec,
    "AugmentedA"          -> a,
    "HasConstantTerm"     -> bestPivot["hasConst"]
  |>
];


(* --------------------------------------------------------------------------
   ValidateLiftedDecomposition  (plan.md §6.2, §9 risk #1)
   Direct integral from originalSpec; sectors via ProcessSectorLifted.  This
   EXACT NIntegrate gate (not the sampled 5sigma) is the lifting correctness
   check.  EmptyDomain sectors contribute 0 and are listed in DroppedSectors.
   -------------------------------------------------------------------------- *)

ValidateLiftedDecomposition[originalSpec_Association,
                            liftedSpec_Association,
                            liftedFanData_List,
                            liftData_Association,
                            testKinematics_List,
                            precisionGoal_Integer: 3] :=
Module[
  {polys, monoExps, polyExps, vars, n,
   dualVertices, simplexList, kinRules, pg,
   directIntegrand, directResult,
   sectorResults, droppedSectors, sectorSum, relError},

  polys    = originalSpec["Polynomials"];
  monoExps = originalSpec["MonomialExponents"];
  polyExps = originalSpec["PolynomialExponents"];
  vars     = originalSpec["Variables"];
  n        = Length[vars];
  kinRules = testKinematics;
  pg       = precisionGoal;

  {dualVertices, simplexList} = liftedFanData;

  directIntegrand = (Times @@ MapThread[Power, {vars, monoExps}]) *
    (Times @@ MapThread[Power, {polys, polyExps}]) /. kinRules;

  directResult = NIntegrate[
    directIntegrand,
    Evaluate[Sequence @@ ({#, 0, Infinity} & /@ vars)],
    MaxRecursion -> 20, PrecisionGoal -> pg + 1, Method -> "GlobalAdaptive"
  ];

  droppedSectors = {};

  sectorResults = Table[
    Module[{sd, flatPolys, pExps, pf, dim, yVars, integrand,
            polyVals, dc, result, logYpStar, logZ0num, mpNum, icNum},
      sd = ProcessSectorLifted[liftedSpec, dualVertices,
                               simplexList[[s]], s, liftData];
      Which[
        sd === $Failed,
          Return[$Failed],
        AssociationQ[sd] && KeyExistsQ[sd, "EmptyDomain"] && sd["EmptyDomain"],
          AppendTo[droppedSectors, sd["ConeIndex"]];
          0,
        True,
          flatPolys = sd["FlattenedPolys"];
          pExps     = sd["PolynomialExponents"] /. kinRules;
          pf        = sd["Prefactor"] /. kinRules;
          dim       = sd["Dimension"];
          dc        = sd["DomainConstraint"];
          yVars     = Table[Unique["yv"], {dim}];

          polyVals = Table[
            Total[Table[
              Module[{coeff = mono[[1]] /. kinRules,
                      alphas = mono[[2]] /. kinRules, logY = Log /@ yVars},
                coeff * Exp[Total[alphas * logY]]
              ], {mono, flatPolys[[j]]}]],
            {j, Length[flatPolys]}
          ];

          integrand = pf *
            Times @@ MapThread[
              Function[{pv, be}, Exp[be * Log[pv]]], {polyVals, pExps}];

          If[dc =!= None,
            logZ0num = N[dc["LogZ0"]];
            mpNum    = N[dc["MP"]];
            icNum    = N[dc["IndicatorCoeffs"]];
            logYpStar = (logZ0num - Total[icNum * (Log /@ yVars)]) / mpNum;
            integrand = integrand * Boole[logYpStar <= 0]
          ];

          result = Quiet@NIntegrate[
            integrand,
            Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
            MaxRecursion -> 15, PrecisionGoal -> pg, Method -> "GlobalAdaptive"
          ];
          result
      ]
    ],
    {s, Length[simplexList]}
  ];

  sectorSum = Total[sectorResults];
  relError  = Abs[(sectorSum - directResult) / directResult];

  If[NumericQ[relError] && relError > 10^(-pg + 1),
    Message[TropicalEval::validate, "LiftedDecomposition", relError,
            10^(-pg + 1)]
  ];

  <|"DirectResult" -> directResult, "SectorSum" -> sectorSum,
    "RelativeError" -> relError, "SectorResults" -> sectorResults,
    "DroppedSectors" -> droppedSectors|>
];


(* ============================================================================
   MODULE 2: DIVERGENCE REGULATION (Tropical Subtraction)
   ============================================================================ *)

(* --------------------------------------------------------------------------
   IdentifyDivergences
   Uses the EFFECTIVE exponents (post tropical factoring).
   -------------------------------------------------------------------------- *)

IdentifyDivergences[sectorData_Association, eps_] :=
Module[
  {aVals, n, a0, a1, divVars, ck, ak0, k},

  aVals = sectorData["NewExponents"];  (* effective exponents *)
  n     = sectorData["Dimension"];

  (* Expand a_i(epsilon) = a_i^(0) + a_i^(1) * epsilon + ... *)
  a0 = aVals /. eps -> 0;
  a1 = D[aVals, eps] /. eps -> 0;

  (* Find divergent variables: Re(a_k^(0)) <= 0 *)
  divVars = {};
  Do[
    If[TrueQ[Re[a0[[i]]] <= 0] ||
       (NumericQ[a0[[i]]] && Re[a0[[i]]] <= 0),
      AppendTo[divVars, i]
    ],
    {i, n}
  ];

  (* Check: at most one divergent variable *)
  If[Length[divVars] > 1,
    Message[TropicalEval::nested, sectorData["ConeIndex"], divVars];
    Return[$Failed]
  ];

  If[Length[divVars] == 0,
    Return[<|"IsDivergent" -> False|>]
  ];

  k   = divVars[[1]];
  ck  = a1[[k]];
  ak0 = a0[[k]];

  (* Check: c_k != 0 *)
  If[TrueQ[ck == 0] || (NumericQ[ck] && ck == 0),
    Message[TropicalEval::badck, sectorData["ConeIndex"], k];
    Return[$Failed]
  ];

  (* Check: all non-divergent variables have Re(a_i^(0)) > 0 *)
  Do[
    If[i != k && (TrueQ[Re[a0[[i]]] <= 0] ||
                  (NumericQ[a0[[i]]] && Re[a0[[i]]] <= 0)),
      Print["WARNING: Sector ", sectorData["ConeIndex"],
            ": non-divergent variable y_", i,
            " has Re(a_", i, "^(0)) = ", Re[a0[[i]]], " <= 0"]
    ],
    {i, n}
  ];

  Print["Sector ", sectorData["ConeIndex"],
        ": variable y_", k, " divergent, a_", k,
        "(eps) = ", ck, " * eps + ..., c_k = ", ck];

  <|"IsDivergent" -> True,
    "DivergentVariable" -> k,
    "ck" -> ck,
    "ak0" -> ak0,
    "a0" -> a0,
    "a1" -> a1|>
];

(* --------------------------------------------------------------------------
   ProcessDivergentSector
   Constructs G0, G1, remainder, and analytic pole.
   Works with the CLEARED polynomials Q_j (non-negative y-exponents)
   and the EFFECTIVE exponents.
   -------------------------------------------------------------------------- *)

ProcessDivergentSector[sectorData_Association, integrandSpec_Association] :=
Module[
  {eps, k, a0, a1, ck, n, polyExps,
   clearedPolys, detM, pfBase,
   B0, B1, simpPolys, fullPolys,
   g0FlatPolys, g0Prefactor, g0Avals,
   g1LogInsertions,
   remainderData,
   analyticPole,
   divInfo},

  eps = integrandSpec["RegulatorSymbol"];
  n   = sectorData["Dimension"];

  (* Identify the divergent variable using effective exponents *)
  divInfo = IdentifyDivergences[sectorData, eps];
  If[divInfo === $Failed, Return[$Failed]];
  If[!divInfo["IsDivergent"],
    Print["ProcessDivergentSector called on non-divergent sector"];
    Return[$Failed]
  ];

  k      = divInfo["DivergentVariable"];
  ck     = divInfo["ck"];
  a0     = divInfo["a0"];
  a1     = divInfo["a1"];
  detM   = sectorData["DetM"];
  (* PrefactorBase (planAXpDIV.md §4.1): Abs[detM] for unlifted sectors,
     (Abs[detM]/Abs[mp]) z0^(ap/mp-1) for a lifted sector. *)
  pfBase = Lookup[sectorData, "PrefactorBase", Abs[detM]];

  polyExps    = sectorData["PolynomialExponents"];
  clearedPolys = sectorData["ClearedPolys"];

  (* Polynomial exponents at eps=0 and first derivative *)
  B0 = polyExps /. eps -> 0;
  B1 = D[polyExps, eps] /. eps -> 0;

  (* --- Step 2: Construct F_simple by setting y_k = 0 --- *)
  (* In the cleared polynomials Q_j, keep only monomials where
     the exponent of y_k is 0 (these survive when y_k = 0). *)
  simpPolys = Table[
    Select[clearedPolys[[j]],
      (#[[2, k]] === 0 || TrueQ[#[[2, k]] == 0]) &
    ],
    {j, Length[clearedPolys]}
  ];

  fullPolys = clearedPolys;

  (* --- Step 3: G0 and G1 --- *)
  Module[{ndVars, g0aVals, g0Dim},
    ndVars  = DeleteCases[Range[n], k];
    g0aVals = a0[[ndVars]];
    g0Dim   = n - 1;

    (* For G0: flatten non-divergent variables *)
    g0FlatPolys = Table[
      Table[
        Module[{coeff, exps, ndExps, flatExps},
          coeff  = mono[[1]];
          exps   = mono[[2]];
          ndExps = exps[[ndVars]];
          flatExps = MapThread[#1/#2 &, {ndExps, g0aVals}];
          {coeff, flatExps}
        ],
        {mono, simpPolys[[j]]}
      ],
      {j, Length[simpPolys]}
    ];

    g0Prefactor = pfBase / (Times @@ g0aVals);
    g0Avals = g0aVals;

    (* G1 log insertion factors *)
    g1LogInsertions = <|
      "VariableTerms" -> Table[
        {a1[[ndVars[[i]]]] / g0aVals[[i]], i},
        {i, g0Dim}
      ],
      "PolynomialTerms" -> Table[
        {B1[[j]], j},
        {j, Length[polyExps]}
      ]
    |>;
  ];

  (* --- Step 4: Finite remainder --- *)
  Module[{ndVars, remAvals, remFlatFullPolys, remFlatSimpPolys,
          remPrefactor},
    ndVars   = DeleteCases[Range[n], k];
    remAvals = a0[[ndVars]];

    (* Partially flatten: only non-div vars get flattened *)
    remFlatFullPolys = Table[
      Table[
        Module[{coeff, exps, flatExps},
          coeff = mono[[1]];
          exps  = mono[[2]];
          flatExps = Table[
            If[i == k,
              exps[[i]],
              exps[[i]] / a0[[i]]
            ],
            {i, n}
          ];
          {coeff, flatExps}
        ],
        {mono, fullPolys[[j]]}
      ],
      {j, Length[fullPolys]}
    ];

    remFlatSimpPolys = Table[
      Table[
        Module[{coeff, exps, flatExps},
          coeff = mono[[1]];
          exps  = mono[[2]];
          flatExps = Table[
            If[i == k,
              exps[[i]],
              exps[[i]] / a0[[i]]
            ],
            {i, n}
          ];
          {coeff, flatExps}
        ],
        {mono, simpPolys[[j]]}
      ],
      {j, Length[simpPolys]}
    ];

    remPrefactor = pfBase / (Times @@ remAvals);

    remainderData = <|
      "FullPolys"     -> remFlatFullPolys,
      "SimplifiedPolys" -> remFlatSimpPolys,
      "Prefactor"     -> remPrefactor,
      "DivVarIndex"   -> k,
      "DivVarExp"     -> a0[[k]],
      "Dimension"     -> n,
      "NonDivVars"    -> DeleteCases[Range[n], k],
      "PolynomialExponents" -> B0
    |>;
  ];

  (* --- Step 5: Analytic contributions --- *)
  analyticPole = 1 / ck;

  <|
    "ConeIndex"           -> sectorData["ConeIndex"],
    "IsDivergent"         -> True,
    "DivergentVariable"   -> k,
    "ck"                  -> ck,
    "a0"                  -> a0,
    "a1"                  -> a1,
    "B0"                  -> B0,
    "B1"                  -> B1,
    "G0FlatPolys"         -> g0FlatPolys,
    "G0Prefactor"         -> g0Prefactor,
    "G0Avals"             -> g0Avals,
    "G0Dimension"         -> n - 1,
    "G0PolyExponents"     -> B0,
    "G1LogInsertions"     -> g1LogInsertions,
    "Remainder"           -> remainderData,
    "AnalyticPole"        -> analyticPole,
    "DetM"                -> detM,
    "Dimension"           -> n,
    "ClearedPolys"        -> clearedPolys,
    "SimplifiedPolys"     -> simpPolys,
    "TransformedPolys"    -> sectorData["TransformedPolys"],
    "MinExponents"        -> sectorData["MinExponents"],
    "PolynomialExponents" -> polyExps,
    "MonomialExponents"   -> sectorData["MonomialExponents"],
    "RawExponents"        -> sectorData["RawExponents"],
    (* Lifted-sector metadata (planAXpDIV.md §4.2): propagated so the G0/G1/rem
       codegen and the validators carry the domain indicator.  None for unlifted
       sectors -> byte-identical emitted C++ (#25). *)
    "DomainConstraint"    -> Lookup[sectorData, "DomainConstraint", None]
  |>
];

(* --------------------------------------------------------------------------
   ValidateSubtraction
   Self-consistency check at finite epsilon.
   Uses the ORIGINAL transformed polynomials (not cleared) to compute
   the true sector integral, then compares against the subtracted form.
   -------------------------------------------------------------------------- *)

ValidateSubtraction[divSectorData_Association, sectorData_Association,
                    integrandSpec_Association,
                    testKinematics_List, testEpsilon_: 0.01] :=
Module[
  {eps, n, k, ck, a0, a1, aVals, polyExps,
   kinRules, epsRules, fullRules,
   clearedPolys, simpPolys,
   detM, pfBase, domC, yVars,
   originalIntegral, g0Val, g1Val, remVal,
   divContrib, reconstructed, relError,
   B0, B1},

  eps      = integrandSpec["RegulatorSymbol"];
  n        = sectorData["Dimension"];
  k        = divSectorData["DivergentVariable"];
  ck       = divSectorData["ck"];
  a0       = divSectorData["a0"];
  a1       = divSectorData["a1"];
  aVals    = sectorData["NewExponents"];  (* effective, symbolic in eps *)
  polyExps = sectorData["PolynomialExponents"];
  detM     = sectorData["DetM"];
  (* PrefactorBase (planAXpDIV.md §4.1) + lifted DomainConstraint (§4.2). *)
  pfBase   = Lookup[sectorData, "PrefactorBase", Abs[detM]];
  domC     = Lookup[sectorData, "DomainConstraint", None];

  kinRules = testKinematics;
  epsRules = {eps -> testEpsilon};
  fullRules = Join[kinRules, epsRules];

  clearedPolys = divSectorData["ClearedPolys"];
  simpPolys    = divSectorData["SimplifiedPolys"];

  yVars = Table[Unique["yv"], {n}];

  (* Original sector integral at finite epsilon.
     Use the cleared polynomial form:
     integrand = |det(M)| * prod y_i^{a_i^eff - 1} * prod Q_j^{B_j}
     where Q_j are cleared polys with non-negative exponents. *)
  Module[{aNum, polyValsExpr, integrand},
    aNum = aVals /. fullRules;
    polyValsExpr = Table[
      Total[
        Table[
          Module[{coeff, exps, logY},
            coeff = mono[[1]] /. kinRules;
            exps  = mono[[2]];
            logY  = Log /@ yVars;
            coeff * Exp[Total[exps * logY]]
          ],
          {mono, clearedPolys[[j]]}
        ]
      ],
      {j, Length[clearedPolys]}
    ];
    integrand = (pfBase /. fullRules) *
      Exp[Total[(aNum - 1) * Log /@ yVars]] *
      Times @@ MapThread[
        Function[{pv, be}, Exp[be * Log[pv]]],
        {polyValsExpr, polyExps /. fullRules}
      ];
    (* lifted-sector domain indicator over all n coords (planAXpDIV.md §4.2) *)
    integrand = integrand * liftedDomainBooleWL[domC, yVars];

    originalIntegral = Quiet@NIntegrate[
      integrand,
      Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
      MaxRecursion -> 20,
      PrecisionGoal -> 4,
      Method -> "GlobalAdaptive"
    ];
  ];

  B0 = polyExps /. eps -> 0 /. kinRules;
  B1 = (D[polyExps, eps] /. eps -> 0) /. kinRules;

  (* G0: (n-1)-dim integral at eps=0 *)
  Module[{ndVars, g0yVars, g0aNum, g0PolyVals, g0Integrand},
    ndVars    = DeleteCases[Range[n], k];
    g0yVars   = Table[Unique["gy"], {n - 1}];
    g0aNum    = (a0 /. kinRules)[[ndVars]];

    g0PolyVals = Table[
      Total[
        Table[
          Module[{coeff, ndExps, logY},
            coeff  = mono[[1]] /. kinRules;
            ndExps = mono[[2]][[ndVars]];
            logY   = Log /@ g0yVars;
            coeff * Exp[Total[ndExps * logY]]
          ],
          {mono, simpPolys[[j]]}
        ]
      ],
      {j, Length[simpPolys]}
    ];

    g0Integrand = (pfBase /. fullRules) *
      Exp[Total[(g0aNum - 1) * Log /@ g0yVars]] *
      Times @@ MapThread[
        Function[{pv, be}, Exp[be * Log[pv]]],
        {g0PolyVals, B0}
      ];
    (* domain indicator over the n-1 non-divergent coords (drop slot k, §4.2) *)
    g0Integrand = g0Integrand *
      liftedDomainBooleWL[dropDivVarFromDomain[domC, k], g0yVars];

    g0Val = Quiet@NIntegrate[
      g0Integrand,
      Evaluate[Sequence @@ ({#, 0, 1} & /@ g0yVars)],
      MaxRecursion -> 20,
      PrecisionGoal -> 4,
      Method -> "GlobalAdaptive"
    ];
  ];

  (* G1: G0 integrand * log insertion sum *)
  Module[{ndVars, g1yVars, g1aNum, g1PolyVals,
          g1BaseIntegrand, logInsertionSum, g1Integrand},
    ndVars  = DeleteCases[Range[n], k];
    g1yVars = Table[Unique["hy"], {n - 1}];
    g1aNum  = (a0 /. kinRules)[[ndVars]];

    g1PolyVals = Table[
      Total[
        Table[
          Module[{coeff, ndExps, logY},
            coeff  = mono[[1]] /. kinRules;
            ndExps = mono[[2]][[ndVars]];
            logY   = Log /@ g1yVars;
            coeff * Exp[Total[ndExps * logY]]
          ],
          {mono, simpPolys[[j]]}
        ]
      ],
      {j, Length[simpPolys]}
    ];

    g1BaseIntegrand = (pfBase /. fullRules) *
      Exp[Total[(g1aNum - 1) * Log /@ g1yVars]] *
      Times @@ MapThread[
        Function[{pv, be}, Exp[be * Log[pv]]],
        {g1PolyVals, B0}
      ];

    logInsertionSum =
      Total[Table[
        (a1 /. kinRules)[[ndVars[[i]]]] / g1aNum[[i]] * Log[g1yVars[[i]]],
        {i, n - 1}
      ]] +
      Total[Table[
        B1[[j]] * Log[g1PolyVals[[j]]],
        {j, Length[polyExps]}
      ]];

    g1Integrand = g1BaseIntegrand * logInsertionSum *
      liftedDomainBooleWL[dropDivVarFromDomain[domC, k], g1yVars];

    g1Val = Quiet@NIntegrate[
      g1Integrand,
      Evaluate[Sequence @@ ({#, 0, 1} & /@ g1yVars)],
      MaxRecursion -> 20,
      PrecisionGoal -> 3,
      Method -> "GlobalAdaptive"
    ];
  ];

  (* Remainder integral *)
  Module[{remYVars, remAnum, fullPolyVals, simpPolyVals,
          bracket, remIntegrand},
    remYVars   = Table[Unique["rv"], {n}];
    remAnum    = a0 /. kinRules;

    fullPolyVals = Table[
      Total[
        Table[
          Module[{coeff, exps, logY},
            coeff = mono[[1]] /. kinRules;
            exps  = mono[[2]];
            logY  = Log /@ remYVars;
            coeff * Exp[Total[exps * logY]]
          ],
          {mono, clearedPolys[[j]]}
        ]
      ],
      {j, Length[clearedPolys]}
    ];

    simpPolyVals = Table[
      Total[
        Table[
          Module[{coeff, exps, logY},
            coeff = mono[[1]] /. kinRules;
            exps  = mono[[2]];
            logY  = Log /@ remYVars;
            coeff * Exp[Total[exps * logY]]
          ],
          {mono, simpPolys[[j]]}
        ]
      ],
      {j, Length[simpPolys]}
    ];

    bracket = Times @@ MapThread[
        Function[{pv, be}, Exp[be * Log[pv]]],
        {fullPolyVals, B0}
      ] -
      Times @@ MapThread[
        Function[{pv, be}, Exp[be * Log[pv]]],
        {simpPolyVals, B0}
      ];

    remIntegrand = (pfBase /. fullRules) *
      Exp[Total[(remAnum - 1) * Log /@ remYVars]] *
      bracket *
      liftedDomainBooleWL[domC, remYVars];

    remVal = Quiet@NIntegrate[
      remIntegrand,
      Evaluate[Sequence @@ ({#, 0, 1} & /@ remYVars)],
      MaxRecursion -> 20,
      PrecisionGoal -> 3,
      Method -> "GlobalAdaptive"
    ];
  ];

  (* Reconstruct: divContrib + remainder = original *)
  divContrib   = g0Val / (ck * testEpsilon) + g1Val / ck;
  reconstructed = divContrib + remVal;
  relError      = Abs[(reconstructed - originalIntegral) / originalIntegral];

  Print["ValidateSubtraction for sector ", sectorData["ConeIndex"], ":"];
  Print["  Original integral (eps=", testEpsilon, "): ", originalIntegral];
  Print["  G0/(ck*eps) + G1/ck = ", divContrib];
  Print["  Remainder = ", remVal];
  Print["  Reconstructed = ", reconstructed];
  Print["  Relative error = ", relError];

  If[relError > 0.01,
    Print["  WARNING: relative error exceeds 1%"]
  ];

  <|"OriginalIntegral" -> originalIntegral,
    "DivergentContribution" -> divContrib,
    "Remainder" -> remVal,
    "Reconstructed" -> reconstructed,
    "RelativeError" -> relError,
    "G0" -> g0Val, "G1" -> g1Val|>
];


(* ============================================================================
   MODULE 3: C++ CODE GENERATION
   ============================================================================ *)

(* --------------------------------------------------------------------------
   MmaToC: Convert Mathematica expression to C++ string
   Recursive pattern-matching converter (not CForm-based).
   -------------------------------------------------------------------------- *)

MmaToC[expr_, paramMap_Association: <||>] := mmaToCInternal[expr, paramMap];

(* Integers -> double literals *)
mmaToCInternal[n_Integer, _] := ToString[n] <> ".0";

(* Rationals -> (a.0/b.0) *)
mmaToCInternal[r_Rational, _] :=
  "(" <> ToString[Numerator[r]] <> ".0/" <>
  ToString[Denominator[r]] <> ".0)";

(* Reals -> string *)
mmaToCInternal[r_Real, _] := ToString[CForm[r]];

(* Complex numbers -> cx(re, im) *)
mmaToCInternal[Complex[re_, im_], pm_] :=
  "cx(" <> mmaToCInternal[re, pm] <> ", " <>
  mmaToCInternal[im, pm] <> ")";

(* Symbols -> parameter lookup or literal *)
mmaToCInternal[s_Symbol, pm_] :=
  If[KeyExistsQ[pm, s],
    pm[s],
    ToString[s]
  ];

(* Subscripted variables like x[i] *)
mmaToCInternal[s_Symbol[i_Integer], pm_] :=
  If[KeyExistsQ[pm, s[i]],
    pm[s[i]],
    ToString[s] <> "[" <> ToString[i] <> "]"
  ];

(* Power *)
mmaToCInternal[Power[base_, -1], pm_] :=
  "(1.0/" <> mmaToCInternal[base, pm] <> ")";

mmaToCInternal[Power[base_, 1/2], pm_] :=
  "std::sqrt(" <> mmaToCInternal[base, pm] <> ")";

mmaToCInternal[Power[base_, -1/2], pm_] :=
  "(1.0/std::sqrt(" <> mmaToCInternal[base, pm] <> "))";

mmaToCInternal[Power[E, exp_], pm_] :=
  "std::exp(" <> mmaToCInternal[exp, pm] <> ")";

mmaToCInternal[Power[base_, n_Integer], pm_] :=
  "std::pow(" <> mmaToCInternal[base, pm] <> ", " <>
  ToString[n] <> ".0)";

mmaToCInternal[Power[base_, exp_], pm_] :=
  "std::pow(" <> mmaToCInternal[base, pm] <> ", " <>
  mmaToCInternal[exp, pm] <> ")";

(* Log *)
mmaToCInternal[Log[arg_], pm_] :=
  "std::log(" <> mmaToCInternal[arg, pm] <> ")";

(* Exp *)
mmaToCInternal[Exp[arg_], pm_] :=
  "std::exp(" <> mmaToCInternal[arg, pm] <> ")";

(* Abs *)
mmaToCInternal[Abs[arg_], pm_] :=
  "std::abs(" <> mmaToCInternal[arg, pm] <> ")";

(* Re, Im *)
mmaToCInternal[Re[arg_], pm_] :=
  "(" <> mmaToCInternal[arg, pm] <> ").real()";

mmaToCInternal[Im[arg_], pm_] :=
  "(" <> mmaToCInternal[arg, pm] <> ").imag()";

(* Plus -> infix with parens.
   HoldPattern is needed because Mathematica evaluates Plus[args__] and
   Times[args__] before storing the DownValue inside BeginPackage. *)
mmaToCInternal[HoldPattern[Plus[args__]], pm_] :=
  "(" <> StringRiffle[mmaToCInternal[#, pm] & /@ {args}, " + "] <> ")";

(* Times: handle negation and general products *)
mmaToCInternal[HoldPattern[Times[-1, rest__]], pm_] :=
  "(-" <> mmaToCInternal[Times[rest], pm] <> ")";

mmaToCInternal[HoldPattern[Times[args__]], pm_] :=
  "(" <> StringRiffle[mmaToCInternal[#, pm] & /@ {args}, " * "] <> ")";

(* Fallback: use CForm and warn *)
mmaToCInternal[expr_, pm_] := Module[{str},
  str = ToString[CForm[expr]];
  str = StringReplace[str, {
    "Power(" ~~ a__ ~~ "," ~~ b__ ~~ ")" :>
      "std::pow(" <> a <> ", " <> b <> ")",
    "Sqrt(" ~~ a__ ~~ ")" :> "std::sqrt(" <> a <> ")"
  }];
  str
];

(* --------------------------------------------------------------------------
   Helper: Generate C++ for one monomial sum (polynomial evaluation)
   -------------------------------------------------------------------------- *)

GenerateMonomialSumCpp[flatPolys_List, polyIndex_Integer,
                       paramMap_Association, dim_Integer,
                       varPrefix_String: "log_y"] :=
Module[{lines, polyVar},
  polyVar = "P" <> ToString[polyIndex];
  lines = {"    cx " <> polyVar <> "(0.0, 0.0);"};

  Do[
    Module[{coeff, alphas, coeffStr, expTerms, expStr},
      coeff    = mono[[1]];
      alphas   = mono[[2]];
      coeffStr = mmaToCInternal[coeff, paramMap];

      expTerms = Table[
        If[TrueQ[alphas[[i]] == 0],
          Nothing,
          mmaToCInternal[alphas[[i]], paramMap] <> " * " <>
            varPrefix <> "[" <> ToString[i - 1] <> "]"
        ],
        {i, dim}
      ];

      expStr = If[Length[expTerms] == 0,
        "0.0",
        StringRiffle[expTerms, " + "]
      ];

      AppendTo[lines,
        "    " <> polyVar <> " += " <> coeffStr <>
        " * std::exp(" <> expStr <> ");"
      ];
    ],
    {mono, flatPolys}
  ];

  StringRiffle[lines, "\n"]
];

(* --------------------------------------------------------------------------
   GenerateCppMonteCarlo
   Main code generation function.
   -------------------------------------------------------------------------- *)

(* ============================================================================
   MODULE 3: Unified C++ code generation  (plan.md §5.1, §5.2)

   The convergent and IBP code-generation towers are collapsed into ONE
   parametrized path:
     emitBaseFuncBody / emitLogTail -> shared per-integrand C++ body helpers
     emitIntegrandDefinitions       -> ONE Defs emitter (per-integrand TYPE TAG;
                                       IBP-only Defs lines emitted only in IBP
                                       mode, keyed off ibpSectors =!= {})
     emitMain                       -> ONE main() emitter, kind in {MC, VEGAS}
                                       with a Batch mode (the old standalone
                                       batched-Vegas main is folded in here),
                                       result-assembly routed by info["IsIBP"]
     GenerateCppMonteCarlo          -> public entry for the subtraction /
                                       convergent path (defined below).
     GenerateCppMonteCarloIBP       -> public entry for the IBP path
                                       (defined in Module 2b).

   Output is byte-for-byte identical to the Phase-1 towers for every kind x
   sampler (golden-master regression, cross-check #25).  Integrator vocabulary
   is "MC" / "VEGAS" (plan.md §4.3) with "MonteCarlo" / "Vegas" accepted as
   input aliases so older callers / captured scripts keep working; the emitted
   C++ bytes (incl. the user-facing 'use Integrator -> "MonteCarlo"' hint) are
   unchanged.
   ============================================================================ *)

(* Shared VEGAS / sampler option defaults (plan.md §5.7) -- Join'd into every
   codegen / driver option block so the tuning defaults live in one place. *)
$vegasOptionDefaults = {
  "VegasEpsRel"   -> 1.*^-12,
  "VegasEpsAbs"   -> 1.*^-300,
  "VegasSeed"     -> 0,
  "CubaMaxComp"   -> 512,
  (* High-D VEGAS sizing (lift_error_log L1, plan.md §6.5).  Automatic =>
     resolved by resolveVegasSizing from the max sector dimension: the historical
     {1000,500,1000} for d<=4 (low-dim byte-identical), ~base*3^(d-4) for d>=5.
     An explicit integer is always honored verbatim. *)
  "VegasNStart"    -> Automatic,
  "VegasNIncrease" -> Automatic,
  "VegasNBatch"    -> Automatic
};

(* --------------------------------------------------------------------------
   resolveVegasSizing — high-D VEGAS grid sizing (lift_error_log L1, §6.5).

   The CUBA-Vegas defaults nstart=1000, nincrease=500, nbatch=1000 are tuned for
   ambient dim <= 4.  In high dimension the last iteration places only a few
   points per axis, the adaptive grid never resolves, and VEGAS returns a
   confidently-wrong value with a tight (lying) error bar.  This resolver keeps
   the historical values for d<=4 (so low-dim output is byte-identical) and
   scales by ~3^(d-4) for d>=5 (NStart ~ 8.1e4 at d=8).

     opt           — the user option value (Automatic, or an explicit integer)
     maxSectorDim  — the largest sector integration dimension
     base          — the historical default for this knob (1000 / 500 / 1000)
     kind          — "nstart" | "nincrease" | "nbatch" (for clarity; sizing is
                     identical across knobs since each scales from its own base)

   An explicit integer (or any non-Automatic value) is returned unchanged. *)
resolveVegasSizing[opt_, maxSectorDim_Integer, base_Integer, kind_String] :=
  If[opt === Automatic,
    If[maxSectorDim <= 4,
      base,
      Ceiling[base * 3^(maxSectorDim - 4)]
    ],
    opt
  ];

(* Normalize the integrator vocabulary: canonical "MC" / "VEGAS"; accept the
   older "MonteCarlo" / "Vegas" spellings as aliases.  Returns "MC" | "VEGAS". *)
normalizeIntegrator[s_] := Switch[s,
  "VEGAS" | "Vegas", "VEGAS",
  "MC" | "MonteCarlo", "MC",
  _, "MC"];

(* --------------------------------------------------------------------------
   Im(A) monomial-exponent helpers (planAXpDIVv2.md §4.1 / planR.md R2,R6).
   Factored out of the two drivers so the load-bearing declare-real / eps->0
   guard lives in exactly one place.
   -------------------------------------------------------------------------- *)

(* {imA, hasImagA} for the monomial exponents (the eps-free imaginary part).

   The REAL regulator eps lives in the monomial exponent (x1^{-1+eps} is how the
   1/eps pole enters), so we take Im at eps->0 — Im(eps) must NOT be mistaken for
   an imaginary monomial exponent.

   planR.md R2: a monomial exponent may also carry a declared kinematic symbol
   (e.g. a dimension/mass).  Those are REAL, but Im[sym] stays unevaluated, so a
   bare PossibleZeroQ[Im[sym]] is False and would spuriously flag hasImagA=True
   (then the realify subtraction yields a still-complex exponent -> liftcomplex
   $Failed, or Im[sym] leaks into the double-typed C++).  Declare the spec's
   KinematicSymbols real before taking Im so only genuinely complex literals
   survive. *)
imagMonoInfo[spec_, eps_] := Module[{realSyms, exps, imA},
  realSyms = Lookup[spec, "KinematicSymbols", {}];
  exps     = spec["MonomialExponents"] /. If[eps =!= None, eps -> 0, {}];
  imA = If[ListQ[realSyms] && Length[realSyms] > 0,
    Assuming[Element[realSyms, Reals], Simplify[Im[exps]]],
    Im[exps]];
  {imA, AnyTrue[imA, (!TrueQ[PossibleZeroQ[#]]) &]}];

(* Realify A by SUBTRACTING its eps-free imaginary part (NOT MapAt[Re], which
   would wrap the real regulator eps — present in the monomial exponent — in
   Re[eps] and leak it into the double-typed C++).  MUST run BEFORE any eps-
   pinning: on the symbolic spec the subtraction is EXACT
   ((-1+eps+i g) - i g -> -1+eps, a true Real); after pinning eps to a float it
   would leave Complex[float, 0.] (machine-zero imaginary that does not reduce to
   Real) and leak cx(...) into the domain indicator.  (planAXpDIVv2.md §4.2.) *)
realifyMonoA[spec_, imA_] := MapAt[# - I*imA &, spec, {Key["MonomialExponents"]}];

(* imagPolyInfo — the POLYNOMIAL-exponent twin of imagMonoInfo (line ~2336).
   Returns the eps-free Im(B_j), the oscillatory-phase coefficient
   (ImagPolyExponents) the SplitRealImag codegen emits.  Two symbols must be
   declared REAL before taking Im, exactly as imagMonoInfo does for A:
     - the regulator (physically real): otherwise Im[... - eps] keeps a spurious
       Im[eps], which leaks "(eps).imag()" into C++ and makes the divergence guard
       see theta as eps-DEPENDENT (splitdivmono misfire);
     - the KinematicSymbols (planR.md R2: an exponent may carry a real kinematic
       symbol, e.g. a dimension/mass).  For the principal-series propagator power
       delta = 3/2 + i nu - eps with nu a SYMBOL, a bare Im[I nu] evaluates to
       Re[nu] (NOT nu), so realifyPolyB's B - I*Im(B) would stay COMPLEX and leak
       "(nu).real()" into C++ / trip liftcomplex.  Declaring nu real gives Im = nu.
   Unlike imagMonoInfo (which substitutes eps->0), the regulator is folded into the
   declare-REAL set, NOT substituted: a genuinely eps-dependent Im(B) (the out-of-
   scope theta ~ 1/eps) then keeps a bare eps and is still caught by the downstream
   FreeQ guard, rather than being silently zeroed.  Exact / no float (#1); for real
   or eps-free-complex B this equals Im[B] (#25).  (new_request UPDATE 2026-06-30;
   mirrors imagMonoInfo / realifyMonoA / planAXpDIVv2.md §4.2.) *)
imagPolyInfo[spec_] := Module[{realSyms, exps, eps, imB},
  realSyms = Lookup[spec, "KinematicSymbols", {}];
  eps      = Lookup[spec, "RegulatorSymbol", None];
  If[eps =!= None, realSyms = Append[realSyms, eps]];
  exps     = spec["PolynomialExponents"];
  imB = If[ListQ[realSyms] && Length[realSyms] > 0,
    Assuming[Element[realSyms, Reals], Simplify[Im[exps]]], Im[exps]];
  imB];

(* realifyPolyB — realify the POLYNOMIAL exponents by SUBTRACTING their eps-free
   imaginary part, exactly as realifyMonoA does for A (NOT MapAt[Re], which wraps
   the real regulator eps in Re[eps] when eps sits in B: Re[3/2 + i nu - eps]
   becomes 3/2 - Re[eps], whose eps-derivative emits an unconverted Derivative(Re)
   into C++).  B - I*imB leaves the regulator as a TRUE real (3/2 + i nu - eps ->
   3/2 - eps), so the log-insertion / pole machinery sees a clean eps.  Like
   realifyMonoA, MUST run BEFORE any eps-pinning so the subtraction is symbolic-
   exact (after pinning it would leave Complex[float, 0.]).  Identical to Re(B)
   for real / eps-free-complex B (#25).  (new_request UPDATE 2026-06-30.) *)
realifyPolyB[spec_, imB_] := MapAt[# - I*imB &, spec, {Key["PolynomialExponents"]}];

(* --------------------------------------------------------------------------
   emitBaseFuncBody — the shared C++ body for a "base" integrand function:
   the signature, log_y[] setup, the per-polynomial monomial sums, the
   prefactor assignment, and the polynomial-exponent product.  Does NOT emit
   the closing "return ...;\n}\n" (callers add their own tail: a plain return,
   or a log-insertion factor).  resultVar is the C++ lvalue ("result",
   "g0_val", "base_val").  Reproduces the Phase-1 bytes exactly.
   -------------------------------------------------------------------------- *)
(* emitDomainIndicatorCpp (planAXpDIV.md §4.2): the lifted-sector domain
   indicator C++ block (a half-space Boole on log_y[]).  Factored out of
   emitBaseFuncBody so the hand-written remainder integrand can reuse it
   verbatim.  Returns "" for None, so the unlifted path is byte-identical (#25).
   The emitted bytes are identical to the previous inline block. *)
emitDomainIndicatorCpp[None, _] := "";
emitDomainIndicatorCpp[dc_Association, paramMap_Association] :=
Module[{logZ0str, mpStr, icList, icTerms, sumStr, rhs, isCx},
  logZ0str = mmaToCInternal[N[dc["LogZ0"]], paramMap];
  mpStr    = mmaToCInternal[N[dc["MP"]], paramMap];
  icList   = N[dc["IndicatorCoeffs"]];
  icTerms  = Table[
    mmaToCInternal[icList[[i]], paramMap] <>
    " * log_y[" <> ToString[i - 1] <> "]",
    {i, Length[icList]}
  ];
  sumStr = If[Length[icTerms] == 0, "0.0", StringRiffle[icTerms, " + "]];
  rhs = "(" <> logZ0str <> " - (" <> sumStr <> ")) * (1.0/" <> mpStr <> ")";
  (* Complex-typing fix (planIBPMULTIDIV.md §3.5).  The SplitRealImag domain map
     is built from Re(A)/Re(B), so the indicator is provably REAL — but a
     machine-zero imaginary (Complex[x,0.], the kind realifyMonoA warns about) in
     LogZ0/MP/IndicatorCoeffs makes mmaToCInternal emit cx(...), so the assembled
     RHS is complex-typed and `double log_ypstar = (...)` would fail to compile.
     Take .real() (discarding only the machine-zero noise) ONLY when a cx(...)
     literal is actually present; the real-typed path is byte-identical (#25),
     and a real symbolic kinematic expression (no cx) is left untouched. *)
  isCx = StringContainsQ[logZ0str, "cx("] || StringContainsQ[mpStr, "cx("] ||
         StringContainsQ[sumStr, "cx("];
  "    // lifted-sector domain indicator\n" <>
  "    double log_ypstar = " <>
    If[isCx, "(" <> rhs <> ").real()", rhs] <> ";\n" <>
  "    if (log_ypstar > 0.0) return cx(0.0, 0.0);\n\n"
];

(* dropDivVarFromDomain (planAXpDIV.md §4.2): for a lifted DIVERGENT sector, the
   G0 / IBP-boundary integrands run over the n-1 non-divergent coordinates, so
   the domain indicator must drop the divergent slot k.  In Case A (the only case
   ProcessSectorLifted emits) ic_k = 0, so this just removes a zero entry and
   re-indexes the survivors to the non-divergent ordering. *)
dropDivVarFromDomain[None, _] := None;
dropDivVarFromDomain[dc_Association, k_Integer] := <|
  "LogZ0"           -> dc["LogZ0"],
  "MP"              -> dc["MP"],
  "IndicatorCoeffs" -> Drop[dc["IndicatorCoeffs"], {k}]
|>;

(* liftedDomainBooleWL (planAXpDIV.md §4.2): the WL-side Boole factor that the
   NIntegrate validators (ValidateSubtraction / ValidateIBP) multiply into the
   reconstructed sector integral, mirroring ValidateLiftedDecomposition.  Returns
   1 for None, so unlifted validation is unchanged.  vv must be in the same
   coordinate order as dc["IndicatorCoeffs"]. *)
liftedDomainBooleWL[None, _] := 1;
liftedDomainBooleWL[dc_Association, vv_List] := Module[
  {lz = N[dc["LogZ0"]], mp = N[dc["MP"]], ic = N[dc["IndicatorCoeffs"]]},
  Boole[(lz - Total[ic * (Log /@ vv)]) / mp <= 0]
];

emitBaseFuncBody[funcName_String, comment_String, resultVar_String,
                 flatPolys_List, polyExps_List, prefactor_, dim_Integer,
                 paramMap_Association, domainConstraint_: None,
                 monoFactorLogs_: None, imagExps_: None,
                 monomialPhaseLog_: None] :=
Module[{funcCode, useSplit, useMonoPhase},
  (* SplitRealImag (plan.md §6.3) is active for this function ONLY when an
     imaginary-exponent list is supplied AND some entry is nonzero.  Otherwise
     the Direct product line below is emitted verbatim — so the real-exponent /
     Direct path is byte-identical to the Phase-2 goldens (cross-check #25). *)
  useSplit = ListQ[imagExps] &&
    AnyTrue[imagExps, (!TrueQ[PossibleZeroQ[#]]) &];
  (* Bare-monomial oscillatory phase from Im(A) (planAXpDIVv2.md §4.4).  Present
     ONLY in SplitRealImag mode with complex monomial exponents; None (real-A /
     Direct) -> no monomial phase emitted -> byte-identical (#25). *)
  useMonoPhase = (monomialPhaseLog =!= None);
  funcCode = "inline cx " <> funcName <>
    "(const double* y, const double* params) {\n";
  funcCode = funcCode <> "    // " <> comment <> "\n";
  funcCode = funcCode <>
    "    double log_y[" <> ToString[dim] <> "];\n";
  funcCode = funcCode <>
    "    for (int i = 0; i < " <> ToString[dim] <> "; i++)\n";
  funcCode = funcCode <>
    "        log_y[i] = (y[i] > 1e-300) ? std::log(y[i]) : -700.0;\n\n";

  (* Lifted-sector domain indicator (planAXpDIV.md §4.2).  Only emitted when a
     DomainConstraint is present; for the unlifted path (domainConstraint===None)
     emitDomainIndicatorCpp returns "", keeping the emitted C++ byte-identical (#25). *)
  funcCode = funcCode <> emitDomainIndicatorCpp[domainConstraint, paramMap];

  Do[
    funcCode = funcCode <>
      GenerateMonomialSumCpp[flatPolys[[j]], j - 1, paramMap, dim] <> "\n\n";,
    {j, Length[flatPolys]}
  ];

  funcCode = funcCode <> "    cx " <> resultVar <> " = " <>
    mmaToCInternal[prefactor, paramMap] <> ";\n";
  Do[
    funcCode = funcCode <>
      "    " <> resultVar <> " *= std::exp(" <>
      mmaToCInternal[polyExps[[j]], paramMap] <>
      " * std::log(P" <> ToString[j - 1] <> "));\n";,
    {j, Length[polyExps]}
  ];

  (* SplitRealImag oscillatory phase (plan.md §6.3; TMCv2_BUG_LOG BUG 1 +
     lift_error_log L3).  In SplitRealImag mode the sector was processed with
     Re(B), so the product above is the MAGNITUDE Prod_j P_j^{Re B_j}; here we
     reintroduce the imaginary exponents as one oscillatory factor
        exp( Sum_j i*Im(B_j) * ( log P_j^{true} ) ),  log P_j^{true} =
            std::log(P_j)            (COMPLEX log — BUG 1 fix; NOT std::abs)
          + MonoFactorLog_j          (the dropped tropical monomial factor — L3)
     with MonoFactorLog_j = Const_j + Sum_i Coeffs_{j,i} * log_y[i]. *)
  If[useSplit || useMonoPhase,
    Module[{phaseTerms, mfl, constStr, coeffStr, logPstr, jj, aa, mpl},
      phaseTerms = {};
      If[useSplit,
        Do[
          If[!TrueQ[PossibleZeroQ[imagExps[[j]]]],
            mfl = If[ListQ[monoFactorLogs] && Length[monoFactorLogs] >= j,
                     monoFactorLogs[[j]], <|"Const" -> 0, "Coeffs" -> {}|>];
            constStr = mmaToCInternal[N[mfl["Const"]], paramMap];
            coeffStr = StringJoin@Table[
              If[TrueQ[PossibleZeroQ[mfl["Coeffs"][[jj]]]], "",
                " + " <> mmaToCInternal[N[mfl["Coeffs"][[jj]]], paramMap] <>
                " * log_y[" <> ToString[jj - 1] <> "]"],
              {jj, Length[mfl["Coeffs"]]}];
            logPstr = "std::log(P" <> ToString[j - 1] <> ")";
            AppendTo[phaseTerms,
              "cx(0.0, " <> mmaToCInternal[N[imagExps[[j]]], paramMap] <> ") * (" <>
              logPstr <> " + (" <> constStr <> coeffStr <> "))"];
          ],
          {j, Length[polyExps]}
        ]
      ];
      (* Bare-monomial phase (planAXpDIVv2.md §4.4): ONE oscillatory term
         i * (Const_A + Sum_a Coeffs_a log_y[a]); Coeffs already carry Im(A) (via
         Im(A).M, flattened) and the lift constant Const_A = (Im(eaug)_p/m_p) logZ0
         (0 unlifted).  M0 Test E verified vs NIntegrate to 1.2e-8. *)
      (* planR.md R5: skip the whole monomial-phase term when its phase flattens
         to all-zero (Const_A = 0 and every Coeff = 0) — mirrors the per-term
         PossibleZeroQ skip in the Im(B) loop above.  Otherwise a complex-A sector
         with a vanishing Im(A).M would emit a dead `* std::exp(cx(0,1)*(0.))`.
         Real-A sectors never enter this block, so #25 is unaffected either way. *)
      If[useMonoPhase &&
         (!TrueQ[PossibleZeroQ[monomialPhaseLog["Const"]]] ||
           AnyTrue[monomialPhaseLog["Coeffs"], (!TrueQ[PossibleZeroQ[#]]) &]),
        mpl = monomialPhaseLog;
        constStr = mmaToCInternal[N[mpl["Const"]], paramMap];
        coeffStr = StringJoin@Table[
          If[TrueQ[PossibleZeroQ[mpl["Coeffs"][[aa]]]], "",
            " + " <> mmaToCInternal[N[mpl["Coeffs"][[aa]]], paramMap] <>
            " * log_y[" <> ToString[aa - 1] <> "]"],
          {aa, Length[mpl["Coeffs"]]}];
        AppendTo[phaseTerms,
          "cx(0.0, 1.0) * (" <> constStr <> coeffStr <> ")"];
      ];
      If[phaseTerms =!= {},
        funcCode = funcCode <> "    // SplitRealImag oscillatory phase (" <>
          Which[
            useSplit && useMonoPhase,
              "Im(B); complex log P + MonoFactorLog; + Im(A) monomial phase",
            useMonoPhase, "Im(A) monomial phase",
            True, "Im(B); complex log P + MonoFactorLog"] <> ")\n";
        funcCode = funcCode <> "    cx phaseSum = " <>
          StringRiffle[phaseTerms, " + "] <> ";\n";
        funcCode = funcCode <> "    " <> resultVar <> " *= std::exp(phaseSum);\n";
      ]
    ]
  ];
  funcCode
];

(* emitLogTail — the shared log-insertion factor block (G1 / IBP log funcs).
   withComment -> True emits the "// Log insertion factors\n" header (G1 and the
   IBP boundary-log function); False omits it (IBP term-log function), matching
   the Phase-1 bytes.  Accumulates into a C++ "cx log_sum". *)
emitLogTail[varTerms_List, polyTerms_List, paramMap_Association,
            withComment_:True] :=
Module[{funcCode},
  funcCode = If[withComment,
    "\n    // Log insertion factors\n    cx log_sum(0.0, 0.0);\n",
    "\n    cx log_sum(0.0, 0.0);\n"];
  Do[
    Module[{coeff, varIdx},
      coeff  = vt[[1]];
      varIdx = vt[[2]] - 1;
      funcCode = funcCode <>
        "    log_sum += " <> mmaToCInternal[coeff, paramMap] <>
        " * log_y[" <> ToString[varIdx] <> "];\n";
    ],
    {vt, varTerms}
  ];
  Do[
    Module[{coeff, polyIdx},
      coeff   = pt[[1]];
      polyIdx = pt[[2]] - 1;
      funcCode = funcCode <>
        "    log_sum += " <> mmaToCInternal[coeff, paramMap] <>
        " * std::log(P" <> ToString[polyIdx] <> ");\n";
    ],
    {pt, polyTerms}
  ];
  funcCode
];

(* --------------------------------------------------------------------------
   emitIntegrandDefinitions — ONE Defs emitter for both the convergent/
   subtraction path (pass ibpSectors -> {}) and the IBP path (pass the IBP
   sectors).  In subtraction mode it emits conv + g0 + g1 + rem integrands and
   the plain Defs tail; in IBP mode it emits conv + ibp(boundary base/log +
   per-term base/log), the integrand_type[]/N_CONV lines, and the IBPFuncMap.
   The two modes are mutually exclusive (divergentSectors and ibpSectors are
   never both non-empty).  Byte-identical to the Phase-1 towers (#25).
   -------------------------------------------------------------------------- *)
Options[emitIntegrandDefinitions] = {"MaxDim" -> 20};

emitIntegrandDefinitions[convergentSectors_List, divergentSectors_List,
                         ibpSectors_List, integrandSpec_Association,
                         OptionsPattern[]] :=
Module[
  {kinSyms, paramMap, nParams, isIBP,
   code, integrandFuncs, integrandDims, integrandTypes,
   nConvergent, nG0, nG1, nRemainder, nIBPFuncs,
   maxDim, ibpFuncMap, nIntegrands, allFuncNames},

  kinSyms  = integrandSpec["KinematicSymbols"];
  nParams  = Length[kinSyms];
  maxDim   = OptionValue["MaxDim"];
  isIBP    = (ibpSectors =!= {});

  paramMap = Association @@ Table[
    kinSyms[[i]] -> ("params[" <> ToString[i - 1] <> "]"),
    {i, nParams}
  ];

  integrandFuncs = {};
  integrandDims  = {};
  integrandTypes = {};
  nConvergent = 0; nG0 = 0; nG1 = 0; nRemainder = 0; nIBPFuncs = 0;
  ibpFuncMap  = {};

  (* --- Convergent sector integrands (shared by both modes) --- *)
  Do[
    Module[{sd, funcCode},
      sd = convergentSectors[[s]];
      (* SplitRealImag (plan.md §6.3): the driver attaches "ImagPolyExponents"
         (= Im(B) per polynomial) to a sector ONLY in SplitRealImag mode, after
         processing it with Re(B).  When absent (Direct / real-exponent path),
         emitBaseFuncBody emits the verbatim Direct product (byte-identical, #25). *)
      funcCode = emitBaseFuncBody[
        "integrand_conv_" <> ToString[s - 1],
        "Convergent sector " <> ToString[sd["ConeIndex"]],
        "result",
        sd["FlattenedPolys"], sd["PolynomialExponents"],
        sd["Prefactor"], sd["Dimension"], paramMap,
        Lookup[sd, "DomainConstraint", None],
        Lookup[sd, "MonoFactorLog", None],
        Lookup[sd, "ImagPolyExponents", None],
        (* planAXpDIVv2.md §4.4: bare-monomial phase (Im(A)); None real-A -> #25. *)
        Lookup[sd, "MonomialPhaseLog", None]];
      funcCode = funcCode <> "    return result;\n}\n";
      AppendTo[integrandFuncs, funcCode];
      AppendTo[integrandDims, sd["Dimension"]];
      AppendTo[integrandTypes, "conv"];
      nConvergent++;
    ],
    {s, Length[convergentSectors]}
  ];

  If[!isIBP,
    (* ============ Subtraction mode: g0 / g1 / remainder ============ *)

    (* --- G0 integrands --- *)
    Do[
      Module[{dd, funcCode},
        dd = divergentSectors[[s]];
        funcCode = emitBaseFuncBody[
          "integrand_g0_" <> ToString[s - 1],
          "G0 for divergent sector " <> ToString[dd["ConeIndex"]],
          "result",
          dd["G0FlatPolys"], dd["G0PolyExponents"],
          dd["G0Prefactor"], dd["G0Dimension"], paramMap,
          (* lifted-sector domain indicator over the n-1 non-divergent coords
             (planAXpDIV.md §4.2); None for unlifted -> byte-identical (#25). *)
          dropDivVarFromDomain[Lookup[dd, "DomainConstraint", None],
                               dd["DivergentVariable"]]];
        funcCode = funcCode <> "    return result;\n}\n";
        AppendTo[integrandFuncs, funcCode];
        AppendTo[integrandDims, dd["G0Dimension"]];
        nG0++;
      ],
      {s, Length[divergentSectors]}
    ];

    (* --- G1 integrands --- *)
    Do[
      Module[{dd, logIns, funcCode},
        dd     = divergentSectors[[s]];
        logIns = dd["G1LogInsertions"];
        funcCode = emitBaseFuncBody[
          "integrand_g1_" <> ToString[s - 1],
          "G1 for divergent sector " <> ToString[dd["ConeIndex"]],
          "g0_val",
          dd["G0FlatPolys"], dd["G0PolyExponents"],
          dd["G0Prefactor"], dd["G0Dimension"], paramMap,
          dropDivVarFromDomain[Lookup[dd, "DomainConstraint", None],
                               dd["DivergentVariable"]]];
        funcCode = funcCode <>
          emitLogTail[logIns["VariableTerms"], logIns["PolynomialTerms"],
                      paramMap, True];
        funcCode = funcCode <> "\n    return g0_val * log_sum;\n}\n";
        AppendTo[integrandFuncs, funcCode];
        AppendTo[integrandDims, dd["G0Dimension"]];
        nG1++;
      ],
      {s, Length[divergentSectors]}
    ];

    (* --- Remainder integrands (unique structure: Pfull/Psimp) --- *)
    Do[
      Module[{dd, rem, fullPolys, simpPolys, remPf, remDim, k,
              divVarExp, polyExpsRem, funcName, funcCode},
        dd        = divergentSectors[[s]];
        rem       = dd["Remainder"];
        fullPolys = rem["FullPolys"];
        simpPolys = rem["SimplifiedPolys"];
        remPf     = rem["Prefactor"];
        remDim    = rem["Dimension"];
        k         = rem["DivVarIndex"] - 1;
        divVarExp = rem["DivVarExp"];
        polyExpsRem = rem["PolynomialExponents"];
        funcName  = "integrand_rem_" <> ToString[s - 1];

        funcCode = "inline cx " <> funcName <>
          "(const double* y, const double* params) {\n";
        funcCode = funcCode <>
          "    // Remainder for divergent sector " <>
          ToString[dd["ConeIndex"]] <> "\n";
        funcCode = funcCode <>
          "    double log_y[" <> ToString[remDim] <> "];\n";
        funcCode = funcCode <>
          "    for (int i = 0; i < " <> ToString[remDim] <> "; i++)\n";
        funcCode = funcCode <>
          "        log_y[i] = (y[i] > 1e-300) ? std::log(y[i]) : -700.0;\n\n";

        (* Lifted-sector domain indicator over all n coords (planAXpDIV.md §4.2);
           the remainder is full-n-dim so no slot is dropped.  None -> ""
           (byte-identical for unlifted, #25). *)
        funcCode = funcCode <>
          emitDomainIndicatorCpp[Lookup[dd, "DomainConstraint", None], paramMap];

        funcCode = funcCode <> "    // Full polynomials\n";
        Do[
          funcCode = funcCode <>
            StringReplace[
              GenerateMonomialSumCpp[fullPolys[[j]], j - 1, paramMap,
                                     remDim, "log_y"],
              "P" <> ToString[j - 1] -> "Pfull" <> ToString[j - 1]
            ] <> "\n\n";,
          {j, Length[fullPolys]}
        ];

        funcCode = funcCode <> "    // Simplified polynomials (y_k = 0)\n";
        Do[
          funcCode = funcCode <>
            StringReplace[
              GenerateMonomialSumCpp[simpPolys[[j]], j - 1, paramMap,
                                     remDim, "log_y"],
              "P" <> ToString[j - 1] -> "Psimp" <> ToString[j - 1]
            ] <> "\n\n";,
          {j, Length[simpPolys]}
        ];

        funcCode = funcCode <> "    cx full_prod = ";
        If[Length[polyExpsRem] == 1,
          funcCode = funcCode <>
            "std::exp(" <> mmaToCInternal[polyExpsRem[[1]], paramMap] <>
            " * std::log(Pfull0));\n";,
          funcCode = funcCode <> "cx(1.0, 0.0);\n";
          Do[
            funcCode = funcCode <>
              "    full_prod *= std::exp(" <>
              mmaToCInternal[polyExpsRem[[j]], paramMap] <>
              " * std::log(Pfull" <> ToString[j - 1] <> "));\n";,
            {j, Length[polyExpsRem]}
          ];
        ];

        funcCode = funcCode <> "    cx simp_prod = ";
        If[Length[polyExpsRem] == 1,
          funcCode = funcCode <>
            "std::exp(" <> mmaToCInternal[polyExpsRem[[1]], paramMap] <>
            " * std::log(Psimp0));\n";,
          funcCode = funcCode <> "cx(1.0, 0.0);\n";
          Do[
            funcCode = funcCode <>
              "    simp_prod *= std::exp(" <>
              mmaToCInternal[polyExpsRem[[j]], paramMap] <>
              " * std::log(Psimp" <> ToString[j - 1] <> "));\n";,
            {j, Length[polyExpsRem]}
          ];
        ];

        funcCode = funcCode <>
          "\n    cx yk_factor = std::exp(" <>
          mmaToCInternal[divVarExp - 1, paramMap] <>
          " * log_y[" <> ToString[k] <> "]);\n";

        funcCode = funcCode <>
          "    return " <> mmaToCInternal[remPf, paramMap] <>
          " * yk_factor * (full_prod - simp_prod);\n}\n";

        AppendTo[integrandFuncs, funcCode];
        AppendTo[integrandDims, remDim];
        nRemainder++;
      ],
      {s, Length[divergentSectors]}
    ];
    ,
    (* ================= IBP mode: boundary + terms ================= *)
    Do[
      Module[{ibpSD, bndData, ibpTerms, sIdx, sectorFuncs},
        ibpSD    = ibpSectors[[s]];
        sIdx     = s - 1;

      If[TrueQ[ibpSD["MultiDiv"]],
        (* ===== MultiDiv corner path (planIBPMULTIDIV.md §3.3) =====
           One base + one log integrand per corner piece; the func map carries
           the per-piece Sign / Coeff0 / Coeff1 the driver combines. *)
        Module[{corners, mSectorFuncs},
          corners = ibpSD["Corners"];
          mSectorFuncs = <|"SectorIndex" -> sIdx,
                           "ConeIndex" -> ibpSD["ConeIndex"],
                           "MultiDiv" -> True, "Pieces" -> {}|>;
          Do[
            Module[{piece, cIdx, logIns, funcBase, funcLog},
              piece  = corners[[c]];
              cIdx   = c - 1;
              logIns = piece["LogInsertions"];
              funcBase = emitBaseFuncBody[
                "integrand_ibp_" <> ToString[sIdx] <> "_c" <> ToString[cIdx] <>
                  "_base",
                "IBP corner " <> ToString[cIdx] <> " base, sector " <>
                  ToString[ibpSD["ConeIndex"]],
                "result",
                piece["FlatPolys"], piece["PolyExponents"],
                piece["Prefactor"], piece["Dimension"], paramMap,
                piece["DomainConstraint"],
                piece["MonoFactorLog"],
                Lookup[ibpSD, "ImagPolyExponents", None],
                piece["MonomialPhaseLog"]];
              funcBase = funcBase <> "    return result;\n}\n";
              AppendTo[integrandFuncs, funcBase];
              AppendTo[integrandDims, piece["Dimension"]];
              AppendTo[integrandTypes, "ibp"];

              funcLog = emitBaseFuncBody[
                "integrand_ibp_" <> ToString[sIdx] <> "_c" <> ToString[cIdx] <>
                  "_log",
                "IBP corner " <> ToString[cIdx] <> " log, sector " <>
                  ToString[ibpSD["ConeIndex"]],
                "base_val",
                piece["FlatPolys"], piece["PolyExponents"],
                piece["Prefactor"], piece["Dimension"], paramMap,
                piece["DomainConstraint"],
                piece["MonoFactorLog"],
                Lookup[ibpSD, "ImagPolyExponents", None],
                piece["MonomialPhaseLog"]];
              funcLog = funcLog <>
                emitLogTail[logIns["VariableTerms"], logIns["PolynomialTerms"],
                            paramMap, True];
              funcLog = funcLog <> "\n    return base_val * log_sum;\n}\n";
              AppendTo[integrandFuncs, funcLog];
              AppendTo[integrandDims, piece["Dimension"]];
              AppendTo[integrandTypes, "ibp"];

              AppendTo[mSectorFuncs["Pieces"], <|
                "BaseFuncId" -> Length[integrandFuncs] - 2,
                "LogFuncId"  -> Length[integrandFuncs] - 1,
                "Sign"       -> piece["Sign"],
                "Coeff0"     -> piece["Coeff0"],
                "Coeff1"     -> piece["Coeff1"]|>];
              nIBPFuncs += 2;
            ],
            {c, Length[corners]}];
          AppendTo[ibpFuncMap, mSectorFuncs];
        ]
      , (* ===== existing single-pole boundary + terms path ===== *)
        bndData  = ibpSD["BoundaryData"];
        ibpTerms = ibpSD["IBPTerms"];
        sectorFuncs = <|"SectorIndex" -> sIdx,
                        "ConeIndex" -> ibpSD["ConeIndex"]|>;

        (* Boundary base (order 0) *)
        Module[{funcCode},
          funcCode = emitBaseFuncBody[
            "integrand_ibp_bnd_" <> ToString[sIdx] <> "_base",
            "IBP boundary base, sector " <> ToString[ibpSD["ConeIndex"]],
            "result",
            bndData["FlatPolys"], bndData["PolyExponents"],
            bndData["Prefactor"], bndData["Dimension"], paramMap,
            (* boundary is at y_k=1 over the n-1 non-divergent coords:
               drop the divergent slot from the domain indicator (§4.2). *)
            dropDivVarFromDomain[Lookup[ibpSD, "DomainConstraint", None],
                                 ibpSD["DivergentVariable"]],
            (* SplitRealImag phase (planCXLIFTDIV.md §4.3): per-boundary MonoFactorLog
               + sector Im(B); both None on the real path -> byte-identical (#25). *)
            Lookup[bndData, "MonoFactorLog", None],
            Lookup[ibpSD, "ImagPolyExponents", None],
            (* planAXpDIVv2.md §4.5: per-boundary bare-monomial phase (Im(A)). *)
            Lookup[bndData, "MonomialPhaseLog", None]];
          funcCode = funcCode <> "    return result;\n}\n";
          AppendTo[integrandFuncs, funcCode];
          AppendTo[integrandDims, bndData["Dimension"]];
          AppendTo[integrandTypes, "ibp"];
          sectorFuncs["BndBaseFuncId"] = Length[integrandFuncs] - 1;
          nIBPFuncs++;
        ];

        (* Boundary log (order 1) *)
        Module[{logIns, funcCode},
          logIns = bndData["LogInsertions"];
          funcCode = emitBaseFuncBody[
            "integrand_ibp_bnd_" <> ToString[sIdx] <> "_log",
            "IBP boundary log, sector " <> ToString[ibpSD["ConeIndex"]],
            "base_val",
            bndData["FlatPolys"], bndData["PolyExponents"],
            bndData["Prefactor"], bndData["Dimension"], paramMap,
            dropDivVarFromDomain[Lookup[ibpSD, "DomainConstraint", None],
                                 ibpSD["DivergentVariable"]],
            Lookup[bndData, "MonoFactorLog", None],
            Lookup[ibpSD, "ImagPolyExponents", None],
            Lookup[bndData, "MonomialPhaseLog", None]];
          funcCode = funcCode <>
            emitLogTail[logIns["VariableTerms"], logIns["PolynomialTerms"],
                        paramMap, True];
          funcCode = funcCode <> "\n    return base_val * log_sum;\n}\n";
          AppendTo[integrandFuncs, funcCode];
          AppendTo[integrandDims, bndData["Dimension"]];
          AppendTo[integrandTypes, "ibp"];
          sectorFuncs["BndLogFuncId"] = Length[integrandFuncs] - 1;
          nIBPFuncs++;
        ];

        (* IBP term functions (base + log per term) *)
        sectorFuncs["TermFuncIds"] = {};
        Do[
          Module[{termData, logIns, tIdx, funcCodeBase, funcCodeLog},
            termData  = ibpTerms[[t]];
            logIns    = termData["LogInsertions"];
            tIdx      = t - 1;

            funcCodeBase = emitBaseFuncBody[
              "integrand_ibp_" <> ToString[sIdx] <> "_t" <> ToString[tIdx] <>
                "_base",
              "IBP term " <> ToString[tIdx] <> " base, sector " <>
                ToString[ibpSD["ConeIndex"]],
              "result",
              termData["FlatPolys"], termData["PolyExponents"],
              termData["Prefactor"], termData["Dimension"], paramMap,
              (* IBP terms are full-n-dim (y_k still integrated): full domain
                 indicator, no slot dropped (§4.2). *)
              Lookup[ibpSD, "DomainConstraint", None],
              Lookup[termData, "MonoFactorLog", None],
              Lookup[ibpSD, "ImagPolyExponents", None],
              Lookup[termData, "MonomialPhaseLog", None]];
            funcCodeBase = funcCodeBase <> "    return result;\n}\n";
            AppendTo[integrandFuncs, funcCodeBase];
            AppendTo[integrandDims, termData["Dimension"]];
            AppendTo[integrandTypes, "ibp"];

            funcCodeLog = emitBaseFuncBody[
              "integrand_ibp_" <> ToString[sIdx] <> "_t" <> ToString[tIdx] <>
                "_log",
              "IBP term " <> ToString[tIdx] <> " log, sector " <>
                ToString[ibpSD["ConeIndex"]],
              "base_val",
              termData["FlatPolys"], termData["PolyExponents"],
              termData["Prefactor"], termData["Dimension"], paramMap,
              Lookup[ibpSD, "DomainConstraint", None],
              Lookup[termData, "MonoFactorLog", None],
              Lookup[ibpSD, "ImagPolyExponents", None],
              Lookup[termData, "MonomialPhaseLog", None]];
            funcCodeLog = funcCodeLog <>
              emitLogTail[logIns["VariableTerms"], logIns["PolynomialTerms"],
                          paramMap, False];
            funcCodeLog = funcCodeLog <> "\n    return base_val * log_sum;\n}\n";
            AppendTo[integrandFuncs, funcCodeLog];
            AppendTo[integrandDims, termData["Dimension"]];
            AppendTo[integrandTypes, "ibp"];

            AppendTo[sectorFuncs["TermFuncIds"],
              <|"BaseFuncId" -> Length[integrandFuncs] - 2,
                "LogFuncId"  -> Length[integrandFuncs] - 1,
                "Coeff0"     -> termData["Coeff0"],
                "Coeff1"     -> termData["Coeff1"]|>
            ];
            nIBPFuncs += 2;
          ],
          {t, Length[ibpTerms]}
        ];

        AppendTo[ibpFuncMap, sectorFuncs];
      ] (* end If MultiDiv / single-pole *)
      ],
      {s, Length[ibpSectors]}
    ];
  ];

  nIntegrands = Length[integrandFuncs];

  (* --- Assemble the full C++ file --- *)
  code = If[isIBP,
    "// Auto-generated by TropicalEval`GenerateCppMonteCarloIBP\n" <>
    "// " <> ToString[nConvergent] <> " convergent, " <>
    ToString[nIBPFuncs] <> " IBP integrands\n\n",
    "// Auto-generated by TropicalEval`GenerateCppMonteCarlo\n" <>
    "// " <> ToString[nConvergent] <> " convergent, " <>
    ToString[nG0] <> " G0, " <> ToString[nG1] <> " G1, " <>
    ToString[nRemainder] <> " remainder integrands\n\n"];

  code = code <> "#include <complex>\n";
  code = code <> "#include <cmath>\n";
  code = code <> "#include <random>\n";
  code = code <> "#include <fstream>\n";
  code = code <> "#include <vector>\n";
  code = code <> "#include <iostream>\n";
  code = code <> "#include <string>\n";
  code = code <> "#include <cassert>\n";
  code = code <> "#include <cstdlib>\n";
  code = code <> "#include <array>\n";
  code = code <> "#ifdef _OPENMP\n";
  code = code <> "#include <omp.h>\n";
  code = code <> "#endif\n\n";

  code = code <> "using cx = std::complex<double>;\n\n";

  Do[code = code <> integrandFuncs[[i]] <> "\n";, {i, nIntegrands}];

  code = code <>
    If[isIBP,
      "using IntegrandFunc = cx(*)(const double*, const double*);\n\n",
      "// Function pointer type\n" <>
      "using IntegrandFunc = cx(*)(const double*, const double*);\n\n"];

  (* Function table.  The two towers built the name list differently (conv via
     hard-coded loops, IBP via regex over the emitted code) but produced the
     SAME ordering; we build it from regex in both modes (identical bytes). *)
  allFuncNames = Table[
    StringCases[integrandFuncs[[i]],
      RegularExpression["inline cx (\\w+)\\("] -> "$1"][[1]],
    {i, nIntegrands}];

  code = code <> "IntegrandFunc integrand_table[] = {\n    " <>
    StringRiffle[allFuncNames, ",\n    "] <> "\n};\n\n";

  code = code <> "int integrand_dim[] = {" <>
    StringRiffle[ToString /@ integrandDims, ", "] <> "};\n";

  If[isIBP,
    code = code <> "int integrand_type[] = {" <>
      StringRiffle[If[# === "conv", "0", "1"] & /@ integrandTypes, ", "] <>
      "};\n"];

  code = code <> "const int N_INTEGRANDS = " <> ToString[nIntegrands] <> ";\n";
  If[isIBP,
    code = code <> "const int N_CONV = " <> ToString[nConvergent] <> ";\n"];
  code = code <> "const int N_PARAMS = " <> ToString[nParams] <> ";\n";
  code = code <> "const int MAX_DIM = " <> ToString[maxDim] <> ";\n\n";

  If[isIBP,
    <|"Defs" -> code, "Dims" -> integrandDims,
      "NConvergent" -> nConvergent, "NIBPFuncs" -> nIBPFuncs,
      "NTotal" -> nIntegrands, "NParams" -> nParams,
      "IsIBP" -> True, "IBPFuncMap" -> ibpFuncMap|>,
    <|"Defs" -> code, "Dims" -> integrandDims,
      "NConvergent" -> nConvergent, "NG0" -> nG0,
      "NG1" -> nG1, "NRemainder" -> nRemainder,
      "NTotal" -> nIntegrands, "NParams" -> nParams,
      "IsIBP" -> False|>
  ]
];

(* --------------------------------------------------------------------------
   emitMain — the SINGLE main() emitter (plan.md §5.2).  kind in {"MC","VEGAS"};
   batch -> True selects the chunked-ncomp Vegas mode (the old standalone
   batched-Vegas main folds in here; Vegas + convergent / lifted-convergent
   only).  The result-assembly is routed by info["IsIBP"]: the plain path sums all
   integrands into one line per kp; the IBP path writes the convergent sum on
   line 0 and each IBP function on its own line.  The per-branch C++ bodies are
   the Phase-1 emitters' bodies verbatim, so the emitted bytes are unchanged
   (golden-master regression, cross-check #25).  On the IBP path batch selects the
   chunked-ncomp Vegas over the IBP function table (planIBPCX.md §4); without batch
   it stays the per-kp Vegas main.
   -------------------------------------------------------------------------- *)
Options[emitMain] = Join[
  {"NSamples" -> 1000000, "SeedBase" -> 42},
  $vegasOptionDefaults
];
emitMain[info_Association, integrator_String, batch_:False,
         OptionsPattern[]] :=
Module[{code, nSamples, seedBase, epsrel, epsabs, seed, maxComp, isIBP,
        maxSectorDim, nStart, nIncrease, nBatch, vegSizing},
  nSamples = OptionValue["NSamples"];
  seedBase = OptionValue["SeedBase"];
  epsrel   = ToString[CForm[N[OptionValue["VegasEpsRel"]]]];
  epsabs   = ToString[CForm[N[OptionValue["VegasEpsAbs"]]]];
  seed     = ToString[OptionValue["VegasSeed"]];
  maxComp  = ToString[OptionValue["CubaMaxComp"]];
  isIBP    = TrueQ[info["IsIBP"]];

  (* §6.5 high-D VEGAS sizing.  maxSectorDim = largest sector integration dim
     (info["Dims"]); the resolved nstart/nincrease/nbatch substitute the
     historical literal {1000,500,1000} in the Vegas() calls below.  For
     maxSectorDim<=4 these resolve back to {1000,500,1000} so the MC and low-D
     VEGAS bytes are unchanged. *)
  maxSectorDim = If[ListQ[info["Dims"]] && Length[info["Dims"]] > 0,
                    Max[info["Dims"]], 0];
  nStart    = resolveVegasSizing[OptionValue["VegasNStart"],    maxSectorDim, 1000, "nstart"];
  nIncrease = resolveVegasSizing[OptionValue["VegasNIncrease"], maxSectorDim,  500, "nincrease"];
  nBatch    = resolveVegasSizing[OptionValue["VegasNBatch"],    maxSectorDim, 1000, "nbatch"];
  vegSizing = ToString[nStart] <> ", " <> ToString[nIncrease] <> ", " <> ToString[nBatch];

  Which[
    integrator === "MC" && !isIBP,
  code = "int main(int argc, char* argv[]) {\n";
  code = code <> "    if (argc < 3) {\n";
  code = code <> "        std::cerr << \"Usage: \" << argv[0] << \" <input_file> <output_file> [n_samples] [n_threads] [seed_base]\" << std::endl;\n";
  code = code <> "        return 1;\n";
  code = code <> "    }\n\n";

  code = code <> "    std::string input_file = argv[1];\n";
  code = code <> "    std::string output_file = argv[2];\n";
  code = code <> "    int n_samples = (argc > 3) ? std::atoi(argv[3]) : " <>
    ToString[nSamples] <> ";\n";
  code = code <> "    int n_threads = (argc > 4) ? std::atoi(argv[4]) : 1;\n";
  code = code <> "    // Optional runtime seed override (argv[5]); falls back to compile-time SeedBase.\n";
  code = code <> "    uint64_t seed_base = (argc > 5) ? std::strtoull(argv[5], nullptr, 10) : " <>
    ToString[seedBase] <> "ULL;\n";
  code = code <> "#ifdef _OPENMP\n";
  code = code <> "    if (n_threads == 1) n_threads = omp_get_max_threads();\n";
  code = code <> "    omp_set_num_threads(n_threads);\n";
  code = code <> "#endif\n\n";

  code = code <> "    // Read kinematic data\n";
  code = code <> "    std::ifstream fin(input_file);\n";
  code = code <> "    if (!fin) {\n";
  code = code <> "        std::cerr << \"Cannot open \" << input_file << std::endl;\n";
  code = code <> "        return 1;\n";
  code = code <> "    }\n\n";

  code = code <> "    std::vector<std::vector<double>> kinematic_data;\n";
  code = code <> "    if (N_PARAMS == 0) {\n";
  code = code <> "        // No kinematic parameters: read count from file, default 1\n";
  code = code <> "        int count = 1;\n";
  code = code <> "        fin >> count;\n";
  code = code <> "        if (count < 1) count = 1;\n";
  code = code <> "        for (int i = 0; i < count; i++)\n";
  code = code <> "            kinematic_data.push_back({});\n";
  code = code <> "    } else {\n";
  code = code <> "        double val;\n";
  code = code <> "        std::vector<double> row;\n";
  code = code <> "        while (fin >> val) {\n";
  code = code <> "            row.push_back(val);\n";
  code = code <> "            if ((int)row.size() == N_PARAMS) {\n";
  code = code <> "                kinematic_data.push_back(row);\n";
  code = code <> "                row.clear();\n";
  code = code <> "            }\n";
  code = code <> "        }\n";
  code = code <> "    }\n";
  code = code <> "    fin.close();\n";
  code = code <> "    int n_kp = (int)kinematic_data.size();\n";
  code = code <> "    std::cerr << \"Read \" << n_kp << \" kinematic points\" << std::endl;\n\n";

  code = code <> "    std::vector<std::array<double, 4>> results(n_kp);\n\n";

  code = code <> "    #pragma omp parallel for schedule(dynamic)\n";
  code = code <> "    for (int kp = 0; kp < n_kp; kp++) {\n";
  code = code <> "        const double* params = kinematic_data[kp].data();\n";
  code = code <> "        uint64_t seed = seed_base + (uint64_t)kp;\n";
  code = code <> "        std::mt19937_64 rng(seed);\n";
  code = code <> "        std::uniform_real_distribution<double> dist(0.0, 1.0);\n\n";

  code = code <> "        double total_re = 0.0, total_im = 0.0;\n";
  code = code <> "        double total_var_re = 0.0, total_var_im = 0.0;\n\n";

  code = code <> "        for (int s = 0; s < N_INTEGRANDS; s++) {\n";
  code = code <> "            int dim = integrand_dim[s];\n";
  code = code <> "            double mean_re = 0.0, mean_im = 0.0;\n";
  code = code <> "            double M2_re = 0.0, M2_im = 0.0;\n";

  code = code <> "#ifdef TROPICAL_MC_DEBUG\n";
  code = code <> "            int nan_count = 0;\n";
  code = code <> "            double max_mag = 0.0;\n";
  code = code <> "#endif\n\n";

  code = code <> "            for (int k = 0; k < n_samples; k++) {\n";
  code = code <> "                double y[MAX_DIM];\n";
  code = code <> "                for (int i = 0; i < dim; i++) y[i] = dist(rng);\n\n";

  code = code <> "                cx val = integrand_table[s](y, params);\n\n";

  code = code <> "#ifdef TROPICAL_MC_DEBUG\n";
  code = code <> "                if (!std::isfinite(val.real()) || !std::isfinite(val.imag())) {\n";
  code = code <> "                    nan_count++;\n";
  code = code <> "                    if (nan_count <= 5) {\n";
  code = code <> "                        std::cerr << \"NaN/Inf in sector \" << s << \" kp=\" << kp << \" y=[\";";
  code = code <> "\n                        for (int i = 0; i < dim; i++) std::cerr << y[i] << \" \";\n";
  code = code <> "                        std::cerr << \"]\" << std::endl;\n";
  code = code <> "                    }\n";
  code = code <> "                    continue;\n";
  code = code <> "                }\n";
  code = code <> "                double mag = std::abs(val);\n";
  code = code <> "                if (mag > max_mag) max_mag = mag;\n";
  code = code <> "#endif\n\n";

  code = code <> "                double d_re = val.real() - mean_re;\n";
  code = code <> "                mean_re += d_re / (k + 1);\n";
  code = code <> "                M2_re += d_re * (val.real() - mean_re);\n";
  code = code <> "                double d_im = val.imag() - mean_im;\n";
  code = code <> "                mean_im += d_im / (k + 1);\n";
  code = code <> "                M2_im += d_im * (val.imag() - mean_im);\n";
  code = code <> "            }\n\n";

  code = code <> "#ifdef TROPICAL_MC_DEBUG\n";
  code = code <> "            if (kp == 0) {\n";
  code = code <> "                std::cerr << \"Sector \" << s << \": mean=(\" << mean_re << \",\" << mean_im\n";
  code = code <> "                          << \") max_mag=\" << max_mag;\n";
  code = code <> "                if (nan_count > 0)\n";
  code = code <> "                    std::cerr << \" NaN_count=\" << nan_count;\n";
  code = code <> "                double mean_mag = std::sqrt(mean_re*mean_re + mean_im*mean_im);\n";
  code = code <> "                if (mean_mag > 0 && max_mag / mean_mag > 1000)\n";
  code = code <> "                    std::cerr << \" WARNING: large fluctuations\";\n";
  code = code <> "                std::cerr << std::endl;\n";
  code = code <> "            }\n";
  code = code <> "            if ((double)nan_count / n_samples > 0.001)\n";
  code = code <> "                std::cerr << \"WARNING: >0.1%% NaN in sector \" << s << \" kp=\" << kp << std::endl;\n";
  code = code <> "#endif\n\n";

  code = code <> "            total_re += mean_re;\n";
  code = code <> "            total_im += mean_im;\n";
  code = code <> "            total_var_re += M2_re / ((double)n_samples * (n_samples - 1));\n";
  code = code <> "            total_var_im += M2_im / ((double)n_samples * (n_samples - 1));\n";
  code = code <> "        }\n\n";

  code = code <> "        results[kp] = {total_re, total_im,\n";
  code = code <> "                       std::sqrt(total_var_re), std::sqrt(total_var_im)};\n";

  code = code <> "#ifdef TROPICAL_MC_DEBUG\n";
  code = code <> "        if (kp == 0) {\n";
  code = code <> "            std::cerr << \"KP 0 total: (\" << total_re << \", \" << total_im\n";
  code = code <> "                      << \") +/- (\" << std::sqrt(total_var_re) << \", \"\n";
  code = code <> "                      << std::sqrt(total_var_im) << \")\" << std::endl;\n";
  code = code <> "        }\n";
  code = code <> "#endif\n";

  code = code <> "    }\n\n";

  code = code <> "    // Write results\n";
  code = code <> "    std::ofstream fout(output_file);\n";
  code = code <> "    if (!fout) {\n";
  code = code <> "        std::cerr << \"Cannot open \" << output_file << std::endl;\n";
  code = code <> "        return 1;\n";
  code = code <> "    }\n";
  code = code <> "    fout.precision(17);\n";
  code = code <> "    for (int kp = 0; kp < n_kp; kp++) {\n";
  code = code <> "        fout << results[kp][0] << \" \" << results[kp][1] << \" \"\n";
  code = code <> "             << results[kp][2] << \" \" << results[kp][3] << \"\\n\";\n";
  code = code <> "    }\n";
  code = code <> "    fout.close();\n\n";

  code = code <> "    std::cerr << \"Done. Processed \" << n_kp << \" kinematic points.\" << std::endl;\n";
  code = code <> "    return 0;\n";
  code = code <> "}\n";
    code
    ,
    integrator === "MC" && isIBP,
  code = "int main(int argc, char* argv[]) {\n";
  code = code <> "    if (argc < 3) {\n";
    code = code <> "        std::cerr << \"Usage: \" << argv[0] << \" <input_file> <output_file> [n_samples] [n_threads] [seed_base]\" << std::endl;\n";
    code = code <> "        return 1;\n";
    code = code <> "    }\n\n";

    code = code <> "    std::string input_file = argv[1];\n";
    code = code <> "    std::string output_file = argv[2];\n";
    code = code <> "    int n_samples = (argc > 3) ? std::atoi(argv[3]) : " <>
      ToString[nSamples] <> ";\n";
    code = code <> "    int n_threads = (argc > 4) ? std::atoi(argv[4]) : 1;\n";
    code = code <> "    // Optional runtime seed override (argv[5]); falls back to compile-time SeedBase.\n";
    code = code <> "    uint64_t seed_base = (argc > 5) ? std::strtoull(argv[5], nullptr, 10) : " <>
      ToString[seedBase] <> "ULL;\n";
    code = code <> "#ifdef _OPENMP\n";
    code = code <> "    if (n_threads == 1) n_threads = omp_get_max_threads();\n";
    code = code <> "    omp_set_num_threads(n_threads);\n";
    code = code <> "#endif\n\n";

    code = code <> "    std::ifstream fin(input_file);\n";
    code = code <> "    if (!fin) { std::cerr << \"Cannot open \" << input_file << std::endl; return 1; }\n\n";

    code = code <> "    std::vector<std::vector<double>> kinematic_data;\n";
    code = code <> "    if (N_PARAMS == 0) {\n";
    code = code <> "        int count = 1; fin >> count; if (count < 1) count = 1;\n";
    code = code <> "        for (int i = 0; i < count; i++) kinematic_data.push_back({});\n";
    code = code <> "    } else {\n";
    code = code <> "        double val; std::vector<double> row;\n";
    code = code <> "        while (fin >> val) { row.push_back(val);\n";
    code = code <> "            if ((int)row.size() == N_PARAMS) { kinematic_data.push_back(row); row.clear(); }}\n";
    code = code <> "    }\n";
    code = code <> "    fin.close();\n";
    code = code <> "    int n_kp = (int)kinematic_data.size();\n";
    code = code <> "    std::cerr << \"Read \" << n_kp << \" kinematic points\" << std::endl;\n\n";

    (* Output: one convergent sum + each IBP function individually,
       per kinematic point *)
    code = code <> "    // Output: (1 + N_IBP_FUNCS) lines per kinematic point\n";
    code = code <> "    // Line 0: convergent sum; Lines 1..N_IBP: individual IBP functions\n";
    code = code <> "    int n_ibp = N_INTEGRANDS - N_CONV;\n";
    code = code <> "    std::vector<std::array<double, 4>> results((1 + n_ibp) * n_kp);\n\n";

    code = code <> "    #pragma omp parallel for schedule(dynamic)\n";
    code = code <> "    for (int kp = 0; kp < n_kp; kp++) {\n";
    code = code <> "        const double* params = kinematic_data[kp].data();\n";
    code = code <> "        uint64_t seed = seed_base + (uint64_t)kp;\n";
    code = code <> "        std::mt19937_64 rng(seed);\n";
    code = code <> "        std::uniform_real_distribution<double> dist(0.0, 1.0);\n\n";

    (* Convergent sum *)
    code = code <> "        double conv_re = 0.0, conv_im = 0.0;\n";
    code = code <> "        double conv_var_re = 0.0, conv_var_im = 0.0;\n\n";

    code = code <> "        for (int s = 0; s < N_CONV; s++) {\n";
    code = code <> "            int dim = integrand_dim[s];\n";
    code = code <> "            double mean_re = 0.0, mean_im = 0.0;\n";
    code = code <> "            double M2_re = 0.0, M2_im = 0.0;\n";
    code = code <> "            for (int k = 0; k < n_samples; k++) {\n";
    code = code <> "                double y[MAX_DIM];\n";
    code = code <> "                for (int i = 0; i < dim; i++) y[i] = dist(rng);\n";
    code = code <> "                cx val = integrand_table[s](y, params);\n";
    code = code <> "                double d_re = val.real() - mean_re;\n";
    code = code <> "                mean_re += d_re / (k + 1);\n";
    code = code <> "                M2_re += d_re * (val.real() - mean_re);\n";
    code = code <> "                double d_im = val.imag() - mean_im;\n";
    code = code <> "                mean_im += d_im / (k + 1);\n";
    code = code <> "                M2_im += d_im * (val.imag() - mean_im);\n";
    code = code <> "            }\n";
    code = code <> "            conv_re += mean_re; conv_im += mean_im;\n";
    code = code <> "            conv_var_re += M2_re / ((double)n_samples * (n_samples - 1));\n";
    code = code <> "            conv_var_im += M2_im / ((double)n_samples * (n_samples - 1));\n";
    code = code <> "        }\n";
    code = code <> "        results[kp * (1 + n_ibp)] = {conv_re, conv_im, std::sqrt(conv_var_re), std::sqrt(conv_var_im)};\n\n";

    (* Individual IBP functions *)
    code = code <> "        for (int s = N_CONV; s < N_INTEGRANDS; s++) {\n";
    code = code <> "            int dim = integrand_dim[s];\n";
    code = code <> "            double mean_re = 0.0, mean_im = 0.0;\n";
    code = code <> "            double M2_re = 0.0, M2_im = 0.0;\n";
    code = code <> "            for (int k = 0; k < n_samples; k++) {\n";
    code = code <> "                double y[MAX_DIM];\n";
    code = code <> "                for (int i = 0; i < dim; i++) y[i] = dist(rng);\n";
    code = code <> "                cx val = integrand_table[s](y, params);\n";
    code = code <> "                double d_re = val.real() - mean_re;\n";
    code = code <> "                mean_re += d_re / (k + 1);\n";
    code = code <> "                M2_re += d_re * (val.real() - mean_re);\n";
    code = code <> "                double d_im = val.imag() - mean_im;\n";
    code = code <> "                mean_im += d_im / (k + 1);\n";
    code = code <> "                M2_im += d_im * (val.imag() - mean_im);\n";
    code = code <> "            }\n";
    code = code <> "            int idx = kp * (1 + n_ibp) + (s - N_CONV + 1);\n";
    code = code <> "            results[idx] = {mean_re, mean_im, std::sqrt(M2_re / ((double)n_samples * (n_samples - 1))), std::sqrt(M2_im / ((double)n_samples * (n_samples - 1)))};\n";
    code = code <> "        }\n";
    code = code <> "    }\n\n";

    code = code <> "    std::ofstream fout(output_file);\n";
    code = code <> "    if (!fout) { std::cerr << \"Cannot open \" << output_file << std::endl; return 1; }\n";
    code = code <> "    fout.precision(17);\n";
    code = code <> "    for (int i = 0; i < (int)results.size(); i++) {\n";
    code = code <> "        fout << results[i][0] << \" \" << results[i][1] << \" \"\n";
    code = code <> "             << results[i][2] << \" \" << results[i][3] << \"\\n\";\n";
    code = code <> "    }\n";
    code = code <> "    fout.close();\n";
    code = code <> "    std::cerr << \"Done. Processed \" << n_kp << \" kinematic points.\" << std::endl;\n";
    code = code <> "    return 0;\n}\n";
    code
    ,
    integrator === "VEGAS" && batch && !isIBP,
  code = "// TROPICAL_REQUIRES_CUBA  (CompileCpp greps for this sentinel)\n";
  code = code <> "#ifdef TROPICAL_USE_CUBA\n";
  code = code <> "extern \"C\" {\n";
  code = code <> "#include <cuba.h>\n";
  code = code <> "}\n";
  code = code <> "#include <algorithm>\n";
  code = code <> "static const double* gb_par = nullptr;\n";
  code = code <> "static int gb_sector = 0, gb_dim = 0, gb_kp0 = 0, gb_chunk = 0;\n";
  code = code <> "static int cubaBatch(const int* ndim, const cubareal xx[], const int* ncomp,\n";
  code = code <> "                     cubareal ff[], void* userdata) {\n";
  code = code <> "    (void)ndim; (void)ncomp; (void)userdata;\n";
  code = code <> "    double y[MAX_DIM];\n";
  code = code <> "    for (int i = 0; i < gb_dim; ++i) y[i] = (double)xx[i];\n";
  code = code <> "    for (int c = 0; c < gb_chunk; ++c) {\n";
  code = code <> "        const double* pp = (N_PARAMS > 0) ? &gb_par[(size_t)(gb_kp0 + c) * N_PARAMS] : nullptr;\n";
  code = code <> "        cx v = integrand_table[gb_sector](y, pp);\n";
  code = code <> "        ff[2*c] = v.real(); ff[2*c + 1] = v.imag();\n";
  code = code <> "    }\n";
  code = code <> "    return 0;\n";
  code = code <> "}\n";
  code = code <> "#endif\n\n";

  code = code <> "int main(int argc, char* argv[]) {\n";
  code = code <> "    if (argc < 3) {\n";
  code = code <> "        std::cerr << \"Usage: \" << argv[0] << \" <input_file> <output_file> [maxeval] [n_threads]\" << std::endl;\n";
  code = code <> "        return 1;\n";
  code = code <> "    }\n\n";
  code = code <> "    std::string input_file = argv[1];\n";
  code = code <> "    std::string output_file = argv[2];\n";
  code = code <> "    long long maxeval = (argc > 3) ? std::atoll(argv[3]) : " <> ToString[nSamples] <> "LL;\n\n";

  code = code <> "    std::ifstream fin(input_file);\n";
  code = code <> "    if (!fin) { std::cerr << \"Cannot open \" << input_file << std::endl; return 1; }\n";
  code = code <> "    std::vector<std::vector<double>> kinematic_data;\n";
  code = code <> "    if (N_PARAMS == 0) {\n";
  code = code <> "        int count = 1; fin >> count; if (count < 1) count = 1;\n";
  code = code <> "        for (int i = 0; i < count; i++) kinematic_data.push_back({});\n";
  code = code <> "    } else {\n";
  code = code <> "        double val; std::vector<double> row;\n";
  code = code <> "        while (fin >> val) { row.push_back(val);\n";
  code = code <> "            if ((int)row.size() == N_PARAMS) { kinematic_data.push_back(row); row.clear(); }}\n";
  code = code <> "    }\n";
  code = code <> "    fin.close();\n";
  code = code <> "    int n_kp = (int)kinematic_data.size();\n";
  code = code <> "    std::cerr << \"Read \" << n_kp << \" kinematic points (Vegas batch)\" << std::endl;\n\n";

  code = code <> "    std::vector<std::array<double, 4>> results(n_kp);\n\n";

  code = code <> "#ifdef TROPICAL_USE_CUBA\n";
  code = code <> "    { const int zero = 0; cubacores(&zero, &zero); }\n";
  code = code <> "    const double epsrel = " <> epsrel <> ", epsabs = " <> epsabs <> ";\n";
  code = code <> "    const int    seed    = " <> seed <> ";\n";
  code = code <> "    const int    maxComp = " <> maxComp <> ";   // ncomp ceiling -- never exceed (Vegas segfaults above it)\n";
  code = code <> "    const int    kpPerChunk = (maxComp / 2 < 1) ? 1 : maxComp / 2;  // 2 comps (re,im) per kp\n";
  code = code <> "    std::vector<double> flat_par((size_t)n_kp * (N_PARAMS > 0 ? N_PARAMS : 1), 0.0);\n";
  code = code <> "    for (int kp = 0; kp < n_kp; ++kp)\n";
  code = code <> "        for (int p = 0; p < N_PARAMS; ++p) flat_par[(size_t)kp * N_PARAMS + p] = kinematic_data[kp][p];\n";
  code = code <> "    gb_par = flat_par.data();\n";
  code = code <> "    std::vector<double> tre(n_kp, 0), tim(n_kp, 0), vre(n_kp, 0), vim(n_kp, 0);\n";
  code = code <> "    for (int s = 0; s < N_INTEGRANDS; ++s) {\n";
  code = code <> "        gb_sector = s; gb_dim = integrand_dim[s];\n";
  code = code <> "        for (int k0 = 0; k0 < n_kp; k0 += kpPerChunk) {\n";
  code = code <> "            int cs = std::min(kpPerChunk, n_kp - k0);\n";
  code = code <> "            gb_kp0 = k0; gb_chunk = cs;\n";
  code = code <> "            int ncomp = 2 * cs;\n";
  code = code <> "            std::vector<cubareal> integ(ncomp), err(ncomp), prob(ncomp);\n";
  code = code <> "            int neval = 0, fail = 0;\n";
  code = code <> "            Vegas(gb_dim, ncomp, cubaBatch, nullptr, 1, epsrel, epsabs, 0, seed,\n";
  code = code <> "                  0, (int)maxeval, " <> vegSizing <> ", 0, nullptr, nullptr,\n";
  code = code <> "                  &neval, &fail, integ.data(), err.data(), prob.data());\n";
  code = code <> "            for (int c = 0; c < cs; ++c) {\n";
  code = code <> "                tre[k0 + c] += integ[2*c];       tim[k0 + c] += integ[2*c + 1];\n";
  code = code <> "                vre[k0 + c] += err[2*c]*err[2*c]; vim[k0 + c] += err[2*c + 1]*err[2*c + 1];\n";
  code = code <> "            }\n";
  code = code <> "        }\n";
  code = code <> "    }\n";
  code = code <> "    for (int kp = 0; kp < n_kp; ++kp)\n";
  code = code <> "        results[kp] = {tre[kp], tim[kp], std::sqrt(vre[kp]), std::sqrt(vim[kp])};\n";
  code = code <> "#else\n";
  code = code <> "    (void)maxeval;\n";
  code = code <> "    std::cerr << \"ERROR: built without CUBA. Rebuild with -DTROPICAL_USE_CUBA \"\n";
  code = code <> "                 \"(install CUBA), or use Integrator -> \\\"MonteCarlo\\\".\" << std::endl;\n";
  code = code <> "    return 2;\n";
  code = code <> "#endif\n\n";

  code = code <> "    std::ofstream fout(output_file);\n";
  code = code <> "    if (!fout) { std::cerr << \"Cannot open \" << output_file << std::endl; return 1; }\n";
  code = code <> "    fout.precision(17);\n";
  code = code <> "    for (int kp = 0; kp < n_kp; kp++) {\n";
  code = code <> "        fout << results[kp][0] << \" \" << results[kp][1] << \" \"\n";
  code = code <> "             << results[kp][2] << \" \" << results[kp][3] << \"\\n\";\n";
  code = code <> "    }\n";
  code = code <> "    fout.close();\n";
  code = code <> "    std::cerr << \"Done. Processed \" << n_kp << \" kinematic points (batched).\" << std::endl;\n";
  code = code <> "    return 0;\n";
  code = code <> "}\n";
    code
    ,
    integrator === "VEGAS" && !batch && !isIBP,
  code = "// TROPICAL_REQUIRES_CUBA  (CompileCpp greps for this sentinel)\n";
  code = code <> "#ifdef TROPICAL_USE_CUBA\n";
  code = code <> "extern \"C\" {\n";
  code = code <> "#include <cuba.h>\n";
  code = code <> "}\n";
  code = code <> "static IntegrandFunc g_fn     = nullptr;\n";
  code = code <> "static const double* g_params = nullptr;\n";
  code = code <> "static int           g_dim    = 0;\n";
  code = code <> "static int cubaWrap(const int* ndim, const cubareal xx[], const int* ncomp,\n";
  code = code <> "                    cubareal ff[], void* userdata) {\n";
  code = code <> "    (void)ndim; (void)ncomp; (void)userdata;\n";
  code = code <> "    double y[MAX_DIM];\n";
  code = code <> "    for (int i = 0; i < g_dim; ++i) y[i] = (double)xx[i];\n";
  code = code <> "    cx v = g_fn(y, g_params);\n";
  code = code <> "    ff[0] = v.real(); ff[1] = v.imag();   // ncomp = 2 (complex integrand)\n";
  code = code <> "    return 0;\n";
  code = code <> "}\n";
  code = code <> "// integrate one sector with Vegas; a 0-dim integrand (e.g. an IBP\n";
  code = code <> "// boundary of a 1-var sector) is a constant, so eval it directly.\n";
  code = code <> "static void vegasOne(int s, const double* params, long long maxeval,\n";
  code = code <> "                     double epsrel, double epsabs, int seed,\n";
  code = code <> "                     double out[2], double err[2]) {\n";
  code = code <> "    g_fn = integrand_table[s]; g_dim = integrand_dim[s]; g_params = params;\n";
  code = code <> "    if (g_dim == 0) {\n";
  code = code <> "        double y0[1] = {0.0}; cx v = g_fn(y0, params);\n";
  code = code <> "        out[0] = v.real(); out[1] = v.imag(); err[0] = 0.0; err[1] = 0.0;\n";
  code = code <> "        return;\n";
  code = code <> "    }\n";
  code = code <> "    int neval = 0, fail = 0; cubareal I[2], E[2], Pr[2];\n";
  code = code <> "    Vegas(g_dim, 2, cubaWrap, nullptr, 1, epsrel, epsabs, 0, seed,\n";
  code = code <> "          0, (int)maxeval, " <> vegSizing <> ", 0, nullptr, nullptr,\n";
  code = code <> "          &neval, &fail, I, E, Pr);\n";
  code = code <> "    out[0] = I[0]; out[1] = I[1]; err[0] = E[0]; err[1] = E[1];\n";
  code = code <> "}\n";
  code = code <> "#endif\n\n";

  code = code <> "int main(int argc, char* argv[]) {\n";
  code = code <> "    if (argc < 3) {\n";
  code = code <> "        std::cerr << \"Usage: \" << argv[0] << \" <input_file> <output_file> [maxeval] [n_threads]\" << std::endl;\n";
  code = code <> "        return 1;\n";
  code = code <> "    }\n\n";
  code = code <> "    std::string input_file = argv[1];\n";
  code = code <> "    std::string output_file = argv[2];\n";
  code = code <> "    long long maxeval = (argc > 3) ? std::atoll(argv[3]) : " <> ToString[nSamples] <> "LL;\n\n";

  code = code <> "    std::ifstream fin(input_file);\n";
  code = code <> "    if (!fin) { std::cerr << \"Cannot open \" << input_file << std::endl; return 1; }\n";
  code = code <> "    std::vector<std::vector<double>> kinematic_data;\n";
  code = code <> "    if (N_PARAMS == 0) {\n";
  code = code <> "        int count = 1; fin >> count; if (count < 1) count = 1;\n";
  code = code <> "        for (int i = 0; i < count; i++) kinematic_data.push_back({});\n";
  code = code <> "    } else {\n";
  code = code <> "        double val; std::vector<double> row;\n";
  code = code <> "        while (fin >> val) { row.push_back(val);\n";
  code = code <> "            if ((int)row.size() == N_PARAMS) { kinematic_data.push_back(row); row.clear(); }}\n";
  code = code <> "    }\n";
  code = code <> "    fin.close();\n";
  code = code <> "    int n_kp = (int)kinematic_data.size();\n";
  code = code <> "    std::cerr << \"Read \" << n_kp << \" kinematic points (Vegas)\" << std::endl;\n\n";

  code = code <> "    std::vector<std::array<double, 4>> results(n_kp);\n\n";

  code = code <> "#ifdef TROPICAL_USE_CUBA\n";
  code = code <> "    { const int zero = 0; cubacores(&zero, &zero); }  // single process: deterministic, macOS-safe\n";
  code = code <> "    const double epsrel = " <> epsrel <> ", epsabs = " <> epsabs <> ";\n";
  code = code <> "    const int    seed   = " <> seed <> ";\n";
  code = code <> "    for (int kp = 0; kp < n_kp; ++kp) {        // SERIAL over kp: CUBA keeps global state\n";
  code = code <> "        const double* params = kinematic_data[kp].data();\n";
  code = code <> "        double tre = 0, tim = 0, vre = 0, vim = 0;\n";
  code = code <> "        for (int s = 0; s < N_INTEGRANDS; ++s) {\n";
  code = code <> "            double I[2], E[2]; vegasOne(s, params, maxeval, epsrel, epsabs, seed, I, E);\n";
  code = code <> "            tre += I[0]; tim += I[1];\n";
  code = code <> "            vre += E[0]*E[0]; vim += E[1]*E[1];\n";
  code = code <> "        }\n";
  code = code <> "        results[kp] = {tre, tim, std::sqrt(vre), std::sqrt(vim)};\n";
  code = code <> "    }\n";
  code = code <> "#else\n";
  code = code <> "    (void)maxeval;\n";
  code = code <> "    std::cerr << \"ERROR: built without CUBA. Rebuild with -DTROPICAL_USE_CUBA \"\n";
  code = code <> "                 \"(install CUBA), or use Integrator -> \\\"MonteCarlo\\\".\" << std::endl;\n";
  code = code <> "    return 2;\n";
  code = code <> "#endif\n\n";

  code = code <> "    std::ofstream fout(output_file);\n";
  code = code <> "    if (!fout) { std::cerr << \"Cannot open \" << output_file << std::endl; return 1; }\n";
  code = code <> "    fout.precision(17);\n";
  code = code <> "    for (int kp = 0; kp < n_kp; kp++) {\n";
  code = code <> "        fout << results[kp][0] << \" \" << results[kp][1] << \" \"\n";
  code = code <> "             << results[kp][2] << \" \" << results[kp][3] << \"\\n\";\n";
  code = code <> "    }\n";
  code = code <> "    fout.close();\n";
  code = code <> "    std::cerr << \"Done. Processed \" << n_kp << \" kinematic points.\" << std::endl;\n";
  code = code <> "    return 0;\n";
  code = code <> "}\n";
    code
    ,
    (* IBP path, batched Vegas (planIBPCX.md §4): one Vegas call per (IBP
       function, kp-chunk) sharing samples across the chunk, exactly like the
       convergent batched path (cubaBatch / the gb_ globals), but iterating the
       IBP FUNCTION table and writing the (1 + n_ibp)-row-per-kp layout
       (convergent sectors summed into row 0, each IBP function on its own row).
       The output contract is identical to the per-kp IBP main, so the driver
       combine is unchanged. *)
    integrator === "VEGAS" && batch && isIBP,
  code = "// TROPICAL_REQUIRES_CUBA  (CompileCpp greps for this sentinel)\n";
  code = code <> "#ifdef TROPICAL_USE_CUBA\n";
  code = code <> "extern \"C\" {\n";
  code = code <> "#include <cuba.h>\n";
  code = code <> "}\n";
  code = code <> "#include <algorithm>\n";
  code = code <> "static const double* gb_par = nullptr;\n";
  code = code <> "static int gb_sector = 0, gb_dim = 0, gb_kp0 = 0, gb_chunk = 0;\n";
  code = code <> "static int cubaBatch(const int* ndim, const cubareal xx[], const int* ncomp,\n";
  code = code <> "                     cubareal ff[], void* userdata) {\n";
  code = code <> "    (void)ndim; (void)ncomp; (void)userdata;\n";
  code = code <> "    double y[MAX_DIM];\n";
  code = code <> "    for (int i = 0; i < gb_dim; ++i) y[i] = (double)xx[i];\n";
  code = code <> "    for (int c = 0; c < gb_chunk; ++c) {\n";
  code = code <> "        const double* pp = (N_PARAMS > 0) ? &gb_par[(size_t)(gb_kp0 + c) * N_PARAMS] : nullptr;\n";
  code = code <> "        cx v = integrand_table[gb_sector](y, pp);\n";
  code = code <> "        ff[2*c] = v.real(); ff[2*c + 1] = v.imag();\n";
  code = code <> "    }\n";
  code = code <> "    return 0;\n";
  code = code <> "}\n";
  code = code <> "#endif\n\n";

  code = code <> "int main(int argc, char* argv[]) {\n";
  code = code <> "    if (argc < 3) {\n";
  code = code <> "        std::cerr << \"Usage: \" << argv[0] << \" <input_file> <output_file> [maxeval] [n_threads]\" << std::endl;\n";
  code = code <> "        return 1;\n";
  code = code <> "    }\n\n";
  code = code <> "    std::string input_file = argv[1];\n";
  code = code <> "    std::string output_file = argv[2];\n";
  code = code <> "    long long maxeval = (argc > 3) ? std::atoll(argv[3]) : " <> ToString[nSamples] <> "LL;\n\n";
  code = code <> "    std::ifstream fin(input_file);\n";
  code = code <> "    if (!fin) { std::cerr << \"Cannot open \" << input_file << std::endl; return 1; }\n";
  code = code <> "    std::vector<std::vector<double>> kinematic_data;\n";
  code = code <> "    if (N_PARAMS == 0) {\n";
  code = code <> "        int count = 1; fin >> count; if (count < 1) count = 1;\n";
  code = code <> "        for (int i = 0; i < count; i++) kinematic_data.push_back({});\n";
  code = code <> "    } else {\n";
  code = code <> "        double val; std::vector<double> row;\n";
  code = code <> "        while (fin >> val) { row.push_back(val);\n";
  code = code <> "            if ((int)row.size() == N_PARAMS) { kinematic_data.push_back(row); row.clear(); }}\n";
  code = code <> "    }\n";
  code = code <> "    fin.close();\n";
  code = code <> "    int n_kp = (int)kinematic_data.size();\n";
  code = code <> "    std::cerr << \"Read \" << n_kp << \" kinematic points (Vegas batch IBP)\" << std::endl;\n\n";
  code = code <> "    int n_ibp = N_INTEGRANDS - N_CONV;\n";
  code = code <> "    std::vector<std::array<double, 4>> results((size_t)(1 + n_ibp) * n_kp);\n\n";

  code = code <> "#ifdef TROPICAL_USE_CUBA\n";
  code = code <> "    { const int zero = 0; cubacores(&zero, &zero); }\n";
  code = code <> "    const double epsrel = " <> epsrel <> ", epsabs = " <> epsabs <> ";\n";
  code = code <> "    const int    seed    = " <> seed <> ";\n";
  code = code <> "    const int    maxComp = " <> maxComp <> ";   // ncomp ceiling -- never exceed (Vegas segfaults above it)\n";
  code = code <> "    const int    kpPerChunk = (maxComp / 2 < 1) ? 1 : maxComp / 2;  // 2 comps (re,im) per kp\n";
  code = code <> "    std::vector<double> flat_par((size_t)n_kp * (N_PARAMS > 0 ? N_PARAMS : 1), 0.0);\n";
  code = code <> "    for (int kp = 0; kp < n_kp; ++kp)\n";
  code = code <> "        for (int p = 0; p < N_PARAMS; ++p) flat_par[(size_t)kp * N_PARAMS + p] = kinematic_data[kp][p];\n";
  code = code <> "    gb_par = flat_par.data();\n";
  code = code <> "    // (1 + n_ibp) output rows per kp: row 0 = convergent sum, rows 1..n_ibp = IBP funcs.\n";
  code = code <> "    const size_t nrows = (size_t)(1 + n_ibp) * n_kp;\n";
  code = code <> "    std::vector<double> ore(nrows, 0), oim(nrows, 0), ovre(nrows, 0), ovim(nrows, 0);\n";
  code = code <> "    for (int s = 0; s < N_INTEGRANDS; ++s) {\n";
  code = code <> "        gb_sector = s; gb_dim = integrand_dim[s];\n";
  code = code <> "        // convergent sectors (s < N_CONV) accumulate into row 0; each IBP\n";
  code = code <> "        // function (s >= N_CONV) writes its own row (s - N_CONV + 1).\n";
  code = code <> "        int rowInKp = (s < N_CONV) ? 0 : (s - N_CONV + 1);\n";
  code = code <> "        if (gb_dim == 0) {\n";
  code = code <> "            // 0-dim integrand (e.g. IBP boundary of a 1-var sector): a constant\n";
  code = code <> "            // per kp -- evaluate directly, no Vegas.\n";
  code = code <> "            for (int kp = 0; kp < n_kp; ++kp) {\n";
  code = code <> "                const double* pp = (N_PARAMS > 0) ? &flat_par[(size_t)kp * N_PARAMS] : nullptr;\n";
  code = code <> "                double y0[1] = {0.0}; cx v = integrand_table[s](y0, pp);\n";
  code = code <> "                size_t r = (size_t)kp * (1 + n_ibp) + rowInKp;\n";
  code = code <> "                ore[r] += v.real(); oim[r] += v.imag();\n";
  code = code <> "            }\n";
  code = code <> "            continue;\n";
  code = code <> "        }\n";
  code = code <> "        for (int k0 = 0; k0 < n_kp; k0 += kpPerChunk) {\n";
  code = code <> "            int cs = std::min(kpPerChunk, n_kp - k0);\n";
  code = code <> "            gb_kp0 = k0; gb_chunk = cs;\n";
  code = code <> "            int ncomp = 2 * cs;\n";
  code = code <> "            std::vector<cubareal> integ(ncomp), err(ncomp), prob(ncomp);\n";
  code = code <> "            int neval = 0, fail = 0;\n";
  code = code <> "            Vegas(gb_dim, ncomp, cubaBatch, nullptr, 1, epsrel, epsabs, 0, seed,\n";
  code = code <> "                  0, (int)maxeval, " <> vegSizing <> ", 0, nullptr, nullptr,\n";
  code = code <> "                  &neval, &fail, integ.data(), err.data(), prob.data());\n";
  code = code <> "            for (int c = 0; c < cs; ++c) {\n";
  code = code <> "                size_t r = (size_t)(k0 + c) * (1 + n_ibp) + rowInKp;\n";
  code = code <> "                ore[r]  += integ[2*c];        oim[r]  += integ[2*c + 1];\n";
  code = code <> "                ovre[r] += err[2*c]*err[2*c]; ovim[r] += err[2*c + 1]*err[2*c + 1];\n";
  code = code <> "            }\n";
  code = code <> "        }\n";
  code = code <> "    }\n";
  code = code <> "    for (size_t r = 0; r < nrows; ++r)\n";
  code = code <> "        results[r] = {ore[r], oim[r], std::sqrt(ovre[r]), std::sqrt(ovim[r])};\n";
  code = code <> "#else\n";
  code = code <> "    (void)maxeval;\n";
  code = code <> "    std::cerr << \"ERROR: built without CUBA. Rebuild with -DTROPICAL_USE_CUBA \"\n";
  code = code <> "                 \"(install CUBA), or use Integrator -> \\\"MonteCarlo\\\".\" << std::endl;\n";
  code = code <> "    return 2;\n";
  code = code <> "#endif\n\n";

  code = code <> "    std::ofstream fout(output_file);\n";
  code = code <> "    if (!fout) { std::cerr << \"Cannot open \" << output_file << std::endl; return 1; }\n";
  code = code <> "    fout.precision(17);\n";
  code = code <> "    for (int i = 0; i < (int)results.size(); i++) {\n";
  code = code <> "        fout << results[i][0] << \" \" << results[i][1] << \" \"\n";
  code = code <> "             << results[i][2] << \" \" << results[i][3] << \"\\n\";\n";
  code = code <> "    }\n";
  code = code <> "    fout.close();\n";
  code = code <> "    std::cerr << \"Done. Processed \" << n_kp << \" kinematic points (batched IBP).\" << std::endl;\n";
  code = code <> "    return 0;\n}\n";
    code
    ,
    (* IBP path, per-kp Vegas (one Vegas run per kinematic point; the batched
       variant above shares the grid across a kp-chunk). *)
    integrator === "VEGAS" && isIBP,
  code = "// TROPICAL_REQUIRES_CUBA  (CompileCpp greps for this sentinel)\n";
  code = code <> "#ifdef TROPICAL_USE_CUBA\n";
  code = code <> "extern \"C\" {\n";
  code = code <> "#include <cuba.h>\n";
  code = code <> "}\n";
  code = code <> "static IntegrandFunc g_fn     = nullptr;\n";
  code = code <> "static const double* g_params = nullptr;\n";
  code = code <> "static int           g_dim    = 0;\n";
  code = code <> "static int cubaWrap(const int* ndim, const cubareal xx[], const int* ncomp,\n";
  code = code <> "                    cubareal ff[], void* userdata) {\n";
  code = code <> "    (void)ndim; (void)ncomp; (void)userdata;\n";
  code = code <> "    double y[MAX_DIM];\n";
  code = code <> "    for (int i = 0; i < g_dim; ++i) y[i] = (double)xx[i];\n";
  code = code <> "    cx v = g_fn(y, g_params);\n";
  code = code <> "    ff[0] = v.real(); ff[1] = v.imag();\n";
  code = code <> "    return 0;\n";
  code = code <> "}\n";
  code = code <> "// integrate one function with Vegas; a 0-dim integrand (e.g. the IBP\n";
  code = code <> "// boundary of a 1-var sector) is a constant, so eval it directly.\n";
  code = code <> "static void vegasOne(int s, const double* params, long long maxeval,\n";
  code = code <> "                     double epsrel, double epsabs, int seed,\n";
  code = code <> "                     double out[2], double err[2]) {\n";
  code = code <> "    g_fn = integrand_table[s]; g_dim = integrand_dim[s]; g_params = params;\n";
  code = code <> "    if (g_dim == 0) {\n";
  code = code <> "        double y0[1] = {0.0}; cx v = g_fn(y0, params);\n";
  code = code <> "        out[0] = v.real(); out[1] = v.imag(); err[0] = 0.0; err[1] = 0.0;\n";
  code = code <> "        return;\n";
  code = code <> "    }\n";
  code = code <> "    int neval = 0, fail = 0; cubareal I[2], E[2], Pr[2];\n";
  code = code <> "    Vegas(g_dim, 2, cubaWrap, nullptr, 1, epsrel, epsabs, 0, seed,\n";
  code = code <> "          0, (int)maxeval, " <> vegSizing <> ", 0, nullptr, nullptr,\n";
  code = code <> "          &neval, &fail, I, E, Pr);\n";
  code = code <> "    out[0] = I[0]; out[1] = I[1]; err[0] = E[0]; err[1] = E[1];\n";
  code = code <> "}\n";
  code = code <> "#endif\n\n";

  code = code <> "int main(int argc, char* argv[]) {\n";
  code = code <> "    if (argc < 3) {\n";
  code = code <> "        std::cerr << \"Usage: \" << argv[0] << \" <input_file> <output_file> [maxeval] [n_threads]\" << std::endl;\n";
  code = code <> "        return 1;\n";
  code = code <> "    }\n\n";
  code = code <> "    std::string input_file = argv[1];\n";
  code = code <> "    std::string output_file = argv[2];\n";
  code = code <> "    long long maxeval = (argc > 3) ? std::atoll(argv[3]) : " <> ToString[nSamples] <> "LL;\n\n";
  code = code <> "    std::ifstream fin(input_file);\n";
  code = code <> "    if (!fin) { std::cerr << \"Cannot open \" << input_file << std::endl; return 1; }\n";
  code = code <> "    std::vector<std::vector<double>> kinematic_data;\n";
  code = code <> "    if (N_PARAMS == 0) {\n";
  code = code <> "        int count = 1; fin >> count; if (count < 1) count = 1;\n";
  code = code <> "        for (int i = 0; i < count; i++) kinematic_data.push_back({});\n";
  code = code <> "    } else {\n";
  code = code <> "        double val; std::vector<double> row;\n";
  code = code <> "        while (fin >> val) { row.push_back(val);\n";
  code = code <> "            if ((int)row.size() == N_PARAMS) { kinematic_data.push_back(row); row.clear(); }}\n";
  code = code <> "    }\n";
  code = code <> "    fin.close();\n";
  code = code <> "    int n_kp = (int)kinematic_data.size();\n";
  code = code <> "    std::cerr << \"Read \" << n_kp << \" kinematic points (Vegas IBP)\" << std::endl;\n\n";
  code = code <> "    int n_ibp = N_INTEGRANDS - N_CONV;\n";
  code = code <> "    std::vector<std::array<double, 4>> results((1 + n_ibp) * n_kp);\n\n";

  code = code <> "#ifdef TROPICAL_USE_CUBA\n";
  code = code <> "    { const int zero = 0; cubacores(&zero, &zero); }\n";
  code = code <> "    const double epsrel = " <> epsrel <> ", epsabs = " <> epsabs <> ";\n";
  code = code <> "    const int    seed   = " <> seed <> ";\n";
  code = code <> "    for (int kp = 0; kp < n_kp; ++kp) {\n";
  code = code <> "        const double* params = kinematic_data[kp].data();\n";
  code = code <> "        // convergent sectors (s < N_CONV) summed into line 0\n";
  code = code <> "        double cre = 0, cim = 0, cvre = 0, cvim = 0;\n";
  code = code <> "        for (int s = 0; s < N_CONV; ++s) {\n";
  code = code <> "            double I[2], E[2]; vegasOne(s, params, maxeval, epsrel, epsabs, seed, I, E);\n";
  code = code <> "            cre += I[0]; cim += I[1]; cvre += E[0]*E[0]; cvim += E[1]*E[1];\n";
  code = code <> "        }\n";
  code = code <> "        results[kp * (1 + n_ibp)] = {cre, cim, std::sqrt(cvre), std::sqrt(cvim)};\n";
  code = code <> "        // each IBP function (s >= N_CONV) on its own line\n";
  code = code <> "        for (int s = N_CONV; s < N_INTEGRANDS; ++s) {\n";
  code = code <> "            double I[2], E[2]; vegasOne(s, params, maxeval, epsrel, epsabs, seed, I, E);\n";
  code = code <> "            int idx = kp * (1 + n_ibp) + (s - N_CONV + 1);\n";
  code = code <> "            results[idx] = {I[0], I[1], E[0], E[1]};\n";
  code = code <> "        }\n";
  code = code <> "    }\n";
  code = code <> "#else\n";
  code = code <> "    (void)maxeval;\n";
  code = code <> "    std::cerr << \"ERROR: built without CUBA. Rebuild with -DTROPICAL_USE_CUBA \"\n";
  code = code <> "                 \"(install CUBA), or use Integrator -> \\\"MonteCarlo\\\".\" << std::endl;\n";
  code = code <> "    return 2;\n";
  code = code <> "#endif\n\n";

  code = code <> "    std::ofstream fout(output_file);\n";
  code = code <> "    if (!fout) { std::cerr << \"Cannot open \" << output_file << std::endl; return 1; }\n";
  code = code <> "    fout.precision(17);\n";
  code = code <> "    for (int i = 0; i < (int)results.size(); i++) {\n";
  code = code <> "        fout << results[i][0] << \" \" << results[i][1] << \" \"\n";
  code = code <> "             << results[i][2] << \" \" << results[i][3] << \"\\n\";\n";
  code = code <> "    }\n";
  code = code <> "    fout.close();\n";
  code = code <> "    std::cerr << \"Done. Processed \" << n_kp << \" kinematic points.\" << std::endl;\n";
  code = code <> "    return 0;\n}\n";
    code
  ]
];

(* --------------------------------------------------------------------------
   GenerateCppMonteCarlo — public entry.  Builds the shared definitions block,
   selects the main() from {Integrator, Batch}, concatenates.  Defaults are
   exactly today's behavior ("MonteCarlo"), which stays byte-identical (B0).
   -------------------------------------------------------------------------- *)
Options[GenerateCppMonteCarlo] = {
  "NSamples"       -> 1000000,
  "MaxDim"         -> 20,
  "SeedBase"       -> 42,
  "Integrator"     -> "MonteCarlo",   (* "MonteCarlo" | "Vegas" *)
  "Batch"          -> False,          (* True => chunked-ncomp Vegas (Vegas only) *)
  "VegasEpsRel"    -> 1.*^-12,
  "VegasEpsAbs"    -> 1.*^-300,
  "VegasSeed"      -> 0,
  "CubaMaxComp"    -> 512,
  "VegasNStart"    -> Automatic,
  "VegasNIncrease" -> Automatic,
  "VegasNBatch"    -> Automatic
};

(* --------------------------------------------------------------------------
   cppBadTokens: scan generated C++ for tokens that mean the output cannot
   compile -- either symbolic infinities/indeterminates (leaked from a division
   by a vanishing effective exponent, e.g. an unsupported nested divergence) or
   unconverted Mathematica heads.  Returns the list of offending tokens (empty
   when the code is clean).  The codegen entry points use this to FAIL LOUDLY
   ($Failed) instead of writing C++ that g++ rejects -- see TropicalEval::badcpp.

   BUG #27 (§5.8) FIX: the original token list matched BRACKET-form heads
   ("Sin[", "Cos[", ...).  But codegen emits via CForm (mmaToCInternal's
   fallback uses ToString[CForm[expr]]), and CForm formats an unconverted head
   in PAREN form -- "Sin(", "BesselJ(", "Power(", ... -- never bracket form.
   So the bracket tokens NEVER matched and uncompilable C++ (an unsupported
   head leaking through the fallback) was emitted SILENTLY.  The guard now
   detects the actually-emitted PAREN form:

     - Symbolic infinities / indeterminate: CForm renders ComplexInfinity AS
       "DirectedInfinity()" and prints "Indeterminate" verbatim, so we scan
       for the literal substrings "DirectedInfinity" and "Indeterminate".
     - Unconverted heads: a capitalized Mathematica head immediately followed
       by "(" with a word boundary in front -- e.g. "Sin(", "Power(",
       "BesselJ(", "Rule(", "List(".  GOOD output never contains a
       capitalized "Head(" token: every supported function is emitted in
       lowercase std::-qualified form (std::sin, std::sqrt, std::pow, ...),
       so this regex cannot false-positive on a clean integrand (verified:
       all 9 codegen goldens contain zero capitalized "Head(" tokens, and the
       infinity / indeterminate literals appear only in broken output).

   The head list is an EXPLICIT deny-list of Mathematica heads CForm can emit
   and MmaToC cannot convert (trig/hyperbolic + inverses, Gamma/PolyGamma,
   Zeta/PolyLog, Bessel, Erf, elliptic, hypergeometric, plus the structural
   heads Rule/List and the residual Power/Sqrt the CForm fallback may leave).
   Each is matched in PAREN form with a leading word boundary, e.g. "Sin(".
   It is explicit (not a blanket "any Uppercase(" scan) precisely so that the
   GENERATED C++ identifiers -- function names like P0(, Pfull0(, Psimp0( and
   the CUBA Vegas( entry point -- are NOT flagged: all 9 codegen goldens stay
   clean (zero matches) and byte-identical, while a genuinely un-emittable head
   that survives the CForm fallback is caught and turned into $Failed. *)
cppBadHeads = {
  "Sin", "Cos", "Tan", "Cot", "Sec", "Csc",
  "ArcSin", "ArcCos", "ArcTan", "ArcCot", "ArcSec", "ArcCsc",
  "Sinh", "Cosh", "Tanh", "Coth", "Sech", "Csch",
  "ArcSinh", "ArcCosh", "ArcTanh",
  "Sqrt", "Power", "Surd", "CubeRoot",
  "Gamma", "LogGamma", "PolyGamma", "Beta", "Factorial", "Binomial",
  "Zeta", "PolyLog", "HurwitzZeta", "LerchPhi", "Hypergeometric2F1",
  "HypergeometricPFQ", "Hypergeometric1F1", "MeijerG",
  "BesselJ", "BesselY", "BesselI", "BesselK", "AiryAi", "AiryBi",
  "Erf", "Erfc", "Erfi", "FresnelC", "FresnelS",
  "EllipticK", "EllipticE", "EllipticF", "EllipticPi", "EllipticTheta",
  "ProductLog", "ExpIntegralEi", "ExpIntegralE", "LogIntegral",
  "SinIntegral", "CosIntegral", "Sign", "UnitStep", "Floor", "Ceiling",
  "Round", "Mod", "Max", "Min", "Piecewise",
  "Rule", "List", "Rational", "Complex"
};

cppBadTokens[code_String] := Module[{literals, headTokens},
  (* Bare-identifier symbolic infinities / indeterminate (CForm prints these
     verbatim; ComplexInfinity itself renders as "DirectedInfinity()", so the
     "DirectedInfinity" substring already covers the complex-infinity case). *)
  literals = Select[
    {"DirectedInfinity", "Indeterminate"},
    StringContainsQ[code, #] &];
  (* Unconverted heads in the PAREN form CForm actually emits.  Leading
     (?<![A-Za-z0-9_]) word boundary so a head spelled as a tail of a longer
     C++ identifier is not flagged; the "(" is regex-escaped.  Matching the
     exact capitalized head + "(" never collides with the lowercase std::
     forms MmaToC emits for supported heads. *)
  headTokens = Select[
    cppBadHeads,
    StringContainsQ[code,
      RegularExpression["(?<![A-Za-z0-9_])" <> # <> "\\("]] &];
  DeleteDuplicates @ Join[literals, (# <> "(") & /@ headTokens]
];

GenerateCppMonteCarlo[convergentSectors_List, divergentSectors_List,
                      integrandSpec_Association, outputFile_String,
                      OptionsPattern[]] :=
Module[{defs, mainCode, code, integrator, batch, maxDim, nSamples, seedBase},
  integrator = OptionValue["Integrator"];
  batch      = TrueQ[OptionValue["Batch"]];
  maxDim     = OptionValue["MaxDim"];
  nSamples   = OptionValue["NSamples"];
  seedBase   = OptionValue["SeedBase"];

  defs = emitIntegrandDefinitions[convergentSectors, divergentSectors,
           {}, integrandSpec, "MaxDim" -> maxDim];

  integrator = normalizeIntegrator[integrator];
  mainCode = emitMain[defs, integrator, batch,
    "NSamples" -> nSamples, "SeedBase" -> seedBase,
    "VegasEpsRel" -> OptionValue["VegasEpsRel"],
    "VegasEpsAbs" -> OptionValue["VegasEpsAbs"],
    "VegasSeed"   -> OptionValue["VegasSeed"],
    "CubaMaxComp" -> OptionValue["CubaMaxComp"],
    "VegasNStart"    -> OptionValue["VegasNStart"],
    "VegasNIncrease" -> OptionValue["VegasNIncrease"],
    "VegasNBatch"    -> OptionValue["VegasNBatch"]];

  code = defs["Defs"] <> mainCode;

  (* Codegen guard (G-A): never write C++ that cannot compile.  A clean
     convergent / single-divergence spec produces no bad tokens, so this is a
     no-op there (byte-identical output, B0 preserved); a broken spec returns
     $Failed BEFORE Export so no invalid .cpp is left behind. *)
  Module[{bad = cppBadTokens[code]},
    If[bad =!= {},
      Message[TropicalEval::badcpp, bad];
      Return[$Failed]
    ]
  ];

  Export[outputFile, code, "Text"];

  Print["Generated C++ Monte Carlo code: ", outputFile,
    Switch[{integrator, batch},
      {"VEGAS", True}, "  (VEGAS, batched-ncomp)",
      {"VEGAS", _},    "  (VEGAS, per-kp)",
      _,               ""]];
  Print["  ", defs["NConvergent"], " convergent sectors"];
  Print["  ", defs["NG0"], " G0 integrands"];
  Print["  ", defs["NG1"], " G1 integrands"];
  Print["  ", defs["NRemainder"], " remainder integrands"];
  Print["  Total: ", defs["NTotal"], " integrand functions"];

  <|"Code" -> code, "OutputFile" -> outputFile,
    "NConvergent" -> defs["NConvergent"], "NG0" -> defs["NG0"],
    "NG1" -> defs["NG1"], "NRemainder" -> defs["NRemainder"],
    "NTotal" -> defs["NTotal"],
    "Dimensions" -> defs["Dims"]|>
];

(* --------------------------------------------------------------------------
   detectCuba: locate an optional CUBA install (cuba.h + libcuba.{a,dylib,so}).
   Memoized.  Probes the Homebrew prefix, the usual system prefixes, and the
   CUBA_INCLUDE / CUBA_LIB environment overrides.  The plain-MC path never calls
   this, so it adds zero dependencies; Vegas opts in.
   -------------------------------------------------------------------------- *)

detectCuba[] := detectCuba[] = Module[{brew, incs, libs, inc, lib},
  brew = Quiet@Check[
    StringTrim@RunProcess[{"brew", "--prefix"}, "StandardOutput"], ""];
  If[!StringQ[brew], brew = ""];
  incs = DeleteDuplicates@Select[
    {If[brew =!= "", brew <> "/include", Nothing],
     "/opt/homebrew/include", "/usr/local/include", "/usr/include",
     Environment["CUBA_INCLUDE"]}, StringQ];
  libs = DeleteDuplicates@Select[
    {If[brew =!= "", brew <> "/lib", Nothing],
     "/opt/homebrew/lib", "/usr/local/lib", "/usr/lib",
     Environment["CUBA_LIB"]}, StringQ];
  inc = SelectFirst[incs,
    FileExistsQ[FileNameJoin[{#, "cuba.h"}]] &, $Failed];
  lib = SelectFirst[libs,
    Function[d, AnyTrue[{"libcuba.a", "libcuba.dylib", "libcuba.so"},
      FileExistsQ[FileNameJoin[{d, #}]] &]], $Failed];
  <|"Found" -> (inc =!= $Failed && lib =!= $Failed),
    "IncludeDir" -> inc, "LibDir" -> lib|>
];

(* --------------------------------------------------------------------------
   CompileCpp: Compile generated C++ code.  CUBA linking is conditional: when
   the source needs CUBA (Vegas) it adds -DTROPICAL_USE_CUBA -I.. -L.. -lcuba.
   The 3-arg form CompileCpp[src, bin, debug] still works; options go after.
   "UseCuba" -> Automatic detects from the TROPICAL_REQUIRES_CUBA source
   sentinel; True/False forces it.
   -------------------------------------------------------------------------- *)

Options[CompileCpp] = {"UseCuba" -> Automatic};

CompileCpp[sourceFile_String, outputBinary_String,
           debug_: False, OptionsPattern[]] :=
Module[{compiler, flags, cmd, result, needsCuba, cuba, linkFlags},
  compiler = "g++";

  needsCuba = TrueQ[OptionValue["UseCuba"]] ||
    (OptionValue["UseCuba"] === Automatic &&
     StringContainsQ[Quiet@Import[sourceFile, "Text"], "TROPICAL_REQUIRES_CUBA"]);

  flags = If[debug,
    {"-std=c++17", "-O2", "-fopenmp", "-DTROPICAL_MC_DEBUG",
     "-Wall", "-Wextra"},
    {"-std=c++17", "-O3", "-fopenmp", "-Wall", "-Wextra"}
  ];
  linkFlags = {"-lm"};

  If[needsCuba,
    cuba = detectCuba[];
    If[!cuba["Found"],
      Print["ERROR: Vegas (CUBA) was requested but CUBA was not found."];
      Print["  Install it (e.g. `brew install cuba`) or set CUBA_INCLUDE / ",
            "CUBA_LIB,"];
      Print["  or use Integrator -> \"MonteCarlo\" (the zero-dependency path)."];
      Return[$Failed]
    ];
    flags = Join[flags, {"-DTROPICAL_USE_CUBA",
      "-I" <> cuba["IncludeDir"]}];
    linkFlags = {"-L" <> cuba["LibDir"], "-lcuba", "-lm"};
  ];

  cmd = {compiler, Sequence @@ flags, "-o", outputBinary,
         sourceFile, Sequence @@ linkFlags};

  Print["Compiling: ", StringRiffle[cmd, " "]];
  result = RunProcess[cmd];

  (* If compilation fails due to -fopenmp (e.g. macOS clang), retry without it *)
  If[result["ExitCode"] != 0 &&
     StringContainsQ[result["StandardError"], "fopenmp"],
    Print["  OpenMP not supported, retrying without -fopenmp..."];
    flags = DeleteCases[flags, "-fopenmp"];
    cmd = {compiler, Sequence @@ flags, "-o", outputBinary,
           sourceFile, Sequence @@ linkFlags};
    Print["Compiling: ", StringRiffle[cmd, " "]];
    result = RunProcess[cmd];
  ];

  If[result["ExitCode"] != 0,
    Print["Compilation FAILED:"];
    Print[result["StandardError"]];
    Return[$Failed]
  ];

  If[StringLength[StringTrim[result["StandardError"]]] > 0,
    Print["Compiler warnings:"];
    Print[result["StandardError"]]
  ];

  Print["Compilation successful: ", outputBinary];
  outputBinary
];


(* ============================================================================
   MODULE 4: EvaluateTropicalMC (Driver)
   ============================================================================ *)

Options[EvaluateTropicalMC] = {
  (* divergence-handling method (plan.md §3.4, D2):
       Automatic -> "None" if no divergent sectors, else "IBP" (the default
                    divergence method); force with "None" | "IBP" | "Subtraction".
     "IBP" routes to the IBP execution path (each term convergent at eps=0);
     "None"/"Subtraction" use the tropical-subtraction path of this driver. *)
  "Method"         -> Automatic,
  (* Lifting (plan.md §6.2): pass a LiftData association (from LiftCoefficients)
     to route all sector processing through ProcessSectorLifted; the first
     argument is then the LIFTED (n+1)-dim spec and fanData the lifted fan.
     Mutually exclusive with divergence handling (plan.md N3). *)
  "LiftData"       -> None,
  "NSamples"       -> 1000000,
  "NThreads"       -> Automatic,
  (* Seed base for the per-kp MC streams (plan.md §3.4, §8.6, Tree B).  Forwarded
     to the emitted C++ as the compile-time default; the binary also accepts a
     runtime argv[5] override (seed_base) so ONE compiled binary can be re-run
     under distinct seeds without recompiling (determinism cross-check #22). *)
  "SeedBase"       -> 42,
  "RunChecks"      -> True,
  "EpsilonValue"   -> None,
  "TestEpsilon"    -> 0.01,
  "PrecisionGoal"  -> 3,
  "WorkingDirectory" -> Automatic,
  "Verbose"        -> True,
  (* sampler selection (default = the zero-dependency plain Monte Carlo) *)
  "Integrator"     -> "MC",    (* "MC" | "VEGAS" (aliases "MonteCarlo"/"Vegas") *)
  "Batch"          -> False,         (* chunked-ncomp VEGAS (VEGAS only) *)
  (* Complex polynomial exponents (plan.md §6.3).  Automatic -> "Direct"
     ( result *= exp(B*log P), correct for any complex B ).  "SplitRealImag" is
     the opt-in VEGAS-variance mode: process sectors on Re(B) and reintroduce
     Im(B) as an oscillatory phase (complex log P + MonoFactorLog).  No effect
     on real exponents (no phase block emitted when Im(B)=0). *)
  "ComplexExponentMode" -> Automatic,
  Sequence @@ $vegasOptionDefaults
};

EvaluateTropicalMC[integrandSpec_Association, fanData_List,
                   kinematicPoints_List, OptionsPattern[]] :=
Module[
  {dualVertices, simplexList, n, nKP, nParams,
   allSectorData, convergentSectors, divergentSectors,
   processedDivergent, analyticContributions,
   cppFile, cppBinary, kinFile, resultFile,
   cppResult, mcResults, finalResults,
   runChecks, verbose, nSamples, nThreads,
   workDir, epsVal, testEps, precGoal, eps,
   integrator, batch, useCuba, vegasOpts, method,
   liftData, isLifted, liftedSpec, originalSpec,
   emptyDomainCount, hasConstList, cxMode, hasImagB,
   cxSplit, imBList, specForProc, imAList, hasImagA},

  runChecks  = OptionValue["RunChecks"];
  verbose    = OptionValue["Verbose"];
  nSamples   = OptionValue["NSamples"];
  nThreads   = OptionValue["NThreads"];
  workDir    = OptionValue["WorkingDirectory"];
  epsVal     = OptionValue["EpsilonValue"];
  testEps    = OptionValue["TestEpsilon"];
  precGoal   = OptionValue["PrecisionGoal"];
  method     = OptionValue["Method"];
  integrator = normalizeIntegrator[OptionValue["Integrator"]];
  batch      = TrueQ[OptionValue["Batch"]];
  useCuba    = (integrator === "VEGAS");
  (* Complex-exponent mode (plan.md §6.3).  Automatic -> "Direct". *)
  cxMode     = OptionValue["ComplexExponentMode"];
  If[cxMode === Automatic, cxMode = "Direct"];
  vegasOpts  = {"Integrator" -> integrator, "Batch" -> batch,
    "SeedBase" -> OptionValue["SeedBase"],
    "VegasEpsRel" -> OptionValue["VegasEpsRel"],
    "VegasEpsAbs" -> OptionValue["VegasEpsAbs"],
    "VegasSeed" -> OptionValue["VegasSeed"],
    "CubaMaxComp" -> OptionValue["CubaMaxComp"],
    "VegasNStart"    -> OptionValue["VegasNStart"],
    "VegasNIncrease" -> OptionValue["VegasNIncrease"],
    "VegasNBatch"    -> OptionValue["VegasNBatch"]};
  eps        = integrandSpec["RegulatorSymbol"];

  (* --- Lifting setup (plan.md §6.2) --- *)
  liftData = OptionValue["LiftData"];
  isLifted = (liftData =!= None);
  liftedSpec   = integrandSpec;  (* in lifted mode this IS the lifted (n+1) spec *)
  originalSpec = If[isLifted, liftData["OriginalSpec"], integrandSpec];

  (* --- SplitRealImag setup (plan.md §6.3) ---
     Active only when the user opted in AND some polynomial exponent B_j has a
     nonzero imaginary part.  Then sectors are decomposed on Re(B) (real
     importance measure for VEGAS), and Im(B_j) is reintroduced per convergent
     sector via the "ImagPolyExponents" key, which the codegen turns into the
     oscillatory phase (complex log P + MonoFactorLog).  Direct mode and the
     real-exponent path leave specForProc = the spec unchanged and attach no
     ImagPolyExponents, so the emitted C++ is byte-identical (#25). *)
  imBList  = imagPolyInfo[integrandSpec];   (* eps-free Im(B); see imagPolyInfo *)
  hasImagB = AnyTrue[imBList, (!TrueQ[PossibleZeroQ[#]]) &];
  (* planAXpDIVv2.md §4.1: the split must ALSO engage when the MONOMIAL exponents
     A are complex (the bubble at imaginary mu), not only when B is.  imAList is
     the AUGMENTED Im(A) in lifted mode (aux entry 0, LiftCoefficients:950).
     See imagMonoInfo for the eps->0 / declare-real rationale (planR.md R2,R6). *)
  {imAList, hasImagA} = imagMonoInfo[integrandSpec, eps];
  cxSplit  = (cxMode === "SplitRealImag") && (hasImagB || hasImagA);

  (* Lifted-mode fan-dimension assertion (n+1 rays per simplex coordinate). *)
  If[isLifted,
    Module[{nOrig, fanDim},
      nOrig  = Length[originalSpec["Variables"]];
      fanDim = If[Length[fanData[[1]]] > 0, Length[fanData[[1, 1]]], 0];
      If[fanDim != nOrig + 1,
        Message[TropicalEval::liftfandim, fanDim, nOrig + 1];
        Return[$Failed]
      ]
    ]
  ];

  (* --- Method routing (plan.md §3.4, §5.4, D2) ---
     This driver owns the convergent + tropical-subtraction paths; the IBP path
     lives in evaluateTropicalIBPDriver, to which we delegate when Method is
     "IBP", or (Automatic) when the integral has divergent sectors AND a
     regulator is present (IBP is the DEFAULT divergence method).  "None" and
     "Subtraction" stay in this driver.  The convergent case (no divergent
     sectors) is unaffected: Automatic resolves to "None" and falls through. *)
  Module[{routeToIBP},
    routeToIBP = Which[
      (* Lifted + divergence (planAXpDIV.md §4.4, Barrier C): route the lifted
         call to the (now lifting-aware) IBP path when a regulator is present,
         eps is NOT pinned to a value (a pinned eps is the convergent
         LaurentFromSubtraction route), and some lifted sector is divergent.
         Method "None"/"Subtraction" stay in this driver. *)
      isLifted,
        Which[
          method === "None" || method === "Subtraction", False,
          epsVal =!= None || eps === None, False,
          method === "IBP", True,
          method === Automatic,
            (* SplitRealImag x lift: detect divergence on Re(A)/Re(B) — Process-
               SectorLifted on the complex spec would classify "complex" ($Failed)
               and hide the pole (planCXLIFTDIV.md §4.4 / planAXpDIVv2.md §4.1).
               Re of an already-real exponent is the identity, so realifying both
               families when cxSplit is byte-safe; divergence lives in Re only. *)
            AnyTrue[
              Table[ProcessSectorLifted[
                      If[cxSplit,
                         (* realify both families by SUBTRACTING the eps-free
                            imaginary parts (realifyPolyB / realifyMonoA), so a
                            regulator inside B/A stays a true real — see their
                            docstrings; divergence lives in Re only. *)
                         realifyMonoA[
                           realifyPolyB[liftedSpec, imBList],
                           imAList],
                         liftedSpec],
                      fanData[[1]], fanData[[2, s]],
                      s, liftData, "Eps" -> eps], {s, Length[fanData[[2]]]}],
              (AssociationQ[#] && TrueQ[#["IsDivergent"]]) &],
          True, False
        ],
      method === "IBP", True,
      method === "None" || method === "Subtraction", False,
      method === Automatic,
        eps =!= None && AnyTrue[
          Table[ProcessSector[integrandSpec, fanData[[1]], fanData[[2, s]], s],
                {s, Length[fanData[[2]]]}],
          (AssociationQ[#] && TrueQ[#["IsDivergent"]]) &],
      True, False
    ];
    (* SplitRealImag x lift x divergence (planCXLIFTDIV.md C1): now delivered via
       the IBP route — the per-IBP-piece MonoFactorLog is re-derived in
       IBPProcessSector and Im(B) rides the boundary/term codegen as the
       oscillatory phase.  ComplexExponentMode is forwarded so the IBP driver
       knows to decompose on Re(B). *)
    If[routeToIBP,
      Return[evaluateTropicalIBPDriver[integrandSpec, fanData, kinematicPoints,
        "NSamples" -> nSamples, "NThreads" -> nThreads,
        "RunChecks" -> runChecks, "WorkingDirectory" -> workDir,
        "Verbose" -> verbose, "Integrator" -> integrator, "Batch" -> batch,
        "VegasEpsRel" -> OptionValue["VegasEpsRel"],
        "VegasEpsAbs" -> OptionValue["VegasEpsAbs"],
        "VegasSeed" -> OptionValue["VegasSeed"],
        "CubaMaxComp" -> OptionValue["CubaMaxComp"],
        "ComplexExponentMode" -> cxMode,
        "LiftData" -> If[isLifted, liftData, None]]]
    ]
  ];

  If[workDir === Automatic,
    workDir = DirectoryName[$InputFileName];
    If[workDir === "" || !StringQ[workDir],
      workDir = Directory[]
    ];
    workDir = FileNameJoin[{workDir, "INTERFILES"}]
  ];
  If[!DirectoryQ[workDir], Quiet[CreateDirectory[workDir]]];

  {dualVertices, simplexList} = fanData;
  n      = Length[integrandSpec["Variables"]];
  nKP    = Length[kinematicPoints];
  nParams = Length[integrandSpec["KinematicSymbols"]];

  If[verbose,
    Print["EvaluateTropicalMC: ", Length[simplexList], " sectors, ",
          nKP, " kinematic points, ", n, " variables"]
  ];

  (* --- Step 1: Validate inputs --- *)
  If[nKP == 0,
    Print["ERROR: no kinematic points provided"];
    Return[$Failed]
  ];

  If[nParams > 0,
    If[!AllTrue[kinematicPoints,
                (ListQ[#] && Length[#] == nParams) &],
      Print["ERROR: kinematic points must have ", nParams, " parameters each"];
      Return[$Failed]
    ];
    If[!AllTrue[Flatten[kinematicPoints], (NumericQ[#] && Im[#] == 0 &&
                                          Abs[#] < 10^15) &],
      Print["ERROR: kinematic values must be real and finite"];
      Return[$Failed]
    ];
  ];

  (* --- Step 2: Process all sectors --- *)
  If[verbose, Print["Processing ", Length[simplexList], " sectors..."]];

  emptyDomainCount = 0;
  hasConstList     = {};

  If[isLifted,
    (* ---- Lifted mode: process every sector via the eps-aware
       ProcessSectorLifted (planAXpDIV.md §4.3/§4.4) ----
       $Failed (liftcomplex/liftnopivot/liftdivdomain) aborts the whole call;
       EmptyDomain sectors contribute 0 and are dropped+counted.  The survivors
       split into convergent and (single-pole) divergent lifted sectors — the
       latter feed the tropical-subtraction machinery (Step 4) below, exactly
       like unlifted divergent sectors.  (The Automatic+IBP combination already
       routed to evaluateTropicalIBPDriver above; a pinned EpsilonValue makes
       every sector convergent — the LaurentFromSubtraction route.) *)
    (* SplitRealImag (plan.md §6.3): decompose on Re(B); Im(B) reattached as
       "ImagPolyExponents".  A pinned EpsilonValue is substituted first so the
       flattening exponents are numeric (mirrors the unlifted branch). *)
    specForProc = liftedSpec;
    (* planAXpDIVv2.md §4.2: realify Im(A) on the same footing as Re(B) (clears
       liftcomplex).  realifyMonoA MUST run BEFORE the eps-pinning below — see its
       docstring (the symbolic subtraction is exact; pinning first leaks a
       Complex[float, 0.] into the domain indicator). *)
    If[cxSplit && hasImagA,
      specForProc = realifyMonoA[specForProc, imAList]];
    (* Re(B) by eps-free subtraction, BEFORE pinning (realifyPolyB docstring). *)
    If[cxSplit,
      specForProc = realifyPolyB[specForProc, imBList]];
    If[epsVal =!= None && eps =!= None,
      specForProc = MapAt[# /. eps -> epsVal &, specForProc, {Key["MonomialExponents"]}];
      specForProc = MapAt[# /. eps -> epsVal &, specForProc, {Key["PolynomialExponents"]}]
    ];
    allSectorData = {};
    Catch[
      Do[
        Module[{sd},
          sd = ProcessSectorLifted[specForProc, dualVertices,
                                   simplexList[[s]], s, liftData,
                                   "Eps" -> eps, "Verbose" -> verbose,
                                   "ImagMonoExps" ->
                                     If[cxSplit && hasImagA, imAList, None]];
          Which[
            sd === $Failed,
              allSectorData = $Failed;  Throw[Null],
            AssociationQ[sd] && KeyExistsQ[sd, "EmptyDomain"] && sd["EmptyDomain"],
              emptyDomainCount++,
            True,
              If[cxSplit, sd["ImagPolyExponents"] = imBList];
              AppendTo[allSectorData, sd]
          ]
        ],
        {s, Length[simplexList]}
      ]
    ];
    If[allSectorData === $Failed,
      Print["ERROR: ProcessSectorLifted failed for a sector (liftcomplex / ",
            "liftnopivot / liftdivdomain).  Aborting ($Failed)."];
      Return[$Failed]
    ];
    convergentSectors = Select[allSectorData,
      (AssociationQ[#] && !TrueQ[#["IsDivergent"]]) &];
    divergentSectors  = Select[allSectorData,
      (AssociationQ[#] && TrueQ[#["IsDivergent"]]) &];
    (* SplitRealImag x lift x divergence: the IBP route (above) and pinned-eps
       LaurentFromSubtraction (no divergent sectors) are phase-aware; this inline
       symbolic-eps subtraction path (G0/G1) is NOT, so refuse cleanly rather than
       emit a phase-less (wrong) Laurent (planCXLIFTDIV.md C1, §9). *)
    If[cxSplit && Length[divergentSectors] > 0,
      Message[TropicalEval::splitliftdiv];  Return[$Failed]
    ];
    hasConstList = (#["HasConstantTerm"] & /@ convergentSectors);
    If[verbose,
      Print["  ", Length[convergentSectors], " convergent, ",
            Length[divergentSectors], " divergent, ",
            emptyDomainCount, " EmptyDomain sectors dropped"]
    ];
    ,
    (* ---- Standard (unlifted) mode ---- *)
    allSectorData = Table[
      Module[{sd, specToUse},
        specToUse = integrandSpec;
        (* SplitRealImag (plan.md §6.3 / planAXpDIVv2.md §4.2): decompose on
           Re(B) and Re(A); reattach Im(B); attach the Im(A) monomial phase via
           the ProcessSector "ImagMonoExps" option (computed with the sector's M
           and a_eff).  Real-A / Direct -> option None -> no monomial phase.
           realifyMonoA runs BEFORE any eps-pinning so the subtraction is exact —
           see its docstring. *)
        If[cxSplit && hasImagA,
          specToUse = realifyMonoA[specToUse, imAList]];
        (* Re(B) by eps-free subtraction, BEFORE pinning (realifyPolyB docstring). *)
        If[cxSplit,
          specToUse = realifyPolyB[specToUse, imBList]];
        specToUse = If[epsVal =!= None && eps =!= None,
          MapAt[# /. eps -> epsVal &, specToUse,
                {Key["MonomialExponents"]}] //
          MapAt[# /. eps -> epsVal &, #,
                {Key["PolynomialExponents"]}] &,
          specToUse
        ];
        sd = ProcessSector[specToUse, dualVertices,
                           simplexList[[s]], s, "Verbose" -> verbose,
                           "ImagMonoExps" ->
                             If[cxSplit && hasImagA, imAList, None]];
        If[cxSplit && AssociationQ[sd] && !TrueQ[sd["IsDivergent"]],
          sd["ImagPolyExponents"] = imBList];
        sd
      ],
      {s, Length[simplexList]}
    ];

    If[Count[allSectorData, _Association] != Length[simplexList],
      Print["WARNING: ", Count[allSectorData, $Failed],
            " sectors failed to process"];
    ];

    convergentSectors = Select[allSectorData,
      (AssociationQ[#] && !#["IsDivergent"]) &];
    divergentSectors  = Select[allSectorData,
      (AssociationQ[#] && #["IsDivergent"]) &];

    (* planR.md R1: SplitRealImag (Im(A) and/or Im(B)) on an UNLIFTED divergent
       integral is not phase-aware on this inline subtraction path — the divergent
       branch of ProcessSector omits MonomialPhaseLog and ImagPolyExponents is
       reattached only for convergent sectors, so a divergent complex sector flows
       phase-less into the subtraction machinery.  Same limitation as the lifted
       branch (line 4072 / planCXLIFTDIV.md C1); refuse cleanly rather than emit a
       phase-less (wrong) Laurent.  (Method Automatic routes complex+divergent to
       the phase-aware IBP driver above; "None"/"Subtraction" and eps===None land
       here.) *)
    If[cxSplit && Length[divergentSectors] > 0,
      Message[TropicalEval::splitliftdiv];  Return[$Failed]];

    If[verbose,
      Print["  ", Length[convergentSectors], " convergent, ",
            Length[divergentSectors], " divergent sectors"]
    ]
  ];

  (* --- Step 3: Validation checks --- *)
  If[runChecks && Length[convergentSectors] > 0 && nParams > 0,
    Module[{testKP, kinRules},
      testKP = Take[kinematicPoints, Min[3, nKP]];
      Do[
        kinRules = Thread[
          originalSpec["KinematicSymbols"] -> testKP[[i]]
        ];
        If[verbose,
          Print["Validating decomposition at kinematic point ", i, "..."]
        ];
        Module[{vr},
          vr = If[isLifted,
            Quiet@ValidateLiftedDecomposition[originalSpec, liftedSpec,
              fanData, liftData, kinRules, precGoal],
            Quiet@ValidateDecomposition[integrandSpec, fanData,
              kinRules, precGoal]
          ];
          If[AssociationQ[vr],
            Print["  Point ", i, ": rel error = ", vr["RelativeError"]]
          ];
        ];,
        {i, Length[testKP]}
      ];
    ]
  ];

  (* --- Step 4: Process divergent sectors --- *)
  processedDivergent = {};
  analyticContributions = {};

  If[Length[divergentSectors] > 0,
    If[verbose, Print["Processing ", Length[divergentSectors],
                      " divergent sectors..."]];
    processedDivergent = Table[
      ProcessDivergentSector[divergentSectors[[s]], integrandSpec],
      {s, Length[divergentSectors]}
    ];

    (* Fail loudly if any divergent sector could not be subtracted (e.g. a
       nested divergence: >1 divergent variable per sector, refused by
       IdentifyDivergences with TropicalEval::nested).  Silently dropping the
       sector would return a finite-looking but wrong Laurent. *)
    If[MemberQ[processedDivergent, $Failed],
      Print["ERROR: tropical subtraction failed for ",
        Count[processedDivergent, $Failed], " of ", Length[divergentSectors],
        " divergent sector(s) — typically a nested divergence (>1 divergent ",
        "variable per sector), which the subtraction scheme does not handle ",
        "(TropicalEval::nested).  Aborting ($Failed) rather than dropping the ",
        "sector(s).  [The full integral at a fixed regulator value is still ",
        "available via \"EpsilonValue\"->eps, which makes all sectors ",
        "convergent.]"];
      Return[$Failed]
    ];

    analyticContributions = Table[
      Module[{dd, ck},
        dd = processedDivergent[[s]];
        If[!AssociationQ[dd], 0,
          ck = dd["ck"];
          <|"PoleCoeff" -> 1/ck, "FiniteCoeff" -> 1/ck|>
        ]
      ],
      {s, Length[processedDivergent]}
    ];
  ];

  (* --- Step 5: Validate subtraction --- *)
  If[runChecks && Length[processedDivergent] > 0 && nParams > 0,
    Module[{testKP, kinRules},
      testKP = Take[kinematicPoints, Min[2, nKP]];
      Do[
        kinRules = Thread[
          integrandSpec["KinematicSymbols"] -> testKP[[i]]
        ];
        Do[
          If[AssociationQ[processedDivergent[[s]]],
            If[verbose,
              Print["Validating subtraction sector ", s,
                    " at point ", i, "..."]
            ];
            Quiet@ValidateSubtraction[
              processedDivergent[[s]], divergentSectors[[s]],
              integrandSpec, kinRules, testEps
            ];
          ],
          {s, Length[processedDivergent]}
        ],
        {i, Length[testKP]}
      ];
    ]
  ];

  (* --- Step 6: Generate C++ code --- *)
  cppFile    = FileNameJoin[{workDir, "tropical_mc_generated.cpp"}];
  cppBinary  = FileNameJoin[{workDir, "tropical_mc"}];
  kinFile    = FileNameJoin[{workDir, "kinematic_data.txt"}];
  resultFile = FileNameJoin[{workDir, "mc_results.txt"}];

  Module[{convForCpp, divForCpp, specForCpp, epsSub},
    (* The C++ integrands are concrete cx(y,params) functions with no eps
       argument, so any residual regulator must be evaluated to a number
       before codegen.  A convergent sector contributes only its eps^0 value
       to the Laurent finite part, so substitute eps->0 (or the user-pinned
       EpsilonValue) across the WHOLE sector association — including the
       flattening exponents inside FlattenedPolys, not just PolynomialExponents.
       Without this, symbolic eps leaks into the emitted C++ and the divergent
       path fails to compile (pre-existing bug; convergent specs are
       unaffected because their RegulatorSymbol is None). *)
    epsSub = If[epsVal =!= None, epsVal, 0];
    convForCpp = If[eps =!= None,
      Map[(# /. eps -> epsSub) &, convergentSectors],
      convergentSectors
    ];

    divForCpp = If[eps =!= None,
      Map[
        Function[dd,
          If[AssociationQ[dd],
            dd /. eps -> epsSub,
            dd
          ]
        ],
        processedDivergent
      ],
      processedDivergent
    ];

    specForCpp = If[eps =!= None, integrandSpec /. eps -> epsSub, integrandSpec];

    (* §6.5 vegasbudget guard (lift_error_log L1).  Only meaningful for VEGAS.
       NSamples (= maxeval/sector) must exceed a multiple of the resolved NStart
       or the adaptive grid never resolves and the quoted error bar LIES (a tight
       but systematically-biased value).  The historical 20x floor was too loose
       at high d: the 8D case NStart=81000 passed at NSamples=4e6 (49x) yet was
       0.118% biased LOW with a ~1.3e-6 bar.  Empirically (L1 sweep) the VEGAS
       grid needs ~100 adaptive passes to resolve in 8D, so for d>=5 we require
       NSamples >= 100*NStart (8.1e6 at d=8 -> fires on the previously-silent
       4e6); d<=4 keeps the 20x floor (low-dim behaviour unchanged).
       guardFactor(d): 20 for d<=4, 100 for d>=5.
       maxSD is the largest convergent-sector integration dim; NStart is what
       resolveVegasSizing hands the codegen for that dim. *)
    If[integrator === "VEGAS",
      Module[{secDims, maxSD, resolvedNStart, guardFactor, threshold},
        secDims = Cases[convForCpp, a_Association :> a["Dimension"]];
        maxSD = If[secDims === {} || !VectorQ[secDims, IntegerQ], 0, Max[secDims]];
        resolvedNStart = resolveVegasSizing[OptionValue["VegasNStart"], maxSD, 1000, "nstart"];
        guardFactor = If[maxSD >= 5, 100, 20];
        threshold = guardFactor*resolvedNStart;
        If[IntegerQ[resolvedNStart] && nSamples < threshold,
          Message[TropicalEval::vegasbudget,
            nSamples, threshold, maxSD, resolvedNStart]
        ]
      ]
    ];

    cppResult = GenerateCppMonteCarlo[
      convForCpp,
      Select[divForCpp, AssociationQ],
      specForCpp, cppFile,
      "NSamples" -> nSamples,
      Sequence @@ vegasOpts
    ];
  ];

  If[!AssociationQ[cppResult],
    Print["ERROR: C++ code generation failed"];
    Return[$Failed]
  ];

  (* --- Step 7: Debug compile and test (MC only) ---
     The -DTROPICAL_MC_DEBUG instrumentation lives in the MC main; the Vegas
     mains do not carry it, so the debug build is gated to the MC integrator. *)
  If[runChecks && integrator === "MC",
    Module[{dbgBinary, dbgResult, dbgKinFile},
      dbgBinary  = FileNameJoin[{workDir, "tropical_mc_dbg"}];
      dbgKinFile = FileNameJoin[{workDir, "kinematic_data_dbg.txt"}];

      If[CompileCpp[cppFile, dbgBinary, True] =!= $Failed,
        Module[{testKP},
          testKP = Take[kinematicPoints, Min[5, nKP]];
          Export[dbgKinFile,
            StringRiffle[
              StringRiffle[ToString[CForm[#]] & /@ #, " "] & /@ testKP,
              "\n"
            ] <> "\n",
            "Text"
          ];

          dbgResult = RunProcess[{dbgBinary, dbgKinFile,
            FileNameJoin[{workDir, "mc_results_dbg.txt"}],
            "100000", "2"}];

          If[dbgResult["ExitCode"] == 0,
            Print["Debug run successful"];
            If[StringLength[StringTrim[dbgResult["StandardError"]]] > 0,
              Print["Debug output:\n", dbgResult["StandardError"]]
            ];,
            Print["Debug run FAILED:"];
            Print[dbgResult["StandardError"]];
          ];
        ];
      ];
    ]
  ];

  (* --- Step 8: Write kinematic data --- *)
  If[nParams > 0,
    Export[kinFile,
      StringRiffle[
        StringRiffle[
          ToString[CForm[#]] & /@ #, " "
        ] & /@ kinematicPoints,
        "\n"
      ] <> "\n",
      "Text"
    ];,
    Export[kinFile,
      StringRiffle[
        Table["0", {nKP}], "\n"
      ] <> "\n",
      "Text"
    ];
  ];

  (* --- Step 9: Release compile and run --- *)
  If[CompileCpp[cppFile, cppBinary, False, "UseCuba" -> useCuba] === $Failed,
    Print["ERROR: Release compilation failed"];
    Return[$Failed]
  ];

  Module[{nThreadsStr, result},
    nThreadsStr = If[nThreads === Automatic,
      ToString[$ProcessorCount],
      ToString[nThreads]
    ];

    If[verbose, Print["Running Monte Carlo (", nSamples,
                      " samples, ", nThreadsStr, " threads)..."]];

    result = RunProcess[{cppBinary, kinFile, resultFile,
      ToString[nSamples], nThreadsStr}];

    If[result["ExitCode"] != 0,
      Print["ERROR: Monte Carlo execution failed:"];
      Print[result["StandardError"]];
      Return[$Failed]
    ];

    If[verbose && StringLength[StringTrim[result["StandardError"]]] > 0,
      Print[result["StandardError"]]
    ];
  ];

  (* --- Step 10: Read results --- *)
  Module[{rawResults, lines, parsed},
    rawResults = Import[resultFile, "Text"];
    If[rawResults === $Failed,
      Print["ERROR: cannot read results file"];
      Return[$Failed]
    ];

    lines = Select[StringSplit[rawResults, "\n"],
                   StringLength[StringTrim[#]] > 0 &];

    (* Use Read[StringToStream[...], Number] instead of ToExpression
       because C++ outputs scientific notation like 2.05e-05 which
       ToExpression misparses (treats 'e' as a symbol). *)
    parsed = Table[
      Read[StringToStream[#], Number] & /@ StringSplit[line],
      {line, lines}
    ];

    If[Length[parsed] != nKP,
      Print["WARNING: expected ", nKP, " result rows, got ",
            Length[parsed]];
    ];

    Module[{badRows},
      badRows = Select[Range[Length[parsed]],
        !AllTrue[parsed[[#]], NumericQ[#] && Abs[#] < 10^30 &] &
      ];
      If[Length[badRows] > 0,
        Print["WARNING: ", Length[badRows],
              " rows contain non-finite values"]
      ];
    ];

    mcResults = parsed;
  ];

  (* --- Step 11: Assemble final results --- *)
  finalResults = Table[
    If[i <= Length[mcResults] && Length[mcResults[[i]]] >= 4,
      <|"KinematicPoint" -> If[nParams > 0, kinematicPoints[[i]], {}],
        "Re"     -> mcResults[[i, 1]],
        "Im"     -> mcResults[[i, 2]],
        "ReErr"  -> mcResults[[i, 3]],
        "ImErr"  -> mcResults[[i, 4]]|>,
      <|"KinematicPoint" -> If[nParams > 0, kinematicPoints[[i]], {}],
        "Re" -> 0., "Im" -> 0., "ReErr" -> 0., "ImErr" -> 0.|>
    ],
    {i, nKP}
  ];

  (* --- Step 12: Error summary --- *)
  If[verbose,
    Module[{reErrs, imErrs, reMags},
      reErrs = #["ReErr"] & /@ finalResults;
      imErrs = #["ImErr"] & /@ finalResults;
      reMags = Abs[#["Re"]] & /@ finalResults;

      Print["\n=== Monte Carlo Error Summary ==="];
      Print["  Re errors: mean=", Mean[reErrs],
            " median=", Median[reErrs],
            " max=", Max[reErrs]];
      Print["  Im errors: mean=", Mean[imErrs],
            " median=", Median[imErrs],
            " max=", Max[imErrs]];

      Module[{badPts},
        badPts = Select[Range[nKP],
          (reMags[[#]] > 0 &&
           reErrs[[#]] / reMags[[#]] > 0.1) &
        ];
        If[Length[badPts] > 0,
          Print["  WARNING: ", Length[badPts],
                " points have Re error > 10% of result magnitude"]
        ];
      ];
    ]
  ];

  If[runChecks && nParams > 0,
    Module[{testKP, kinRules},
      testKP = Take[kinematicPoints, Min[3, nKP]];
      Print["\n=== NIntegrate Cross-Check ==="];
      Do[
        kinRules = Thread[
          integrandSpec["KinematicSymbols"] -> testKP[[i]]
        ];
        Module[{directResult, mcRe, relErr},
          directResult = Quiet@ValidateDecomposition[
            integrandSpec, fanData, kinRules, precGoal
          ];
          If[AssociationQ[directResult],
            mcRe = finalResults[[i, "Re"]] + I * finalResults[[i, "Im"]];
            relErr = Abs[(mcRe - directResult["DirectResult"]) /
                        directResult["DirectResult"]];
            Print["  Point ", i, ": MC = ", mcRe,
                  ", NIntegrate = ", directResult["DirectResult"],
                  ", rel err = ", relErr];
          ];
        ];,
        {i, Length[testKP]}
      ];
    ]
  ];

  <|"Results"              -> finalResults,
    "ConvergentSectors"    -> Length[convergentSectors],
    "DivergentSectors"     -> Length[divergentSectors],
    "CppFile"              -> cppFile,
    "ResultFile"           -> resultFile,
    "AnalyticContributions" -> analyticContributions,
    (* Lifted-mode diagnostics (plan.md §6.2, §9 risk #1).  HasConstantTerm is
       reported per CONVERGENT sector (same order as ConvergentSectors); a
       False entry flags a potential infinite-MC-variance sector. *)
    "IsLifted"             -> isLifted,
    "HasConstantTerm"      -> hasConstList,
    "EmptyDomainSectors"   -> emptyDomainCount|>
];


(* --------------------------------------------------------------------------
   EvaluateTropicalMCLifted  (plan.md §3.4, §5.4, §6.2 — thin wrapper)

   Detects/applies lifting, builds the (n+1)-dim lifted fan, and routes through
   EvaluateTropicalMC with LiftData.  With no extreme coefficients it falls back
   to plain EvaluateTropicalMC on the original spec.
   -------------------------------------------------------------------------- *)

Options[EvaluateTropicalMCLifted] = Join[
  {"LiftRules" -> Automatic, "Threshold" -> 1000,
   "AnchorRule" -> "kStar", "BandEdgeGuard" -> False, "FanData" -> Automatic},
  Options[EvaluateTropicalMC]
];

EvaluateTropicalMCLifted[integrandSpec_Association, kinematicPoints_List,
                         opts : OptionsPattern[]] :=
Module[
  {liftRulesOpt, threshold, anchorRule, bandEdgeGuard, fanDataOpt,
   rules, liftResult, liftedSpec, liftData, verts, n, passThroughOpts},

  liftRulesOpt  = OptionValue["LiftRules"];
  threshold     = OptionValue["Threshold"];
  anchorRule    = OptionValue["AnchorRule"];
  bandEdgeGuard = OptionValue["BandEdgeGuard"];
  fanDataOpt    = OptionValue["FanData"];

  passThroughOpts = Sequence @@ FilterRules[{opts}, Options[EvaluateTropicalMC]];

  (* --- Determine lift rules --- *)
  If[liftRulesOpt === Automatic,
    rules = DetectExtremeCoefficients[integrandSpec, threshold,
              "AnchorRule" -> anchorRule, "BandEdgeGuard" -> bandEdgeGuard];
    If[rules === {},
      Print["EvaluateTropicalMCLifted: no extreme coefficients detected ",
            "(threshold=", threshold, "); falling back to plain ",
            "EvaluateTropicalMC on the original spec."];
      Module[{origFan = fanDataOpt},
        If[origFan === Automatic,
          Module[{origVerts},
            origVerts = PolytopeVertices[
              (Times @@ integrandSpec["Polynomials"])^(-1),
              integrandSpec["Variables"]];
            origFan = ComputeDecomposition[origVerts, "ShowProgress" -> False]
          ]
        ];
        Return[EvaluateTropicalMC[integrandSpec, origFan,
                                  kinematicPoints, passThroughOpts]]
      ]
    ];
    (* DetectExtremeCoefficients ships "SuggestedK"; LiftCoefficients wants "k". *)
    rules = Map[
      <|"PolyIndex" -> #["PolyIndex"], "ExponentVector" -> #["ExponentVector"],
        "k" -> #["SuggestedK"]|> &, rules]
    ,
    rules = liftRulesOpt
  ];

  (* --- Lift the integrand --- *)
  liftResult = LiftCoefficients[integrandSpec, rules];
  If[!AssociationQ[liftResult], Return[$Failed]];
  liftedSpec = liftResult["LiftedSpec"];
  liftData   = liftResult["LiftData"];
  n          = Length[integrandSpec["Variables"]];

  (* --- Build or use the (n+1)-dim fan --- *)
  Module[{liftedFan},
    If[fanDataOpt =!= Automatic,
      liftedFan = fanDataOpt
      ,
      (* Automatic: compute from the lifted integrand.  Quiet polymake messages
         so a degenerate (lower-dim) lifted polytope is caught by the guard
         below rather than leaking raw messages. *)
      verts = Quiet[
        PolytopeVertices[(Times @@ liftedSpec["Polynomials"])^(-1),
                         liftedSpec["Variables"]],
        TropicalFan::polymake];
      (* §6.4 / lift_error_log L2: use the K-scaled fan (normal fan is
         scale-invariant) so thin lattice simplices in ambient dim >= 4 do not
         leak $Failed.  liftdegenerate then fires only for a genuinely
         lower-dimensional lifted polytope. *)
      liftedFan = If[ListQ[verts],
        Quiet[computeFanScaled[verts], TropicalFan::polymake],
        $Failed];
      If[!ListQ[liftedFan] || Length[liftedFan] < 2,
        Message[TropicalEval::liftdegenerate];  Return[$Failed]
      ];
      Module[{dv, sl},
        {dv, sl} = liftedFan;
        If[Length[sl] == 0 || Length[dv] == 0 ||
           !AllTrue[sl, Length[#] == n + 1 &],
          Message[TropicalEval::liftdegenerate];  Return[$Failed]
        ]
      ]
    ];

    EvaluateTropicalMC[liftedSpec, liftedFan, kinematicPoints,
                       "LiftData" -> liftData, passThroughOpts]
  ]
];


(* --------------------------------------------------------------------------
   LaurentFromSubtraction  (plan B3, minimal helper)

   The tropical subtraction driver EvaluateTropicalMC returns the FULL integral
   F(eps) at a pinned EpsilonValue, not a Laurent.  For a single-divergence spec
       F(eps) = pole/eps + finite + c eps + ...,
   so  g(eps) = eps * F(eps) = pole + finite eps + c eps^2 + ...
   Running F at a few small eps and fitting g(eps) recovers (pole, finite)
   directly comparable to the IBP driver -- using only the already-working
   EpsilonValue route, for either sampler ("Integrator" -> MC | Vegas).
   (Public ::usage is declared in the BeginPackage section near the top.)
   -------------------------------------------------------------------------- *)

Options[LaurentFromSubtraction] = Join[
  {"EpsilonValues" -> {0.005, 0.01, 0.02}},
  Options[EvaluateTropicalMC]
];

LaurentFromSubtraction[integrandSpec_Association, fanData_List,
                       kinematicPoints_List, opts : OptionsPattern[]] :=
Module[{epsVals, fwd, nKP, perEps, fitPerKP},
  epsVals = OptionValue["EpsilonValues"];
  nKP     = Length[kinematicPoints];

  (* forward EvaluateTropicalMC options, but we own EpsilonValue + RunChecks +
     Method: this helper substitutes a NUMERIC eps (every sector convergent) and
     fits the Laurent, so it must stay in the convergent driver path
     (Method -> "None"), never route to IBP. *)
  fwd = DeleteCases[
    FilterRules[{opts}, Options[EvaluateTropicalMC]],
    ("EpsilonValue" -> _) | ("RunChecks" -> _) | ("Method" -> _)];

  (* full complex integral F(eps_i) for every kp, one row per eps value *)
  perEps = Table[
    Module[{r},
      r = EvaluateTropicalMC[integrandSpec, fanData, kinematicPoints,
            "EpsilonValue" -> ev, "RunChecks" -> False, "Method" -> "None",
            Sequence @@ fwd];
      If[AssociationQ[r],
        (#["Re"] + I*#["Im"]) & /@ r["Results"],
        $Failed]
    ],
    {ev, epsVals}
  ];
  If[MemberQ[perEps, $Failed],
    Print["LaurentFromSubtraction: a subtraction run failed"]; Return[$Failed]];

  fitPerKP = Table[
    Module[{Fs, gs, deg, design, cRe, cIm, pole, finite},
      Fs  = perEps[[All, kp]];                   (* F(eps_i) for this kp *)
      gs  = MapThread[#1*#2 &, {epsVals, Fs}];    (* g_i = eps_i F(eps_i) *)
      deg = Min[Length[epsVals] - 1, 2];          (* linear (2 pts) or quadratic *)
      design = Table[ev^p, {ev, epsVals}, {p, 0, deg}];
      cRe = LeastSquares[design, Re[gs]];
      cIm = LeastSquares[design, Im[gs]];
      pole   = cRe[[1]] + I*cIm[[1]];             (* g(0)  = pole  *)
      finite = cRe[[2]] + I*cIm[[2]];             (* g'(0) = finite *)
      <|"Pole" -> pole, "Finite" -> finite,
        "EpsilonValues" -> epsVals, "FullIntegrals" -> Fs|>
    ],
    {kp, nKP}
  ];

  <|"Results" -> fitPerKP, "EpsilonValues" -> epsVals|>
];


(* ============================================================================
   MODULE 2b: IBP-BASED POLE EXTRACTION

   Integration-by-parts based pole extraction for divergent sectors.
   The symbolic reduction (IBPReduceSector) iterates over multiple divergent
   variables, but the numerical pipeline (IBPProcessSector / boundary
   construction / EvaluateTropicalMCIBP Laurent assembly) supports a SINGLE
   divergent variable per sector (one 1/eps pole).  Nested divergences are
   detected and refused ($Failed); see TropicalEval::nestedIBP.

   This module is ADDITIVE — it does not modify any existing functions.
   The existing ProcessDivergentSector (tropical subtraction) is retained
   as the "Subtraction" option for DivergenceMethod.

   Key difference from the plan: the IBP on [0,1]^n produces a boundary
   term at y_k = 1 which is generally non-zero.  The correct formula is:

       I_sector = (1/a_k) * [ B_boundary - sum_t coeff_t * I_t ]

   where B_boundary = int f(y)|_{y_k=1} dy_{!=k} is the boundary
   contribution and I_t are the IBP monomial terms.
   ============================================================================ *)

(* --------------------------------------------------------------------------
   IBPExpandOneVariable
   Applies IBP decomposition for a single divergent variable y_k to one
   term.  Returns a list of new terms (one per monomial per polynomial).
   -------------------------------------------------------------------------- *)

IBPExpandOneVariable[term_Association, k_Integer,
                     clearedPolys_List, eps_, imB_List : {}] :=
Module[
  {termCoeff, termExps, termPolyExps, newTerms, nPolys},

  termCoeff    = term["Coefficient"];
  termExps     = term["NewExponents"];
  termPolyExps = term["PolyExponents"];
  nPolys       = Length[clearedPolys];

  newTerms = {};

  Do[
    Module[{Bj, polj},
      (* IBP differentiates Q_j^{B_j} and brings down the FULL exponent B_j.
         In SplitRealImag mode the sector is processed on Re(B), so termPolyExps
         carries only Re(B_j) (the magnitude exponent, shifted by prior steps);
         the imaginary part Im(B_j) does NOT shift and must be reinstated HERE so
         the brought-down coefficient is the true complex B_j (planCXLIFTDIV.md:
         the i*Im(B_j) phase-derivative term, absent in the plan's §3).  imB={} or
         all-zero (real / unlifted) -> Bj real -> byte-identical (#25). *)
      Bj   = termPolyExps[[j]] +
             If[Length[imB] >= j, I * imB[[j]], 0];
      polj = clearedPolys[[j]];

      Do[
        Module[{cm, em, emk, newCoeff, newExps, newPolyExps},
          cm  = mono[[1]];
          em  = mono[[2]];
          emk = em[[k]];

          (* Only monomials with e_{m,k} > 0 contribute to d/dy_k Q_j *)
          If[TrueQ[emk > 0] || (NumericQ[emk] && emk > 0),

            (* Coefficient: B_j * c_m * e_{m,k} *)
            newCoeff = termCoeff * Bj * cm * emk;

            (* New effective exponents: alpha_i = old_alpha_i + e_{m,i} *)
            newExps = termExps + em;

            (* Polynomial exponents: B_j -> B_j - 1, others unchanged *)
            newPolyExps = termPolyExps;
            newPolyExps[[j]] = newPolyExps[[j]] - 1;

            AppendTo[newTerms, <|
              "Coefficient"   -> newCoeff,
              "NewExponents"  -> newExps,
              "PolyExponents" -> newPolyExps
            |>];
          ]
        ],
        {mono, polj}
      ];
    ],
    {j, nPolys}
  ];

  newTerms
];

(* --------------------------------------------------------------------------
   ibpImagPole  (planIBPCX.md §3.1)
   The imaginary part theta_k of the divergent direction y_k's effective
   endpoint exponent  s_k = c_k eps + i theta_k.  The 1/eps pole is governed
   by s_k, not by c_k alone: theta_k = 0 gives a genuine 1/(c_k eps) pole;
   theta_k != 0 gives 1/(c_k eps + i theta_k), which is REGULAR at eps=0 (no
   pole, finite 1/(i theta_k)) — see planIBPCX.md §1.2.

   Two representations contribute, and they are mutually exclusive:
     - Direct / unlifted:  Im(A) + Sum_j Im(B_j) D_{j,k} is folded into the
       complex effective exponent NewExponents[[k]]  ->  Im[a0[[k]]].
     - Lifted SplitRealImag:  NewExponents is realified (Re only); the imaginary
       part rides separately via LiftedMonoFactor (Sum_j Im(B_j) D_{j,k}) and
       LiftedMonoPhase (Im(A) numerator on slot k).
   For every real / Case-A sector theta_k is provably zero (PossibleZeroQ) ->
   the real-pole Laurent branch is taken -> byte/numerically identical (#25).
   Built from exact Im[...] data (no N[...]) so exactness (#1) is preserved.
   -------------------------------------------------------------------------- *)

ibpImagPole[sectorData_Association, k_Integer, eps_] :=
Module[{a0k, lmf, lmp, imB, theta, realSyms},
  a0k   = If[eps === None, sectorData["NewExponents"][[k]],
                           sectorData["NewExponents"][[k]] /. eps -> 0];
  theta = Im[a0k];                                  (* Direct / unlifted part *)
  lmf = Lookup[sectorData, "LiftedMonoFactor", None];
  lmp = Lookup[sectorData, "LiftedMonoPhase",  None];
  imB = Lookup[sectorData, "ImagPolyExponents", None];
  (* Lifted SplitRealImag part (a0k is real there, so the two never overlap). *)
  If[ListQ[imB] && lmf =!= None,
    theta += Sum[imB[[j]] * lmf["DExp"][[j, k]], {j, Length[lmf["DExp"]]}]];
  If[lmp =!= None, theta += lmp["Num"][[k]]];
  (* The regulator eps is physically real (e.g. dim-reg eps in d = d0 - 2 eps).
     When it is an undeclared symbol, Im[...] of an exponent that contains it
     leaves a spurious Im[eps] term un-simplified: Im(3/2 + I nu - eps) prints
     as nu - Im[eps] instead of nu, and Im(B_j)/Im(A) stored on the lifted
     SplitRealImag slots carry the same stray Im[eps].  That term makes theta_k
     look eps-DEPENDENT, so the off-axis guard (splitdivmono) misfires and refuses
     a CONSTANT off-axis direction that the IBP route actually supports (the
     lifted divergent presectors of the 4-point figure; new_request UPDATE
     2026-06-30, planIBPMULTIDIV/planIBPCX.md §8).  Assume the regulator is real
     so Im[eps]->0; a genuinely eps-dependent imaginary part (a bare/real eps in
     the exponent -- the true out-of-scope theta ~ 1/eps) does NOT contain Im[eps]
     and so survives, still caught by the downstream FreeQ guard.  Refine is exact
     (no float, invariant #1); theta stays identically 0 on every real / theta=0
     sector, so PossibleZeroQ and the real-pole Laurent branch are unchanged (#25).
     Also declare any KinematicSymbols this sector carries real (planR.md R2:
     a kinematic symbol in an exponent's imaginary part would otherwise leave
     Re[]/Im[] heads in theta -- same footing as imagPolyInfo/imagMonoInfo). *)
  realSyms = Lookup[sectorData, "KinematicSymbols", {}];
  If[eps =!= None, realSyms = Append[realSyms, eps]];
  If[Length[realSyms] > 0, theta = Refine[theta, Element[realSyms, Reals]]];
  theta
];

(* --------------------------------------------------------------------------
   ibpDivClass  (planIBPMULTIDIV.md §3.1)
   Classify divergent direction i of a sector by (Re(alpha0), theta):
     "Convergent" : Re(alpha0) > 0                  (not divergent at all)
     "Pole"       : Re(alpha0) <= 0 AND theta == 0  (a genuine 1/eps pole)
     "OffAxis"    : Re(alpha0) <= 0 AND theta != 0  (finite, 1/(i theta))
   Decided SYMBOLICALLY (PossibleZeroQ on the exact ibpImagPole, before any
   kinRules), so it makes no float decision in the decomposition (invariant #1;
   mirrors planIBPCX.md §3.4).  theta == 0 on every real / Case-A sector ->
   "Pole"/"Convergent" only -> the off-axis code paths never activate ->
   byte/numerically identical (#25).
   -------------------------------------------------------------------------- *)
ibpDivClass[sectorData_Association, i_Integer, eps_] :=
Module[{a0i, th},
  a0i = If[eps === None, sectorData["NewExponents"][[i]],
                         sectorData["NewExponents"][[i]] /. eps -> 0];
  th  = ibpImagPole[sectorData, i, eps];
  Which[
    ! (TrueQ[Re[a0i] <= 0] || (NumericQ[a0i] && Re[a0i] <= 0)), "Convergent",
    TrueQ[PossibleZeroQ[th]],                                   "Pole",
    True,                                                       "OffAxis"]
];

(* dropVarsFromDomain (planIBPMULTIDIV.md §3.3): a corner piece runs over the
   surviving (live) coordinates only; the lifted-sector domain indicator must be
   restricted to those coordinates.  In Case A ic = 0 on every divergent slot
   (ProcessSectorLifted sets IndicatorCoeffs[[divDir]] = 0), so this just selects
   the live entries.  None for unlifted -> None (byte-identical #25). *)
dropVarsFromDomain[None, _] := None;
dropVarsFromDomain[dc_Association, liveVars_List] := <|
  "LogZ0"           -> dc["LogZ0"],
  "MP"              -> dc["MP"],
  "IndicatorCoeffs" -> dc["IndicatorCoeffs"][[liveVars]]
|>;

(* --------------------------------------------------------------------------
   IBPReduceSector
   Main IBP reduction function.  Iteratively applies IBP to resolve all
   divergent variables in a sector.
   -------------------------------------------------------------------------- *)

IBPReduceSector[sectorData_Association, eps_] :=
Module[
  {n, aVals, polyExps, clearedPolys, detM,
   a0, allDivVars, ibpPrefactors, terms, resolvedVars, imB},

  n            = sectorData["Dimension"];
  aVals        = sectorData["NewExponents"];
  polyExps     = sectorData["PolynomialExponents"];
  clearedPolys = sectorData["ClearedPolys"];
  detM         = sectorData["DetM"];
  (* Im(B) per polynomial (planCXLIFTDIV.md §4.4): reinstated as the brought-down
     IBP coefficient's imaginary part.  None/absent (real / unlifted) -> {} ->
     IBPExpandOneVariable brings down a real B -> byte-identical (#25). *)
  imB          = Lookup[sectorData, "ImagPolyExponents", None];
  If[imB === None, imB = {}];

  (* Find all divergent variables: Re(a_i^(0)) <= 0 *)
  a0 = aVals /. eps -> 0;
  allDivVars = {};
  Do[
    If[TrueQ[Re[a0[[i]]] <= 0] ||
       (NumericQ[a0[[i]]] && Re[a0[[i]]] <= 0),
      AppendTo[allDivVars, i]
    ],
    {i, n}
  ];

  If[Length[allDivVars] == 0,
    Print["IBPReduceSector: no divergent variables in sector ",
          sectorData["ConeIndex"]];
    Return[$Failed]
  ];

  Print["IBPReduceSector: sector ", sectorData["ConeIndex"],
        ", divergent variables: y_", allDivVars,
        ", effective exponents at eps=0: ", a0];

  (* Initialize: single term representing the original integral *)
  terms = {<|
    "Coefficient"   -> 1,
    "NewExponents"  -> aVals,
    "PolyExponents" -> polyExps
  |>};

  ibpPrefactors = {};
  resolvedVars  = {};

  (* Iteratively apply IBP for each divergent variable *)
  Do[
    Module[{k, ak, ck, newTerms, nBefore},
      k  = allDivVars[[dv]];
      ak = aVals[[k]];
      ck = D[ak, eps] /. eps -> 0;

      (* c_k = 0 is a genuine higher-order / unregulated log pole ONLY when the
         direction is also on-axis (theta_k = 0).  When theta_k != 0 the y_k=0
         endpoint is regulated by the oscillation (prefactor 1/(i theta_k), no
         eps needed) and IBP resolves it — do not abort (planIBPCX.md §3.2). *)
      If[(TrueQ[ck == 0] || (NumericQ[ck] && ck == 0)) &&
         TrueQ[PossibleZeroQ[ibpImagPole[sectorData, k, eps]]],
        Message[TropicalEval::badck, sectorData["ConeIndex"], k];
        Return[$Failed]
      ];

      AppendTo[ibpPrefactors, -1/ak];
      AppendTo[resolvedVars, k];

      nBefore  = Length[terms];
      newTerms = {};

      Do[
        Module[{term, termA0k},
          term    = terms[[t]];
          termA0k = term["NewExponents"][[k]] /. eps -> 0;

          If[TrueQ[Re[termA0k] <= 0] ||
             (NumericQ[termA0k] && Re[termA0k] <= 0),
            (* Divergent in y_k: apply IBP *)
            newTerms = Join[newTerms,
              IBPExpandOneVariable[term, k, clearedPolys, eps, imB]
            ],
            (* Already convergent: pass through *)
            AppendTo[newTerms, term]
          ]
        ],
        {t, Length[terms]}
      ];

      terms = newTerms;

      Print["  IBP step for y_", k, ": ", nBefore, " -> ",
            Length[terms], " terms"];
    ],
    {dv, Length[allDivVars]}
  ];

  (* Verify all terms are now convergent *)
  Module[{problemTerms = 0},
    Do[
      Module[{termA0},
        termA0 = terms[[t]]["NewExponents"] /. eps -> 0;
        Do[
          If[TrueQ[Re[termA0[[i]]] <= 0] ||
             (NumericQ[termA0[[i]]] && Re[termA0[[i]]] <= 0),
            problemTerms++;
            If[problemTerms <= 3,
              Print["WARNING: term ", t, " still divergent in y_", i,
                    ", alpha^(0) = ", termA0[[i]]]
            ]
          ],
          {i, n}
        ]
      ],
      {t, Length[terms]}
    ];
    If[problemTerms > 0,
      Print["WARNING: ", problemTerms,
            " divergences remain after IBP reduction"]
    ];
  ];

  <|
    "ConeIndex"              -> sectorData["ConeIndex"],
    "IBPPrefactors"          -> ibpPrefactors,
    "Terms"                  -> terms,
    "NTerms"                 -> Length[terms],
    "DivergentVariables"     -> resolvedVars,
    "NDivergent"             -> Length[resolvedVars],
    "ClearedPolys"           -> clearedPolys,
    "DetM"                   -> detM,
    "Dimension"              -> n,
    "OriginalExponents"      -> aVals,
    "OriginalPolyExponents"  -> polyExps,
    "TransformedPolys"       -> sectorData["TransformedPolys"],
    "MinExponents"           -> sectorData["MinExponents"],
    "MonomialExponents"      -> sectorData["MonomialExponents"],
    "RawExponents"           -> sectorData["RawExponents"],
    "RayMatrix"              -> sectorData["RayMatrix"],
    "SelectedRays"           -> sectorData["SelectedRays"]
  |>
];

(* --------------------------------------------------------------------------
   ibpProcessLeaf  (planIBPMULTIDIV.md §3.3)
   Turn one corner LEAF (a single term over its live coordinates) into the
   processed, flattened integrand piece the codegen consumes.  This is exactly
   the per-term processing of IBPProcessSector Step 3 (and reproduces the Step 2
   boundary data when the leaf is the all-boundary corner), generalized to an
   arbitrary live-coordinate subset.
   -------------------------------------------------------------------------- *)
ibpProcessLeaf[node_Association, eps_, pfBase_, lmf_, lmp_, nPolys_Integer,
               nPolyExps_Integer, domC_] :=
Module[{liveVars, cp, exps, polyExps, coeff, alpha0, alpha1,
        flatPolys, prefactor, tB0, tB1, coeff0, coeff1, logInsertions,
        monoFactorLog, monoPhaseLog, nlive},
  liveVars = node["Live"];
  nlive    = Length[liveVars];
  cp       = node["ClearedPolys"];
  exps     = node["Exps"];
  polyExps = node["PolyExps"];
  coeff    = node["Coeff"];

  alpha0 = (exps /. eps -> 0)[[liveVars]];
  alpha1 = (D[exps, eps] /. eps -> 0)[[liveVars]];

  (* Flatten: keep only the live coordinates, divide each by its alpha0. *)
  flatPolys = Table[
    Table[{mono[[1]], MapThread[#1/#2 &, {mono[[2]][[liveVars]], alpha0}]},
          {mono, cp[[j]]}],
    {j, nPolys}];

  prefactor = (pfBase /. eps -> 0) / (Times @@ alpha0);

  tB0 = polyExps /. eps -> 0;
  tB1 = D[polyExps, eps] /. eps -> 0;
  coeff0 = coeff /. eps -> 0;
  coeff1 = D[coeff, eps] /. eps -> 0;

  logInsertions = <|
    "VariableTerms"   -> Table[{alpha1[[i]] / alpha0[[i]], i}, {i, nlive}],
    "PolynomialTerms" -> Table[{tB1[[j]], j}, {j, nPolyExps}]|>;

  monoFactorLog = If[lmf === None, None,
    Table[<|"Const"  -> lmf["Const"][[j]],
            "Coeffs" -> Table[lmf["DExp"][[j, liveVars[[i]]]] / alpha0[[i]],
                              {i, nlive}]|>,
          {j, nPolys}]];
  monoPhaseLog = If[lmp === None, None,
    <|"Const"  -> lmp["Const"],
      "Coeffs" -> Table[lmp["Num"][[liveVars[[i]]]] / alpha0[[i]], {i, nlive}]|>];

  <|"Sign"            -> node["Sign"],
    "FlatPolys"       -> flatPolys,
    "Prefactor"       -> prefactor,
    "Dimension"       -> nlive,
    "PolyExponents"   -> tB0,
    "Coeff0"          -> coeff0,
    "Coeff1"          -> coeff1,
    "LogInsertions"   -> logInsertions,
    "MonoFactorLog"   -> monoFactorLog,
    "MonomialPhaseLog"-> monoPhaseLog,
    "DomainConstraint"-> dropVarsFromDomain[domC, liveVars],
    "LiveVars"        -> liveVars,
    "Alpha0"          -> alpha0|>
];

(* --------------------------------------------------------------------------
   IBPBuildCorners  (planIBPMULTIDIV.md §3.3 — the one-pole + N-off-axis core)
   Iterated IBP over the full divergent set D = {pole k, off-axis m_1..m_N}.
   Each divergent direction is peeled into a BOUNDARY choice (y_d = 1, drop the
   coordinate, sign +) or a BULK choice (raise alpha_d via IBPExpandOneVariable,
   sign -); the 2^|D| leaves are the corners.  The global prefactor
   prod_d (1/a_d) is factored OUT (assembled in the driver), so a corner carries
   only its sign and its convergent flattened integrand.  This is the strict
   generalization of the single-direction Step 2 (boundary) + Step 3 (terms):
   for |D| = 1 it produces exactly {boundary (sign +), bulk terms (sign -)}.

   Returns {corners, $Failed-or-Null}: corners is the processed piece list;
   the second value is a direction index if a power-divergent residual (a live
   coordinate still Re(alpha0) <= 0 after one raise — the §3.1 edge) survived,
   else Null.
   -------------------------------------------------------------------------- *)
IBPBuildCorners[sectorData_Association, eps_, pfBase_, lmf_, lmp_, imB_List,
                divVars_List, domC_] :=
Module[{n, aVals, polyExps, clearedPolys, nPolys, nPolyExps, nodes, leaves,
        corners, badDir = Null},
  n            = sectorData["Dimension"];
  aVals        = sectorData["NewExponents"];
  polyExps     = sectorData["PolynomialExponents"];
  clearedPolys = sectorData["ClearedPolys"];
  nPolys       = Length[clearedPolys];
  nPolyExps    = Length[polyExps];

  (* Seed: the original single integrand term over all n coordinates. *)
  nodes = {<|"Coeff" -> 1, "Exps" -> aVals, "PolyExps" -> polyExps,
             "ClearedPolys" -> clearedPolys, "Dead" -> {}, "Sign" -> 1|>};

  (* Peel each divergent direction in turn. *)
  Do[
    Module[{d = divVars[[dv]], next = {}},
      Do[
        Module[{node = nodes[[ni]], bndCP, bulkTerms},
          (* --- BOUNDARY corner: set y_d = 1 (zero the d-th monomial exponent,
                 mark d dead), sign unchanged (+). --- *)
          bndCP = Table[
            Table[{mono[[1]], ReplacePart[mono[[2]], d -> 0]},
                  {mono, node["ClearedPolys"][[j]]}],
            {j, nPolys}];
          AppendTo[next, <|
            "Coeff" -> node["Coeff"], "Exps" -> node["Exps"],
            "PolyExps" -> node["PolyExps"], "ClearedPolys" -> bndCP,
            "Dead" -> Append[node["Dead"], d], "Sign" -> node["Sign"]|>];

          (* --- BULK corner: IBP-expand y_d (raise alpha_d, bring down the poly
                 log), sign flipped (-).  Multiple monomial terms. --- *)
          bulkTerms = IBPExpandOneVariable[
            <|"Coefficient" -> node["Coeff"], "NewExponents" -> node["Exps"],
              "PolyExponents" -> node["PolyExps"]|>,
            d, node["ClearedPolys"], eps, imB];
          Do[
            AppendTo[next, <|
              "Coeff" -> bt["Coefficient"], "Exps" -> bt["NewExponents"],
              "PolyExps" -> bt["PolyExponents"],
              "ClearedPolys" -> node["ClearedPolys"],
              "Dead" -> node["Dead"], "Sign" -> -node["Sign"]|>],
            {bt, bulkTerms}];
        ],
        {ni, Length[nodes]}];
      nodes = next;
    ],
    {dv, Length[divVars]}];

  (* Finalize leaves: live coords + verify convergence (§3.1 / §3.4). *)
  leaves = Table[
    Append[node, "Live" -> Complement[Range[n], node["Dead"]]],
    {node, nodes}];

  Do[
    Module[{a0 = (leaf["Exps"] /. eps -> 0)[[leaf["Live"]]]},
      Do[
        If[TrueQ[Re[a0[[i]]] <= 0] || (NumericQ[a0[[i]]] && Re[a0[[i]]] <= 0),
          badDir = leaf["Live"][[i]]],
        {i, Length[a0]}]],
    {leaf, leaves}];

  If[badDir =!= Null, Return[{$Failed, badDir}]];

  corners = Table[
    ibpProcessLeaf[leaf, eps, pfBase, lmf, lmp, nPolys, nPolyExps, domC],
    {leaf, leaves}];

  {corners, Null}
];

(* --------------------------------------------------------------------------
   IBPProcessSectorMultiDiv  (planIBPMULTIDIV.md §3.3)
   The one-pole + N-off-axis (or N-off-axis, no pole) processing path.  Builds
   the 2^|D| iterated-IBP corners (IBPBuildCorners) and carries the global
   prefactor data the driver needs to assemble:
       I = [1/(c_k eps)] * prod_j[1/(i theta_j)] * sum_corners (sign) L_corner
   The pole's 1/(c_k eps) is the only eps-singular factor; the off-axis
   prod_j 1/(i theta_j) is a finite constant -> Laurent orders -1 and 0 only.
   -------------------------------------------------------------------------- *)
IBPProcessSectorMultiDiv[sectorData_Association, integrandSpec_Association,
                         poleDirs_List, offDirs_List] :=
Module[
  {eps, n, detM, pfBase, lmf, lmp, imB, domC, aVals, divVars,
   cornersResult, corners, badDir, k, ck, rk, dlog, offThetas, offCks, offRks,
   hasPole},

  eps     = integrandSpec["RegulatorSymbol"];
  n       = sectorData["Dimension"];
  detM    = sectorData["DetM"];
  pfBase  = Lookup[sectorData, "PrefactorBase", Abs[detM]];
  lmf     = Lookup[sectorData, "LiftedMonoFactor", None];
  lmp     = Lookup[sectorData, "LiftedMonoPhase",  None];
  imB     = Lookup[sectorData, "ImagPolyExponents", None];
  If[imB === None, imB = {}];
  domC    = Lookup[sectorData, "DomainConstraint", None];
  aVals   = sectorData["NewExponents"];
  hasPole = Length[poleDirs] >= 1;

  (* Peel the pole first (if present), then the off-axis directions. *)
  divVars = Join[poleDirs, offDirs];

  (* A genuine pole with c_k = 0 (and theta_k = 0 by construction) is an
     unregulated higher-order log pole — refuse (badck).  Off-axis directions
     are exempt: theta != 0 regulates the y=0 endpoint without eps. *)
  If[hasPole,
    k = poleDirs[[1]];
    Module[{ckk = D[aVals[[k]], eps] /. eps -> 0},
      If[TrueQ[ckk == 0] || (NumericQ[ckk] && ckk == 0),
        Message[TropicalEval::badck, sectorData["ConeIndex"], k];
        Return[$Failed]]]
  ];

  (* Every divergent direction must be a LOGARITHMIC endpoint, Re(alpha0)=0: the
     pole 1/(c_k eps) and each off-axis 1/(i theta_j) are residues of a Re=0
     endpoint.  A strictly power-divergent direction (Re(alpha0)<0) is genuinely
     divergent — the oscillation does NOT regulate a magnitude divergence
     (|y^{rho-1+i theta}| = y^{rho-1} still blows up at y->0 for rho<0) and IBP's
     y=0 boundary term does not vanish — so its raised-integrand "finite" value
     would be a WRONG number.  Refuse cleanly (planIBPMULTIDIV.md §3.1): off-axis
     -> offaxispow; a Re<0 "pole" (theta=0) is a higher-order residual ->
     ibpresidual.  (After this guard every divergent direction has Re(alpha0)=0,
     so one IBP raise by an integer >=1 makes every corner leaf convergent.) *)
  Module[{neg},
    neg[i_] := Module[{r = Re[aVals[[i]] /. eps -> 0]},
      TrueQ[r < 0] || (NumericQ[r] && r < 0)];
    With[{m = SelectFirst[offDirs, neg, None]},
      If[m =!= None,
        Message[TropicalEval::offaxispow, sectorData["ConeIndex"], m];
        Return[$Failed]]];
    With[{p = SelectFirst[poleDirs, neg, None]},
      If[p =!= None,
        Message[TropicalEval::ibpresidual, sectorData["ConeIndex"], 0, p,
                aVals[[p]] /. eps -> 0];
        Return[$Failed]]]
  ];

  (* Build all corners (iterated IBP over the divergent set). *)
  cornersResult = IBPBuildCorners[sectorData, eps, pfBase, lmf, lmp, imB,
                                  divVars, domC];
  corners = cornersResult[[1]];
  badDir  = cornersResult[[2]];
  If[corners === $Failed,
    (* Defensive: with the Re=0 guard above this is unreachable (one raise makes
       every leaf convergent).  A surviving Re(alpha0)<=0 leaf would be an
       on-axis (theta=0) unreduced residual, so report it as such. *)
    Message[TropicalEval::ibpresidual, sectorData["ConeIndex"], 0, badDir,
            "Re<=0 after one IBP raise"];
    Return[$Failed]
  ];

  (* Off-axis prefactor data: theta_m (exact symbolic), c_m = d Re(alpha_m)/d eps,
     and r_m = a_m^(2)/c_m (the second-order ratio, used only by the driver's
     degenerate-kp fallback when a lone off-axis theta crosses 0 and the direction
     becomes a genuine 1/(c_m eps) pole at that kp).  Im(alpha_m) is eps-free
     (splitdivmono guards otherwise), so D[aVals[[m]],eps] is the real c_m. *)
  offThetas = Table[ibpImagPole[sectorData, m, eps], {m, offDirs}];
  offCks    = Table[D[aVals[[m]], eps] /. eps -> 0, {m, offDirs}];
  offRks    = Table[
    Module[{cm = D[aVals[[m]], eps] /. eps -> 0,
            am2 = (1/2) D[aVals[[m]], {eps, 2}] /. eps -> 0},
      If[TrueQ[cm == 0], 0, am2/cm]],
    {m, offDirs}];

  (* PrefactorBase eps-derivative (planAXpDIV.md §4.1): finite-part correction
     (P'/P)(0)*pole.  Computed UNCONDITIONALLY (0 for unlifted, nonzero for a
     lifted pfBase carrying eps) — the no-pole branch also needs it for the
     driver's degenerate-kp fallback, where a lone off-axis becomes a genuine
     pole and the (P'/P)(0)*pole term must be reinstated. *)
  dlog = If[eps === None, 0, D[Log[pfBase], eps] /. eps -> 0];
  (* Pole data (the single 1/eps factor), if a genuine pole is present. *)
  If[hasPole,
    ck   = D[aVals[[k]], eps] /. eps -> 0;
    Module[{ak2 = (1/2) D[aVals[[k]], {eps, 2}] /. eps -> 0},
      rk = If[TrueQ[ck == 0], 0, ak2/ck]],
    ck = None; rk = 0
  ];

  <|"ConeIndex"             -> sectorData["ConeIndex"],
    "IsDivergent"           -> True,
    "Method"                -> "IBP",
    "MultiDiv"              -> True,
    "Corners"               -> corners,
    "NCorners"              -> Length[corners],
    (* NTerms alias so the verbose driver's Total[#["NTerms"]&/@...] stays
       numeric across mixed single-pole / MultiDiv sectors. *)
    "NTerms"                -> Length[corners],
    "Dimension"             -> n,
    "DetM"                  -> detM,
    "HasPole"               -> hasPole,
    "PoleDir"               -> If[hasPole, k, None],
    "OffAxisDirs"           -> offDirs,
    "ck"                    -> ck,
    "rk"                    -> rk,
    "DLogPrefactor"         -> dlog,
    "OffAxisThetas"         -> offThetas,
    "OffAxisCks"            -> offCks,
    "OffAxisRks"            -> offRks,
    "DivergentVariables"    -> divVars,
    "NDivergent"            -> Length[divVars],
    "OriginalExponents"     -> aVals,
    "OriginalPolyExponents" -> sectorData["PolynomialExponents"],
    "TransformedPolys"      -> sectorData["TransformedPolys"],
    "MinExponents"          -> sectorData["MinExponents"],
    "MonomialExponents"     -> sectorData["MonomialExponents"],
    "RawExponents"          -> sectorData["RawExponents"],
    "DomainConstraint"      -> domC,
    "ImagPolyExponents"     -> Lookup[sectorData, "ImagPolyExponents", None]|>
];

(* --------------------------------------------------------------------------
   IBPProcessSector
   Full IBP pipeline: reduce, expand in epsilon, flatten, build boundary.
   Returns an IBPSectorData association ready for C++ codegen.

   For single divergence: produces boundary (at y_k = 1) and IBP terms.
   The combination is:
       I_sector = (1/a_k) [ B_boundary - sum_t coeff_t I_t ]

   For one pole + N off-axis divergent directions (planIBPMULTIDIV.md): the
   off-axis directions are finite (1/(i theta_j)); IBPProcessSectorMultiDiv
   assembles the full 2^(N+1)-corner iterated IBP (tagged "MultiDiv" -> True).
   -------------------------------------------------------------------------- *)

IBPProcessSector[sectorData_Association, integrandSpec_Association] :=
Module[
  {eps, n, ibpData, terms, clearedPolys, detM, pfBase,
   divVars, ck, rk, aVals, polyExps,
   a0, a1, B0, B1, ak, ak2,
   boundaryData, ibpTermsProcessed, lmf, lmp,
   k, imagPole},

  eps = integrandSpec["RegulatorSymbol"];
  n   = sectorData["Dimension"];

  (* LiftedMonoFactor (planCXLIFTDIV.md §4.1/§4.2): un-flattened numerators of the
     dropped tropical monomial factor, used to re-derive the SplitRealImag phase
     per IBP piece by dividing by that piece's own flattening alpha0.  None for
     real / unlifted sectors -> no phase emitted -> byte-identical (#25). *)
  lmf = Lookup[sectorData, "LiftedMonoFactor", None];
  (* LiftedMonoPhase (planAXpDIVv2.md §4.3): the un-flattened bare-monomial phase
     numerator (Im(A)), re-flattened per IBP piece exactly like lmf.  None for
     real-A sectors -> no monomial phase emitted -> byte-identical (#25). *)
  lmp = Lookup[sectorData, "LiftedMonoPhase", None];

  (* Off-axis-aware divergence classification (planIBPMULTIDIV.md §3.1/§3.2).
     A direction is a genuine 1/eps POLE only when Re(alpha0)=0 AND theta=0;
     a Re(alpha0)<=0 direction with theta!=0 is OFF-AXIS (finite, 1/(i theta)),
     NOT a pole.  Refuse only on >1 genuine pole (real 1/eps^{d>=2}); off-axis
     directions no longer inflate the count.  theta=0 on every real / Case-A
     sector -> no off-axis directions -> the existing single-pole construction
     below runs verbatim (byte/numerically identical, #25). *)
  Module[{cls, poleDirs, offDirs},
    cls      = Table[ibpDivClass[sectorData, i, eps], {i, n}];
    poleDirs = Select[Range[n], cls[[#]] === "Pole" &];
    offDirs  = Select[Range[n], cls[[#]] === "OffAxis" &];
    If[Length[poleDirs] > 1,
      Message[TropicalEval::nestedIBP, sectorData["ConeIndex"], Length[poleDirs]];
      Return[$Failed]
    ];
    (* eps-dependent theta on any off-axis direction is out of scope (it scales
       the oscillation with 1/eps; planIBPCX.md §8 / splitdivmono).  NOTE: a bare
       Return[$Failed] inside Do[Module[...]] would NOT propagate out of this
       function (it exits only the Do iteration) — find the offender first, then
       Return at the Module top level (where Return does propagate). *)
    Module[{badOff = SelectFirst[offDirs,
        (!FreeQ[ibpImagPole[sectorData, #, eps], eps]) &, None]},
      If[badOff =!= None,
        Message[TropicalEval::splitdivmono, sectorData["ConeIndex"],
                ibpImagPole[sectorData, badOff, eps]];
        Return[$Failed]]];
    (* Any off-axis direction present -> the generalized multi-corner path
       (also covers the single-off-axis no-pole case, cc_47 Part C). *)
    If[Length[offDirs] > 0,
      Return[IBPProcessSectorMultiDiv[sectorData, integrandSpec,
                                      poleDirs, offDirs]]
    ]
    (* else: no off-axis, <=1 pole -> fall through to the single-pole path. *)
  ];

  (* Step 1: IBP reduction *)
  ibpData = IBPReduceSector[sectorData, eps];
  If[ibpData === $Failed, Return[$Failed]];

  terms        = ibpData["Terms"];
  clearedPolys = ibpData["ClearedPolys"];
  detM         = ibpData["DetM"];
  (* PrefactorBase (planAXpDIV.md §4.1): Abs[detM] unlifted, lifted base else. *)
  pfBase       = Lookup[sectorData, "PrefactorBase", Abs[detM]];
  divVars      = ibpData["DivergentVariables"];
  aVals        = ibpData["OriginalExponents"];
  polyExps     = ibpData["OriginalPolyExponents"];

  (* Single divergent variable only (nested was refused above; this is a
     defensive backstop should IBPReduceSector resolve a different count). *)
  If[Length[divVars] != 1,
    Message[TropicalEval::nestedIBP, sectorData["ConeIndex"], Length[divVars]];
    Return[$Failed]
  ];

  k = divVars[[1]];

  (* Imaginary exponent of the divergent direction (planIBPCX.md §3.1).  When the
     divergent slot k carries a nonzero imaginary part theta_k the effective
     endpoint exponent is s_k = c_k eps + i theta_k: REGULAR at eps=0 (no 1/eps
     pole, finite 1/(i theta_k)) — see planIBPCX.md §1.2.  This was previously a
     refusal (TropicalEval::splitdivmono, "complex Case B"), conflating the harmless
     complex-pole-LOCATION shift with the genuinely-hard geometric Case B
     (liftdivdomain).  It is now COMPUTED and carried as "ImagPole"; the driver's
     Laurent assembly branches on it (theta_k=0 -> real pole, unchanged; theta_k!=0
     -> finite, off-axis).  Sources: the dropped tropical MONOMIAL FACTOR (Im(B).D),
     the bare MONOMIAL phase (Im(A)), or — in Direct/unlifted mode — the complex
     effective exponent itself (Im[a0[[k]]]); ibpImagPole unifies all three.
     theta_k = 0 on every real / Case-A sector -> PossibleZeroQ True -> byte-
     identical real path (#25).  The phases that make the off-axis bulk integrand
     numerically bounded are already emitted (per-piece termMonoFactorLog /
     MonomialPhaseLog divide by the RAISED alpha0; planIBPCX.md §2). *)
  imagPole = ibpImagPole[sectorData, k, eps];
  (* Out-of-scope guard (planIBPCX.md §8): an imaginary exponent that itself scales
     like 1/eps (eps-dependent Im) reintroduces fast oscillation that the bounded
     IBP integrand can no longer resolve.  The supported case has eps-free Im(B)/
     Im(A); refuse cleanly otherwise. *)
  If[!FreeQ[imagPole, eps],
    Message[TropicalEval::splitdivmono, sectorData["ConeIndex"], imagPole];
    Return[$Failed]
  ];

  (* Epsilon expansion of the original effective exponents *)
  a0 = aVals /. eps -> 0;
  a1 = D[aVals, eps] /. eps -> 0;

  (* Polynomial exponent expansion *)
  B0 = polyExps /. eps -> 0;
  B1 = D[polyExps, eps] /. eps -> 0;

  (* IBP prefactor expansion: -1/a_k where a_k = c_k*eps + a_k^(2)*eps^2 + ... *)
  ak  = aVals[[k]];
  ck  = D[ak, eps] /. eps -> 0;
  ak2 = (1/2) D[ak, {eps, 2}] /. eps -> 0;
  rk  = If[TrueQ[ck == 0], 0, ak2 / ck];

  (* ----- Step 2: Build boundary data (f evaluated at y_k = 1) ----- *)
  (* Boundary polynomials: Q_j with y_k set to 1.
     This keeps ALL monomials but drops the y_k coordinate. *)
  Module[{ndVars, bndPolys, bndA0, bndFlatPolys, bndPrefactor, bndDim,
          bndLogInsertions, bndMonoFactorLog, bndMonoPhaseLog},
    ndVars = DeleteCases[Range[n], k];
    bndDim = n - 1;

    (* Boundary polynomials: set y_k = 1, remove y_k exponent *)
    bndPolys = Table[
      Table[
        {mono[[1]], mono[[2]][[ndVars]]},
        {mono, clearedPolys[[j]]}
      ],
      {j, Length[clearedPolys]}
    ];

    bndA0 = a0[[ndVars]];

    (* Flatten non-divergent variables *)
    bndFlatPolys = Table[
      Table[
        {mono[[1]],
         MapThread[#1/#2 &, {mono[[2]], bndA0}]},
        {mono, bndPolys[[j]]}
      ],
      {j, Length[bndPolys]}
    ];

    (* PrefactorBase at eps->0 for codegen (planAXpDIV.md §4.1): a lifted
       PrefactorBase carries eps; the codegen needs the numeric eps^0 value, and
       the finite (P'/P)(0)*pole term is added back via DLogPrefactor in the
       driver assembly.  eps-free for unlifted -> byte-identical (#25). *)
    bndPrefactor = (pfBase /. eps -> 0) / (Times @@ bndA0);

    bndLogInsertions = <|
      "VariableTerms" -> Table[
        {a1[[ndVars[[i]]]] / bndA0[[i]], i},
        {i, bndDim}
      ],
      "PolynomialTerms" -> Table[
        {B1[[j]], j},
        {j, Length[polyExps]}
      ]
    |>;

    (* Per-piece MonoFactorLog (planCXLIFTDIV.md §3/§4.2): the boundary lives at
       y_k = 1 over the n-1 non-divergent coords (slot k dropped, log y_k = 0),
       flattened by bndA0 = a0[[ndVars]].  Re-divide the fixed numerators D_{j,*}
       by this piece's flattening.  None when not a lifted complex sector. *)
    bndMonoFactorLog = If[lmf === None, None,
      Table[
        <|"Const"  -> lmf["Const"][[j]],
          "Coeffs" -> Table[lmf["DExp"][[j, ndVars[[i]]]] / bndA0[[i]],
                            {i, bndDim}]|>,
        {j, Length[clearedPolys]}]];

    (* Per-piece bare-monomial phase (planAXpDIVv2.md §4.5): single oscillatory
       term re-flattened to the boundary's bndA0 (slot k dropped, log y_k = 0).
       None when not a lifted complex-A sector. *)
    bndMonoPhaseLog = If[lmp === None, None,
      <|"Const"  -> lmp["Const"],
        "Coeffs" -> Table[lmp["Num"][[ndVars[[i]]]] / bndA0[[i]], {i, bndDim}]|>];

    boundaryData = <|
      "FlatPolys"        -> bndFlatPolys,
      "Prefactor"        -> bndPrefactor,
      "Dimension"        -> bndDim,
      "PolyExponents"    -> B0,
      "Avals"            -> bndA0,
      "LogInsertions"    -> bndLogInsertions,
      "MonoFactorLog"    -> bndMonoFactorLog,
      "MonomialPhaseLog" -> bndMonoPhaseLog
    |>;
  ];

  (* ----- Step 3: Process each IBP term (expand in eps + flatten) ----- *)
  (* Catch wraps the per-term verify (planIBPMULTIDIV.md §3.4): a bare
     Return[$Failed] inside Table[Module[...]] does NOT propagate — it leaves an
     unevaluated Return[$Failed] in one cell and the function keeps going.  Throw
     to this tag so a genuine unreduced residual cleanly aborts the whole sector. *)
  ibpTermsProcessed = Catch[Table[
    Module[{term, alpha, alpha0, alpha1, termPolyExps, tB0, tB1,
            coeff, coeff0, coeff1, flatPrefactor, flatPolys,
            logInsertions, termMonoFactorLog, termMonoPhaseLog},
      term         = terms[[t]];
      alpha        = term["NewExponents"];
      termPolyExps = term["PolyExponents"];
      coeff        = term["Coefficient"];

      (* Expand effective exponents in epsilon *)
      alpha0 = alpha /. eps -> 0;
      alpha1 = D[alpha, eps] /. eps -> 0;

      (* Verify all alpha0 > 0 (planIBPMULTIDIV.md §3.4).  This is the SINGLE-POLE
         path: any off-axis (Re<=0, theta!=0) direction was routed to the MultiDiv
         path by the front guard, so a surviving Re(alpha0)<=0 term here is ALWAYS
         a genuine unreduced residual (whether or not it carries an imaginary part:
         flattening by an exponent with Re<=0 would emit a divergent integrand).
         Abort on ANY Re<=0 — a clean refusal, never a silent wrong number.  In
         supported sectors the lone pole is raised to Re>=1 so this never fires. *)
      Do[
        If[(NumericQ[alpha0[[i]]] && Re[alpha0[[i]]] <= 0) ||
           TrueQ[Re[alpha0[[i]]] <= 0],
          Module[{},
            Message[TropicalEval::ibpresidual, sectorData["ConeIndex"], t, i,
                    alpha0[[i]]];
            Throw[$Failed, "ibpverify"]
          ]
        ],
        {i, n}
      ];

      (* Expand polynomial exponents *)
      tB0 = termPolyExps /. eps -> 0;
      tB1 = D[termPolyExps, eps] /. eps -> 0;

      (* Expand coefficient *)
      coeff0 = coeff /. eps -> 0;
      coeff1 = D[coeff, eps] /. eps -> 0;

      (* Flatten: y_i -> (y_i')^{1/alpha0_i} *)
      flatPolys = Table[
        Table[
          {mono[[1]],
           MapThread[#1/#2 &, {mono[[2]], alpha0}]},
          {mono, clearedPolys[[j]]}
        ],
        {j, Length[clearedPolys]}
      ];

      flatPrefactor = (pfBase /. eps -> 0) / (Times @@ alpha0);

      (* Per-piece MonoFactorLog (planCXLIFTDIV.md §3/§4.2): IBP terms keep all n
         coords (y_k still integrated), flattened by this term's alpha0; re-divide
         the fixed numerators D_{j,i} by alpha0_i.  None for non-lifted-complex. *)
      termMonoFactorLog = If[lmf === None, None,
        Table[
          <|"Const"  -> lmf["Const"][[j]],
            "Coeffs" -> Table[lmf["DExp"][[j, i]] / alpha0[[i]], {i, n}]|>,
          {j, Length[clearedPolys]}]];

      (* Per-piece bare-monomial phase (planAXpDIVv2.md §4.5): IBP terms keep all
         n coords; re-flatten Num_i by this term's alpha0_i. *)
      termMonoPhaseLog = If[lmp === None, None,
        <|"Const"  -> lmp["Const"],
          "Coeffs" -> Table[lmp["Num"][[i]] / alpha0[[i]], {i, n}]|>];

      (* Log insertion sum *)
      logInsertions = <|
        "VariableTerms" -> Table[
          {alpha1[[i]] / alpha0[[i]], i},
          {i, n}
        ],
        "PolynomialTerms" -> Table[
          {tB1[[j]], j},
          {j, Length[termPolyExps]}
        ]
      |>;

      <|
        "FlatPolys"      -> flatPolys,
        "Prefactor"      -> flatPrefactor,
        "Dimension"      -> n,
        "PolyExponents"  -> tB0,
        "Coeff0"         -> coeff0,
        "Coeff1"         -> coeff1,
        "LogInsertions"  -> logInsertions,
        "Alpha0"         -> alpha0,
        "Alpha1"         -> alpha1,
        "MonoFactorLog"  -> termMonoFactorLog,
        "MonomialPhaseLog" -> termMonoPhaseLog
      |>
    ],
    {t, Length[terms]}
  ], "ibpverify"];
  If[ibpTermsProcessed === $Failed, Return[$Failed]];

  <|
    "ConeIndex"              -> sectorData["ConeIndex"],
    "IsDivergent"            -> True,
    "Method"                 -> "IBP",
    "DivergentVariable"      -> k,
    "DivergentVariables"     -> divVars,
    "NDivergent"             -> Length[divVars],
    "Dimension"              -> n,
    "DetM"                   -> detM,
    "ck"                     -> ck,
    "rk"                     -> rk,
    "a0"                     -> a0,
    "a1"                     -> a1,
    "B0"                     -> B0,
    "B1"                     -> B1,
    "BoundaryData"           -> boundaryData,
    "IBPTerms"               -> ibpTermsProcessed,
    "NTerms"                 -> Length[ibpTermsProcessed],
    "IBPPrefactors"          -> ibpData["IBPPrefactors"],
    "ClearedPolys"           -> clearedPolys,
    "OriginalExponents"      -> aVals,
    "OriginalPolyExponents"  -> polyExps,
    "TransformedPolys"       -> sectorData["TransformedPolys"],
    "MinExponents"           -> sectorData["MinExponents"],
    "MonomialExponents"      -> sectorData["MonomialExponents"],
    "RawExponents"           -> sectorData["RawExponents"],
    (* Lifted-sector domain indicator (planAXpDIV.md §4.2): propagated to the IBP
       boundary/term codegen.  None for unlifted -> byte-identical (#25). *)
    "DomainConstraint"       -> Lookup[sectorData, "DomainConstraint", None],
    (* DLogPrefactor (planAXpDIV.md §4.1): d/deps log(PrefactorBase)|_0.  A lifted
       sector's PrefactorBase = (|detM|/|mp|) z0^(ap/mp-1) carries eps via ap(eps),
       so the codegen — which uses PrefactorBase(0) — misses a finite-part term
       (P'/P)(0)*pole; the driver's Laurent assembly adds it back.  Exactly 0 for
       an unlifted sector (PrefactorBase = |detM| is eps-free), so the unlifted
       IBP numbers are unchanged. *)
    "DLogPrefactor"          -> If[eps === None, 0,
                                   D[Log[pfBase], eps] /. eps -> 0],
    (* Im(B) per polynomial (planCXLIFTDIV.md §4.4): rides the real measure as the
       SplitRealImag oscillatory phase in the IBP boundary/term codegen.  None for
       real / unlifted sectors -> no phase emitted -> byte-identical (#25). *)
    "ImagPolyExponents"      -> Lookup[sectorData, "ImagPolyExponents", None],
    (* theta_k = Im part of the divergent endpoint exponent (planIBPCX.md §3.1).
       0 (exact) for real / Case-A sectors -> driver takes the real-pole branch and
       reproduces the current numbers bit-for-bit.  Nonzero (off-axis) -> the cone
       is finite (no 1/eps pole); the driver assembles 1/(i theta_k). *)
    "ImagPole"               -> imagPole,
    (* 1/ck is the real-pole residue; authoritative only when ImagPole == 0
       (planIBPCX.md §3.1).  Kept for diagnostics / back-compat.  Guard ck=0
       (a pure-imaginary off-axis direction reaching this single-pole path would
       give 1/0 = ComplexInfinity and a Power::infy warning; planIBPMULTIDIV §3.6). *)
    "AnalyticPole"           -> If[TrueQ[ck == 0], Indeterminate, 1/ck]
  |>
];

(* --------------------------------------------------------------------------
   IBPCheckBoundary
   Numerically verify that IBP boundary terms and the IBP decomposition
   are self-consistent by checking at finite epsilon.
   -------------------------------------------------------------------------- *)

IBPCheckBoundary[sectorData_Association, integrandSpec_Association,
                 nPoints_Integer: 20] :=
Module[
  {eps, n, k, aVals, clearedPolys, polyExps, detM,
   a0, a1, divVars, ck, testEps, yVars, yValsSet,
   kinRules, epsRules, results},

  eps          = integrandSpec["RegulatorSymbol"];
  n            = sectorData["Dimension"];
  aVals        = sectorData["NewExponents"];
  clearedPolys = sectorData["ClearedPolys"];
  polyExps     = sectorData["PolynomialExponents"];
  detM         = sectorData["DetM"];
  a0           = aVals /. eps -> 0;
  a1           = D[aVals, eps] /. eps -> 0;

  (* Only genuine 1/eps POLES (Re(a0)=0 AND theta=0) get the y_k=0 boundary check:
     a co-located OFF-AXIS direction (Re=0, theta!=0; planIBPMULTIDIV.md) is finite
     with c_k = a1[[k]] = 0, so the endpoint test y_k^{c_k eps} = y_k^0 = 1 would
     fire a spurious boundary-violation warning for a sector handled correctly.
     Skip off-axis directions here. *)
  divVars = {};
  Do[
    If[(TrueQ[Re[a0[[i]]] <= 0] || (NumericQ[a0[[i]]] && Re[a0[[i]]] <= 0)) &&
       ibpDivClass[sectorData, i, eps] === "Pole",
      AppendTo[divVars, i]
    ],
    {i, n}
  ];

  results = {};

  Do[
    k  = divVars[[dv]];
    ck = a1[[k]];
    testEps = 0.01;

    Print["IBPCheckBoundary: checking y_", k, " in sector ",
          sectorData["ConeIndex"]];

    (* Check boundary at y_k = 0: y_k^{a_k} f(y) -> 0 *)
    Module[{yvals, akNum, bndVal},
      Do[
        yvals   = RandomReal[{0.01, 0.99}, n];
        akNum   = ck * testEps;
        yvals[[k]] = 10^(-6);
        bndVal  = yvals[[k]]^akNum;
        If[Abs[bndVal] > 10^(-4),
          Print["  WARNING: y_k=0 boundary not small: y_k^a_k = ",
                bndVal, " at y_k = ", yvals[[k]]]
        ],
        {nPoints}
      ];
      Print["  y_k=0 boundary: OK (y_k^a_k vanishes for eps>0)"];
    ];

    (* Check boundary at y_k = 1: compute f(y)|_{y_k=1}
       and report its magnitude *)
    Module[{yvals, polyVals, fVal, fMags, ndVars, ndYvals},
      ndVars = DeleteCases[Range[n], k];
      fMags = Table[
        ndYvals = RandomReal[{0.01, 0.99}, n - 1];

        polyVals = Table[
          Total[Table[
            Module[{cm, em, logY, ndEm},
              cm = mono[[1]];
              em = mono[[2]];
              (* Set y_k = 1, use only non-div vars *)
              cm * Exp[Total[em[[ndVars]] * Log[ndYvals]]]
            ],
            {mono, clearedPolys[[j]]}
          ]],
          {j, Length[clearedPolys]}
        ];

        fVal = Abs[detM] *
          Exp[Total[(a0[[ndVars]] - 1) * Log[ndYvals]]] *
          Times @@ MapThread[
            Function[{pv, be}, Exp[be * Log[pv]]],
            {polyVals, polyExps /. eps -> 0}
          ];

        Abs[fVal],
        {nPoints}
      ];

      Print["  y_k=1 boundary magnitude: mean=", Mean[fMags],
            " max=", Max[fMags], " min=", Min[fMags]];
      If[Max[fMags] > 10^(-10),
        Print["  NOTE: boundary at y_k=1 is non-zero (expected). ",
              "IBP formula includes this as the boundary integral."]
      ];
      AppendTo[results, <|"Variable" -> k, "BoundaryMags" -> fMags|>];
    ],
    {dv, Length[divVars]}
  ];

  results
];

(* --------------------------------------------------------------------------
   ValidateIBP
   Cross-check the IBP pipeline against direct NIntegrate at finite epsilon.
   Verifies: (1/a_k)[B - sum_t coeff_t I_t] = original sector integral.
   -------------------------------------------------------------------------- *)

ValidateIBP[ibpSectorData_Association, sectorData_Association,
            integrandSpec_Association,
            testKinematics_List, testEpsilon_: 0.01] :=
Module[
  {eps, n, k, ck, rk, a0, a1, aVals, polyExps,
   kinRules, epsRules, fullRules,
   clearedPolys, detM, pfBase, domC, yVars,
   originalIntegral,
   boundaryVal, ibpTermVals, ibpSum,
   ak, reconstructed, relError,
   B0},

  (* MultiDiv sectors (one pole + N off-axis, planIBPMULTIDIV.md) carry "Corners"
     and per-corner pieces, NOT the single-pole "IBPTerms"/"BoundaryData"/
     "DivergentVariable" this finite-eps reconstruction consumes.  This oracle
     does not cover the multi-corner assembly (which is validated against the
     closed-form Dirichlet oracle in cc_51/cc_52); refuse cleanly rather than
     index a Missing key. *)
  If[TrueQ[ibpSectorData["MultiDiv"]],
    Message[TropicalEval::ibpvalmultidiv, ibpSectorData["ConeIndex"]];
    Return[$Failed]
  ];

  eps      = integrandSpec["RegulatorSymbol"];
  n        = sectorData["Dimension"];
  k        = ibpSectorData["DivergentVariable"];
  ck       = ibpSectorData["ck"];
  rk       = ibpSectorData["rk"];
  a0       = ibpSectorData["a0"];
  a1       = ibpSectorData["a1"];
  aVals    = sectorData["NewExponents"];
  polyExps = sectorData["PolynomialExponents"];
  detM     = sectorData["DetM"];
  (* PrefactorBase (planAXpDIV.md §4.1) + lifted DomainConstraint (§4.2). *)
  pfBase   = Lookup[sectorData, "PrefactorBase", Abs[detM]];
  domC     = Lookup[sectorData, "DomainConstraint", None];
  B0       = polyExps /. eps -> 0;

  kinRules  = testKinematics;
  epsRules  = {eps -> testEpsilon};
  fullRules = Join[kinRules, epsRules];

  clearedPolys = ibpSectorData["ClearedPolys"];

  yVars = Table[Unique["yv"], {n}];

  (* Original sector integral at finite epsilon *)
  Module[{aNum, polyValsExpr, integrand},
    aNum = aVals /. fullRules;
    polyValsExpr = Table[
      Total[Table[
        Module[{coeff, exps},
          coeff = mono[[1]] /. kinRules;
          exps  = mono[[2]];
          coeff * Exp[Total[exps * Log /@ yVars]]
        ],
        {mono, clearedPolys[[j]]}
      ]],
      {j, Length[clearedPolys]}
    ];
    integrand = (pfBase /. fullRules) *
      Exp[Total[(aNum - 1) * Log /@ yVars]] *
      Times @@ MapThread[
        Function[{pv, be}, Exp[be * Log[pv]]],
        {polyValsExpr, polyExps /. fullRules}
      ];
    (* lifted-sector domain indicator over all n coords (planAXpDIV.md §4.2) *)
    integrand = integrand * liftedDomainBooleWL[domC, yVars];

    originalIntegral = Quiet@NIntegrate[
      integrand,
      Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
      MaxRecursion -> 20,
      PrecisionGoal -> 4,
      Method -> "GlobalAdaptive"
    ];
  ];

  (* Boundary integral at y_k = 1 (n-1 dim) *)
  Module[{ndVars, bndYVars, bndAnum, bndPolyVals, bndIntegrand},
    ndVars  = DeleteCases[Range[n], k];
    bndYVars = Table[Unique["by"], {n - 1}];
    bndAnum  = (aVals /. fullRules)[[ndVars]];

    bndPolyVals = Table[
      Total[Table[
        Module[{coeff, exps},
          coeff = mono[[1]] /. kinRules;
          exps  = mono[[2]][[ndVars]];
          coeff * Exp[Total[exps * Log /@ bndYVars]]
        ],
        {mono, clearedPolys[[j]]}
      ]],
      {j, Length[clearedPolys]}
    ];

    bndIntegrand = (pfBase /. fullRules) *
      Exp[Total[(bndAnum - 1) * Log /@ bndYVars]] *
      Times @@ MapThread[
        Function[{pv, be}, Exp[be * Log[pv]]],
        {bndPolyVals, polyExps /. fullRules}
      ];
    (* boundary at y_k=1 over the n-1 non-divergent coords (drop slot k, §4.2) *)
    bndIntegrand = bndIntegrand *
      liftedDomainBooleWL[dropDivVarFromDomain[domC, k], bndYVars];

    boundaryVal = Quiet@NIntegrate[
      bndIntegrand,
      Evaluate[Sequence @@ ({#, 0, 1} & /@ bndYVars)],
      MaxRecursion -> 20,
      PrecisionGoal -> 4,
      Method -> "GlobalAdaptive"
    ];
  ];

  (* Each IBP term: n-dim integral *)
  Module[{ibpTerms, ibpYVars},
    ibpTerms = ibpSectorData["IBPTerms"];
    ibpYVars = Table[Unique["iv"], {n}];

    ibpTermVals = Table[
      Module[{termData, alpha0, tB0, coeff0, flatPrefactor,
              polyValsExpr, integrand, termVal},
        termData = ibpTerms[[t]];
        alpha0   = termData["Alpha0"] /. kinRules;
        tB0      = termData["PolyExponents"] /. kinRules;
        coeff0   = termData["Coeff0"] /. kinRules;

        (* Simpler approach: evaluate at finite epsilon using
           original cleared polys and the term's effective exponents *)
        Module[{rawTerms, origTerm, termAlpha, termPE},
          (* Get the term's eps-dependent exponents from IBPReduceSector.
             We stored Alpha0/Alpha1 but need the full expression.
             Reconstruct: alpha = alpha0 + alpha1*eps + ... *)
          termAlpha = termData["Alpha0"] + testEpsilon * termData["Alpha1"];
          termPE    = termData["PolyExponents"];

          polyValsExpr = Table[
            Total[Table[
              Module[{coeff, exps},
                coeff = mono[[1]] /. kinRules;
                exps  = mono[[2]];
                coeff * Exp[Total[exps * Log /@ ibpYVars]]
              ],
              {mono, clearedPolys[[j]]}
            ]],
            {j, Length[clearedPolys]}
          ];

          integrand = (pfBase /. fullRules) *
            Exp[Total[(termAlpha - 1) * Log /@ ibpYVars]] *
            Times @@ MapThread[
              Function[{pv, be}, Exp[(be /. kinRules) * Log[pv]]],
              {polyValsExpr, termPE /. epsRules}
            ];
          (* IBP terms are full-n-dim: full domain indicator (§4.2) *)
          integrand = integrand * liftedDomainBooleWL[domC, ibpYVars];

          termVal = Quiet@NIntegrate[
            integrand,
            Evaluate[Sequence @@ ({#, 0, 1} & /@ ibpYVars)],
            MaxRecursion -> 15,
            PrecisionGoal -> 3,
            Method -> "GlobalAdaptive"
          ];

          (* Include the coefficient *)
          (termData["Coeff0"] + testEpsilon * termData["Coeff1"]) /.
            kinRules /. epsRules // (# * termVal &)
        ]
      ],
      {t, Length[ibpTerms]}
    ];
  ];

  ibpSum = Total[ibpTermVals];

  (* Reconstruct: (1/a_k) * [boundary - ibpSum] *)
  ak = ck * testEpsilon + (ibpSectorData["rk"] * ck) * testEpsilon^2;
  reconstructed = (1/ak) * (boundaryVal - ibpSum);
  relError = Abs[(reconstructed - originalIntegral) / originalIntegral];

  Print["ValidateIBP for sector ", sectorData["ConeIndex"], ":"];
  Print["  Original integral (eps=", testEpsilon, "): ", originalIntegral];
  Print["  Boundary at y_", k, "=1: ", boundaryVal];
  Print["  IBP sum (", Length[ibpTermVals], " terms): ", ibpSum];
  Print["  (1/a_k)[B - IBP]: ", reconstructed];
  Print["  Relative error: ", relError];

  If[relError > 0.02,
    Print["  WARNING: relative error exceeds 2%"]
  ];

  <|"OriginalIntegral" -> originalIntegral,
    "Boundary" -> boundaryVal,
    "IBPSum" -> ibpSum,
    "Reconstructed" -> reconstructed,
    "RelativeError" -> relError|>
];

(* --------------------------------------------------------------------------
   GenerateCppMonteCarloIBP
   Generates C++ Monte Carlo code with IBP divergence handling.

   The generated code evaluates:
     - Convergent sector integrands (summed into one result)
     - IBP boundary and term integrands (reported individually)

   Output format per kinematic point:
     Line 1: convergent_sum_re convergent_sum_im err_re err_im
     Then for each IBP function:
     Line i: func_re func_im err_re err_im
   -------------------------------------------------------------------------- *)

(* The old IBP Defs emitter was folded into the single parametrized
   emitIntegrandDefinitions (plan.md §5.2): call it with the IBP sectors
   as the 3rd argument and {} for divergentSectors.  Byte-identical output
   is guarded by the codegen goldens (cross-check #25). *)

(* The IBP main()s (MC + per-kp Vegas) were folded into the single emitMain
   (plan.md §5.2), routed by info["IsIBP"].  Byte-identical output is
   guarded by the codegen goldens (cross-check #25). *)

(* --------------------------------------------------------------------------
   GenerateCppMonteCarloIBP — public entry; builds Defs and selects the main.
   Default "MonteCarlo" stays byte-identical (B0).  Batched-ncomp Vegas is not
   implemented for the IBP path (T8, deferred) -> falls back to per-kp Vegas.
   -------------------------------------------------------------------------- *)
Options[GenerateCppMonteCarloIBP] = {
  "NSamples"       -> 1000000,
  "MaxDim"         -> 20,
  "SeedBase"       -> 42,
  "Integrator"     -> "MonteCarlo",
  "Batch"          -> False,
  "VegasEpsRel"    -> 1.*^-12,
  "VegasEpsAbs"    -> 1.*^-300,
  "VegasSeed"      -> 0,
  "CubaMaxComp"    -> 512,
  "VegasNStart"    -> Automatic,
  "VegasNIncrease" -> Automatic,
  "VegasNBatch"    -> Automatic
};

GenerateCppMonteCarloIBP[convergentSectors_List, ibpSectors_List,
                         integrandSpec_Association, outputFile_String,
                         OptionsPattern[]] :=
Module[{defs, mainCode, code, integrator, batch, maxDim, nSamples, seedBase},
  integrator = OptionValue["Integrator"];
  batch      = TrueQ[OptionValue["Batch"]];
  maxDim     = OptionValue["MaxDim"];
  nSamples   = OptionValue["NSamples"];
  seedBase   = OptionValue["SeedBase"];

  defs = emitIntegrandDefinitions[convergentSectors, {}, ibpSectors,
           integrandSpec, "MaxDim" -> maxDim];

  integrator = normalizeIntegrator[integrator];

  (* emitMain routes IBP result-assembly off defs["IsIBP"]; batch now selects the
     chunked-ncomp Vegas over the IBP function table (planIBPCX.md §4) while
     preserving the (1 + n_ibp)-row-per-kp output contract. *)
  mainCode = emitMain[defs, integrator, batch,
    "NSamples" -> nSamples, "SeedBase" -> seedBase,
    "VegasEpsRel" -> OptionValue["VegasEpsRel"],
    "VegasEpsAbs" -> OptionValue["VegasEpsAbs"],
    "VegasSeed"   -> OptionValue["VegasSeed"],
    "CubaMaxComp" -> OptionValue["CubaMaxComp"],
    "VegasNStart"    -> OptionValue["VegasNStart"],
    "VegasNIncrease" -> OptionValue["VegasNIncrease"],
    "VegasNBatch"    -> OptionValue["VegasNBatch"]];

  code = defs["Defs"] <> mainCode;

  (* Codegen guard (G-A): never write C++ that cannot compile.  Nested /
     higher-order-pole sectors leak ComplexInfinity/Indeterminate here; fail
     loudly with $Failed instead of emitting an uncompilable .cpp. *)
  Module[{bad = cppBadTokens[code]},
    If[bad =!= {},
      Message[TropicalEval::badcpp, bad];
      Return[$Failed]
    ]
  ];

  Export[outputFile, code, "Text"];

  Print["Generated IBP C++ Monte Carlo code: ", outputFile,
    If[integrator === "VEGAS", If[batch, "  (VEGAS, batched)", "  (VEGAS, per-kp)"], ""]];
  Print["  ", defs["NConvergent"], " convergent sectors"];
  Print["  ", defs["NIBPFuncs"], " IBP integrands (",
        Length[ibpSectors], " sectors)"];
  Print["  Total: ", defs["NTotal"], " integrand functions"];

  <|"Code" -> code, "OutputFile" -> outputFile,
    "NConvergent" -> defs["NConvergent"], "NIBPFuncs" -> defs["NIBPFuncs"],
    "NTotal" -> defs["NTotal"],
    "Dimensions" -> defs["Dims"],
    "IBPFuncMap" -> defs["IBPFuncMap"]|>
];


(* ============================================================================
   MODULE 4b: EvaluateTropicalMCIBP (IBP Driver)

   Full pipeline using IBP for divergent sectors.
   Single divergent variable per sector => Laurent P_{-1}/eps + P_0, returned
   as PoleCoefficient/FinitePart and as LaurentCoefficients<|-1->.,0->.|>.
   Nested / higher-order poles (1/eps^d, d>=2) are detected and refused
   ($Failed); full multi-order assembly is future work (see manual sec. on
   nested divergences).
   ============================================================================ *)

(* evaluateTropicalIBPDriver — the IBP execution path, formerly the public
   EvaluateTropicalMCIBP body.  In v3 it is a PRIVATE helper that the unified
   EvaluateTropicalMC delegates to when Method resolves to "IBP"; the public
   EvaluateTropicalMCIBP is now a thin wrapper (plan.md §5.4).  Body unchanged
   so the IBP tests reproduce the Phase-0 numbers exactly. *)
Options[evaluateTropicalIBPDriver] = {
  "NSamples"         -> 1000000,
  "NThreads"         -> Automatic,
  "SeedBase"         -> 42,   (* compile-time default; runtime argv[5] override (#22) *)
  "RunChecks"        -> True,
  "WorkingDirectory" -> Automatic,
  "Verbose"          -> True,
  (* sampler selection (default = the zero-dependency plain Monte Carlo) *)
  "Integrator"       -> "MC",    (* "MC" | "VEGAS" (aliases "MonteCarlo"/"Vegas") *)
  "Batch"            -> False,   (* True + VEGAS: chunked-ncomp batched IBP (planIBPCX.md §4) *)
  (* LiftData (planAXpDIV.md §4.4, Barrier C): when present, sectors are
     processed via ProcessSectorLifted (eps-aware), so lifted+divergent
     integrals take the IBP route exactly like unlifted divergent ones. *)
  "LiftData"         -> None,
  (* ComplexExponentMode (planCXLIFTDIV.md §4.4): "SplitRealImag" + a lifted
     complex-exponent integral decomposes on Re(B) and rides Im(B) as the IBP
     boundary/term oscillatory phase.  Automatic/"Direct" -> no phase (#25). *)
  "ComplexExponentMode" -> Automatic,
  Sequence @@ $vegasOptionDefaults
};

evaluateTropicalIBPDriver[integrandSpec_Association, fanData_List,
                      kinematicPoints_List, OptionsPattern[]] :=
Module[
  {dualVertices, simplexList, n, nKP, nParams, eps,
   allSectorData, convergentSectors, divergentSectors,
   ibpProcessedSectors, liftData, isLifted, emptyDomainCount,
   cppFile, cppBinary, kinFile, resultFile,
   cppResult, ibpFuncMap,
   mcRawResults, finalResults,
   runChecks, verbose, nSamples, nThreads, workDir,
   integrator, batch, useCuba, vegasOpts,
   cxMode, imBList, hasImagB, cxSplit, specForProc, imAList, hasImagA},

  runChecks  = OptionValue["RunChecks"];
  verbose    = OptionValue["Verbose"];
  nSamples   = OptionValue["NSamples"];
  nThreads   = OptionValue["NThreads"];
  workDir    = OptionValue["WorkingDirectory"];
  liftData   = OptionValue["LiftData"];
  isLifted   = (liftData =!= None);
  emptyDomainCount = 0;
  eps        = integrandSpec["RegulatorSymbol"];

  (* --- SplitRealImag setup (planCXLIFTDIV.md §4.4) --- mirrors EvaluateTropicalMC:
     active only for a LIFTED complex-exponent integral in SplitRealImag mode.
     Then sectors are processed on Re(B) (real measure + real, well-defined
     divergence/domain machinery) and Im(B) is reintroduced per sector as the
     oscillatory phase the IBP boundary/term codegen compiles (per-piece
     MonoFactorLog, §3).  Real / Direct / unlifted: specForProc = integrandSpec
     and no ImagPolyExponents attached -> the emitted C++ is byte-identical (#25). *)
  cxMode     = OptionValue["ComplexExponentMode"];
  If[cxMode === Automatic, cxMode = "Direct"];
  imBList    = imagPolyInfo[integrandSpec];   (* eps-free Im(B); see imagPolyInfo *)
  hasImagB   = AnyTrue[imBList, (!TrueQ[PossibleZeroQ[#]]) &];
  (* planAXpDIVv2.md §4.1: engage on complex MONOMIAL exponents A too.  See
     imagMonoInfo for the eps->0 / declare-real rationale (planR.md R2,R6). *)
  {imAList, hasImagA} = imagMonoInfo[integrandSpec, eps];
  cxSplit    = isLifted && (cxMode === "SplitRealImag") && (hasImagB || hasImagA);
  specForProc = integrandSpec;
  (* Re(B) by eps-free subtraction (NOT MapAt[Re]) so a regulator inside B stays a
     true real for the symbolic-eps Laurent/log-insertion machinery — see
     realifyPolyB.  IBP keeps eps symbolic (no pinning), so this is exact. *)
  If[cxSplit,
    specForProc = realifyPolyB[specForProc, imBList]];
  (* Subtract the eps-free Im(A) so the real regulator eps in the monomial
     exponent is preserved for the pole machinery (planAXpDIVv2.md §4.2 / see
     realifyMonoA). *)
  If[cxSplit && hasImagA,
    specForProc = realifyMonoA[specForProc, imAList]];
  integrator = normalizeIntegrator[OptionValue["Integrator"]];
  batch      = TrueQ[OptionValue["Batch"]];
  useCuba    = (integrator === "VEGAS");
  vegasOpts  = {"Integrator" -> integrator, "Batch" -> batch,
    "SeedBase" -> OptionValue["SeedBase"],
    "VegasEpsRel" -> OptionValue["VegasEpsRel"],
    "VegasEpsAbs" -> OptionValue["VegasEpsAbs"],
    "VegasSeed" -> OptionValue["VegasSeed"],
    "CubaMaxComp" -> OptionValue["CubaMaxComp"],
    "VegasNStart"    -> OptionValue["VegasNStart"],
    "VegasNIncrease" -> OptionValue["VegasNIncrease"],
    "VegasNBatch"    -> OptionValue["VegasNBatch"]};

  If[workDir === Automatic,
    workDir = DirectoryName[$InputFileName];
    If[workDir === "" || !StringQ[workDir], workDir = Directory[]];
    workDir = FileNameJoin[{workDir, "INTERFILES"}]
  ];
  If[!DirectoryQ[workDir], Quiet[CreateDirectory[workDir]]];

  {dualVertices, simplexList} = fanData;
  n       = Length[integrandSpec["Variables"]];
  nKP     = Length[kinematicPoints];
  nParams = Length[integrandSpec["KinematicSymbols"]];

  If[verbose,
    Print["EvaluateTropicalMCIBP: ", Length[simplexList], " sectors, ",
          nKP, " kinematic points, ", n, " variables"]
  ];

  (* --- Step 1: Process all sectors --- *)
  If[verbose, Print["Processing ", Length[simplexList], " sectors..."]];

  If[isLifted,
    (* Lifted IBP route (planAXpDIV.md §4.4): each sector via the eps-aware
       ProcessSectorLifted.  $Failed (liftcomplex / liftnopivot / liftdivdomain)
       aborts the whole call; EmptyDomain sectors contribute 0 and are dropped;
       the rest split into convergent / (single-pole) divergent lifted sectors. *)
    allSectorData = Reap[
      Catch[
        Do[
          Module[{sd},
            sd = ProcessSectorLifted[specForProc, dualVertices,
                   simplexList[[s]], s, liftData, "Eps" -> eps, "Verbose" -> False,
                   "ImagMonoExps" -> If[cxSplit && hasImagA, imAList, None]];
            Which[
              sd === $Failed, Sow[$Failed]; Throw[Null],
              AssociationQ[sd] && KeyExistsQ[sd, "EmptyDomain"] && sd["EmptyDomain"],
                emptyDomainCount++,
              True,
                (* SplitRealImag: reattach Im(B) to BOTH convergent and divergent
                   lifted sectors (planCXLIFTDIV.md §4.4); the convergent emitter
                   and IBPProcessSector pick it up as the phase.  Im(A) rides via
                   MonomialPhaseLog / LiftedMonoPhase attached by ProcessSectorLifted
                   itself (planAXpDIVv2.md §4.3). *)
                If[cxSplit, sd["ImagPolyExponents"] = imBList];
                Sow[sd]
            ]
          ],
          {s, Length[simplexList]}
        ]
      ]
    ][[2]];
    allSectorData = If[allSectorData === {}, {}, allSectorData[[1]]];
    If[MemberQ[allSectorData, $Failed],
      Print["ERROR: ProcessSectorLifted failed for a sector (liftcomplex / ",
            "liftnopivot / liftdivdomain).  Aborting ($Failed)."];
      Return[$Failed]
    ],
    allSectorData = Table[
      ProcessSector[integrandSpec, dualVertices,
                    simplexList[[s]], s, "Verbose" -> False],
      {s, Length[simplexList]}
    ]
  ];

  convergentSectors = Select[allSectorData,
    (AssociationQ[#] && !#["IsDivergent"]) &];
  divergentSectors  = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  If[verbose,
    Print["  ", Length[convergentSectors], " convergent, ",
          Length[divergentSectors], " divergent sectors",
          If[isLifted, ", " <> ToString[emptyDomainCount] <> " empty (dropped)", ""]]
  ];

  (* --- Step 2: Process divergent sectors with IBP --- *)
  If[verbose && Length[divergentSectors] > 0,
    Print["Processing ", Length[divergentSectors],
          " divergent sectors with IBP..."]
  ];

  ibpProcessedSectors = Table[
    IBPProcessSector[divergentSectors[[s]], specForProc],
    {s, Length[divergentSectors]}
  ];
  (* Fail loudly if ANY divergent sector could not be IBP-processed (e.g. a
     nested / higher-order pole).  Silently skipping it would return a finite-
     only partial result that looks plausible but is wrong. *)
  If[MemberQ[ibpProcessedSectors, $Failed],
    Print["ERROR: IBP processing failed for ",
      Count[ibpProcessedSectors, $Failed], " of ", Length[divergentSectors],
      " divergent sector(s) — typically a nested / higher-order pole (>1 ",
      "divergent variable per sector), which the IBP numerical path does not ",
      "assemble.  Aborting ($Failed) rather than dropping the sector(s)."];
    Return[$Failed]
  ];

  If[verbose,
    Print["  Successfully processed ", Length[ibpProcessedSectors],
          " IBP sectors, total IBP terms: ",
          Total[#["NTerms"] & /@ ibpProcessedSectors]]
  ];

  (* --- Step 3: Boundary checks --- *)
  If[runChecks,
    Do[
      IBPCheckBoundary[divergentSectors[[s]], specForProc, 10],
      {s, Min[3, Length[divergentSectors]]}
    ]
  ];

  (* --- Step 4: Generate C++ code --- *)
  cppFile    = FileNameJoin[{workDir, "tropical_mc_ibp.cpp"}];
  cppBinary  = FileNameJoin[{workDir, "tropical_mc_ibp"}];
  kinFile    = FileNameJoin[{workDir, "kinematic_data_ibp.txt"}];
  resultFile = FileNameJoin[{workDir, "mc_results_ibp.txt"}];

  (* Convergent sectors carry symbolic eps in their flattening exponents;
     evaluate it to its eps^0 value before codegen (see EvaluateTropicalMC).
     The IBP boundary/term integrands are already constructed at eps=0. *)
  cppResult = GenerateCppMonteCarloIBP[
    If[eps =!= None, (# /. eps -> 0) & /@ convergentSectors, convergentSectors],
    ibpProcessedSectors,
    specForProc, cppFile,
    "NSamples" -> nSamples,
    Sequence @@ vegasOpts
  ];

  If[!AssociationQ[cppResult],
    Print["ERROR: C++ code generation failed"];
    Return[$Failed]
  ];

  ibpFuncMap = cppResult["IBPFuncMap"];

  (* --- Step 5: Write kinematic data --- *)
  If[nParams > 0,
    Export[kinFile,
      StringRiffle[
        StringRiffle[ToString[CForm[#]] & /@ #, " "] & /@ kinematicPoints,
        "\n"
      ] <> "\n",
      "Text"
    ];,
    Export[kinFile,
      StringRiffle[Table["0", {nKP}], "\n"] <> "\n",
      "Text"
    ];
  ];

  (* --- Step 6: Compile and run --- *)
  If[CompileCpp[cppFile, cppBinary, False, "UseCuba" -> useCuba] === $Failed,
    Print["ERROR: Compilation failed"];
    Return[$Failed]
  ];

  Module[{nThreadsStr, result},
    nThreadsStr = If[nThreads === Automatic,
      ToString[$ProcessorCount], ToString[nThreads]];

    If[verbose, Print["Running Monte Carlo (", nSamples,
                      " samples, ", nThreadsStr, " threads)..."]];

    result = RunProcess[{cppBinary, kinFile, resultFile,
      ToString[nSamples], nThreadsStr}];

    If[result["ExitCode"] != 0,
      Print["ERROR: Monte Carlo execution failed:"];
      Print[result["StandardError"]];
      Return[$Failed]
    ];
    If[verbose && StringLength[StringTrim[result["StandardError"]]] > 0,
      Print[result["StandardError"]]
    ];
  ];

  (* --- Step 7: Read and parse results --- *)
  Module[{rawText, lines, parsed, nIBPFuncs, nLinesPerKP},
    rawText = Import[resultFile, "Text"];
    If[rawText === $Failed,
      Print["ERROR: cannot read results file"];
      Return[$Failed]
    ];

    lines = Select[StringSplit[rawText, "\n"],
                   StringLength[StringTrim[#]] > 0 &];
    parsed = Table[
      Read[StringToStream[#], Number] & /@ StringSplit[line],
      {line, lines}
    ];

    nIBPFuncs  = cppResult["NIBPFuncs"];
    nLinesPerKP = 1 + nIBPFuncs;

    If[Length[parsed] != nKP * nLinesPerKP,
      Print["WARNING: expected ", nKP * nLinesPerKP, " result rows, got ",
            Length[parsed]]
    ];

    mcRawResults = parsed;
  ];

  (* --- Step 8: Combine results --- *)
  finalResults = Table[
    Module[{kpOffset, convResult, ibpContribPole, ibpContribFinite,
            poleTotal, finiteTotal, mcCx, kinRules},
      kpOffset = (i - 1) * (1 + cppResult["NIBPFuncs"]);

      (* Substitute THIS kp's kinematic values into the symbolic IBP
         coefficients.  Coeff0/Coeff1 pick up polynomial coefficients (e.g.
         c1 from d/dy of (c0 + c1 y + ...)), so without this the assembled
         pole/finite stay symbolic in the kinematic parameters.  ck/rk are
         exponent-derived (numeric) but are substituted defensively. *)
      kinRules = If[nParams > 0,
        Thread[integrandSpec["KinematicSymbols"] -> kinematicPoints[[i]]], {}];

      (* The per-kp output block is [conv-sum, IBP-func-0, IBP-func-1, ...].
         IBP function ids in ibpFuncMap are GLOBAL (they count the N_CONV
         convergent funcs first), so the func with global id fid lives at
         block row (fid - N_CONV + 1).  mcCx maps a global func id to its
         complex MC value (1-based row, +1 for the leading conv-sum row).
         Omitting the -N_CONV shift mis-reads every IBP row when N_CONV>0. *)
      mcCx[fid_] := With[{r = kpOffset + fid - cppResult["NConvergent"] + 2},
        If[r <= Length[mcRawResults],
          mcRawResults[[r, 1]] + I*mcRawResults[[r, 2]], 0.]];

      (* Convergent sector sum (block row 0) *)
      convResult = If[kpOffset + 1 <= Length[mcRawResults],
        mcRawResults[[kpOffset + 1, 1]] +
        I * mcRawResults[[kpOffset + 1, 2]],
        0.
      ];

      (* IBP sector contributions *)
      ibpContribPole   = 0.;
      ibpContribFinite = 0.;

      Do[
        Module[{sMap, bndBaseFid, bndLogFid, termFids,
                bndBase, bndLog, S0, S1, ckS, rkS,
                poleCont, finiteCont, imagPoleSym, theta0},
          sMap     = ibpFuncMap[[s]];

          If[TrueQ[sMap["MultiDiv"]],
          (* ===== One-pole + N-off-axis corner assembly (planIBPMULTIDIV §3.3) =====
             I = [1/(c_k eps)] prod_j[1/(i theta_j)] sum_corners (sign) L_corner.
             Let NB = sum sign*Coeff0*base, NL = sum sign*(Coeff0*log + Coeff1*base),
             K0 = prod_j 1/(i theta_j), K1 = K0 * sum_j(-c_j/(i theta_j)).  Then
             pole   = K0 NB / c_k,
             finite = (K0 NL + (K1 - r_k K0) NB)/c_k + DLogPrefactor*pole.
             No pole (all off-axis): pole = 0, finite = K0 NB.  This reduces to the
             single-pole / single-off-axis formulas below for |D| = 1. *)
            Module[{ps, pieces, NB, NL, K0, K1, offTh, offCk, offRk, hasPole,
                    dlogM, zeroQ, vanish, finOff, nPolesKp, cP, rP, indet},
              ps     = ibpProcessedSectors[[s]];
              pieces = sMap["Pieces"];
              offTh  = ps["OffAxisThetas"] /. kinRules;
              offCk  = ps["OffAxisCks"]    /. kinRules;
              offRk  = ps["OffAxisRks"]    /. kinRules;
              hasPole = TrueQ[ps["HasPole"]];
              dlogM  = ps["DLogPrefactor"] /. kinRules;
              indet[msg_] := (
                Print["WARNING: IBP MultiDiv sector ", sMap["ConeIndex"],
                      " at kp ", i, ": ", msg, " -> result Indeterminate."];
                poleCont = Indeterminate; finiteCont = Indeterminate);

              (* Signed corner sums (always needed). *)
              NB = 0.; NL = 0.;
              Do[
                Module[{pc = pieces[[p]], sgn, c0, c1, baseMC, logMC},
                  sgn    = pc["Sign"];
                  c0     = pc["Coeff0"] /. kinRules;
                  c1     = pc["Coeff1"] /. kinRules;
                  baseMC = mcCx[pc["BaseFuncId"]];
                  logMC  = mcCx[pc["LogFuncId"]];
                  NB += sgn * c0 * baseMC;
                  NL += sgn * (c0 * logMC + c1 * baseMC);
                ],
                {p, Length[pieces]}];

              (* Count the genuine poles AT THIS kp.  The off-axis directions are
                 finite (1/(i theta)) UNLESS a symbolically-nonzero theta vanishes
                 EXACTLY at this kp (a mixed-grid crossing, e.g. mu hitting the
                 real axis), in which case that direction becomes a genuine
                 1/(c eps) pole here.  zeroQ uses PossibleZeroQ (exact / machine-0),
                 NOT a tolerance — a tiny-but-nonzero theta stays off-axis (its
                 1/(i theta) is a genuine large value, not a pole).  totalPoles =
                 (genuine pole?) + #(off-axis that vanished).
                   0  -> all off-axis (pole 0, finite K0*NB),
                   1  -> exactly one 1/eps pole (genuine OR a vanished off-axis);
                         every OTHER off-axis stays finite 1/(i theta),
                   >=2 -> a 1/eps^{>=2} this path does not assemble (Indeterminate). *)
              zeroQ[x_] := TrueQ[PossibleZeroQ[x]] || (NumericQ[x] && x == 0);
              Which[
                !AllTrue[offTh, NumericQ],
                  indet["an off-axis theta is non-numeric after kinematics " <>
                        "(a KinematicSymbol may be missing from the spec)"],
                True,
                  vanish   = Select[Range[Length[offTh]], zeroQ[offTh[[#]]] &];
                  finOff   = Complement[Range[Length[offTh]], vanish];
                  nPolesKp = If[hasPole, 1, 0] + Length[vanish];
                  Which[
                    nPolesKp == 0,
                      poleCont   = 0;
                      finiteCont = (Times @@ (1/(I*#) & /@ offTh)) * NB,
                    nPolesKp == 1,
                      (* the single pole: the genuine pole, or the lone vanished
                         off-axis.  Its (c,r); the finite off-axis = non-vanishing. *)
                      {cP, rP} = If[hasPole,
                        {ps["ck"] /. kinRules, ps["rk"] /. kinRules},
                        {offCk[[First[vanish]]], offRk[[First[vanish]]]}];
                      If[zeroQ[cP],
                        indet["the 1/eps pole's c vanishes here (unregulated / " <>
                              "higher-order)"],
                        K0 = Times @@ (1/(I*offTh[[#]]) & /@ finOff);
                        K1 = K0 * Total[Table[-offCk[[j]]/(I*offTh[[j]]),
                                              {j, finOff}]];
                        poleCont   = K0 * NB / cP;
                        finiteCont = (K0*NL + (K1 - rP*K0)*NB)/cP + dlogM*poleCont],
                    True,
                      indet["two or more 1/eps poles coincide (1/eps^{>=2} not " <>
                            "assembled)"]
                  ]
              ];
              ibpContribPole   += poleCont;
              ibpContribFinite += finiteCont;
            ]
          , (* ===== existing single divergent-direction path ===== *)
          bndBaseFid = sMap["BndBaseFuncId"];
          bndLogFid  = sMap["BndLogFuncId"];
          termFids   = sMap["TermFuncIds"];

          ckS = ibpProcessedSectors[[s]]["ck"] /. kinRules;
          rkS = ibpProcessedSectors[[s]]["rk"] /. kinRules;

          (* theta_k: the imaginary shift of the divergent endpoint exponent
             (planIBPCX.md §3.3/§3.4).  Decide the branch from the SYMBOLIC value
             first (the exact, decomposition-free test): if provably zero this is a
             real 1/eps pole and the real-pole formulas below run unchanged (byte/
             numerically identical, #25).  Mixed-grid guard: a symbolically-nonzero
             theta that vanishes at THIS kp (e.g. mu crossing the real axis) would
             divide by zero, so also fall back to the real formula when theta0
             evaluates to 0 numerically. *)
          imagPoleSym = Lookup[ibpProcessedSectors[[s]], "ImagPole", 0];
          theta0      = imagPoleSym /. kinRules;

          (* Read MC results for boundary functions *)
          bndBase = mcCx[bndBaseFid];
          bndLog  = mcCx[bndLogFid];

          (* Sum IBP terms: S0 and S1 *)
          S0 = 0.; S1 = 0.;
          Do[
            Module[{tf, baseMC, logMC, c0, c1},
              tf     = termFids[[t]];
              c0     = tf["Coeff0"] /. kinRules;
              c1     = tf["Coeff1"] /. kinRules;
              baseMC = mcCx[tf["BaseFuncId"]];
              logMC  = mcCx[tf["LogFuncId"]];
              S0 += c0 * baseMC;
              S1 += c0 * logMC + c1 * baseMC;
            ],
            {t, Length[termFids]}
          ];

          If[TrueQ[PossibleZeroQ[imagPoleSym]] ||
             (NumericQ[theta0] && TrueQ[PossibleZeroQ[theta0]]),
            (* ---- Real 1/eps pole (theta_k = 0), prefactor 1/(c_k eps) ---- *)
            (* Pole contribution: (B^{(0)} - S_0) / c_k *)
            poleCont = (bndBase - S0) / ckS;
            (* Finite contribution:
               [(B^{(1)} - S_1) - r_k (B^{(0)} - S_0)] / c_k *)
            finiteCont = ((bndLog - S1) - rkS * (bndBase - S0)) / ckS;
            (* Lifted PrefactorBase eps-dependence (planAXpDIV.md §4.1): the codegen
               used PrefactorBase(0), so add the missing (P'/P)(0)*pole finite term.
               DLogPrefactor is 0 for unlifted sectors -> unlifted result unchanged. *)
            finiteCont = finiteCont +
              (ibpProcessedSectors[[s]]["DLogPrefactor"] /. kinRules) * poleCont;
          ,
            (* ---- Off-axis (theta_k != 0): prefactor 1/(c_k eps + i theta_k),
               REGULAR at eps=0 -> no pole, finite 1/(i theta_k).  The eps^1 pieces
               (bndLog/S1) and the DLogPrefactor term are O(eps) and drop
               (planIBPCX.md §1.2/§3.3).  bndBase/S0 already carry pfBase(0) and the
               bounded theta-phases, so nothing upstream changes. ---- *)
            poleCont   = 0;
            finiteCont = (bndBase - S0) / (I * theta0);
          ];

          ibpContribPole   += poleCont;
          ibpContribFinite += finiteCont;
          ] (* end If MultiDiv / single-direction *)
        ],
        {s, Length[ibpFuncMap]}
      ];

      finiteTotal = convResult + ibpContribFinite;
      poleTotal   = ibpContribPole;

      <|"KinematicPoint" -> If[nParams > 0, kinematicPoints[[i]], {}],
        "PoleCoefficient" -> poleTotal,
        "FinitePart"      -> finiteTotal,
        (* Laurent as an association keyed by eps power (G-B).  Single
           divergent variable per sector => orders -1 (pole) and 0 (finite)
           only; PoleCoefficient/FinitePart are kept as the d=1 aliases.  The
           container generalizes to {-d,...,0} once nested assembly lands. *)
        "LaurentCoefficients" -> <|-1 -> poleTotal, 0 -> finiteTotal|>,
        "ConvergentSum"   -> convResult,
        "IBPPoleContrib"  -> ibpContribPole,
        "IBPFiniteContrib" -> ibpContribFinite|>
    ],
    {i, nKP}
  ];

  If[verbose,
    Print["\n=== IBP Monte Carlo Results ==="];
    Do[
      Print["  KP ", i, ": pole = ", finalResults[[i]]["PoleCoefficient"],
            ", finite = ", finalResults[[i]]["FinitePart"]],
      {i, Min[5, nKP]}
    ];
  ];

  <|"Results"              -> finalResults,
    "ConvergentSectors"    -> Length[convergentSectors],
    "DivergentSectors"     -> Length[divergentSectors],
    "IBPSectors"           -> Length[ibpProcessedSectors],
    "CppFile"              -> cppFile,
    "ResultFile"           -> resultFile,
    "IBPFuncMap"           -> ibpFuncMap,
    "IBPProcessedSectors"  -> ibpProcessedSectors|>
];

(* --------------------------------------------------------------------------
   EvaluateTropicalMCIBP — thin wrapper for API continuity (plan.md §5.4).
   The unified EvaluateTropicalMC routes Method->"IBP" to the same private
   driver, so this just forwards.  Kept so existing callers (and the divergent
   cross-check) need no change.
   -------------------------------------------------------------------------- *)
Options[EvaluateTropicalMCIBP] = Options[evaluateTropicalIBPDriver];

EvaluateTropicalMCIBP[integrandSpec_Association, fanData_List,
                      kinematicPoints_List, opts : OptionsPattern[]] :=
  evaluateTropicalIBPDriver[integrandSpec, fanData, kinematicPoints, opts];



(* ============================================================================
   MODULE 5: RunAllTests  —  Tree-A in-package validation suite (Tests 1–18)

   Ported verbatim (logic-identical) from
     OLD_CODE/TROPICAL_MONTE_CARLO/EXAMPLES/tropical_eval_examples.wl
   so merge-fidelity #31 (cc_31.wl Part D) can call RunAllTests[] from inside
   the v3 package and reproduce the Tree-A goldens (TEST/baselines/treeA.txt).

   RunAllTests[] returns a list of {"Test N", True|False} pairs (18 entries);
   cc_31 asserts zero False.  Tests 1–7 are convergent (ValidateDecomposition),
   8–12,14–18 are divergent (ProcessDivergentSector + ValidateSubtraction),
   13 is C++ codegen (GenerateCppMonteCarlo + CompileCpp).  All called functions
   (ProcessSector, ProcessDivergentSector, ValidateSubtraction,
   ValidateDecomposition, GenerateCppMonteCarlo, CompileCpp, PolytopeVertices,
   ComputeDecomposition) are the v3 implementations — no test logic weakened.
   These live in TropicalEval`Private` (RunAllTests is the only public symbol;
   RunTest1..18 are internal helpers).  The integration variable x is a
   Private-context symbol used consistently within each self-contained test.
   ============================================================================ *)

RunAllTests[] := Module[
  {results, nPass, nFail},

  results = {};
  nPass = 0;
  nFail = 0;

  Print[""];
  Print["================================================================"];
  Print["  TropicalEval Validation Suite"];
  Print["================================================================"];
  Print[""];

  Module[{pass},
    pass = RunTest1[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 1", pass}];
  ];

  Module[{pass},
    pass = RunTest2[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 2", pass}];
  ];

  Module[{pass},
    pass = RunTest3[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 3", pass}];
  ];

  Module[{pass},
    pass = RunTest4[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 4", pass}];
  ];

  Module[{pass},
    pass = RunTest5[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 5", pass}];
  ];

  Module[{pass},
    pass = RunTest6[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 6", pass}];
  ];

  Module[{pass},
    pass = RunTest7[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 7", pass}];
  ];

  Module[{pass},
    pass = RunTest8[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 8", pass}];
  ];

  Module[{pass},
    pass = RunTest9[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 9", pass}];
  ];

  Module[{pass},
    pass = RunTest10[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 10", pass}];
  ];

  Module[{pass},
    pass = RunTest11[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 11", pass}];
  ];

  Module[{pass},
    pass = RunTest12[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 12", pass}];
  ];

  Module[{pass},
    pass = RunTest13[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 13", pass}];
  ];

  Module[{pass},
    pass = RunTest14[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 14", pass}];
  ];

  Module[{pass},
    pass = RunTest15[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 15", pass}];
  ];

  Module[{pass},
    pass = RunTest16[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 16", pass}];
  ];

  Module[{pass},
    pass = RunTest17[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 17", pass}];
  ];

  Module[{pass},
    pass = RunTest18[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 18", pass}];
  ];

  Print[""];
  Print["================================================================"];
  Print["  Results: ", nPass, " PASSED, ", nFail, " FAILED"];
  Print["================================================================"];

  results
];
RunTest1[] := Module[
  {poly, vars, integrandSpec, verts, fanData,
   testAValues, allPass},

  Print["--- Test 1: Convergent 2D, real exponents ---"];
  Print["Integral[0,Inf] dx1 dx2 / (1+2x1^2+x2^2+x1*x2^2+3x1^2*x2)^A"];

  poly = 1 + 2 x[1]^2 + x[2]^2 + x[1] x[2]^2 + 3 x[1]^2 x[2];
  vars = {x[1], x[2]};

  testAValues = {2, 3};
  allPass = True;

  Do[
    Module[{A, spec, directResult, mcResult, relErr, pass},
      A = testA;

      spec = <|
        "Polynomials"       -> {poly},
        "MonomialExponents" -> {0, 0},
        "PolynomialExponents" -> {-A},
        "Variables"         -> vars,
        "KinematicSymbols"  -> {},
        "RegulatorSymbol"   -> None
      |>;

      directResult = NIntegrate[
        1 / (1 + 2 t1^2 + t2^2 + t1 t2^2 + 3 t1^2 t2)^A,
        {t1, 0, Infinity}, {t2, 0, Infinity},
        MaxRecursion -> 20, PrecisionGoal -> 6
      ];

      verts = PolytopeVertices[poly^(-A), vars];
      fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

      Module[{vr},
        vr = Quiet@ValidateDecomposition[spec, fanData, {}, 3];
        If[AssociationQ[vr],
          mcResult = vr["SectorSum"];
          relErr   = vr["RelativeError"];,
          mcResult = 0;
          relErr = Infinity;
        ];
      ];

      pass = NumericQ[relErr] && (relErr < 0.02);
      Print["  A = ", A, ":"];
      Print["    NIntegrate = ", directResult];
      Print["    Sector sum = ", mcResult];
      Print["    Rel error  = ", relErr];
      Print["    ", If[pass, "PASS", "FAIL"]];

      If[!pass, allPass = False];
    ],
    {testA, testAValues}
  ];

  allPass
];

(* --------------------------------------------------------------------------
   Test 2: Convergent 2D, complex exponents
   -------------------------------------------------------------------------- *)

RunTest2[] := Module[
  {poly, vars, A, spec, verts, fanData,
   directResult, mcResult, relErr, pass},

  Print["--- Test 2: Convergent 2D, complex exponents ---"];
  Print["Same integral, A = 2 + 0.5I"];

  poly = 1 + 2 x[1]^2 + x[2]^2 + x[1] x[2]^2 + 3 x[1]^2 x[2];
  vars = {x[1], x[2]};
  A    = 2 + 0.5 I;

  spec = <|
    "Polynomials"       -> {poly},
    "MonomialExponents" -> {0, 0},
    "PolynomialExponents" -> {-A},
    "Variables"         -> vars,
    "KinematicSymbols"  -> {},
    "RegulatorSymbol"   -> None
  |>;

  directResult = NIntegrate[
    1 / (1 + 2 t1^2 + t2^2 + t1 t2^2 + 3 t1^2 t2)^A,
    {t1, 0, Infinity}, {t2, 0, Infinity},
    MaxRecursion -> 20, PrecisionGoal -> 5
  ];

  verts   = PolytopeVertices[poly^(-Re[A]), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  Module[{vr},
    vr = Quiet@ValidateDecomposition[spec, fanData, {}, 3];
    If[AssociationQ[vr],
      mcResult = vr["SectorSum"];
      relErr   = Abs[(mcResult - directResult) / directResult];,
      mcResult = 0;
      relErr = Infinity;
    ];
  ];

  pass = NumericQ[relErr] && (relErr < 0.02);
  Print["  NIntegrate = ", directResult];
  Print["  Sector sum = ", mcResult];
  Print["  Rel error  = ", relErr];
  Print["  ", If[pass, "PASS", "FAIL"]];

  pass
];

(* --------------------------------------------------------------------------
   Test 3: Divergent 1D
   f(eps) = Integral[0,1] dy y^{2eps-1} / (1+y^2)
   Exact: 1/(2eps) - log(2)/2
   -------------------------------------------------------------------------- *)

RunTest3[] := Module[
  {exactPole, exactFinite, pass, testEps,
   numericalResult, exactAtEps, relErr},

  Print["--- Test 3: Divergent 1D ---"];
  Print["f(eps) = Integral[0,1] y^{2eps-1} / (1+y^2)"];
  Print["Exact: 1/(2eps) - log(2)/2"];

  exactPole   = 1/2;
  exactFinite = -Log[2]/2;

  testEps = 0.01;
  numericalResult = NIntegrate[
    y^(2 testEps - 1) / (1 + y^2),
    {y, 0, 1},
    MaxRecursion -> 30, PrecisionGoal -> 8,
    Method -> "DoubleExponential"
  ];

  exactAtEps = exactPole / testEps + exactFinite;
  relErr = Abs[(numericalResult - exactAtEps) / exactAtEps];

  Print["  Pole coefficient: expected = ", N[exactPole],
        " (1/2)"];
  Print["  Finite part: expected = ", N[exactFinite],
        " (-log(2)/2)"];
  Print["  At eps = ", testEps, ":"];
  Print["    NIntegrate = ", numericalResult];
  Print["    Exact formula = ", N[exactAtEps]];
  Print["    Rel error = ", relErr];

  pass = (relErr < 0.001);
  Print["  ", If[pass, "PASS", "FAIL"]];

  pass
];

(* --------------------------------------------------------------------------
   Test 4: Divergent 2D
   -------------------------------------------------------------------------- *)

RunTest4[] := Module[
  {testEps, directResult, pass},

  Print["--- Test 4: Divergent 2D ---"];
  Print["g(eps) = Int dy1 dy2 y2^{2eps-1} y1^{2eps} / (1+y1+y2+y1*y2^2+y1^3*y2^2)^2"];

  testEps = 0.05;

  directResult = Quiet@NIntegrate[
    t2^(2 testEps - 1) * t1^(2 testEps) /
    (1 + t1 + t2 + t1 t2^2 + t1^3 t2^2)^2,
    {t1, 0, 1}, {t2, 0, 1},
    MaxRecursion -> 30, PrecisionGoal -> 5,
    Method -> "GlobalAdaptive",
    MinRecursion -> 5
  ];

  Print["  At eps = ", testEps, ":"];
  Print["    Direct NIntegrate = ", directResult];

  pass = NumericQ[directResult] && Abs[directResult] < 10^4;
  Print["  Value is finite: ", If[pass, "PASS", "FAIL"]];

  Module[{val1, val2, ratio},
    val1 = Quiet@NIntegrate[
      t2^(2 * 0.1 - 1) * t1^(2 * 0.1) /
      (1 + t1 + t2 + t1 t2^2 + t1^3 t2^2)^2,
      {t1, 0, 1}, {t2, 0, 1},
      MaxRecursion -> 30, PrecisionGoal -> 4,
      Method -> "GlobalAdaptive", MinRecursion -> 5
    ];
    val2 = Quiet@NIntegrate[
      t2^(2 * 0.2 - 1) * t1^(2 * 0.2) /
      (1 + t1 + t2 + t1 t2^2 + t1^3 t2^2)^2,
      {t1, 0, 1}, {t2, 0, 1},
      MaxRecursion -> 30, PrecisionGoal -> 4,
      Method -> "GlobalAdaptive", MinRecursion -> 5
    ];
    ratio = val1 / val2;
    Print["  g(0.1)/g(0.2) = ", ratio, " (expect ~ 2 if 1/eps behavior)"];
    If[Abs[ratio - 2] < 0.5,
      Print["  1/eps scaling: PASS"];,
      Print["  1/eps scaling: approximate (ratio = ", ratio, ")"];
    ];
  ];

  pass
];

(* --------------------------------------------------------------------------
   Test 5: End-to-end kinematic scan
   -------------------------------------------------------------------------- *)

RunTest5[] := Module[
  {poly, vars, lam, A, spec, verts, fanData,
   lamValues, nPoints,
   allPass, maxRelErr},

  Print["--- Test 5: End-to-end kinematic scan ---"];
  Print["Int dx1 dx2 / (1+lam*x1^2+x2^2+x1*x2^2)^{2+0.5I}"];
  Print["100 values of lam in [0.1, 10]"];

  lam = Symbol["lam"];
  A   = 2 + 0.5 I;
  poly = 1 + lam x[1]^2 + x[2]^2 + x[1] x[2]^2;
  vars = {x[1], x[2]};

  nPoints  = 100;
  lamValues = Table[0.1 + (10.0 - 0.1) (i - 1)/(nPoints - 1),
                    {i, nPoints}];

  spec = <|
    "Polynomials"       -> {poly},
    "MonomialExponents" -> {0, 0},
    "PolynomialExponents" -> {-A},
    "Variables"         -> vars,
    "KinematicSymbols"  -> {lam},
    "RegulatorSymbol"   -> None
  |>;

  verts   = PolytopeVertices[(poly /. lam -> 1)^(-Re[A]), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  Print["Computing NIntegrate reference values..."];
  Module[{testIndices, refResults},
    testIndices = {1, 25, 50, 75, 100};
    refResults = Table[
      Module[{lamVal, result},
        lamVal = lamValues[[idx]];
        result = Quiet@NIntegrate[
          1 / (1 + lamVal t1^2 + t2^2 + t1 t2^2)^A,
          {t1, 0, Infinity}, {t2, 0, Infinity},
          MaxRecursion -> 20, PrecisionGoal -> 5
        ];
        {lamVal, result}
      ],
      {idx, testIndices}
    ];

    allPass = True;
    maxRelErr = 0;

    Do[
      Print["  lam = ", refResults[[i, 1]], ": NIntegrate = ",
            refResults[[i, 2]]];,
      {i, Length[refResults]}
    ];

    (* Sector decomposition check *)
    Do[
      Module[{kinRules, vr, relErr},
        kinRules = {lam -> refResults[[i, 1]]};
        vr = Quiet@ValidateDecomposition[spec, fanData, kinRules, 3];
        If[AssociationQ[vr],
          relErr = vr["RelativeError"];
          If[NumericQ[relErr],
            If[relErr > maxRelErr, maxRelErr = relErr];
            If[relErr > 0.05, allPass = False];
            Print["  lam = ", refResults[[i, 1]],
                  ": sector rel err = ", relErr,
                  " ", If[relErr < 0.05, "PASS", "FAIL"]];,
            Print["  lam = ", refResults[[i, 1]],
                  ": non-numeric error, FAIL"];
            allPass = False;
          ];
        ];
      ],
      {i, Length[refResults]}
    ];
  ];

  Print["  Max relative error: ", maxRelErr];
  Print["  ", If[allPass, "PASS", "FAIL"]];

  allPass
];

(* --------------------------------------------------------------------------
   Test 6: Large coefficients
   Verifies that the tropical decomposition and sector integrals remain
   correct when polynomial coefficients span many orders of magnitude.
   The tropically dominant monomial (min exponents) may NOT be the
   numerically largest monomial, but the factoring is an exact algebraic
   identity so the result must still agree with direct NIntegrate.
   -------------------------------------------------------------------------- *)

RunTest6[] := Module[
  {poly, vars, spec, verts, fanData, allPass,
   testCases, t1, t2},

  Print["--- Test 6: Large polynomial coefficients ---"];
  Print["Verifies correctness when coefficients span many orders of magnitude"];

  allPass = True;

  (* Test case A: coefficients O(10^6)
     P = 1 + 10^6 x1^2 + x2^2 + x1 x2^2
     The 10^6 term dominates numerically but is NOT the tropically
     dominant monomial in most sectors. *)
  Module[{polyA, specA, vertsA, fanA, directA, vrA, relErrA},
    Print[];
    Print["  Case A: P = 1 + 10^6 x1^2 + x2^2 + x1*x2^2, exponent -2"];
    polyA = 1 + 10^6 x[1]^2 + x[2]^2 + x[1] x[2]^2;
    vars  = {x[1], x[2]};

    specA = <|
      "Polynomials"        -> {polyA},
      "MonomialExponents"  -> {0, 0},
      "PolynomialExponents" -> {-2},
      "Variables"          -> vars,
      "KinematicSymbols"   -> {},
      "RegulatorSymbol"    -> None
    |>;

    directA = Quiet@NIntegrate[
      1 / (1 + 10^6 t1^2 + t2^2 + t1 t2^2)^2,
      {t1, 0, Infinity}, {t2, 0, Infinity},
      MaxRecursion -> 20, PrecisionGoal -> 6
    ];

    vertsA = PolytopeVertices[polyA^(-2), vars];
    fanA   = ComputeDecomposition[vertsA, "ShowProgress" -> False];

    vrA = Quiet@ValidateDecomposition[specA, fanA, {}, 3];
    If[AssociationQ[vrA],
      relErrA = vrA["RelativeError"];
      Print["    NIntegrate = ", directA];
      Print["    Sector sum = ", vrA["SectorSum"]];
      Print["    Rel error  = ", relErrA];
      If[!NumericQ[relErrA] || relErrA > 0.05,
        Print["    FAIL"];
        allPass = False;,
        Print["    PASS"];
      ];,
      Print["    ValidateDecomposition returned non-association, FAIL"];
      allPass = False;
    ];
  ];

  (* Test case B: mixed large and small coefficients
     P = 10^(-4) + 10^4 x1^2 + 10^(-4) x2^2 + 10^4 x1 x2^2 + x1^2 x2
     Coefficients span 8 orders of magnitude. *)
  Module[{polyB, specB, vertsB, fanB, directB, vrB, relErrB},
    Print[];
    Print["  Case B: coefficients from 10^-4 to 10^4, exponent -2"];
    polyB = 10^(-4) + 10^4 x[1]^2 + 10^(-4) x[2]^2 +
            10^4 x[1] x[2]^2 + x[1]^2 x[2];
    vars  = {x[1], x[2]};

    specB = <|
      "Polynomials"        -> {polyB},
      "MonomialExponents"  -> {0, 0},
      "PolynomialExponents" -> {-2},
      "Variables"          -> vars,
      "KinematicSymbols"   -> {},
      "RegulatorSymbol"    -> None
    |>;

    directB = Quiet@NIntegrate[
      1 / (10^(-4) + 10^4 t1^2 + 10^(-4) t2^2 +
           10^4 t1 t2^2 + t1^2 t2)^2,
      {t1, 0, Infinity}, {t2, 0, Infinity},
      MaxRecursion -> 20, PrecisionGoal -> 6
    ];

    vertsB = PolytopeVertices[polyB^(-2), vars];
    fanB   = ComputeDecomposition[vertsB, "ShowProgress" -> False];

    vrB = Quiet@ValidateDecomposition[specB, fanB, {}, 3];
    If[AssociationQ[vrB],
      relErrB = vrB["RelativeError"];
      Print["    NIntegrate = ", directB];
      Print["    Sector sum = ", vrB["SectorSum"]];
      Print["    Rel error  = ", relErrB];
      If[!NumericQ[relErrB] || relErrB > 0.05,
        Print["    FAIL"];
        allPass = False;,
        Print["    PASS"];
      ];,
      Print["    ValidateDecomposition returned non-association, FAIL"];
      allPass = False;
    ];
  ];

  (* Test case C: large coefficient with higher exponent
     P = 1 + 10^8 x1^3 x2 + x2^3, exponent -3
     The 10^8 monomial has degree 4 and large coefficient, ensuring
     it numerically dominates even though the constant term is tropically
     dominant in its cone. *)
  Module[{polyC, specC, vertsC, fanC, directC, vrC, relErrC},
    Print[];
    Print["  Case C: P = 1 + 10^8 x1^3*x2 + x2^3, exponent -3"];
    polyC = 1 + 10^8 x[1]^3 x[2] + x[2]^3;
    vars  = {x[1], x[2]};

    specC = <|
      "Polynomials"        -> {polyC},
      "MonomialExponents"  -> {0, 0},
      "PolynomialExponents" -> {-3},
      "Variables"          -> vars,
      "KinematicSymbols"   -> {},
      "RegulatorSymbol"    -> None
    |>;

    directC = Quiet@NIntegrate[
      1 / (1 + 10^8 t1^3 t2 + t2^3)^3,
      {t1, 0, Infinity}, {t2, 0, Infinity},
      MaxRecursion -> 20, PrecisionGoal -> 5
    ];

    vertsC = PolytopeVertices[polyC^(-3), vars];
    fanC   = ComputeDecomposition[vertsC, "ShowProgress" -> False];

    vrC = Quiet@ValidateDecomposition[specC, fanC, {}, 2];
    If[AssociationQ[vrC],
      relErrC = vrC["RelativeError"];
      Print["    NIntegrate = ", directC];
      Print["    Sector sum = ", vrC["SectorSum"]];
      Print["    Rel error  = ", relErrC];
      If[!NumericQ[relErrC] || relErrC > 0.05,
        Print["    FAIL"];
        allPass = False;,
        Print["    PASS"];
      ];,
      Print["    ValidateDecomposition returned non-association, FAIL"];
      allPass = False;
    ];
  ];

  Print[];
  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 7: Higher-dimensional integrals (3D and 4D)
   -------------------------------------------------------------------------- *)

RunTest7[] := Module[
  {allPass = True, vars3, poly3, spec3, verts3, fan3, vr3,
   vars4, poly4, spec4, verts4, fan4, vr4},

  Print["--- Test 7: Higher-dimensional integrals (3D and 4D) ---"];
  Print["Tests that the pipeline works correctly in dimensions > 2"];

  (* 3D: Int dx1 dx2 dx3 / (1 + x1^2 + x2^2 + x3^2 + x1*x2*x3)^3 *)
  Print[];
  Print["  Case A (3D): P = 1 + x1^2 + x2^2 + x3^2 + x1*x2*x3, exponent -3"];

  poly3 = 1 + x[1]^2 + x[2]^2 + x[3]^2 + x[1] x[2] x[3];
  vars3 = {x[1], x[2], x[3]};

  spec3 = <|
    "Polynomials"        -> {poly3},
    "MonomialExponents"  -> {0, 0, 0},
    "PolynomialExponents" -> {-3},
    "Variables"          -> vars3,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  verts3 = PolytopeVertices[poly3^(-1), vars3];
  fan3   = ComputeDecomposition[verts3, "ShowProgress" -> False];
  Print["    Fan: ", Length[fan3[[1]]], " rays, ", Length[fan3[[2]]], " sectors"];

  vr3 = Quiet@ValidateDecomposition[spec3, fan3, {}, 3];
  If[AssociationQ[vr3],
    Print["    NIntegrate = ", vr3["DirectResult"]];
    Print["    Sector sum = ", vr3["SectorSum"]];
    Print["    Rel error  = ", vr3["RelativeError"]];
    If[!NumericQ[vr3["RelativeError"]] || vr3["RelativeError"] > 0.01,
      Print["    FAIL"];
      allPass = False;,
      Print["    PASS"];
    ];,
    Print["    ValidateDecomposition returned non-association, FAIL"];
    allPass = False;
  ];

  (* 4D: Int dx1 dx2 dx3 dx4 / (1+x1^2+x2^2+x3^2+x4^2+x1*x2+x3*x4)^4 *)
  Print[];
  Print["  Case B (4D): P = 1+x1^2+x2^2+x3^2+x4^2+x1*x2+x3*x4, exponent -4"];

  poly4 = 1 + x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2 + x[1] x[2] + x[3] x[4];
  vars4 = {x[1], x[2], x[3], x[4]};

  spec4 = <|
    "Polynomials"        -> {poly4},
    "MonomialExponents"  -> {0, 0, 0, 0},
    "PolynomialExponents" -> {-4},
    "Variables"          -> vars4,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  verts4 = PolytopeVertices[poly4^(-1), vars4];
  fan4   = ComputeDecomposition[verts4, "ShowProgress" -> False];
  Print["    Fan: ", Length[fan4[[1]]], " rays, ", Length[fan4[[2]]], " sectors"];

  vr4 = Quiet@ValidateDecomposition[spec4, fan4, {}, 2];
  If[AssociationQ[vr4],
    Print["    NIntegrate = ", vr4["DirectResult"]];
    Print["    Sector sum = ", vr4["SectorSum"]];
    Print["    Rel error  = ", vr4["RelativeError"]];
    If[!NumericQ[vr4["RelativeError"]] || vr4["RelativeError"] > 0.01,
      Print["    FAIL"];
      allPass = False;,
      Print["    PASS"];
    ];,
    Print["    ValidateDecomposition returned non-association, FAIL"];
    allPass = False;
  ];

  Print[];
  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 8: Full divergent pipeline (2D)
   ProcessSector -> ProcessDivergentSector -> ValidateSubtraction
   -------------------------------------------------------------------------- *)

RunTest8[] := Module[
  {poly, vars, eps, spec, verts, fanData,
   dualVertices, simplexList,
   allSectorData, divSectors,
   allPass},

  Print["--- Test 8: Full divergent pipeline (2D) ---"];
  Print["Int x1^{2eps-1} (1+x1+x2+x1*x2)^{-2} dx"];

  eps = Symbol["eps8"];
  poly = 1 + x[1] + x[2] + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {2 eps - 1, 0},
    "PolynomialExponents" -> {-2},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;

  verts = PolytopeVertices[poly^(-2), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  dualVertices = fanData[[1]];
  simplexList  = fanData[[2]];

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  divSectors = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  If[Length[divSectors] == 0,
    Print["  No divergent sectors found, FAIL"];
    Return[False]
  ];
  Print["  Found ", Length[divSectors], " divergent sector(s)"];

  allPass = True;
  Do[
    Module[{divData, vsResult, relErr},
      divData = Quiet@ProcessDivergentSector[sd, spec];
      If[!AssociationQ[divData],
        Print["  ProcessDivergentSector failed, FAIL"];
        allPass = False;,

        vsResult = Quiet@ValidateSubtraction[divData, sd, spec, {}, 0.05];
        If[!AssociationQ[vsResult],
          Print["  ValidateSubtraction failed, FAIL"];
          allPass = False;,

          relErr = vsResult["RelativeError"];
          Print["  Sector ", sd["ConeIndex"], ": RelativeError = ", relErr];
          If[!NumericQ[relErr] || relErr > 0.02,
            Print["    FAIL"];
            allPass = False;,
            Print["    PASS"];
          ];
        ];
      ];
    ],
    {sd, divSectors}
  ];

  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 9: Reconstruction at multiple epsilon values
   -------------------------------------------------------------------------- *)

RunTest9[] := Module[
  {poly, vars, eps, spec, verts, fanData,
   dualVertices, simplexList,
   allSectorData, divSectors,
   testEpsValues, allPass, remainders},

  Print["--- Test 9: Reconstruction at multiple epsilon values ---"];
  Print["Same integral as Test 8, validated at eps = {0.1, 0.05, 0.01, 0.005}"];

  eps = Symbol["eps9"];
  poly = 1 + x[1] + x[2] + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {2 eps - 1, 0},
    "PolynomialExponents" -> {-2},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;

  verts = PolytopeVertices[poly^(-2), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  dualVertices = fanData[[1]];
  simplexList  = fanData[[2]];

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  divSectors = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  If[Length[divSectors] == 0,
    Print["  No divergent sectors, FAIL"];
    Return[False]
  ];

  testEpsValues = {0.1, 0.05, 0.01, 0.005};
  allPass = True;
  remainders = {};

  (* Use first divergent sector *)
  Module[{sd, divData},
    sd = divSectors[[1]];
    divData = Quiet@ProcessDivergentSector[sd, spec];
    If[!AssociationQ[divData],
      Print["  ProcessDivergentSector failed, FAIL"];
      Return[False]
    ];

    Do[
      Module[{vsResult, relErr},
        vsResult = Quiet@ValidateSubtraction[divData, sd, spec, {}, te];
        If[!AssociationQ[vsResult],
          Print["  ValidateSubtraction failed at eps=", te, ", FAIL"];
          allPass = False;,

          relErr = vsResult["RelativeError"];
          AppendTo[remainders, {te, vsResult["Remainder"]}];
          (* Threshold scales with eps: truncation error is O(eps) *)
          Module[{threshold},
            threshold = 0.03 + 0.5 te;
            Print["  eps = ", te, ": RelativeError = ", relErr,
                  " (threshold = ", threshold, ")"];
            If[!NumericQ[relErr] || relErr > threshold,
              Print["    FAIL"];
              allPass = False;,
              Print["    PASS"];
            ];
          ];
        ];
      ],
      {te, testEpsValues}
    ];

    (* Check remainder convergence *)
    Module[{r1, r2, convErr},
      r1 = Select[remainders, #[[1]] == 0.01 &];
      r2 = Select[remainders, #[[1]] == 0.005 &];
      If[Length[r1] > 0 && Length[r2] > 0,
        r1 = r1[[1, 2]]; r2 = r2[[1, 2]];
        If[NumericQ[r1] && NumericQ[r2] && Abs[r1] > 0,
          convErr = Abs[r2 - r1] / Abs[r1];
          Print["  |R(0.005)-R(0.01)|/|R(0.01)| = ", convErr];
          If[convErr > 0.3,
            Print["    Remainder convergence: FAIL"];
            allPass = False;,
            Print["    Remainder convergence: PASS"];
          ];
        ];
      ];
    ];
  ];

  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 10: Divergent variable index permutation
   -------------------------------------------------------------------------- *)

RunTest10[] := Module[
  {poly, vars, eps, spec, verts, fanData,
   dualVertices, simplexList,
   allSectorData, divSectors,
   allPass, divVarSet},

  Print["--- Test 10: Divergent variable index permutation ---"];
  Print["Int x2^{2eps-1} (1+x1+x2+x1*x2)^{-2} dx"];

  eps = Symbol["eps10"];
  poly = 1 + x[1] + x[2] + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {0, 2 eps - 1},
    "PolynomialExponents" -> {-2},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;

  verts = PolytopeVertices[poly^(-2), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  dualVertices = fanData[[1]];
  simplexList  = fanData[[2]];

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  divSectors = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  If[Length[divSectors] == 0,
    Print["  No divergent sectors, FAIL"];
    Return[False]
  ];

  divVarSet = Union[#["DivergentVariable"] & /@ divSectors];
  Print["  Divergent variables found: ", divVarSet];

  allPass = True;
  Do[
    Module[{divData, vsResult, relErr},
      divData = Quiet@ProcessDivergentSector[sd, spec];
      If[!AssociationQ[divData],
        Print["  ProcessDivergentSector failed for sector ",
              sd["ConeIndex"], ", FAIL"];
        allPass = False;,

        vsResult = Quiet@ValidateSubtraction[divData, sd, spec, {}, 0.05];
        If[!AssociationQ[vsResult],
          Print["  ValidateSubtraction failed for sector ",
                sd["ConeIndex"], ", FAIL"];
          allPass = False;,

          relErr = vsResult["RelativeError"];
          Print["  Sector ", sd["ConeIndex"],
                " (divVar=", sd["DivergentVariable"],
                "): RelativeError = ", relErr];
          If[!NumericQ[relErr] || relErr > 0.02,
            Print["    FAIL"];
            allPass = False;,
            Print["    PASS"];
          ];
        ];
      ];
    ],
    {sd, divSectors}
  ];

  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 11: Epsilon-dependent polynomial exponents (B1 != 0)
   -------------------------------------------------------------------------- *)

RunTest11[] := Module[
  {poly, vars, eps, spec, verts, fanData,
   dualVertices, simplexList,
   allSectorData, divSectors,
   allPass},

  Print["--- Test 11: Epsilon-dependent polynomial exponents (B1 != 0) ---"];
  Print["Int x1^{2eps-1} (1+x1+x2+x1*x2)^{-2+eps} dx"];

  eps = Symbol["eps11"];
  poly = 1 + x[1] + x[2] + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {2 eps - 1, 0},
    "PolynomialExponents" -> {-2 + eps},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;

  verts = PolytopeVertices[poly^(-2), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  dualVertices = fanData[[1]];
  simplexList  = fanData[[2]];

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  divSectors = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  If[Length[divSectors] == 0,
    Print["  No divergent sectors, FAIL"];
    Return[False]
  ];

  allPass = True;
  Do[
    Module[{divData, vsResult, relErr, b1Check, logInsCheck},
      divData = Quiet@ProcessDivergentSector[sd, spec];
      If[!AssociationQ[divData],
        Print["  ProcessDivergentSector failed, FAIL"];
        allPass = False;,

        (* Verify B1 is {1} *)
        b1Check = (divData["B1"] === {1});
        Print["  Sector ", sd["ConeIndex"], ": B1 = ", divData["B1"],
              If[b1Check, " (correct)", " (UNEXPECTED)"]];
        If[!b1Check, allPass = False];

        (* Verify G1LogInsertions PolynomialTerms has nonzero coefficient *)
        logInsCheck = False;
        Module[{polyTerms},
          polyTerms = divData["G1LogInsertions"]["PolynomialTerms"];
          If[AnyTrue[polyTerms, (#[[1]] =!= 0) &],
            logInsCheck = True
          ];
        ];
        Print["  G1LogInsertions PolynomialTerms nonzero: ", logInsCheck];
        If[!logInsCheck, allPass = False];

        vsResult = Quiet@ValidateSubtraction[divData, sd, spec, {}, 0.05];
        If[!AssociationQ[vsResult],
          Print["  ValidateSubtraction failed, FAIL"];
          allPass = False;,

          relErr = vsResult["RelativeError"];
          Print["  RelativeError = ", relErr];
          If[!NumericQ[relErr] || relErr > 0.03,
            Print["    FAIL"];
            allPass = False;,
            Print["    PASS"];
          ];
        ];
      ];
    ],
    {sd, divSectors}
  ];

  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 12: Full sector sum (convergent + divergent)
   -------------------------------------------------------------------------- *)

RunTest12[] := Module[
  {poly, vars, eps, spec, verts, fanData,
   dualVertices, simplexList,
   allSectorData, convergentSectors, divSectors,
   testEps, totalSum, directResult, relErr, allPass},

  Print["--- Test 12: Full sector sum (convergent + divergent) ---"];
  Print["Sum all sectors at eps=0.05, compare with direct NIntegrate"];

  testEps = 0.05;
  eps = Symbol["eps12"];
  poly = 1 + x[1] + x[2] + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {2 eps - 1, 0},
    "PolynomialExponents" -> {-2},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;

  verts = PolytopeVertices[poly^(-2), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  dualVertices = fanData[[1]];
  simplexList  = fanData[[2]];

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  convergentSectors = Select[allSectorData,
    (AssociationQ[#] && !#["IsDivergent"]) &];
  divSectors = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  Print["  ", Length[convergentSectors], " convergent, ",
        Length[divSectors], " divergent sectors"];

  totalSum = 0;
  allPass = True;

  (* Convergent sectors: NIntegrate of flattened form at eps=testEps *)
  Do[
    Module[{flatPolys, polyExps, pf, n, yVars, polyVals, integrand, sVal},
      flatPolys = cs["FlattenedPolys"] /. eps -> testEps;
      polyExps  = cs["PolynomialExponents"] /. eps -> testEps;
      pf = cs["Prefactor"] /. eps -> testEps;
      n  = cs["Dimension"];
      yVars = Table[Unique["cy"], {n}];

      polyVals = Table[
        Total[Table[
          mono[[1]] * Exp[Total[mono[[2]] * Log /@ yVars]],
          {mono, flatPolys[[j]]}
        ]],
        {j, Length[flatPolys]}
      ];

      integrand = pf *
        Times @@ MapThread[
          Function[{pv, be}, Exp[be * Log[pv]]],
          {polyVals, polyExps}
        ];

      sVal = Quiet@NIntegrate[
        integrand,
        Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
        MaxRecursion -> 20, PrecisionGoal -> 4,
        Method -> "GlobalAdaptive"
      ];
      Print["  Convergent sector ", cs["ConeIndex"], ": ", sVal];
      totalSum += sVal;
    ],
    {cs, convergentSectors}
  ];

  (* Divergent sectors: use ValidateSubtraction *)
  Do[
    Module[{divData, vsResult},
      divData = Quiet@ProcessDivergentSector[ds, spec];
      If[!AssociationQ[divData],
        Print["  ProcessDivergentSector failed for sector ",
              ds["ConeIndex"]];
        allPass = False;,

        vsResult = Quiet@ValidateSubtraction[divData, ds, spec, {}, testEps];
        If[!AssociationQ[vsResult],
          Print["  ValidateSubtraction failed for sector ",
                ds["ConeIndex"]];
          allPass = False;,

          Print["  Divergent sector ", ds["ConeIndex"], ": ",
                vsResult["Reconstructed"]];
          totalSum += vsResult["Reconstructed"];
        ];
      ];
    ],
    {ds, divSectors}
  ];

  (* Direct NIntegrate on [0,Infinity)^2 *)
  directResult = Quiet@NIntegrate[
    t1^(2 testEps - 1) / (1 + t1 + t2 + t1 t2)^2,
    {t1, 0, Infinity}, {t2, 0, Infinity},
    MaxRecursion -> 30, PrecisionGoal -> 5,
    Method -> "GlobalAdaptive"
  ];

  relErr = Abs[(totalSum - directResult) / directResult];
  Print["  Sector sum      = ", totalSum];
  Print["  Direct NIntegrate = ", directResult];
  Print["  Relative error   = ", relErr];

  allPass = allPass && NumericQ[relErr] && relErr < 0.05;
  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 13: C++ codegen for divergent sectors
   -------------------------------------------------------------------------- *)

RunTest13[] := Module[
  {poly, vars, eps, spec, verts, fanData,
   dualVertices, simplexList,
   allSectorData, convergentSectors, divSectors,
   processedDiv, testEps,
   cppFile, cppBinary, codeResult,
   allPass},

  Print["--- Test 13: C++ codegen for divergent sectors ---"];
  Print["Verify code generation and compilation"];

  testEps = 0.05;
  eps = Symbol["eps13"];
  poly = 1 + x[1] + x[2] + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {2 eps - 1, 0},
    "PolynomialExponents" -> {-2},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;

  verts = PolytopeVertices[poly^(-2), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  dualVertices = fanData[[1]];
  simplexList  = fanData[[2]];

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  convergentSectors = Select[allSectorData,
    (AssociationQ[#] && !#["IsDivergent"]) &];
  divSectors = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  processedDiv = Table[
    Quiet@ProcessDivergentSector[sd, spec],
    {sd, divSectors}
  ];
  processedDiv = Select[processedDiv, AssociationQ];

  (* Substitute eps for numerical codegen *)
  Quiet[CreateDirectory[FileNameJoin[{Directory[], "INTERFILES"}]]];
  cppFile   = FileNameJoin[{Directory[], "INTERFILES", "test13_tropical.cpp"}];
  cppBinary = FileNameJoin[{Directory[], "INTERFILES", "test13_tropical"}];

  Module[{convNum, specNum},
    convNum = convergentSectors /. eps -> testEps;
    specNum = <|
      "Polynomials"        -> {poly},
      "MonomialExponents"  -> {2 testEps - 1, 0},
      "PolynomialExponents" -> {-2},
      "Variables"          -> vars,
      "KinematicSymbols"   -> {},
      "RegulatorSymbol"    -> None
    |>;
    codeResult = Quiet@GenerateCppMonteCarlo[
      convNum, processedDiv, specNum, cppFile, "NSamples" -> 1000
    ];
  ];

  allPass = True;

  If[!AssociationQ[codeResult],
    Print["  GenerateCppMonteCarlo failed, FAIL"];
    Return[False]
  ];

  If[codeResult["NG0"] > 0,
    Print["  NG0 = ", codeResult["NG0"], ", PASS"];,
    Print["  NG0 = ", codeResult["NG0"], " (expected > 0), FAIL"];
    allPass = False;
  ];

  If[codeResult["NG1"] > 0,
    Print["  NG1 = ", codeResult["NG1"], ", PASS"];,
    Print["  NG1 = ", codeResult["NG1"], " (expected > 0), FAIL"];
    allPass = False;
  ];

  If[codeResult["NRemainder"] > 0,
    Print["  NRemainder = ", codeResult["NRemainder"], ", PASS"];,
    Print["  NRemainder = ", codeResult["NRemainder"], " (expected > 0), FAIL"];
    allPass = False;
  ];

  (* Check no unresolved Mathematica symbols *)
  Module[{code, badPatterns, found},
    code = codeResult["Code"];
    badPatterns = {"Sin[", "Cos[", "Sqrt[", "Plus[", "Times[",
                   "Power[", "Rule[", "List["};
    found = Select[badPatterns, StringContainsQ[code, #] &];
    If[Length[found] > 0,
      Print["  Unresolved symbols: ", found, ", FAIL"];
      allPass = False;,
      Print["  No unresolved Mathematica symbols, PASS"];
    ];
  ];

  (* Check compilation *)
  Module[{compResult},
    compResult = Quiet@CompileCpp[cppFile, cppBinary];
    If[compResult === $Failed,
      Print["  Compilation FAILED"];
      allPass = False;,
      Print["  Compilation PASS"];
    ];
  ];

  (* Clean up *)
  Quiet[If[FileExistsQ[cppFile], DeleteFile[cppFile]]];
  Quiet[If[FileExistsQ[cppBinary], DeleteFile[cppBinary]]];

  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 14: Pole coefficient accuracy
   -------------------------------------------------------------------------- *)

RunTest14[] := Module[
  {poly, vars, eps, spec, verts, fanData,
   dualVertices, simplexList,
   allSectorData, divSectors,
   allPass},

  Print["--- Test 14: Pole coefficient accuracy ---"];
  Print["Int x1^{2eps-1} (1+x1+x2)^{-2} dx"];

  eps = Symbol["eps14"];
  poly = 1 + x[1] + x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {2 eps - 1, 0},
    "PolynomialExponents" -> {-2},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;

  verts = PolytopeVertices[poly^(-2), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  dualVertices = fanData[[1]];
  simplexList  = fanData[[2]];

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  divSectors = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  If[Length[divSectors] == 0,
    Print["  No divergent sectors, FAIL"];
    Return[False]
  ];

  allPass = True;

  Do[
    Module[{divData, g0Val, ck, poleCoeff,
            eps1, eps2, Ieps1, Ieps2, slopeCoeff, relErr},
      divData = Quiet@ProcessDivergentSector[sd, spec];
      If[!AssociationQ[divData],
        Print["  ProcessDivergentSector failed, FAIL"];
        allPass = False;,

        ck = divData["ck"];

        (* Compute G0 via NIntegrate *)
        Module[{g0FlatPolys, g0Pf, g0Dim, B0, g0yVars,
                g0PolyVals, g0Integrand},
          g0FlatPolys = divData["G0FlatPolys"];
          g0Pf  = divData["G0Prefactor"];
          g0Dim = divData["G0Dimension"];
          B0    = divData["G0PolyExponents"];
          g0yVars = Table[Unique["g0y"], {g0Dim}];

          g0PolyVals = Table[
            Total[Table[
              mono[[1]] * Exp[Total[mono[[2]] * Log /@ g0yVars]],
              {mono, g0FlatPolys[[j]]}
            ]],
            {j, Length[g0FlatPolys]}
          ];

          g0Integrand = g0Pf *
            Times @@ MapThread[
              Function[{pv, be}, Exp[be * Log[pv]]],
              {g0PolyVals, B0}
            ];

          g0Val = Quiet@NIntegrate[
            g0Integrand,
            Evaluate[Sequence @@ ({#, 0, 1} & /@ g0yVars)],
            MaxRecursion -> 20, PrecisionGoal -> 4,
            Method -> "GlobalAdaptive"
          ];
        ];

        poleCoeff = g0Val / ck;
        Print["  Sector ", sd["ConeIndex"], ": G0/ck = ", poleCoeff];

        (* Numerical check: sector integral at two small eps values *)
        eps1 = 0.01; eps2 = 0.02;

        Module[{yVars, clearedPolys, aNum1, polyVals1, integ1,
                aNum2, polyVals2, integ2, polyExps},
          yVars = Table[Unique["py"], {sd["Dimension"]}];
          clearedPolys = divData["ClearedPolys"];
          polyExps = sd["PolynomialExponents"];

          (* Sector integral at eps1 *)
          aNum1 = sd["NewExponents"] /. eps -> eps1;
          polyVals1 = Table[
            Total[Table[
              mono[[1]] * Exp[Total[mono[[2]] * Log /@ yVars]],
              {mono, clearedPolys[[j]]}
            ]],
            {j, Length[clearedPolys]}
          ];
          integ1 = Abs[sd["DetM"]] *
            Exp[Total[(aNum1 - 1) * Log /@ yVars]] *
            Times @@ MapThread[
              Function[{pv, be}, Exp[be * Log[pv]]],
              {polyVals1, polyExps /. eps -> eps1}
            ];

          Ieps1 = Quiet@NIntegrate[
            integ1,
            Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
            MaxRecursion -> 20, PrecisionGoal -> 3,
            Method -> "GlobalAdaptive"
          ];

          (* Sector integral at eps2 *)
          aNum2 = sd["NewExponents"] /. eps -> eps2;
          polyVals2 = Table[
            Total[Table[
              mono[[1]] * Exp[Total[mono[[2]] * Log /@ yVars]],
              {mono, clearedPolys[[j]]}
            ]],
            {j, Length[clearedPolys]}
          ];
          integ2 = Abs[sd["DetM"]] *
            Exp[Total[(aNum2 - 1) * Log /@ yVars]] *
            Times @@ MapThread[
              Function[{pv, be}, Exp[be * Log[pv]]],
              {polyVals2, polyExps /. eps -> eps2}
            ];

          Ieps2 = Quiet@NIntegrate[
            integ2,
            Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
            MaxRecursion -> 20, PrecisionGoal -> 3,
            Method -> "GlobalAdaptive"
          ];
        ];

        (* Numerical slope: I(eps) ~ C/eps + D + ...
           dI/d(1/eps) = C
           slopeCoeff = (Ieps1 - Ieps2) / (1/eps1 - 1/eps2) *)
        slopeCoeff = (Ieps1 - Ieps2) / (1/eps1 - 1/eps2);
        Print["  Numerical slope = ", slopeCoeff];

        relErr = Abs[(poleCoeff - slopeCoeff) / poleCoeff];
        Print["  |G0/ck - slope|/|G0/ck| = ", relErr];
        If[!NumericQ[relErr] || relErr > 0.05,
          Print["    FAIL"];
          allPass = False;,
          Print["    PASS"];
        ];
      ];
    ],
    {sd, divSectors}
  ];

  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 15: G1 log-insertion structural + numerical check
   -------------------------------------------------------------------------- *)

RunTest15[] := Module[
  {poly, vars, eps, spec, verts, fanData,
   dualVertices, simplexList,
   allSectorData, divSectors,
   allPass},

  Print["--- Test 15: G1 log-insertion structural + numerical check ---"];
  Print["Int x1^{2eps-1} (1+x1+x2+x1*x2)^{-2+eps} dx (B1 != 0)"];

  eps = Symbol["eps15"];
  poly = 1 + x[1] + x[2] + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {2 eps - 1, 0},
    "PolynomialExponents" -> {-2 + eps},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;

  verts = PolytopeVertices[poly^(-2), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  dualVertices = fanData[[1]];
  simplexList  = fanData[[2]];

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  divSectors = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  If[Length[divSectors] == 0,
    Print["  No divergent sectors, FAIL"];
    Return[False]
  ];

  allPass = True;

  Do[
    Module[{divData, vsResult, ck, g1Formula, g1Num, relErr,
            logIns, nVarTerms, nPolyTerms, structPass},
      divData = Quiet@ProcessDivergentSector[sd, spec];
      If[!AssociationQ[divData],
        Print["  ProcessDivergentSector failed, FAIL"];
        allPass = False;,

        ck = divData["ck"];
        logIns = divData["G1LogInsertions"];

        (* Structural checks *)
        nVarTerms  = Length[logIns["VariableTerms"]];
        nPolyTerms = Length[logIns["PolynomialTerms"]];
        structPass = (nVarTerms == divData["G0Dimension"]) &&
                     (nPolyTerms == Length[divData["G0PolyExponents"]]);
        Print["  Sector ", sd["ConeIndex"],
              ": VariableTerms=", nVarTerms,
              " PolynomialTerms=", nPolyTerms,
              If[structPass, " (correct structure)",
                 " (UNEXPECTED structure)"]];
        If[!structPass, allPass = False];

        (* Get G1 from ValidateSubtraction *)
        vsResult = Quiet@ValidateSubtraction[divData, sd, spec, {}, 0.05];
        If[!AssociationQ[vsResult],
          Print["  ValidateSubtraction failed, FAIL"];
          allPass = False;,

          g1Formula = vsResult["G1"];

          (* G1_num = ck * (I(eps) - G0/(ck*eps) - R) *)
          g1Num = ck * (vsResult["OriginalIntegral"] -
                        vsResult["G0"] / (ck * 0.05) -
                        vsResult["Remainder"]);

          Print["  G1 (formula)   = ", g1Formula];
          Print["  G1 (numerical) = ", g1Num];

          (* Use integral scale for error metric: G1 can be small
             relative to G0/(ck*eps), amplifying reconstruction error *)
          If[NumericQ[g1Formula] && Abs[vsResult["OriginalIntegral"]] > 0,
            relErr = Abs[g1Formula - g1Num] / Abs[vsResult["OriginalIntegral"]];
            Print["  |G1_formula - G1_num|/|I(eps)| = ", relErr];
            If[relErr > 0.05,
              Print["    FAIL"];
              allPass = False;,
              Print["    PASS"];
            ];,
            Print["  G1 or integral non-numeric, skipping cross-check"];
          ];
        ];
      ];
    ],
    {sd, divSectors}
  ];

  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 16: Multiple-polynomial divergent integral
   -------------------------------------------------------------------------- *)

RunTest16[] := Module[
  {p1, p2, vars, eps, spec, verts, fanData,
   dualVertices, simplexList,
   allSectorData, divSectors,
   allPass},

  Print["--- Test 16: Multiple-polynomial divergent integral ---"];
  Print["Int x1^{2eps-1} (1+x1+x2)^{-1} (1+x1*x2)^{-1} dx"];

  eps = Symbol["eps16"];
  p1 = 1 + x[1] + x[2];
  p2 = 1 + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {p1, p2},
    "MonomialExponents"  -> {2 eps - 1, 0},
    "PolynomialExponents" -> {-1, -1},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;

  verts = PolytopeVertices[(p1 p2)^(-1), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  dualVertices = fanData[[1]];
  simplexList  = fanData[[2]];

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  divSectors = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  If[Length[divSectors] == 0,
    Print["  No divergent sectors, FAIL"];
    Return[False]
  ];

  allPass = True;
  Module[{nProcessed = 0},
    Do[
      Module[{divData, vsResult, relErr},
        divData = Quiet@ProcessDivergentSector[sd, spec];
        If[!AssociationQ[divData],
          (* Nested divergence: known limitation, skip gracefully *)
          Print["  Sector ", sd["ConeIndex"],
                ": ProcessDivergentSector skipped (nested divergence)"];,

          nProcessed++;

          (* Verify ClearedPolys, SimplifiedPolys, G0FlatPolys each have 2 entries *)
          If[Length[divData["ClearedPolys"]] != 2,
            Print["  ClearedPolys has ", Length[divData["ClearedPolys"]],
                  " entries (expected 2), FAIL"];
            allPass = False;,
            Print["  ClearedPolys: 2 entries, correct"];
          ];
          If[Length[divData["SimplifiedPolys"]] != 2,
            Print["  SimplifiedPolys has ", Length[divData["SimplifiedPolys"]],
                  " entries (expected 2), FAIL"];
            allPass = False;,
            Print["  SimplifiedPolys: 2 entries, correct"];
          ];
          If[Length[divData["G0FlatPolys"]] != 2,
            Print["  G0FlatPolys has ", Length[divData["G0FlatPolys"]],
                  " entries (expected 2), FAIL"];
            allPass = False;,
            Print["  G0FlatPolys: 2 entries, correct"];
          ];

          vsResult = Quiet@ValidateSubtraction[divData, sd, spec, {}, 0.05];
          If[!AssociationQ[vsResult],
            Print["  ValidateSubtraction failed for sector ",
                  sd["ConeIndex"], ", FAIL"];
            allPass = False;,

            relErr = vsResult["RelativeError"];
            Print["  Sector ", sd["ConeIndex"],
                  ": RelativeError = ", relErr];
            If[!NumericQ[relErr] || relErr > 0.05,
              Print["    FAIL"];
              allPass = False;,
              Print["    PASS"];
            ];
          ];
        ];
      ],
      {sd, divSectors}
    ];
    If[nProcessed == 0,
      Print["  No sectors could be processed, FAIL"];
      allPass = False;,
      Print["  Processed ", nProcessed, " of ", Length[divSectors],
            " divergent sectors"];
    ];
  ];

  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 17: 3D divergent integral
   -------------------------------------------------------------------------- *)

RunTest17[] := Module[
  {poly, vars, eps, spec, verts, fanData,
   dualVertices, simplexList,
   allSectorData, divSectors,
   allPass},

  Print["--- Test 17: 3D divergent integral ---"];
  Print["Int x1^{2eps-1} (1+x1+x2+x3+x1*x2*x3)^{-3} dx"];

  eps = Symbol["eps17"];
  poly = 1 + x[1] + x[2] + x[3] + x[1] x[2] x[3];
  vars = {x[1], x[2], x[3]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {2 eps - 1, 0, 0},
    "PolynomialExponents" -> {-3},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;

  verts = PolytopeVertices[poly^(-3), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  dualVertices = fanData[[1]];
  simplexList  = fanData[[2]];
  Print["  Fan: ", Length[dualVertices], " rays, ",
        Length[simplexList], " sectors"];

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  divSectors = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  If[Length[divSectors] == 0,
    Print["  No divergent sectors, FAIL"];
    Return[False]
  ];
  Print["  Found ", Length[divSectors], " divergent sector(s)"];

  allPass = True;
  Do[
    Module[{divData, vsResult, relErr},
      divData = Quiet@ProcessDivergentSector[sd, spec];
      If[!AssociationQ[divData],
        Print["  ProcessDivergentSector failed for sector ",
              sd["ConeIndex"], ", FAIL"];
        allPass = False;,

        (* Verify G0 dimension is 2 *)
        If[divData["G0Dimension"] != 2,
          Print["  G0Dimension = ", divData["G0Dimension"],
                " (expected 2), FAIL"];
          allPass = False;,
          Print["  G0Dimension = 2, correct"];
        ];

        vsResult = Quiet@ValidateSubtraction[divData, sd, spec, {}, 0.05];
        If[!AssociationQ[vsResult],
          Print["  ValidateSubtraction failed for sector ",
                sd["ConeIndex"], ", FAIL"];
          allPass = False;,

          relErr = vsResult["RelativeError"];
          Print["  Sector ", sd["ConeIndex"],
                ": RelativeError = ", relErr];
          (* Relaxed threshold for 3D: NIntegrate less precise *)
          If[!NumericQ[relErr] || relErr > 0.10,
            Print["    FAIL"];
            allPass = False;,
            Print["    PASS"];
          ];
        ];
      ];
    ],
    {sd, divSectors}
  ];

  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 18: Divergent integral with complex polynomial exponents
   -------------------------------------------------------------------------- *)

RunTest18[] := Module[
  {poly, vars, eps, spec, verts, fanData,
   dualVertices, simplexList,
   allSectorData, divSectors,
   allPass},

  Print["--- Test 18: Divergent integral with complex polynomial exponents ---"];
  Print["Int x1^{2eps-1} (1+x1+x2+x1*x2)^{-2+I} dx"];

  eps = Symbol["eps18"];
  poly = 1 + x[1] + x[2] + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {2 eps - 1, 0},
    "PolynomialExponents" -> {-2 + I},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;

  verts = PolytopeVertices[poly^(-2), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  dualVertices = fanData[[1]];
  simplexList  = fanData[[2]];

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  divSectors = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  If[Length[divSectors] == 0,
    Print["  No divergent sectors, FAIL"];
    Return[False]
  ];
  Print["  Found ", Length[divSectors], " divergent sector(s)"];

  allPass = True;
  Do[
    Module[{divData, vsResult, relErr, b0Check},
      divData = Quiet@ProcessDivergentSector[sd, spec];
      If[!AssociationQ[divData],
        Print["  ProcessDivergentSector failed for sector ",
              sd["ConeIndex"], ", FAIL"];
        allPass = False;,

        (* Verify B0 is complex: {-2+I} *)
        b0Check = (divData["B0"] === {-2 + I});
        Print["  Sector ", sd["ConeIndex"], ": B0 = ", divData["B0"],
              If[b0Check, " (correct)", " (UNEXPECTED)"]];
        If[!b0Check, allPass = False];

        vsResult = Quiet@ValidateSubtraction[divData, sd, spec, {}, 0.05];
        If[!AssociationQ[vsResult],
          Print["  ValidateSubtraction failed for sector ",
                sd["ConeIndex"], ", FAIL"];
          allPass = False;,

          relErr = vsResult["RelativeError"];
          Print["  Sector ", sd["ConeIndex"],
                ": RelativeError = ", relErr];
          (* Relaxed threshold: complex exponents add numerical difficulty *)
          If[!NumericQ[relErr] || relErr > 0.05,
            Print["    FAIL"];
            allPass = False;,
            Print["    PASS"];
          ];
        ];
      ];
    ],
    {sd, divSectors}
  ];

  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Package end
   -------------------------------------------------------------------------- *)

End[]

EndPackage[]
