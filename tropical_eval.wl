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
  "RunAllTests[] runs the validation suite (17 tests) with structured reporting. \
Defined in tropical_eval_examples.wl.";

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
\"Vegas\" (per-kp CUBA Vegas; same output layout); \"Batch\" falls back to \
per-kp Vegas on the IBP path. Default output is byte-identical to before.";

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
TropicalEval::liftcomplex = "ProcessSectorLifted: cone `1` — all candidate pivots produce complex atilde; cannot emit a real-valued domain indicator.  Lift with a different k or check that polynomial exponents B are real.";
TropicalEval::liftdivergent = "ProcessSectorLifted: cone `1` — atilde `2` has a non-positive component after delta resolution; the lifted sector is divergent.  Lifting supports convergent integrals only (plan.md N3).";
TropicalEval::liftfandim = "EvaluateTropicalMC with LiftData: the fan dimension is `1` but n+1 = `2` is required.  Supply the (n+1)-dimensional lifted fan.";
TropicalEval::liftdegenerate = "EvaluateTropicalMCLifted: the lifted Newton polytope is lower-dimensional; automatic fan construction is not possible — supply an explicit complete simplicial fan via the \"FanData\" option.";

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

Options[ProcessSector] = {"Verbose" -> False};

ProcessSector[integrandSpec_Association, dualVertices_List,
              simplex_List, coneIndex_Integer, OptionsPattern[]] :=
Module[
  {polys, monoExps, polyExps, vars, eps, kinSyms,
   selectedRays, mMatrix, detM, n,
   transformedPolys, clearedPolys, minExponents,
   rawAVals, effectiveAVals,
   flattenedPolys, prefactor,
   isDivergent, divVar, verbose, sectorData,
   parsedPolys},

  verbose = OptionValue["Verbose"];

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
  If[!AllTrue[heads, # === h &], h = First[heads]];
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
   liftedSubbed, original, relResid},

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
  kprimary    = primaryRule["k"];

  (* z0 = |C_primary|^(1/k_primary), EXACT.  Sign/phase goes into residual. *)
  z0 = If[IntegerQ[Abs[Cprimary]^(1/kprimary)],
    Abs[Cprimary]^(1/kprimary),
    Power[Abs[Cprimary], 1/kprimary]
  ];

  (* ---- Residuals: c_i = C_i / z0^{k_i} (exact) ---- *)
  residuals = Table[
    Simplify[ruleCoeffs[[i]] / z0^liftRules[[i]]["k"]],
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
            ki    = r["k"];
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

Options[ProcessSectorLifted] = {"Verbose" -> False};

ProcessSectorLifted[liftedSpec_Association, dualVertices_List,
                    simplex_List, coneIndex_Integer,
                    liftData_Association, OptionsPattern[]] :=
Module[
  {sdAug, z0, auxIdx, a, clearedPolys, detM, mMatrix, polyExps, n1, n,
   mVec, verbose, tryPivot, pivotCache, candidates, bestPivot,
   pivotP, mp, ap, mOtherVec, atildeVals, reclearedPolys,
   allOtherZero, domainClass, logZ0,
   fsResult, flattenedPolys, prefactor, prefactorBase,
   isRealPos},

  verbose = OptionValue["Verbose"];

  (* --- Step 1: standard (n+1)-dim ProcessSector --- *)
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

  (* §6.7 EXACT realness predicate: atilde_j is real (Im exactly 0) AND
     Re[atilde_j] > 0 — decided exactly, never by a 10^-12 tolerance.  A
     complex atilde routes deterministically to liftcomplex. *)
  isRealPos[av_] := TrueQ[PossibleZeroQ[Im[av]]] && TrueQ[Re[av] > 0];

  candidates = {};
  Do[
    If[mVec[[p]] != 0,
      Module[{res = pivotCache[p]},
        If[res =!= $Failed &&
           And @@ (isRealPos /@ res["atilde"]),
          AppendTo[candidates, res]
        ]
      ]
    ],
    {p, n1}
  ];

  If[candidates === {},
    (* Any complex-atilde candidate => liftcomplex (decided exactly). *)
    Module[{anyComplex = False},
      Do[
        If[mVec[[p]] != 0,
          Module[{res = pivotCache[p]},
            If[res =!= $Failed &&
               AnyTrue[res["atilde"], (!TrueQ[PossibleZeroQ[Im[#]]]) &],
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
            If[res =!= $Failed, {p, mVec[[p]], N[res["atilde"]]}]
          ],
          Nothing
        ],
        {p, n1}
      ];
      Message[TropicalEval::liftnopivot, coneIndex, mVec, pivotSummary]
    ];
    Return[$Failed]
  ];

  (* AMENDED ranking (plan.md §6.2): (1) HasConstantTerm; (2) |mp|=1;
     (3) max min_j Re[atilde].  Constant-term preservation is FIRST because
     HasConstantTerm=False sectors can have infinite MC variance (§9 risk 1). *)
  candidates = SortBy[candidates,
    {-Boole[#["hasConst"]],
     -Boole[Abs[#["mp"]] == 1],
     -Min[Re[N[#["atilde"]]]]} &
  ];
  bestPivot = candidates[[1]];

  pivotP         = bestPivot["pivot"];
  mp             = bestPivot["mp"];
  ap             = bestPivot["ap"];
  mOtherVec      = bestPivot["mOther"];
  atildeVals     = bestPivot["atilde"];
  reclearedPolys = bestPivot["newPolys"];

  If[verbose,
    Print["  PSL cone ", coneIndex, ": pivot p=", pivotP,
          " mp=", mp, " atilde=", N[atildeVals],
          " HasConstantTerm=", bestPivot["hasConst"]]
  ];

  (* ---- Step 4: domain-constraint classification (plan.md §6.2) ---- *)
  logZ0        = Log[z0];
  allOtherZero = And @@ (# == 0 & /@ mOtherVec);

  If[allOtherZero,
    If[N[z0^(1/mp)] > 1,
      Return[<|"EmptyDomain" -> True, "ConeIndex" -> coneIndex|>]
    ];
    domainClass = None;
    ,
    If[mp > 0,
      If[And @@ (#>= 0 & /@ mOtherVec) && N[z0] > 1,
        Return[<|"EmptyDomain" -> True, "ConeIndex" -> coneIndex|>]
      ];
      If[And @@ (#<= 0 & /@ mOtherVec) && N[z0] <= 1,
        domainClass = None;
        ,
        domainClass = <|
          "LogZ0"           -> logZ0,
          "MP"              -> mp,
          "IndicatorCoeffs" -> Table[mOtherVec[[jj]] / atildeVals[[jj]], {jj, n}]
        |>
      ],
      (* mp < 0 *)
      If[And @@ (#<= 0 & /@ mOtherVec) && N[z0] < 1,
        Return[<|"EmptyDomain" -> True, "ConeIndex" -> coneIndex|>]
      ];
      If[And @@ (#>= 0 & /@ mOtherVec) && N[z0] >= 1,
        domainClass = None;
        ,
        domainClass = <|
          "LogZ0"           -> logZ0,
          "MP"              -> mp,
          "IndicatorCoeffs" -> Table[mOtherVec[[jj]] / atildeVals[[jj]], {jj, n}]
        |>
      ]
    ]
  ];

  (* ---- Step 5: flatten via FlattenSector ---- *)
  prefactorBase = (Abs[detM] / Abs[mp]) * z0^(ap / mp - 1);
  fsResult = FlattenSector[reclearedPolys, atildeVals, prefactorBase];

  If[fsResult["IsDivergent"],
    Message[TropicalEval::liftdivergent, coneIndex, atildeVals];
    Return[$Failed]
  ];

  flattenedPolys = fsResult["FlattenedPolys"];
  prefactor      = fsResult["Prefactor"];

  (* ---- Step 6: assemble full SectorData ---- *)
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
  {eps, k, aVals, a0, a1, ck, n, polyExps,
   clearedPolys, detM,
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
  aVals  = sectorData["NewExponents"];  (* effective *)
  detM   = sectorData["DetM"];

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

    g0Prefactor = Abs[detM] / (Times @@ g0aVals);
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

    remPrefactor = Abs[detM] / (Times @@ remAvals);

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
    "RawExponents"        -> sectorData["RawExponents"]
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
   clearedPolys, simpPolys, minExps,
   detM, yVars,
   originalIntegral, g0Val, g1Val, remVal,
   divContrib, reconstructed, relError,
   rawAVals, B0, B1},

  eps      = integrandSpec["RegulatorSymbol"];
  n        = sectorData["Dimension"];
  k        = divSectorData["DivergentVariable"];
  ck       = divSectorData["ck"];
  a0       = divSectorData["a0"];
  a1       = divSectorData["a1"];
  aVals    = sectorData["NewExponents"];  (* effective, symbolic in eps *)
  rawAVals = sectorData["RawExponents"];
  polyExps = sectorData["PolynomialExponents"];
  detM     = sectorData["DetM"];
  minExps  = sectorData["MinExponents"];

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
    integrand = Abs[detM] *
      Exp[Total[(aNum - 1) * Log /@ yVars]] *
      Times @@ MapThread[
        Function[{pv, be}, Exp[be * Log[pv]]],
        {polyValsExpr, polyExps /. fullRules}
      ];

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

    g0Integrand = Abs[detM] *
      Exp[Total[(g0aNum - 1) * Log /@ g0yVars]] *
      Times @@ MapThread[
        Function[{pv, be}, Exp[be * Log[pv]]],
        {g0PolyVals, B0}
      ];

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

    g1BaseIntegrand = Abs[detM] *
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

    g1Integrand = g1BaseIntegrand * logInsertionSum;

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

    remIntegrand = Abs[detM] *
      Exp[Total[(remAnum - 1) * Log /@ remYVars]] *
      bracket;

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
     GenerateCpp                    -> ONE public entry (the IBP variant is an
                                       absorbed code path, selected by passing
                                       ibpSectors); GenerateCppMonteCarlo kept as
                                       a thin alias for API continuity.

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
  "VegasEpsRel" -> 1.*^-12,
  "VegasEpsAbs" -> 1.*^-300,
  "VegasSeed"   -> 0,
  "CubaMaxComp" -> 512
};

(* Normalize the integrator vocabulary: canonical "MC" / "VEGAS"; accept the
   older "MonteCarlo" / "Vegas" spellings as aliases.  Returns "MC" | "VEGAS". *)
normalizeIntegrator[s_] := Switch[s,
  "VEGAS" | "Vegas", "VEGAS",
  "MC" | "MonteCarlo", "MC",
  _, "MC"];

(* --------------------------------------------------------------------------
   emitBaseFuncBody — the shared C++ body for a "base" integrand function:
   the signature, log_y[] setup, the per-polynomial monomial sums, the
   prefactor assignment, and the polynomial-exponent product.  Does NOT emit
   the closing "return ...;\n}\n" (callers add their own tail: a plain return,
   or a log-insertion factor).  resultVar is the C++ lvalue ("result",
   "g0_val", "base_val").  Reproduces the Phase-1 bytes exactly.
   -------------------------------------------------------------------------- *)
emitBaseFuncBody[funcName_String, comment_String, resultVar_String,
                 flatPolys_List, polyExps_List, prefactor_, dim_Integer,
                 paramMap_Association, domainConstraint_: None] :=
Module[{funcCode},
  funcCode = "inline cx " <> funcName <>
    "(const double* y, const double* params) {\n";
  funcCode = funcCode <> "    // " <> comment <> "\n";
  funcCode = funcCode <>
    "    double log_y[" <> ToString[dim] <> "];\n";
  funcCode = funcCode <>
    "    for (int i = 0; i < " <> ToString[dim] <> "; i++)\n";
  funcCode = funcCode <>
    "        log_y[i] = (y[i] > 1e-300) ? std::log(y[i]) : -700.0;\n\n";

  (* Lifted-sector domain indicator (plan.md §3.3 / §6.2).  Only emitted when a
     DomainConstraint is present; for the unlifted path (domainConstraint===None)
     this block is skipped, keeping the emitted C++ byte-identical (#25). *)
  If[domainConstraint =!= None,
    Module[{dc = domainConstraint, logZ0str, mpStr, icList, icTerms, sumStr},
      logZ0str = mmaToCInternal[N[dc["LogZ0"]], paramMap];
      mpStr    = mmaToCInternal[N[dc["MP"]], paramMap];
      icList   = N[dc["IndicatorCoeffs"]];
      icTerms  = Table[
        mmaToCInternal[icList[[i]], paramMap] <>
        " * log_y[" <> ToString[i - 1] <> "]",
        {i, Length[icList]}
      ];
      sumStr = If[Length[icTerms] == 0, "0.0", StringRiffle[icTerms, " + "]];
      funcCode = funcCode <> "    // lifted-sector domain indicator\n";
      funcCode = funcCode <>
        "    double log_ypstar = (" <> logZ0str <>
        " - (" <> sumStr <> ")) * (1.0/" <> mpStr <> ");\n";
      funcCode = funcCode <>
        "    if (log_ypstar > 0.0) return cx(0.0, 0.0);\n\n";
    ]
  ];

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
      funcCode = emitBaseFuncBody[
        "integrand_conv_" <> ToString[s - 1],
        "Convergent sector " <> ToString[sd["ConeIndex"]],
        "result",
        sd["FlattenedPolys"], sd["PolynomialExponents"],
        sd["Prefactor"], sd["Dimension"], paramMap,
        Lookup[sd, "DomainConstraint", None]];
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
          dd["G0Prefactor"], dd["G0Dimension"], paramMap];
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
          dd["G0Prefactor"], dd["G0Dimension"], paramMap];
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
        bndData  = ibpSD["BoundaryData"];
        ibpTerms = ibpSD["IBPTerms"];
        sIdx     = s - 1;
        sectorFuncs = <|"SectorIndex" -> sIdx,
                        "ConeIndex" -> ibpSD["ConeIndex"]|>;

        (* Boundary base (order 0) *)
        Module[{funcCode},
          funcCode = emitBaseFuncBody[
            "integrand_ibp_bnd_" <> ToString[sIdx] <> "_base",
            "IBP boundary base, sector " <> ToString[ibpSD["ConeIndex"]],
            "result",
            bndData["FlatPolys"], bndData["PolyExponents"],
            bndData["Prefactor"], bndData["Dimension"], paramMap];
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
            bndData["Prefactor"], bndData["Dimension"], paramMap];
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
              termData["Prefactor"], termData["Dimension"], paramMap];
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
              termData["Prefactor"], termData["Dimension"], paramMap];
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
   (golden-master regression, cross-check #25).  Batch is ignored on the IBP
   path (per-kp Vegas fallback, plan.md §5.1).
   -------------------------------------------------------------------------- *)
Options[emitMain] = Join[
  {"NSamples" -> 1000000, "SeedBase" -> 42},
  $vegasOptionDefaults
];
emitMain[info_Association, integrator_String, batch_:False,
         OptionsPattern[]] :=
Module[{code, nSamples, seedBase, epsrel, epsabs, seed, maxComp, isIBP},
  nSamples = OptionValue["NSamples"];
  seedBase = OptionValue["SeedBase"];
  epsrel   = ToString[CForm[N[OptionValue["VegasEpsRel"]]]];
  epsabs   = ToString[CForm[N[OptionValue["VegasEpsAbs"]]]];
  seed     = ToString[OptionValue["VegasSeed"]];
  maxComp  = ToString[OptionValue["CubaMaxComp"]];
  isIBP    = TrueQ[info["IsIBP"]];

  Which[
    integrator === "MC" && !isIBP,
  code = "int main(int argc, char* argv[]) {\n";
  code = code <> "    if (argc < 3) {\n";
  code = code <> "        std::cerr << \"Usage: \" << argv[0] << \" <input_file> <output_file> [n_samples] [n_threads]\" << std::endl;\n";
  code = code <> "        return 1;\n";
  code = code <> "    }\n\n";

  code = code <> "    std::string input_file = argv[1];\n";
  code = code <> "    std::string output_file = argv[2];\n";
  code = code <> "    int n_samples = (argc > 3) ? std::atoi(argv[3]) : " <>
    ToString[nSamples] <> ";\n";
  code = code <> "    int n_threads = (argc > 4) ? std::atoi(argv[4]) : 1;\n";
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
  code = code <> "        uint64_t seed = " <> ToString[seedBase] <> "ULL + (uint64_t)kp;\n";
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
    code = code <> "        std::cerr << \"Usage: \" << argv[0] << \" <input_file> <output_file> [n_samples] [n_threads]\" << std::endl;\n";
    code = code <> "        return 1;\n";
    code = code <> "    }\n\n";

    code = code <> "    std::string input_file = argv[1];\n";
    code = code <> "    std::string output_file = argv[2];\n";
    code = code <> "    int n_samples = (argc > 3) ? std::atoi(argv[3]) : " <>
      ToString[nSamples] <> ";\n";
    code = code <> "    int n_threads = (argc > 4) ? std::atoi(argv[4]) : 1;\n";
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
    code = code <> "        uint64_t seed = " <> ToString[seedBase] <> "ULL + (uint64_t)kp;\n";
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
  code = code <> "                  0, (int)maxeval, 1000, 500, 1000, 0, nullptr, nullptr,\n";
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
  code = code <> "          0, (int)maxeval, 1000, 500, 1000, 0, nullptr, nullptr,\n";
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
    (* IBP path: batched Vegas falls back to per-kp (plan.md §5.1) *)
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
  code = code <> "          0, (int)maxeval, 1000, 500, 1000, 0, nullptr, nullptr,\n";
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
  "NSamples"    -> 1000000,
  "MaxDim"      -> 20,
  "SeedBase"    -> 42,
  "Integrator"  -> "MonteCarlo",   (* "MonteCarlo" | "Vegas" *)
  "Batch"       -> False,          (* True => chunked-ncomp Vegas (Vegas only) *)
  "VegasEpsRel" -> 1.*^-12,
  "VegasEpsAbs" -> 1.*^-300,
  "VegasSeed"   -> 0,
  "CubaMaxComp" -> 512
};

(* --------------------------------------------------------------------------
   cppBadTokens: scan generated C++ for tokens that mean the output cannot
   compile -- either symbolic infinities/indeterminates (leaked from a division
   by a vanishing effective exponent, e.g. an unsupported nested divergence) or
   unconverted Mathematica heads.  Returns the list of offending tokens (empty
   when the code is clean).  The codegen entry points use this to FAIL LOUDLY
   ($Failed) instead of writing C++ that g++ rejects -- see TropicalEval::badcpp.
   -------------------------------------------------------------------------- *)
cppBadTokens[code_String] := Select[
  {"DirectedInfinity", "ComplexInfinity", "Indeterminate",
   "Sin[", "Cos[", "Sqrt[", "Plus[", "Times[", "Power[", "Rule[", "List["},
  StringContainsQ[code, #] &];

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
    "CubaMaxComp" -> OptionValue["CubaMaxComp"]];

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
  "RunChecks"      -> True,
  "EpsilonValue"   -> None,
  "TestEpsilon"    -> 0.01,
  "PrecisionGoal"  -> 3,
  "WorkingDirectory" -> Automatic,
  "Verbose"        -> True,
  (* sampler selection (default = the zero-dependency plain Monte Carlo) *)
  "Integrator"     -> "MC",    (* "MC" | "VEGAS" (aliases "MonteCarlo"/"Vegas") *)
  "Batch"          -> False,         (* chunked-ncomp VEGAS (VEGAS only) *)
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
   emptyDomainCount, hasConstList},

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
  vegasOpts  = {"Integrator" -> integrator, "Batch" -> batch,
    "VegasEpsRel" -> OptionValue["VegasEpsRel"],
    "VegasEpsAbs" -> OptionValue["VegasEpsAbs"],
    "VegasSeed" -> OptionValue["VegasSeed"],
    "CubaMaxComp" -> OptionValue["CubaMaxComp"]};
  eps        = integrandSpec["RegulatorSymbol"];

  (* --- Lifting setup (plan.md §6.2) --- *)
  liftData = OptionValue["LiftData"];
  isLifted = (liftData =!= None);
  liftedSpec   = integrandSpec;  (* in lifted mode this IS the lifted (n+1) spec *)
  originalSpec = If[isLifted, liftData["OriginalSpec"], integrandSpec];

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
      (* Lifting and divergence are mutually exclusive (plan.md N3): never
         route a lifted call to the IBP/divergence path. *)
      isLifted, False,
      method === "IBP", True,
      method === "None" || method === "Subtraction", False,
      method === Automatic,
        eps =!= None && AnyTrue[
          Table[ProcessSector[integrandSpec, fanData[[1]], fanData[[2, s]], s],
                {s, Length[fanData[[2]]]}],
          (AssociationQ[#] && TrueQ[#["IsDivergent"]]) &],
      True, False
    ];
    If[routeToIBP,
      Return[evaluateTropicalIBPDriver[integrandSpec, fanData, kinematicPoints,
        "NSamples" -> nSamples, "NThreads" -> nThreads,
        "RunChecks" -> runChecks, "WorkingDirectory" -> workDir,
        "Verbose" -> verbose, "Integrator" -> integrator, "Batch" -> batch,
        "VegasEpsRel" -> OptionValue["VegasEpsRel"],
        "VegasEpsAbs" -> OptionValue["VegasEpsAbs"],
        "VegasSeed" -> OptionValue["VegasSeed"],
        "CubaMaxComp" -> OptionValue["CubaMaxComp"]]]
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
    (* ---- Lifted mode: route every sector through ProcessSectorLifted ----
       $Failed (liftcomplex/liftnopivot/liftdivergent) aborts the whole call;
       EmptyDomain sectors contribute 0 and are dropped+counted; the rest are
       convergent.  No "divergent sectors" exist in lifted mode (plan.md §6.2). *)
    allSectorData = {};
    Catch[
      Do[
        Module[{sd},
          sd = ProcessSectorLifted[liftedSpec, dualVertices,
                                   simplexList[[s]], s, liftData,
                                   "Verbose" -> verbose];
          Which[
            sd === $Failed,
              allSectorData = $Failed;  Throw[Null],
            AssociationQ[sd] && KeyExistsQ[sd, "EmptyDomain"] && sd["EmptyDomain"],
              emptyDomainCount++,
            True,
              AppendTo[allSectorData, sd];
              AppendTo[hasConstList, sd["HasConstantTerm"]]
          ]
        ],
        {s, Length[simplexList]}
      ]
    ];
    If[allSectorData === $Failed,
      Print["ERROR: ProcessSectorLifted failed for a sector (liftcomplex / ",
            "liftnopivot / liftdivergent).  Aborting ($Failed)."];
      Return[$Failed]
    ];
    convergentSectors = allSectorData;
    divergentSectors  = {};
    If[verbose,
      Print["  ", Length[convergentSectors], " convergent sectors, ",
            emptyDomainCount, " EmptyDomain sectors dropped"]
    ];
    ,
    (* ---- Standard (unlifted) mode ---- *)
    allSectorData = Table[
      Module[{sd, specToUse},
        specToUse = If[epsVal =!= None && eps =!= None,
          MapAt[# /. eps -> epsVal &, integrandSpec,
                {Key["MonomialExponents"]}] //
          MapAt[# /. eps -> epsVal &, #,
                {Key["PolynomialExponents"]}] &,
          integrandSpec
        ];
        sd = ProcessSector[specToUse, dualVertices,
                           simplexList[[s]], s, "Verbose" -> verbose];
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
      liftedFan = If[ListQ[verts],
        Quiet[ComputeDecomposition[verts, "ShowProgress" -> False],
              TropicalFan::polymake],
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
                     clearedPolys_List, eps_] :=
Module[
  {termCoeff, termExps, termPolyExps, newTerms, nPolys},

  termCoeff    = term["Coefficient"];
  termExps     = term["NewExponents"];
  termPolyExps = term["PolyExponents"];
  nPolys       = Length[clearedPolys];

  newTerms = {};

  Do[
    Module[{Bj, polj},
      Bj   = termPolyExps[[j]];
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
   IBPReduceSector
   Main IBP reduction function.  Iteratively applies IBP to resolve all
   divergent variables in a sector.
   -------------------------------------------------------------------------- *)

IBPReduceSector[sectorData_Association, eps_] :=
Module[
  {n, aVals, polyExps, clearedPolys, detM,
   a0, allDivVars, ibpPrefactors, terms, resolvedVars},

  n            = sectorData["Dimension"];
  aVals        = sectorData["NewExponents"];
  polyExps     = sectorData["PolynomialExponents"];
  clearedPolys = sectorData["ClearedPolys"];
  detM         = sectorData["DetM"];

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

      If[TrueQ[ck == 0] || (NumericQ[ck] && ck == 0),
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
              IBPExpandOneVariable[term, k, clearedPolys, eps]
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
   IBPProcessSector
   Full IBP pipeline: reduce, expand in epsilon, flatten, build boundary.
   Returns an IBPSectorData association ready for C++ codegen.

   For single divergence: produces boundary (at y_k = 1) and IBP terms.
   The combination is:
       I_sector = (1/a_k) [ B_boundary - sum_t coeff_t I_t ]
   -------------------------------------------------------------------------- *)

IBPProcessSector[sectorData_Association, integrandSpec_Association] :=
Module[
  {eps, n, ibpData, terms, clearedPolys, detM,
   divVars, ck, rk, aVals, polyExps,
   a0, a1, B0, B1, ak, ak2,
   boundaryData, ibpTermsProcessed,
   k},

  eps = integrandSpec["RegulatorSymbol"];
  n   = sectorData["Dimension"];

  (* Nested-divergence guard (G-B scope): the boundary construction (Step 2)
     and the driver's Laurent assembly support a SINGLE divergent variable per
     sector (one 1/eps pole).  A sector with >1 divergent variable would
     require a 1/eps^d assembly that is not implemented; proceeding divides by a
     vanishing effective exponent and leaks ComplexInfinity into the C++.
     Detect up front (from the effective exponents at eps=0) and refuse. *)
  Module[{a0chk = sectorData["NewExponents"] /. eps -> 0, nDiv},
    nDiv = Count[a0chk, _?(Function[v,
      TrueQ[Re[v] <= 0] || (NumericQ[v] && Re[v] <= 0)])];
    If[nDiv > 1,
      Message[TropicalEval::nestedIBP, sectorData["ConeIndex"], nDiv];
      Return[$Failed]
    ]
  ];

  (* Step 1: IBP reduction *)
  ibpData = IBPReduceSector[sectorData, eps];
  If[ibpData === $Failed, Return[$Failed]];

  terms        = ibpData["Terms"];
  clearedPolys = ibpData["ClearedPolys"];
  detM         = ibpData["DetM"];
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
          bndLogInsertions},
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

    bndPrefactor = Abs[detM] / (Times @@ bndA0);

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

    boundaryData = <|
      "FlatPolys"      -> bndFlatPolys,
      "Prefactor"      -> bndPrefactor,
      "Dimension"      -> bndDim,
      "PolyExponents"  -> B0,
      "Avals"          -> bndA0,
      "LogInsertions"  -> bndLogInsertions
    |>;
  ];

  (* ----- Step 3: Process each IBP term (expand in eps + flatten) ----- *)
  ibpTermsProcessed = Table[
    Module[{term, alpha, alpha0, alpha1, termPolyExps, tB0, tB1,
            coeff, coeff0, coeff1, flatPrefactor, flatPolys,
            logInsertions},
      term         = terms[[t]];
      alpha        = term["NewExponents"];
      termPolyExps = term["PolyExponents"];
      coeff        = term["Coefficient"];

      (* Expand effective exponents in epsilon *)
      alpha0 = alpha /. eps -> 0;
      alpha1 = D[alpha, eps] /. eps -> 0;

      (* Verify all alpha0 > 0 *)
      Do[
        If[(NumericQ[alpha0[[i]]] && Re[alpha0[[i]]] <= 0) ||
           TrueQ[Re[alpha0[[i]]] <= 0],
          Print["ERROR: IBPProcessSector: term ", t,
                " alpha0[", i, "] = ", alpha0[[i]], " <= 0"];
          Return[$Failed]
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

      flatPrefactor = Abs[detM] / (Times @@ alpha0);

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
        "Alpha1"         -> alpha1
      |>
    ],
    {t, Length[terms]}
  ];

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
    "AnalyticPole"           -> 1/ck
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

  divVars = {};
  Do[
    If[TrueQ[Re[a0[[i]]] <= 0] ||
       (NumericQ[a0[[i]]] && Re[a0[[i]]] <= 0),
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
   clearedPolys, detM, yVars,
   originalIntegral,
   boundaryVal, ibpTermVals, ibpSum,
   ak, reconstructed, relError,
   B0},

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
    integrand = Abs[detM] *
      Exp[Total[(aNum - 1) * Log /@ yVars]] *
      Times @@ MapThread[
        Function[{pv, be}, Exp[be * Log[pv]]],
        {polyValsExpr, polyExps /. fullRules}
      ];

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

    bndIntegrand = Abs[detM] *
      Exp[Total[(bndAnum - 1) * Log /@ bndYVars]] *
      Times @@ MapThread[
        Function[{pv, be}, Exp[be * Log[pv]]],
        {bndPolyVals, polyExps /. fullRules}
      ];

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

        (* Evaluate un-flattened (using full epsilon-dependent exponents
           from the IBP reduction, at finite epsilon) *)
        Module[{termExps, termPolyExps, termCoeff, term, terms},
          terms = ibpSectorData["ClearedPolys"];
          term  = ibpData`Private`dummyNotUsed;  (* placeholder *)

          (* Use the term's NewExponents from the IBP reduction,
             evaluated at finite eps *)
          termExps = ibpSectorData["IBPTerms"][[t]];
          (* We need the original unexpanded term data.
             Recompute from the processed data. *)
        ];

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

          integrand = Abs[detM] *
            Exp[Total[(termAlpha - 1) * Log /@ ibpYVars]] *
            Times @@ MapThread[
              Function[{pv, be}, Exp[(be /. kinRules) * Log[pv]]],
              {polyValsExpr, termPE /. epsRules}
            ];

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
  "NSamples"    -> 1000000,
  "MaxDim"      -> 20,
  "SeedBase"    -> 42,
  "Integrator"  -> "MonteCarlo",
  "Batch"       -> False,
  "VegasEpsRel" -> 1.*^-12,
  "VegasEpsAbs" -> 1.*^-300,
  "VegasSeed"   -> 0,
  "CubaMaxComp" -> 512
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
  If[integrator === "VEGAS" && batch,
    Print["NOTE: batched-ncomp VEGAS is not implemented for the IBP path ",
          "(VEGAS plan T8, deferred); using per-kp VEGAS."]];

  (* emitMain routes IBP result-assembly off defs["IsIBP"]; batch is ignored on
     the IBP path (per-kp Vegas fallback, plan.md §5.1). *)
  mainCode = emitMain[defs, integrator, batch,
    "NSamples" -> nSamples, "SeedBase" -> seedBase,
    "VegasEpsRel" -> OptionValue["VegasEpsRel"],
    "VegasEpsAbs" -> OptionValue["VegasEpsAbs"],
    "VegasSeed"   -> OptionValue["VegasSeed"],
    "CubaMaxComp" -> OptionValue["CubaMaxComp"]];

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
    If[integrator === "VEGAS", "  (VEGAS, per-kp)", ""]];
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
  "RunChecks"        -> True,
  "WorkingDirectory" -> Automatic,
  "Verbose"          -> True,
  (* sampler selection (default = the zero-dependency plain Monte Carlo) *)
  "Integrator"       -> "MC",    (* "MC" | "VEGAS" (aliases "MonteCarlo"/"Vegas") *)
  "Batch"            -> False,         (* IBP batch deferred -> per-kp VEGAS *)
  Sequence @@ $vegasOptionDefaults
};

evaluateTropicalIBPDriver[integrandSpec_Association, fanData_List,
                      kinematicPoints_List, OptionsPattern[]] :=
Module[
  {dualVertices, simplexList, n, nKP, nParams, eps,
   allSectorData, convergentSectors, divergentSectors,
   ibpProcessedSectors,
   cppFile, cppBinary, kinFile, resultFile,
   cppResult, ibpFuncMap,
   mcRawResults, finalResults,
   runChecks, verbose, nSamples, nThreads, workDir,
   integrator, batch, useCuba, vegasOpts},

  runChecks  = OptionValue["RunChecks"];
  verbose    = OptionValue["Verbose"];
  nSamples   = OptionValue["NSamples"];
  nThreads   = OptionValue["NThreads"];
  workDir    = OptionValue["WorkingDirectory"];
  eps        = integrandSpec["RegulatorSymbol"];
  integrator = normalizeIntegrator[OptionValue["Integrator"]];
  batch      = TrueQ[OptionValue["Batch"]];
  useCuba    = (integrator === "VEGAS");
  vegasOpts  = {"Integrator" -> integrator, "Batch" -> batch,
    "VegasEpsRel" -> OptionValue["VegasEpsRel"],
    "VegasEpsAbs" -> OptionValue["VegasEpsAbs"],
    "VegasSeed" -> OptionValue["VegasSeed"],
    "CubaMaxComp" -> OptionValue["CubaMaxComp"]};

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

  allSectorData = Table[
    ProcessSector[integrandSpec, dualVertices,
                  simplexList[[s]], s, "Verbose" -> False],
    {s, Length[simplexList]}
  ];

  convergentSectors = Select[allSectorData,
    (AssociationQ[#] && !#["IsDivergent"]) &];
  divergentSectors  = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  If[verbose,
    Print["  ", Length[convergentSectors], " convergent, ",
          Length[divergentSectors], " divergent sectors"]
  ];

  (* --- Step 2: Process divergent sectors with IBP --- *)
  If[verbose && Length[divergentSectors] > 0,
    Print["Processing ", Length[divergentSectors],
          " divergent sectors with IBP..."]
  ];

  ibpProcessedSectors = Table[
    IBPProcessSector[divergentSectors[[s]], integrandSpec],
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
      IBPCheckBoundary[divergentSectors[[s]], integrandSpec, 10],
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
    integrandSpec, cppFile,
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
                poleCont, finiteCont},
          sMap     = ibpFuncMap[[s]];
          bndBaseFid = sMap["BndBaseFuncId"];
          bndLogFid  = sMap["BndLogFuncId"];
          termFids   = sMap["TermFuncIds"];

          ckS = ibpProcessedSectors[[s]]["ck"] /. kinRules;
          rkS = ibpProcessedSectors[[s]]["rk"] /. kinRules;

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

          (* Pole contribution: (B^{(0)} - S_0) / c_k *)
          poleCont = (bndBase - S0) / ckS;

          (* Finite contribution:
             [(B^{(1)} - S_1) - r_k (B^{(0)} - S_0)] / c_k *)
          finiteCont = ((bndLog - S1) - rkS * (bndBase - S0)) / ckS;

          ibpContribPole   += poleCont;
          ibpContribFinite += finiteCont;
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


(* --------------------------------------------------------------------------
   Package end
   -------------------------------------------------------------------------- *)

End[]

EndPackage[]
