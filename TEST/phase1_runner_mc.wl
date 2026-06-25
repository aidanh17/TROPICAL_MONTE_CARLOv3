(* ============================================================================
   TEST/phase1_runner_mc.wl    --  Phase 1.1 reproduction: MC sampler
   Loads v3's tropical_eval.wl (verbatim copy of Tree A).
   Runs RunAllTests[] (Tests 1-18) + IBP Tests 19-22 + test_divergent_crosscheck.
   Uses INTERFILES/phase1_mc/ to avoid clobbering Tree A INTERFILES.
   ============================================================================ *)

v3Root = DirectoryName[$InputFileName] // FileNameJoin[{#, ".."}] &;
SetDirectory[v3Root];

Print["======================================================================="];
Print["PHASE 1 MC RUNNER -- started: ", DateString[]];
Print["v3 root: ", Directory[]];
Print["======================================================================="];
Print[""];

(* Step 1: load v3's package *)
Get["tropical_eval.wl"];
Print["v3 tropical_eval.wl loaded OK"];

(* Step 2: load test definitions from OLD_CODE EXAMPLES, but keep our
   SetDirectory on v3 root.  The examples file calls SetDirectory[parent of
   $InputFileName]; we restore ours immediately after. *)
exDir = FileNameJoin[{Directory[], "OLD_CODE", "TROPICAL_MONTE_CARLO",
                       "EXAMPLES"}];

(* The examples script SetDirectories to its parent (v3Root in old code).
   To make C++ output go into v3 INTERFILES rather than OLD_CODE INTERFILES,
   we temporarily create symlinks or (simpler) set an override variable.
   The examples file uses Directory[]+"/INTERFILES" verbatim, so we must
   set Directory[] to a location whose INTERFILES we control.
   Approach: keep Directory[] = v3Root; the examples file's SetDirectory
   call will also land on v3Root (since DirectoryName of exDir/../ = v3Root).
   Wait: DirectoryName[exDir <> "/tropical_eval_examples.wl"] = exDir,
   and exDir <> "/" <> ".." = OLD_CODE/TROPICAL_MONTE_CARLO.  So the examples
   file's SetDirectory would move us to OLD_CODE/TROPICAL_MONTE_CARLO.
   Fix: We read and evaluate the examples definitions in a way that avoids
   the problematic SetDirectory.  We strip that line by reading the file as
   text and prepending a SetDirectory to v3Root, then evaluating. *)

exFile = FileNameJoin[{exDir, "tropical_eval_examples.wl"}];
exText = Import[exFile, "Text"];

(* Remove the first two lines (SetDirectory + Get) of the examples file --
   they would redirect us to OLD_CODE and reload OLD_CODE's tropical_eval.wl.
   We have already loaded v3's copy. *)
exLines = StringSplit[exText, "\n"];
(* Find the line after "Get[..." which is the package load - skip those *)
startIdx = 1;
Do[
  If[StringMatchQ[exLines[[i]], "*SetDirectory*"] ||
     StringMatchQ[exLines[[i]], "*Get[*tropical_eval*"],
    startIdx = i + 1],
  {i, 1, Min[25, Length[exLines]]}
];
exStripped = StringJoin[Riffle[exLines[[startIdx ;;]], "\n"]];

(* Evaluate the definitions with v3 root as working directory *)
ToExpression[exStripped];
Print["OLD_CODE examples definitions loaded (RunAllTests etc. now defined)"];
Print["Current directory after loading: ", Directory[]];
SetDirectory[v3Root];  (* restore in case the examples file changed it *)

Print[""];
Print["======================================================================="];
Print["  RUNNING RunAllTests[] (Tests 1-18)"];
Print["======================================================================="];
Print[""];

allResults18 = RunAllTests[];

Print[""];
Print["RunAllTests[] results: ", allResults18];
Print[""];

(* ============================================================
   IBP Tests 19-22
   ============================================================ *)
Print["======================================================================="];
Print["  RUNNING IBP TESTS (19-22)"];
Print["======================================================================="];
Print[""];

ibpResults = {};
Do[
  Module[{fn, pass},
    fn = Symbol["RunTest" <> ToString[n]];
    Print["--- Test ", n, " ---"];
    pass = fn[];
    AppendTo[ibpResults, {"Test " <> ToString[n], pass}]
  ],
  {n, 19, 22}
];

Print[""];
Print["IBP test results: ", ibpResults];
Print[""];

(* ============================================================
   test_divergent_crosscheck: load just the definitions (no reload of eval.wl)
   ============================================================ *)
Print["======================================================================="];
Print["  RUNNING test_divergent_crosscheck.wl"];
Print["======================================================================="];
Print[""];

ccFile = FileNameJoin[{Directory[], "OLD_CODE", "TROPICAL_MONTE_CARLO",
                        "EXAMPLES", "test_divergent_crosscheck.wl"}];
ccText = Import[ccFile, "Text"];
ccLines = StringSplit[ccText, "\n"];
(* Strip the SetDirectory and Get lines (first ~3 lines) that reload eval.wl *)
ccStart = 1;
Do[
  If[StringMatchQ[ccLines[[i]], "*SetDirectory*"] ||
     StringMatchQ[ccLines[[i]], "*Get[*tropical_eval*"],
    ccStart = i + 1],
  {i, 1, Min[10, Length[ccLines]]}
];
(* Also strip the final RunDivergentCrossCheck[] auto-call -- we call it ourselves *)
ccEnd = Length[ccLines];
Do[
  If[StringMatchQ[ccLines[[i]], "If[!TrueQ*RunDivergentCrossCheck*"],
    ccEnd = i - 1],
  {i, Max[1, Length[ccLines]-5], Length[ccLines]}
];
ccStripped = StringJoin[Riffle[ccLines[[ccStart ;; ccEnd]], "\n"]];
ToExpression[ccStripped];
SetDirectory[v3Root];  (* restore after potential directory change *)

ccPass = RunDivergentCrossCheck[];
Print[""];
Print["crosscheck result: ", If[TrueQ[ccPass], "PASS", "FAIL"]];

Print[""];
Print["======================================================================="];
Print["PHASE 1 MC RUNNER COMPLETE -- ended: ", DateString[]];
Print["======================================================================="];
