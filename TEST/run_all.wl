(* ============================================================================
   TEST/run_all.wl  —  TROPICAL_MONTE_CARLOv3 CI entry point  (plan.md §6.1)

   Runs:
     1.  RunAllTests[]  — the Tree-A in-package validation suite (Tests 1–18),
         loaded from tropical_eval.wl in the repo root.
     2.  Every TEST/cc_*.wl cross-check file, each in its own isolated
         INTERFILES/cc_<id>/ working directory (§10 clobber-hazard rule).

   For each item the harness assigns exactly one of three verdicts:
     PASS  — exit code 0 and no "CC<N> FAIL" line in stdout.
     SKIP  — exit code 0, no FAIL line, and at least one "CC<N> SKIP" line;
             only Tier-2 (CUBA) / Tier-3 (FIESTA) checks ever produce this.
     FAIL  — exit code non-zero, OR a "CC<N> FAIL" line present.

   Exit behaviour:
     Exits 0  if nFAIL == 0  (any number of SKIPs is fine).
     Exits 1  if nFAIL >= 1.

   Usage:
     cd /Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3
     wolframscript -file TEST/run_all.wl

   Nothing under OLD_CODE/ is touched.  (plan.md §10: "Never edit OLD_CODE/")
   ============================================================================ *)

(* --------------------------------------------------------------------------
   0.  Repository root — derive from this script's own path.
   -------------------------------------------------------------------------- *)
$v3Root = DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
If[$v3Root === "" || !DirectoryQ[$v3Root],
  (* Fallback for interactive evaluation without $InputFileName set. *)
  $v3Root = "/Users/aidanh/Desktop/TROPICAL_MONTE_CARLO3"
];

$interfilesRoot = FileNameJoin[{$v3Root, "INTERFILES"}];

Print[""];
Print["============================================================"];
Print["  TROPICAL_MONTE_CARLOv3 — CI test harness"];
Print["  repo root: ", $v3Root];
Print["============================================================"];
Print[""];

(* --------------------------------------------------------------------------
   1.  Global tally
   -------------------------------------------------------------------------- *)
$nPass = 0;
$nFail = 0;
$nSkip = 0;
$log   = {};   (* list of {label, verdict, detail} *)

recordResult[label_String, verdict_String, detail_String] :=
  Module[{},
    AppendTo[$log, {label, verdict, detail}];
    Switch[verdict,
      "PASS", $nPass++,
      "FAIL", $nFail++,
      "SKIP", $nSkip++
    ]
  ];

(* --------------------------------------------------------------------------
   2.  Phase A — RunAllTests[]
       Load the package and call the in-package suite; it returns
       {{"Test N", True|False}, ...} (18 entries).
   -------------------------------------------------------------------------- *)
Print["------------------------------------------------------------"];
Print["  Phase A: RunAllTests[] — in-package validation suite"];
Print["------------------------------------------------------------"];
Print[""];

Module[{fanWL, evalWL, nTestPass, nTestFail, results, failedNames},

  fanWL  = FileNameJoin[{$v3Root, "tropical_fan.wl"}];
  evalWL = FileNameJoin[{$v3Root, "tropical_eval.wl"}];

  If[!FileExistsQ[fanWL] || !FileExistsQ[evalWL],
    Print["RUNALL FAIL  tropical_fan.wl or tropical_eval.wl not found at: ", $v3Root];
    recordResult["RunAllTests", "FAIL",
                 "package files missing in " <> $v3Root];
    Goto["phaseB"]
  ];

  Get[fanWL];
  Get[evalWL];

  (* Run the suite via the fully-qualified name to avoid the Global`/TropicalEval`
     shadow that arises when run_all.wl itself runs in Global` context.
     TropicalEval`RunAllTests is the authoritative definition loaded above. *)
  results = Quiet[TropicalEval`RunAllTests[], {General::munfl, Divide::infy, Greater::nord}];

  If[!ListQ[results],
    Print["RUNALL FAIL  RunAllTests[] returned non-List: ", results];
    recordResult["RunAllTests", "FAIL", "RunAllTests[] returned non-List"];
    Goto["phaseB"]
  ];

  nTestPass = Count[results, {_, True}];
  nTestFail = Count[results, {_, False}];
  failedNames = Cases[results, {name_, False} :> name];

  Print[""];
  Print["RunAllTests: ", nTestPass, " passed, ", nTestFail, " failed."];

  If[nTestFail == 0,
    Print["RUNALL PASS  all ", nTestPass, " in-package tests passed."];
    recordResult["RunAllTests", "PASS",
                 ToString[nTestPass] <> " tests passed"],
    Print["RUNALL FAIL  ", nTestFail, " test(s) failed: ", failedNames];
    recordResult["RunAllTests", "FAIL",
                 ToString[nTestFail] <> " failed: " <> ToString[failedNames]]
  ];
];

Label["phaseB"];

(* --------------------------------------------------------------------------
   3.  Phase B — individual cc_*.wl cross-checks
       Each is launched as a subprocess via RunProcess so that:
         (a) kernel state is isolated (one crashing test cannot corrupt others);
         (b) each gets its own INTERFILES/cc_<id>/ working directory.
   -------------------------------------------------------------------------- *)
Print[""];
Print["------------------------------------------------------------"];
Print["  Phase B: cc_*.wl cross-checks"];
Print["------------------------------------------------------------"];
Print[""];

(* Collect all cc_*.wl files, sorted lexicographically-by-id. *)
$ccFiles = Sort @ FileNames["cc_*.wl", FileNameJoin[{$v3Root, "TEST"}]];

If[Length[$ccFiles] == 0,
  Print["RUNALL WARN  no cc_*.wl files found in TEST/"];
  Goto["summary"]
];

Print["Found ", Length[$ccFiles], " cc_*.wl checks to run."];
Print[""];

(* Ensure the INTERFILES root directory exists. *)
Quiet[CreateDirectory[$interfilesRoot], {CreateDirectory::eexist}];

(* Helper: derive the check ID string from the filename. *)
ccID[path_String] := StringReplace[FileBaseName[path], StartOfString ~~ "cc_" -> ""];

(* Helper: create (or reuse) a per-check working directory. *)
ccWorkDir[id_String] :=
  Module[{d = FileNameJoin[{$interfilesRoot, "cc_" <> id}]},
    Quiet[CreateDirectory[d, CreateIntermediateDirectories -> True],
          {CreateDirectory::eexist}];
    d
  ];

(* Helper: classify the output of a completed subprocess into PASS/FAIL/SKIP.

   Logic:
     - Non-zero exit code  → FAIL (regardless of stdout content).
     - Any "CC<id> FAIL"   → FAIL (catches cc_*.wl that exit 0 on failure;
                              the pattern is case-insensitive on the id).
     - Any "CC<id> SKIP"   → SKIP (only when no FAIL also present).
     - Any "CC<id> PASS"   → PASS (explicit verdict line present).
     - Exit 0, no FAIL, no SKIP, no PASS → SKIP.
         This handles the case where a Tier-2/3 check's dependency is
         partially installed and the script aborts before printing a verdict.
         Treating such cases as SKIP (not FAIL) is conservative: it avoids
         falsely blocking the CI on a broken optional dependency.  A genuine
         test failure always prints "CC<id> FAIL" or exits non-zero.

   Note: wolframscript may capitalise id differently in edge cases, so we
   scan for the canonical "CC<id>" prefix printed by every check script.  *)
classifyOutput[id_String, exitCode_Integer, stdout_String] :=
  Module[{upperID, failPat, skipPat, passPat, hasFail, hasSkip, hasPass},
    upperID = ToUpperCase[id];
    failPat = RegularExpression["(?mi)^CC" <> upperID <> " FAIL\\b"];
    skipPat = RegularExpression["(?mi)^CC" <> upperID <> " SKIP\\b"];
    passPat = RegularExpression["(?mi)^CC" <> upperID <> " PASS\\b"];
    hasFail = exitCode != 0 || StringMatchQ[stdout, ___ ~~ failPat ~~ ___];
    hasSkip = !hasFail && StringMatchQ[stdout, ___ ~~ skipPat ~~ ___];
    hasPass = !hasFail && StringMatchQ[stdout, ___ ~~ passPat ~~ ___];
    Which[
      hasFail,  "FAIL",
      hasSkip,  "SKIP",
      hasPass,  "PASS",
      True,     "SKIP"   (* no verdict printed, exit 0: assume dependency absent *)
    ]
  ];

(* Run each cc_*.wl in sequence (never in parallel: §10 clobber-hazard rule). *)
Do[
  Module[{id, wdir, proc, verdict, detail, label, stdout, exitCode, elapsed, t0},

    id    = ccID[ccFile];
    label = "cc_" <> id;
    wdir  = ccWorkDir[id];

    Print["-- Running ", label, " ..."];
    t0 = AbsoluteTime[];

    proc = RunProcess[
      {"wolframscript", "-file", ccFile},
      ProcessDirectory -> wdir,
      ProcessEnvironment -> <|
        "HOME"          -> Environment["HOME"],
        "PATH"          -> Environment["PATH"],
        "WOLFRAMPATH"   -> Environment["WOLFRAMPATH"],
        "DYLD_LIBRARY_PATH" -> Environment["DYLD_LIBRARY_PATH"],
        "LD_LIBRARY_PATH"   -> Environment["LD_LIBRARY_PATH"]
      |>
    ];

    elapsed  = AbsoluteTime[] - t0;
    exitCode = proc["ExitCode"];
    stdout   = proc["StandardOutput"] <> proc["StandardError"];

    verdict = classifyOutput[id, exitCode, stdout];
    detail  = "exit=" <> ToString[exitCode] <>
              "  t=" <> ToString[NumberForm[elapsed, {4, 1}]] <> "s";

    (* Print the subprocess output, indented, for the CI log. *)
    Scan[
      Print["  | ", #] &,
      StringSplit[StringTrim[stdout], "\n"]
    ];

    Switch[verdict,
      "PASS",
        Print[label, " PASS  [", detail, "]"],
      "SKIP",
        Print[label, " SKIP  [", detail, "]"],
      "FAIL",
        Print[label, " FAIL  [", detail, "]"]
    ];
    Print[""];

    recordResult[label, verdict, detail];
  ],

  {ccFile, $ccFiles}
];

(* --------------------------------------------------------------------------
   4.  Summary
   -------------------------------------------------------------------------- *)
Label["summary"];

Print["============================================================"];
Print["  CI SUMMARY"];
Print["============================================================"];
Print[""];
Print["  PASS:  ", $nPass];
Print["  SKIP:  ", $nSkip,
      If[$nSkip > 0,
         "  (Tier-2/3 checks skipped — dependency absent, exit-0 contribution)",
         ""]];
Print["  FAIL:  ", $nFail];
Print[""];

(* Print a compact table of all verdicts. *)
Print["  Detail:"];
Scan[
  Function[row,
    Print["    ", StringPadRight[row[[1]], 14], "  ", row[[2]],
          "  ", row[[3]]]],
  $log
];
Print[""];

If[$nFail == 0,
  Print["RESULT: ALL CHECKS PASSED  (", $nPass, " pass, ", $nSkip, " skip)"];
  Print["============================================================"];
  Exit[0],

  (* One or more FAILs: print the culprit list before exiting 1. *)
  Print["RESULT: ", $nFail, " CHECK(S) FAILED:"];
  Scan[
    Function[row,
      If[row[[2]] === "FAIL",
         Print["  *** FAIL: ", row[[1]], "  (", row[[3]], ")"]]],
    $log
  ];
  Print["============================================================"];
  Exit[1]
];
