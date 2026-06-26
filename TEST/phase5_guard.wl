(* Guard arithmetic unit test (no fan build): confirms the strengthened
   threshold fires on the previously-silent 8D 4e6 case and not on 4e7,
   and is unchanged at low d. Mirrors the in-driver guardFactor logic. *)
SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
Get[FileNameJoin[{Directory[], "tropical_eval.wl"}]];

(* replicate the driver guard predicate *)
fires[nSamp_, maxSD_] := Module[{ns, gf, thr},
  ns = TropicalEval`Private`resolveVegasSizing[Automatic, maxSD, 1000, "nstart"];
  gf = If[maxSD >= 5, 100, 20];
  thr = gf*ns;
  {ns, thr, nSamp < thr}];

Print["d=8 NStart=", fires[4000000, 8][[1]]];
Print["d=8 4e6 (prev-silent) threshold=", fires[4000000, 8][[2]],
      " FIRES=", fires[4000000, 8][[3]]];
Print["d=8 4e7 threshold=",  fires[40000000, 8][[2]],
      " FIRES=", fires[40000000, 8][[3]]];
Print["d=8 8.1e6 (=100*NStart) FIRES=", fires[8100001, 8][[3]]];
Print["d=4 1e6 (low-d, 20x) FIRES=", fires[1000000, 4][[3]],
      " (NStart=", fires[1000000, 4][[1]], ", thr=", fires[1000000, 4][[2]], ")"];
Print["d=2 1e6 (low-d) FIRES=", fires[1000000, 2][[3]]];

(* The decisive checks for the gate: *)
g1 = fires[4000000, 8][[3]] === True;     (* prev-silent case MUST fire *)
g2 = fires[40000000, 8][[3]] === False;   (* sufficient budget passes *)
g3 = fires[1000000, 2][[3]] === False;    (* low-d unchanged *)
g4 = fires[1000000, 4][[2]] === 20*1000;  (* low-d uses 20x of base 1000 *)
Print["GUARD_PASS=", g1 && g2 && g3 && g4];
