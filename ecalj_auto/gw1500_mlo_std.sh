#!/bin/bash
# The standard MLO recipe of GW1500 for one material (2026-10-05; user: one recipe that everyone can reproduce, as simple as
# possible). <src> holds PlotBand/ of gw1500_rerun.sh plus efermi.lmf of the QSGW run (as for gw1500_mlo.sh).
#
#   gw1500_mlo_std.sh <src> <out> [np]
#
#   <out>/b1     baseline 1: the [mlo] of gwinit (gw1500_mlo.sh)
#   <out>/b2     only when b1 is not PASS: + EH2 s,p of the cations (the mlo_lm2 rows of gwinit; MLO_LM2=1)
#   <out>/b2all  only when b1 is not PASS: + EH2 s,p of every atom but the transition metals, 4f and 5f (MLO_LM2=all)
#   <out>/best.txt   "<variant> <grade>": the best grade (PASS > OK > FAIL, gw1500_mlo_grade.py); on a tie the simpler one
# The model matrices are removed after each variant (the recipe makes them again). Env: ECALJ_BIN (default ~/bin).
src=$1; out=$2; np=${3:-4}; H=$(dirname "$(readlink -f "$0")")
mkdir -p "$out"
run(){  # variant, MLO_LM2
  rm -rf "$out/$1"; cp -r "$src" "$out/$1"
  MLO_LM2=$2 "$H/gw1500_mlo.sh" "$out/$1" $np > "$out/$1.line" 2>&1
  ( cd "$out/$1" && rm -rf __* HamRsMLO PROCAR.* rst.* sigm sigm.* save.* )
  g=$(cat "$out/$1/grade.json" 2>/dev/null | python3 -c "import json,sys; print(json.load(sys.stdin)['grade'])" 2>/dev/null)
  echo "${g:-NONE}"
}
g1=$(run b1 "")
best="b1 $g1"
if [ "$g1" != PASS ]; then
  g2=$(run b2 1); g3=$(run b2all all)
  for vg in "b2 $g2" "b2all $g3"; do
    set -- $vg $best
    r(){ case $1 in PASS) echo 0;; OK) echo 1;; FAIL) echo 2;; *) echo 3;; esac; }
    [ $(r $2) -lt $(r $4) ] && best="$1 $2"
  done
fi
echo "$best" > "$out/best.txt"
echo "$(basename "$out") $best"
