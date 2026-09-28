#!/bin/bash
# Run the testecalj targets of Samples group by group and write one summary (2026-09-28; user: "Samples のテストをしっかり").
#
#   samples_tests.sh [--gpu] [-np N] [-np2 M] [group ...]      default: every group, in the order below
#
# Env: BIN (default ~/bin: the installed testecalj and binaries), LOG (default <tree>/samples_tests_<host>_<date>).
# Groups (the tracked dirs with a test.py; Samples/README.md "Verification status"):
#   install  TestInstall --all                      eps      EPS/*            procar  PROCAR/*
#   mlo      MLOsamples/*                           mloqsgw  MLOQSGW/*        afsym   Legacy/AFsymmetry/*
#   bench    BenchmarkTest/* (with -np2 1: two GW ranks on one 32 GB GPU run out of memory)
#   heavy    TestInstall cugase2_gwsc222 nio_gwsc444 pdo_gwsc443 gas_gwsc666
#   magnon   Legacy/Magnon/*  (work dirs of 4-5 GB)
# A target is a subdirectory with test.py (the <target>_work copies that testecalj makes are skipped).
set -u
ROOT=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
BIN=${BIN:-$HOME/bin}
NP=8; NP2=""; GPU=""; GRP=()   # (not GROUPS: a bash variable that ignores assignments)
while [ $# -gt 0 ]; do
  case $1 in
    --gpu) GPU=--gpu ;;
    -np) NP=$2; shift ;;
    -np2) NP2=$2; shift ;;
    *) GRP+=("$1") ;;
  esac; shift
done
[ ${#GRP[@]} = 0 ] && GRP=(install eps procar mlo mloqsgw afsym bench heavy magnon)
LOG=${LOG:-$ROOT/samples_tests_$(hostname -s)_$(date +%Y%m%d_%H%M)}
mkdir -p $LOG; SUM=$LOG/summary.txt
echo "# samples_tests.sh $(date '+%F %T') on $(hostname -s): tree $ROOT ($(head -1 $ROOT/SRC/.ecalj_rev 2>/dev/null || git -C $ROOT rev-parse --short HEAD 2>/dev/null)), BIN=$BIN, -np $NP ${NP2:+-np2 $NP2} $GPU" | tee $SUM

targets(){ # <dir>: the subdirectories with test.py
  ( cd $ROOT/Samples/$1 && for d in */; do d=${d%/}; [ -f $d/test.py ] && [[ $d != *_work ]] && echo $d; done )
}
run(){ # <group> <dir under Samples> <testecalj arguments...>
  local g=$1 d=$2; shift 2
  local t0=$(date +%s) st n np nf
  echo "=== $g: Samples/$d: $*" >> $SUM
  ( cd $ROOT/Samples/$d && $BIN/testecalj "$@" -np $NP $GPU ${NP2:+-np2 $NP2} ) > $LOG/$g.log 2>&1
  # testecalj prints the summary of all targets so far, with a status line, after each target: read only the last one
  # (2026-09-28: counting every PASSED line of the log gave 832 for the 64 checks of TestInstall --all, and a failure in a
  # later target was hidden behind the "OK! ALL PASSED" of the earlier ones)
  st=$(grep -E "OK! ALL PASSED|FAILED at some tests" $LOG/$g.log | tail -1)
  case "$st" in "OK! ALL PASSED"*) st=PASSED ;; "FAILED at some tests"*) st=FAILED ;; *) st=STOPPED ;; esac
  n=$(grep -n "^=== END test at" $LOG/$g.log | tail -1 | cut -d: -f1)
  np=$(tail -n +${n:-1} $LOG/$g.log | grep -c "^PASSED!")
  nf=$(tail -n +${n:-1} $LOG/$g.log | grep -E "^FAILED" | grep -v -c "^FAILED at some tests")
  printf "%-8s %-7s %4d checks passed, %3d failed  %6d s\n" $g $st $np $nf $(( $(date +%s)-t0 )) | tee -a $SUM
  [ $nf -gt 0 ] && tail -n +${n:-1} $LOG/$g.log | grep -E "^FAILED" | grep -v "^FAILED at some tests" | cut -c1-200 | sed 's/^/         /' | tee -a $SUM
  [ $st = STOPPED ] && echo "         no summary: see $g.log" | tee -a $SUM
  return 0
}
for g in "${GRP[@]}"; do
  case $g in
    install) run install TestInstall --all ;;
    eps)     run eps EPS $(targets EPS) ;;
    procar)  run procar PROCAR $(targets PROCAR) ;;
    mlo)     run mlo MLOsamples $(targets MLOsamples) ;;
    mloqsgw) run mloqsgw MLOQSGW $(targets MLOQSGW) ;;
    afsym)   run afsym Legacy/AFsymmetry $(targets Legacy/AFsymmetry) ;;
    bench)   NP2=1 run bench BenchmarkTest $(targets BenchmarkTest) ;;
    heavy)   run heavy TestInstall cugase2_gwsc222 nio_gwsc444 pdo_gwsc443 gas_gwsc666 ;;
    magnon)  run magnon Legacy/Magnon $(targets Legacy/Magnon) ;;
    *) echo "unknown group $g" | tee -a $SUM ;;
  esac
done
echo "# done $(date '+%F %T')" | tee -a $SUM
