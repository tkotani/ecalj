#!/bin/bash
# Run the testecalj targets of Samples group by group and write one summary (2026-09-28; user: "Samples のテストをしっかり").
#
#   samples_tests.sh [--gpu] [--mp] [--run-args=ARGS] [-np N] [-np2 M] [group ...]      default: every group, in the order below
#   (--mp and --run-args go to testecalj: e.g. --gpu --mp is tf32, --gpu --run-args=--prec=fp32 is fp32; ForDevelopers section 5)
#
# Env: BIN (default ~/bin: the installed testecalj and binaries), LOG (default <tree>/samples_tests_<host>_<date>).
# Groups (the tracked dirs with a test.py; Samples/README.md "Verification status"):
#   install  TestInstall --all                      gwall    TestInstall --gwall (the GW targets only)
#   eps      EPS/*                                  procar   PROCAR/*
#   mlo      MLOsamples/*                           mloqsgw  MLOQSGW/*        afsym   AFsymmetry/*
#   samples  FermiSurface Doping HomoGas BoltzTraP SLAB SOC LDAU EffectiveMass Relax kBT/scanT: one summary line per directory
#            (the samples rebuilt from Legacy in 2026-09, and the temperature scan; a directory that is not there is skipped)
#   bench    BenchmarkTest/* (with --gpu: -np2 1, two GW ranks on one 32 GB GPU run out of memory)
#   heavy    TestInstall cugase2_gwsc222 nio_gwsc444 pdo_gwsc443 gas_gwsc666
#   magnon   Magnon/*  (work dirs of 4-5 GB)
#   inputs   every Samples/**/ctrlg.<sname>.toml outside the *_work dirs: lmchk reads it (the [struc] rules of the loader), and
#            the rules of ctrlg_update.py (pylib/ctrlg_rules.py: no nbas/nspec, t_tetrakbt present, no SmearX0) leave it unchanged
# A target is a subdirectory with test.py (the <target>_work copies that testecalj makes are skipped).
set -u
ROOT=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
BIN=${BIN:-$HOME/bin}
NP=8; NP2=""; GPU=""; MP=""; RUNARGS=""; GRP=()   # (not GROUPS: a bash variable that ignores assignments)
while [ $# -gt 0 ]; do
  case $1 in
    --gpu) GPU=--gpu ;;
    --mp) MP=--mp ;;
    --run-args=*) RUNARGS=${1#--run-args=} ;;
    -np) NP=$2; shift ;;
    -np2) NP2=$2; shift ;;
    *) GRP+=("$1") ;;
  esac; shift
done
[ ${#GRP[@]} = 0 ] && GRP=(inputs install eps procar mlo mloqsgw afsym samples bench heavy magnon)
# --gpu without -np2: one GW rank per visible GPU (2026-09-28: 8 GW ranks on the one 32 GB GPU of kr7 ran out of memory in the heavy group)
if [ -n "$GPU" ] && [ -z "$NP2" ]; then
  if [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then NP2=$(echo $CUDA_VISIBLE_DEVICES | tr "," "\n" | grep -c .); else NP2=$(nvidia-smi -L 2>/dev/null | grep -c "^GPU"); fi
  [ "${NP2:-0}" -ge 1 ] 2>/dev/null || NP2=1
fi
LOG=${LOG:-$ROOT/samples_tests_$(hostname -s)_$(date +%Y%m%d_%H%M)}
mkdir -p $LOG; SUM=$LOG/summary.txt
echo "# samples_tests.sh $(date '+%F %T') on $(hostname -s): tree $ROOT ($(head -1 $ROOT/SRC/.ecalj_rev 2>/dev/null || git -C $ROOT rev-parse --short HEAD 2>/dev/null)), BIN=$BIN, -np $NP ${NP2:+-np2 $NP2} $GPU $MP ${RUNARGS:+--run-args=$RUNARGS}" | tee $SUM

targets(){ # <dir>: the subdirectories with test.py
  ( cd $ROOT/Samples/$1 && for d in */; do d=${d%/}; [ -f $d/test.py ] && [[ $d != *_work ]] && echo $d; done )
}
run(){ # <group> <dir under Samples> <testecalj arguments...>
  local g=$1 d=$2; shift 2
  local t0=$(date +%s) st n np nf
  echo "=== $g: Samples/$d: $*" >> $SUM
  ( cd $ROOT/Samples/$d && $BIN/testecalj "$@" -np $NP $GPU $MP ${NP2:+-np2 $NP2} ${RUNARGS:+--run-args=$RUNARGS} ) > $LOG/$g.log 2>&1
  # testecalj prints the summary of all targets so far, with a status line, after each target: read only the last one
  # (2026-09-28: counting every PASSED line of the log gave 832 for the 64 checks of TestInstall --all, and a failure in a
  # later target was hidden behind the "OK! ALL PASSED" of the earlier ones)
  st=$(grep -E "OK! ALL PASSED|FAILED at some tests" $LOG/$g.log | tail -1)
  case "$st" in "OK! ALL PASSED"*) st=PASSED ;; "FAILED at some tests"*) st=FAILED ;; *) st=STOPPED ;; esac
  # a target whose test stopped (runprogs "Error exit!") leaves a START without an END and no summary of its own, so the
  # last status line is that of the target before it (2026-09-28: PROCAR/Ni2MnGa on mic counted as passed)
  local ns ne; ns=$(grep -c "^=== START test at" $LOG/$g.log); ne=$(grep -c "^=== END test at" $LOG/$g.log)
  [ $ne -lt $ns ] && st=STOPPED
  n=$(grep -n "^=== END test at" $LOG/$g.log | tail -1 | cut -d: -f1)
  np=$(tail -n +${n:-1} $LOG/$g.log | grep -c "^PASSED!")
  nf=$(tail -n +${n:-1} $LOG/$g.log | grep -E "^FAILED" | grep -v -c "^FAILED at some tests")
  printf "%-8s %-7s %4d checks passed, %3d failed  %6d s\n" $g $st $np $nf $(( $(date +%s)-t0 )) | tee -a $SUM
  [ $nf -gt 0 ] && tail -n +${n:-1} $LOG/$g.log | grep -E "^FAILED" | grep -v "^FAILED at some tests" | cut -c1-200 | sed 's/^/         /' | tee -a $SUM
  [ $st = STOPPED ] && echo "         stopped in $(grep "^=== START test at" $LOG/$g.log | tail -1 | sed 's|.*/||') (no summary of its own): see $g.log" | tee -a $SUM
  return 0
}
run_inputs(){ # the inputs group (see the header)
  local t0=$(date +%s) np=0 nf=0 f d s w rc rule
  echo "=== inputs: every ctrlg.<sname>.toml of Samples: lmchk reads it, and the rules of ctrlg_update.py leave it unchanged" >> $SUM
  while read f; do
    d=$(dirname $f); s=$(basename $f .toml); s=${s#ctrlg.}
    w=$LOG/inputs_work/$(echo ${d#$ROOT/Samples/} | tr / _); rm -rf $w; mkdir -p $w; cp $f $w/
    ( cd $w && mpirun -np 1 $BIN/lmchk $s > llmchk 2>&1 < /dev/null ); rc=$?   # (mpirun reads stdin: it ate the file list)
    rule=$(python3 -c "import sys; sys.path.insert(0, sys.argv[1]); from pylib.ctrlg_rules import update
t = open(sys.argv[2]).read(); n, m = update(t, sys.argv[2]); print('same' if n == t else m)" $ROOT/SRC/exec $f 2>&1 | tail -1)
    if [ $rc = 0 ] && [ "$rule" = same ]; then np=$((np+1))
    else nf=$((nf+1)); echo "         FAILED ${f#$ROOT/}: lmchk rc=$rc ($(grep -m1 -i -E "error|abort|stop|rx:" $w/llmchk | cut -c1-80)); rules: $rule" | tee -a $SUM; fi
  done < <(find $ROOT/Samples -name 'ctrlg.*.toml' -not -path '*_work*' | sort)
  printf "%-8s %-7s %4d checks passed, %3d failed  %6d s\n" inputs $([ $nf = 0 ] && echo PASSED || echo FAILED) $np $nf $(( $(date +%s)-t0 )) | tee -a $SUM
}
for g in "${GRP[@]}"; do
  case $g in
    inputs)  run_inputs ;;
    install) run install TestInstall --all ;;
    gwall)   run gwall TestInstall --gwall ;;
    eps)     run eps EPS $(targets EPS) ;;
    procar)  run procar PROCAR $(targets PROCAR) ;;
    mlo)     run mlo MLOsamples $(targets MLOsamples) ;;
    mloqsgw) run mloqsgw MLOQSGW $(targets MLOQSGW) ;;
    afsym)   run afsym AFsymmetry $(targets AFsymmetry) ;;
    samples) for d in FermiSurface Doping HomoGas BoltzTraP SLAB SOC LDAU EffectiveMass Relax kBT/scanT; do
               [ -d $ROOT/Samples/$d ] && [ -n "$(targets $d)" ] && run ${d//\//_} $d $(targets $d)
             done ;;
    bench)   if [ -n "$GPU" ]; then NP2=1 run bench BenchmarkTest $(targets BenchmarkTest)   # one GW rank even with two GPUs (README)
             else run bench BenchmarkTest $(targets BenchmarkTest); fi ;;  # CPU: -np ranks (2026-09-30: with -np2 1 the GW
                                                                         # programs ran on one core, on mic over 26 hours)
    heavy)   run heavy TestInstall cugase2_gwsc222 nio_gwsc444 pdo_gwsc443 gas_gwsc666 ;;
    magnon)  run magnon Magnon $(targets Magnon) ;;
    *) echo "unknown group $g" | tee -a $SUM ;;
  esac
done
echo "# done $(date '+%F %T')" | tee -a $SUM
