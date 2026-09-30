#!/bin/bash
# 2026-10-01 mlocheck: one material, DFT -> band -> MLO -> plot.  usage: run_one.sh <Material> <np>
mat=$1; np=$2; B=$HOME/bin
cd $(dirname $0)/$mat || exit 1
s=$(ls ctrlg.*.toml | sed 's/^ctrlg\.//; s/\.toml$//')
export OMP_NUM_THREADS=1
st() { echo "$(date '+%F %T') $mat $*" >> ../status.log; }
st START np=$np
mpirun -np 1 $B/lmfa $s > llmfa 2>&1 || { st FAIL lmfa; exit 1; }
timeout 4h mpirun -np $np $B/lmf $s > llmf 2>&1 || { st FAIL lmf rc=$?; exit 1; }
st lmf done "$(grep -E '^[ ]*c |^x ' llmf | tail -1)"
python3 $B/getsyml $s --nobzview > lgetsyml 2>&1 || { st FAIL getsyml; exit 1; }
timeout 2h $B/job_band $s -np $np --nognuplot > ljob_band 2>&1 || { st FAIL job_band; exit 1; }
ls bnd001.spin1 >/dev/null 2>&1 || { st FAIL job_band nobnd; exit 1; }
timeout 4h $B/job_mlo $s -np $np --nognuplot > ljob_mlo 2>&1 || { st FAIL job_mlo; exit 1; }
ls band_MLO_spin1.dat >/dev/null 2>&1 || { st FAIL job_mlo noMLOband; exit 1; }
python3 $B/mlo_bandplot.py . -o mlo_$mat.png --label $mat > lplot 2>&1 || st WARN plot
st DONE
