#!/bin/bash
# 2026-09-30/10-01: GW1500, the 1389 materials that May left GOOD (1210) or NOTCONV (179), from scratch in fp32 with the
# present code (user 2026-09-30 23:0x: "the TF32 of May may be suspect; check with the new one"). Same settings as the rerun
# of the 147 (gw1500_rerun.sh, MAXITER 10, TOL 0.1); gwscconv of ecfac2c6b (judges a metal by its eigenvalues).
# Fortran binaries of b81da2342, as for the 147. Queue (user: problem materials first): 1 suspect GOOD (gap below LDA,
# oscillation over 0.5 eV after the first step) 18, 2 NOTCONV 179, 3 GOOD drifting by 0.08 eV or more in the last 3, 63,
# 4 the other GOOD 1128; fewest atoms first within each. Column 5 is the gap of May (TF32) for comparison.
# Cores 64-111 (run2m and NaCoO2 use 16-63); ranks unbound. About 6 days with 4 workers (estimate from the 147).
R=/mnt/data1/gw1500_rerun
export RUN_DIR=$R/run3 BIN=$HOME/bin_frozen_b81da2342m TOOLS_BIN=$HOME/bin_dev POSCAR_DIR=$HOME/ecalj/ecalj_auto/INPUT/gw1500/POSCARALL
export ENVSH=/mnt/data1/LiTi2O4_kbt_runs/liti_src_full9/env.sh NP=12 GPUS=0,1 PREC=fp32 LIMIT=28800 MAXITER=10
export OMPI_MCA_hwloc_base_binding_policy=none
mkdir -p $R/run3
[ -e $R/run3_queue.txt ] || cp $R/run3_queue.orig.txt $R/run3_queue.txt
i=0
for w in R1 R2 R3 R4; do
  a=$((64+12*i)); b=$((a+11)); i=$((i+1))
  nohup setsid taskset -c $a-$b bash $R/gw1500_rerun_v3.sh $R/run3_queue.txt $w > $R/run3/worker_$w.out 2>&1 < /dev/null &
done
disown -a
