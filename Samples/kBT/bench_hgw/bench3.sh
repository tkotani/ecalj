#!/bin/bash
# hgw timing bench on the LiTi2O4 6^3 state (copy of liti_mlo_v9 after iteration 10).
#   bench3.sh <tag> <np> [VAR=value ...]  -> res_<tag>/ (lgw, SE* files, stdout.*.hgw), line in bench.log
#   np = number of MPI ranks = number of GPUs (one rank per GPU)
#   BENCH_BIN=<dir> selects the binaries (default ~/bin); BENCH_ARGS adds hgw options; BENCH_PRE wraps the launch;
#   BENCH_PREC the precision options (default --use_fp32; tf32 = "--use_fp32 --sigma_tf32");
#   GPUS the GPUs (default 0; kt1's GPU 1 is for the user's jobs, GPUS=0,1 only when allowed; the runs of 2026-09-27 used 0,1);
#   BENCH_DIR the state to run on (default the 6^3 bench_hgw666; bench_hgw999 for 9^3)
set -u
cd ${BENCH_DIR:-/mnt/data1/LiTi2O4_kbt_runs/bench_hgw666}
source env.sh
tag=$1; np=$2; shift 2
BIN=${BENCH_BIN:-$HOME/bin}
R=res_$tag; rm -rf $R; rm -f stdout.*.hgw; mkdir -p $R
t0=$(date +%s.%N)
env "$@" CUDA_VISIBLE_DEVICES=${GPUS:-0} ${BENCH_PRE:-} mpirun -np $np $BIN/${BENCH_EXE:-hgw_mp_gpu} liti2o4 --jobgw=1 ${BENCH_PREC:---use_fp32} --ntqxx --mlo ${BENCH_ARGS:-} > $R/lgw 2>&1
rc=$?
t1=$(date +%s.%N)
echo "$(date +%m-%d_%H:%M) $tag np=$np bin=$BIN env=[$*] args=[${BENCH_ARGS:-}] rc=$rc wall=$(echo "$t1-$t0"|bc)" | tee -a bench.log
cp SECU SEC2U SEXU SEX2U SEXcoreU XCU $R/ 2>/dev/null; mv stdout.*.hgw $R/ 2>/dev/null
