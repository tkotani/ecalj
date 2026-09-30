#!/bin/bash
# 2026-09-30 23:2x: start the workers of the GW1500 recheck (run3, kt1_recheck_0930.sh's settings) on the cores as they come free.
# kt1 has cores 0-63 only (SMT off); taskset 64- fails. Now: tests 0-15, run2m workers M1-M3 16-51, NaCoO2 52-63.
#   tests done         -> run3 worker R1 on 0-11
#   run2m worker Mi    -> run3 worker R(i+1) on its cores
#   NaCoO2 done        -> mp-1179832 (Rb8; lmf --jobgw=1 ran out of memory with 12 ranks beside 3 other workers) with 6 ranks
#                         on 52-63, into run2m; then run3 worker R5 there
R=/mnt/data1/gw1500_rerun
export BIN=$HOME/bin_frozen_b81da2342m TOOLS_BIN=$HOME/bin_dev POSCAR_DIR=$HOME/ecalj/ecalj_auto/INPUT/gw1500/POSCARALL
export ENVSH=/mnt/data1/LiTi2O4_kbt_runs/liti_src_full9/env.sh GPUS=0,1 PREC=fp32 LIMIT=28800
export OMPI_MCA_hwloc_base_binding_policy=none
start_run3(){  # $1 worker $2 cores
  RUN_DIR=$R/run3 NP=12 MAXITER=10 nohup setsid taskset -c $2 bash $R/gw1500_rerun_v3.sh $R/run3_queue.txt $1 \
    > $R/run3/worker_$1.out 2>&1 < /dev/null &
  echo "$(date '+%F %T') chain: started run3 worker $1 on cores $2" >> $R/run3/chain.log
}
mkdir -p $R/run3
( while [ ! -e /mnt/data1/ecalj_test0930d/tests_done_0930d.txt ]; do sleep 120; done; start_run3 R1 0-11 ) &
( while pgrep -f "rerun_v3.sh $R/run2m_queue.txt M1\$" > /dev/null; do sleep 120; done; start_run3 R2 16-27 ) &
( while pgrep -f "rerun_v3.sh $R/run2m_queue.txt M2\$" > /dev/null; do sleep 120; done; start_run3 R3 28-39 ) &
( while pgrep -f "rerun_v3.sh $R/run2m_queue.txt M3\$" > /dev/null; do sleep 120; done; start_run3 R4 40-51 ) &
( while pgrep -f "gwscconv .*mp-18921" > /dev/null; do sleep 120; done
  printf 'mp-1179832\tRb8\t8\n' > $R/run2m_rb8_queue.txt
  echo "$(date '+%F %T') chain: Rb8 with 6 ranks on 52-63" >> $R/run3/chain.log
  RUN_DIR=$R/run2m NP=6 MAXITER=15 taskset -c 52-63 bash $R/gw1500_rerun_v3.sh $R/run2m_rb8_queue.txt B1 > $R/run2m/worker_B1.out 2>&1 < /dev/null
  start_run3 R5 52-63 ) &
disown -a
