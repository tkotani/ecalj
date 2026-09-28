#!/bin/bash
# The LDA band of a run started from LDA (2026-09-28): job_band on the LDA rst that gwsc keeps in <run>/LDA, along the
# Gamma-X path of run_gwsc10.sh (211 points) -> <run>/bnd_lda.dat (all bands, as bench_hgw999/bnd_lda.dat).
# For the LDA row of the figures (update_six_rows.sh).  Runs with the same sites share one LDA band; a run with other
# sites (the empty spheres of qmlo_k9_tf32_es_i15) needs its own.  CPU only (job_band -np 8), a few minutes for 9^3.
#   usage: draw_lda_band.sh <tag> <bindir>        (RUNS_DIR as in run_gwsc10.sh)
set -u
TAG=$1; BIN=$2
S=${RUNS_DIR:-/mnt/data1/LiTi2O4_kbt_runs}; T=liti2o4; D=$S/$TAG; L=$D/LDA; W=$D/ldaband
[ -f $L/rst.$T ] || { echo "no LDA rst in $L"; exit 1; }
rm -rf $W; mkdir -p $W
cp -f $L/ctrlg.$T.toml $L/rst.$T $D/env.sh $D/__atm.$T $W/ || exit 1
cp -f $L/atmpnu.*.$T $W/ 2>/dev/null
printf "211 0.0 0.0 0.0 0.0 -1.0 0.0 GAMMA X\n0 !terminator\n" > $W/syml.$T    # the path of run_gwsc10.sh
( cd $W && source env.sh && CUDA_VISIBLE_DEVICES= $BIN/job_band $T -np 8 --nognuplot > lb_path.log 2>&1 )
if [ -s $W/bnd001.spin1 ]; then
  cp $W/bnd001.spin1 $D/bnd_lda.dat && rm -rf $W && echo "$D/bnd_lda.dat written $(date '+%m-%d %H:%M')"
else
  echo "job_band FAILED: see $W/lb_path.log"; exit 1
fi
