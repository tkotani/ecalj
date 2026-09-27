#!/bin/bash
# Continue a run_gwsc10.sh run (2026-09-28): iterations FROM..TO, one `gwsc 1` each, and after each THE MLO band
# (draw_mloband.sh) along Gamma-X, so that the bands can be followed iteration by iteration (SIGM_BAND=1: the conventional
# sigm band by job_band too).  gwsc 1 in the run directory continues from the last QPU.<N>run.
#   usage: cont_gwsc.sh <tag> <from> <to> <bindir>
#   PREC (tf32 | fp32 | fp64, default fp32), GPUS (default 0; GPUS=0,1 for two), RUNS_DIR as in run_gwsc10.sh
# Output in <RUNS_DIR>/<tag>: bnd_mlo_iter<N>.dat, bndPMT_iter<N>.dat, snap/iter<N>/ (the state after iteration N),
# steps.log (one line per iteration, then the band lines).
set -u
TAG=$1; FROM=$2; TO=$3; BIN=$4
S=${RUNS_DIR:-/mnt/data1/LiTi2O4_kbt_runs}; T=liti2o4; D=$S/$TAG
HERE=$(dirname "$(readlink -f "$0")")
DRAW=$S/draw_mloband_qmlo.sh; [ -f $DRAW ] || DRAW=$HERE/draw_mloband.sh
QPATH=${QPATH:-$S/qplist_path_GX211.dat}    # the 211 points of Gamma-X (see draw_iter_bands.sh)
GPUS=${GPUS:-0}; NP2=$(echo $GPUS | tr ',' '\n' | wc -l)
cd $D || exit 1
source env.sh
export ECALJ_MLO_MIX=1
L=$D/steps.log
say(){ echo "$(date '+%m-%d %H:%M') $*" >> $L; }
B=$(dirname "$(readlink -f $BIN/lmf)")
say "continue $FROM..$TO bin=$BIN rev=$(head -1 $B/../.ecalj_rev 2>/dev/null || head -1 $B/FROZEN_REV 2>/dev/null) prec=${PREC:-fp32} GPUS=$GPUS"
for it in $(seq $FROM $TO); do
  t0=$(date +%s)
  CUDA_VISIBLE_DEVICES=$GPUS $BIN/gwsc 1 -np 60 -np2 $NP2 --gpu --prec=${PREC:-fp32} --ntqxx --mlo $T > gwsc_it$it.log 2>&1; rc=$?
  say "iter $it rc=$rc secs=$(( $(date +%s)-t0 )) mloON=$(grep -c 'MLO Sigma interpolation ON' llmf 2>/dev/null) $(grep ehf llmf 2>/dev/null | tail -1 | tr -s ' ')"
  [ $rc -eq 0 ] || { say ABORT; exit 1; }
  SD=$D/snap/iter$it; mkdir -p $SD
  for f in rst.$T sigm sigm.$T QMLO_SigRs QMLO_z HamRsMLO ctrlg.$T.toml env.sh __atm.$T efermi.lmf __mixsig __QMLO_mixsig; do
    cp -f $f $SD/ 2>/dev/null
  done
  printf "211 0.0 0.0 0.0 0.0 -1.0 0.0 GAMMA X\n0 !terminator\n" > $SD/syml.$T    # the path of run_gwsc10.sh
  if [ "${SIGM_BAND:-0}" = 1 ] || [ ! -s $QPATH ]; then
    ( cd $SD && source env.sh && CUDA_VISIBLE_DEVICES= $BIN/job_band $T -np 8 --nognuplot > lb_path.log 2>&1 )
    [ -s $SD/bnd001.spin1 ] && cp $SD/bnd001.spin1 $D/bndPMT_iter$it.dat
    sb="sigm $( [ -s $SD/bnd001.spin1 ] && echo ok || echo FAILED )"
  else    # the MLO band only (2026-09-28, user: only the MLO band is compared); qplist.dat as in draw_iter_bands.sh
    python3 -c 'import sys; ef=float(open(sys.argv[1]).read().split()[0].replace("D","E")); print(f"{ef:16.10f} ! ef or estaticav (efermi.lmf)"); sys.stdout.write(open(sys.argv[2]).read())' \
      $SD/efermi.lmf $QPATH > $SD/qplist.dat
    sb="sigm not drawn"
  fi
  mb=$(MLOBAND_BUILD=$B bash $DRAW $SD $D/mlobandwork_it$it $D/bnd_mlo_iter$it.dat 2>&1 | tail -1)
  say "   bands after iter $it: $sb; MLO: $mb"
done
say "continue done"; touch $S/$TAG.cont$TO.done
