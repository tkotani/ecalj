#!/bin/bash
# Draw the bands of iterations FROM..TO of a gwsc run afterwards, from what gwsc kept in QSGW.<N>run (2026-09-28): the
# conventional sigm band (job_band) and THE MLO band (draw_mloband.sh) along Gamma-X, the same two that cont_gwsc.sh
# draws after each iteration.  QSGW.<N>run holds the state at the end of iteration N (rst of the SCF with the new Sigma).
# The MLO band needs QMLO_SigRs, QMLO_z and efermi.lmf there, which gwsc keeps since 2026-09-28; for an older run only the
# sigm band is drawn.  HamRsMLO and __atm are taken from the run directory (--mlofreeze keeps one HamRsMLO for the chain).
# CPU only (job_band -np 8, draw_mloband.sh -np 8); 9^3 about 2 minutes per iteration.
# The MLO band only by default (user 2026-09-28: only the MLO band is compared); SIGM_BAND=1 draws the sigm band too.
# Without job_band, qplist.dat is made here: the E_F of efermi.lmf, then the path of QPATH (lines 2- of a qplist.dat that
# job_band wrote for the same path).  mlo takes the E_F of its energy window from that first line; in the gwsc loop it reads
# efermi.lmf itself, and job_band's E_F differs from it by 0.5 meV (9^3 iteration 10).
#   usage: draw_iter_bands.sh <tag> <from> <to> <bindir>        (RUNS_DIR as in run_gwsc10.sh)
# Output in <RUNS_DIR>/<tag>: bnd_mlo_iter<N>.dat (bndPMT_iter<N>.dat with SIGM_BAND=1); a line per iteration in steps.log.
set -u
TAG=$1; FROM=$2; TO=$3; BIN=$4
S=${RUNS_DIR:-/mnt/data1/LiTi2O4_kbt_runs}; T=liti2o4; D=$S/$TAG
HERE=$(dirname "$(readlink -f "$0")")
DRAW=$S/draw_mloband_qmlo.sh; [ -f $DRAW ] || DRAW=$HERE/draw_mloband.sh
QPATH=${QPATH:-$S/qplist_path_GX211.dat}    # the 211 points of Gamma-X of run_gwsc10.sh
B=$(dirname "$(readlink -f $BIN/lmf)")
cd $D || exit 1
say(){ echo "$(date '+%m-%d %H:%M') $*" >> $D/steps.log; }
for it in $(seq $FROM $TO); do
  R=$D/QSGW.${it}run; W=$D/iterbands/iter$it
  [ -f $R/rst.$T ] && [ -f $R/sigm.$T ] || { say "   bands of iter $it: no rst/sigm in $R"; continue; }
  rm -rf $W; mkdir -p $W
  cp -f $R/ctrlg.$T.toml $R/rst.$T $D/env.sh $D/__atm.$T $W/
  cp -f $R/sigm.$T $W/sigm; ln -s sigm $W/sigm.$T                 # as in the run directory: sigm.<target> -> sigm
  for f in efermi.lmf QMLO_SigRs QMLO_z; do [ -f $R/$f ] && cp -f $R/$f $W/; done   # copies: the kept files must not change
  printf "211 0.0 0.0 0.0 0.0 -1.0 0.0 GAMMA X\n0 !terminator\n" > $W/syml.$T    # the path of run_gwsc10.sh
  mlo_ok=0; [ -e $W/QMLO_SigRs ] && [ -e $W/QMLO_z ] && [ -e $W/efermi.lmf ] && mlo_ok=1
  if [ "${SIGM_BAND:-0}" = 1 ] || [ $mlo_ok = 0 ] || [ ! -s $QPATH ]; then
    ( cd $W && source env.sh && CUDA_VISIBLE_DEVICES= $BIN/job_band $T -np 8 --nognuplot > lb_path.log 2>&1 )
    [ -s $W/bnd001.spin1 ] && cp $W/bnd001.spin1 $D/bndPMT_iter$it.dat
    sb="sigm $( [ -s $W/bnd001.spin1 ] && echo ok || echo FAILED )"
  else
    python3 -c 'import sys; ef=float(open(sys.argv[1]).read().split()[0].replace("D","E")); print(f"{ef:16.10f} ! ef or estaticav (efermi.lmf)"); sys.stdout.write(open(sys.argv[2]).read())' \
      $W/efermi.lmf $QPATH > $W/qplist.dat
    sb="sigm not drawn"
  fi
  if [ $mlo_ok = 1 ]; then
    ln -s $D/HamRsMLO $W/HamRsMLO      # 0.9 GB for 9^3; only draw_mloband.sh reads it, from its own copy
    mb=$(MLOBAND_BUILD=$B bash $DRAW $W $W/mlobandwork $D/bnd_mlo_iter$it.dat 2>&1 | tail -1)
  else
    mb="not drawn: QSGW.${it}run has no QMLO_SigRs/QMLO_z/efermi.lmf (gwsc before 2026-09-28)"
  fi
  say "   bands of iter $it (from QSGW.${it}run): $sb; MLO: $mb"
  rm -rf $W
done
rmdir $D/iterbands 2>/dev/null
