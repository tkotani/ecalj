#!/bin/bash
# LiTi2O4 MLO-QSGW with the QMLO_* code (2026-09-27): ONE `gwsc 10` from LDA, to reproduce the chains
# that ran `gwsc 1` ten times with the old file names (liti_mlo_v9 = 6^3, liti_mlo_k9 = 9^3).
#   usage: run_gwsc10.sh <tag> <input dir> <bindir>
#   <input dir>: ctrlg.liti2o4.toml + env.sh (liti_src_full9 = 6^3, liti_src_full9_k9 = 9^3 on kt1; the ctrlg is
#                Samples/kBT/LiTi2O4/input/qmlo, where the README says how to make env.sh and the 9^3 input)
#   <bindir>   : where gwsc and the executables are (~/bin -> ~/ecalj, ~/bin_dev -> ~/ecalj_dev). For a run of
#                hours use a copy that a rebuild cannot change: cp -rL ~/bin_dev ~/bin_frozen_<rev>, and write the
#                revision into <copy>/FROZEN_REV (the executables find the copied .so through RUNPATH $ORIGIN)
#   RUNS_DIR   : where the run directory <tag> is made (default /mnt/data1/LiTi2O4_kbt_runs, kt1)
#   PREC       : precision of the GW programs, tf32 | fp32 (default) | fp64 (gwsc --prec)
#   GWSC_EXTRA : more gwsc options, e.g. --prec-final=fp32:2 (the last two iterations in fp32)
# gwsc keeps llmf.<N>run and QPU.<N>run of every iteration; the state after iteration 10 is copied to
# snap/final and its bands are drawn there: the conventional sigm band (job_band, also makes qplist.dat)
# and THE MLO band (draw_mloband_qmlo.sh = Samples/kBT/LiTi2O4/draw_mloband.sh).
set -u
TAG=$1; SRC=$2; BIN=$3
S=${RUNS_DIR:-/mnt/data1/LiTi2O4_kbt_runs}; T=liti2o4; D=$S/$TAG
HERE=$(dirname "$(readlink -f "$0")")
DRAW=$S/draw_mloband_qmlo.sh; [ -f $DRAW ] || DRAW=$HERE/draw_mloband.sh     # the same script (kt1 keeps a copy in $S)
[ -e $D ] && { echo "$D exists; refusing to overwrite"; exit 1; }
mkdir -p $D; cd $D
cp $SRC/ctrlg.$T.toml $SRC/env.sh .
source env.sh
L=$D/steps.log; : > $L
say(){ echo "$(date '+%m-%d %H:%M') $*" >> $L; }
B=$(dirname "$(readlink -f $BIN/lmf)")
say "tag=$TAG src=$SRC bin=$BIN build=$B rev=$(head -1 $B/../.ecalj_rev 2>/dev/null || head -1 $B/FROZEN_REV 2>/dev/null)"
say "binary: __QMLO_zNew in libecaljF_mp_gpu.so = $(strings $B/libecaljF_mp_gpu.so 2>/dev/null | grep -c __QMLO_zNew) (0 = a binary from before the QMLO_* names)"
say "settings: $(grep -hE '^nkabc|^n1n2n3|mlo_nkabc|^pwmode|^mixbeta' ctrlg.$T.toml | tr -s ' ' | tr '\n' ';')"
export ECALJ_MLO_MIX=1
say "precision: ${PREC:-fp32}  extra gwsc options: [${GWSC_EXTRA:-}]"
say "ECALJ_MLO_MIX=$ECALJ_MLO_MIX  GPUs: $(nvidia-smi --query-compute-apps=pid,process_name --format=csv,noheader | wc -l) processes before start"
t0=$(date +%s)
CUDA_VISIBLE_DEVICES=0,1 $BIN/gwsc 10 -np 60 -np2 2 --gpu --prec=${PREC:-fp32} --ntqxx --mlo ${GWSC_EXTRA:-} $T > gwsc10.log 2>&1; rc=$?
say "gwsc 10 rc=$rc secs=$(( $(date +%s)-t0 ))"
for it in $(seq 1 10); do
  f=llmf.${it}run
  say "iter $it mloON=$(grep -c 'MLO Sigma interpolation ON' $f 2>/dev/null) $(grep ehf $f 2>/dev/null | tail -1 | tr -s ' ') $(grep 'gap' $f 2>/dev/null | tail -1 | tr -s ' ')"
done
[ $rc -eq 0 ] || { say ABORT; exit 1; }
SD=$D/snap/final; mkdir -p $SD
for f in rst.$T sigm sigm.$T QMLO_SigRs QMLO_z HamRsMLO ctrlg.$T.toml env.sh __atm.$T efermi.lmf __mixsig __QMLO_mixsig; do
  cp -f $f $SD/ 2>/dev/null
done
printf "211 0.0 0.0 0.0 0.0 -1.0 0.0 GAMMA X\n0 !terminator\n" > $SD/syml.$T    # the path of the earlier chains
( cd $SD && source env.sh && CUDA_VISIBLE_DEVICES= $BIN/job_band $T -np 8 NoGnuplot > lb_path.log 2>&1 )
[ -s $SD/bnd001.spin1 ] && cp $SD/bnd001.spin1 $D/bndPMT_final.dat
mb=$(MLOBAND_BUILD=$B bash $DRAW $SD $D/mlobandwork $D/bnd_mlo_final.dat 2>&1 | tail -1)
say "bands: sigm drawing $( [ -s $SD/bnd001.spin1 ] && echo ok || echo FAILED ); MLO band: $mb"
say done; touch $S/$TAG.done
