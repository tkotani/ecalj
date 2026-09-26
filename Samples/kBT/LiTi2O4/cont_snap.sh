#!/bin/bash
# Continue an existing run_snap.sh chain: iterations $2..$3, snapshot + --mlo band each.
set -u
TAG=$1; FROM=$2; TO=$3; MLO=${4:-}
S=/mnt/data1/LiTi2O4_kbt_runs; TB=$S/tbin_frozen; T=liti2o4; D=$S/$TAG
cd $D; source env.sh; export PATH=$TB:$PATH
L=$D/steps.log
say(){ echo "$(date '+%m-%d %H:%M') $*" >> $L; }
say "continuing iterations $FROM..$TO"
for it in $(seq $FROM $TO); do
  t0=$(date +%s)
  CUDA_VISIBLE_DEVICES=0,1 $TB/gwsc 1 -np 60 -np2 2 --gpu --mp --fp32 --ntqxx $MLO $T > g$it.log 2>&1; rc=$?
  say "iter $it rc=$rc secs=$(( $(date +%s)-t0 )) mloON=$(grep -c 'MLO Sigma interpolation ON' llmf 2>/dev/null) gwdrv=$(grep -c 'MLO Sigma interpolation ON' llmfgw01 2>/dev/null) zslot=$(ls QMLO_z 2>/dev/null | wc -l) $(grep ehf llmf | tail -1 | tr -s ' ')"
  if [ $rc -ne 0 ]; then say ABORT; break; fi
  SD=$D/snap/iter$it; mkdir -p $SD
  ## __mixsig / __QMLO_mixsig: the Anderson histories.  Losing them on a restart damps
  ## Sigma by beta against x_0 = 0 once, which reads as a spurious jump in the band.
  for f in rst.$T sigm sigm.$T QMLO_SigRs QMLO_z HamRsMLO ctrlg.$T.toml ctrl.$T env.sh __atm.$T efermi.lmf syml.$T __mixsig __QMLO_mixsig; do cp -f $f $SD/ 2>/dev/null; done
  tail -400 llmf > $SD/llmf.tail 2>/dev/null; cp -f $L $SD/steps.log 2>/dev/null
  ( cd $SD && source env.sh && export PATH=$TB:$PATH && CUDA_VISIBLE_DEVICES= $TB/job_band $T -np 8 NoGnuplot --mlo > lb.log 2>&1 )
  if [ -f $SD/bnd001.spin1 ]; then cp $SD/bnd001.spin1 $D/bnd_iter$it.dat
    say "   band ok (--mlo, mloON=$(grep -c 'MLO Sigma interpolation ON' $SD/llmf_band 2>/dev/null))"
  else say "   BAND FAILED"; fi
done
say "continue done"; touch $S/$TAG.cont.done
