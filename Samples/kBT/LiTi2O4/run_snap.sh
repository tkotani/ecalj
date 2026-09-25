#!/bin/bash
# MLO-QSGW chain with a SNAPSHOT after every iteration, so any iteration can be
# restarted or re-examined later.  2026-09-25: the LiTi2O4 state was lost when kt1
# went down because only logs and bands had been kept.
#
#   usage: run_snap.sh <tag> <niter> [--mlo]        (add --mlofreeze with --mlo)
#
# Each snapshot holds everything needed to resume that iteration:
#   rst.<T>  sigm  SigRsMLO  bnd_iter<N>.dat  llmf(tail)  steps.log
# HamRsMLO is written once (it is the frozen model) into snap/HamRsMLO.
set -u
TAG=$1; NITER=$2; MLO=${3:-}
S=/mnt/data1/LiTi2O4_kbt_runs
TB=$S/tbin_frozen
SRC=$S/liti_src
D=$S/$TAG; T=liti2o4
rm -rf $D; mkdir -p $D/snap; cd $D
cp $SRC/ctrlg.$T.toml $SRC/env.sh . 2>/dev/null
cp $SRC/ctrl.$T . 2>/dev/null
source env.sh; export PATH=$TB:$PATH
L=$D/steps.log; : > $L
say(){ echo "$(date '+%m-%d %H:%M') $*" >> $L; }

say "tag=$TAG niter=$NITER opt='$MLO'"
say "binary: ZmloRef in libecaljF.so = $(strings /home/takao/ecalj/SRC/build_nvfortran/libecaljF.so 2>/dev/null | grep -c ZmloRef) (0 = z^MLO rebuilt per lmf run)"
say "stale: sigm=$(ls sigm* 2>/dev/null|wc -l) mix=$(ls __mix* 2>/dev/null|wc -l) mlo=$(ls HamRsMLO SigRsMLO ZmloRef 2>/dev/null|wc -l)  [must be 0]"
say "settings: $(grep -hE '^nkabc|^n1n2n3|mlo_nkabc|^pwmode|^nit' ctrlg.$T.toml | tr -s ' ' | tr '\n' ';')"

CUDA_VISIBLE_DEVICES=0,1 mpirun -np 1 $TB/lmfa $T > llmfa 2>&1; say "lmfa rc=$?"
CUDA_VISIBLE_DEVICES=0,1 mpirun -np 8 $TB/lmf $T > llmf_lda 2>&1
say "LDA rc=$? $(tail -3 llmf_lda | grep ehf | tr -s ' ')"
printf "211 0.0 0.0 0.0 0.0 -1.0 0.0 GAMMA X\n0 !terminator\n" > syml.$T
CUDA_VISIBLE_DEVICES= $TB/job_band $T -np 8 NoGnuplot > lb0.log 2>&1
cp bnd001.spin1 bnd_lda.dat; say "LDA band rc=$?"
mkdir -p snap/lda && cp rst.$T snap/lda/ 2>/dev/null && cp bnd_lda.dat snap/lda/

for it in $(seq 1 $NITER); do
  t0=$(date +%s)
  CUDA_VISIBLE_DEVICES=0,1 $TB/gwsc 1 -np 60 -np2 2 --gpu --mp --fp32 --ntqxx $MLO $T > g$it.log 2>&1; rc=$?
  say "iter $it rc=$rc secs=$(( $(date +%s)-t0 )) mloON=$(grep -c 'MLO Sigma interpolation ON' llmf 2>/dev/null) gwdrv=$(grep -c 'MLO Sigma interpolation ON' llmfgw01 2>/dev/null) $(grep ehf llmf | tail -1 | tr -s ' ')"
  if [ $rc -ne 0 ]; then say ABORT; break; fi
  ## --- snapshot + band, both inside snap/iter<N>.  The band run must NOT happen in the
  ## chain directory: job_band's first stage is `lmf --quit=band`, which rewrites
  ## efermi.lmf (m_bndfp.f90:274) -- and the MLO window now reads that file.  Doing it in
  ## its own directory keeps the chain untouched and leaves the whole state kept.
  ## --mlo is essential: job_band forwards unknown args to lmf (parse_known_args), and
  ## without it c0_mlo is false, sigmlo_init returns at once and the band is drawn
  ## through the conventional sigm instead of SigRsMLO.
  SD=$D/snap/iter$it; mkdir -p $SD
  for f in rst.$T sigm sigm.$T SigRsMLO HamRsMLO ctrlg.$T.toml ctrl.$T env.sh __atm.$T \
           efermi.lmf syml.$T; do cp -f $f $SD/ 2>/dev/null; done
  cp -f ZmloSig.* $SD/ 2>/dev/null
  tail -400 llmf > $SD/llmf.tail 2>/dev/null
  cp -f $L $SD/steps.log 2>/dev/null
  ## Band in the MLO space, not the PMT space.  H^MLO(R) and Sigma^MLO(R) are both
  ## matrices over the MLO channels and therefore q-periodic, so the reduced eigenproblem
  ## can be solved at ANY k with no lift to PMT and no z^MLO(k) -- which is the one thing
  ## that cannot be interpolated.  Drawing in PMT space instead mixes the stored chi~ at
  ## the mesh points with a rebuilt one in between, and that mixture, not the method, is
  ## what made the 2026-09-25 band plots look rough.  --mlofreeze keeps HamRsMLO so the
  ## channel list stays the one SigRsMLO is a matrix over.
  ( cd $SD && source env.sh && export PATH=$TB:$PATH && CUDA_VISIBLE_DEVICES= \
      $TB/job_mlo $T -np 8 --mlofreeze --mlo > lb.log 2>&1 )
  BND=$(ls $SD/band_MLO_spin1.dat $SD/bnd001.spin1 2>/dev/null | head -1)
  if [ -n "$BND" ]; then
    cp $BND $D/bnd_iter$it.dat
    say "   band ok (MLO space, in snap/iter$it; $(basename $BND)); snapshot $(du -sh $SD 2>/dev/null | cut -f1)"
  else
    say "   BAND FAILED: $(tail -2 $SD/lb.log 2>/dev/null | tr '\n' ' ')"
  fi
done
say done; touch $S/$TAG.done
