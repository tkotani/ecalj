#!/bin/bash
# MLO-QSGW chain with a SNAPSHOT after every iteration, so any iteration can be
# restarted or re-examined later.  2026-09-25: the LiTi2O4 state was lost when kt1
# went down because only logs and bands had been kept.
#
#   usage: run_snap.sh <tag> <niter> [--mlo]        (add --mlofreeze with --mlo)
#
# Each snapshot holds everything needed to resume that iteration:
#   rst.<T>  sigm  QMLO_SigRs  QMLO_z  bnd_iter<N>.dat  llmf(tail)  steps.log
# HamRsMLO is written once (it is the frozen model) into snap/HamRsMLO.
set -u
TAG=$1; NITER=$2; MLO=${3:-}
S=/mnt/data1/LiTi2O4_kbt_runs
TB=$S/tbin_frozen
SRC=${LITI_SRC:-$S/liti_src}   # override to run the same chain from a different ctrlg
D=$S/$TAG; T=liti2o4
rm -rf $D; mkdir -p $D/snap; cd $D
cp $SRC/ctrlg.$T.toml $SRC/env.sh . 2>/dev/null
cp $SRC/ctrl.$T . 2>/dev/null
source env.sh; export PATH=$TB:$PATH
L=$D/steps.log; : > $L
say(){ echo "$(date '+%m-%d %H:%M') $*" >> $L; }

say "tag=$TAG niter=$NITER opt='$MLO'"
say "binary: __QMLO_zNew in libecaljF.so = $(strings /home/takao/ecalj/SRC/build_nvfortran/libecaljF.so 2>/dev/null | grep -c __QMLO_zNew) (0 = a binary from before the QMLO_* names)"
say "stale: sigm=$(ls sigm* 2>/dev/null|wc -l) mix=$(ls *mix* *MIX* 2>/dev/null|wc -l) mlo=$(ls HamRsMLO QMLO_SigRs QMLO_z 2>/dev/null|wc -l)  [must be 0]"
## any mixing history at all -- mix*, __mix*, __mixm, mixsigma -- means this is not a
## clean start.  Refuse rather than let an inherited history steer iteration 1.
if ls *mix* *MIX* >/dev/null 2>&1; then say "ABORT: mixing files present at start: $(ls *mix* *MIX* 2>/dev/null | tr '\n' ' ')"; exit 1; fi
## ECALJ_MLO_MIX=1 makes hqpe_sc Anderson-mix Sigma^MLO with [gw] mixbeta, the same way
## sigm is mixed.  Default is off (the MLO route otherwise runs at an effective beta=1),
## so record which way this chain ran -- the two are not comparable.
say "ECALJ_MLO_MIX=${ECALJ_MLO_MIX:-unset}  (1 = Sigma^MLO is Anderson-mixed at [gw] mixbeta)"
say "settings: $(grep -hE '^nkabc|^n1n2n3|mlo_nkabc|^pwmode|^nit' ctrlg.$T.toml | tr -s ' ' | tr '\n' ';')"

CUDA_VISIBLE_DEVICES=0,1 mpirun -np 1 $TB/lmfa $T > llmfa 2>&1; say "lmfa rc=$?"
CUDA_VISIBLE_DEVICES=0,1 mpirun -np 8 $TB/lmf $T > llmf_lda 2>&1
say "LDA rc=$? $(tail -3 llmf_lda | grep ehf | tr -s ' ')"
printf "211 0.0 0.0 0.0 0.0 -1.0 0.0 GAMMA X\n0 !terminator\n" > syml.$T
CUDA_VISIBLE_DEVICES= $TB/job_band $T -np 8 NoGnuplot > lb0.log 2>&1
cp bnd001.spin1 bnd_lda.dat; say "LDA band rc=$?"
mkdir -p snap/lda && cp rst.$T snap/lda/ 2>/dev/null && cp bnd_lda.dat snap/lda/
## The LDA lmf leaves its own density-mixing history.  gwsc redoes the LDA at iteration 1
## (there is no sigm yet) but only removes __mixm*<ext>*, so clear every mixing file here
## and say what was there.
say "after LDA, mixing files: $(ls *mix* *MIX* 2>/dev/null | tr '\n' ' ')  -> removed"
rm -f *mix* *MIX*

for it in $(seq 1 $NITER); do
  t0=$(date +%s)
  CUDA_VISIBLE_DEVICES=0,1 $TB/gwsc 1 -np 60 -np2 2 --gpu --mp --fp32 --ntqxx $MLO $T > g$it.log 2>&1; rc=$?
  say "iter $it rc=$rc secs=$(( $(date +%s)-t0 )) mloON=$(grep -c 'MLO Sigma interpolation ON' llmf 2>/dev/null) gwdrv=$(grep -c 'MLO Sigma interpolation ON' llmfgw01 2>/dev/null) mlomix=$(grep -c 'Sigma^MLO mixing section' lqpe 2>/dev/null) nofile=$(grep -c 'No mixing file' lqpe 2>/dev/null) mixfiles=[$(ls *mix* 2>/dev/null | tr '\n' ' ')] $(grep ehf llmf | tail -1 | tr -s ' ')"
  if [ $rc -ne 0 ]; then say ABORT; break; fi
  ## --- snapshot + band, both inside snap/iter<N>.  The band run must NOT happen in the
  ## chain directory: job_band's first stage is `lmf --quit=band`, which rewrites
  ## efermi.lmf (m_bndfp.f90:274) -- and the MLO window now reads that file.  Doing it in
  ## its own directory keeps the chain untouched and leaves the whole state kept.
  ## --mlo is essential: job_band forwards unknown args to lmf (parse_known_args), and
  ## without it c0_mlo is false, sigmlo_init returns at once and the band is drawn
  ## through the conventional sigm instead of QMLO_SigRs.
  SD=$D/snap/iter$it; mkdir -p $SD
  ## __mixsig / __QMLO_mixsig are the Anderson histories of sigm and of Sigma^MLO.  Without
  ## them a restart from this snapshot re-enters mixsigma with x_0 = 0, i.e. it damps the
  ## whole Sigma by beta once, which looks like a spurious jump in the band.
  for f in rst.$T sigm sigm.$T QMLO_SigRs QMLO_z HamRsMLO ctrlg.$T.toml ctrl.$T env.sh __atm.$T \
           efermi.lmf syml.$T __mixsig __QMLO_mixsig; do cp -f $f $SD/ 2>/dev/null; done
  tail -400 llmf > $SD/llmf.tail 2>/dev/null
  cp -f $L $SD/steps.log 2>/dev/null
  ## Two bands per iteration, both drawn inside/next to snap/iter<N>:
  ##  - snap/iter<N>/bnd001.spin1 : PMT band through the conventional sigm (job_band, no --mlo).
  ##    Not the band of the SCF state -- the SCF runs on Sigma^MLO -- but the same drawing as
  ##    the conventional chain, and job_band also writes qplist.dat, the k path mlo needs.
  ##  - $D/mloband/bnd_iter<N>.dat : THE MLO band, draw_mloband.sh = MLO model of the SCF
  ##    Hamiltonian (written at the mesh points with the stored chi~, QMLO_SigRs set aside,
  ##    unfrozen mlo).  NOT job_mlo --mlofreeze: that adds the current QMLO_SigRs (iteration-N
  ##    chi~) to the frozen HamRsMLO (step-0d LDA H, step-0d chi~) -- two bases, 0.3-0.5 eV
  ##    off.  It used to crash (zMLO overrun, e9f633d79); fixed, it would silently succeed
  ##    with that wrong band, so it is not run here at all.
  ( cd $SD && source env.sh && export PATH=$TB:$PATH && CUDA_VISIBLE_DEVICES= \
      $TB/job_band $T -np 8 NoGnuplot > lb_path.log 2>&1 )
  [ -s $SD/bnd001.spin1 ] && cp $SD/bnd001.spin1 $D/bndPMT_iter$it.dat
  mkdir -p $D/mloband
  mb=$(bash $S/draw_mloband.sh $SD $D/mlobandwork/iter$it $D/mloband/bnd_iter$it.dat 2>&1 | tail -1)
  say "   bands: sigm drawing $( [ -s $SD/bnd001.spin1 ] && echo ok || echo FAILED ); MLO band: $mb; snapshot $(du -sh $SD 2>/dev/null | cut -f1)"
done
say done; touch $S/$TAG.done
