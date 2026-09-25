#!/bin/bash
# Redraw the earlier conventional 9^3 chain on the 211-point Gamma->X line the MLO chains use.
# Source: n999_P_smearx0_nofilter_mix05/band/iter<N> (15 iterations, each a full state:
# rst, sigm, ctrlg, efermi.lmf).  That chain differs from liti_mlo_k9 in two settings:
# nkabc = 16^3 (the GW mesh is 9^3, so Sigma is interpolated into the density) and
# pwmode = 1 (|G|).  Drawn exactly as refbands.sh drew the 6^3 reference: job_band, no --mlo,
# CPU only (CUDA_VISIBLE_DEVICES empty) so the running chain keeps both GPUs.
set -u
S=/mnt/data1/LiTi2O4_kbt_runs; TB=$S/tbin_frozen; T=liti2o4
SRC=$S/n999_P_smearx0_nofilter_mix05/band; OUTD=$S/refbands999_211; mkdir -p $OUTD
source $S/liti_src/env.sh; export PATH=$TB:$PATH; export CUDA_VISIBLE_DEVICES=
for n in $(seq 1 15); do
  tag=iter$n; W=$S/_rb9_$tag; rm -rf $W; cp -a $SRC/$tag $W || { echo "$tag missing"; continue; }
  cd $W; rm -f bnd*.spin1
  printf "211 0.0 0.0 0.0 0.0 -1.0 0.0 GAMMA X\n0 !terminator\n" > syml.$T
  t0=$(date +%s); $TB/job_band $T -np 8 NoGnuplot > jb.log 2>&1
  if [ -f bnd001.spin1 ]; then cp bnd001.spin1 $OUTD/bnd_$tag.dat; echo "$(date +%H:%M) $tag ok $(( $(date +%s)-t0 ))s $(grep -c '' bnd001.spin1) lines"
  else echo "$(date +%H:%M) $tag FAILED $(tail -2 jb.log | tr '\n' ' ')"; fi
  cd $S; rm -rf $W
done
echo DONE
