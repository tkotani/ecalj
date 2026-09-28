#!/bin/bash
# The six patterns {9^3, 6^3} x {tf32, fp32, fp64} side by side, rows = LDA and iterations 1..22 (user 2026-09-28):
# (update_six_rows.sh) fetch the MLO bands from kt1, fill each panel from the column's own run (the current code), and where that run has
# not got there yet, from what exists already, marked [09-26 code] / [09-27 code].  -> liti2o4_six_rows.png
set -u
SP=${WORK:-/tmp/liti2o4_six}     # local work: the fetched band files and the panel links
S=${RUNS_DIR:-/mnt/data1/LiTi2O4_kbt_runs}   # on kt1
C=$SP/six_cache; R=$SP/rows_six
mkdir -p $C
get(){ # <cache name> <remote glob...>: copy the band files of a run
  local n=$1; shift; mkdir -p $C/$n
  for g in "$@"; do rsync -q -t -e 'ssh -o BatchMode=yes' "kt1:$S/$g" $C/$n/ 2>/dev/null; done   # -t: keep the dates
}
for t in qmlo_k9_tf32_i15 qmlo_k9_fp32n qmlo_k9_fp64_i15 qmlo_k6_tf32_i15 qmlo_k6_fp32_i15 qmlo_k6_fp64_i15; do get $t "$t/bnd_mlo_iter*.dat"; done
get qmlo_k9_tf32n "qmlo_k9_tf32n/bnd_mlo_iter*.dat" "qmlo_k9_tf32n/bnd_mlo_final.dat"
get qmlo_k6_tf32n "qmlo_k6_tf32n/bnd_mlo_final.dat"
get qmlo_k6_gwsc10 "qmlo_k6_gwsc10/bnd_mlo_iter9.dat" "qmlo_k6_gwsc10/bnd_mlo_final.dat"
get liti_mlo_k9 "liti_mlo_k9/mloband/bnd_iter*.dat"
get liti_mlo_v9 "liti_mlo_v9/mloband/bnd_iter*.dat"
get lda9 "bench_hgw999/bnd_lda.dat"; get lda6 "liti_mlo_v9/bnd_lda.dat"
# When each iteration was computed: the date of llmf.<N>run (written at the end of iteration N) of every run
for t in qmlo_k9_tf32_i15 qmlo_k9_tf32n qmlo_k9_fp32n liti_mlo_k9 qmlo_k9_fp64_i15 qmlo_k6_tf32_i15 qmlo_k6_tf32n qmlo_k6_fp32_i15 \
         qmlo_k6_gwsc10 liti_mlo_v9 qmlo_k6_fp64_i15; do
  ssh -o BatchMode=yes kt1 "cd $S/$t 2>/dev/null && for f in llmf.*run; do [ -e \$f ] && echo \${f#llmf.} \$(date -r \$f '+%m-%d %H:%M'); done" 2>/dev/null \
    | sed 's/run / /' > $C/$t.times
done
rm -rf $R; mkdir -p $R/k9tf32 $R/k9fp32 $R/k9fp64 $R/k6tf32 $R/k6fp32 $R/k6fp64
put(){ # <column> <iteration> <file> <run> [label]: the first candidate that exists wins
  local col=$1 i=$2 f=$3 run=$4 lab=${5:-}
  local dst=$R/$col/bnd_iter$i.dat; [ $i = 0 ] && dst=$R/$col/bnd_lda.dat
  [ -e $dst ] && return 0; [ -s $f ] || return 1
  ln -s $f $dst; [ -n "$lab" ] && echo "$lab" > $dst.label
  local when=$(awk -v i=$i '$1==i {print $2, $3}' $C/$run.times 2>/dev/null)
  echo "${when:-$(date -r $f '+%m-%d %H:%M')}" > $dst.stamp    # the end of the iteration; else the date of the file
  local cd=$(basename $(dirname $f)) rd                            # where the band file lies on kt1 (kept in the npz)
  case $cd in liti_mlo_k9|liti_mlo_v9) rd=$cd/mloband;; lda9) rd=bench_hgw999;; lda6) rd=liti_mlo_v9;; *) rd=$cd;; esac
  echo "kt1:$S/$rd/$(basename $f)" > $dst.src; return 0
}
for i in $(seq 1 30); do   # iterations up to 30 (the tf32 continuation)
  put k9tf32 $i $C/qmlo_k9_tf32_i15/bnd_mlo_iter$i.dat qmlo_k9_tf32_i15
  put k9tf32 $i $C/qmlo_k9_tf32n/bnd_mlo_iter$i.dat qmlo_k9_tf32n "gwsc10 run"
  [ $i = 10 ] && put k9tf32 10 $C/qmlo_k9_tf32n/bnd_mlo_final.dat qmlo_k9_tf32n "gwsc10 run"
  put k9fp32 $i $C/qmlo_k9_fp32n/bnd_mlo_iter$i.dat qmlo_k9_fp32n
  put k9fp32 $i $C/liti_mlo_k9/bnd_iter$i.dat liti_mlo_k9 "09-26 code"
  put k9fp64 $i $C/qmlo_k9_fp64_i15/bnd_mlo_iter$i.dat qmlo_k9_fp64_i15
  put k6tf32 $i $C/qmlo_k6_tf32_i15/bnd_mlo_iter$i.dat qmlo_k6_tf32_i15
  [ $i = 10 ] && put k6tf32 10 $C/qmlo_k6_tf32n/bnd_mlo_final.dat qmlo_k6_tf32n "gwsc10 run"
  put k6fp32 $i $C/qmlo_k6_fp32_i15/bnd_mlo_iter$i.dat qmlo_k6_fp32_i15
  put k6fp32 $i $C/qmlo_k6_gwsc10/bnd_mlo_iter$i.dat qmlo_k6_gwsc10 "09-27 code"
  [ $i = 10 ] && put k6fp32 10 $C/qmlo_k6_gwsc10/bnd_mlo_final.dat qmlo_k6_gwsc10 "09-27 code"
  put k6fp32 $i $C/liti_mlo_v9/bnd_iter$i.dat liti_mlo_v9 "09-26 code"
  put k6fp64 $i $C/qmlo_k6_fp64_i15/bnd_mlo_iter$i.dat qmlo_k6_fp64_i15
done
for c in k9tf32 k9fp32 k9fp64; do put $c 0 $C/lda9/bnd_lda.dat none "job_band"; done    # LDA: the date of the file
for c in k6tf32 k6fp32 k6fp64; do put $c 0 $C/lda6/bnd_lda.dat none "job_band"; done
rows=$(ls $R/*/bnd_iter*.dat 2>/dev/null | sed 's/.*iter//; s/.dat//' | sort -n | uniq | tr '\n' ',' | sed 's/,$//')
cd "$(dirname "$(readlink -f "$0")")"
HILITE=1 EMPTY=frame MESH=9,6 MESHCOLS=9,9,9,6,6,6 ROWS=0,$rows DPI=90 DATAOUT=liti2o4_six_rows_data.npz python3 mlo_rows.py $R liti2o4_six_rows.png \
  k9tf32,k9fp32,k9fp64,k6tf32,k6fp32,k6fp64 \
  "9^3 tf32  MLO BAND,9^3 fp32  MLO BAND,9^3 fp64  MLO BAND,6^3 tf32  MLO BAND,6^3 fp32  MLO BAND,6^3 fp64  MLO BAND" \
  'tab:green,tab:green,tab:green,tab:green,tab:green,tab:green' noref
# 12 columns (user 2026-09-28): the t2g and, beside it, the eg bands b45-52 with b53 (the band whose part near x = 0.38 ran away
# in 9^3 tf32 iterations 18-22) of every pattern; an eg column is a link to the t2g column (mlo_rows.py: a name ending in _eg)
for c in k9tf32 k9fp32 k9fp64 k6tf32 k6fp32 k6fp64; do ln -sfn $c $R/${c}_eg; done
HILITE=1 EMPTY=frame MESH=9,6 MESHCOLS=9,9,9,9,9,9,6,6,6,6,6,6 ROWS=0,$rows DPI=80 COLW=3.6 TITLEFS=7 DATAOUT=liti2o4_six_rows12_data.npz \
  python3 mlo_rows.py $R liti2o4_six_rows12.png \
  k9tf32,k9tf32_eg,k9fp32,k9fp32_eg,k9fp64,k9fp64_eg,k6tf32,k6tf32_eg,k6fp32,k6fp32_eg,k6fp64,k6fp64_eg \
  "9^3 tf32 t2g,9^3 tf32 2.6-8 eV,9^3 fp32 t2g,9^3 fp32 2.6-8 eV,9^3 fp64 t2g,9^3 fp64 2.6-8 eV,6^3 tf32 t2g,6^3 tf32 2.6-8 eV,6^3 fp32 t2g,6^3 fp32 2.6-8 eV,6^3 fp64 t2g,6^3 fp64 2.6-8 eV" \
  'tab:green,tab:green,tab:green,tab:green,tab:green,tab:green,tab:green,tab:green,tab:green,tab:green,tab:green,tab:green' noref
# the occupied t2g pair only (OCC=1, 0 < x < 0.5), same panels (its numbers are in the npz above: the two lowest of the t2g)
HILITE=1 EMPTY=frame OCC=1 MESH=9,6 MESHCOLS=9,9,9,6,6,6 ROWS=$rows DPI=90 python3 mlo_rows.py $R liti2o4_six_occ.png \
  k9tf32,k9fp32,k9fp64,k6tf32,k6fp32,k6fp64 \
  "9^3 tf32  MLO BAND,9^3 fp32  MLO BAND,9^3 fp64  MLO BAND,6^3 tf32  MLO BAND,6^3 fp32  MLO BAND,6^3 fp64  MLO BAND" \
  'tab:green,tab:green,tab:green,tab:green,tab:green,tab:green' noref
# how much the MLO bands move per iteration (max |E(N)-E(N-1)|, within one run), with its numbers
python3 plot_mlo_conv_six.py $R liti2o4_six_conv.png liti2o4_six_conv_data.npz
for c in k9tf32 k9fp32 k9fp64 k6tf32 k6fp32 k6fp64; do
  echo "$c: $(ls $R/$c | grep -v -E 'label|stamp|src' | sed 's/bnd_iter//; s/bnd_lda/LDA/; s/.dat//' | sort -n | tr '\n' ' ')"
done
