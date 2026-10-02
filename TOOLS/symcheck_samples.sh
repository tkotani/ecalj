#!/bin/bash
# Step S0 of ecaljdoc/MD/symmetry_spglib.md (2026-10-02): for every ctrlg.<sname>.toml of Samples (outside *_work), run lmchk and
# symfind.py --check in a scratch copy, and list whether ecalj and spglib find the same operations.
#   TOOLS/symcheck_samples.sh <out dir> [symprec]      -> <out dir>/summary.txt (one line per input), and the work dirs
OUT=${1:?usage: symcheck_samples.sh <out dir> [symprec]}; SP=${2:-1e-5}
ROOT=$(cd $(dirname $0)/.. && pwd)
mkdir -p $OUT; : > $OUT/summary.txt
export OMP_NUM_THREADS=1
find $ROOT/Samples -name 'ctrlg.*.toml' -not -path '*_work/*' | sort | while read f; do
  rel=${f#$ROOT/Samples/}; tag=$(dirname $rel | tr '/' '_'); s=$(basename $f .toml); s=${s#ctrlg.}
  w=$OUT/${tag}__$s; mkdir -p $w; cp $f $w/   # one dir per input (2026-10-02 05:23 fix: two ctrlg in one dir overwrote llmchk)
  ( cd $w && mpirun -np 1 lmchk $s > llmchk 2>&1 < /dev/null; rc=$?
    if [ $rc -ne 0 ]; then echo "$tag $s LMCHK_FAILED rc=$rc"; exit; fi
    r=$(python3 $ROOT/SRC/exec/symfind.py $s --symprec $SP --check llmchk 2>&1 < /dev/null | grep -E 'symcheck:|symfind:' | tail -1)
    echo "$tag $r" ) >> $OUT/summary.txt
done
echo "# done $(date '+%F %T')" >> $OUT/summary.txt
