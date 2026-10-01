#!/bin/bash
# 2026-10-02: lmchk writes symmetry.<sname>.json with the vendored spglib (C); compare its operations, as a set, with those of
# symfind.py (spglib from Python) for every ctrlg.<sname>.toml of Samples.   symcheck_spglib_c.sh <out dir>  -> <out dir>/summary.txt
OUT=${1:?usage: symcheck_spglib_c.sh <out dir>}
ROOT=$(cd $(dirname $0)/.. && pwd)
mkdir -p $OUT; : > $OUT/summary.txt
export OMP_NUM_THREADS=1
find $ROOT/Samples -name 'ctrlg.*.toml' -not -path '*_work/*' | sort | while read f; do
  rel=${f#$ROOT/Samples/}; tag=$(dirname $rel | tr '/' '_'); s=$(basename $f .toml); s=${s#ctrlg.}
  w=$OUT/${tag}__$s; rm -rf $w; mkdir -p $w; cp $f $w/
  ( cd $w && mpirun -np 1 lmchk $s > llmchk 2>&1 < /dev/null; rc=$?
    if [ $rc -ne 0 ]; then echo "$tag $s LMCHK_FAILED rc=$rc"; exit; fi
    if [ ! -f symmetry.$s.json ]; then echo "$tag $s NO_JSON (symgrp: $(grep -m1 '^symgrp' ctrlg.$s.toml | cut -c1-40))"; exit; fi
    python3 $ROOT/SRC/exec/symfind.py $s --out py.$s.json > /dev/null 2>&1 < /dev/null
    python3 - $s <<'PY'
import json,sys,numpy as np
s=sys.argv[1]
a=json.load(open(f'symmetry.{s}.json')); b=json.load(open(f'py.{s}.json'))
S=lambda d:{(tuple(np.array(o['rotation']).ravel()),tuple(np.round(np.array(o['translation'])%1,6)%1),o['time_reversal']) for o in d['operations']}
print(s, a['n_operations'], a['spacegroup']['international'], 'SAME' if S(a)==S(b) else 'DIFFERENT')
PY
  ) | sed "s/^/$tag /" >> $OUT/summary.txt
done
echo "# done $(date '+%F %T')" >> $OUT/summary.txt
