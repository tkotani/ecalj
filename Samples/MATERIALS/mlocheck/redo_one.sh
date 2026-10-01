#!/bin/bash
# 2026-10-01 15:xx: the MLO step only, with the new default (semicore LO above EF-17 eV added; weight over the species)
mat=$1; np=$2; B=$HOME/bin
cd $(dirname $0)/$mat || exit 1
s=$(ls ctrlg.*.toml | sed 's/^ctrlg\.//; s/\.toml$//'); export OMP_NUM_THREADS=1
st() { echo "$(date '+%F %T') $mat $*" >> ../status_redo.log; }
so=$(python3 -c "import tomllib;print(tomllib.load(open('ctrlg.$s.toml','rb')).get('ham',{}).get('so',0))")
if [ "$so" = 1 ]; then
  timeout 4h $B/job_mlo_soc $s -np $np --NoGnuplot > ljob_mlo_soc 2>&1 || { st FAIL job_mlo_soc; exit 1; }
  [ -s band_MLO_spin2.dat ] || rm -f band_MLO_spin2.dat
else
  timeout 4h $B/job_mlo $s -np $np --nognuplot > ljob_mlo 2>&1 || { st FAIL job_mlo; exit 1; }
fi
[ -s band_MLO_spin1.dat ] || { st FAIL noMLOband; exit 1; }
python3 $B/mlo_bandplot.py . -o mlo_$mat.png --label $mat > lplot 2>&1
st DONE $(grep -o 'HamRsMTO= *[0-9]*' lmlo | tail -1)
