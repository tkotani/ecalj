#!/bin/bash
# 2026-10-01 15:xx: the baseline (auto semicore) plus EH2 s,p on the cations (atoms that are not N O F P S Cl As Se Br Sb Te I);
# the MLO step only, on the DFT of ~/work/mlocheck_auto_20261001/<Material>
m=$1; np=$2; B=$HOME/bin; W=$(cd $(dirname $0); pwd)
src=/home/takao/work/mlocheck_auto_20261001/$m
rm -rf $W/$m; mkdir -p $W/$m
rsync -a --exclude='__HamiltonianPMT' --exclude='PROCAR*' --exclude='__amlo.data' --exclude='band_MLO*' --exclude='lmlo' \
      --exclude='HamRsMLO' --exclude='mlo_*.png' --exclude='ljob_mlo*' --exclude='lwriteham' $src/ $W/$m/
cd $W/$m || exit 1
s=$(ls ctrlg.*.toml | sed 's/^ctrlg\.//; s/\.toml$//'); export OMP_NUM_THREADS=1
st() { echo "$(date '+%F %T') $m $*" >> $W/status.log; }
python3 - $s <<'PY'
import sys, tomllib
s = sys.argv[1]; f = f'ctrlg.{s}.toml'; t = open(f).read(); d = tomllib.loads(t)
z = {x['atom']: x['z'] for x in d['spec']}
anions = {7, 8, 9, 15, 16, 17, 33, 34, 35, 51, 52, 53}
rows = [f"{i} {x['atom']:4s} 1 2 3 4" for i, x in enumerate(d['site'], 1) if z[x['atom']] not in anions and z[x['atom']] > 0]
t = t.replace('mlo_lm = """\n', 'mlo_lm2 = """\n' + '\n'.join(rows) + '\n"""\nmlo_lm = """\n', 1)
open(f, 'w').write(t)
PY
so=$(python3 -c "import tomllib;print(tomllib.load(open('ctrlg.$s.toml','rb')).get('ham',{}).get('so',0))")
if [ "$so" = 1 ]; then timeout 4h $B/job_mlo_soc $s -np $np --NoGnuplot > ljob_mlo_soc 2>&1 || { st FAIL job_mlo_soc; exit 1; }
  [ -s band_MLO_spin2.dat ] || rm -f band_MLO_spin2.dat
else timeout 4h $B/job_mlo $s -np $np --nognuplot > ljob_mlo 2>&1 || { st FAIL job_mlo $(grep -h -m1 -E "completeness loss|Error|error" lmlo ljob_mlo 2>/dev/null | cut -c1-80); exit 1; }; fi
[ -s band_MLO_spin1.dat ] || { st FAIL noMLOband; exit 1; }
python3 $B/mlo_bandplot.py . -o mlo_${m}_eh2cat.png --label "$m + EH2 s,p on cations" > lplot 2>&1
st DONE $(grep -o 'HamRsMTO= *[0-9]*' lmlo | tail -1)
