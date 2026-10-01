#!/bin/bash
# MLO against the Wannier functions (maximally localized, genMLWFx): cRPA U, J and the spread, Ni d and SrVO3 t2g.
# Runs only at the commit tagged last-wannier (the Wannier path was removed after it). 2026-10-02.
#   ./run.sh [-np N]      -> work_<sname>_{mlo,wan}/ and the table on stdout
set -e
NP=${2:-4}; [ "$1" = "-np" ] || NP=4
H=$(cd $(dirname $0) && pwd)
export OMP_NUM_THREADS=1
for s in ni srvo3; do
  for m in mlo wan; do
    w=$H/work_${s}_$m; rm -rf "${w:?}"; mkdir -p $w; cp $H/$s/* $w/
    ( cd $w && mpirun -np 1 lmfa $s > llmfa && mpirun -np $NP lmf $s > llmf && getsyml $s --nobzview > lgetsyml 2>&1 \
      && job_band $s -np $NP --NoGnuplot > ljob_band 2>&1
      if [ $m = mlo ]; then job_mlo $s -np $NP --NoGnuplot > ljob_mlo 2>&1 && job_mloW $s -np $NP --crpa > ljob_mloW 2>&1 \
                             && mlo_spread.py $s -np $NP > lspread 2>&1
      else genMLWFx $s -np $NP > lgenMLWFx 2>&1; fi )
  done
done
cd $H
python3 ujk.py work_ni_mlo work_ni_wan work_srvo3_mlo work_srvo3_wan
for s in ni srvo3; do tail -8 work_${s}_mlo/lspread; python3 wan_spread.py work_${s}_wan; done
