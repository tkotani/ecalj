#!/bin/bash
# 2026-10-01 mlocheck, so=1 materials: DFT(so=1) -> band(so=1) -> job_mlo_soc (SOC as perturbation on the MLO model) -> plot
mat=$1; np=$2; B=$HOME/bin
cd $(dirname $0)/$mat || exit 1
s=$(ls ctrlg.*.toml | sed 's/^ctrlg\.//; s/\.toml$//')
export OMP_NUM_THREADS=1
st() { echo "$(date '+%F %T') $mat $*" >> ../status.log; }
st START soc np=$np
if [ ! -e bnd001.spin1 ]; then
  mpirun -np 1 $B/lmfa $s > llmfa 2>&1 || { st FAIL lmfa; exit 1; }
  timeout 4h mpirun -np $np $B/lmf $s > llmf 2>&1 || { st FAIL lmf rc=$?; exit 1; }
  st lmf done "$(grep -E '^[ ]*c |^x ' llmf | tail -1)"
  python3 $B/getsyml $s --nobzview > lgetsyml 2>&1 || { st FAIL getsyml; exit 1; }
  timeout 2h $B/job_band $s -np $np --nognuplot > ljob_band 2>&1 || { st FAIL job_band; exit 1; }
fi
timeout 4h $B/job_mlo_soc $s -np $np --NoGnuplot > ljob_mlo_soc 2>&1 || { st FAIL job_mlo_soc; exit 1; }
ls band_MLO_spin1.soc.dat >/dev/null 2>&1 || { st FAIL job_mlo_soc noMLOband; exit 1; }
st DONE soc
