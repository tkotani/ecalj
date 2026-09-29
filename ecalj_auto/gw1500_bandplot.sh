#!/bin/bash
# Band plots of the materials that gw1500_rerun.sh finished (2026-09-30).
#   gw1500_bandplot.sh <run dir> [<mpid> ...]      default: every material logged CONVERGED or MAXITER in <run dir>/rerun.log
# Env: BIN (frozen binaries, required), TOOLS_BIN (getsyml; default ~/bin_dev), NP (default 8), ENVSH
# In <run dir>/<mpid>/PlotBand: ctrlg, rst, sigm and atmpnu.* of the run, then lmfa (the atom file), getsyml and job_band.
# A material whose PlotBand already holds bnd001.spin1 is skipped.
set -u
R=$(realpath "$1"); shift
: "${BIN:?set BIN}"; TOOLS_BIN=${TOOLS_BIN:-$HOME/bin_dev}; NP=${NP:-8}
[ -n "${ENVSH:-}" ] && source "$ENVSH"
export PATH=$BIN:$PATH CUDA_VISIBLE_DEVICES= OMP_NUM_THREADS=1
export OMPI_MCA_hwloc_base_binding_policy=${OMPI_MCA_hwloc_base_binding_policy:-none}   # ranks unbound: they stay in the taskset of the caller (see gw1500_rerun.sh)
list=("$@"); [ ${#list[@]} -eq 0 ] && list=($(awk '$5=="CONVERGED" || $5=="MAXITER" {print $4}' "$R/rerun.log"))
for m in "${list[@]}"; do
  D=$R/$m; [ -s "$D/sigm" ] && [ -s "$D/rst.$m" ] || { echo "$m: no sigm or rst"; continue; }
  P=$D/PlotBand; [ -s "$P/bnd001.spin1" ] && continue
  rm -rf "${P:?}"; mkdir -p "$P"; cp -f "$D/ctrlg.$m.toml" "$D/rst.$m" "$D/sigm" "$D"/atmpnu.*.$m "$P/"
  ( cd "$P" && ln -sf sigm sigm.$m \
    && mpirun -np 1 "$BIN/lmfa" $m > llmfa 2>&1 < /dev/null \
    && "$TOOLS_BIN/getsyml" $m --nobzview > lgetsyml 2>&1 < /dev/null \
    && "$BIN/job_band" $m -np $NP --NoGnuplot > ljob_band 2>&1 < /dev/null )
  rm -rf "$P"/__* 2>/dev/null
  echo "$(date '+%F %T') $m bands: $(ls "$P" | grep -c 'bnd.*spin') files, gap $(grep -h 'gap =' "$P/llmf_ef" 2>/dev/null | tail -1 | sed 's/.*Ry = *//')"
done
