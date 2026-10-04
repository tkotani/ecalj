#!/bin/bash
# Run GW1500 materials again from scratch with the present ecalj (2026-09-30; GW1500_status.md).
# 2026-10-02 23:15 (the GW1500 database run): T_TETRAKBT (K) replaces t_tetrakbt of the template in ctrlg (user: -300 with
#    t_sigmaw 300), and the LDA bands are drawn too, in PlotBand_LDA/ from LDA/ (rst of the LDA step).
#
#   gw1500_rerun.sh <queue file> <worker name>
#
# The queue is a text file, one mpid per line (further columns are ignored).  A worker takes the first line under flock,
# so several workers can share one queue; start as many as the machine carries (the GPU programs of gwsc take the GPU
# locks of pylib/run_cmd.py, so the workers share the GPUs without a scheduler).
# Give each worker its own cores by taskset.  The ranks are left unbound (OMPI_MCA_hwloc_base_binding_policy=none below):
# OpenMPI otherwise binds rank i to core i of the machine whatever the taskset of the worker is, so the ranks of every
# worker sat on cores 0-11, and every 1-rank GPU program on core 0 (2026-09-30 02:49: a CPU job whose rank 0 was spinning
# in the NVIDIA driver on core 0 kept them from running, and kt1 stood still for 17 minutes).
#
# Env:
#   RUN_DIR     where <mpid>/ is made (required; a new place, never the directories of the old runs)
#   BIN         frozen binaries: a copy (cp -rL) of the bindir, so that a rebuild cannot change a run of hours (required)
#   TOOLS_BIN   bindir with the Python tools that live beside their modules: vasp2ctrl, getsyml (default ~/bin_dev)
#   POSCAR_DIR  default <this directory>/INPUT/gw1500/POSCARALL
#   NP          CPU ranks of lmf and of the CPU steps (default 16);  GPUS  the GPUs the worker may use (default 0,1)
#   PREC        fp32 (default) | tf32 | fp64;  MAXITER 10;  TOL 0.1 (eV);  CONV_QP 0.03 (eV, metals);  SSIG 0.8 (QSGW80)
#   GWSCCONV    default $BIN/gwscconv. It calls the gwsc of its own directory, so give a bindir that has both
#   LIMIT       seconds per material (default 21600); a material over the limit is logged TIMEOUT
#   ES_RMIN     if set (a.u.): empty spheres at the voids larger than this, before LDA (gw1500_addes.py; 2026-10-04: bands tied
#               to the vacuum level need them). The log line gets ES=<n>
#   GPU_OPT     default --gpu; empty for a machine without GPU
#   MLO         1 (default): make the MLO model after the band plot (gw1500_mlo.sh); 0: skip
#   T_TETRAKBT  if set, t_tetrakbt of the generated ctrlg (K; < 0: T=0 tetrahedron with a Gaussian smoothing of Im chi0)
#   ENVSH       a file to source first (PATH and LD_LIBRARY_PATH of the compiler and MPI)
#
# Per material, in $RUN_DIR/<mpid>:
#   POSCAR -> vasp2ctrl -> ctrlgenToml.py --ssig=$SSIG  (nkabc 8 8 8, n1n2n3 4 4 4, the [gw] template: t_tetrakbt 300, t_sigmaw 300)
#   gwscconv (LDA, then one QSGW iteration at a time until the gap of the last 3 iterations stays within TOL twice;
#   a metal, with no gap, until the eigenvalues within 5 eV of E_F change by less than CONV_QP twice, 2026-09-30)
#   band plot in PlotBand/ (getsyml, job_band), then the work files __* SEBK STDOUT are removed
#   MLO (MLO=1, default; 2026-10-04): the standard MLO model on the QSGW80 bands in PlotBand/ (gw1500_mlo.sh with the
#   binaries of BIN; result PlotBand/bandcheck.json, band_MLO_spin1.dat)
# Log: $RUN_DIR/rerun.log, one line per material:
#   <date time> <worker> <mpid> <verdict> iter=<n> gapLDA=<eV> gap=<eV> <seconds>s dqp=<eV>
#   verdict: CONVERGED | CONVERGED_METAL (no gap; judged by the eigenvalues) | MAXITER (no crash, not converged) |
#            NOGAP (gwscconv rc=3, only with --no-metal in GWSCCONV_OPT) | FAIL(rc) | TIMEOUT
#   dqp: the last max change of the eigenvalues within 5 eV of E_F between two iterations (gwscconv, 2026-09-30)
set -u
Q=$(realpath "$1"); WK=$2
HERE=$(dirname "$(readlink -f "$0")")
: "${RUN_DIR:?set RUN_DIR}" "${BIN:?set BIN}"
TOOLS_BIN=${TOOLS_BIN:-$HOME/bin_dev}; POSCAR_DIR=${POSCAR_DIR:-$HERE/INPUT/gw1500/POSCARALL}
NP=${NP:-16}; GPUS=${GPUS:-0,1}; PREC=${PREC:-fp32}; MAXITER=${MAXITER:-10}; TOL=${TOL:-0.1}; SSIG=${SSIG:-0.8}; LIMIT=${LIMIT:-21600}
CONV_QP=${CONV_QP:-0.03}; GWSCCONV=${GWSCCONV:-$BIN/gwscconv}; MLO=${MLO:-1}
[ -n "${ENVSH:-}" ] && source "$ENVSH"
export PATH=$BIN:$PATH CUDA_VISIBLE_DEVICES=$GPUS OMP_NUM_THREADS=1
export OMPI_MCA_hwloc_base_binding_policy=${OMPI_MCA_hwloc_base_binding_policy:-none}
mkdir -p "$RUN_DIR"; LOG=$RUN_DIR/rerun.log
say(){ echo "$(date '+%F %T') $WK $*" >> "$LOG"; }

pop(){  # take the first mpid of the queue
  ( flock 9
    m=$(awk 'NF && $1 !~ /^#/ {print $1; exit}' "$Q")
    [ -n "$m" ] && sed -i "0,/^$m\([[:space:]]\|\$\)/{/^$m\([[:space:]]\|\$\)/d}" "$Q"
    echo "$m" ) 9> "$Q.lock"
}

while :; do
  [ -e "$RUN_DIR/STOP" ] && { say "- STOP file found, worker ends"; break; }
  m=$(pop); [ -z "$m" ] && { say "- queue empty, worker ends"; break; }
  D=$RUN_DIR/$m
  [ -e "$D" ] && { say "$m SKIP ($D exists)"; continue; }
  [ -f "$POSCAR_DIR/POSCAR.$m" ] || { say "$m FAIL(no POSCAR)"; continue; }
  mkdir -p "$D"; cd "$D" || continue
  t0=$(date +%s)
  cp "$POSCAR_DIR/POSCAR.$m" POSCAR
  "$TOOLS_BIN/vasp2ctrl" POSCAR > lvasp2ctrl 2>&1 && cp ctrls.POSCAR.vasp2ctrl ctrls.$m
  "$BIN/ctrlgenToml.py" --ssig=$SSIG $m > lctrlgen 2>&1
  [ -s ctrlg.$m.toml ] || { say "$m FAIL(input) $(tail -1 lctrlgen | cut -c1-80)"; continue; }
  nes=
  if [ -n "${ES_RMIN:-}" ]; then
    python3 "$HERE/gw1500_addes.py" $m $ES_RMIN "$BIN" > laddes 2>&1 || { say "$m FAIL(ES) $(tail -1 laddes | cut -c1-80)"; continue; }
    nes=" ES=$(tail -1 laddes | awk '{print $2}')"
  fi
  if [ -n "${T_TETRAKBT:-}" ]; then
    sed -i "s/^t_tetrakbt *= *[-+0-9.eE]*/t_tetrakbt    = $T_TETRAKBT/" ctrlg.$m.toml
    grep -Eq "^t_tetrakbt *= *$T_TETRAKBT( |$)" ctrlg.$m.toml || { say "$m FAIL(t_tetrakbt not set)"; continue; }
  fi
  timeout -k 60 $LIMIT "$GWSCCONV" -np $NP -np2 1 ${GPU_OPT---gpu} --prec=$PREC $m --conv-tol $TOL --conv-qp $CONV_QP --max-iter $MAXITER \
    ${GWSCCONV_OPT:-} > osgw.conv.out 2>&1
  rc=$?
  it=$(ls -d QSGW.*run 2>/dev/null | sed 's/QSGW\.//; s/run//' | sort -n | tail -1)
  glda=$(grep -h 'gap =' llmf_lda 2>/dev/null | tail -1 | sed 's/.*Ry = *//; s/ *eV.*//')
  g=$(grep -h 'gap =' llmf 2>/dev/null | tail -1 | sed 's/.*Ry = *//; s/ *eV.*//')
  dqp=$(grep -ho 'max QP change[^=]*= *[0-9][0-9.]*' osgw.conv.out 2>/dev/null | tail -1 | sed 's/.*= *//')
  if   grep -q 'gwscconv: CONVERGED.*(metal' osgw.conv.out; then v=CONVERGED_METAL
  elif grep -q 'gwscconv: CONVERGED' osgw.conv.out; then v=CONVERGED
  elif [ $rc -eq 124 ] || [ $rc -eq 137 ]; then v=TIMEOUT
  elif [ $rc -eq 3 ]; then v=NOGAP
  elif [ $rc -ne 0 ]; then v="FAIL($rc)"
  else v=MAXITER; fi
  if [ "$v" = CONVERGED ] || [ "$v" = CONVERGED_METAL ] || [ "$v" = MAXITER ]; then
    # atmpnu.* too: the input has ham.readp = true, so lmf reads the pnu from them (2026-09-30 01:10: without them the
    # band step ended in a segmentation fault; the band plots of the first run were made afterwards by gw1500_bandplot.sh)
    mkdir -p PlotBand && cp -f ctrlg.$m.toml rst.$m sigm __atm.$m atmpnu.*.$m PlotBand/ 2>/dev/null
    ( cd PlotBand && ln -sf sigm sigm.$m && "$TOOLS_BIN/getsyml" $m --nobzview > lgetsyml 2>&1 \
      && "$BIN/job_band" $m -np $NP --NoGnuplot > ljob_band 2>&1 )
    cp -f efermi.lmf PlotBand/ 2>/dev/null
    [ "$MLO" = 1 ] && [ -s PlotBand/bnd001.spin1 ] && ECALJ_BIN=$BIN "$HERE/gw1500_mlo.sh" PlotBand $NP > PlotBand/lmlo_recipe 2>&1
    if [ -s LDA/rst.$m ]; then
      mkdir -p PlotBand_LDA && cp -f LDA/ctrlg.$m.toml LDA/rst.$m LDA/atmpnu.*.$m PlotBand_LDA/ 2>/dev/null
      cp -f __atm.$m PlotBand_LDA/ 2>/dev/null
      ( cd PlotBand_LDA && "$TOOLS_BIN/getsyml" $m --nobzview > lgetsyml 2>&1 \
        && "$BIN/job_band" $m -np $NP --NoGnuplot > ljob_band 2>&1 )
    fi
  fi
  rm -rf __* SEBK STDOUT PlotBand/__* PlotBand_LDA/__* 2>/dev/null
  say "$m $v iter=${it:-0} gapLDA=${glda:-none} gap=${g:-none} $(( $(date +%s)-t0 ))s dqp=${dqp:-none}${nes}"
done
