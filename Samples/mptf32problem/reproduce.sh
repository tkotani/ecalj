#!/bin/bash
# Reproduce the TF32-vs-FP32 QSGW problem on mp-8196 (AgNO3).
#
# Runs one QSGW iteration (from LDA) in three GEMM-precision modes and prints
# sigma_m, the band sum sev, convergence flag, and the SEc of a few states:
#   cpu   : CPU  (no GPU)                -- reference truth
#   tf32  : --gpu --mp                   -- TF32 tensor cores (current default)
#   fp32  : --gpu --mp --fp32            -- true FP32 in the GEMMs (the fix)
#
# Expected (see reference/summary.txt):
#   cpu/fp32 : sigma_m ~ 529, sev ~ -204.7 eV, converged ('c').
#   tf32     : sigma_m ~ 55000, sev ~ -5059 eV, diverges ('i') -> QSGW fails.
# AgNO3's NO3- gives a strong local field (eps_head ~1536 vs eps_wLFC ~3.5,
# condition number ~400); TF32's ~1e-3 error is amplified by the dielectric
# inversion into ~40% error on W/SEc. FP32 (~1e-7) keeps single-precision
# storage but restores correctness (only ~7% slower than TF32).
#
# Usage:  ./reproduce.sh [NP] [NP2]
#   NP  = MPI ranks for lmf/lmfa   (default 8)
#   NP2 = GPU ranks for hgw/hsfp0  (default 2)
# Needs ctrlgenToml.py, gwsc and the GPU build (hgw_mp_gpu etc.) on PATH.
set -u
NP=${1:-8}
NP2=${2:-2}
HERE=$(cd "$(dirname "$0")" && pwd)
EXT=mp-8196

run_one() {            # $1=tag ; $2.. = gwsc precision flags
  local tag=$1; shift
  local d="$HERE/run_$tag"
  rm -rf "$d"; mkdir -p "$d"; cd "$d"
  cp "$HERE/ctrls.$EXT" .
  ctrlgenToml.py --ssig=0.8 "$EXT" > l_ctrlgen 2>&1 \
      || { echo "$tag: ctrlgenToml FAILED (see $d/l_ctrlgen)"; return 1; }
  gwsc "$@" 1 "$EXT" > l_gwsc 2>&1 || true   # tf32 is expected to fail
  local sm sev cv
  sm=$(grep -a "sum check of sigma_m" lqpe 2>/dev/null | tail -1 | grep -oE "[0-9]+\.[0-9]+" | tail -1)
  sev=$(tail -1 save.$EXT 2>/dev/null | grep -oE "sev\(eV\)=[-0-9.]+" | cut -d= -f2)
  cv=$(tail -1 save.$EXT 2>/dev/null | cut -c1)
  printf "%-5s  sigma_m=%-18s  sev=%-13s  conv=%-3s\n" "$tag" "${sm:-NA}" "${sev:-NA}" "${cv:-NA}"
  # SEc of states 4,5 (the ones TF32 destroys): columns = ... state SEx SExcore SEc ...
  awk '$1=="0.00000" && ($4=="4"||$4=="5"){printf "         state %s: SEx=%s SEc=%s\n",$4,$5,$7}' QPU 2>/dev/null
}

echo "=== mp-8196 (AgNO3) QSGW 1-iter: GEMM precision comparison (NP=$NP NP2=$NP2) ==="
run_one cpu  -np "$NP"
run_one tf32 -np "$NP" -np2 "$NP2" --gpu --mp
run_one fp32 -np "$NP" -np2 "$NP2" --gpu --mp --fp32
echo ""
echo "Compare with reference/summary.txt and reference/QPU.{cpu,tf32,fp32}.head"
