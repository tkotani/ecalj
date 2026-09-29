#!/bin/bash
# Temperature scan of QSGW: the same input run at several electron temperatures (Samples/kBT/README.md, 2026-09-30).
#
#   scan_T.sh <input ctrlg.<sname>.toml> <work dir> <niter> <np> <setting> [<setting> ...]
#
#   <setting> = <label>:<t_tetrakbt>:<t_sigmaw>     e.g.  T0:0:300  T1000:1000:1000  G1000:-1000:1000
#       t_tetrakbt (K): 0 = chi0 at T=0, T > 0 = finite-temperature tetrahedron, -T = Gaussian on Im chi0 (ecaljdoc kBT, table 1)
#       t_sigmaw   (K): Fermi-Dirac width of the levels in Sigma
#   Env: BIN (default: the directory of gwsc found in PATH), GWSC_OPT (more gwsc options, e.g. "--gpu --prec=tf32 -np2 1"),
#        SYML (a syml.<sname> file to draw the bands along; default: getsyml makes one)
#
# For each setting: <work dir>/<label>/ gets the input with the two keys replaced, then gwsc <niter> from LDA and job_band.
# Kept in <work dir>/<label>/: everything gwsc leaves.  collect_T.py reads llmf.<N>run, QPU.<N>run, EFERMI* and bnd*.spin*.
set -u
IN=$(realpath "$1"); W=$(realpath -m "$2"); NIT=$3; NP=$4; shift 4
S=$(basename "$IN" .toml); S=${S#ctrlg.}
BIN=${BIN:-$(dirname "$(command -v gwsc)")}
for set in "$@"; do
  IFS=: read -r lab tt ts <<< "$set"
  D=$W/$lab; rm -rf "${D:?}"; mkdir -p "$D"; cd "$D" || exit 1
  # replace the two keys (both must be present in the input: t_tetrakbt is required since 2026-09-28)
  grep -q '^t_tetrakbt' "$IN" && grep -q '^t_sigmaw' "$IN" || { echo "$IN: t_tetrakbt and t_sigmaw must both be written in [gw]"; exit 1; }
  sed -e "s/^t_tetrakbt .*/t_tetrakbt    = $tt     # (K) set by scan_T.sh ($lab)/" \
      -e "s/^t_sigmaw .*/t_sigmaw      = $ts     # (K) set by scan_T.sh ($lab)/" "$IN" > ctrlg.$S.toml
  t0=$(date +%s)
  $BIN/gwsc $NIT -np $NP ${GWSC_OPT:-} $S > lgwsc 2>&1; rc=$?
  if [ -n "${SYML:-}" ]; then cp "$SYML" syml.$S; else $BIN/getsyml $S > lgetsyml 2>&1; fi
  $BIN/job_band $S -np $NP NoGnuplot > ljob_band 2>&1; rb=$?
  echo "$lab t_tetrakbt=$tt t_sigmaw=$ts gwsc rc=$rc job_band rc=$rb $(( $(date +%s)-t0 )) s  $(grep -h 'gap' llmf 2>/dev/null | tail -1 | tr -s ' ')" | tee -a "$W/scan_T.log"
done
