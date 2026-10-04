#!/bin/bash
# GW1500 MLO: the standard model of ecalj on the QSGW80 result of one material (2026-10-04, user: the recipe must be
# reproducible by anyone; ecaljdoc manual/mlo.md "GW1500 の標準処方").
#
#   gw1500_mlo.sh <dir> [np]
#
# <dir> holds what gw1500_rerun.sh leaves in PlotBand/ (ctrlg, rst, sigm, atmpnu,
# syml, the QSGW80 bands bnd*.spin1 of job_band) plus efermi.lmf of the last lmf of the QSGW run (the window of the model
# is measured from its band edges). Steps:
#   1. [mlo] of the ctrlg rewritten with the present gwinit rules (gw1500_mlo_regen.py)
#   2. job_mlo (lmf --writeham --mlo, mlo): mlo_method 4, mlo_delta = mlo_w = 2 eV, Lowdin-orthogonalized MLOs. The
#      ctrlg has [ham] ssig 0.8, so the model is that of the QSGW80 Hamiltonian
#   3. mlo_bandcheck.py: the model against the QSGW80 bands on the path -> bandcheck.json, lbandcheck
# Writes mlo_version.txt (ecalj revision of the binaries) and one line to stdout: <m> DONE|FAIL ...
# Env: ECALJ_BIN (default ~/bin); MLO_LM2=1: baseline 2 (the mlo_lm2 rows of gwinit taken, gw1500_mlo_regen.py --lm2);
#      MLO_LM2=all: EH2 s,p on every atom but the transition metals, 4f and 5f (--lm2all, a test of 2026-10-04).
d=$1; np=${2:-4}; H=$(dirname "$(readlink -f "$0")"); B=${ECALJ_BIN:-$HOME/bin}
cd "$d" || { echo "$d FAIL nodir"; exit 1; }
m=$(ls ctrlg.*.toml 2>/dev/null | head -1 | sed 's/^ctrlg\.//; s/\.toml$//')
[ -n "$m" ] || { echo "$d FAIL no ctrlg"; exit 1; }
for f in ctrlg.$m.toml rst.$m sigm efermi.lmf syml.$m bnd001.spin1; do [ -e $f ] || { echo "$m FAIL missing $f"; exit 1; }; done
export OMP_NUM_THREADS=1
src=$(dirname "$(readlink -f "$B/lmf")"); rev=$(git -C "$src" rev-parse --short=9 HEAD 2>/dev/null || echo unknown)
dirty=$(git -C "$src" status --porcelain --untracked-files=no -- . 2>/dev/null | head -1)
echo "ecalj $rev${dirty:+ (modified)} binaries $(readlink -f "$B/lmf") $(date '+%F %T')" > mlo_version.txt
case "${MLO_LM2:-}" in 1) o=--lm2;; all) o=--lm2all;; *) o=;; esac
python3 "$H/gw1500_mlo_regen.py" . $m $o > lregen 2>&1 || { echo "$m FAIL regen"; exit 1; }
s=$(date +%s)
timeout 2h "$B/job_mlo" $m -np $np --nognuplot > ljob_mlo 2>&1; rc=$?
[ -s band_MLO_spin1.dat ] || { echo "$m FAIL job_mlo rc=$rc $(grep -m1 -i -E 'error|abort|stop' ljob_mlo lmlo lwriteham 2>/dev/null | cut -c1-100)"; exit 1; }
python3 "$B/mlo_bandcheck.py" . --json bandcheck.json > lbandcheck 2>&1
echo "$m DONE $(( $(date +%s) - s ))s $(grep -m1 -o 'CHECK [A-Z]*' lbandcheck)"
