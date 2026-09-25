#!/bin/bash
# Wait for iteration N of a chain, copy its logs at once, wait for its band, check it.
#   usage: watch_iter.sh <chain dir> <N> [expect_mix(1|0)]
D=$1; N=$2; MIX=${3:-1}
cd "$D" || exit 1
until grep -q "iter $N rc=" steps.log || grep -q ABORT steps.log; do sleep 10; done
C=chk/iter$N; mkdir -p $C; cp -f lqpe lmlo_sigr llmf llmfgw01 $C/ 2>/dev/null
t=0; until [ -f snap/iter$N/bnd001.spin1 ] || grep -qE "iter $N .*BAND|band ok" <(sed -n "/iter $N rc/,\$p" steps.log) || [ $t -gt 60 ]; do sleep 5; t=$((t+1)); done
bash "$(dirname "$0")/check_iter.sh" "$D" "$N" "$MIX" "$C"
