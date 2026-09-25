#!/bin/bash
# Per-iteration health check for a run_snap.sh / cont_snap.sh MLO-QSGW chain.
#   usage: check_iter.sh <chain dir> <N> [expect_mix(1|0)] [log dir]
# [log dir] holds lqpe/lmlo_sigr/llmf/llmfgw01 copied the moment "iter N rc=" appeared
# (watch_iter.sh does that); reading the live files later would see iteration N+1.
# Prints the facts, then one line per ANOMALY.  Known-normal behaviour is listed
# below so that it is NOT reported -- mistaking normal for abnormal wastes as much
# time as the reverse:
#   - "BAND FAILED: ... mlo --mlofreeze --mlo": the MLO stage of job_band dies on
#     LiTi2O4 (A2').  lmf --band runs before it, so bnd001.spin1 is still produced.
#   - gwdrv=0 at N=1: the GW driver of iteration 1 has no SigRsMLO yet.
#   - x_0 from __SigmMLO.q.prev = F at N=1: nothing to read yet.
#   - mloON=0 in the band's llmf_band: bands are drawn through sigm.
#   - E_HF moving by eV in early iterations.
D=$1; N=$2; MIX=${3:-1}; L=${4:-.}
## Expectations measured on the 6^3 chain.  They depend on the mesh and the MPI layout, so a
## different setup must pass its own (e.g. 9^3: CHECK_MAXSECS=15000).  A hard-coded 6^3 value
## would call every 9^3 iteration a hang -- normal reported as abnormal.
EXP_MLOON=${CHECK_MLOON:-16}; EXP_GWDRV=${CHECK_GWDRV:-36}; MAXSECS=${CHECK_MAXSECS:-2000}
cd "$D" || { echo "ANOMALY: no chain dir $D"; exit 1; }
S=snap/iter$N
line=$(grep "iter $N rc=" steps.log | tail -1)
echo "steps : $line"
rc=$(echo "$line" | grep -oE "rc=[0-9]+" | cut -d= -f2)
secs=$(echo "$line" | grep -oE "secs=[0-9]+" | cut -d= -f2)
mloON=$(echo "$line" | grep -oE "mloON=[0-9]+" | cut -d= -f2)
gwdrv=$(echo "$line" | grep -oE "gwdrv=[0-9]+" | cut -d= -f2)
mlomix=$(echo "$line" | grep -oE "mlomix=[0-9]+" | cut -d= -f2)
[ -z "$line" ] && echo "ANOMALY: no steps.log line for iter $N"
[ "${rc:-1}" != 0 ] && echo "ANOMALY: rc=$rc"
grep -q ABORT steps.log && echo "ANOMALY: ABORT in steps.log"
[ "${mloON:-0}" != "$EXP_MLOON" ] && echo "ANOMALY: mloON=$mloON (expect $EXP_MLOON)"
[ "$N" -ge 2 ] && [ "${gwdrv:-0}" != "$EXP_GWDRV" ] && echo "ANOMALY: gwdrv=$gwdrv (expect $EXP_GWDRV from N=2)"
[ "${secs:-0}" -gt "$MAXSECS" ] && echo "ANOMALY: secs=$secs (> $MAXSECS; possible hang)"
if [ "$MIX" = 1 ]; then
  [ "${mlomix:-0}" != 1 ] && echo "ANOMALY: mlomix=$mlomix (expect 1)"
  x0=$(grep -h "x_0 from __SigmMLO.q.prev" $L/lqpe 2>/dev/null | grep -oE "= [TF]" | tr -d '= ')
  out=$(grep -h "SigmMLO_out" $L/lqpe 2>/dev/null | grep -oE "[0-9]+\.[0-9]+$")
  mix=$(grep -h "SigmMLO_mixed" $L/lqpe 2>/dev/null | grep -oE "[0-9]+\.[0-9]+$")
  echo "mix   : x_0 from .prev=$x0  |S_out|=$out  |S_mixed|=$mix"
  [ "$N" -ge 2 ] && [ "$x0" != T ] && echo "ANOMALY: x_0 not read from __SigmMLO.q.prev at N=$N"
fi
prom=$(grep -c "promoted ZmloNew" $L/lmlo_sigr 2>/dev/null)
echo "slots : promoted=$prom"
[ "${prom:-0}" -lt 1 ] && echo "ANOMALY: chi~ slot not promoted in lmlo_sigr"
miss=$(cat llmf $L/llmfgw01 2>/dev/null | grep -c "outside ZmloSig\|MISS=[1-9]")
[ "$miss" -gt 0 ] && echo "ANOMALY: $miss chi~ misses"
nan=$(cat $L/lqpe llmf $L/lmlo_sigr 2>/dev/null | grep -ciE "\bnan\b|infinity")
[ "$nan" -gt 0 ] && echo "ANOMALY: $nan NaN/Inf lines in lqpe/llmf/lmlo_sigr"
lastit=$(grep -oE "^ +it +[0-9]+ +of *[0-9]+" $L/llmf | tail -1)
dq=$(grep -oE "RMS DQ= *[0-9.E+-]+" $L/llmf | tail -1)
echo "scf   : $lastit   $dq"
echo "$lastit" | grep -qE "it +300 +of" && echo "ANOMALY: lmf SCF hit nit=300 (not converged)"
if [ -f "$S/bnd001.spin1" ]; then
  echo "band  : $S/bnd001.spin1 ok  mloON(band)=$(grep -c 'MLO Sigma interpolation ON' $S/llmf_band 2>/dev/null)"
else
  echo "ANOMALY: no $S/bnd001.spin1"
fi
