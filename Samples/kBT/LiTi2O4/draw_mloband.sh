#!/bin/bash
# Draw THE MLO band of a chain snapshot: the MLO model of the Hamiltonian the SCF solved.
#   usage: draw_mloband.sh <snapshot dir> <work dir> <out .dat>
#
# Why not the other two routes (2026-09-26, liti_mlo_v9 iter 3):
#  - job_mlo --mlofreeze adds QMLO_SigRs (iteration-N chi~) to the FROZEN HamRsMLO (step-0d LDA H in
#    the step-0d chi~): two bases, 0.3-0.5 eV off.
#  - lmf --band --mlo looks chi~ up by the literal q; path points are symmetry-equivalent to the
#    stored irreducible q but not equal, so chi~ is rebuilt from the Sigma-free H and Sigma is read
#    in the wrong basis.  Valid only at Gamma.
# This route:
#  1) lmf --writeham --noinv --mlo: H = H_LDA[rho_N] + senex(Sigma^MLO_N) on the SCF mesh, with chi~
#     taken from QMLO_z.  ECALJ_MLO_ALLOW_REBUILD is NOT set, so a single miss aborts.
#  2) QMLO_SigRs set aside (the H from step 1 already contains Sigma), HamRsMLO removed,
#     mlo --mlo WITHOUT --mlofreeze: build the MLO model of that H and solve it on qplist.dat.
# Needs the zMLO fix (c0e17f444): mlo with fewer ranks than q points overran zMLO before it.
S=$(realpath "$1"); W=$(realpath -m "$2"); OUT=$(realpath -m "$3"); T=liti2o4   # absolute: we cd below
B=${MLOBAND_BUILD:-/home/takao/ecalj_dev/SRC/build_nvfortran}   # a build that knows the QMLO_* files (4fb6df3ba or later);
                                                              # kt1's development tree (2026-09-28; was the production ~/ecalj)
MPI=/opt/nvidia/hpc_sdk/Linux_x86_64/26.1/comm_libs/13.1/hpcx/latest/ompi/bin/mpirun
unset ECALJ_MLO_ALLOW_REBUILD ECALJ_MLO_NOSIG
export CUDA_VISIBLE_DEVICES=
mkdir -p "$(dirname "$OUT")"; rm -f "$OUT"
rm -rf "$W"; mkdir -p "$W"
for f in ctrlg.$T.toml rst.$T __atm.$T efermi.lmf syml.$T qplist.dat sigm sigm.$T QMLO_SigRs QMLO_z HamRsMLO; do
  cp -p "$S/$f" "$W/" || { echo "FAIL: $S/$f missing"; exit 1; }
done
cd "$W"
$MPI -np 8 $B/lmf $T --writeham --mkprocar --noinv --mlo > lwriteham 2>&1 || { echo "FAIL writeham: $(tail -2 lwriteham | tr '\n' ' ')"; exit 1; }
on=$(grep -c 'MLO Sigma interpolation ON' lwriteham); miss=$(grep -c 'chi~ MISS' lwriteham)
[ "$on" -gt 0 ] && [ "$miss" -eq 0 ] || { echo "FAIL: mloON=$on miss=$miss"; exit 1; }
mv QMLO_SigRs QMLO_SigRs.hidden; rm -f HamRsMLO
$MPI -np 8 $B/mlo $T --mlo > lmlo 2>&1 || { echo "FAIL mlo: $(tail -2 lmlo | tr '\n' ' ')"; exit 1; }
[ -f band_MLO_spin1.dat ] || { echo "FAIL: no band_MLO_spin1.dat"; exit 1; }
# convert to the bnd001.spin1 layout mlo_rows.py reads: "i x E-EF[eV]" per point, blank line per band
python3 - "$OUT" <<'PY' || { echo "FAIL: conversion"; exit 1; }
import sys, numpy as np
ef=float(open('efermi.lmf').read().split()[0].replace('D','E'))
m=np.loadtxt('band_MLO_spin1.dat'); nb=int(m[:,3].max()); nk=len(m)//nb
x=m[:,0].reshape(nk,nb)[:,0]; e=13.605698*(m[:,1].reshape(nk,nb)-ef)
with open(sys.argv[1],'w') as f:
    f.write(f'# MLO band (MLO model of the SCF H), EF={ef} Ry from efermi.lmf\n')
    for b in range(nb):
        for i in range(nk): f.write(f'{i+1:5d} {x[i]:.8f} {e[i,b]:.10f}\n')
        f.write('\n')
PY
[ -s "$OUT" ] || { echo "FAIL: $OUT empty"; exit 1; }
echo "ok: mloON=$on miss=$miss nbands=$(grep -c '^$' "$OUT") -> $OUT"
