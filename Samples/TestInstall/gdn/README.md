# gdn — GdN (rocksalt), LDA+U on Gd 4f + total DOS

GdN in the rocksalt structure (Gd + N, nspin=2) with LDA+U on the Gd 4f
shell: `idu=[0,0,0,2]`, `uh=[0,0,0,0.515]` Ry (around-mean-field type 2).
The spin moment is mmom = 7 (Gd f^7). By design the run is a short fixed
2-iteration demo (`nit=2` inherited from the legacy ctrl) — it exercises
the LDA+U machinery (dmats, Etot(LDA+U) term) without aiming at full
self-consistency. Total DOS via `job_tdos`.

## Run / test

From `Samples/TestInstall/`:

    testecalj -np 12 gdn

or by hand in a work copy of this directory:

    lmfa gdn
    mpirun -np 12 lmf gdn
    job_tdos gdn -np 12       # -> dos.tot.gdn (+ .glt/.pdf via gnuplot)

Checked quantities: `out.lmf.gdn` (etot / ehf / ehk / mmom / Etot(LDA+U)
via test1_check) and `dos.tot.gdn` (abs tol 1e-3).

## Notes

- `ctrl.gdn` is the legacy input this sample was converted from
  (`Legacy2toml.py gdn`). The conversion resolves the legacy CONST
  variables (ALAT=a*fv) and takes the `%else` branch of `%ifdef coref`
  with `%ifdef ldau` enabled — i.e. the LDA+U variant without a 4f
  core hole.
- To converge for real work, raise the iteration count:
  `--ctrlg:iter.nit=50`.
