# fe — bcc Fe, spin-polarized LDA + total DOS

Ferromagnetic bcc Fe (1 atom/cell, nspin=2, tetrahedron METAL=3,
nkabc=10). Self-consistency runs from scratch (`lmfa` start) to
conv=1e-4 Ry; the converged moment is mmom ~ -2.30 mu_B (sign follows
the MMOM=-2.2 initial condition). Total DOS is produced with `job_tdos`.

## Run / test

From `Samples/TestInstall/`:

    testecalj -np 12 fe

or by hand in a work copy of this directory:

    lmfa fe
    mpirun -np 12 lmf fe
    job_tdos fe -np 12        # -> dos.tot.fe (+ .glt/.pdf via gnuplot)

Checked quantities: `out.lmf.fe` (etot / ehf / ehk / mmom / RMS DQ via
test1_check) and `dos.tot.fe` (abs tol 1e-3).

## Notes

- `ctrl.fe` is the legacy input this sample was converted from
  (`Legacy2toml.py fe` regenerates `ctrlg.fe.toml`). The legacy CONST
  variables (a=2.87 AA, R=0.99*sqrt(3)a/4) are baked into the TOML.
- `CLSinput` is kept for the core-level spectroscopy demo: `lmf fe --cls`
  runs with it (see `--cls` in `ecalj_cmdopts.list`); the legacy
  `lmdos --cls` post-processing flow is not part of the test.
- k-mesh can be overridden at run time: `--ctrlg:bz.nkabc=[6,6,6]`.
