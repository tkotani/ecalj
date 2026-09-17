# Si666gwsc — MLO bands on top of converged QSGW (Si, 6x6x6)

Bands and MLO (muffin-tin localized orbitals) for Si generated on top of
a **converged QSGW self-energy**. Unlike TestInstall/si_gwsc (which runs
the gwsc cycle itself), this sample ships the converged state:

- `sigm.si` — converged QSGW static self-energy (6x6x6, read by lmf
  automatically since `[ham] rdsig` is on)
- `rst.si`  — matching charge density restart

so the test only runs `job_band` + `job_mlo` (fast) and checks the
resulting `bnd00*.spin1` and `band_MLO_spin1.dat` against references.

## Run / test

From `Samples/MLOsamples/`:

    testecalj -np 12 Si666gwsc

To regenerate `sigm.si`/`rst.si` from scratch: run `gwsc` with
`ctrlg.si.toml` (n1n2n3=[6,6,6]) to self-consistency,
then copy sigm.si/rst.si here.
