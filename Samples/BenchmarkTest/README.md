# BenchmarkTest — one QSGW iteration of InAs/GaSb superlattices

Two targets for timing and for a check of the GPU build: `inas2gasb2` (8 atoms) and `inas4gasb4` (16 atoms),
`gwsc 0` from LDA. The test compares the `QPU` of the run with the reference `QPU.1run` (`dqpu`, 0.011 eV).

```bash
cd Samples/BenchmarkTest
testecalj inas2gasb2 inas4gasb4 -np 8 --gpu -np2 1      # one GPU of 32 GB: -np2 1
```

**Table 1**. Times with one RTX 5090 and 8 CPU ranks (kr7, nvfortran 26.1, double precision, 2026-09-30)

| target | atoms | time |
|---|---|---|
| `inas2gasb2` | 8 | 56 min |
| `inas4gasb4` | 16 | 2 h 8 min |

On CPU alone the two targets take about a day (8 ranks): `TOOLS/samples_tests.sh` runs the group `bench` last.

The references were made on 2026-09-30 (kr7, GPU, double precision). Those of 2026-05-31 differed from the present
results by 0.03 eV in SEx and SEc of the band-edge states (they cancel) and by 0.01-0.02 eV in the quasiparticle
energy shift: the systems have a gap of 0.14 eV in LDA, and the occupation of the band-edge states follows the
Fermi-Dirac kernel of `t_sigmaw` since 2026-09-20 (ecaljdoc [kBT](https://ecalj.github.io/ecaljdoc/manual/kBT), section 3).

`ISSP/` holds job scripts for the ISSP supercomputers.
