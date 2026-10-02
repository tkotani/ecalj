# GetStarted with AI

**For people**: give the address of this page (or of the ecalj repository, https://github.com/tkotani/ecalj) to an AI assistant
— Claude Code is what ecalj is developed with — and ask, for example, *"Read this and guide me through ecalj"*. The AI then
explains what ecalj is, asks what you want to do, and walks you through installation, a first calculation and the tests,
in your language.

**For the AI**: the rest of this page is written for you. It tells you what to read, in which order, and how to guide the person.
(Written 2026-10-03. When this page and the code disagree, the code and `Changes.md` of the repository are right.)

## 1. ecalj in five lines

1. A first-principles electronic-structure package: DFT (LDA/GGA, LDA+U) with the PMT basis (augmented plane waves + muffin-tin orbitals), and
   **QSGW** (quasiparticle self-consistent GW), which gives band structures, band gaps and effective masses close to experiment
   without adjustable parameters.
2. On top of QSGW: dielectric functions, spin fluctuations and magnons, self-energy spectra, impact ionization, and model
   Hamiltonians by MLO (muffin-tin-orbital-based localized orbitals, with U, J and cRPA).
3. Input: one TOML file per material, `ctrlg.<sname>.toml`, generated from a POSCAR (VASP format) by `vasp2ctrl` and `ctrlgenToml.py`.
   Most settings are automatic.
4. Runs on CPUs (MPI) and NVIDIA GPUs (`--gpu`; the precision of the GW products is chosen by `--prec=tf32|fp32|fp64`).
5. Open source (AGPLv3), developed by T. Kotani and co-workers, largely together with Claude Code.

## 2. What to read, and in which order

Read only as deep as the person's question needs. Pages of this site are given by their names; open them as `./<name>.md`
next to this page.

**Step 1, this site (the manual)**

| order | page | why |
| --- | --- | --- |
| 1 | [What's new](./whatsnew.md) | the recent changes and the big topics; what is current |
| 2 | [MainDocument](./README_tutorial.md) | what ecalj can do, the overview of QSGW, and **GetStarted** (GaAs, step by step: POSCAR → ctrlg → LDA → band plot → QSGW) |
| 3 | [Install](../install/install.md) | `InstallAll.py` and the install tests |
| 4 | [Samples](./samples.md) | which sample shows what (each sample directory has a README) |
| 5 | [TOML migration](./toml_migration.md) | the input file, the override `--ctrlg:<path>=<value>`, the change from the old `ctrl` + `GWinput` |
| 6 | by topic | [lmf](./lmf.md) (DFT, symmetry), [gwsc](./gwsc.md) (QSGW), [GPU version](./ecaljgpu.md), [kBT](./kBT.md) (temperature and smearing), [MLO](./mlo.md) (models, U/J/cRPA, magnons), [MLO-gwsc](./mlo_gwsc.md), [Dielectric function](./optical.md), [spectrum](./spectrum.md), [UsageDetailed](./UsageDetailed.md) (antiferromagnets, SOC, LDA+U), [Command-line options](./cmdopts.md) |

**Step 2, the repository** (when the person has cloned it, or through GitHub)

- `README.md` (three links), then `CLAUDE.md`: the entry for an AI working on the code. It loads the rules and the map of the
  documents. `Changes.md`: every change, newest first (the canonical record; [What's new](./whatsnew.md) is its summary).
- `Samples/<dir>/README.md`: what each sample is, how to run it, what its test compares.
- The development notes are in the directory `MD/` of the repository (not on this site). Start from `MD/README.md`: the reading order,
  then jumps by topic (finite temperature, MLO, magnons, GPU, tests). Read them only when the person develops ecalj or asks why
  something is as it is.
- The structure of the code (main programs → modules): `MD/module_map.md`.

## 3. How to guide a person

Explain in the person's language. One step at a time; show the commands; let the person run them and report back.

1. **Ask first** (three short questions):
   - What do you want to do? (try ecalj / band gap or band structure of a material / dielectric function / magnons / a model
     Hamiltonian / finite temperature / large-scale automatic runs / develop the code)
   - Your background: DFT? GW? (decides how much theory to explain)
   - Your computer: OS, number of cores, NVIDIA GPU or not, compilers (gfortran, ifx, nvfortran). Is ecalj already installed?
2. **Give the overview** in a few lines (section 1), adjusted to the answers. For a newcomer, say what QSGW buys over DFT (gaps without
   adjustable parameters) and what it costs (minutes to hours instead of seconds).
3. **Route by the goal** (table 1), then do it together.
4. **Check the results** with the person: the gap printed by lmf (`gap = ... eV` in `llmf`), the band plot files, the iteration
   history of QSGW. Compare with the README of the sample or with known values.
5. **When something fails**: read the end of the log (`llmf`, `lgwsc`, the `l*` files of the step that failed), compare with the
   README of the nearest sample, and look the error up in `Changes.md` of the repository.

*Table 1*. Goal → where to start

| goal | start with | then |
| --- | --- | --- |
| try ecalj | [Install](../install/install.md), then [GetStarted](./README_tutorial.md) with GaAs | QSGW of GaAs with `gwsc` ([gwsc](./gwsc.md)) |
| band gap / bands of a material | GetStarted with the person's POSCAR | `gwsc` for QSGW; `--ssig=0.8` in `ctrlgenToml.py` for QSGW80 |
| GPU | [GPU version](./ecaljgpu.md) | `--prec=tf32` is fast and accurate enough for gaps |
| dielectric function | [Dielectric function](./optical.md), samples `Samples/EPS` | |
| magnons, U/J, cRPA, models | [MLO](./mlo.md) §6, sample `Samples/Magnon/Fe_mlo_magnon` | |
| finite temperature | [kBT](./kBT.md): two keys only, `t_tetrakbt` and `t_sigmaw` | sample `Samples/kBT/scanT` |
| antiferromagnets, SOC, LDA+U | [UsageDetailed](./UsageDetailed.md), [lmf](./lmf.md) | samples `Samples/AFsymmetry`, `Samples/SOC`, `Samples/LDAU` |
| many materials automatically | `ecalj_auto/` of the repository (`README.md` there) | |
| develop the code | `CLAUDE.md` and `MD/README.md` of the repository | |

## 4. Tests

After `InstallAll.py` (it runs the minimum tests), the full tests are run by `testecalj` and, group by group, by
`TOOLS/samples_tests.sh` of the repository:

```bash
cd ecalj/Samples/TestInstall && testecalj -np 8 --all          # the install tests (CPU)
ecalj/TOOLS/samples_tests.sh -np 8                               # every group of the samples (CPU)
ecalj/TOOLS/samples_tests.sh --gpu -np 8 install mlo heavy       # some groups, GPU
```

Each group prints one line `<group> PASSED <n> checks passed, 0 failed <seconds> s`; the summary is `summary.txt` in the log
directory. Count the checks from the last summary only (the log repeats the earlier lines).

*Table 2*. Time of each test group on the machines of the developers (2026-10-03, from their `summary.txt`: the fastest passed run with
the standard options, `-np 8`; t14 and mic on CPU, kt1 and kr7 with `--gpu`. The records of mic, kt1 and kr7 are of 2026-09-30, when
some groups had other targets, e.g. magnon with four samples instead of one. ✗: no passed run, the references of bench are old.)

| test group | t14 gfortran (CPU, 16 cores) | mic ifx (CPU, 32 cores) | kt1 (RTX 5090 ×2) | kr7 (RTX 5090) |
| --- | --- | --- | --- | --- |
| inputs | 2 min | 56 s | 58 s | 56 s |
| install (TestInstall, all) | 6 min | 4 min | 4 min | 4 min |
| eps | 2 min | 15 s | 53 s | 12 s |
| procar | 47 s | 19 s | 22 s | 20 s |
| mlo | 11 min | 5 min | 4 min | 4 min |
| mloqsgw | 2 min | 54 s | 54 s | 47 s |
| afsym | 55 s | 12 s | 9 s | 7 s |
| affix | 2 min | — | — | — |
| FermiSurface | 9 s | 6 s | 4 s | 5 s |
| Doping | 22 s | 14 s | 9 s | 10 s |
| HomoGas | 5 s | 3 s | 3 s | 3 s |
| BoltzTraP | 6 s | 4 s | 3 s | 3 s |
| SLAB | 32 s | 18 s | 11 s | 12 s |
| SOC | 4 min | 2 min | 1 min | 1 min |
| LDAU | 2 min | 1 min | 40 s | 42 s |
| EffectiveMass | 2 min | 42 s | 33 s | 26 s |
| Relax | 4 min | 2 min | 2 min | 2 min |
| DOS | 57 s | 32 s | 22 s | 21 s |
| IIR | 45 s | 35 s | 21 s | 19 s |
| AtomDimer | 27 min | 11 min | 7 min | 7 min |
| kBT_scanT | 28 s | 15 s | 27 s | 16 s |
| magnon | 2 min | 11 min | 8 min | 8 min |
| heavy (four large GW tests) | — | 25 min | 7 min | 14 min |
| bench (BenchmarkTest) | — | 25.9 h ✗ | — | 2.8 h ✗ |
| **total without bench and heavy** | **71 min** | **41 min** | **32 min** | **31 min** |

So: everything except bench and heavy takes 30–70 minutes; bench adds about 3 hours on a GPU. On a machine busy with other jobs
the tests slow down a lot (on kt1, install took over 2 hours next to production GW runs).

## 5. Points to get right when you explain

- QSGW gives band structures and gaps (eigenvalues and wave functions), **not total energies**.
- The input is `ctrlg.<sname>.toml` (since 2026-05). The old `ctrl.<sname>` + `GWinput` are not read; `Legacy2toml.py` converts them.
- Temperature and smearing are set by **two keys only**, `t_tetrakbt` (χ0) and `t_sigmaw` (Σ) ([kBT](./kBT.md)).
- `ctrlgenToml.py` writes `scaledsigma = 1.0` (full QSGW) by default; QSGW80 (`--ssig=0.8`) mixes 80 % of the QSGW self-energy and
  corrects the average overestimate of QSGW gaps.
- The gap printed by lmf is taken on the k mesh; when the band extremum is off the mesh (layered or hexagonal materials, K points),
  the gap along the band path is the better estimate.
- Symmetry is found by spglib (`symgrp = "find"`, the default); generators in `symgrp` lower the symmetry.
- MLO are orthonormalized (Löwdin) localized orbitals; models, U/J, cRPA and magnons use them ([MLO](./mlo.md) §6).
- GPU precision: `--prec=tf32` changes gaps by about 1 meV against fp32; the old option `--mp` (all products in TF32, before 2026-09)
  could shift gaps by 0.1–0.3 eV.
