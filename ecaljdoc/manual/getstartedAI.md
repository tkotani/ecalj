# GetStarted with AI

::: tip How to use this page (for people) / このページの使い方
**1. With an AI chat (no install).** Open an AI assistant that can read web pages (Claude, for example) and paste:

`Read https://ecalj.github.io/ecaljdoc/manual/getstartedAI and guide me through ecalj.`

日本語なら: `https://ecalj.github.io/ecaljdoc/manual/getstartedAI を読んで、ecalj を日本語で案内してください。`

**2. With Claude Code on your computer (to run calculations together; recommended).**

```bash
git clone https://github.com/tkotani/ecalj
cd ecalj
claude        # then say: "Read ecaljdoc/manual/getstartedAI.md and guide me through ecalj"
```

The AI asks what you want to do, explains ecalj briefly, and walks you through installation, a first calculation (GaAs) and the
tests, in your language. ecalj is developed together with Claude Code.

**Note**: many pages of this manual and of the repository mix Japanese and English. Ask your AI to translate them for you.
:::

**For the AI**: the rest of this page is written for you. It tells you what to read, in which order, and how to guide the person.
Many pages here and in the repository are partly in Japanese: translate them into the person's language when you explain.
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

## 6. Exploring the repository without cloning it

Every file of ecalj can be read on GitHub. For an AI that fetches web pages, the raw text is at

```text
https://raw.githubusercontent.com/tkotani/ecalj/main/<path>        e.g. .../main/Changes.md, .../main/Samples/README.md
https://github.com/tkotani/ecalj/tree/main/<directory>             the listing of a directory
```

*Table 3*. The top of the repository

| path | what is there |
| --- | --- |
| `README.md`, `CLAUDE.md`, `Changes.md` | entry (three links and the publishing steps), the entry for an AI working on the code, every change (newest first) |
| `InstallAll.py` | build and install (`--fc` with `gfortran`, `ifx` or `nvfortran`; `--gpu`; `--bindir`) |
| `SRC/subroutines/` | the Fortran library `libecaljF.so` (modules `m_*.f90`; the programs are thin entries) |
| `SRC/main/` | the entries of the Fortran programs (table 4) |
| `SRC/exec/` | the Python and shell drivers (`gwsc`, `job_*`, `ctrlgenToml.py`, `testecalj`, ...) and `pylib/` |
| `Samples/` | worked examples and tests, one directory per topic, each with a README (`Samples/README.md` lists them) |
| `ecaljdoc/` | the source of this site |
| `ecalj_auto/` | automatic runs of many materials (GW1500), the database builder `gw1500db_build.py` |
| `TOOLS/` | maintenance scripts: `samples_tests.sh` (tests by group), `publish_*.sh`, `doclinks.py`, `module_map.py` |
| `MD/` | development notes (not on this site; `MD/README.md` first) |
| `GetSyml/`, `StructureTool/` | k paths for band plots, structure conversion |

## 7. Programs and drivers

*Table 4*. What runs what (installed into the bin directory by `InstallAll.py`)

| command | what it does |
| --- | --- |
| `vasp2ctrl POSCAR` | POSCAR → `ctrls.POSCAR.vasp2ctrl` (structure in ecalj form; copy it to `ctrls.<sname>`) |
| `ctrlgenToml.py <sname>` | `ctrls.<sname>` → the input `ctrlg.<sname>.toml` (basis, k meshes, GW settings; `--ssig=0.8` for QSGW80) |
| `lmfa <sname>` | free atoms (the starting density and the atom files) |
| `lmf <sname>` | the DFT self-consistency (also the QSGW one-body part when `sigm` exists); `lmchk <sname>` checks the structure |
| `getsyml <sname>` | the k path for band plots |
| `job_band <sname> -np N` | band plot (files `bnd*.spin*`, `bandplot.isp1.glt`); `job_tdos`, `job_pdos`: densities of states |
| `gwsc <n> -np N <sname>` | `n` QSGW iterations: lmf → GW inputs → W and Σ (`hgw`) → `hqpe_sc` (QP energies, `sigm`) → lmf |
| `gwscconv -np N <sname> --conv-tol 0.1` | QSGW repeated until the gap (or, for metals, the levels near E_F) stops changing |
| `job_eps <sname>` | dielectric function (RPA) |
| `job_mlo <sname>`, `job_mlo_soc`, `job_mloW <sname> [--crpa]`, `job_mlo_magnon` | MLO model, with SOC, U/J/cRPA, magnons |
| `testecalj [--all] [targets]` | the tests (in `Samples/TestInstall`, or in a sample directory) |
| `ctrlg_update.py`, `Legacy2toml.py` | update an older `ctrlg`, convert the old `ctrl` + `GWinput` |

Every program prints its options with `--help`; [Command-line options](./cmdopts.md) lists them.

## 8. The input file and the output files

`ctrlg.<sname>.toml` has the sections `[struc]` (lattice), `[[site]]` (atoms), `[[spec]]` (species: basis, radii, LDA+U),
`[bz]` (k mesh, metal or not), `[iter]` (mixing), `[ham]` (functional, spin, SOC, `scaledsigma`, plane waves), `[options]`,
`[gw]` (GW meshes and cut-offs, `t_tetrakbt`, `t_sigmaw`), `[mlo]`, and at the end `[product_basis]` (made automatically).
Any key can be overridden on the command line: `lmf <sname> --ctrlg:bz.nkabc=[8,8,8]` ([TOML migration](./toml_migration.md)).

*Table 5*. Output files to look at

| file | what |
| --- | --- |
| `llmf` (stdout of lmf in gwsc), `log.<sname>` | the course of the run; the gap line `VBmax= ... CBmin= ... gap = ... eV`; `efermi.lmf` the Fermi level |
| `rst.<sname>` | the density (restart) |
| `sigm`, `sigm.<sname>` | the QSGW self-energy (static, Hermitian) used by lmf |
| `QPU`, `QPD` | QP energies of the last iteration (spin up, down) |
| `bnd*.spin1` | band lines of `job_band` |
| `dos.tot.<sname>*` | total DOS |
| `QSGW.<n>run/` | the files of each iteration (`gwsc`, `gwscconv`) |

## 9. Words

PMT: the basis, augmented plane waves plus muffin-tin orbitals. QSGW: quasiparticle self-consistent GW; the static
Hermitian part of Σ replaces the exchange-correlation potential and is iterated to self-consistency. QSGW80: 80 % of the QSGW Σ
and 20 % of LDA (`scaledsigma = 0.8`). MLO: localized orbitals made from the muffin-tin orbitals, orthonormalized (Löwdin),
for models. cRPA: the constrained RPA for the screened U. ESM: the effective screening medium for slabs. `<sname>`: the name of
a material, the extension of its files (`ctrlg.gaas.toml` → `gaas`).

## 10. Data, and where to look when something goes wrong

- The GW1500 database: QSGW80 band gaps, bands and DOS of 1546 materials of Materials Project,
  https://github.com/tkotani/DOSnpSupplement/tree/main/QSGW80_2026 (its README gives the conditions of each value and the checks).
  The 2025 tables there are the supplement of arXiv:2507.19189.
- An old `ctrl.<sname>` or `GWinput`: not read any more; convert with `Legacy2toml.py`. An old `ctrlg` with `SmearX0` stops;
  `ctrlg_update.py` rewrites it. `t_tetrakbt` must be in `[gw]`.
- Numbers off from a README: check the k mesh, `scaledsigma`, and the number of iterations first; then `Changes.md` for a change
  of defaults since the README was written.
- Questions and reports: the GitHub repository, https://github.com/tkotani/ecalj
