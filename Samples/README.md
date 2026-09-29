# ecalj/Samples — tests and example calculations

Every sample reads one input file, `ctrlg.<sname>.toml`. A sample directory that holds a `test.py` is a test target:
`testecalj <target> -np 8`, run in the directory that contains the target, copies it to `<target>_work`, runs it and
compares the results with the reference files of the sample. `TOOLS/samples_tests.sh` runs the groups of table 1 one
after another and writes one summary (how to read it: ecaljdoc
[ForDevelopers](https://ecalj.github.io/ecaljdoc/manual/ForDevelopers), section 5).

**Table 1**. Sample directories

| directory | what it shows | targets | group of `samples_tests.sh` |
|---|---|---|---|
| [GetStarted/](GetStarted/README.md) | the seed of the [tutorial](https://ecalj.github.io/ecaljdoc/manual/README_tutorial): GaAs from `ctrls.gaas` | — | `inputs` |
| [TestInstall/](TestInstall/README_testecalj.md) | installation tests: LDA (`c`, `co`, `cu`, `fe`, `gdn` with LDA+U, `te` with relaxation, `felz`, `gasls` with SOC, ...), GW and QSGW (`si_gwsc`, `gas_gwsc`, `nio_gwsc`, `fe_gwsc`, `fe_kbt` at 3000 K, `si_gw_lmfh`, ...), dielectric functions, spin susceptibility, cRPA | 26 in `--all` (64 checks); 4 heavy GW targets run by name | `install`, `gwall`, `heavy` |
| [EPS/](EPS/) | dielectric function, `job_eps`: Cu, Ag, GaAs | 3 | `eps` |
| [PROCAR/](PROCAR/README.md) | fat bands and orbital weights: MgO, Ni2MnGa | 2 | `procar` |
| [MLOsamples/](MLOsamples/README.md) | MLO (muffin-tin localized orbitals, the model Hamiltonian that replaces Wannier functions): semiconductors, metals, SOC, 4f, screened W by `job_mloW` | 25 | `mlo` |
| [MLOQSGW/](MLOQSGW/README.md) | QSGW with the self-energy interpolated in the MLO representation (`gwsc --mlo`): GaAs, NiO | 2 | `mloqsgw` |
| [AFsymmetry/](AFsymmetry/README.md) | antiferromagnetic symmetry (`symgrpaf`): NiO and NiSe in LDA, NiO in QSGW | 3 | `afsym` |
| [Magnon/](Magnon/README.md) | magnon spectra through Wannier functions (`job_magnon`): Fe, Ni, FeCo, Fe in the simple cubic cell | 4 | `magnon` |
| [kBT/](kBT/README.md) | electron temperature and smearing in GW (`t_tetrakbt`, `t_sigmaw`): temperature scans of Si and Fe; MLO-QSGW of LiTi2O4; the research log | `scanT/Si` | `samples` |
| [FermiSurface/](FermiSurface/Cu/README.md) | Fermi surface of Cu for xcrysden (`job_fermisurface`) | `Cu` | `samples` |
| [Doping/](Doping/Si/README.md) | doping of Si: fractional nuclear charge, background charge `zbak`, fixed spin moment | `Si` | `samples` |
| [HomoGas/](HomoGas/es/README.md) | Lindhard function of the electron gas by the tetrahedron method, against the analytic form | `es` | `samples` |
| [BoltzTraP/](BoltzTraP/Si/README.md) | input files of BoltzTraP written by `lmf --boltztrap` | `Si` | `samples` |
| [SLAB/](SLAB/Cu001/README.md) | Cu(001) slab of four layers: bands, DOS, work function | `Cu001` | `samples` |
| [SOC/](SOC/FePt_MAE/README.md) | magnetic anisotropy energy by the force theorem (SOC for two spin axes on one potential): FePt, MnGa | `FePt_MAE`, `MnGa_MAE` | `samples` |
| [LDAU/](LDAU/ReN/README.md) | LDA+U on the 4f shell with spin-orbit coupling (`so = 2`) and starting occupations `occnum.<sname>`: GdN, PrN | `ReN` | `samples` |
| [EffectiveMass/](EffectiveMass/GaAs/README.md) | effective masses of GaAs from the QSGW bands with SOC (the mass mode of `syml.<sname>`, fit by `massfit.py`) | `GaAs` | `samples` |
| [Relax/](Relax/LaGaO3/README.md) | relaxation of the atomic positions: LaGaO3 | `LaGaO3` | `samples` |
| [BenchmarkTest/](BenchmarkTest/) | one QSGW iteration of InAs/GaSb superlattices (16 and 32 atoms), for GPU | 2 | `bench` |
| [mptf32problem/](mptf32problem/README.md) | AgNO3 and other cases where the old TF32 mode of the GPU build failed (GW1500) | — | — |
| [Legacy/](Legacy/) | older examples that are not rebuilt (table 2) | — | `inputs` |

The group `inputs` takes every `ctrlg.<sname>.toml` of this tree: `lmchk` must read it, and the rules of
`ctrlg_update.py` must leave it unchanged.

## Running

```bash
cd Samples/TestInstall && testecalj --all -np 8           # the installation tests (4 minutes)
cd Samples/AFsymmetry  && testecalj NiO NiO_gwsc -np 8    # targets by name
TOOLS/samples_tests.sh                                    # every group (hours: bench, heavy and magnon are long)
TOOLS/samples_tests.sh --gpu inputs install samples       # some groups, GPU build
```

Notes:

- The magnon targets make work directories of 1-5 GB each.
- `BenchmarkTest` with one GPU of 32 GB: `-np2 1` (`samples_tests.sh --gpu` does it).
- `TestInstall/yh3fcc_gwsc666` does not run with the present code (the radial solver fails on the extended local
  orbitals of Y); it is kept with its old references and is not in any group.
- An input of your own in the legacy form (`ctrl.<sname>` + `GWinput`) is converted by `Legacy2toml.py <sname>`
  (ecaljdoc [toml_migration](https://ecalj.github.io/ecaljdoc/manual/toml_migration)).

## Legacy

**Table 2**. Directories under `Legacy/`. Their inputs are converted to `ctrlg.<sname>.toml` (and are in the group
`inputs`), but the job scripts and READMEs are of the time before 2026-05 and are not tests.

| directory | content | state |
|---|---|---|
| `AHC/` | anomalous Hall conductivity of bcc Fe (`job_AHC`, `hx0ahc.py`) | to be rebuilt; `hx0ahc.py` loses a single extra option |
| `IIR/` | impact ionization rate from Im Sigma (diamond, 12x12x12: 200 minutes) | long |
| `CMDsample/` | tutorial inputs of 2019 (Si, InAs, ZnS, GaN, Fe, NiO, BaTiO3 QSGW) | the commands of its README are gone |
| `MATERIALS/La2CuO4`, `InAsGaSb/`, `Samples_ISSP/`, `bga2o3deformation/` | QSGW inputs and results of larger systems | long |
| `ReNcub/`, `mass_fit_test/`, `LaGaO3_relax/` | the older runs behind `LDAU/`, `EffectiveMass/`, `Relax/` | kept for their results |
| `FermiSurface/`, `Si_doping_sample/`, `TETRAHEDRON_HomoGas/`, `TETRAHEDRON_HomoGas_test/`, `BOLZTRAP/`, `SLAB/`, `SOCAXIS/` | the originals of `FermiSurface/`, `Doping/`, `HomoGas/`, `BoltzTraP/`, `SLAB/`, `SOC/` of table 1 | rebuilt; the originals are kept for now |
| `GdNldau/`, `MATERIALS/erasldau`, `MATERIALS/pdo_gwsc443`, `MATERIALS/yh3fcc_gwsc666`, `SOC/` | covered by `TestInstall` (`gdn`, `eras`, `pdo_gwsc443`, `felz`) and `MLOsamples/FeSoc` | duplicates |
| `TestHomoDimerAtom/`, `UUmatSOC/`, `AHCSOCtest/` | scripts of 2012 in Python 2; work files without input; a copy of the AHC README | not usable |
| `superlattice/` | generator of strained zincblende superlattice POSCARs | a tool, no calculation |
