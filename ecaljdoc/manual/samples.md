# Samples — what's where in `ecalj/Samples/`

This page is a thin index over the [`ecalj/Samples/`](https://github.com/tkotani/ecalj/tree/main/Samples)
tree.  Each row points at the per-directory README that owns the
authoritative description; the role / what-to-start-with column is
meant to help you pick the right entry point without opening every
sub-tree.

## Top-level layout (since 2026-05)

```
ecalj/Samples/
├── README.md          (https://github.com/tkotani/ecalj/blob/main/Samples/README.md)
├── GetStarted/        minimal seeds for the README_tutorial walk-through
├── TestInstall/       testecalj-runnable (install validation)
├── EPS/               testecalj-runnable
├── PROCAR/            testecalj-runnable
├── MLOsamples/        testecalj-runnable
├── MLOQSGW/           testecalj-runnable (QSGW with gwsc --mlo)
├── AFsymmetry/        testecalj-runnable (antiferromagnetic symmetry, LDA and QSGW)
├── Magnon/            testecalj-runnable (magnon spectra, job_mlo_magnon)
├── kBT/               electron temperature in GW: scanT/ (testecalj-runnable), LiTi2O4, research log
├── FermiSurface/ Doping/ HomoGas/ BoltzTraP/ SLAB/ SOC/ LDAU/ EffectiveMass/ Relax/ DOS/ IIR/
│                      short examples rebuilt from Legacy/, each testecalj-runnable
├── MATERIALS/         LDA and MLO samples of 65 materials; inputs and results of QSGW of larger systems: not tests
├── BenchmarkTest/     testecalj-runnable (QSGW of InAs/GaSb superlattices, for GPU)
├── mptf32problem/     GW1500 failure corpus: TF32 vs FP32 on AgNO3 etc.
└── AtomDimer/         a molecule and an atom in a box (N2): bond length and binding energy
```

The inputs of all these directories are `ctrlg.<sname>.toml`, what `lmf` /
`gwsc` / `job_eps` / ... read out-of-the-box.  The groups of test targets
are run one after another by `TOOLS/samples_tests.sh` of ecalj.
The targets of every directory and the group of `samples_tests.sh` that runs
them are in table 1 of ecalj
[`Samples/README.md`](https://github.com/tkotani/ecalj/blob/main/Samples/README.md).
For your own directory with the legacy `ctrl.<sname>` + `GWinput`, see
[TOML migration](./toml_migration) for the conversion flow.

## TOML-migrated samples (testecalj-runnable)

| dir | physics / role | per-dir README | doc page |
|---|---|---|---|
| **GetStarted/GaAs** | minimal worked example for the [tutorial](./README_tutorial#getstarted) — GaAs zinc-blende, 2 atoms, ships `ctrls.gaas` + `ctrlg.gaas.toml` | [GetStarted/README.md](https://github.com/tkotani/ecalj/blob/main/Samples/GetStarted/README.md) / [GaAs README](https://github.com/tkotani/ecalj/blob/main/Samples/GetStarted/GaAs/README.md) | [./README_tutorial](./README_tutorial) |
| **EPS/EPS_Cu** | dielectric ε(q,ω), `job_eps`, FCC metal | [EPS_Cu/](https://github.com/tkotani/ecalj/tree/main/Samples/EPS/EPS_Cu) (no README; see test.py) | [Dielectric function: job_eps](./optical) |
| **EPS/EPS_GaAs** | ε(q,ω), GaAs zinc-blende semiconductor | [EPS_GaAs/](https://github.com/tkotani/ecalj/tree/main/Samples/EPS/EPS_GaAs) (no README; see test.py) | [./optical](./optical) |
| **EPS/EPS_Ag** | ε(q,ω), Ag (FCC, no LFC) | [EPS_Ag/](https://github.com/tkotani/ecalj/tree/main/Samples/EPS/EPS_Ag) (no README; see test.py) | [./optical](./optical) |
| **MLOsamples/** (25 test targets) | MTO Localized Orbitals (Wannier replacement); SOC variants; on-site W via `job_mloW` | [MLOsamples/README.md](https://github.com/tkotani/ecalj/blob/main/Samples/MLOsamples/README.md) | README_SOC.md (`MD/mlo_notes/mlo_soc_memo.md`) |
| **PROCAR/MgO_PROCAR** | fat-band weight; O-2p projection on MgO bands | [MgO_PROCAR/](https://github.com/tkotani/ecalj/tree/main/Samples/PROCAR/MgO_PROCAR) (no README; see test.py) | [./UsageDetailed § PROCAR mode](./UsageDetailed#procar-mode) |
| **PROCAR/Ni2MnGa_L21_PROCAR** | per-atom fat band; FM Heusler; ships converged `rst.ni2mnga` | [Ni2MnGa_L21_PROCAR/](https://github.com/tkotani/ecalj/tree/main/Samples/PROCAR/Ni2MnGa_L21_PROCAR) (no README; see test.py) | [./UsageDetailed § PROCAR mode](./UsageDetailed#procar-mode) |
| **TestInstall/** (26 targets in `testecalj --all`, 66 checks) | install validation: ground-state, GW (`gwsc`, `gw_lmfh`), finite-T (`fe_kbt`), eps (`job_eps`), ChiPM, cRPA of the MLO model. The heavy GW targets `cugase2_gwsc222`, `nio_gwsc444`, `pdo_gwsc443`, `gas_gwsc666` are run one by one | (per dir; driven by `testecalj --all`) | [./gwsc](./gwsc), [./optical](./optical) |
| **MLOQSGW/** (`GaAs`, `NiO`) | QSGW with the self-energy interpolated in the MLO representation (`gwsc --mlo`) | [MLOQSGW/README.md](https://github.com/tkotani/ecalj/blob/main/Samples/MLOQSGW/README.md) | [./mlo_gwsc](./mlo_gwsc) |
| **BenchmarkTest/** (`inas2gasb2`, `inas4gasb4`) | one QSGW iteration of InAs/GaSb superlattices against the stored `QPU.1run`; with one 32 GB GPU use `-np2 1` | (no README; see test.py) | [./ecaljgpu](./ecaljgpu) |
| **AFsymmetry/** (`NiO`, `NiSe`, `NiO_gwsc`) | antiferromagnetic symmetry `symgrpaf` with `af = 1 / -1` on the sites: LDA, and QSGW (the GW programs compute both spins) | [AFsymmetry/README.md](https://github.com/tkotani/ecalj/blob/main/Samples/AFsymmetry/README.md) | [./lmf](./lmf) |
| **Magnon/** (`Fe_mlo_magnon`) | magnon spectra through the MLO model (`job_mlo_magnon`, bcc Fe), compared with the Wannier functions of the version before 2026-10-02 (tag `last-wannier`) | [Magnon/README.md](https://github.com/tkotani/ecalj/blob/main/Samples/Magnon/README.md) | [./UsageDetailed](./UsageDetailed) |
| **kBT/scanT/** (`Si`) | the same input at several electron temperatures (`scan_T.sh`, `collect_T.py`): gap of Si, bands and moment of bcc Fe | [kBT/scanT/README.md](https://github.com/tkotani/ecalj/blob/main/Samples/kBT/scanT/README.md) | [./kBT](./kBT) section 5 |
| **FermiSurface/Cu**, **Doping/Si**, **HomoGas/es**, **BoltzTraP/Si**, **SLAB/Cu001**, **SOC/** (`FePt_MAE`, `MnGa_MAE`), **LDAU/ReN**, **EffectiveMass/** (`GaAs`, `CdS`, `GaN`), **Relax/LaGaO3**, **DOS/** (`ZnS`, `Fe`), **IIR/C**, **AtomDimer/N2** | short examples (each below a few minutes): Fermi surface for xcrysden, doping by fractional Z or background charge, Lindhard function against the analytic form, input of BoltzTraP, slab with work function, magnetic anisotropy energy, LDA+U, effective masses with SOC, relaxation of atomic positions, total and partial DOS, impact ionization rate (Im of the self-energy), a molecule in a box | README.md in each directory | [./UsageDetailed](./UsageDetailed) |

## Not testecalj targets

| dir | physics / role | per-dir README | doc page |
|---|---|---|---|
| **MATERIALS/** | LDA and MLO samples of 65 materials (62 simple materials, La2CuO4, InAs/GaSb, BaTiO3; how to choose the MLO model: [mlo](./mlo) §9). Also QSGW of the larger systems: La2CuO4 (20 iterations, 7 atoms), InAs/GaSb superlattices of 16 and 40 atoms, BaTiO3 (three stages), inputs and results kept, long runs | [MATERIALS/README.md](https://github.com/tkotani/ecalj/blob/main/Samples/MATERIALS/README.md) | [./gwsc](./gwsc) |
| **kBT/LiTi2O4** and the rest of kBT/ | LiTi₂O₄ (metallic spinel): the input of the MLO-QSGW runs (`LiTi2O4/input/qmlo`: `t_tetrakbt = -992.4`, `t_sigmaw = 1000`) with the scripts and figures, the finite-T QSGW runs (`t_tetrakbt` > 0), the control experiment on bcc Fe, the reports on the GPU precision, and the research log. **Inputs and results only** — the runs are too expensive for a regression test (the test targets at finite temperature are `TestInstall/fe_kbt` and `kBT/scanT/Si`) | [kBT/README.md](https://github.com/tkotani/ecalj/blob/main/Samples/kBT/README.md) | [./kBT](./kBT) |

### `Samples/MLOsamples/` quick map (25 testecalj cases)

The dedicated [MLOsamples/README.md](https://github.com/tkotani/ecalj/blob/main/Samples/MLOsamples/README.md)
table has roles + reference data per case.  Highlights:

- `GaAsSoc`, `FeSoc`, `FeMgOSoc` — SOC-as-perturbation tests (`job_mlo_soc`; compared with the DFT bands with SOC that it draws, [mlo](./mlo) §4)
- `C`, `C.sp`, `Cu`, `SrTiO3`, `GaAs`, `Si666gwsc` (QSGW) — non-magnetic semiconductors / metals
- `Ag`, `Al`, `NaCl`, `SiC`, `CdTe`, `ZnO`, `TiO2` — Materials Project structures with the default settings
- `Fe`, `FeCo`, `FeMgO` — magnetic metals (no SOC)
- `Al2O3_Cr`, `GdCo5`, `GdION`, `SmP`, `NiO666lda`, `RuO2` — LDA+U / 4f / AFM
- `Fe`'s test.py drives `job_mloW` and verifies on-site V, W−V vs
  inline reference values — see
  [MLOsamples/README.md § Regression test (testecalj Fe)](https://github.com/tkotani/ecalj/blob/main/Samples/MLOsamples/README.md#regression-test-testecalj-fe).

### `Samples/TestInstall/` what each tests

Driven by `testecalj --all -np 8` from inside the directory.  Coverage
table (one-line role per dir):

| dir | tool exercised |
|---|---|
| `c`, `te`, `zrt`, `co`, `cr3si6`, `felz`, `gasls`, `eras`, `crn`, `cu`, `na`, `fe`, `gdn` | LDA / lmf install validation (`fe`: spin-polarized DOS, `gdn`: LDA+U) |
| `copt` | LDA + force calculation |
| `gas_eps_lmfh`, `gas_epsPP_lmfh` | dielectric ε(q,ω) (with / without LFC) on GaAs |
| `fe_epsPP_lmfh_chipm` | spin susceptibility χ⁻⁻ on bcc Fe |
| `si_gw_lmfh`, `gas_pw_gw_lmfh`, `si_gwsc`, `gas_gwsc`, `nio_gwsc`, `fe_gwsc` | one-shot GW vs full QSGW (`gwsc`) |
| `fe_kbt` | QSGW with the finite-temperature tetrahedron method (`t_tetrakbt` > 0) on bcc Fe |
| `ni_crpa`, `srvo3_crpa` | constrained RPA |

### EPS samples — what differs across the three

Detailed worked descriptions live in [./optical](./optical).  Quick
distinguishing features:

| dir | sname | binary | BZ mesh | QforEPS | what's specific |
|---|---|---|---|---|---|
| `EPS_Cu` | `cu` | `job_eps` | 12³ (lmf) / 10³ (gw) | 3 small q (5e-4 / 1e-3 / 2e-3 a.u.) | FCC metal, has band plot scripts |
| `EPS_GaAs` | `gaas` | `job_eps` | 4³ / 8³ | same Cu-style 3-q | zinc-blende semiconductor |
| `EPS_Ag` | `ag` | `job_eps` | 8³ / 8³ | Cu-style | FCC metal, no LFC; reference data 2026-05-08 generated |

All three carry committed reference EPS files (`EPS000{1,2,3}.nlfc.dat.{inter,intra}bandonly`)
and the test.py compares them with skipcond on the intra-band 1/(ω+iη)
divergence.  `testecalj EPS_Cu / EPS_GaAs / EPS_Ag -np 4` each yields
6/6 PASSED.

## Legacy/ — gone

The older examples of `Samples/Legacy/` were rebuilt as the directories above (2026-09-29 and 30) and the directory
was removed from the repository on 2026-09-30 (ecalj history before commit `932c59a6d`). The last one,
`TestHomoDimerAtom` (dimers and atoms in a box, 2012), became `AtomDimer/N2`.

## Where each Samples-related description in ecaljdoc lives

To minimise duplication across pages, each topic is owned by one
page; this Samples page only points:

| topic | owner page |
|---|---|
| TOML migration / `Legacy2toml.py` workflow | [./toml_migration](./toml_migration) |
| `ctrlg.<sname>.toml` schema and the legacy ↔ TOML key map | [./lmf](./lmf) |
| keys of `[gw]` (the former `GWinput`) | [./gwinput](./gwinput) |
| finite temperature, the keys `t_tetrakbt` / `t_sigmaw` | [./kBT](./kBT) |
| MLO-QSGW (`gwsc --mlo`) | [./mlo_gwsc](./mlo_gwsc) |
| QSGW with `gwsc` end-to-end | [./gwsc](./gwsc) |
| Dielectric ε(q,ω) (`job_eps`) | [./optical](./optical) |
| Fat-band weights / `PROCAR` (`MgO_PROCAR`, `Ni2MnGa_L21_PROCAR`) | [./UsageDetailed § PROCAR mode](./UsageDetailed#procar-mode) |
| MLO bands + on-site W (`job_mloW`) | [MLOsamples/README.md](https://github.com/tkotani/ecalj/blob/main/Samples/MLOsamples/README.md) |
| `ecalj_auto` batch runner | `MD/auto.md` |
| End-to-end tutorial (POSCAR → QSGW) | [./README_tutorial](./README_tutorial) |

If you want to add a new Samples directory, the conventions live in
`MD/developer.md` (`testecalj`, `comp.py`, naming).
