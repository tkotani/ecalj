# ecalj/Samples — methodological tests and example calculations

As of 2026-05 the Fortran binaries (`lmf`, `lmfa`, `lmchk`, `gwsc`,
`hsfp0`, ...) read structured TOML only:

- `ctrlG.<sname>.toml` — merged ctrl + GW driver sections + PB cut-offs
- `PB.toml` — per-atom product basis tables (sname-free)

The four directories below are already migrated and pass `testecalj`.
The rest are kept under [Legacy/](Legacy/) and still use legacy
`ctrl.<sname>` + `GWinput`; run `Legacy2toml.py <sname>` inside each
working dir to convert before invoking `lmf`/`gwsc`/etc.

## TOML-ified samples (testecalj-driven)

| dir | role | what's inside |
|---|---|---|
| [MLOsamples/](MLOsamples/) | MuffinTin Localized Orbitals (Wannier replacement) | 17 samples — semiconductors, magnetic metals, multilayers, 4f systems. See [MLOsamples/README.md](MLOsamples/README.md). |
| [TestInstall/](TestInstall/) | install validation suite | 23 samples — ground-state, GW (gwsc), eps (eps_lmfh, epsPP_lmfh), magnetic susceptibility (chipm), cRPA. Driven by `testecalj --all`. |
| [EPS/](EPS/) | dielectric function ε(q,ω) | 3 samples — `EPS_Cu`, `EPS_GaAs`, `EPS_Ag`. epsPP0 with no LFC; small q probe + intra/inter band split. |
| [PROCAR/](PROCAR/) | fat-band weight / orbital projection | 2 samples — `MgO_PROCAR` (O-2p weight on bands), `Ni2MnGa_L21_PROCAR` (per-atom fat band, FM Heusler). |

## Legacy (under Legacy/, awaiting TOML migration)

These still depend on the legacy `ctrl.<sname>` / `GWinput` text format.
They will not run under current `lmf` until converted via
`Legacy2toml.py <sname>` (or the equivalent in the directory's
instructions).

### with test.py (testecalj-runnable once migrated)

| dir | role |
|---|---|
| [Legacy/Magnon/](Legacy/Magnon/) | magnon spectra — Fe, FeCo, Ni, Fe_bcc_in_sc (epsPP_magnon_chipm). Heavy: 4-5 GB work-dirs. |
| [Legacy/AFsymmetry/](Legacy/AFsymmetry/) | antiferromagnetic symmetry tests — NiO (rocksalt AFM), NiSe (NiAs-type AFM). |

### example / 教材 (no automated test)

| dir | role |
|---|---|
| [Legacy/AHC/](Legacy/AHC/) | Anomalous Hall Conductivity — bcc Fe, with/without SOC. |
| [Legacy/AHCSOCtest/](Legacy/AHCSOCtest/) | AHC + SOC variants. |
| [Legacy/BOLZTRAP/](Legacy/BOLZTRAP/) | generate Boltztrap input files (`--boltztrap` option). |
| [Legacy/CMDsample/](Legacy/CMDsample/) | command-line invocation samples — Si, InAs, ZnS GGA/LDA. |
| [Legacy/FermiSurface/](Legacy/FermiSurface/) | Fermi surface plotting (`job_fermisurface`) — Cu reference. |
| [Legacy/GdNldau/](Legacy/GdNldau/) | GdN with LDA+U (4f localized). |
| [Legacy/IIR/](Legacy/IIR/) | impact ionization rate calculation. |
| [Legacy/InAsGaSb/](Legacy/InAsGaSb/) | InAs/GaSb superlattice. |
| [Legacy/LaGaO3_relax/](Legacy/LaGaO3_relax/) | atomic position relaxation (LaGaO3, 20 atoms/cell). |
| [Legacy/MATERIALS/](Legacy/MATERIALS/) | per-material starter inputs (large catalog). |
| [Legacy/ReNcub/](Legacy/ReNcub/) | rare-earth nitride cubic structure — systematic series. |
| [Legacy/SLAB/](Legacy/SLAB/) | slab geometry calculations. |
| [Legacy/SOC/](Legacy/SOC/) | spin-orbit coupling tests. |
| [Legacy/SOCAXIS/](Legacy/SOCAXIS/) | SOC axis-dependent MAE — Hsoc(100) vs Hsoc(110). |
| [Legacy/Samples_ISSP/](Legacy/Samples_ISSP/) | ISSP-supplied tutorial samples. |
| [Legacy/Si_doping_sample/](Legacy/Si_doping_sample/) | Si doping example. |
| [Legacy/TETRAHEDRON_HomoGas/](Legacy/TETRAHEDRON_HomoGas/) | dielectric function of homogeneous electron gas (Lindhard). |
| [Legacy/TETRAHEDRON_HomoGas_test/](Legacy/TETRAHEDRON_HomoGas_test/) | parameter-swept variant of above. |
| [Legacy/TestHomoDimerAtom/](Legacy/TestHomoDimerAtom/) | homonuclear dimer + atom in supercell. |
| [Legacy/UUmatSOC/](Legacy/UUmatSOC/) | U-U matrix elements with SOC. |
| [Legacy/bga2o3deformation/](Legacy/bga2o3deformation/) | β-Ga2O3 under uniaxial deformation. |
| [Legacy/mass_fit_test/](Legacy/mass_fit_test/) | effective mass fitting from band curvature. |
| [Legacy/superlattice/](Legacy/superlattice/) | superlattice geometry generator. |
| [Legacy/BK/](Legacy/BK/) | local backup directory (no specific test). |

## Running the tests

```bash
# all install tests (TOML-ified)
cd Samples/TestInstall && testecalj --all -np 8

# specific test
testecalj <dirname> -np <N>

# convert legacy inputs in cwd to TOML
Legacy2toml.py <sname>
```

See [MLOsamples/README.md](MLOsamples/README.md) for the MLO-specific
documentation and screened-W (job_mloW) flow.
