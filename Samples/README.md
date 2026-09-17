# ecalj/Samples — methodological tests and example calculations

As of 2026-05 the Fortran binaries (`lmf`, `lmfa`, `lmchk`, `gwsc`,
`hsfp0`, ...) read structured TOML only:

- `ctrlg.<sname>.toml` — ctrl + GW driver sections; `[product_basis]` at the end carries the cut-offs and the per-atom tables (nlx/valence/core)

The five directories below are already migrated and pass `testecalj`.
The rest are kept under [Legacy/](Legacy/) and still use legacy
`ctrl.<sname>` + `GWinput`; run `Legacy2toml.py <sname>` inside each
working dir to convert before invoking `lmf`/`gwsc`/etc.

## TOML-ified samples (testecalj-driven)

| dir | role | what's inside |
|---|---|---|
| [GetStarted/](GetStarted/) | minimal seeds for the [ecaljdoc tutorial](https://ecalj.github.io/ecaljdoc/manual/README_tutorial#getstarted) | `GaAs/` ships `ctrls.gaas` + `ctrlg.gaas.toml`. See [GetStarted/README.md](GetStarted/README.md). |
| [MLOsamples/](MLOsamples/) | MuffinTin Localized Orbitals (Wannier replacement) | 17 samples — semiconductors, magnetic metals, multilayers, 4f systems. See [MLOsamples/README.md](MLOsamples/README.md). |
| [TestInstall/](TestInstall/) | install validation suite | 25 samples — ground-state (incl. `fe` spin-pol DOS, `gdn` LDA+U), GW (gwsc), eps (eps_lmfh, epsPP_lmfh), magnetic susceptibility (chipm), cRPA. Driven by `testecalj --all`. Extra heavy GW targets (`cugase2_gwsc222`, `nio_gwsc444`, `pdo_gwsc443`, `gas_gwsc666`) have their own `test.py` and run individually: `testecalj -np 60 <dir>`. |
| [EPS/](EPS/) | dielectric function ε(q,ω) | 3 samples — `EPS_Cu`, `EPS_GaAs`, `EPS_Ag`. epsPP0 with no LFC; small q probe + intra/inter band split. |
| [PROCAR/](PROCAR/) | fat-band weight / orbital projection | 2 samples — `MgO_PROCAR` (O-2p weight on bands), `Ni2MnGa_L21_PROCAR` (per-atom fat band, FM Heusler). |

## Verification status (2026-08-18/19 sweep)

All TOML-ified samples above were executed end-to-end and their
`test.py` checks pass:

- TestInstall `--all` (25 targets incl. new `fe`, `gdn`) — ALL PASSED
- TestInstall heavy GW extras — `cugase2_gwsc222`, `nio_gwsc444`,
  `pdo_gwsc443`, `gas_gwsc666` PASSED with references regenerated
  2026-08 (the 2025-10 references predated the 2026-06 real-axis
  binning fix and the TOML migration).
- EPS 3/3, PROCAR 2/2, MLOsamples 17/17 — ALL PASSED
- BenchmarkTest 2/2 (`inas2gasb2`, `inas4gasb4`) — PASSED against the
  stored QPU.1run references. On a single 32GB GPU run with `-np2 1`
  (two GW ranks on one GPU run out of memory for these sizes).
- GetStarted/GaAs — tutorial seed verified (lmfa + lmf converge).
- Legacy — every dir with a usable `ctrl.<sname>` (32 targets) is now
  TOML-ified and passes an lmfa smoke test; `AFsymmetry/NiO`+`NiSe`
  are full testecalj targets again (fresh rst + references).

Known issues:

- `Legacy/Magnon/*`: RESOLVED 2026-08-19 — the hmagnon crash was a
  bootstrap-order bug (cmdopt registry read before MPI__Initialize
  populated it, so --geteta/--sp1/--sp2 were silently ignored); fixed
  in main_hmagnon/main_mlo_magnon/main_hhomogas. All four magnon
  tests now PASS with references regenerated 2026-08 (kt1, np=60).
  Note TrKpm/TrRpm values near magnon resonances are pole-position
  sensitive; cross-compiler comparisons may need loose tolerances
  there.
- `idu>=10` (LDA+U with the 10-offset mode) combined with `symgrpaf`
  crashes zhev_tk4 (nev/=nevout) — reproduced with NiSe; its sample
  now runs without the (U=0, no-op) idu block.

- `TestInstall/yh3fcc_gwsc666`: broken with the current code — the
  radial solver runs away (e~20 Ry) on the Y `pz=[4.9,4.9]` extended
  local orbitals and lmf aborts on a huge allocation. Needs a
  physics-level revisit of the basis setup; kept with its (stale)
  references for now.
- `MLOsamples/Fe_job_mloW`: work-in-progress dir, no test.py.
- PROCAR band weights require np-safe concatenation (fixed 2026-08 in
  `pylib merge_files` / MgO test.py); if you concatenate
  `PROCAR.UP.*` by hand, sort the rank suffix numerically.

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
