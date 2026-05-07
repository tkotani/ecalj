# MTO-based Localized Orbital (MLO) — sample tests and how to run

`mlo` builds a localized real-space tight-binding-style Hamiltonian on top of
ecalj's PMT (MTO + APW) basis. The samples in this directory exercise four
representative cases. The SOC-as-perturbation variant (`job_mlo_soc`) is
documented separately in [README_SOC.md](README_SOC.md).

## Quick start

```bash
cd ~/ecalj/Samples/MLOsamples
# SOC test set (job_mlo + job_mlo_soc)
testecalj -np 8 GaAsSoc       # semiconductor, with SOC
testecalj -np 8 FeSoc         # ferromagnetic Fe (nspin=2), with SOC
testecalj -np 8 FeMgOSoc      # FeMgO multilayer (nspin=2), with SOC

# Plain MLO non-SOC
testecalj -np 8 Si666gwsc     # Si QSGW, non-magnetic, MLO + DFT band check
testecalj -np 8 Al2O3_Cr      # Cr-doped Al2O3 QSGW80, mlo_emax=7 example
testecalj -np 8 C             # diamond C
testecalj -np 8 C.sp          # graphite-like C (sp-only)
testecalj -np 8 Cu            # fcc Cu
testecalj -np 8 SrTiO3        # perovskite, non-magnetic
testecalj -np 8 Fe            # bcc Fe (nspin=2, no SOC)
testecalj -np 8 FeCo          # FeCo alloy (nspin=2)
testecalj -np 8 FeMgO         # FeMgO multilayer (nspin=2, no SOC)
testecalj -np 8 GdCo5         # GdCo5 (Gd 4f, nspin=2)
testecalj -np 8 GdION         # Gd ion QSGW
testecalj -np 8 NiO666lda     # NiO LDA AFM
testecalj -np 8 RuO2          # RuO2 QSGW
testecalj -np 8 SmP           # SmP (Sm 4f + so=2)
```

Each test runs the full pipeline (lmf → mlo, plus the 4-step
`job_mlo_soc` flow on `*Soc` samples), then numerically diffs the produced
`band_MLO_spin{1,2}*.dat` against committed references with tolerance
`7.4e-5 Ry ≈ 0.001 eV`. PASS line:

    max deviation = 0.0  tolerance = 7.4e-05

## Sample summary

| dir | system | nspin | SOC test | reference |
|---|---|---|---|---|
| `Si666gwsc`  | Si QSGW (6×6×6)        | 1 | – | `bnd00{1..6}.spin1`, `band_MLO_spin1.dat` |
| `GaAsSoc`    | GaAs QSGW              | 1 | ✓ | `band_MLO_spin1.dat`, `band_MLO_spin1.soc.dat` |
| `FeSoc`      | bcc Fe DFT             | 2 | ✓ | `band_MLO_spin{1,2}.dat`, `band_MLO_spin1.soc.dat` |
| `FeMgOSoc`   | FeMgO multilayer DFT   | 2 | ✓ | `band_MLO_spin{1,2}.dat`, `band_MLO_spin1.soc.dat` |
| `Al2O3_Cr`   | Cr-doped Al2O3 QSGW80  | 2 | – | `band_MLO_spin{1,2}.dat` (mlo_emax=7) |
| `C`          | diamond C              | 1 | – | `band_MLO_spin1.dat` |
| `C.sp`       | C graphite-like (sp)   | 1 | – | `band_MLO_spin1.dat` |
| `Cu`         | fcc Cu                 | 1 | – | `band_MLO_spin1.dat` |
| `SrTiO3`     | SrTiO3 perovskite      | 1 | – | `band_MLO_spin1.dat` |
| `Fe`         | bcc Fe DFT (no SOC)    | 2 | – | `band_MLO_spin{1,2}.dat` |
| `FeCo`       | FeCo alloy             | 2 | – | `band_MLO_spin{1,2}.dat` |
| `FeMgO`      | FeMgO multilayer       | 2 | – | `band_MLO_spin{1,2}.dat` |
| `GdCo5`      | GdCo5 (4f magnetic)    | 2 | – | `band_MLO_spin{1,2}.dat` |
| `GdION`      | Gd ion QSGW            | 2 | – | `band_MLO_spin{1,2}.dat` |
| `NiO666lda`  | NiO LDA AFM            | 2 | – | `band_MLO_spin{1,2}.dat` |
| `RuO2`       | RuO2 QSGW              | 2 | – | `band_MLO_spin{1,2}.dat` |
| `SmP`        | SmP (4f, so=2 Lz·Sz)   | 2 | – | `band_MLO_spin{1,2}.dat` |

For `Si666gwsc` the MLO test was added on top of an existing `job_band`
DFT-band check, so testecalj runs both. `Al2O3_Cr` is the canonical
example for `mlo_emax > 0` (Cr 3d states sit ~7 eV above E_F).

## Manual run

```bash
cd <target>_work          # produced by testecalj
gnuplot -p bandplot_MLO.isp1.glt        # red = MLO bands, line = DFT bands
gnuplot -p bandplot_MLO.isp2.glt        # spin-down (FeSoc / FeMgOSoc only)
```

For SOC samples, the same `bandplot_MLO.isp1.glt` shows the 2N-spinor MLO
overlay after `job_mlo_soc` finishes (Fermi level reset to `efermi_soc`).
See `README_SOC.md` for what to expect.

To run only the MLO step manually:

```bash
job_mlo gaas -np 8                       # non-SOC
job_mlo_soc gaas -np 8                   # SOC as perturbation
```

## What gets shown in the plot

Both DFT and MLO bands plotted on the same panel (`Energy − E_F`, eV).

- `Si666gwsc`: dispersive sp bands. MLO traces the lower 4 valence + a few
  conduction bands faithfully along Γ–X–U–Γ–L–W–X.
- `GaAsSoc`: VBM around Γ shows the As-p triplet. With `job_mlo_soc`,
  the Δ_SO ≈ 0.34 eV split-off band is reproduced.
- `FeSoc`: Fe 3d manifold ±2 eV around E_F, both spins.
- `FeMgOSoc`: dense band structure of the Fe-MgO interface region.

## Settings to know (in `GWinput`)

```
mlo_method 0       # 0 = standard, 1 = emax-cut, 2 = MTO-energy cut
mlo_emax   0       # eV relative to E_F (semiconductor: 0 ≈ at E_F).
                   # User-set 0 is honoured literally (sentinel default for
                   # "key absent" → emax*Ry, applied at runtime).
                   # Al2O3_Cr-style insulators: try 7 eV.
<Worb>             # which lm orbitals per atom enter the MLO basis
1 Ga 1 2 3 4 5 6 7 8 9
2 As 1 2 3 4 5 6 7 8 9
</Worb>
```

The `<Worb>` indices map to `(l, m)`: `1 → s; 2..4 → p; 5..9 → d; 10..16 → f`
(see `bandplot.isp1.glt` header for the full table).

## cRPA W on MLO (separate flow)

```bash
job_mloW fe -np 8
```

Produces:
- `Coulomb_v.*` — bare V on MLO basis
- `Screening_W-v.*` — cRPA-screened W − V on MLO basis

What actually matters for downstream analysis:
- **On-site (one-atom) blocks** of V and W − V (the U/J parameters in
  Hubbard-style models live here).
- **W at ω = 0** (static screened interaction).

The full ω-dependence is stored on disk (frequency mesh from
`HistBin_dw` / `HistBin_ratio`), but for typical use cases only the
ω = 0 slice and the on-site sub-block need post-processing. cRPA
mode may need minor follow-up; flag if problems show up.

## Recent fixes (2026-05)

The MLO + TOML migration introduced two regressions, fixed in commits
`10e642915` / `c11663856` / `86e58d0c4` / `802878e82`:

- **`mlo_emax 0`**: TOML reader treated 0 as "use default", overriding
  user intent. Now sentinel-based — explicit 0 is honoured.
- **`<Worb>` duplicate blocks**: `gwinput2toml.py` last-write-wins silently
  dropped real data when GWinput had a real block plus a commented
  example. Now keeps the first.
- **`job_mlo_soc` -v overrides**: `-vnspin=2 -vso=0` form is silently
  ignored under TOML-only mode; replaced with `-v[ham.nspin]=2`,
  `-v[ham.so]=0/1`.
- **gfortran 13.3 miscompile workaround** in `m_HamPMT.f90` (legacy
  reader kept as dead code in the else branch — removing it triggers
  a partial-array codegen bug).

See `README_SOC.md` for SOC-specific implementation notes.
