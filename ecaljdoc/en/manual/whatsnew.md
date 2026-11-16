# What's new (summary)

> English version of the page [manual/whatsnew](../../manual/whatsnew.md). The pages at the top of ecaljdoc are the reference; this one was made from manual/whatsnew.md of commit 54d72993d (source blob 73b38f44f054).

A summary of the new features of ecalj and of the fixes that change results, newest first. The detailed record (all changes) is ecalj's
[Changes.md](https://github.com/tkotani/ecalj/blob/main/Changes.md). Each item is explained on the linked page.

## 2026-11: old ctrlg files stop in GW (how to fix)

::: warning Change of usage (2026-11-14)
`t_sigmaw`, `HistBin_dw`, `HistBin_ratio`, `QforEPSau` and `EMAXforGW` of `[gw]` are now required.
With an old `ctrlg.<sname>.toml` that lacks them, the GW programs (`gwsc` etc.) say so and stop.
One line adds them. The values added are the former code defaults, so **the results do not change**.

```bash
ctrlg_update.py ctrlg.<sname>.toml      # the original is kept as ctrlg.<sname>.toml.bak_update
```

The old `esmr` and `t_sigmakbt` are no longer read; `ctrlg_update.py` converts them to `t_sigmaw`.
The meaning of the keys is in [GW settings](./gwinput.md), the temperature keys in [kBT](./kBT.md).
:::

## Recent major topics (2026-09 to 10)

### Temperature and smoothing are given only by `t_tetrakbt` and `t_sigmaw`; results of the LiTi₂O₄ iterations (tf32, fp32, fp64)

- **Temperature and smoothing are now specified only by the two keys `t_tetrakbt` and `t_sigmaw` in `[gw]`** (unit K) → [kBT](./kBT.md) Table 4
  - `t_tetrakbt` (χ0 side, required): 0 is the T = 0 tetrahedron method (no smoothing); a positive value T is the finite-temperature tetrahedron method at temperature T;
    a negative value −T smooths the Im χ0 made by the T = 0 tetrahedron method with a Gaussian of the same width as the Fermi–Dirac function at temperature T → [kBT](./kBT.md) §3
  - `t_sigmaw` (Σ side, required from 2026-11-14): the width of the Fermi–Dirac kernel that smooths the intermediate levels of Σx and Σc → [kBT](./kBT.md) §4
  - The former keys `SmearX0`, `esmr`, `tetrakbt = true`, `t_sigmakbt` are not used (if `SmearX0` remains, the program stops; `ctrlg_update.py` rewrites it.
    The history of the keys is in Table 2 of [kBT](./kBT.md) §3.5)
  - `wcsmear` (default true) is not a temperature but a switch to use the same kernel also for the pole positions in the real-axis pole term of Σc. Keep the default.
    It smooths the cause of rough bands in metals that depend on the iteration and the k mesh (the real-axis pole term of Σc picks up the plasmon pole of W at small q) → [kBT](./kBT.md) §5
- **MLO-QSGW of LiTi₂O₄ (metallic spinel, 14 atoms), 40 iterations from LDA** (6³ and 9³, `gwsc --mlo`). Changing the precision of the matrix products on the GPU (`gwsc --prec=tf32|fp32|fp64`)
  changes the MLO bands by a few meV or less (at 10 iterations on 6³, tf32 vs fp32 is 0.8 meV for t2g; in the first iteration, fp64 vs fp32 is 1 meV), one order of magnitude smaller than
  the differences due to iterations or meshes (tens of meV). The fast tf32 can be used (one iteration on 6³: tf32 3.7 min, fp64 17 min; kt1, RTX 5090 × 2). Stopping at 10 iterations is midway through convergence; 30 to 40 iterations are needed
  → [Summary of LiTi₂O₄](https://github.com/tkotani/ecalj/blob/main/Samples/kBT/LiTi2O4/README.md) (Tables 3 to 5, Figs. 1 and 2), [MLO-gwsc](./mlo_gwsc.md), [GPU version](./ecaljgpu.md)

### MLO and magnons

- **MLO now uses Löwdin-orthogonalized functions (projected Wannier functions) as the standard** (2026-10-02). They keep the orbital names and symmetry, and the model is only $H(\mathbf R)$ in the orthogonal basis.
  $U$, $J$, cRPA (`job_mloW`) and magnons (`job_mlo_magnon`) also use this basis; the Wannier-function (maximally localized) path was removed → [MLO](./mlo.md) §6
- **Magnons** (bcc Fe, FeCo, Ni): they no longer depend on the MLO window and have come closer to the earlier calculations with Wannier functions. For bcc Fe they almost coincide for q ≥ 0.4; for q ≤ 0.3
  MLO is about 20 % higher. The scaling factor of $W$ for the Goldstone condition is 1.24 for Fe, 1.28 for FeCo, 1.75 for Ni
  → [Sample of Fe magnons](https://github.com/tkotani/ecalj/blob/main/Samples/Magnon/Fe_mlo_magnon/README.md) (Table 2, Fig. 1)
- **cRPA $U$**: for the Ni d, MLO 2.90 eV, Wannier 3.78 eV (due to the different choice of the subspace; the RPA $U$ with all screening included is 1.43 and 1.58 eV, close to each other) → [MLO](./mlo.md) §6

## 2026-10-03

- **[GetStarted with AI](./getstartedAI.md) (New!)**: give the URL of this page to an AI (Claude Code etc.), and it guides you interactively through a summary of ecalj, what it can do, the reading order,
  installation, the first calculation and the tests (with a table of times per machine). The page is written on the assumption that an AI reads it

## 2026-10-02

- **MLO are now Löwdin-orthogonalized functions**. They are projected Wannier functions of the MLO subspace and keep the orbital names (t₂g, e_g, etc.) and symmetry.
  The model is held only as $H(\mathbf R)$ in the orthogonal basis; spin-orbit coupling, $U$, $J$, cRPA (`job_mloW`), magnons (`job_mlo_magnon`) and the $\Sigma$ of MLO-QSGW also use this basis.
  The bands on the mesh do not change, and the interpolation off the mesh improves near $E_F$ (over 63 materials, the median of the maximum error 0.019 → 0.011 eV).
  Magnons no longer depend on the MLO window (`mlo_delta`, `mlo_w`) and have come closer to the earlier calculations with Wannier functions.
  → [MLO](./mlo.md) §6 "Standard of MLO", Table M8
  - **The old `HamRsMLO` is not read**. Rerun `job_mlo`. The earlier non-orthogonalized model is `job_mlo <sname> --mlo_raw` (only for comparison)
  - Criterion 3 of the model check (`mlo_bandcheck.py`) is now "spikes in a band" (the model band jumps at one point; a breakdown of linear dependence) → [MLO](./mlo.md) §9 Eq. (12)
- **Symmetry from spglib** (when `symgrp = "find"`; spglib 2.6.0 is bundled, Python is not needed). The operations found are written to
  `symmetry.<sname>.json` in the working directory and read from there afterwards. Pure translations of supercells and the magnetic symmetry of antiferromagnets (`af` in `[[site]]`) are included.
  Inputs that give generators in `symgrp` work as before. The old search (`ECALJ_SYMFIND=ecalj`) was removed (all ways of writing the input still work as they are) → SYMGRP in [lmf](./lmf.md)
  - The order of the operations differs from before, so outputs that write the operations in their order (BoltzTraP `.struct`, etc.) change their order
- **For antiferromagnetic symmetry, write `symgrpaf = "find"`** and spglib finds it from the `af` marks (no need to write generators). When using `symgrpaf`, use `pwmode = 11`
  → Antiferro symmetry in [UsageDetailed](./UsageDetailed.md)
- **The Wannier-function (maximally localized) path, the anomalous Hall conductivity (AHC) and `lmfham2` were removed**. The model is MLO only. cRPA is `job_mloW <sname> --crpa`,
  the spread is `mlo_spread.py`, magnons are `job_mlo_magnon`. The version before the removal is the git tag `last-wannier` → [MLO](./mlo.md) §6
- MLO are normalized in real space ($U$ and $J$ become values for orbitals of norm 1; the $W$ of the Fe d moved by 1 to 6 %) → [MLO](./mlo.md) §6 Eq. (7a)
- Fixed error: `lmchk --getwsr` (estimate of MT radii) stopped for spin-polarized inputs

## 2026-10-01

- **The defaults of the MLO model were reorganized** (criteria 1, 2, 3. Criterion 1 is EH and automatic local orbitals; criterion 2 adds EH2 for the s,p of cations; criterion 3 adds empty spheres in the voids. [MLO](./mlo.md) §9 Table M3). Local orbitals for semicore states are chosen automatically from the band positions (there is no input key). gwinit writes the f of 4f atoms and
  `mlo_lm2` with `!` (removing it gives criterion 2; this fits MgTe, ZnTe, CdTe, AlN, etc., which have large voids) → [MLO](./mlo.md) §1, §9
- **Model check** `mlo_bandcheck.py` (run at the end of `job_mlo` and `job_mlo_soc`; prints CHECK PASS/FAIL) → [MLO](./mlo.md) §9
- `job_mlo_soc` also plots the DFT bands with spin-orbit coupling for comparison → [MLO](./mlo.md) §4
- Put the Materials Project API key in `MaterialProject.key` at the top of ecalj (ignored by git) or in the environment variable `MP_API_KEY`

## 2026-09-30

- **Fixes that change results**:
  - With `pwmode = 11` (default), the density matrix of LDA+U, the orbital moments with spin-orbit coupling, and `--cls` read the eigenvectors with the dimension of another k (since 2026-03-30).
    LDA+U SCF with `pwmode = 11` during this period is wrong
  - `idu = 10 + mode` in `[[spec]]` (e.g. `idu = 12` for 4f) became LDA+U without subtracting the double-counting term when there is no `sigm` (since 2023-09)
  - QSGW with antiferromagnetic symmetry (`symgrpaf`) stopped → now it computes both spins and runs
- `--ctrlg:<path>=<value>` now also adds keys that are not in the file and makes them effective (before, the run went on with the default values) → [Input file](./toml_migration.md)
- `gwscconv`: for metals and semimetals, convergence is judged by the change of eigenvalues → [gwsc](./gwsc.md)
- The samples were reorganized (magnetic anisotropy with SOC, LDA+U, effective masses, structural relaxation, DOS, IIR, molecules in a box, scans of electronic temperature, etc.) → [Samples](./samples.md)

## Before 2026-09 (major topics)

| Period | Topic | Description |
| --- | --- | --- |
| 2026-09 | MLO-QSGW (QSGW that interpolates the self-energy in the MLO representation, `gwsc --mlo`) | [MLO-gwsc](./mlo_gwsc.md) |
| 2026-09 | Finite temperature (`t_tetrakbt` for χ0, the level width `t_sigmaw` for Σ, `wcsmear`) | [kBT](./kBT.md) |
| 2026-09 | GPU acceleration of GW (switching precision `--prec=tf32\|fp32\|fp64`, method chosen automatically from a table) | [GPU version](./ecaljgpu.md) |
| 2026-09 | Automatic MLO window (`mlo_method = 4`), treatment of local orbitals | [MLO](./mlo.md) |
| 2026-06 | Finite-temperature tetrahedron method, GPU for `gw_lmfh` | [kBT](./kBT.md), [gwsc](./gwsc.md) |
| 2026-05 | A single input `ctrlg.<sname>.toml` (TOML), overrides by `--ctrlg:` | [Input file](./toml_migration.md) |
