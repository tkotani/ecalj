# What's new (summary)

> English version of the page [manual/whatsnew](../../manual/whatsnew.md). The pages at the top of ecaljdoc are the reference; this one was made from manual/whatsnew.md of commit 0b4c05734 (source blob 6261a5fd4927).

A summary of new features of ecalj and of fixes that change results. Newest first. The detailed record (all changes) is ecalj's
[Changes.md](https://github.com/tkotani/ecalj/blob/main/Changes.md). Each item is explained on the linked page.

## Recent major topics (2026-09 to 10)

### 1. Temperature and smearing are only `t_tetrakbt` and `t_sigmaw`; results of the LiTi₂O₄ iterations (tf32, fp32, fp64)

- **Temperature and smearing are now specified only by the two `[gw]` keys `t_tetrakbt` and `t_sigmaw`** (unit K) → [kBT](./kBT.md) Table 4
  - `t_tetrakbt` (χ0 side, required): 0 is the T = 0 tetrahedron method (no smearing); a positive value T is the finite-temperature tetrahedron method at temperature T;
    a negative value −T smears the Im χ0 made by the T = 0 tetrahedron method with a Gaussian of the same width as the Fermi–Dirac function at temperature T → [kBT](./kBT.md) §2
  - `t_sigmaw` (Σ side, <del class="clfix">default 1000</del><span class="clfix">required since 2026-11-14</span>): width of the Fermi–Dirac kernel that smears the intermediate levels of Σx and Σc → [kBT](./kBT.md) §3
  - The former keys `SmearX0`, `esmr`, `tetrakbt = true`, `t_sigmakbt` are not used (if `SmearX0` remains, the program stops; `ctrlg_update.py` rewrites it.
    History of the keys in Table 2 of [kBT](./kBT.md) §2.5)
  - `wcsmear` (default true) is not a temperature but a switch to use the same kernel also for the pole positions in the real-axis pole term of Σc. Keep the default.
    It smooths the cause of rough bands in metals depending on iteration or k mesh (the real-axis pole term of Σc picks up the plasmon pole of W at small q) → [kBT](./kBT.md) §3.5
- **MLO-QSGW of LiTi₂O₄ (metallic spinel, 14 atoms) for 40 iterations from LDA** (6³ and 9³, `gwsc --mlo`). Changing the precision of the GPU matrix products (`gwsc --prec=tf32|fp32|fp64`),
  the difference of the MLO bands is a few meV or less (at 10 iterations on 6³, tf32 vs fp32 0.8 meV for t2g; in the 1st iteration, fp64 vs fp32 1 meV), one order of magnitude smaller than
  the differences due to iteration and mesh (tens of meV). You may run with the fast tf32 (1 iteration: 6³ tf32 3.7 min, fp64 17 min; kt1, RTX 5090 × 2). Stopping at 10 iterations is midway through convergence; 30 to 40 iterations are needed
  → [Summary of LiTi₂O₄](https://github.com/tkotani/ecalj/blob/main/Samples/kBT/LiTi2O4/README.md) (Tables 3 to 5, Figs. 1 and 2), [MLO-gwsc](./mlo_gwsc.md), [GPU version](./ecaljgpu.md)

### 2. MLO and magnons

- **MLO now uses Löwdin-orthogonalized functions (projected Wannier functions) as the standard** (2026-10-02). They keep orbital names and symmetry, and the model is only $H(\mathbf R)$ in the orthogonal basis.
  $U$, $J$, cRPA (`job_mloW`) and magnons (`job_mlo_magnon`) also use this basis; the path through Wannier functions (maximally localized) was removed → [MLO](./mlo.md) §6
- **Magnons** (bcc Fe, FeCo, Ni): no longer depend on the MLO window, and are closer to the previous calculations with Wannier functions. For bcc Fe they nearly coincide for q ≥ 0.4, and for q ≤ 0.3
  MLO is about 20 % higher. The scaling factors of $W$ for the Goldstone condition are Fe 1.24, FeCo 1.28, Ni 1.75
  → [Fe magnon sample](https://github.com/tkotani/ecalj/blob/main/Samples/Magnon/Fe_mlo_magnon/README.md) (Table 2, Fig. 1)
- **cRPA $U$**: for Ni d, MLO 2.90 eV, Wannier 3.78 eV (difference in the choice of subspace; the fully screened RPA $U$ is 1.43 and 1.58 eV, close to each other) → [MLO](./mlo.md) §6

## 2026-10-03

- **[GetStarted with AI](./getstartedAI.md) (New!)**: give the URL of this page to an AI (Claude Code etc.), and it guides you interactively through a summary of ecalj, what it can do, the reading order,
  installation, the first calculation and tests (with a table of times per machine). The page is written to be read by an AI

## 2026-10-02

- **MLO became Löwdin-orthogonalized functions**. These are projected Wannier functions of the MLO subspace, keeping orbital names (t₂g, e_g, etc.) and symmetry.
  The model is held only as $H(\mathbf R)$ in the orthogonal basis; spin-orbit, $U$, $J$, cRPA (`job_mloW`), magnons (`job_mlo_magnon`), and the $\Sigma$ of MLO-QSGW also use this basis.
  Bands on the mesh are unchanged, and interpolation off the mesh improves near $E_F$ (median of the maximum error over 63 materials 0.019 → 0.011 eV).
  Magnons no longer depend on the MLO window (`mlo_delta`, `mlo_w`) and are closer to the previous calculations with Wannier functions.
  → [MLO](./mlo.md) §6 "The MLO standard", Table M8
  - **The old `HamRsMLO` is not read**. Rerun `job_mlo`. The previous non-orthogonalized model is `job_mlo <sname> --mlo_raw` (only for comparison)
  - Criterion 3 of the model check (`mlo_bandcheck.py`) became "band spikes" (a model band jumps at one point; breakdown of linear independence) → [MLO](./mlo.md) §9 Eq. (12)
- **Symmetry from spglib** (when `symgrp = "find"`; spglib 2.6.0 bundled, Python not needed). The operations found are written to
  `symmetry.<sname>.json` in the working directory and read from then on. Pure translations of supercells and magnetic symmetry of antiferromagnets (`af` of `[[site]]`) are included.
  Inputs with generators written in `symgrp` work as before. The old search (`ECALJ_SYMFIND=ecalj`) was removed (all input styles can still be used as they are) → SYMGRP in [lmf](./lmf.md)
  - The order of operations differs from before, so outputs that write the order as is (BoltzTraP `.struct` etc.) change order
- **Antiferromagnetic symmetry: write `symgrpaf = "find"`** and spglib determines it from the `af` marks (no need to write generators). When using `symgrpaf`, use `pwmode = 11`
  → Antiferro symmetry in [UsageDetailed](./UsageDetailed.md)
- **The path through Wannier functions (maximally localized), the anomalous Hall conductivity (AHC), and `lmfham2` were removed**. The only model is MLO. cRPA is `job_mloW <sname> --crpa`,
  spreads `mlo_spread.py`, magnons `job_mlo_magnon`. The version before the removal is git tag `last-wannier` → [MLO](./mlo.md) §6
- MLO normalized in real space ($U$, $J$ become values for orbitals of unit norm; $W$ of Fe d moved by 1 to 6 %) → [MLO](./mlo.md) §6 Eq. (7a)
- Bug fixed: `lmchk --getwsr` (estimate of MT radii) stopped for spin-polarized input

## 2026-10-01

- **Defaults of the MLO model were organized** (criteria 1, 2, 3<span class="clfix">: criterion 1 is the EH plus the automatic local orbitals, criterion 2 adds EH2 on the s,p of the cations, criterion 3 adds empty spheres in the voids; [MLO](./mlo.md) §9 Table M3</span>). Semicore local orbitals are selected automatically from the band positions (no input key). gwinit writes f for 4f atoms and
  `mlo_lm2` with `!` (removing it gives criterion 2, which fits MgTe, ZnTe, CdTe, AlN etc. with large voids) → [MLO](./mlo.md) §1, §9
- **Model check** `mlo_bandcheck.py` (run at the end by `job_mlo` and `job_mlo_soc`, prints CHECK PASS/FAIL) → [MLO](./mlo.md) §9
- `job_mlo_soc` also plots the DFT bands with spin-orbit for comparison → [MLO](./mlo.md) §4
- Put the Materials Project API key in `MaterialProject.key` at the top of ecalj (ignored by git) or in the environment variable `MP_API_KEY`

## 2026-09-30

- **Fixes that change results**:
  - With `pwmode = 11` (default), the density matrix of LDA+U, the orbital moment with spin-orbit, and `--cls` read eigenvectors with the dimension of a different k (since 2026-03-30).
    LDA+U SCF with `pwmode = 11` in this period is wrong
  - `idu = 10 + mode` of `[[spec]]` (`idu = 12` for 4f etc.) was LDA+U without subtracting the double-counting term when there is no `sigm` (since 2023-09)
  - QSGW with antiferromagnetic symmetry (`symgrpaf`) stopped → now works by computing both spins
- `--ctrlg:<path>=<value>` now also adds and applies keys not in the file (previously the run went on with the default value) → [Input file](./toml_migration.md)
- `gwscconv`: for metals and semimetals, convergence is judged by the change of eigenvalues → [gwsc](./gwsc.md)
- Samples were reorganized (SOC magnetic anisotropy, LDA+U, effective mass, structural relaxation, DOS, IIR, molecule in a box, electron-temperature scan, etc.) → [Samples](./samples.md)

## Up to 2026-09 (major topics)

| Period | Topic | Explanation |
| --- | --- | --- |
| 2026-09 | MLO-QSGW (QSGW interpolating the self-energy in the MLO representation, `gwsc --mlo`) | [MLO-gwsc](./mlo_gwsc.md) |
| 2026-09 | Finite temperature (`t_tetrakbt` for χ0, level width `t_sigmaw` for Σ, `wcsmear`) | [kBT](./kBT.md) |
| 2026-09 | GPU acceleration of GW (precision switch `--prec=tf32\|fp32\|fp64`, method chosen automatically from a table) | [GPU version](./ecaljgpu.md) |
| 2026-09 | Automatic MLO window (`mlo_method = 4`), treatment of local orbitals | [MLO](./mlo.md) |
| 2026-06 | Finite-temperature tetrahedron method, GPU for `gw_lmfh` | [kBT](./kBT.md), [gwsc](./gwsc.md) |
| 2026-05 | Input unified into `ctrlg.<sname>.toml` (TOML), override with `--ctrlg:` | [Input file](./toml_migration.md) |
