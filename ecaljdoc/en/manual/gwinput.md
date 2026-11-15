# GW settings (the `[gw]` and related keys) {#gwinput-legacy}

> English version of the page [manual/gwinput](../../manual/gwinput.md). The pages at the top of ecaljdoc are the reference; this one was made from manual/gwinput.md of commit 296542c3a (source blob 40bd24ce291b).

The sections `[gw]`, `[mlo]`, `[blocks]` and `[product_basis]` of `ctrlg.<sname>.toml` set the conditions of a GW/QSGW calculation. This page explains their keys.
For the whole file and the conversion from the old `GWinput`, see [Input files](./toml_migration).

## Automatic generation {#mkgwinput}
- The GW settings (`[gw]`, `[mlo]`, `[blocks]`, `[product_basis]` in `ctrlg.<sname>.toml`) are generated automatically, and most parameters can be used at their default values.

When you start from `ctrls.<sname>` (`<sname>` is the material name),
```bash
ctrlgenToml.py <sname>
```
writes `ctrlg.<sname>.toml` including the GW settings (internally it runs `lmfa`, `lmf --jobgw=0` and `gwinit`). To add them to a `ctrlg.<sname>.toml` without the GW sections, use
```bash
ctrlgenToml.py <sname> --addgw   # stops if [gw] exists
```
([Input files](./toml_migration)).

`t_tetrakbt`, `t_sigmaw`, `HistBin_dw`, `HistBin_ratio`, `QforEPSau` and `EMAXforGW` are required; if they are missing in `[gw]`, the programs stop (`t_tetrakbt` since 2026-09-28, the others since 2026-11-14).
For old files, add them with `ctrlg_update.py ctrlg.<sname>.toml`. What it adds are the default values of the earlier code, so the results do not change.

The generated `[gw]` has `n1n2n3 = [4, 4, 4]`, `t_sigmaw = 300`, `t_tetrakbt = 300`. Rewrite them for your calculation.
For insulators and semiconductors, the gap is much larger than $k_BT$, so a temperature left at 300 K does not affect the results ([kBT](./kBT) §7).

## Format
- Written in TOML syntax. From `#` to the end of the line is a comment.
- The form is `key = value`; the order does not matter, but do not write the same key more than once. Logical values are `true` / `false`,
  lists are `[4, 4, 4]`, reals are `1e-5` (`1d-5` cannot be written), and multi-line tables (such as `QforGW`) are enclosed in `"""`.

## Main keys {#parameters-in-gwinput}
The generated `[gw]` describes the meaning of the parameters; below we describe the important variables and the ones not described there.

### `n1n2n3`
Number of k points in the $GW$ self-energy calculation. It can be set differently from `[bz] nkabc` (the k points of `lmf`) (in `gwsc --mlo` it must be the same as `nkabc`; otherwise the program stops). Note that the computational cost is proportional to the square of this number of k points.
- **type** : integer list

We set number of k points for self-energy in `[gw]`. These can be 1/2 or 2/3 of k points specified in `[bz] nkabc`, in order to reduce computational time.
The 6x6x6 k points is feasible setting for ZB structure (2 atoms per cell). That is, 6x6x6x2 =432 \sim `k points \times atom number`.
For example, when we try 8 atoms per cell, 4x4x4 or 3x3x3
is fine because 4x4x4x8 \sim 400.
For metallic systems, larger is fine, but limited by computational
time. See takao kotani's papers, for example, https://doi.org/10.1103/PhysRevB.93.075125
Good news is that we don't need to use so many k points as in `[bz] nkabc`.
In cases, 4x4x4 for ZB is not so bad --- roughly speaking, this is a low limit for
publication. This means 3x3x3x8 (for 8 atom case) is not so bad. (3x3x3x8 > 4x4x4x2).
See our examinations shown in  https://doi.org/10.7566/JPSJ.83.094711 and
http://doi.org/10.7567/JJAP.55.051201

### `KeepEigen`  
Whether to keep the wave functions in memory. The wave functions are written to the files `CPHI` and `GEIG`. If their size is comparable to or larger than the memory available to one MPI process, always set `false`. Set `false` when the number of $k$ points or orbitals is large.
- **type**: boolean
- **default**: true

### `zmel_batch_gb`
Batch size (in GB, per MPI rank) for building the matrix elements (zmel) at a time. A value of about one third to one quarter of the memory available to one MPI process is good. It is used both in the self-energy calculation and in the calculation of $W$ (the polarization function). For large systems a larger value is faster, but if the calculation dies from lack of memory, make it smaller.
- **type** : float
- **default** : 0.4 (values below 0.4 on CPU and below 2.0 on GPU are raised to 0.4 and 2.0, respectively)

If the value is too small for the system, the calculation may stop with `sxcf_fal2_count.sc. Too small memory for nmbatch. Enlarge zmel_batch_gb in ctrlg.<sname>.toml [gw]`. In that case make it larger.

### `zmel_max_size`
The old name of `zmel_batch_gb`. It is still read as `zmel_batch_gb`. `MEMnmbatch` in old inputs is not read (write `zmel_batch_gb`).

### `KeepPpb`
Whether to keep in memory the matrix elements between the product basis in the MT and the wave-function basis in the MT (the variable `ppb`) for all symmetry operations. Usually the default `false` is fine (then only the ones for the symmetry operation in use are kept). With `true`, the memory use grows by the number of symmetry operations.
- **type**: boolean
- **default**: false

### `t_tetrakbt`  *(required, 2026-09-28)*
Selects how $\chi_0$ is smeared, by the sign of a temperature (K). If not written, the program stops.
- `0`: no smearing (T=0 tetrahedron). For insulators and semiconductors the template value 300 may be kept (it has no effect if gap ≫ $k_BT$. [kBT](kBT) §7).
- `T > 0`: finite-temperature tetrahedron. The occupations in $\chi_0$ become Fermi-Dirac at temperature $T$, and the Fermi level is `EFERMI_kbt` (also on the $\Sigma$ side). [kBT](kBT) §3.
- `-T`: keeps the T=0 tetrahedron and smears $\mathrm{Im}\chi_0$ along the frequency with a Gaussian. The width is the same as the Fermi-Dirac at temperature $T$ (standard deviation $\pi k_BT/\sqrt3$, $-1000$ K = 0.00574 Ha = 0.156 eV). Use it for metals where the plasmon pole of $W$ is rough (the practical setting for LiTi₂O₄ is $-992.4$. See the settings for LiTi₂O₄ at the top of [kBT](kBT)). It also affects $W(i\omega)$ and the static $W(0)$.
The gwinit template writes 300.
- **type**: float (K)
- **default**: none (required)

Keys appearing in old inputs are handled as in Table 1. `ctrlg_update.py ctrlg.<sname>.toml` rewrites `SmearX0` (Ha, the standard deviation of the
Gaussian for $\mathrm{Im}\chi_0$) to `t_tetrakbt = -T` of the same width (0.0057 Ha → $-992.4$ K), and adds
`t_tetrakbt = 0` to files without `t_tetrakbt`. The calculation does not change. `Legacy2toml.py` does the same rewriting in the conversion. For the history of the keys, see Table 2 in [kBT](kBT) §3.5.

**Table 1**. Old keys in `[gw]` and how they are handled now

| Old key | Current handling |
|---|---|
| `SmearX0`, `SmearX0q0` | The program stops if written. Write `t_tetrakbt = -T` |
| `tetrakbt` (logical) | Not read. Finite temperature is `t_tetrakbt > 0` |
| `esmr` (Ry, Gaussian) | Not read (since 2026-11-14). `ctrlg_update.py` converts it to `t_sigmaw` with the same standard deviation and writes that (before, it was converted and read when `t_sigmaw` was absent) |
| `t_sigmakbt` (K) | Not read (since 2026-11-14). `ctrlg_update.py` writes the value as it is as `t_sigmaw` |
| `GaussianFilterX0`, `GaussSmear`, `delta`, `dw`, `omg_c`, `WgtQ0P`, `MEMnmbatch` | Not read (no effect even if written) |
| `zmel_max_size` | Read as `zmel_batch_gb` |

### `t_sigmaw`  *(2026-09-20, formerly `esmr` / `t_sigmakbt`)*
Width of the Fermi-Dirac kernel that smears the intermediate levels (the poles of G) in the self-energy $\Sigma = iGW$. Unit K.
The same kernel $-\partial f/\partial\varepsilon$ ($k_BT$ = `t_sigmaw`) is used for all of: the occupation weights in $\Sigma_x = iGv$, and the window weights and the pole positions (`wcsmear`)
of the real-axis pole term of $\Sigma_c = iG(W-v)$.
**This is a numerical width, not a complete finite-temperature self-energy** (Bose factors etc. are not included).
The Fermi level is `EFERMI_kbt` if `t_tetrakbt > 0` ($\chi_0$ also at finite temperature); otherwise it is determined by counting the levels on the k mesh (Table 5 in [kBT](kBT) §1).
The margin of the candidate window for intermediate states and of the real-axis $\omega$ mesh is $E_F \pm 15\,k_BT$. The same kernel also enters the step of the
half residue in the imaginary-axis integral of $\Sigma_c$ (contract: the occupation weights, the window weights of the pole term, and the step on the imaginary axis, these 3 places use the same kernel. [manual/kBT](kBT) §6).

- Required (since 2026-11-14). The template is 300 K. The default of the earlier code was 1000 K, a width close to the old `esmr = 0.01` Ry (Gaussian $\sigma$)
  (the same standard deviation is 872 K; the std of FD is $\pi k_BT/\sqrt3 = 1.81\,k_BT$). The old `esmr = 0.003` Ry (the value in TestInstall, e.g. Si) corresponds to 262 K.
- Insulators: if the width is much smaller than the band gap, it does not affect the results (the difference between 262 K and 1000 K is a few meV for Si).
- Metals: in principle the limit $n_1n_2n_3 \to \infty$, `t_sigmaw` $\to 0$, but a finite width is needed for convergence.
  300 to 1000 K together with a fine k mesh. Choose per material (1000 K for LiTi₂O₄).
- The old `esmr` and `t_sigmakbt` are not read. `ctrlg_update.py` writes `t_sigmaw = esmr[Ry]·13.6057/(1.81·8.617e-5)` (for `t_sigmakbt`, the value as it is).
- **type**: float (K)
- **default**: none (required since 2026-11-14; the template is 300, the default of the earlier code was 1000)

### `wcsmear`  *(2026-09-19, default true since 2026-09-20)*
In the real-axis pole term of $\Sigma_c$, applies the Fermi-Dirac kernel of the intermediate levels (`t_sigmaw`) also to the pole position $\omega_\varepsilon$
and integrates $W_c(\omega)$ with the kernel. `false` gives the former one-point evaluation at the mean energy (3-point Lagrange interpolation,
PRB76,165106 Eq.58). In systems where $W_c$ is smooth, both are nearly the same (a few meV up to 77 meV at most for Si).
Required for metals with a sharp plasmon pole (LiTi₂O₄ at ≲ 1500 K). The cost of the real-axis pole term becomes several times larger (+7 to 26 % in total).
- **type**: boolean
- **default**: true

> Parallelization of the tetrahedron weights (`tetwt5`): since 2026-09-22 **MPI k parallelization** (each rank of the W-build computes only its own k range).
> Starting hgw with `-np2` with more ranks than GPUs (GPUs shared round-robin) makes it faster accordingly. It matters when $\chi_0$ is at finite temperature
> (`t_tetrakbt > 0`, 150 s per q at 9³). The OpenMP version (`omp_tetwt`) tried during 09-19 to 22 was removed.

### `unit_2pioa`
  `unit_2pioa = false`   # false --> a.u.; true --> unit of QpGcut_* are in 2*pi/alat

  Usually `false` is fine.

###  `alpha_OffG`
  `alpha_OffG = 1.0`   # (a.u.) Used in auxially function in the offset-Gamma method.

  Not need to change.

### `QpGcut_psi`

Cutoff ~~energy (unit Ryd)~~ of the IPW basis for the wave functions. Note that the basis functions of the wave functions differ from those of the PMT method used in the `lmf` calculation. The wave functions in the interstitial region are re-expanded in IPWs. IPW means interstitial plane wave.
- **type**: float
- **default** : 4.0
`|q+G| < QpGcut_psi` in units of 1/a.u.
The number of them is written in the output of qg4gw (lqg4gw).
To convert to energy (Ry), square it (**2).

The IPW cutoff for the wave functions is also related to `[ham] pwemax`. For pwemax, 2 or 3 is often tried
(pwemax is an energy cutoff in Ry, $|\mathbf q+\mathbf G|^2 \lt {\rm pwemax}$ ($|\mathbf q+\mathbf G|$ in 1/a.u.).)

### `QpGcut_cou`
Cutoff ~~energy (unit Ryd)~~ of the IPW part of the product basis 
`|q+G| < QpGcut_cou` in units of 1/a.u.
To convert to energy (Ry), square it (**2).
- **type**: float
- **default** : 3.0

The number of them is written in the output of qg4gw (lqg4gw).

### `emax_sigm`
Maximum energy, measured from the Fermi level, of the quasiparticle states whose self-energy is calculated: unit Ryd. In QSGW the self-energy is used to construct the potential, so making it smaller also lowers the accuracy of the potential.
- **type**: float
- **default** : 3.0

## `HistBin_dw`, `HistBin_ratio` — the real-axis mesh of the polarization function {#change-histbin-for-high-resolution-energy-mesh-near-ef-for-metal}
This defines histogram bins to accumulate weight (imaginary part) of the polarization functions. `getfreq` in `m_freq.f90` makes this mesh.

Both are required (since 2026-11-14). The template (the values `ctrlgenToml.py` writes into `[gw]` through `gwinit`) is `HistBin_dw = 1e-4`, `HistBin_ratio = 1.05`.
Until 2026-11-14, if the keys were not written, the code defaults `1e-5`, `1.03` were used. `ctrlg_update.py` adds these values to old files, so the results do not change.

History of the template values:

| Period | `HistBin_dw` | `HistBin_ratio` |
| --- | --- | --- |
| until 2025-08-19 | 2e-3 | 1.08 |
| from 2025-08-22 | 1e-4 | 1.05 |

The finer the mesh (smaller `HistBin_dw`, `HistBin_ratio` closer to 1), the higher the energy resolution near $\omega\sim 0$, but the heavier the calculation.
Usually the GW results do not depend much on this choice (see the explanation below). Make it finer, for example, when you plot $W(\omega)$ in detail near $\omega\sim 0$.

* HistBin_dw and HistBin_ratio specify bins along the real frequency axis which we accumulate imaginary part weight of polarization functions. The bins are (see frhis in m_freq.f90) 
$\omega_i = b\,(e^{a(i-1)}-1)$, $a={\rm HistBin\_ratio}-1$, $b={\rm HistBin\_dw}/a$
where $[\omega_i, \omega_{i+1}]$ $(i = 1, 2, ...,\verb#nwhis#+1)$ is the i-th bin. HistBin_dw is bin width at ω = 0.
The ratio ω(i + 1)/ω(i) for large ω is $e^{a}\approx$ HistBin_ratio. This choice of getting coarser at high energy is because we think W (ω) around ω ∼ 0 gives most important contribution to the GW approximation. If histogram bins are too wide, dielectric function can be less accurate, but results may be not so much affected.
In GW calculation, the Plasmon pole is important. It is determined not by the Drude weight; it usually gives small contribution to the Plasmon pole (for example, Si is described well by a Plasmon pole model, but Si has no Fermi surface). We expect that the GW results are not so sensitive to the choice of HistBin_dw, HistBin_ratio usually. We may use fine mesh when we plot quantities such as W (ω) near ω = 0.
<!-- The ecalj gives ¯W (ω = 0) ∼ 0 for metal; where ¯W (ω) is the effective interaction averaged in
the Γ-cell [?]. And ¯W (ω get closer to v for larger ω. -->

## `niw`
 Number of frequencies along Im axis. Used for integration along imaginary axis. Probably 10 is goog enough. 
See routines `wintzsg.f90`. 

The integration points are $i \omega'(n)= i( 1/x(n) -1)$, where $x(n)$ is the
usual Gaussian-integration points for the interval $[0,1]$. In addition, we give the special analytical treatment for the peaky part at $\omega'= 0$. Out tests shows niw=6 for Si is good enough for 0.01 eV accuracy. The convergence as for niw is quite good. This integration scheme has been developed by Ferdi Aryasetiawan. The number of points should be the one of 6,10,12,16,20,24,32,40,or 48, because we use a subroutine gauss in mate.f90 prepared by Ferdi.

## Q0PChoice
`Q0Pchoice` sets the size of the offset-Gamma points Q0P (the factor `deltaq_scale` on half of the mesh spacing).
`Q0Pchoice = 1` (default) gives `deltaq_scale` = 0.1 (close to the $q\to0$ limit), `Q0Pchoice = 2` gives $1/\sqrt3$ (mean value in a sphere).
`deltaq_scale = x` ($x>0$) in `[gw]` sets the factor directly and overrides `Q0Pchoice`. See lqg4gw for the Q0P in use.

### `deltaw`
real (a.u.) only for one-shot case.
deltaw is the interval for the numerical derivative ∂Σ(ω)/∂ω.  We calculate $\Sigma(\omega),\Sigma(\omega\pm {\rm deltaw})$.
From these values, we can calculate two Z (or second-derivative of Σ(ω)). 

## `<QforGW>` section
This is only for one-shot/spectrum function mode. 
In `[gw]`, `QforGW` gives the q vectors, one per line between `"""` (in the unit of `2pi/alat`. See BZ.html and syml file to check):
```toml
QforGW = """
 0.0 0.0 0.0
 0.5 0.0 0.0
"""
```
If no `QforGW` (or `QforGWIBZ = true`), all irreducible q points are used.

## EMAXforGW
one-shot GW. Set this (eV,above the Fermi energy) to specify up to which bands you calculate. 
Required (since 2026-11-14). The template is 15. 99999 means no upper limit (the default of the earlier code, and the value `ctrlg_update.py` adds).
(we have EMINforGW).

## `--readQforGW` option
We can use any q in the QforGW section
Because of the shifted-mesh method, we can use any q point which is not on the mesh points.
If on the mesh point, we need to prepare less number of eigenfunctions for GW calculations.

## `<QforEPS>`
This is for dielectric function mode. Set in `[gw]` as
```toml
QforEPS = """
 0 0 0.00050
 0 0 0.00100
 0 0 0.00200
"""
```
Because of the shifted-mesh method, we can use any q point which is not on the mesh points.
### QforEPSau 
Whether q is in the unit of `2pi/alat` or in the unit of a.u.(=`bohr^{-1}`). I think `QforEPSau = true` might be convenient.
Required (since 2026-11-14). The template is true. The default of the earlier code was false (the value `ctrlg_update.py` adds).

## Product basis `[product_basis]` {#product-basis}

The product basis (the union of local products in the MTs and plane waves in the interstitial region) is the key to how ecalj represents the screened Coulomb interaction $W$
([The product basis that expands products of wave functions](../theory/gw_product_basis)).
It is specified in the **`[product_basis]` section** of `ctrlg.<sname>.toml` (always at the end of the file) (Table PB1).

Table PB1. Keys of `[product_basis]`

| Key | Content |
|---|---|
| `pb_tolerance`, `pb_lcutmx` | Global limits. These are the ones you may adjust |
| `nlx`, `valence`, `core` | Per-atom tables of radial functions. `gwinit` writes them from `[[spec]]` (through `ctrlgenToml.py` / `Legacy2toml.py`). You rarely edit them by hand |

(Until 2026-09, the three tables were in a separate `PB.<sname>.toml`. The present programs do not read that file; they stop and point to `ctrlg_absorb.py <sname>`.
`ctrlg_absorb.py` adds those tables to `[product_basis]` and keeps the original file as `.bk`.)

The product basis originates from F. Aryasetiawan and O. Gunnarsson ([PRB 49, 16214](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.49.16214)),
and the mixed product basis ([Solid State Commun. 2002](https://doi.org/10.1016/S0038-1098(02)00028-5)) is its successor.

### Limits: `pb_tolerance`, `pb_lcutmx`

```toml
[product_basis]
pb_tolerance = [0.001]   # cut off of linear dependency (default 1e-3)
                         #   larger => fewer product basis. 1d-2 is OK
                         #   when you want a cheaper calculation.
pb_lcutmx    = [4, 4]    # per-site l-cutoff for the product basis.
                         #   4 needed when valence 3d exists (Ni, Ga, ...);
                         #   2 OK for oxygen; 6 needed for 4f / 5f atoms.
                         #   site order = [[site]] order; see lmchk SiteInfo.
```

- `pb_tolerance`: criterion for removing nearly linearly dependent products. The larger it is, the fewer product basis functions. 1e-2 is acceptable when you want a cheaper calculation.
- `pb_lcutmx`: per-site upper limit of $l$ of the product basis in the MT. To reduce computational time, set 2 for s, p atoms such as oxygen (probably the same for other elements of the same period).
  Atoms with 3d in the valence (Ni, Ga, etc.) need 4, and 4f / 5f atoms need 6. Example:
  ```toml
  [product_basis]
  pb_lcutmx = [4, 4, 4, 2, 2, 2]
  ```
  The order of the atoms is the order of the `[[site]]` tables (`lmchk` shows it, and it is also in `SiteInfo.lmchk`).

After changing them, rerun `lmf --jobgw=1` (or the next step of `gwsc`). `hbasfp0` uses these values and writes the size of the resulting product basis to `lbas` / `lbasC`.

### Per-atom tables (`nlx`, `valence`, `core`)

They follow the limits, in the same `[product_basis]`. Every row begins with the site number `iatom` (in the `[[site]]` order).
`gwinit` writes them from `rsmh`/`eh`/`pz` of `[[spec]]` and labels each row. The explanation below is also written in the file as comments, and the three tables are separated by `# ----` rules.

```toml
[product_basis]
pb_tolerance = [1e-3]
pb_lcutmx    = [4]
# ----------------------------------------------------------------
# nlx [iatom, l, nnvv, nnc]: one row per (iatom, l) ...
nlx = [
  [1, 0, 2, 3],
  [1, 1, 2, 2],
  [1, 2, 2, 0],
  [1, 3, 2, 0],
  [1, 4, 2, 0],
]
# ----------------------------------------------------------------
# valence [iatom, l, n, occ, unocc]: one row per valence radial function ...
valence = [
  [1, 0, 1, 1, 1],   # 4s_phi
  [1, 0, 2, 0, 0],   # 4s_phidot
  [1, 1, 1, 1, 1],   # 4p_phi
  [1, 1, 2, 0, 0],   # 4p_phidot
  [1, 2, 1, 1, 1],   # 3d_phi
  [1, 2, 2, 0, 0],   # 3d_phidot
  [1, 3, 1, 0, 1],   # 4f_phi
  [1, 3, 2, 0, 0],   # 4f_phidot
  [1, 4, 1, 0, 0],   # 5g_phi
  [1, 4, 2, 0, 0],   # 5g_phidot
]
# ----------------------------------------------------------------
# core [iatom, l, n, occ, unocc, forX0, forSxc]: one row per core level ...
core = [
  [1, 0, 1, 0, 0, 0, 0],   # 1S
  [1, 0, 2, 0, 0, 0, 0],   # 2S
  [1, 0, 3, 0, 0, 0, 0],   # 3S
  [1, 1, 1, 0, 0, 0, 0],   # 2P
  [1, 1, 2, 0, 0, 0, 0],   # 3P
]
```

Meaning (the same as the old `<PRODUCT_BASIS>` block; only the container has changed):

- `nlx` row `[iatom, l, nnvv, nnc]`: how many radial functions ecalj has for each `l` of that site. `nnvv = 2` is `{phi, phidot}`, `nnvv = 3` adds the local orbital `phiz`. `nnc` is the number of core levels of that `l`. These are derived values, so do not edit them.
- `valence` row `[iatom, l, n, occ, unocc]`: one per radial function. `n = 1` is phi, `2` is phidot, `3` is phiz.
  `occ` / `unocc` are 0/1 switches for including the function in the "occupied" / "unoccupied" set. The product basis is built from the products (occupied set) × (unoccupied set) and thinned by `pb_tolerance`.
  `phidot` is usually excluded (`0 0`) for speed. Including it makes the product basis more accurate.
- `core` row `[iatom, l, n, occ, unocc, forX0, forSxc]`: one per core level (`n = 1` is the deepest). `forX0` / `forSxc` come from the CORE1/CORE2 distinction of PRB 76, 165106 (2007). They are historical; in present calculations they are all 0 in every core row.

`Legacy2toml.py` converts the `<PRODUCT_BASIS>` block of a pre-TOML `GWinput` into this `[product_basis]` section.

## `Smaller IPW related part to reduce computational time`
To reduce compuational time, we may reduce MPB coming from IPS, reduce the size of IPW for psi, reduce emax_sigm and pwemax as follows.
   + QpGcut_cou is fort the Interstitial plane wave (IPW) for MPB. 
   + QpGcut_psi is for expantion of eigenfunctions.
   + emax_sigm is the upper cutoff (relative to the Fermi energy) to calculate self energy.
   + `[ham] pwemax` is the APW basis cutoff for the eigenfunciton.

To reduce computational time, we may use
```toml
[gw]
QpGcut_psi = 3.0
QpGcut_cou = 2.5
emax_sigm  = 2.0
```
In addition, we may use
```toml
[ham]
pwemax = 2
```
This is for the APWs in the band calculation.

We sometimes use this setting. It is better to check how the numerical results are affected. 
