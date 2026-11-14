# GWinput (legacy)

> English version of the page [manual/gwinput](../../manual/gwinput.md). The pages at the top of ecaljdoc are the reference; this one was made from manual/gwinput.md of commit 85cdee981 (source blob c626e5938cfc).

> ⚠️ The Fortran programs do not read the `GWinput` text file. The settings of this page are keys of `ctrlg.<sname>.toml` (sections `[gw]` / `[mlo]` / `[blocks]` / `[product_basis]`; the last one also holds the per-atom tables), written in TOML syntax (`key = value`). For a directory with the legacy `ctrl.<sname>` + `GWinput`, run `Legacy2toml.py <sname>` to convert. See [Input file](./toml_migration) for the key mapping.

The `[gw]` section of `ctrlg.<sname>.toml` sets the computatioanl conditions for  GW/QSGW calculations.

## Automatic generation by mkGWinput
- The GW settings (`[gw]`, `[mlo]`, `[blocks]`, `[product_basis]` of `ctrlg.<sname>.toml`) can be generated automatically, and most parameters can be used with their default settings.

When starting from `ctrls.<sname>` (`<sname>` is the material name),
```bash
ctrlgenToml.py <sname>
```
writes `ctrlg.<sname>.toml` including the GW settings. To add the GW settings to an existing `ctrlg.<sname>.toml`,
```bash
mkGWinput <sname>                # regenerate the GW settings; the original file is kept as ctrlg.<sname>.toml.bakup
ctrlgenToml.py <sname> --addgw   # add to a file without [gw] (stops if [gw] exists)
```
> [!INFO]
> Both internally run executables of the ecalj package (`lmfa`, `lmf --jobgw=0`, `gwinit`) with `mpirun`.
> The executables in the same directory as the script (the binary install directory of the ecalj package) are used.

The generated `[gw]` has `n1n2n3 = [4, 4, 4]`, `t_sigmaw = 300`, `t_tetrakbt = 300`. Rewrite them for your calculation
(`t_tetrakbt = 0` if you do not smear $\chi_0$ for insulators and semiconductors<span class="clfix">; the template value 300 does not affect the result either when the gap is much larger than $k_BT$ ([kBT](kBT) §4)</span>).


## Format
- Written in TOML syntax. From `#` to the end of the line is a comment
- The form is `key = value`; the order does not matter, but do not write the same key more than once. Logical values are `true` / `false`,
  lists are `[4, 4, 4]`, reals are `1e-5` (`1d-5` cannot be written), and multi-line tables (`QforGW` etc.) are enclosed in `"""`

## Parameters in GWinput
The meaning of the parameters is written in the generated `[gw]`; below we describe important variables and variables without a description there.

### `n1n2n3`
Number of k points in the $GW$ self-energy calculation. Points different from `[bz] nkabc` (the k points of `lmf`) can be specified (with `gwsc --mlo` it must be the same as `nkabc`; otherwise it stops). Note that the computational cost is proportional to the square of this number of k points.
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
Whether to keep the wave functions in memory. The wave functions are written to the files `CPHI` and `GEIG`. If they are comparable to or larger than the memory available to one MPI process, always set `false`. Set `false` when the number of $k$ points or orbitals is large.
- **type**: boolean
- **default**: true

### `zmel_batch_gb`
Batch size for building the matrix elements (zmel) at once, in GB (per MPI rank). About one third to one quarter of the memory available to one MPI process is a good value. It is used both in the self-energy calculation and in the calculation of $W$ (polarization function). For large systems a larger value is faster, but reduce it if the calculation crashes from lack of memory.
- **type** : float
- **default** : 0.4 (values below 0.4 on CPU and below 2.0 on GPU are raised to 0.4 and 2.0, respectively)

If the value is too small for the system, the calculation may stop with `sxcf_fal2_count.sc. Too small memory for nmbatch. Enlarge zmel_batch_gb in ctrlg.<sname>.toml [gw]`. In that case enlarge it.

### `zmel_max_size`
Old name of `zmel_batch_gb`. It is still read as `zmel_batch_gb`. `MEMnmbatch` in old inputs is not read (write `zmel_batch_gb`).

### `KeepPpb`
Whether to keep in memory, for all symmetry operations, the matrix elements between the product basis in the MT and the wave-function basis in the MT (the `ppb` variable). Usually the default `false` (then only the symmetry operations in use are kept) is fine. With `true` the memory usage grows by the number of symmetry operations.
- **type**: boolean
- **default**: false

### `t_tetrakbt`  *(required, 2026-09-28)*
Selects how $\chi_0$ is smeared, by the sign of a temperature (K). It stops if not written.
- `0`: no smearing (T=0 tetrahedron). <del class="clfix">Insulators and semiconductors.</del><span class="clfix">For insulators and semiconductors the template value 300 is also fine (no effect if the gap ≫ $k_BT$; [kBT](kBT) §4).</span>
- `T > 0`: finite-temperature tetrahedron. The occupations in $\chi_0$ become Fermi-Dirac at temperature $T$, and the Fermi level is `EFERMI_kbt` (also on the $\Sigma$ side). [kBT](kBT) §2.
- `-T`: keeps the T=0 tetrahedron and smears $\mathrm{Im}\chi_0$ along frequency with a Gaussian. The width is the same as Fermi-Dirac at temperature $T$ (standard deviation $\pi k_BT/\sqrt3$, $-1000$ K = 0.00574 Ha = 0.156 eV). Used when the plasmon pole of $W$ is rough in metals (the practical setting for LiTi₂O₄ is $-992.4$; <del class="clfix">"Practical settings" in [kBT](kBT)</del><span class="clfix">"Settings for LiTi₂O₄" at the top of [kBT](kBT)</span>). It also affects $W(i\omega)$ and the static $W(0)$.
The gwinit template writes 300.
- **type**: float (K)
- **default**: none (required)

Keys appearing in old inputs are handled as in Table 1. `ctrlg_update.py ctrlg.<sname>.toml` rewrites `SmearX0` (Ha, standard deviation of the
Gaussian for $\mathrm{Im}\chi_0$) into `t_tetrakbt = -T` of the same width (0.0057 Ha → $-992.4$ K), and adds `t_tetrakbt = 0` to files without
`t_tetrakbt`. The calculation does not change. `Legacy2toml.py` does the same rewriting during conversion. For the history of the keys, see Table 2 in [kBT](kBT) §2.5.

**Table 1**. Old keys of `[gw]` and how they are handled now

| Old key | Current handling |
|---|---|
| `SmearX0`, `SmearX0q0` | stops if written. Write `t_tetrakbt = -T` |
| `tetrakbt` (logical) | not read. Finite temperature is `t_tetrakbt > 0` |
| `esmr` (Ry, Gaussian) | if `t_sigmaw` is absent, read after conversion to `t_sigmaw` of the same standard deviation |
| `t_sigmakbt` (K) | if `t_sigmaw` is absent, its value becomes `t_sigmaw` as is |
| `GaussianFilterX0`, `GaussSmear`, `delta`, `dw`, `omg_c`, `WgtQ0P`, `MEMnmbatch` | not read (no effect even if written) |
| `zmel_max_size` | read as `zmel_batch_gb` |

### `t_sigmaw`  *(2026-09-20, formerly `esmr` / `t_sigmakbt`)*
Width of the Fermi-Dirac kernel that smears the intermediate levels (poles of G) in the self-energy <del class="clfix">$\Sigma = G(W)$</del><span class="clfix">$\Sigma = iGW$</span>. Unit K.
The same kernel $-\partial f/\partial\varepsilon$ ($k_BT$ = `t_sigmaw`) is used for all of: the occupation weights of <del class="clfix">$\Sigma_x = Gv$</del><span class="clfix">$\Sigma_x = iGv$</span>, and the window weights and pole positions (`wcsmear`) of the real-axis pole term of <del class="clfix">$\Sigma_c = G(W-v)$</del><span class="clfix">$\Sigma_c = iG(W-v)$</span>.
**This is a numerical width, not a complete finite-temperature self-energy** (Bose factors etc. are not included).
The Fermi level is `EFERMI_kbt` if `t_tetrakbt > 0` ($\chi_0$ also at finite temperature); otherwise it is determined by counting levels on the k mesh (Table 5 in [kBT](kBT) §0).
The candidate window of intermediate states and the margin of the real-axis $\omega$ mesh are $E_F \pm 15\,k_BT$. The same kernel also enters the
half-residue step of the imaginary-axis integral of $\Sigma_c$ (contract: the occupation weights, the window weights of the pole term and the imaginary-axis step, all 3 use the same kernel. [manual/kBT](kBT) §3.6).

- The default 1000 K is a width close to the old `esmr = 0.01` Ry (Gaussian $\sigma$) (the same standard deviation is 872 K; the std of FD is
  $\pi k_BT/\sqrt3 = 1.81\,k_BT$). The old `esmr = 0.003` Ry (the value for Si etc. in TestInstall) corresponds to 262 K.
- Insulators: if the width is sufficiently smaller than the band gap, it does not affect the result (for Si the difference between 262 K and 1000 K is a few meV).
- Metals: in principle the limit $n_1n_2n_3 \to \infty$, `t_sigmaw` $\to 0$, but a finite width is needed for convergence.
  300–1000 K together with a fine k mesh. Choose per material (1000 K for LiTi₂O₄).
- If `esmr` remains in ctrlg, it is read after conversion `t_sigmaw = esmr[Ry]·13.6057/(1.81·8.617e-5)`
  (`t_sigmaw` takes precedence if present). `t_sigmakbt` is copied as is to `t_sigmaw`.
- **type**: float (K)
- **default**: 1000.0 (the gwinit template writes 300)

### `wcsmear`  *(2026-09-19, default true from 2026-09-20)*
In the real-axis pole term of $\Sigma_c$, the Fermi-Dirac kernel of the intermediate levels (`t_sigmaw`) is also applied to the pole position $\omega_\varepsilon$,
and $W_c(\omega)$ is integrated with the kernel. `false` gives the conventional one-point evaluation at the mean energy (3-point Lagrange interpolation,
PRB76,165106 Eq.58). In systems where $W_c$ is smooth both are almost the same (a few meV up to 77 meV at most for Si).
Required for metals with sharp plasmon poles (LiTi₂O₄ ≲ 1500 K). The cost of the real-axis pole term grows several times (+7–26 % in total).
- **type**: boolean
- **default**: true

> Parallelization of the tetrahedron weights (`tetwt5`): from 2026-09-22 **MPI k parallelization** (each rank of the W-build computes only its own k range).
> Launching hgw with `-np2` with more ranks than GPUs (GPUs shared round-robin) makes it correspondingly faster. It matters when χ₀ is at finite temperature
> (`t_tetrakbt > 0`, 150 s per q for 9³). The OpenMP (`omp_tetwt`) tried 09-19 to 22 was removed.

### `unit_2pioa`
  `unit_2pioa = false`   # false --> a.u.; true --> unit of QpGcut_* are in 2*pi/alat

  Usually `false` is fine.

###  `alpha_OffG`
  `alpha_OffG = 1.0`   # (a.u.) Used in auxially function in the offset-Gamma method.

  Not need to change.

### `QpGcut_psi`

Cutoff ~~energy: unit Ryd~~ of the IPW basis of the wave functions. Note that the basis functions of the wave functions differ from the PMT method used in the `lmf` calculation. The wave functions in the interstitial region are re-expanded and described in IPW. IPW means interstitial plane wave.
- **type**: float
- **default** : 4.0
`|q+G| < QpGcut_psi` in units of 1/a.u.
The number is written in the output of qg4gw (lqg4gw).
To convert to energy (Ry), take **2.

The IPW cutoff of the wave functions is also related to `[ham] pwemax`. For pwemax, 2 or 3 is often tried
<del class="clfix">(its unit is an energy cutoff, i.e. **2, in Ryd).</del><span class="clfix">(pwemax is an energy cutoff in Ry: $|\mathbf q+\mathbf G|^2 \lt {\rm pwemax}$, with $|\mathbf q+\mathbf G|$ in 1/a.u.)</span>

### `QpGcut_cou`
Cutoff ~~energy: unit Ryd~~ of the IPW part of the product basis 
`|q+G| < QpGcut_cou` in units of 1/a.u.
To convert to energy (Ry), take **2.
- **type**: float
- **default** : 3.0

The number is written in the output of qg4gw (lqg4gw).

### `emax_sigm`
Maximum energy, measured from the Fermi level, of the quasiparticle states whose self-energy is calculated: unit Ryd. In QSGW the self-energy is used to build the potential, so making it small also lowers the accuracy of the potential.
- **type**: float
- **default** : 3.0


## Change `HistBin` for high resolution energy mesh near Ef for metal 
This defines histogram bins to accumulate weight (imaginary part) of the polarization functions.

For metal, it might be better to test
```toml
HistBin_dw    = 1e-5   # 1e-5 is fine mesh (good for metal?)  (a.u.) BinWidth along real axis at omega=0.
HistBin_ratio = 1.03   # 1.03 maybe safer. frhis(iw)= b*(exp(a*(iw-1))-1), where a=ratio-1.0 and dw=b*a
```
This mesh is made in getfreq of m_freq.f90.

* HistBin_dw and HistBin_ratio specify <del class="clfix">real space bins</del><span class="clfix">bins along the real frequency axis</span> which we accumulate imaginary part weight of polarization functions. The bins are (see frhis in m_freq.f90) 
<del class="clfix">$\omega_i = {{dw}} ∗ (\exp({\rm ratio} ∗ (i − 1)) − 1)$</del><span class="clfix">$\omega_i = b\,(e^{a(i-1)}-1)$, $a={\rm HistBin\_ratio}-1$, $b={\rm HistBin\_dw}/a$</span>
where $[\omega_i, \omega_{i+1}]$ $(i = 1, 2, ...,\verb#nwhis#+1)$ is the i-th bin. HistBin_dw is bin width at ω = 0.
The ratio ω(i + 1)/ω(i) for large ω is <del class="clfix">exp(a)=HistBin_ratio</del><span class="clfix">$e^{a}\approx$ HistBin_ratio</span>. This choice of getting coarser at high energy is because we think W (ω) around ω ∼ 0 gives most important contribution to the GW approximation. If histogram bins are too wide, dielectric function can be less accurate, but results may be not so much affected.
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



## `Smaller lcutmx to reduce computatinal time`
  To reduce the computational time, we reduce number of MPB (mixed product basis).
  One is lcutoff of Product Basis (PB) within MT. Use lcutmx=2 for oxygen or something (s,p block atoms). 
  Thus it  is like
```toml
[product_basis]
pb_lcutmx = [4, 4, 4, 2, 2, 2]   # maximum l-cutoff for the product basis, one per atom. 4 is required for atoms with valence d
```
in `ctrlg.<sname>.toml`. The atoms are in the order of the `[[site]]` tables (`lmchk` shows it, also in `SiteInfo.lmchk`).

* NOTE: we need lcutmx =6 is requied for 4f/5f atoms, while 2 is enough for Oxygen (probably other species on the same row).


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
