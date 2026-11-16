# kBT — finite temperature and smearing (`t_tetrakbt` / `t_sigmaw` / `wcsmear`)

> English version of the page [manual/kBT](../../manual/kBT.md). The pages at the top of ecaljdoc are the reference; this one was made from manual/kBT.md of commit 4f74d9132 (source blob 6db99d0d6c18).

This page explains the `[gw]` keys that set the electronic temperature and the numerical smearing in GW and QSGW. There are separate keys for the $W$ side ($\chi_0$) and the $\Sigma$ side ($G$) (Table 4).
`t_sigmaw` on the $\Sigma$ side is a numerical width that is always in effect; `t_tetrakbt` on the $\chi_0$ side is required: if positive it is the physical electronic temperature,
if negative it is the width (written as a temperature) of the Gaussian that smears $\mathrm{Im}\chi_0$, and if 0 there is no smearing.

This page describes the present specification. For the history of the keys, needed when reading old inputs and records, see Table 2 in §3.5. The course of events in 2026-06 to 09 and the explanation of the calculations of that time are in
the development record (MD/kBT/history.md in the repository, not on the site).

**Table 4**. Keys of `[gw]` (temperatures in K)

| Key | Where it acts | Default |
|---|---|---|
| `t_tetrakbt` (required) | Selects the smearing of $\chi_0$ by its sign (0: no smearing, $T>0$: finite-temperature tetrahedron, $-T$: Gaussian. Table 1 in §3) | None (stops if not given) |
| `t_sigmaw` (required) | Width of the Fermi-Dirac kernel that smears the intermediate levels in $\Sigma_x=Gv$ and $\Sigma_c=G(W-v)$ (§4) | None (required since 2026-11-14. The default of the code before that was 1000) |
| `wcsmear` | In the real-axis pole term of $\Sigma_c$, integrates $W_c$ using the same kernel also for the position of the pole (§5) | true |
| `chi0_skip_window` = [emin, emax] (eV) | Removes from $\chi_0$ the band pairs whose both ends fall in this window (cRPA-like. Under test) | off |

The template of gwinit (also used by `ctrlgenToml.py`) writes `t_tetrakbt = 300` and `t_sigmaw = 300`.

**Which to use**

- Insulators and semiconductors: the template as it is (§7).
- To see the effect of the electronic temperature: `t_tetrakbt = t_sigmaw = T` (Si in §8.1, Fe in §8.2).
- Metals whose bands become rough depending on the iteration or the k mesh (plasmon poles of $W$, §5): `wcsmear` (default) and `t_sigmaw = 1000`.
  To smear $\chi_0$ while keeping T=0, make `t_tetrakbt` negative.

**Setting for LiTi₂O₄**: $\chi_0$ is the T=0 tetrahedron with $\mathrm{Im}\chi_0$ smeared by a Gaussian, `t_tetrakbt = -992.4`
(the same width as Fermi-Dirac at 992.4 K, standard deviation 0.155 eV), `t_sigmaw = 1000`, `wcsmear = true`, `mixbeta = 0.5`.
With this, $6^3$ converges in 8 iterations (the change of the QP energies from iteration 7 to 8 is at most 10 meV, 2 meV on average), and the one-shot difference between $6^3$ and $9^3$ in `dSEnoZ` (the column of `QPU`: the diagonal elements of $\langle\Sigma-V_{\rm xc}\rangle$ without multiplying by $Z$. [gwsc](./gwsc)) is at most 13 meV.
The MLO-QSGW calculation (§8.3, ecalj `Samples/kBT/LiTi2O4`) uses the same setting.

---

## What changes between T=0 and finite temperature (overview)

The temperatures on the $\chi_0$ side and the $\Sigma$ side are independent keys. If you use `t_tetrakbt > 0`, it is natural to set the same T in `t_sigmaw` on the $\Sigma$ side
(not required). The setting for LiTi₂O₄ keeps $\chi_0$ at T=0 and smears it with the Gaussian of `t_tetrakbt < 0` ("Setting for LiTi₂O₄" at the top).

### $\chi_0$ side (hx0fp0 / W-build of hgw)

The differences between the three cases of `t_tetrakbt` (0, $T>0$, $-T$) are in Table 1 and Eqs. (3) to (8) of §3. In all cases `chi0_skip_window = [emin, emax]`
(eV; removes from $\chi_0$ the band pairs whose both ends fall in the window. Under test) can be used.

### $\Sigma$ side (hsfp0_sc / $\Sigma_c$ of hgw, $\Sigma_x$ of hsfp0)

On the $\Sigma$ side there is no distinction between "T=0 and finite temperature". The kernel is always Fermi-Dirac, with the width `t_sigmaw`.
The only thing tied to the temperature on the $\chi_0$ side is the Fermi level.

**Table 5**. What on the $\Sigma$ side is and is not tied to the temperature of $\chi_0$

| | `t_tetrakbt` ≤ 0 | `t_tetrakbt > 0` |
|---|---|---|
| Fermi level | Determined by counting the levels on the k mesh smeared by a Gaussian of width $k_BT$ (`t_sigmaw`) (`efsimplef2ax`. Determined separately from the `EFERMI` used by $\chi_0$) | `EFERMI_kbt` (stops if absent) |
| Occupation kernel of the levels (weight of $\Sigma_x$, window weight of the pole term) | Fermi-Dirac, width $k_BT$ = `t_sigmaw` | Same as left |
| Position at which the pole term picks up $W_c$ | `wcsmear = true` (default): integrates $W_c$ with the FD kernel. `false`: the mean $\bar\omega$ (3-point interpolation) | Same as left |
| Candidate window of intermediate states, margin of the $\omega$ mesh | $E_F \pm 15\,k_BT$ (`sig_window`) | Same as left |
| Imaginary-axis integral, $W_c$ itself, mixing of sigm (`mixbeta`) | Common | Common |

- `t_sigmaw` is a numerical width, not a full finite-temperature self-energy (Bose term, §11.7).
  The standard deviation of the kernel $-\partial f/\partial\varepsilon$ is $\pi k_BT/\sqrt3 \simeq 1.81\,k_BT$ (0.156 eV at 1000 K).
- `wcsmear` is not a temperature switch but a switch to "use the same distribution as the occupation kernel also for the position of the pole" (§5).
- On the $\chi_0$ side the $\theta$ is truly sharp at T=0, so a width appears only when `t_tetrakbt` is made nonzero.

---

## Why it is needed

There are two uses.

1. **Physics of the electronic temperature**. At $T>0$ the occupation near the Fermi level becomes the Fermi-Dirac distribution, and in semiconductors carriers thermally excited across the gap
   join the screening. This is handled by `t_tetrakbt = T` on the $\chi_0$ side (§3.3) and `t_sigmaw = T` on the $\Sigma$ side (§4).
   In Si the gap shrinks by 0.16 eV at an electronic temperature of 3000 K (§8.1).
2. **Numerical stabilization of QSGW for metals**. In metals $W_c(\omega)$ at small $q$ has sharp plasmon poles, and if the real-axis pole term of $\Sigma_c$
   picks them up at a single point, the value changes strongly with the k mesh and at every iteration (§5). What smears this are `wcsmear` and `t_sigmaw` on the $\Sigma$ side,
   and `t_tetrakbt < 0` on the $\chi_0$ side (Gaussian, §3.4). The temperature here is a numerical width, not a physical temperature.

---

## $\chi_0$ side — `t_tetrakbt`

`[gw] t_tetrakbt` (K) is **required** (stops if not given). Its sign selects the smearing of $\chi_0$ (Table 1).
The template of gwinit writes `t_tetrakbt = 300`. The name and meaning of the key changed several times, so when reading old inputs and records see Table 2 in §3.5.

**Table 1**. The three cases of `t_tetrakbt`

| `t_tetrakbt` | Occupation difference $f_a-f_b$ | $\mathrm{Im}\chi_0$ in the $\omega$ direction | Fermi level of $\chi_0$ ($\Sigma$: Table 5 in §1) | Use |
|---|---|---|---|---|
| 0 | Sharp step, Eq. (3) | As it is | `EFERMI` | Insulators and semiconductors |
| $T>0$ | Fermi-Dirac at temperature $T$, Eqs. (4)(5) | As it is | `EFERMI_kbt`, Eq. (6) | Physics of the electronic temperature (§2, §8.1) |
| $-T$ | Sharp step, Eq. (3) | Smeared by a Gaussian, Eqs. (7)(8) | `EFERMI` | Metals where the plasmon poles of $W$ make things rough (§5) |

```toml
[gw]
t_tetrakbt = 0        # (K) required. 0: no smearing. T > 0: finite-temperature tetrahedron at temperature T. -T: smear Im chi0 with a Gaussian
```

### Common — histogram of $\mathrm{Im}\chi_0$ and the Hilbert transform

In all cases, $\chi_0$ is built by first making $\mathrm{Im}\chi_0$ as a **histogram** for $\omega>0$, and from it the real part on the real axis and
the values on the imaginary axis are made by the Hilbert transform. Bin $j$ is $[\omega_j,\omega_{j+1}]$; its center is written $\bar\omega_j$ and its width $\Delta_j=\omega_{j+1}-\omega_j$
(`frhis`. An exponential mesh determined by `HistBin_dw` and `HistBin_ratio`, wider at larger $\omega$). The weight of a bin, accumulated by the tetrahedron method over the transitions from a level $e_a$ at $\mathbf k$
to a level $e_b$ at $\mathbf k+\mathbf q$, is

$$
h_j \propto \sum_{\mathbf k}\sum_{a,b}\int_{\omega_j}^{\omega_{j+1}}\!d\omega\;
\bigl|M_{ab}\bigr|^2\,\bigl(f_a-f_b\bigr)\,\delta\bigl(\omega-(e_b-e_a)\bigr),
\qquad f_a=f(e_a)
\tag{1}
$$

($M_{ab}$ is the matrix element with the product basis; factors and the sum over spins are omitted). The tetrahedron method linearly interpolates the levels and matrix elements
inside each tetrahedron, and within that approximation integrates the $\delta$ function over the bin exactly. $h_j/\Delta_j$ is the mean of $\mathrm{Im}\chi_0$ in bin $j$.
Since $\mathrm{Im}\chi_0$ is an odd function of $\omega$ ($\mathrm{Im}\chi_0(-\omega)=-\mathrm{Im}\chi_0(\omega)$), only $\omega>0$ is kept, and
the bin means are connected piecewise linearly between the centers, and

$$
\chi_0(z)=\frac1\pi\int_{-\infty}^{\infty}\!d\omega'\,\frac{\mathrm{Im}\chi_0(\omega')}{\omega'-z}
\tag{2}
$$

is integrated analytically bin by bin (`hilbertmat`, `dpsion_init` and `dpsion_chiq_*` in `dpsion5.f90`). On the real axis, the real part at $z=\omega+i0$
is $\mathrm{Re}\chi_0(\omega)$, and for the imaginary part the bin means are used as they are. On the imaginary axis $z=i\nu$, and since the function is odd the values are real.
The three cases differ only in the occupation difference $f_a-f_b$ of Eq. (1) and in the treatment of $h_j$ before it is passed to Eq. (2).

### `t_tetrakbt = 0` — T=0 tetrahedron

The occupation is a step at the Fermi level $E_F$ (`EFERMI`), and the occupation difference in Eq. (1) is

$$
f_a-f_b\;\to\;\theta(E_F-e_a)\,\theta(e_b-E_F)
\tag{3}
$$

(`lindtet6`. Transitions from an occupied $a$ to an empty $b$, falling at $\omega=e_b-e_a>0$).

### `t_tetrakbt = T > 0` — finite-temperature tetrahedron (method B′)

Method B′ is the name used during development; it refers to the method that averages the weights of the T=0 tetrahedron over the Fermi level with the thermal kernel, as in Eq. (5).

The finite-temperature occupation difference equals the step of Eq. (3) averaged over the Fermi level $E$ with the thermal kernel $-f'$. With $f(x)=1/(e^{x/k_BT}+1)$,

$$
\int dE\,\bigl(-f'(E-\mu)\bigr)\,\theta(E-e_a)\,\theta(e_b-E)=f(e_a-\mu)-f(e_b-\mu)\qquad(e_a<e_b)
\tag{4}
$$

Since Eq. (1) is linear in $f_a-f_b$, the finite-temperature bin weight is the T=0 weight $h^{(0)}_j(E)$ with the Fermi level placed at $E$,
averaged over $E$. Putting $E=\mu+2k_BT\,t$ gives $(-f')\,dE=\tfrac12\,\mathrm{sech}^2t\,dt$, and this is integrated with 20 points:

$$
h_j(T)=\sum_{i=1}^{20}c_i\,h^{(0)}_j\bigl(\mu+2k_BT\,t_i\bigr),\qquad
c_i=\frac{w_i\,\tfrac12\mathrm{sech}^2t_i}{\sum_{i'}w_{i'}\,\tfrac12\mathrm{sech}^2t_{i'}}
\tag{5}
$$

The nodes $t_i$ and weights $w_i$ are 5-point Gauss-Legendre (`gausq_fd`) on each of the 4 intervals $[-6,-1],[-1,0],[0,1],[1,6]$ into which
$t\in[-6,6]$ (that is, $|E-\mu|\le12\,k_BT$) is divided. The sum of the $c_i$ is normalized to 1, so the sum of the occupations is kept even though the outside of the range is cut off.
When at the 4 vertices of a tetrahedron both $a$ and $b$ lie together outside the window $|e-\mu|\le12\,k_BT$ on one side (the occupation is flat over the window), the result is the T=0 value
itself, so the 20 evaluations are skipped (`lindtet6_kbt`, `tetwt5.f90`). Eq. (5) only exchanges the order of the $\mathbf k$ integral and the $E$ integral,
and the only error is that of the numerical $E$ integration (§10.1, §10.2). The numerator obtained is the Adler–Wiser $f_a-f_b$, not $f_a(1-f_b)$
(this is the correct one for a response function. What is dropped on the $\Sigma$ side is in §10.4).

$\mu$ = `EFERMI_kbt` is obtained by `heftet`, which smears the tetrahedron number of states $N(E)$ with the same thermal kernel and

$$
\sum_{i=1}^{20}c_i\,N\bigl(\mu+2k_BT\,t_i\bigr)=N_{\rm val}
\tag{6}
$$

solves this by bisection, writing it to the file `EFERMI_kbt` (`fermi_kbt_tetra`, `main_heftet.f90`). As $T\to0$ it returns to the tetrahedron $E_F$.
When the solution is an interval in the gap (the left-hand side is flat at $N_{\rm val}$ across the gap), the middle of the interval is taken (item 2 of §12).
When `t_tetrakbt > 0`, both $\chi_0$ and $\Sigma$ use this, and the program stops if the file is absent (§1, §4). When `t_tetrakbt` ≤ 0, $\chi_0$ uses `EFERMI` of the tetrahedron method of `heftet`, and $\Sigma$ uses the Fermi level counted on the k mesh (Table 5 in §1).

Check in the output (standard output of `heftet`):

```
tetrakbt_init: T[K], kbt[Ry], kbt[eV]  2000  0.12667E-01 0.17234E+00
```

### `t_tetrakbt = -T` — smearing $\mathrm{Im}\chi_0$ with a Gaussian

The occupation and the Fermi level stay at T=0 (Eq. (3), `EFERMI`), and the bin weights $h_j$ are smeared along $\omega$ with a Gaussian before being passed to Eq. (2).
The standard deviation of the Gaussian is matched to that of the Fermi-Dirac kernel $-f'$ at temperature $T$:

$$
\sigma=\Bigl[\int dx\,x^2\bigl(-f'(x)\bigr)\Bigr]^{1/2}=\frac{\pi k_BT}{\sqrt3}\simeq1.81\,k_BT
\tag{7}
$$

($\sigma=0.1563$ eV $=0.005744$ Ha at $T=1000$ K). The smearing is

$$
\tilde h_i=\sum_j G_{ij}\,h_j,\qquad
G_{ij}=\frac{\bigl[K(\bar\omega_i-\bar\omega_j)-K(\bar\omega_i+\bar\omega_j)\bigr]\,\Delta_i}{\sum_{i'}K(\bar\omega_{i'}-\bar\omega_j)\,\Delta_{i'}},\qquad
K(x)=e^{-x^2/2\sigma^2}
\tag{8}
$$

where all sums run over the bins with $\omega>0$ (`gaussianfilterhis` and `smearx0_apply`, `dpsion5.f90`). Column $j$ redistributes the weight of bin $j$
to the surrounding bins with the Gaussian:

- $\Delta_i$ in the numerator (width of the destination bin): even on an exponential mesh with varying widths, this amounts to convolving the mean $h/\Delta$ with the Gaussian
- The second term $-K(\bar\omega_i+\bar\omega_j)$ in the numerator: extension to an odd function. The part that spills into negative $\omega$ is folded back with its sign changed.
  $\tilde h$ goes to 0 as $\omega\to0$
- Denominator: the sum of the symmetric kernel over $\omega>0$. For $\bar\omega_j\gg\sigma$ it becomes $\sqrt{2\pi}\,\sigma$, and Eq. (8) becomes the continuous form

$$
\widetilde{\mathrm{Im}\chi_0}(\omega)=\int_0^\infty\!d\omega'\,\bigl[g_\sigma(\omega-\omega')-g_\sigma(\omega+\omega')\bigr]\,\mathrm{Im}\chi_0(\omega'),
\qquad g_\sigma(x)=\frac{e^{-x^2/2\sigma^2}}{\sqrt{2\pi}\,\sigma}
\tag{9}
$$

  that is, the convolution of $\mathrm{Im}\chi_0$ extended as an odd function with $g_\sigma$. Eq. (9) conserves the f-sum $\int_0^\infty\omega\,\mathrm{Im}\chi_0\,d\omega$
  (the sum of the weights over $\omega>0$ decreases by the part folded back near 0). For $\bar\omega_j\lesssim\sigma$ only part of the Gaussian lies at $\omega>0$,
  so the denominator is smaller than $\sqrt{2\pi}\,\sigma$ and Eq. (8) deviates from Eq. (9)
- When both $\pm\omega$ are kept (`npm = 2`, e.g. $\chi^{+-}$ for magnons): the symmetric kernel without folding is applied separately to $\omega>0$ and $\omega<0$

Since the real part and the imaginary-axis values are made by Eq. (2) from the smeared $\tilde h_j$, the real-axis $W(\omega)$, the imaginary-axis $W(i\nu)$, and the static $W(0)$ all come from the same
smeared $\mathrm{Im}\chi_0$.

$-T$ is a convention to write the Gaussian width as a temperature, and is different from the finite-temperature $\chi_0$ (Eq. (5)). Eq. (5) smears the occupations of the levels,
while Eq. (8) smears the distribution of the transition energies $e_b-e_a$.

The input value is converted by `m_GWinput` to $\sigma$ (Ha; the internal variable name is `SmearX0`), and the following line appears in the output:

```
 SmearX0 (Ha) =   5.7438E-03  : Gaussian smearing of Im chi0 along omega ([gw] t_tetrakbt < 0)
```

### History of the keys

**Table 2**. History of the inputs for the smearing of $\chi_0$ and $\Sigma$

| Period | Smearing of $\chi_0$ | Width of the levels in $\Sigma$ (§4) |
|---|---|---|
| to 2026-06 | `GaussianFilterX0` (GWinput, Gaussian on $\mathrm{Im}\chi_0$) | `esmr` (Ry, Gaussian) |
| 2026-06 to 09-18 | Finite temperature by the logical key `tetrakbt = true` and `t_tetrakbt` (K). Gaussian by `SmearX0` (Ha, $\sigma$; renamed from `GaussianFilterX0`) | `esmr`, or `t_sigmakbt` (K, Fermi-Dirac) |
| 09-19 to 09-27 | Finite temperature by `t_tetrakbt > 0` alone (default 0 = no smearing). `SmearX0` exclusive with `t_tetrakbt > 0` from 09-20 | `t_sigmaw` (K, Fermi-Dirac) from 09-20. `esmr` and `t_sigmakbt` are read with conversion |
| 09-28 to 11-13 | `t_tetrakbt` is required, with three cases by its sign (Table 1). `SmearX0` is not read (stops if given) | `t_sigmaw` (if absent, `esmr` or `t_sigmakbt` is read with conversion; if those are absent too, 1000 K) |
| 11-14 onward | Same as above | `t_sigmaw` is required. `esmr` and `t_sigmakbt` are not read (`ctrlg_update.py` converts them and writes `t_sigmaw`) |

- `SmearX0 = s` and `t_tetrakbt = -T` give the same calculation when $s\,[\mathrm{Ha}]=\pi k_BT/\sqrt3$ by Eq. (7) (0.0057 Ha and $-992.4$ K).
  `ctrlg_update.py ctrlg.<sname>.toml` rewrites with this conversion, and adds 0 (the former default) to files without `t_tetrakbt`
- `esmr` = $e$ (Ry, standard deviation of the Gaussian) and `t_sigmaw` = $T$ give the same width when $e=\pi k_BT/\sqrt3$ (0.003 Ry and 262 K, 0.01 Ry and 872 K. These values arise because the code rounds $\pi/\sqrt3$ to 1.81 in the conversion; without rounding they are 261 K and 870 K).
  From 2026-11-14, `esmr` and `t_sigmakbt` are not read. `ctrlg_update.py` adds `t_sigmaw`, converting `esmr` with this relation and taking `t_sigmakbt` as it is
- The template of gwinit wrote `tetrakbt = true` and `t_tetrakbt = 300` from 2026-06 to 09-16 (removed on 09-17). Since 2026-09-28, when `t_tetrakbt` became required, it again writes `t_tetrakbt = 300` (`main_gwinit.f90`).
  From 09-28 it writes `t_tetrakbt = 300` and `t_sigmaw = 300`

---

## $\Sigma$ side — `t_sigmaw`

The intermediate levels of $\Sigma$ (the poles of $G$) are smeared by a Fermi-Dirac kernel. Its width is `t_sigmaw` (K), and
the same kernel is used for all of the occupation weight of $\Sigma_x = Gv$, the window weight of the real-axis pole term of $\Sigma_c = G(W-v)$,
and the distribution of the pole position in `wcsmear`.

```toml
[gw]
t_sigmaw   = 1000.0   # (K) width of the FD kernel that smears the levels of Sigma (required. The default up to 2026-11-13 was 1000)
t_tetrakbt = 1000.0   # if chi0 is also to be at finite temperature. The Fermi level of Sigma becomes EFERMI_kbt
```

The window weight is the Fermi–Dirac cumulative distribution

$$
\Phi(x)=\frac{1}{1+e^{-x/k_BT}},
\qquad
w(e_l,e_h;\varepsilon_k)=\Phi(e_h-\varepsilon_k)-\Phi(e_l-\varepsilon_k)
\tag{10}
$$

**The same kernel also enters the imaginary-axis integral of $\Sigma_c$**
— the pole term and the imaginary-axis integral are the two halves of one contour, and unless the smearing of the levels enters both with the same kernel,
$\Sigma_c$ breaks. Since this is an important point, it is described separately in [§6](#kbt-contour-smearing).
The candidate window of intermediate states and the margin of the real-axis $\omega$ mesh
are the tail of the kernel $E_F \pm 15\,k_BT$ (`sig_window`), and widen automatically when `t_sigmaw` is raised (§10.3).

### Choosing the temperature

- Insulators: no effect if sufficiently smaller than the gap (in earlier records, the difference between 262 K and 1000 K in Si was a few meV. Not checked with the present code).
- Metals: in principle the limit $n_1n_2n_3 \to \infty$, `t_sigmaw` → 0, but a width comparable to the level spacing on the Fermi surface
  of the k mesh is needed. 300 to 1000 K, per material. 1000 K for LiTi₂O₄ $9^3$.
- When `t_tetrakbt` is used, it is natural to set the same temperature (not required).

### Check in the output

```
    t_sigmaw =   1000.0 K  (Sigma kernel width kBT =    0.006333 Ry)
 sigmaw_setup: Sigma kernel Fermi-Dirac, t_sigmaw[K]=  1000.0  ef: counted on the k mesh (efsimplef2ax); chi0 at T=0
```

If `t_tetrakbt > 0`, `ef<-EFERMI_kbt=...`. In this case the program aborts if `EFERMI_kbt` (written by `heftet`) is absent.

---

## What makes metals rough, and `wcsmear`

In QSGW for metals, the bands can become rough by tens of meV depending on the k mesh and the iteration, and a one-shot GW from LDA can produce spikes of the size of eV
near $E_F$ (for LiTi₂O₄ at 1000 K or below, the bottom of the O 2p band is rough by 30 to 100 meV at $9^3$, and one-shot GW gives spikes of −2.5 to −3.4 eV ($6^3$, 500 K and 1000 K)). The cause is as follows
(the record of the investigation is in the ecalj development record (MD/research_kotani_log.md in the repository, not on the site), 2026-09-18 to 19).

- **The culprit is the real-axis pole term of $\Sigma_c$** (cutting it changes the anomalous off-diagonal ⟨48|Σc|33⟩ from −13 to −0.35 eV. ⟨48|Σc|33⟩ is the off-diagonal element $\langle 48|\Sigma_c(\varepsilon_{48})|33\rangle$ between band 48 ($E_F$+3.7 eV) and band 33 at $\mathbf q=(-\tfrac16,\tfrac16,\tfrac16)$ in one-shot GW from LDA for LiTi₂O₄ $6^3$ at 1000 K. The same below).
  The imaginary-axis integral, the mixing of sigm, the head at offset-Γ, and the $\omega$ interpolation method are unrelated.
- What the pole term hits is the **t2g plasmon pole** in $W_c(q,\omega)$ at $q$ of the first shell
  ($\omega_p \simeq 1.78$ eV, width ≈ 0.1 eV. Checked directly with `--dumpW`). The pole term, for each intermediate state
  $n'$, picks up $W_c(\mathbf k,\, \varepsilon_{\mathbf q n} - \varepsilon_{\mathbf q-\mathbf k, n'})$ at
  a single point, so a term whose $\varepsilon$ falls within ±0.1 eV of the pole receives a value of $W_c$ of about ±50 a.u. (at the pole $\mathrm{Re}\,W_c$ jumps from −54 to +4 a.u. to the next bin),
  and crosses to the other side of the pole whenever the eigenvalues move in the iteration (seesaw).
- It vanishes monotonically with temperature: ⟨48|Σc|33⟩ = −27 / −13 / −3.2 / −1.9 eV (500 / 1000 / 1500 / 2000 K).
  This is why the runs at 2000 K and above had no trouble.

**Remedy `wcsmear = true`** (default): the approach of PRB 76, 165106 Eq. 58 (`wcsmear = false`) uses the smearing kernel $g$ of the intermediate levels (the Fermi-Dirac
$-\partial f/\partial\varepsilon$, width `t_sigmaw`) only for the window weight $w_{n'} = \int g$, and evaluates $W_c$ at a single point,
the mean energy $\bar\omega$. As a discretization of the $k$ integral it is consistent to use the same kernel for the weight and for the position of the pole, and

$$
\Sigma^{\rm pole}_{nn''} = \sum_{\mathbf k n'} s\,M^{*}_{n'n}
\left[\int_{e_l}^{e_h} de\; g(e-\varepsilon')\,W_c(\mathbf k,\,\omega-e)\right] M_{n'n''}
\tag{11}
$$

the $W_c$ side is integrated (Eq. (11)). If $W_c$ is smooth within the width of the kernel this agrees with the single-point evaluation (at most 77 meV in Si,
a few meV for the low-lying states), and the pole is smeared by the width of the kernel (Table 6).

**Table 6**. Effect of `wcsmear`

| LiTi₂O₄ 1000 K | `wcsmear = false` | `wcsmear = true` |
|---|---|---|
| $6^3$ one-shot, spike on Γ–L | −3.4 eV | Disappears |
| $6^3$ one-shot, roughness of O 2p | 12.2 meV | 7.2 meV |
| $9^3$ iter 2→3, roughness of O 2p (no mixing) | 47.8 meV | 12.7 meV |
| Cost of the real-axis pole term | — | ×6.7 (overall $6^3$ +7 %, $9^3$ +26 %) |

```toml
[gw]
t_sigmaw   = 1000.0   # default value up to 2026-11-13. wcsmear = true is the default
t_tetrakbt = 1000.0   # if the chi0 side is also to be at finite temperature
```

Notes: (1) Conceptually, this carries the approximation "treat a pole of G as a spectral function of width $k_BT$" through the whole pole term,
while the edges of the contour (E_F and ω) remain sharp. (2) Using the Fermi-Dirac $-\partial f/\partial\varepsilon$ at finite temperature
as the distribution of the pole position is a heuristic regularization (the exact finite-temperature formula only replaces θ → f with the levels kept sharp).
The proper width should come from the $k$ integral (tetrahedron). (3) It is implemented in the QSGW path (`hsfp0_sc` / `hgw`)
and in the pole term of one-shot GW `hsfp0` (`sxcf_fal2`, `gw_lmfh`).

**Limitation**: `wcsmear` smears the roughness of the states that "pick up the pole across it" (O 2p, the first-shot spikes), but
it cannot remove the off-diagonal $\tfrac12[\Sigma(\varepsilon_i)+\Sigma(\varepsilon_j)]_{ij}$ of **pairs with a plasmon pole between the two levels**, such as a level just below $E_F$ and a level at +3 eV
(for a pair spanning the interval where $\Sigma(\omega)$ changes abruptly,
the premise of the static form breaks down). For LiTi₂O₄ $9^3$ at 1000 K, spikes of this type remain near $E_F$
from iteration 2 on. As long as the RPA $W$ has sharp undamped poles, this is a limit of applicability of static QSGW;
they do not appear at 1500 K and above (the poles are Landau damped).

`wcsmear = false` returns to the single-point evaluation (differences of a few meV in insulators and tens of meV in metals). A method that makes the weights of the pole term by the same tetrahedron
integration as $\chi_0$ is not implemented.

---

## Contour decomposition of $\Sigma_c$ and consistency of the level smearing {#kbt-contour-smearing}

The kernel of `t_sigmaw` is not only a matter of "occupation numbers". In the contour decomposition that splits $\Sigma_c = G(W-v)$ into the imaginary-axis integral and the real-axis pole term,
**unless the kernel smearing the intermediate levels enters both terms in the same form, the steps do not cancel for levels near $\omega$,
and an error of the order of $|M|^2 W_c(0)$ (of the order of eV) appears**. When changing the kernel or trying another kernel, always come back here.

### Formula for sharp levels

For one intermediate level $\varepsilon'$, with $w_e \equiv (\omega-\varepsilon')/2$ (Hartree), and omitting the matrix elements $|M|^2$ of $W_c$
(the skeleton of PRB 76, 165106 Eq. 57–58):

$$
\Sigma_c^{(\varepsilon')}(\omega)
= \underbrace{-\frac{1}{\pi}\int_0^\infty d\omega'\,\frac{w_e\,W_c(i\omega')}{w_e^2+\omega'^2}}_{I(w_e)\ \text{imaginary-axis integral}}
\;+\;\underbrace{\theta\bigl[(\varepsilon'-\omega)(\varepsilon'-E_F)<0\bigr]\; s\,W_c(|w_e|)}_{P(w_e)\ \text{real-axis pole term}}
\tag{12}
$$

($s=\pm1$ for $\omega \gtrless E_F$, $W_c$ in $P$ is the real-axis value, and $\theta$ means "the level lies between $\omega$ and $E_F$").
As $w_e \to 0$, $I$ has a **step** $I \to -\tfrac12\,\mathrm{sign}(w_e)\,W_c(0)$, which exactly cancels the step of $P$ jumping $0 \to W_c(0)$ at the edge of the window of the pole term
($\varepsilon'=\omega$), so the sum is a continuous function of $\omega$.
The code splits $W_c(i\omega')$ into $W_c(0)\,e^{-(u_a\omega')^2}$ and the rest, and writes the former, which contains the step, analytically as

$$
I(w_e) = -\frac1\pi\int_0^\infty d\omega'\,\frac{w_e\,[W_c(i\omega')-W_c(0)e^{-(u_a\omega')^2}]}{w_e^2+\omega'^2}
\;-\;\frac12\,\mathrm{sign}(w_e)\,e^{a^2}\mathrm{erfc}(a)\,W_c(0),\qquad a=u_a|w_e|
\tag{13}
$$

(`wintz_npm`, the branch for core levels in m_sxcf_sc). The rest is smooth, so `niw`-point Gauss–Legendre suffices.

### The case of the Gaussian kernel (code up to 2026-09-19, `esmr`)

For valence levels, instead of (13), the following was used, which removes a piece of a two-dimensional Gaussian from the integrand:

$$
\frac{w_e}{w_e^2+\omega'^2}\;\to\;\frac{w_e\bigl(1-e^{-(w_e^2+\omega'^2)/2\sigma^2}\bigr)}{w_e^2+\omega'^2},
\qquad \sigma = \texttt{esmr}/2\ [\mathrm{Ha}]
\tag{14}
$$

The analytic integral of the removed piece is $+\tfrac12\mathrm{sign}(w_e)\,e^{a^2}\mathrm{erfc}\bigl(\sqrt{a^2+w_e^2/2\sigma^2}\bigr)W_c(0)$, and
together with the step term of (13) in the limit $u_a \to 0$ it gives

$$
-\tfrac12\,\mathrm{sign}(w_e) + \tfrac12\,\mathrm{sign}(w_e)\,\mathrm{erfc}\!\Bigl(\frac{|w_e|}{\sqrt2\sigma}\Bigr)
= -\tfrac12\,\mathrm{erf}\!\Bigl(\frac{w_e}{\sqrt2\sigma}\Bigr)
\tag{15}
$$

— that is, **the step on the imaginary-axis side smeared by a Gaussian (standard deviation $\sigma$)**. The window weight (10) on the pole-term side was also the Gaussian cumulative
$\Phi_G$ with the same $\sigma$, so the two sides were consistent and formed "$\Sigma_c$ with levels smeared by a Gaussian".
In other words, `sig` of PRB 76, 165106 Eq. 57 ($\sigma$ of Eq. (14)) was not a numerical regularization but **the level smearing itself**.

### What happens if the kernel is changed only in the pole term

If the weight of the pole term is changed to the Fermi–Dirac $\Phi_{FD}$ while the imaginary-axis side is left as (14), then for each level
$[\Phi_{FD}-\Phi_G](\varepsilon'-\omega)\times|M|^2 W_c(0)$ remains uncancelled. It matters only for levels with $|\varepsilon'-\omega|$ of
a few $k_BT$, but since the head $W_c(0)$ at $q\to0$ multiplies it, its size can be of the order of eV.
Measured (NiO 2³, O 2s pair at the L point, 70 meV difference in LDA, `t_sigmaw = 262 K`):

**Table 7**. $\Sigma_c$ when the kernel is changed on one side only (NiO 2³, pair of O 2s at the L point)

| SEc [eV] (state 1 / 2) | Gaussian on both sides | FD only in the pole term | FD on both sides (present) wcsmear=false | FD on both sides, wcsmear=true 30 / 262 / 473 K |
|---|---|---|---|---|
| | 10.523 / 10.541 | 9.909 / 10.992 | **10.523 / 10.541** (QSGW eigenvalues also agree) | 10.495 / 10.503 / 10.510 ; 10.524 / 10.528 / 10.533 |

With FD only in the pole term the error was ±0.6 eV (anticorrelated), 0.47 eV in the QSGW eigenvalues, and moreover strongly dependent on the width
(10.54/10.48/10.31/9.91/9.72 eV at 30/100/150/262/473 K). When both sides are made consistent, the width dependence almost disappears.

### Present formula — averaging with the FD kernel on both sides

$\Sigma_c$ with levels smeared by the kernel $g(x) = -\partial f/\partial x$ (width $k_BT$ = `t_sigmaw`) is averaged

$$
\bar\Sigma_c(\omega) = \int dx\,g(x)\,\Sigma_c^{(\varepsilon'+x)}(\omega)
= \int dx\,g(x)\,I\!\Bigl(\tfrac{\omega-\varepsilon'-x}{2}\Bigr)
+ \int dx\,g(x)\,P\!\Bigl(\tfrac{\omega-\varepsilon'-x}{2}\Bigr)
\tag{16}
$$

**in both terms**. The second term is `wcsmear` (window weight $\Phi_{FD}$, $W_c$ integrated with the kernel, §5).
The first term uses the sharp formula (13), and only the discontinuous $-\tfrac12\mathrm{sign}(w_e)$ is replaced analytically by

$$
\int dx\,g(x)\Bigl[-\tfrac12\,\mathrm{sign}(\omega-\varepsilon'-x)\Bigr]
= -\tfrac12\bigl[2\,\Phi_{FD}(\omega-\varepsilon')-1\bigr],\qquad
\Phi_{FD}(y) = \frac{1}{1+e^{-y/k_BT}}
\tag{17}
$$

while the continuous rest ($-\tfrac12\mathrm{sign}(w_e)\,(e^{a^2}\mathrm{erfc}(a)-1)$ and the numerical part) is averaged
with the midpoint rule at $x_j = k_BT\ln\frac{u_j}{1-u_j}$, $u_j = (j-\tfrac12)/40$
($u$ is uniform because the FD cumulative distribution is the measure). For levels with $|\omega-\varepsilon'| > 30\,k_BT$ the integrand is smooth, so
the sharp formula is used at one point. Core levels ($k_BT = 0$) stay with (13). `sig` is gone.
Implementation: the `CorrelationSelfEnergyImagAxis` block in `m_sxcf_sc.f90` (QSGW), and `wintzsg_npm` in `wintzsg.f90` (one-shot GW).

Notes:
- With `wcsmear = false`, the second term becomes "only the weight averaged with the kernel, $W_c$ at one point at the mean energy", so it is not the exact
  average of (16) (the difference is $O(k_BT^2\,W_c'')$). The cancellation of the steps is determined by the weights alone, so it is not broken.
- The separation $W_c(i\omega') \simeq W_c(0)e^{-(u_a\omega')^2}$ and the accuracy of the `niw` points are as before.
- The cost of the quadrature is $40\times$`niw` scalar operations per (itp, it) pair (on the CPU side, before the transfer to the GPU). A few seconds per rank even for 9³ LiTi₂O₄.
- **Rule when changing the kernel**: the three places, the occupation weight (`wfacx`), the window weight of the pole term (`wfacx2` / `pole_weights`), and the imaginary-axis step (17),
  must refer to the same $\Phi$. To go back to Gaussian, revert all three places as in (14)–(15).

The course of events and the numbers are in the development record (MD/research_kotani_log.md in the repository, not on the site), 2026-09-20 12:17 / 13:10.

---

## Choosing the width and the temperature

- **Insulators and semiconductors**: if the gap is sufficiently larger than $k_BT$, the width does not affect the results. The template value of 300 K is fine
  (Si: the difference of the gap between `t_tetrakbt` = 0 and 300 K is below 0.001 meV, and −1 meV at 1000 K. Table 9 in §8.1).
- **`t_sigmaw` for metals**: in principle the limit $n_1n_2n_3\to\infty$, `t_sigmaw` → 0 is what is wanted, but with a finite k mesh
  a width at least comparable to the level spacing that the mesh has near the Fermi surface is needed. A guide is

$$
k_BT \gtrsim \frac{v_F \cdot (\text{size of the BZ})}{N}
\tag{18}
$$

  ($N$ is the number of mesh points in one direction). 1000 K ($k_BT = 0.086$ eV) for LiTi₂O₄ $6^3$ and $9^3$.
- **`t_tetrakbt`**: to see the physics of the electronic temperature, make it positive and set `t_sigmaw` to the same temperature. To only smear the poles of $W$ numerically, make it negative
  (Gaussian). In this case the Fermi level and the occupations stay at T=0.
- Check that the quantity you want is converged with respect to the width, and if possible that two k meshes agree with the same width.
- The electronic temperature enters only $\chi_0$ and $\Sigma$. The density made by `lmf` keeps the T=0 occupations (tetrahedron method), and the lattice temperature (phonons) does not enter either.

---

## Samples — `Samples/kBT/`

They are in `Samples/kBT/` of ecalj (Table 8). All numbers and figures, and how to run them, are in the README of each directory.

**Table 8**. Samples

| Directory | Contents | Time |
|---|---|---|
| `scanT/` | Temperature scan. `scan_T.sh` runs QSGW on the same input changing only `t_tetrakbt` and `t_sigmaw`, and `collect_T.py` makes tables and figures. Si and GaAs (§8.1), and bcc Fe and fcc Cu (§8.2) | 2 to 10 min per setting |
| `LiTi2O4/` | MLO-QSGW of the metallic spinel LiTi₂O₄ (§8.3) | 21 min per iteration for $9^3$ (RTX 5090 × 2) |
| `scanT/Si`, `Samples/TestInstall/fe_kbt` | Tests of `testecalj` (Si and Fe at 3000 K, 1 iteration) | 20 to 40 s |

In all samples, the electronic temperature enters only $\chi_0$ and $\Sigma$ (§7). `scanT` has a coarse k mesh and few iterations, so
its numbers are for seeing the change with temperature, not converged values.

### Si and GaAs — gap and electronic temperature

`nkabc` = `n1n2n3` = $4^3$, 5 QSGW iterations from LDA. `t_tetrakbt = t_sigmaw = T`.

**Table 9**. Gap of Si (minimum on the k mesh of the SCF; after 5 QSGW iterations)

| `t_tetrakbt` (K) | `t_sigmaw` (K) | Gap (eV) | Change from $T=0$ (meV) |
|---|---|---|---|
| 0 | 300 | 1.223 | 0 |
| 300 | 300 | 1.223 | 0 |
| 1000 | 1000 | 1.222 | −1 |
| 2000 | 2000 | 1.188 | −35 |
| 3000 | 3000 | 1.066 | −157 |
| 5000 | 5000 | 0.913 | −310 |
| 3000 | 300 | 1.090 | −133 |
| 0 | 3000 | 1.176 | −46 |
| −3000 | 3000 | 1.173 | −50 |

![Gap of Si and electronic temperature](../../manual/kBT/si_T_gap.png)

**Figure 1**. Left: gap of Si versus electronic temperature (`t_tetrakbt = t_sigmaw = T`). Right: gap at each iteration

- At 300 K it is the same as at $T=0$ (difference below 0.001 meV). The shrinking becomes visible from 2000 K, and the gap shrinks by 0.16 eV at 3000 K.
- The shape of the valence band hardly changes; the conduction band moves down as a whole relative to the valence band.
- The shrinking comes mainly from the $\chi_0$ side (of the −157 meV at 3000 K, −133 meV when only $\chi_0$ is at 3000 K, and −46 meV with only $\Sigma$).
  Carriers thermally excited across the gap join the screening, and $W$ becomes smaller.
- With `t_tetrakbt = -3000` (Gaussian), the result is almost the same as the case with only $\Sigma$ (−50 meV). The Gaussian does not change the occupations, so the screening by carriers does not enter (§3.4).
- GaAs ($4^3$, 5 QSGW iterations, direct gap 1.415 eV) shows the same trend: no change at 300 K, −4 meV at 1000 K, −20 meV at 2000 K, −95 meV at 3000 K,
  and −267 meV at 5000 K. At 3000 K, −83 meV with $\chi_0$ only and −39 meV with $\Sigma$ only (numbers and figures in Table 3 of `Samples/kBT/scanT/README.md`).

### bcc Fe and fcc Cu — bands and magnetic moment of metals

`nkabc` = $7^3$, `n1n2n3` = $5^3$, spin polarized, 3 QSGW iterations from LDA.

**Table 10**. bcc Fe (after 3 QSGW iterations). The exchange splitting is the difference between spin 2 and spin 1 of the d bands at Γ (bands 5 and 6)

| `t_tetrakbt` (K) | `t_sigmaw` (K) | Magnetic moment ($\mu_B$) | `EFERMI_kbt` − `EFERMI` (eV) | Exchange splitting (eV) |
|---|---|---|---|---|
| 0 | 300 | 2.429 | — | 2.619 |
| 300 | 300 | 2.426 | +0.006 | 2.604 |
| 1000 | 1000 | 2.443 | +0.045 | 2.688 |
| 3000 | 3000 | 2.420 | +0.272 | 2.492 |
| 5000 | 5000 | 2.252 | +0.279 | 1.842 |
| 3000 | 300 | 2.441 | +0.258 | 2.676 |
| 0 | 3000 | 2.409 | — | 2.450 |
| −1000 | 1000 | 2.443 | — | 2.685 |

- At 300 K it is almost the same as at $T=0$ (the difference of the bands at the symmetry points is at most 14 meV).
- Up to 1000 K the bands move by at most 63 meV, not monotonically with temperature. With the $5^3$ mesh, the bands move by this much depending on which states near $E_F$
  become fractionally occupied (item 1 of §12).
- From 3000 K the exchange splitting and the moment decrease. In metals the $\Sigma$ side matters: the change of the exchange splitting at 3000 K, −0.13 eV, is close to the case with only $\Sigma$ at 3000 K (−0.17 eV),
  while with only $\chi_0$ it is +0.06 eV. This is the opposite of Si.
- `t_tetrakbt = -1000` (Gaussian) and `t_tetrakbt = 1000` (finite temperature) give almost the same results at this temperature (the difference of the bands is at most 7 meV).
- For nonmagnetic fcc Cu (`nkabc` = $8^3$, `n1n2n3` = $6^3$, 3 QSGW iterations), the bands move by 11 meV or less up to 1000 K. From 3000 K the d bands approach $E_F$
  (the d band at Γ moves by about +61 meV at 3000 K and +130 meV at 5000 K). +46 meV with $\chi_0$ only and +24 meV with $\Sigma$ only, so both matter
  (Table 4 of `Samples/kBT/scanT/README.md`).

### LiTi₂O₄ — MLO-QSGW

The metallic spinel LiTi₂O₄ (14 atoms) was run from LDA up to 40 iterations with the "Setting for LiTi₂O₄" at the top, by QSGW with $\Sigma$ interpolated by MLOs ([mlo_gwsc](./mlo_gwsc)).
The summary is in `Samples/kBT/LiTi2O4/README.md` of ecalj.

- The difference of the MLO bands due to the precision of the GPU matrix products (tf32, fp32, fp64) is a few meV or less.
- For $6^3$ the bands are smooth up to 10 iterations; from 20 iterations a shoulder appears in the occupied bands, which remains even at 40 iterations.
- For $9^3$ the wiggles between mesh points decrease only slowly, but become fairly smooth by 40 iterations. The interpolating curve through the mesh points looks slightly wavy.

---

## Where it is implemented

**Table 3**. Implementation of finite temperature and smearing (`SRC/subroutines/`)

| File | Role |
|---|---|
| `m_GWinput.f90` | Reads `t_tetrakbt` (required, Table 1), `t_sigmaw`, and `wcsmear` of `[gw]`. Converts `t_tetrakbt < 0` to $\sigma$ (internal `SmearX0`, Ha) (Eq. (7)) |
| `tetwt5.f90` | Tetrahedron weights of $\chi_0$. `lindtet6` (Eq. (3)) and `lindtet6_kbt` (Eq. (5)) |
| `fpiint.f90` | `gausq_fd`: the 20 points of the thermal kernel (Eqs. (5)(6)) |
| `main_heftet.f90` | `fermi_kbt_tetra`: `EFERMI_kbt` (Eq. (6)) |
| `m_tetwt.f90` | Uses `lindtet6_kbt` and `EFERMI_kbt` if `t_tetrakbt > 0` |
| `dpsion5.f90`, `hilbertmat.f90` | Hilbert transform (Eq. (2)) and Gaussian (`gaussianfilterhis`, `smearx0_apply`, Eq. (8)) |
| `genallcf_mod.f90` | `sigmakbt_setup`: Fermi level of $\Sigma$ (`EFERMI_kbt` if `t_tetrakbt > 0`; stops if absent) |
| `wfacx.f90` | Fermi–Dirac kernel on the $\Sigma$ side (`fd_cdf`, `fd_iav`. Width `t_sigmaw`, Eq. (10)) |
| `m_sxcf_sc.f90`, `sxcf_fal2.f90` | Real-axis pole term of $\Sigma_c$ and `wcsmear` (§5, §6) |

---

## Verification and limitations of the implementation

What was checked by rereading the implementation, and the limitations that remain. The contour identity on the $\Sigma$ side (imaginary axis + pole term) is in
[§6](#kbt-contour-smearing).

### The identity is exact (checked)

Method B′ of `lindtet6_kbt` (`tetwt5.f90`) uses Eq. (4) (§3.3) **inside the k integral**. Since the $E$ integral and the $k$ integral
commute, it is exact. The numerator obtained is the Adler–Wiser $f_a-f_b$, not $f_a(1-f_b)$
— this is the correct one for a response function.

Checked individually:

- **Gate** (at the top of `lindtet6_kbt`) — all the cases in which a/b are flat outside the window were traced,
  and all agree with the value at $T=0$
- **node zero-skip** (the loop over $E$ in `lindtet6_kbt`) — exactly the condition under which `lindtet6` is identically 0.
  `knorm` is normalized over all nodes, so skipping is consistent
- (4) holds even if `efermia` ≠ `efermib` (it becomes $f$ at the respective Fermi levels)
- **`EFERMI_kbt`** (`fermi_kbt_tetra`, `main_heftet.f90`, Eq. (6))
  — the NOS convolved with the thermal kernel, then bisection. As $T\to0$ it returns exactly to the tetrahedron $E_F$.
  The χ0 side also reads it (`m_tetwt`)

### Limitation — quadrature of the $E$ integral

The $E$ integral of Eqs. (5)(6) uses 20 points: 5-point Gauss-Legendre on each of the 4 intervals $[-6,-1],[-1,0],[0,1],[1,6]$ into which $t\in[-6,6]$ is divided
(`gausq_fd`). The node closest to $E_F$ is at $|t|=0.047$ ($|E-E_F|=0.094\,k_BT$). With $\Delta$ the width of the band inside a tetrahedron,
Table 11 gives the measured worst-case error of the occupation factor.

**Table 11**. Quadrature of the $E$ integral and the worst-case error of the occupation factor

| Quadrature | Nodes | $\Delta=0.4k_BT$ | $1k_BT$ | $2k_BT$ | $4k_BT$ | $8k_BT$ |
|---|---|---|---|---|---|---|
| **Present implementation**: 4 intervals $[-6,-1,0,1,6]$ × GL5 | 20 | 2.9e-2 | 1.6e-2 | 4.8e-3 | 2.6e-3 | 1.3e-3 |
| 6 intervals × GL4 | 24 | 2.9e-2 | 1.1e-2 | 2.8e-3 | 1.7e-3 | 9.6e-4 |
| 2 intervals $[-6,0,6]$ × GL10 | 20 | 5.5e-2 | 1.9e-2 | 1.1e-2 | 3.8e-3 | 1.9e-3 |
| 1 interval $[-6,6]$ × GL20 (code up to 2026-09-26) | 20 | 1.7e-1 | 9.7e-2 | 2.3e-2 | 1.3e-2 | 6.8e-3 |

For bands with large dispersion the error is small. When a flat band, a band edge, or a van Hove singularity is at $E_F$, the occupation factor is off by 0.01 to 0.03.
Increasing the number of nodes reduces the error only as $1/N$; how the intervals are cut matters more. The proper way is to cut the intervals at the kinks of the integrand (the 8 energies $e_a$, $e_b$ at the 4 corners of the tetrahedron),
but this is not implemented. The error of `EFERMI_kbt` for Fe 3000 K $5^3$ (relative to a 600-point quadrature) is 1.2 meV.

### Window of states on the $\Sigma$ side

The candidate window of intermediate states and the margin of the real-axis $\omega$ mesh are determined by the tail of the kernel $E_F\pm15\,k_BT$ ($k_BT$ is `t_sigmaw`; `sig_window`).
The occupation cut off is $f(15)=3\times10^{-7}$, and the window widens when `t_sigmaw` is raised. The center of the window and the center of the weights are the same Fermi level
(`EFERMI_kbt` if `t_tetrakbt > 0`). In the code up to 2026-09-19 the window was determined by `esmr`, and at high temperatures occupations were cut off
(the window was $E_F\pm10\,$`esmr`; for example, at 2000 K with `esmr` = 0.002 Ry, 17 % of the occupation is cut off. The development record is MD/kBT/history.md §4 in the repository, not on the site).

### Scope — what $\Sigma$ needs is the "sum", not the "difference"

Making $\chi_0$ finite-temperature alone does not make $\Sigma$ consistent. The reason can be stated in one line.

With the particle‑hole pair of the intermediate state written $(2,3)$ and $\omega'=\varepsilon_2-\varepsilon_3>0$,

**Table 12**. Weights of a particle‑hole pair

| | Weight |
|---|---|
| **Creating** a pair | $f_3(1-f_2)$ |
| **Destroying** an existing pair | $f_2(1-f_3)$ |

What $\mathrm{Im}\,\chi_0$ contains is the **difference** $f_3(1-f_2)-f_2(1-f_3)=f_3-f_2$, but
the weight of the intermediate states of $\Sigma$ is the **sum** (both have the same denominator).
The present implementation takes the difference from $\mathrm{Im}W$ and multiplies by the outer $f_j$, so
**the part $f_2(1-f_3)$ is dropped.** From detailed balance $f_2(1-f_3)/f_3(1-f_2)=e^{-\beta\omega'}$,

$$
f_2(1-f_3)=\bigl(f_3-f_2\bigr)\,n_B(\omega')
\tag{19}
$$

so this is what is usually called the Bose factor. **However, its origin is the Fermi
occupation numbers**, not a bosonic heat bath put in from outside. For the derivation see
[§11](#kbt-real-time). The physics is simple: at $T>0$
particle‑hole pairs are already thermally excited, and electrons can absorb them; that is all.

**Size.** $f_2(1-f_3)$ is exponentially small for $\omega'\gg k_BT$, so only
the particle‑hole continuum with $\omega'\lesssim k_BT$ matters. For a metal with
$\mathrm{Im}W_c\simeq c\,\omega'$,

$$
\int_0^\infty\!\! d\omega'\,\mathrm{Im}W_c(\omega')\,n_B(\omega')
= c\,(k_BT)^2\frac{\pi^2}{6}
\tag{20}
$$

At 2000 K, $(k_BT)^2=0.029$ eV², which for $\omega-\varepsilon_j\sim1$ eV gives a few meV to a few tens of meV.

**This is not a defect but a definition.** The ecalj scheme is "an effective one-body Hamiltonian made with the occupation numbers at temperature $T$",
the same standpoint as Mermin's finite-temperature DFT, which closes with the Fermi occupation numbers alone.
QSGW makes a **static Hermitian** one-body Hamiltonian, not a dynamical $\Sigma$,
and partly because what is dropped is mainly on the detailed-balance (= lifetime) side,
the effect is much lighter than when computing the spectral function of one-shot GW.
It is simply that (19) is needed if one wants the $\Sigma$ of equilibrium many-body perturbation theory.

### Core exchange (CoreEx)

The exchange with the core (`hsfp0_sc --job=3`) places the Fermi level below the bottom of the valence states even at finite temperature, so that only core states enter the sum.
The core levels are much deeper than $k_BT$, so the result is the T=0 value regardless of the temperature.

**Caution on finite-temperature results of the code from 2026-06-26 to 09-27**: the code in this period had two errors in CoreEx, and
in systems where the number of core states $n_\mathrm{ctot}$ is smaller than the number of occupied bands (such as compounds of light elements whose only core is 1s), part of the valence states
entered the exchange with the core. This can be judged by comparing `CoreEx mode: ef nspin nctot=` in `lsxC` with `valn=` (the number of valence electrons)
(note that systems with semicore states in the valence have fewer core states). Fe ($n_\mathrm{ctot}=9$) and LiTi₂O₄ ($n_\mathrm{ctot}=46$) are not affected.
The errors were that the Fermi level of CoreEx was overwritten by `EFERMI_kbt`, so that occupied valence states entered the sum, and an index shift in `sxs_ekc` (both fixed on 2026-09-27). The development record is MD/kBT/history.md §5 in the repository, not on the site.

---

## Appendix — formulation without imaginary time {#kbt-real-time}

Where $n_B$ comes from becomes clear without going through the Matsubara formalism. Rather,
that way it is easier to see that **everything comes from the Fermi occupation numbers**.
Textbooks often describe the closed time path and equilibrium separately,
so here it is written through in one piece.

Below, $\hbar=1$, and energies are measured from $\mu$.
$f(\varepsilon)=1/(e^{\beta\varepsilon}+1)$, $n_B(\omega)=1/(e^{\beta\omega}-1)$.

### Why the $T=0$ argument cannot be used as it is

Real-time perturbation theory at $T=0$ relies on Gell-Mann–Low. When the interaction is switched on adiabatically,
$U(\infty,-\infty)|\Phi_0\rangle = e^{i\theta}|\Phi_0\rangle$, that is,
**the ground state returns to itself up to a phase**, so the $+\infty$ side can be
exchanged with the $-\infty$ side, and

$$
\langle\Psi_0|T\{\cdots\}|\Psi_0\rangle
=\frac{\langle\Phi_0|T\{S\cdots\}|\Phi_0\rangle}{\langle\Phi_0|S|\Phi_0\rangle}
\tag{21}
$$

it can be written as a one-way time-ordered product.

At $T>0$ this cannot be used. $\rho_0=e^{-\beta(H_0-\mu N)}/Z_0$ is not an eigenstate, and
even under adiabatic evolution **the weights remain those of the noninteracting system** ($e^{-\beta E_n^{(0)}}$, not
$e^{-\beta E_n}$), so the $+\infty$ side cannot be identified with the $-\infty$ side.
The essence of the problem is that the counterpart of the denominator of (21) cannot be written.

### Closed time path — going and coming back

So **the state at $+\infty$ is not used at all**. If one takes the path $C$ that goes
$-\infty\to+\infty\to-\infty$ and back, then as long as no operators are inserted

$$
\mathrm{Tr}\bigl[\rho_0\,U(-\infty,+\infty)\,U(+\infty,-\infty)\bigr]
=\mathrm{Tr}\,\rho_0 = 1
\tag{22}
$$

holds **identically**. No denominator is needed. This is all there is to the closed time path (Keldysh).

With the path-ordered product $T_C$ (the ordering in which the return branch comes after the forward branch),

$$
G(1,2)=-i\,\mathrm{Tr}\Bigl[\rho_0\,T_C\Bigl\{
e^{-i\int_C dt\,H_1(t)}\;\psi(1)\psi^\dagger(2)\Bigr\}\Bigr]
\tag{23}
$$

The fields are in the interaction picture of $H_0$. It is in real time, not $\tau$, and the $T$ product has merely become the $T_C$ product.
Four components arise depending on which branch $1,2$ are placed on, but only 2 are independent, and
below

$$
G^<(1,2)=+i\langle\psi^\dagger(2)\psi(1)\rangle,\qquad
G^>(1,2)=-i\langle\psi(1)\psi^\dagger(2)\rangle
\tag{24}
$$

are used.

### Wick's theorem — the only input is $f$

$\rho_0$ is Gaussian with respect to $H_0$, so Wick's theorem holds as it is,
and the contractions become free path propagators. For a level of energy $\varepsilon_j$,

$$
g_j^<(\omega)=2\pi i\,f_j\,\delta(\omega-\varepsilon_j),\qquad
g_j^>(\omega)=-2\pi i\,(1-f_j)\,\delta(\omega-\varepsilon_j)
\tag{25}
$$

**This is the only place where temperature enters, and what enters is only the Fermi distribution.**
There is no Bose distribution anywhere at this stage.

### Writing GW on the path

With path arguments, the form is the same as at $T=0$:

$$
\chi_0(1,2)=-i\,G(1,2)G(2,1),\qquad
W=v+v\chi_0 W,\qquad
\Sigma(1,2)=i\,G(1,2)W(2,1)
\tag{26}
$$

To reduce to real-time components, Langreth's rules are used. For a product with the same arguments
$C(1,2)=A(1,2)B(2,1)$,

$$
C^{\gtrless}(1,2)=A^{\gtrless}(1,2)\,B^{\lessgtr}(2,1)
\tag{27}
$$

The key point is that **$\gtrless$ is flipped**.

### The $\gtrless$ components of $\chi_0$ — this is where $n_B$ appears

From (26)(27) and (25), using time-translation symmetry,

$$
\chi_0^>(\omega)=-2\pi i\sum_{23} f_3(1-f_2)\,\delta(\omega-\varepsilon_2+\varepsilon_3)
$$
$$
\chi_0^<(\omega)=-2\pi i\sum_{23} f_2(1-f_3)\,\delta(\omega-\varepsilon_2+\varepsilon_3)
\tag{28}
$$

$\chi_0^>$ is the process that **creates a pair**, and $\chi_0^<$ the process that **destroys an existing pair**.
For $\omega>0$, $\varepsilon_2>\varepsilon_3$.

The spectral function is the difference,

$$
\chi_0^>-\chi_0^< \;\propto\; f_3(1-f_2)-f_2(1-f_3)=f_3-f_2
\;=\;2i\,\mathrm{Im}\chi_0^R
\tag{29}
$$

This is the quantity that the $T=0$ code computes. On the other hand, the ratio in (28), simply by substituting
$f_j=1/(e^{\beta\varepsilon_j}+1)$, becomes

$$
\frac{\chi_0^<(\omega)}{\chi_0^>(\omega)}
=\frac{f_2(1-f_3)}{f_3(1-f_2)}=e^{-\beta\omega}
\tag{30}
$$

and using $1/(1-e^{-\beta\omega})=1+n_B(\omega)$,

$$
\chi_0^>=\bigl(1+n_B\bigr)\bigl(\chi_0^>-\chi_0^<\bigr),\qquad
\chi_0^<=n_B\,\bigl(\chi_0^>-\chi_0^<\bigr)
\tag{31}
$$

**$n_B$ appears here for the first time. It is not an assumption but a consequence of the Fermi distribution in (25)**.
The Bose distribution appears even though no bosonic heat bath was put in from outside, because
$\chi_0$ is a bosonic correlation function ((30) is exactly the bosonic KMS condition).

Dressing with RPA does not change the ratio. Since $W^{\gtrless}=\epsilon^{-1,R}\,v\chi_0^{\gtrless}v\,
\epsilon^{-1,A}$ with $\epsilon^{-1,A}=(\epsilon^{-1,R})^\dagger$,
(30) is inherited as it is by $W^{\gtrless}$.

### The $\gtrless$ components of $\Sigma$

From (26)(27), $\Sigma^{\gtrless}(t)=i\,G^{\gtrless}(t)\,W^{\lessgtr}(-t)$,
and in frequency

$$
\Sigma^{\gtrless}(\omega)=i\sum_j\int\!\frac{d\nu}{2\pi}\,
g_j^{\gtrless}(\nu)\,W^{\lessgtr}(\nu-\omega)
\tag{32}
$$

Inserting (25), and with $B(\Omega)\equiv W^>(\Omega)-W^<(\Omega)=2i\,\mathrm{Im}W^R(\Omega)$,
$W^>=(1+n_B)B$, $W^<=n_B B$,
from $2i\,\mathrm{Im}\Sigma^R=\Sigma^>-\Sigma^<$

$$
\mathrm{Im}\,\Sigma^R(\omega)=\sum_j
\mathrm{Im}W^R(\varepsilon_j-\omega)\;
\bigl[\,n_B(\varepsilon_j-\omega)+f_j\,\bigr]
\tag{33}
$$

Splitting by the sign of $\varepsilon_j-\omega$ and using $n_B(-\Omega)=-(1+n_B(\Omega))$,
one obtains the familiar emission factor $1-f_j+n_B$ and absorption factor $f_j+n_B$. The real part is
determined by Kramers–Kronig.

Check of $T\to0$: $n_B(\Omega>0)\to0$, $n_B(\Omega<0)\to-1$, $f_j\to\theta(-\varepsilon_j)$.
For a quasiparticle with $\omega>0$: if $j$ is occupied, $[-1+1]=0$ (Pauli blocking);
if $j$ is empty and $0<\varepsilon_j<\omega$, it contributes $[-1+0]=-1$ (decay); and if $\varepsilon_j>\omega$, $[0+0]=0$ (no final state). This returns to the conventional $T=0$ GW.

### What the present implementation keeps and what it drops

What ecalj computes is $\mathrm{Im}W^R$, that is, the **difference** in (29).
When building $\Sigma$, only the outer $f_j$ is made finite-temperature, and the $n_B$ of (33) is not included.
That is,

$$
\Sigma^{\text{ecalj}}:\quad n_B\to 0,\qquad f_j\to f_j(T)
$$

This breaks the detailed balance of $\Sigma$

$$
\Sigma^<(\omega)=-e^{-\beta\omega}\,\Sigma^>(\omega)
\tag{34}
$$

(from (32) and (30), (34) holds only when **both** are kept).
Therefore the $\Sigma$ obtained is, strictly speaking, not that of an equilibrium state at any temperature $T$.

On the other hand, this is a light matter in QSGW. What QSGW makes is a **static Hermitian** one-body
Hamiltonian, and $\mathrm{Im}\Sigma$ is discarded from the start.
That (34) is broken is mainly a matter of the lifetime side, and the effect on the
$\mathrm{Re}\,\Sigma(\varepsilon_j)$ used by QSGW stays at the $O((k_BT)^2)$ of (20).

### Note on initial correlations

In 8.1–8.2 it was said that the interaction is switched on adiabatically from $\rho_0$ (thermal equilibrium without interaction), but
this **drops the initial correlations**. Strictly, one uses the Kadanoff–Baym path, in which the vertical imaginary-time branch
$[t_0,\,t_0-i\beta]$ is added to the path, and the initial correlations are
carried there (Danielewicz). In the equilibrium and adiabatic case the vertical branch separates from the real-time part,
and the result agrees with the Matsubara formalism.

Here $\tau$ shows up **to specify the initial condition**;
the perturbation expansion itself remains in real time. It has nothing to do with whether $n_B$ appears —
as seen in 8.5, that comes only from the Fermi distribution in (25).

---

## Remaining issues

1. **Verification of the choice of `t_sigmaw` in metals**. A few states within $\pm k_BT$ of $E_F$ become fractionally occupied, and their SEx and SEc move by the size of eV and
   cancel each other. On a coarse mesh this is sensitive to "which states cross $E_F$", so the convergence with respect to the mesh ($5^3\to9^3$) needs to be pinned down numerically
   (temperature dependence of Fe in §8.2).
2. **`EFERMI_kbt` for insulators**. When the solution is an interval in the gap, the middle of the interval is returned. Physically it should return the intrinsic $E_F$ determined by the asymmetry of the density of states,
   and the result changes with the relation between the window of `fermi_kbt_tetra` ($E_F\pm60\,k_BT$) and the reach of the thermal kernel ($\pm12\,k_BT$) (GaAs at 300 K: 0.15345 Ry; the middle at T=0 is 0.15221 Ry).
3. **Bose term of $O((k_BT)^2)$** (§10.4, §11.7). $f_2(1-f_3)=(f_3-f_2)\,n_B$ is dropped. Whether it matters at 3000 K has not been estimated.
4. **Cutting the $E$ integral at the kinks of the integrand** (§10.2).
5. **Off-diagonal elements of pairs with a plasmon pole between the two levels** (limitation in §5). Not removed by `wcsmear`.
6. **The density of `lmf` keeps the T=0 occupations** (§7). The contribution of carriers thermally excited at high temperature to the density is not included.
7. **LiTi₂O₄ $9^3$ at 3000 K**: there is no calculation that follows the iterations.
8. **The Fermi level of $\Sigma$ when `t_tetrakbt` ≤ 0** is determined by Gaussian counting (Table 5 in §1), while the kernel smearing the levels is Fermi-Dirac.
   They are not made the same kernel.

---

## About gwkbt-dev

Everything on this page is in `main`. There is no need to switch branches.
The separate `gwkbt-dev` is a half-finished attempt to build the finite-temperature $\Sigma$ with a different design (gwkbt Stage A/B), and is not in a usable state.
It is different from `t_sigmaw` on this page.

---

## Related

- [QSGW calculation (gwsc)](./gwsc)
- Samples `ecalj/Samples/kBT/` (§8)
- Development record (MD/kBT/history.md in the repository, not on the site): the course of events in 2026-06 to 09, explanations of the calculations of that time, errors fixed
- [The `[gw]` section of `ctrlg.<sname>.toml`](./lmf#file-structure-sections)
- For QP energies of one-shot GW, evaluation on arbitrary k lines, and diagnostics of $W$ and the real-axis integral, see
  the ecalj development record (MD/past_log.md §5 in the repository, not on the site) (the essentials of the 2026-06 guide `FiniteT_and_QPE_HOWTO.md`; this page develops its §1)
