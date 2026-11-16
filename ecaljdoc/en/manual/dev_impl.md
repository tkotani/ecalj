# Implementation notes (details of the GW implementation) {#dev-impl}

> English version of the page [manual/dev_impl](../../manual/dev_impl.md). The pages at the top of ecaljdoc are the reference; this one was made from manual/dev_impl.md of commit 4f74d9132 (source blob 26a10313902e).

::: info About this page
This page is made from the sections, still useful, of the old developer memo `ecaljdetails.tex` (T. Kotani, last version 2024-11; the original is `ecaljdetails/` of the ecalj commit `3284d0b`, `git show 3284d0b:ecaljdetails/ecaljdetails.tex`),
reviewed against the present code and M. Obata et al. (2025).
The equation numbers (I1)– were given when copied; the numbers of the removed equations are missing. `[key]` is a key of the references (`ecaljdetails/ecaljrefs_detail.bib` of the same commit).
For the theory, see [Band calculation method of ecalj](../theory/pmt), [QSGW](../theory/qsgw) and the references there. The mathematical notes are in the [Appendix](../theory/appendix).
:::


The main reference for the formulas of the GW part is [Obata2025]: M. Obata, T. Kotani, T. Oda, "Efficient implementation of the quasiparticle self-consistent $GW$ method on GPU", Comput. Mater. Sci. **260**, 114190 (2025), [doi:10.1016/j.commatsci.2025.114190](https://doi.org/10.1016/j.commatsci.2025.114190), section 2. The equations below that carry "[Obata2025] Eq. (n)" are taken from it.

## Representation of eigenfunctions

In the PMT method [kotani_fusion_2010],
the valence eigenfunctions for a given $H^0$ are represented
in the linear combinations of the Bloch-summed MTOs
$\chi^{ {\bf k} }_{ {\bf R}L{j} }({\bf r})$ and the APWs $\chi^{ {\bf k} }_{\bf G}({\bf r})$;
$$\begin{aligned}
{\Psi_{ {\bf k}n} }({\bf r}) = \sum_{ {\bf R}L{j} } z^{ {{\bf k} }n}_{ {\bf R}L{j} }
\chi^{ {\bf k} }_{ {\bf R}L{j} }({ {\bf r} })+ \sum_{ {\bf G} } z^{ {\bf k}n}_{ {\bf G} }
\chi^{ {\bf k} }_{\bf G}({\bf r}),
\end{aligned}
\tag{I1}
$$
where we use indexes of the wave vector ${\bf k}$, band index $n$, and
reciprocal lattice vector ${\bf G}$. The MTOs in the primitive cell are
specified by the index of MT site ${\bf R}$,
angular momentum $L=(l,m)$, and $j$ for radial functions.
As for core eigenfunctions, we calculate them under the condition
that they are restricted within MTs.
Then we take into accounts the contributions of the cores to the exchange part defined in Eq. (I15) in the following. But not to the correlation part.
All cores are treated as "core1" of Ref.[kotani07a] (the `core` rows of `[product_basis]` are all zero).

###  MTO part

Within MTs, the Bloch sum of the MTO,
$\chi^{ {\bf k} }_{ {\bf R}L{j} }({\bf r})$, is
expressed by a linear combination of atomic like orbitals
$A_{ {\bf R}{Lj} }({\bf r})
\equiv \{ 
\phi_{ {\bf R}{Lj} }(r),\dot{\phi}_{ {\bf R}{Lj} }(r),\phi^z_{ {\bf R}{Lj} }(r) \}
\times Y_L$. ($\phi^z$ means local orbital).
These radial functions are solutions of
the radial Schrödinger equations(or their energy derivatives)
within ${\bf R}$.
The MTO basis is specified by smHankel functions which
contains two parameters ($E=-|\kappa|^2,R_{\rm sm}$).

$A_{ {\bf R}{Lj} }({\bf r})$ makes orthonormalized basis for each MT ${\bf R}$.
Then the MTO including tail part can be written as
$$\begin{aligned}
\chi^{ {\bf k} }_{ {\bf R}{Lj} }({\bf r}) & =  \sum_{ {\bf R}{Lj} } C^{ {\bf k} }_{ {\bf R}{Lj} } A^{\bf k}_{ {\bf R}{Lj} }({\bf r})   \ \ \ {\rm \ if \ {\bf r} \in any \ MT}  \\
        & =    H^{\kappa,R_s,{\bf k} }_{ {\bf R}{L} }({\bf r}) \ \ \ {\rm otherwise},
\end{aligned}
\tag{I2}
$$
where we use the Bloch sums,
$$\begin{aligned}
A^{\bf k}_{ {\bf R}{Lj} }({\bf r}) & \equiv  \sum_{\bf T} A_{ {\bf R}{Lj} }({\bf r-{\bf R}-{\bf T} }) \exp(i {\bf kT}), \\
H^{ {\bf k}s}({\bf r})   & \equiv  \sum_{\bf T} H_{s}({\bf r-{\bf R}-{\bf T} })
 \exp(i {\bf kT}).
\end{aligned}
\tag{I3}
$$
Here the smoothe Hankel functions $H^{ {\bf k}s}({\bf r})$ are
the envelope functions of MTOs.

The MTOs have the symmetry of the angular momentum $L$.
This is also the casee for any atomic-like localized functions.
All such kinds of functions follow the same mapping formulas by space-group operations.
See Sec.([§](#spacegroup-rotation)).

###  APW part

The APW $\chi^{ {\bf k} }_{\bf G}({\bf r})$ are
given as a linear combination of atomic like orbitals
$A_{ {\bf R}{Lu} }({\bf r})\equiv 
\{ \phi_{ {\bf R}{l}u}(r) Y_L(\hat{\bf r}),
    \dot{\phi}_{ {\bf R}{l}u}(r) Y_L(\hat{\bf r}) \}$
within MTs, and just the usual plane waves
within the interstitial region. Here
$\phi_{ {\bf R}{l}u}(r)$ and $\dot{\phi}_{ {\bf R}{l}u}(r)$
denote two solutions of
the radial Schrödinger equations for each $l$, with the logarithmic derivative at the MT radius given.
With `[ham] readp = true` and `pnufix = true` (the defaults), the logarithmic derivatives are those of the free atom of `lmfa` (`atmpnu.*`) and are kept fixed
(with `pnufix = false`, `pnunew.f90` sets them at the center of gravity of the occupied PDOS; [lmf § readp](./lmf#readp-option)). $\dot{\phi}$ means energy derivatives
(or something similar).
$u$ is the composite index
to diffrenciate $\phi$ and $\dot{\phi}$.
${\bf R}$ is the index to specify MTs in the primitive cell.
The APW basis is specified by $s \equiv { {\bf R}{j}L}$, where
$L\equiv(l,m)$ is the angular momentum index,
and $j$ is the additional index (principle quantum number or so).
$A_{ {\bf R}{L}u}({\bf r})$ makes normalized-orthogonal basis in
each MT ${\bf R}$.
The APW can be written as
$$\begin{aligned}
\chi^{\bf k+G}({\bf r}) & =  \sum_{a u} C^{\bf k+G}_{a u} A^{\bf k}_{a u}({\bf r})   \ \ \ {\rm \ if \ {\bf r} \in any \ MT}  \\
        & =    \exp(i ({\bf k+G}){\bf r}) \ \ \ {\rm otherwise},
\end{aligned}
\tag{I4}
$$
where we use the Bloch sums,
$$\begin{aligned}
A^{\bf k}_{ {\bf R}u}({\bf r}) & \equiv  \sum_{\bf T} 
A_{ {\bf R}u}({\bf r-R-T}) \exp(i {\bf kT}),
\end{aligned}
\tag{I5}
$$

The number of ${\bf G}$ is limited by the condition
$|{\bf k+G}|\lt  {\tt QpGcut\_psi}$ (IPWpsi).
The coefficients $\alpha^{ {\bf k}n}_{au}$ can be calculated as
$$\begin{aligned}
\alpha^{ {\bf k}n}_{au} = \sum_{\bf G} C^{ {\bf k+G} }_{a u} z^{\bf k+G}_n.
\end{aligned}
\tag{I6}
$$

## Re-expansion of eigenfunctions: CPHI and GEIG

To perform the $GW$ calculation, we first
have to prepare all eigenfunctions (and eigenvalues) for
given setting of BZ mesh. Then the eigenfunctions are
represented as follows; we
re-expand ${\Psi_{ {\bf k}n} }({\bf r})$ in Eq. (I1)
as the sum of the augmentation parts in MTs and the
PW parts in the interstitial region.
$$\begin{aligned}
{\Psi_{ {\bf k}n} }({\bf r})
= \sum_{ {\bf R}u}  \alpha^{ {{\bf k} }n}_{ {\bf R}u} \varphi^{\bf k}_{ {\bf R}u}({\bf r})
 + \sum_{\bf G}  \beta^{ {{\bf k} }n}_{\bf G} P^{\bf k}_{\bf G}({\bf r}),
\end{aligned}
\tag{I7}
$$
where the interstitial plane wave (IPW) is defined as
$$\begin{aligned}
P^{\bf k}_{\bf G}({\bf r}) =
\begin{cases}
0    & \text{if {\bf r} } \in \text{any MT} \\
\exp(i ({\bf k+G})\cdot{\bf r})& \text{otherwise}
\end{cases}
\end{aligned}
\tag{I8}
$$
and $\varphi^{\bf k}_{R u}({\bf r})$ are Bloch sums of the atomic functions
$\varphi_{R u}({\bf r})$ defined within the MT at $R$,
$$\begin{aligned}
\varphi^{\bf k}_{R u}({\bf r}) & \equiv  \sum_{\bf T} \varphi_{R u}({\bf r-R-T}) \exp(i {\bf k\cdot{}T}).
\end{aligned}
\tag{I9}
$$
**T** and **G** are lattice translation vectors in real and reciprocal spaces, respectively.
Eq. (I7) is Eq. (17) in Ref.[kotani_quasiparticle_2014].
`lmf --jobgw=1` (`sugw.f90`) writes $\alpha_{ {\bf R}u}^{ {\bf k}n}$ to `__CPHI`, $\beta_{\bf G}^{ {\bf k}n}$ to `__GEIG` and the eigenvalues to `__EValue`. They are read by `readcphi` and `readgeig` in `readeigen.f90` (used e.g. in `m_zmel.f90`).

We need $\Psi_{ {\bf k}n}$ at the ${\bf k}$ points generated by **qg4gw** (the regular mesh and the offset-Γ points).
Then $P_{\bf G}^{\bf k}({\bf r})$ is just specified by
${\bf G}$, which is also generated by **qg4gw**.
$\varphi_{ {\bf R}u}^{\bf k}({\bf r})$ is specified by
radial functions. It is in `__PHIVC`, written by `sugw.f90` and read in `main_hbasfp0.f90`. Number of radial functions are
`ncore(ic)+nrad(ic)`
(dependes on $l,n,\sigma$, but not on $m$).
`lmxa` is given per species (`[[spec]] lmxa`) and may differ; the arrays are dimensioned by the maximum `lmxax` (`sugw.f90`).

In the GW calculation of ecalj, the eigenfunctions enter only through the matrix elements ([Obata2025] Eq. (8))
$$
Z^{\mathbf q-\mathbf k n',\mathbf q n}_{\mathbf k\mu} = \langle E^{\mathbf k}_\mu \Psi_{\mathbf q-\mathbf k n'} | \Psi_{\mathbf q n}\rangle ,
\tag{I11}
$$
where $E^{\mathbf k}_\mu$ is the MPB orthonormalized and diagonal to the Coulomb interaction (Eq. (I14) below).
The polarization function uses them with the labels $Z^{\mathbf k n,\mathbf q+\mathbf k n'}_{\mathbf q\mu}=\langle E^{\mathbf q}_\mu\Psi_{\mathbf k n}|\Psi_{\mathbf q+\mathbf k n'}\rangle$.
They are made on the fly by `build_zmel` and stored in `zmel` of `m_zmel.f90` ([Matrix elements used in GW](../theory/gw_zmel)).

Coefficients of Eq. (I7) (here is MTO part only)
are calculated as
$$\begin{aligned}
 &\alpha^{ {\bf k} n}_{au} = \sum_s C^{ {\bf k}s}_{a u} z^{ {\bf k}n}_s \\
 &\beta^{ {\bf k}n}_{\bf G} = \sum_{ {\bf G}'s} 
\langle P^{\bf k}_{\bf G}|P^{\bf k}_{ {\bf G}'}\rangle^{-1}
\langle P^{\bf k}_{ {\bf G}'}|H^{ {\bf k}s}\rangle z^{ {\bf k}n}_s,
\end{aligned}
\tag{I12}
$$
where the number of ${\bf G}$ is limited by the condition
$|{\bf k+G}|\lt  {\tt QpGcut\_psi}$;
${\bf G}'$ is limited by the `gmax` of `lmf` (the overlaps are made by `pwmat`, called in `sugw.f90`).
The basis $s$ includes the APWs as well as the MTOs. The coefficients are then orthonormalized by Gram-Schmidt (skipped with `--skipGS`).

## Overview of GW calculation

In the $GW$ calculation, we need not only the basis set for
eigenfunctions, but also the basis set for expanding the product of eigenfunctions.
The basis is called the mixed product basis (MPB) $\{M^{\bf k}_I({\bf r}) \}$ first introduced in Ref.[kotani_all-electron_2002] by Kotani. The MPB consists of the product basis (PB) within MTs [aryasetiawan_product-basis_1994]
and the IPW in the interstitial region.
Since $\{M^{\bf k}_I({\bf r}) \}$ contains IPWs which are not orthogonal, we have the overlap matrix and the dual basis
$$\begin{aligned}
 & O^{\bf k}_{IJ} = \langle M^{\bf k}_{I} |  M^{\bf k}_J \rangle, \\
 & |\tilde{M}^{\bf k}_{I} \rangle \equiv \sum_{J}
   |M^{\bf k}_{J} \rangle (O^{\bf k})^{-1}_{JI} .
\end{aligned}
\tag{I13}
$$

With $v_{IJ}^{\bf k}= \langle M^{\bf k}_{I} |v|  M^{\bf k}_J \rangle$, the generalized eigenvalue problem
$\sum_J (v_{IJ}^{\bf k}- v_\mu({\bf k}) O^{\bf k}_{IJ} ) z^{\bf k}_{\mu J}= 0$
gives the eigenvalues $v_\mu({\bf k})$ and the coefficients $z^{\bf k}_{\mu J}$ (`zcousq` in the code) of the orthonormal MPB
$|E^{\bf k}_\mu\rangle=\sum_J z^{\bf k}_{\mu J}|M^{\bf k}_J\rangle$, which diagonalizes the Coulomb interaction ([Obata2025] Eq. (7)),
$$\begin{aligned}
v({\bf k})=\sum_{\mu} | E^{ {\bf k} }_\mu \rangle {v_\mu({\bf k})}
\langle E^{ {\bf k} }_\mu |.
\end{aligned}
\tag{I14}
$$

**Two bases in the code.** The code keeps both the MPB $M^{\bf k}_I$ and the orthonormal $E^{\bf k}_\mu$, and changes between them with $z^{\bf k}$ (`m2e_prod_basis` of `m_zmel.f90`, set from `zcousq` by `set_m2e_prod_basis`).
`build_zmel` first makes $\langle M^{\bf k}_I\Psi_{\mathbf q-\mathbf k n'}|\Psi_{\mathbf q n}\rangle$ (`zmelt`), and returns it as it is with `is_m_basis=.true.`,
or as $Z^{\mathbf q-\mathbf k n',\mathbf q n}_{\mathbf k\mu}=\sum_I (z^{\bf k}_{\mu I})^*\langle M^{\bf k}_I\Psi_{\mathbf q-\mathbf k n'}|\Psi_{\mathbf q n}\rangle$ of Eq. (I11) with `is_m_basis=.false.`.
In `hgw` (`main_hgw.f90`), the large sums over $\mathbf k$ and bands are done in the M basis, and only the $\mathrm{ngb}\times\mathrm{ngb}$ matrices are transformed at each $\mathbf q$ and $\omega$:
$$
\begin{aligned}
\Pi^{E}_{\mu\nu} &= \sum_{IJ} (z_{\mu I})^*\,\Pi^{M}_{IJ}\, z_{\nu J}, \\
W^{{\rm c},M}_{IJ} &= \sum_{\mu\nu} z_{\mu I}\, W^{{\rm c},E}_{\mu\nu}\, (z_{\nu J})^*, \\
W^{\rm c} &= \sum_{IJ}|M_I\rangle W^{{\rm c},M}_{IJ}\langle M_J| .
\end{aligned}
\tag{I14a}
$$
- $\Pi^M$ is accumulated from `zmelt` (`x0kf_zxq` with `is_m_basis=.true.`, called from `m_screened_coulomb.f90`), turned into $\Pi^E$, and Eqs. (I16) and (I19) are evaluated in the E basis (`WVRllwR`, `WVIllwI` of `m_llw.f90` with `is_x0_m_basis`, `is_wc_m_basis` true).
  $W^{\rm c}$ is then stored in the M basis. At $\mathbf q=0$ the transformation to the M basis is done after the offset-$\Gamma$ correction (`w0w0i` with `is_wc_m_basis=.true.`).
- $\Sigma^{\rm c}$ (Eq. (I20)) uses `zmelt` in the M basis with $W^{{\rm c},M}$ (`m_sxcf_sc.f90`, `build_zmel(..., is_m_basis=.true.)`).
- $\Sigma^{\rm x}$ (Eq. (I15)) uses the E basis, where $v$ is diagonal (`build_zmel(..., is_m_basis=.false.)`). For the core states (core exchange) each atom's block of the MPB is transformed alone, since `zmelt` of a core state vanishes outside its MT.
- `hx0fp0` (`main_hx0fp0.f90`) works in the E basis throughout (`is_x0_m_basis=.false.`).

The paper [Obata2025] writes everything in the E basis; Eqs. (I14a) give the same $\Pi^E$, $W^{\rm c}$ and $\Sigma$.

For the all-electron full-potential $GW$ approximation,
Eq. (I14) is introduced in Ref.[friedrich_efficient_2010].
This corresponds to the representation in the plane wave expansion
$v({\bf k}+{\bf G},{\bf k}+{\bf G}')=\frac{4 \pi \delta_{ {\bf G}{\bf G}'} }{|{\bf k}+{\bf G}|^2}$.
$\mu=1$ corresponds to the largest eigenvalue of $v_{\mu}$, and
$v_{\mu=1}({\bf k})$ is $\sim \frac{4 \pi}{|{\bf k}|^2}$, which is related to
the divergent term discussed in Sec.([§](#offset-gamma)).

With $\langle A| B\rangle =\int d^3r A^*({\bf r}) B({\bf r})$,
the exchange self-energy is ([Obata2025] Eq. (9))
$$\begin{aligned}
\Sigma^{\rm x}_{nm}({\bf q})&=
\langle \Psi_{ {\bf q}n} |\Sigma^{\rm x} |\Psi_{ {\bf q}m} \rangle \\
&=-\sum^{\rm BZ}_{ {\bf k} }  \sum^{\rm  occ}_{n'} \sum_\mu
[Z^{\mathbf q-\mathbf k n',\mathbf q n}_{\mathbf k\mu}]^*\, v_{\mu}({\bf k})\,
Z^{\mathbf q-\mathbf k n',\mathbf q m}_{\mathbf k\mu}.
\end{aligned}
\tag{I15}
$$
Replacing the occupied valence states $\Psi_{\mathbf q-\mathbf k n'}$ with the core states gives the core exchange $\Sigma^{\rm Cx}$.
Its MPB (MPB$_{\rm C}$, written $M_{{\rm C}I}$ and $E_{{\rm C}\mu}$ here) differs from that for the valence, since it represents products of the core and valence eigenfunctions.
`gwsc` makes it in a separate step before the valence part: `hbasfp0 --job=3` makes the product basis $M_{\rm C}$ (`genallcf_v3` with `incwfin=-2`: the core states times the valence radial functions),
`hvccfp0 --job=3` makes the Coulomb matrix in it and its eigenvectors ($E_{\rm C}$, in `__Vcoud.*`, which are overwritten by the valence `hvccfp0 --job=0` afterwards),
and `hsfp0_sc --job=3` evaluates Eq. (I15) with the core states in the $E_{\rm C}$ basis.
Thus `gwsc` computes the exchange twice in each iteration: the core exchange $\Sigma^{\rm Cx}$ here (in $E_{\rm C}$), and the valence exchange $\Sigma^{\rm x}$ inside `hgw` (in $E$) after `hbasfp0 --job=0` and `hvccfp0 --job=0`.
The core exchange can be computed without the product basis, since the core states are within the MTs (a task for the future).

The screened Coulomb interaction is given in the symmetric form ([Obata2025] Eqs. (10), (11))
$$
\begin{aligned}
\epsilon_{\mu\nu}(\mathbf q,\omega) &= \delta_{\mu\nu} - \sqrt{v_\mu(\mathbf q)}\,\Pi_{\mu\nu}(\mathbf q,\omega)\sqrt{v_\nu(\mathbf q)}, \\
W_{\mu\nu}(\mathbf q,\omega) &= \sqrt{v_\mu(\mathbf q)}\,\epsilon^{-1}_{\mu\nu}(\mathbf q,\omega)\sqrt{v_\nu(\mathbf q)},
\end{aligned}
\tag{I16}
$$
where the Lindhard polarization function is ([Obata2025] Eq. (12))
$$\begin{aligned}
\Pi_{\mu \nu}({\bf q},\omega)
 &=
\sum^{\rm BZ}_{ {\bf k} } \sum^{\rm occ}_{n } \sum^{\rm unocc}_{n'}
\frac{Z^{\mathbf k n,\mathbf q+\mathbf k n'}_{\mathbf q\mu}[Z^{\mathbf k n,\mathbf q+\mathbf k n'}_{\mathbf q\nu}]^*}
{\omega-(\varepsilon_{ {\bf q+k} n'}-\varepsilon_{ {\bf k}n})+i \delta} \\
 &+ \sum^{\rm BZ}_{ {\bf k} } \sum^{\rm  unocc}_{n } \sum^{\rm occ}_{n'}
\frac{Z^{\mathbf k n,\mathbf q+\mathbf k n'}_{\mathbf q\mu}[Z^{\mathbf k n,\mathbf q+\mathbf k n'}_{\mathbf q\nu}]^*}
{-\omega-(\varepsilon_{ {\bf k}n}-\varepsilon_{ {\bf q+k} n'})+i \delta}.
\end{aligned}
\tag{I17}
$$
With time-reversal symmetry ($\Psi^*_{\mathbf k n}$ is the degenerate eigenfunction of $\Psi_{-\mathbf k n}$), the second term equals the first one with $n\leftrightarrow n'$, and ([Obata2025] Eq. (14))
$$\begin{aligned}
\Pi_{\mu \nu}({\bf q},\omega)
 &=\sum^{\rm BZ}_{ {\bf k} }  \sum^{\rm  occ}_{n} \sum^{\rm  unocc}_{n'}
Z^{\mathbf k n,\mathbf q+\mathbf k n'}_{\mathbf q\mu}[Z^{\mathbf k n,\mathbf q+\mathbf k n'}_{\mathbf q\nu}]^* \\
 &\quad\times\Big(\frac{1}{\omega-\varepsilon_{ {\bf q+k}n'}+\varepsilon_{ {\bf k}n}+i \delta} \\
 &\qquad -\frac{1}{\omega+\varepsilon_{ {\bf q+k}n'}-\varepsilon_{ {\bf k}n}-i \delta}\Big).
\end{aligned}
\tag{I18}
$$
This definition is the same as that in the Fetter-Walecka textbook.
At $\omega=0$ it is a negative definite matrix, and $\Pi_{\mu\nu}(\mathbf q,\omega)=\Pi_{\mu\nu}(\mathbf q,-\omega)$.
Without time-reversal symmetry, the negative-energy part of Eq. (I17) is accumulated as well (Eq. (I18) may be used as an approximation).

We first accumulate the anti-Hermitian part of $\Pi$ on histogram bins $\omega_j=[\omega_{{\rm L}j},\omega_{{\rm R}j}]$ of the real axis
by the tetrahedron method [rath_generalized_1975], and then obtain the Hermitian part by the Hilbert transformation ([Obata2025] Eqs. (15)–(18)):
$$
\begin{aligned}
&\mathrm{Im}[\Pi_{\mu\nu}(\mathbf q,\omega_j)] = \frac{\pi}{i}\sum^{\rm BZ}_{\mathbf k}\sum^{\rm occ}_{n}\sum^{\rm unocc}_{n'} \\
&\qquad\qquad Z^{\mathbf k n,\mathbf q+\mathbf k n'}_{\mathbf q\mu}\,[Z^{\mathbf k n,\mathbf q+\mathbf k n'}_{\mathbf q\nu}]^*\, w^{\mathbf q,\mathbf k}_{n'n}(\omega_j), \\
&w^{\mathbf q,\mathbf k}_{n'n}(\omega_j) = \sum_{\mathbf k\in{\rm BZtet}(i_{\rm tet})}\frac14\, w(n',n,\mathbf q,\mathbf k,i_{\rm tet},\omega_j), \\
&w(n',n,\mathbf q,\mathbf k,i_{\rm tet},\omega_j) \\
&\qquad\qquad = \int_{{\rm BZtet}(i_{\rm tet})} d^3k \int_{\omega_{{\rm L}j}}^{\omega_{{\rm R}j}} \delta(\omega-\varepsilon_{\mathbf q+\mathbf k n'}+\varepsilon_{\mathbf k n})\,d\omega, \\
&\mathrm{Re}[\Pi_{\mu\nu}(\mathbf q,\omega_i)] = \sum_j K_{ij}\,\mathrm{Im}[\Pi_{\mu\nu}(\mathbf q,\omega_j)].
\end{aligned}
\tag{I18a}
$$
The bins are a logarithmic mesh, dense near $\omega=0$ (`[gw] HistBin_dw`, `HistBin_ratio`; [kotani_quasiparticle_2007]).
$K_{ij}$ is the Hilbert-transformation matrix, a discretized $1/(\omega_i-\omega_j)$ that takes the size of the bins into account.
The tetrahedron weights are made by `gettetwt` (`m_tetwt.f90`).
This is more efficient and safer than calculating the real part directly.

The extended irreducible zone (EIBZ) symmetrization of Ref.[friedrich_efficient_2010] is not used at present: `eibzgen.f90` is there but not called,
`[gw] EIBZmode` is read but not used, and `main_hx0fp0.f90` runs with `ngrpx=1` (no symmetry operation). It may be used again ([IBZ and EIBZ](#ibz-and-eibz-scheme)).

The correlation part $W^{\rm c}=W-v$ is ([Obata2025] Eq. (19))
$$
W^{\rm c}_{\mu\nu}(\mathbf k,\omega) = \sqrt{v_\mu(\mathbf k)}\,(\epsilon^{-1}_{\mu\nu}(\mathbf k,\omega)-\delta_{\mu\nu})\sqrt{v_\nu(\mathbf k)},
\tag{I19}
$$
and the correlation self-energy is ([Obata2025] Eq. (20))
$$
\begin{aligned}
\Sigma^{\rm c}_{nm}({\bf q},\omega) &= \sum_{\mathbf k,n'}\sum_{\mu\nu} \int_{-\infty}^\infty \frac{i\,d\omega'}{2\pi}\,
\frac{e^{-i\delta\omega'}}{\omega-\omega'-\varepsilon_{\mathbf q-\mathbf k n'}\pm i \delta} \\
&\qquad\times [Z^{\mathbf q-\mathbf k n',\mathbf q n}_{\mathbf k\mu}]^*\,W^{\rm c}_{\mu\nu}(\mathbf k,\omega')\,Z^{\mathbf q-\mathbf k n',\mathbf q m}_{\mathbf k\nu},
\end{aligned}
\tag{I20}
$$
with $-i\delta$ for the occupied states $\mathbf q-\mathbf k n'$ and $+i\delta$ for the unoccupied ones.
In QSGW we need the Hermitian part of $\Sigma^{\rm c}_{nm}(\mathbf q,\varepsilon_{\mathbf q n})$ to obtain $V^{\rm xc}$ by [Obata2025] Eq. (3) (Eq. (Q3) of [QSGW](../theory/qsgw)).
The contour of $\omega'$ is deformed, and $\Sigma^{\rm c}$ is the sum of the imaginary-axis part $I^{\rm c}$ and the pole part $R^{\rm c}$ ([Obata2025] Eqs. (21)–(24)):
$$
\begin{aligned}
\Sigma^{\rm c}_{nm}(\mathbf q,\varepsilon_{\mathbf q n}) &= I^{\rm c}_{nm}(\mathbf q,\varepsilon_{\mathbf q n}) + R^{\rm c}_{nm}(\mathbf q,\varepsilon_{\mathbf q n}), \\
I^{\rm c}_{nm}(\mathbf q,\varepsilon_{\mathbf q n}) &= -\sum_{\mathbf k,n'}\sum_{\mu\nu}\int_{-\infty}^{\infty}\frac{d\omega}{2\pi}
 \frac{1}{\varepsilon_{\mathbf q n}-i\omega-\varepsilon_{\mathbf q-\mathbf k n'}} \\
 &\qquad\times [Z^{\mathbf q-\mathbf k n',\mathbf q n}_{\mathbf k\mu}]^*\,W^{\rm c}_{\mu\nu}(\mathbf k,i\omega)\,Z^{\mathbf q-\mathbf k n',\mathbf q m}_{\mathbf k\nu} \\
 &\simeq -\sum_{\mathbf k,n'}\sum_{i}\sum_{\mu\nu} w^{\rm I}_i(\varepsilon_{\mathbf q n}-\varepsilon_{\mathbf q-\mathbf k n'}) \\
 &\qquad\times [Z^{\mathbf q-\mathbf k n',\mathbf q n}_{\mathbf k\mu}]^*\,W^{\rm c}_{\mu\nu}(\mathbf k,i\omega^{\rm I}_i)\,Z^{\mathbf q-\mathbf k n',\mathbf q m}_{\mathbf k\nu}, \\
R^{\rm c}_{nm}(\mathbf q,\varepsilon_{\mathbf q n}) &= \sum_{\mathbf k,n'}\sum_{\mu\nu}[Z^{\mathbf q-\mathbf k n',\mathbf q n}_{\mathbf k\mu}]^* \\
 &\qquad\times W^{\rm c}_{\mu\nu}(\mathbf k,\varepsilon_{\mathbf q n}-\varepsilon_{\mathbf q-\mathbf k n'})\,Z^{\mathbf q-\mathbf k n',\mathbf q m}_{\mathbf k\nu} \\
 &\quad\times\big[\theta(\varepsilon_{\mathbf q n}-E_{\rm F})\,\theta(\varepsilon_{\mathbf q n}-\varepsilon_{\mathbf q-\mathbf k n'}) \\
 &\qquad\qquad\times\theta(\varepsilon_{\mathbf q-\mathbf k n'}-E_{\rm F}) \\
 &\qquad -\theta(E_{\rm F}-\varepsilon_{\mathbf q n})\,\theta(\varepsilon_{\mathbf q-\mathbf k n'}-\varepsilon_{\mathbf q n}) \\
 &\qquad\qquad\times\theta(E_{\rm F}-\varepsilon_{\mathbf q-\mathbf k n'})\big] \\
 &\simeq \sum_{\mathbf k,n'}\sum_{\omega^{\rm R}_i}\sum_{\mu\nu} w^{\rm R}_i(\varepsilon_{\mathbf q n}-E_{\rm F},\varepsilon_{\mathbf q-\mathbf k n'}-E_{\rm F}) \\
 &\qquad\times[Z^{\mathbf q-\mathbf k n',\mathbf q n}_{\mathbf k\mu}]^* \\
 &\qquad\times W^{\rm c}_{\mu\nu}(\mathbf k,\omega^{\rm R}_i)\,Z^{\mathbf q-\mathbf k n',\mathbf q m}_{\mathbf k\nu}.
\end{aligned}
\tag{I20a}
$$
The weights $w^{\rm I}_i(\Delta E)$ at the imaginary-axis points $i\omega^{\rm I}_i$ follow Aryasetiawan and Gunnarsson;
since $W^{\rm c}(\omega)=W^{\rm c}(-\omega)$ only $\omega^{\rm I}_i\ge0$ are used, typically about 10 points.
$\omega^{\rm R}_i$ are the boundaries of the histogram bins of Eq. (I18a), and $w^{\rm R}_i$ are the three-point interpolation coefficients of $W^{\rm c}$ at $\varepsilon_{\mathbf q n}-\varepsilon_{\mathbf q-\mathbf k n'}$, including the step functions.
Only $R^{\rm c}$ gives the lifetime (impact ionization rate); in the static $\Sigma^{\rm c}$ of QSGW, $W^{\rm c}$ in $R^{\rm c}$ is Hermitized.
In the present code the step functions are smeared by the Fermi-Dirac kernel of `[gw] t_sigmaw` (K) instead of the infinitesimal $\delta$,
and by default (`[gw] wcsmear = true`) $W^{\rm c}(\omega)$ itself is integrated over that kernel; the three-point interpolation of [Obata2025] is `wcsmear = false` ([kBT](./kBT)). The integration over the kernel (`wcsmear = true`) is an improvement made after [Obata2025]; it removes the sensitivity to sharp plasmon poles of $W^{\rm c}(\omega)$.

The first key point, given in Sec.([§](#offset-gamma)),
is the improved offset-$\Gamma$ method, which treats
the divergence of $W^{\rm c}({\bf k}\to 0,\omega)$ in Eq. (I20).
For this purpose, we define the non-divergent
effective interaction $\overline{W^{\rm c} }({\bf k}=0,\omega)$ instead of
${W^{\rm c} }({\bf k}=0,\omega)$.
Then we can take a simple discrete sum for both expressions of
Eqs.(I15) and (I20).

The second point in Sec.([§](#sigma-interpolation)) is
how to perform
an interpolation to give $V^{\rm xc}_{\bf q}$ at any ${\bf q}$ in the whole BZ,
from $V^{\rm xc}_{\bf q}$ calculated only at limited numbers of ${\bf q}$ points.
This is required in the offset-$\Gamma$ method shown in Sec.([§](#offset-gamma)),
that is, we have to calculate eigenfunctions at some ${\bf q}$ points near ${\bf q}=0$.
For the interpolation, we expand the static
nonlocal potential $V^{\rm xc}$ in [Obata2025] Eq. (3) in highly localized
MTOs in real space. Thus such MTOs are used for two purposes:
one as the basis of the eigenfunctions; and two as the basis of
expanding $V^{\rm xc}$. The interpolation procedure of $V^{\rm xc}_{\bf k}({\bf r},{\bf r}')$
becomes stabler and simpler than the
complicated interpolation procedure in
Ref.[kotani_quasiparticle_2007].
This is because we now use highly localized MTOs. `sigm.<sname>` holds $\Sigma-V_{\rm xc}$ as a matrix in the MTO basis at the irreducible $\mathbf q$;
`gwsc --mlo` (in development) holds it in the MLO representation instead ([MLO-gwsc](./mlo_gwsc)).
In the planewave-based QSGW method by Hamann and Vanderbilt [hamann09],
they expand $V^{\rm xc}$ in the maximally localized Wannier functions instead of MTOs.

In practical implementation, the LDA or GGA exchange-correlation
potential $V^{\rm xc}_{\rm LDA}$ is used to perform efficient numerical calculations.
That is, it is used in order to generate core eigenfunctions as well as radial functions within MTs
(in this paper, we use the subscript LDA even when we use GGA. "LDA/GGA" means LDA or GGA).
The difference $V^{\rm xc}-V^{\rm xc}_{\rm LDA}$ is used
for the interpolation in the BZ (explained in Sec.([§](#sigma-interpolation))),
because this difference is numerically small as long as
$V^{\rm xc}_{\rm LDA}$ roughly gives an approximation to $V^{\rm xc}$.
This procedure utilizing $V^{\rm xc}_{\rm LDA}$ to perform efficient numerical calculations
give a very weak dependence to the final numerical results in practice,
although the results formally does not depend on the LDA/GGA exchange-correlation functions.

## ${\bf q}$ and ${\bf G}$ vector generation. qg4gw

To see the theory, read Section 3.2 in [kotani_quasiparticle_2014].
${\bf k}$ points are the regular mesh points used for the integration
in the BZ. Regular mesh points are given as
$$\begin{aligned}
{\bf k}_{i_1,i2,i3}=(i1/N1)*{\bf P}_1+(i2/N2)*{\bf P}_2+(i3/N3)*{\bf P}_3
\end{aligned}
\tag{I21}
$$
, where ${\bf P}_i$ is the primitive vector $\frac{2\pi}{\rm alat} \times$ `qlat(:,i)`, the reciprocal lattice vectors ${\bf Q}_i$ (`[gw] chi_RegQbz = true`, the default, includes ${\bf k}=0$; with `false` the mesh is shifted by half a step).
Quantities $f({\bf k})$ (periodic in BZ) can
be integrated just as a sum on the regular mesh points.
However, when $f({\bf k}=0)$ is divergent, we have careful treatments.

`mkQG2` in `mkqg.f90` is the main part of **qg4gw** (`main_qg4gw.f90`). The job `--job=N` (`iq0pin`) is 1 (generate the offset-Γ points Q0P; `gwsc`), 2 (read Q0P; `job_eps` uses `--job=2 --n1n2n3eps`), 3 (epsPP), 4 (χ+-; ngp = ngc = 0), 201 (Q0P plus the Q0P along x, y, z), 10002 (skip the Γ cell), 20002 (Γ cell only), 40001 (`QforEPSL` for magnons). With `--n1n2n3eps`, the mesh `[gw] n1n2n3eps` is used if given.
The purpose is to generate the required q and G vectors at the regular mesh points and at the offset-Γ points.
It writes the BZ data `__BZDATA` and, for each q, the G vectors with $|{\bf q}+{\bf G}|$ below `[gw] QpGcut_psi` (for the eigenfunctions, `__QGpsi`) and below `QpGcut_cou` (for the IPW of the MPB, `__QGcou`).

### Make $G$ vectors: getgv2

To get ${\bf G}$ vectors, we use `getgv2` in `getgv2.f90`, whose head is

```fortran
subroutine getgv2(alat,plat,qlat,q, QpGcut,job, ng,ngvec) ! Set up a list of recip vectors within cutoff |q+G| < QpGcut a.u.
  ! job==1 -> return ng (number of G ) only
  ! job==2 -> return ng and ngvec
  ! True G is given by  G(1:3,1:ng) = 2*pi/alat * matmul(qlat * ngvec(1:3,1:ng))
```

We can use this to get $\{ {\bf G}\}$ for given $q$.
The algorism of `getgv2` is a little complicated.
We first give the lower and upper limits `n1min` $\le n_1 \le$ `n1max` (with a margin of 1),
where ${\bf G}= n_1 Q_1 + n_2 Q_2 + n_3 Q_3$, $n_2$ and $n_3$ as well.

Let us consider the three dimensional space of ${\bf x}= {\bf q}+{\bf G}$.
For given ${\bf q}$, allowed ${\bf G}$ make a set of regular mesh points $\{ {\bf q}+{\bf G}\}$.
The purpose of `getgv2` is picking up only mesh points satisfying
$|{\bf q}+{\bf G}|\lt $`QpGcut` among these mesh points.

At first, we can calculate allowed range of $n_1$ for given maximum of
$|{\bf q}+{\bf G}|$ (=`QpGcut`). Note $|{\bf x}|=|{\bf q}+{\bf G}|$=`QpGcut`
gives a sphere; we have to pick up mesh points within the sphere.
When we spefcify $n_1$, we have a plane
(allowing $n_2,n_3$ can take any values).
The range is determined by the condition that the sphere
$|{\bf q}+{\bf G}|$=`QpGcut` cross the plane specified by $n_1$
(exactly speaking, such $n_1$ is real number).
The vector normal to the plane is the external product $Q_2 \times Q_3$.

After we get the range of $n_1$, as well as $n_2,n_3$ ,
we simply test whether ${\bf q}+{\bf G}$ for $(n_1,n_2,n_3)$ is allowed or not.

### Make $G$ vectors: shortn3

Find shortest vector in modulo of $\{Q_i\}$. That is, pull back ${\bf q}$ in the 1st BZ. (`shortn3`, in `m_shortn3.f90`)
Caution; it can be not unique when ${\bf q}$ is on the BZ boundary;
then we need to know all ${\bf q}$ and degeneracy.

## Mixed Product basis {#mixed-product-basis}

The mixed product basis consists of two types of basis sets, that is
the product basis and the IPW:
$\{M^{\bf k}_I({\bf r}) \}\equiv 
\{ P^{\bf k}_{\bf G}({\bf r}), B_{ {\bf R}\mu}^{\bf k}({\bf r})\}$,
where the index $I\equiv \{ {\bf G},{\bf R}\mu\}$
classifies the members of the basis.
The PB $B^{\bf k}_{ {\bf R}\mu}({\bf r})$ is defined as
$$\begin{aligned}
  B^{\bf k}_{ {\bf R}\mu}({\bf r}) & =  
          \sum_{\bf T} B_{ {\bf R}\mu}({\bf r-R}-{\bf T}) e^{i {\bf k\cdot T} },
\end{aligned}
\tag{I22}
$$
where $B_{ {\bf R}\mu}({\bf r})$ is made from the products of radial functions.
$B_{ {\bf R}\mu}({\bf r})$ is real and
zero for outside of MT, $|{\bf r}| \gt  R$ (See Sec.([§](#product-basis-hbasfp0))).
We set up $\{B_{ {\bf R}\mu}({\bf r})\}$ so that they are orthonormalized;
$$\begin{aligned}
  \int_{|{\bf r}| < R} B_{ {\bf R}\mu}({\bf r}) B_{ {\bf R}\mu'}({\bf r}) d^3r =
  \delta_{\mu \mu'}.
  \end{aligned}
\tag{I23}
$$
In addition, it is trivial that $\{B_{ {\bf R}\mu}({\bf r})\}$ and
$\{ P^{\bf k}_{\bf G}({\bf r})\}$ are orthogonal.
The overlap matrix
$$\begin{aligned}
 O^{\bf k}_{IJ}= 
  \int_\Omega \{M_I^{\bf k}({\bf r})\}^* {M}_J^{\bf k}({\bf r}) d^3r
\end{aligned}
\tag{I24}
$$
is therefore the unit matrix in the PB block and is not diagonal only in the IPW block, since the IPWs $\{ P^{\bf k}_{\bf G}({\bf r})\}$ are not orthogonal.
Thus it is convenient to define the dual of ${M}^{\bf k}_{I}({\bf r})$ as
$\tilde{M}^{\bf k}_{I}({\bf r})$ in the manner of Eq. (I13).

Functions made from the product of eigenfunctions can be
virtually completely expanded in the basis of ${M}^{\bf k}_{I}({\bf r})$
in this manner;
$$\begin{aligned}
\begin{cases}
  \displaystyle 
  F^{\bf k}({\bf r}) = \sum_I M_I^{\bf k} ({\bf r}) F_I({\bf k}) & \cr
  \displaystyle
  F_I({\bf k}) = 
  \int_\Omega \{\tilde{M}_I^{\bf k}({\bf r})\}^* 
  F^{\bf k}({\bf r}) d^3r.
\end{cases}
  \end{aligned}
\tag{I25}
$$

The GW part uses both $M^{\bf k}_I$ and the Coulomb-orthonormalized $E^{\bf k}_\mu$ of Eq. (I14); see Eq. (I14a) for which step uses which.

### Product basis (hbasfp0) {#product-basis-hbasfp0}

We denote the radial function of atom $a$ as
$$\begin{aligned}
   u_{apl\sigma}(r) = r \phi_{apl\sigma}(r),
\end{aligned}
\tag{I26}
$$
where the index $p$ takes 1 for $\phi$, 2 for $\dot{\phi}$,
3 for local orbital as well. We do not allow $m$ dependence
($m$ is $m$ of $L=(l,m)$.) for the radial functions.

In addition, $p$ also takes the indexes of the core functions: the core and valence radial functions are read together from `__PHIVC`
(written by `lmf --jobgw=1`, `sugw.f90`) in `main_hbasfp0.f90`.
For each atom and spin, the cores come first (`ncore`, with their $l$ and $n$ in `lcindx`, `ncindx`), then the valence radial functions (`nrad`, `lindx_r`, `nindx_r`),
so that the valence $n$ is shifted by the number of cores of that $l$ (`nc_max`).
The orthonormalized radial functions $u_{apl\sigma}(r)$ are stored in `phitoto`, and the raw ones in `phitotr` (the same for the cores).

Note that the *true* radial function is
$\phi_{apl\sigma}(r) = u_{apl\sigma}(r)/r$.
Normalization is
$1=\int_0^{R_a} \{ u_{apl\sigma}(r) \}^2 dr = 
   \int_0^{R_a} \{ \phi_{apl\sigma}(r) \}^2 r^2dr$.
The function $u_{apl\sigma}(r)$ is stored in `phitoto`.

When producing the product functions,
we use spin-averaged function `phiav` given as
$$\begin{aligned}
   u_{apl}(r) = \frac{1}{N_{\rm spin} }\sum_{\sigma}
   u_{apl\sigma}(r).
\end{aligned}
\tag{I27}
$$
(See `subroutine basnfp_v2` in `basnfp.f90`; it averages the raw `phitotr`, with the sign of the spin-down function aligned to the spin-up one).
From them, we make the product functions `rprod`
$$\begin{aligned}
  \tilde{b}_{al\nu}(r) = \frac{1}{r}{u_{a p l }(r) u_{a p' l'}(r)}
  = r \phi_{a p l}(r) \phi_{a p' l'}(r),
\end{aligned}
\tag{I28}
$$
where the index $l$ runs $|l-l'| \le l \le |l+l'|$;
$\nu$ is the index of the combination $(p,p')$.
Note the *true* product functions are given as
$$\begin{aligned}
  \tilde{B}_{al\nu}(r) = \frac{1}{r}\tilde{b}_{al\nu}(r).
\end{aligned}
\tag{I29}
$$
This relation is as same as $\phi_{apl}(r)=u_{apl}(r)/r$.

Then we calculate the overlap matrix `ovmt`,
$$\begin{aligned}
  O_{\nu_1\nu_2} &= \int_0^{R_a} 
  \tilde{B}_{al\nu_1}(r) \tilde{B}_{al\nu_2}(r) r^2 dr \\
  &= \int_0^{R_a}  \phi_{a p_1 l_1}(r) \phi_{a p'_1 l'_1}(r)
                  \phi_{a p_2 l_2}(r) \phi_{a p'_2 l'_2}(r) r^2 dr
\end{aligned}
\tag{I30}
$$
and solve the eigenvalue problem of the overlap matrix,
$Oz_{\nu}=\epsilon_{\nu} z_{\nu}$, by `call rss(..)`.
(See `basnfp_v2`.)

After neglecting eigenvectors
$z_{\nu}$ with eigenvalues $\epsilon_{\nu} \lt  {\rm tolerance} \sim 10^{-4}$ (given for each $l$ by `[product_basis] pb_tolerance` ([GW input](./gwinput#product-basis)); the last value is used for the larger $l$),
we finally have the optimal product functions as the linear combinations
of the product functions as
$$\begin{aligned}
   b_{al\nu}(r) = \frac{1}{\sqrt{\epsilon_{\nu} }}
   \sum_{\nu'}\tilde{b}_{al\nu'}(r)z_{\nu'\nu},
\end{aligned}
\tag{I31}
$$
which are stored in `rprodx` and written into `__BASFP<ic>`
and used in the successive Coulomb matrix routine `hvccfp0`.
Of course, *true* product function is
$B_{al\nu}(r) = b_{al\nu}(r)/r$.

We check the normalization of the optimal product function
in standard output
(See `lbasC` and `lbas` written by `gwsc`):

    ...
    Use rs diagonalization for real symmetric 
    Diag ibx ovv=  1 0.9999999999999930D+00 eb=   0.2716113799D-01 nod=   2
    Diag ibx ovv=  2 0.9999999999999980D+00 eb=   0.4993303381D-01 nod=   3
    Diag ibx ovv=  3 0.1000000000000001D+01 eb=   0.1467546915D+00 nod=   3
    Diag ibx ovv=  4 0.9999999999999996D+00 eb=   0.4415639258D+01 nod=   0

    ...

In `basnfp`, we calculate all the required
radial integrations $\lt \phi \phi B\gt $`=ppbrd`;
$$\begin{aligned}
   \langle\phi \phi B\rangle &= \int_0^{R_a}
   \phi_{ap_1 l_1}(r) \phi_{ap_2 l_2}(r) B_{al\nu}(r) r^2 dr \\
   &= \int_0^{R_a} \frac{1}{r}
     u_{ap_1l_1}(r) u_{ap_2l_2}(r) b_{al\nu}(r) dr,
\end{aligned}
\tag{I32}
$$
which are stored into `__PPBRD_V2_<ic>`.
At `call rdpp(...` in `rdpp.f90` (module `m_rdpp`, used by `m_zmel`),
we allocate and read `ppbrd`.

In addition, we read the \"rotated Clebsch-Gordan coefficient\"
$C(L,L_1,L_2,g)$ `cgr(lm,lm1,lm2,ng)`
, where $g$ is the index for space group
coefficient ( rotated by point group symmetries).

## The Coulomb matrix (hvccfp0)

We have to calculate
$v_{IJ}^{\bf k}= \langle M^{\bf k}_{I} |v|  M^{\bf k}_J \rangle$, which
appears right after Eq. (I13).
Our `hvccfp0` (`main_hvccfp0.f90`) can handle the
case $\frac{\exp(- |\kappa| |{\bf r}_1-{\bf r}_2|)}{|{\bf r}_1-{\bf r}_2|}$.
The energy variable $E$ is given by $E=\kappa^2$, where
imaginary part of $\kappa$ is defined to be positive for negative $E$.
This $E=$`eee=screenfac()` in `main_hvccfp0.f90`; `screenfac()` in `switch.f90` returns $-$`TFscreen`$^2$.
The default is `[gw] TFscreen` $=\sqrt{10^{-5}}$, that is, $E=-10^{-5}$, nearly the bare Coulomb interaction (a too small value gives a negative eigenvalue of the Coulomb matrix at $q=0$).
For example, we can set $|\kappa|=0.1$ (a.u.) by `TFscreen = 0.1` in `[gw]` of `ctrlg.<sname>.toml` if necessary.

The MPB is made of IPWcou ($|{\bf q}+{\bf G}|$`<QpGcut_cou`) and PB,
that is, $\{M^{\bf k}_J\}=\{P_{\bf G}^{\bf k}({\bf r}), B^{\bf k}_{ {\bf R}{liL}\mu}({\bf r})\}$.
For given ${\bf k}$ (not explicitly shown in cases for simplicity),
IPWcou is specified just by the ${\bf G}$ vector.
PB can be specified just by the radial functions for each $l$
(because we have neither $m$ nor spin dependence);
PB is generated in the manner of Sec.([§](#mixed-product-basis)).
For $GW$ calculation, we need the Coulomb matrix elements
$\langle B|v|B \rangle,\langle B|v|P_{\bf G}\rangle,$ and $\langle P|v|P \rangle$.
`vcoulq_4` (`mkjp.f90`), called in `main_hvccfp0.f90`, handles these calculations.

### Spherical Bessel and related functions

We use a notation such that $X_L({\bf r})=X_l(r) Y_L(\hat{r})$; note that
their radial part is dependent only on $l$.

For given energy $E$, Bessel functions
$J_L(E,{\bf r}),H_L(E,{\bf r})$ are
(here we introduce $\sqrt{E}=\kappa$,
where $\sqrt{E}=\kappa=i|\kappa|$ for $E\lt 0$),
$$\begin{aligned}
 &{J}_L(E,{\bf r})=
j_l(\sqrt{E} r)Y_L(\hat{r}) =
\begin{cases}
j_l(i|\kappa| r)Y_L(\hat{r}) \ {\rm for} \ E<0, \\
j_l(\kappa r)Y_L(\hat{r})    \ {\rm for \ E>0},
\end{cases}
 \\
 &{H}_L(E,{\bf r})=h_l(i|\kappa| r)Y_L(\hat{r})
\ {\rm for \ E<0}, 
\end{aligned}
\tag{I33}
$$
where $j_l(z)$ and $h_l(z)$ are usual spherical Bessel and Hankel functions
which behaves
$$\begin{aligned}
 &j_l(z) \sim \frac{z^l}{(2l+1)!!}  \\
 &h_l(z) \sim \frac{-i(2l-1)!!}{z^{l+1} }
\end{aligned}
\tag{I34}
$$
at $z \sim 0$.
For convenience, we define
the **Methfessel's Bessel functions (convension)**
$\bar{J}_L=\bar{J}_lY_L$ and $\bar{H}_L=\bar{H}_lY_L$, where
$$\begin{aligned}
 &\bar{J}_l(E,{\bf r})=j_l(i|\kappa| r) / (i |\kappa|)^l  \\
 &\bar{H}_l(E,{\bf r})=h_l(i|\kappa| r) i (i |\kappa|)^{l+1}. 
%&\bar{J}_l(E,\bfr)=(2l+1)!! j_l(i|\kappa| r) / (i |\kappa|)^l  \\
%&\bar{H}_l(E,\bfr)=\frac{1}{(2l-1)!!} h_l(i|\kappa| r) i (i |\kappa|)^{l+1}. 
\end{aligned}
\tag{I35}
$$
(memo: for example, $7!!=7\cdot5\cdot3\cdot1$ `;function fac2m(i)`).
These are real functions for $E\le 0$.
See note for `genjh` in `mkjp.f90`.
At $E=0$, this is reduced to be
$$\begin{aligned}
 &\bar{J}_l(E=0,{\bf r})=r^l/(2l+1)!!,  \\
 &\bar{H}_l(E=0,{\bf r})=r^{-l-1} (2l-1)!!.
\end{aligned}
\tag{I36}
$$

(We use very small negetive $E=-10^{-5}$ as default instead of $E=0$
for Coulomb matrix generator `hvccfp0`
to avoid numerical troubles (see default `screenfac` in `switch.f90`)).

The source codes to define bessl is a little confusing because
of some convensions are mixed up\...(note at the beginning of `besslr.f90`).

Main one is `call bessl(ex2,lx,phi(0:lx),psi(0:lx))` in `mkjp.f90` (`genjh`).
This is defined in `bessl(ex2,lmax,phi,psi)` in `besslr.f90` (the module `m_bessl` has the same `bessl`, used in `vcoulq_4`).
`lx` is the upper limit of $l$.
For given `ex2`$=E\times x^2$,
this return spherical Bessel functions $j_l(\kappa x)$
in the following manner.


$$\begin{aligned}
&E<0,\ \kappa=i|\kappa|,\ |\kappa|=\sqrt{|E|}: \\
&\qquad \bar{J}_l(E,{\bf r})=j_l(i|\kappa|r)/(i|\kappa|)^l={\tt phi(l)}\,r^l, \\
&\qquad \bar{H}_l(E,{\bf r})=h_l(i|\kappa|r)\,i (i|\kappa|)^{l+1}={\tt psi(l)}/r^{l+1}, \\
&E>0,\ \kappa=\sqrt{E}: \\
&\qquad \bar{J}_l(E,{\bf r})=j_l(\kappa r)/\kappa^l = {\tt phi(l)}\,r^l, \\
&\qquad \bar{N}_l(E,{\bf r}) = {\tt psi(l)}/r^{l+1}\quad (\text{Neumann function})
\end{aligned}
\tag{I37}
$$

for $0\le l \le{\tt lmax}$. These definition gives `phi`$=1/(2l+1)!!$
and `psi`$=(2l-1)!!$ at $E \to 0$.
Here $n_l(r)$ is the spherical Neumann functions.
That is, `psi` is for Hankel function
for negative $E$, and for Neumann function for positive $E$.
See `hansmr` in `m_hansr.f90`, for example. `psi` is not used so often.

For convenience, we sometimes use
$$\begin{aligned}
 &R^{+l}=\bar{J}_l(E,r)(2l+1)!!={\tt phi(l)}r^l(2l+1)!!,   \\
 &R^{-l}=\bar{H}_l(E,r)/(2l-1)!!={\tt psi(l)}r^{-l-1}(2l-1)!!.
\end{aligned}
\tag{I38}
$$

Here note $R^{+k} \rightarrow r^k$ for $E\to 0$,and
$R^{-k} \rightarrow r^{-k-1}$ for $E\to 0$.
These $rR^{+k},rR^{-k}$ are `rkpr,rkmr` generated in subroutine `genjh`.

### Green function

(Readers can skip this subsection).
We explain the free-space Green's function
$G({\bf r}-{\bf r}',E)$ here.
Let us start from $G({\bf r}-{\bf r}',E)$, which satisfies
$$\begin{aligned}
(E+i\delta+\nabla^2)G({\bf r}-{\bf r}',E)=\delta({\bf r}-{\bf r}').
\end{aligned}
\tag{I41}
$$
($+i\delta$ is to specify the boundary condition along time axis. This pick up the retarded Green's function).
Roughly speaking, this is $(\omega -H) G=1$.
Its Fourier transform is easily written as
$G({\bf k},E)=1/(E-|{\bf k}|^2+ i \delta)$.
We apply back Fourier transformation to this, and get $G({\bf r}-{\bf r}',E)$.
It is
$$\begin{aligned}
G({\bf r}-{\bf r}',E)=-\frac{1}{4 \pi} 
\frac{e^{i \kappa|{\bf r}-{\bf r}'|} }{|{\bf r}-{\bf r}'|}, \end{aligned}
\tag{I42}
$$
where $\kappa \equiv \sqrt{E}$; Imaginary part of $\kappa$ is positive
for $E\lt 0$, that is, $\kappa = i |\kappa|= i \sqrt{-E}$ for negative $E$.
Eq. (I42) is nothing but the solution of Helmholtz
differential equation; it reduces to the usual Poisson equation at
$E=0$. For $E\lt 0$, we have Thomas-Fermi type function;
the numerator of Eq. (I42) is $\exp(i \kappa |{\bf r}-{\bf r}'|)=\exp(- |\kappa| |{\bf r}-{\bf r}'|)$.

### Used formulas {#used-formulas}

$$\begin{aligned}
\frac{1}{|{\bf r}-{\bf r}'|}
=4 \pi \sum_K \frac{r^k_<}{r^{k+1}_>}\frac{1}{2k+1}
Y^*_K(\hat{ {\bf r} })Y_K(\hat{ {\bf r} }') \end{aligned}
\tag{I43}
$$
See Appendix A in Ref.[svane_evaluation_1986].
This is generalized to be
$$\begin{aligned}
\frac{e^{-|\kappa| |{\bf r}-{\bf r}'|} }{|{\bf r}-{\bf r}'|}
&=4 \pi \sum_K R^{+k}(r_<) R^{-k}(r_>)\frac{1}{2k+1}
Y^*_K(\hat{ {\bf r} })Y_K(\hat{ {\bf r} }') \\
&=4\pi \sum_K 
\bar{J}_k(E,r_<) \bar{H}_k(E,r_>)\, Y^*_K(\hat{ {\bf r} })Y_K(\hat{ {\bf r} }'),
\end{aligned}
\tag{I44}
$$
The definition of $\bar{J},R^{+l}$ are in Eq. (I35) and Eq. (I38).
Here $E\lt 0$ and $\kappa= i|\kappa|=i\sqrt{-E}$.

$$\begin{aligned}
 &\exp( i {\bf k} {\bf r}) = 4 \pi \sum_L i^l 
j_l(|{\bf k}|r)  Y^*_L( {\bf \widehat{ k} }) Y_L(\widehat{ {\bf r} }) 
\\
 &\frac{2l+1}{4 \pi} P_l(\cos\Theta) = \sum_m Y^*_L( {\bf \widehat{r} }_1) 
Y_L({\bf \widehat{r} }_2) \ \ \ \ \ \ 
[\cos \Theta = {\bf \widehat{r} }_1 \cdot {\bf \widehat{r} }_2].
\end{aligned}
\tag{I45}
$$

The overlap of the IPWs (`mkppovl2` and `matgg2` in `pwmat.f90`) depends only on ${\bf G}'-{\bf G}$ and not on ${\bf k}$:
$$\begin{aligned}
\langle P^{\bf k}_{\bf G}|P^{\bf k}_{ {\bf G}'} \rangle
&= \Omega \delta_{ {\bf G},{\bf G}'} \\ &\quad - \sum_{a} e^{ i ({\bf G}'-{\bf G})\cdot {\bf R}_a}\,\frac{4\pi}{g^3}\big(\sin(gR_a)-gR_a\cos(gR_a)\big),
\end{aligned}
\tag{I46}
$$
where $g=|{\bf G}'-{\bf G}|$, and the term of $a$ is $4\pi R_a^3/3$ for $g=0$.

### Hankel function and Structure constant

We can expand
$v({\bf r},{\bf r}')=\frac{e^{-|\kappa||{\bf r}-{\bf r}'|} }{|{\bf r}-{\bf r}'|}$
in the one-center expansion Eq. (I44).

For the Hankel function in Eq. (I44), we use
the off-center expansion theorem of the
Hankel function Eq. (I47), that is,
a Hankel whose center is at ${\bf X}\equiv{\bf R}+{\bf T}$, $H_L({\bf r}-{\bf X})$,
can be expanded in the Bessel functions whose center is at ${\bf X}'$;
$$\begin{aligned}
\bar{H}_L(E,{\bf r}-{\bf X})=\frac{1}{4\pi} \sum_{L'}\bar{J}_{L'}(E,{\bf r}-{\bf X}') S_{ {\bf X}'{L'},{\bf X}{L} }\ ,
\end{aligned}
\tag{I47}
$$
where $E\lt 0$ (Hankel functions of negative energy).
Here $E$-dependence of $S_{ {\bf X}'{L},{\bf X}{L'} }$ is not explicitly shown.
Note the difference between $\bar{J}_L$ and $J_L$ ($\bar{H}_L$, as well).

Thus, for $({\bf R}',{\bf T}')\neq({\bf R},{\bf T})$,
we have two-center expansion;
$$
\begin{aligned}
&\frac{e^{-|\kappa||{\bf r}+{\bf R}+{\bf T}-({\bf r}'+{\bf R}'+{\bf T}')|} }
{|{\bf r}+{\bf R}+{\bf T}-({\bf r}'+{\bf R}'+{\bf T}')|} \\
&\qquad = \sum_L \sum_{L'} \bar{J}_L(E,{\bf r})\,
S_{ {\bf R}+{\bf T}{L},{\bf R}'+{\bf T}'{L'} }\, \bar{J}_{L'}(E,{\bf r}')
\end{aligned}
\tag{I48}
$$
for $({\bf R}',{\bf T}')\neq({\bf R},{\bf T})$.

The Bloch sum of $S_{ {L},{\bf X}{L'} }$
gives the structure constant of ${\bf k}$ as
$$\begin{aligned}
S^{\bf k}_{ {\bf R}{L},{\bf R}'{L'} }=
\sum_{ {\bf T}'} S_{ {\bf R}{L},{\bf R}'+{\bf T}'{L'} }\exp(i {\bf k}{\bf T}').
\end{aligned}
\tag{I49}
$$
We use $E=-10^{-5}$ by default (`screenfac`), nearly the bare Coulomb interaction.

See the top of `strxq` defined in `m_hvccfp0_util.f90` (called in `main_hvccfp0.f90`; for $E\lt 0$ all pairs at once by `strxq_all`).
This routine is for the one-center expansion of
usual Bloch summed Hankels. (not for smooth Hankels).
This result is finally converted to be the
Bloch sum of the structure constant
$S^{ {\bf k} }_{ {\bf R}{L},{\bf R}'{L'} }$ used in Eq. (I50).

`strxq` returns the structure constant without the factor $4\pi$ (its expansion theorem is $H_{RL}(E,r)=\sum J_{R'L'}(E,r)\,S_{R'L',RL}$, with $S_{R'L',RL}=4\pi\sum_{L''}C_{LL'L''}(-1)^l(-E)^{(l+l'-l'')/2}H_{L''}(E,R-R')$ in its notation, where $R$ includes $T$).
`main_hvccfp0.f90` multiplies it by $4\pi$ (`strx = fpi*s`), and `strx` is $S^{\bf k}_{{\bf R}L,{\bf R}'L'}$ of Eqs. (I48)–(I50).

### $\langle B|v|B \rangle$ part {#b-v-b-part}

Let us start from $\langle B|v|B \rangle$ part.
For this calculation, we need structure constant, and a few types of
radial integrals.
With the Bloch-summed structure constant `strx`
$=S^{\bf k}_{RL,R'L'}$, we have
$$
\begin{aligned}
v^{\bf k}({\bf r},{\bf r}')
&=\sum_{ {\bf T}'} v({\bf r},{\bf r}'+{\bf T}') e^{i {\bf k}{\bf T}'} \\
&=e^2 \sum_{L,L'} \bar{J}_L(E,{\bf r}-{\bf R})\,
S^{ {\bf k} }_{ {\bf R}{L},{\bf R}'{L'} }\, \bar{J}_{L'}(E,{\bf r}'-{\bf R}'),
\end{aligned}
\tag{I50}
$$
for ${\bf R}\ne {\bf R}'$.
Note that we use $e^2=1$ (a.u.) in `hvccfp0`.
(See note at the beginning of `mkjp.f90`. And note the normalization check
at the end of `hvccfp0`; $\langle \exp(i {\bf q}{\bf r}) |v |\exp(i {\bf q}{\bf r})\rangle=4\pi \Omega/|{\bf q}|^2$, where $\Omega$ is the cell volume.)

This `strx` is used in `call vcoulq_4`.
`nlx1` means `(l+1)**2` for $R$
`(=ibas1)`, `nlx2` as well.

Except the contribution for $({\bf R},{\bf T})=({\bf R}',{\bf T}')$,
we can evaluate
$\langle B^{\bf k}_{ {\bf R}{L}\mu}({\bf r})|v^{\bf k}({\bf r},{\bf r}')|B^{\bf k}_{ {\bf R}'{L'}\mu'}({\bf r}') \rangle$, from the\
`rojb` integrals $\rho^l(B_{ {\bf R}{l}\mu})$ as
$$
\begin{aligned}
\rho^l(B_{ {\bf R}{l}\mu})
&=\int_0^{R} r\bar{J}_l(E,r) \ r B_{ {\bf R}{l}\mu}(r)  dr \\
&= \frac{1}{(2l+1)!!}\int_0^{R} {\tt rkpr(r) \ rprodx(r)} \ dr
\end{aligned}
\tag{I51}
$$
Here $B_{ {\bf R}{l}\mu}(r)$ is the radial part of $B_{ {\bf R}{L}\mu}({\bf r})$.
This `rojb` integrals are calculated in
the subroutine `mkjb_4` in `mkjp.f90`.
We use radial functions
`rprodx` $=r B_{ {\bf R}{l}\mu}(r)$,
`rkpr`$=r\bar{J}_l(E,r)(2l+1)!!$,
and `rkmr`$=r\bar{H}_l(E,r)/(2l-1)!!$.
(Here `rkpr` and `rkmr` are propotional to
$r^l$ and $r^{-l-1}$ for $E =0$.)

The contribution from $({\bf R},{\bf T})=({\bf R}',{\bf T}')$
should be added. This is ${\bf k}$-independent, and
given by the\
`sgbb` integral, which is also calculated in `mkjb_4`.
$$
\begin{aligned}
\sigma^l(B_{ {\bf R}{l}\mu},B_{ {\bf R}{l}\nu})
&={4\pi}\int_0^{R}\!\!\int_0^{R} r_<\bar{J}_l(E,r_<)\, r_>\bar{H}_l(E,r_>) \\
&\qquad\times r B_{ {\bf R}{l}\mu}(r)\, r' B_{ {\bf R}{l}\nu}(r')\, dr dr'  \\
&=\frac{4\pi}{2l+1}\int_0^{R}\!\!\int_0^{R} {\tt rkpr(r_<) \ rkmr(r_>)} \\
&\qquad\times {\tt rprodx(n1,r)}\ {\tt rprodx(n2, r')}\, dr dr'.
\end{aligned}
\tag{I52}
$$

With the integrals `rojb` and `sgbb`, we can calculate
$\langle B|v|B \rangle$ in `vcoulq_4` as follows
(`nbloch` means the total number of PB; the loop over the pairs of PB in `vcoulq_4`):

```fortran
vcoul(ibl1,ibl2) = rojb(n1, l1, ibas1) *strx(lm1,ibas1,lm2,ibas2) *rojb(n2, l2, ibas2)    ! offsite Coulomb
if(ibas1==ibas2.AND.lm1==lm2) vcoul(ibl1,ibl2) = vcoul(ibl1,ibl2) + sgbb(n1,n2,l1, ibas1) ! sigma-type onsite parts
```

### $RL$ expansion of $|P^{\bf k}_{ {\bf G}'} \rangle$

To evaluate
$\langle P^{\bf k}_{\bf G}|v| P^{\bf k}_{ {\bf G}'} \rangle$,
we can use
$$\begin{aligned}
\bar{P}^{\bf k}_{ {\bf G} } \equiv 
\left( 1- \sum_{ {\bf R}{L} }^{l\le l_{\rm Pmax} } {P}_{ {\bf R}{L} } \right) 
e^{i ({\bf k+ G}){\bf r} }, \end{aligned}
\tag{I53}
$$
in the place of ${P}^{\bf k}_{ {\bf G} }$
as long as we use large enough $l_{\rm Pmax}$.
Here ${P}_{ {\bf R}{L} }$ denotes the projection operator
to extract the component of ${ {\bf R}{L} }$ contribution.
In fact, we use large enough $l_{\rm Pmax}$;
$l_{\rm Pmax} = 2 \times l_{\rm max}$`=2*LMXA` (`lxx = 2*lmxax` in `basnfp.f90` and `hvccfp0`),
where $l_{\rm max}$ denotes the maximum angular momentum
for the expansion of eigenfunctions within MT (maximum
$l$ cutoff for $\alpha^{ {{\bf k} }n}_{ {\bf R}u}$ in Eq. (I1)).
In the defalut setting, we use $l_{\rm Pmax}=8$
since we use $l_{\rm max}=4$.

The matrix elements $\langle P_1 |v| P_2 \rangle$
can be calculated as
$$
\begin{aligned}
\langle P_1(\phi\phi) |v| P_2(\phi\phi) \rangle
&=\sum_{\bf G_{1'} G_{1''} G_{2'} G_{2''} } 
\langle P_1(\phi\phi) | {P_{1'} } \rangle  
\langle {P}_{1'}| {P}_{1''} \rangle^{-1} \\
&\qquad\times\langle {P}_{1''} |v| {P_{2''} } \rangle  
\langle {P}_{2''}| {P}_{2'} \rangle^{-1} 
\langle {P_{2'} } | P_2(\phi\phi) \rangle,
\end{aligned}
\tag{I54}
$$
where $1\equiv ({\bf k},{\bf G_1})$ and so on.
Here $P_1(\phi\phi)$ indicates
that an IPW made from a product
of IPWs.The matrix elements
$\langle {P}_{2''}| {P}_{2'} \rangle^{-1}$
(stored into `__PPOvlp`; the inverse is taken where it is used)
and $\langle {P_{2'} } | P_2(\phi\phi) \rangle$
(stored into `__PPOvlpGG`)
are given at `rdata4gw` (a block in `sugw.f90`, that is, `lmf --jobgw=1`).

### Overlap matrix of the IPWs (`__PPOvlp*` files)

The overlap matrices of the IPWs are written by `lmf --jobgw=1` (the block `rdata4gw` in `sugw.f90`) and read by `getppx2` of the module `m_read_ppovl` (`rppovl.f90`), called in `build_zmel` (`m_zmel.f90`) for the IPW part of the matrix elements `<IPW psi|psi>`.
There are two kinds:

- $\langle P^{\bf k}_{\bf G}|P^{\bf k}_{{\bf G}'}\rangle$ for the IPWs of the MPB (`QpGcut_cou`), for each $q$ of the IBZ and the offset-Γ points (`qibze`):
  `__PPOvlp` (direct access, `ppovl(ngc,ngc)`), `__PPOvlpG` (the G vectors `ngvecc`) and `__PPOvlp.info` (sizes).
  The matrix depends only on ${\bf G}'-{\bf G}$ (Eq. (I46)), but it is stored for each $q$ since the set of G (`ngc`) depends on $q$.
- `__PPOvlpGG`: $\langle P_{\bf G}|P_{{\bf G}'}\rangle$ as a function of ${\bf G}'-{\bf G}$ (`ggg`, on the set `nvggg`), used for $\langle G({\rm eigenfunction})\,G({\rm eigenfunction})|G({\rm Coulomb})\rangle$
  (a product of the IPWs of two eigenfunctions is expanded in the IPWs of the MPB). It is for $k=0$ only. The ranges are set by
  `QpGcutggg = (2+0.01)*QpGcut_psi + QpGcut_cou + 2π/alat*dQpG` and `QpGcutgcgp = (1+0.01)*QpGcut_psi + QpGcut_cou + 2*2π/alat*dQQ`
  (the triangle inequality; `dQpG`, `dQQ` enlarge the ranges for the offset-Γ points).

### $\langle {P}^{\bf k}_{\bf G}|v| {P}^{\bf k}_{ {\bf G}'} \rangle$ {#p-v-p-part}

To evaluate Eq. (I54), we need to know
its main part $\langle {P}^{\bf k}_{\bf G}|v| {P}^{\bf k}_{ {\bf G}'} \rangle$. It is written as
$$
\begin{aligned}
\langle {P}^{\bf k}_{\bf G}|v| {P}^{\bf k}_{ {\bf G}'} \rangle
&\approx \langle \bar{P}^{\bf k}_{\bf G}|v| \bar{P}^{\bf k}_{ {\bf G}'} \rangle \\
&= \langle e^{i({\bf k}+{\bf G}){\bf r}} |v| e^{i({\bf k}+{\bf G}'){\bf r}} \rangle
-\sum_{ {\bf R}{L} } \langle P^{ {\bf k}+{\bf G} }_{ {\bf R}{L} }|v| e^{i({\bf k}+{\bf G}'){\bf r}} \rangle  \\
&\quad - \sum_{ {\bf R}'{L'} } \langle e^{i({\bf k}+{\bf G}){\bf r}} |v| P^{ {\bf k}+{\bf G}'}_{ {\bf R}'{L'} }\rangle \\
&\quad +\sum_{ {\bf R}{L} } \sum_{ {\bf R}'{L'} } 
\langle P^{ {\bf k}+{\bf G} }_{ {\bf R}{L} }|v|P^{ {\bf k}+{\bf G}'}_{ {\bf R}'{L'} } \rangle,
\end{aligned}
\tag{I55}
$$
where $P^{ {\bf k}+{\bf G} }_{ {\bf R}{L} }$ denotes
the projection of PW to ${ {\bf R}{L} }$,
That is,
$P^{ {\bf k}+{\bf G} }_{ {\bf R}{L} } \equiv 
{P}_{ {\bf R}{L} } e^{i ({\bf k+ G}){\bf r} }$.

##### The first term

The first term in the right-hand side of Eq. (I55) is
$$\begin{aligned}
\langle \exp(i({\bf k}+{\bf G}){\bf r}) |v| \exp(i({\bf k}+{\bf G}'){\bf r}) \rangle
=\frac{4 \pi \Omega}{|{\bf k}+{\bf G}|^2+|E|} \delta_{ {\bf G}{\bf G}'},
\end{aligned}
\tag{I56}
$$
because we simply use $\exp(i({\bf k}+{\bf G}){\bf r})$
(no prefactor for normalization) for IPW.
Here $E$`=eee` is negative (or (almost) zero).
This is coded by a line in `vcoulq_4` (`mkjp.f90`)

    vcoul(nbloch + ig1, nbloch + ig1) = vcoul(nbloch + ig1, nbloch + ig1) + fpivol/(absqg2(ig1) - eee)

##### The second and third term

In the second term,
we can replace $v$ with $\frac{4 \pi \Omega}{|{\bf k}+{\bf G}'|^2+|E|}$
since $v$ is diagonal for $|\exp(i({\bf k}+{\bf G}'){\bf r}) \rangle$, the third term as well.
Without $v$, we have
$$
\begin{aligned}
\langle e^{i({\bf k}+{\bf G}){\bf r}}| P^{ {\bf k}+{\bf G}'}_{ {\bf R}'{L'} }\rangle
&=\sum_{ {\bf R}{L} } \langle P^{ {\bf k}+{\bf G} }_{ {\bf R}{L} }| P^{ {\bf k}+{\bf G}'}_{ {\bf R}'{L'} }\rangle  \\
&=\sum_{ {\bf R}{L} } \big({\tt pjyl\_}({\bf k}+{\bf G},L)\, e^{i({\bf k}+{\bf G}){\bf R}}\big)^* \\
&\qquad\times R^{JJ}(|{\bf k}+{\bf G}|,|{\bf k}+{\bf G}'|,l) \\
&\qquad\times {\tt pjyl\_}({\bf k}+{\bf G}',L)\, e^{i({\bf k}+{\bf G}'){\bf R}},
\end{aligned}
\tag{I57}
$$
where we use
$$\begin{aligned}
 &P^{ {\bf k}+{\bf G} }_{ {\bf R}{L} }({\bf r})=4 \pi i^l j_l(|{ {\bf k}+{\bf G} }|r) 
Y_L(\widehat{ {\bf k}+{\bf G} }) Y_L(\hat{ {\bf r} })\exp(i({\bf k}+{\bf G}){\bf R}) \\
 &= {\tt pjyl\_}({\bf k}+{\bf G},L) \bar{J}_l(|{ {\bf k}+{\bf G} }|r) Y_L(\hat{ {\bf r} })
\exp(i({\bf k}+{\bf G}){\bf R}), 
\end{aligned}
\tag{I58}
$$
where ${\bf r}$ is measured from the center ${\bf R}$.
Here we use `pjyl_` defined as
$$\begin{aligned}
{\tt pjyl\_}({\bf k}+{\bf G},L) 
=4 \pi i^l |{\bf k}+{\bf G}|^l Y_L(\widehat{ {\bf k}+{\bf G} }) 
\end{aligned}
\tag{I59}
$$
(recall the definition of $\bar{J}_l$.
In codes, `cy(lm)*yl(lm)`
$=Y_L(\widehat{ {\bf k}+{\bf G} })$.)
Search `pjyl_` in `mkjp.f90`.
The Bessel functions appear here in the expansion of PW; see Eq. (I45).
$R^{JJ}(|{\bf k}+{\bf G}|,|{\bf k}+{\bf G}'|,l)$
is given as
$$\begin{aligned}
R^{JJ}(|{\bf k}+{\bf G}|,|{\bf k}+{\bf G}'|,l)
= \int_0^R  r^2 \bar{J}_l(|{ {\bf k}+{\bf G} }|r)
\bar{J}_l(|{ {\bf k}+{\bf G}'}|r) dr ,
\end{aligned}
\tag{I60}
$$
which can be calculated by the wronskian (`wronkj`) by the formula
$$
\begin{aligned}
R^{JJ}(\kappa_{\rm A},\kappa_{\rm B},l)
&=\int_0^R \!\!r^2 \bar{J}_l(\kappa_{\rm A} r) \bar{J}_l(\kappa_{\rm B} r) dr \\
&=\frac{R^2}{\kappa_{\rm A}^2-\kappa_{\rm B}^2}\Big(
\bar{J}_l(\kappa_{\rm A} r) \frac{d \bar{J}_l(\kappa_{\rm B} r)}{dr} \\
&\qquad -\frac{d \bar{J}_l(\kappa_{\rm A} r)}{dr}\bar{J}_l(\kappa_{\rm B} r)
\Big)\Big|_{r=R} = {\tt -fjj}
\end{aligned}
\tag{I61}
$$

(`-fjj` is used in `mkjp.f90`.)
In `vcoulq_4` the second, third and the on-site part of the fourth term of Eq. (I55) are added together for each group of atoms with the same radial mesh (`eee = 0`; for `eee < 0` they come from `vcoul_termb`):

```fortran
vcoul(nbloch+ig1,nbloch+ig2) = vcoul(nbloch+ig1,nbloch+ig2) &
  + cPhi * sum( dconjg(pjyl_(1:lm2x,ig1)) * pjyl_(1:lm2x,ig2) &
    * ((fpi/(absqg2(ig1)-eee)+fpi/(absqg2(ig2)-eee))*fjj(llx(1:lm2x)) + radsig(llx(1:lm2x))) )
```

where `cPhi` is the sum of the phase factors $e^{-i({\bf k}+{\bf G}_1){\bf R}}e^{i({\bf k}+{\bf G}_2){\bf R}}$ over the atoms of the group,
`fpi` is $4\pi$, `absqg2(ig)` is $|{\bf k}+{\bf G}|^2$, `-eee` is $|E|$, and `radsig` is the $\sigma$-type integral of Eq. (I62) (`sigintpp`).

##### The forth term

The last term of Eq. (I55) can be calculated essentially
the same manner with $\langle B| v| B\rangle$, where we use the
Bessel function instead of $B_{ {\bf R}{l}\mu}(r)$ appeared in Sec.([§](#b-v-b-part)).
Then we define integrals
`rojp` and `sgpp` defined as (in `mkjp.f90`; `sgpp` is not stored as an array now, see `radsig` above);
$$
\begin{aligned}
{\tt rojp}({\bf R}L)&= {\tt pjyl\_}\, \exp(i ({\bf q}+{\bf G}) {\bf R})\,  \rho^l(\bar{J}_l), \\
{\tt sgpp}({\bf R}L,{\bf G},{\bf G}')&= {\tt pjyl\_^*(ig1)} \exp(-i ({\bf q}+{\bf G}) {\bf R}) \\
&\qquad\times {\tt pjyl\_(ig2)} \exp(i ({\bf q}+{\bf G}') {\bf R})\, {\tt radsig}, \\
{\tt radsig}&=\sigma^l(\bar{J}_l(|{\bf q}+{\bf G}|r),\bar{J}_l(|{\bf q}+{\bf G}'|r)).
\end{aligned}
\tag{I62}
$$
Thus, `rojp` and `sgpp` made of coefficients for explansion
and radial integral.
The off-site part of the fourth term is $\sum$ `rojp`$^*$ `strx` `rojp` (`zmm` in `vcoulq_4`).

### $\langle P|B \rangle$ part

$$
\begin{aligned}
\langle {P}^{\bf k}_{ {\bf G} } |v| B^{\bf k}_{ {\bf R}{L}\mu} \rangle
&=\langle  e^{i({\bf k}+{\bf G}){\bf r}} |v|B^{\bf k}_{ {\bf R}{L}\mu} \rangle \\
&\quad -\sum_{ {\bf R}'{L'} } \langle P^{ {\bf k}+{\bf G} }_{ {\bf R}'{L'} }|v| B^{\bf k}_{ {\bf R}{L}\mu} \rangle  \\
&={\tt fouvb}({\bf G},{\bf R}L \mu) - {\tt sgpb}({\bf G},{\bf R}L \mu) \\
&\quad -\sum_{ {\bf R}'L'} {\tt rojp}^*\, {\tt strx}\, {\tt rojb}
\end{aligned}
\tag{I63}
$$
The first term is stored in `fouvb(ngc,nxx,nlxx,nbas)`, made by `mkjp_4` (`mkjp.f90`) together with `rojp` and `sgpb`, and used in `vcoulq_4`.
Since $v$ is diagonal to PWs,
we can evaluate this in the similar manner of Eq. (I57).
The second term can be evaluated from
`rojb` and `sgpb` in the same manner of last section: the on-site part ${\bf R}'={\bf R}$ is `sgpb`, and the off-site part is $\sum$ `rojp`$^*$ `strx` `rojb` (`zmm` in the block `PvB2` of `vcoulq_4`).

## Offset-$\Gamma$ method; W(k=0) averaged in the $\Gamma$ cell. {#offset-gamma}

The offset-$\Gamma$ method, originally
invented for Ref.[kotani_all-electron_2002] by Kotani
(it is described in Ref.[kotani_quasiparticle_2007]),
was a key to perform accurate $GW$
calculation in our papers.
It is for the integration of ${\bf k}$
in Eqs.(I15) and (I20), where
we have the integrands that diverge at ${\bf k}\to 0$.
The original offset $\Gamma$ method works well
for highly symmetric systems; however,
it may be problematic to apply to less symmetric systems,
because the anisotropic divergence of the integrands
may not be treated accurately.

Here we show an improved offset-$\Gamma$ method, which treats
the anisotropy of $W({\bf k},\omega)$ accurately.
In the followings, we use expression $W({\bf k})$
for simplicity (omit subscripts and $\omega$)
instead of $W_{\mu \nu}({\bf k},\omega)$,
since we are concerned with the ${\bf k}$ integral here.

Let us give a formula for calculating
$\int_{\rm BZ} f({\bf k}) d^3k$ using a discrete sum on ${\bf k}$-mesh,
where $f({\bf k})=G({\bf q}-{\bf k}) \times W({\bf k})$.
For the ${\bf k}$-mesh, we use
$$
{\bf k}(i_1,i_2,i_3) = 2 \pi \Big(\frac{i_1}{N_1} {\bf b}_1 
+ \frac{i_2}{N_2} {\bf b}_2 + \frac{i_3}{N_3} {\bf b}_3\Big),
\tag{I64}
$$
where ${\bf b}_1, {\bf b}_2$, and ${\bf b}_3$ are the primitive reciprocal
vectors (the same as the Eq.(47) in Ref.[kotani_quasiparticle_2007]).
The 1st BZ is divided into $N=N_1 \times N_2 \times N_3$
microcells ($i_1=0,1,... N_1-1$, and also the same for $i_2$ and $i_3$.).
The microcell including the $\Gamma$ point is called the
$\Gamma$ cell [freysoldt_dielectric_2007].
The main problem is how to evaluate the contribution of the $\Gamma$ cell.
The divergent part of $f({\bf k})$ behaves $\approx$ (analytic function of ${\bf k}$)
$/({\bf k}^{\rm T}{\bf L}{\bf k})$,
where ${\bf k}^{\rm T}$ denotes the transpose of ${\bf k}$; ${\bf L}$
is a $3\times3$ Hermitian matrix [friedrich_efficient_2010].
We neglect an odd part of ${\bf k}$ in the above (analytic function of ${\bf k}$)
because it has no contribution to the integral around ${\bf k}=0$.
Thus it is sufficient to consider the integral for $f({\bf k})$ whose divergent parts behave
as $f({\bf k})=\sum_L \frac{f_L Y_L(\widehat{ {\bf k} })}{|{\bf k}|^2}$ at ${\bf k}\to 0$,
where $l$ of $L\equiv(l,m)$ is restricted to be even numbers.
We evaluate the integral using the formula
$$\begin{aligned}
\int_{\rm BZ} f({\bf k}) d^3k \approx \frac{1}{N}\sum^{ {\bf k}\ne 0} f({\bf k})
 + \sum_L f_L w_L + \frac{1}{N} \tilde{f},
\end{aligned}
\tag{I65}
$$
which is introduced in Ref.[freysoldt_dielectric_2007].
Here the weight $w_L$ is determined in a manner as follows,
so as to take into account the contributions of
the divergent part of $f({\bf k})$ at ${\bf k}\to 0$ in the $\Gamma$ cell.
$\tilde{f}$ is the constant part of $f({\bf k})$ at ${\bf k}\to 0$.

To determine $w_L$, we can use the following procedure instead of that given
in Ref.[freysoldt_dielectric_2007].
We first introduce the auxiliary function
$$\begin{aligned}
F_L({\bf k}) =\sum_{ {\bf G} } \frac{\exp(- \alpha
 |{\bf k}-{\bf G}|^2)Y_L(\widehat{ {\bf k}-{\bf G} })}{|{\bf k}-{\bf G}|^2}.
\end{aligned}
\tag{I66}
$$
This is a generalization of an auxiliary function used in the
offset-$\Gamma$ method (then we only used $F_{00}$ [kotani_quasiparticle_2007]).
We usually take the $\alpha \to 0$ limit, or a sufficiently small $\alpha$ instead (`[gw] alpha_offG`, read in `mkqg.f90`; `m_q0p.f90`).
Let us apply Eq. (I65) to $F_L({\bf k})$.
Then we can evaluate the left-hand side of Eq. (I65)
exactly (the exact values are zero except for $L=(0,0)$).
On the other hand, the first and
third terms on the right-hand side of Eq. (I65)
can be evaluated numerically.
In addition, we know that $f_{L'}$ for $F_L({\bf k})$ is unity for $L'=L$,
and zero otherwise.
Thus we can determine $w_L$ in Eq. (I65)
so that Eq. (I65) is exactly satisfied for $F_L({\bf k})$ for any $L$.

Let us apply Eq. (I65) to $f({\bf k})=G({\bf q}-{\bf k}) \times W({\bf k})$.
Then we perform an approximation taking only the most divergent term
in $W({\bf k})$ in addition to its analytic part. That is, we use
$$\begin{aligned}
W_{\mu \nu}({\bf k}) \sim \widetilde{W}_{\mu \nu}({\bf 0}) +
 \frac{4 \pi}{ {\bf k}^{\rm T} {\bf L} {\bf k} }
 \delta_{1\mu}\delta_{1\nu}
\end{aligned}
\tag{I67}
$$
at ${\bf k}\to 0$.
$\widetilde{W}_{\mu \nu}({\bf 0})=0$ for $\mu=1$ or $\nu=1$.
See Eq.(36) in Ref.[friedrich_efficient_2010] to know what is
neglected in the approximation of Eq. (I67).

Then we finally obtain
$$\begin{aligned}
\int_{\rm BZ} d^3k G({\bf q}-{\bf k}) W({\bf k}) 
\approx \overline{\sum G({\bf q}-{\bf k}) W({\bf k})},\end{aligned}
\tag{I68}
$$
where its right-hand side is defined as
$$\begin{aligned}
 &\overline{\sum G({\bf q}-{\bf k}) W({\bf k})}  \\
 &\equiv 
\frac{1}{N} \sum_{ {\bf k}\ne 0} G({\bf q}-{\bf k}) W({\bf k}) + \frac{1}{N} G({\bf q})\overline{W}({\bf 0}),\\
 &\overline{W}({\bf 0})\equiv N\sum w_L W_L + \widetilde{W}({\bf 0}).
\end{aligned}
\tag{I69}
$$
Here $\overline{W}({\bf 0})$ is an
average of $W$ in the $\Gamma$ cell.
With this $\overline{W}({\bf 0})$, we can
evaluate integrals just as the sum on the discrete ${\bf k}$-mesh.
When the matrix ${\bf L}$ is given
(a method of calculating ${\bf L}$ is given in the next paragraph),
the non-analytic (but non-divergent) function
${ {\bf k}^{\rm T}{\bf L}{\bf k} }/|{\bf k}^2|$
is expanded in the spherical harmonics. Then
$W_L$ is calculated for a given ${\bf L}$ in the manner shown in Ref.[friedrich_efficient_2010].
We can evaluate the accuracy of integrals with a discrete ${\bf k}$-mesh
in combination with the approximation of Eq. (I67)
by calculations while changing the size of the ${\bf k}$-mesh.\

The code uses only $w_L$ for $L=(0,0)$ ($w_L=0$ for $L\ne(0,0)$; `lxklm = 0` in `m_q0p.f90`, which writes the weights `wklm`).
With all $w_L$, the real part of $W(\omega)$ became non-monotonic (against causality), which came from the $w_L$ of $L\ne(0,0)$.
In the case of metals at $\omega=0$, ${\bf L}$ is very large (virtually infinite), and there is no $1/|{\bf k}|^2$ behavior.

The remaining problem is how to calculate
the matrix ${\bf L}$ in Eq. (I67); there are two possible ways.
One is the ${\bf k}\cdot {\bf p}$ method (perturbation) used in Ref.[friedrich_efficient_2010];
the other is the numerical method to calculate ${\bf L}$ at some
${\bf k}$ points near ${\bf k}=0$. Here we use the latter method.
Because of the point-group symmetry of the system, ${\bf L}$ can be expressed by the
linear combination of invariant tensors $\mu_{ij}^g$ for the symmetry of the unit cell,
$$\begin{aligned}
L_{ij}(\omega)=  \sum_{g=1}^{N_g} a_g(\omega) \mu_{ij}^g,
\end{aligned}
\tag{I70}
$$
where $g$ is the index of the invariant tensor.
The number of $g$'s $N_g$, can be
from one (cubic symmetry) through six (no symmetry).
It is possible to determine the coefficient $a_g(\omega)$
from the dielectric functions ${\hat{ {\bf k} }_{0i} }^{\rm T}{\bf L}\hat{ {\bf k} }_{0i}$
calculated at $\{ {\bf k}_{0i}\}$ points around
${\bf k}=0$, where $\{ {\bf k}_{0i}; i=1,N_g\}$ is a set of the offset-$\Gamma$ points.
The offset-$\Gamma$ points are chosen so that the conversion matrix
from $\hat{ {\bf k} }_{0i}^{\rm T}{\bf L}(\omega)\hat{ {\bf k} }_{0i}$ to $a_g(\omega)$
is not numerically degenerated. The length $|{\bf k}_{0i}|$ (`[gw] Q0Pchoice`, or `deltaq_scale` given directly)
can be chosen to be sufficiently small, but avoiding numerical error
as the average of $W({\bf k})$ in the $\Gamma$ cell.
The improved offset-$\Gamma$ method shown here can be applicable even to metal cases, as long as
$\hat{ {\bf k} }_{0i}^{\rm T}{\bf L}(\omega)\hat{ {\bf k} }_{0i}$ contains the contribution
of intraband transition.

## Data structure of the tetrahedron method in ecalj {#tetrahedron-data-structure}

The tetrahedron weights of the polarization function, Eqs. (I18a), are made by `gettetwt` of the module `m_tetwt` (`m_tetwt.f90`):

```fortran
call gettetwt(q, iq, isp_k, isp_kq, ekxx1, ekxx2, nband=nband, ikbz_in=k_lo, fkbz_in=k_hi)   ! in x0kf_v4h.f90
```

- `q`, `iq`: the vector ${\bf q}$ and its index.
- `isp_k`, `isp_kq` ($\alpha$, $\beta$): the spins of ${\bf k}$ and ${\bf k}+{\bf q}$ (different for the transverse spin susceptibility).
- `ekxx1(nband,nqbz)`, `ekxx2(nband,nqbz)`: the eigenvalues at ${\bf k}$ and ${\bf k}+{\bf q}$ for all ${\bf k}$ of the BZ (`readeval`).
- `ikbz_in`, `fkbz_in`: the range of ${\bf k}$ to be done (the ${\bf k}$ loop is split for memory and for the parallelization).
- The other inputs come from modules: the histogram bins `frhis(1:nwhis+1)` and `npm` of `m_freq` (`npm = 2` also accumulates the negative-frequency part, needed without time-reversal symmetry),
  the tetrahedra (`ntetf`, `idtetf`, `ib1bz`) of `m_read_bzdata`, the core levels of `m_core_state`, $E_{\rm F}$ of `m_ReadEfermi`, and `ebmx`, `nbmx`, `mtet` of `m_readgwinput`.
  With `[gw] t_tetrakbt`, the Fermi level of the finite-temperature occupation is used ([kBT](./kBT)).

The output is in the module variables:

- `nbnb(k,jpm)`: the number of pairs $(n,n')$ with non-zero weight at ${\bf k}$; `n1b(ibib,k,jpm)` (occupied, cores after the valence) and `n2b(ibib,k,jpm)` (unoccupied) are their band indexes.
- `ihw(ibib,k,jpm)`, `nhw(ibib,k,jpm)`: the first histogram bin and the number of bins with non-zero weight for the pair;
  `jhw(ibib,k,jpm)` points to the weights `whw(jhw : jhw+nhw-1)`.

They are used in `x0kf_v4h.f90` as (simplified)

```fortran
do jpm = 1, npm
  do ibib = 1, nbnb(k,jpm)
    it  = n1b(ibib,k,jpm)                      ! n  at k
    itp = n2b(ibib,k,jpm)                      ! n' at k+q
    ! M(I,J) = conjg(zmel(I,it,itp)) * zmel(J,it,itp)   (zmel from build_zmel)
    do iw = ihw(ibib,k,jpm), ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1
      imagweight = whw(jhw(ibib,k,jpm)+iw-ihw(ibib,k,jpm))
      ! Im Pi(:,:,iw,jpm) += M(:,:)*imagweight
    enddo
  enddo
enddo
```

which is the sum over $n,n'$ and ${\bf k}$ of Eq. (I18a).

## Fourier transformation of non-local quantity

The static self-energy in the MTO basis is Fourier transformed in `m_rdsigm2` (`rdsigm2.f90`): `rdsigm2` expands `sigm.<sname>` to the full BZ, `fftz3` gives the real-space matrix `hrr`, and `bloch2` transforms it back at any ${\bf k}$.
This is for the three dimensional FFT.
The usual FT is by
$$\begin{aligned}
f({\bf T})=\sum_{\bf k}f({\bf k}) \exp(i {\bf k}{\bf T}),
\end{aligned}
\tag{I71}
$$
where $\{ {\bf k}\}$ is are on the regular mesh points.
The total number of its members is $N_1 \times N_2 \times N_3$.
(the number is the same as that of $\{ {\bf T}\}$).
Note that we have periodicity both in ${\bf k}$ points
and in ${\bf T}$ points. Because of the periodicity, the range of $\{ {\bf k}\}$
is not unique, the range of $\{ {\bf T}\}$ as well.

Let us think about non-local quantity which is dependent on ${\bf T}-{\bf T}'$. Then we have
$$\begin{aligned}
f({\bf R}{\bf T},{\bf R}'{\bf T}')=\sum_{\bf k}
f_{ {\bf R}{\bf R}'}({\bf k}) \exp(i {\bf k}({\bf T}-{\bf T}')),
\end{aligned}
\tag{I72}
$$
In practical calculations (static version of self-energy treated by
`bloch2`),
we first calculate $f_{ {\bf R}{\bf R}'}({\bf k})$ on ${\bf k}$ of regular mesh points. Then we
need to obtain its real-space representation
$f({\bf R}{\bf T},{\bf R}'{\bf T}')$. Because of the periodicity,
we have ambiguity for the choice of possible
$|{\bf T}-{\bf T}'|$. If we introduce $\bar{\bf T}={\bf T}-{\bf T}'$,
Eq. (I72) is written as
$f({\bf R}\bar{ {\bf T} },{\bf R}'0)=\sum_{\bf k}
f_{ {\bf R}{\bf R}'}({\bf k}) \exp(i {\bf k}\bar{ {\bf T} })$ because of translational symmetry.

A reasonable choice is that we allow $\bar{ {\bf T} }$
which satisfy $|{\bf R}-{\bf R}'+\bar{ {\bf T} }|\le\eta_{\rm FTmax}$.
Here we should choose $\eta_{\rm FTmax}$ so that
the number of alllowed $\{\bar{\bf T}\}$ is $N_1 \times N_2 \times N_3$.
Howerever, it can be not possible, because of
deneneracy; for the largest value of $|{\bf R}-{\bf R}'+\bar{ {\bf T} }|$ in the
alllowed $\{\bar{\bf T}\}$, we may
have some of $\bar{ {\bf T} }$ (we say degenerated).
Then we need to give fractional weight for such $\bar{ {\bf T} }$.

In the code, the list of $\bar{\bf T}$ is made automatically: for each pair of sites $({\bf R},{\bf R}')$ and each ${\bf T}$ of the $N_1\times N_2\times N_3$ cell,
`m_gennlat_sig` (`m_gennlat.f90`) keeps the shortest equivalent vectors $|{\bf R}-{\bf R}'+\bar{\bf T}|$ (`nlat`, with the range `nlatS:nlatE`),
and `bloch2` averages the phase factors over the degenerate ones (equal fractional weights).

## Interpolation of the self-energy in the Brillouin zone {#sigma-interpolation}

Here we show
an interpolation procedure for giving
$V^{\rm xc}_{ {\bf k} }$ at any ${\bf k}$,
from $V^{\rm xc}_{ {\bf k} }$ calculated only at the regular mesh points ${\bf k}(i_1,i_2,i_3)$.
This interpolation is used for the offset-$\Gamma$ method that
requires $W(\omega)$ at $\{ {\bf k}_{0i}\}$;
to calculate this $W(\omega)$, we need
eigenfunctions and eigenvalues not only
at the regular mesh points ${\bf k}(i_1,i_2,i_3)$
but also at ${\bf k}(i_1,i_2,i_3)+{\bf k}_{0i}$.
This interpolation is also useful for plotting
energy bands.
A key point of the interpolation is that $V^{\rm xc}$
is expanded in real space in highly localized MTOs as follows.

At the end of a QSGW iteration (`hqpe_sc`), we obtain the matrix elements
$\langle {\Psi_{ {\bf k}n} }| {\it \Delta} V^{\rm xc}_{\bf k}| {\Psi_{ {\bf k}m} }\rangle$ on the regular mesh
points of ${\bf k}$, where
${\it \Delta} V^{\rm xc}_{\bf k}=V^{\rm xc}_{\bf k}-V^{\rm xc,LDA}_{\bf k}$.
Then it is converted to the representation in the APW and MTO bases as
$$\begin{aligned}
\langle \chi^{ {\bf k} }_a| {\it \Delta} V^{\rm xc}_{\bf k}| \chi^{ {\bf k} }_b \rangle
= \sum_{n,m} \left(z^{-1}\right)^*_{an}\langle {\Psi_{ {\bf k}n} }|{\it \Delta} V^{\rm xc}_{\bf k}|
{\Psi_{ {\bf k}m} }\rangle z_{bm}^{-1},  \\
\end{aligned}
\tag{I73}
$$
where we use the simplified basis index $a$, which is the index for specifying a basis
(${\bf R}L{j}$ for MTO or ${ {\bf G} }$ for APW).
Thus $\chi^{ {\bf k} }_a$ denotes the APWs or MTOs in Eq. (I1);
$z_{na}$ (${\bf k}$ is omitted for simplicity)
denotes the coefficients of the eigenfunctions at ${\bf k}$, that is,
$z^{ {{\bf k} }n}_{ {\bf R}L{j} }$ and $z^{ {\bf k}n}_{ {\bf G} }$ in Eq. (I1) together.
This $z_{an}$ is identified as a conversion matrix that connects
eigenfunctions (band index $n$) and the APW and MTO bases (basis index $a$).

To obtain real-space representation of ${\it \Delta} V^{\rm xc}$, we need a representation
expanded in the basis that consist of the Bloch-summed localized orbitals,
which are periodic for ${\bf k}$ in the BZ. However, this is not the case for the APWs in
Eq. (I73). To overcome this problem, we use an approximation
in which we only take the matrix elements
related to MTOs, that is, the elements
$\langle \chi^{ {\bf k} }_a| {\it \Delta} V^{\rm xc}_{\bf k}| \chi^{ {\bf k} }_b \rangle$
where $a$ and $b$ specify MTOs.
This means that the part of ${\it \Delta} V^{\rm xc}$ related to APWs is projected onto the basis of MTOs.
This approximation can be reasonable as long as the main part of ${\it \Delta} V^{\rm xc}$ can be
well expanded in MTOs,
although we need numerical tests to confirm the accuracy.
Then we obtain a real-space representation of
${\it \Delta} V^{\rm xc}$ expanded in MTOs from the MTO part of
$\langle \chi^{ {\bf k} }_a| {\it \Delta} V^{\rm xc}_{\bf k}| \chi^{ {\bf k} }_b \rangle$ by Fourier transformation.
Then we can have interpolated ${\it \Delta} V^{\rm xc}$ at any ${\bf k}$ by inverse Fourier transformation.
Since we use highly localized MTOs,
this interpolation is more stable
than the previous one in FP-LMTO-QSGW [kotani_quasiparticle_2007].
The complicated interpolation procedure
given in Sec.II-G in Ref.[kotani_quasiparticle_2007] is no
longer necessary.

To reduce the computational time,
we calculate the matrix elements $\langle {\Psi_{ {\bf k}n} }|{\it \Delta} V^{\rm xc}_{\bf k}| {\Psi_{ {\bf k}m} }\rangle$ only up to the states whose
eigenvalues are less than $E_{\rm MAX}^\Sigma$ (`[gw] emax_sigm`).
Then the high energy parts of the matrix elements are assumed
to be diagonal, where their values are given by an energy-weighted average (`ESEAVR` in `hqpe_sc`)
of calculated diagonal elements.

## Crystal symmetry for GW calculation

**Space-group operations.** `lmf` finds the space-group operations (with spglib by default; [lmf: symmetry with spglib](./lmf#spglib-symmetry)) and,
with `--jobgw`, writes the point-group parts `symops(3,3,ig), ig=1,ngrp` and the translations `ag` into `__HAMindex` (`m_hamindexW_init` in `m_hamindex.f90`, called in `main_lmf.f90`).
For example, GaAs has 24 operations.

The GW programs read them through the module `m_hamindex` (`symops`, `ngrp`, `invgx`).
`mptauof_zmel` in `m_zmel.f90` calls `mptauof` (`m_symderive.f90`), which gives, for each operation, the atom to which each atom is mapped (`miat`) and the lattice translation (`tiat`).
See the comments at the beginning of `mptauof`.

The operations can be checked by `lmchk` in advance. For the ${\bf k}$ points, the inversion is added because of the time-reversal symmetry
($|\psi_{-{\bf k}}({\bf r})|^2=|\psi_{\bf k}({\bf r})|^2$ for a real local potential); this needs attention when there is an orbital moment (LDA+U, SOC).

The symmetry is controlled by `symgrp` at the top level of `ctrlg.<sname>.toml`: `"find"` (the default) uses the symmetry of the crystal structure,
and generators can be given to use a lower symmetry on purpose (for example, NiO is deformed along (111) in reality, but we may calculate it without the deformation).
**Antiferromagnetic symmetry**: write `af = ±k` on the pair of sites in `[[site]]`, and `symgrpaf = "find"`; the operations that map the up-spin site to the down-spin one
(with the spin inversion) are added. Use it with `pwmode = 11`.

## Spherical Harmonics and Real harmonics in ecalj

In ecalj, fortran codes are mainly based on the real harmonics $y_{lm}(\hat{\bf r})$,
instead of the usual sperical (complex) harmonics $Y_{lm}(\hat{\bf r})$.
The coefficients of eigenfunctions and so on are ordered as, e.g.
$(m=-2, m=-1, m=0, m=1,m=2)$ for $l=2$.
For example, `LMXA=4`, we have `(4+1)**2=25` harmonics,
ordered as
$y_{00},y_{-11},y_{01},y_{11},y_{-22},y_{-12},... ,y_{22},y_{-33},...y_{33},
y_{-44},...y_{44}$.
The real harmonics (times $r^l$) are made by `sylm` (`sylm.f90`) with the normalization `cy` of `sylmnc`; for $l=1$ the order is $y, z, x$, and for $l=2$ it is $xy, yz, 3z^2-r^2, xz, x^2-y^2$ (up to the normalization). $m\gt 0$ is the $\cos m\phi$ and $m\lt 0$ the $\sin |m|\phi$ type.

$y_{lm}(\hat{\bf r})$ is defined from $Y_{lm}(\hat{\bf r})$.
Here $\hat{\bf r}=(\theta, \phi)$.

$$\begin{aligned}
 y_{l0}(\hat{\bf r}) 
  & \equiv  Y_{l0}(\hat{\bf r}). \\
 y_{lm}(\hat{\bf r}) 
  & \equiv  \frac{1}{\sqrt{2} }
           [ (-1)^m Y_{lm}(\hat{\bf r}) + Y_{l-m}(\hat{\bf r}) ]. \\
 y_{l-m}(\hat{\bf r})
  & \equiv  \frac{1}{\sqrt{2}i}
           [ (-1)^m Y_{lm}(\hat{\bf r}) - Y_{l-m}(\hat{\bf r}) ].
\end{aligned}
\tag{I87}
$$
, where $m\gt 0$. Equivalently,
$$\begin{aligned}
 Y_{l0}(\hat{\bf r}) 
  & \equiv  y_{l0}(\hat{\bf r}). \\
 Y_{lm}(\hat{\bf r}) 
  & \equiv  \frac{(-1)^m}{\sqrt{2} }
           [ y_{lm}(\hat{\bf r}) + iy_{l-m}(\hat{\bf r}) ]. \\
 Y_{l-m}(\hat{\bf r})
  & \equiv  \frac{1}{\sqrt{2} }
           [ y_{lm}(\hat{\bf r}) - iy_{l-m}(\hat{\bf r}) ].
\end{aligned}
\tag{I88}
$$

The definition of $Y_{lm}(\hat{\bf r})$ is
$$\begin{aligned}
 &Y_{lm}(\theta, \phi)
=(-1)^m \left[ \frac{(2l+1)(l-m)!}{4 \pi (l+m)!} \right]^{\frac{1}{2} } P^m_l(\cos(\theta)) e^{i m \phi}, \\
 &P^m_l(x) = \frac{(1-x^2)^{m/2} }{2^l l!}\frac{d^{l+m} \ \ }{dx^{l+m} } (x^2-1)^l
\end{aligned}
\tag{I89}
$$

We take these definitions from
(1)A.R.Edmonds, Angular Momentum in quantum Mechanics,
Princeton University Press, 1960,
(2)M.E.Rose, Elementary Theory of angular Momentum,
John Wiley & Sons, INC. 1957. The definition of spherical harmonics is the same in these books.

## Rotation of eigenfunctions and MPB by the space-group operations {#spacegroup-rotation}

### space-group rotation of PMT eigenfunction

We sometimes have to map eigenfuncions by space group operations.

An space-group operation is specified by $(g,{\it \Delta})$, which contains
a $3 \times 3$ space-rotaiton matrix $g$ together with a translation vector
${\it \Delta}$. Recall that ${\it \Delta}= 0$ for the symmorphic cases.
This $(g,{\it \Delta})$ makes a rotation ${\bf r}\to {\bf r}'$ as
$$\begin{aligned}
 &{\bf r}'=g({\bf r})+ {\it \Delta},
\end{aligned}
\tag{I96}
$$
which is equilalent to
$$\begin{aligned}
 &{\bf r}=g^{-1}({\bf r}')+ {\it \Delta}^{-1},\\
 &{\it \Delta}^{-1}=-g^{-1}({\it \Delta})
\end{aligned}
\tag{I97}
$$
Further, we calculate which atomic site is mapped to which
atomic site for give $(g,{\it \Delta})$. It is given as
$$\begin{aligned}
 &g({\bf R})+{\it \Delta}={\bf R}'+{\it \Delta}{\bf T}_{\bf R}.
\end{aligned}
\tag{I98}
$$
In `ecalj`, the array `miat(ipos=1,npos)` gives ${\bf R}'$ in
the primitve cell.
`tiat(3,ipos)` gives the translation of multiples of primitive cell vectors
${\it \Delta}{\bf T}_{\bf R}$. Here ${\bf R}$ is the atom in the unit cell.

The rotation of an function $F \to F'$ by $(g,{\it \Delta})$
is defined so that $F'({\bf r}')=F({\bf r})$. Here $F'({\bf r}')$ is written as $g[F]({\bf r}')$.
Thus we have
$$\begin{aligned}
g[F]({\bf r})=F(g^{-1}({\bf r})+{\it \Delta}^{-1})
\end{aligned}
\tag{I99}
$$

Let us apply this to the bloch sum function
$A^{\bf k}_{ {\bf R}u L}({\bf r}) \equiv \sum_{\bf T} 
A_{ {\bf R}u L}({\bf r-R-T}) \exp(i {\bf kT})$ in Eq. (I3).
(Here we get the angular momentum index $L$. $u$ is for radial part of index.)
Then we have
$$
\begin{aligned}
g[A^{\bf k}_{ {\bf R}u L}]({\bf r})&=A^{\bf k}_{ {\bf R}u L}(g^{-1}({\bf r})+{\it \Delta}^{-1}) \\
&=\sum_{\bf T}A_{ {\bf R}u L}(g^{-1}({\bf r})+{\it \Delta}^{-1} -{\bf R}-{\bf T}) \exp(i {\bf kT})  \\
&=\sum_{\bf T}A_{ {\bf R}u L}\big(g^{-1}({\bf r}+g({\it \Delta}^{-1}) -g({\bf R})-g({\bf T}))\big) \\
&\qquad\times \exp(i {\bf kT})  \\
&=\sum_{\bf T}\sum_{L'}  A_{ {\bf R}u L'} ({\bf r}+g({\it \Delta}^{-1}) -g({\bf R})-g({\bf T})) \\
&\qquad\times D_{L'L}(g) \exp(i {\bf kT})  \\
&=\sum_{\bf T}\sum_{L'}  A_{ {\bf R}' u L'}({\bf r}-{\bf R}'-{\it \Delta}{\bf T}_{\bf R}-g({\bf T})) \\
&\qquad\times D_{L'L}(g)\exp(i {\bf kT})  \\
&=\sum_{ {\bf T}''} \sum_{L'}  A_{ {\bf R}' u L'}({\bf r}-{\bf R}'-{\bf T}'')\, D_{L'L}(g) \\
&\qquad\times \exp(i g({\bf k}){ {\bf T} }'') \\
&\qquad\times \exp(-i g({\bf k}){\it \Delta}{\bf T}_{\bf R})  \\
&= \sum_{L'}  A^{g({\bf k})}_{ {\bf R}' u L'}({\bf r})\, D_{L'L}(g) \exp(-i g({\bf k}){\it \Delta}{\bf T}_{\bf R}),
\end{aligned}
\tag{I100}
$$
where we use $A_{ {\bf R}' u L}({\bf r})=A_{ {\bf R}u L}({\bf r})$ since $R$ and $R'$
sites should be equilalent sites. We use cancellation ${\it \Delta}=-g({\it \Delta}^{-1})$.
$D_{L'L}(g)$ is the rotaion matrix of the real spherical harmonics for $g$.

### Code: `rotevec`

The routines are in `rotwave.f90`: `rotevec` rotates the coefficients of an eigenfunction in the PMT basis (MTO and APW parts), and `rotmatMTO` gives the rotation matrix of the MTO part.

With this Eq. (I100), coeffecients of eigenfunctions for the MTO part of
eigenfunctions are rotated (mapped) by the space group operation $g$.
It is coded in `subroutine rotevec`(eigenfunction rotation routine) as

    ...
     phase = [(exp(-img2pi*sum(qtarget*tiat(:,ibas,igg))),ibas=1,nbas)]
     do iorb=1,norbmto !orbital-blocks are specified by ibas, l, and k.
       ibas = ibastab(iorb)
       l   = ltab(iorb)
       k   = ktab(iorb)
       init1 = offl(iorb)+1
       iend1 = offl(iorb)+2*l+1
       init2 = offlrev(miat(ibas,igg),l,k)+1
       iend2 = offlrev(miat(ibas,igg),l,k)+2*l+1
       evecout(init2:iend2,:)= matmul(dlmm(-l:l,-l:l,l,igg),evec(init1:iend1,:))*phase(ibas)
     enddo
    ...

Here

> `ndimh = nlmto + napw`: the size of PMT Hamiltonian\
> `evec(1:nlmto,1:nband)` : the coeffecients of eigenfunctions on MTO\
> `igg`: index specifying the space-group operation $g$\
> `qtarget`: qtarget=$g({\bf k})$\
> `ibas`: site index within the primitive cell.\
> `iorb` is the index of MTO part of block ($2l+1$ elements for $m=-l,-l+1,... ,l$).\
> `l` is the angularmomentum index. l of L.\
> `k` is the radial index.\
> `miat`: site index ibas is mapped to `miat(ibas,igg)`\
> `tiat`: ${\it \Delta}{\bf T}_{\bf R}$\
> `offl`: offset for the `iorb` block\
> `offlrev`: offset for the (ibas,l,k) block\
> The Block `init1:iend1` is mapped to the block `init2:iend2`.

In order to call `rotevec`, we have to call `readhamindex()` so as to set up variables
`ltab,ktab,offl,offlrev,...` in `m_hamindex` module in advance.

Orbital index tables (`ibas_tbl`, `l_tbl`, `k_tbl`, ...) are read from a file `@MNLA_CPHI`
(written by `sugw.f90`) by `readmnla_cphi` (`readeigen.f90`).
The last column in `@MNLA_CPHI` contains the orbital-block index.
The orbital-blocks are specified by $l$ and $k$, its size is $2l+1$.
The number of $i$ is exactly the number of lines of `MNLA_CPHI`.\

In `subroutine rotevec`(eigenfunction rotation routine) we also have part for the APW.
Formula is
$$
\begin{aligned}
g[\exp(i ({\bf k}+{\bf G}){\bf r})]  
&=\exp\big(i ({\bf k}+{\bf G}) (g^{-1}({\bf r})+{\it \Delta}^{-1})\big) \\
&=\exp\big(i (g({\bf k})+g({\bf G})) ({\bf r}+g({\it \Delta}^{-1}))\big)  \\
&=\exp\big(i g({\bf k}+{\bf G}) {\bf r}\big) \exp\big(i ({\bf k}+{\bf G}) {\it \Delta}^{-1}\big)  \\
&= \exp\Big({i \big(\bar{g}({\bf k}) + g({\bf G})+ (g({\bf k})-\bar{g}({\bf k}))\big) {\bf r} }\Big) \\
&\qquad\times\exp\big(-i\, g({\bf k}+ {\bf G}){\it \Delta} \big).
\end{aligned}
\tag{I101}
$$

(Our definition of PW (IPW) has 'no normalization factor'.)
Here $g({\bf k})$ is pulled back to be $\bar{g}({\bf k})$ in the 1stBZ.
Their difference is $g({\bf k})-\bar{g}({\bf k})$.
We define `barG` as $\bar{G}=g({\bf G})+ (g({\bf k})-\bar{g}({\bf k}))$.
Thus it is coded as\
(here we use ${\bf q}$ instead of ${\bf k}$, sorry for confusion).

    ...
      igg: index for space group operation.
    nlmto: the size of MTO part.

      ikt  = getikt(q)       !index for q 
      ikt2 = getikt(qtarget) !index for bargq = \bar{g}(q) 
      do ig = 1,napw_in !number of APW 
        qpg  = q + matmul( qlat(:,:),igv2(:,ig,ikt))  ! q+G
        qpgr = matmul(symops(:,:,igg),qpg)            !g(q+G)
        nnn= nint(matmul(platt,qpgr-qtarget)) !integer sets for barG= g(G)+ g(q)-bargq
        ig2 = ngvecprev(nnn(1),nnn(2),nnn(3))   ! index for barG
        evecout(nlmto+ig2,:)= evec(nlmto+ig,:) * exp( -img2pi*sum(qpgr*shtvg(:,igg)) )
      enddo
    ...

### space-group rotation of MTIPW expansion of eigenfunction

In ecalj, we re-expand eigenfuncitons in the manner of MTIPW form, that is,
MTpart+IPWpart as in Eq. (I7) in the GW calculations.
In the form of MTIPW, we have to map eigenfunctions by space-group operations.
As in the manner of `rotevec`, we have `rotmto` and `rotipw` (`rotwave.f90`).
These are for rotation MTpart and IPW parts; they are used in `readeigen.f90` to give `__CPHI` and `__GEIG` at the ${\bf k}$ points that are not stored, from the stored equivalent ones.

### space-group rotation of Mixed Product Basis

We need all the matrix elements $\langle \Psi_{ {\bf q+k}n'}| \Psi_{ {\bf k}n} E^{\bf q}_\nu \rangle$.
for given inputs of ${\bf q},{\bf k},n,n'$. Here $E^{\bf q}_\nu$ is the Coulomb-diagonalized mixed product basis
[friedrich_efficient_2010].

This is a memo for the EIBZ scheme, which is not used now ([IBZ and EIBZ](#ibz-and-eibz-scheme)); `m_rotMPB2` and `x0kf_v4hz_symmetrize` are not in the present code, while `rotmto2` and `rotipw2` remain in `rotwave.f90`.
In the procedure of EIBZ (See Eq. (I95) around), we have to rotate MPB.
Module `module m_rotMPB2 (in m_rotMPB.F)` is for generating the rotation matrix `zrotm`.
With `zrotm`, we rotate MPB in `call X0kf_v4hz_symmetrize`, where
we multiply the conversion matrix `zcousq` for conversion from MPB to $\{E^{\bf q}_\nu\}$ as\
`  rcxq(:,:,iw,jpm) = matmul(zcousqc,matmul(rcxq_core,zcousq))` (see code in `x0kf_v4hz_symmetrize`).

`zrotm` is generated by calling `rotmto2,rotipw2`.
These are essentially similar with `rotmto,rotipw`, but the basis set in the MT region,
and the basis set for the interstitial region are different from the case of eigenfunctions.
Thus we supply `ibas_tbl` and so on as

          use m_pbindex,only: norbt, ibas_tbl,l_tbl,k_tbl,offset_tbl,offset_rev_tbl,
         &     max_ibas_tbl,max_l_tbl,max_k_tbl,max_offset_tbl

in `module m_rotMPB2`.

## SOC matrix for the MLO

The Wannier-function route (`hsocmat`) was removed. The SOC matrix is now used for the MLO (`job_mlo_soc`; [MLO §4](./mlo#mlo-soc)):

1. `lmf --writeham --socmatrix` (scalar-relativistic, `[ham] so = 0`) writes the Hamiltonian in the PMT basis (`__HamiltonianPMT`) and, separately,
   the SOC matrix $\langle F_{{\bf k}i}|H_{\rm SOC}|F_{{\bf k}j}\rangle$ (`hammhso(ndimh,ndimh,3)` made by `aughsoc`; the three are the $\uparrow\uparrow$, $\downarrow\downarrow$ and $\uparrow\downarrow$ blocks, per orbital without spinor doubling) into `__HamiltonianPMTsoc` (`m_bandcal.f90`).
   Here $F_{{\bf k}i}$ are the basis functions of the PMT method, and $\Psi_{{\bf k}n\sigma}=\sum_i F_{{\bf k}i}\alpha^{\bf k}(i,n,\sigma)$.
2. `mlo --socmatrix` transforms it to the MLO basis together with the Hamiltonian (`hammhsor` in `HamRsMLO`, read in `m_mlo_ham.f90`),
   and adds it as a perturbation to the MLO Hamiltonian (spinor dimension `nspinor = 2`).
3. The Fermi energy is taken from a full-SOC `lmf` run (`[ham] so = 1`, `efermi_soc`).

The SOC splitting of the MLO bands should be checked against the bands of the full calculation with SOC (`job_mlo_soc` draws both).

## Harris-Foulkner energy and Kohn-Sham energy

In LDA/GGA, on the way to self-consistency, input density and output
density are not the same (not self-consitent).
Thus we define two total energy for given input density $n^{\rm in}$,
the Harris-Foulkner energy $E_{\rm Harris}$ and the Hohenberg-Kohn energy $E_{\rm HoKohn}$
$$
\begin{aligned}
E_{\rm Harris}& =  E_{\rm k}^{\rm core} + E_{\rm B} \\
 &\quad - V_{\rm es}[n^{\rm Zc}+n^{\rm in},{\bf R}_a] \cdot n^{\rm in}
- V_{\rm xc}[n^{\rm c}+n^{\rm in}] \cdot n^{\rm in} \\
 &\quad + E_{\rm es}[n^{\rm Zc}+n^{\rm in},{\bf R}_a] + E_{\rm xc}[n^{\rm c}+n^{\rm in}],\\
E_{\rm HoKohn}& =  E_{\rm k}^{\rm core} + E_{\rm B} \\
 &\quad - V_{\rm es}[n^{\rm Zc}+n^{\rm in},{\bf R}_a] \cdot n^{\rm out}
- V_{\rm xc}[n^{\rm c}+n^{\rm in}] \cdot n^{\rm out} \\
 &\quad + E_{\rm es}[n^{\rm Zc}+n^{\rm out},{\bf R}_a] + E_{\rm xc}[n^{\rm c}+n^{\rm out}], \\
E_{\rm B} & =  \sum_p^{\rm occupied} \alpha_{p}^{i*} \langle F_i|H^{\rm in}|F_j \rangle \alpha_p^j,
\end{aligned}
\tag{I80}
$$
See Eqs.(B.1),(B.2)in Ref.[kotani_formulation_2015].
These are slightly wrong because Eq.(31) for $E_{\rm xc}$ is wrong; it is
not the functional of total charges (inlcuding nucleus) $n^{\rm Zc}+n^{\rm in}$,
but the total electron density $n^{\rm c}+n^{\rm in}$.
Based on the 3-component formalism in Ref.[kotani_formulation_2015],
$E_{\rm es}$ have explicit dependendence on atomic positions ${\bf R}_a$, but $E_{\rm xc}$ does not.
The dependence is via the mutipole transformation in Eq.(14).
Search `ham_ehf,ham_ehk` in `m_bndfp.f90`.
The $E_{\rm Harris}$ is given at `call m_mkehkf_etot1(sev, eharris)`, while
$E_{\rm HoKohn}$ is given by
`call m_mkehkf_etot2(ekinval, eksham)` (`mkehkf.f90`).
$E_{\rm Harris}$ and $E_{\rm HoKohn}$ are shown in `save.*` file.

The band energy $E_{\rm B}$ is stored in the module variable `sev` in `m_bndfp.f90`;
this is shown as `sev` in `save.*`. The kinetic energy
$E_{\rm k}^{\rm core}+E_{\rm B} 
- V_{\rm es}[n^{\rm Zc}+n^{\rm in},{\bf R}_a] \cdot n^{\rm out}
- V_{\rm xc}[n^{\rm c}+n^{\rm in}] \cdot n^{\rm out}$ contained in
the $E_{\rm HoKohn}$ is calculated as `ekinval` by `mkekin` (`m_bndfp_util.f90`).

In LDA/GGA calculations by `lmf`, `save.*` file contains a line per iteration (written in `m_lmfp.f90`; the energies are relative to the sum of the reference energies `eref` of the species over the atoms, `ham_seref` of `m_lmfinit.f90`).

    c ehf(eV)=-15730.0982239 ehk(eV)=-15730.0982222 sev(eV)=-15.7737393

shows $E_{\rm Harris}=$-15730.0982239 eV and $E_{\rm HoKohn}=$-15730.0982222 eV, as well as the valence
band energy $=-15.7737393$ eV.
In principle, $E_{\rm Harris}$ and $E_{\rm HoKohn}$ should be exactly the same when converged;
the difference is the numerical error. `c` at the begining of line means "converged".
${\tt h}$ means the 1st iteration from `atm.*` file (superposition of atomic density); the letters are `svflg='cxhiC'` in `m_lmfp.f90` (`x`: the maximum number of iterations reached, `i`: an intermediate iteration, `C`: converged relaxation).

## IBZ and EIBZ scheme {#ibz-and-eibz-scheme}

The IBZ part below is what the code uses: the self-energy is summed over ${\bf k}$ of the IBZ and the rotations (`kx` and `irot` in `m_sxcf_sc.f90`, with the table `irkip`). The EIBZ part is a memo; it is not used now (`eibzgen.f90` is there but not called, `[gw] EIBZmode` is read but not used, and `m_eibz`, `x0kf_v4hz_symmetrize`, `zsecsym` and `rotwvigg` are not in the present code). It may be used again.

We describe the EIBZ procedure given in III.G. in Ref.[friedrich_efficient_2010], for the polarization function and for the self-energy.
[IBZ]{.underline}\
Recall IBZ first. We have a set of space-group operation
$S_A=\{A_i |i=1,... N_A\}$ (except translation by reciprocal vectors).
$A_i$ is represented by space rotation and translation.
This $S_A$ can include time-reversal operation. ${\bf q}'=A_i({\bf q})$ belongs to BZ.
(Simple multiplication of rotation matrix to ${\bf q}$ may (or may not) need pulling back to (1st) BZ).
For a set of mesh points "$\{ {\bf q}\}$ in the BZ", $A_i$ gives an one-to-one mapping between them.
"$\{ {\bf q}\}$ in IBZ" is a sub set of 1stBZ.

Consider a set of functions dependent of ${\bf q}$ as $\{g_{ {\bf q} }({\bf r}, n)\}$.
We can calculate $X$ defined as a sum of $F[g_{\bf q}]$ as
$$\begin{aligned}
X=\sum_{ {\bf q}\in {\rm BZ} } F[g_{\bf q}]
= \sum_{A_i \in S_A } 
\sum_{ {\bf q}\in {\rm IBZ} } \frac{N_{\bf q} }{N_A} F[g_{A_i({\bf q})}] \end{aligned}
\tag{I91}
$$
Here IBZ is a set of ${\bf q}$ from which we can generate all the mesh points in the BZ.
$N_{\bf q}$ is the number of mesh points generated from ${\bf q}$ by $S_A$ (the number of the set ${\bf q}^*$).
For ${\bf q}=0$, $N_{\bf q}=1$. For general ${\bf q}$ (lowest symmetric points), $N_{\bf q}=N_A$.

If $A_i[g_{\bf q}]= g_{A_i({\bf q})}({\bf r},n)$ is satisfied, we have
$$
\begin{aligned}
X=\sum_{ {\bf q}\in {\rm BZ} } F[g_{\bf q}]
&= \sum_{A_i \in S_A } \sum_{ {\bf q}\in {\rm IBZ} } \frac{N_{\bf q} }{N_A} F[g_{A_i({\bf q})}] \\
&= \sum_{A_i \in S_A } \sum_{ {\bf q}\in {\rm IBZ} } \frac{N_{\bf q} }{N_A} F[A_i[g_{\bf q}]] .
\end{aligned}
\tag{I92}
$$
This shows how we evaluate $X$ only from $\{g_{ {\bf q} }({\bf r}, n)|{\bf q}\in {\rm IBZ}\}$.

(NOTE: eigenfunctions calculated by diagonalization satisfy
$A_i[g_{\bf q}]= g_{A_i({\bf q})}({\bf r},n)$ except phase factor. We consider cases that the phase factor is irrelevant.)

 \
[EIBZ]{.underline}\
We consider sub group of $S_A$ with keeping given ${\bf k}$.
This sub group is written as $S^{\bf k}_A=\{A^{\bf k}_i |i=1,2,...,N_A^{\bf k}\}$.
For this sub group, we apply the same logic of IBZ, resulting
$$\begin{aligned}
X=\sum_{ {\bf q}\in {\rm BZ} } F[g_{\bf q}] = 
\sum_{A^{\bf k}_i \in S_A^{\bf k} } \sum_{ {\bf q}\in {\rm EIBZ({\bf k})} }
\frac{N_{\bf q}^{\bf k} }{N_A^{\bf k} } F[A^{\bf k}_i[g_{\bf q}]] \end{aligned}
\tag{I93}
$$
Here EIBZ$({\bf k})$ is a set of ${\bf q}$ from which we can generate all the mesh points in the BZ.
$N^{\bf k}_{\bf q}$ is the number of mesh points
generated from $N_{\bf q}$ by $S^{\bf k}_A$. Note that ${\bf k}$ is fixed in Eq. (I93).
Since $S_A^{\bf k}\subset S_A$, IBZ $\subset$ EIBZ(${\bf k}$). That is, equivalent $q$ poinst in IBZ
might be diffrenciated in EIBZ(${\bf k}$).

See Eq.(50) in Ref.[friedrich_efficient_2010].
`call Seteibz -->call eibzgen` returns $S_A^{\bf k}$ stored in the modele `m_eibz`.

Note normazization is
$$\begin{aligned}
N=\sum_{\bf q}1 = 
\frac{1}{N_A^{\bf k} } \sum_{A^{\bf k}_i \in S_A^{\bf k} } 
\sum_{ {\bf q}\in {\rm EIBZ}({\bf k})}N_{\bf q}^{\bf k}=
\sum_{ {\bf q}\in {\rm EIBZ}({\bf k})}N_{\bf q}^{\bf k}\end{aligned}
\tag{I94}
$$

[For poralization]{.underline}\
We apply Eq. (I93) to $F^{\bf k}[g_{\bf q}]$, resulting
$$\begin{aligned}
X({\bf k})
=\frac{1}{N_A^{\bf k} } \sum_{A^{\bf k}_i \in S_A^{\bf k} } \left( 
\sum_{ {\bf q}\in {\rm EIBZ}({\bf k})}N_{\bf q}^{\bf k}F^{\bf k}[A^{\bf k}_i[g_{\bf q}]] \right).
\end{aligned}
\tag{I95}
$$
This means that we have to calculate weighted average only for
$g_{\bf q}\in {\rm EIBZ}({\bf k})$, and then finally symmetrized (symmetrization is in
`call x0kf_v4hz_symmetrize`
in `hx0fp0.sc.m.F,hx0fp0.m.F` or newly developing `hrcxq`).

We apply Eq. (I95) to Eq.(22) in Ref.[friedrich_efficient_2010]
resulting Eq.(51) in Ref.[friedrich_efficient_2010].
Then we need to rotate MPB for given space group rotation $A^{\bf k}_i$.
Note that $A^{\bf k}_i$ keeps ${\bf k}$ of MPB. This is in the next section.

Note that `call eibzgen` stores not only $S_A^{\bf q}$(`eibzsym`),
but also the weight $N_{\bf q}^{\bf k}$ (`nwgt`) in the module `m_eibz`.
(role of ${\bf q}$ and ${\bf k}$ may be turn around in `m_eibz`)

 \
[For self-energy]{.underline}\
In `hsfp0.sc.m.F`(not in `hsfp0.m.F`), we have symmetrization
as Eq.(52) in Ref.[friedrich_efficient_2010].
But used formula is a little different.

We use Eq. (I94) as it is.
As in the same manner of poralization function, we call `eibzgen`
at first. To obtain weighted sum of self-energy at ${\bf q}$,
we calculate weigted sum for ${\bf k}\in {\rm EIBZ({\bf q})}$.
(note the sum is controlled by `nrkip` which is the weights copied
from `nwgt` (MPI parallelization is used). When `nrkip=0`, we
skip corresponding do loop.)

Symmetrizer is `call zsecsym`, called from the main program of `hsfp0.sc.m.F`.
Its core part is

           zsect(ii1:ie1,ii2:ie2)= zsect(ii1:ie1,ii2:ie2) 
         &   + matmul( dconjg(transpose(rmatjj(1:ne1,1:ne1,iblk1))), 
         &     matmul(zsec(ii1:ie1,ii2:ie2,iqxx),
         &     rmatjj(1:ne2,1:ne2,iblk2)) )

. This is the self-energy rotated by the rotation matrix `rmatjj`,
which is generated by `call rotwvigg` in it.
(NOTE: we carefully treat degenerated bands to keep symmetry well.
It makes block decomposition(due to degeneracy) of self-energy matrix.
A block is `zsect(ii1:ie1,ii2:ie2)`, where band indeces
`ii1:ie1` are degenerated bands).

## Electron-phonon coupling in MPB (suspended) {#electron-phonon}

This work was suspended. The static $\epsilon^{-1}$ for phonons is written by `hgw` into `__W4PHONON.<iq>` when `w4pmode` is on, that is, when the offset-Γ points along $x,y,z$ were added by `qg4gw --job=201` (`ixyz` of `m_q0p.f90`; `m_llw.f90`).

To calculate derivative for electron-phonon (EP) coupling, we need
$$
\begin{aligned}
&\frac{\partial}{\partial {\bf r}_i} \left(\bar{J}_l(E,{\bf r}) Y_L(\hat{r})\right)\Big|_{ {\bf r}=0}=
\frac{1}{3} \sqrt{\frac{3}{4\pi} } \\
&\qquad (\text{only for } l=1;\ \text{zero otherwise}).
\end{aligned}
\tag{I39}
$$
For real spherical harmonics, cases are:
$i=y$ and $L=(1,-1)$; $i=z$ and $L=(1,0)$; $i=x$ and $L=(1,1)$.
Because we use MPB, it is necessary to evaluate the bare matrix element of the Coulomb interaction between nucleus and
$\langle B\rangle$ as
$$
\langle B|\frac{\partial v({\bf r}-{\bf R})}{\partial {\bf R}_i}\rangle
\tag{I40}
$$
