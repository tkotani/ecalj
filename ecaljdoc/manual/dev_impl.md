# 実装メモ（GW の実装の詳細） {#dev-impl}

::: warning 写したばかり（2026-11-15）
このページは、開発者向けの古いメモ `ecaljdetails.tex`（T. Kotani、最終版 2024-11。原本は ecalj の `ecaljdoc/ecaljdetails/`）から、
いまも役に立つ節を写したものです（英語のまま）。まだ見直していません。いまのコードに見つからない名前には朱の印を付けました。
式番号 (I1)〜 は写したときに付け直したもの、`[key]` は文献のキー（`ecaljdoc/ecaljdetails/ecaljrefs_detail.bib`）です。
理論の説明は [ecalj のバンド計算法](../theory/pmt)、[QSGW](../theory/qsgw) と、そこにある文献を見てください。
:::

[[toc]]

### Representation of eigenfunctions

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
Then we take into accounts the contributions of the cores to the exchange part defined in Eq. ((I15)) in the following. But not to the correlation part.
(caution: we now usally apply "core1 treatment" give in
Ref.[kotani07a] for all cores. Rarely use core2).

####  MTO part

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
The MTO basis is specified by smHankel functions which
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

####  APW part

The APW $\chi^{ {\bf k} }_{\bf G}({\bf r})$ are
given as a linear combination of atomic like orbitals
$A_{ {\bf R}{Lu} }({\bf r})\equiv 
\{ \phi_{ {\bf R}{l}u}(r) Y_L(\hat{\bf r}),
    \dot{\phi}_{ {\bf R}{l}u}(r) Y_L(\hat{\bf r}) \}$
within MTs, and just the usual plane waves
within the interstitial region. Here
$\phi_{ {\bf R}{l}u}(r)$ and $\dot{\phi}_{ {\bf R}{l}u}(r)$
denote two solutions of
the radial Schrödinger equations at an energy `enu`
for each $l$ (an usual choice of `enu` is the center of
gravity of occupied PDOS). $\dot{\phi}$ means energy derivatives
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

### Re-expansion of eigenfunctions: CPHI and GEIG

To perform the $GW$ calculation, we first
have to prepare all eigenfunctions (and eigenvalues) for
given setting of BZ mesh. Then the eigenfunctions are
represented as follows; we
re-expand ${\Psi_{ {\bf k}n} }({\bf r})$ in Eq. ((I1))
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
We explain how they can be represented in codes
in Sec. (ref sec:crystalstructure in the original).

We expand the eigenfunctions as the sum of the augmentation parts in MTs and the PW parts in the interstitial region
$$\begin{aligned}
\Psi_{ {\bf k}n} = \sum_{ {\bf R}u}\alpha_{ {\bf R}u}^{ {\bf k}n}\phi_{ {\bf R}u}^{\bf k}({\bf r}) + \sum_{\bf G}\beta_{\bf G}^{ {\bf k}n}P_{\bf G}^{\bf k}({\bf r})
\end{aligned}
\tag{I10}
$$
This is Eq.(17) in Ref.[kotani_quasiparticle_2014].
Here, Files CPHI contains the information of $\alpha_{ {\bf R}u}^{ {\bf k}n}$ and GEIG contains $\beta_{\bf G}^{ {\bf k}n}$.
We use subroutines `readcphi` and `readeig` to
read them; see `m_zeml.F`<span class="clfix">(not in the present SRC, 2026-11-15)</span> for example.
(In future, we may start from better representation
based on the 3 component formalism in Ref.[kotani2015pmt].)

We need $\Psi_{ {\bf k}n}$ for given ${\bf q}$ points
(in this text, we mix up ${\bf q}$ and ${\bf k}$\... Sorry.).
Then $P_{\bf G}^{\bf k}({\bf r})$ is just specified by
${\bf G}$, which is generated by **qg4gw**.
$\phi_{ {\bf R}u}^{\bf k}({\bf r})$ is specified by
radial functions. It is contained in PHIVC read in hbasfp0.m.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span>. Number of radial functions are
`ncore(ic)+nrad(ic)`
(dependes on $l,n,\sigma$, but not on $m$).
For simplicity, maximum of $l$ is fixed by LMXA in ctrl file.
It must be the same for all MTs.

In the GW calculation of ecalj, important matrix elements related
to the eigenfunction is only the matrix element as
$$\begin{aligned}
\langle E^{\bf q}_\mu \Psi_{ {\bf k}n} |\Psi_{ {\bf q+k}n'} \rangle,
\end{aligned}
\tag{I11}
$$
, where $E^{\bf q}_\mu$ is the MPB (an unitary transformation of MPB).
The information of eigenfunctions are used to calculate this matrix elements, which is read by `get_zmel`.

Coefficients of Eq. ((I7)) (here is MTO part only)
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
${\bf G}'$ is by $|{\bf k+G}'|\lt  {\tt  QpGcutHakel}$.

### Overview of GW calculation

In the $GW$ calculation, we need not only the basis set for
eigenfunctions, but also the basis set for expanding the product of eigenfunctions.
The basis is called the mixed product basis (MPB) $\{M^{\bf k}_I({\bf r}) \}$ first introduced in Ref.[kotani_all-electron_2002] by Kotani. The MPB consists of the product basis (PB) within MTs [aryasetiawan_product-basis_1994]
and the IPW in the interstitial region.
Since $\{M^{\bf k}_I({\bf r}) \}$ contains IPWs which are not orthogonal,
we define dual for $\{M^{\bf k}_I({\bf r}) \}$ as
$$\begin{aligned}
 & |\tilde{M}^{\bf k}_{I} \rangle \equiv \sum_{I'}
   |M^{\bf k}_{I'} \rangle (O^{\bf k})^{-1}_{I'I} \, , \\
 & O^{\bf k}_{I'I} = \langle M^{\bf k}_{I'} |  M^{\bf k}_I \rangle.
\end{aligned}
\tag{I13}
$$
From $v_{IJ}^{\bf k}= \langle M^{\bf k}_{I} |v|  M^{\bf k}_J \rangle$,
we calculate the eigenfunction for the generalized eigenvalue problem defined by
$\sum_J (v_{IJ}^{\bf k}- v^{\bf k}_\mu O^{ {\bf k} }_{IJ} ) w_{\mu J}^{\bf k}= 0$, where
$v_\mu({\bf k})$ are the eigenvalues of the Coulomb interaction matrix.
Then we have the Coulomb interaction represented by matrix elements as
$$\begin{aligned}
v({\bf k})=\sum_{\mu} | E^{ {\bf k} }_\mu \rangle {v_\mu({\bf k})} 
\langle E^{ {\bf k} }_\mu |,
\end{aligned}
\tag{I14}
$$
where we define a new MPB
$|E^{\bf k}_\mu({\bf r})\rangle=\sum_J |M^{\bf k}_J\rangle w^{\bf k}_{\mu J}$,
which is orthonormal and is diagonal to the Coulomb interaction $v({\bf k})$.
For the all-electron full-potential $GW$ approximation,
Eq. ((I14)) is introduced in Ref.[friedrich_efficient_2010].
This corresponds to the representation in the plane wave expansion
$v({\bf k}+{\bf G},{\bf k}+{\bf G}')=\frac{4 \pi \delta_{ {\bf G}{\bf G}'} }{|{\bf k}+{\bf G}|^2}$.
$\mu=1$ corresponds to the largest eigenvalue of $v_{\mu}$, and
$v_{\mu=1}$ is $\sim \frac{4 \pi e^2}{|{\bf k}|^2}$, which is related to
the divergent term discussed in Sec.([§](#offset-gamma)).

With the definition of
$\langle A| B\rangle =\int d^3r A^*({\bf r}) B({\bf r})$,
the exchange part of $\Sigma(\omega)$ is written as
$$\begin{aligned}
\Sigma^{\rm x}_{nm}({\bf q})=
\langle {\Psi_{ {\bf q}n} }|\Sigma_{\rm x} |{\Psi_{ {\bf q}m} }\rangle
=-\sum^{\rm BZ}_{ {\bf k} }  \sum^{\rm  occ}_{n'}
\langle {\Psi_{ {\bf q}n} }| \Psi_{ {\bf q-k}n'} E_\mu^{\bf k} \rangle
v_{\mu}({\bf k})
\langle E_\mu^{\bf k}\Psi_{ {\bf q-k}n'} | {\Psi_{ {\bf q}m} }\rangle.
\end{aligned}
\tag{I15}
$$

The screened Coulomb interaction $W(\omega)$ is calculated from
$$W = \epsilon^{-1} v = \left(1-v \Pi\right)^{-1} v,
\tag{I16}
$$
where the Lindhard polarization function $\Pi(\omega)$
is written as
$$\begin{aligned}
\Pi_{\mu \nu}({\bf q},\omega)
 &=
\sum^{\rm BZ}_{ {\bf k} } \sum^{\rm occ}_{n } \sum^{\rm unocc}_{n'}
\frac{
\langle E^{\bf q}_\mu {\Psi_{ {\bf k}n} }|\Psi_{ {\bf q+k}n'} \rangle
\langle \Psi_{ {\bf q+k}n'}| {\Psi_{ {\bf k}n} }E^{\bf q}_\nu \rangle
}{\omega-(\varepsilon_{ {\bf q+k} n'}-\varepsilon_{ {\bf k}n})+i \delta} \\
 &+ \sum^{\rm BZ}_{ {\bf k} } \sum^{\rm  unocc}_{n } \sum^{\rm occ}_{n'}
\frac{
\langle E^{\bf q}_\mu \Psi_{ {\bf k}n} |\Psi_{ {\bf q+k}n'} \rangle
\langle \Psi_{ {\bf q+k}n'}| \Psi_{ {\bf k}n} E^{\bf q}_\nu \rangle
}{-\omega-(\varepsilon_{ {\bf k}n}-\varepsilon_{ {\bf q+k} n'})+i \delta}.
\end{aligned}
\tag{I17}
$$
When time-reversal symmetry is assumed (if $\Psi_{ {\bf k}n}$ is the eigenfunction,
$\Psi^*_{ {\bf k}n}$ is also the eigenfuntion with the same eigenvalue
$\varepsilon_{ {\bf k}n}$),
$\Pi(\omega)$ is simplified to be
$$\begin{aligned}
\Pi_{\mu \nu}({\bf q},\omega)
 &=\sum^{\rm BZ}_{ {\bf k} }  \sum^{\rm  occ}_{n} \sum^{\rm  unocc}_{n'}
\langle E^{\bf q}_\mu \Psi_{ {\bf k}n} |\Psi_{ {\bf q+k}n'} \rangle
\langle \Psi_{ {\bf q+k}n'}| \Psi_{ {\bf k}n} E^{\bf q}_\nu \rangle  \\
 & \times
\left(\frac{1}{\omega-\varepsilon_{ {\bf q+k}n'}+\varepsilon_{ {\bf k}n}+i \delta}
-\frac{1}{\omega+\varepsilon_{ {\bf q+k}n'}-\varepsilon_{ {\bf k}n}-i \delta}\right). \end{aligned}
\tag{I18}
$$
This definition is the same as that in the Fetter-Walecka text book.
At $\omega=0$, this is negative definite matrix.
Im part at any $\omega$ is negative definite matrix. This is even function of $\omega$.

To evaluate Eq. ((I17)) or Eq. ((I18)),
we first accumulate the imaginary parts (anti-Hermitian part) of $\Pi_{\mu \nu}({\bf q},\omega)$
along bins of histograms on the real axis $\omega$
by the tetrahedron technique [rath_generalized_1975],
and then determine the real part via the Hilbert transformation.
The bins are dense near the Fermi energy and coarser at higher energy as described in
Ref.[kotani_quasiparticle_2007].
This procedure is not only more efficient but also safer than the methods
of calculating the real part directly.
We also use the extended irreducible zone (EIBZ)
symmetrization procedure described in Ref.[friedrich_efficient_2010].

The correlation part of the screened Coulomb interaction $W^c(\omega)=W(\omega)-v$,
which is calculated from $v$ and $\Pi(\omega)$, is given as
$$\begin{aligned}
W^{\rm c}({\bf k},\omega)=\sum_{\mu\nu} | E^{ {\bf k} }_\mu \rangle {W^{\rm c}_{\mu\nu}({\bf k},\omega)} 
\langle E^{ {\bf k} }_\mu |.
\end{aligned}
\tag{I19}
$$
With this $W^{\rm c}({\bf k},\omega)$, we have the correlation part of the self-energy as
$$\begin{aligned}

\Sigma^{\rm c}_{n,n'}({\bf q},\omega)= \sum_{ {\bf k},m} \int_{-\infty}^\infty  \!d \omega' \sum_{\mu,\nu} 
\frac{\langle{\Psi_{ {\bf q}n} }|{\Psi_{ {\bf q}-{\bf k}m} }E^{ {\bf k} }_\mu \rangle 
W^{\rm c}_{\mu\nu}({\bf k},\omega')\langle E^{ {\bf k} }_\nu {\Psi_{ {\bf q}-{\bf k}m} }
|{\Psi_{ {\bf q}n'} }\rangle e^{-i \delta \omega'} }
{\omega-\omega'-\epsilon_{ {\bf q}-{\bf k}m}\pm i \delta}.
\end{aligned}
\tag{I20}
$$
Here, we use $-i \delta$ for occupied states of
${ {\bf q}\!-\!{\bf k}m}$, and $+i \delta$ for unoccupied states.
In QSGW, we have to calculate the Hermitian part of $\Sigma_{nn'}({\bf q},\epsilon_{ {\bf q}n})$,
to obtain $V^{\rm xc}_{\bf q}$ using Eq. ((ref eq:vxceq in the original)).

There are two key points to handle the $GW$ procedure given above.
The first key point, given in Sec.([§](#offset-gamma)),
is the improved offset-$\Gamma$ method, which treats
the divergence of $W^{\rm c}({\bf k}\to 0,\omega)$ in Eq. ((I20)).
For this purpose, we define the non-divergent
effective interaction $\overline{W^{\rm c} }({\bf k}=0,\omega)$ instead of
${W^{\rm c} }({\bf k}=0,\omega)$.
Then we can take a simple discrete sum for both expressions of
Eqs.((I15)) and ((I20)).

The second point in Sec.([§](#sigma-interpolation)) is
how to perform
an interpolation to give $V^{\rm xc}_{\bf q}$ at any ${\bf q}$ in the whole BZ,
from $V^{\rm xc}_{\bf q}$ calculated only at limited numbers of ${\bf q}$ points.
This is required in the offset-$\Gamma$ method shown in Sec.([§](#offset-gamma)),
that is, we have to calculate eigenfunctions at some ${\bf q}$ points near ${\bf q}=0$.
For the interpolation, we expand the static
nonlocal potential $V^{\rm xc}$ in Eq. ((ref eq:vxceq in the original)) in highly localized
MTOs in real space. Thus such MTOs are used for two purposes:
one as the basis of the eigenfunctions; and two as the basis of
expanding $V^{\rm xc}$. The interpolation procedure of $V^{\rm xc}_{\bf k}({\bf r},{\bf r}')$
becomes stabler and simpler than the
complicated interpolation procedure in
Ref.[kotani_quasiparticle_2007].
This is because we now use highly localized MTOs.
In the planewave-based QSGW method by Hamann and Vanderbilt [hamann09],
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
give a very weak dependence to the final numerical results in practice as seen in Sec.(ref sec:numtest in the original),
although the results formally does not depend on the LDA/GGA exchange-correlation functions.

## ${\bf q}$ and ${\bf G}$ vector generation. qg4gw

*\... xxxxx under construction xxxxx\...\
*
To see the theory, read Section 3.2 in in [kotani_quasiparticle_2014]
${\bf k}$ points are the regular mesh points used for the integration
in the BZ. Regular mesh points are given as
$$\begin{aligned}
{\bf k}_{i_1,i2,i3}=(i1/N1)*{\bf P}_1+(i2/N2)*{\bf P}_2+(i3/N3)*{\bf P}_3
\end{aligned}
\tag{I21}
$$
, where ${\bf P}_i$ is the primitive vector $\frac{2\pi}{\rm alat} \times$
Quantities $f({\bf k})$ (periodic in BZ) can
be integrated just as a sum on the regular mesh points.
However, when $f({\bf k}=0)$ is divergent, we have careful treatments.

`qlat`(:,i1), where qlat=${\bf Q}_i$ is the reciprocal vectors.
See Sec.3.2 .

mkqg.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span> is the main part of **qg4gw** (fpgw/main/qg4gw.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span>).
iq0pin is the job control of mkqg. iq0pin=101 is only for backward compatibility.
The purpose is generates required q and G vectors.
Not only regular mesh points, or q along symmetry lines,
but also offset Gamma points.

iq0pin=1: normal mode. Generate q and G for regular mesh point and Q0P
points (offset Gamma points).

iq0pin: input to qg4gw


    ncindx, lcindx

    getkeyvalue
    phi

    radial mesh
    CLASS
    symgg
    core, radial functions

#### Make $G$ vectors: getgv2

To get ${\bf G}$ vectors, we use an algorithm in `fpgw/getgv2.F`, whose head
is

          subroutine getgv2(alat,plat,qlat,q, QpGcut,job,
         o                 ng, ngvec)
    !! == Set up a list of recip vectors within cutoff |Q+G| < QpGcut a.u. ==
    !! job==1 -> return ng (number of G ) and imx(as ngvec(1,1));mar2012takao add imx.
    !! job==2 -> return ng and ngvec
    !! True G is given as
    !!    G(1:3,1:ng) = 2*pi/alat * matmul(qlat * ngvec(1:3,1:ng))
    !! NOTE: we need some geometorial consideration for this routine.
    !!   Consiser ellipsoid. Takao need to give more detailed explanation...
    !! -------------------------------------------------------------

We can use this to get $\{ {\bf G}\}$ for given $q$.
The algorism of `getgv2` is a little complicated.
We first gives the upper and lower limits `n1max` $\le n_1 \le$ `n1min`,
where ${\bf G}= n_1 Q_1 + n_2 Q_2 + n_3 Q_3$, $n_2$ and $n_3$ as well.

*\... xxxxx under construction xxxxx\...\
*


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

#### Make $G$ vectors: shortn3

Find shortest vector in modulo of $\{Q_i\}$. That is, pull back ${\bf q}$ in the 1st BZ.
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
Thus only the elements of overlap matrix is
$$\begin{aligned}
 O^{\bf k}_{IJ}= 
  \int_\Omega \{M_I^{\bf k}({\bf r})\}^* {M}_J^{\bf k}({\bf r}) d^3r,
  ,
\end{aligned}
\tag{I24}
$$
because $\{ P^{\bf k}_{\bf G}({\bf r})\}$
are not orthogonal.
Thus it is convenient to define the dual of ${M}^{\bf k}_{I}({\bf r})$ as
$\tilde{M}^{\bf k}_{I}({\bf r})$ in the manner of Eq. ((I13)).

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

---- NOTE: ----
Now we use Coulomb matrix orthogonalized basis $\{E^{\bf k}_\mu\}$ as
discussed in Eq. ((I14)).

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
In addition, $p$ can taks indexes to specify core functions:
we combine core and valence functions (these
are stored in `PHIVC`, which read in `hbasfp0.m.F`).
Here is a part of copy to read `PHIVC` in hvccfp0.m.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span>

>           ifphi  = iopen('PHIVC', 0,-1,0)! augmentation wave and core
>           read(ifphi) nbas, nradmx, ncoremx
>           allocate(  ncindx(ncoremx,nbas),
>          &           lcindx(ncoremx,nbas),
>          &           nrad(nbas), nindx_r(1:nradmx,1:nbas),
>          &           lindx_r(1:nradmx,1:nbas),
>          &        aa(nbas),bb(nbas),zz(nbas), rr(nrx,nbas), nrofi(nbas) ,
>          &        phitoto(nrx,0:nl-1,nn,nbas,nsp),
>          &        phitotr(nrx,0:nl-1,nn,nbas,nsp),
>          &        nc_max(0:nl-1,nbas),ncore(nbas) )
>           read(ifphi) nrad(1:nbas)
>           read(ifphi) nindx_r(1:nradmx,1:nbas),lindx_r(1:nradmx,1:nbas)
>           nc_max=0
>           do ibas=1,nbas
>             write(6,*)' --- read PHIVC of ibas=',ibas
>             ic = ibas
>             read(ifphi) ncore(ic), ncoremx                            !core
>             read(ifphi) ncindx(1:ncoremx,ibas),lcindx(1:ncoremx,ibas) !core
>             read(ifphi) icx,zz(ic),nrofi(ic),aa(ic),bb(ic)
>             if(ic/=icx) then
>               write(6,*) 'ic icx=',ic,icx
>               call rx( 'hbasfp0: ic/=icx')
>             endif
>             read(ifphi) rr(1:nrofi(ic),ic)
>             do isp = 1, nsp
>               write(6,*)'---  isp nrad ncore(ic)=',isp, nrad(ic),ncore(ic)
>               do icore = 1, ncore(ic)
>                 l =  lcindx(icore,ic)
>                 n =  ncindx(icore,ic)
>                 read(ifphi) phitoto(1:nrofi(ic),l,n, ic,isp)!core orthogonal
>                 phitotr(1:nrofi(ic),l,n, ic,isp)=   !we set core raw= core orthgonal
>          &      phitoto(1:nrofi(ic),l,n, ic,isp)       
>                 if(n>nc_max(l,ic)) nc_max(l,ic)=n
>               enddo
>               do irad = 1, nrad(ic)
>                 l = lindx_r (irad,ic)
>                 n = nindx_r (irad,ic) + nc_max(l,ic)
>                 read(ifphi) phitoto(1:nrofi(ic),l,n, ic,isp) !valence orthogonal
>                 read(ifphi) phitotr(1:nrofi(ic),l,n, ic,isp) !valence raw
>               enddo
>             enddo
>           enddo
>
> \[note: The orthonomalized radial functions $u_{apl\sigma}(r)$
> are stored in `phitoto`; we also have the un-orthonormalized ones in `phitotr`.\]

Note that the *true* radial function is
$\phi_{apl\sigma}(r) = u_{apl\sigma}(r)/r$.
Normalization is
$1=\int_0^{R_a} \{ u_{apl\sigma}(r) \}^2 dr = 
   \int_0^{R_a} \{ \phi_{apl\sigma}(r) \}^2 r^2dr$.
The function $u_{apl\sigma}(r)$ is stored in `phitot`.

When producing the product functions,
we use spin-averaged function `phiav` given as
$$\begin{aligned}
   u_{apl}(r) = \frac{1}{N_{\rm spin} }\sum_{\sigma}
   u_{apl\sigma}(r).
\end{aligned}
\tag{I27}
$$
(See `subroutine basnfp_v2`).
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
  O_{\nu_1\nu_2} = \int_0^{R_a} 
  \tilde{B}_{al\nu_1}(r) \tilde{B}_{al\nu_2}(r) r^2 dr
  = \int_0^{R_a}  \phi_{a p_1 l_1}(r) \phi_{a p'_1 l'_1}(r)
                  \phi_{a p_2 l_2}(r) \phi_{a p'_2 l'_2}(r) r^2 dr
\end{aligned}
\tag{I30}
$$
and solve the eigenvalue problem of the overlap matrix,
$Oz_{\nu}=\epsilon_{\nu} z_{\nu}$, by `call rs(..)`.
(See `basnfp_v2`.)

After neglecting eigenvectors
$z_{\nu}$ with eigenvalues $\epsilon_{\nu} \lt  {\rm tolerance} \sim 10^{-4}$ (given in GWinput,
we finally have the optimal product functions as the linear combinations
of the product functions as
$$\begin{aligned}
   b_{al\nu}(r) = \frac{1}{\sqrt{\epsilon_{\nu} }}
   \sum_{\nu'}\tilde{b}_{al\nu'}(r)z_{\nu'\nu},
\end{aligned}
\tag{I31}
$$
which are stored in `rprodx` and written into BASFP\*
and used in the successive Coulomb matrix routine hvccfp0.m.f.
Of course, *true* product function is
$B_{al\nu}(r) = b_{al\nu}(r)/r$.

We check the normalization of the optimal product function
in standard output
(See `lbasC` and `lbas` when you did **gw_lmf**):

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
   <\phi \phi B> = \int_0^{R_a}
   \phi_{ap_1 l_1}(r) \phi_{ap_2 l_2}(r) B_{al\nu}(r) r^2 dr
   = \int_0^{R_a} \frac{1}{r}
     u_{ap_1l_1}(r) u_{ap_2l_2}(r) b_{al\nu}(r) dr,
\end{aligned}
\tag{I32}
$$
which are stored into PPBRD\*.
At `call rdpp(...` in hx0fp0.m.f hxfp0.m.f,
we allocate and read `ppbrd`.

In addition, we read the \"rotated Clebsh-Gordon coefficient\"
$C(L,L_1,L_2,g)$ `cgr(lm,lm1,lm2,ng)`
, where $g$ is the index for space group
coefficient ( rotated by point group symmetries).

## The Coulomb matrix (hvccfp0.m.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span>)

We have to calculate
$v_{IJ}^{\bf k}= \langle M^{\bf k}_{I} |v|  M^{\bf k}_J \rangle$, which
appears right after Eq. ((I13)).
Our `hvccfp0.m.F` can handle the
case $\frac{\exp(- |\kappa| |{\bf r}_1-{\bf r}_2|)}{|{\bf r}_1-{\bf r}_2|}$.
The default is the bare Coulomb interaction, that is, $\kappa=0$.
The energy variable $E$ is given by $E=\kappa^2$, where
imaginary part of $\kappa$ is defined to be positive for negative $E$.
This $E=$`eee=screenfac()` in hvccfp0.m.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span>.
For example, we can set $|\kappa|=0.1$ (a.u.) by a line
\"`TFscreen 0.1`\" in GWinput if necessary. Look for `TFscreen` in `swithes.F`<span class="clfix">(not in the present SRC, 2026-11-15)</span>.

The MPB is made of IPWcou ($|{\bf q}+{\bf G}|$`<QpGcut_cou`) and PB,
that is, $\{M^{\bf k}_J\}=\{P_{\bf G}^{\bf k}({\bf r}), B^{\bf k}_{ {\bf R}{liL}\mu}({\bf r})\}$.
For given ${\bf k}$ (not explicitly shown in cases for simplicity),
IPWcou is specified just by the ${\bf G}$ vector.
PB can be specified just by the radial functions for each $l$
PB is generated in the manner of Sec.([§](#mixed-product-basis)) .
(because we have neither $m$ nor spin dependence).
For $GW$ calculation, we need the Coulomb matrix elements
$\langle B|v|B \rangle,\langle B|v|P_{\bf G}\rangle,$ and $\langle B|v|P \rangle$.
`hvccfp0.m.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span>-->vcoulq_4` handles these calculations.

### Spherical Bessel and related functions

We use a notation such that $X_L({\bf r})=X_l(r) Y_L(\hat{r})$; note that
their radial part is dependent only on $l$.
The ordinary definition of $L$-dependent spherical
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
See note for `genjh` in `mkjp.F`<span class="clfix">(not in the present SRC, 2026-11-15)</span>.
At $E=0$, this is reduced to be
$$\begin{aligned}
 &\bar{J}_l(E=0,{\bf r})=r^l/(2l+1)!!,  \\
 &\bar{H}_l(E=0,{\bf r})=r^{-l-1} (2l-1)!!.
\end{aligned}
\tag{I36}
$$

(We use very small negetive $E=-10^{-5}$ as default instead of $E=0$
for Coulomb matrix generator `hvccfp0.m.F`
to avoid numerical troubles (see default `screenfac` in switch.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span>)).

The source codes to define bessl is a little confusing because
of some convensions are mixed up\...(note at the beginning of besslr.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span>).

Main one is `call bessel(ex2,lx,phi(0:lx),psi(0:lx))` in mkjp.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span>.
This is defined in `bessl(ex2,lmax,phi,psi)` in `besslr.F`
(it calls `besslr` with `loka=F`).
`lx` is the upper limit of $l$.
For given `ex2`$=E\times x^2$,
this return spherical Bessel functions $j_l(\kappa x)$
in the following manner.
$$\begin{aligned}
\left.
\begin{array}{ll}
\bar{J}_l(E,{\bf r})=j_l(i|\kappa|r)/(i|\kappa|)^l={\tt phi(l)}r^l\\
\bar{H}_l(E,{\bf r})=h_l(i|\kappa|r){i (i|\kappa|)}^{l+1}={\tt psi(l)}/r^{l+1} 
\end{array}
\right\}& \ {\rm for}\ E<0, \ \kappa=i|\kappa|,
|\kappa|=\sqrt{|E|}  \\
\left.
\begin{array}{ll}
\bar{J}_l(E,{\bf r})=j_l(\kappa r)/(\kappa )^l = {\tt phi(l)}r^l  \\
\bar{N}_l(\kappa r) = {\tt psi(l)}\frac{1}{(\kappa)^{l+1} }={\tt psi(l)}/r^{l+1}
\end{array} 
\right\}& \ {\rm for}\ E>0, \ \kappa=\sqrt{E}
\end{aligned}
\tag{I37}
$$
for $0\le l \le{\tt lmax}$. These definition gives `phi`$=1/(2l+1)!!$
and `psi`$=(2l-1)!!$ at $E \to 0$.
Here $n_l(r)$ is the spherical Neumann functions.
That is, `psi` is for Hankel function
for negative $E$, and for Neumann function for positive $E$.
See `hansmr.F`, for example. `psi` is not used so often.
We have `lm7K/subs/besslr.F`, it is similar.

For convenience, we sometimes use
$$\begin{aligned}
 &R^{+l}=\bar{J}_l(E,r)(2l+1)!!={\tt psi(l)}r^l(2l+1)!!,   \\
 &R^{-l}=\bar{H}_l(E,r)/(2l-1)!!={\tt psi(l)}r^{-l-1}(2l-1)!!.
\end{aligned}
\tag{I38}
$$
Here note $R^{+k} \rightarrow r^k$ for $E\to 0$,and
$R^{-k} \rightarrow r^{-k-1}$ for $E\to 0$.
These $rR^{+k},rR^{-k}$ are `rkpr,rkmr` generated in subroutine `genjh`.\
xxxx We may need to simplify our treatment of bessel functions in
future\... xxxx

##### Electron-phonon coupling in MPB

To calculate derivative for electron-phonon (EP) coupling, we need
$$\begin{aligned}
 &\frac{\partial}{\partial {\bf r}_i} \left(\bar{J}_l(E,{\bf r})
                    Y_L(\hat{r})\right)\big|_{ {\bf r}=0}=
\frac{1}{3} \sqrt{\frac{3}{4\pi} }, \ \ \ \ {\rm only \ for} \ L=1,\ {\rm zero otherwize}.
\end{aligned}
\tag{I39}
$$
For real spherical harmonics, cases are:
$i=y$ and $L=(1,-1)$; $i=z$ and $L=(1,0)$; $i=x$ and $L=(1,1)$.
Because we use MPB, it is necessary to evaluate the
the bare matrix element of the Coulomb interaction between nucleus and
$\langle B\rangle$ as
$$\begin{aligned}
 &\langle B|\frac{\partial v({\bf r}-{\bf R})}{\partial {\bf R}_i}\rangle
\end{aligned}
\tag{I40}
$$

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
Eq. ((I42)) is nothing but the solution of Helmhorz
differential equation; it reduces to the usual Poisson equation at
$E=0$. For $E\lt 0$, we have Thomas-Fermi type function;
the numerator of Eq. ((I42)) is $\exp(i \kappa |{\bf r}-{\bf r}'|)=\exp(- |\kappa| |{\bf r}-{\bf r}'|)$.

### Used formulas {#used-formulas}

$$\begin{aligned}
\frac{1}{|{\bf r}-{\bf r}'|}
=4 \pi \sum_K \frac{r^k_<}{r^{k+1}_>}\frac{1}{2k+1}
Y^*_K(\hat{ {\bf r} })Y_K(\hat{ {\bf r} }) \end{aligned}
\tag{I43}
$$
See Appendix A in Ref.[svane_evaluation_1986].
This is generalized to be
$$\begin{aligned}
\frac{e^{-|\kappa| |{\bf r}-{\bf r}'|} }{|{\bf r}-{\bf r}'|}
=4 \pi \sum_K R^{+k}(r_<) R^{-k}(r_>)\frac{1}{2k+1}
Y^*_K(\hat{ {\bf r} })Y_K(\hat{ {\bf r} }) 
=4\pi \sum_L 
\bar{J}_L(E,r_<) \bar{H}_L(E,r_>), %Y^*_L(\hat{\bfr})Y_L(\hat{\bfr})
\end{aligned}
\tag{I44}
$$
The definition of $\bar{J},R^{+l}$ are in Eq. ((I35)) and Eq. ((I38)).
Here $E\lt 0$ and $\kappa= i|\kappa|=\sqrt{-E}$.

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

$$\begin{aligned}
\langle P^{\bf k}_{\bf G}|P^{\bf k}_{ {\bf G}'} \rangle
= \Omega \delta_{ {\bf G},{\bf G}'} -  
\sum_{a,L} \exp( i ({\bf G}-{\bf G}') {\bf R_a}) \times Y_L(\widehat{ {\bf k+G} }') 
Y_L( {\bf \widehat{ k+G} })  \\
\times \int_0^{R_a} j_l(|{\bf k+G}|r)j_l(|{\bf k+G}'|r) 4 \pi^2 r^2 dr,
\end{aligned}
\tag{I46}
$$

### Hankel function and Structure constant

We can expand
$v({\bf r},{\bf r}')=\frac{e^{-|\kappa||{\bf r}-{\bf r}'|} }{|{\bf r}-{\bf r}'|}$
in the one-center expansion Eq. ((I44)).

For the Hankel function in Eq. ((I44)), we use
the off-center expansion theorem of the
Hankel function Eq. ((I47)), that is,
a Hankel whose center is at ${\bf X}\equiv{\bf R}+{\bf T}$, $H_L({\bf r}-{\bf X})$,
can be expanded in the Bessel functions whose center is at ${\bf X}'$;
$$\begin{aligned}
\bar{H}_L(E,{\bf r}-{\bf X})=\frac{1}{4\pi} \sum_{L'}\bar{J}_{L'}(E,{\bf r}-{\bf X}') S_{ {\bf X}'{L'},{\bf X}{L} }\ ,
\end{aligned}
\tag{I47}
$$
where the Hankel function for negative energy $E$.
Here $E$-dependence of $S_{ {\bf X}'{L},{\bf X}{L'} }$ is not explicitly shown.
Note the difference between $\bar{J}_L$ and $J_L$ ($\bar{H}_L$, as well).

Thus, for $({\bf R}',{\bf T}')\neq({\bf R},{\bf T})$,
we have two-center expansion;
$$\begin{aligned}
 &
\frac{e^{-|\kappa||{\bf r}+{\bf R}+{\bf T}-({\bf r}'+{\bf R}'+{\bf T}')|} }
{|{\bf r}+{\bf R}+{\bf T}-({\bf r}'+{\bf R}'+{\bf T}')|}
%= -4 \pi G(\bfr+\bfR+\bfT-(\bfr'+\bfR'+\bfT'),E)
% \\
%&= -4 \pi \sum_L (-i\kappa) \bar{J}_L(E,\hat{\bfr}) 
%\bar{H}_L(E,\bfr'+\bfR'+\bfT'-\bfR-\bfT) 
= \sum_L \sum_{L'}
{\bar{J}_L(E,{\bf r})}
%{(2l+1)!!} 
S_{ {\bf R}+{\bf T}{L},{\bf R}'+{\bf T}'{L'} } {\bar{J}_{L'}(E,{\bf r}')} 
%{(2l'+1)!!}
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
Usually we use bare Coulomb at $E=0$.\
[ No 4$\pi$ factor in the definition of $S^{ {\bf k} }_{ {\bf R}{L},{\bf R}'{L'} }$]{.underline}\
See the top of `strxq` defined in strxq.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span> (called in hvccfp0.m.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span>).
This routine is for the one-center expansion of
usual Bloch summed Hankels. (not for smooth Hankels).
This result is finally converted to be the
Bloch sum of the structure constant
$S^{ {\bf k} }_{ {\bf R}{L},{\bf R}'{L'} }$ used in Eq. ((I50)).

> ====================\
> NOTE in strxq.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span>; it says
>
>     -------------------------------------
>     Cr Expansion Theorem: H_{RL}(r) = H_L(r-R)
>     Cr H_{RL}(E,r) = J_{R'L'}(E,r) * S_{R'L',RL}
>     Cr S_R'L',RL = 4 pi Sum_l" C_{LL'L"} (-1)^l (-E)^(l+l'-l")/2 H_L"(E,R-R')
>     -----
>     CAUTION!: We use R to denote MT position in the primitive cell; 
>     thus this R is R+T in our notation.
>     --------------------------------

Thus it seeminly $4 \pi$ factor is missing in (I48).
However, righ after the subroutine `strxq` is called, we have added
$4\pi$ factor as `strx(1:nlx1,ibas1,1:nlx2,ibas2) = fpi*s` in
`hvcc0.m.F`<span class="clfix">(not in the present SRC, 2026-11-15)</span>. Thus Eq. (I48)
don't include $4\pi$ factors.

### $\langle B|v|B \rangle$ part {#b-v-b-part}

Let us start from $\langle B|v|B \rangle$ part.
For this calculation, we need structure constant, and a few types of
radial integrals.
With the Bloch-summed structure constant `strx`
$=S^{\bf k}_{RL,R'L'}$, we have
$$\begin{aligned}
v^{\bf k}({\bf r},{\bf r}')
=\sum_{ {\bf T}'} v({\bf r},{\bf r}'+{\bf T}') e^{i {\bf k}{\bf T}'}
=e^2 \sum_{L,L'} 
{\bar{J}_L(E,{\bf r}-{\bf R})}
%{(2l+1)!!} 
S^{ {\bf k} }_{ {\bf R}{L},{\bf R}'{L'} }
{\bar{J}_{L'}(E,{\bf r}'-{\bf R}')},
%{(2l+1)!!}.
\end{aligned}
\tag{I50}
$$
for ${\bf R}\ne {\bf R}'$.
Note that we use $e^2=1$ (a.u.) in hvccfp0.m.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span>.
(See note at the beginning of mkjp.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span>. And note the normalization check
at the end of hvccfp0.m.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span>; $\langle \exp(i {\bf q}{\bf r}) |v |\exp(i {\bf q}{\bf r})\rangle=4\pi \Omega/|{\bf q}|^2$, where $\Omega$ is the cell volume.)

`strx calculated by "call strxq@L806:hvccfp0.m.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span>"`
means $S^{\bf k}_{ {\bf R}{L},{\bf R}'L'}$. Note that $4 \pi$ is mulipled right
after `call strxq` to obtain .
This `strx` is used in `call vcoulq_4`.
`nlx1` means `(l+1)**2` for $R$
`(=ibas1)`, `nlx2` as well.

Except the contribution for $({\bf R},{\bf T})=({\bf R}',{\bf T}')$,
we can evaluate
$\langle B^{\bf k}_{ {\bf R}{L}\mu}({\bf r})|v^{\bf k}({\bf r},{\bf r}')|B^{\bf k}_{ {\bf R}'{L'}\mu'}({\bf r}') \rangle$, from the\
`rojb` integrals $\rho^l(B_{ {\bf R}{l}\mu})$ as
$$\begin{aligned}
 &\rho^l(B_{ {\bf R}{l}\mu})=
%\frac{1}{(2l+1)!!} 
\int_0^{R} r\bar{J}_l(E,r) \ r B_{ {\bf R}{l}\mu}(r)  dr
= \frac{1}{(2l+1)!!}\int_0^{R} {\tt rkpr(r) \ rprodx(r)} \ dr
\end{aligned}
\tag{I51}
$$
Here $B_{ {\bf R}{l}\mu}(r)$ is the radial part of $B_{ {\bf R}{L}\mu}({\bf r})$.
This `rojb` integrals are calculated in
the subrouitne `mkjb_4` in `mkjp.F`<span class="clfix">(not in the present SRC, 2026-11-15)</span>.
We use radial functions
`rprodox` $=r B_{ {\bf R}{l}\mu}(r)$,
`rkpr`$=r\bar{J}_l(E,r)(2l+1)!!$,
and `rkmr`$=r\bar{H}_l(E,r)/(2l-1)!!$.
(Here `rkpr` and `rkmr` are propotional to
$r^l$ and $r^{-l-1}$ for $E =0$.)

The contribution from $({\bf R},{\bf T})=({\bf R}',{\bf T}')$
should be added. This is ${\bf k}$-independent, and
given by the\
`sgbb` integral, which is also calculated in `mkjb_4`.
$$\begin{aligned}
 &\sigma^l(B_{ {\bf R}{l}\mu},B_{ {\bf R}{l}\nu})=
{4\pi}
\int_0^{R}\int_0^{R} (r_<)\bar{J}_l(E,r_<) (r_>)\bar{H}_l(E,r_>)
r B_{ {\bf R}{l}\mu}(r) r' B_{ {\bf R}{l}\nu}(r') dr dr'  \\
 &=\frac{4\pi}{2l+1}
\int_0^{R}\int_0^{R} {\tt rkpr(r_<) \ rkmr(r_>) \ rprodx(n1,r) \ rprodx(n2, r')} 
\ dr dr'.
\end{aligned}
\tag{I52}
$$

With the integrals `rojb` and `sgbb`, we can calculate
$\langle B|v|B \rangle$ in `vcoulq_4` as follows
(`nbloch` means the total number of PB);

>           do ibl1= 1, nbloch
>             ibas1= ibasbl(ibl1)
>             n1   = nbl (ibl1)
>             l1   = lbl (ibl1)
>             m1   = mbl (ibl1)
>             lm1  = lmbl(ibl1)
>             do ibl2= 1, ibl1
>               ibas2= ibasbl(ibl2)
>               n2   = nbl (ibl2)
>               l2   = lbl (ibl2)
>               m2   = mbl (ibl2)
>               lm2  = lmbl(ibl2)
>               vcoul(ibl1,ibl2) =
>          &     rojb(n1, l1, ibas1) *strx(lm1,ibas1,lm2,ibas2)
>          &    *rojb(n2, l2, ibas2)
>               if(ibas1==ibas2 .and. lm1==lm2) then
>                 vcoul(ibl1,ibl2) = vcoul(ibl1,ibl2) + sgbb(n1,n2,l1, ibas1)
>                 ! sigma-type contribution. onsite coulomb
>               endif
>             enddo
>           enddo

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
$l_{\rm Pmax} = 2 \times l_{\rm max}$`=2*LMXA`,
where $l_{\rm max}$ denotes the maximum angular momentum
for the expansion of eigenfunctions within MT (maximum
$l$ cutoff for $\alpha^{ {{\bf k} }n}_{ {\bf R}u}$ in Eq. ((I1))).
In the defalut setting, we use $l_{\rm Pmax}=8$
since we use $l_{\rm max}=4$.

The matrix elements $\langle P_1 |v| P_2 \rangle$
can be calculated as
$$\begin{aligned}
\langle P_1(phiphi) |v| P_2(phiphi) \rangle=
\sum_{\bf G_{1'} G_{1''} G_{2'} G_{2''} } 
\langle P_1(phiphi) | {P_{1'} } \rangle  
\langle {P}_{1'}| {P}_{1''} \rangle^{-1} 
\langle {P}_{1''} |v| {P_{2''} } \rangle  
\langle {P}_{2''}| {P}_{2'} \rangle^{-1} 
\langle {P_{2'} } | P_2(phiphi) \rangle,  
\end{aligned}
\tag{I54}
$$
where $1\equiv ({\bf k},{\bf G_1})$ and so on.
Here $P_1(phiphi)$ indidates
that an IPW made from a product
of IPWs.The matrix elements
$\langle {P}_{2''}| {P}_{2'} \rangle^{-1}$
(stored into PPOVLG,PPOVLI)
and $\langle {P_{2'} } | P_2(phiphi) \rangle$
(stored into PPOVLGG)
are given at `rdata4gw`.

### Overlap matrix of PPOVL\* files

PPOVL\* files contains the
overlap matrix of IPWs
$\langle {P}^{\bf k}_{ {\bf G} }|{P}^{\bf k}_{ {\bf G'} }\rangle$
.
We have two types of $\langle {P}^{\bf k}_{ {\bf G} }|{P}^{\bf k}_{ {\bf G'} }\rangle$.
One is for the Coulomb matrix (PPOVLG,PPOVLI).
The other is for generating
$\langle G({\rm eigenfun.}) G({\rm eigenfun.})| G({\rm cou}) \rangle$.
(a product of IPWs of eigenfunctions can be
expanded by IPWs for Coulomb matrix).

The overlap matrix elements
$\langle {P}^{\bf k}_{ {\bf G} }|{P}^{\bf k}_{ {\bf G'} }\rangle$
are generated in `rdata4gw`.
These are read and allocated in the module
`m_read_ppovl (rppovl.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span>)` when we call
`getppx2`. We have \"call getppx2\" in the
subroutine melpln2t in ppbafp.fal.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span>.
The melpln2t is for generating the matrix element of
`<IPW psi |psi>`.

    (1)PPOVLG + PPOVLI:
       For q in qibze(1:3,1:nqnumt) (=IBZ + Q0P points), 
       number of IPWcou =ngc can be dependent of q.
       We have <k+G|k+G'>= ppovl(ngc,ngc) 
       PPOVLG: G vectors as ngvecc(1:ngc).
       PPOVLI: ppovl^-1(ngc,ngc). Inverse of PPOVL0
       (PPOVL0 is unused now. It is divided into PPOVLG and PPOVLI).

       In principle the matrix element itself is k-independent,
       (just the difference of G vectors due to periodicity).
       But, for convenience, we generate them separately for each k.

    (2)PPOVLGG: 
      This is used for <Gphi Gphi|Gc>.
      ppovl(nggg,ngcgp) for nvggg,nvgcgp
      Range of G for nvggg is |Gc+Gp+Gp|< |Gcou|+ |Gphi|+ |Gphi| 
      (triangle inequality.)
      This is only for k=0 (Thus we remove k-dependece).
      ngcgp 
        QpGcutggg = (2d0+1d-2)*QpGcut_psi+QpGcut_cou+ 2d0*pi/alat*dQpG 
        QpGcutgcgp= (1d0+1d-2)*QpGcut_psi+QpGcut_cou+ 2d0* 2d0*pi/alat*dQQ
      dQpG, dQQ is to enlarge range related to Q0P points.

### $\langle {P}^{\bf k}_{\bf G}|v| {P}^{\bf k}_{ {\bf G}'} \rangle$

To evaluate Eq. ((I54)), we need to know
its main part $\langle {P}^{\bf k}_{\bf G}|v| {P}^{\bf k}_{ {\bf G}'} \rangle$. It is written as
$$\begin{aligned}
\langle {P}^{\bf k}_{\bf G}|v| {P}^{\bf k}_{ {\bf G}'} \rangle
\approx
\langle \bar{P}^{\bf k}_{\bf G}|v| \bar{P}^{\bf k}_{ {\bf G}'} \rangle
& = \langle \exp(i({\bf k}+{\bf G}){\bf r}) |v| \exp(i({\bf k}+{\bf G}'){\bf r}) \rangle
-\sum_{ {\bf R}{L} } 
\langle P^{ {\bf k}+{\bf G} }_{ {\bf R}{L} }|v| \exp(i({\bf k}+{\bf G}'){\bf r}) \rangle  \\
& - \sum_{ {\bf R}'{L'} } 
\langle \exp(i({\bf k}+{\bf G}'){\bf r}) |v| P^{ {\bf k}+{\bf G}'}_{ {\bf R}'{L'} }\rangle
+\sum_{ {\bf R}{L} } \sum_{ {\bf R}'{L'} } 
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

The first term in the right-hand side of Eq. ((I55)) is
$$\begin{aligned}
\langle \exp(i({\bf k}+{\bf G}){\bf r}) |v| \exp(i({\bf k}+{\bf G}'){\bf r}) \rangle
=\frac{4 \pi \Omega}{|{\bf k}+{\bf G}|^2+|E|} \delta_{ {\bf G}{\bf G}'},
\end{aligned}
\tag{I56}
$$
because we simply use $\exp(i({\bf k}+{\bf G}){\bf r})$
(no prefactor for normalization) for IPW.
Here $E$`=eee` is negative (or (almost) zero).
This is coded by a line

      if(ig1==ig2) vcoul(ipl1,ipl2) = fpivol/(absqg2(ig1) -eee)
    subr:vroulq_4@L281:mkjp.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span>.

##### The second and third term

In the second term,
we can replace $v$ with $\frac{4 \pi \Omega}{|{\bf k}+{\bf G}'|^2+|E|}$
since $v$ is diagonal for $|\exp(i({\bf k}+{\bf G}'){\bf r}) \rangle$, the third term as well.
Without $v$, we have
$$\begin{aligned}
 &\langle \exp(i({\bf k}+{\bf G}){\bf r})| P^{ {\bf k}+{\bf G}'}_{ {\bf R}'{L'} }\rangle
=\sum_{ {\bf R}{L} } \langle 
P^{ {\bf k}+{\bf G} }_{ {\bf R}{L} }|
P^{ {\bf k}+{\bf G}'}_{ {\bf R}'{L'} }\rangle  \\
 &=\sum_{ {\bf R}{L} } 
({\tt pjyl\_}({\bf k}+{\bf G},L) \exp(i({\bf k}+{\bf G}){\bf R}))^* 
\times R^{JJ}(|{\bf k}+{\bf G}|,|{\bf k}+{\bf G}'|,l)  \\
 &\times {\tt pjyl\_}({\bf k}+{\bf G}',L) \exp(i({\bf k}+{\bf G}'){\bf R}), \end{aligned}
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
Search `pjyl_` in `mkjp.F`<span class="clfix">(not in the present SRC, 2026-11-15)</span>.
The Bessel functions appear here in the expansion of PW; see Eq. ((I45)).
$R^{JJ}(|{\bf k}+{\bf G}|,|{\bf k}+{\bf G}'|,l)$
is given as
$$\begin{aligned}
R^{JJ}(|{\bf k}+{\bf G}|,|{\bf k}+{\bf G}'|,l)
= \int_0^R  r^2 \bar{J}_l(|{ {\bf k}+{\bf G} }|r)
\bar{J}_l(|{ {\bf k}+{\bf G}'}|r) dr ,
\end{aligned}
\tag{I60}
$$
which can be calculated by the wronskian (`wronskj`) by the formula
$$\begin{aligned}
R^{JJ}(\kappa_{\rm A},\kappa_{\rm B},l)=
\int_0^R \!\!r^2 \bar{J}_l(\kappa_{\rm A} r)
\bar{J}_l(\kappa_{\rm B} r) dr
=R^2\frac{
\bar{J}_l(\kappa_{\rm A} r) \frac{d \bar{J}_l(\kappa_{\rm B} r)}{dr}
-\frac{d \bar{J}_l(\kappa_{\rm A} r)}{dr}\bar{J}_l(\kappa_{\rm B} r)
}{\kappa_{\rm A}^2-\kappa_{\rm B}^2}\Big|_{r=R}= {\tt -fjj}
\end{aligned}
\tag{I61}
$$
(`-fjj` is used in `mkjp.F`<span class="clfix">(not in the present SRC, 2026-11-15)</span>.
In codes, the contributions to the second and third terms of Eq. ((I55)) due to ${\bf R}{L}$ components are given as (simplified for illustration)

      fouvp_ig1_ig2 = fpi/(absqg2(ig1)-eee)                  &
             * dconjg(pjyl_(lm2,ig1)*phase(ig1,ibas2))      &  
             * (-fjj(l)) * pjyl_(lm2,ig2)*phase(ig2,ibas2)
      fouvp_ig2_ig1 = fpi/(absqg2(ig2)-eee)                  & 
             * dconjg(pjyl_(lm2,ig2)*phase(ig2,ibas2))      &
             * (-fjj(l)) * pjyl_(lm2,ig1)*phase(ig1,ibas2)

Look for the keyword `fourvp` in `mkjp.F`<span class="clfix">(not in the present SRC, 2026-11-15)</span>.
Correspondences are

    fpi --> 4 pi
    absqg2(ig1) --> |q+G1|**2
    -eee        --> |E|
    lm2         --> L
    ibas2       --> R

##### The forth term

The last term of Eq. ((I55)) can be calculated essentially
the same manner with $\langle B| v| B\rangle$, where we use the
Bessel function instead of $B_{ {\bf R}{l}\mu}(r)$ appeared in Sec.([§](#b-v-b-part)).
Then we define integrals
`rojp` and `sgpp` defined as (in `fpgw/gwsrc/mkjp.F`);
$$\begin{aligned}
{\tt rojp}({\bf R}L)& = {\tt pjyl\_} \ 
\exp(i ({\bf q}+{\bf G}) {\bf R})  \rho^l(\bar{J}_l) \\
{\tt sgpp}({\bf R}L,{\bf G},{\bf G}')& = {\tt pjyl\_^*(ig1)} \exp(-i ({\bf q}+{\bf G}) {\bf R})
{\tt pjyl\_(ig2)} \exp(i ({\bf q}+{\bf G}') {\bf R}) 
 \\
 &  \times {\tt  radsig},  \\
 & \ {\rm where\ }
{\tt radsig}=\sigma^l(\bar{J}_l(|{\bf q}+{\bf G}|r),\bar{J}_l(|{\bf q}+{\bf G}'|r)).
\end{aligned}
\tag{I62}
$$
Thus, `rojp` and `sgpp` made of coefficients for explansion
and radial integral.
Search `sgpp_ig1_ig2`<span class="clfix">(not in the present SRC, 2026-11-15)</span> in `mkjp.F`<span class="clfix">(not in the present SRC, 2026-11-15)</span>.

### $\langle P|B \rangle$ part

$$\begin{aligned}
 &\langle {P}^{\bf k}_{ {\bf G} } |v| B^{\bf k}_{ {\bf R}{L}\mu}({\bf r}) \rangle
=\langle  \exp(i({\bf k}+{\bf G}){\bf r}) 
|v|B^{\bf k}_{ {\bf R}{L}\mu}({\bf r}) \rangle
-\sum_{ {\bf R}'{L'} } 
\langle P^{ {\bf k}+{\bf G} }_{ {\bf R}{L} }|v|
B^{\bf k}_{ {\bf R}{L}\mu}({\bf r})  \rangle  \\
 &={\tt fouvb}({\bf G},{\bf R}L \mu)
-\sum_{ {\bf R}'L'} {\tt rojp * strx * rojb}
\end{aligned}
\tag{I63}
$$
The first term is stored in `fouvb(ngc,nxx,nlxx,nbas)` allocated at
`L824:hvccfp0.m.F`.
Search `fouvb` in `mkjp.F`<span class="clfix">(not in the present SRC, 2026-11-15)</span>.
Since $v$ is diagonal to PWs,
we can evaluate this in the similar manner of Eq. ((I57)).
The second term can be evaluated from
`rojb` and `sgpb` in the same manner of last section.

*\... xxxxx under construction xxxxx\...\
*(we will detail a little more\...)

## Offset-$\Gamma$ method; W(k=0) averaged in the $\Gamma$ cell. {#offset-gamma}

The offset-$\Gamma$ method, originally
invented for Ref.[kotani_all-electron_2002] by Kotani
(it is described in Ref.[kotani_quasiparticle_2007]),
was a key to perform accurate $GW$
calculation in our papers.
It is for the integration of ${\bf k}$
in Eqs.((I15)) and ((I20)), where
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
$$\begin{aligned}
{\bf k}(i_1,i_2,i_3) & =  2 \pi (\frac{i_1}{N_1} {\bf b}_1 
+ \frac{i_2}{N_2} {\bf b}_2 + \frac{i_3}{N_3} {\bf b}_3),
%\\
%\sum_{\bf k}^{\rm BZ} & \approx  \frac{1}{N_1N_2N_3} \sum_{i_1,i_2,i_3},

\end{aligned}
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
We usually take the $\alpha \to 0$ limit, or a sufficiently small $\alpha$ instead.
Let us apply Eq. ((I65)) to $F_L({\bf k})$.
Then we can evaluate the left-hand side of Eq. ((I65))
exactly (the exact values are zero except for $L=(0,0)$).
On the other hand, the first and
third terms on the right-hand side of Eq. ((I65))
can be evaluated numerically.
In addition, we know that $f_{L'}$ for $F_L({\bf k})$ is unity for $L'=L$,
and zero otherwise.
Thus we can determine $w_L$ in Eq. ((I65))
so that Eq. ((I65)) is exactly satisfied for $F_L({\bf k})$ for any $L$.

Let us apply Eq. ((I65)) to $f({\bf k})=G({\bf q}-{\bf k}) \times W({\bf k})$.
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
neglected in the approximation of Eq. ((I67)).

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
in combination with the approximation of Eq. ((I67))
by calculations while changing the size of the ${\bf k}$-mesh.\

> **\[NOTE 2016-03-18:\]**
> We now use only $w_L$ for $L=(0,0)$.
> T.Kotani found strange behavior for $W(\omega)$ for the Wannier
> functions (for La2CuO4) when we use all $w_L$.
> Then I observed non-monotonic behavior of its real part; recall that
> we should see monotonic behavior as long as causality is satisfied.
> I found it is originated from non zero $w_L$ for $L\ne(0,0)$.
> Thus we now use the formula \"$w_L=0$ for $L\ne(0,0)$\" in
> Eq. ((I65)).In the case of metal at $\omega=0$, we have very large
> ${\bf L}$ (virtually infinite); thus no $1/|{\bf k}|^2$ behavior.\

The remaining problem is how to calculate
the matrix ${\bf L}$ in Eq. ((I67)); there are two possible ways.
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
from $\hat{ {\bf k} }_{0i}^{\rm T}{\bf L}(\omega)\hat{ {\bf k}^{0i} }$ to $a_g(\omega)$
is not numerically degenerated. The length $|{\bf k}^{0i}|$
can be chosen to be sufficiently enough, but avoiding numerical error
as the average of $W({\bf k})$ in the $\Gamma$ cell.
The improved offset-$\Gamma$ method shown here can be applicable even to metal cases, as long as
$\hat{ {\bf k} }_{0i}^{\rm T}{\bf L}(\omega)\hat{ {\bf k}_{0i} }$ contains the contribution
of intraband transition.

> **\[NOTE 2016-Nov:\]**
> T.Kotani think evaluation of $L_ij(\omega)$ can be more effectively
> performed as for the convergence on number of ${\bf k}$ points.
> At least, without local field correction, we only need to calculate head part.
> For the head part, we may use large number of ${\bf k}$ points.

## Fourier transformation of non-local quantity

We have \"`call bloch`\" in `lm7K/fp/bndfp.F`.
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
`bloch` called in `fp/bndfp.F`),
we first calculate $f_{ {\bf R}{\bf R}'}({\bf k})$ on ${\bf k}$ of regular mesh points. Then we
need to obtain its real-space representation
$f({\bf R}{\bf T},{\bf R}'{\bf T}')$. Because of the periodicity,
we have ambiguity for the choice of possible
$|{\bf T}-{\bf T}'|$. If we introduce $\bar{\bf T}={\bf T}-{\bf T}'$,
Eq. ((I72)) is written as
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

To get a list of $\bar{ {\bf T} }$,
we need to collect them satisfying
$|{\bf R}-{\bf R}'+\bar{\bf T}|\lt \eta_{\rm FTmax}$.
$\eta_{\rm FTmax}$ should be automatically chosen.
However, in the \"`bloch`\" subroutine, this it too primitive yet(aug2015);
we need to specify possible upper limit of
\"range of allowed pairs\" $({\bf R}\bar{ {\bf T} },{\bf R}'0)$ by hand.
(`RSRNGE` in `ctrl` file). This should be fixed in future.
In the current version `iaxs (=``sham%iv_a_oiaxs``)` contains such pair table.
It is generated by `call hft2rs` in `call seneinterp` in `bndfp.F`, I think.
We will have to replace \"`bloch`\" with better version. Pair table must be generated
in a simple manner( with the technique of
`getgv2` (`getgv2` is given in `fpgw/gwsrc/getgv2.F` and `lm7K/subs/pairs.F`.

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

At the end of step (IV) in Sec.(ref sec:qsgwth in the original), we obtain the matrix elements
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
Thus $\chi^{ {\bf k} }_a$ denotes the APWs or MTOs in Eq. ((I1));
$z_{na}$ (${\bf k}$ is omitted for simplicity)
denotes the coefficients of the eigenfunctions at ${\bf k}$, that is,
$z^{ {{\bf k} }n}_{ {\bf R}L{j} }$ and $z^{ {\bf k}n}_{ {\bf G} }$ in Eq. ((I1)) together.
This $z_{an}$ is identified as a conversion matrix that connects
eigenfunctions (band index $n$) and the APW and MTO bases (basis index $a$).

To obtain real-space representation of ${\it \Delta} V^{\rm xc}$, we need a representation
expanded in the basis that consist of the Bloch-summed localized orbitals,
which are periodic for ${\bf k}$ in the BZ. However, this is not the case for the APWs in
Eq. ((I73)). To overcome this problem, we use an approximation
in which we only take the matrix elements
related to MTOs, that is, the elements
$\langle \chi^{ {\bf k} }_a| {\it \Delta} V^{\rm xc}_{\bf k}| \chi^{ {\bf k} }_b \rangle$
where $a$ and $b$ specify MTOs.
This means that the part of ${\it \Delta} V^{\rm xc}$ related to APWs is projected onto the basis of MTOs.
This approximation can be reasonable as long as the main part of ${\it \Delta} V^{\rm xc}$ can be
well expanded in MTOs,
although we need numerical tests to confirm the accuracy as shown in Sec.(ref sec:numtest in the original).
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
eigenvalues are less than $E_{\rm MAX}^\Sigma$.
Then the high energy parts of the matrix elements are assumed
to be diagonal, where their values are given by a simple average
of calculated diagonal elements.

## Magnon project {#magnon-project}

Generalized Lindhard polarization function (or Kernel) is defined as
(c.f.Eq. ((I17)));
$$\begin{aligned}
 &-K^{\alpha \beta}=({\bf r}_1,{\bf r}_2;{\bf r}_3,{\bf r}_4;\omega) =
-iG_\alpha(1,3)G_\beta(4,2)=  \\
 &=
\sum^{\rm BZ}_{ {\bf k} }  \sum^{\rm occ}_{n } 
\sum^{\rm BZ}_{ {\bf k}'} \sum^{\rm unocc}_{n'}
\frac{
\Psi_{ {\bf k}n \beta}^*({\bf r}_2) \Psi_{ {\bf k}n \beta}({ {\bf r} }_4) 
\Psi_{ {\bf k}'n' \alpha}({\bf r}_1) \Psi^*_{ {\bf k'}n'\alpha}({\bf r}_3)
}{\omega-(\varepsilon_{ {\bf k}' n'\alpha}-\varepsilon_{ {\bf k}n \beta})+i \delta} \\
 &
+ \sum^{\rm BZ}_{ {\bf k} } \sum^{\rm  unocc}_{n } 
  \sum^{\rm BZ}_{ {\bf k}'} \sum^{\rm occ}_{n'}
\frac{
\Psi^*_{ {\bf k}n \beta}({\bf r}_2)\Psi_{ {\bf k}n \beta}({\bf r}_4) 
\Psi_{ {\bf k}'n' \alpha}({\bf r}_1) \Psi^*_{ {\bf k}' n' \alpha}({\bf r}_3)
}
{-\omega-(\varepsilon_{ {\bf k}n \beta}-\varepsilon_{ {{\bf k}'} n' \alpha})+i \delta}.
\end{aligned}
\tag{I74}
$$
**(--- NEED diagram FIGURE HERE for $K_{\alpha\beta}$ ---)**\
Here $t=t_1=t_2$ and $t'=t_3=t_4$.
To understand this, note
the correspondence between real-time and $\omega$ space;
$$\begin{aligned}
\frac{1}{\omega-\varepsilon +i \delta} &\leftrightarrow &
i\theta(t-t') \exp (-i \varepsilon (t-t')) \\
\frac{1}{-\omega-\varepsilon +i \delta} &\leftrightarrow &
i\theta(t'-t) \exp (-i \varepsilon (t'-t)).
\end{aligned}
\tag{I75}
$$
It is easy to make the product $G_\alpha(1,3)G_\beta(4,2)$ in real space
and real time represantaiton.
This Eq. ((I76)) is a general time-ordered linear response function for a non-interacting system. Note that this is reduced to be Eq. ((I17)) for ${\bf r}_1={\bf r}_2$
and ${\bf r}_3={\bf r}_4$.
This is the same as $K^{\alpha \beta}$
in Eq.(15) in Ref.[sasioglu_wannier-function_2010] by Sasioglu.
By the Fourier transformation of Eq. ((I76)),
we have
$$\begin{aligned}
-K^{\alpha \beta}({\bf q},\omega)& =  
\sum^{\rm BZ}_{ {\bf k} }  \sum^{\rm occ}_{n}\sum^{\rm unocc}_{n'}  
\frac{
\Psi_{ {\bf k}n \beta}^*({\bf r}_2) \Psi_{ {\bf k}n \beta}({ {\bf r} }_4) 
\Psi_{ {{\bf k}+{\bf q} }n' \alpha}({\bf r}_1) \Psi^*_{ {{\bf k}+{\bf q} }n'\alpha}({\bf r}_3)
}{\omega-(\varepsilon_{ {\bf k}+{\bf q}n'\alpha}-\varepsilon_{ {\bf k}n \beta})+i \delta} \\
& +  \sum^{\rm BZ}_{ {\bf k} } \sum^{\rm  unocc}_{n} \sum^{\rm occ}_{n'}
\frac{
\Psi^*_{ {\bf k}n \beta}({\bf r}_2)\Psi_{ {\bf k}n \beta}({\bf r}_4) 
\Psi_{ {\bf k}+{\bf q}n' \alpha}({\bf r}_1) \Psi^*_{ {\bf k}+{\bf q}n' \alpha}({\bf r}_3)
}
{-\omega-(\varepsilon_{ {\bf k}n \beta}-\varepsilon_{ {{\bf k}+{\bf q} } n' \alpha})+i \delta}.
\end{aligned}
\tag{I76}
$$
When we expand eigenfunctions by atom-centered localized functions as
$$\begin{aligned}
\Psi_{ {\bf k}n}^{\bf k}({\bf r}) = \sum_{ {\bf R}i} a^\alpha_{ {\bf R}i} w_{ {\bf R}i \alpha}^{\bf k}({\bf r}),
{\rm \bf (correct? \ Notation \ check)}
\end{aligned}
\tag{I77}
$$
in the restricted Model Hilbert space,
$K^{\alpha \beta}$ is represented as $K^{\alpha
\beta}_{ {\bf R}ij,{\bf R}'kl}({\bf q},\omega)$, where ${\bf R}i,{\bf R}j,{\bf R}'
k,{\bf R}'l$ are orbital indexes in the unit cell.
As we have $W=\langle {\bf R}i {\bf R}j|W|{\bf R}' k {\bf R}' l\rangle$ in the mRPA method
(or cRPA), we can calculate $\chi_{+-}={K}/(1-WK)$ in the model space.
**(need detailed equation with indexes!)**\
[**Spectrum funciton and Hilbert transformation**]{.underline}\
Imaginary part of $K^{\alpha \beta}({\bf q},\omega)$ is obtained just by
replacement $1/(\omega- \epsilon+i\delta) \rightarrow \pi \delta(\omega -\epsilon)$.
That is, it is given as
$$\begin{aligned}
{\rm Im}[K]= \sum_{\bf k}\sum_n \sum_{n'} M({\bf k},n,\beta;{\bf k}+{\bf q},n',\alpha)
\delta(\omega-(\varepsilon_{ {\bf k}+{\bf q}n'\alpha}-\varepsilon_{ {\bf k}n \beta})),
\end{aligned}
\tag{I78}
$$
where we define the matrix element
$M=\Psi_{ {\bf k}n \beta}^*({\bf r}_2) \Psi_{ {\bf k}n \beta}({ {\bf r} }_4) 
\Psi_{ {{\bf k}+{\bf q} }n' \alpha}({\bf r}_1) \Psi^*_{ {{\bf k}+{\bf q} }n'\alpha}({\bf r}_3)$.
By the replacement $\delta(\omega-\varepsilon) \rightarrow 1/(\omega-\varepsilon
+i \delta)$ by the Hilbert transformation, we can recover real part from the Im\[$K$\].
$\omega$ can be positive or negative; positive is for spin excitation $\beta
\to \alpha$, negative is for $\alpha \to \beta$.

In the tetrahedron method of ecalj, we can calculate Eq. ((I78)) as follows;
$$\begin{aligned}
{\rm Im}[K]\,([\omega_i,\omega_{i+1}])= 
\overline{\sum_{\bf k} } \sum_n \sum_{n'} M({\bf k},n,\beta;{\bf k}+{\bf q},n',\alpha)
W_{\rm tet}(i,{\bf k},n,{\bf k}+{\bf q},n'),
\end{aligned}
\tag{I79}
$$
where $\overline{\sum_{\bf k} }$ is for discritized ${\bf k}$ points set by ecalj.
${\rm Im}[K]\,([\omega_i,\omega_{i+1}])$ means the imaginary part (weight)
in the interval $[\omega_i,\omega_{i+1}]$.\
I think the t2g-eg separated motion is one of the interesting theme.
(need to write a little more\...)\
**Test of magnon calculation (transversal spin fluctuation)**\
Look into `ecalj/MATERIALS/Fe_magnon` and `Ni_magnon`<span class="clfix">(not in the present SRC, 2026-11-15)</span>. Run a script.
We have results of spin fluctuations.

### Data structure of the tetrahedron method in ecalj {#tetrahedron-data-structure}

At first, you can consider only the `mtet=F` case.
To make it safer, check `mtet=F` (if it is True, error exit by 'call rx').
In `m_tetw.F`<span class="clfix">(not in the present SRC, 2026-11-15)</span> we have `gettetwt` called from `hx0fp0.m.F` as

              call gettetwt(q,iq,is,isf,nwgt(:,iq),frhis,nwhis,npm,
         i     qbas,ginv, ef, nqibz, nband,ekxx1,ekxx2, nctot,ecore,
         i     nqbz,qbz,nqbzw,qbzw,  ntetf,idtetf,ib1bz,
         i     nbmx,ebmx,mtet,eibzmode) !nov2016

. All arguments are input.
This returns tetrahedron weight $W_{\rm tet}(i,{\bf k},n,{\bf k}+{\bf q},n')$
in the common data area of `m_tetw`<span class="clfix">(not in the present SRC, 2026-11-15)</span>.\

-   `q`: ${\bf q}$ vector and index of ${\bf q}$ vector

-   `iq`: unused for `mtet=F`. It may be safer to set -9999 for safe.

-   `is,isf`: $\alpha$ and $\beta$

-   Eigenvalues are supplied by

              do kx = 1, nqbz
                call readeval(qbz(:,kx),   is,  ekxx1(1:nband, kx) ) 
                call readeval(q+qbz(:,kx), isf, ekxx2(1:nband, kx) )
              enddo

    at `hx0fp0.m.F`

-   `nwgt`: This is for EIBZ mode.
    We skip microtetrahedrons when all of `nwgt(kx0:kx3)`
    are zero (here ,kx1,kx2,kx3@ are four corners of ${\bf k}$ vectors).
    Thus we calculate tetrahedron weigh only for ${\bf k}$ points for `nwgw(k)=1`.
    In EIBZ mode, we will symmetrize finally obtained correct $K^{\alpha\beta}$.
    (takao: I like to check this code again\...)\

-   `frhis(1:1+nwhis),nwhis`: Histogram bins $[\omega_i,\omega_{i+1}]$
    are @. We set `frhis(1)=0d0`.
    See `m_freq.F`. We calculate values on $K^{\alpha,\beta}(\omega)$ where
    $\omega$ is specified by `freq_r` (this is center of the bins.
    See `freq_r` file.)

-   If @, we calculate negative frequency part. For negative $\omega$,
    Histogram bins $[\omega_i,\omega_{i+1}]$
    are @.
    (Opposite spin flip excitation to positive energy.)

-   We can set `ncore=0` (genallcf_v3@).

-   For BZ data (m_read_bzdata@), you can set them
    by read_bzdata@.

-   (3,3,1:ngrp)@ is the rotation part of space group operation.
    SYMOPS file is created by `echo 0|lmfgw`, and read by qg4gw.
    Then symops information is written into HAMindex file.

Output data sets are
real(8),allocatable :: whw(:)

          integer,allocatable:: ihw(:,:,:),nhw(:,:,:),jhw(:,:,:),ibjb(:,:,:,:)
          integer:: nbnbx,nhwtot
          integer,allocatable :: n1b(:,:,:),n2b(:,:,:),nbnb(:,:)

These are used as in `x0kf_v4h.F`. Look into it. Here is a simplified version for showing how to use these data set.

            do 25 jpm  = 1, npm !
            do 25 ibib = 1, nbnb(k,jpm)
              n1b(ibib,k,jpm) !band index for k
              n2b(ibib,k,jpm) !band index for q+k
              it  =  n1b(ibib,k,jpm)                ! index for n  for q
              itp =  1+  n2b(ibib,k,jpm)-n2bminimum ! index for n' for q+k
              ... M(:,:) is calculated here ...  
                    call get_zmelt2 --> Get zmel=<MPB psi_k|psi_{k+q}>
                 zmel=dconjg(zmel)  ! -->  zmel= <psi_{k+q}| psi_k MPB_q>
                 M(igb1,igb2) = dconjg(zmel(igb1,it,itp))* zmel(igb2, it,itp)
                 Here M = <M_igb1 psi_it | psi_itp> < psi_itp | psi_it M_ibg2 >
              do iw = ihw(ibib,k,jpm),ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1 
                 imagweight = whw(jhw(ibib,k,jpm)+iw-ihw(ibib,k,jpm))
                 Kmatrix(:,:,iw,jpm) = Kmatrix(:,:,iw,jpm) + M(:,:)*imagweight 
              enddo          
       25  continue

We do $\sum_n \sum_{n'}$ in Eq. ((I78)) as follows;
before this `do 25`, we set ${\bf k},{\bf k}+{\bf q}$.
$k$ is index for ${\bf k}$ vector.
`n2bminimum` is introduced (this is for unoccupied states for `npm=1`.
$W_{\rm tet}$ is given as , `iw` is
the index for the histogram bin.

### Requirement of the spin symmetic Hamiltonian

(this idea is moved to another text modelSX.pdf.)

I think we can set up a model spin-symmetric screened exchange energy,
which is consistent with the linear response theory.
Thus we have two ways, one is calculate Heisenberg's $J$ by the
difference of total energy.
The other is by the linear response theory as in the magnon project.

## Harris-Foulkner energy and Kohn-Sham energy

In LDA/GGA, on the way to self-consistency, input density and output
density are not the same (not self-consitent).
Thus we define two total energy for given input density $n^{\rm in}$,
the Harris-Foulkner energy $E_{\rm Harris}$ and the Hohenberg-Kohn energy $E_{\rm HoKohn}$
$$\begin{aligned}
E_{\rm Harris}& =  E_{\rm k}^{\rm core} + E_{\rm B} 
- V_{\rm es}[n^{\rm Zc}+n^{\rm in},{\bf R}_a] \cdot n^{\rm in}
- V_{\rm xc}[n^{\rm c}+n^{\rm in}] \cdot n^{\rm in} \\
 &+ E_{\rm es}[n^{\rm Zc}+n^{\rm in},{\bf R}_a] + E_{\rm xc}[n^{\rm c}+n^{\rm in}],\\
E_{\rm HoKohn}& =  E_{\rm k}^{\rm core} + E_{\rm B} 
- V_{\rm es}[n^{\rm Zc}+n^{\rm in},{\bf R}_a] \cdot n^{\rm out}
- V_{\rm xc}[n^{\rm c}+n^{\rm in}] \cdot n^{\rm out} \\
 &+ E_{\rm es}[n^{\rm Zc}+n^{\rm out},{\bf R}_a] + E_{\rm xc}[n^{\rm c}+n^{\rm out}], \\
E_{\rm B} & =  \sum_p^{\rm occupied}
\alpha_{p}^{i*} 
\langle F_i|H^{\rm in}|F_j \rangle 
\alpha_p^j,
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
Search `ham_ehf,ham_ehk` in ecalj/lm7K/fp/bndfp.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span>.
The $E_{\rm Harris}$ is given at `call m_mkehkf_etot1(sev, eharris)`, while
$E_{\rm HoKohn}$ is given by\
`call m_mkehkf_etot2(sev,sumtv, eksham)`.
$E_{\rm Harris}$ and $E_{\rm HoKohn}$ are shown in `save.*` file.

The band energy $E_{\rm B}$ is stored in the module variable in bndfp.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span>;
this is shown as `sev` in `save.*`. The kinetic energy
$E_{\rm k}^{\rm core}+E_{\rm B} 
- V_{\rm es}[n^{\rm Zc}+n^{\rm in},{\bf R}_a] \cdot n^{\rm out}
- V_{\rm xc}[n^{\rm c}+n^{\rm in}] \cdot n^{\rm out}$ contained in
the $E_{\rm HoKohn}$ is calculated as `sumtv` in bndfp.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span>-mkekin.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span>

In LDA/GGA calculations by lmf-MPIK, `save.*` file contains a line per iteration.

    c ehf(eV)=-15730.0982239 ehk(eV)=-15730.0982222 sev(eV)=-15.7737393

shows $E_{\rm Harris}=$-15730.0982239 eV and $E_{\rm HoKohn}=$-15730.0982222 eV, as well as the valence
band energy $=-15.7737393$ eV.
In principle, $E_{\rm Harris}$ and $E_{\rm HoKohn}$ should be exactly the same when converged;
the difference is the numerical error. `c` at the begining of line means "converged".
${\tt h}$ means the 1st iteration from `atm.*` file (superposition of atomic density).

## Block inversion used for dielectric functions and downfolding

See Christph's and Pick's paper
$$\begin{aligned}
\left(\begin{array}{cc} P & Q \\ R & S \\ \end{array} \right) 
\left(\begin{array}{cc} W & -WQS^{-1} \\ -S^{-1}RW & S^{-1}+S^{-1}RWQS^{-1} \\ \end{array} \right) 
=\left(\begin{array}{cc} 1 & 0 \\ 0 & 1 \\ \end{array} \right), 
\end{aligned}
\tag{I81}
$$
where $P$ and $S$ are square matrices, and
$$\begin{aligned}
W=(P-QS^{-1}R)^{-1}. \end{aligned}
\tag{I82}
$$

We refer $X=-WQS^{-1}$ and $Y=-S^{-1}RW$.

Proof:
$$\begin{aligned}
 &{\rm (1,1) \ component}= PW-QS^{-1}RW=(P-QS^{-1}R)W = 1,\\
 &{\rm (1,2) \ component}= -PWQS^{-1}+QS^{-1}+QS^{-1}RWQS^{-1}\\ 
 &=(-PW+1+QS^{-1}RW)QS^{-1}=(-(P-QS^{-1}R)W+1)QS^{-1}=0\\
 &{\rm (2,1) \ component}= RW-SS^{-1}RW=0,\\
 &{\rm (2,2) \ component}= -RWQS^{-1}+SS^{-1}+SS^{-1}RWQS^{-1}=1
\end{aligned}
\tag{I83}
$$
(I think this proof is a little too complicated).

## Downfolding

Downfolding is a general concept which often appears in physics of varieties of contexts. It is based on the Block inversion Eq. ((I81)).
This appears for the inversion of one-body problem, or divide the Fock space
for many-body theory.
(divide the Fock space into two Hilbert spaces; one-particle excited states and
states with more than one-particles).

For exaple, the Green function is the inversion of the matrix $\omega-H$,
where $H$ is divided into to
$$\begin{aligned}
\left(\begin{array}{cc} H_{11} & H_{12} \\ H_{21}& H_{22} \\ \end{array} \right).
\end{aligned}
\tag{I84}
$$

If we have $H_{12}=H_{21}=0$, we have $G^0_{11}=1/(\omega -H_{11})$
as the main part of Green function. This is completely separated from
from the residual part $G_{22}^0=1/(\omega -H_{22})$.

When $H_{12}$ and $H_{21}$ are non-zero,
we have to take into their effect by perturbation, or by the block inversion
of Eq. ((I81)). Then we have
$$\begin{aligned}
\left(\begin{array}{cc} G_{11} & G_{12} \\ G_{21} & G_{22} \\ \end{array} \right)
=
\left(\begin{array}{cc} G_{11} & -G_{11}H_{12}G^0_{22} 
                    \\ -G^0_{22}H_{21}G_{11} & G^0_{22}+G^0_{22}H_{21}G_{11}H_{12}G^0_{22} \\ \end{array} \right),
\end{aligned}
\tag{I85}
$$
where $W$ in Eq. ((I81)) is $G_{11}$. Note that the $\omega$ dependence is
in $G^0_{11}$ and $G^0_{22}$. Here $G_{11}$ given as (See Eq. ((I82))).
$$\begin{aligned}
G_{11}=\frac{1}{\omega - H_{11} - H_{12}G^0_{22}H_{21} }.
\end{aligned}
\tag{I86}
$$
Then $H_{12}G^0_{22}H_{21}$ is identified as the self-energy.

There are possible cases about how to choose Hilbert space $\Omega_1$ and
$\Omega_2$ corresponding to the division Eq. ((I84)).

-   In the one-body problem, for example, we take $\Omega_1$ as low energy part
    and $\Omega_2$ as high energy part. or 3d parts and others.

-   In a case of many-body theory,
    $\Omega_1$ is the one-particle Fock space.
    Then $G_{11}$ is the one-body propagetor, and $G_{22}$ is many-body (two-
    and more particles-) propagetor.

-   We may take $\Omega_1$ as the model space, $\Omega_2$ as the residual space.
    Then we make like to identify $H_{11}$ as the Hubbard Hamiltonian.
    However, it is not so simple since $H_{11}$ contains screening effect
    due to the degree of freedom of $\Omega_2$. We have to consider a little more
    complicated downfolding procedure.

Warn: not be confused with the division of many-body Hamiltonian $H$ into
$H_0 + (H-H_0)$ for perturbation.

### Causality and analytic property

The equation $(\omega-H) G(\omega)=1$ is not well defined.
We need to consider original equation of motion in real time as
$(i \frac{\partial }{\partial t}-H) G(t-t')=\delta(t-t')$ where
we have to take into account the bondary condition of retarded response
responding to the source term of impulse at the right-hand side.

Causality means \"cause gives result\".
This is represented by the step function, for example, as
$\theta(t-t') \exp(i \omega_0 (t-t'))$,
whose FT gives $1/(\omega-\omega_0-i \delta)$. Thus the position of pole
(upper or lower plane) is important to determine the direction of time (real time representation).

Sum rule is related to the causality but a little different.
For example, sum rule for Imaginary part of $G_{11}$ is controlled
only by the beheavior $G_{11}(\omega)$ at $|\omega| \to \infty$. Thus, as long as
$H_{12}G^0_{22}H_{21} \to 0$ for $|\omega| \to \infty$,
the sum rule is satisfied.

## Spherical Harmonics and Real harmonics in ecalj

In ecalj, fortran codes are mainly based on the real harmonics $y_{lm}(\hat{\bf r})$,
instead of the usual sperical (complex) harmonics $Y_{lm}(\hat{\bf r})$.
The coefficients of eigenfunctions and so on are ordered as, e.g.
$(m=-2, m=-1, m=0, m=1,m=2)$ for $l=2$.
For example, `LMXA=4`, we have `(4+1)**2=25` harmonics,
ordered as
$y_{00},y_{-11},y_{01},y_{11},y_{-22},y_{-12},... ,y_{22},y_{-33},...y_{33},
y_{-44},...y_{44}$.

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

----------------------------------------------------\
The definition of $Y_{lm}(\hat{\bf r})$ are
$$\begin{aligned}
 &Y_{lm}(\theta, \phi)
=(-1)^m \left[ \frac{(2l+1)(l-m)!}{4 \pi (l+m)!} \right]^{\frac{1}{2} } P^m_l(\cos(\theta)) e^{i m \phi}, \\
 &P^m_l(x) = \frac{(1-x^2)^{m/2} }{2^l l!}\frac{d^{l+m} \ \ }{dx^{l+m} } (x^2-1)^l
\end{aligned}
\tag{I89}
$$

We take these definitions from\
(1)A.R.Edmonds, Angular Momentum in quantum Mechanics,
Princeton University Press, 1960,\
(2)M.E.Rose, Elementary Theory of angular Momentum,
John Wiley & Sons, INC. 1957,\
if necessary. The definition of spherical hermonics are the same in these books.

### SOC matrix for the Wannier functions

We calculate SOC matrix `hammsoc`<span class="clfix">(not in the present SRC, 2026-11-15)</span>
$\langle F_{ {\bf k}i \sigma} |H_{\rm SOC}|F_{ {\bf k}j \sigma} \rangle$ in `sugw.F`\
(`lmf-MPIK --jobgw --job --socmatrix` mode around the end of `genMLFhso`)
for all ${\bf k}$ points in the Brillowin zone. Here
$F_{ {\bf k}i \sigma}$ means the basis functions in the PMT method.
Eigenfunction is given as
$$\begin{aligned}
\Psi_{ {\bf k}n \sigma}=\sum_i F_{ {\bf k}i \sigma} \alpha^{\bf k}(i,n,\sigma) .
\end{aligned}
\tag{I90}
$$

Thus the required SOC matrix
$\langle W_{lm\sigma}({\bf r}) | H_{\rm SOC} |W_{lm'\sigma'}({\bf r}) \rangle$
(`hammso`<span class="clfix">(not in the present SRC, 2026-11-15)</span> in `hsocmat.F`<span class="clfix">(not in the present SRC, 2026-11-15)</span>) is calcualted in `hsocmat`<span class="clfix">(not in the present SRC, 2026-11-15)</span> from

1.  $\alpha^{\bf k}(i,n,\sigma)$ (`evec(1:nz,iband,iqbz,isp)` in `hsocmat.F`<span class="clfix">(not in the present SRC, 2026-11-15)</span>)

2.  $\langle F_{ {\bf k}i\sigma} |H_{\rm SOC}|F_{ {\bf k}j\sigma'} \rangle$
    (`hso(1:nz,1:nz,iqbz,ispc)` in `hsocmat.F`<span class="clfix">(not in the present SRC, 2026-11-15)</span>. `ispc=3` means spin off-diagonal)

3.  $\langle \psi_{ {\bf k}n\sigma}|G_{l\sigma}(r)y_{lm}(\hat{ {\bf r} })\rangle$
    (`dnk(iband,iwf,iqbz,isp)` in `hsocmat.F`<span class="clfix">(not in the present SRC, 2026-11-15)</span>)

See `hsocmat.F`<span class="clfix">(not in the present SRC, 2026-11-15)</span>. We have `evecw=matmul(evec,dnk)`.
From `hammso`<span class="clfix">(not in the present SRC, 2026-11-15)</span>, we calculate the trace of the
square of $H_{\rm SOC}$ as ${\rm Tr}(H_{\rm SOC}^2)$ in the space spanned by the Wannier functions.
By construction the trace is not dependent on $y_{lm}$ or on $Y_{lm}$ because projected fourteen Wannier functions (when $4f$) spans the same Hilbert space. The trace is equilannt to the sum of square of eigenvalues
of $H_{\rm SOC}$.

${\rm Tr}(H_{\rm SOC}^2)$, as well as
$\langle W_{lm \sigma}| H_{\rm SOC} |W_{lm'\sigma'} \rangle$,
can be Wannier dependent because the Wannier is not unique.
But we expect the dependence is not so large.

Important check point to evaluate the size of SOC is whether we can reproduce the SOC splitting in the original bands
(not the bands of the 14-space model but the bands of full QSGW).

## Crystal symmetry for GW calculation

[space group operation (space group rotation):]{.underline}\
In the preparation state of GW calculation, we have
3x3 matrices `symops(3,3,ig),ig=1,ngrp`
for point group part of the space group operations.
It is written in `call m_hamindex_init,m_hamindex.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span>:L128`) at `lmv7.F` (lmfgw-MPIK job=0 mode).
For example, in the case of GaAs, it shows 24 matrices in it.

In GW main routines such as `hsfp0.sc.m.F`,
we do `call genallcf_v3`. This set `symgg (=symops)` in
the module `m_genallcf_v3` as well as `invg` pointing the .
In the module `m_zmel.F`, we do `call mptauof`.
This set all required symmetry operations (find accompanied translation vector
and atom-site mapping information)
from point group operations. It gives informations of space group symmetry
(which atom is mapped to which atom, and so on).
See document at the begining of `subroutine mptauof`.

In `lmv7.F`(main routine), we have `m_hamindex` called for lmfgw mode=0.
This gives a file Hamindex, containing Hamindex.

The space group symmetry operation can be checked by
"lmchk" in advance. But be careful. It also includes inversion symmetry,
which means $|\psi_{-{\bf k} }({\bf r})|^2=|\psi_{\bf k}({\bf r})|^2$
(this is because of time-reversal symmetry as for real local $V({\bf r})$ (DFT case)).
Need attention as for LDA+U case of orbital moments appears (no time-reversal case).

The given crystal symmetry can be controlled by `SYMGRP` category
in `GWinput`. You can supply generators of space group (see
explanation in GWinput), or set `find`. If it is `find`,
we assume electronic structure has the symmetry of crystal structure.
Note that symmetry of electronic structure can be lower than the crystal symmetry.
For example, NiO is deformed along (111) direction in realitiy,
but we may like to use structure without such deformation for calcualtions.
In such a case, we may need to supply generators for lower symmetry.\
[AF symmetry]{.underline}:\
Recently, we add a new option `SYMGRP_AF`<span class="clfix">(not in the present SRC, 2026-11-15)</span>.
We have to set

    /home/usr2/h70252a/ecaljAFtest/NdTEST
    ~/SWJ/Nd2CuO4
    ...
        ATOM=Ndup  POS=  0.0000  {a}/2   {c}/2-{znd}*{c}  AF=1
        ATOM=Nddn  POS=  {a}/2   0.0000  {znd}*{c}       AF=-1
    SYMGRPAF r2x
    SYMGRP r4z::(-1/2,1/2,0) 
    ...

Here, AF=1,AF=-1 give AF pairs.
SYMGRP is to keep z axis. SYMGRPAF is an antiferro magnetic symmetry operation
(we add r2x+spin inversion), then the ATOM AF=1 is mapped to AF=-1.
Ask to T.Kotani as for AF symmetry setting, because it is relatively new.

## IBZ and EIBZ scheme

We include EIBZ procedure given in III.G. in Ref.[friedrich_efficient_2010].
It is for the polarization function 0fp0,hx0fp0_sc@ and
for self-energy \_sc@.\
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
$$\begin{aligned}
X=\sum_{ {\bf q}\in {\rm BZ} } F[g_{\bf q}]
= \sum_{A_i \in S_A } \sum_{ {\bf q}\in {\rm IBZ} } \frac{N_{\bf q} }{N_A} F[g_{A_i({\bf q})}] 
= \sum_{A_i \in S_A } \sum_{ {\bf q}\in {\rm IBZ} } \frac{N_{\bf q} }{N_A} F[A_i[g_{\bf q}]] .
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
generated from $N_{\bf q}$ by $S^{\bf k}_A$. Note that ${\bf k}$ is fixed in Eq. ((I93)).
Since $S_A^{\bf k}\subset S_A$, IBZ $\subset$ EIBZ(${\bf k}$). That is, equivalent $q$ poinst in IBZ
might be diffrenciated in EIBZ(${\bf k}$).

See Eq.(50) in Ref.[friedrich_efficient_2010].
`call Seteibz -->call eibzgen` returns $S_A^{\bf k}$ stored in the modele `m_eibz`<span class="clfix">(not in the present SRC, 2026-11-15)</span>.

Note normazization is
$$\begin{aligned}
N=\sum_{\bf q}1 = 
\frac{1}{N_A^{\bf k} } \sum_{A^{\bf k}_i \in S_A^{\bf k} } 
\sum_{ {\bf q}\in {\rm EIBZ}({\bf k})}N_{\bf q}^{\bf k}=
\sum_{ {\bf q}\in {\rm EIBZ}({\bf k})}N_{\bf q}^{\bf k}\end{aligned}
\tag{I94}
$$

[For poralization]{.underline}\
We apply Eq. ((I93)) to $F^{\bf k}[g_{\bf q}]$, resulting
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
in `hx0fp0.sc.m.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span>,hx0fp0.m.F` or newly developing `hrcxq`.

We apply Eq. ((I95)) to Eq.(22) in Ref.[friedrich_efficient_2010]
resulting Eq.(51) in Ref.[friedrich_efficient_2010].
Then we need to rotate MPB for given space group rotation $A^{\bf k}_i$.
Note that $A^{\bf k}_i$ keeps ${\bf k}$ of MPB. This is in the next section.

Note that `call eibzgen` stores not only $S_A^{\bf q}$(`eibzsym`),
but also the weight $N_{\bf q}^{\bf k}$ (`nwgt`) in the module `m_eibz`<span class="clfix">(not in the present SRC, 2026-11-15)</span>.
(role of ${\bf q}$ and ${\bf k}$ may be turn around in `m_eibz`<span class="clfix">(not in the present SRC, 2026-11-15)</span>)

 \
[For self-energy]{.underline}\
In `hsfp0.sc.m.F`(not in `hsfp0.m.F`), we have symmetrization
as Eq.(52) in Ref.[friedrich_efficient_2010].
But used formula is a little different.

We use Eq. ((I94)) as it is.
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

## Rotation of eigenfuncitons and MPB by the space-group operations {#spacegroup-rotation}

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
A_{ {\bf R}u L}({\bf r-R-T}) \exp(i {\bf kT})$ in Eq. ((I3)).
(Here we get the angular momentum index $L$. $u$ is for radial part of index.)
Then we have
$$\begin{aligned}
 &g[A^{\bf k}_{ {\bf R}u L}]({\bf r})=A^{\bf k}_{ {\bf R}u L}(g^{-1}({\bf r})+{\it \Delta}^{-1})
=\sum_{\bf T}A_{ {\bf R}u L}(g^{-1}({\bf r})+{\it \Delta}^{-1} -{\bf R}-{\bf T}) \exp(i {\bf kT})  \\
 &=\sum_{\bf T}A_{ {\bf R}u L}(g^{-1}({\bf r}+g({\it \Delta}^{-1}) -g({\bf R})-g({\bf T})))
 \exp(i {\bf kT})  \\
 &=\sum_{\bf T}\sum_{L'}  A_{ {\bf R}u L'} ({\bf r}+g({\it \Delta}^{-1}) -g({\bf R})-g({\bf T}))D_{L'L}(g) \exp(i {\bf kT})  \\
 &=\sum_{\bf T}\sum_{L'}  A_{ {\bf R}' u L'}({\bf r}-{\bf R}'-{\it \Delta}{\bf T}_{\bf R}-g({\bf T})) D_{L'L}(g)\exp(i {\bf kT})  \\
 &=\sum_{ {\bf T}''} \sum_{L'}  A_{ {\bf R}' u L'}({\bf r}-{\bf R}'-{\bf T}'') D_{L'L}(g)
 \exp(i g({\bf k}){ {\bf T} }'') \exp(-i g({\bf k}){\it \Delta}{\bf T}_{\bf R})  \\
 &= \sum_{L'}  A^{g({\bf k})}_{ {\bf R}' u L'}({\bf r}) D_{L'L}(g) \exp(-i g({\bf k}){\it \Delta}{\bf T}_{\bf R}),
\end{aligned}
\tag{I100}
$$
where we use $A_{ {\bf R}' u L}({\bf r})=A_{ {\bf R}u L}({\bf r})$ since $R$ and $R'$
sites should be equilalent sites. We use cancellation ${\it \Delta}=-g({\it \Delta}^{-1})$.
$D_{L'L}(g)$ is the rotaion matrix of the real spherical harmonics for $g$.

### space-group rotation of PMT eigenfunction

We have eigenfunction rotation routine `subroutine rotwvigg` called from `zsecsym.F`<span class="clfix">(not in the present SRC, 2026-11-15)</span>.


With this Eq. ((I100)), coeffecients of eigenfunctions for the MTO part of
eigenfunctions are rotated (mapped) by the space group operation $g$.
It is coded in `subroutine rotwvigg`(eigenfunction rotation routine) as

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

In order to call rotwvigg, we have to call `readhamindex()` so as to set up variables
`ltab,ktab,offl,oflrev,...` in `m_hamindex` module in advance.

Orbital index tables (`_tbl`<span class="clfix">(not in the present SRC, 2026-11-15)</span> is read from a file `MNLA_CPHI`
(See Sec.(ref sec:atmnla in the original)) by `readmnla_cphi`.
The last column in `MNLA_CPHI` contains the orbital-block index.
The orbital-blocks are specified by $l$ and $k$, its size is $2l+1$.
The number of $i$ is exactly the number of lines of `MNLA_CPHI`.\

In `subroutine rotwvigg`(eigenfunction rotation routine) we also have part for the APW.
Formula is
$$\begin{aligned}
 &g[\exp(i ({\bf k}+{\bf G}){\bf r})]  
=\exp(i ({\bf k}+{\bf G}) (g^{-1}({\bf r})+{\it \Delta}^{-1})) 
=\exp(i (g({\bf k})+g({\bf G})) ({\bf r}+g({\it \Delta}^{-1}))  \\
 &=\exp(i (g({\bf k}+{\bf G}) {\bf r}) \exp(i ({\bf k}+{\bf G}) {\it \Delta}^{-1})  \\
 &= \exp({i (\bar{g}({\bf k}) + g({\bf G})+ (g({\bf k})-\bar{g}({\bf k}))) {\bf r} }) \exp({-i ({\bf k}+ {\bf G}){\it \Delta} }).
\end{aligned}
\tag{I101}
$$
(Our definition of PW (IPW) has 'no nomalization factor'.)
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
MTpart+IPWpart as in Eq. ((I7)) in the GW calculations.
In the form of MTIPW, we have to map eigenfunctions by space-group operations.
As in the manner of `rotwvigg`, we have `rotmto` and `rotipw`.
These are for rotation MTpart and IPW parts.

### space-group rotation of Mixed Product Basis

We need all the matrix elements $\langle \Psi_{ {\bf q+k}n'}| \Psi_{ {\bf k}n} E^{\bf q}_\nu \rangle$.
for given inputs of ${\bf q},{\bf k},n,n'$. Here $E^{\bf q}_\nu$ is the Coulomb-diagonalized mixed product basis
[friedrich_efficient_2010].

In the procedure of EIBZ (See Eq. ((I95)) around), we have to rotate MPB.\
Module `module m_rotMPB2 (in m_rotMPB.F<span class="clfix">(the source is now *.f90; check the name, 2026-11-15)</span>)` is for generating the rotation matrix `zrotm`.
With `zrotm`, we rotate MPB in `call X0kf_v4hz_symmetrize`, where
we multiple conversion matrix `\zcousq` for conversion from MPB to $\{\}E^{\bf q}_\nu\}$ as\
`  rcxq(:,:,iw,jpm) = matmul(zcousqc,matmul(rcxq_core,zcousq))` (see code in `x0kf_v4hz_symmetrize`<span class="clfix">(not in the present SRC, 2026-11-15)</span>).

`zrotm` is generated by calling `rotmto2,rotipw2`.
These are essentially similar with `rotmto,rotipw`, but the basis set in the MT region,
and the basis set for the interstitial region are different from the case of eigenfunctions.
Thus we supply `ibas_tbl` and so on as

          use m_pbindex,only: norbt, ibas_tbl,l_tbl,k_tbl,offset_tbl,offset_rev_tbl,
         &     max_ibas_tbl,max_l_tbl,max_k_tbl,max_offset_tbl

in `module m_rotMPB2`.
