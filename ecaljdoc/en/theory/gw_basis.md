# Representation of wave functions in GW

> English version of the page [theory/gw_basis](../../theory/gw_basis.md). The pages at the top of ecaljdoc are the reference; this one was made from theory/gw_basis.md of commit 296542c3a (source blob 0958e9992c41).

::: warning Under construction
This page is a draft. For the definitions of the symbols and the formulas, [Band-structure method of ecalj](./pmt) and section 2 of the main reference of [QSGW](./qsgw) (Obata, Kotani, Oda 2025) are newer and more accurate.
:::

## Representation of one-particle wave functions
The wave functions are expanded in two types of basis functions as follows:
$$
\ket{Ψ_{𝒌n}} = \sum_{𝑹u} α_{𝑹u}^{𝒌n} \ket{φ^𝒌_{𝑹u}} + \sum_{𝑮} β_{𝑮}^{𝒌n} \ket{P^𝒌_{𝑮}},
\tag{B1}
$$
where
$\ket{φ_{𝑹u}}$ is a function confined in the MT sphere at $𝑹$ (a radial function times a spherical harmonic; zero outside the MT), and $\ket{P^𝒌_{𝑮}}$ is an interstitial plane wave (zero inside the MTs). The eigenfunctions of `lmf`, expanded in MTOs and APWs, are re-expanded in these two kinds of functions for the GW calculation.

### Detail of $\ket{φ^𝒌_{𝑹u}}$

The $\ket{φ^𝒌_{𝑹u}}$ is the Bloch sum of the atomic orbital $\ket{φ_{𝑹u}}$:
$$
φ^𝒌_{𝑹u}(𝒓) = \sum_{𝑻} φ_{𝑹u}(𝒓 - 𝑹 - 𝑻) \exp(i𝒌⋅𝑻)
\tag{B2}
$$
The radial part of $φ_{𝑹u}$ is one of $φ_{lν}(r)$, $\dot φ_{lν}(r)$ and $φ^{z}_{l}(r)$ (the functions used in the augmentation of the smooth Hankel functions); $u$ labels the radial function and $lm$.
<!-- ::: details $\braket{φ^𝒌_{𝑹u}|φ^𝒌_{𝑹'u'}}$ -->
<!---->
<!-- ::: -->

The index of the atom is omitted for simplicity. Each $l$ has 2 or 3 reference energy of $ϵ_{lν}$.
$φ_{lν}(r)$ is the $l$ channel orbital of solution on radial Schrödinger equation at energy $ϵ_{lν}$
,$\dot φ_{lν}(r)$ is energy derivative of $φ_{lν}(r)$, $φ^{z}_{l}(r)$ is a local orbital, also a solution of the radial Schrödinger equation, at an energy far from $ϵ_{lν}$.
With `[ham].readp = true` and `[ham].pnufix = true` (both written by `ctrlgenToml.py`), the boundary conditions of $φ_{lν}$ at the MT radius (the logarithmic derivatives) are read from the results of `lmfa` (`readp`) and kept fixed (`pnufix`), so the basis set is fixed during the QSGW cycle.

::: info About the radial equation
- The radial functions $φ_{lν}$ and $\dot φ_{lν}$ are made by solving the radial equation in the spherical part of the LDA (GGA) potential in the MT. The same holds in the QSGW cycle: the nonlocal QSGW potential does not enter the radial equation.
$$
\left[ -\frac{1}{2} \frac{d^2}{dr^2} + V_\mathrm{eff}(r) - \frac{Z}{r} + \frac{l(l+1)}{2r^2}  -ϵ_{lν} \right] r φ_{lν}(r) = 0
\tag{Bq1}
$$
- The electron density that determines the potential is the sum over the occupied states of the wave functions (the radial equation uses its spherical part in the MT):
$$
n(𝒓) = \sum_{𝒌n}^{\mathrm{occ}} |Ψ_{𝒌n}(𝒓)|^2
\tag{Bq2}
$$
- Boundary conditions: the solution regular at the origin is taken, with the logarithmic derivative ($P$) given at the MT radius. $ϵ_{lν}$ is the energy that matches this logarithmic derivative, not an eigenvalue of the radial equation (a solution vanishing at infinity). The logarithmic derivatives are determined by the spherical free-atom calculation of `lmfa` and written to `atmpnu.*`; `lmf` reads them (`[ham].readp = true`) and keeps them fixed during the cycle (`[ham].pnufix = true`). See §2 (note) and §6 of [Band-structure method of ecalj](./pmt).
- The local orbital $φ^{z}_{l}$ is also determined by its logarithmic derivative ($P_z$) at the MT radius; $P_z$ is also read from `atmpnu.*` written by `lmfa` and kept fixed. (With `pnufix = false`, `pnunew` updates $P$; an $l$ without a local orbital is solved at the energy of the center of gravity of the occupied states (the center of the partial DOS). This was used before, but was unstable.)
:::

### Detail of $\ket{P^𝒌_{𝑮}}$
The interstitial plane wave is expressed as follows:
$$
\begin{align}
{P^𝒌_{𝑮}}(𝒓) = 
\begin{cases}
0    & \text{if}~ 𝒓 \in \text{any MT} \\
  \exp(i (𝒌 + 𝑮)⋅𝒓)& \text{otherwise}
\end{cases}
\tag{B3}
\end{align}
$$
Since the MT regions are cut out, these functions are not orthogonal; their overlap matrix is
$$
\braket{P^𝒌_{𝑮}|P^𝒌_{𝑮'}} = \Omega\,\delta_{𝑮,𝑮'} - \sum_{𝑹_I} \int_{MT_I} e^{i(𝒓+𝑹_I) ⋅(𝑮'-𝑮)} d𝒓
\tag{B4}
$$
where $\Omega$ is the volume of the unit cell, and
$$
\int_{MT_I} e^{i(𝒓+𝑹_I) ⋅(𝑮'-𝑮)} d𝒓
 = \frac{4\pi e^{i𝑹_I ⋅(𝑮'-𝑮)}}{|𝑮'-𝑮|} \int_0^{R_I} dr\, r \sin(|𝑮'-𝑮|r)
\tag{B5}
$$

### Formula of inner product
This is for checking the orthonormality of the wave functions, and helps understanding.
$$
\braket{Ψ_{𝒌n}|Ψ_{𝒌m}} = \sum_{𝑹uu'} α_{𝑹u}^{𝒌n*} α_{𝑹u'}^{𝒌m}\braket{φ^𝒌_{𝑹u}|φ^𝒌_{𝑹u'}} + \sum_{𝑮𝑮'} β_{𝑮}^{𝒌n*} β_{𝑮'}^{𝒌m}\braket{P^𝒌_{𝑮}|P^𝒌_{𝑮'}}
\tag{B6}
$$
The cross terms vanish, since $φ_{𝑹u}$ is zero outside the MTs and $P^𝒌_{𝑮}$ is zero inside.
