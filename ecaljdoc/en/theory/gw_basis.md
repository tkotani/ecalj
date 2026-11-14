# Basis set in the 𝘎𝘞 calculation

> English version of the page [theory/gw_basis](../../theory/gw_basis.md). The pages at the top of ecaljdoc are the reference; this one was made from theory/gw_basis.md of the tree of c9f9a3003 (source blob cfcb556370fb).

::: warning Under construction
This page is a draft. The definitions and equations in [Band-structure method of ecalj](./pmt) and in section 2 of the main reference on the [QSGW](./qsgw) page (Obata, Kotani, Oda 2025) are newer and more accurate.
:::

## Representation of one-particle wave functions
The wave functions are expanded in two types of basis functions as follows:
$$
\ket{Ψ_{𝒌n}} = \sum_{𝑹u} α_{𝑹u}^{𝒌n} \ket{φ^𝒌_{𝑹u}} + \sum_{𝑮} β_{𝑮}^{𝒌n} \ket{P^𝒌_{𝑮}},
$$
where
$\ket{φ_{𝑹u}}$ is an atomic MTO basis function (an augmented smooth Hankel function) and $\ket{P^𝒌_{𝑮}}$ is an interstitial plane wave (zero inside the MTs).

### Detail of $\ket{φ^𝒌_{𝑹u}}$

The $\ket{φ^𝒌_{𝑹u}}$ is the Bloch sum of the atomic orbital $\ket{φ_{𝑹u}}$:
$$
φ^𝒌_{𝑹u}(𝒓) = \sum_{𝑻} φ_{𝑹u}(𝒓 - 𝑹 - 𝑻) \exp(i𝒌⋅𝑻)
$$
The radial part of atomic orbital are constructed from $φ_{lν}(r)$, $\dot φ_{lν}(r)$, and $φ^{z}_{l}(r)$ by the augmentation of smooth Hankel functions.
<!-- ::: details $\braket{φ^𝒌_{𝑹u}|φ^𝒌_{𝑹'u'}}$ -->
<!---->
<!-- ::: -->

The index of the atom is omitted for simplicity. Each $l$ has 2 or 3 reference energy of $ϵ_{lν}$.
$φ_{lν}(r)$ is the $l$ channel orbital of solution on radial Schrödinger equation at energy $ϵ_{lν}$
,$\dot φ_{lν}(r)$ is energy derivative of $φ_{lν}(r)$, $φ^{z}_{l}(r)$ is a local orbital, also a solution of the radial Schrödinger equation, at an energy far from $ϵ_{lν}$.
With `[ham].readp = true` and `[ham].pnufix = true` (both written by `ctrlgenToml.py`), the boundary conditions of $φ_{lν}$ at the MT radius (the logarithmic derivatives) are read from the results of `lmfa` (`readp`) and kept fixed (`pnufix`), so the basis set is fixed during the QSGW cycle.

::: warning ?? About radial equation??
- Is LDA used?
$$
\left[ -\frac{1}{2} \frac{d^2}{dr^2} + V_\mathrm{eff}[n(r)] - \frac{Z}{r} + \frac{l(l+1)}{2r^2}  -ϵ_{nlν} \right] φ_{nlν}(r) = 0
$$

$$
n(r) = 2\sum_{nl} (2l+1) |φ_{nl}(r)|^2/r^2
$$

 - In the $n(r)$ in radial equation obtained from self-consistent calculation?
 - What is the boundary conditions of radial equation.
 - $ϵ_{lν}$ is not the solution of radial equation, the $φ_{lν}(r)$ does not satisfy the some boundary conditions, for example $φ_{lν}(r) → 0 (r → 0)$
 - how to set the reference energy in the case of $φ^{z}_{l}(r)$?
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
\end{align}
$$
Since the MT regions are cut out, these functions are not orthogonal; their overlap matrix is
$$
\braket{P^𝒌_{𝑮}|P^𝒌_{𝑮'}} = \delta_{𝑮,𝑮'} - \sum_{𝑹_I} \int_{MT_I} e^{i(𝒓+𝑹_I) ⋅(𝑮'-𝑮)} d𝒓
$$
where,
$$
\int_{MT_I} e^{i(𝒓+𝑹_I) ⋅(𝑮'-𝑮)} d𝒓
 = \frac{4\pi e^{i𝑹_I ⋅(𝑮'-𝑮)}}{|𝑮'-𝑮|} \int_0^{R_I} dr r \sin(|𝑮'-𝑮|r)
$$

### Formula of inner product
This is for checking the orthonormality of the wave functions, and helps understanding.
$$
\braket{Ψ_{𝒌n}|Ψ_{𝒌m}} = \sum_{𝑹uu'} α_{𝑹u}^{𝒌n*} α_{𝑹u'}^{𝒌m}\braket{φ^𝒌_{𝑹u}|φ^𝒌_{𝑹u'}} + \sum_{𝑮𝑮'} β_{𝑮}^{𝒌n*} β_{𝑮'}^{𝒌m}\braket{P^𝒌_{𝑮}|P^𝒌_{𝑮'}}
$$
The cross terms vanish.
