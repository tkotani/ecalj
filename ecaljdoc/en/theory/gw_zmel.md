# Matrix elements used in GW

> English version of the page [theory/gw_zmel](../../theory/gw_zmel.md). The pages at the top of ecaljdoc are the reference; this one was made from theory/gw_zmel.md of commit 83026e73d (source blob de47a883c56c).

::: warning Under construction
This page is a draft. For the definitions of the symbols and the formulas, [Band-structure method of ecalj](./pmt) and section 2 of the main reference of [QSGW](./qsgw) (Obata, Kotani, Oda 2025) are newer and more accurate.
:::

## Calculation of the matrix elements

$$
\braket{E^𝒒_μ Ψ_{𝒌n} |Ψ_{𝒒+𝒌n'}} = \sum_{I} z^{𝒒*}_{μI} \braket{M^𝒒_I Ψ_{𝒌n} |Ψ_{𝒒+𝒌n'}}
\tag{Z1}
$$
<!-- 基底関数は,MT とMT外完全に空間的に分離されているので, MT内(原子基底)とMT外(IPW)の寄与に分離できる。 -->
::: info
Regarding $\ket{E^𝒒_μ}$, it is a linear combination of $\ket{M^𝒒_I}$, so it is not separated into inside or outside of MT for a certain $μ$.
:::

In the code these matrix elements are called `zmel` (module `m_zmel`, `SRC/subroutines/m_zmel.f90`): `zmel` is $\braket{E^𝒒_μ Ψ_{𝒌n} |Ψ_{𝒒+𝒌n'}}$ above.
Since the product basis $M$ is not orthogonal, $E$ (orthonormal, made from $M$ by the eigenvectors $z$ of the Coulomb matrix normalized with $O_{IJ}=\braket{M^𝒒_I|M^𝒒_J}$) is used.
The code multiplies $\braket{M^𝒒_I Ψ_{𝒌n} |Ψ_{𝒒+𝒌n'}}$ by $z$ (`m2e_prod_basis` = `zcousq`). Before 2025-10-10 it went through the expansion coefficients in $M$,
$\sum_J (O^{-1})_{IJ} \braket{M^𝒒_J Ψ|Ψ}$, and then multiplied by `ppovlz` $= O z$ (the `ZO^-1` of older comments); the result is the same.

For the case where MPB is a function within MT:
$$
\braket{M^𝒒_I Ψ_{𝒌n} |Ψ_{𝒒+𝒌n'}} = 
\sum_{𝑹uu'} α_{𝑹u}^{𝒌n*} α_{𝑹u'}^{𝒒+𝒌n'} \braket{M^𝒒_I φ^{𝒌}_{𝑹u}| φ^{𝒒+𝒌}_{𝑹u'}}
\tag{Z2}
$$

$
\braket{M^𝒒_I φ^{𝒌}_{𝑹u}| φ^{𝒒+𝒌}_{𝑹u'}}
$ : variable name PPB

In the case of IPW:
$$
\braket{M^𝒒_𝑮  Ψ_{𝒌n} |Ψ_{𝒒+𝒌n'}} = 
\sum_{𝑮_1 𝑮_2} β_{𝑮_1}^{𝒌n*} β_{𝑮_2}^{𝒒+𝒌n'} \braket{M^𝒒_𝑮 P^{𝒌}_{𝑮_1}| P^{𝒒+𝒌}_{𝑮_2}}_\mathrm{I}
\tag{Z3}
$$

## Inside the MT

## IPW part 
The calculation for the region outside MT is obtained by calculating the MT part and subtracting it, with $𝑲 \equiv 𝑮_2-𝑮_1-𝑮$:
$$
\begin{align}
\braket{M^𝒒_𝑮 P^{𝒌}_{𝑮_1}| P^{𝒒+𝒌}_{𝑮_2}}_\mathrm{I} &= \int_{\mathrm I} d𝒓 e^{i𝑲⋅𝒓 } \\
&= \Omega\,\delta_{𝑲,0} - \sum_{𝑹} \int_\mathrm{MT} d𝒓_I e^{i𝑲⋅(𝒓_I+𝑹) } \\
&≡ 𝒪_{𝑮+𝑮_1,𝑮_2}
\tag{Z4}
\end{align}
$$
Here $\Omega$ is the volume of the unit cell (as in Eq. (B4) of [Representation of wave functions in GW](./gw_basis)), and $𝒓=𝑹+𝒓_I$.
The integral over $𝒓_I$ is written in spherical coordinates because the integration region is a sphere:

$$
\begin{align}
\int_\mathrm{MT} d𝒓_I e^{i𝑲⋅(𝒓_I+𝑹)} &=  
\frac{4\pi e^{i𝑲⋅𝑹}}{|𝑲|} \\
&\quad \times \int_0^{R} dr_I r_I\sin(|𝑲|r_I)
\tag{Z5}
\end{align}
$$
::: info
It is also acceptable to expand the plane wave with spherical Bessel functions and consider only the contribution of $l=0$ as the result of the integration.
:::

## Total
$$
\begin{align}
\braket{E^𝒒_μ Ψ_{𝒌n} |Ψ_{𝒒+𝒌n'}} &= \sum_{I \in \mathrm{MT}} z^{𝒒*}_{μI} \sum_{𝑹uu'} α_{𝑹u}^{𝒌n*} α_{𝑹u'}^{𝒒+𝒌n'} \braket{M^𝒒_I φ^{𝒌}_{𝑹u}| φ^{𝒒+𝒌}_{𝑹u'}} \\
 &+ \sum_{𝑮} z^{𝒒*}_{μ𝑮} \sum_{𝑮_1 𝑮_2}𝒪_{𝑮+𝑮_1,𝑮_2} β_{𝑮_1}^{𝒌n*} β_{𝑮_2}^{𝒒+𝒌n'} 
\tag{Z6}
\end{align}
$$
