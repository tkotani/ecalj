# Product basis for expanding products of wave functions

> English version of the page [theory/gw_product_basis](../../theory/gw_product_basis.md). The pages at the top of ecaljdoc are the reference; this one was made from theory/gw_product_basis.md of commit be745e319 (source blob 3867242f623a).

::: warning Under construction
This page is a draft. For the definitions of the symbols and the formulas, [Band-structure method of ecalj](./pmt) and section 2 of the main reference of [QSGW](./qsgw) (Obata, Kotani, Oda 2025) are newer and more accurate.
:::

## Product basis $M$
The basis set to represent products of one-particle wave functions, with a reduced dimension.
$$
\ket{M^𝒌_I}  = \{ \ket{φ^{𝒌_1*}_{𝑹u} φ^{𝒌_2}_{𝑹u'}},  \ket{P^{𝒌_1*}_{𝑮_1} P^{𝒌_2}_{𝑮_2}} \}
 \equiv \{ \ket{M^𝒌_I}_\mathrm{MT}, \ket{M^𝒌_{𝑮}}_\mathrm{IPW} \}
\tag{PB1}
$$
,where $𝒌 = 𝒌_2 - 𝒌_1, 𝑮 = 𝑮_2 - 𝑮_1$.
The dimension is reduced as follows. A product of IPWs is again an IPW, $P^{𝒌_1*}_{𝑮_1}P^{𝒌_2}_{𝑮_2} = P^{𝒌}_{𝑮}$, so the IPW part of $M$ is the set of IPWs up to a cutoff of $|𝒌+𝑮|$. In the MT part, the products of radial functions that are nearly linearly dependent are removed.

::: warning
In this notation, $\ket{φ^{𝒌_1*}_{𝑹u} φ^{𝒌_2}_{𝑹u'}}$ DOES NOT indicate the direct product between $\ket{φ^{𝒌_1*}_{𝑹u}}$ and $\ket{φ^{𝒌_2}_{𝑹u'}}$.
$\braket{ 𝒓| φ^{𝒌_1*}_{𝑹u} φ^{𝒌_2}_{𝑹u'}} = φ^{𝒌_1*}_{𝑹u}(𝒓)φ^{𝒌_2}_{𝑹u'}(𝒓)$
:::
## Product basis $E$

The MPB $M$ introduced above is not orthogonal. Therefore, we introduce an orthonormal basis.
This basis also diagonalizes the Coulomb matrix. By doing so, the calculation of the exchange self-energy becomes easier.
The new product basis $E$ is a linear combination of $M$:
$$
\ket{E^𝒒_μ} = \sum_{I} z^𝒒_{μI}\ket{M^𝒒_I}.
\tag{PB2}
$$
Then $E$ satisfies
$$
\braket{E^𝒒_μ|E^𝒒_ν} = δ_{μν}, \qquad
v(𝒒) \ket{E^𝒒_μ} =v_μ(𝒒)\ket{E^𝒒_μ}
\tag{PB3}
$$
where $v(𝒒)$ is the Coulomb matrix and $v_μ(𝒒)$ its eigenvalue. The coefficients $z^𝒒_{μI}$ and $v_μ(𝒒)$ are obtained by solving the generalized eigenvalue problem
$$
\sum_J (v_{IJ}^{𝒒} - v_{\mu}(𝒒) O^{𝒒}_{IJ} ) z^{𝒒}_{\mu J} = 0,
\tag{PB4}
$$
where $v_{IJ}^{𝒒}$ is the Coulomb matrix in $M$, namely $\braket{M^{𝒒}_{I}|v|{M^{𝒒}_J}}$, and $O^{𝒒}_{IJ} = \braket{M^{𝒒}_{I}|M^{𝒒}_J}$ is the overlap matrix.
By using $E$, Coulomb interaction operator is represented as follows:
$$
v(𝒒)=\sum_{μ} \ket{E^𝒒_μ}v_μ(𝒒)\bra{E^𝒒_μ}.
\tag{PB5}
$$

::: warning about $\braket{M^{𝒒}_{I}|v|{M^{𝒒}_J}}$
Since $v$ is a non-local function, the calculation of this matrix element includes cross terms of $M_\text{MT}$ and $M_\text{IPW}$.
:::

