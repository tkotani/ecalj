# Product basis set in the 𝘎𝘞 calculation

## Product basis $M$
The basis set to represent product of one-particle wave functions.
This is the way to reduce the dimension of product.
$$
\ket{M^𝒌_I}  = \{ \ket{φ^{𝒌_1*}_{𝑹u} φ^{𝒌_2}_{𝑹u'}},  \ket{P^{𝒌_1*}_{𝑮_1} P^{𝒌_2}_{𝑮_2}} \}
 \equiv \{ \ket{M^𝒌_I}_\mathrm{MT}, \ket{M^𝒌_{𝑮}}_\mathrm{IPW} \}
$$
,where $𝒌 = 𝒌_2 - 𝒌_1, 𝑮 = 𝑮_2 - 𝑮_1$.
Some of pairs of ($𝒌_1$,$𝒌_2$) or ($𝑮_1$,$𝑮_2$) which has the same $𝒌$ or $𝑮$ are got rid of in the set of $M$.

::: warning
In this notation, $\ket{φ^{𝒌_1*}_{𝑹u} φ^{𝒌_2}_{𝑹u'}}$ DOES NOT indicate the direct product between $\ket{φ^{𝒌_1*}_{𝑹u}}$ and $\ket{φ^{𝒌_2}_{𝑹u'}}$.
$\braket{ 𝒓| φ^{𝒌_1*}_{𝑹u} φ^{𝒌_2}_{𝑹u'}} = φ^{𝒌_1*}_{𝑹u}(𝒓)φ^{𝒌_2}_{𝑹u'}(𝒓)$
:::
## Product basis $E$

The MPB $M$ introduced above does not satisfy orthogonality. Therefore, we introduce an orthogonal basis.
This basis also diagonalizes the Coulomb matrix. By doing so, the calculation of the exchange self-energy becomes easier.
In equations, new product basis $E$ is represent by the liner combination of $M$, i.e.:
$$
\ket{E^𝒒_μ} = \sum_{I} z^𝒒_{μI}\ket{M^𝒒_I}.
$$
Then, $E$ is satisfy the following releations.
$$
\braket{E^𝒒_μ|E^𝒒_ν} = δ_{μν},
v(𝒒) \ket{E^𝒒_μ} =v_μ(𝒒)\ket{E^𝒒_μ}
$$
where $v(𝒒)$ is Coulomb matrix and $v_μ(𝒒)$ is eigen value. The coefficient $z^𝒒_{μI}$ and $v_μ(𝒒)$ are obtained by soliving following generalized eigenvalue equation:
$$
\sum_J (v_{IJ}^{𝒒} - v_{\mu}(𝒒) O^{𝒒}_{IJ} ) z^{𝒒}_{\mu J} = 0,
$$
where $v_{IJ}^{𝒒}$  is Coulomb matrix represented by $M$, namely, $\braket{M^{𝒒}_{I}|v|{M^{𝒒}_J}}$. 
By using $E$, Coulomb interaction operator is represented as follows:
$$
v(𝒒)=\sum_{μ} \ket{E^𝒒_μ}v_μ(𝒒)\bra{E^𝒒_μ}.
$$

::: warning about $\braket{M^{𝒒}_{I}|v|{M^{𝒒}_J}}$
Since $v$ is a non-local function, the calculation of this matrix element includes cross terms of $M_\text{MT}$ and $M_\text{IPW}$.
:::

