# Matrix elements

::: warning 工事中
このページは書きかけです。記号の定義と式は、[ecalj のバンド計算法](./pmt) と、[QSGW](./qsgw) の第一の文献（Obata, Kotani, Oda 2025）の §2 のほうが新しく正確です。
:::

## Calculation of the matrix elements

$$
\braket{E^𝒒_μ Ψ_{𝒌n} |Ψ_{𝒒+𝒌n'}} = \sum_{I} z^{𝒒*}_{μI} \braket{M^𝒒_I Ψ_{𝒌n} |Ψ_{𝒒+𝒌n'}}
$$
<!-- 基底関数は,MT とMT外完全に空間的に分離されているので, MT内(原子基底)とMT外(IPW)の寄与に分離できる。 -->
::: info
Regarding $\ket{E^𝒒_μ}$, it is a linear combination of $\ket{M^𝒒_I}$, so it is not separated into inside or outside of MT for a certain $μ$.
:::

In the code these matrix elements are called `zmel` (module `m_zmel`, `SRC/subroutines/m_zmel.f90`).

For the case where MPB is a function within MT:
$$
\braket{M^𝒒_I Ψ_{𝒌n} |Ψ_{𝒒+𝒌n'}} = 
\sum_{𝑹uu'} α_{𝑹u}^{𝒌n*} α_{𝑹u'}^{𝒒+𝒌n'} \braket{M^𝒒_I φ^{𝒌}_{𝑹u}| φ^{𝒒+𝒌}_{𝑹u'}}
$$

$
\braket{M^𝒒_I φ^{𝒌}_{𝑹u}| φ^{𝒒+𝒌}_{𝑹u'}}
$ : variable name PPB

In the case of IPW:
$$
\braket{M^𝒒_𝑮  Ψ_{𝒌n} |Ψ_{𝒒+𝒌n'}} = 
\sum_{𝑮_1 𝑮_2} β_{𝑮_1}^{𝒌n*} β_{𝑮_2}^{𝒒+𝒌n'} \braket{M^𝒒_𝑮 P^{𝒌}_{𝑮_1}| P^{𝒒+𝒌}_{𝑮_2}}_\mathrm{I}
$$

## inside the MT

## IPW part 
The calculation for the region outside MT is obtained by calculating the MT part and subtracting it:
$$
\begin{align}
\braket{M^𝒒_𝑮 P^{𝒌}_{𝑮_1}| P^{𝒒+𝒌}_{𝑮_2}}_\mathrm{I} &= \int_{\mathrm I} d𝒓 e^{i(-𝑮-𝑮_1+𝑮_2)⋅𝒓 } \\
&= \delta_{-𝑮-𝑮_1+𝑮_2} - \sum_{𝑹} \int_\mathrm{MT} d𝒓_I e^{i(-𝑮-𝑮_1+𝑮_2)⋅(𝒓_I+𝑹) } \\
&≡ 𝒪_{𝑮+𝑮_1,𝑮_2}
\end{align}
$$
It should be noted (there might be a Fourier transform factor). $𝒓=𝑹+𝒓_I$.
The integral of $𝒓=𝒓_I$ is written in spherical coordinates because the integration region is a sphere:

$$
\begin{align}
\int_\mathrm{MT} d𝒓_I e^{i(-𝑮-𝑮_1+𝑮_2)⋅(𝒓_I+𝑹)} =  
\frac{4\pi e^{i(-𝑮-𝑮_1+𝑮_2)⋅𝑹}}{|-𝑮-𝑮_1+𝑮_2|} \int_0^{R} dr_I r_I\sin(|-𝑮-𝑮_1+𝑮_2|r_I)
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
\end{align}
$$
