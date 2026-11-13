# Basics 2

> English version of the page [theory/basic2](../../theory/basic2.md). The pages at the top of ecaljdoc are the reference; this one was made from theory/basic2.md of the tree of c9f9a3003 (source blob e86d73963d75).

`Under construction`

- On the relation among time-ordered, advanced and retarded functions
- On the Hilbert transform / Kramers-Kronig relation


The Green function and the physical quantities that go with it are time-ordered, and their analyticity as complex functions is poor.
Time-ordered quantities are used because Wick's theorem applies and they are useful in perturbation theory, but sometimes it is better to use functions with good analyticity.
For example, the Kramers-Kronig relation lets one construct the imaginary (real) part from the real (imaginary) part, and it can also be used as a technique for efficient numerical calculation; but the time-ordered Green function, which is regular neither in the upper nor in the lower half plane, cannot be used with the KK relation as it is.
The KK relation also serves as a barometer for checking numerical accuracy, so we examine the relation to the retarded and advanced functions, which have good analyticity.

## Hilbert's formula
$$
\lim_{𝜂 → 0 }\int_{-∞}^{∞} \frac{𝜒(ω')}{ω' - ω ± i 𝜂} dω' = 𝒫 \int_{-∞}^{∞} \frac{𝜒(ω')}{ω' - ω } dω' ∓ iπ𝜒(ω)
$$

## Relation among advanced, retarded and time-ordered functions
$\mathrm{R}$ denotes retarded and $\mathrm{A}$ advanced.

### Green function
Without time ordering, the propagator of the anticommutator of the field operators is the retarded (advanced) Green function.

$$ \begin{align}
iG(𝒓_1t_1, 𝒓_2t_2) &= \braket{\Psi_0|T[\hat 𝜓_\mathrm{H}(𝒓_1,t_1)\hat 𝜓^†_\mathrm{H}(𝒓_2,t_2)]|\Psi_0}  \\
iG^\mathrm{R}(𝒓_1t_1, 𝒓_2t_2) &=  \braket{\Psi_0|\{\hat 𝜓_\mathrm{H}(𝒓_1,t_1),\hat 𝜓^†_\mathrm{H}(𝒓_2,t_2)\}|\Psi_0} 𝜃(t_1-t_2) \\
iG^\mathrm{A}(𝒓_1t_1, 𝒓_2t_2) &= -\braket{\Psi_0|\{\hat 𝜓_\mathrm{H}(𝒓_1,t_1),\hat 𝜓^†_\mathrm{H}(𝒓_2,t_2)\}|\Psi_0} 𝜃(t_2-t_1) \\
\end{align} $$
The Fourier representation is
$$ \begin{align}
G(𝒓_1, 𝒓_2;𝜔)  &=  \sum_{n'}^{\mathrm{unocc}} \frac{𝜙_{n'}(𝒓_1)𝜙^*_{n'}(𝒓_2)}{𝜔-𝜀_{n'}/\hbar + i𝛿} +
    \sum_{n}^{\mathrm{occ}} \frac{𝜙_{n}(𝒓_1)𝜙^*_{n}(𝒓_2)}{𝜔 - 𝜀_{n}/\hbar -i𝛿} = G_{\mathrm{ele}} + G_{\mathrm{hole}} \\
G^\mathrm{R}(𝒓_1, 𝒓_2;𝜔)  &=  \sum_{n'}^{\mathrm{unocc}} \frac{𝜙_{n'}(𝒓_1)𝜙^*_{n'}(𝒓_2)}{𝜔-𝜀_{n'}/\hbar + i𝛿} +
    \sum_{n}^{\mathrm{occ}} \frac{𝜙_{n}(𝒓_1)𝜙^*_{n}(𝒓_2)}{𝜔 - 𝜀_{n}/\hbar +i𝛿} = G_{\mathrm{ele}} + G^†_{\mathrm{hole}} \\
G^\mathrm{A}(𝒓_1, 𝒓_2;𝜔)  &=  \sum_{n'}^{\mathrm{unocc}} \frac{𝜙_{n'}(𝒓_1)𝜙^*_{n'}(𝒓_2)}{𝜔-𝜀_{n'}/\hbar - i𝛿} +
    \sum_{n}^{\mathrm{occ}} \frac{𝜙_{n}(𝒓_1)𝜙^*_{n}(𝒓_2)}{𝜔 - 𝜀_{n}/\hbar -i𝛿} = G^†_{\mathrm{ele}} + G_{\mathrm{hole}}
\end{align} $$
Here $†$ means taking the complex conjugate and exchanging $𝒓_1, 𝒓_2$.
For $G^\mathrm{R}$ and $G^\mathrm{A}$ the poles in the $ω$ plane are separated into the lower half plane (retarded) and the upper half plane (advanced), and they satisfy the Kramers-Kronig (KK) relation
$$
\Re G^\mathrm{R,A} (𝒓_1, 𝒓_2;𝜔) = ± 𝒫 \int_{-∞}^{∞} \frac{dω'}{π} \frac{ \Im G^{\mathrm{R,A}}(𝒓_1, 𝒓_2;𝜔') }{ω' - ω}
$$
.

### Density response
The density response function has time-reversal symmetry.
$$ \begin{align}
iD(𝒓_1t_1, 𝒓_2t_2) &= \braket{T[\hat {𝛿𝜌}(𝒓_1,t_1)\hat {𝛿𝜌}(𝒓_2,t_2)]} = G(𝒓_1t_1,𝒓_2t_2)G(𝒓_2t_2,𝒓_1t_1) \\
iD^\mathrm{R}(𝒓_1t_1, 𝒓_2t_2) &= \braket{[\hat {𝛿𝜌}(𝒓_1,t_1),\hat {𝛿𝜌}(𝒓_2,t_2)]} 𝜃(t_1-t_2) \\
iD^\mathrm{A}(𝒓_1t_1, 𝒓_2t_2) &= \braket{[\hat {𝛿𝜌}(𝒓_1,t_1),\hat {𝛿𝜌}(𝒓_2,t_2)]} 𝜃(t_2-t_1)
\end{align} $$
> [!INFO]
> A problem at the end of a chapter of Fetter and Walecka:
>  $\Re D(ω) = \Re D^\mathrm{R}(ω)$ and $\Im D(ω) \mathrm{sgn}(ω) = \Im D^\mathrm{R}(ω)$


### $GW$ self-energy
$$ \begin{align}
𝛴(𝒓_1t_1, 𝒓_2t_2) &=iG(𝒓_1t_1, 𝒓_2t_2)W(𝒓_1t_1, 𝒓_2t_2) \\
                  &= iG_\mathrm{elec}(𝒓_1t_1, 𝒓_2t_2)W(𝒓_1t_1, 𝒓_2t_2) 𝜃(t_1-t_2)
                   + iG_\mathrm{hole}(𝒓_1t_1, 𝒓_2t_2)W(𝒓_1t_1, 𝒓_2t_2) 𝜃(t_2-t_1)
\end{align} $$
$$ \begin{align}
𝛴_{nn}(𝜔) &= (G_\mathrm{ele}W)(𝜔)_{nn} + (G_\mathrm{hole}W)(𝜔)_{nn} = 𝛴_\mathrm{ele}^{nn}(𝜔)+𝛴_\mathrm{hole}^{nn}(𝜔) \\
𝛴^\mathrm{R}_{nn}(𝜔) &= 𝛴_\mathrm{ele}^{nn}(𝜔)+𝛴_\mathrm{hole}^{nn}(𝜔) \\
\end{align} $$

## Hilbert transform

### Numerical calculation
There is a principal-value integral, but basically an ordinary numerical integration method is used, avoiding treating the divergent point numerically.
With the data points at $ω_{C_i}$, let $ω_{L_i}$ and $ω_{R_i}$ be the midpoints to the neighbouring data points.
- Zeroth approximation
$$ \begin{align}
F(ω) &= 𝒫 \int_{-∞}^{∞} dω' \frac{f(ω')}{ω' -ω} = \sum_{i} 𝒫 \int_{ω_{R_i}}^{ω_{L_i}} dω' \frac{f(ω')}{ ω' -ω} \\
& \simeq \sum_{i} f(ω_{C_i})𝒫 \int_{ω_{R_i}}^{ω_{L_i}} dω' \frac{1}{ω'-ω}
= \sum_{i} f(ω_{C_i}) \log{\frac{|ω_{L_i} -ω|} {| ω_{R_i}-ω|}}
\end{align} $$
 Usually the real (imaginary) part is computed on the same mesh points as the imaginary (real) part.
