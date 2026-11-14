# Basic expressions

> English version of the page [theory/basic](../../theory/basic.md). The pages at the top of ecaljdoc are the reference; this one was made from theory/basic.md of the tree of c9f9a3003 (source blob 0c3198c06f9f).

On the Green function at zero temperature

## References
- Fetter Walecka

## Definition of the (one-particle) Green function

$$
\begin{align}
iG(𝒓_1t_1, 𝒓_2t_2) &= \braket{\Psi_0|T[\hat 𝜓_\mathrm{H}(𝒓_1,t_1)\hat 𝜓^†_\mathrm{H}(𝒓_2,t_2)]|\Psi_0} \\
&=\braket{\Psi_0|\hat 𝜓_\mathrm{H}(𝒓_1,t_1)\hat 𝜓^†_\mathrm{H}(𝒓_2,t_2)|\Psi_0} 𝜃(t_1-t_2) -\braket{\Psi_0|\hat 𝜓^†_\mathrm{H}(𝒓_2,t_2)\hat 𝜓_\mathrm{H}(𝒓_1,t_1)|\Psi_0} 𝜃(t_2-t_1)
\end{align}
$$

- $\hat 𝜓_\mathrm{H}(𝒓,t)$ : the field operator in the Heisenberg picture, $\hat 𝜓_\mathrm{H}(𝒓,t) = e^{i\hat Ht/\hbar}\hat 𝜓(𝒓)e^{-i\hat Ht/\hbar}$.
- $\hat 𝜓^†_\mathrm{H}(𝒓,t) = e^{i\hat Ht/\hbar}\hat 𝜓^†(𝒓)e^{-i\hat Ht/\hbar}$.
- The field operators of fermions satisfy the anticommutation relation $\{𝜓_\mathrm{H}(𝒓_1,t),𝜓_\mathrm{H}^†(𝒓_2,t)\} = 𝛿(𝒓_1-𝒓_2)$. Note that this is the equal-time anticommutation relation.
- $T$ is the time-ordered product; it reorders the operators so that later times stand to the left and earlier times to the right.
- $𝜃(t)$: the step function, 1 for $t>0$ and 0 for $t<0$.
- $\Psi_0$: the $N$-particle ground state in the Heisenberg picture, $\hat H\ket{\Psi_0} = E \ket{\Psi_0}$

## Green function of a non-interacting system

### Hamiltonian
We take the Hamiltonian as $\hat H -μ\hat N$ and set the ground-state energy to 0. $𝜇$: the Fermi energy; $\hat N = \int d𝒓 𝜓^†(𝒓)𝜓(𝒓)$: the particle-number operator; with $N$ the number of particles in the ground state, $\hat N \ket{\Psi_0} = N \ket{\Psi_0}$.<span class="clfix">Below, the $\hat H$ of the Heisenberg picture $e^{i\hat Ht/\hbar}\cdots e^{-i\hat Ht/\hbar}$ is also read as this $\hat H-𝜇\hat N$ (for the density operator $\hat 𝜌_\mathrm{H}$ it makes no difference, since $\hat N$ commutes with $\hat 𝜌$).</span>
Without interaction the Hamiltonian is a sum over one-particle states,
$$
\hat H - \mu \hat N =  \sum_n 𝜀_n \hat c^†_{n}\hat c_{n}
$$
.

> [!WARNING]
>  Note that the energy reference (zero) is the Fermi level.

Here $𝜙_{n}(𝒓)$ are the eigenfunctions of the one-particle Hamiltonian (they form a complete set), and $\hat c_n$ are the coefficients of the expansion of the field operator in them.
$$
\hat 𝜓(𝒓) = \sum_n \hat c_n 𝜙_{n}(𝒓)
$$
In the ground state without interaction, the orbitals of energy $𝜀_n$ are occupied up to the Fermi level following the Fermi distribution.

### Real-time representation
The two processes, an electron created in an unoccupied state by $\hat 𝜓^\dagger \ket{\Psi_0}$, or a hole created in an occupied state by $\hat 𝜓 \ket{\Psi_0}$ (an occupied electron annihilated), correspond to the first and the second terms of the Green function.
$$
\begin{align}
\braket{\Psi_0|\hat 𝜓_\mathrm{H}(𝒓_1,t_1)\hat 𝜓^†_\mathrm{H}(𝒓_2,t_2)|\Psi_0} &= 
\braket{\Psi_0|\hat 𝜓(𝒓_1)e^{-i(\hat H- 𝜇 \hat N)(t_1-t_2)/\hbar}  \hat 𝜓^†(𝒓_2)|\Psi_0} \\
&= \sum_{n'}^{\mathrm{unocc}} 𝜙_{n'}(𝒓_1)𝜙^*_{n'}(𝒓_2)e^{-i𝜀_{n'}(t_1-t_2)/\hbar} :\ \text{propagator of an electron from creation to annihilation} \\
\braket{\Psi_0|\hat 𝜓^†_\mathrm{H}(𝒓_2,t_2)\hat 𝜓_\mathrm{H}(𝒓_1,t_1)|\Psi_0} &= 
\braket{\Psi_0|\hat 𝜓^†(𝒓_2)e^{-i(\hat H- 𝜇 \hat N)(t_2-t_1)/\hbar} \hat 𝜓(𝒓_1)|\Psi_0} \\
&= \sum_{n}^{\mathrm{occ}} 𝜙_{n}(𝒓_1)𝜙^*_{n}(𝒓_2)e^{i𝜀_{n}(t_2-t_1)/\hbar}:\ \text{propagator of a hole from creation to annihilation}
\end{align}
$$
- Creation and annihilation act on the same one-particle state $n$; otherwise the term vanishes.
- A state with a created hole is also an excited state; its excitation energy from the ground state is $-𝜀_n>0$.
::: details Derivation

<div class="clfix">

Since the eigenvalue of $\hat H-𝜇\hat N$ in the ground state is set to 0, $e^{-i(\hat H-𝜇\hat N)t_2/\hbar}\ket{\Psi_0}=\ket{\Psi_0}$.
Of $\hat 𝜓^†(𝒓_2)\ket{\Psi_0}=\sum_{n'}𝜙^*_{n'}(𝒓_2)\hat c^†_{n'}\ket{\Psi_0}$ only the unoccupied $n'$ survive, and $(\hat H-𝜇\hat N)\hat c^†_{n'}\ket{\Psi_0}=𝜀_{n'}\hat c^†_{n'}\ket{\Psi_0}$.
Multiplying from the left by $\bra{\Psi_0}\hat 𝜓(𝒓_1)=\sum_{m}𝜙_{m}(𝒓_1)\bra{\Psi_0}\hat c_{m}$ and using $\braket{\Psi_0|\hat c_{m}\hat c^†_{n'}|\Psi_0}=𝛿_{mn'}$ gives the first expression.
The second follows in the same way: in $\hat 𝜓(𝒓_1)\ket{\Psi_0}$ only the occupied $n$ survive, and $(\hat H-𝜇\hat N)\hat c_{n}\ket{\Psi_0}=-𝜀_{n}\hat c_{n}\ket{\Psi_0}$.

</div>

:::

Hence,
$$
\begin{align}
iG(𝒓_1t_1,𝒓_2t_2) =& \sum_{n'}^{\mathrm{unocc}} 𝜙_{n'}(𝒓_1)𝜙^*_{n'}(𝒓_2)e^{-i𝜀_{n'}(t_1-t_2)/\hbar} 𝜃(t_1-t_2) \\
  &- \sum_{n}^{\mathrm{occ}} 𝜙_{n}(𝒓_1)𝜙^*_{n}(𝒓_2)e^{i𝜀_{n}(t_2-t_1)/\hbar}𝜃(t_2-t_1) \\
  & \equiv iG_\mathrm{ele}(𝒓_1t_1,𝒓_2t_2)+ iG_\mathrm{hole}(𝒓_1t_1,𝒓_2t_2)
\end{align}
$$
In the last line we define the electron and hole contributions separately. This helps in understanding products of the Green function with something else.
> [!INFO]
> The Green function is easier to picture in real space and real time.

### Frequency representation
Since the Green function is a function of $t_1-t_2$, its Fourier transform with respect to it is also used. It is called the Lehmann representation or the spectral representation.
As the Fourier representation of the step function we have
$$
𝜃(t) = - \int \frac{d𝜔'}{2𝜋i} \frac{e^{-i𝜔't}}{𝜔'+i𝛿}
$$
Here $𝛿$ is a positive infinitesimal. The step function is expressed by whether or not the residue in the lower half of the $𝜔'$ plane is picked up.

::: info
The Fourier transform has an ambiguity in the normalization factor; in this text we take the following.
$$
\begin{align}
f(𝜔) &= \int dt f(t) e^{i𝜔t} \\
f(t) &= \frac{1}{2\pi}\int d𝜔 f(𝜔) e^{-i𝜔t} \\
𝛿(𝜔) &= \frac{1}{2\pi} \int dt e^{i𝜔t}
\end{align}
$$
:::

Consider the Fourier transform of $𝜃(\pm t)e^{-iEt/\hbar}$ (for $𝜃(-t)$, replace $e^{-i𝜔't}$ in $𝜃(t)$ by $e^{i𝜔't}$):
$$
\begin{align}
\int dt 𝜃(\pm t)e^{-iEt/\hbar}e^{i𝜔t} = - \int \frac{d𝜔'}{2𝜋i} \int dt \frac{e^{-i(\pm 𝜔'+ E/\hbar-𝜔)t}}{𝜔'+ i𝛿} \\
= - \int \frac{d𝜔'}{i}\frac{1}{𝜔'+ i𝛿}𝛿(\pm 𝜔'-𝜔+ E/\hbar) = i\frac{1}{\pm (𝜔 - E/\hbar) + i𝛿}
\end{align}
$$
Then,

$$
\begin{align}
G(𝒓_1,𝒓_2;𝜔) &= \int d(t_1-t_2) G(𝒓_1t_1,𝒓_2t_2) e^{i𝜔(t_1-t_2)} \\
&=  \sum_{n'}^{\mathrm{unocc}} \frac{𝜙_{n'}(𝒓_1)𝜙^*_{n'}(𝒓_2)}{𝜔-𝜀_{n'}/\hbar + i𝛿} +
    \sum_{n}^{\mathrm{occ}} \frac{𝜙_{n}(𝒓_1)𝜙^*_{n}(𝒓_2)}{𝜔 - 𝜀_{n}/\hbar -i𝛿}
\end{align}
$$
is obtained. In the $𝜔$ plane extended to complex numbers, the poles of the unoccupied states (electron excitations) lie in the lower half plane at $\Re 𝜔 >0$, and the poles of the occupied states (hole excitations) in the upper half plane at $\Re 𝜔 <0$.

![Poles of the Green function](../../theory/figures/green_poles.png)

<span class="clfix">"occupied stats" and "unoccupied stats" in the figure should read "occupied states" and "unoccupied states". The script that drew the figure is lost, so the figure has not been redrawn.</span>

> [!INFO]
> The $i$ in the definition of the real-time representation was attached so that the expression in the frequency representation becomes clean.

### Density of states

$$
\lim_{𝛿 →0}\frac{1}{x \pm i 𝛿} = \mathcal{P} \frac{1}{x} \mp i𝜋𝛿(x)
$$
Using the relation above,
$$
\begin{align}
\int d𝒓 G(𝒓,𝒓;𝜔) =  \sum_{n} \frac{1}{𝜔-𝜀_{n}/\hbar \mp i𝛿} = \sum_n \left( \mathcal{P} \frac{1}{𝜔-𝜀_n/\hbar} \pm i 𝜋 𝛿(𝜔-𝜀_n/\hbar)\right)
\end{align}
$$
<span class="clfix">(upper signs: occupied, $𝜀_n\lt 0$; lower signs: unoccupied, $𝜀_n>0$),</span> and formally the density of states <del class="clfix">$D(𝜀) = \sum_n 𝛿(𝜀-𝜀_n)$</del><span class="clfix">$N(𝜀) = \sum_n 𝛿(𝜀-𝜀_n)$ (to distinguish it from the polarization function $D$)</span> can be written as

<div class="clfix-old">

$$
\begin{align}
D(𝜀) = \pm \frac{1}{𝜋}\Im \text{Tr}G(𝜀/\hbar)\quad \text{where } (+:𝜀<0,\ -:𝜀>0)
\end{align}
$$
(are $\hbar$ and the signs correct??)

</div>

<div class="clfix">

$$
\begin{align}
N(𝜀) = \pm \frac{1}{𝜋\hbar}\Im \text{Tr}G(𝜀/\hbar)\quad \text{where } (+:𝜀<0,\ -:𝜀>0)
\end{align}
$$
The factor $1/\hbar$ comes from $𝛿(𝜀/\hbar-𝜀_n/\hbar)=\hbar\,𝛿(𝜀-𝜀_n)$. The signs follow the double signs above.

</div>

## Time-ordered polarization function

With the density operator in the Heisenberg picture $\hat 𝜌_\mathrm{H}(𝒓,t)$, $\hat 𝜌_\mathrm{H}(𝒓,t) = e^{i\hat Ht/\hbar}\hat 𝜓^†(𝒓)\hat 𝜓(𝒓)e^{-i\hat Ht/\hbar}$.<span class="clfix">(previously written with a stray `^` after $\hat 𝜓$)</span>
We introduce the operator of the deviation from the expectation value of the density, $\hat {𝛿\rho} = \hat \rho - \braket{\hat \rho}$.
The time-ordered polarization function $D(𝒓_1t_1,𝒓_2t_2)$ is then defined as the time-ordered product of densities (density correlation):
$$
\begin{align}
iD(𝒓_1t_1,𝒓_2t_2) = 
\braket{T[\hat {𝛿𝜌}(𝒓_1,t_1)\hat {𝛿𝜌}(𝒓_2,t_2)]} = 
\braket{\hat   {𝛿𝜌}(𝒓_1,t_1)\hat {𝛿𝜌}(𝒓_2,t_2)}𝜃(t_1-t_2)  \\
+\braket{\hat  {𝛿𝜌}(𝒓_2,t_2)\hat {𝛿𝜌}(𝒓_1,t_1)}𝜃(t_2-t_1)
\end{align}
$$

## Time-ordered polarization function of a non-interacting system
In a non-interacting system the normal-ordered (N-product) terms have zero ground-state expectation value, so by Wick's theorem only the full contractions remain:
$$
\begin{align}
iD(𝒓_1t_1,𝒓_2t_2) = 
\braket{T[\hat {𝛿𝜌}(𝒓_1,t_1)\hat {𝛿𝜌}(𝒓_2,t_2)]} = G(𝒓_1t_1,𝒓_2t_2)G(𝒓_2t_2,𝒓_1t_1)
\end{align}
$$
Recalling that $G(𝒓_1t_1,𝒓_2t_2) = G_\mathrm{ele}(𝒓_1t_1,𝒓_2t_2) + G_\mathrm{hole}(𝒓_1t_1,𝒓_2t_2)$ can be separated, with $G_\mathrm{ele}$ remaining for $t_1 > t_2$ and $G_\mathrm{hole}$ for $t_1 < t_2$,
$$
\begin{align}
iD(𝒓_1t_1,𝒓_2t_2) =& G_\mathrm{ele}(𝒓_1t_1,𝒓_2t_2)G_\mathrm{hole}(𝒓_2t_2,𝒓_1t_1) \\
                  &+ G_\mathrm{hole}(𝒓_1t_1,𝒓_2t_2)G_\mathrm{ele}(𝒓_2t_2,𝒓_1t_1) \\
               =& \sum_{n'}^\mathrm{unocc} \sum_{n}^\mathrm{occ} 𝜙_{n'}(𝒓_1) 𝜙^*_{n'}(𝒓_2) 𝜙_{n}(𝒓_2) 𝜙^*_{n}(𝒓_1)e^{-i(𝜀_{n'}-𝜀_{n})(t_1-t_2)/\hbar} \theta(t_1-t_2) \\
                & +\sum_{n}^\mathrm{occ} \sum_{n'}^\mathrm{unocc} 𝜙_{n}(𝒓_1) 𝜙^*_{n}(𝒓_2) 𝜙_{n'}(𝒓_2) 𝜙^*_{n'}(𝒓_1)e^{-i(𝜀_{n'}-𝜀_{n})(t_2-t_1)/\hbar} \theta(t_2-t_1)
\end{align}
$$
Its Fourier transform, in the same way as above, is
$$
\begin{align}
D(𝒓_1,𝒓_2;𝜔) =& \sum_{n'}^\mathrm{unocc} \sum_{n}^\mathrm{occ} 
                   \frac{𝜙_{n'}(𝒓_1) 𝜙^*_{n'}(𝒓_2) 𝜙_{n}(𝒓_2) 𝜙^*_{n}(𝒓_1)}{  𝜔-(𝜀_{n'}-𝜀_n)/\hbar + i𝛿} \\
              & +\sum_{n}^\mathrm{occ} \sum_{n'}^\mathrm{unocc} 
                   \frac{𝜙_{n}(𝒓_1) 𝜙^*_{n}(𝒓_2) 𝜙_{n'}(𝒓_2) 𝜙^*_{n'}(𝒓_1)}{-(𝜔-(𝜀_{n}-𝜀_{n'})/\hbar) + i𝛿} \\
            =& \sum_{n'}^\mathrm{unocc} \sum_{n}^\mathrm{occ} \left[
                   \frac{𝜙_{n'}(𝒓_1) 𝜙^*_{n'}(𝒓_2) 𝜙_{n}(𝒓_2) 𝜙^*_{n}(𝒓_1)}{  𝜔-(𝜀_{n'}-𝜀_n)/\hbar + i𝛿}
                  -\frac{𝜙_{n}(𝒓_1) 𝜙^*_{n}(𝒓_2) 𝜙_{n'}(𝒓_2) 𝜙^*_{n'}(𝒓_1)}{  𝜔+(𝜀_{n'}-𝜀_{n})/\hbar - i𝛿} \right]
\end{align}
$$
.
Let $E_\mu(𝒓)$ be an orthonormalized basis in which products of wave functions are expanded as $𝜙^*_n(𝒓)𝜙_m(𝒓) = \sum_{𝜇} Z_{nm𝜇}E_𝜇(𝒓)$, with expansion coefficients $Z_{nm𝜇}$
(satisfying $\braket{E_𝜇|E_𝜈} = 𝛿_{𝜇𝜈}$). The expansion coefficients of $D(𝜔)$ in $E_{𝜇}$ are then:
$$
\begin{align}
D_{𝜇𝜈}(𝜔) &= \int d𝒓_1 d𝒓_2 E_{𝜇}^*(𝒓_1) D(𝒓_1,𝒓_2;𝜔)  E_{𝜈}(𝒓_2) \\
         &= \sum_{n'}^\mathrm{unocc} \sum_{n}^\mathrm{occ} \left[
               \frac{Z_{nn'𝜇} Z^*_{nn'𝜈}}{  𝜔-(𝜀_{n'}-𝜀_n)/\hbar + i𝛿}
              -\frac{Z_{n'n𝜇} Z^*_{n'n𝜈}}{  𝜔+(𝜀_{n'}-𝜀_{n})/\hbar - i𝛿} \right]
\end{align}
$$
<!-- Because it is time-symmetrized, it contains both the retarded (poles in the lower $ω$ plane) and the advanced (poles in the upper $ω$ plane) parts -->


<!-- ::: details Derivation -->
<!-- $$ -->
<!-- \begin{align} -->
<!-- \braket{T[\hat {𝛿𝜌}(𝒓_1,t_1)\hat {𝛿𝜌}(𝒓_2,t_2)]} &=  -->
<!--    \braket{T[\hat 𝜌_\mathrm{H}(𝒓_1,t_1) \hat 𝜌_\mathrm{H}(𝒓_2,t_2)]} -->
<!--  - \braket{T[\hat 𝜌_\mathrm{H}(𝒓_1,t_1) \braket{\hat \rho_\mathrm{H}(𝒓_2,t_2)}]}  \\ -->
<!-- &- \braket{T[\hat 𝜌_\mathrm{H}(𝒓_2,t_2) \braket{\hat \rho_\mathrm{H}(𝒓_1,t_1)}]}  -->
<!-- + \braket{T[\braket{\hat \rho_\mathrm{H}(𝒓_1,t_1)} \braket{\hat \rho_\mathrm{H}(𝒓_2,t_2)}]} \\ -->
<!-- &= \braket{T[\hat 𝜌_\mathrm{H}(𝒓_1,t_1) \hat 𝜌_\mathrm{H}(𝒓_2,t_2)]}+{\braket{\hat \rho_\mathrm{H}(𝒓_1,t_1)} \braket{\hat \rho_\mathrm{H}(𝒓_2,t_2)}} \\ -->
<!-- \end{align} -->
<!-- $$ -->
<!---->
<!-- Since $\hat \rho_\mathrm{H}(𝒓,t) = \hat 𝜓^\dagger_\mathrm{H}(𝒓,t) \hat 𝜓_\mathrm{H}(𝒓,t)$,  -->
<!-- \braket{T[\hat 𝜓^†_\mathrm{H}(𝒓_1,t_1)𝜓_\mathrm{H}(𝒓_1,t_1)𝜓^†_\mathrm{H}(𝒓_1,t_1)𝜓_\mathrm{H}(𝒓_2,t_2)]} =  -->
<!-- ::: -->
