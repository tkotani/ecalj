# Appendix: mathematical notes {#appendix}

> English version of the page [theory/appendix](../../theory/appendix.md). The pages at the top of ecaljdoc are the reference; this one was made from theory/appendix.md of commit 1dd77bb12 (source blob 1fdc16fa378e).

Mathematical tools used in the other pages.

## Block inversion {#block-inversion}

Used for the dielectric function and for downfolding ([§2](#downfolding)). With square matrices $P$ and $S$,

$$
\begin{pmatrix} P & Q \\ R & S \end{pmatrix}
\begin{pmatrix} W & -WQS^{-1} \\ -S^{-1}RW & S^{-1}+S^{-1}RWQS^{-1} \end{pmatrix}
=\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix},
\tag{A1}
$$

$$
W=(P-QS^{-1}R)^{-1}.
\tag{A2}
$$

We write the off-diagonal blocks as $X=-WQS^{-1}$ and $Y=-S^{-1}RW$.

To check it, compute the four blocks of the product:

$$
\begin{aligned}
(1,1)&: \ PW-QS^{-1}RW=(P-QS^{-1}R)W = 1,\\
(1,2)&: \ -PWQS^{-1}+QS^{-1}+QS^{-1}RWQS^{-1}\\
&\qquad=\big(1-(P-QS^{-1}R)W\big)QS^{-1}=0,\\
(2,1)&: \ RW-SS^{-1}RW=0,\\
(2,2)&: \ -RWQS^{-1}+SS^{-1}+SS^{-1}RWQS^{-1}=1.
\end{aligned}
\tag{A3}
$$

## Downfolding {#downfolding}

Downfolding is a general idea that appears in many contexts of physics, based on the block inversion (A1).
It appears when a one-body problem is divided into two spaces, or when the Fock space of a many-body problem is divided into two
(one-particle excited states, and states with more particles).

For example, the Green function is the inverse of $\omega-H$, and we divide $H$ as

$$
H=\begin{pmatrix} H_{11} & H_{12} \\ H_{21}& H_{22} \end{pmatrix}.
\tag{A4}
$$

If $H_{12}=H_{21}=0$, the main part $G^0_{11}=1/(\omega -H_{11})$ is completely separated from the rest $G_{22}^0=1/(\omega -H_{22})$.

When $H_{12}$ and $H_{21}$ are not zero, their effect is taken into account by perturbation, or by the block inversion (A1). $W$ of (A1) is $G_{11}$, and

$$
\begin{pmatrix} G_{11} & G_{12} \\ G_{21} & G_{22} \end{pmatrix}
=
\begin{pmatrix} G_{11} & -G_{11}H_{12}G^0_{22} \\ -G^0_{22}H_{21}G_{11} & G^0_{22}+G^0_{22}H_{21}G_{11}H_{12}G^0_{22} \end{pmatrix},
\tag{A5}
$$

$$
G_{11}=\frac{1}{\omega - H_{11} - H_{12}G^0_{22}H_{21} }
\tag{A6}
$$

(from (A2)). The $\omega$ dependence enters through $G^0_{22}$. $H_{12}G^0_{22}H_{21}$ is the self-energy.

The spaces $\Omega_1$ and $\Omega_2$ of the division (A4) can be chosen, for example, as follows.

- In a one-body problem, $\Omega_1$ is the low-energy part and $\Omega_2$ the high-energy part; or the 3d part and the others.
- In a many-body problem, $\Omega_1$ is the one-particle Fock space. Then $G_{11}$ is the one-body propagator, and $G_{22}$ the propagator of two or more particles.
- $\Omega_1$ is the model space and $\Omega_2$ the rest. One would like to regard $H_{11}$ as the Hubbard Hamiltonian, but it is not so simple:
  $H_{11}$ should include the screening by the degrees of freedom of $\Omega_2$, which needs a more involved downfolding.

This is different from dividing the many-body Hamiltonian as $H=H_0 + (H-H_0)$ for perturbation.

## Causality and analytic property {#causality}

The equation $(\omega-H) G(\omega)=1$ alone does not determine $G$. We go back to the equation of motion in real time,

$$
\Big(i \frac{\partial }{\partial t}-H\Big) G(t-t')=\delta(t-t'),
\tag{A7}
$$

with the boundary condition that the response comes after the impulse (the source term) on the right-hand side ([Basic I](./basic), impulse response).

Note that this response function need not be a number. The same holds for $G(x,x',t-t')$ with $H$ a matrix $H(x,x')$ and the right-hand side $\delta(x-x')\delta(t-t')$.
Finally, both sides may be multiplied from the left by the same matrix $C$. In matrix notation,

$$
C\Big(i \frac{\partial }{\partial t}-H\Big) G(t-t')=C\,\delta(t-t').
\tag{A8}
$$

Here $H$, $C$ and $G$ are matrices with the index $x$; the components of $G$ are $G(x,x',t-t')$.

In fact, such forms are used in many-body perturbation theory.
In short, a Green function (a linear response) can be defined for any inhomogeneous differential equation whose right-hand side is a source term $\delta(t-t')$ at the point $x'$.
The matrix $C$ on the left can be regarded as the sensitivity of a sensor placed at the point $x$.
Abstractly, $G$ is the integral kernel that gives the particular solution of a general homogeneous linear equation to which a source term is added.
Once this kernel is known, for the equation

$$
\Big(i \frac{\partial }{\partial t}-H\Big) R(t)=S(t)
\tag{A9}
$$

that maps a source $S$ to a response $R$, the solution multiplied by $C$ from the left is written at once as the integral

$$
C R(t)=\int dt'\, C\, G(t-t')\, S(t')
\tag{A10}
$$

(matrix notation; in components $R(x,t)$ and $S(x,t)$, and the right-hand side contains the sum over $x'$.)
By causality $G(t-t')=0$ for $t<t'$, so the integral runs only over $t'\le t$.
In a word, (A9) has the form $AR=S$ with $A=i\frac{\partial}{\partial t}-H$, and $A^{-1}$ is the Green function. Note, however, that the inverse matrix can be taken only after the boundary condition (initial condition) of the retarded type is imposed.
The Green function here is the retarded one. In many-body perturbation theory the complicated T-product appears, but it can be regarded as a combination of retarded functions (the theory can be built without the T-product, with only the Green function of electrons (particles) going forward in time and the Green function of holes going forward in time).

Causality means "the cause gives the result", which is represented by the step function. For example, the Fourier transform of
$-i\theta(t-t') e^{-i \omega_0 (t-t')}$ is $1/(\omega-\omega_0+i \delta)$ (the retarded Green function, the same convention as [Basic I](./basic)).
Whether the pole is in the upper or lower half plane determines the direction of time.

The sum rule is related to causality but is a little different. For example, the sum rule for the imaginary part of $G_{11}$ is determined only by
the behavior of $G_{11}(\omega)$ at $|\omega| \to \infty$. Thus the sum rule is satisfied as long as $H_{12}G^0_{22}H_{21} \to 0$ for $|\omega| \to \infty$.
