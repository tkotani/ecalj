# Basics of the electrical conductivity calculation {#conductivity}

> English version of the page [theory/conductivity](../../theory/conductivity.md). The pages at the top of ecaljdoc are the reference; this one was made from theory/conductivity.md of commit 6678ecd65 (source blob 1998761b20a8).

::: info Overview
As a **one-body problem**, the electrical resistance due to impurities is derived from linear response (the Kubo formula).
The key point is that the average over the impurity configurations is taken of a **product of Green functions** (averaging each one first and then multiplying is not enough).
The sign convention of the Green function is the same as in [Basic I](./basic); the retarded Green function is $G^R(E)=1/(E-H+i\delta)$ and the advanced one $G^A(E)=1/(E-H-i\delta)$.
For retarded response functions see also the [Appendix](./appendix#causality).
:::

## Setting {#setting}

Here we consider a system of atoms (impurities) placed at random in vacuum.
There are $N_i$ identical atoms in the volume $\Omega$, and their positions ${\bf R}_j$ are distributed uniformly and independently of each other.
Non-interacting electrons (mass $m$, $\varepsilon_{\bf k}=\hbar^2k^2/2m$) feel the potential of these atoms,

$$
V({\bf r})=\sum_{j} u({\bf r}-{\bf R}_j)
\tag{C1}
$$

($u$ is the potential of one atom). The number density of the impurities is $n_i=N_i/\Omega$, and the Fourier transform of $u$ is $u({\bf q})=\int d^3r\, e^{-i{\bf q}\cdot{\bf r}}u({\bf r})$.
What is observed is the average $\langle\cdots\rangle_{\rm imp}$ over the configurations $\{ {\bf R}_j \}$ of the impurities.

The discussion here carries over directly to impurities in a solid: replace the plane waves $|{\bf k}\rangle$ and $\varepsilon_{\bf k}$ by the Bloch states and band energies of the host crystal, and $u$ by the change of the potential caused by the substituted impurity.

**What is treated here is a one-body problem.** There is no interaction between the electrons; each electron moves independently in the static potential (C1) of the impurities.
The Hamiltonian is the one-body $H=\hat{\bf p}^2/2m+V({\bf r})$, and $G^R(E)=1/(E-H+i\delta)$ below is its one-body Green function.
Once a configuration of the impurities is fixed, the problem is nothing but solving the one-body Schrödinger equation.
The self-energy and the vertex correction that appear below look like many-body diagrams, but they all come from the average over the impurity configurations (not from interactions between the electrons).

## Kubo formula {#kubo}

### General linear response {#linear-response}

Suppose that a time-dependent external field $F(t)$ acts on the system through an operator $\hat B$:

$$
\hat H(t)=\hat H_0-\hat B\,F(t).
\tag{C2}
$$

The operator $\hat A$ with which the response is observed may itself contain $F$ (the current operator used below does).
Taking the expectation value $\langle\hat A\rangle(t)={\rm Tr}[\rho(t)\hat A(t)]$ to first order in $F$ gives the sum of the term in which the change of $\hat A$ itself, $(\partial\hat A/\partial F)F$, is averaged with the equilibrium density matrix $\rho_0$, and the term from the change $\delta\rho$ of the density matrix.
The latter is obtained by expanding, in the interaction picture, $\delta\rho_{\rm I}(t)=(i/\hbar)\int_{-\infty}^{t}dt'\,[\hat B_{\rm I}(t'),\rho_0]\,F(t')$ to first order. The result is

$$
\begin{aligned}
&\delta\langle\hat A\rangle(t)=\Big\langle\frac{\partial\hat A}{\partial F}\Big\rangle F(t)+\int dt'\,\chi^R_{AB}(t-t')\,F(t'),\\
&\chi^R_{AB}(t)=\frac{i}{\hbar}\,\theta(t)\,\big\langle\big[\big(\hat A\big|_{F=0}\big)_{\rm H}(t),\hat B\big]\big\rangle
\end{aligned}
\tag{C3}
$$

where $\hat A_{\rm H}(t)=e^{i\hat H_0t/\hbar}\hat A\,e^{-i\hat H_0t/\hbar}$, $\langle\cdots\rangle$ is the thermal average for $\hat H_0$. The commutator in the second term contains $\hat A\big|_{F=0}$, of zeroth order in $F$.
If $\hat A$ does not contain $F$, the first term is absent.
Because of $\theta(t)$, $\chi^R$ is retarded: the response comes only after the cause $F$ (see the [Appendix](./appendix#causality)).

In a one-body system, with the one-body eigenstates $|n\rangle$ of $\hat H_0$ (energy $E_n$, occupation $f_n$) and one-body operators $\hat A$, $\hat B$, the Fourier transform is

$$
\chi^R_{AB}(\omega)=\int dt\,e^{i\omega t}\chi^R_{AB}(t)
=\sum_{nm}(f_n-f_m)\,\frac{\langle n|\hat A|m\rangle\langle m|\hat B|n\rangle}{E_m-E_n-\hbar\omega-i\delta}
\tag{C4}
$$

(using $\langle[\hat c^\dagger_n\hat c_m,\hat c^\dagger_m\hat c_n]\rangle=f_n-f_m$).

> Note why a commutator appears in Eq. (C3) and elsewhere. This can be understood from the fact that differentiating $\hat A_{\rm H}(t)$ with respect to $t$ gives a commutator, $d\hat A_{\rm H}/dt=(i/\hbar)[\hat H_0,\hat A]_{\rm H}$.
> Time evolution is a unitary rotation in Hilbert space, and an infinitesimal rotation $e^{i\epsilon\hat G}$ changes an operator by $e^{i\epsilon\hat G}\hat A\,e^{-i\epsilon\hat G}\simeq\hat A+i\epsilon[\hat G,\hat A]$.
> When the perturbation $-\hat B F$ is added, an infinitesimal rotation by $\hat B$ is added to the time evolution, so the change of $\hat A$ is the commutator with $\hat B$.
> This is the same as how a $3\times 3$ Hermitian (real symmetric) matrix $M$ transforms under an infinitesimal spatial rotation $R=e^{i\epsilon L}\simeq 1+i\epsilon L$ ($L$ Hermitian): $RMR^{\dagger}\simeq M+i\epsilon[L,M]$.

### Response to an electric field: the Kubo formula {#response-to-field}

A uniform electric field is introduced through a vector potential $A_x(t)$ as $E_x=-\partial A_x/\partial t$ (this avoids the position operator even with periodic boundary conditions).
With the electron charge $-e$, the kinetic energy in the Hamiltonian is $\sum_i(\hat{\bf p}_i+e{\bf A})^2/2m$, and its expansion in ${\bf A}$ is

$$
\begin{aligned}
&\hat H(t)=\hat H_0-\hat J^0_x\,A_x(t)+\frac{e^2}{2m}\hat N A_x(t)^2,\qquad \hat J^0_x=-e\sum_i\frac{\hat p_{x,i}}{m},\\
&\hat J_x\equiv-\frac{\partial\hat H}{\partial A_x}=\hat J^0_x-\frac{e^2}{m}\hat N A_x
\end{aligned}
\tag{C5}
$$

($\hat N$ is the operator of the number of electrons; the sum is over the electrons).
$\hat J_x$ is the total current operator (the sum of the charge $-e$ times the velocity $(\hat{\bf p}_i+e{\bf A})/m$), and $\hat J^0_x$ is $\hat J_x$ at ${\bf A}=0$.

Comparing with Eqs. (C2) and (C3), $\hat A=\hat J_x$, $\hat B=\hat J^0_x$ and $F=A_x$.
Substituting these into Eq. (C3),

$$
\begin{aligned}
&\delta\langle\hat J_x\rangle(t)=\Big\langle\frac{\partial\hat J_x}{\partial A_x}\Big\rangle A_x(t)+\int dt'\,\chi^R_{JJ}(t-t')\,A_x(t'),\\
&\chi^R_{JJ}(t)=\frac{i}{\hbar}\,\theta(t)\,\big\langle\big[\big(\hat J_x\big|_{A_x=0}\big)_{\rm H}(t),\hat J^0_x\big]\big\rangle
=\frac{i}{\hbar}\,\theta(t)\,\big\langle[\hat J^0_{x,{\rm H}}(t),\hat J^0_x]\big\rangle
\end{aligned}
\tag{C6}
$$

(the last equality uses $\hat J_x\big|_{A_x=0}=\hat J^0_x$ of Eq. (C5)).
From Eq. (C5), $\partial\hat J_x/\partial A_x=-(e^2/m)\hat N$, and its expectation value is $-(e^2/m)\langle\hat N\rangle=-(e^2/m)n\Omega$ (the volume of the system is $\Omega$ and the electron density $n$). Therefore

$$
\delta\langle\hat J_x\rangle(t)=-\frac{e^2}{m}n\Omega\,A_x(t)+\int dt'\,\chi^R_{JJ}(t-t')\,A_x(t')
\tag{C7}
$$

The first term is the diamagnetic part and the second the paramagnetic part; both come from the linear-response formula (C3).

For the $e^{-i\omega t}$ component the convolution in the second term becomes the product $\chi^R_{JJ}(\omega)A_x(\omega)$.
The current density is $\delta\langle\hat J_x\rangle/\Omega$, and $A_x(\omega)=E_x(\omega)/i\omega$, so

$$
\sigma(\omega)=\frac{\delta\langle\hat J_x\rangle(\omega)}{\Omega\,E_x(\omega)}=\frac{i}{\omega}\Big[\frac{ne^2}{m}-\frac{1}{\Omega}\,\chi^R_{JJ}(\omega)\Big]
\tag{C8}
$$

This is the result of the Kubo formula. What fixes the dissipation (Joule heat) is the real part of $\sigma$. In the diamagnetic term the phases of the electric field and the current are orthogonal, so it gives no dissipation (it affects only the imaginary part of $\sigma$; multiplying by the imaginary unit $i$ means a phase delay of 90 degrees. This is an inductance effect, not due to a magnetic field but a kinetic inductance from the inertia of the electrons; the energy is stored as the kinetic energy of the electrons).

> In Kusunose's book (H. Kusunose, *Electron Theory of Spin-Orbit-Coupled Physics* (Kodansha, 2019, in Japanese), §1.3, Eqs. (1.106)–(1.109)), the electric field is put in as $-\Omega\,{\bf P}\cdot{\bf E}$ (${\bf P}$ the polarization), and the conductivity is written using partial integration in time (Eqs. (1.107) and (1.108) of the book; there $\Pi_{\mu\nu}=\chi^R_{JJ}$ and $V=\Omega$). Side by side with Eq. (C8),
>
> $$
> \begin{aligned}
> \sigma(\omega)&=\frac{\chi^R_{JJ}(\omega)-ne^2\Omega/m}{i\omega\Omega}&&\text{(Eq. (C8))}\\
> \sigma(\omega)&=\frac{\chi^R_{JJ}(\omega)-\chi^R_{JJ}(0)}{i\omega\Omega}&&\text{(Eq. (1.109) of the book)}
> \end{aligned}
> \tag{C9}
> $$
>
> and the only difference is in the second term. This is because the first term on the right-hand side of Eq. (C6), $-(ne^2\Omega/m)A_x(t)$, is deliberately pushed into the form of the second term, $-\chi^R_{JJ}(\omega=0)\,A_x(t)$.
> $\chi^R_{JJ}(\omega=0)\,A_x$ is the current due to the change of the density matrix $\rho$ when $A_x$ is switched on infinitely slowly (adiabatically). This is seen by retracing the derivation of the linear-response formula. Switching on a static $A_x$ as $e^{\delta t'}$, and setting $F(t')=e^{\delta t'}A_x$ in the formula for $\delta\rho_{\rm I}$ used to derive Eq. (C3), at $t=0$,
>
> $$
> \begin{aligned}
> \delta\rho&=\frac{i}{\hbar}\int_{-\infty}^{0}dt'\,e^{\delta t'}\,[\hat J^0_{x,{\rm H}}(t'),\rho_0]\,A_x,\\
> \langle n|\delta\rho|m\rangle&=(f_m-f_n)\,\frac{\langle n|\hat J^0_x|m\rangle}{E_n-E_m-i\delta}\,A_x,\\
> {\rm Tr}[\delta\rho\,\hat J^0_x]&=\sum_{nm}(f_n-f_m)\,\frac{|\langle n|\hat J^0_x|m\rangle|^2}{E_m-E_n-i\delta}\,A_x=\chi^R_{JJ}(0)\,A_x
> \end{aligned}
> \tag{C10}
> $$
>
> (with $\hbar\delta$ written again as $\delta$). The $i\delta$ is the trace of the slow switching on, and the third line is Eq. (C4) at $\omega=0$.
> That is, the book's form amounts to assuming that under a static $A_x$ the current from the change of $\rho$ exactly cancels the diamagnetic current, so that the total current is 0. In a system where a current remains under a static $A_x$, such as a superconductor, the book's form cannot be used, and Eq. (C8) is used.

### The direct-current limit: the Kubo–Greenwood formula {#kubo-greenwood}

Taking the real part of Eq. (C8), the diamagnetic term drops out and

$$
\begin{aligned}
{\rm Re}\,\sigma(\omega)&=\frac{1}{\omega\Omega}\,{\rm Im}\,\chi^R_{JJ}(\omega)\\
&=\frac{1}{\omega\Omega}\,{\rm Im}\sum_{nm}(f_n-f_m)\,\frac{|\langle n|\hat J^0_x|m\rangle|^2}{E_m-E_n-\hbar\omega-i\delta}\\
&=\frac{\pi}{\omega\Omega}\sum_{nm}(f_n-f_m)\,|\langle n|\hat J^0_x|m\rangle|^2\,\delta(E_m-E_n-\hbar\omega)
\end{aligned}
\tag{C11}
$$

In the second line Eq. (C4) is inserted with $\hat A=\hat B=\hat J^0_x$, and in the third ${\rm Im}\,[1/(x-i\delta)]=\pi\delta(x)$ is used.
For $\omega\to 0$ the delta function gives $E_m=E_n+\hbar\omega$, so $f_n-f_m\simeq -f'(E_n)\,\hbar\omega$ and $\omega$ cancels.
At low temperature $-f'(E)=\delta(E-E_F)$. The one-body matrix element of $\hat J^0_x$ is $-e\langle n|\hat v_x|m\rangle$, and taking the spin into account by multiplying by 2 ($n$, $m$ are orbital states),

$$
\sigma_{xx}=\frac{2\pi\hbar e^2}{\Omega}\sum_{nm}|\langle n|\hat v_x|m\rangle|^2\,\delta(E_F-E_n)\,\delta(E_F-E_m)
\tag{C12}
$$

Here the delta functions are regarded as the limit $\delta\to+0$ of the Lorentzian $\delta_\delta(x)=(1/\pi)\,\delta/(x^2+\delta^2)$ of width $\delta$ (the subscript $\delta$ is the $\delta$ in the denominator of Eq. (C11)). A consideration keeping $\delta$ finite is needed; it is done later in Eqs. (C19) and (C20).
This is fixed only by one-body quantities at the Fermi level $E_F$.
Rewriting it with $\delta(E-H)=-{\rm Im}\,G^R(E)/\pi$ gives the Kubo–Greenwood formula

$$
\sigma_{xx}=\frac{2\hbar e^2}{\pi\Omega}\,
{\rm Tr}\big[\hat v_x\,{\rm Im}\,G^R(E_F)\,\hat v_x\,{\rm Im}\,G^R(E_F)\big]
\tag{C13}
$$

$\hat v_x$ is the velocity operator, ${\rm Tr}$ the sum over one-body (orbital) states, and the factor 2 is the spin.

> For a better understanding, we look at the relation to the viewpoint of many-body perturbation theory. In a one-body system the ring diagram (a product of two $G$, closed with $J$ at the end points) gives the exact result.
> Its linear response function is written with the time-ordered Green function $G^T(E)=\sum_m|m\rangle\langle m|/\big(E-E_m+i\delta\,{\rm sgn}(E_m-\mu)\big)$ as
>
> $$
> \chi^T_{JJ}(\omega)=\frac{i}{2\pi}\int dE\,{\rm Tr}\big[J\,G^T(E)\,J\,G^T(E+\hbar\omega)\big]
> \tag{C14}
> $$
>
> (at zero temperature; $J=-e\hat v_x$ is a one-body operator). Taking the residues in the $E$ integral, only the pairs of occupied $n$ and unoccupied $m$ remain. The poles are at the same places as in Eq. (C4); only the sign of $i\delta$ on the $\omega\lt 0$ side differs.
>
> However, replacing $G^T$ in Eq. (C14) simply by $G^R$ does not make it retarded. $G^T$ has the form of $G^A$ for the occupied states and of $G^R$ for the unoccupied states, so the rewriting into the retarded form is done on the response function $\chi$ (moving all the poles in $\omega$ to the lower half-plane, $\omega\to\omega+i\delta$).

> We rewrite Eq. (C4) with Green functions. Writing $J_{nm}=\langle n|J|m\rangle$ and $f_n=\int dE\,f(E)\,\delta(E-E_n)$, and setting $E=E_n$ in the sum with $f_n$ and $E=E_m$ in the sum with $f_m$,
>
> $$
> \begin{aligned}
> \chi^R_{JJ}(\omega)&=\sum_{nm}(f_n-f_m)\,\frac{J_{nm}J_{mn}}{E_m-E_n-\hbar\omega-i\delta}\\
> &=-\int dE\,f(E)\,{\rm Tr}\big[J\,\delta(E-H)\,J\,G^R(E+\hbar\omega)\\
> &\qquad+J\,G^A(E-\hbar\omega)\,J\,\delta(E-H)\big]
> \end{aligned}
> \tag{C15}
> $$
>
> (using $1/(E_m-E-\hbar\omega-i\delta)=-G^R_m(E+\hbar\omega)$ and $1/(E-E_n-\hbar\omega-i\delta)=G^A_n(E-\hbar\omega)$; the subscripts are components in the eigenstates).
> In the first line of Eq. (C15), the energies with $\hbar\omega=E_m-E_n$ are the excitation energies (excitation from $n$ to $m$). $\chi^R_{JJ}$ has the retarded form with respect to these excitations (the poles $\hbar\omega=E_m-E_n-i\delta$ lie in the lower half-plane of $\omega$). The last step of Eq. (C15) replaces this excitation energy $\hbar\omega=E_m-E_n$ by the quasiparticle energies of the retarded Green function.
> The poles at negative frequency (the pairs with $E_m\lt E_n$) represent the response to a perturbation whose phase rotates in the opposite direction (time still runs forward). This is seen by rewriting the first line of Eq. (C15) in real time:
>
> $$
> \chi^R_{JJ}(t)=\frac{i}{\hbar}\,\theta(t)\sum_{nm}(f_n-f_m)\,J_{nm}J_{mn}\,e^{-i(E_m-E_n)t/\hbar-\delta t/\hbar}
> \tag{C16}
> $$
>
> Each term has the form $\theta(t)\,e^{\mp i|E_m-E_n|t/\hbar-\delta t/\hbar}$, rotating with positive frequency for $E_m\gt E_n$ and with negative frequency for $E_m\lt E_n$, but both carry $\theta(t)$ and respond only for $t\gt 0$. The Fourier transform $\int dt\,e^{i\omega t}$ gives back the first line of Eq. (C15).
> As in Eqs. (C15) and (C16), a response function can always be written as a sum of $1/(\varepsilon-\hbar\omega-i\delta)$ (in the frequency form) or a step function $\times\,e^{-i\varepsilon t/\hbar}$ (in the real-time form), multiplied by matrix elements and occupation factors ($\varepsilon=E_m-E_n$ is the excitation energy).

> Inserting $\delta(E-H)=-\big(G^R(E)-G^A(E)\big)/2\pi i$ into Eq. (C15), at any temperature
>
> $$
> \begin{aligned}
> \chi^R_{JJ}(\omega)&=\frac{1}{2\pi i}\int dE\,f(E)\\
> &\quad\times{\rm Tr}\big[J\,\big(G^R(E)-G^A(E)\big)\,J\,G^R(E+\hbar\omega)\\
> &\qquad+J\,G^A(E-\hbar\omega)\,J\,\big(G^R(E)-G^A(E)\big)\big]
> \end{aligned}
> \tag{C17}
> $$

> The difference between $G^R(E)=(E-H+i\delta)^{-1}$ and $G^A(E)=(E-H-i\delta)^{-1}$ is whether the poles (the eigenvalues of $H$) are placed below or above the real axis. Their difference $G^R-G^A=-2\pi i\,\delta(E-H)$ is the delta function made by those poles.
> Writing $G^R={\rm Re}\,G+i\,{\rm Im}\,G$ and $G^A={\rm Re}\,G-i\,{\rm Im}\,G$ (${\rm Re}\,G$ and ${\rm Im}\,G$ are both Hermitian; ${\rm Im}\,G=-\pi\delta(E-H)$ is the part from the poles only), the imaginary part is a combination of ${\rm Re}\,G\cdot{\rm Im}\,G$ and ${\rm Im}\,G\cdot{\rm Im}\,G$.
> The traces of the terms with ${\rm Re}\,G$ are real and do not contribute to the imaginary part; only the products of ${\rm Im}\,G$ (the poles) remain.

> Shifting $E\to E+\hbar\omega$ in the second term,
>
> $$
> \begin{aligned}
> {\rm Re}\,\sigma(\omega)&=\frac{1}{\omega\Omega}\,{\rm Im}\,\chi^R_{JJ}(\omega)\\
> &=\frac{1}{\pi\omega\Omega}\int dE\,\big[f(E)-f(E+\hbar\omega)\big]\\
> &\quad\times{\rm Tr}\big[J\,{\rm Im}\,G^R(E)\,J\,{\rm Im}\,G^R(E+\hbar\omega)\big]
> \end{aligned}
> \tag{C18}
> $$
>
> This holds for any $\omega$ and temperature (for $\omega\neq 0$ it gives the optical conductivity). In the limit $\omega\to 0$ at low temperature it agrees with Eq. (C12).

The sum (integral) over the product of delta functions in Eq. (C12) may be hard to read, so it is explained here.
Keeping $\delta$ finite, the delta functions in Eq. (C12) remain $\delta_\delta$.
$\sum_n\delta_\delta(E-E_n)$ is the density of states $D(E)$ (the number of states per energy in the whole system).

Let us consider this a little. If the matrix elements $|\langle n|\hat v_x|m\rangle|^2$ are uniform near $E_F$ (average $\overline{|v_x|^2}$),

$$
\begin{aligned}
\sigma_{xx}&\simeq\frac{2\pi\hbar e^2}{\Omega}\,\overline{|v_x|^2}\,\Big[\sum_n\delta_\delta(E_F-E_n)\Big]^2\\
&=\frac{2\pi\hbar e^2}{\Omega}\,\overline{|v_x|^2}\,D(E_F)^2
\end{aligned}
\tag{C19}
$$

independent of $\delta$. When the matrix elements are only diagonal, only the terms $n=m$ remain, and with $\int dx\,\delta_\delta(x)^2=1/(2\pi\delta)$

$$
\begin{aligned}
\sigma_{xx}&\simeq\frac{2\pi\hbar e^2}{\Omega}\sum_n|\langle n|\hat v_x|n\rangle|^2\,\delta_\delta(E_F-E_n)^2\\
&=\frac{\hbar e^2}{\Omega}\,\overline{v_x^2}\,\frac{D(E_F)}{\delta}
\end{aligned}
\tag{C20}
$$

which diverges for $\delta\to 0$ (without scattering the current does not decay). Roughly, $\delta$ may be regarded as the width $\hbar/2\tau$ from the lifetime $\tau$ of quasiparticles due to phonons and so on.
In practice the matrix elements are not only diagonal, but the number of states $m$ connected to $n$ is finite (by completeness $\sum_m|\langle n|\hat v_x|m\rangle|^2=\langle n|\hat v_x^2|n\rangle$ is finite), so Eq. (C20) may be regarded as correct.

Eq. (C20)'s "the electrons on the Fermi surface contribute to the conduction" is connected to Drude's "all electrons are accelerated for the time $\tau$" as follows.
The following factor on the right-hand side of Eq. (C20), for one band and per spin, becomes by integration by parts in $k_x$, using $v_x\,\delta(E_F-\varepsilon_{\bf k})=-\hbar^{-1}\partial_{k_x}\theta(E_F-\varepsilon_{\bf k})$ ($\theta$ is the step function of the occupation, $v_x=\hbar^{-1}\partial\varepsilon_{\bf k}/\partial k_x$; the boundary term vanishes by the periodicity of the Brillouin zone),

$$
\begin{aligned}
\frac{\overline{v_x^2}\,D(E_F)}{\Omega}&=\int\frac{d^3k}{(2\pi)^3}\,v_x^2\,\delta(E_F-\varepsilon_{\bf k})\\
&=\frac{1}{\hbar}\int\frac{d^3k}{(2\pi)^3}\,\theta(E_F-\varepsilon_{\bf k})\,\frac{\partial v_x}{\partial k_x}\\
&=\int_{\rm occ}\frac{d^3k}{(2\pi)^3}\,\frac{1}{\hbar^2}\frac{\partial^2\varepsilon_{\bf k}}{\partial k_x^2}
\end{aligned}
\tag{C21}
$$

One $v_x$ of $v_x^2$ combines with the delta function into the derivative of $\theta$ (differentiating $\theta$ in $k_x$ brings out $\partial\varepsilon_{\bf k}/\partial k_x=\hbar v_x$), and the integration by parts moves the derivative onto the other $v_x$, giving $\partial v_x/\partial k_x=\hbar^{-1}\partial^2\varepsilon_{\bf k}/\partial k_x^2$ (the inverse effective mass).
The left-hand side is a Fermi-surface quantity, and the right-hand side is the sum of the inverse effective masses over all occupied states.
For free electrons $\partial^2\varepsilon_{\bf k}/\partial k_x^2=\hbar^2/m$ is constant, so the right-hand side is $n/2m$ ($n$ is the electron density including both spins, and $n/2$ is per spin).
Therefore, with $\delta=\hbar/2\tau$ in Eq. (C20) and using Eq. (C21),

$$
\sigma_{xx}=2e^2\tau\,\frac{\overline{v_x^2}\,D(E_F)}{\Omega}=\frac{ne^2\tau}{m}
\tag{C22}
$$

which is [the Drude formula](#drude).
This approximation of a constant effective mass on the Fermi surface should hold well for Na, Cu, Al and so on.

> The structure of Eq. (C21), an integral over the occupied region connected to an integral over the Fermi surface by integration by parts, is the same as in the calculation of the anomalous Hall conductivity with the Berry curvature (the integral of the Berry curvature over the occupied states is rewritten as an integral of the Berry phase on the Fermi surface).


## Average over the impurity configurations {#impurity-average}

Eq. (C13) is for one fixed configuration of the impurities. Writing ${\rm Im}\,G^R=(G^R-G^A)/2i$, and since $G^A=(G^R)^\dagger$ makes ${\rm Tr}[\hat v_xG^R\hat v_xG^A]$ real and ${\rm Tr}[\hat v_xG^A\hat v_xG^A]$ the complex conjugate of ${\rm Tr}[\hat v_xG^R\hat v_xG^R]$, the first line of Eq. (C23) holds for each configuration (no approximation is used). The observed conductivity is the average of Eq. (C13) over the impurity configurations, which gives the second line (all $G$ at $E=E_F$):

$$
\begin{aligned}
{\rm Tr}\big[\hat v_x\,{\rm Im}\,G^R\,\hat v_x\,{\rm Im}\,G^R\big]&=\frac{1}{2}\Big\{ {\rm Tr}\big[\hat v_x\,G^R\,\hat v_x\,G^A\big]\\
&\qquad-{\rm Re}\,{\rm Tr}\big[\hat v_x\,G^R\,\hat v_x\,G^R\big]\Big\},\\
\sigma_{xx}&=\frac{\hbar e^2}{\pi\Omega}\Big\{\Big\langle{\rm Tr}\big[\hat v_x\,G^R\,\hat v_x\,G^A\big]\Big\rangle_{\rm imp}\\
&\qquad-{\rm Re}\,\Big\langle{\rm Tr}\big[\hat v_x\,G^R\,\hat v_x\,G^R\big]\Big\rangle_{\rm imp}\Big\}
\end{aligned}
\tag{C23}
$$

The two terms are expanded together below (that the second term is small is seen in the section on the [Drude formula](#drude)).
To compute the average of the product of Green functions, we go back to the original one-body problem. $G^R(E)=(E-H_0-V+i\delta)^{-1}$ is the Green function of the one-body Hamiltonian $H=H_0+V$ for one fixed configuration of the impurities ($H_0$ is the kinetic energy and $V$ is Eq. (C1)). It is expanded in $V$, and the average of the product is decomposed into products of averages of pairs of $V$.

::: tip Preparation: the average when the impurities are placed at random
**Before the average** (for one fixed configuration): the matrix element between plane waves $|{\bf k}\rangle$ is $\langle{\bf k}|V|{\bf k}'\rangle=\Omega^{-1}u({\bf k}-{\bf k}')\sum_j e^{-i({\bf k}-{\bf k}')\cdot{\bf R}_j}$.

**After the average**: since the impurity positions are random, the phases $e^{-i({\bf k}-{\bf k}')\cdot{\bf R}_j}$ cancel on averaging (random phases), the first-order average vanishes for ${\bf k}\neq{\bf k}'$, and at second order only the pairs on the same impurity $j$, whose phases cancel, remain:

$$
\begin{aligned}
&\big\langle \langle{\bf k}|V|{\bf k}'\rangle\big\rangle_{\rm imp}=n_i\,u(0)\,\delta_{ {\bf k}{\bf k}'},\\
&\big\langle \langle{\bf k}|V|{\bf k}'\rangle\langle{\bf k}'|V|{\bf k}\rangle\big\rangle_{\rm imp}=\frac{n_i}{\Omega}\,|u({\bf k}-{\bf k}')|^2
\end{aligned}
\tag{C24}
$$

:::

With this we expand the two terms of Eq. (C23). Writing $G^{\pm}=(E-H_0-V\pm i\delta)^{-1}$ ($G^+=G^R$, $G^-=G^A$), we consider ${\rm Tr}[\hat v_xG^{\pm}\hat v_xG^{\pm'}]$ (the two signs are taken independently; $(+,-)$ gives the first term, $(+,+)$ the second). In plane waves $|{\bf k}\rangle$ the operator $\hat v_x$ is diagonal ($v_x({\bf k})=\hbar k_x/m$), so the trace is

$$
\begin{aligned}
{\rm Tr}\big[\hat v_x G^{\pm}\hat v_x G^{\pm'}\big]
&=\sum_{ {\bf k}{\bf k}'}v_x({\bf k})\,G^{\pm}_{ {\bf k}{\bf k}'}\,v_x({\bf k}')\,G^{\pm'}_{ {\bf k}'{\bf k}}
\end{aligned}
\tag{C25}
$$

($G_{ {\bf k}{\bf k}'}=\langle{\bf k}|G|{\bf k}'\rangle$; $E=E_F$ is not written). Expanding $G^{\pm}$ in $V$ gives

$$
\begin{aligned}
G^{\pm}_{ {\bf k}{\bf k}'}={}&G^{\pm}_0({\bf k})\,\delta_{ {\bf k}{\bf k}'}+G^{\pm}_0({\bf k})\,V_{ {\bf k}{\bf k}'}\,G^{\pm}_0({\bf k}')\\
&+\sum_{ {\bf k}_1}G^{\pm}_0({\bf k})\,V_{ {\bf k}{\bf k}_1}\,G^{\pm}_0({\bf k}_1)\,V_{ {\bf k}_1{\bf k}'}\,G^{\pm}_0({\bf k}')+\cdots
\end{aligned}
\tag{C26}
$$

where $G^{\pm}_0({\bf k})=1/(E-\varepsilon_{\bf k}\pm i\delta)$, $\varepsilon_{\bf k}=\hbar^2k^2/2m$ and $V_{ {\bf k}{\bf k}'}=\langle{\bf k}|V|{\bf k}'\rangle$.
We put this into Eq. (C25) and average. Of the terms first order in $V$ only the component ${\bf k}={\bf k}'$, $n_iu(0)$, remains; it is a uniform energy shift and is absorbed in $E_F$. The second-order terms fall into three groups according to where the two $V$ are, and with Eq. (C24)

$$
\begin{aligned}
&\Big\langle{\rm Tr}\big[\hat v_x G^{\pm}\hat v_x G^{\pm'}\big]\Big\rangle_{\rm imp}\\
&\quad=\sum_{\bf k}v_x({\bf k})^2\,G^{\pm}_0({\bf k})\,G^{\pm'}_0({\bf k})\\
&\qquad+\color{#1f6fd0}{\sum_{\bf k}v_x({\bf k})^2\,G^{\pm}_0({\bf k})\,\Sigma^{\pm}_2({\bf k})\,G^{\pm}_0({\bf k})\,G^{\pm'}_0({\bf k})}\\
&\qquad+\color{#1f6fd0}{\sum_{\bf k}v_x({\bf k})^2\,G^{\pm}_0({\bf k})\,G^{\pm'}_0({\bf k})\,\Sigma^{\pm'}_2({\bf k})\,G^{\pm'}_0({\bf k})}\\
&\qquad+\color{#d26a00}{\sum_{ {\bf k}{\bf k}'}v_x({\bf k})\,G^{\pm}_0({\bf k})\,G^{\pm'}_0({\bf k})}\\
&\qquad\quad\color{#d26a00}{\times\frac{n_i}{\Omega}\,|u({\bf k}-{\bf k}')|^2\,G^{\pm}_0({\bf k}')\,G^{\pm'}_0({\bf k}')\,v_x({\bf k}')}\\
&\qquad+O(V^4)
\end{aligned}
\tag{C27}
$$

where $\Sigma^{\pm}_2({\bf k})=\sum_{ {\bf k}_1}(n_i/\Omega)\,|u({\bf k}-{\bf k}_1)|^2\,G^{\pm}_0({\bf k}_1)$.
The second and third terms of Eq. (C27) (<span style="color:#1f6fd0">blue</span>) pair the two $V$ within one $G$ and have one $\Sigma_2$ inserted (the average within one $G$ remains only for ${\bf k}={\bf k}'$ because of the random phases, so it is diagonal in ${\bf k}$). At higher orders in $V$, two, three, ... $\Sigma_2$ are inserted. Keeping only the pairings that do not cross each other (Born approximation), they form a geometric series,

$$
\begin{aligned}
\bar G^{\pm}({\bf k},E)&=G^{\pm}_0+G^{\pm}_0\Sigma^{\pm}G^{\pm}_0+G^{\pm}_0\Sigma^{\pm}G^{\pm}_0\Sigma^{\pm}G^{\pm}_0+\cdots\\
&=\frac{1}{E-\varepsilon_{\bf k}-\Sigma^{\pm}({\bf k},E)}
\end{aligned}
\tag{C28}
$$

($\Sigma^{\pm}=\Sigma^{\pm}_2$). $\bar G$ is translationally invariant and diagonal in ${\bf k}$.
The fourth term (<span style="color:#d26a00">orange</span>) pairs a $V$ in $G^{\pm}$ with a $V$ in $G^{\pm'}$; it is the correlation of scattering by the same impurity, which remains because the product of the two $G$ is averaged (it does not appear if each $G$ is averaged first and then multiplied). Dropping the terms of this type, Eq. (C23) becomes

$$
\begin{aligned}
\sigma_{xx}=\frac{\hbar e^2}{\pi\Omega}\sum_{\bf k}v_x({\bf k})^2\Big\{&\bar G^R({\bf k},E_F)\,\bar G^A({\bf k},E_F)\\
&-{\rm Re}\,\bar G^R({\bf k},E_F)^2\Big\}
\end{aligned}
\tag{C29}
$$

Below, $\Sigma$ is obtained (next section) and the two terms of Eq. (C29) are evaluated (section on the [Drude formula](#drude)). The terms of the type of the fourth term are included in the section on the [vertex correction](#vertex-correction). Pairings that cross each other give the weak-localization correction, which is not treated here.

## Summing the second and third terms: the lifetime and the Drude formula {#self-energy-sum}

### The first term on the right-hand side of Eq. (C29): the lifetime and the Drude formula {#averaged-g}

We now obtain the $\Sigma^{\pm}$ of Eq. (C28) explicitly. Turning the sum in $\Sigma_2$ into an integral, the self-energy is

$$
\Sigma^R({\bf k},E)=n_i\int\frac{d^3k'}{(2\pi)^3}\,
\frac{|u({\bf k}-{\bf k}')|^2}{E-\varepsilon_{\bf k'}+i\delta}
\tag{C30}
$$

and its imaginary part gives a time $\tau$ (the real part is an energy shift, absorbed in $E_F$):

$$
{\rm Im}\,\Sigma^R=-\frac{\hbar}{2\tau},
\tag{C31}
$$

Eq. (C28) is the expression of the Green function that includes the impurity scattering fully as a one-body problem (within the Born approximation). Its poles are at a distance $\hbar/2\tau$ from the real axis, and $G$ decays in time as $e^{-t/2\tau}$ (the probability $|G|^2$ as $e^{-t/\tau}$). So this $\tau$ can be regarded as the quasiparticle lifetime.
At $E=E_F$

$$
\frac{1}{\tau}=\frac{2\pi}{\hbar}\,n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\,\delta(E_F-\varepsilon_{\bf k'})
\tag{C32}
$$

$1/\tau$ is just the scattering rate by Fermi's golden rule (see the correspondence between the golden rule and the imaginary part of the self-energy in [Basic I](./basic#td-vs-ti-perturbation)).

Putting this into Eq. (C29), $\bar G^R({\bf k},E_F)=1/(E_F-\varepsilon_{\bf k}+i\hbar/2\tau)$ and $\bar G^A$ is its complex conjugate. In the first term, $\bar G^R\bar G^A=1/\big[(E_F-\varepsilon_{\bf k})^2+(\hbar/2\tau)^2\big]$ has a sharp peak of width $\hbar/2\tau$ at $\varepsilon_{\bf k}=E_F$, and the sum over ${\bf k}$ is proportional to $\tau$. The second term $(\bar G^R)^2$ does not behave so, because its poles are on the same side of the real axis. They are evaluated explicitly below.

<a id="drude"></a>

In Eq. (C29), $\bar G$ contains the width $\hbar/2\tau$ from the self-energy.
When $\tau$ is long ($\hbar/\tau\ll E_F$),

$$
\begin{aligned}
\bar G^R({\bf k},E)\,\bar G^A({\bf k},E)
&=\frac{1}{(E-\varepsilon_{\bf k})^2+(\hbar/2\tau)^2}\\
&\simeq\frac{2\pi\tau}{\hbar}\,\delta(E-\varepsilon_{\bf k})
\end{aligned}
\tag{C33}
$$

which corresponds to the Lorentzian of Eq. (C20) with the width $\delta$ set to $\hbar/2\tau$. So the first term of Eq. (C29) becomes

$$
\sigma_{xx}=2e^2\tau\,\frac{\overline{v_x^2}\,D(E_F)}{\Omega}
\tag{C34}
$$

the same Drude formula as Eq. (C22) ($ne^2\tau/m$ for free electrons).

### The second term on the right-hand side of Eq. (C29) {#rr-term}

We estimate the second term of Eq. (C29) precisely. Writing $g(\varepsilon)=\sum_{\bf k}v_x({\bf k})^2\,\delta(\varepsilon-\varepsilon_{\bf k})$ and $\gamma=\hbar/2\tau$,

$$
\begin{aligned}
\sum_{\bf k}v_x({\bf k})^2\,\bar G^R({\bf k},E_F)^2
&=\int d\varepsilon\,\frac{g(\varepsilon)}{(E_F-\varepsilon+i\gamma)^2}\\
&=-\int d\varepsilon\,\frac{g'(\varepsilon)}{E_F-\varepsilon+i\gamma}
\end{aligned}
\tag{C35}
$$

The second equality is an integration by parts with $(E_F-\varepsilon+i\gamma)^{-2}=\partial_\varepsilon(E_F-\varepsilon+i\gamma)^{-1}$ (the boundary terms vanish). The right-hand side stays finite for $\gamma\to0$, and its real part is $-{\rm P}\!\int d\varepsilon\,g'(\varepsilon)/(E_F-\varepsilon)$, independent of $\tau$. Since $g$ varies on an energy scale of about $E_F$, its size is of the order of $g(E_F)$.
The first term, on the other hand, written with $g$ using Eq. (C33), is

$$
\begin{aligned}
\sum_{\bf k}v_x({\bf k})^2\,\bar G^R({\bf k},E_F)\,\bar G^A({\bf k},E_F)&=\int d\varepsilon\,\frac{g(\varepsilon)}{(E_F-\varepsilon)^2+\gamma^2}\\
&\simeq\frac{\pi}{\gamma}\,g(E_F)=\frac{2\pi\tau}{\hbar}\,g(E_F)
\end{aligned}
\tag{C36}
$$

proportional to $\tau$. The ratio of Eqs. (C35) and (C36) is of the order of $\hbar/(E_F\tau)$, small for weak scattering ($\hbar/\tau\ll E_F$), so the second term was dropped.

### Summary {#lifetime-summary}

In summary, in the approximation dropping the terms of the type of the fourth term (orange), the first term of Eq. (C29) gives the Drude formula (C34) with the lifetime $\tau$, and the second term is smaller by a factor of the order of $\hbar/(E_F\tau)$.
A naive $\delta\to 0$ gave zero resistance, but the scattering width due to the impurities plays the role of $\delta$, and the impurities give a resistance through the lifetime $\tau$. That is, we obtain the result of the naive Drude picture, in which electrons are accelerated by the field and their lifetime is set by the time $\tau$ until they are scattered by the impurities.

## Vertex correction: the transport relaxation time {#vertex-correction}

The product average is taken by summing the terms of the type of the fourth term of Eq. (C27) (orange) to all orders (the ladder; the approximation that goes with the Born self-energy (C30)).
The current vertex $\Gamma_x({\bf k})$ satisfies

$$
\begin{aligned}
\Gamma_x({\bf k})={}&v_x({\bf k})+n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\\
&\times\bar G^R({\bf k}',E_F)\bar G^A({\bf k}',E_F)\,\Gamma_x({\bf k}')
\end{aligned}
\tag{C37}
$$

Putting $\Gamma_x({\bf k})=\gamma\,v_x({\bf k})$ and using (C33), the average over the direction of ${\bf k}'$ turns $v_x({\bf k}')$ into $v_x({\bf k})\cos\theta_{\bf kk'}$, so that
$\gamma=1+\gamma\,(1-\tau/\tau_{\rm tr})$, that is $\gamma=\tau_{\rm tr}/\tau$. Here

$$
\frac{1}{\tau_{\rm tr}}=\frac{2\pi}{\hbar}\,n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\,
\big(1-\cos\theta_{\bf kk'}\big)\,\delta(E_F-\varepsilon_{\bf k'})
\tag{C38}
$$

is the transport relaxation time. Replacing one $v_x$ in Eq. (C34) by $\Gamma_x$ gives

$$
\sigma_{xx}=2e^2\tau_{\rm tr}\,\frac{\overline{v_x^2}\,D(E_F)}{\Omega}
\tag{C39}
$$

($ne^2\tau_{\rm tr}/m$ for free electrons).

## Two relaxation times: the lifetime $\tau$ and the transport relaxation time $\tau_{\rm tr}$ {#two-taus}

What is important is the difference between the lifetime $\tau$ (Eq. (C32)) and the transport relaxation time $\tau_{\rm tr}$ (Eq. (C38)). Side by side, <strong>the only difference is the weight $1-\cos\theta_{\bf kk'}$ in the integrand.</strong>

| | quasiparticle lifetime $\tau$ (Eq. (C32)) | transport relaxation time $\tau_{\rm tr}$ (Eq. (C38)) |
| --- | --- | --- |
| other names | single-particle relaxation time, quantum lifetime | momentum relaxation time, transport lifetime |
| weight of the scattering | $1$ | $1-\cos\theta_{\bf kk'}$ |
| where it comes from | the average of a single $G$ (self-energy; second and third terms of Eq. (C27)) | the average of the product of $G$ (sum of terms of the type of the fourth term of Eq. (C27)) |
| time of what | the time an electron stays in the state ${\bf k}$; scattering in any direction loses ${\bf k}$ | the time in which the current (the direction of the momentum) is lost |
| where it appears | the width $\hbar/\tau$ of the spectral function | the conductivity $\sigma_{xx}$ (Eq. (C39)) |
| in the theory of scattering cross sections | total cross section $\sigma_{\rm tot}$ | transport (momentum-transfer) cross section $\sigma_{\rm tr}$ |

In terms of the theory of scattering cross sections: with the differential cross section $d\sigma/d\Omega$ of one impurity (in the Born approximation $d\sigma/d\Omega=(m/2\pi\hbar^2)^2|u({\bf q})|^2$, $|{\bf q}|=2k_F\sin(\theta/2)$), Eqs. (C32) and (C38) are written as

$$
\begin{aligned}
&\frac{1}{\tau}=n_i\,v_F\,\sigma_{\rm tot},\qquad \sigma_{\rm tot}=\int d\Omega\,\frac{d\sigma}{d\Omega},\\
&\frac{1}{\tau_{\rm tr}}=n_i\,v_F\,\sigma_{\rm tr},\qquad \sigma_{\rm tr}=\int d\Omega\,(1-\cos\theta)\,\frac{d\sigma}{d\Omega}
\end{aligned}
\tag{C40}
$$

The total cross section $\sigma_{\rm tot}$ counts scattering in any direction and corresponds to the rate at which particles leave a beam. The transport cross section $\sigma_{\rm tr}$ is weighted by the momentum transfer $\hbar k_F(1-\cos\theta)$, and it is the one that appears in the resistance and diffusion. $\tau$ appears in the width of the spectral function (for example the lifetime obtained from the Dingle temperature of de Haas–van Alphen oscillations, also called the quantum lifetime). Beyond the Born approximation, $|u|^2$ is replaced by the T matrix of one impurity, and the correspondence holds as it is.

- **Forward scattering** ($\theta\simeq0$) moves the electron to another state ${\bf k}'$, so it contributes to $1/\tau$, but it hardly changes the direction of the momentum, so it hardly contributes to $1/\tau_{\rm tr}$.
- **Backward scattering** ($\theta\simeq\pi$) reverses the momentum, so it contributes to $1/\tau_{\rm tr}$ with twice the weight.
- Hence, from $1-\cos\theta\le2$, $\tau_{\rm tr}\ge\tau/2$, and the stronger the forward scattering, the longer $\tau_{\rm tr}$ is than $\tau$.
- If $u$ is point-like ($u({\bf q})$ independent of ${\bf q}$), the scattering does not depend on the direction, the average of $\cos\theta$ is zero, and $\tau_{\rm tr}=\tau$ (the vertex correction vanishes). The more $u({\bf q})$ is concentrated at small ${\bf q}$, as for screened Coulomb impurities, the stronger the forward scattering, and $\tau_{\rm tr}\gg\tau$ is possible.

If each $G$ is averaged first and then multiplied, as in Eq. (C29), $\sigma_{xx}$ contains $\tau$ (Eq. (C34)). To obtain the correct $\tau_{\rm tr}$, the average of the product of $G$ (the vertex correction) is needed. This is what it means that the average is taken of the product of $G$.
Eq. (C39) with $\tau_{\rm tr}$ is the same as the result of the Boltzmann equation.

<hr style="margin-top:4em;margin-bottom:2.5em">

This concludes the discussion of impurities. The following are notes for calculations such as the anomalous Hall effect.

## Note: integrals over the Brillouin zone and its topology $T^3$ {#bz-topology}

Linear response computes integrals over the Brillouin zone. We therefore review the geometry (topology) of the Brillouin zone. Once it is understood on the two-dimensional torus, it can be grasped intuitively to a fair extent. This is the key point in discussing topological physical quantities; jumping straight into the calculation makes things hard to follow.

${\bf k}$ and ${\bf k}+{\bf G}$ (${\bf G}$ a reciprocal lattice vector) represent the same point, so the Brillouin zone, with its opposite faces glued together, is a three-dimensional torus $T^3$.
$T^2$ (the two-dimensional torus) is a two-dimensional space like the surface of a doughnut. $T^3$ (the three-dimensional torus) is the three-dimensional version of $T^2$, a kind of closed three-dimensional space without boundary.
On $T^2$ there are two kinds of loops ($T^1$) that cannot be shrunk to a point (loops that can be deformed continuously into each other are regarded as the same): calling the large circle of the doughnut the "ring" and the small circle of its cross-section the "tube", one goes along the ring and one goes around the tube. Other loops are combinations of these two, for example a path that goes twice around the tube and then once along the ring.
On $T^3$ there are three such kinds of loops ($T^1$), winding in the $k_x$, $k_y$, $k_z$ directions, and in addition three kinds of closed surfaces ($T^2$) that cannot be shrunk to a point. These are closed surfaces without boundary (two-dimensional tori), specified by $k_x=$ const and so on.
For example, the cross-section $T^2(k_z)$ at $k_z=$ const wraps the Brillouin zone in the $k_x$ and $k_y$ directions. An integral over such a surface cannot be rewritten as an integral over a boundary.
A region that fits inside the Brillouin zone (like the pocket $S$ below), on the other hand, has a boundary, and an integral over it can be rewritten by Stokes' theorem as an integral over the boundary.

Why do we need to think about this at all? Because in linear-response calculations, problems can arise (as a topological response). The linear-response integrals, in the form of sums over states, do not depend on how the phases are chosen. When they are rewritten, by partial integration, as integrals of the Berry curvature or over the Fermi surface, however, one uses the fact that the phase can be chosen smoothly. On a closed surface $T^2$ in $T^3$, this may not be possible.

Thinking of it in simplified form, it suffices to focus on the case in which a phase factor $e^{i\theta(x,y,z)}$ is specified at each point $(x,y,z)$ of $T^3$. It may also be given in the form amplitude × phase, but the amplitude is not essential and can be considered later.

The phase is assumed to vary continuously and smoothly. However, singular lines (vortex lines; the limit of the phase depends on the direction of approach, so the phase is undefined there, and it goes around by an integer multiple of $2\pi$ on a loop around them; for example $\theta=\arg(x+iy)$ on the $z$ axis) and singular points where bands are degenerate are also considered in some cases. In three dimensions, a singularity of the phase is generically a line.

The number of times $e^{i\theta}$ goes around the unit circle along a loop $L$ is the winding number

$$
w_L=\frac{1}{2\pi}\oint_L\nabla\theta\cdot d{\bf r}
\tag{C41}
$$

(an integer). It does not change under continuous deformation, except when a singular line is crossed. So a loop that can be shrunk to a point has winding number 0. In ordinary three-dimensional space, a nonzero winding number means that the loop winds around a singular line.
On $T^3$, however, the loops ($T^1$) winding in the $x$, $y$, $z$ directions cannot be shrunk, so winding numbers are possible without any singular line. For example, in a periodic box with edges $a$, $b$, $c$,

$$
e^{i\theta}=e^{i{\bf G}\cdot{\bf r}},\qquad {\bf G}=2\pi\Big(\frac{w_x}{a},\frac{w_y}{b},\frac{w_z}{c}\Big)
\tag{C42}
$$

Advancing $x$ by $a$ advances the phase by $2\pi w_x$, and $e^{i\theta}$ comes back to itself (likewise for $y$ and $z$). So, as a function on the periodic box ($T^3$), it is smooth wherever one looks. Even so, going once around in the $x$ direction turns the phase by $2\pi w_x$, so it has the winding numbers $(w_x,w_y,w_z)$. ${\bf G}$ is a reciprocal lattice vector, and the winding numbers are its components (in units of $2\pi/a$, $2\pi/b$, $2\pi/c$). This is the essence of topological quantization. Maps $e^{i\theta(x,y,z)}$ from $T^3$ to the unit circle are classified by these winding numbers; only maps with the same winding numbers can be connected by continuous deformation.
For wave functions (reading $(x,y,z)$ as ${\bf k}$), the winding number of a gauge transformation $e^{i\phi({\bf k})}$ corresponds to this (the $w$ of Eq. (C47) below). Moreover, the phase may not be choosable as a single smooth function over a whole closed $T^2$; the winding number of the mismatch is the Chern number (Eqs. (C44) and (C45) below). With a Fermi surface, the topology of the Fermi surface also matters (for example, a closed pocket Fermi surface cannot carry a Chern number unless a band degeneracy is inside), but we do not go into it here.

## Note: Fermi-surface and Fermi-sea terms {#fermi-sea-surface}

The static conductivity $\sigma_{ij}$ of the Kubo formula is sometimes split into a **Fermi-surface term**, with the weight $-\partial f/\partial E$, and a **Fermi-sea term**, integrated over all occupied states with the weight $f(E)$ (Středa's split; often used in the literature on the anomalous Hall effect).
The Fermi-sea term is antisymmetric in $i\leftrightarrow j$, so it does not contribute to $\sigma_{xx}$. The first term of Eq. (C23) is the $G^RG^A$ part of the Fermi-surface term.

This split, however, is not unique, because integration by parts rewrites one into the other.
In Eq. (C21), the integral over the occupied region (Drude's $n/m$) and $\overline{v_x^2}D(E_F)$ on the Fermi surface were the same quantity.
For the anomalous Hall conductivity too, the integral of the Berry curvature over the occupied states can be rewritten as the Berry phase on the Fermi surface (F. D. M. Haldane, Phys. Rev. Lett. **93**, 206602 (2004)).
So the distinction "only the electrons on the Fermi surface contribute, or all electrons" has no physical meaning.

The exception is the integer part. We see this with formulas in the following setting.
In a three-dimensional crystal (orthorhombic for simplicity, with lattice constants $a$, $b$, $c$), band 1 is completely filled (no Fermi surface), and band 2 is occupied only in a region $S$ of the Brillouin zone (an electron pocket). The Fermi surface is the closed surface $\partial S$, the boundary of $S$.
The Berry connection and Berry curvature of the periodic part $|u_{n{\bf k}}\rangle$ of the Bloch functions are

$$
\begin{aligned}
&{\bf a}_n({\bf k})=i\langle u_{n{\bf k}}|\nabla_{\bf k}u_{n{\bf k}}\rangle,\qquad
{\bf\Omega}_n({\bf k})=\nabla_{\bf k}\times{\bf a}_n({\bf k}),\\
&\Omega_n^z=-2\,{\rm Im}\sum_{m\neq n}\frac{\langle u_{n{\bf k}}|\partial_{k_x}H_{\bf k}|u_{m{\bf k}}\rangle\langle u_{m{\bf k}}|\partial_{k_y}H_{\bf k}|u_{n{\bf k}}\rangle}{(E_{n{\bf k}}-E_{m{\bf k}})^2}
\end{aligned}
\tag{C43}
$$

($H_{\bf k}=e^{-i{\bf k}\cdot{\bf r}}He^{i{\bf k}\cdot{\bf r}}$). The second line is the form as a sum over states, which follows directly from the Kubo formula.
The intrinsic anomalous Hall conductivity (the sign depends on the orientation convention of ${\bf\Omega}$), split into cross-sections at each $k_z$, is

$$
\begin{aligned}
\sigma_{xy}&=\frac{e^2}{\hbar}\sum_n\int_{\rm occ}\frac{d^3k}{(2\pi)^3}\,\Omega_n^z\\
&=\frac{e^2}{h}\,\frac{1}{2\pi}\int dk_z\,\Big[\frac{1}{2\pi}\int_{T^2(k_z)}\Omega_1^z\,d^2k+\frac{1}{2\pi}\int_{S(k_z)}\Omega_2^z\,d^2k\Big]
\end{aligned}
\tag{C44}
$$

where $T^2(k_z)$ is the cross-section of the Brillouin zone at fixed $k_z$ and $S(k_z)$ the cross-section of the pocket in that plane.

The integrand $\Omega_n^z$ does not depend on the choice of phases, in either term. Changing the phase as $|u_{n{\bf k}}\rangle\to e^{i\phi({\bf k})}|u_{n{\bf k}}\rangle$ changes the Berry connection as ${\bf a}_n\to{\bf a}_n-\nabla_{\bf k}\phi$, but $\nabla_{\bf k}\times\nabla_{\bf k}\phi=0$. In the second line of Eq. (C43), too, the phase of each state cancels in the numerator. Only the Berry connection ${\bf a}_n$ depends on the choice of phases.
When bands are degenerate, there are choices that mix the degenerate bands, so $\Omega_n$ is not determined band by band. The sum over the occupied bands, $\sum_n\Omega_n$, is determined.

**First term (band 1)**: The integral on the closed surface $T^2(k_z)$, $C_1^z=\frac{1}{2\pi}\int_{T^2(k_z)}\Omega_1^z\,d^2k$, is a Chern number, an integer. Band 1 has a gap, so it does not depend on $k_z$. Since $T^2(k_z)$ is closed, this cannot be rewritten as an integral over a boundary.
The $k_z$ integral gives $2\pi/c$, so the first term is $(e^2/h)\,C_1^z/c$. Together with the Chern numbers $C_1^x$, $C_1^y$ of the surfaces wrapping the Brillouin zone in the other directions,

$$
\sigma_{ij}^{(1)}=\frac{e^2}{2\pi h}\,\epsilon_{ijl}\,G_l,\qquad
{\bf G}=2\pi\Big(\frac{C_1^x}{a},\frac{C_1^y}{b},\frac{C_1^z}{c}\Big)
\tag{C45}
$$

where ${\bf G}$ is a reciprocal lattice vector (B. I. Halperin, Jpn. J. Appl. Phys. **26**, Suppl. 26-3, 1913 (1987)).

**Second term (band 2)**: In each cross-section, with a gauge smooth on $S(k_z)$, Stokes' theorem gives

$$
\int_{S(k_z)}\Omega_2^z\,d^2k=\oint_{\partial S(k_z)}{\bf a}_2\cdot d{\bf k}\equiv\gamma_2(k_z)
\tag{C46}
$$

The loop $\partial S(k_z)$ is the Fermi surface cut at $k_z$ (the cyclotron orbit in a magnetic field along $z$), and $\gamma_2(k_z)$ is the Berry phase on it.
From the states on the Fermi surface alone, however, $\gamma_2(k_z)$ is fixed only modulo $2\pi$. A gauge transformation $|u_{2{\bf k}}\rangle\to e^{i\phi({\bf k})}|u_{2{\bf k}}\rangle$ gives ${\bf a}_2\to{\bf a}_2-\nabla_{\bf k}\phi$, so that

$$
\gamma_2(k_z)\to\gamma_2(k_z)-2\pi w
\tag{C47}
$$

($w$ is the integer number of times $\phi$ winds by $2\pi$ around the loop). Which $w$ is right (which gauge connects smoothly inside $S(k_z)$) cannot be decided without looking at the states inside the Fermi surface.

In summary, the information on the Fermi surface alone (the Fermi surface and the Berry phases on it) fixes $\sigma_{ij}$ only up to a term of the form (C45) (a reciprocal lattice vector ${\bf G}$). That integer part is information of the Fermi sea, such as band 1, which has no Fermi surface (Haldane's conclusion).
As an example, consider two systems A and B in which band 2 (the Fermi surface and the states on it) is exactly the same and only band 1 differs. For instance, under a metal, without coupling, put either an ordinary insulator (A, ${\bf G}=0$) or Chern-insulator layers stacked along $z$ with period $c$ (B, $C_1^z=1$). Eq. (C44) gives

$$
\sigma_{xy}^{\rm B}-\sigma_{xy}^{\rm A}=\frac{e^2}{hc}
\tag{C48}
$$

and this difference cannot be told from the Fermi surface. The extreme case is a three-dimensional quantum Hall insulator, in which $S$ is empty (no Fermi surface) but $\sigma_{xy}=C_1^z e^2/hc$.

A weak coupling only changes the Fermi surface and the Berry phase continuously; the integer does not change unless the gap closes.

In the bulk formula with periodic boundary conditions, the boundary term of the integration by parts vanishes by the periodicity of the Brillouin zone, so only the integer part remains in the Fermi sea.
$\sigma_{xx}$ has nothing corresponding to this integer, so the integration by parts in Eq. (C21) moves it completely onto the Fermi surface.

