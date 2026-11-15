# Basics of the electrical conductivity calculation {#conductivity}

> English version of the page [theory/conductivity](../../theory/conductivity.md). The pages at the top of ecaljdoc are the reference; this one was made from theory/conductivity.md of commit f54cf0be6 (source blob 2aef425e709e).

::: info Draft
As a **one-body problem**, the electrical resistance due to impurities is derived from linear response (the Kubo formula).
The key point is that the average over the impurity configurations is taken of a **product of Green functions** (averaging each one first and then multiplying is not enough).
The sign convention of the Green function is the same as in [Basic I](./basic); the retarded Green function is $G^R(E)=1/(E-H+i\delta)$.
For retarded response functions see also the [Appendix](./appendix#causality).
:::

## Setting {#setting}

Non-interacting electrons (mass $m$, $\varepsilon_{\bf k}=\hbar^2k^2/2m$) feel the potential of impurities placed at random in the volume $\Omega$,

$$
V({\bf r})=\sum_{j} u({\bf r}-{\bf R}_j)
\tag{C1}
$$

The number density of the impurities is $n_i$, and the Fourier transform of $u$ is $u({\bf q})=\int d^3r\, e^{-i{\bf q}\cdot{\bf r}}u({\bf r})$.
What is observed is the average $\langle\cdots\rangle_{\rm imp}$ over the configurations $\{ {\bf R}_j \}$ of the impurities.

**What is treated here is a one-body problem.** There is no interaction between the electrons; each electron moves independently in the static potential (C1) of the impurities.
The Hamiltonian is the one-body $H=\hat{\bf p}^2/2m+V({\bf r})$, and $G^R(E)=1/(E-H+i\delta)$ below is its one-body Green function (the resolvent).
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

### The direct-current limit: the Kubo–Greenwood formula {#kubo-greenwood}

Taking the real part of Eq. (C8), the diamagnetic term drops out and

$$
\begin{aligned}
{\rm Re}\,\sigma(\omega)&=\frac{1}{\omega\Omega}\,{\rm Im}\,\chi^R_{JJ}(\omega)\\
&=\frac{1}{\omega\Omega}\,{\rm Im}\sum_{nm}(f_n-f_m)\,\frac{|\langle n|\hat J^0_x|m\rangle|^2}{E_m-E_n-\hbar\omega-i\delta}\\
&=\frac{\pi}{\omega\Omega}\sum_{nm}(f_n-f_m)\,|\langle n|\hat J^0_x|m\rangle|^2\,\delta(E_m-E_n-\hbar\omega)
\end{aligned}
\tag{C9}
$$

In the second line Eq. (C4) is inserted with $\hat A=\hat B=\hat J^0_x$, and in the third ${\rm Im}\,[1/(x-i\delta)]=\pi\delta(x)$ is used.
For $\omega\to 0$ the delta function gives $E_m=E_n+\hbar\omega$, so $f_n-f_m\simeq -f'(E_n)\,\hbar\omega$ and $\omega$ cancels.
At low temperature $-f'(E)=\delta(E-E_F)$. The one-body matrix element of $\hat J^0_x$ is $-e\langle n|\hat v_x|m\rangle$, and taking the spin into account by multiplying by 2 ($n$, $m$ are orbital states),

$$
\sigma_{xx}=\frac{2\pi\hbar e^2}{\Omega}\sum_{nm}|\langle n|\hat v_x|m\rangle|^2\,\delta(E_F-E_n)\,\delta(E_F-E_m)
\tag{C10}
$$

This is fixed only by one-body quantities at the Fermi level $E_F$.
Rewriting it with $\delta(E-H)=-{\rm Im}\,G^R(E)/\pi$ and averaging over the impurity configurations gives the Kubo–Greenwood formula

$$
\sigma_{xx}=\frac{2\hbar e^2}{\pi\Omega}\,
\Big\langle {\rm Tr}\big[\hat v_x\,{\rm Im}\,G^R(E_F)\,\hat v_x\,{\rm Im}\,G^R(E_F)\big]\Big\rangle_{\rm imp}
\tag{C11}
$$

$\hat v_x$ is the velocity operator, ${\rm Tr}$ the sum over one-body (orbital) states, and the factor 2 is the spin.

Alternatively, one can start from the definition of $\chi^R_{JJ}$ and write it with Green functions from the beginning.
The Fourier transform of $\chi^R_{JJ}$ of Eq. (C6) is Eq. (C4) with $\hat A=\hat B=\hat J^0_x$.
Insert $1=\int dE\,\delta(E-E_n)$ in the $f_n$ term and use $\sum_m|m\rangle\langle m|/(E_m-E-\hbar\omega-i\delta)=-G^R(E+\hbar\omega)$.
Treating the $f_m$ term in the same way ($E_m\to E$, $\sum_n|n\rangle\langle n|/(E-\hbar\omega-E_n-i\delta)=G^A(E-\hbar\omega)$) and using $\delta(E-H)=-{\rm Im}\,G^R(E)/\pi$,

$$
\begin{aligned}
\chi^R_{JJ}(\omega)=\frac{1}{\pi}\int dE\,f(E)\,
{\rm Tr}\big[&\hat J^0_x\,{\rm Im}\,G^R(E)\,\hat J^0_x\,G^R(E+\hbar\omega)\\
&+\hat J^0_x\,G^A(E-\hbar\omega)\,\hat J^0_x\,{\rm Im}\,G^R(E)\big]
\end{aligned}
\tag{C12}
$$

Here ${\rm Re}\,G=(G^R+G^A)/2$ and ${\rm Im}\,G=(G^R-G^A)/2i$ are both Hermitian operators, and $G^R(E+\hbar\omega)={\rm Re}\,G+i\,{\rm Im}\,G$, $G^A(E-\hbar\omega)={\rm Re}\,G-i\,{\rm Im}\,G$.
The traces of the terms with ${\rm Re}\,G$ are real, so they do not contribute to the imaginary part of $\chi^R_{JJ}$. The imaginary part is given only by the products of ${\rm Im}\,G$, that is, by the poles.
Shifting $E\to E+\hbar\omega$ in the second term,

$$
\begin{aligned}
{\rm Re}\,\sigma(\omega)&=\frac{1}{\omega\Omega}\,{\rm Im}\,\chi^R_{JJ}(\omega)\\
&=\frac{1}{\pi\omega\Omega}\int dE\,\big[f(E)-f(E+\hbar\omega)\big]\\
&\quad\times{\rm Tr}\big[\hat J^0_x\,{\rm Im}\,G^R(E)\,\hat J^0_x\,{\rm Im}\,G^R(E+\hbar\omega)\big]
\end{aligned}
\tag{C13}
$$

This holds for any $\omega$ and temperature (for $\omega\neq 0$ it gives the optical conductivity).

### The $G^RG^A$ terms {#grga}

Writing ${\rm Im}\,G^R=(G^R-G^A)/2i$ and expanding gives terms with $G^RG^A$ and terms with $G^RG^R$, $G^AG^A$.
For weak impurities the $G^RG^A$ terms dominate:

$$
\sigma_{xx}\simeq\frac{\hbar e^2}{\pi\Omega}\,{\rm Re}\,
\Big\langle {\rm Tr}\big[\hat v_x\,G^R(E_F)\,\hat v_x\,G^A(E_F)\big]\Big\rangle_{\rm imp}
\tag{C14}
$$

(the $G^RG^R$ terms only give corrections of higher order in $1/\tau$).

**The key point**: the average is taken of the **product** of $G^R$ and $G^A$:

$$
\big\langle \hat v_x G^R \hat v_x G^A\big\rangle_{\rm imp}\neq
\hat v_x\langle G^R\rangle_{\rm imp}\,\hat v_x\langle G^A\rangle_{\rm imp}.
\tag{C15}
$$

A correlation remains between $G^R$ and $G^A$ being scattered by the same impurity (in diagrams, impurity lines connecting the two lines). This is the vertex correction.

## Averaged single Green function: the lifetime {#averaged-g}

Averaging a single $G^R$ restores the translational symmetry, and it becomes diagonal in the wave vector ${\bf k}$:

$$
\langle G^R\rangle_{\rm imp}({\bf k},E)=\frac{1}{E-\varepsilon_{\bf k}-\Sigma^R({\bf k},E)}.
\tag{C16}
$$

The self-energy to second order in the impurity potential (Born approximation) is

$$
\Sigma^R({\bf k},E)=n_i\int\frac{d^3k'}{(2\pi)^3}\,
\frac{|u({\bf k}-{\bf k}')|^2}{E-\varepsilon_{\bf k'}+i\delta}
\tag{C17}
$$

and its imaginary part gives the lifetime $\tau$ of the quasiparticle (the real part is an energy shift, absorbed in $E_F$):

$$
\begin{aligned}
&{\rm Im}\,\Sigma^R=-\frac{\hbar}{2\tau},\\
&\frac{1}{\tau}=\frac{2\pi}{\hbar}\,n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\,\delta(E-\varepsilon_{\bf k'}).
\end{aligned}
\tag{C18}
$$

$1/\tau$ is just the scattering rate by Fermi's golden rule (see the correspondence between the golden rule and the imaginary part of the self-energy in [Basic I](./basic#td-vs-ti-perturbation)).

## Without the vertex correction: the Drude formula {#drude}

Replacing the average in Eq. (C14) by the product of the single averages (C16) (dropping the vertex correction) gives a sum over ${\bf k}$.
When $\tau$ is long ($\hbar/\tau\ll E_F$),

$$
\begin{aligned}
\langle G^R\rangle({\bf k},E)\,\langle G^A\rangle({\bf k},E)
&=\frac{1}{(E-\varepsilon_{\bf k})^2+(\hbar/2\tau)^2}\\
&\simeq\frac{2\pi\tau}{\hbar}\,\delta(E-\varepsilon_{\bf k})
\end{aligned}
\tag{C19}
$$

so, with the angular average $v_F^2/3$ of $v_x^2$ and the density of states per spin $N_1(E_F)$ (per unit volume),

$$
\sigma_{xx}=\frac{\hbar e^2}{\pi}\cdot\frac{2\pi\tau}{\hbar}\cdot\frac{v_F^2}{3}\,N_1(E_F)
=\frac{2e^2\tau v_F^2N_1(E_F)}{3}=\frac{ne^2\tau}{m}
\tag{C20}
$$

The last step uses $N_1(E_F)=mk_F/(2\pi^2\hbar^2)$, $n=k_F^3/(3\pi^2)$ and $v_F=\hbar k_F/m$ of free electrons. This is the Drude formula.

## Vertex correction: the transport relaxation time {#vertex-correction}

The product average (C15) is taken up to the ladder diagrams, in which impurity lines connect the two $G$ (the approximation that goes with the Born self-energy (C17)).
The current vertex $\Gamma_x({\bf k})$ satisfies

$$
\begin{aligned}
\Gamma_x({\bf k})={}&v_x({\bf k})+n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\\
&\times\langle G^R\rangle({\bf k}',E_F)\langle G^A\rangle({\bf k}',E_F)\,\Gamma_x({\bf k}')
\end{aligned}
\tag{C21}
$$

Putting $\Gamma_x({\bf k})=\gamma\,v_x({\bf k})$ and using (C19), the average over the direction of ${\bf k}'$ turns $v_x({\bf k}')$ into $v_x({\bf k})\cos\theta_{\bf kk'}$, so that
$\gamma=1+\gamma\,(1-\tau/\tau_{\rm tr})$, that is $\gamma=\tau_{\rm tr}/\tau$. Here

$$
\frac{1}{\tau_{\rm tr}}=\frac{2\pi}{\hbar}\,n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\,
\big(1-\cos\theta_{\bf kk'}\big)\,\delta(E_F-\varepsilon_{\bf k'})
\tag{C22}
$$

is the transport relaxation time. Replacing one $v_x$ in Eq. (C20) by $\Gamma_x$ gives

$$
\sigma_{xx}=\frac{ne^2\tau_{\rm tr}}{m}
\tag{C23}
$$

Forward scattering ($\cos\theta\simeq 1$) hardly reduces the current, so what matters for the resistance is the scattering weighted by $1-\cos\theta$.
The lifetime $\tau$ (fixed by the average of a single $G$) and $\tau_{\rm tr}$, which fixes the resistance (by the average of the product), differ in general.
If $u$ is point-like ($u({\bf q})$ independent of ${\bf q}$), the average of $\cos\theta$ is zero, the vertex correction vanishes, and $\tau_{\rm tr}=\tau$.
This result is the same as that of the Boltzmann equation.

Adding diagrams with crossed impurity lines, instead of the ladder, gives the weak-localization correction (not treated here).
