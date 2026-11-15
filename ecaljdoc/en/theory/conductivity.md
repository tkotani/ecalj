# Basics of the electrical conductivity calculation {#conductivity}

> English version of the page [theory/conductivity](../../theory/conductivity.md). The pages at the top of ecaljdoc are the reference; this one was made from theory/conductivity.md of commit bf1d12508 (source blob c5da84811423).

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

Expanding the state (or the density matrix) in the interaction picture to first order in $F$, the change of the expectation value of another operator $\hat A$ is

$$
\begin{aligned}
&\delta\langle\hat A\rangle(t)=\int_{-\infty}^{\infty}dt'\,\chi^R_{AB}(t-t')\,F(t'),\\
&\chi^R_{AB}(t)=\frac{i}{\hbar}\,\theta(t)\,\big\langle[\hat A_{\rm H}(t),\hat B]\big\rangle
\end{aligned}
\tag{C3}
$$

where $\hat A_{\rm H}(t)=e^{i\hat H_0t/\hbar}\hat A\,e^{-i\hat H_0t/\hbar}$ and $\langle\cdots\rangle$ is the thermal average for $\hat H_0$.
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
Writing the response function in (C3) as $\chi^R_{JJ}$, the current density is ($n$ is the electron density)

$$
j_x(t)=\frac{1}{\Omega}\,\langle\hat J_x\rangle(t)
=-\frac{ne^2}{m}A_x(t)+\frac{1}{\Omega}\int dt'\,\chi^R_{JJ}(t-t')\,A_x(t')
\tag{C6}
$$

$-ne^2A_x/m$ is the diamagnetic part and the $\chi^R_{JJ}$ term the paramagnetic part.

For the $e^{-i\omega t}$ component the convolution in the paramagnetic term becomes the product $\chi^R_{JJ}(\omega)A_x(\omega)$ ($\chi^R_{JJ}(\omega)$ is given by Eq. (C4)).
Moreover $A_x(\omega)=E_x(\omega)/i\omega$, so

$$
\sigma(\omega)=\frac{j_x(\omega)}{E_x(\omega)}=\frac{i}{\omega}\Big[\frac{ne^2}{m}-\frac{1}{\Omega}\,\chi^R_{JJ}(\omega)\Big].
\tag{C7}
$$

This is the Kubo formula. The real part, which fixes the dissipation (Joule heat), is ${\rm Re}\,\sigma(\omega)={\rm Im}\,\chi^R_{JJ}(\omega)/(\omega\Omega)$; the diamagnetic term affects only the imaginary part.

### The direct-current limit: the Kubo–Greenwood formula {#kubo-greenwood}

Inserting Eq. (C4) and using ${\rm Im}\,[1/(x-i\delta)]=\pi\delta(x)$,

$$
{\rm Re}\,\sigma(\omega)=\frac{\pi}{\omega\Omega}\sum_{nm}(f_n-f_m)\,|\langle n|\hat J^0_x|m\rangle|^2\,\delta(E_m-E_n-\hbar\omega).
\tag{C8}
$$

For $\omega\to 0$ the delta function gives $E_m=E_n+\hbar\omega$, so $f_n-f_m\simeq -f'(E_n)\,\hbar\omega$ and $\omega$ cancels.
At low temperature $-f'(E)=\delta(E-E_F)$. The one-body matrix element of $\hat J^0_x$ is $-e\langle n|\hat v_x|m\rangle$, and the spin sum gives a factor 2 ($n$, $m$ are orbital states), so
$\sigma_{xx}=(2\pi\hbar e^2/\Omega)\sum_{nm}|\langle n|\hat v_x|m\rangle|^2\delta(E_F-E_n)\delta(E_F-E_m)$.
This is fixed only by one-body quantities at the Fermi level $E_F$.
Rewriting it with $\delta(E-H)=-{\rm Im}\,G^R(E)/\pi$ and averaging over the impurity configurations gives the Kubo–Greenwood formula

$$
\sigma_{xx}=\frac{2\hbar e^2}{\pi\Omega}\,
\Big\langle {\rm Tr}\big[\hat v_x\,{\rm Im}\,G^R(E_F)\,\hat v_x\,{\rm Im}\,G^R(E_F)\big]\Big\rangle_{\rm imp}.
\tag{C9}
$$

$\hat v_x$ is the velocity operator, ${\rm Tr}$ the sum over one-body (orbital) states, and the factor 2 is the spin.

### The $G^RG^A$ terms {#grga}


Writing ${\rm Im}\,G^R=(G^R-G^A)/2i$ and expanding gives terms with $G^RG^A$ and terms with $G^RG^R$, $G^AG^A$.
For weak impurities the $G^RG^A$ terms dominate:

$$
\sigma_{xx}\simeq\frac{\hbar e^2}{\pi\Omega}\,{\rm Re}\,
\Big\langle {\rm Tr}\big[\hat v_x\,G^R(E_F)\,\hat v_x\,G^A(E_F)\big]\Big\rangle_{\rm imp}
\tag{C10}
$$

(the $G^RG^R$ terms only give corrections of higher order in $1/\tau$).

**The key point**: the average is taken of the **product** of $G^R$ and $G^A$:

$$
\big\langle \hat v_x G^R \hat v_x G^A\big\rangle_{\rm imp}\neq
\hat v_x\langle G^R\rangle_{\rm imp}\,\hat v_x\langle G^A\rangle_{\rm imp}.
\tag{C11}
$$

A correlation remains between $G^R$ and $G^A$ being scattered by the same impurity (in diagrams, impurity lines connecting the two lines). This is the vertex correction.

## Averaged single Green function: the lifetime {#averaged-g}

Averaging a single $G^R$ restores the translational symmetry, and it becomes diagonal in the wave vector ${\bf k}$:

$$
\langle G^R\rangle_{\rm imp}({\bf k},E)=\frac{1}{E-\varepsilon_{\bf k}-\Sigma^R({\bf k},E)}.
\tag{C12}
$$

The self-energy to second order in the impurity potential (Born approximation) is

$$
\Sigma^R({\bf k},E)=n_i\int\frac{d^3k'}{(2\pi)^3}\,
\frac{|u({\bf k}-{\bf k}')|^2}{E-\varepsilon_{\bf k'}+i\delta}
\tag{C13}
$$

and its imaginary part gives the lifetime $\tau$ of the quasiparticle (the real part is an energy shift, absorbed in $E_F$):

$$
\begin{aligned}
&{\rm Im}\,\Sigma^R=-\frac{\hbar}{2\tau},\\
&\frac{1}{\tau}=\frac{2\pi}{\hbar}\,n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\,\delta(E-\varepsilon_{\bf k'}).
\end{aligned}
\tag{C14}
$$

$1/\tau$ is just the scattering rate by Fermi's golden rule (see the correspondence between the golden rule and the imaginary part of the self-energy in [Basic I](./basic#td-vs-ti-perturbation)).

## Without the vertex correction: the Drude formula {#drude}

Replacing the average in Eq. (C10) by the product of the single averages (C12) (dropping the vertex correction) gives a sum over ${\bf k}$.
When $\tau$ is long ($\hbar/\tau\ll E_F$),

$$
\begin{aligned}
\langle G^R\rangle({\bf k},E)\,\langle G^A\rangle({\bf k},E)
&=\frac{1}{(E-\varepsilon_{\bf k})^2+(\hbar/2\tau)^2}\\
&\simeq\frac{2\pi\tau}{\hbar}\,\delta(E-\varepsilon_{\bf k})
\end{aligned}
\tag{C15}
$$

so, with the angular average $v_F^2/3$ of $v_x^2$ and the density of states per spin $N_1(E_F)$ (per unit volume),

$$
\sigma_{xx}=\frac{\hbar e^2}{\pi}\cdot\frac{2\pi\tau}{\hbar}\cdot\frac{v_F^2}{3}\,N_1(E_F)
=\frac{2e^2\tau v_F^2N_1(E_F)}{3}=\frac{ne^2\tau}{m}
\tag{C16}
$$

The last step uses $N_1(E_F)=mk_F/(2\pi^2\hbar^2)$, $n=k_F^3/(3\pi^2)$ and $v_F=\hbar k_F/m$ of free electrons. This is the Drude formula.

## Vertex correction: the transport relaxation time {#vertex-correction}

The product average (C11) is taken up to the ladder diagrams, in which impurity lines connect the two $G$ (the approximation that goes with the Born self-energy (C13)).
The current vertex $\Gamma_x({\bf k})$ satisfies

$$
\begin{aligned}
\Gamma_x({\bf k})={}&v_x({\bf k})+n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\\
&\times\langle G^R\rangle({\bf k}',E_F)\langle G^A\rangle({\bf k}',E_F)\,\Gamma_x({\bf k}')
\end{aligned}
\tag{C17}
$$

Putting $\Gamma_x({\bf k})=\gamma\,v_x({\bf k})$ and using (C15), the average over the direction of ${\bf k}'$ turns $v_x({\bf k}')$ into $v_x({\bf k})\cos\theta_{\bf kk'}$, so that
$\gamma=1+\gamma\,(1-\tau/\tau_{\rm tr})$, that is $\gamma=\tau_{\rm tr}/\tau$. Here

$$
\frac{1}{\tau_{\rm tr}}=\frac{2\pi}{\hbar}\,n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\,
\big(1-\cos\theta_{\bf kk'}\big)\,\delta(E_F-\varepsilon_{\bf k'})
\tag{C18}
$$

is the transport relaxation time. Replacing one $v_x$ in Eq. (C16) by $\Gamma_x$ gives

$$
\sigma_{xx}=\frac{ne^2\tau_{\rm tr}}{m}
\tag{C19}
$$

Forward scattering ($\cos\theta\simeq 1$) hardly reduces the current, so what matters for the resistance is the scattering weighted by $1-\cos\theta$.
The lifetime $\tau$ (fixed by the average of a single $G$) and $\tau_{\rm tr}$, which fixes the resistance (by the average of the product), differ in general.
If $u$ is point-like ($u({\bf q})$ independent of ${\bf q}$), the average of $\cos\theta$ is zero, the vertex correction vanishes, and $\tau_{\rm tr}=\tau$.
This result is the same as that of the Boltzmann equation.

Adding diagrams with crossed impurity lines, instead of the ladder, gives the weak-localization correction (not treated here).
