# Basics of the electrical conductivity calculation {#conductivity}

> English version of the page [theory/conductivity](../../theory/conductivity.md). The pages at the top of ecaljdoc are the reference; this one was made from theory/conductivity.md of commit dc084e8c4 (source blob 1832343159a5).

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

Applying an electric field ${\bf E}(t)$ as a perturbation and taking the response of the current density to first order, the conductivity is written by the retarded current-current correlation function (the Kubo formula).
For non-interacting electrons, the direct-current ($\omega\to 0$) conductivity at low temperature is written only by one-body quantities at the Fermi level $E_F$ (the Kubo–Greenwood formula):

$$
\sigma_{xx}=\frac{2\hbar e^2}{\pi\Omega}\,
\Big\langle {\rm Tr}\big[\hat v_x\,{\rm Im}\,G^R(E_F)\,\hat v_x\,{\rm Im}\,G^R(E_F)\big]\Big\rangle_{\rm imp}.
\tag{C2}
$$

$\hat v_x$ is the velocity operator, ${\rm Tr}$ the sum over one-body states, and the factor 2 is the spin.
Inserting ${\rm Im}\,G^R(E)=-\pi\delta(E-H)$ gives the sum over eigenstates $n$, $\sigma_{xx}=(2\pi\hbar e^2/\Omega)\sum_{nm}|\langle n|\hat v_x|m\rangle|^2\delta(E_F-E_n)\delta(E_F-E_m)$.

Writing ${\rm Im}\,G^R=(G^R-G^A)/2i$ and expanding gives terms with $G^RG^A$ and terms with $G^RG^R$, $G^AG^A$.
For weak impurities the $G^RG^A$ terms dominate:

$$
\sigma_{xx}\simeq\frac{\hbar e^2}{\pi\Omega}\,{\rm Re}\,
\Big\langle {\rm Tr}\big[\hat v_x\,G^R(E_F)\,\hat v_x\,G^A(E_F)\big]\Big\rangle_{\rm imp}
\tag{C3}
$$

(the $G^RG^R$ terms only give corrections of higher order in $1/\tau$).

**The key point**: the average is taken of the **product** of $G^R$ and $G^A$:

$$
\big\langle \hat v_x G^R \hat v_x G^A\big\rangle_{\rm imp}\neq
\hat v_x\langle G^R\rangle_{\rm imp}\,\hat v_x\langle G^A\rangle_{\rm imp}.
\tag{C4}
$$

A correlation remains between $G^R$ and $G^A$ being scattered by the same impurity (in diagrams, impurity lines connecting the two lines). This is the vertex correction.

## Averaged single Green function: the lifetime {#averaged-g}

Averaging a single $G^R$ restores the translational symmetry, and it becomes diagonal in the wave vector ${\bf k}$:

$$
\langle G^R\rangle_{\rm imp}({\bf k},E)=\frac{1}{E-\varepsilon_{\bf k}-\Sigma^R({\bf k},E)}.
\tag{C5}
$$

The self-energy to second order in the impurity potential (Born approximation) is

$$
\Sigma^R({\bf k},E)=n_i\int\frac{d^3k'}{(2\pi)^3}\,
\frac{|u({\bf k}-{\bf k}')|^2}{E-\varepsilon_{\bf k'}+i\delta}
\tag{C6}
$$

and its imaginary part gives the lifetime $\tau$ of the quasiparticle (the real part is an energy shift, absorbed in $E_F$):

$$
\begin{aligned}
&{\rm Im}\,\Sigma^R=-\frac{\hbar}{2\tau},\\
&\frac{1}{\tau}=\frac{2\pi}{\hbar}\,n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\,\delta(E-\varepsilon_{\bf k'}).
\end{aligned}
\tag{C7}
$$

$1/\tau$ is just the scattering rate by Fermi's golden rule (see the correspondence between the golden rule and the imaginary part of the self-energy in [Basic I](./basic#td-vs-ti-perturbation)).

## Without the vertex correction: the Drude formula {#drude}

Replacing the average in Eq. (C3) by the product of the single averages (C5) (dropping the vertex correction) gives a sum over ${\bf k}$.
When $\tau$ is long ($\hbar/\tau\ll E_F$),

$$
\begin{aligned}
\langle G^R\rangle({\bf k},E)\,\langle G^A\rangle({\bf k},E)
&=\frac{1}{(E-\varepsilon_{\bf k})^2+(\hbar/2\tau)^2}\\
&\simeq\frac{2\pi\tau}{\hbar}\,\delta(E-\varepsilon_{\bf k})
\end{aligned}
\tag{C8}
$$

so, with the angular average $v_F^2/3$ of $v_x^2$ and the density of states per spin $N_1(E_F)$ (per unit volume),

$$
\sigma_{xx}=\frac{\hbar e^2}{\pi}\cdot\frac{2\pi\tau}{\hbar}\cdot\frac{v_F^2}{3}\,N_1(E_F)
=\frac{2e^2\tau v_F^2N_1(E_F)}{3}=\frac{ne^2\tau}{m}
\tag{C9}
$$

The last step uses $N_1(E_F)=mk_F/(2\pi^2\hbar^2)$, $n=k_F^3/(3\pi^2)$ and $v_F=\hbar k_F/m$ of free electrons. This is the Drude formula.

## Vertex correction: the transport relaxation time {#vertex-correction}

The product average (C4) is taken up to the ladder diagrams, in which impurity lines connect the two $G$ (the approximation that goes with the Born self-energy (C6)).
The current vertex $\Gamma_x({\bf k})$ satisfies

$$
\begin{aligned}
\Gamma_x({\bf k})={}&v_x({\bf k})+n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\\
&\times\langle G^R\rangle({\bf k}',E_F)\langle G^A\rangle({\bf k}',E_F)\,\Gamma_x({\bf k}')
\end{aligned}
\tag{C10}
$$

Putting $\Gamma_x({\bf k})=\gamma\,v_x({\bf k})$ and using (C8), the average over the direction of ${\bf k}'$ turns $v_x({\bf k}')$ into $v_x({\bf k})\cos\theta_{\bf kk'}$, so that
$\gamma=1+\gamma\,(1-\tau/\tau_{\rm tr})$, that is $\gamma=\tau_{\rm tr}/\tau$. Here

$$
\frac{1}{\tau_{\rm tr}}=\frac{2\pi}{\hbar}\,n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\,
\big(1-\cos\theta_{\bf kk'}\big)\,\delta(E_F-\varepsilon_{\bf k'})
\tag{C11}
$$

is the transport relaxation time. Replacing one $v_x$ in Eq. (C9) by $\Gamma_x$ gives

$$
\sigma_{xx}=\frac{ne^2\tau_{\rm tr}}{m}
\tag{C12}
$$

Forward scattering ($\cos\theta\simeq 1$) hardly reduces the current, so what matters for the resistance is the scattering weighted by $1-\cos\theta$.
The lifetime $\tau$ (fixed by the average of a single $G$) and $\tau_{\rm tr}$, which fixes the resistance (by the average of the product), differ in general.
If $u$ is point-like ($u({\bf q})$ independent of ${\bf q}$), the average of $\cos\theta$ is zero, the vertex correction vanishes, and $\tau_{\rm tr}=\tau$.
This result is the same as that of the Boltzmann equation.

Adding diagrams with crossed impurity lines, instead of the ladder, gives the weak-localization correction (not treated here).
