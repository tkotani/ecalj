# A finite sample with surfaces and the Thouless pump {#surface-and-pump}

> English version of the page [unfinished/surface_and_pump](../../unfinished/surface_and_pump.md). The pages at the top of ecaljdoc are the reference; this one was made from unfinished/surface_and_pump.md of commit 92b331c91 (source blob d2d852fe4).

::: info Unfinished
The part moved out of the notes of [Basics of the electrical conductivity calculation](../theory/conductivity), not yet sorted out. Equation numbers (C…) refer to that page. The equations of this page are numbered (S1), (S2).
:::

## A finite sample with surfaces {#slab-edge}

In a finite sample with surfaces, the integer part can also be regarded as carried by edge states (states localized at the surface) passing through $E_F$. We look at an example below.

**One sheet**: First consider a single two-dimensional Chern-insulator sheet spread in the $xy$ plane (Chern number $C$). Its Brillouin zone is the $T^2$ of $(k_x,k_y)$.
Cutting it at $x=0$ into the half-plane $x\gt 0$, $k_x$ is no longer a good quantum number, and the states are labeled only by $k_y$ along the edge ($T^1$).
On the edge there are $C$ states propagating in one direction only (chiral edge states); they disperse in $k_y$, cross the gap and pass through $E_F$.
The Hall conductance of the sheet is $\sigma_{xy}=C\,e^2/h$, and this Hall current is carried by the edge states. This is the bulk-edge correspondence.

> **Why there are $C$ edge states**: regard $k_y$ as the time $t$ of the [Thouless pump](#thouless-pump).
> At fixed $k_y$, the problem is a one-dimensional insulator in $x$. With periodic boundary conditions, its occupied bands define the Zak phase $\gamma(k_y)$ and the polarization in $x$, $P(k_y)=-e\gamma(k_y)/2\pi$ (Eq. (S1) with $t$ read as $k_y$).
> When $k_y$ goes once around the Brillouin zone, $\gamma$ advances by $2\pi C$ (the Chern number is the winding number of the Zak phase; Eq. (S2)). So one cycle of $k_y$ carries exactly $C$ electrons in $x$.
> In the half-plane with an edge, the carried electrons can only accumulate at the edge (or leave it). Since the Hamiltonian returns to itself after one cycle of $k_y$, $C$ levels must have crossed $E_F$ in between. The bulk has a gap at $E_F$, so the crossing levels are localized at the edge, and their branches $E(k_y)$ cross the gap.
> What is guaranteed is the net number (branches crossing upward minus those crossing downward); details of the surface can add or remove pairs of upward and downward branches, but the net $C$ does not change unless the bulk gap closes. The direction of crossing (the slope of $E(k_y)$) is the direction of the group velocity $v_y$, so the net motion is in one direction.

## The Thouless pump {#thouless-pump}

**Thouless pump** (D. J. Thouless, Phys. Rev. B **27**, 6083 (1983)): The Hamiltonian of a one-dimensional insulator (lattice constant $a$) is varied slowly with period $T$ ($H(t+T)=H(t)$). The gap between the occupied and unoccupied bands stays open during the whole cycle.
To each time $t$ we assign the phase obtained by integrating the Berry connection of the periodic parts $|u_{nk}(t)\rangle$ of the occupied Bloch functions over the $T^1$ of $k$ (the one-dimensional Brillouin zone), the Zak phase. The polarization is written with it (R. D. King-Smith and D. Vanderbilt, Phys. Rev. B **47**, 1651 (1993)):

$$
\gamma(t)=\sum_{n\in{\rm occ}}\oint_{T^1}i\langle u_{nk}(t)|\partial_k u_{nk}(t)\rangle\,dk,\qquad
P(t)=-\frac{e}{2\pi}\,\gamma(t)
\tag{S1}
$$

(the electron charge is $-e$). $\gamma$ changes by integer multiples of $2\pi$ under gauge transformations (Eq. (C47)), so $P$ is fixed only modulo $e$. This corresponds to the crystal being the same when an electron is moved to the neighboring unit cell.
Under an adiabatic variation, the current is $dP/dt$, and the charge passing a unit-cell boundary in one cycle is

$$
Q=\int_0^T\frac{dP}{dt}\,dt=-\frac{e}{2\pi}\big[\gamma(T)-\gamma(0)\big]=-eC
\tag{S2}
$$

($\gamma(t)$ is followed continuously in $t$). Since $H(T)=H(0)$, $e^{i\gamma(T)}=e^{i\gamma(0)}$, and $C$ is the winding number (C41) of the map $t\mapsto e^{i\gamma(t)}$ from the time loop $T^1$ to the unit circle. The charge carried in one cycle is quantized in integer multiples of $e$.

No singularity appears here. The gap stays open, and $|u_{nk}(t)\rangle$ is smooth in $(k,t)$. The integer comes only from $t$ being periodic, so that the loop cannot be shrunk.
$C$ is also the Chern number on the $T^2$ formed by $(k,t)$, a quantity of the same kind as $C_1^z$ in Eq. (C45) with $k_y$ replaced by $t$. Seen in a larger parameter space, it counts how many times the loop traced by the Hamiltonian goes around a gap-closing point lying outside the loop (in the Rice–Mele model, for example, $C=\pm 1$ when the circle traced by the dimerization and the sublattice potential difference encloses the point where both vanish).

In a finite chain with ends, a charge $-eC$ is carried from one end to the other in one cycle. Since the bulk has a gap, this happens by edge-state levels crossing the gap ($C$ of them per cycle). This is the one-dimensional version of the bulk-edge correspondence. For the two-dimensional Hall conductivity too, making the sample a cylinder and threading a flux lets $k_y$ play the role of $t$, and the same argument holds (Laughlin's argument).
