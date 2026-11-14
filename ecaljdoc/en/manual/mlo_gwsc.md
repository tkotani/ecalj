# MLO-gwsc — QSGW interpolating the self-energy in the MLO representation (under development)

> English version of the page [manual/mlo_gwsc](../../manual/mlo_gwsc.md). The pages at the top of ecaljdoc are the reference; this one was made from manual/mlo_gwsc.md of commit 8e38c105c (source blob a5f53889fe6c).

> **Status (2026-09-30)**: LiTi₂O₄ at 6³ and 9³ was run from LDA up to 40 iterations (§3.5). The difference in the MLO bands due to the precision of the GPU matrix products (tf32, fp32, fp64) is
> a few meV or less, so the fast tf32 may be used. At 10 iterations it is still converging; at 6³ a shoulder appears in the occupied bands from iteration 20, and at 9³ the bands become fairly smooth at 40 iterations.
> In bands drawn in the MLO representation, the humps between mesh points that appeared with the `sigm` interpolation (interpolation ringing) are greatly reduced (§3.2).

> By default nothing changes (active only when `gwsc --mlo` is given). Combination with antiferromagnetic symmetry (`symgrpaf`) has not been checked (§5).
> The record of progress is <del class="clfix">ecalj's `MD/research_kotani_log.md` (2026-09-24 to 29); design details are in `MD/kBT/sigma_mlo_design.md` §9–13</del><span class="clfix">a development note (MD/research_kotani_log.md in the repository, 2026-09-24 to 29, not on this site); design details are in a development note (MD/mlo_notes/sigma_mlo_design.md §9–13 in the repository, not on this site)</span>.

This page has 3 parts: **§1 Theory**, **§2 Settings and running**, **§3 Example: LiTi₂O₄**.

---

## 1. Theory

### 1.1 What it solves

The static self-energy $\Sigma$ of QSGW is computed on the GW $q$ mesh, and at an arbitrary $k$ it is used after **interpolation in k space** from the values at the mesh points (Fourier interpolation: transform to real space and back).
Conventionally this is brought to real space as a matrix in the **MTO basis** (`sigm`) and brought back by a Bloch sum.
The MTO basis is non-orthogonal, EH and EH2 are nearly linearly dependent, and the tails are long, so $\Sigma^{\rm MTO}(R)$ may not decay at long distance
(in LiTi₂O₄ the Ti 3d block levels off at 0.3–0.6 meV up to the edge of the BvK cell). In k space this means it is not flat. That is,
humps of tens of meV, and when large hundreds of meV, appear in the bands **only "between" the mesh points**.
Making the mesh finer does not make the humps smaller; in LiTi₂O₄ they actually became larger (§3.2).

MLO-gwsc turns the same $\Sigma$ into a matrix represented in a small subspace $\{\tilde\chi_\alpha\}$ of **MLOs (localized orbitals based on MTOs)** and interpolates it in k space.
The MLOs are rebuilt at every iteration; $\Sigma$ is written with those MLOs and read back with the same MLOs (MLO Sigma interpolation, §1.4).

### 1.2 How the MLO $\tilde\chi$ is made

At each $k$, from the eigenstates $\Psi^{\rm PMT}_i(k)$ of the PMT Hamiltonian (eigenvalues $\varepsilon^{\rm PMT}_i$) and the eigenstates $\Psi^{\rm MTO}_j(k)$ obtained by diagonalizing only the MTO block (eigenvalues $\varepsilon^{\rm MTO}_j$), the projection

$$\hat P(k)=\sum_{ij}|\Psi^{\rm PMT}_i\rangle\,\bar\theta_{ij}\,\langle\Psi^{\rm PMT}_i|\Psi^{\rm MTO}_j\rangle\langle\Psi^{\rm MTO}_j| \tag{1}$$

is built and applied to the MTO seeds $\chi^{\rm MTO}_\alpha$ selected by `mlo_lm`: 
$$|\tilde\chi_\alpha(k)\rangle=\hat P(k)|\chi^{\rm MTO}_\alpha(k)\rangle \tag{2}$$
(`Hreduction` in `m_hreduction.f90`).
The window function $\bar\theta_{ij}$ is that of `mlo_method = 4` of [MLO](./mlo),
$$\bar\theta_{ij} =f\big((\varepsilon^{\rm PMT}_i-\varepsilon^{\rm cut}_j)/w\big),\qquad
\varepsilon^{\rm cut}_j=\max(\varepsilon_{\rm CBM}+\Delta,\ \varepsilon^{\rm MTO}_j) \tag{3}$$
($\varepsilon_{\rm CBM}$ is the conduction-band bottom of the SCF of each iteration; for metals it is almost $E_F$). $\Delta$ = `mlo_delta` and $w$ = `mlo_w`, both 2 eV by default, and $f$ is the Fermi function. All samples of the new MLO model use this form.
Since it only contains the energy function $\bar\theta_{ij}$, and the only feature of the system needed is $\varepsilon_{\rm CBM}$, it should roughly work even for a system containing, e.g., two slabs of Fe and Si.
$\Psi$ enters only in the form $|\Psi\rangle\langle\Psi|$, so there is no gauge (phase) arbitrariness.

We write the coefficients of $\tilde\chi$ in the PMT basis as $z^{\rm MLO}(k)$: 
$$|\tilde\chi_\alpha(k)\rangle=\sum_m|\chi^{\rm PMT}_m(k)\rangle z_{m\alpha}(k) \tag{4}$$

$\tilde\chi$ is not a stored table but **the prescription applied to the Hamiltonian at that moment**; it changes when $H$ changes.
From iteration 2 on it is built from $H$ including $\Sigma$. The standard MLO method also builds from `__HamiltonianPMT` ($H$ after adding senex), so this is the same.
In any case, MLOs can be built as long as the Hamiltonian (including the basis functions) is available on the regular mesh points.

* In the end, the MLO method landed on a "very simple approach". In short, $\bar\theta_{ij}$, which firmly extracts the low-energy part, just entered the complete-set expansion of Eq. (1).
  That is, "if a good localized basis (like MTOs) is chosen, a projection that respects the low-energy side models the system cleanly enough".
  There is no need to worry about symmetry or about convergence. For slabs etc. MLOs corresponding to empty spheres may need to be brought in, but this is the same for other methods.
* Why $\varepsilon^{\rm cut}_j$ in Eq. (3) is not set below $\varepsilon^{\rm MTO}_j$: if bands are drawn with the MTO basis only, without the APW basis, the energies come out
  higher by the amount the basis is lacking. The PMT eigenstates corresponding to MTO state $j$ lie below $\varepsilon^{\rm MTO}_j$, so
  a PMT eigenstate is taken in if its energy is at or below $\varepsilon^{\rm MTO}_j$.

### 1.3 Holding $\Sigma$ in MLOs and returning to PMT

What GW produces is $\Sigma^\psi_{ij}(q)=\langle\psi_i|\hat\Sigma|\psi_j\rangle$ with band indices. This is transferred to MLOs:

$$\Sigma^{\rm MLO}=c^\dagger\,\Sigma^\psi\,c,\qquad c=\langle\psi|\tilde\chi\rangle \tag{5}$$

This is literally $\langle\tilde\chi_\alpha|\hat\Sigma|\tilde\chi_\beta\rangle$ (<del class="clfix">states outside the window $\bar\theta$ are dropped, so it is exact only inside the window</del><span class="clfix">$c_i$ is nonzero only for bands where $\bar\theta$ is nonzero. It is exact if $\Sigma^\psi$ covers all bands with non-negligible $\bar\theta$, including the tail of $f$ (`ntq` large enough); otherwise the tail is cut off</span>).
The PMT basis is non-orthogonal; we write its overlap matrix as $S^{\rm PMT}_{mn}(k)=\langle\chi^{\rm PMT}_m(k)|\chi^{\rm PMT}_n(k)\rangle$.
Expanding the eigenstates as $|\psi_i\rangle=\sum_m|\chi^{\rm PMT}_m\rangle\psi_{mi}$ ($\psi$ is the coefficient matrix), Eq. (4) gives $c=\psi^\dagger S^{\rm PMT}z^{\rm MLO}$.
$\Sigma^{\rm MLO}(q)$ is brought to real space as $\Sigma^{\rm MLO}(R)$ (`QMLO_SigRs`), interpolated in k space, and at an arbitrary $k$ returned to the PMT basis and added to $H$ (`getsenex`).
For general (non-orthogonal) MLOs the inverse of the overlap is needed on both sides (from 2026-10-02 the MLOs are Löwdin-orthogonalized functions, so $O=1$ and Eq. (7) is $A\,\Sigma^{\rm MLO}A^\dagger$. [mlo](mlo) §6):

$$A=S^{\rm PMT}z^{\rm MLO}=\langle\chi^{\rm PMT}|\tilde\chi\rangle,\qquad O=z^{{\rm MLO}\dagger}A=\langle\tilde\chi|\tilde\chi\rangle \tag{6}$$

$$\mathrm{senex}=A\,O^{-1}\,\Sigma^{\rm MLO}\,O^{-1}A^\dagger=\langle\chi^{\rm PMT}|\hat P_{\rm MLO}\hat\Sigma\hat P_{\rm MLO}|\chi^{\rm PMT}\rangle \tag{7}$$

$\hat P_{\rm MLO}=|\tilde\chi\rangle O^{-1}\langle\tilde\chi|$ is the projector onto the MLO subspace (different from $\hat P(k)$ of Eq. (1)). **Only $\Sigma^{\rm MLO}$ is interpolated in k space.**

> **Note**: the conventional `sigm` is a quantity projected onto the MTO block with a pseudo-inverse, i.e. a form with $O^{-1}$ folded in on both sides.
> The magnitude of $\Sigma(R)$ cannot be compared directly with $\Sigma^{\rm MLO}(R)$.

### 1.4 MLO Sigma interpolation — $\Sigma$ is read with the $\tilde\chi$ that wrote it

The MLOs are rebuilt at every iteration; the $\Sigma$ of that iteration is turned into the matrix $\Sigma^{\rm MLO}(k)$ in the MLO representation and **interpolated in k space**
(values at an arbitrary $k$ are made from the values at the mesh points), and the same MLOs are used when returning to the PMT basis.
We call this **MLO Sigma interpolation**. This form is determined by the following 2 facts.

1. **What can be interpolated in k space is the matrix in the MLO representation; $z^{\rm MLO}(k)$ cannot be interpolated.** Both rows and columns of $\Sigma^{\rm MLO}_{\alpha\beta}(k)$ are MLO indices,
   which point to the same channels at every $k$, so it can be interpolated in k space. On the other hand $z^{\rm MLO}(k)$ has rows in the PMT basis.
   With `pwmode = 11` ($|q+G|$ cutoff) the set of APWs itself changes with $k$, so the same row number points to a different function at each $k$, and
   it cannot be interpolated in k space ($A$ of Eq. (6) likewise). For the same reason $H^{\rm PMT}(k)$ cannot be interpolated either.
2. **The numbers of $\Sigma^{\rm MLO}$ carry the stamp of the $\tilde\chi$ that wrote them.** At iteration $n-1$ they are made in a single chain
   $z_{n-1}\to c_{n-1}=\psi^\dagger S^{\rm PMT} z_{n-1}\to\Sigma_{n-1}=c_{n-1}^\dagger\Sigma^\psi c_{n-1}$, so
   when returning by Eq. (7) **the same $z_{n-1}$ must always** be used. Reading with another $z$ reinterprets the same numbers as matrix elements between different functions,
   and changes the operator itself.

Therefore **the $z$ made on the writing side is saved for all $k$ that will be read later** ($\tilde\chi$ changes when $H$ changes, so rebuilding it on the reading side gives something different).
Within one iteration both "$z_{n-1}$ for reading" and "$z_n$ for writing" are needed at the same time, so there are 2 places to keep them.

| File | Content |
|---|---|
| `QMLO_z` | the $\tilde\chi$ with which the current `QMLO_SigRs` was written (for reading)|
| `__QMLO_zNew` | the $\tilde\chi$ made from the $H$ of this iteration (for writing)|

Each is a single file with a fixed-length record ($q$, spin, `ndimh` and $z$) per GW $(q,$ spin$)$.
In stage a' each process writes the records of its own $q$ with MPI-IO; `getsenex` first reads only the index and then reads only the records of the $k$ it needs.
**Right after** `mlo` writes `QMLO_SigRs`, `__QMLO_zNew` is renamed (promoted) to `QMLO_z`, so $\Sigma$ and $\tilde\chi$ always form a pair.
To keep "all $k$ that will be read later" within the saved points, **`nkabc` = `n1n2n3` = `mlo_nkabc`** (SCF $k$ ⊂ GW $q$) is required.
`gwsc` checks this before the calculation and stops if they differ.
Attempting to read at a $k$ that is not saved also stops by default (since that would read $\Sigma$ in another basis).

### 1.5 Flow of one iteration (`gwsc 1 --mlo`)

Only once, at the start of the chain, `lmf --writeham --mlo` and `mlo --mlo` build the **MLO index** in `HamRsMLO`.
The index is which MTO channels are used as seeds (`ndimMTO` of them). $\Sigma^{\rm MLO}$ is a matrix on those channels, so the index is fixed throughout the chain (`--mlofreeze`).
The reference of the window ($E_F$ and the conduction edge) is not included in the index; stage a' reads it each iteration from `efermi.lmf` of that iteration's SCF (from <span class="clfix">commit 5118934d3 on the evening of </span>2026-09-26; §2.4).
On the other hand, the coefficients $z$ are rebuilt in ① at every iteration and follow the Hamiltonian.

| # | Run | About $\tilde\chi$ |
|---|---|---|
| ① | `lmf --jobgw=1` | `getsenex` reads `QMLO_SigRs` ($\Sigma_{n-1}$) with $z_{n-1}$ of `QMLO_z` → diagonalization → $\psi$ (passed to GW). In the following stage a', $z_n$ is made from $H$ including $\Sigma$ and written to `__QMLO_zNew`, and $c_n=\psi^\dagger S^{\rm PMT} z_n$ to `__cmlo.data` |
| ② | `heftet`, `hbasfp0`, `hvccfp0`, `hsfp0_sc`, `hgw` | $\Sigma^\psi_{ij}(q)$. Unrelated to MLO |
| ③ | `hqpe_sc` | Eq. (5) gives $\Sigma^{\rm MLO}_n$ → mixing (§1.6) → `__QMLO_Sig`. The conventional `sigm` is also written alongside |
| ④ | `mlo --mlofreeze --mlo` | symmetrization → full BZ → FFT → `QMLO_SigRs` (Σ only; MLO indices from `HamRsMLO`)|
| ⑤ | promotion (inside the `mlo` of ④)| `__QMLO_zNew` → `QMLO_z`. Here the "current $\tilde\chi$" switches to $z_n$ |
| ⑥ | `lmf` (SCF)| `getsenex` reads $\Sigma_n$ with $z_n$ |

### 1.6 Mixing

In the MLO path, what the SCF reads is `QMLO_SigRs` (made from `__QMLO_Sig`), not `sigm`.
`sigm` is Anderson-mixed in `hqpe_sc`, but unless $\Sigma^{\rm MLO}$ is also mixed, **only the MLO chain runs with an effective β=1**,
and in LiTi₂O₄ it overshot greatly at the 3rd iteration (§3.3).
$\Sigma^{\rm MLO}$ gets the same Anderson mixing as `sigm` (β = `[gw] mixbeta`) (default from 2026-09-30; `ECALJ_MLO_MIX=0` turns it off.
In earlier code it was mixed only when `ECALJ_MLO_MIX=1` was given).

- The history is `__QMLO_mixsig` (`__mixsig` on the `sigm` side). The input of the mixing $x_0$ (the $\Sigma^{\rm MLO}$ actually used in that iteration) is
  the Bloch sum of `QMLO_SigRs` at the irreducible $q$ (the same matrix that `getsenex` put into $H$). It plays the same role as `sigm` in the `sigm` mixing, so the history may be discarded
- **The first time (from LDA) starts with a factor $\beta$**: there is no previous $\Sigma$, so mixing starts from $x_0=0$ ($\Sigma=0$), and iteration 1 is $\beta\Sigma^{\rm out}$.
  Both `sigm` and $\Sigma^{\rm MLO}$ are scaled by $\beta$; with the default $\beta=1$ the full amount enters ($\beta=0.5$ for LiTi₂O₄. [mixbeta](./lmf#mixbeta-and-mixpriorit))

### 1.7 How to draw the MLO bands

To see the state of an intermediate iteration at an arbitrary $k$ while keeping $\Sigma^{\rm MLO}$, build and solve **the MLO model of the Hamiltonian solved by the SCF**
(ecalj's `Samples/kBT/LiTi2O4/draw_mloband.sh`).

1. On a snapshot, run `lmf --writeham --mkprocar --noinv --mlo` and, with the saved `QMLO_z`, write out at the mesh points the same
   $H=H^{\rm LDA}[\rho_n]+\mathrm{senex}(\Sigma^{\rm MLO}_n)$ as the SCF (it stops if $\tilde\chi$ is not found at even one point)
2. After moving `QMLO_SigRs` aside ($\Sigma$ is already in the $H$ of step 1, so do not add it twice), build the MLO model with `mlo --mlo` without freezing and solve it along the path

At the mesh points it is the SCF $H$ itself, and in between it is the Fourier interpolation of the standard MLO. **$\tilde\chi$ is not rebuilt between the mesh points.**
The following 3 ways do not give MLO bands.

| Path | What happens |
|---|---|
| `job_band` (without `--mlo`)| drawn with the conventional `sigm` interpolation. It combines the density of the MLO chain with the conventional $\Sigma$, not the $\Sigma^{\rm MLO}$ that drove the SCF (usable for comparison)|
| `job_band --mlo` | the search for $\tilde\chi$ compares $k$ literally, so it is not found even when a $k$ on the path is symmetry-equivalent to a saved irreducible point. It is rebuilt from $H$ without $\Sigma$, so it is invalid except at Γ |
| `job_mlo --mlofreeze --mlo` | the current `QMLO_SigRs` (current $\tilde\chi$) is added to the frozen `HamRsMLO` (the LDA $H$ and $\tilde\chi$ at the start of the chain), so 2 bases are mixed and it deviates by 0.3–0.5 eV |

`job_band` writes the $E_F$ for drawing to `efermi.lmf.job_band` and does not touch `efermi.lmf`, from which stage a' reads the window reference (the same for `job_dos` and `job_fermisurface`.
It can also be given with `--efermi=`). Drawing inside the chain's directory does not move the window.

---

## 2. Settings and running

### 2.1 `ctrlg.<sname>.toml`

The settings used for LiTi₂O₄ at 6³ (for 9³ only the 3 mesh entries were changed to `[9, 9, 9]`):

```toml
[bz]
nkabc  = [6, 6, 6]          # same as [gw] n1n2n3 (required)
[ham]
pwmode = 11                 # |q+G| cutoff (default)
[gw]
n1n2n3  = [6, 6, 6]
mixbeta = 0.5               # mixing of sigm and Sigma^MLO
t_sigmaw   = 1000.0         # (K) Fermi-Dirac width of the intermediate levels of Sigma (default)
t_tetrakbt = -992.4         # (K) negative: chi0 stays T=0 and Im chi0 is smeared by a Gaussian as wide as FD at about 1000 K (former SmearX0 = 0.0057 Ha). Smears the poles of W of metals (<del class="clfix">"Practical settings" of [kBT](./kBT)</del><span class="clfix">"Settings for LiTi₂O₄" at the top of [kBT](./kBT), and its §2.4</span>)
[mlo]
mlo_method = 4              # default
mlo_delta  = 2.0            # (eV) default
mlo_w      = 2.0            # (eV) default
mlo_nkabc  = [6, 6, 6]      # same as n1n2n3 (required)
mlo_lm = """
1 Li   1 2 3 4 5 6 7 8 9
2 Li   1 2 3 4 5 6 7 8 9
3 Ti   1 2 3 4 5 6 7 8 9
...
14 O   1 2 3 4 5 6 7 8 9
"""
```

- The `gwinit` default of `mlo_lm` is s,p (lm 1–4) for Z ≤ 10, s,p,d (1–9) from Na on, and also f (10–16) for atoms with 4f valence (from 2026-10-01. [MLO](./mlo) §1).
  Local orbitals of semicore states enter automatically by the band position (table in MLO §1).
  For LiTi₂O₄ all atoms were given s+p+d (126 orbitals). With the default it would be 76 (s,p for Li and O). Which of 76 and 126 is better is not settled (§5)
- **In principle at most one MLO per (site, lm).** Adding `mlo_lm2` (EH2) on the same lm collapses to the same function under the projection and makes $O$ nearly singular
  (exception: s,p of cations with a conduction-band bottom spread into the voids: criterion 2 of [MLO](./mlo) §9; not tried in MLO-QSGW).
  <span class="clfix">With the raw MLOs used until 2026-10-02, </span>Eq. (7) sandwiche<del class="clfix">s</del><span class="clfix">d</span> $O^{-1}$ twice, so the effect <del class="clfix">is</del><span class="clfix">was</span> squared (for Ni d of NiO the condition number goes $1.3\times10^2\to7.4\times10^4$)<span class="clfix">. Now the nearly singular $O$ enters through the $O^{-1/2}$ of the Löwdin orthogonalization</span>

### 2.2 Running

```bash
gwsc 10 -np 60 -np2 2 --gpu --prec=fp32 --ntqxx --mlo liti2o4                     # 10 iterations
gwsc 10 -np 60 -np2 2 --gpu --prec=tf32 --prec-final=fp32:2 --ntqxx --mlo liti2o4 # tf32 in between, fp32 for the last 2 iterations
```

`gwsc N --mlo` may be run for N iterations at once (it is not strictly the same as running `gwsc 1` N times, but once converged the difference in ehf is 0.1–0.2 meV; research log 2026-09-25 00:06).
For the GPU precision `--prec=tf32|fp32|fp64` and the automatic selection of the matrix-product method, see [ecaljgpu](./ecaljgpu) and [gwsc](./gwsc) (`--prec=fp32` is the same as the old `--mp --fp32`).
One `hgw` run for LiTi₂O₄ takes, at 6³, 367 s with fp32 and 173 s with tf32 (kt1, RTX 5090 × 2, evening of 2026-09-27), and at 9³, 2990 s with fp32 and 1975 s with tf32 (version of the morning of the same day).
tf32 makes only the last matrix product of Σc take inputs with a 10-bit mantissa (FP16 on RTX 5090<span class="clfix">, chosen by the per-GPU table; FP16 also has only a 5-bit exponent, so its range is narrow. [ecaljgpu](./ecaljgpu)</span>); Re Σc within $E_F$ ±1 eV is within 1.3 meV of double precision (<del class="clfix">`MD/ForDevelopers.md` §11.2</del><span class="clfix">development note MD/ForDevelopers.md §11.2 in the repository, not on this site</span>).

- ecalj's `Samples/kBT/LiTi2O4/run_gwsc10.sh <tag> <input> <bindir>` runs `gwsc 10` and, in the final state, draws the `sigm` plot and the MLO bands.
  `run_snap.sh <tag> <niter> --mlo` repeats `gwsc 1` and keeps a snapshot of every iteration (the full set needed to restart) and the 2 kinds of bands
  (it stops if even one mixing file remains after LDA). `run_gwsc10.sh` can be given `RUNS_DIR` (working location) and `GPUS` (default one GPU, GPU 0);
  the input is `Samples/kBT/LiTi2O4/input/qmlo`. `run_snap.sh` has kt1 paths hard-coded, so rewrite them when using it

### 2.3 Checks at each iteration

What to look at (in the `run_snap.sh` chain, `Samples/kBT/LiTi2O4/check_iter.sh` examines the logs saved by `watch_iter.sh` at the end of each iteration).
When run at once with `gwsc N`, `gwsc` copies `lgw`, `lqpe`, `llmfgw01`, `lmlo_sigr`, `llmf` (`llmf.<N>run`) to `QSGW.<N>run/` at the end of each iteration
(from 2026-09-28), so afterwards grep there as in the table below (`check_iter.sh` assumes the format of `steps.log` of run_snap):

| Item | Where to look | Normal |
|---|---|---|
| whether the SCF entered the MLO path | number of `MLO Sigma interpolation ON` lines in `llmf` | number of processes holding $k$ points $\lceil n_k/\lceil n_k/n_p\rceil\rceil$ |
| whether the GW driver entered the MLO path | the same line in `llmfgw01` | number of processes holding GW $q$ points (`QPLIST.jobgw1`). **0 in iteration 1** (no `QMLO_SigRs` yet)|
| promotion | `m_HamPMT: promoted __QMLO_zNew -> QMLO_z` in `lmlo_sigr` | once per iteration |
| mixing | `x_0 from QMLO_SigRs =` and `SigmMLO_out`, `SigmMLO_mixed` in `lqpe` | in iteration 1, 2 lines of `... No mixing file`, x_0 = F, mixed/out = β. From 2 on, x_0 = T |
| missed $\tilde\chi$ | `m_sigmlo: WARNING chi~ MISS` | does not appear (if it does, it stops by default)|
| SCF | number of iterations in `llmf` | not reaching the upper limit of `nit` |

The line counts **are determined by the mesh and the number of processes**, so do not apply values for another mesh.
With `-np 60`, 6³ has 16 irreducible points, giving 16, and 72 GW points, giving 36. 9³ gives 35 for 35 points and 50 for 200 points.
If **fewer** than this, some processes did not enter the MLO path and fell back to the conventional `sigm`.

### 2.4 Files used, and what to keep

Marks: **●** specific to MLO-QSGW (appears only with `gwsc --mlo`), **○** MLO method in general (also appears with `job_mlo`), no mark = common with ordinary QSGW.
Keep: **result** = keep as a calculation result, **restart** = needed to continue the chain, **—** = work file used only within an iteration.
① to ⑥ are the stages of §1.5.

| File | Mark | Keep | Content |
|---|---|---|---|
| `ctrlg.<sname>.toml` | | result, restart | settings |
| `rst.<sname>` | | result, restart | density (written by `lmf` of each SCF)|
| `efermi.lmf` | | result, restart | $E_F$ and conduction-band bottom of the last SCF. Read every iteration by stage a' as the window reference |
| `sigm`, `QPU`, `QPD` | | result (`sigm` also restart)| $\Sigma$ in the conventional representation (plots for comparison, $x_0$ of the `sigm` mixing), table of QP energies |
| `__mixsig` | | restart | history of the `sigm` mixing |
| `QMLO_SigRs` | ● | result, restart | $\Sigma^{\rm MLO}(R)$. The $\Sigma$ of MLO-QSGW itself (written by ④). The mixing $x_0$ is also made from it |
| `QMLO_z` | ● | result, restart | the $\tilde\chi$ with which `QMLO_SigRs` was written (a single file). **Always paired with `QMLO_SigRs`** |
| `HamRsMLO` | ○ | result, restart | MLO index (built at the start of the chain, then fixed). Deleting it in the middle makes it rebuilt and no longer consistent with `QMLO_SigRs` |
| `band_MLO_spin1.dat` | ○ | result | MLO bands (§1.7)|
| `llmf`, `llmfgw01`, `lqpe`, `lmlo_sigr` | | result | logs (used for the checks of §2.3)|
| `__QMLO_mixsig` | ● | — | history of the $\Sigma^{\rm MLO}$ mixing. Deleting it only breaks the Anderson extrapolation once ($x_0$ is made from `QMLO_SigRs`)|
| `__QMLO_Sig` | ● | — | $\Sigma^{\rm MLO}(q)$ written by ③ and read by ④ |
| `__QMLO_zNew` | ● | — | the current $\tilde\chi$ written by stage a' of ① (a single file, MPI-IO)|
| `__QMLO_zdump_*` | ● | — | for diagnostics (`ECALJ_ZMLO_DUMP=1`)|
| `__cmlo.*`, `__HamiltonianGW(.info)` | ○ | — | $c_n$ and the PMT $H$ written by ①. In the chain only ① writes `__cmlo` (in `job_mloW` and `job_mlo_magnon`, `mlo --mlo` writes it from `__HamiltonianGW`)|
| `HamiltonianPMTInfo` | ○ | — | read every iteration by the `mlo` of ④ (lattice, k mesh, table of pairs). If absent, `gwsc` rebuilds it at startup with `lmf --writeham` |
| `__HamiltonianPMT(.info)` | ○ | — | read only by the `mlo --mlo` that builds the index (at the start of the chain). ④ does not read it (2026-09-27)|
| `__amlo.*` | ○ | — | written by the `mlo --mlo` that builds the index, but nothing reads it (④ does not write it)|

Files beginning with `__` may all be deleted between calls of `gwsc` (deleting `__mixsig` and `__QMLO_mixsig` only breaks the Anderson extrapolation once).
`sigm` is not used by the SCF as a value, but without it both `gwsc` and `lmf` restart from LDA, so it is required for restarting.
The snapshot of `run_snap.sh` is most of the "result" files plus `sigm.<sname>`, `__atm.<sname>`, `syml.<sname>`, etc.

### 2.5 Options and environment variables

| | Effect |
|---|---|
| `--mlo` | switch for the whole MLO path. Without it nothing MLO-related runs |
| `--mlofreeze` | uses the index of the existing `HamRsMLO` and rewrites only `QMLO_SigRs` (`gwsc` passes it to the `mlo` within the iteration)|
| `ECALJ_MLO_MIX=0` | turns off the Anderson mixing of $\Sigma^{\rm MLO}$ (§1.6; mixed by default)|
| `ECALJ_MLO_ALLOW_REBUILD=1` | allows rebuilding from the current $H$ at a $k$ not in `QMLO_z` (stops by default). Bases get mixed, so for diagnostics |
| `ECALJ_MLO_NOSIG=1` | ignores `QMLO_z` and rebuilds from the current $H$ at all $k$. For diagnostics |
| `ECALJ_MLO_NOCACHE=1` | stops using both `QMLO_z` and reuse, rebuilding at every call (old behavior). For diagnostics |
| `ECALJ_ZMLO_DUMP=1` | writes the $z$ of stage a' to `__QMLO_zdump_sugw` and the $z$ of `getsenex` to `__QMLO_zdump_getsenex`. They should agree between a' and the SCF of ⑥ in the same iteration |
| `ECALJ_SIGMLO_RT=1` | checks at the mesh points whether the Bloch sum of $\Sigma^{\rm MLO}(R)$ reproduces `__QMLO_Sig` |
| `ECALJ_SIGMLO_CHECK=1` | prints, at each $k$, the maximum difference between the senex of the MLO path and the senex of the conventional `sigm` path |

---

## 3. Example: LiTi₂O₄

### 3.1 Settings

| Item | Value |
|---|---|
| MLO | **126 orbitals = 14 atoms (Li₂Ti₄O₈) × (s+p+d)** |
| Mesh | `nkabc` = `n1n2n3` = `mlo_nkabc` = 6³ and 9³ |
| Mixing | `[gw] mixbeta = 0.5` (both `sigm` and $\Sigma^{\rm MLO}$)|
| Basis | `pwmode = 11` |
| Window | `mlo_method = 4`, $\Delta = w = 2.0$ eV (default)|
| Machine | kt1 (RTX 5090 × 2), `gwsc 1 --mp --fp32` 10 times (code of 2026-09-26<span class="clfix">, before the window fix of that evening</span>)|
| Window reference | fixed to the $E_F$ at the start of the chain (LDA) (the current code uses the value of each iteration's SCF, <del class="clfix">§1.2</del><span class="clfix">§1.5</span>)|
| One iteration | 6³: 950–1000 s, 9³: about 5300 s (89 min). `hgw` takes 86 % (6³) and 94 % (9³); processing for MLO is about 1 %. After the speedup of 2026-09-27, 555 s for 6³ with fp32 |
| Reference for comparison | conventional QSGW (MTO) 6³, `pwmode = 1`, β = 0.5, 10 iterations (only the LDA row uses `pwmode = 11`)|

### 3.2 Bands

![Evolution of the t2g bands of LiTi2O4 over iterations (conventional 6³ / MLO-QSGW 6³ / MLO-QSGW 9³)](../../manual/mlo_gwsc/liti2o4_5col_all.png)

t2g (bands 33–44), Γ→X, 211 points. Columns from the left: **conventional QSGW 6³** / **6³ MLO-QSGW drawn with `sigm`** / **6³ MLO bands** /
**9³ MLO-QSGW drawn with `sigm`** / **9³ MLO bands**; rows from the top: LDA, iterations 1–10.
Red × are the mesh points of $\Sigma$ (points where the interpolation is exact). Green is the lowest band, and blue is the band that bulges in the middle of Γ→X;
it is tracked across crossings, and colored only in iterations where the bulge is 50 meV or more.

- **The hump of the bulging band (blue) does not appear in the MLO bands.** When the same state is drawn with the `sigm` interpolation (columns 2 and 4), humps form between the mesh points.
  At 6³ it swings strongly at iteration 3, at 9³ at iterations 3–6, and settles from 7 on.
  It never appears in the 6³ MLO bands, and in the 9³ MLO bands it does not appear from iteration 5 on
- **The humps are larger at 9³** (`sigm` plot). In conventional QSGW too, making the mesh denser enlarged the ripples of the occupied t2g.
  Our view is that $\Sigma^{\rm MTO}(R)$ does not decay at long distance, so the larger the periodic cell, the more of that component is picked up and interpolated
- The lowest band (green) shows step-like shoulders in the early iterations (the step spacing is close to the mesh spacing).
  These appear both in the `sigm` plot and in the MLO bands
- In the conventional QSGW column, blue appears even at iterations 9–10. It does not appear in the 2 columns of MLO-QSGW at the same 6³

### 3.3 Convergence

![Band change per iteration](../../manual/mlo_gwsc/liti2o4_conv.png)

Change of t2g per iteration (left max, right rms). 9³ solid lines, 6³ dashed lines.

| rms [meV] | 7 | 8 | 9 | 10 |
|---|---|---|---|---|
| **9³ MLO bands** | 7.4 | 3.2 | 5.2 | **3.2** |
| **6³ MLO bands** | 6.9 | 2.9 | 6.2 | **3.9** |
| conventional 6³ | 6.7 | 6.6 | 4.7 | 6.5 |

- Both meshes come down to the same level as the conventional one (rms 3–7 meV). In all 10 iterations all checks of §2.3 passed
- Without mixing $\Sigma^{\rm MLO}$ (equivalent to β=1), rms jumps to 120–150 meV at the 3rd iteration (35 meV for the conventional one).
  With mixing it decreases in the same way as the conventional one (31 meV at the 3rd iteration)
- **The movement of the bands when the mesh goes from 6³ to 9³** (iteration 10) is rms 13 meV for the MLO bands and 27 meV for the `sigm` plot.
  **Interpolation in the MLO representation is less sensitive to the mesh**

### 3.4 Difference from the conventional method

The difference between the MLO bands at iteration 10 and conventional 6³ is rms 46 meV (mean +30) at 6³ and rms 49 meV (mean +34) at 9³.
At 6³ it is not a uniform shift: the 2 occupied bands are 7–13 meV lower and the unoccupied t2g are 22–61 meV higher. The total t2g width is almost the same, 1505 vs 1516 meV.
**They differ by 43 meV even on the mesh points**, so this is a difference of state, not of interpolation.
However, the reference for comparison uses `pwmode = 1` and the MLO side `pwmode = 11`. The `pwmode` difference in LDA (rms 5.6 meV) alone cannot explain it, but
the possibility that it was amplified by QSGW remains. **Whether it is a difference of the $\Sigma$ representation or of `pwmode` cannot be separated without running a conventional chain with `pwmode = 11`** (not done).

### 3.5 Up to 40 iterations (2026-09-29)

With the settings of §3.1 (but with the code from 2026-09-27 on: window reference is the value of each iteration's SCF, $\chi_0$ with `t_tetrakbt = -992.4`, `t_sigmaw = 1000`),
we ran from LDA `gwsc 10 --prec=tf32`, then continued one iteration at a time up to 40 iterations. The table of settings and all figures and numbers are in ecalj's `Samples/kBT/LiTi2O4/README.md`.

![MLO bands of t2g of LiTi2O4, 6³ and 9³, iterations 10, 20, 30, 40](../../manual/mlo_gwsc/liti2o4_summary_t2g.png)

**Fig. 3.5-1**. MLO bands of t2g (Γ→X). Left 6³, right 9³. From the top: LDA, iterations 10, 20, 30, 40<span class="clfix">. Below, $x$ is the position along Γ→X (Γ = 0, X = 1)</span>

**Table 3.5-1**. Maximum change of the MLO bands from the previous iteration (meV), and a measure of the waviness of the lowest band (maximum residual of a fit of the correction from LDA with cos up to 3rd order)

| | Iteration | change of t2g | change of the 2 occupied bands | waviness |
|---|---|---|---|---|
| 6³ | 10 | — | — | 12.7 |
| 6³ | 20 | 13.6 | 8.7 | 26.0 |
| 6³ | 40 | 3.3 | 0.8 | 26.7 |
| 9³ | 10 | — | — | 17.2 |
| 9³ | 20 | 9.6 | 2.1 | 14.4 |
| 9³ | 40 | 5.1 | 0.2 | 12.9 |

- **Precision**: the difference between the MLO bands of tf32 and fp32 with the same input is at most 0.8 meV at iteration 10 of 6³ and at most 1.2 meV at iterations 10–12 of 9³ (t2g).
  fp64 and fp32 differ by at most 1 meV at the 1st iteration of 6³. The time of one iteration at 9³ is 21 min with tf32 and 47 min with fp32
- **6³**: smooth at iteration 10, but from 20 a shoulder forms at $x$ = 0.35–0.45 of the lowest band and remains at 40. The bulge of the band in the middle of t2g
  is largest at iteration 22 with 131 meV, and decreases to 49 meV at 40 (still decreasing by 2–3 meV per iteration)
- **9³**: at iterations 10–15 the lowest band just outside Γ is 6–12 meV below Γ, but this becomes 1.5 meV at 40.
  The interpolation curve connecting the mesh points still looks slightly wavy at 40 iterations
- At the mesh points common to 6³ and 9³ (Γ and $x$ = 2/3), the difference of the lowest band at 40 iterations is 4 meV or less. Between mesh points it is at most 37 meV

---

## 4. Pitfalls

| Mistake | Why it is wrong / what is correct |
|---|---|
| not mixing $\Sigma^{\rm MLO}$ (`ECALJ_MLO_MIX=0`; with code up to 2026-09-29, forgetting `ECALJ_MLO_MIX=1`) | overshoots in the early iterations (§3.3) |
| taking the `job_band` plot for the MLO result | it is the conventional `sigm` interpolation. Draw the MLO bands by the method of §1.7 |
| drawing MLO bands with `job_band --mlo` or `job_mlo --mlofreeze` | invalid, as in the table of §1.7 |
| applying check values for another mesh (line counts of the marks, time taken) | normal is judged abnormal. The values are determined by the mesh and the number of processes (§2.3)|

## 5. Unresolved and future work

- Run a **conventional chain with `pwmode = 11`** for 10 iterations and separate the difference of §3.4 into the difference of representation and of `pwmode`
- **Antiferromagnetic symmetry (`symgrpaf`)**: from 2026-09-30, QSGW works in the form where the GW programs compute both spins (ecalj `Samples/AFsymmetry`).
  Combination with `--mlo` has not been checked
- The figures and numbers of §3.2–§3.4 are from the code of 2026-09-26<span class="clfix"> before the window fix of that evening</span> (window reference is the $E_F$ at the start of the chain). The MLO bands of the current code after 9 iterations differ from them by
  rms 8 meV (max 23 meV), an almost uniform lowering of 7–10 meV (because the window reference moved by 0.55 eV, they diverge from the 2nd iteration). §3.5 is the current code
- `--prec-final=fp32:N` (fp32 only at the end) has not yet been run and checked (the difference between tf32 and fp32 is about 1 meV, as in §3.5)
- It is not known whether the shoulder at 6³ (§3.5) is a property of the converged solution or disappears with more iterations
- Effect of the number of MLOs: 76 (`gwinit` default), 126, and all MTOs (154) have not been compared with $\Sigma^{\rm MLO}$ mixing
- Whether adding empty spheres reduces oscillations and stabilizes (unverified)
- Extend the search for $\tilde\chi$ to symmetry-equivalent $k$ (then `job_band --mlo` becomes usable at the mesh points).
  To draw PMT bands via MLO after convergence, rebuild $z$ at each $k$ with $H$ including $\Sigma$ and make it self-consistent
- Speed: the remaining 80% of `hgw` at 9³ is the matrix products of Σc, which already run near the single-precision capability of this GPU ([ecaljgpu](./ecaljgpu)).
  Building $\Sigma^{\rm MLO}$ directly from the vertex $\langle\tilde\chi_n|\psi_{i'}M_I\rangle$ would make one side of the matrix product $n_{\rm MLO}$, but
  for LiTi₂O₄ ($n_{\rm MLO}$ = 126, `ntq` = 159) it reduces it only by about 10%
