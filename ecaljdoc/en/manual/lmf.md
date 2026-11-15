# `lmf,lmfa,lmchk`

> English version of the page [manual/lmf](../../manual/lmf.md). The pages at the top of ecaljdoc are the reference; this one was made from manual/lmf.md of commit b55836332 (source blob 43ef5dca1bfc).

The input file is `ctrlg.<sname>.toml` (`ctrl.<sname>` in the old format).
`<sname>` is a **positional argument** passed to the program (for example, `lmf si` reads `ctrlg.si.toml`). On this page,
`<sname>` denotes the material name chosen with `ctrlgenToml.py` or `Legacy2toml.py`.

## lmfa: spherically symmetric atom calculation
Example:
```
lmfa si
```
It finishes in a few seconds.
```
lmfa si |grep conf
```
shows the electron configuration. `lmfa` computes the spherically symmetric electron density of the atoms contained in the unit cell.
From their superposition it makes the initial electron density for `lmf`.
`lmfa` also writes the logarithmic derivatives of the radial functions at the MT radii to `atmpnu.*.<sname>` ([`[ham].readp`](#readp-option)).

## lmchk: checking the crystal symmetry
`lmchk` displays the space-group symmetry that `lmf` recognizes. It also determines the MT radii. The notation for space-group operations is in [SYMGRP](#symgrp).

### MT radii determined by lmchk
The ideal sphere radius is the one that best approximates a potential that is spherically symmetric inside the MT and flat outside.
`lmchk` contains one method to determine reasonable initial values of the MT radii.
In all directions from an atom toward relatively nearby atoms, it computes the (electrostatic) potential made by overlapping the free-atom densities,
and takes the MT radius up to the first maximum of the potential that appears in any direction.
This is a fairly reasonable estimate, in the sense that the potential inside becomes nearly spherically symmetric.
For a perfectly symmetric bond, the potential maximum falls exactly at the bond midpoint, so the two spheres just touch and the potentials are also equal.
To obtain these radii, run `lmchk --getwsr`. `ctrlgenToml.py` does this automatically.
* `ctrlgenToml.py` also sets an upper limit on the radii of some atoms (`r_upper_limit`). The limits and the reason for them (spurious levels) are in [MT radii](./mtradius).

* For molecules (dimers in particular), very small MT radii must be used. We have confirmed that PMT works well even for such systems, but it needs more investigation.

## lmf: solving the one-body problem
Performs LDA/GGA and LDA+U band calculations. The QSGW potential term can also be added through `sigm.<sname>`.

Example:
```
mpirun -np 4 lmf si [options] > llmf
```
This is for the case where `ctrlg.si.toml` exists.

### Options
* Options to stop midway: `--quit=show`, `--quit=ham`, `--quit=mkpot`, `--quit=dmat`, `--quit=band` (see [command-line options](./cmdopts)).
* `--tdos`: total density of states. Usually `job_tdos` is used.
* `--ctrlg:<section.key>=<value>`: changes a key for one run only. See [command-line options](./cmdopts#run-time-ctrlg-overrides).

There are also options for electron-density plots, BoltzTraP, and so on (not documented yet).

### Input files

- `syml.<sname>`: symmetry lines for band plots (usually `job_band` calls `lmf` internally)
- `sigm.<sname>`: part of the QSGW calculation
- `atmpnu.{1,2,3}.<sname>`: made by `lmfa` (when `[ham] readp = true`)

#### `[ham].readp` and `[ham].pnufix` {#readp-option}
With `readp = true`, `lmf` reads the logarithmic derivatives at the MT radii (valence $P$ and local-orbital $P_z$) that `lmfa` wrote to `atmpnu.*`.
With `pnufix = true`, it keeps them fixed during the iterations. Both are defaults written by `ctrlgenToml.py`. Hence the basis does not change during the iterations.
For details, see "On the radial equation" in [Representation of wave functions used in GW](../theory/gw_basis).

### Output files
* `rst.<sname>`
  The restart file, containing the electron density and related information.
  If `ctrlg.<sname>.toml`, `atmpnu*`, `sigm*`, and `rst*` are present, the `lmf` calculation can be restarted
  (LDA+U also uses `dmats*`).

* `save.<sname>`
  A minimal history of `lmf`, `lmfa`, and `lmchk`; it contains the total energy at each self-consistency iteration of `lmf` ([below](#save-sname)).

* `__mixm.<sname>`
  The file for electron-density mixing. It holds pairs of input and output densities from previous iterations and is used for Anderson or Broyden mixing to speed up convergence to self-consistency.
  `lmf` discards this file from a previous run at the start (the option `--keepmixm` keeps it). This is so that it is not used for a new calculation, such as one with a changed lattice constant.
  If the history of a run that ended in a nonmagnetic solution remains, a new spin-polarized calculation is pulled into the nonmagnetic solution at the first Broyden step (Fe: 2.13 → 0.02 μB).
  Only `rst.<sname>` is carried over. When changing `nspin` on top of an existing `rst`, delete `rst.<sname>` first.

**Main source**: [`SRC/subroutines/main_lmf.f90`](https://github.com/tkotani/ecalj/blob/main/SRC/subroutines/main_lmf.f90)

## ctrlg.&lt;sname&gt;.toml — the current input format (2026-05–)

Since 2026-05, `lmf`, `lmfa`, and `lmchk` read **only TOML**.
A single file `ctrlg.<sname>.toml` contains everything that used to be in `ctrl.<sname>`, plus the GW settings and the product-basis section that used to be in `GWinput`.
At the end of the file there are per-atom product-basis tables. In an old working directory, make it with `Legacy2toml.py <sname>`.

The `[product_basis]` tables at the end of the file are for GW; `gwinit` writes them and they are usually not edited by hand. See [GW settings: product basis](./gwinput#product-basis).

#### Console output {#reading-the-output}
The console output is now mainly for debugging.
Still, some values have to be read from the console output (we intend to fix this).
From the console output you can check the progress of convergence, band energies, the Fermi energy, and whether `sigm` and `rst` were read correctly.

#### `save.<sname>` and `log.<sname>` {#save-sname}
* `save.<sname>` keeps the history of invocations of `lmf`, `lmchk`, and `lmfa`. In addition, it writes one line with the total energy for each `lmf` iteration.
  The first character of each line indicates the status. c: converged, x: reached the iteration limit without converging, i: in the middle of the iterations,
  h: first iteration (Harris energy from the superposition of free-atom densities), C: converged structural relaxation (`svflg='cxhiC'` in `m_lmfp.f90`).
* The `--ctrlg:` overrides used are also recorded.

* Two total energies are printed, Kohn-Sham and Harris-Foulkes. [The two should be nearly the same](#q-should-the-harris-foulkes-and-hohenberg-kohn-sham-functionals-agree-at-self-consistency),
  but in large systems they may differ somewhat. Use either one.
* In QSGW mode (`lmf` using the `sigm` file), the total energy has no meaning. The `lmf` iterations inside `gwsc` also print total energies, but look at them only as a measure of convergence.

* Log file
  `log.<sname>`, made by `lmf`, `lmchk`, and `lmfa`, is currently for debugging.

### File structure (sections) {#file-structure-sections}

```toml
# Top-level scalars
symgrp  = "find"        # space-group spec: 'find' = auto-detect
verbose = 35            # console verbosity
time    = [0, 0]        # CPU timing log: [depth, on-the-fly]

[struc]   alat / plat                              (no nbas / nspec: the counts are those of the tables below)
[[site]]  atom / pos / xpos / af / relax           (one table per atom)
[[spec]]  atom / z / r / rsmh / eh / lmx / ...     (one table per species)
[bz]      nkabc / metal / tetra / npts / ...
[iter]    nit / mix / b / conv / convc / umix / tolu
[ham]     nspin / rel / so / phispinsym / xcfun / gmax / pwmode / pwemax / oveps / ...
[esm]     boundary / origin / shiftmode / zb / potential / field   (slabs only)
[gw]      n1n2n3 / QpGcut_psi / HistBin_dw / iSigMode / niw / t_sigmaw / t_tetrakbt / ...
          QforEPS / QforGW (multi-line q lists)
[mlo]     mlo_method / mlo_delta / mlo_w / mlo_nkabc (required for job_mlo) / mlo_lm (formerly Worb)
[blocks]  QPNT, QforEPSL, hrotr (raw multi-line blocks with no better home)
[product_basis]   pb_tolerance / pb_lcutmx + nlx / valence / core   (always the last section)
```

Do not write `nbas` / `nspec` in `[struc]`. See [Input files](./toml_migration#an-older-ctrlg-sname-toml-ctrlg-update-py).

### Example: fcc Cu (1 atom, nonmagnetic)

The whole of `ctrlg.cu.toml` in `Samples/EPS/EPS_Cu/` is shown below. The `#` comments at the ends of the lines give the role of each key.
Copy it as a starting point for an elemental fcc metal and adapt `[spec]`, `[bz]`, and `[ham]` to each case.
Since this file was converted from an old input of the sample, some values differ from those written by the current `ctrlgenToml.py` or from the defaults
(for example, `pwmode = 1` here while the default is 11; `[iter] mix = "A2"`, `nit = 30`; `[gw] n1n2n3`, `QpGcut_*`, `t_sigmaw = 262`, `t_tetrakbt = 0`).
For the current template values, see [GW settings](./gwinput) and the comments in a newly made `ctrlg.<sname>.toml`.

```toml
# Top-level scalars
symgrp  = "find"  # 'find' = auto-detect space group from lattice + species
verbose = 35      # 31 default, 35 verbose, 41+ debug
time    = [0, 0]  # CPU timing log: [depth, on-the-fly]

[struc]
alat  = 6.798     # lattice constant (a.u.)
plat  = [[0.0, 0.5, 0.5],   # primitive vectors in units of alat
         [0.5, 0.0, 0.5],
         [0.5, 0.5, 0.0]]

[[site]]
atom = "Cu"
pos  = [0.0, 0.0, 0.0]

[[spec]]
atom   = "Cu"
z      = 29
r      = 2.33    # MT radius (a.u.)
rsmh   = [1.17, 1.17, 1.17, 1.17]   # smoothing radii per l
eh     = [-1.0, -1.0, -1.0, -1.0]   # tail energies
rsmh2  = [1.17, 1.17, 1.17]         # second basis group
eh2    = [-2.0, -2.0, -2.0]
lmx    = 3       # basis l-cutoff
lmxa   = 4       # augmentation l-cutoff
nmcore = 1       # 1: spin-averaged core; 0: spin-polarized core
kmxa   = 5       # radial expansion order at tail sites

[bz]
nkabc  = [12, 12, 12]   # k-mesh divisions; 1 entry => isotropic
metal  = 3              # 0 insulator, 3 tetra+broadening (safe metal)
savdos = true           # write dos.tot.* file
npts   = 2001
dosmax = 1.5

[iter]
nit   = 30        # max self-consistency iterations
mix   = "A2"      # B3 = Broyden hist=3 (default); A2 = Anderson hist=2
conv  = 1e-05     # max dE between iterations (Ry)
convc = 1e-05     # max d(rho_out - rho_in) (Ry)
b     = 0.2       # mixing ratio

[ham]
nspin       = 1   # 1 nonmag, 2 spin-polarized (then spec.mmom seeds)
rel         = true
so          = 0   # 0 none, 1 L.S, 2 Lz.Sz
gmax        = 12.0   # real-space mesh G-cutoff (Ry); alt: ftmesh
frzwf       = false
xcfun       = 1   # 1=VWN, 2=Barth-Hedin, 103=PBE-GGA
forces      = 0
scaledsigma = 1.0   # QSGW mixing: 1.0 full, 0.8 = QSGW80
oveps       = 0.0
pwmode      = 1   # 0=MTO 1=PMT 2=APW 11=PMT(|q+G|, default) 12=LAPW
pwemax      = 3.0  # APW cutoff (Ry)

# === GW driver (used by lmf --jobgw=N and gwsc) ===
[gw]
n1n2n3        = [10, 10, 10]   # GW BZ mesh (typically smaller than [bz].nkabc)
QpGcut_psi    = 3.0     # |q+G| cutoff for eigenfunctions
QpGcut_cou    = 2.0     # |q+G| cutoff for Coulomb / W
unit_2pioa    = false   # false: a.u.;  true: 2*pi/alat
HistBin_dw    = 1e-05   # bin width on real-axis at omega=0
HistBin_ratio = 1.03
iSigMode      = 3       # self-energy mode (3 = standard QSGW)
niw           = 10      # # of imag-axis frequencies
deltaw        = 0.02
t_sigmaw      = 262     # (K) Fermi-Dirac width of the Sigma levels
t_tetrakbt    = 0       # (K) required. 0: no smearing of chi0; T > 0: finite-T tetrahedron; -T: Gaussian on Im chi0
QforEPS = """           # q points for eps (a.u. if QforEPSau = true)
 0 0 0.00050
 0 0 0.00100
 0 0 0.00200
"""
QforGW = """            # q points for the one-shot GW driver (gw_lmfh)
 0.0 0.0 0.0
 0.1 0.0 0.0
"""

# === MLO (muffin-tin based localized orbitals; see manual/mlo) ===
[mlo]
mlo_method = 4          # theta = sigma((eps - ecut_j)/mlo_w), ecut_j = max(CBM + mlo_delta, eps^MTO_j)
mlo_delta  = 2.0        # (eV) how far above the band edge the model must be accurate
mlo_w      = 2.0        # (eV) width of the fall-off; the knob if the residual is too large
mlo_lm     = """        # which lm channels of which atom make the model (1=s, 2-4=p, 5-9=d, 10-16=f)
  1 Cu   5 6 7 8 9
"""

[product_basis]
pb_tolerance = [0.001]  # drop near-linear-dep products (default 1e-3)
pb_lcutmx    = [4]      # max l-cutoff per atom
# ----------------------------------------------------------------
# nlx [iatom, l, nnvv, nnc] ... (per-atom tables written by gwinit; see manual/gwinput)
nlx = [
  [1, 0, 2, 3],
  [1, 1, 2, 2],
  [1, 2, 2, 0],
  [1, 3, 2, 0],
  [1, 4, 2, 0],
]
# ----------------------------------------------------------------
valence = [
  [1, 0, 1, 1, 1],   # 4s_phi
  [1, 0, 2, 0, 0],   # 4s_phidot
  # ...
]
# ----------------------------------------------------------------
core = [
  [1, 0, 1, 0, 0, 0, 0],   # 1S
  # ...
]
```

(The q-point lists `QforEPS` / `QforGW` are placed at the end of `[gw]` as multi-line strings, as in the example above.) That is the whole input.

The layout of the generated file (blank lines, `# ----` rules, order of sections) is described in [Input files](./toml_migration#file-layout).

#### Differences by physics (starting from Cu)

- **Spin polarization**: set `[ham].nspin = 2` and write `mmom` in each `[[spec]]`. See [Other explanations](./UsageDetailed#spin-polarized-case-without-soc). Start from a clean directory (`__mixm.<sname>` above).
- **SOC**: `[ham].nspin = 2` with `[ham].so = 1` (L·S) or `2` (Lz·Sz only). See [Other explanations](./UsageDetailed#spin-orbit-coupling) and `Samples/MLOsamples/FeSoc/`.
- **Insulators and semiconductors** (Si, GaAs): `[bz].metal = 0`, `[bz].tetra = true`. `[bz].nkabc` can be smallish, and `[ham].pwemax` of about 3 is enough.
- **LDA+U**: see [LDA+U](#lda-u-toml) below.
- **Relaxation of atomic positions**: see [Forces and structural relaxation](./relax).
- **GW BZ mesh**: `[gw].n1n2n3` ≤ `[bz].nkabc`. See [GW settings](./gwinput).

### Run-time overrides with `--ctrlg:`

To change a key for one run only, add `--ctrlg:<section.key>=<value>` (for example `lmf si --ctrlg:bz.nkabc=[8,8,8]`). The file is not rewritten. Details of the syntax, the handling of keys not in the file, and the retired old syntax (such as `-v<NAME>=<VAL>`; using it stops the program) are in [command-line options](./cmdopts#run-time-ctrlg-overrides).

The correspondence between the keys of the old `ctrl.<sname>` and the TOML paths is in [Input files](./toml_migration#legacy-toml-map).
The `[esm]` section for ESM (slabs with vacuum) is on the [ESM](./esm) page.

### Convergence of lmf

If `lmf` itself does not converge, make `[iter] b` smaller, such as 0.015, 0.01, 0.005. Mixing of the self-energy in `gwsc` (`[gw].mixbeta`, `mixpriorit`) is described in [gwsc](./gwsc#mixbeta-and-mixpriorit).

### Checking the MTOs in use
The MTOs in use are shown near the beginning of the console output (for GaAs).
 ```
 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol= 1.00E-06
 spec      l    rsm    eh     gmax    last term    cutoff
  Ga       0*   1.13  -1.00   6.579    1.19E-06    1459
  Ga       1*   1.13  -1.00   7.028    1.26E-06    1807
  Ga       2*   1.13  -1.00   7.475    1.09E-06    2109
  Ga       3    1.13  -1.00   7.920    1.06E-06    2637
  Ga       0*   1.13  -2.00   6.579    1.19E-06    1459
  Ga       1*   1.13  -2.00   7.028    1.26E-06    1807
  Ga       2    1.13  -2.00   7.475    1.09E-06    2109
  As       0*   1.18  -1.00   6.300    2.13E-06    1243
  As       1*   1.18  -1.00   6.720    1.26E-06    1471
  As       2*   1.18  -1.00   7.140    1.37E-06    1837
  As       3    1.18  -1.00   7.558    1.05E-06    2229
  As       0*   1.18  -2.00   6.300    2.13E-06    1243
  As       1*   1.18  -2.00   6.720    1.26E-06    1471
  As       2    1.18  -2.00   7.140    1.37E-06    1837
 ```

#### Local orbitals (`pz`) and the `ecore>evalence` message
When the core spillout is not small (larger than about 0.05), it is better to use `pz` (local orbitals) and treat that core as valence.
Bi 4d is an example; set `pz = [0.0, 0.0, 4.9]` or similar.

The message `ecore>evalence` means that a core level is higher than a valence level. For safety, this is not allowed.
Compare the ECORE file with the valence levels in the log or console output.

### Molecules
`ctrlgenToml.py <sname> --systype=molecule` makes a template for a molecule in a box (Γ point only, `[bz] nkabc = [1, 1, 1]`).
For molecules, write the following in `[bz]`.
```toml
tetra = false
n = -1        # Negative is the Fermi distribution function. w gives temperature.
w = 0.001     # w=0.001 corresponds to T=157K as shown in console
```
When there is a magnetic moment, `fsmom` (n_up − n_down) is also needed (`fsmommethod = 1`).
Sample: `Samples/AtomDimer/N2` (N₂ and an N atom in a box; bond length and binding energy).
The template written by `ctrlgenToml.py --systype=molecule` has `tetra = true`, `w = 0.01`, which differ from the values above (those of the sample `Samples/AtomDimer/N2`).
Which is better has not been checked yet (2026-11-14).
An example of the PMT method applied to molecules from H₂ to Kr₂: T. Kotani and H. Kino, J. Phys. Soc. Jpn. **82**, 124714 (2013),
[doi:10.7566/JPSJ.82.124714](https://doi.org/10.7566/JPSJ.82.124714)
([pdf](../../presentations/KotaniKino2013PMTMolecule.pdf)).

---

## LDA+U {#lda-u-toml}

LDA+U is selected per species by writing the vectors `idu`, `uh`, `jh` in each `[[spec]]` table. The index is `[s, p, d, f]`.

```toml
[[spec]]
atom = "Er"
z    = 68
r    = 3.00
# ... rsmh, eh, ... ...
idu  = [0, 0, 2, 2]          # 0 = none, 1 = AMF, 2 = FLL.  +10 disables when sigm.* exists
uh   = [0.0, 0.0, 0.10, 0.632]    # Hubbard U (Ry)
jh   = [0.0, 0.0, 0.00, 0.055]    # exchange J (Ry)
```

Values of `idu[l]`:
- `0` — no U for this l
- `1` — **AMF** (Around-Mean-Field), Petukhov, PRB 67, 153106 (2003)
- `2` — **FLL** (Fully Localized Limit), Liechtenstein, PRB 52, R5467 (1995)
- `12` / `21` / ... — adding 10 to the same value means "do not use LDA+U when `sigm.<sname>` exists"
  (a QSGW calculation with `sigm` becomes plain QSGW).

`dmats.<sname>` (density matrix) is read and written automatically. It is the LDA+U counterpart of `rst.<sname>`. It is made by the first LDA calculation without U.
Writing `occnum.<sname>` by hand (one line of `2l+1` numbers per spin) lets you specify the initial occupations. An ErAs example and pitfalls of the spin configuration are in
[old notes on LDA+U](./lmf_legacy_ctrl#lda-u).

The density matrix is mixed **separately** from the charge mixing, by linear mixing.

```toml
[iter]
umix = 0.2          # dmat mixing beta (0..1)
tolu = 0.0          # rms tol where dmat mixing stops (0 = always mix)
```

Check the density matrix once with `lmf <sname> --quit=dmat`. It stops after the initial setup and prints the matrix without writing the restart file.

## Core hole {#core-hole-toml}

By partially changing the core occupation, Slater's transition-state idea is realized (useful for XAS or for removing 4f). It is set per species.

```toml
[[spec]]
atom    = "Nd"
z       = 60
r       = 3.00
p       = [0.0, 0.0, 0.0, 5.2]    # treat 4f as core, 5f as valence
c-hole  = "4f"                     # core channel to perturb
c-hq    = [-11.0]                  # excess electron charge in that channel
# c-hq    = [-7.0, 7.0]            # optional 2nd entry = spin moment
```

Charge accounting:
```
Q(spin1) = full_core_count + c-hq[1]/2 + c-hq[2]/2
Q(spin2) = full_core_count + c-hq[1]/2 - c-hq[2]/2
```
A negative `c-hq[1]` removes electrons (a hole). The second entry `c-hq[2]` polarizes the core
(for example, the Gd 4f technique: `c-hq = [-7, 7]` fixes the 4f core at 7 up and 0 down, and 5f stays valence).

Check charge neutrality with `lmfa <sname> | grep -i 'add core hole'`. Examples for Gd, Nd, and CrN are in
[old explanation of core holes](./lmf_legacy_ctrl#core-hole-options).

## Local orbitals {#local-orbitals-toml}

The local orbital `pz` is the counterpart of the continuous principal quantum number `p`. Both are per-species, per-`l` vectors in `[[spec]]`.

```toml
[[spec]]
atom = "Ga"
# p[0..3] = s,p,d,f valence principal quantum numbers (continuous);
#          fractional part 0.5 = mid-band, ->0 = bottom, ->1 = top.
p   = [0.0, 0.0, 4.0, 0.0]       # 4d "valence-like"
# pz[0..3] = semicore (lower n) or high-lying (higher n) local-orbital prinq;
#           fractional part ~0.9 = semicore, ~0.3 = high-lying.
pz  = [0.0, 0.0, 3.9, 0.0]       # 3d treated as semicore local orbital
# add 10 to pz[l] to attach a smooth-Hankel tail (extended local orbital).
```

Notes:
- The integer parts of `pz` and `p` determine the principal quantum numbers of the radial functions (for semicore `int(pz) = int(p) − 1`, for high-lying levels `int(pz) = int(p) + 1`).
- The fractional part is typically `0.9` for semicore and `0.3` for high-lying levels. The smaller the fractional part, the lower the linearization energy.
- `pz[l] + 10` (that is, writing `pz = [0, 0, 13.9, 0]`) attaches a smooth-Hankel tail (extended local orbital). Use it when a deep state spills out into the interstitial region.
- Check the configuration with `lmfa <sname> | grep -i conf`. It shows which states are treated as core, valence, and local orbitals.

Common pitfalls:
- When an ordinary (non-extended) local orbital is used for a semicore state, its representation outside the MT sphere relies on the valence envelope functions. Make at least one `rsmh` of that l small (about 1) and `eh` fairly negative (about −1).
- If the fractional part of `pz` is too close to 1.0, matching of the extended local orbital may fail ("`mtchre: failed to match phi'/phi`"). Make it smaller (for example `0.98 → 0.85`).
- The `p` with the larger principal quantum number is fixed automatically, and only the deepest state can float (controlled by `idmod`).

The old discussion of `mtchre` for smooth Hankels is in [old notes on local orbitals](./lmf_legacy_ctrl#old-note-for-local-orbitals).

## SYMGRP / symgrp — specifying the space group {#symgrp}

In TOML, the space group is specified by the top-level string `symgrp`.

```toml
symgrp     = "find"           # auto-detect from lattice + species (default)
symgrpaf   = "find"           # AF symmetry (optional): spglib finds it from the af labels; generators only for the old search
```

`symgrp` takes the string `"find"` (automatic detection) or a whitespace-separated list of generators. The generators are the same as in the old `SYMGRP`.

- `R<n>X`, `R<n>Y`, `R<n>Z`, `R<n>D` — `n`-fold rotation about X, Y, Z, (1,1,1)
- `M<axis>` — mirror
- `I` — inversion
- `<axis> = X | Y | Z | D | (nx,ny,nz)` — a Cartesian axis, or a direction given directly
- `*` takes a product, e.g. `I*R4X`
- Translation: `r6z:(0,0,0.4778973)` (Cartesian), `r6z::(0,0,1/3)` (fractional, multiples of `plat`)

Example (cubic):
```toml
symgrp = "R4X MX R3D"        # 4-fold around X, mirror in X, 3-fold around (1,1,1) ⇒ 48 symops
```

When generators are written together with `find` (the old form, still accepted), spglib finds the whole group and the generators must be among its operations (otherwise `lmf` stops). The generators do not change the group.
```toml
symgrp = "R4Z*I MX*I R3D find"   # zinc blende: the same 24 operations as symgrp = "find"
```

To turn off all symmetry:
```toml
symgrp = "e"
```

`lmchk` prints the recognized symmetry operations at the beginning of the console output. Use this to check that the symmetry was read as intended.
When the potential is local, inversion is implicitly included in DFT calculations. See the note about inversion at the top of the console output.

### Symmetry from spglib: `symmetry.<sname>.json` (2026-10) {#spglib-symmetry}

With `symgrp = "find"` (the default), `lmf`, `lmchk`, and the MLO programs find the space group with spglib. spglib is built into ecalj (`SRC/external/spglib`, version 2.6.0; Python is not needed).
The operations $x' = Rx + t$ (in fractional coordinates of `plat`, with integer $R$) are written to `symmetry.<sname>.json` in the working directory, together with the space-group number and symbol and the structure they belong to.
The next run reads this file. If the structure differs (alat, plat, species, `af` labels, positions, including `--ctrlg:` overrides), the operations are found again and the file is rewritten.
The console output tells which happened: `mksym: N operations, read from ...` or `found by spglib, <symbol> (<number>), written to ...`.

* For a supercell (a cell larger than the primitive cell, e.g. the 8-atom cubic cell of Si), the pure translations are also included in the operations. For Si8 there are 192. The total energy, the MLO bands, and the QSGW quasiparticle energies agree with those obtained with the operations of the primitive cell.
* Antiferromagnets: writing `af = ±k` on the sites of an antiferromagnetic pair lets spglib's magnetic symmetry give the operations. Operations accompanied by time reversal exchange up and down spins. With `symgrpaf = "find"`, these operations put the calculation into antiferromagnetic mode (2026-10-02).
  Generators in `symgrpaf` (the old form, e.g. `"i:( 1 1 1 )"` for NiO) are also accepted, and they must be among these antiferromagnetic operations. Use it together with `pwmode = 11`.
  The procedure for NiO and others is in [Other explanations](./UsageDetailed#antiferro-symmetry-without-soc).
* With a `symgrp` that lists generators (to lower the symmetry on purpose, e.g. `"r4z"` for a 4f state with orbital order), the group generated by those generators is used as before, and no file is written. In this case `symgrpaf` needs generators (`"find"` cannot be used).
* The old ecalj search (searching for lattice operations and the operations that keep the crystal invariant; `ECALJ_SYMFIND=ecalj`) was removed on 2026-10-02. Only spglib does the search. All the input forms above are still accepted.
* `symfind.py <sname>` (the Python spglib) makes the same file outside `lmf` as well.

## Q: Should the Harris-Foulkes and Hohenberg-Kohn-Sham functionals agree at self-consistency? {#q-should-the-harris-foulkes-and-hohenberg-kohn-sham-functionals-agree-at-self-consistency}
(by Mark van Schilfgaarde)
A: Yes. If they do not agree, either self-consistency has not been reached, or the tolerance for evaluating matrix elements is too loose. A common cause is that the integration mesh for the interstitial region is not fine enough.
The following is part of the output when a mesh of 14×14×14 divisions is used for the interstitial density.

```
 GVLIST: gmax = 4.702 a.u. created 973 vectors of 2744 (35%)
         mesh has 14 x 14 x 14 divisions; length 0.588, 0.588, 0.668
 SGVSYM: 170 symmetry stars found for 973 reciprocal lattice vectors
...

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol= 1.00E-05
 spec      l    rsm    eh     gmax    last term    cutoff
  Te       0*   1.30  -0.10   5.220    9.30E-05     973*
  Te       1*   1.30  -0.10   5.597    4.36E-04     973*
  Te       2    1.30  -0.10   5.976    2.04E-03     973*
  Te       0*   1.40  -1.00   4.847    2.11E-05     973*
  Te       1    1.40  -1.00   5.182    9.90E-05     973*
  E        0*   1.50  -0.30   4.524    1.06E-05     847
  E        1    1.50  -0.30   4.823    2.01E-05     973*
```

This table shows how well each basis orbital is converged in the plane-wave expansion used to make the matrix elements of the interstitial potential.
Only the E s orbital is converged below the specified tolerance ($10^{-5}$). Convergence of the Te d orbital is particularly poor, but since Te d hardly participates in bonding, it does not affect the energy much.
In this case (3 Te atoms), the difference between the HF and HK total energies is about 10⁻⁴ Ry. The following table shows the convergence of the HF energy with the number of divisions, and the difference between the HK and HF energies.
```
  n       ehf         ehf-ehk
  14 -14561.140624    0.000089
  15 -14561.133966    0.000017
  16 -14561.132309    0.000005
  20 -14561.131170    0.000000
```

## Old documents

Most of the old explanation of `Category_Token` and the old notes on physics (APW tests, the discussion of `mtchre` for smooth Hankels, old SPEC notes, old notes on the tetrahedron method) have been moved to a separate page.

- [Explanation of the old `ctrl.<sname>` (Category_Token)](./lmf_legacy_ctrl)

Look up the *meaning* of an old key on that page, and convert it to TOML with [correspondence between old keys and TOML paths](./toml_migration#legacy-toml-map).
