# gwsc: a script to run QSGW calculations

> English version of the page [manual/gwsc](../../manual/gwsc.md). The pages at the top of ecaljdoc are the reference; this one was made from manual/gwsc.md of commit 296542c3a (source blob db175b30078f).

gwsc is the script to run QSGW calculations.
A QSGW calculation is run by calling several Fortran executables.
The input is read from `ctrlg.<sname>.toml` (GW settings are in `[gw]`; the keys are explained in [gwinput](./gwinput)).

## Usage
```
gwsc -np NP [-np2 NP2] [--gpu] [--prec=tf32|fp32|fp64] [--prec-final=<prec>[:N]] [--no-tetwt-helper] [--mp] [--fp32] [--mlo] [--dumpW] [--fullstdo] nloop extension [Options]
```

### `-np NP`
Number of MPI processes.

### `-np2 NP2`
Specify only when using the GPU version. Number of MPI processes for the calculations run on GPUs.
Usually set to the number of available GPUs.

### `--gpu`
Specify only when using the GPU version.

### `--prec=tf32|fp32|fp64` (2026-09)
Precision of the GW part. `fp64` is the double-precision version; `fp32` is the mixed-precision version with single-precision products in FP32 (same as `--mp --fp32`);
`tf32` is `fp32` except that only the last product of Σc uses a 10-bit mantissa for the inputs (TF32 or FP16, chosen by a per-GPU table). The legacy `--mp` alone is the same as `tf32`.
Which method computes products and inverse matrices on the GPU is chosen automatically by a table measured at installation (`<bindir>/ecalj_linalg_policy.toml`).
Accuracy and speed are compared in [GPU manual § choosing the precision](./ecaljgpu#gpu-precision).
When given together with `--mp` and `--fp32`, they must not contradict (a contradiction stops the run).

### `--prec-final=<prec>[:N]` (2026-09)
Use a different precision only for the last N iterations of the run (default 1). Example: `gwsc 10 --prec=tf32 --prec-final=fp32:2 ...` runs 8 iterations with tf32 and the last 2 iterations with fp32.

### `--no-tetwt-helper` (2026-09)
With `--gpu`, do not use the helper run that computes the tetrahedron weights on the CPU beforehand (`hgw --tetwt_write`).
It is used by default ([GPU manual § tetrahedron weights](./ecaljgpu#gpu-tetwt-cpu)).

### `--mp`
Use the mixed-precision version. Since 2026-09-27 it means the same as `--prec=tf32`. The earlier meaning, and the cases where TF32 broke $W$, are in [GPU manual § mixed precision](./ecaljgpu#gpu-mp-fp32).

### `--fp32` (2026-06)
Used with `--mp`; same meaning as `--prec=fp32` (the products of Σc also in FP32). The GW1500 batch of ecalj_auto (`gwscconv --gpu --mp --fp32 --conv-tol 0.1`) uses this by default.

### `--mlo`
QSGW in which the self-energy is interpolated in the MLO representation (MLO-QSGW, under development). `[bz] nkabc` and `[mlo] mlo_nkabc` must equal `[gw] n1n2n3` (it stops otherwise). See [MLO-gwsc](./mlo_gwsc).

### `nloop`
Number of QSGW iterations.

### `extension`
The `<sname>` of `ctrlg.<sname>.toml`.

### `Options`
Additional options.
Additional options are passed as run-time arguments to all executables. A list of additional options follows.
TOML overrides for lmf (`--ctrlg:ham.so=1` etc., of the form `--ctrlg:<toml-path>=<value>`) are also written here.

------
#### Memory usage and MPI distribution
The screened Coulomb interaction $W$ is held in memory inside `hgw` (`hgw_gpu`, `hgw_mp_gpu` with `--gpu`) and is not written to files
(`--dumpW` writes it to `__WVR.<iq>`, `__WVI.<iq>` for analysis). The distribution of MPI ranks (groups of $q$ points, and the parallelization over $\omega$ and $k$ within them) is
determined automatically from the available memory and shown at `MPI layout:` in `lgw`. Tuning is done not with command-line options but with keys in `[gw]` of `ctrlg.<sname>.toml`.

- `KeepWV`: keep $W$ in GPU memory when calculating the self-energy (correlation part). Default `true` for the GPU version, `false` for the CPU version. With `false`, $W$ is sent from main memory frequency by frequency (less GPU memory but slower)
- `zmel_batch_gb`: size of the batch in which matrix elements are made at once ([gwinput](./gwinput))
- `mpi_worker_exch`, `mpi_worker_corr`: number of ranks assigned to one group of $q$ points in the exchange and correlation calculations. 0 (default) is automatic. Give a divisor of the total number of ranks

There are no options `--keepwv`, `--nwpara=X` or `--nb=X` (giving them stops the program; [Command-line options](./cmdopts)).

#### `--tetwtk`
If given, the tetrahedron weights between the coupled states are not kept in memory when calculating the polarization function. Use it when memory runs short in calculations with many $k$ points.

#### `--skipGS`
Used by `lmf --jobgw=1`.
In the GW calculation, the wave functions obtained by the DFT calculation (`lmf`) are re-expanded in basis functions different from those of `lmf`. Gram-Schmidt orthonormalization is applied to the re-expanded wave functions. Specify this option to skip that orthonormalization.

Usually it need not be given, but when the `lmf --jobgw=1` calculation is slow, giving it can be expected to speed up the calculation.
* In ecalj the dielectric function can be calculated at finite q; then cancellation between numerator and denominator occurs, so the orthogonality of the wave functions must be accurate, and --skipGS should not be used in that case.

#### `--normcheck`
Used by `lmf --jobgw=1`. Use it to check the orthonormality of the wave functions used in the GW calculation.
`normchk.<sname>` looks like
```
> head -20 normchk.si
#       IPW         IPW(diag)   Onsite(tot)   Onsite(phi)      Total
      0.436015      0.805123      0.563972      0.562573      0.999988
      0.339134      0.620353      0.660515      0.656881      0.999649
      0.339133      0.620353      0.660516      0.656882      0.999649
      0.339133      0.620353      0.660516      0.656882      0.999649
      0.507738      0.648515      0.492040      0.487673      0.999778
...
```
The rightmost value should be 1, but because of the re-expansion it used to become considerably smaller (around 0.8) at high energies (further down; not seen here). Recently Gram-Schmidt orthonormalization is applied by default, so orthonormality is kept to about 8 digits or better.

The first line (corresponding to 1st band of 1st q point) means that total normalization almost unity = 0.999988 = 0.436015 + 0.563972. 

### `--ntqxx`
This fix the number of bands to calculate self energy at the first iteration for each $\bf q$ point in the IBZ.
In principle, the number is determined by the energy window: at each $\bf q$ point and spin, the bands up to `emax_sigm` above the Fermi energy (and at most `nband_sigm` bands) are used, and this count is made again at every iteration. With `--ntqxx`, the counts of the first iteration are written to the file `NTQXX` and read in the later iterations, so that the number of bands of $\Sigma$ (`sigm`, `QPU`) stays fixed through the iterations. Delete `NTQXX` to count again.

## Cautions
* gwsc checks QPU.[number]run and recognizes that number QSGW iterations have finished.
(To run from the initial state, delete all `*run*` directories and files.)

* A QSGW.[number]run directory stores the results of the number-th QSGW iteration, `rst`, `sigm` (plus `atmpnu`, `ctrlg.<sname>.toml`, `QPU`, `efermi.lmf`, and the logs of that iteration); band plots etc. can be made from them.

* Check progress and memory usage in the logs of each step (`lgw`, `lsxC`, `lvcc` etc.); the GW programs write the screen output of rank 0 there.
  Output per rank and how errors appear: `--fullstdo` in [Command-line options](./cmdopts).



---

## Mixing of the iterations: mixbeta and mixpriorit {#mixbeta-and-mixpriorit}

If `gwsc` does not converge stably, try `[gw].mixbeta = 0.5` in `ctrlg.<sname>.toml`.
(If `lmf` itself does not converge, use a smaller `[iter] b`, such as 0.015, 0.01, 0.005.)

What they do: at the end of every iteration `hqpe_sc` mixes the static self-energy `sigm`.
With $x_0$ the `sigm` the iteration ran with and $f$ the one GW has just made, the next `sigm` is

$$x_{\rm new} = x^* + \beta\,(f^* - x^*) \tag{M1}$$

where $x^*$, $f^*$ are the Anderson combination of the current and up to `mixpriorit` previous iterations (history in `__mixsig`).
$\beta$ = `mixbeta` (default 1.0), `mixpriorit` default 3. With no previous iteration, $x^*=x_0$ and $f^*=f$, i.e. plain linear mixing.

**The first iteration from LDA.** There is no previous `sigm` and no history, so $x_0=0$ (no self-energy, i.e. LDA) and Eq. (M1) gives $\beta f$:
the first step takes only $\beta$ of the GW self-energy. With the default $\beta=1$ it takes all of it; with $\beta=0.5$, half.
The second iteration is again linear mixing (one history entry); the Anderson extrapolation starts at the third.
The converged result does not depend on $\beta$; only the path does.
Table M1 shows `Samples/TestInstall/si_gwsc`, one `gwsc 1` from LDA (2026-09-26).

Table M1. `mixbeta` and the first iteration

| | `sigm` before → after mixing | band gap |
|---|---|---|
| LDA | — | 0.571 eV |
| $\beta=1.0$ (default) | 6.0488 → 6.0488 | 1.293 eV |
| $\beta=0.5$ | 6.0488 → 3.0244 | 0.934 eV (half the shift) |

The mixing of $\Sigma^{\rm MLO}$ in MLO-gwsc (`gwsc --mlo`) is in §1.6 of [MLO-gwsc](./mlo_gwsc).

## Other scripts 

`cleargw`: clean up temporary files

`gwscconv`: a wrapper that runs `gwsc` until automatic convergence. `--conv-tol <eV>` (e.g. `0.1`)
sets the stopping condition of the QSGW iterations. **Convergence criterion (revised 2026-06)**: stop when the gap fluctuation over the last 3 iterations
is below `conv-tol` twice in a row (prevents stopping on a single accidental convergence).
**Metals (2026-09-30)**: iterations where `lmf` gives no gap are judged by the change of eigenvalues.
Among the states with |e| < 5 eV in `eQP` (relative to E_F) of `QPU.<n>run` and `QPU.<n-1>run`, stop when the maximum change is
below `--conv-qp <eV>` (default 0.03) for 2 iterations in a row. To stop where there is no gap as before (exit 3), use
`--no-metal`. `--gpu --mp --fp32` given to it are passed on to `gwsc` as they are.
The GW1500 batch of ecalj_auto also calls this script.

`qsgw_status` (2026-06): a progress dashboard of the running `gwsc` (or gwscconv)
iterations. Use it in the run directory as follows.
```bash
qsgw_status               # show once
qsgw_status --last 21     # given the final iter number, also computes the ETA of completion
watch -n 30 qsgw_status   # resident monitor
```
It shows the elapsed time per iteration, the pace (median) and the
ETA from the mtime of `QPU.*run` / `sigm.iter*`; the current step (`hgw` etc.) and the elapsed time in the step from the end of the `gwsc` log;
and the sev difference of `save.<target>` and `rmsdel` of `lqpe` (σ convergence indicator).
The target name, the log and the iteration anchors are detected automatically. To avoid disturbing the run, it only
reads output files and does no heavy processing.

`gw_lmfh`: The one-shot GW calculation. Lifetime(impact ionization rate) of QPs.

`job_eps`: dielectric function. No local field corrections ([optical](./optical))
 
(`job_eps --lcf` : Dielectric function with local-field corrections. computationally expensive.)


<!-- \item
{\bf epsPP\_lmfh\_chipm} : non-interacting spin susceptibility. 
One-degree of freedom like Rigid moment approx.
After it ends, you need to do \verb#calj_nlfc_metal# and/or \verb#calj_summary_mat#
to get the full spin susceptibility. -->

`job_mloW` : the MLO model and the matrix elements of v, W and the cRPA W (`--crpa`) in it ([MLO](./mlo) section 6). (`genMLWFx`, the Wannier functions, was removed 2026-10-02.)


## Files used in gwsc
Temporary files are with `__*`. Thus we can delete __* (or use `cleargw`) after you finish `gwsc/job_eps` and so on.

A small test of `gwsc` is the install test `si_gwsc` ([Test framework](./dev_testecalj)).
In the `si_gwsc_work/` that `testecalj` makes, `ls -rlt` roughly shows which step generates which files.

The console output of `gwsc 1 si -np 2` for the same input is in the [main document](./README_tutorial#flow-of-qsgw-calculation-with-the-script-gwsc).
`lmf --jobgw=0` and `qg4gw` are run once, before the iterations; the other steps are repeated in every iteration.

### How to perform paper-quality QSGW calculations with minimum costs {#paper-quality-qsgw}
We expect that the accuracy of band gaps can be  ~0.1eV (or larger for larger band gap materials). 
In cases, it is easy but in cases not so easy (magnetic metals requires many number of iterations). 
So, it is better to use your own "simple criterion".
"Not stick to convergence so much. Just stick to Reproducibility."

As for the band gaps, 0.1 eV is a kind of limiting accuracy for usual semiconductors with 0~3 eV band gap materials.
For Si or something simple, you may be able to have better accuracy. So, just say ~0.1eV level of accuracy.

In principle, QSGW results should not depend on VWN or GGA. But there is minor dependence
(frozen core density, core eigenfunctions, radial basis functions, and a slight numerical reason; it affects the band gap of GaAs by about 0.02 eV).
Use VWN (`[ham] xcfun = 1`) as standard. Such technical things affect the band gap at the 0.05 eV level of error.

#### LDA calculation
We need to confirm LDA-level of calculations first.
`ctrlg.<sname>.toml` is generated just from ctrls.* (crystal structure file).
For calculation of QSGW, use large enough `[bz] nkabc`, so as to avoid
convergence check on them.

#### QSGW: how to check convergence
QSGW iteration cycle by gwsc contains (1) and (2).
- (1) One-body self-consistent calculation (where we add `sigm` = Σ − V_xc^LDA to the one-body potential) to determine the one-body Hamiltonian $H_0$.
- (2) For given $H_0$, we calculate the `sigm` file.

Big iteration cycle of QSGW is made from (1)+(2).
With (1), we have small iteration cycle of one-body calculation with keeping given sigm.

In `save.*`, we see a line per iteration of (1) with the total energy (which is not meaningful in the QSGW mode with the sigm file).
The meaning of the first letter of each line is in [lmf: `save.<sname>`](./lmf#save-sname); an `x` on the way is no problem as long as we finally see `c` in (1).

The command `grep '[cx] ' save.*` gives an indicator for
going to be converged or not. Or you can take `grep gap llmf.*run`.

Another way:
`dqpu QPU.3run QPU.6run` 
is to compare two QPU files which contains QP energies.
(note: QP energies shown are calculated just at the beginning of iteration).

For insulators, (I think), comparing band gap for each iteration 
is good enough to check convergence. But for metal, it is better to plot energy bands
for some of final iterations, and overlapped(`cd QSGW.*run` and run `job_band`).

Another way is `grep rms lqpe*`. This gives rmsdel. Difference of self-energy
(at least we see it is getting smaller for initial first cycles). 

Avoid unbalanced k points for an anisotropic cell (Comput. Phys. Commun. 176 (2007) 1-13); the offset-Γ points are set by
`[gw] Q0Pchoice` ([GW settings](./gwinput#q0pchoice)).

One-shot QSGW (not one-shot GW) can be useful in cases.
As it contains off-diagonal part, we can resolve band tanglement problem in Ge (no band gap).

#### How to make 80%QSGW +20% LDA, and SO setting
Note that the `sigm` file contains $V_{\rm xc}^{\rm QSGW}-V_{\rm xc}^{\rm LDA}$.
If `sigm` exists, `lmf` reads it, and runs self-consistent calculations
with adding `sigm` to the one-body potential.

Justifications of QSGW80 are:
* [JJAP 55, 051201](https://iopscience.iop.org/article/10.7567/JJAP.55.051201/pdf)
* [PRB 101, 205120](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.101.205120)
* [PRB 108, 165104](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.108.165104)

##### 1. QSGW80(NoSC)
For practical prediction of band structure, such as band gap and so
on, it may be better to use 80% QSGW +20% LDA procedure when you
make band plot. After you have the `rst` and `sigm` files determined self-consistently, run
```
job_band gaas -np 4 --ctrlg:ham.scaledsigma=0.80
```
This overrides `[ham] scaledsigma` in `ctrlg.gaas.toml` at run time, giving the QSGW80nosc result in Table II.

##### 2. QSGW80(NoSC)+SO
80%QSGW+20%LDA with `[ham].so = 1` (L·S method; needs `nspin = 2`).
If you like to include the L·S method:
```
mpirun lmf gaas -np 4 --ctrlg:ham.scaledsigma=0.80 --ctrlg:ham.so=1 --ctrlg:ham.nspin=2
```
This procedure makes self-consistency with keeping the sigm file. This may/(or may not) required. If you expect large orbital moment this procedure may be needed.
```
job_band gaas -np 4 --ctrlg:ham.scaledsigma=0.80 --ctrlg:ham.so=1 --ctrlg:ham.nspin=2
```
NOTE: nspin=2 is required for so=1. rst and sigm are expanded for nspin=2 (you can not run with nspin=1, after rst and sigm are expanded).
The general procedure for SOC bands is in [Other explanations](./UsageDetailed#spin-orbit-coupling).

##### 3. QSGW80
With `[ham] scaledsigma = 0.80` (or `gwsc ... --ctrlg:ham.scaledsigma=0.80`), you can run QSGW calculation in gwsc.
Then you have self-consistent results of QSGW80.
You can simultaneously use the setting so=2 (Lz·Sz scheme).
Be careful for z-direction and setting of `symgrp` (so as to keep the z-axis), for so=2.
If you like to get results of QSGW80+SO, you need to set so=1 after self-consistency of `sigm` is attained.

##### 4. Example of GaAs
Good example to check band gap, and SO splitting at top of valence of Gamma point for ZB structure as GaAs.
The keys are `scaledsigma`, `so` and `nspin` in `[ham]` of `ctrlg.gaas.toml`; look at them by
```
>grep -E "^(scaledsigma|so|nspin) " ctrlg.gaas.toml
```
They are overridden for a run by `--ctrlg:ham.scaledsigma=0.80`, `--ctrlg:ham.so=1`, `--ctrlg:ham.nspin=2` as in 1. and 2.
(the keys need not be written in the file for this).

#### Papers
It is instructive to reproduce samples in [Deguchi paper](../../presentations/deguchi2016.pdf). 
We can set up templates for your calculations. Ask us.
We have a latest paper at [arXiv:2506.03477](https://arxiv.org/abs/2506.03477) for GPU version.
It shows some details of computational steps in ecalj.
In Japanese, pages by Dr.Gomi at [gomisai.blog75.fc2.com](http://gomisai.blog75.fc2.com/blog-entry-675.html) and [Qiita](https://qiita.com/takaokotani/items/9bdf5f1551000771dc48).

### `lmfa`, `lmf`
`atmpnu*`, `rst.<sname>`, `__mixm.<sname>` and so on are in [lmf: output files](./lmf). Only those for GW are given here.

#### `__atm.<sname>`
This contains electron densities of spherical atoms.

#### QPLIST*.chk
QPLIST.jobgw1.chk (jobgw1 means --jobgw=1 for lmf) is for human, containing q for irr=1.
QPLIST.lmf.chk  (no jobgw option) is for human. Irreducible q points for lmf self-consistent calculations.

#### __HAMindex
q points table and so on for generating Hamiltonian

### `lmf --jobgw=0`
####  NLAindx.chk
This is for human. It shows index to expand eigenfunctions in MTs.

#### __HAMindex0
Generated by `m_hamindex0_init` called in main_lmf.f90.
Index of MTOs, space-group symmetries and so on.

#### QBZ.chk
q points mesh. Just for human reading.

### `qg4gw`

#### __QGpsi, __QGcou
q+G of the interstitial plane wave (IPW). Type lqg4gw, which shows
```
--- Max number of G for psi, G for Cou= 116 36
 iq=       1 q=  0.000000  0.000000  0.000000 ngp ngc=    111    29 irr.= 1 <--R
 iq=       2 q= -0.250000 -0.250000  0.750000 ngp ngc=    106    34 irr.= 1 <--R
 iq=       3 q= -0.250000  0.750000 -0.250000 ngp ngc=    106    34 irr.= 0 <--R
 iq=       4 q= -0.500000  0.500000  0.500000 ngp ngc=    116    36 irr.= 1 <--R
 iq=       5 q=  0.750000 -0.250000 -0.250000 ngp ngc=    106    34 irr.= 0 <--R
 iq=       6 q=  0.500000 -0.500000  0.500000 ngp ngc=    116    36 irr.= 0 <--R
 iq=       7 q=  0.500000  0.500000 -0.500000 ngp ngc=    116    36 irr.= 0 <--R
 iq=       8 q=  0.250000  0.250000  0.250000 ngp ngc=    116    36 irr.= 1 <--R
 iq=       9 q= -0.012500 -0.012500  0.037500 ngp ngc=    111    29 irr.= 1 <--Q0P
 iq=      10 q= -0.262500 -0.262500  0.787500 ngp ngc=    106    34 irr.= 1 <--Q0P+R
 iq=      11 q= -0.262500  0.737500 -0.212500 ngp ngc=    106    34 irr.= 1 <--Q0P+R
 iq=      12 q= -0.512500  0.487500  0.537500 ngp ngc=    116    36 irr.= 1 <--Q0P+R
 iq=      13 q=  0.737500 -0.262500 -0.212500 ngp ngc=    106    34 irr.= 0 <--Q0P+R
 iq=      14 q=  0.487500 -0.512500  0.537500 ngp ngc=    116    36 irr.= 0 <--Q0P+R
 iq=      15 q=  0.487500  0.487500 -0.462500 ngp ngc=    116    36 irr.= 1 <--Q0P+R
 iq=      16 q=  0.237500  0.237500  0.287500 ngp ngc=    116    36 irr.= 1 <--Q0P+R
 iq=      17 q= -0.012500  0.012500  0.012500 ngp ngc=    111    29 irr.= 1 <--Q0P
 iq=      18 q= -0.262500 -0.237500  0.762500 ngp ngc=    106    34 irr.= 1 <--Q0P+R
 iq=      19 q= -0.262500  0.762500 -0.237500 ngp ngc=    106    34 irr.= 0 <--Q0P+R
 iq=      20 q= -0.512500  0.512500  0.512500 ngp ngc=    116    36 irr.= 1 <--Q0P+R
 iq=      21 q=  0.737500 -0.237500 -0.237500 ngp ngc=    106    34 irr.= 1 <--Q0P+R
 iq=      22 q=  0.487500 -0.487500  0.512500 ngp ngc=    116    36 irr.= 1 <--Q0P+R
 iq=      23 q=  0.487500  0.512500 -0.487500 ngp ngc=    116    36 irr.= 0 <--Q0P+R
 iq=      24 q=  0.237500  0.262500  0.262500 ngp ngc=    116    36 irr.= 1 <--Q0P+R
  OK! End of qg4gw 
```
Here, we have regular mesh points specified by `<--R`. Q0P is the offset Gamma points shown in Q0P. 
irr=1 shows the irreducible q points at which we have to generate eigenfunctions.
ngp is the number of IPW for the expansion of eigenfunctions (controlled by `QpGcut_psi`).
ngc is for IPW for the MPB (controlled by `QpGcut_cou`).

#### __BZDATA 
__BZDATA contains info on regular mesh points, and offset Gamma (Q0P). Data for tetrahedron division. 

### `lmf --jobgw=1`
Generate all the following data to perform GW.
See at subroutines/sugw.f90.

#### @MNLA_core.chk 
human readable: core index

#### @MNLA_CPHI 
human readable (but program use this):
Eigenfunctions expanded within MT

#### hbe.d.chk 
human readable: size file for check

####  __PHIVC
radial functions.

#### __MTOindex
MTO index

#### __VxcEvec*
 Coefficients of eigenfunctions and  eigenvalues for $\langle F_i|H^0|F_j\rangle$ in the basis of PMT$\{F_i\}$.

#### __GEIG, __CPHI, __EValue
__GEIG: Coefficients of eigenfunctions. IPW part
__CPHI: Coefficients of eigenfunctions. MT  part
__EValue: eigenvalue

#### __PPOvlp, __PPOvlpG, __PPOvlpGG
overlap matrix of IPW.

#### __VXCFP
XC term in LDA (only diagonal part). A part of vxcevec
Used at hsfp0 but not essential (only for convenience of presentantion).

### `heftet --job=1`

#### EFERMI
The Fermi energy in the tetrahedron method

#### EFERMI_kbt
Written when `[gw] t_tetrakbt` > 0: the Fermi energy at that temperature, used for $\chi_0$ and $\Sigma$ ([kBT](./kBT) §3.3).


### hbasfp0 (we have --job=3 for core and --job=0 for valence)
#### __BASFP* 
 Product basis functions
#### __PPBRD*
 Radial integrals on each MT, symbolically written as  $\int \phi(r) \phi(r) B(r) dr$

### `hvccfp0` 
We call two times one for core, and the other for valence.

#### __Vcoud*
the Coulomb matrix (eigenvalues and eigenfunctions of the Coulomb matrix in the expansion of MPB)

### `hgw`
`hgw --jobgw=1` calculates the valence exchange, the screened Coulomb interaction $W$ and the valence correlation in one run.
$W-v$ in the expansion of mixed product basis (along the real axis and along the imag axis) is kept in memory and handed to the
calculation of the correlation $q$ by $q$; it is not written to files.
Note that we use a technique to define $W-v$ at ${\bf q}=0$ as an average of the Gamma cell.

#### __WV.d
Size of the dielectric function

#### __WVR.*, __WVI.*
Only with `gwsc --dumpW` (for analysis): $W-v$ along the real axis and along the imag axis, a file for each $q$.

#### freq_r
human readable: energy mesh to accumrate imaginary parts of W-v.


### `hsfp0_sc --job=3` (core exchange) and `hgw` (valence exchange and valence correlation)
These are moved to SEBK at the end of gwsc iteration cycle.

#### SEXcoreU, SEXcore2U
core exchange by `hsfp0_sc --job=3` . SEXcoreU is diagonal part for human check but not used so much. 
We have *D for down spin (isp=2)as well.
#### SEXU, SEX2U
valence exchange by `hgw` .       SEXU is diagonal part for human check but not used so much. 
#### SECU, SEC2U
valence correlation by `hgw` .     SECU is diagonal part for human check but not used so much. 
#### XCU
LDA exchange correlation

### `lqpe`

#### QPU, QPD
human readable format. 
decomposion of self-energy for diagonal elements.QPD is for isp=2.

An example of one-shot GW by `gw_lmfh si` (small size calculation in TestInstall) is:
```
...
           q               state  SEx   SExcore SEc    vxc    dSE  dSEnoZ  eLDA    eQP  eQPnoZ   eHF  Z    FWHM=2Z*Simg  ReS(elda)
  0.00000  0.00000  0.00000  1  -16.91  -1.85   6.62 -12.47   0.22   0.33 -12.24 -12.03 -11.92 -18.54 0.66   1.25708    -12.14031
  0.00000  0.00000  0.00000  2  -13.87  -1.96   2.81 -13.61   0.47   0.59  -0.31   0.16   0.28  -2.53 0.80   0.00000    -13.02308
  0.00000  0.00000  0.00000  3  -13.87  -1.96   2.81 -13.61   0.47   0.59  -0.31   0.16   0.28  -2.53 0.80   0.00000    -13.02308
  0.00000  0.00000  0.00000  4  -13.87  -1.96   2.81 -13.61   0.47   0.59  -0.31   0.16   0.28  -2.53 0.80   0.00000    -13.02308
  0.00000  0.00000  0.00000  5   -4.60  -1.41  -4.27 -11.81   1.19   1.52   2.23   3.42   3.75   8.03 0.78  -0.02515    -10.28546
  0.00000  0.00000  0.00000  6   -4.60  -1.41  -4.27 -11.81   1.19   1.52   2.23   3.42   3.75   8.03 0.78  -0.02515    -10.28546
  0.00000  0.00000  0.00000  7   -4.60  -1.41  -4.27 -11.81   1.19   1.52   2.23   3.42   3.75   8.03 0.78  -0.02515    -10.28546
  0.00000  0.00000  0.00000  8   -5.11  -3.79  -5.14 -15.23   0.91   1.20   2.95   3.86   4.14   9.28 0.76  -0.07397    -14.03341

  0.50000  0.00000  0.00000  1  -16.80  -1.91   5.97 -12.64  -0.06  -0.10 -11.18 -11.24 -11.27 -17.25 0.67   0.84187    -12.73584
  0.50000  0.00000  0.00000  2  -13.77  -2.37   2.91 -13.84   0.47   0.60  -3.84  -3.37  -3.24  -6.14 0.78   0.04107    -13.24015
  0.50000  0.00000  0.00000  3  -13.57  -1.74   3.01 -12.76   0.37   0.47  -2.20  -1.84  -1.74  -4.75 0.78   0.00000    -12.29249
  0.50000  0.00000  0.00000  4  -13.57  -1.74   3.01 -12.76   0.37   0.47  -2.20  -1.84  -1.74  -4.75 0.78   0.00000    -12.29249
  0.50000  0.00000  0.00000  5   -4.27  -1.19  -4.08 -10.97   1.17   1.43   0.76   1.92   2.19   6.27 0.81  -0.00000     -9.53378
  0.50000  0.00000  0.00000  6   -3.83  -1.10  -4.38 -10.98   1.33   1.68   2.98   4.31   4.65   9.03 0.79  -0.01203     -9.30307
  0.50000  0.00000  0.00000  7   -4.73  -1.66  -4.56 -12.65   1.34   1.70   5.45   6.79   7.14  11.70 0.79  -0.08743    -10.95431
  0.50000  0.00000  0.00000  8   -4.73  -1.66  -4.56 -12.65   1.34   1.70   5.45   6.79   7.14  11.70 0.79  -0.08743    -10.95431

  1.00000  0.00000  0.00000  1  -15.60  -2.13   4.43 -13.20  -0.08  -0.11  -8.13  -8.21  -8.24 -12.67 0.78   0.61936    -13.30329
  1.00000  0.00000  0.00000  2  -15.60  -2.13   4.43 -13.20  -0.08  -0.11  -8.13  -8.21  -8.24 -12.67 0.78   0.61936    -13.30329
  1.00000  0.00000  0.00000  3  -13.66  -1.70   3.17 -12.58   0.30   0.39  -3.16  -2.86  -2.77  -5.94 0.77   0.07961    -12.19216
  1.00000  0.00000  0.00000  4  -13.66  -1.70   3.17 -12.58   0.30   0.39  -3.16  -2.86  -2.77  -5.94 0.77   0.07961    -12.19216
  1.00000  0.00000  0.00000  5   -3.97  -0.91  -3.95 -10.33   1.22   1.50   0.31   1.53   1.81   5.76 0.81  -0.00000     -8.82976
  1.00000  0.00000  0.00000  6   -3.97  -0.91  -3.95 -10.33   1.22   1.50   0.31   1.53   1.81   5.76 0.81  -0.00000     -8.82976
  1.00000  0.00000  0.00000  7   -3.59  -2.37  -5.91 -13.53   1.20   1.66   9.81  11.01  11.47  17.37 0.72  -0.40499    -11.86796
  1.00000  0.00000  0.00000  8   -3.59  -2.37  -5.91 -13.53   1.20   1.66   9.81  11.01  11.47  17.37 0.72  -0.40499    -11.86796
```
All of the unit of energy is in eV. Detailed value of  `eLDA` is in `TOTE.UP`.
For insulators, the Fermi energy is at the middle of the band gap (the energies are measured from it). For metals, one shot GW can be problematic
if we consider the self-consistency of the Fermi energy.

- q: q vector
- state: band index n for valence
- SEx: $\langle\Psi_{{\bf k}n}|\Sigma_{\rm x}^{\rm valence}({\bf r},{\bf r}^{\prime})|\Psi_{{\bf k}n}\rangle$
- SExcore: $\langle\Psi_{{\bf k}n}|\Sigma_{\rm x}^{\rm core}({\bf r},{\bf r}^{\prime})|\Psi_{{\bf k}n}\rangle$
- SEc: $\langle\Psi_{\mathbf{k}n}|\Sigma_{\rm c}^{\rm valence}({\bf r},{\bf r}^{\prime}, \epsilon_n({\bf k}))|\Psi_{\mathbf{k}n}\rangle$
- vxc: LDA exchange correlation energy $\langle \Psi_{{\bf k}n}|V_{\rm xc}^{\rm LDA}([n_{\rm total}],{\bf r})|\Psi_{{\bf k}n}\rangle$
- dSE: $Z_{n{\bf k}}\times$ dSEnoZ
- dSEnoZ: $\langle\Psi_{{\bf k}n}|\Sigma_{\rm x}^{\rm core}+\Sigma_{\rm xc}^{\rm valence}(\epsilon_n({\bf k}))|\Psi_{{\bf k}n}\rangle - \langle\Psi_{{\bf k}n}|V_{\rm xc}^{\rm LDA}|\Psi_{{\bf k}n}\rangle$ = SEx + SExcore + SEc − vxc
- eLDA: LDA eigenvalue $\epsilon_n({\bf k})$
- eQP: QP energy $\epsilon_n({\bf k})$ + dSE
- eQPnoZ: QP energy without $Z$, $\epsilon_n({\bf k})$ + dSEnoZ
- eHF: HF energy of the first iteration, $\epsilon_n({\bf k})$ + SEx + SExcore − vxc
- Z: Z factor $Z_{n{\bf k}}$
- FWHM=2Z*S_{img}: quasiparticle lifetime $2 Z_{n{\bf k}} \times {\rm Im}\langle\Psi_{{\bf k}n}|\Sigma_{\rm c}^{\rm valence}(\epsilon_n({\bf k}))|\Psi_{{\bf k}n}\rangle$
- ReS(elda): ${\rm Re} \langle\Psi_{{\bf k}n}|\Sigma_{\rm x}^{\rm core}+\Sigma_{\rm xc}^{\rm valence}(\epsilon_n({\bf k}))|\Psi_{{\bf k}n}\rangle$

* NOTE: QPU for `gwsc` is a little different. No Z and no life time shown;
  the shown eQP is just the eigenvalues of the starting point of `lmf`. To compare QP energies, look at the column dSEnoZ.

#### TOTE.UP
Numerical detailed values of QPU. TOTE.DN for QPD.
In one-shot GW `gw_lmfh`, TOTE.UP contains LDA and QP energies (two kinds of QP energies, `QP` and `QPnoZ`).

#### QPU, TOTE in one-shot GW (`gw_lmfh`)
In `gw_lmfh` (one-shot GW), `hqpe` gives `TOTE.UP` and `QPU` files (`TOTE.DN` and `QPD` as well). They contain the same values.
One-shot GW is also used to calculate the impact ionization rate: in the sample `Samples/IIR/C` (diamond, a coarse mesh), `gw_lmfh` writes the widths to `QPU_life`.

#### __mixsig
Mixing file for `sigm.<sname>`.

#### `sigm.<sname>`
$\Sigma-V_{\rm xc}$ (the static Hermitian QSGW self-energy minus the LDA $V_{\rm xc}$, in Ry) as a matrix in the MTO basis at each irreducible $\bf q$ point, mixed with the previous iterations (`__mixsig`). `lmf` Bloch-sums it to real space and interpolates it to any $\bf k$ ([MLO-gwsc](./mlo_gwsc) §1.1).

### The Fermi energies in GW {#the-fermi-energies-in-gw}
We use two kinds of Fermi energy $E_{\rm FERMI}^{\rm smear}$ and $E_{\rm FERMI}^{\rm tetra}$.
This is because of numerical techniques.

One is $E_{\rm FERMI}^{\rm tetra}$, contained in EFERMI, for the polarization function since we use tetrahedron method for it.

The other is $E_{\rm FERMI}^{\rm smear}$ for $G \times W$ since we use smeared pole technique for $G$
(the width of the smearing is `[gw] t_sigmaw`).
$E_{\rm FERMI}^{\rm smear}$ is given at the beginning of $G \times W$ calculations. 
See the line `ef     esmr   valn=` in `lgw` and `lsxC`.

With `[gw] t_tetrakbt` > 0 (the polarization function at finite temperature), both the polarization function and
$G \times W$ (except the exchange of the cores) use the Fermi energy at that temperature, contained in EFERMI_kbt. See [kBT](./kBT).

### Product basis {#product-basis}
The product-basis settings (the cut-offs of `[product_basis]` and the per-atom tables) are in [GW settings: product basis](./gwinput#product-basis).

## MEMO
* We need to explain how to set Gamma-cell averaged $\tilde{W}({\bf q}=0,\omega)$.
* The old theory notes for developers, ecaljdetails (2015 to 2022, offset-Γ method, Coulomb matrix, EIBZ etc.), were removed on 2026-10-02. The original is `ecaljdetails/ecaljdetails.tex` in ecaljdoc commit `3284d0b`.
