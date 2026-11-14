# gwsc: a script to run QSGW calculation

> English version of the page [manual/gwsc](../../manual/gwsc.md). The pages at the top of ecaljdoc are the reference; this one was made from manual/gwsc.md of commit 034c2b739 (source blob ecf993d644ea).

> ⚠️ `gwsc` reads `ctrlg.<sname>.toml` only, and stops when the file is missing. A directory with the legacy `ctrl.<sname>` + `GWinput` is converted by `Legacy2toml.py <sname>` (run once per directory). See [TOML migration](./toml_migration) and the [Samples/TestInstall/*_gwsc/](https://github.com/tkotani/ecalj/tree/main/Samples/TestInstall) directories as templates.

gwsc is the script to run QSGW calculations.
A QSGW calculation is run by calling several Fortran executables.
The input is read from `ctrlg.<sname>.toml` (GW settings are in `[gw]`; the keys are explained in [gwinput](./gwinput.md)).

## Usage
**Usage**: <del class="clfix">` gwsc -np NP [-np2 NP2] [--gpu] [--prec=tf32|fp32|fp64] [--mp] [--fp32] nloop extension [Options]`</del><span class="clfix">` gwsc -np NP [-np2 NP2] [--gpu] [--prec=tf32|fp32|fp64] [--prec-final=<prec>[:N]] [--no-tetwt-helper] [--mp] [--fp32] [--mlo] [--dumpW] [--fullstdo] nloop extension [Options]`</span>

### `-np NP`
Number of MPI processes.

### `-np2 NP2`
Specify only when using the GPU version. Number of MPI processes for the calculations run on GPUs.
Usually set to the number of available GPUs.

### `--gpu`
Specify only when using the GPU version.

### `--prec=tf32|fp32|fp64` (2026-09)
Precision of the GW part. `fp64` is the double-precision version; `fp32` is the mixed-precision version with single-precision products in FP32 (same as `--mp --fp32`); `tf32` is `fp32` except that
only the last product of Σc uses a 10-bit mantissa for the inputs (TF32 or FP16, chosen by a per-GPU table) (for LiTi₂O₄ 6³, Re Σc within E_F ±1 eV agrees with double precision within 1.3 meV, and `hgw` on RTX 5090 takes about half the time of fp32: 173 vs 367 s, version of the night of 2026-09-27). The legacy `--mp` alone is the same as `tf32`
(until 2026-09-27 all products were TF32; the error was 6 times larger and the time shrank by only 10 %).
Which method computes products and inverse matrices on the GPU is chosen automatically by a table measured at installation (`<bindir>/ecalj_linalg_policy.toml`).
Details in [GPU manual § choice of precision and automatic choice of matrix operations](./ecaljgpu#choosing-the-precision-and-automatic-selection-of-matrix-operations-2026-09).
When given together with `--mp` and `--fp32`, they must not contradict (a contradiction stops the run).

### `--prec-final=<prec>[:N]` (2026-09)
Use a different precision only for the last N iterations of the run (default 1). Example: `gwsc 10 --prec=tf32 --prec-final=fp32:2 ...` runs 8 iterations with tf32 and the last 2 iterations with fp32.

### `--no-tetwt-helper` (2026-09)
With `--gpu`, do not use the helper run that computes the tetrahedron weights on the CPU beforehand (`hgw --tetwt_write`).
It is used by default ([GPU manual § tetrahedron weights](./ecaljgpu#computing-the-tetrahedron-weights-first-on-the-cpu)).

### `--mp`
Use the mixed-precision version. Since 2026-09-27 it means the same as `--prec=tf32` (only the last product of Σc with a 10-bit input mantissa, other single-precision products in FP32).
Before that, all single-precision products on the GPU were TF32, and for ill-conditioned dielectric matrices (heavy elements + molecular anions NO3 / N3 / ClO etc.) the TF32 error was
amplified in $W$ / $\Sigma^c$. For the choice of precision see `--prec` above and [GPU manual § choice of precision](./ecaljgpu#choosing-the-precision-and-automatic-selection-of-matrix-operations-2026-09).

### `--fp32` (2026-06)
Used with `--mp`; same meaning as `--prec=fp32` (the products of Σc also in FP32). The GW1500 batch of ecalj_auto (`gwscconv --gpu --mp --fp32 --conv-tol 0.1`)
uses this by default. Reproductions of cases broken by TF32 are in [`Samples/mptf32problem/`](https://github.com/tkotani/ecalj/tree/main/Samples/mptf32problem).


<div class="clfix">

### `--mlo`
QSGW in which the self-energy is interpolated in the MLO representation (MLO-QSGW, under development). `[bz] nkabc` and `[mlo] mlo_nkabc` must equal `[gw] n1n2n3` (it stops otherwise). See [MLO-gwsc](./mlo_gwsc).

</div>

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
- `zmel_batch_gb`: size of the batch in which matrix elements are made at once ([gwinput](./gwinput.md))
- `mpi_worker_exch`, `mpi_worker_corr`: number of ranks assigned to one group of $q$ points in the exchange and correlation calculations. 0 (default) is automatic. Give a divisor of the total number of ranks

There are no options `--keepwv` or `--nwpara=X` (giving them stops the program). `--nb=X` is accepted but not used in the `gwsc` calculation.

#### `--tetwtk`
If given, the tetrahedron weights between the coupled states are not kept in memory when calculating the polarization function. Use it when memory runs short in calculations with many $k$ points.

#### `--skipGS`
Used by `lmf --jobgw=1`.
In the GW calculation, the wave functions obtained by the DFT calculation (`lmf`) are re-expanded in basis functions different from those of `lmf`. Gram-Schmidt orthonormalization is applied to the re-expanded wave functions. Specify this option to skip that orthonormalization.

Usually it need not be given, but when the `lmf --jobgw=1` calculation is slow, giving it can be expected to speed up the calculation.
* In ecalj the dielectric function can be calculated at finite q; then cancellation between numerator and denominator occurs, so the orthogonality of the wave functions must be accurate, and --skipGS should not be used in that case.

#### `--normcheck`
Used by `lmf --jobgw=1`. Use it to check the orthonormality of the wave functions used in the GW calculation.
normchk.fobar looks like
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
In principle, the number is determined by<span class="clfix"> the energy window: at each $\bf q$ point and spin, the bands up to `emax_sigm` above the Fermi energy (and at most `nband_sigm` bands) are used, and this count is made again at every iteration. With `--ntqxx`, the counts of the first iteration are written to the file `NTQXX` and read in the later iterations, so that the number of bands of $\Sigma$ (`sigm`, `QPU`) stays fixed through the iterations. Delete `NTQXX` to count again.</span>

## Cautions
* gwsc checks QPU.[number]run and recognizes that number QSGW iterations have finished.
(To run from the initial state, delete all `*run*` directories and files.)

* A QSGW.[number]run directory stores the results of the number-th QSGW iteration, `rst`, `sigm` (plus `atmpnu`, `ctrlg.<sname>.toml`, `QPU`, `efermi.lmf`, and the logs of that iteration); band plots etc. can be made from them.

* Check progress and memory usage in the logs of each step (`lgw`, `lsxC`, `lvcc` etc.). The GW programs write the screen output of rank 0 there,
  and discard the output of other ranks (since 2026-09-27). When per-rank output is needed, use `--fullstdo` (it goes to `STDOUT/stdout.<rank>.<prog>`).
  Errors from any rank appear in the standard error of `gwsc` with the rank number.



---

# Other scripts 

`cleargw`: clean up temporary files

`gwscconv`: a wrapper that runs `gwsc` until automatic convergence. `--conv-tol <eV>` (e.g. `0.1`)
sets the stopping condition of the QSGW iterations. **Convergence criterion (revised 2026-06)**: stop when the gap fluctuation over the last 3 iterations
is below `conv-tol` twice in a row (prevents stopping on a single accidental convergence).
**Metals (2026-09-30)**: iterations where `lmf` gives no gap are judged by the change of eigenvalues.
Among the states with |e| < 5 eV in `eQP` (relative to E_F) of `QPU.<n>run` and `QPU.<n-1>run`, stop when the maximum change is
below `--conv-qp <eV>` (default 0.03) for 2 iterations in a row. To stop where there is no gap as before (exit 3), use
`--no-metal`. `--gpu --mp --fp32` given to it are passed on to `gwsc` as they are.
The GW1500 batch of ecalj_auto also calls this script (<del class="clfix">`MD/auto.md`</del><span class="clfix">development record (MD/auto.md in the repository, not on this site)</span>).

`qsgw_status` (2026-06): a progress dashboard of the running `gwsc` (or gwscconv)
iterations. Use it in the run directory as follows:
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

`gw_lmfh`: The one-shot <del class="clfix">\GW</del><span class="clfix">GW</span> calculation. Lifetime(impact ionization rate) of QPs.

`job_eps`: dielectric function. No local field corrections ([optical](./optical))
 
(`job_eps --lcf` : Dielectric function with local-field corrections. computationally expensive.)


<!-- \item
{\bf epsPP\_lmfh\_chipm} : non-interacting spin susceptibility. 
One-degree of freedom like Rigid moment approx.
After it ends, you need to do \verb#calj_nlfc_metal# and/or \verb#calj_summary_mat#
to get the full spin susceptibility. -->

`job_mloW` : the MLO model and the matrix elements of v, W and the cRPA W (`--crpa`) in it ([MLO](./mlo) section 6). (`genMLWFx`, the Wannier functions, was removed 2026-10-02.)


# Files used in gwsc
Temporary files are with `__*`. Thus we can delete __* (or use `cleargw`) after you finish `gwsc/job_eps` and so on.

To repeat a small test for gwsc:

```bash
cd ~/ecalj/Samples/TestInstall
testecalj si_gwsc -np 8
```

This is one of the install tests. `testecalj` creates
`si_gwsc_work/` next to `si_gwsc/`, copies inputs in, and runs there.
After that, `ls -rlt si_gwsc_work/` roughly shows which step generates
which files.

The console output of `gwsc 1 si -np 2` for the same input is
```text
===== Ititial band structure ====== 
--> No sigm. LDA caculation for eigenfunctions 
00:00:00.007   mpirun -np 1 /home/takao/bin/lmfa si stdout='llmfa'  Elap. 1.5s
00:00:01.481   mpirun -np 2 /home/takao/bin/lmf si --ctrlg:iter.b=0.5 stdout='llmf_lda'  Elap. 6.1s
lmf successful with b=0.5
Removing __mixm.si
00:00:07.563   mpirun -np 2 /home/takao/bin/lmf si --jobgw=0 stdout='llmfgw00'  Elap. 1.8s
00:00:09.335   mpirun -np 1 /home/takao/bin/qg4gw si --job=1 stdout='lqg4gw'  Elap. 1.0s
===== QSGW iteration start iter 1 === (GW precision fp64)
00:00:10.381   mpirun -np 2 /home/takao/bin/lmf si --jobgw=1 stdout='llmfgw01'  Elap. 3.7s
00:00:14.042   mpirun -np 1 /home/takao/bin/heftet si --job=1 stdout='leftet'  Elap. 1.1s
00:00:15.193   mpirun -np 1 /home/takao/bin/hbasfp0 si --job=3 stdout='lbasC'  Elap. 1.4s
00:00:16.570   mpirun -np 2 /home/takao/bin/hvccfp0 si --job=3 stdout='lvccC'  Elap. 1.4s
00:00:17.966   mpirun -np 2 /home/takao/bin/hsfp0_sc si --job=3 stdout='lsxC'  Elap. 1.9s
00:00:19.862   mpirun -np 1 /home/takao/bin/hbasfp0 si --job=0 stdout='lbas'  Elap. 2.1s
00:00:21.973   mpirun -np 2 /home/takao/bin/hvccfp0 si --job=0 stdout='lvcc'  Elap. 2.5s
00:00:24.472   mpirun -np 2 /home/takao/bin/hgw si --jobgw=1 stdout='lgw'  Elap. 4.8s
00:00:29.276   mpirun -np 1 /home/takao/bin/hqpe_sc si stdout='lqpe'  Elap. 1.5s
Using cached b-value 0.5 for /home/takao/bin/lmf
00:00:30.814   mpirun -np 2 /home/takao/bin/lmf si --ctrlg:iter.b=0.5 stdout='llmf'  Elap. 8.0s
lmf successful with b=0.5
===== QSGW iteration end   iter 1 ===
```
`lmf --jobgw=0` and `qg4gw` are run once, before the iterations; the other steps are repeated in every iteration.

## `lmfa`
### atmpnu*
Logalithmic derivative at MT boundaries generated by `lmfa`. This is
used by `lmf` when `[ham].readp = true` in `ctrlg.<sname>.toml`.

### `__atm.<sname>`
This contains electron densities of spherical atoms.

## `lmf`
### `rst.<sname>`
This contains self-consistent electron density revised by each iteration of lmf

### QPLIST*.chk
QPLIST.jobgw1.chk (jobgw1 means --jobgw=1 for lmf) is for human, containing q for irr=1.
QPLIST.lmf.chk  (no jobgw option) is for human. Irreducible q points for lmf self-consistent calculations.

### __HAMindex
q points table and so on for generating Hamiltonian

### `__mixm.<sname>`
mixing file for lda

## `lmf --jobgw=0`
###  NLAindx.chk
This is for human. It shows index to expand eigenfunctions in MTs.

### __HAMindex0
Generated by `m_hamindex0_init` called in main_lmf.f90.
Index of MTOs, space-group symmetries and so on.

### QBZ.chk
q points mesh. Just for human reading.

## `qg4gw`

### __QGpsi, __QGcou
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
ngp is the number of IPW for the expansion of eigenfunctions (controlled by QpGcut_phi). 
ngc is for IPW for the MPB (controlled by QpGcut_cou).

### __BZDATA 
__BZDATA contains info on regular mesh points, and offset Gamma (Q0P). Data for tetrahedron division. 

## `lmf --jobgw=1`
Generate all the following data to perform GW.
See at subroutines/sugw.f90.

### @MNLA_core.chk 
human readable: core index

### @MNLA_CPHI 
human readable (but program use this):
Eigenfunctions expanded within MT

### hbe.d.chk 
human readable: size file for check

###  __PHIVC
radial functions.

### __MTOindex
MTO index

### __VxcEvec*
 Coefficients of eigenfunctions and  eigenvalues for <del class="clfix">$\langle F_i|H^0|\psi F_j$</del><span class="clfix">$\langle F_i|H^0|F_j\rangle$</span> in the basis of PMT$\{F_i\}$.

### __GEIG,__CPHI,__EValue
__GEIG: Coefficients of eigenfunctions. IPW part
__CPHI: Coefficients of eigenfunctions. MT  part
__EValue: eigenvalue

### __PPOvlp, __PPOvlpG, __PPOvlpGG
overlap matrix of IPW.

### __VXCFP
XC term in LDA (only diagonal part). A part of vxcevec
Used at hsfp0 but not essential (only for convenience of presentantion).

## `heftet --job=1`

### EFERMI
The Fermi energy in the tetrahedron method

### EFERMI_kbt
Written when `[gw] t_tetrakbt` > 0: the Fermi energy at that temperature, used for $\chi_0$ and $\Sigma$ ([kBT](./kBT) §2.3).


## hbasfp0 (we have --job=3 for core and --job=0 for valence)
### __BASFP* 
 Product basis functions
### __PPBRD*
 Radial integrals on each MT, symbolically written as  $\int \phi(r) \phi(r) B(r) dr$

## `hvccfp0` 
We call two times one for core, and the other for valence.

### __Vcoud*
the Coulomb matrix (eigenvalues and eigenfunctions of the Coulomb matrix in the expansion of MPB)

## `hgw`
`hgw --jobgw=1` calculates the valence exchange, the screened Coulomb interaction $W$ and the valence correlation in one run.
$W-v$ in the expansion of mixed product basis (along the real axis and along the imag axis) is kept in memory and handed to the
calculation of the correlation $q$ by $q$; it is not written to files.
Note that we use a technique to define $W-v$ at ${\bf q}=0$ as an average of the Gamma cell.

### __WV.d
Size of the dielectric function

### __WVR.*, __WVI.*
Only with `gwsc --dumpW` (for analysis): $W-v$ along the real axis and along the imag axis, a file for each $q$.

### freq_r
human readable: energy mesh to accumrate imaginary parts of W-v.


## `hsfp0_sc --job=3` (core exchange) and `hgw` (valence exchange and valence correlation)
These are moved to SEBK at the end of gwsc iteration cycle.

### SEXcoreU,SEXcore2U :  
core exchange by `hsfp0_sc --job=3` . SEXcoreU is diagonal part for human check but not used so much. 
We have *D for down spin (isp=2)as well.
### SEXU,SEX2U 
valence exchange by `hgw` .       SEXU is diagonal part for human check but not used so much. 
### SECU,SEC2U:  
valence correlation by `hgw` .     SECU is diagonal part for human check but not used so much. 
### XCU
LDA exchange correlation
## `lqpe`

### QPU, QPD
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
All of the unit of energy is in eV. Detailed value of  <del class="clfix">eLDA} is in {TOTE.UP}</del><span class="clfix">`eLDA` is in `TOTE.UP`</span>.
For insulators, the Fermi energy is at the middle of <del class="clfix">band</del><span class="clfix">the band gap (the energies are measured from it)</span>. For metals, one shot GW can be problematic
if we consider the self-consistency of the Fermi energy.

q  : q vector

state:  Band index n for valence.

SEx: $= \langle\Psi_{{\bf k}n}|\Sigma_{\rm x}^{\rm valence}({\bf r},{\bf r}^{\prime})|\Psi_{{\bf k}n}\rangle$

SExcore: $= \langle\Psi_{{\bf k}n}|\Sigma_{\rm x}^{\rm core}({\bf r},{\bf r}^{\prime})|\Psi_{{\bf k}n}\rangle$

SEc: <del class="clfix">$= \langle\Psi_{{\bf k}n}|\Sigma_{\rm c}^{\rm valence}({\bf r},{\bf r}^{\prime}, \epsilon_n({\bf k})|\Psi_{{\bf k}n}\rangle$</del><span class="clfix">$= \langle\Psi_{\mathbf{k}n}|\Sigma_{\rm c}^{\rm valence}({\bf r},{\bf r}^{\prime}, \epsilon_n({\bf k}))|\Psi_{\mathbf{k}n}\rangle$</span>

vxc: LDA exchange correlation energy.$\langle \Psi_{{\bf k}n}|V_{\rm xc}^{\rm LDA}([n_{\rm total}],{\bf r})|\Psi_{{\bf k}n}\rangle$

dSE: $=Z_{n{\bf k}}\times$ dSEnoZ

dSEnoZ: $=\langle\Psi_{{\bf k}n}|\Sigma_{\rm x}^{\rm core}({\bf r},{\bf r}^{\prime})+\Sigma_{\rm xc}^{\rm valence}({\bf r},{\bf r}^{\prime},\epsilon_n({\bf k}))|\Psi_{{\bf k}n}\rangle - \langle\Psi_{{\bf k}n}|V_{\rm xc}^{\rm LDA}([n_{\rm total}],{\bf r})|\Psi_{{\bf k}n}\rangle$= SEx + SExcore + SEc - vxc

eLDA: LDA eigenvalues. $\epsilon_n({\bf k})$

eQP: QP energy.  $\epsilon_n({\bf k})$+dSE

eQPnoZ: QP energy without $Z$. $\epsilon_n({\bf k})$+dSEnoZ

eHF: HF energy of 1st iteration. $\epsilon_n({\bf k})$+SEx + SExcore -vxc
Z: Z factor. $Z_{n{\bf k}}$

FWHM=2Z*S_{img}: Quasi-particle life time. $2 Z_{n{\bf k}} \times {\rm Im}\langle\Psi_{{\bf k}n}|\Sigma_{\rm c}^{\rm valence}({\bf r},{\bf r}^{\prime},\epsilon_n({\bf k}))|\Psi_{{\bf k}n}\rangle$  

ReS(elda):${\rm Re} \langle\Psi_{{\bf k}n}|\Sigma_{\rm x}^{\rm core}({\bf r},{\bf r}^{\prime})+\Sigma_{\rm xc}^{\rm valence}({\bf r},{\bf r}^{\prime},\epsilon_n({\bf k}))|\Psi_{{\bf k}n}\rangle$ 

* NOTE: QPU for `gwsc` is a little different. No Z and no life time shown.
  Shown eQP is just the eigenvalues of starting point of lmf. 

### TOTE.UP
numerical detailed values of QPU. TOTE.DN for QPD

In one-shot GW `gw_lmfh`, TOTE.UP contains LDA and QP energies. 
It contains two kind of QP energies <del class="clfix">{\tt QP QPnoZ}</del><span class="clfix">`QP` and `QPnoZ`</span>.

### __mixsig
mixing file for `sigm.<sname>`

### `sigm.<sname>`
<del class="clfix">self-energy file in the expansion of eigenfunctions of $H^0$.</del><span class="clfix">$\Sigma-V_{\rm xc}$ (the static Hermitian QSGW self-energy minus the LDA $V_{\rm xc}$, in Ry) as a matrix in the MTO basis at each irreducible $\bf q$ point, mixed with the previous iterations (`__mixsig`). `lmf` Bloch-sums it to real space and interpolates it to any $\bf k$ ([MLO-gwsc](./mlo_gwsc) §1.1).</span>


# Product basis

The product basis (the union of MT-local products and interstitial
plane waves) is the heart of how ecalj represents the screened Coulomb
interaction `W`. It is specified in the **`[product_basis]` section of
`ctrlg.<sname>.toml`**, always the last section of the file:

| key | what |
|---|---|
| `pb_tolerance`, `pb_lcutmx` | global cut-offs — the part you may tune |
| `nlx`, `valence`, `core` | per-atom radial-function tables, written by `gwinit` (via `ctrlgenToml.py` / `Legacy2toml.py`) from `[[spec]]`; almost never hand-edited |

(Until 2026-09 the three tables were a separate `PB.<sname>.toml`. The
binaries do not read such a file: they abort and point at
`ctrlg_absorb.py <sname>`, which appends its tables to `[product_basis]`
and keeps the original as `.bk`.)

> Product basis is originally given by F.Aryasetiawan and O.Gunnarsson
> ([PRB 49, 16214](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.49.16214)).
> The mixed product basis
> ([Solid State Commun. 2002](https://doi.org/10.1016/S0038-1098(02)00028-5))
> is the successor of original.

## `[product_basis]` cut-offs (in `ctrlg.<sname>.toml`)

```toml
[product_basis]
pb_tolerance = [0.001]   # cut off of linear dependency (default 1e-3)
                         #   larger ⇒ fewer product basis. 1d-2 is OK
                         #   when you want a cheaper calculation.
pb_lcutmx    = [4, 4]    # per-site l-cutoff for the product basis.
                         #   4 needed when valence 3d exists (Ni, Ga, ...);
                         #   2 OK for oxygen; 6 needed for 4f / 5f atoms.
                         #   site order = [[site]] order; see lmchk SiteInfo.
```

After editing, re-run `lmf --jobgw=1` (or the next `gwsc` step) — the
`hbasfp0` driver consumes these values and reports the resulting basis
sizes in `lbas` / `lbasC`.

## The per-atom tables (`nlx`, `valence`, `core`)

They follow the cut-offs in the same `[product_basis]` table; every row
starts with the site index `iatom` (in `[[site]]` order). `gwinit` writes
them from `rsmh`/`eh`/`pz` of `[[spec]]` and adds a label per row; the
explanation below is written into the file itself as comments, and the
three tables are divided by `# ----` rules.

```toml
[product_basis]
pb_tolerance = [1e-3]
pb_lcutmx    = [4]
# ----------------------------------------------------------------
# nlx [iatom, l, nnvv, nnc]: one row per (iatom, l) ...
nlx = [
  [1, 0, 2, 3],
  [1, 1, 2, 2],
  [1, 2, 2, 0],
  [1, 3, 2, 0],
  [1, 4, 2, 0],
]

# ----------------------------------------------------------------
# valence [iatom, l, n, occ, unocc]: one row per valence radial function ...
valence = [
  [1, 0, 1, 1, 1],   # 4s_phi
  [1, 0, 2, 0, 0],   # 4s_phidot
  [1, 1, 1, 1, 1],   # 4p_phi
  [1, 1, 2, 0, 0],   # 4p_phidot
  [1, 2, 1, 1, 1],   # 3d_phi
  [1, 2, 2, 0, 0],   # 3d_phidot
  [1, 3, 1, 0, 1],   # 4f_phi
  [1, 3, 2, 0, 0],   # 4f_phidot
  [1, 4, 1, 0, 0],   # 5g_phi
  [1, 4, 2, 0, 0],   # 5g_phidot
]

# ----------------------------------------------------------------
# core [iatom, l, n, occ, unocc, forX0, forSxc]: one row per core level ...
core = [
  [1, 0, 1, 0, 0, 0, 0],   # 1S
  [1, 0, 2, 0, 0, 0, 0],   # 2S
  [1, 0, 3, 0, 0, 0, 0],   # 3S
  [1, 1, 1, 0, 0, 0, 0],   # 2P
  [1, 1, 2, 0, 0, 0, 0],   # 3P
]
```

Semantics (same as the legacy `<PRODUCT_BASIS>` block — only the
container changed):

- `nlx` row `[iatom, l, nnvv, nnc]` — how many radial functions ecalj has
  on each `l` channel of that site: `nnvv = 2` is `{phi, phidot}`,
  `nnvv = 3` adds the local orbital `phiz`; `nnc` is the number of core
  levels at that `l`. Derived, do not edit.
- `valence` row `[iatom, l, n, occ, unocc]` — one per radial function,
  `n = 1` phi, `2` phidot, `3` phiz. `occ` / `unocc` are 0/1 switches that
  put the function into the "occupied" / "unoccupied" group; the product
  basis is formed from products (occupied group) × (unoccupied group) and
  pruned by `pb_tolerance`. `phidot` is usually left out (`0 0`) for speed;
  including it gives a more accurate basis.
- `core` row `[iatom, l, n, occ, unocc, forX0, forSxc]` — one per core
  level (`n = 1` deepest). `forX0` / `forSxc` came from the CORE1/CORE2
  distinction of PRB 76, 165106 (2007). Historical: modern runs leave
  every core row all-zero.

The site index `iatom` follows the order of `[[site]]` in `ctrlg.<sname>.toml`;
`lmchk` prints the order at the top of its console output.

> **Legacy input.** `Legacy2toml.py` translates the pre-TOML
> `<PRODUCT_BASIS>` block of `GWinput` into the `[product_basis]` section
> above.


# MEMO
* We need to explain how to set Gamma-cell averaged $\tilde{W}({\bf q}=0,\omega)$.
* The old theory notes for developers, ecaljdetails (2015 to 2022, offset-Γ method, Coulomb matrix, EIBZ etc.), were removed on 2026-10-02. The points still valid are in <del class="clfix">ecalj's `MD/past_log.md` §14</del><span class="clfix">the development record (MD/past_log.md §14 in the repository, not on this site)</span>; the original is `ecaljdetails/ecaljdetails.tex` in ecaljdoc commit `3284d0b`.
