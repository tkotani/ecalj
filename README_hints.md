## SOME NOTES

### -- Runtime option to lmf lmfa lmchk 
ecalj made from band structure part (lmf), and GW part (plus some additional functionalities such as Wannier). 

lmfa (spherical atom for initial condition) and lmchk (crystal structure check) are by single core.

> usage: mpirun -np 8 lmf [--OPTION] [-vxxx] [ext]

[ext] is for ctrl.[ext]. -vxxx is often used to give value. 

 --h, --help         Print this message, and quit
  --pr#1         Set the verbosity 
 --time=#1[,#2] Print timing info to # levels (#1=summary; #2=on-the-fly). Try --time=5,5

 -vnam=expr     Define numerical variable "nam"; set to result of 'expr'
                Eg. -vnk1=3 replace {nk1} in ctrl with 3.  

--quit=dmatu 
    Quit after initial setup. Convenient for check.
--quit=band 
    Quit after band calculation.

Search --quit option in ```SRC/*/*.f90```


See some details of ctrl file at at
```ecalj/Document/help_lmf.org, help_lmfjobgw.org, and help_lmfa.org```,
but not maintained well. 

###  -- lmchk  
Note that ordering of atoms in the cell are shown by 
 >lmchk si

or so. The atoms ordering specified by ctrl do not mean atom id in the calculation.


###  -- how to perform paper-quarilty QSGW calculations with minimum costs. 

The accuracy of band gaps can be  ~0.1eV or larger for larger band gap materials.. 
In cases, it is easy, but in cases not so easy. So, it is better to use your own "simple criterion".
"Not stick to convergence so much. Just stick to Reproducibility."

Caution: For 4f and probably also for 5f systems, some special
care is required; just defaults ctrlgenM1 do not work.; 
See [HowToSet4f_GdQSGW4.pdf](Document/HowToSet4f_GdQSGW4.pdf)

Except 4f systems, use default setting (just change k points).
However, we may need to reduce computaiional time. Following are hints.

#### LDA calculation
 We need to confirm LDA-level of calculations first.
 The ctrl file is generated just from ctrls.* (crystal structure file)
 by ctrlgen.*. However, we pay attention to MMOM (initial magnetic moment), and k
 points NKABC or nk1,nk2,nk3.

 For calculation of GW, use large enough NKABC, so as to avoid
 convergence check on them. 

#### Setting in GWinput.
#####  Set # of k points
  We use different number of k points for self-energy.
  The 6x6x6 k points is good setting for ZB structure (2 atoms per cell).
  It is better to use this level of k points.
  In other words, 6x6x6x2 \sim 432 \sim (k points \times atom number)
  should be used for calculaitons.
  For example, when we try 8 atoms per cell, we can use 4x4x4 or 3x3x3
  is fine because 4x4x4x8 \sim 400.
  For metallic systems, larger is fine, but limited by computational
  time. See takao kotani's papers, for example, https://doi.org/10.1103/PhysRevB.93.075125
  In my observation, good news is that we don't need to use so many k
  points as in one-body part (NKABC in ctrl).
  4x4x4 for ZB is not so bad --- roughly speaking, a lower limit for
  publication probably. This means 3x3x3x8 (for 8 atom case) is not so bad. (3x3x3x8 > 4x4x4x2).
  See an examination 
  https://doi.org/10.7566/JPSJ.83.094711
  http://doi.org/10.7567/JJAP.55.051201
##### Reduce # of lcutmx
  In my experience, this is effective to reduce computaional time.
  To reduce the computaitonal time, we reduce number of MPB  (mixed product basis).
  One is lcut off of PB within MT. Use 2 for oxygen or something (s,p block atoms). Thus it
  is like
  > lcutmx(atom) = maximum l-cutoff for the product basis.  =4 is required for atoms with valence d, like N
   4 4 4 2 2 2 

   in a section of GWinput
  (NOTE: we know that lcutmx =6 is requied for 4f systems. See above)

##### Reduce IPW of MPB, Reduce IPW of psi and Reduce emax_sigm, and pwemax.
   QpGcut_cou is the Interstitial plane wave (IPW) for MPB. 
   QpGcut_psi is for expantion of eigenfunctions.
   emax_sigm is the upper cutoff (relative to the Fermi energy) to calculate self energy.
   pwemax (in ctrl) is the APW basis cutoff for the eigenfunciton.
   To reduce computational time, we may use

   >QpGcut_psi 3.0
   QpGcut_cou 2.5
   emax_sigm 2.0
   pwemax=2 (in ctrl file).

   We sometimes use this setting as long as the numerical
   results are affected little (check this with small number of k
   points).

##### high resolution energy mesh near Ef for metal. 
So it may be better to set
HistBin_dw    1d-5 ! 1d-5 is fine mesh (good for metal?) !(a.u.) BinWidth along real axis at omega=0.
HistBin_ratio 1.03 ! 1.03 maybe safer. frhis(iw)= b*(exp(a*(iw-1))-1), where a=ratio-1.0 and dw=b*a
GaussianFilterX0 0.01 !a.u. this is smearing fo X0 !this setting might be irrelevant.




### -- QSGW: how to check convergence
QSGW calculation contains (1) and (2)
  (1). One-body self-consistent calculation 
      (where we add sigm = Sigma-Vxc^LDA to one-body potential).
      H_0 is determined.
  (2). For given H0, we calculate sigm file.

Big iteration cycle of QSGW is made from (1)+(2).
(gwsc script. not run_arg is a subroutine of bash script) 
With (1), we have small iteration cycle of one-body calculaiton with keeping given sigm.

In save.*, we see total energy (but not the total energy in the QSGW
mode), a line per each iteration of (1). A line "c ..." is the final
iteration cycle of (1)."x ..." is unconverged 
(but no problem as long as we finally see "c ...").

The command "grep '[cx] ' save.*" gives an indicator for 
going to be converged or not.
Or you can take "grep gap llmf.*run" (see it bottom.)

Another way:
$~/ecalj/TestInstall/bin/diffnum QPU.3run QPU.6run 
is to compare two QPU files which contains QP energies.
(note: QP energies shown are calculated just at the begininig of iteration).

For insulater, (I think), comparing band gap for each iteration 
is good enough to check onvergence. But for metal, it is better to plot energy bands
for some of final iterations, and overlapped (cd RUN.ITER* and run
job_band).

Another way is
>grep rms lqpe*

This gives rmsdel. Diffence of self-energy
(at least we see it is getting smaller for initial first cycles). 

### -- How to make 80%QSGW +20% LDA, and SO setting
  Note that sigm file contains Vxc^QSGW-Vxc^LDA.
  If sigm exists, lmf read it, and run self-consistent calculations
  with adding sigm to the one-body potential.

  See TableII in 
  https://iopscience.iop.org/article/10.7567/JJAP.55.051201/pdf

###### 1. QSGW80(NoSC)
  For practical prediction of band structure, such as band gap and so
  on, it may be better to use 80% QSGW +20% LDA procedure when you
  make band plot. 
  After, you have rst and sigm files determined self-consistently
  Run 
  ```job_band gaas -np 4 -vssig=0.80 ```
  (check ssig is defined and cited as
      ScaledSigma={ssig} 
   in the ctrl file).
  This gives a result of QSGW80nosc in the TableII.

###### 2. QSGW80(Nosc)+SO 
   80%QSGW+20%LDA with SO=1 (L.S method).
   If you like to include L.S method 

```mpirun lmf gaas -np 4 -vssig=0.80 -vso=1 -vnspin=2```
    this procedure makes self-consistency with keeping the sigm file. This may/(or may not) required.
     If you expect large obital moment this procedure may be needed.)

```job_band gaas -np 4 -vssig=0.80 -vso=1 -vnspin=2```

  NOTE: nspin=2 is required for so=1
  rst and sigm are expanded for npsin=2 (you can not run nspin=2, after rst and sigm are expanded).
  
###### 3. QSGW80
  With ssig=0.80, you can run QSGW calculaiton in gwsc.
  Then you have self-consistent results of QSGW80.
  You can simultaneously use the setting so=2 (Lz.Sz scheme).
  Be careful for z-direction and setting of SYMOPS (so as to keep the z-axis), for so=2.
  If you like to get results of QSGW80+SO, you need to set so=1 after self-consistent of 
  sigm atteined.

###### 4. Example of GaAs   
  Good example to check band gap, and SO splitting at top of valence of Gamma point for
  ZB structure as GaAs.
  Before run it, make sure your ctrl file include variables ssig, so,
  nspin by   
  >grep ssig ctrl.gaas
  >grep so   ctrl.gaas
  >grep nspin ctrl.gaas

  to know the variable ssigm is defined and used as
  ScaledSigma={ssig}, NSPIN={nspin}.
  For -vso=1 work, you also need to so is defined and SO={so} is set.


### -- ecalj/MATERIALS/ 
At ~/ecalj/MATERIALS/, run ./job_materials.py
  It shows a help with a list of materials.
  It contains samples of simple materials.
  It performs LDA calculations and generates GWinput for materials.
 This job_materials.py works as follows for given material names.
 Step 1. Generate ctrls.* file for Materials.ctrls.database. (names are in DATASECTION:)
 Step 2. Generate ctrl by ctrlgenM1.py
 Make directtory such as Si/

We have deguchi paper 
https://sci-hub.tw/https://doi.org/10.7567/JJAP.55.051201
All calculation is by the default setting in QSGW on the PMT method.
   No empty spheres. EH=-1,EH=-2, MT radius is -3% untouching.
   RSMH=RSMH2=R/2


### -- ecalj/Samples

We have samples of 
 
* Magnon: magnon 
   now with MaxlocWannier. going to move to MLO
* MLO: new localized basis
  Maximally localized Wannier function and cRPA interaction
* CuepsPP 
Dielectric function. We will improved this 

* MgO_PROCAR : PROCAR generation sample
 Run jobprocar. This gives *.eps file which shows Fat band picture.
 PROCAR (vasp format) is generated and analysed by a script BandWeight.py.

* MLWF_samples
  Wannier function generator and cRPA 
  wannier90 method implemented in ecalj and cRPA. 
  (a cRPA method by Juelich group).
  See Samples_MLWF/README.
   ~/ecalj/README_wannier.md This is going to be replaced with README_MLO.md

* mass_fit_test
  Effective mass calculation. See README. probably not maitained\dots

* IIR
  impact ionization rate

In addition, we have Boltztorap, FermiSurface, LaGaO3_relax(in LDA/GGA), ReN(cubic) samples\dots
We are going to organize samples, but not clean yet. Ask us.

### ecalj/ecalj_auto
massive calculations for given sets of POSCAR files.