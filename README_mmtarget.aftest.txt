QSGW+{B}
Memo for NiSe 02jul2021

=== Install ===
Please git clone it to your directory.
If you are not so familiar with github,
move old ecalj to ecalj_backup or something.
Then run 'git clone' from bitbucket.


== QSGW+{B} ===
We can run QSGW+{B} calculation. We need to follow instruction below.
Here B is the alternative magnetic field to keep the magnetic moment
with keeping the antiferro order. Antiferro version of the Fixed moment method.

Current code is only for crystal structure of Antiferro-magnetic NiSe
(because of minor reason not woking for NiO).

If you set 0.1 in the file 'mmtarget.aftest.',
lmf-MPIK try to find out U values (J=0),
to keep the size of magnetic moment 0.1 (and -0.1 at down sites).

When you succeeded QSGW+{B} calcultion, you can get the size of B field,
so as to keep the specified magnetic moments, given in mmtarget.aftest.
After iteration cycle become stable (converged),
you can plot band structures, and so on.

=== Setting ===
1. ctrl.nise:
---
AF magnetic pair setting (AF=
  ATOM=Niup POS=     0.00       0.000   0.000 AF=1
  ATOM=Nidn POS=     0.00       0.000   0.500 AF=-1
In addtion,
  SYMGRPAF  i:(0,0,1/2)
(this means inversion with (0,0,1/2) translation gives AF symmetry).
---
IDU=0 0 2
UH=0 0 0
JH=0 0 0
 for Niup and for Nidn (FLL mode of LDA+U).
 You can supply UH and JH if you like to run QSGW+U+{B} or something.
 But we did not test it yet.
---
ScaledSigma=0.8
 to rmove overshoot of QSGW.
 https://iopscience.iop.org/article/10.7567/JJAP.55.051201
---
BZ_FSMOMMETHOD=0
BZ_FSMOM=0.0
 These are to enforce no total magnetic moment. Not necessary but may make stable calculation
---
HAM_FRZWF=T
 This fix the boundary conditions within MuffinTin to calculate radial wave functions.

---
SPEC_ATOM_P, SPEC_ATOM_PZ
 Since we set HAM_FRZWF, we need to set them reasonable values. Use values in kotani's example.
 

2. GWinput:
Settings to make calculations smoother.
--------
mixbeta 0.3            <-- self-energy mixing for stabilize convergence.
GaussianFilterX0 0.05  <-- smoothing the density of state for the polarization fuction chi o
esmr   0.03            <-- smoothing the pole of G in G \times W.

You can use my setting
Try n1n2n3 4 4 4, at first.
If you like try 5 5 5 to check stability---but 4 4 4 is probably not so bad.

3. mmtarget.aftest:
----------------
Put a file 'mmtarget.aftest' which is specifing the magnetic moment. (sigle line).
When a file mmtarget.aftest exist, lmf-MPIK works in the "AFTEST" mode.
lmf-MPIK is adding bias field so as to keep the magnetic moment given in mmtrget.aftest.
Use gwsc_sym (included) instead of gwsc. Probably this is necessary to perform numerically stable
converence.

4. Required commandline arguments.(this is automatic when you use gwsc_sym calling run_arg2 in it).
--phispinsym : spin symmetrized SPEC_ATOM_P,PZ.
               This is probably not necessary but use this makes things safer.
--v0fix  : fix the radial potentail to generate radial functions.
           At the beginig of lmf-MPIK, it generates v0pot.xxx files containing the potential.
           When v0pot* file exist, lmf-MPIK read v0pot* files.

If you use gwsc_sym, which callsed run_arg2,
it automatically includes --phispinsym and --v0fix when runnnig gwsc.

-------------------------
WARN!!!: You have to use run_arg2 even in job_band and so on.
test it please.
-------------------------


=== Convergence check ===================
We have two convergence.
As you know, QSGW made from two iteration cycle.
Big iteration is GW cycle (*run).
lmf-MPIK makes an inner iteration cycle
(lmf cyle; each line in save.nise. It is shown in llmf.*run).

* Check of lmf cycle (inner cycle for given sigm file):
>grep -n -A2 -B2 'From last ' llmf
is to monitor iteration cycle of lmf-MPIK.
For example, it shows
...
 2653-   it  8  of 30    ehf=   -1195.278537   ehk=   -1195.278537
 2654: From last iter    ehf=   -1195.278540   ehk=   -1195.278537
 2655- diffe(q)=  0.000003 (0.000013)    tol= 0.000010 (0.000010)   more=T
 2656-i ehf=-1195.2785371 ehk=-1195.2785374
...
Here diffe(q)= 0.000003 (0.000013) show energy and density difference from previous lmf iteration.
When they become smaller than tol. It is judged as converged. (CONV and CONVC in ctrl).
In NiSe, diffe is slight larger than cliteron and going to run until NIT (upper limit).
But no problem in QSGW usually
(for purpose to get band structure--- we may use a little loose CONV CONVC probably).

* magnetic moment check of lmf cycle (inner cycle for give sigm file):
We need to check the magnetic moments
are really controlled to be the given number in mmtarget.aftest in lmf.
To check it, we need to do 
> grep -A4 'true mm' llmf
(or grep -A4 'true mm' llmf_lda')
Then read "true mm". If you set 0.2 in mmtarget.aftest, you will see
---------
mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    8.748453    2.309224    6.439229      0.200175    0.022252    0.177924
   2    8.748458    2.309611    6.438847     -0.199945   -0.021957   -0.177988
---------
You see true mm are almost 0.2 (and -0.2) for ibas=1 (Niup) and ibas=2 (Nidn).

*mmagfield.aftest
This contains magnetic field to controll the size of magnetic moment.
Right after lmf-MPIK finished, mmagfield.aftest is essentially deleted (only single line).

* magnetic moment check of RUN cycle.
Read and check the script hmmom (inlcuded), and run
./hhmom
It shows like 
...
xxxxxxxxx
run ix=,160 ,uhval: 37 0.002732 0.199157 -0.200829 -0.005411 0.004808
 run ix=,161 ,uhval: 72 0.003801 0.200180 -0.199685 -0.003650 0.003424
 run ix=,162 ,uhval: 21 0.005352 0.199403 -0.200113 -0.003975 0.003251
 run ix=,163 ,uhval: 45 0.007362 0.199974 -0.199989 -0.002940 0.002785
 run ix=,164 ,uhval: 30 0.008642 0.199926 -0.200228 -0.002124 0.002276
 run ix=,165 ,uhval: 24 0.009568 0.200270 -0.199745 -0.002666 0.002538
 run ix=,166 ,uhval: 58 0.008485 0.200110 -0.199914 -0.002824 0.002748
 run ix=,167 ,uhval: 100 0.008020 0.202850 -0.197137 -0.001750 0.002638
 run ix=,168 ,uhval: 36 0.008810 0.199701 -0.200316 -0.001170 0.001748
 run ix=,169 ,uhval: 80 0.007699 0.201160 -0.198861 -0.001836 0.002055
...
Check the content of the script hhmom.
Here ix is the number of RUN cycle (so I run so many RUN cycle 'gwsc_sym 180 -np 32 nise')
See the line ix=,169. 0.007699 shows the U value to keep magnetic moment to be 0.2.
0.201160 -0.198861 -0.001836 0.002055 are magnetic moment as Niup, Nidn, Se1, Se2 (last resulst of lmf cycle).
So, manetic moments are well controlled to be ~0.2.
But we see UH value (U value) is changing from 0.002732 to 0.007699.
This implies instability of bias field so as to keep the magnetic moment=0.2.

grep mmaftest llmf, shows current converging status (see hmmom).


* grep 'true mm' llmf -A4, show how lmf-MPI is going to be converged
  as for magnetic moments for given sigm.nise

* save.nise check:
save.nise is not showing total energy in QSGW; just take
it as an indicator for convergence. See save.nise. Take "grep  '[xc] ' save.nise".


* QPU check:
Finally, in principle, we need to confirm QPU.*run becomes almost stable.
