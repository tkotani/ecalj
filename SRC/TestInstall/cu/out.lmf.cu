INFO: Ubuntu 20.04.4 LTS \n \l
INFO: GNU Fortran (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0
INFO: -O2 -g -fimplicit-none -finit-integer=NaN -finit-real=NaN -JOBJ.gfortran -IOBJ.gfortran
INFO: MATH: -lmkl_rt
INFO: git: commit 913e769c0a5a77a2254ce7ce7011c5bc7fb5168a
INFO:    : Author: Takao Kotani <takaokotani@gmail.com>
INFO:    : Date:   Mon Feb 13 19:43:59 2023 +0900
INFO: linked at Tue Feb 14 12:44:43 JST 2023
=== START LFMA ===
 mpisize=           1
 m_lmfinit:program LMFA
 cmdl for python=/home/takao/ecalj/SRC/TestInstall/bin/ctrl2ctrlp.py  --no-iactiv cu -vnk=8<ctrl.cu >ctrlp.cu
rval2: STRUC_NSPEC             requ n= 1 val= 1.00000000
rval2: STRUC_NBAS              requ n= 1 val= 1.00000000
rval2: HAM_NSPIN               defa n= 1 val= 1.00000000
rval2: IO_VERBOS               defa n= 1 val= 31.00000000
rval2: IO_TIM                  defa n= 1 val= 1.00000000
rval2: STRUC_ALAT              ---- n= 1 val= 6.79800000
rval2: STRUC_DALAT             ---- n= 1 val= 0.00000000
rval2: STRUC_PLAT              requ n= 9 val= 0.00000000  0.50000000  0.50000000  0.50000000  0.00000000  0.50000000  0.50000000  0.50000000  0.00000000
rval2: OPTIONS_HF              defa n= 1 val= 0.00000000
rval2: HAM_REL                 defa n= 1 val= 1.00000000
rval2: HAM_SO                  defa n= 1 val= 0.00000000
rval2: HAM_SOCAXIS             defa n= 3 val= 0.00000000  0.00000000  1.00000000
rval2: HAM_GMAX                defa n= 1 val= 9.00000000
rval2: HAM_FTMESH              defa n= 3 val= 10.00000000  10.00000000  10.00000000
rval2: HAM_TOL                 defa n= 1 val= 0.00000100
rval2: HAM_FRZWF               defa n= 1 val= 0.00000000
rval2: HAM_XCFUN               defa n= 1 val= 2.00000000
rval2: HAM_FORCES              defa n= 1 val= 0.00000000
rval2: HAM_RDSIG               defa n= 1 val= 1.00000000
rval2: HAM_ScaledSigma         defa n= 1 val= 1.00000000
rval2: HAM_EWALD               defa n= 1 val= 0.00000000
rval2: HAM_OVEPS               defa n= 1 val= 0.00000000
rval2: HAM_PWMODE              defa n= 1 val= 0.00000000
rval2: HAM_PWEMAX              defa n= 1 val= 3.00000000
rval2: HAM_READP               defa n= 1 val= 0.00000000
rval2: HAM_V0FIX               defa n= 1 val= 0.00000000
rval2: HAM_PNUFIX              defa n= 1 val= 0.00000000
rval2: SPEC_Z@1                ---- n= 1 val= 29.00000000
rval2: SPEC_R@1                ---- n= 1 val= 2.31127105
rval2: SPEC_R/W@1              ---- n= 0 val= 
rval2: SPEC_R/A@1              ---- n= 0 val= 
rval2: SPEC_A@1                defa n= 1 val= 0.02500000
rval2: SPEC_NR@1               defa n= 1 val= 0.00000000
rval2: SPEC_RSMH@1             ---- n= 3 val= 2.50000000  2.50000000  1.00000000
rval2: SPEC_EH@1               requ n= 3 val= -0.01000000 -0.01000000 -0.01000000
rval2: SPEC_RSMH2@1            ---- n= 5 val= 0.00000000  0.00000000  0.00000000  0.00000000  0.00000000
rval2: SPEC_EH2@1              requ n= 0 val= 
rval2: SPEC_LMX@1              defa n= 1 val= 3.00000000
rval2: SPEC_LMXA@1             defa n= 1 val= 3.00000000
rval2: SPEC_LMXL@1             defa n= 1 val= 3.00000000
rval2: SPEC_P@1                ---- n= 4 val= 4.65000000  4.34000000  3.87000000  4.11000000
rval2: SPEC_Q@1                ---- n= 0 val= 
rval2: SPEC_NMCORE@1           defa n= 1 val= 0.00000000
rval2: SPEC_PZ@1               ---- n= 3 val= 5.50000000  5.50000000  4.50000000
rval2: SPEC_LFOCA@1            defa n= 1 val= 1.00000000
rval2: SPEC_KMXA@1             defa n= 1 val= 4.00000000
rval2: SPEC_RSMA@1             defa n= 1 val= 0.92450842
rval2: SPEC_IDMOD@1            ---- n= 5 val= 0.00000000  0.00000000  0.00000000  1.00000000  1.00000000
rval2: SPEC_FRZWF@1            defa n= 1 val= 0.00000000
rval2: SPEC_IDU@1              ---- n= 0 val= 
rval2: SPEC_UH@1               ---- n= 0 val= 
rval2: SPEC_JH@1               ---- n= 0 val= 
rval2: SPEC_C-HQ@1             defa n= 2 val= -1.00000000  0.00000000
rval2: SPEC_EREF1              defa n= 1 val= 0.00000000
rval2: SITE_POS@1              ---- n= 3 val= 0.00000000  0.00000000  0.00000000
rval2: SITE_RELAX@1            defa n= 3 val= 1.00000000  1.00000000  1.00000000
rval2: SITE_AF@1               defa n= 1 val= 0.00000000
rval2: STR_RMAXS               ---- n= 0 val= 
rval2: STR_RMAX                ---- n= 0 val= 
rval2: STR_MXNBR               defa n= 1 val= 0.00000000
rval2: BZ_NKABC                ---- n= 1 val= 8.00000000
rval2: BZ_BZJOB                ---- n= 1 val= 1.00000000
rval2: BZ_METAL                defa n= 1 val= 2.00000000
rval2: BZ_TETRA                defa n= 1 val= 1.00000000
rval2: BZ_N                    defa n= 1 val= 0.00000000
rval2: BZ_W                    defa n= 1 val= 0.00200000
rval2: BZ_ZBAK                 defa n= 1 val= 0.00000000
rval2: BZ_SAVDOS               defa n= 1 val= 1.00000000
rval2: BZ_NPTS                 defa n= 1 val= 1001.00000000
rval2: BZ_DOSMAX               defa n= 1 val= 2.93992268
rval2: BZ_EFMAX                defa n= 1 val= 5.00000000
rval2: BZ_NEVMX                defa n= 1 val= 0.00000000
rval2: BZ_FSMOM                defa n= 1 val= -99999.00000000
rval2: BZ_FSMOMMETHOD          defa n= 1 val= 0.00000000
rval2: EWALD_AS                defa n= 1 val= 2.00000000
rval2: EWALD_TOL               defa n= 1 val= 0.00000000
rval2: EWALD_NKDMX             defa n= 1 val= 600.00000000
rval2: ITER_NIT                defa n= 1 val= 12.00000000
rval2: ITER_NRMIX              defa n= 1 val= 80.00000000
rval2: ITER_CONV               defa n= 1 val= 0.00001000
rval2: ITER_CONVC              defa n= 1 val= 0.00001000
rval2: ITER_UMIX               defa n= 1 val= 0.50000000
rval2: ITER_TOLU               defa n= 1 val= 0.00000000
rval2: DYN_MODE                defa n= 1 val= 0.00000000
rval2: DYN_NIT                 defa n= 1 val= 1.00000000
rval2: DYN_HESS                defa n= 1 val= 1.00000000
rval2: DYN_XTOL                defa n= 1 val= 0.00100000
rval2: DYN_GTOL                defa n= 1 val= 0.00000000
rval2: DYN_STEP                defa n= 1 val= 0.01500000
rval2: DYN_NKILL               defa n= 1 val= 0.00000000
mixing parameters: A/B nmix wt: 0 -1 1.000000  1.000000  0.000000 beta elin wc killj=  1.000000 -1.000000 -1
 ===> for --jobgw, pwmode is switched to be  0

 ... Species  1
  bndfp (warning): no sigm file found ... LDA calculation only
pnuall: j isp pnu= 1 1 4.650000  4.340000  3.870000  4.110000
pnzall: j isp  pz= 1 1 5.500000  5.500000  4.500000  0.000000


mmm === MTO setting ===
mmm ispec lmxb lpz nkapii nkaphh=    1    2    1    1    1
mmm rsmh1    1  2.50  2.50  1.00
mmm   eh1    1 -0.01 -0.01 -0.01
mmm pz       1  5.50  5.50  4.50
mmm lh       2  2
xxx isp pz= 1 5.500000  5.500000  4.500000  0.000000
 goto freats

conf:------------------------------------------------------
conf:SPEC_ATOM= A : --- Table for atomic configuration ---
conf:  isp  l  int(P) int(P)z    Qval     Qcore   CoreConf
conf:    1  0       4  5         1.000    6.000 => 1,2,3,
conf:    1  1       4  5         0.000   12.000 => 2,3,
conf:    1  2       3  4        10.000    0.000 => 
conf:    1  3       4  0         0.000    0.000 => 
usedQ=     1.000     0.000    10.000     0.000
conf: Species  A        Z=  29.00 Qc=  18.000 R=  2.311271 Q=  0.000000 nsp= 1 mom=  0.000000
conf: rmt rmax a=  2.311271  48.805862  0.025000 nrmt nr= 393 515
 goto atomc xxx
 atomsc nmcore=           0

 end of atomsc xxxxx
 vsum=  -130.79144076069792                1
sumev= -4.333254 etot= -3304.416258 eref=  0.000000 etot-eref= -3304.416258

 Free-atom wavefunctions:
 valence:      eval       node at      max at       c.t.p.   rho(r>rmt)       pnu
   4s      -0.36411         0.890       2.256       3.582     0.643062       4.761  0
   5s      -0.00028         3.669      10.794      19.873     0.990448       5.848  1
   4p      -0.06295         0.975       3.484       7.414     0.901829       4.561  0
   5p       0.00796         6.760      30.414      48.806*    0.999240       5.593  1
   3d      -0.39691         0.000       0.600       3.429     0.056076       3.888  0
   4d       0.01308         1.868      33.290      48.806*    0.999995       4.148  1
   4f       0.01948         0.000      35.393      48.806*    1.000000       4.137  0

 core:        ecore       node at      max at       c.t.p.   rho(r>rmt)
   1s    -649.07634         0.000       0.034       0.069     0.000000
   2s     -77.91382         0.070       0.197       0.308     0.000000
   2p     -67.32532         0.000       0.158       0.335     0.000000
   3s      -8.39248         0.288       0.614       0.895     0.000141
   3p      -5.29682         0.260       0.619       1.078     0.000727
 tailsm: init

 tailsm: fit tails to 6 smoothed hankels, rmt= 2.31127, rsm= 1.15564
  ---E:energies of smHankels. C:fitting coeeficient for core tail. ---
 E:    -1.00000    -2.00000    -4.00000    -6.00000    -9.00000   -15.00000
 C:    -0.07160    10.75053  -187.49213  1222.02349 -4717.78530 21166.80769
        r          rho         fit         diff
    2.311271    0.017797    0.017766    0.000031
    2.967767    0.005662    0.005658    0.000005
    3.810725    0.001517    0.001518   -0.000001
    4.893104    0.000305    0.000305   -0.000000
    6.282906    0.000041    0.000041   -0.000001
    8.067448    0.000003    0.000003    0.000000
    q(fit):     1.203836    rms diff:   0.000016
    fit: r>rmt  1.203836   r<rmt  3.442816   qtot  4.646652
    rho: r>rmt  1.203836   r<rmt  9.796164   qtot 11.000000
 tailsm:  fit tails to        6 functions with

 rsm=  0.11556D+01 rms error=  0.16285D-04
conf: Core rhoc(rmt)= 0.003922 spillout= 0.004646
 Fit with Hankel e=-24.082483 coeff=764.352513
      r            rhoc          fit
    2.311271    0.02095279    0.02095279
    2.429779    0.01229068    0.01231367
    2.753317    0.00285262    0.00285190
    3.119934    0.00054243    0.00053465
    3.535366    0.00008235    0.00007888
    4.006112    0.00000969    0.00000887
    4.539536    0.00000085    0.00000073
    5.143985    0.00000005    0.00000004
 end of freats: spid nmcore=A                  0
Sum of reference energies:                      0.000000000000
OK! end of LMFA ======================
INFO: Ubuntu 20.04.4 LTS \n \l
INFO: GNU Fortran (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0
INFO: -O2 -g -fimplicit-none -finit-integer=NaN -finit-real=NaN -JOBJ.gfortran -IOBJ.gfortran
INFO: MATH: -lmkl_rt
INFO: git: commit 913e769c0a5a77a2254ce7ce7011c5bc7fb5168a
INFO:    : Author: Takao Kotani <takaokotani@gmail.com>
INFO:    : Date:   Mon Feb 13 19:43:59 2023 +0900
INFO: linked at Tue Feb 14 12:44:43 JST 2023
 m_lmfinit:program LMF
===START LMF with   ===
 mpisize=           4
 m_lmfinit:program LMF
 m_lmfinit:program LMF
 m_lmfinit:program LMF
 cmdl for python=/home/takao/ecalj/SRC/TestInstall/bin/ctrl2ctrlp.py  --no-iactiv cu -vnk=8 -vbigbas=f<ctrl.cu >ctrlp.cu
rval2: STRUC_NSPEC             requ n= 1 val= 1.00000000
rval2: STRUC_NBAS              requ n= 1 val= 1.00000000
rval2: HAM_NSPIN               defa n= 1 val= 1.00000000
rval2: IO_VERBOS               defa n= 1 val= 31.00000000
rval2: IO_TIM                  defa n= 1 val= 1.00000000
rval2: STRUC_ALAT              ---- n= 1 val= 6.79800000
rval2: STRUC_DALAT             ---- n= 1 val= 0.00000000
rval2: STRUC_PLAT              requ n= 9 val= 0.00000000  0.50000000  0.50000000  0.50000000  0.00000000  0.50000000  0.50000000  0.50000000  0.00000000
rval2: OPTIONS_HF              defa n= 1 val= 0.00000000
rval2: HAM_REL                 defa n= 1 val= 1.00000000
rval2: HAM_SO                  defa n= 1 val= 0.00000000
rval2: HAM_SOCAXIS             defa n= 3 val= 0.00000000  0.00000000  1.00000000
rval2: HAM_GMAX                defa n= 1 val= 9.00000000
rval2: HAM_FTMESH              defa n= 3 val= 10.00000000  10.00000000  10.00000000
rval2: HAM_TOL                 defa n= 1 val= 0.00000100
rval2: HAM_FRZWF               defa n= 1 val= 0.00000000
rval2: HAM_XCFUN               defa n= 1 val= 2.00000000
rval2: HAM_FORCES              defa n= 1 val= 0.00000000
rval2: HAM_RDSIG               defa n= 1 val= 1.00000000
rval2: HAM_ScaledSigma         defa n= 1 val= 1.00000000
rval2: HAM_EWALD               defa n= 1 val= 0.00000000
rval2: HAM_OVEPS               defa n= 1 val= 0.00000000
rval2: HAM_PWMODE              defa n= 1 val= 0.00000000
rval2: HAM_PWEMAX              defa n= 1 val= 3.00000000
rval2: HAM_READP               defa n= 1 val= 0.00000000
mixing parameters: A/B nmix wt: 0 -1 1.000000  1.000000  0.000000 beta elin wc killj=  1.000000 -1.000000 -1

 ... Species  1
rval2: HAM_V0FIX               defa n= 1 val= 0.00000000
rval2: HAM_PNUFIX              defa n= 1 val= 0.00000000
rval2: SPEC_Z@1                ---- n= 1 val= 29.00000000
rval2: SPEC_R@1                ---- n= 1 val= 2.31127105
rval2: SPEC_R/W@1              ---- n= 0 val= 
rval2: SPEC_R/A@1              ---- n= 0 val= 
rval2: SPEC_A@1                defa n= 1 val= 0.02500000
rval2: SPEC_NR@1               defa n= 1 val= 0.00000000
rval2: SPEC_RSMH@1             ---- n= 3 val= 2.50000000  2.50000000  1.00000000
rval2: SPEC_EH@1               requ n= 3 val= -0.01000000 -0.01000000 -0.01000000
mixing parameters: A/B nmix wt: 0 -1 1.000000  1.000000  0.000000 beta elin wc killj=  1.000000 -1.000000 -1
mixing parameters: A/B nmix wt: 0 -1 1.000000  1.000000  0.000000 beta elin wc killj=  1.000000 -1.000000 -1

 ... Species  1

 ... Species  1
rval2: SPEC_RSMH2@1            ---- n= 5 val= 0.00000000  0.00000000  0.00000000  0.00000000  0.00000000
rval2: SPEC_EH2@1              requ n= 0 val= 
rval2: SPEC_LMX@1              defa n= 1 val= 3.00000000
rval2: SPEC_LMXA@1             defa n= 1 val= 3.00000000
rval2: SPEC_LMXL@1             defa n= 1 val= 3.00000000
rval2: SPEC_P@1                ---- n= 4 val= 4.65000000  4.34000000  3.87000000  4.11000000
rval2: SPEC_Q@1                ---- n= 0 val= 
rval2: SPEC_NMCORE@1           defa n= 1 val= 0.00000000
rval2: SPEC_PZ@1               ---- n= 3 val= 5.50000000  5.50000000  4.50000000
rval2: SPEC_LFOCA@1            defa n= 1 val= 1.00000000
rval2: SPEC_KMXA@1             defa n= 1 val= 4.00000000
rval2: SPEC_RSMA@1             defa n= 1 val= 0.92450842
rval2: SPEC_IDMOD@1            ---- n= 5 val= 0.00000000  0.00000000  0.00000000  1.00000000  1.00000000
rval2: SPEC_FRZWF@1            defa n= 1 val= 0.00000000
rval2: SPEC_IDU@1              ---- n= 0 val= 
rval2: SPEC_UH@1               ---- n= 0 val= 
rval2: SPEC_JH@1               ---- n= 0 val= 
rval2: SPEC_C-HQ@1             defa n= 2 val= -1.00000000  0.00000000
rval2: SPEC_EREF1              defa n= 1 val= 0.00000000
rval2: SITE_POS@1              ---- n= 3 val= 0.00000000  0.00000000  0.00000000
rval2: SITE_RELAX@1            defa n= 3 val= 1.00000000  1.00000000  1.00000000
rval2: SITE_AF@1               defa n= 1 val= 0.00000000
rval2: STR_RMAXS               ---- n= 0 val= 
rval2: STR_RMAX                ---- n= 0 val= 
rval2: STR_MXNBR               defa n= 1 val= 0.00000000
rval2: BZ_NKABC                ---- n= 1 val= 8.00000000
rval2: BZ_BZJOB                ---- n= 1 val= 1.00000000
rval2: BZ_METAL                defa n= 1 val= 2.00000000
rval2: BZ_TETRA                defa n= 1 val= 1.00000000
rval2: BZ_N                    defa n= 1 val= 0.00000000
rval2: BZ_W                    defa n= 1 val= 0.00200000
rval2: BZ_ZBAK                 defa n= 1 val= 0.00000000
rval2: BZ_SAVDOS               defa n= 1 val= 1.00000000
rval2: BZ_NPTS                 defa n= 1 val= 1001.00000000
rval2: BZ_DOSMAX               defa n= 1 val= 2.93992268
rval2: BZ_EFMAX                defa n= 1 val= 5.00000000
rval2: BZ_NEVMX                defa n= 1 val= 0.00000000
rval2: BZ_FSMOM                defa n= 1 val= -99999.00000000
rval2: BZ_FSMOMMETHOD          defa n= 1 val= 0.00000000
rval2: EWALD_AS                defa n= 1 val= 2.00000000
rval2: EWALD_TOL               defa n= 1 val= 0.00000000
rval2: EWALD_NKDMX             defa n= 1 val= 600.00000000
rval2: ITER_NIT                defa n= 1 val= 12.00000000
rval2: ITER_NRMIX              defa n= 1 val= 80.00000000
rval2: ITER_CONV               defa n= 1 val= 0.00001000
rval2: ITER_CONVC              defa n= 1 val= 0.00001000
rval2: ITER_UMIX               defa n= 1 val= 0.50000000
rval2: ITER_TOLU               defa n= 1 val= 0.00000000
rval2: DYN_MODE                defa n= 1 val= 0.00000000
rval2: DYN_NIT                 defa n= 1 val= 1.00000000
rval2: DYN_HESS                defa n= 1 val= 1.00000000
rval2: DYN_XTOL                defa n= 1 val= 0.00100000
rval2: DYN_GTOL                defa n= 1 val= 0.00000000
rval2: DYN_STEP                defa n= 1 val= 0.01500000
rval2: DYN_NKILL               defa n= 1 val= 0.00000000
mixing parameters: A/B nmix wt: 0 -1 1.000000  1.000000  0.000000 beta elin wc killj=  1.000000 -1.000000 -1
 ===> for --jobgw, pwmode is switched to be  0

 ... Species  1
  bndfp (warning): no sigm file found ... LDA calculation only
pnuall: j isp pnu= 1 1 4.650000  4.340000  3.870000  4.110000
pnzall: j isp  pz= 1 1 5.500000  5.500000  4.500000  0.000000
 imx=           3           3           3
 imx=           3           3           3


mmm === MTO setting ===
mmm ispec lmxb lpz nkapii nkaphh=    1    2    1    1    1
mmm rsmh1    1  2.50  2.50  1.00
mmm   eh1    1 -0.01 -0.01 -0.01
mmm pz       1  5.50  5.50  4.50
mmm lh       2  2
 imx=           3           3           3
 imx=           4           4           4
 imx=           4           4           4

                Plat                                  Qlat
   0.000000   0.500000   0.500000       -1.000000   1.000000   1.000000
   0.500000   0.000000   0.500000        1.000000  -1.000000   1.000000
   0.500000   0.500000   0.000000        1.000000   1.000000  -1.000000
 imx=           4           4           4
  Cell vol=   78.538660
 imx=           3           3           3
 imx=           4           4           4

 LATTC:  as= 2.000   tol= 1.00E-12   alat= 6.79800   awald= 0.467
         r1=  2.252   nkd= 201      q1=  6.871   nkg= 331
SGROUP:  1 symmetry operations from 0
 SYMLAT: Bravais system is cubic       with 48 symmetry operations.
 SYMCRY: crystal invariant under  48 symmetry operations for tol=  0.000100
 ig  group op
   1  i*i
   2  i
   3  r3(1,1,-1)
   4  i*r3(1,1,-1)
   5  r3(-1,-1,1)
   6  i*r3(-1,-1,1)
   7  r3d
   8  i*r3d
   9  r3(-1,-1,-1)
  10  i*r3(-1,-1,-1)
  11  r2x
  12  mx
  13  r4x
  14  i*r4x
  15  r4(-1,0,0)
  16  i*r4(-1,0,0)
  17  r3(1,-1,-1)
  18  i*r3(1,-1,-1)
  19  r3(-1,1,1)
  20  i*r3(-1,1,1)
  21  r2(1,1,0)
  22  m(1,1,0)
  23  r2(1,0,-1)
  24  m(1,0,-1)
  25  r2y
  26  my
  27  r4y
  28  i*r4y
  29  r4(0,-1,0)
  30  i*r4(0,-1,0)
  31  r2(0,1,-1)
  32  m(0,1,-1)
  33  r2z
  34  mz
  35  r4(0,0,-1)
  36  i*r4(0,0,-1)
  37  r4z
  38  i*r4z
  39  r3(-1,1,-1)
  40  i*r3(-1,1,-1)
  41  r3(1,-1,1)
  42  i*r3(1,-1,1)
  43  r2(1,0,1)
  44  m(1,0,1)
  45  r2(1,-1,0)
  46  m(1,-1,0)
  47  r2(0,1,1)
  48  m(0,1,1)
 nnnnnn         729         889
GROUPG: the following are sufficient to generate the space group:
 Generators:trans(cart)  = i*r3(1,1,-1) r4x
 Generators::trans(frac) = i*r3(1,1,-1) r4x
MKSYM: found  48  space group operations
SPLCLS: ibas iclass ispec label(ispec)
 SPLCLS     1    1    1     A
 BZMESH:     60 irreducible QP from    8   8   8 shift=TTT
 TETIRR: sorting     3072 tetrahedra ...
     264 inequivalent tetrahedron=
 nnnnnn         729         889
MSHSIZ: mesh has 10 x 10 x 10 divisions; length =     0.481     0.481     0.481
      generated from gmax (a.u.)=      9.0000: 889 vectors of 1000 (88%)
 SGVSYM: 38 symmetry stars found for 861 reciprocal lattice vectors

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol= 1.00E-06
 spec      l    rsm    eh     gmax    last term    cutoff
  A        0*   2.50  -0.01   2.974    2.30E-05      27 
  A        1    2.50  -0.01   3.093    1.29E-06      51 
  A        2    1.00  -0.01   8.508    1.16E-06     813 
 m_qplistinit:start
 nnnnnn         729         889
 nnnnnn         729         889

 iors  : read rst restart file (binary mesh density)
 iors  : empty file ... nothing read

rdovfa: read and overlap free-atom densities (mesh density) ...
 rdovfa: expected A,       read A        with rmt=  2.3113  mesh   393  0.025
  ovlpfa: overlap smooth part of FA densities
 site 1 spec 1 pos 0.0000  0.0000  0.0000 Qsmooth 4.6466523386126539 mom 4.6466523386126539
 total smooth Q =  4.6466523386126539

 Free atom and overlapped crystal site charges:
   ib    true(FA)    smooth(FA)  true(OV)    smooth(OV)    local
    1    9.796164    3.442816   10.275300    3.921952    6.353348

 Smooth charge on mesh:            4.646652
 Sum of local charges:             6.353348
 Total valence charge:            11.000000
 Sum of core charges:             18.000000
 Sum of nuclear charges:         -29.000000
 Homogeneous background:           0.000000
 Deviation from neutrality:       -0.000000

 Basis, after reading restart file
 site spec        pos (Cartesian coordinates)         pos (multiples of plat)
   1         A  0.000000   0.000000   0.000000    0.000000   0.000000   0.000000

--- BNDFP:  begin iteration 1 of 12
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.555104
   smooth rhoves     11.022231   charge     4.646652
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -3.843799 -5.010453 -0.851784
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
     potential shift to crystal energy zero:   -0.000099
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -12.157495      -177.337532      -189.495027
   rhoval*ves            -46.690633      -115.324376      -162.015010
   psnuc*ves              68.735095    -12976.662436    -12907.927341
   utot                   11.022231     -6545.993406     -6534.971175
   rho*exc                -3.843799      -126.414298      -130.258096
   rho*vxc                -5.010453      -167.409316      -172.419769
   valence chg             4.646652         6.353348        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of    60 k=  0.0625  0.0625  0.0625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     2 of    60 k= -0.0625  0.1875  0.1875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     3 of    60 k= -0.1875  0.3125  0.3125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     4 of    60 k= -0.3125  0.4375  0.4375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     5 of    60 k= -0.4375  0.5625  0.5625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     6 of    60 k= -0.5625  0.6875  0.6875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     7 of    60 k= -0.6875  0.8125  0.8125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     8 of    60 k= -0.8125  0.9375  0.9375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     9 of    60 k=  0.0625  0.0625  0.3125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    10 of    60 k= -0.0625  0.1875  0.4375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    11 of    60 k= -0.1875  0.3125  0.5625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    12 of    60 k= -0.3125  0.4375  0.6875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    13 of    60 k= -0.4375  0.5625  0.8125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    14 of    60 k= -0.5625  0.6875  0.9375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    15 of    60 k= -0.6875  0.8125  1.0625 ndimh = nmto+napw =    18   18    0
 ... Done MPI k-loop: elapsed time=   0.0382

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.144466;  11.000000 electrons
         Sum occ. bands:   -0.854464, incl. Bloechl correction: -0.006586
 bndfp:Generating TDOS: efermi=  0.144466  dos window emin emax=  -0.672856  3.084388


 m_bandcal_2nd: to fill eigenfunctions**2 up to Efermi

 mkrout:  Qtrue      sm,loc       local
   1    9.927753    3.113496    6.814257
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -0.854464 Vin*nin=    -189.495027 Ek=Eb-Vin*nin=     188.640563
 Ek(core)=    3171.756639 Exc=    -130.258096 Ees=   -6534.971175 Eharris=   -3304.832069

 mkekin:
   nout*Vin = smpart,onsite,total=:     -6.361301   -168.223682   -174.584983
    E_B(band energy sum)=   -0.854464  E_B-nout*Vin=  173.730519

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.677372
   smooth rhoves     13.178923   charge     4.185743
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -3.054120 -3.974966 -0.866699
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=F
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -7.030715      -175.114087      -182.144802
   rhoval*ves            -50.183969      -106.215150      -156.399119
   psnuc*ves              76.541815    -12962.871088    -12886.329272
   utot                   13.178923     -6534.543119     -6521.364196
   rho*exc                -3.054120      -125.587137      -128.641257
   rho*vxc                -3.974966      -166.302309      -170.277275
   valence chg             4.185743         6.814257        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:      0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=      173.730519 Ekcore=      3171.756639 Ektot    =     3345.487159
 Exc=    -128.641257 Ees   =     -6521.364196 EKohnSham=    -3304.518294
 mixrealsmooth= T
 wgtsmooth=   3.1622776601683791E-002
 mixrho: sought 8 iter from file mixm ; read 0 RMS DQ= 4.34E-2
 AMIX: nmix=0 mmix=8  nelts=  2074  beta=1.00000  tm= 5.00000  rmsdel= 4.34D-02
 mixrho: add corrections to qcell smrho = -0.23414D-07 -0.29812D-09

 iors  : write rst restart file (binary mesh density)

   it  1  of 12    ehf=   -3304.832069   ehk=   -3304.518294
h ehf(eV)=-44964.884163 ehk(eV)=-44960.615006 sev(eV)=-11.625667

--- BNDFP:  begin iteration 2 of 12
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.677372
   smooth rhoves     13.178923   charge     4.185743
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -3.054120 -3.974966 -0.866699
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
     potential shift to crystal energy zero:    0.000553
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -7.030715      -175.114087      -182.144801
   rhoval*ves            -50.183969      -106.215150      -156.399119
   psnuc*ves              76.541815    -12962.871088    -12886.329272
   utot                   13.178923     -6534.543119     -6521.364196
   rho*exc                -3.054120      -125.587137      -128.641257
   rho*vxc                -3.974966      -166.302309      -170.277275
   valence chg             4.185743         6.814257        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of    60 k=  0.0625  0.0625  0.0625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     2 of    60 k= -0.0625  0.1875  0.1875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     3 of    60 k= -0.1875  0.3125  0.3125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     4 of    60 k= -0.3125  0.4375  0.4375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     5 of    60 k= -0.4375  0.5625  0.5625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     6 of    60 k= -0.5625  0.6875  0.6875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     7 of    60 k= -0.6875  0.8125  0.8125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     8 of    60 k= -0.8125  0.9375  0.9375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     9 of    60 k=  0.0625  0.0625  0.3125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    10 of    60 k= -0.0625  0.1875  0.4375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    11 of    60 k= -0.1875  0.3125  0.5625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    12 of    60 k= -0.3125  0.4375  0.6875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    13 of    60 k= -0.4375  0.5625  0.8125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    14 of    60 k= -0.5625  0.6875  0.9375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    15 of    60 k= -0.6875  0.8125  1.0625 ndimh = nmto+napw =    18   18    0
 ... Done MPI k-loop: elapsed time=   0.0588

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.257389;  11.000000 electrons
         Sum occ. bands:   -9.368175, incl. Bloechl correction: -0.013546
 bndfp:Generating TDOS: efermi= -0.257389  dos window emin emax=  -0.960038  2.682534


 mkrout:  Qtrue      sm,loc       local
   1   10.453470    1.685884    8.767586
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -9.368175 Vin*nin=    -182.144801 Ek=Eb-Vin*nin=     172.776627
 Ek(core)=    3171.756639 Exc=    -128.641257 Ees=   -6521.364196 Eharris=   -3305.472187

 mkekin:
   nout*Vin = smpart,onsite,total=:     -3.908741   -233.344935   -237.253675
    E_B(band energy sum)=   -9.368175  E_B-nout*Vin=  227.885500

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.359167
   smooth rhoves      3.889627   charge     2.232413
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -1.459173 -1.896418 -0.723230
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=F
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -2.683250      -210.677503      -213.360752
   rhoval*ves            -33.661409      -148.825746      -182.487155
   psnuc*ves              41.440663    -12997.718880    -12956.278216
   utot                    3.889627     -6573.272313     -6569.382686
   rho*exc                -1.459173      -131.947804      -133.406976
   rho*vxc                -1.896418      -174.700211      -176.596629
   valence chg             2.232413         8.767586        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=      227.885500 Ekcore=      3171.756639 Ektot    =     3399.642140
 Exc=    -133.406976 Ees   =     -6569.382686 EKohnSham=    -3303.147522
 wgtsmooth=   3.1622776601683791E-002
 mixrho: sought 8 iter from file mixm ; read 1 RMS DQ= 1.23E-1  last it= 4.34E-2
 AMIX: nmix=1 mmix=8  nelts=  2074  beta=1.00000  tm= 5.00000  rmsdel= 1.23D-01
   tj: 0.82062
 mixrho: add corrections to qcell smrho = -0.96492D-08 -0.12286D-09

 iors  : write rst restart file (binary mesh density)

   it  2  of 12    ehf=   -3305.472187   ehk=   -3303.147522
 From last iter    ehf=   -3304.832069   ehk=   -3304.518294
 diffe(q)= -0.640118 (0.123166)    tol= 0.000010 (0.000010)   more=T
i ehf(eV)=-44973.593479 ehk(eV)=-44941.964555 sev(eV)=-127.461515

--- BNDFP:  begin iteration 3 of 12
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.620291
   smooth rhoves     11.107548   charge     3.835349
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.748140 -3.575952 -0.844926
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
     potential shift to crystal energy zero:    0.000503
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -6.127022      -182.314718      -188.441740
   rhoval*ves            -48.030181      -113.790712      -161.820893
   psnuc*ves              70.245276    -12969.122178    -12898.876902
   utot                   11.107548     -6541.456445     -6530.348897
   rho*exc                -2.748140      -126.718973      -129.467113
   rho*vxc                -3.575952      -167.796505      -171.372457
   valence chg             3.835349         7.164651        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of    60 k=  0.0625  0.0625  0.0625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     2 of    60 k= -0.0625  0.1875  0.1875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     3 of    60 k= -0.1875  0.3125  0.3125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     4 of    60 k= -0.3125  0.4375  0.4375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     5 of    60 k= -0.4375  0.5625  0.5625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     6 of    60 k= -0.5625  0.6875  0.6875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     7 of    60 k= -0.6875  0.8125  0.8125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     8 of    60 k= -0.8125  0.9375  0.9375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     9 of    60 k=  0.0625  0.0625  0.3125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    10 of    60 k= -0.0625  0.1875  0.4375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    11 of    60 k= -0.1875  0.3125  0.5625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    12 of    60 k= -0.3125  0.4375  0.6875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    13 of    60 k= -0.4375  0.5625  0.8125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    14 of    60 k= -0.5625  0.6875  0.9375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    15 of    60 k= -0.6875  0.8125  1.0625 ndimh = nmto+napw =    18   18    0
 ... Done MPI k-loop: elapsed time=   0.0534

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.148169;  11.000000 electrons
         Sum occ. bands:   -5.251649, incl. Bloechl correction: -0.012268
 bndfp:Generating TDOS: efermi= -0.148169  dos window emin emax=  -0.774071  2.791753


 mkrout:  Qtrue      sm,loc       local
   1   10.297519    2.200427    8.097092
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -5.251649 Vin*nin=    -188.441740 Ek=Eb-Vin*nin=     183.190091
 Ek(core)=    3171.756639 Exc=    -129.467113 Ees=   -6530.348897 Eharris=   -3304.869280

 mkekin:
   nout*Vin = smpart,onsite,total=:     -4.778849   -205.317019   -210.095869
    E_B(band energy sum)=   -5.251649  E_B-nout*Vin=  204.844220

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.463277
   smooth rhoves      6.414066   charge     2.902909
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -1.980933 -2.576057 -0.778176
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=F
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -4.013607      -196.297023      -200.310630
   rhoval*ves            -40.418130      -131.103281      -171.521411
   psnuc*ves              53.246261    -12980.824930    -12927.578668
   utot                    6.414066     -6555.964106     -6549.550040
   rho*exc                -1.980933      -129.503698      -131.484631
   rho*vxc                -2.576057      -171.471048      -174.047106
   valence chg             2.902909         8.097092        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:      0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=      204.844220 Ekcore=      3171.756639 Ektot    =     3376.600859
 Exc=    -131.484631 Ees   =     -6549.550040 EKohnSham=    -3304.433811
 wgtsmooth=   3.1622776601683791E-002
 mixrho: sought 8 iter from file mixm ; read 2 RMS DQ= 5.12E-2  last it= 1.23E-1
 AMIX: nmix=2 mmix=8  nelts=  2074  beta=1.00000  tm= 5.00000  rmsdel= 5.12D-02
   tj:-0.76791  -0.08165
 mixrho: add corrections to qcell smrho = -0.31521D-06 -0.40134D-08

 iors  : write rst restart file (binary mesh density)

   it  3  of 12    ehf=   -3304.869280   ehk=   -3304.433811
 From last iter    ehf=   -3305.472187   ehk=   -3303.147522
 diffe(q)=  0.602907 (0.051213)    tol= 0.000010 (0.000010)   more=T
i ehf(eV)=-44965.390444 ehk(eV)=-44959.465546 sev(eV)=-71.452882

--- BNDFP:  begin iteration 4 of 12
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.525743
   smooth rhoves      8.256180   charge     3.313051
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.319647 -3.017521 -0.807475
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
     potential shift to crystal energy zero:    0.000461
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -4.941480      -185.839184      -190.780664
   rhoval*ves            -43.897574      -119.382195      -163.279769
   psnuc*ves              60.409935    -12969.317726    -12908.907791
   utot                    8.256180     -6544.349961     -6536.093780
   rho*exc                -2.319647      -127.961918      -130.281565
   rho*vxc                -3.017521      -169.433789      -172.451310
   valence chg             3.313051         7.686949        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of    60 k=  0.0625  0.0625  0.0625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     2 of    60 k= -0.0625  0.1875  0.1875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     3 of    60 k= -0.1875  0.3125  0.3125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     4 of    60 k= -0.3125  0.4375  0.4375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     5 of    60 k= -0.4375  0.5625  0.5625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     6 of    60 k= -0.5625  0.6875  0.6875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     7 of    60 k= -0.6875  0.8125  0.8125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     8 of    60 k= -0.8125  0.9375  0.9375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     9 of    60 k=  0.0625  0.0625  0.3125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    10 of    60 k= -0.0625  0.1875  0.4375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    11 of    60 k= -0.1875  0.3125  0.5625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    12 of    60 k= -0.3125  0.4375  0.6875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    13 of    60 k= -0.4375  0.5625  0.8125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    14 of    60 k= -0.5625  0.6875  0.9375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    15 of    60 k= -0.6875  0.8125  1.0625 ndimh = nmto+napw =    18   18    0
 ... Done MPI k-loop: elapsed time=   0.0537

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.132346;  11.000000 electrons
         Sum occ. bands:   -0.993219, incl. Bloechl correction: -0.006627
 bndfp:Generating TDOS: efermi=  0.132346  dos window emin emax=  -0.684035  3.072269


 mkrout:  Qtrue      sm,loc       local
   1    9.930892    3.113709    6.817183
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -0.993219 Vin*nin=    -190.780664 Ek=Eb-Vin*nin=     189.787444
 Ek(core)=    3171.756639 Exc=    -130.281565 Ees=   -6536.093780 Eharris=   -3304.831261

 mkekin:
   nout*Vin = smpart,onsite,total=:     -6.042093   -168.378044   -174.420138
    E_B(band energy sum)=   -0.993219  E_B-nout*Vin=  173.426918

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.676372
   smooth rhoves     13.152472   charge     4.182817
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -3.052152 -3.972407 -0.866429
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=F
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -7.029083      -174.837177      -181.866260
   rhoval*ves            -50.154645      -105.981094      -156.135739
   psnuc*ves              76.459589    -12962.467179    -12886.007589
   utot                   13.152472     -6534.224136     -6521.071664
   rho*exc                -3.052152      -125.574326      -128.626478
   rho*vxc                -3.972407      -166.285214      -170.257621
   valence chg             4.182817         6.817183        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:      0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=      173.426918 Ekcore=      3171.756639 Ektot    =     3345.183558
 Exc=    -128.626478 Ees   =     -6521.071664 EKohnSham=    -3304.514584
 wgtsmooth=   3.1622776601683791E-002
 mixrho: sought 8 iter from file mixm ; read 3 RMS DQ= 4.22E-2  last it= 5.12E-2
 AMIX: nmix=3 mmix=8  nelts=  2074  beta=1.00000  tm= 5.00000  rmsdel= 4.22D-02
   tj: 0.74166  -0.16547   0.00023
 mixrho: add corrections to qcell smrho = -0.17530D-06 -0.22320D-08

 iors  : write rst restart file (binary mesh density)

   it  4  of 12    ehf=   -3304.831261   ehk=   -3304.514584
 From last iter    ehf=   -3304.869280   ehk=   -3304.433811
 diffe(q)=  0.038018 (0.042227)    tol= 0.000010 (0.000010)   more=T
i ehf(eV)=-44964.873174 ehk(eV)=-44960.564532 sev(eV)=-13.513546

--- BNDFP:  begin iteration 5 of 12
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.570816
   smooth rhoves      9.541427   charge     3.556294
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.515298 -3.272452 -0.825708
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
     potential shift to crystal energy zero:    0.000477
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -5.474529      -185.001328      -190.475858
   rhoval*ves            -45.954938      -117.383196      -163.338134
   psnuc*ves              65.037792    -12970.249397    -12905.211605
   utot                    9.541427     -6543.816296     -6534.274870
   rho*exc                -2.515298      -127.432934      -129.948232
   rho*vxc                -3.272452      -168.737527      -172.009979
   valence chg             3.556294         7.443706        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of    60 k=  0.0625  0.0625  0.0625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     2 of    60 k= -0.0625  0.1875  0.1875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     3 of    60 k= -0.1875  0.3125  0.3125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     4 of    60 k= -0.3125  0.4375  0.4375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     5 of    60 k= -0.4375  0.5625  0.5625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     6 of    60 k= -0.5625  0.6875  0.6875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     7 of    60 k= -0.6875  0.8125  0.8125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     8 of    60 k= -0.8125  0.9375  0.9375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     9 of    60 k=  0.0625  0.0625  0.3125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    10 of    60 k= -0.0625  0.1875  0.4375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    11 of    60 k= -0.1875  0.3125  0.5625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    12 of    60 k= -0.3125  0.4375  0.6875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    13 of    60 k= -0.4375  0.5625  0.8125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    14 of    60 k= -0.5625  0.6875  0.9375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    15 of    60 k= -0.6875  0.8125  1.0625 ndimh = nmto+napw =    18   18    0
 ... Done MPI k-loop: elapsed time=   0.0535

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.014879;  11.000000 electrons
         Sum occ. bands:   -2.770224, incl. Bloechl correction: -0.009077
 bndfp:Generating TDOS: efermi= -0.014879  dos window emin emax=  -0.725887  2.925044


 mkrout:  Qtrue      sm,loc       local
   1   10.105683    2.685743    7.419940
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.770224 Vin*nin=    -190.475858 Ek=Eb-Vin*nin=     187.705634
 Ek(core)=    3171.756639 Exc=    -129.948232 Ees=   -6534.274870 Eharris=   -3304.760829

 mkekin:
   nout*Vin = smpart,onsite,total=:     -5.495500   -184.335642   -189.831142
    E_B(band energy sum)=   -2.770224  E_B-nout*Vin=  187.060918

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.575486
   smooth rhoves      9.678781   charge     3.580060
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.534251 -3.297145 -0.827512
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=F
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -5.522086      -184.602389      -190.124476
   rhoval*ves            -46.157297      -116.896962      -163.054259
   psnuc*ves              65.514859    -12969.845858    -12904.330999
   utot                    9.678781     -6543.371410     -6533.692629
   rho*exc                -2.534251      -127.351062      -129.885313
   rho*vxc                -3.297145      -168.629433      -171.926578
   valence chg             3.580060         7.419940        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=      187.060918 Ekcore=      3171.756639 Ektot    =     3358.817557
 Exc=    -129.885313 Ees   =     -6533.692629 EKohnSham=    -3304.760384
 wgtsmooth=   3.1622776601683791E-002
 mixrho: sought 8 iter from file mixm ; read 4 RMS DQ= 1.53E-3  last it= 4.22E-2
 AMIX: nmix=4 mmix=8  nelts=  2074  beta=1.00000  tm= 5.00000  rmsdel= 1.53D-03
   tj: 0.42366   0.74179  -0.16550   0.00023
 mixrho: add corrections to qcell smrho = -0.17533D-06 -0.22324D-08

 iors  : write rst restart file (binary mesh density)

   it  5  of 12    ehf=   -3304.760829   ehk=   -3304.760384
 From last iter    ehf=   -3304.831261   ehk=   -3304.514584
 diffe(q)=  0.070433 (0.001529)    tol= 0.000010 (0.000010)   more=T
i ehf(eV)=-44963.914882 ehk(eV)=-44963.908839 sev(eV)=-37.691115

--- BNDFP:  begin iteration 6 of 12
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.570815
   smooth rhoves      9.541421   charge     3.556293
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.515298 -3.272452 -0.825708
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
     potential shift to crystal energy zero:    0.000477
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -5.474530      -185.001335      -190.475864
   rhoval*ves            -45.954929      -117.383208      -163.338137
   psnuc*ves              65.037772    -12970.249413    -12905.211641
   utot                    9.541421     -6543.816310     -6534.274889
   rho*exc                -2.515298      -127.432937      -129.948235
   rho*vxc                -3.272452      -168.737531      -172.009983
   valence chg             3.556293         7.443707        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of    60 k=  0.0625  0.0625  0.0625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     2 of    60 k= -0.0625  0.1875  0.1875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     3 of    60 k= -0.1875  0.3125  0.3125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     4 of    60 k= -0.3125  0.4375  0.4375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     5 of    60 k= -0.4375  0.5625  0.5625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     6 of    60 k= -0.5625  0.6875  0.6875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     7 of    60 k= -0.6875  0.8125  0.8125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     8 of    60 k= -0.8125  0.9375  0.9375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     9 of    60 k=  0.0625  0.0625  0.3125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    10 of    60 k= -0.0625  0.1875  0.4375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    11 of    60 k= -0.1875  0.3125  0.5625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    12 of    60 k= -0.3125  0.4375  0.6875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    13 of    60 k= -0.4375  0.5625  0.8125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    14 of    60 k= -0.5625  0.6875  0.9375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    15 of    60 k= -0.6875  0.8125  1.0625 ndimh = nmto+napw =    18   18    0
 ... Done MPI k-loop: elapsed time=   0.0534

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.014878;  11.000000 electrons
         Sum occ. bands:   -2.770209, incl. Bloechl correction: -0.009077
 bndfp:Generating TDOS: efermi= -0.014878  dos window emin emax=  -0.725888  2.925045


 mkrout:  Qtrue      sm,loc       local
   1   10.108512    2.682701    7.425811
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.770209 Vin*nin=    -190.475864 Ek=Eb-Vin*nin=     187.705656
 Ek(core)=    3171.756639 Exc=    -129.948235 Ees=   -6534.274889 Eharris=   -3304.760829

 mkekin:
   nout*Vin = smpart,onsite,total=:     -5.491139   -184.470980   -189.962118
    E_B(band energy sum)=   -2.770209  E_B-nout*Vin=  187.191909

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.574217
   smooth rhoves      9.644003   charge     3.574189
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.529637 -3.291135 -0.827058
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=F
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -5.510705      -184.677648      -190.188353
   rhoval*ves            -46.105318      -116.998090      -163.103408
   psnuc*ves              65.393325    -12969.910108    -12904.516783
   utot                    9.644003     -6543.454099     -6533.810095
   rho*exc                -2.529637      -127.369377      -129.899013
   rho*vxc                -3.291135      -168.653599      -171.944734
   valence chg             3.574189         7.425811        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:      0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=      187.191909 Ekcore=      3171.756639 Ektot    =     3358.948549
 Exc=    -129.899013 Ees   =     -6533.810095 EKohnSham=    -3304.760560
 wgtsmooth=   3.1622776601683791E-002
 mixrho: sought 8 iter from file mixm ; read 5 RMS DQ= 1.21E-3  last it= 1.53E-3
 AMIX: condition of normal eqns >100000. Reducing nmix to 4
 AMIX: condition of normal eqns >100000. Reducing nmix to 3
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=  2074  beta=1.00000  tm= 5.00000  rmsdel= 1.21D-03
   tj:-3.74145
 mixrho: add corrections to qcell smrho = -0.21809D-06 -0.27768D-08

 iors  : write rst restart file (binary mesh density)

   it  6  of 12    ehf=   -3304.760829   ehk=   -3304.760560
 From last iter    ehf=   -3304.760829   ehk=   -3304.760384
 diffe(q)= -0.000000 (0.001212)    tol= 0.000010 (0.000010)   more=T
i ehf(eV)=-44963.914886 ehk(eV)=-44963.911226 sev(eV)=-37.690905

--- BNDFP:  begin iteration 7 of 12
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.569469
   smooth rhoves      9.514455   charge     3.552221
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.512407 -3.268691 -0.825351
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
     potential shift to crystal energy zero:    0.000478
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -5.468164      -184.958377      -190.426541
   rhoval*ves            -45.909698      -117.376720      -163.286417
   psnuc*ves              64.938608    -12970.150494    -12905.211886
   utot                    9.514455     -6543.763607     -6534.249152
   rho*exc                -2.512407      -127.437916      -129.950323
   rho*vxc                -3.268691      -168.744037      -172.012728
   valence chg             3.552221         7.447779        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of    60 k=  0.0625  0.0625  0.0625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     2 of    60 k= -0.0625  0.1875  0.1875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     3 of    60 k= -0.1875  0.3125  0.3125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     4 of    60 k= -0.3125  0.4375  0.4375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     5 of    60 k= -0.4375  0.5625  0.5625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     6 of    60 k= -0.5625  0.6875  0.6875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     7 of    60 k= -0.6875  0.8125  0.8125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     8 of    60 k= -0.8125  0.9375  0.9375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     9 of    60 k=  0.0625  0.0625  0.3125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    10 of    60 k= -0.0625  0.1875  0.4375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    11 of    60 k= -0.1875  0.3125  0.5625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    12 of    60 k= -0.3125  0.4375  0.6875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    13 of    60 k= -0.4375  0.5625  0.8125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    14 of    60 k= -0.5625  0.6875  0.9375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    15 of    60 k= -0.6875  0.8125  1.0625 ndimh = nmto+napw =    18   18    0
 ... Done MPI k-loop: elapsed time=   0.0563

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.012923;  11.000000 electrons
         Sum occ. bands:   -2.744580, incl. Bloechl correction: -0.009050
 bndfp:Generating TDOS: efermi= -0.012923  dos window emin emax=  -0.724778  2.927000


 mkrout:  Qtrue      sm,loc       local
   1   10.106951    2.686361    7.420590
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.744580 Vin*nin=    -190.426541 Ek=Eb-Vin*nin=     187.681961
 Ek(core)=    3171.756639 Exc=    -129.950323 Ees=   -6534.249152 Eharris=   -3304.760874

 mkekin:
   nout*Vin = smpart,onsite,total=:     -5.491609   -184.340431   -189.832040
    E_B(band energy sum)=   -2.744580  E_B-nout*Vin=  187.087460

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.575100
   smooth rhoves      9.671887   charge     3.579410
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.534005 -3.296829 -0.827416
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=F
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -5.522995      -184.608906      -190.131901
   rhoval*ves            -46.144898      -116.913084      -163.057981
   psnuc*ves              65.488671    -12969.862299    -12904.373628
   utot                    9.671887     -6543.387692     -6533.715805
   rho*exc                -2.534005      -127.354724      -129.888729
   rho*vxc                -3.296829      -168.634273      -171.931102
   valence chg             3.579410         7.420590        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:      0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=      187.087460 Ekcore=      3171.756639 Ektot    =     3358.844100
 Exc=    -129.888729 Ees   =     -6533.715805 EKohnSham=    -3304.760434
 wgtsmooth=   3.1622776601683791E-002
 mixrho: sought 8 iter from file mixm ; read 6 RMS DQ= 1.46E-3  last it= 1.21E-3
 AMIX: condition of normal eqns >100000. Reducing nmix to 5
 AMIX: condition of normal eqns >100000. Reducing nmix to 4
 AMIX: condition of normal eqns >100000. Reducing nmix to 3
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: nmix=2 mmix=8  nelts=  2074  beta=1.00000  tm= 5.00000  rmsdel= 1.46D-03
   tj: 4.56404  -4.53114
 mixrho: add corrections to qcell smrho = -0.25579D-06 -0.32569D-08

 iors  : write rst restart file (binary mesh density)

   it  7  of 12    ehf=   -3304.760874   ehk=   -3304.760434
 From last iter    ehf=   -3304.760829   ehk=   -3304.760560
 diffe(q)= -0.000045 (0.001460)    tol= 0.000010 (0.000010)   more=T
i ehf(eV)=-44963.915500 ehk(eV)=-44963.909510 sev(eV)=-37.342205

--- BNDFP:  begin iteration 8 of 12
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.569322
   smooth rhoves      9.513971   charge     3.552634
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.512992 -3.269456 -0.825339
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
     potential shift to crystal energy zero:    0.000478
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -5.471039      -184.951308      -190.422347
   rhoval*ves            -45.906901      -117.374348      -163.281249
   psnuc*ves              64.934844    -12970.154995    -12905.220151
   utot                    9.513971     -6543.764671     -6534.250700
   rho*exc                -2.512992      -127.438209      -129.951201
   rho*vxc                -3.269456      -168.744431      -172.013888
   valence chg             3.552634         7.447366        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of    60 k=  0.0625  0.0625  0.0625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     2 of    60 k= -0.0625  0.1875  0.1875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     3 of    60 k= -0.1875  0.3125  0.3125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     4 of    60 k= -0.3125  0.4375  0.4375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     5 of    60 k= -0.4375  0.5625  0.5625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     6 of    60 k= -0.5625  0.6875  0.6875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     7 of    60 k= -0.6875  0.8125  0.8125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     8 of    60 k= -0.8125  0.9375  0.9375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     9 of    60 k=  0.0625  0.0625  0.3125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    10 of    60 k= -0.0625  0.1875  0.4375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    11 of    60 k= -0.1875  0.3125  0.5625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    12 of    60 k= -0.3125  0.4375  0.6875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    13 of    60 k= -0.4375  0.5625  0.8125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    14 of    60 k= -0.5625  0.6875  0.9375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    15 of    60 k= -0.6875  0.8125  1.0625 ndimh = nmto+napw =    18   18    0
 ... Done MPI k-loop: elapsed time=   0.0534

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.012436;  11.000000 electrons
         Sum occ. bands:   -2.737978, incl. Bloechl correction: -0.009042
 bndfp:Generating TDOS: efermi= -0.012436  dos window emin emax=  -0.724564  2.927487


 mkrout:  Qtrue      sm,loc       local
   1   10.106420    2.687604    7.418816
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.737978 Vin*nin=    -190.422347 Ek=Eb-Vin*nin=     187.684369
 Ek(core)=    3171.756639 Exc=    -129.951201 Ees=   -6534.250700 Eharris=   -3304.760892

 mkekin:
   nout*Vin = smpart,onsite,total=:     -5.493426   -184.292480   -189.785906
    E_B(band energy sum)=   -2.737978  E_B-nout*Vin=  187.047928

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.575401
   smooth rhoves      9.681379   charge     3.581184
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.535490 -3.298765 -0.827538
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=F
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -5.527168      -184.582274      -190.109442
   rhoval*ves            -46.158337      -116.881180      -163.039518
   psnuc*ves              65.521095    -12969.841484    -12904.320389
   utot                    9.681379     -6543.361332     -6533.679953
   rho*exc                -2.535490      -127.349506      -129.884997
   rho*vxc                -3.298765      -168.627389      -171.926154
   valence chg             3.581184         7.418816        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:      0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=      187.047928 Ekcore=      3171.756639 Ektot    =     3358.804568
 Exc=    -129.884997 Ees   =     -6533.679953 EKohnSham=    -3304.760382
 wgtsmooth=   3.1622776601683791E-002
 mixrho: sought 8 iter from file mixm ; read 7 RMS DQ= 1.55E-3  last it= 1.46E-3
 AMIX: condition of normal eqns >100000. Reducing nmix to 6
 AMIX: condition of normal eqns >100000. Reducing nmix to 5
 AMIX: condition of normal eqns >100000. Reducing nmix to 4
 AMIX: condition of normal eqns >100000. Reducing nmix to 3
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: Reducing nmix to  1: t_j exceeds tm: tj= 13.04540   0.80700
 AMIX: Reducing nmix to  0: t_j exceeds tm: tj= 15.86880
 AMIX: nmix=0 mmix=8  nelts=  2074  beta=1.00000  tm= 5.00000  rmsdel= 1.55D-03
 mixrho: add corrections to qcell smrho = -0.44015D-07 -0.56043D-09

 iors  : write rst restart file (binary mesh density)

   it  8  of 12    ehf=   -3304.760892   ehk=   -3304.760382
 From last iter    ehf=   -3304.760874   ehk=   -3304.760434
 diffe(q)= -0.000018 (0.001555)    tol= 0.000010 (0.000010)   more=T
i ehf(eV)=-44963.915751 ehk(eV)=-44963.908806 sev(eV)=-37.252375

--- BNDFP:  begin iteration 9 of 12
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.575401
   smooth rhoves      9.681379   charge     3.581184
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.535490 -3.298765 -0.827538
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
     potential shift to crystal energy zero:    0.000477
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -5.527168      -184.582274      -190.109442
   rhoval*ves            -46.158337      -116.881180      -163.039518
   psnuc*ves              65.521095    -12969.841484    -12904.320389
   utot                    9.681379     -6543.361332     -6533.679953
   rho*exc                -2.535490      -127.349506      -129.884997
   rho*vxc                -3.298765      -168.627389      -171.926154
   valence chg             3.581184         7.418816        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of    60 k=  0.0625  0.0625  0.0625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     2 of    60 k= -0.0625  0.1875  0.1875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     3 of    60 k= -0.1875  0.3125  0.3125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     4 of    60 k= -0.3125  0.4375  0.4375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     5 of    60 k= -0.4375  0.5625  0.5625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     6 of    60 k= -0.5625  0.6875  0.6875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     7 of    60 k= -0.6875  0.8125  0.8125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     8 of    60 k= -0.8125  0.9375  0.9375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     9 of    60 k=  0.0625  0.0625  0.3125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    10 of    60 k= -0.0625  0.1875  0.4375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    11 of    60 k= -0.1875  0.3125  0.5625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    12 of    60 k= -0.3125  0.4375  0.6875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    13 of    60 k= -0.4375  0.5625  0.8125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    14 of    60 k= -0.5625  0.6875  0.9375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    15 of    60 k= -0.6875  0.8125  1.0625 ndimh = nmto+napw =    18   18    0
 ... Done MPI k-loop: elapsed time=   0.0539

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.034777;  11.000000 electrons
         Sum occ. bands:   -3.063027, incl. Bloechl correction: -0.009520
 bndfp:Generating TDOS: efermi= -0.034777  dos window emin emax=  -0.731137  2.905146


 mkrout:  Qtrue      sm,loc       local
   1   10.135745    2.613886    7.521859
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -3.063027 Vin*nin=    -190.109442 Ek=Eb-Vin*nin=     187.046415
 Ek(core)=    3171.756639 Exc=    -129.884997 Ees=   -6533.679953 Eharris=   -3304.761895

 mkekin:
   nout*Vin = smpart,onsite,total=:     -5.381167   -187.286421   -192.667589
    E_B(band energy sum)=   -3.063027  E_B-nout*Vin=  189.604562

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.558286
   smooth rhoves      9.141783   charge     3.478141
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.449201 -3.186286 -0.820466
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=F
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -5.284195      -186.337251      -191.621446
   rhoval*ves            -45.371656      -118.932242      -164.303898
   psnuc*ves              63.655221    -12971.353052    -12907.697831
   utot                    9.141783     -6545.142647     -6536.000864
   rho*exc                -2.449201      -127.667281      -130.116483
   rho*vxc                -3.186286      -169.046789      -172.233074
   valence chg             3.478141         7.521859        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:      0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=      189.604562 Ekcore=      3171.756639 Ektot    =     3361.361202
 Exc=    -130.116483 Ees   =     -6536.000864 EKohnSham=    -3304.756145
 wgtsmooth=   3.1622776601683791E-002
 mixrho: sought 8 iter from file mixm ; read 8 RMS DQ= 5.99E-3  last it= 1.55E-3
 AMIX: condition of normal eqns >100000. Reducing nmix to 7
 AMIX: condition of normal eqns >100000. Reducing nmix to 6
 AMIX: condition of normal eqns >100000. Reducing nmix to 5
 AMIX: condition of normal eqns >100000. Reducing nmix to 4
 AMIX: condition of normal eqns >100000. Reducing nmix to 3
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=  2074  beta=1.00000  tm= 5.00000  rmsdel= 5.99D-03
   tj: 0.79470
 mixrho: add corrections to qcell smrho = -0.43174D-07 -0.54971D-09

 iors  : write rst restart file (binary mesh density)

   it  9  of 12    ehf=   -3304.761895   ehk=   -3304.756145
 From last iter    ehf=   -3304.760892   ehk=   -3304.760382
 diffe(q)= -0.001003 (0.005985)    tol= 0.000010 (0.000010)   more=T
i ehf(eV)=-44963.929396 ehk(eV)=-44963.851160 sev(eV)=-41.674926

--- BNDFP:  begin iteration 10 of 12
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.571887
   smooth rhoves      9.569334   charge     3.560029
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.517715 -3.275594 -0.826098
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
     potential shift to crystal energy zero:    0.000474
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -5.476913      -184.944863      -190.421777
   rhoval*ves            -45.999357      -117.301784      -163.301141
   psnuc*ves              65.138025    -12970.151813    -12905.013789
   utot                    9.569334     -6543.726799     -6534.157465
   rho*exc                -2.517715      -127.414723      -129.932438
   rho*vxc                -3.275594      -168.713462      -171.989056
   valence chg             3.560029         7.439971        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of    60 k=  0.0625  0.0625  0.0625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     2 of    60 k= -0.0625  0.1875  0.1875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     3 of    60 k= -0.1875  0.3125  0.3125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     4 of    60 k= -0.3125  0.4375  0.4375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     5 of    60 k= -0.4375  0.5625  0.5625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     6 of    60 k= -0.5625  0.6875  0.6875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     7 of    60 k= -0.6875  0.8125  0.8125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     8 of    60 k= -0.8125  0.9375  0.9375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     9 of    60 k=  0.0625  0.0625  0.3125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    10 of    60 k= -0.0625  0.1875  0.4375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    11 of    60 k= -0.1875  0.3125  0.5625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    12 of    60 k= -0.3125  0.4375  0.6875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    13 of    60 k= -0.4375  0.5625  0.8125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    14 of    60 k= -0.5625  0.6875  0.9375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    15 of    60 k= -0.6875  0.8125  1.0625 ndimh = nmto+napw =    18   18    0
 ... Done MPI k-loop: elapsed time=   0.0538

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.020384;  11.000000 electrons
         Sum occ. bands:   -2.849249, incl. Bloechl correction: -0.009197
 bndfp:Generating TDOS: efermi= -0.020384  dos window emin emax=  -0.727388  2.919539


 mkrout:  Qtrue      sm,loc       local
   1   10.116160    2.664206    7.451954
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.849249 Vin*nin=    -190.421777 Ek=Eb-Vin*nin=     187.572528
 Ek(core)=    3171.756639 Exc=    -129.932438 Ees=   -6534.157465 Eharris=   -3304.760736

 mkekin:
   nout*Vin = smpart,onsite,total=:     -5.460979   -185.213400   -190.674380
    E_B(band energy sum)=   -2.849249  E_B-nout*Vin=  187.825131

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.569827
   smooth rhoves      9.505257   charge     3.548046
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.507745 -3.262598 -0.825264
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=F
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -5.449183      -185.110725      -190.559908
   rhoval*ves            -45.906731      -117.506230      -163.412962
   psnuc*ves              64.917246    -12970.275212    -12905.357965
   utot                    9.505257     -6543.890721     -6534.385463
   rho*exc                -2.507745      -127.449231      -129.956976
   rho*vxc                -3.262598      -168.758984      -172.021582
   valence chg             3.548046         7.451954        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:      0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=      187.825131 Ekcore=      3171.756639 Ektot    =     3359.581771
 Exc=    -129.956976 Ees   =     -6534.385463 EKohnSham=    -3304.760669
 wgtsmooth=   3.1622776601683791E-002
 mixrho: sought 8 iter from file mixm ; read 8 RMS DQ= 6.21E-4  last it= 5.99E-3
 AMIX: condition of normal eqns >100000. Reducing nmix to 7
 AMIX: condition of normal eqns >100000. Reducing nmix to 6
 AMIX: condition of normal eqns >100000. Reducing nmix to 5
 AMIX: condition of normal eqns >100000. Reducing nmix to 4
 AMIX: condition of normal eqns >100000. Reducing nmix to 3
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: nmix=2 mmix=8  nelts=  2074  beta=1.00000  tm= 5.00000  rmsdel= 6.21D-04
   tj:-0.01447   0.25005
 mixrho: add corrections to qcell smrho = -0.49584D-07 -0.63133D-09

 iors  : write rst restart file (binary mesh density)

   it 10  of 12    ehf=   -3304.760736   ehk=   -3304.760669
 From last iter    ehf=   -3304.761895   ehk=   -3304.756145
 diffe(q)=  0.001159 (0.000621)    tol= 0.000010 (0.000010)   more=T
i ehf(eV)=-44963.913622 ehk(eV)=-44963.912709 sev(eV)=-38.766309

--- BNDFP:  begin iteration 11 of 12
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.571388
   smooth rhoves      9.554467   charge     3.557343
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.515525 -3.272740 -0.825903
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
     potential shift to crystal energy zero:    0.000474
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -5.471046      -184.960983      -190.432028
   rhoval*ves            -45.977561      -117.329257      -163.306818
   psnuc*ves              65.086496    -12970.151166    -12905.064670
   utot                    9.554467     -6543.740211     -6534.185744
   rho*exc                -2.515525      -127.421139      -129.936665
   rho*vxc                -3.272740      -168.721913      -171.994654
   valence chg             3.557343         7.442657        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of    60 k=  0.0625  0.0625  0.0625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     2 of    60 k= -0.0625  0.1875  0.1875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     3 of    60 k= -0.1875  0.3125  0.3125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     4 of    60 k= -0.3125  0.4375  0.4375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     5 of    60 k= -0.4375  0.5625  0.5625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     6 of    60 k= -0.5625  0.6875  0.6875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     7 of    60 k= -0.6875  0.8125  0.8125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     8 of    60 k= -0.8125  0.9375  0.9375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     9 of    60 k=  0.0625  0.0625  0.3125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    10 of    60 k= -0.0625  0.1875  0.4375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    11 of    60 k= -0.1875  0.3125  0.5625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    12 of    60 k= -0.3125  0.4375  0.6875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    13 of    60 k= -0.4375  0.5625  0.8125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    14 of    60 k= -0.5625  0.6875  0.9375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    15 of    60 k= -0.6875  0.8125  1.0625 ndimh = nmto+napw =    18   18    0
 ... Done MPI k-loop: elapsed time=   0.0534

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018829;  11.000000 electrons
         Sum occ. bands:   -2.826994, incl. Bloechl correction: -0.009164
 bndfp:Generating TDOS: efermi= -0.018829  dos window emin emax=  -0.726910  2.921094


 mkrout:  Qtrue      sm,loc       local
   1   10.113805    2.669437    7.444368
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.826994 Vin*nin=    -190.432028 Ek=Eb-Vin*nin=     187.605034
 Ek(core)=    3171.756639 Exc=    -129.936665 Ees=   -6534.185744 Eharris=   -3304.760735

 mkekin:
   nout*Vin = smpart,onsite,total=:     -5.468248   -185.001362   -190.469610
    E_B(band energy sum)=   -2.826994  E_B-nout*Vin=  187.642616

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.571128
   smooth rhoves      9.545786   charge     3.555632
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.514056 -3.270825 -0.825792
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=F
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -5.466730      -184.987683      -190.454414
   rhoval*ves            -45.965378      -117.360179      -163.325557
   psnuc*ves              65.056950    -12970.171219    -12905.114269
   utot                    9.545786     -6543.765699     -6534.219913
   rho*exc                -2.514056      -127.426020      -129.940076
   rho*vxc                -3.270825      -168.728352      -171.999177
   valence chg             3.555632         7.444368        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:      0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=      187.642616 Ekcore=      3171.756639 Ektot    =     3359.399256
 Exc=    -129.940076 Ees   =     -6534.219913 EKohnSham=    -3304.760733
 wgtsmooth=   3.1622776601683791E-002
 mixrho: sought 8 iter from file mixm ; read 8 RMS DQ= 9.08E-5  last it= 6.21E-4
 AMIX: condition of normal eqns >100000. Reducing nmix to 7
 AMIX: condition of normal eqns >100000. Reducing nmix to 6
 AMIX: condition of normal eqns >100000. Reducing nmix to 5
 AMIX: condition of normal eqns >100000. Reducing nmix to 4
 AMIX: condition of normal eqns >100000. Reducing nmix to 3
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: nmix=2 mmix=8  nelts=  2074  beta=1.00000  tm= 5.00000  rmsdel= 9.08D-05
   tj:-0.08146  -0.00810
 mixrho: add corrections to qcell smrho = -0.45336D-07 -0.57724D-09

 iors  : write rst restart file (binary mesh density)

   it 11  of 12    ehf=   -3304.760735   ehk=   -3304.760733
 From last iter    ehf=   -3304.760736   ehk=   -3304.760669
 diffe(q)=  0.000001 (0.000091)    tol= 0.000010 (0.000010)   more=T
i ehf(eV)=-44963.913605 ehk(eV)=-44963.913586 sev(eV)=-38.463515

--- BNDFP:  begin iteration 12 of 12
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.571338
   smooth rhoves      9.552399   charge     3.556877
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.515097 -3.272182 -0.825878
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
     potential shift to crystal energy zero:    0.000474
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -5.469649      -184.966658      -190.436308
   rhoval*ves            -45.974883      -117.335565      -163.310448
   psnuc*ves              65.079682    -12970.153177    -12905.073496
   utot                    9.552399     -6543.744371     -6534.191972
   rho*exc                -2.515097      -127.422176      -129.937274
   rho*vxc                -3.272182      -168.723279      -171.995461
   valence chg             3.556877         7.443123        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of    60 k=  0.0625  0.0625  0.0625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     2 of    60 k= -0.0625  0.1875  0.1875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     3 of    60 k= -0.1875  0.3125  0.3125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     4 of    60 k= -0.3125  0.4375  0.4375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     5 of    60 k= -0.4375  0.5625  0.5625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     6 of    60 k= -0.5625  0.6875  0.6875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     7 of    60 k= -0.6875  0.8125  0.8125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     8 of    60 k= -0.8125  0.9375  0.9375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt     9 of    60 k=  0.0625  0.0625  0.3125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    10 of    60 k= -0.0625  0.1875  0.4375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    11 of    60 k= -0.1875  0.3125  0.5625 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    12 of    60 k= -0.3125  0.4375  0.6875 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    13 of    60 k= -0.4375  0.5625  0.8125 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    14 of    60 k= -0.5625  0.6875  0.9375 ndimh = nmto+napw =    18   18    0
 bndfp: kpt    15 of    60 k= -0.6875  0.8125  1.0625 ndimh = nmto+napw =    18   18    0
 ... Done MPI k-loop: elapsed time=   0.0547

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018653;  11.000000 electrons
         Sum occ. bands:   -2.824438, incl. Bloechl correction: -0.009160
 bndfp:Generating TDOS: efermi= -0.018653  dos window emin emax=  -0.726872  2.921270


 mkrout:  Qtrue      sm,loc       local
   1   10.113517    2.670106    7.443410
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.824438 Vin*nin=    -190.436308 Ek=Eb-Vin*nin=     187.611870
 Ek(core)=    3171.756639 Exc=    -129.937274 Ees=   -6534.191972 Eharris=   -3304.760736

 mkekin:
   nout*Vin = smpart,onsite,total=:     -5.469239   -184.973444   -190.442683
    E_B(band energy sum)=   -2.824438  E_B-nout*Vin=  187.618245

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.571291
   smooth rhoves      9.550887   charge     3.556590
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.514856 -3.271867 -0.825858
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=F
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -5.468965      -184.970990      -190.439955
   rhoval*ves            -45.972719      -117.340737      -163.313456
   psnuc*ves              65.074493    -12970.156537    -12905.082043
   utot                    9.550887     -6543.748637     -6534.197750
   rho*exc                -2.514856      -127.423015      -129.937871
   rho*vxc                -3.271867      -168.724386      -171.996253
   valence chg             3.556590         7.443410        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:      0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=      187.618245 Ekcore=      3171.756639 Ektot    =     3359.374885
 Exc=    -129.937871 Ees   =     -6534.197750 EKohnSham=    -3304.760736
 wgtsmooth=   3.1622776601683791E-002
 mixrho: sought 8 iter from file mixm ; read 8 RMS DQ= 1.54E-5  last it= 9.08E-5
 AMIX: condition of normal eqns >100000. Reducing nmix to 7
 AMIX: condition of normal eqns >100000. Reducing nmix to 6
 AMIX: condition of normal eqns >100000. Reducing nmix to 5
 AMIX: condition of normal eqns >100000. Reducing nmix to 4
 AMIX: condition of normal eqns >100000. Reducing nmix to 3
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=  2074  beta=1.00000  tm= 5.00000  rmsdel= 1.54D-05
   tj:-0.20380
 mixrho: add corrections to qcell smrho = -0.44968D-07 -0.57256D-09

 iors  : write rst restart file (binary mesh density)

   it 12  of 12    ehf=   -3304.760736   ehk=   -3304.760736
 From last iter    ehf=   -3304.760735   ehk=   -3304.760733
 diffe(q)= -0.000002 (0.000015)    tol= 0.000010 (0.000010)   more=F
x ehf(eV)=-44963.913629 ehk(eV)=-44963.913626 sev(eV)=-38.428734
Exit 0 procid= 0 OK! end of LMF ======================
INFO: Ubuntu 20.04.4 LTS \n \l
INFO: GNU Fortran (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0
INFO: -O2 -g -fimplicit-none -finit-integer=NaN -finit-real=NaN -JOBJ.gfortran -IOBJ.gfortran
INFO: MATH: -lmkl_rt
INFO: git: commit 913e769c0a5a77a2254ce7ce7011c5bc7fb5168a
INFO:    : Author: Takao Kotani <takaokotani@gmail.com>
INFO:    : Date:   Mon Feb 13 19:43:59 2023 +0900
INFO: linked at Tue Feb 14 12:44:43 JST 2023
===START LMF with   ===
 mpisize=           4
 m_lmfinit:program LMF
 m_lmfinit:program LMF
 m_lmfinit:program LMF
 m_lmfinit:program LMF
 cmdl for python=/home/takao/ecalj/SRC/TestInstall/bin/ctrl2ctrlp.py  --no-iactiv cu -vnk=8 -vbigbas=t -vmetal=3 -vrsm2=1.3 -vrsmd1x=1 -vlmx=4 -vpwmode=0 -voveps=0d-7<ctrl.cu >ctrlp.cu
rval2: STRUC_NSPEC             requ n= 1 val= 1.00000000
rval2: STRUC_NBAS              requ n= 1 val= 1.00000000
rval2: HAM_NSPIN               defa n= 1 val= 1.00000000
rval2: IO_VERBOS               defa n= 1 val= 31.00000000
rval2: IO_TIM                  defa n= 1 val= 1.00000000
rval2: STRUC_ALAT              ---- n= 1 val= 6.79800000
rval2: STRUC_DALAT             ---- n= 1 val= 0.00000000
rval2: STRUC_PLAT              requ n= 9 val= 0.00000000  0.50000000  0.50000000  0.50000000  0.00000000  0.50000000  0.50000000  0.50000000  0.00000000
rval2: OPTIONS_HF              defa n= 1 val= 0.00000000
rval2: HAM_REL                 defa n= 1 val= 1.00000000
rval2: HAM_SO                  defa n= 1 val= 0.00000000
rval2: HAM_SOCAXIS             defa n= 3 val= 0.00000000  0.00000000  1.00000000
rval2: HAM_GMAX                defa n= 1 val= 9.00000000
rval2: HAM_FTMESH              defa n= 3 val= 10.00000000  10.00000000  10.00000000
rval2: HAM_TOL                 defa n= 1 val= 0.00000100
rval2: HAM_FRZWF               defa n= 1 val= 0.00000000
rval2: HAM_XCFUN               defa n= 1 val= 2.00000000
mixing parameters: A/B nmix wt: 0 -1 1.000000  1.000000  0.000000 beta elin wc killj=  1.000000 -1.000000 -1
rval2: HAM_FORCES              defa n= 1 val= 0.00000000

 ... Species  1
rval2: HAM_RDSIG               defa n= 1 val= 1.00000000
rval2: HAM_ScaledSigma         defa n= 1 val= 1.00000000
rval2: HAM_EWALD               defa n= 1 val= 0.00000000
rval2: HAM_OVEPS               defa n= 1 val= 0.00000000
rval2: HAM_PWMODE              defa n= 1 val= 0.00000000
rval2: HAM_PWEMAX              defa n= 1 val= 3.00000000
rval2: HAM_READP               defa n= 1 val= 0.00000000
rval2: HAM_V0FIX               defa n= 1 val= 0.00000000
rval2: HAM_PNUFIX              defa n= 1 val= 0.00000000
rval2: SPEC_Z@1                ---- n= 1 val= 29.00000000
rval2: SPEC_R@1                ---- n= 1 val= 2.31127105
rval2: SPEC_R/W@1              ---- n= 0 val= 
rval2: SPEC_R/A@1              ---- n= 0 val= 
rval2: SPEC_A@1                defa n= 1 val= 0.02500000
rval2: SPEC_NR@1               defa n= 1 val= 0.00000000
mixing parameters: A/B nmix wt: 0 -1 1.000000  1.000000  0.000000 beta elin wc killj=  1.000000 -1.000000 -1

 ... Species  1
rval2: SPEC_RSMH@1             ---- n= 3 val= 2.50000000  2.50000000  1.00000000
rval2: SPEC_EH@1               requ n= 3 val= -0.01000000 -0.01000000 -0.01000000
rval2: SPEC_RSMH2@1            ---- n= 5 val= 1.30000000  0.00000000  1.00000000  1.30000000  0.00000000
mixing parameters: A/B nmix wt: 0 -1 1.000000  1.000000  0.000000 beta elin wc killj=  1.000000 -1.000000 -1

 ... Species  1
rval2: SPEC_EH2@1              requ n= 4 val= -1.00000000 -1.00000000 -1.00000000 -0.01000000
rval2: SPEC_LMX@1              defa n= 1 val= 3.00000000
rval2: SPEC_LMXA@1             defa n= 1 val= 4.00000000
rval2: SPEC_LMXL@1             defa n= 1 val= 4.00000000
rval2: SPEC_P@1                ---- n= 4 val= 4.65000000  4.34000000  3.87000000  4.11000000
rval2: SPEC_Q@1                ---- n= 0 val= 
rval2: SPEC_NMCORE@1           defa n= 1 val= 0.00000000
rval2: SPEC_PZ@1               ---- n= 3 val= 5.50000000  5.50000000  4.50000000
rval2: SPEC_LFOCA@1            defa n= 1 val= 1.00000000
rval2: SPEC_KMXA@1             defa n= 1 val= 4.00000000
rval2: SPEC_RSMA@1             defa n= 1 val= 0.92450842
rval2: SPEC_IDMOD@1            ---- n= 5 val= 0.00000000  0.00000000  0.00000000  1.00000000  1.00000000
rval2: SPEC_FRZWF@1            defa n= 1 val= 0.00000000
rval2: SPEC_IDU@1              ---- n= 0 val= 
rval2: SPEC_UH@1               ---- n= 0 val= 
rval2: SPEC_JH@1               ---- n= 0 val= 
rval2: SPEC_C-HQ@1             defa n= 2 val= -1.00000000  0.00000000
rval2: SPEC_EREF1              defa n= 1 val= 0.00000000
rval2: SITE_POS@1              ---- n= 3 val= 0.00000000  0.00000000  0.00000000
rval2: SITE_RELAX@1            defa n= 3 val= 1.00000000  1.00000000  1.00000000
rval2: SITE_AF@1               defa n= 1 val= 0.00000000
rval2: STR_RMAXS               ---- n= 0 val= 
rval2: STR_RMAX                ---- n= 0 val= 
rval2: STR_MXNBR               defa n= 1 val= 0.00000000
rval2: BZ_NKABC                ---- n= 1 val= 8.00000000
rval2: BZ_BZJOB                ---- n= 1 val= 1.00000000
rval2: BZ_METAL                defa n= 1 val= 3.00000000
rval2: BZ_TETRA                defa n= 1 val= 1.00000000
rval2: BZ_N                    defa n= 1 val= 0.00000000
rval2: BZ_W                    defa n= 1 val= 0.00200000
rval2: BZ_ZBAK                 defa n= 1 val= 0.00000000
rval2: BZ_SAVDOS               defa n= 1 val= 1.00000000
rval2: BZ_NPTS                 defa n= 1 val= 1001.00000000
rval2: BZ_DOSMAX               defa n= 1 val= 2.93992268
rval2: BZ_EFMAX                defa n= 1 val= 5.00000000
rval2: BZ_NEVMX                defa n= 1 val= 0.00000000
rval2: BZ_FSMOM                defa n= 1 val= -99999.00000000
rval2: BZ_FSMOMMETHOD          defa n= 1 val= 0.00000000
rval2: EWALD_AS                defa n= 1 val= 2.00000000
rval2: EWALD_TOL               defa n= 1 val= 0.00000000
rval2: EWALD_NKDMX             defa n= 1 val= 600.00000000
rval2: ITER_NIT                defa n= 1 val= 12.00000000
rval2: ITER_NRMIX              defa n= 1 val= 80.00000000
rval2: ITER_CONV               defa n= 1 val= 0.00001000
rval2: ITER_CONVC              defa n= 1 val= 0.00001000
rval2: ITER_UMIX               defa n= 1 val= 0.50000000
rval2: ITER_TOLU               defa n= 1 val= 0.00000000
rval2: DYN_MODE                defa n= 1 val= 0.00000000
rval2: DYN_NIT                 defa n= 1 val= 1.00000000
rval2: DYN_HESS                defa n= 1 val= 1.00000000
rval2: DYN_XTOL                defa n= 1 val= 0.00100000
rval2: DYN_GTOL                defa n= 1 val= 0.00000000
rval2: DYN_STEP                defa n= 1 val= 0.01500000
rval2: DYN_NKILL               defa n= 1 val= 0.00000000
mixing parameters: A/B nmix wt: 0 -1 1.000000  1.000000  0.000000 beta elin wc killj=  1.000000 -1.000000 -1
 ===> for --jobgw, pwmode is switched to be  0

 ... Species  1
  bndfp (warning): no sigm file found ... LDA calculation only
pnuall: j isp pnu= 1 1 4.650000  4.340000  3.870000  4.110000  5.100000
pnzall: j isp  pz= 1 1 5.500000  5.500000  4.500000  0.000000  0.000000


mmm === MTO setting ===
mmm ispec lmxb lpz nkapii nkaphh=    1    3    1    2    2
mmm rsmh1    1  2.50  2.50  1.00
mmm   eh1    1 -0.01 -0.01 -0.01
mmm rsmh2    1  1.30  0.00  1.00  1.30
 imx=           3           3           3
 imx=           3           3           3
mmm  eh2     1 -1.00 -1.00 -1.00 -0.01
 imx=           3           3           3
mmm pz       1  5.50  5.50  4.50
mmm lh       2  3  2
 imx=           4           4           4
 imx=           4           4           4

                Plat                                  Qlat
 imx=           4           4           4
   0.000000   0.500000   0.500000       -1.000000   1.000000   1.000000
   0.500000   0.000000   0.500000        1.000000  -1.000000   1.000000
   0.500000   0.500000   0.000000        1.000000   1.000000  -1.000000
  Cell vol=   78.538660
 imx=           3           3           3
 imx=           4           4           4

 LATTC:  as= 2.000   tol= 1.00E-12   alat= 6.79800   awald= 0.467
         r1=  2.252   nkd= 201      q1=  6.871   nkg= 331
SGROUP:  1 symmetry operations from 0
 SYMLAT: Bravais system is cubic       with 48 symmetry operations.
 SYMCRY: crystal invariant under  48 symmetry operations for tol=  0.000100
 ig  group op
   1  i*i
   2  i
   3  r3(1,1,-1)
   4  i*r3(1,1,-1)
   5  r3(-1,-1,1)
   6  i*r3(-1,-1,1)
   7  r3d
   8  i*r3d
   9  r3(-1,-1,-1)
  10  i*r3(-1,-1,-1)
  11  r2x
  12  mx
  13  r4x
  14  i*r4x
  15  r4(-1,0,0)
  16  i*r4(-1,0,0)
  17  r3(1,-1,-1)
  18  i*r3(1,-1,-1)
  19  r3(-1,1,1)
  20  i*r3(-1,1,1)
  21  r2(1,1,0)
  22  m(1,1,0)
  23  r2(1,0,-1)
  24  m(1,0,-1)
  25  r2y
  26  my
  27  r4y
  28  i*r4y
  29  r4(0,-1,0)
  30  i*r4(0,-1,0)
  31  r2(0,1,-1)
  32  m(0,1,-1)
  33  r2z
  34  mz
  35  r4(0,0,-1)
  36  i*r4(0,0,-1)
  37  r4z
  38  i*r4z
  39  r3(-1,1,-1)
  40  i*r3(-1,1,-1)
  41  r3(1,-1,1)
  42  i*r3(1,-1,1)
  43  r2(1,0,1)
  44  m(1,0,1)
  45  r2(1,-1,0)
  46  m(1,-1,0)
  47  r2(0,1,1)
  48  m(0,1,1)
 nnnnnn         729         889
 nnnnnn         729         889
GROUPG: the following are sufficient to generate the space group:
 Generators:trans(cart)  = i*r3(1,1,-1) r4x
 Generators::trans(frac) = i*r3(1,1,-1) r4x
MKSYM: found  48  space group operations
SPLCLS: ibas iclass ispec label(ispec)
 SPLCLS     1    1    1     A
 BZMESH:     60 irreducible QP from    8   8   8 shift=TTT
 TETIRR: sorting     3072 tetrahedra ...
 nnnnnn         729         889
     264 inequivalent tetrahedron=
 nnnnnn         729         889
MSHSIZ: mesh has 10 x 10 x 10 divisions; length =     0.481     0.481     0.481
      generated from gmax (a.u.)=      9.0000: 889 vectors of 1000 (88%)
 SGVSYM: 38 symmetry stars found for 861 reciprocal lattice vectors

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol= 1.00E-06
 spec      l    rsm    eh     gmax    last term    cutoff
  A        0*   2.50  -0.01   2.974    2.30E-05      27 
  A        1    2.50  -0.01   3.093    1.29E-06      51 
  A        2    1.00  -0.01   8.508    1.16E-06     813 
  A        0    1.30  -1.00   5.718    2.28E-06     259 
  A        2    1.00  -1.00   8.508    1.16E-06     813 
  A        3    1.30  -0.01   6.806    2.09E-06     411 
 m_qplistinit:start
lmv7: Read rst version ID=  2.00

 iors  : read rst restart file (binary mesh density)
          use from  restart file:use window, pnu,
          ignore in restart file:
         site   1:A       :file pnu is  4.66  4.39  3.88  4.11  5.10
         site   1:A       :file pz  is  5.50  5.50  4.50  0.00  0.00
         site   1, species A       : augmentation lmax changed from 3 to 4
         site   1, species A       : inflate local density from nlm= 16 to 25

 Basis, after reading restart file
 site spec        pos (Cartesian coordinates)         pos (multiples of plat)
   1         A  0.000000   0.000000   0.000000    0.000000   0.000000   0.000000

--- BNDFP:  begin iteration 1 of 12
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.571324
   smooth rhoves      9.551927   charge     3.556785
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.515019 -3.272080 -0.825872
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
     potential shift to crystal energy zero:    0.000474
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -5.469421      -184.967587      -190.437008
   rhoval*ves            -45.974215      -117.336774      -163.310989
   psnuc*ves              65.078069    -12970.153544    -12905.075476
   utot                    9.551927     -6543.745159     -6534.193233
   rho*exc                -2.515019      -127.422403      -129.937422
   rho*vxc                -3.272080      -168.723578      -171.995658
   valence chg             3.556785         7.443215        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of    60 k=  0.0625  0.0625  0.0625 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     2 of    60 k= -0.0625  0.1875  0.1875 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     3 of    60 k= -0.1875  0.3125  0.3125 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     4 of    60 k= -0.3125  0.4375  0.4375 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     5 of    60 k= -0.4375  0.5625  0.5625 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     6 of    60 k= -0.5625  0.6875  0.6875 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     7 of    60 k= -0.6875  0.8125  0.8125 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     8 of    60 k= -0.8125  0.9375  0.9375 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     9 of    60 k=  0.0625  0.0625  0.3125 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    10 of    60 k= -0.0625  0.1875  0.4375 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    11 of    60 k= -0.1875  0.3125  0.5625 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    12 of    60 k= -0.3125  0.4375  0.6875 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    13 of    60 k= -0.4375  0.5625  0.8125 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    14 of    60 k= -0.5625  0.6875  0.9375 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    15 of    60 k= -0.6875  0.8125  1.0625 ndimh = nmto+napw =    31   31    0
 ... Done MPI k-loop: elapsed time=   0.1900

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018883;  11.000000 electrons
         Sum occ. bands:   -2.826235, incl. Bloechl correction: -0.009153
 bndfp:Generating TDOS: efermi= -0.018883  dos window emin emax=  -0.727728  2.921040


 m_bandcal_2nd: to fill eigenfunctions**2 up to Efermi

 mkrout:  Qtrue      sm,loc       local
   1   10.128173    2.833459    7.294714
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.826235 Vin*nin=    -190.437008 Ek=Eb-Vin*nin=     187.610773
 Ek(core)=    3171.756639 Exc=    -129.937422 Ees=   -6534.193233 Eharris=   -3304.763242

 mkekin:
   nout*Vin = smpart,onsite,total=:     -5.930174   -184.355225   -190.285399
    E_B(band energy sum)=   -2.826235  E_B-nout*Vin=  187.459164

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.581329
   smooth rhoves     10.109465   charge     3.705286
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.657898 -3.458572 -0.832968
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=F
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -6.010278      -184.303254      -190.313532
   rhoval*ves            -46.616840      -116.580931      -163.197771
   psnuc*ves              66.835769    -12971.741116    -12904.905347
   utot                   10.109465     -6544.161023     -6534.051559
   rho*exc                -2.657898      -127.269517      -129.927415
   rho*vxc                -3.458572      -168.523784      -171.982356
   valence chg             3.705286         7.294714        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:      0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=      187.459164 Ekcore=      3171.756639 Ektot    =     3359.215803
 Exc=    -129.927415 Ees   =     -6534.051559 EKohnSham=    -3304.763170
 mixrealsmooth= T
 wgtsmooth=   3.1622776601683791E-002
 mixrho: sought 8 iter from file mixm ; read 0 RMS DQ= 4.32E-3
 AMIX: nmix=0 mmix=8  nelts=  2236  beta=1.00000  tm= 5.00000  rmsdel= 4.32D-03
 mixrho: add corrections to qcell smrho = -0.28235D-07 -0.35950D-09

 iors  : write rst restart file (binary mesh density)

   it  1  of 12    ehf=   -3304.763242   ehk=   -3304.763170
i ehf(eV)=-44963.947720 ehk(eV)=-44963.946745 sev(eV)=-38.453188

--- BNDFP:  begin iteration 2 of 12
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.581329
   smooth rhoves     10.109465   charge     3.705286
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.657898 -3.458572 -0.832968
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
     potential shift to crystal energy zero:    0.000042
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -6.010278      -184.303254      -190.313532
   rhoval*ves            -46.616840      -116.580931      -163.197771
   psnuc*ves              66.835769    -12971.741116    -12904.905347
   utot                   10.109465     -6544.161023     -6534.051559
   rho*exc                -2.657898      -127.269517      -129.927415
   rho*vxc                -3.458572      -168.523784      -171.982356
   valence chg             3.705286         7.294714        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of    60 k=  0.0625  0.0625  0.0625 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     2 of    60 k= -0.0625  0.1875  0.1875 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     3 of    60 k= -0.1875  0.3125  0.3125 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     4 of    60 k= -0.3125  0.4375  0.4375 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     5 of    60 k= -0.4375  0.5625  0.5625 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     6 of    60 k= -0.5625  0.6875  0.6875 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     7 of    60 k= -0.6875  0.8125  0.8125 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     8 of    60 k= -0.8125  0.9375  0.9375 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     9 of    60 k=  0.0625  0.0625  0.3125 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    10 of    60 k= -0.0625  0.1875  0.4375 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    11 of    60 k= -0.1875  0.3125  0.5625 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    12 of    60 k= -0.3125  0.4375  0.6875 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    13 of    60 k= -0.4375  0.5625  0.8125 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    14 of    60 k= -0.5625  0.6875  0.9375 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    15 of    60 k= -0.6875  0.8125  1.0625 ndimh = nmto+napw =    31   31    0
 ... Done MPI k-loop: elapsed time=   0.1866

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.020802;  11.000000 electrons
         Sum occ. bands:   -2.854375, incl. Bloechl correction: -0.009202
 bndfp:Generating TDOS: efermi= -0.020802  dos window emin emax=  -0.727847  2.919121


 m_bandcal_2nd: to fill eigenfunctions**2 up to Efermi

 mkrout:  Qtrue      sm,loc       local
   1   10.131649    2.824829    7.306820
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.854375 Vin*nin=    -190.313532 Ek=Eb-Vin*nin=     187.459157
 Ek(core)=    3171.756639 Exc=    -129.927415 Ees=   -6534.051559 Eharris=   -3304.763177

 mkekin:
   nout*Vin = smpart,onsite,total=:     -5.992931   -184.685653   -190.678585
    E_B(band energy sum)=   -2.854375  E_B-nout*Vin=  187.824209

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.579288
   smooth rhoves     10.043955   charge     3.693180
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.647653 -3.445216 -0.832150
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=F
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -5.980099      -184.560698      -190.540798
   rhoval*ves            -46.528529      -116.862883      -163.391412
   psnuc*ves              66.616438    -12971.995011    -12905.378572
   utot                   10.043955     -6544.428947     -6534.384992
   rho*exc                -2.647653      -127.311275      -129.958928
   rho*vxc                -3.445216      -168.578929      -172.024145
   valence chg             3.693180         7.306820        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:      0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=      187.824209 Ekcore=      3171.756639 Ektot    =     3359.580849
 Exc=    -129.958928 Ees   =     -6534.384992 EKohnSham=    -3304.763072
 wgtsmooth=   3.1622776601683791E-002
 mixrho: sought 8 iter from file mixm ; read 1 RMS DQ= 8.01E-4  last it= 4.32E-3
 AMIX: nmix=1 mmix=8  nelts=  2236  beta=1.00000  tm= 5.00000  rmsdel= 8.01D-04
   tj: 0.08296
 mixrho: add corrections to qcell smrho = -0.26733D-07 -0.34038D-09

 iors  : write rst restart file (binary mesh density)

   it  2  of 12    ehf=   -3304.763177   ehk=   -3304.763072
 From last iter    ehf=   -3304.763242   ehk=   -3304.763170
 diffe(q)=  0.000065 (0.000801)    tol= 0.000010 (0.000010)   more=T
i ehf(eV)=-44963.946839 ehk(eV)=-44963.945406 sev(eV)=-38.836061

--- BNDFP:  begin iteration 3 of 12
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.579458
   smooth rhoves     10.049381   charge     3.694185
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.648503 -3.446324 -0.832218
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
     potential shift to crystal energy zero:    0.000043
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -5.982600      -184.539360      -190.521960
   rhoval*ves            -46.535871      -116.839494      -163.375366
   psnuc*ves              66.634633    -12971.973948    -12905.339314
   utot                   10.049381     -6544.406721     -6534.357340
   rho*exc                -2.648503      -127.307811      -129.956313
   rho*vxc                -3.446324      -168.574354      -172.020678
   valence chg             3.694185         7.305815        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of    60 k=  0.0625  0.0625  0.0625 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     2 of    60 k= -0.0625  0.1875  0.1875 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     3 of    60 k= -0.1875  0.3125  0.3125 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     4 of    60 k= -0.3125  0.4375  0.4375 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     5 of    60 k= -0.4375  0.5625  0.5625 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     6 of    60 k= -0.5625  0.6875  0.6875 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     7 of    60 k= -0.6875  0.8125  0.8125 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     8 of    60 k= -0.8125  0.9375  0.9375 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     9 of    60 k=  0.0625  0.0625  0.3125 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    10 of    60 k= -0.0625  0.1875  0.4375 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    11 of    60 k= -0.1875  0.3125  0.5625 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    12 of    60 k= -0.3125  0.4375  0.6875 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    13 of    60 k= -0.4375  0.5625  0.8125 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    14 of    60 k= -0.5625  0.6875  0.9375 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    15 of    60 k= -0.6875  0.8125  1.0625 ndimh = nmto+napw =    31   31    0
 ... Done MPI k-loop: elapsed time=   0.1792

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.011989;  11.000000 electrons
         Sum occ. bands:   -2.728375, incl. Bloechl correction: -0.009007
 bndfp:Generating TDOS: efermi= -0.011989  dos window emin emax=  -0.725680  2.927933


 m_bandcal_2nd: to fill eigenfunctions**2 up to Efermi

 mkrout:  Qtrue      sm,loc       local
   1   10.118671    2.849894    7.268776
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.728375 Vin*nin=    -190.521960 Ek=Eb-Vin*nin=     187.793585
 Ek(core)=    3171.756639 Exc=    -129.956313 Ees=   -6534.357340 Eharris=   -3304.763429

 mkekin:
   nout*Vin = smpart,onsite,total=:     -6.025972   -183.426093   -189.452065
    E_B(band energy sum)=   -2.728375  E_B-nout*Vin=  186.723690

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.586225
   smooth rhoves     10.259462   charge     3.731224
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.679136 -3.486250 -0.834837
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=F
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -6.068332      -183.815807      -189.884140
   rhoval*ves            -46.821311      -116.019404      -162.840715
   psnuc*ves              67.340236    -12971.265450    -12903.925214
   utot                   10.259462     -6543.642427     -6533.382964
   rho*exc                -2.679136      -127.180654      -129.859790
   rho*vxc                -3.486250      -168.406450      -171.892700
   valence chg             3.731224         7.268776        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:      0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=      186.723690 Ekcore=      3171.756639 Ektot    =     3358.480329
 Exc=    -129.859790 Ees   =     -6533.382964 EKohnSham=    -3304.762425
 wgtsmooth=   3.1622776601683791E-002
 mixrho: sought 8 iter from file mixm ; read 2 RMS DQ= 2.36E-3  last it= 8.01E-4
 AMIX: nmix=2 mmix=8  nelts=  2236  beta=1.00000  tm= 5.00000  rmsdel= 2.36D-03
   tj: 0.74720  -0.00311
 mixrho: add corrections to qcell smrho = -0.27292D-07 -0.34750D-09

 iors  : write rst restart file (binary mesh density)

   it  3  of 12    ehf=   -3304.763429   ehk=   -3304.762425
 From last iter    ehf=   -3304.763177   ehk=   -3304.763072
 diffe(q)= -0.000251 (0.002357)    tol= 0.000010 (0.000010)   more=T
i ehf(eV)=-44963.950256 ehk(eV)=-44963.936608 sev(eV)=-37.121727

--- BNDFP:  begin iteration 4 of 12
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.581057
   smooth rhoves     10.098682   charge     3.702879
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.655668 -3.455663 -0.832837
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
     potential shift to crystal energy zero:    0.000043
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -6.002540      -184.371330      -190.373869
   rhoval*ves            -46.603620      -116.647902      -163.251522
   psnuc*ves              66.800984    -12971.809097    -12905.008113
   utot                   10.098682     -6544.228499     -6534.129817
   rho*exc                -2.655668      -127.277971      -129.933639
   rho*vxc                -3.455663      -168.534952      -171.990614
   valence chg             3.702879         7.297121        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of    60 k=  0.0625  0.0625  0.0625 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     2 of    60 k= -0.0625  0.1875  0.1875 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     3 of    60 k= -0.1875  0.3125  0.3125 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     4 of    60 k= -0.3125  0.4375  0.4375 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     5 of    60 k= -0.4375  0.5625  0.5625 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     6 of    60 k= -0.5625  0.6875  0.6875 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     7 of    60 k= -0.6875  0.8125  0.8125 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     8 of    60 k= -0.8125  0.9375  0.9375 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     9 of    60 k=  0.0625  0.0625  0.3125 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    10 of    60 k= -0.0625  0.1875  0.4375 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    11 of    60 k= -0.1875  0.3125  0.5625 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    12 of    60 k= -0.3125  0.4375  0.6875 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    13 of    60 k= -0.4375  0.5625  0.8125 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    14 of    60 k= -0.5625  0.6875  0.9375 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    15 of    60 k= -0.6875  0.8125  1.0625 ndimh = nmto+napw =    31   31    0
 ... Done MPI k-loop: elapsed time=   0.1831

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.019147;  11.000000 electrons
         Sum occ. bands:   -2.830219, incl. Bloechl correction: -0.009164
 bndfp:Generating TDOS: efermi= -0.019147  dos window emin emax=  -0.727515  2.920776


 m_bandcal_2nd: to fill eigenfunctions**2 up to Efermi

 mkrout:  Qtrue      sm,loc       local
   1   10.129011    2.830104    7.298907
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.830219 Vin*nin=    -190.373869 Ek=Eb-Vin*nin=     187.543651
 Ek(core)=    3171.756639 Exc=    -129.933639 Ees=   -6534.129817 Eharris=   -3304.763166

 mkekin:
   nout*Vin = smpart,onsite,total=:     -6.000745   -184.420037   -190.420782
    E_B(band energy sum)=   -2.830219  E_B-nout*Vin=  187.590564

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.580714
   smooth rhoves     10.088310   charge     3.701093
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.654215 -3.453769 -0.832707
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=F
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -5.998622      -184.401380      -190.400002
   rhoval*ves            -46.589274      -116.683572      -163.272846
   psnuc*ves              66.765894    -12971.837529    -12905.071635
   utot                   10.088310     -6544.260550     -6534.172240
   rho*exc                -2.654215      -127.283912      -129.938126
   rho*vxc                -3.453769      -168.542794      -171.996563
   valence chg             3.701093         7.298907        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:      0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=      187.590564 Ekcore=      3171.756639 Ektot    =     3359.347203
 Exc=    -129.938126 Ees   =     -6534.172240 EKohnSham=    -3304.763163
 wgtsmooth=   3.1622776601683791E-002
 mixrho: sought 8 iter from file mixm ; read 3 RMS DQ= 1.06E-4  last it= 2.36E-3
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: nmix=2 mmix=8  nelts=  2236  beta=1.00000  tm= 5.00000  rmsdel= 1.06D-04
   tj: 0.08940   0.16482
 mixrho: add corrections to qcell smrho = -0.24059D-07 -0.30633D-09

 iors  : write rst restart file (binary mesh density)

   it  4  of 12    ehf=   -3304.763166   ehk=   -3304.763163
 From last iter    ehf=   -3304.763429   ehk=   -3304.762425
 diffe(q)=  0.000263 (0.000106)    tol= 0.000010 (0.000010)   more=T
i ehf(eV)=-44963.946681 ehk(eV)=-44963.946648 sev(eV)=-38.507388

--- BNDFP:  begin iteration 5 of 12
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.580972
   smooth rhoves     10.096230   charge     3.702483
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.655358 -3.455259 -0.832807
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
     potential shift to crystal energy zero:    0.000043
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -6.001787      -184.375434      -190.377221
   rhoval*ves            -46.600150      -116.653745      -163.253896
   psnuc*ves              66.792610    -12971.812338    -12905.019728
   utot                   10.096230     -6544.233042     -6534.136812
   rho*exc                -2.655358      -127.279187      -129.934545
   rho*vxc                -3.455259      -168.536557      -171.991816
   valence chg             3.702483         7.297517        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of    60 k=  0.0625  0.0625  0.0625 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     2 of    60 k= -0.0625  0.1875  0.1875 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     3 of    60 k= -0.1875  0.3125  0.3125 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     4 of    60 k= -0.3125  0.4375  0.4375 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     5 of    60 k= -0.4375  0.5625  0.5625 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     6 of    60 k= -0.5625  0.6875  0.6875 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     7 of    60 k= -0.6875  0.8125  0.8125 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     8 of    60 k= -0.8125  0.9375  0.9375 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     9 of    60 k=  0.0625  0.0625  0.3125 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    10 of    60 k= -0.0625  0.1875  0.4375 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    11 of    60 k= -0.1875  0.3125  0.5625 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    12 of    60 k= -0.3125  0.4375  0.6875 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    13 of    60 k= -0.4375  0.5625  0.8125 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    14 of    60 k= -0.5625  0.6875  0.9375 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    15 of    60 k= -0.6875  0.8125  1.0625 ndimh = nmto+napw =    31   31    0
 ... Done MPI k-loop: elapsed time=   0.1776

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018829;  11.000000 electrons
         Sum occ. bands:   -2.825670, incl. Bloechl correction: -0.009157
 bndfp:Generating TDOS: efermi= -0.018829  dos window emin emax=  -0.727421  2.921093


 m_bandcal_2nd: to fill eigenfunctions**2 up to Efermi

 mkrout:  Qtrue      sm,loc       local
   1   10.128558    2.830941    7.297617
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.825670 Vin*nin=    -190.377221 Ek=Eb-Vin*nin=     187.551551
 Ek(core)=    3171.756639 Exc=    -129.934545 Ees=   -6534.136812 Eharris=   -3304.763167

 mkekin:
   nout*Vin = smpart,onsite,total=:     -6.001692   -184.378524   -190.380217
    E_B(band energy sum)=   -2.825670  E_B-nout*Vin=  187.554547

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.580952
   smooth rhoves     10.095643   charge     3.702383
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.655277 -3.455153 -0.832799
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=F
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -6.001572      -184.377386      -190.378958
   rhoval*ves            -46.599339      -116.655998      -163.255337
   psnuc*ves              66.790625    -12971.814350    -12905.023724
   utot                   10.095643     -6544.235174     -6534.139530
   rho*exc                -2.655277      -127.279546      -129.934822
   rho*vxc                -3.455153      -168.537030      -171.992183
   valence chg             3.702383         7.297617        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:      0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=      187.554547 Ekcore=      3171.756639 Ektot    =     3359.311186
 Exc=    -129.934822 Ees   =     -6534.139530 EKohnSham=    -3304.763167
 wgtsmooth=   3.1622776601683791E-002
 mixrho: sought 8 iter from file mixm ; read 4 RMS DQ= 6.58E-6  last it= 1.06E-4
 AMIX: condition of normal eqns >100000. Reducing nmix to 3
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=  2236  beta=1.00000  tm= 5.00000  rmsdel= 6.58D-06
   tj:-0.06603
 mixrho: add corrections to qcell smrho = -0.26040D-07 -0.33156D-09

 iors  : write rst restart file (binary mesh density)

   it  5  of 12    ehf=   -3304.763167   ehk=   -3304.763167
 From last iter    ehf=   -3304.763166   ehk=   -3304.763163
 diffe(q)= -0.000001 (0.000007)    tol= 0.000010 (0.000010)   more=F
c ehf(eV)=-44963.946694 ehk(eV)=-44963.946691 sev(eV)=-38.445497
Exit 0 procid= 0 OK! end of LMF ======================
INFO: Ubuntu 20.04.4 LTS \n \l
INFO: GNU Fortran (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0
INFO: -O2 -g -fimplicit-none -finit-integer=NaN -finit-real=NaN -JOBJ.gfortran -IOBJ.gfortran
INFO: MATH: -lmkl_rt
INFO: git: commit 913e769c0a5a77a2254ce7ce7011c5bc7fb5168a
INFO:    : Author: Takao Kotani <takaokotani@gmail.com>
INFO:    : Date:   Mon Feb 13 19:43:59 2023 +0900
 m_lmfinit:program LMF
 m_lmfinit:program LMF
 m_lmfinit:program LMF
INFO: linked at Tue Feb 14 12:44:43 JST 2023
===START LMF with   ===
 mpisize=           4
 m_lmfinit:program LMF
 cmdl for python=/home/takao/ecalj/SRC/TestInstall/bin/ctrl2ctrlp.py  --no-iactiv cu -vnk=8 -vbigbas=t -vmetal=3 -vrsm2=1.3 -vrsmd1x=1 -vlmx=4 -vpwmode=0 -voveps=0d-7 --band:fn=syml<ctrl.cu >ctrlp.cu
rval2: STRUC_NSPEC             requ n= 1 val= 1.00000000
rval2: STRUC_NBAS              requ n= 1 val= 1.00000000
rval2: HAM_NSPIN               defa n= 1 val= 1.00000000
rval2: IO_VERBOS               defa n= 1 val= 31.00000000
rval2: IO_TIM                  defa n= 1 val= 1.00000000
rval2: STRUC_ALAT              ---- n= 1 val= 6.79800000
rval2: STRUC_DALAT             ---- n= 1 val= 0.00000000
rval2: STRUC_PLAT              requ n= 9 val= 0.00000000  0.50000000  0.50000000  0.50000000  0.00000000  0.50000000  0.50000000  0.50000000  0.00000000
rval2: OPTIONS_HF              defa n= 1 val= 0.00000000
rval2: HAM_REL                 defa n= 1 val= 1.00000000
rval2: HAM_SO                  defa n= 1 val= 0.00000000
rval2: HAM_SOCAXIS             defa n= 3 val= 0.00000000  0.00000000  1.00000000
rval2: HAM_GMAX                defa n= 1 val= 9.00000000
rval2: HAM_FTMESH              defa n= 3 val= 10.00000000  10.00000000  10.00000000
rval2: HAM_TOL                 defa n= 1 val= 0.00000100
mixing parameters: A/B nmix wt: 0 -1 1.000000  1.000000  0.000000 beta elin wc killj=  1.000000 -1.000000 -1
rval2: HAM_FRZWF               defa n= 1 val= 0.00000000
mixing parameters: A/B nmix wt: 0 -1 1.000000  1.000000  0.000000 beta elin wc killj=  1.000000 -1.000000 -1

rval2: HAM_XCFUN               defa n= 1 val= 2.00000000

 ... Species  1
 ... Species  1
rval2: HAM_FORCES              defa n= 1 val= 0.00000000
rval2: HAM_RDSIG               defa n= 1 val= 1.00000000
rval2: HAM_ScaledSigma         defa n= 1 val= 1.00000000
rval2: HAM_EWALD               defa n= 1 val= 0.00000000
rval2: HAM_OVEPS               defa n= 1 val= 0.00000000
rval2: HAM_PWMODE              defa n= 1 val= 0.00000000
rval2: HAM_PWEMAX              defa n= 1 val= 3.00000000
rval2: HAM_READP               defa n= 1 val= 0.00000000
rval2: HAM_V0FIX               defa n= 1 val= 0.00000000
rval2: HAM_PNUFIX              defa n= 1 val= 0.00000000
rval2: SPEC_Z@1                ---- n= 1 val= 29.00000000
rval2: SPEC_R@1                ---- n= 1 val= 2.31127105
rval2: SPEC_R/W@1              ---- n= 0 val= 
rval2: SPEC_R/A@1              ---- n= 0 val= 
rval2: SPEC_A@1                defa n= 1 val= 0.02500000
rval2: SPEC_NR@1               defa n= 1 val= 0.00000000
rval2: SPEC_RSMH@1             ---- n= 3 val= 2.50000000  2.50000000  1.00000000
rval2: SPEC_EH@1               requ n= 3 val= -0.01000000 -0.01000000 -0.01000000
rval2: SPEC_RSMH2@1            ---- n= 5 val= 1.30000000  0.00000000  1.00000000  1.30000000  0.00000000
rval2: SPEC_EH2@1              requ n= 4 val= -1.00000000 -1.00000000 -1.00000000 -0.01000000
rval2: SPEC_LMX@1              defa n= 1 val= 3.00000000
rval2: SPEC_LMXA@1             defa n= 1 val= 4.00000000
mixing parameters: A/B nmix wt: 0 -1 1.000000  1.000000  0.000000 beta elin wc killj=  1.000000 -1.000000 -1
rval2: SPEC_LMXL@1             defa n= 1 val= 4.00000000

 ... Species  1
rval2: SPEC_P@1                ---- n= 4 val= 4.65000000  4.34000000  3.87000000  4.11000000
rval2: SPEC_Q@1                ---- n= 0 val= 
rval2: SPEC_NMCORE@1           defa n= 1 val= 0.00000000
rval2: SPEC_PZ@1               ---- n= 3 val= 5.50000000  5.50000000  4.50000000
rval2: SPEC_LFOCA@1            defa n= 1 val= 1.00000000
rval2: SPEC_KMXA@1             defa n= 1 val= 4.00000000
rval2: SPEC_RSMA@1             defa n= 1 val= 0.92450842
rval2: SPEC_IDMOD@1            ---- n= 5 val= 0.00000000  0.00000000  0.00000000  1.00000000  1.00000000
rval2: SPEC_FRZWF@1            defa n= 1 val= 0.00000000
rval2: SPEC_IDU@1              ---- n= 0 val= 
rval2: SPEC_UH@1               ---- n= 0 val= 
rval2: SPEC_JH@1               ---- n= 0 val= 
rval2: SPEC_C-HQ@1             defa n= 2 val= -1.00000000  0.00000000
rval2: SPEC_EREF1              defa n= 1 val= 0.00000000
rval2: SITE_POS@1              ---- n= 3 val= 0.00000000  0.00000000  0.00000000
rval2: SITE_RELAX@1            defa n= 3 val= 1.00000000  1.00000000  1.00000000
rval2: SITE_AF@1               defa n= 1 val= 0.00000000
rval2: STR_RMAXS               ---- n= 0 val= 
rval2: STR_RMAX                ---- n= 0 val= 
rval2: STR_MXNBR               defa n= 1 val= 0.00000000
rval2: BZ_NKABC                ---- n= 1 val= 8.00000000
rval2: BZ_BZJOB                ---- n= 1 val= 1.00000000
rval2: BZ_METAL                defa n= 1 val= 3.00000000
rval2: BZ_TETRA                defa n= 1 val= 1.00000000
rval2: BZ_N                    defa n= 1 val= 0.00000000
rval2: BZ_W                    defa n= 1 val= 0.00200000
rval2: BZ_ZBAK                 defa n= 1 val= 0.00000000
rval2: BZ_SAVDOS               defa n= 1 val= 1.00000000
rval2: BZ_NPTS                 defa n= 1 val= 1001.00000000
rval2: BZ_DOSMAX               defa n= 1 val= 2.93992268
rval2: BZ_EFMAX                defa n= 1 val= 5.00000000
rval2: BZ_NEVMX                defa n= 1 val= 0.00000000
rval2: BZ_FSMOM                defa n= 1 val= -99999.00000000
rval2: BZ_FSMOMMETHOD          defa n= 1 val= 0.00000000
rval2: EWALD_AS                defa n= 1 val= 2.00000000
rval2: EWALD_TOL               defa n= 1 val= 0.00000000
rval2: EWALD_NKDMX             defa n= 1 val= 600.00000000
rval2: ITER_NIT                defa n= 1 val= 12.00000000
rval2: ITER_NRMIX              defa n= 1 val= 80.00000000
rval2: ITER_CONV               defa n= 1 val= 0.00001000
rval2: ITER_CONVC              defa n= 1 val= 0.00001000
rval2: ITER_UMIX               defa n= 1 val= 0.50000000
rval2: ITER_TOLU               defa n= 1 val= 0.00000000
rval2: DYN_MODE                defa n= 1 val= 0.00000000
rval2: DYN_NIT                 defa n= 1 val= 1.00000000
rval2: DYN_HESS                defa n= 1 val= 1.00000000
rval2: DYN_XTOL                defa n= 1 val= 0.00100000
rval2: DYN_GTOL                defa n= 1 val= 0.00000000
rval2: DYN_STEP                defa n= 1 val= 0.01500000
rval2: DYN_NKILL               defa n= 1 val= 0.00000000
mixing parameters: A/B nmix wt: 0 -1 1.000000  1.000000  0.000000 beta elin wc killj=  1.000000 -1.000000 -1
 ===> for --jobgw, pwmode is switched to be  0

 ... Species  1
  bndfp (warning): no sigm file found ... LDA calculation only
pnuall: j isp pnu= 1 1 4.650000  4.340000  3.870000  4.110000  5.100000
pnzall: j isp  pz= 1 1 5.500000  5.500000  4.500000  0.000000  0.000000


mmm === MTO setting ===
mmm ispec lmxb lpz nkapii nkaphh=    1    3    1    2    2
mmm rsmh1    1  2.50  2.50  1.00
mmm   eh1    1 -0.01 -0.01 -0.01
 imx=           3           3           3
 imx=           3           3           3
mmm rsmh2    1  1.30  0.00  1.00  1.30
mmm  eh2     1 -1.00 -1.00 -1.00 -0.01
mmm pz       1  5.50  5.50  4.50
mmm lh       2  3  2
 imx=           3           3           3
 imx=           4           4           4
 imx=           4           4           4

                Plat                                  Qlat
 imx=           4           4           4
   0.000000   0.500000   0.500000       -1.000000   1.000000   1.000000
   0.500000   0.000000   0.500000        1.000000  -1.000000   1.000000
   0.500000   0.500000   0.000000        1.000000   1.000000  -1.000000
  Cell vol=   78.538660
 imx=           3           3           3
 imx=           4           4           4

 LATTC:  as= 2.000   tol= 1.00E-12   alat= 6.79800   awald= 0.467
         r1=  2.252   nkd= 201      q1=  6.871   nkg= 331
SGROUP:  1 symmetry operations from 0
 SYMLAT: Bravais system is cubic       with 48 symmetry operations.
 SYMCRY: crystal invariant under  48 symmetry operations for tol=  0.000100
 ig  group op
   1  i*i
   2  i
   3  r3(1,1,-1)
   4  i*r3(1,1,-1)
   5  r3(-1,-1,1)
   6  i*r3(-1,-1,1)
   7  r3d
   8  i*r3d
   9  r3(-1,-1,-1)
  10  i*r3(-1,-1,-1)
  11  r2x
  12  mx
  13  r4x
  14  i*r4x
  15  r4(-1,0,0)
  16  i*r4(-1,0,0)
  17  r3(1,-1,-1)
  18  i*r3(1,-1,-1)
  19  r3(-1,1,1)
  20  i*r3(-1,1,1)
  21  r2(1,1,0)
  22  m(1,1,0)
  23  r2(1,0,-1)
  24  m(1,0,-1)
  25  r2y
  26  my
  27  r4y
  28  i*r4y
  29  r4(0,-1,0)
  30  i*r4(0,-1,0)
  31  r2(0,1,-1)
  32  m(0,1,-1)
  33  r2z
  34  mz
  35  r4(0,0,-1)
  36  i*r4(0,0,-1)
  37  r4z
  38  i*r4z
  39  r3(-1,1,-1)
  40  i*r3(-1,1,-1)
  41  r3(1,-1,1)
  42  i*r3(1,-1,1)
  43  r2(1,0,1)
  44  m(1,0,1)
  45  r2(1,-1,0)
  46  m(1,-1,0)
  47  r2(0,1,1)
  48  m(0,1,1)
 nnnnnn         729         889
 nnnnnn         729         889
GROUPG: the following are sufficient to generate the space group:
 Generators:trans(cart)  = i*r3(1,1,-1) r4x
 Generators::trans(frac) = i*r3(1,1,-1) r4x
MKSYM: found  48  space group operations
SPLCLS: ibas iclass ispec label(ispec)
 SPLCLS     1    1    1     A
 BZMESH:     60 irreducible QP from    8   8   8 shift=TTT
 TETIRR: sorting     3072 tetrahedra ...
     264 inequivalent tetrahedron=
 nnnnnn         729         889
 nnnnnn         729         889
MSHSIZ: mesh has 10 x 10 x 10 divisions; length =     0.481     0.481     0.481
      generated from gmax (a.u.)=      9.0000: 889 vectors of 1000 (88%)
 SGVSYM: 38 symmetry stars found for 861 reciprocal lattice vectors

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol= 1.00E-06
 spec      l    rsm    eh     gmax    last term    cutoff
  A        0*   2.50  -0.01   2.974    2.30E-05      27 
  A        1    2.50  -0.01   3.093    1.29E-06      51 
  A        2    1.00  -0.01   8.508    1.16E-06     813 
  A        0    1.30  -1.00   5.718    2.28E-06     259 
  A        2    1.00  -1.00   8.508    1.16E-06     813 
  A        3    1.30  -0.01   6.806    2.09E-06     411 
 m_qplistinit:start
  --- Readin syml file --- 
   41   0.5000   0.5000   0.5000    0.0000   0.0000   0.0000 L Gamma
   41   0.0000   0.0000   0.0000    1.0000   0.0000   0.0000 Gamma X
   21   1.0000   0.0000   0.0000    1.0000   0.5000   0.0000 X W
   41   1.0000   0.5000   0.0000    0.0000   0.0000   0.0000 W Gamma
nsyml nkp=    4  144
 -------- qplist --------           4
    1   0.500   0.500   0.500  <-- isyml= 001
    2   0.487   0.487   0.487 
    3   0.475   0.475   0.475 
    4   0.463   0.463   0.463 
    5   0.450   0.450   0.450 
    6   0.438   0.438   0.438 
    7   0.425   0.425   0.425 
    8   0.412   0.412   0.412 
    9   0.400   0.400   0.400 
   10   0.388   0.388   0.388 
   11   0.375   0.375   0.375 
   12   0.362   0.362   0.362 
   13   0.350   0.350   0.350 
   14   0.338   0.338   0.338 
   15   0.325   0.325   0.325 
   16   0.312   0.312   0.312 
   17   0.300   0.300   0.300 
   18   0.287   0.287   0.287 
   19   0.275   0.275   0.275 
   20   0.263   0.263   0.263 
   21   0.250   0.250   0.250 
   22   0.237   0.237   0.237 
   23   0.225   0.225   0.225 
   24   0.213   0.213   0.213 
   25   0.200   0.200   0.200 
   26   0.188   0.188   0.188 
   27   0.175   0.175   0.175 
   28   0.162   0.162   0.162 
   29   0.150   0.150   0.150 
   30   0.138   0.138   0.138 
   31   0.125   0.125   0.125 
   32   0.112   0.112   0.112 
   33   0.100   0.100   0.100 
   34   0.088   0.088   0.088 
   35   0.075   0.075   0.075 
   36   0.062   0.062   0.062 
   37   0.050   0.050   0.050 
   38   0.037   0.037   0.037 
   39   0.025   0.025   0.025 
   40   0.013   0.013   0.013 
   41   0.000   0.000   0.000 
   42   0.000   0.000   0.000  <-- isyml= 002
   43   0.025   0.000   0.000 
   44   0.050   0.000   0.000 
   45   0.075   0.000   0.000 
   46   0.100   0.000   0.000 
   47   0.125   0.000   0.000 
   48   0.150   0.000   0.000 
   49   0.175   0.000   0.000 
   50   0.200   0.000   0.000 
   51   0.225   0.000   0.000 
   52   0.250   0.000   0.000 
   53   0.275   0.000   0.000 
   54   0.300   0.000   0.000 
   55   0.325   0.000   0.000 
   56   0.350   0.000   0.000 
   57   0.375   0.000   0.000 
   58   0.400   0.000   0.000 
   59   0.425   0.000   0.000 
   60   0.450   0.000   0.000 
   61   0.475   0.000   0.000 
   62   0.500   0.000   0.000 
   63   0.525   0.000   0.000 
   64   0.550   0.000   0.000 
   65   0.575   0.000   0.000 
   66   0.600   0.000   0.000 
   67   0.625   0.000   0.000 
   68   0.650   0.000   0.000 
   69   0.675   0.000   0.000 
   70   0.700   0.000   0.000 
   71   0.725   0.000   0.000 
   72   0.750   0.000   0.000 
   73   0.775   0.000   0.000 
   74   0.800   0.000   0.000 
   75   0.825   0.000   0.000 
   76   0.850   0.000   0.000 
   77   0.875   0.000   0.000 
   78   0.900   0.000   0.000 
   79   0.925   0.000   0.000 
   80   0.950   0.000   0.000 
   81   0.975   0.000   0.000 
   82   1.000   0.000   0.000 
   83   1.000   0.000   0.000  <-- isyml= 003
   84   1.000   0.025   0.000 
   85   1.000   0.050   0.000 
   86   1.000   0.075   0.000 
   87   1.000   0.100   0.000 
   88   1.000   0.125   0.000 
   89   1.000   0.150   0.000 
   90   1.000   0.175   0.000 
   91   1.000   0.200   0.000 
   92   1.000   0.225   0.000 
   93   1.000   0.250   0.000 
   94   1.000   0.275   0.000 
   95   1.000   0.300   0.000 
   96   1.000   0.325   0.000 
   97   1.000   0.350   0.000 
   98   1.000   0.375   0.000 
   99   1.000   0.400   0.000 
  100   1.000   0.425   0.000 
  101   1.000   0.450   0.000 
  102   1.000   0.475   0.000 
  103   1.000   0.500   0.000 
  104   1.000   0.500   0.000  <-- isyml= 004
  105   0.975   0.487   0.000 
  106   0.950   0.475   0.000 
  107   0.925   0.463   0.000 
  108   0.900   0.450   0.000 
  109   0.875   0.438   0.000 
  110   0.850   0.425   0.000 
  111   0.825   0.412   0.000 
  112   0.800   0.400   0.000 
  113   0.775   0.388   0.000 
  114   0.750   0.375   0.000 
  115   0.725   0.362   0.000 
  116   0.700   0.350   0.000 
  117   0.675   0.338   0.000 
  118   0.650   0.325   0.000 
  119   0.625   0.312   0.000 
  120   0.600   0.300   0.000 
  121   0.575   0.287   0.000 
  122   0.550   0.275   0.000 
  123   0.525   0.263   0.000 
  124   0.500   0.250   0.000 
  125   0.475   0.237   0.000 
  126   0.450   0.225   0.000 
  127   0.425   0.213   0.000 
  128   0.400   0.200   0.000 
  129   0.375   0.188   0.000 
  130   0.350   0.175   0.000 
  131   0.325   0.162   0.000 
  132   0.300   0.150   0.000 
  133   0.275   0.138   0.000 
  134   0.250   0.125   0.000 
  135   0.225   0.112   0.000 
  136   0.200   0.100   0.000 
  137   0.175   0.088   0.000 
  138   0.150   0.075   0.000 
  139   0.125   0.062   0.000 
  140   0.100   0.050   0.000 
  141   0.075   0.037   0.000 
  142   0.050   0.025   0.000 
  143   0.025   0.013   0.000 
  144   0.000   0.000   0.000 
lmv7: Read rst version ID=  2.00

 iors  : read rst restart file (binary mesh density)
          use from  restart file:use window, pnu,
          ignore in restart file:
         site   1:A       :file pnu is  4.66  4.38  3.88  4.11  5.10
         site   1:A       :file pz  is  5.50  5.50  4.50  0.00  0.00

 Basis, after reading restart file
 site spec        pos (Cartesian coordinates)         pos (multiples of plat)
   1         A  0.000000   0.000000   0.000000    0.000000   0.000000   0.000000

--- BNDFP:  begin iteration 1 of 1
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.580968
   smooth rhoves     10.096128   charge     3.702468
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.655347 -3.455244 -0.832805
  locpot:
   site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
    sm core charge in MT=  0.263001 =total-spillout=  0.267647 -  0.004646
     potential shift to crystal energy zero:    0.000043
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -6.001767      -184.375801      -190.377568
   rhoval*ves            -46.600003      -116.654177      -163.254180
   psnuc*ves              66.792258    -12971.812819    -12905.020561
   utot                   10.096128     -6544.233498     -6534.137371
   rho*exc                -2.655347      -127.279257      -129.934604
   rho*vxc                -3.455244      -168.536649      -171.991894
   valence chg             3.702468         7.297532        11.000000
   core charge            18.000000        -0.000000        18.000000
   Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
   hom background     0.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of   144 k=  0.5000  0.5000  0.5000 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     2 of   144 k=  0.4875  0.4875  0.4875 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     3 of   144 k=  0.4750  0.4750  0.4750 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     4 of   144 k=  0.4625  0.4625  0.4625 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     5 of   144 k=  0.4500  0.4500  0.4500 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     6 of   144 k=  0.4375  0.4375  0.4375 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     7 of   144 k=  0.4250  0.4250  0.4250 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     8 of   144 k=  0.4125  0.4125  0.4125 ndimh = nmto+napw =    31   31    0
 bndfp: kpt     9 of   144 k=  0.4000  0.4000  0.4000 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    10 of   144 k=  0.3875  0.3875  0.3875 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    11 of   144 k=  0.3750  0.3750  0.3750 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    12 of   144 k=  0.3625  0.3625  0.3625 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    13 of   144 k=  0.3500  0.3500  0.3500 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    14 of   144 k=  0.3375  0.3375  0.3375 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    15 of   144 k=  0.3250  0.3250  0.3250 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    16 of   144 k=  0.3125  0.3125  0.3125 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    17 of   144 k=  0.3000  0.3000  0.3000 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    18 of   144 k=  0.2875  0.2875  0.2875 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    19 of   144 k=  0.2750  0.2750  0.2750 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    20 of   144 k=  0.2625  0.2625  0.2625 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    21 of   144 k=  0.2500  0.2500  0.2500 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    22 of   144 k=  0.2375  0.2375  0.2375 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    23 of   144 k=  0.2250  0.2250  0.2250 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    24 of   144 k=  0.2125  0.2125  0.2125 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    25 of   144 k=  0.2000  0.2000  0.2000 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    26 of   144 k=  0.1875  0.1875  0.1875 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    27 of   144 k=  0.1750  0.1750  0.1750 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    28 of   144 k=  0.1625  0.1625  0.1625 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    29 of   144 k=  0.1500  0.1500  0.1500 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    30 of   144 k=  0.1375  0.1375  0.1375 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    31 of   144 k=  0.1250  0.1250  0.1250 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    32 of   144 k=  0.1125  0.1125  0.1125 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    33 of   144 k=  0.1000  0.1000  0.1000 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    34 of   144 k=  0.0875  0.0875  0.0875 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    35 of   144 k=  0.0750  0.0750  0.0750 ndimh = nmto+napw =    31   31    0
 bndfp: kpt    36 of   144 k=  0.0625  0.0625  0.0625 ndimh = nmto+napw =    31   31    0
 ... Done MPI k-loop: elapsed time=   0.4312
  Writing bands to bands file for gnuplot ...
 bndfp: kpt    1 of  144 k jsp=  0.50000  0.50000  0.50000 1 nev=   31
 bndfp: kpt    2 of  144 k jsp=  0.48750  0.48750  0.48750 1 nev=   31
 bndfp: kpt    3 of  144 k jsp=  0.47500  0.47500  0.47500 1 nev=   31
 bndfp: kpt    4 of  144 k jsp=  0.46250  0.46250  0.46250 1 nev=   31
 bndfp: kpt    5 of  144 k jsp=  0.45000  0.45000  0.45000 1 nev=   31
 bndfp: kpt    6 of  144 k jsp=  0.43750  0.43750  0.43750 1 nev=   31
 bndfp: kpt    7 of  144 k jsp=  0.42500  0.42500  0.42500 1 nev=   31
 bndfp: kpt    8 of  144 k jsp=  0.41250  0.41250  0.41250 1 nev=   31
 bndfp: kpt    9 of  144 k jsp=  0.40000  0.40000  0.40000 1 nev=   31
 bndfp: kpt   10 of  144 k jsp=  0.38750  0.38750  0.38750 1 nev=   31
 bndfp: kpt   11 of  144 k jsp=  0.37500  0.37500  0.37500 1 nev=   31
 bndfp: kpt   12 of  144 k jsp=  0.36250  0.36250  0.36250 1 nev=   31
 bndfp: kpt   13 of  144 k jsp=  0.35000  0.35000  0.35000 1 nev=   31
 bndfp: kpt   14 of  144 k jsp=  0.33750  0.33750  0.33750 1 nev=   31
 bndfp: kpt   15 of  144 k jsp=  0.32500  0.32500  0.32500 1 nev=   31
 bndfp: kpt   16 of  144 k jsp=  0.31250  0.31250  0.31250 1 nev=   31
 bndfp: kpt   17 of  144 k jsp=  0.30000  0.30000  0.30000 1 nev=   31
 bndfp: kpt   18 of  144 k jsp=  0.28750  0.28750  0.28750 1 nev=   31
 bndfp: kpt   19 of  144 k jsp=  0.27500  0.27500  0.27500 1 nev=   31
 bndfp: kpt   20 of  144 k jsp=  0.26250  0.26250  0.26250 1 nev=   31
 bndfp: kpt   21 of  144 k jsp=  0.25000  0.25000  0.25000 1 nev=   31
 bndfp: kpt   22 of  144 k jsp=  0.23750  0.23750  0.23750 1 nev=   31
 bndfp: kpt   23 of  144 k jsp=  0.22500  0.22500  0.22500 1 nev=   31
 bndfp: kpt   24 of  144 k jsp=  0.21250  0.21250  0.21250 1 nev=   31
 bndfp: kpt   25 of  144 k jsp=  0.20000  0.20000  0.20000 1 nev=   31
 bndfp: kpt   26 of  144 k jsp=  0.18750  0.18750  0.18750 1 nev=   31
 bndfp: kpt   27 of  144 k jsp=  0.17500  0.17500  0.17500 1 nev=   31
 bndfp: kpt   28 of  144 k jsp=  0.16250  0.16250  0.16250 1 nev=   31
 bndfp: kpt   29 of  144 k jsp=  0.15000  0.15000  0.15000 1 nev=   31
 bndfp: kpt   30 of  144 k jsp=  0.13750  0.13750  0.13750 1 nev=   31
 bndfp: kpt   31 of  144 k jsp=  0.12500  0.12500  0.12500 1 nev=   31
 bndfp: kpt   32 of  144 k jsp=  0.11250  0.11250  0.11250 1 nev=   31
 bndfp: kpt   33 of  144 k jsp=  0.10000  0.10000  0.10000 1 nev=   31
 bndfp: kpt   34 of  144 k jsp=  0.08750  0.08750  0.08750 1 nev=   31
 bndfp: kpt   35 of  144 k jsp=  0.07500  0.07500  0.07500 1 nev=   31
 bndfp: kpt   36 of  144 k jsp=  0.06250  0.06250  0.06250 1 nev=   31
 bndfp: kpt   37 of  144 k jsp=  0.05000  0.05000  0.05000 1 nev=   31
 bndfp: kpt   38 of  144 k jsp=  0.03750  0.03750  0.03750 1 nev=   31
 bndfp: kpt   39 of  144 k jsp=  0.02500  0.02500  0.02500 1 nev=   31
 bndfp: kpt   40 of  144 k jsp=  0.01250  0.01250  0.01250 1 nev=   31
 bndfp: kpt   41 of  144 k jsp=  0.00000  0.00000  0.00000 1 nev=   31
 bndfp: kpt   42 of  144 k jsp=  0.00000  0.00000  0.00000 1 nev=   31
 bndfp: kpt   43 of  144 k jsp=  0.02500  0.00000  0.00000 1 nev=   31
 bndfp: kpt   44 of  144 k jsp=  0.05000  0.00000  0.00000 1 nev=   31
 bndfp: kpt   45 of  144 k jsp=  0.07500  0.00000  0.00000 1 nev=   31
 bndfp: kpt   46 of  144 k jsp=  0.10000  0.00000  0.00000 1 nev=   31
 bndfp: kpt   47 of  144 k jsp=  0.12500  0.00000  0.00000 1 nev=   31
 bndfp: kpt   48 of  144 k jsp=  0.15000  0.00000  0.00000 1 nev=   31
 bndfp: kpt   49 of  144 k jsp=  0.17500  0.00000  0.00000 1 nev=   31
 bndfp: kpt   50 of  144 k jsp=  0.20000  0.00000  0.00000 1 nev=   31
 bndfp: kpt   51 of  144 k jsp=  0.22500  0.00000  0.00000 1 nev=   31
 bndfp: kpt   52 of  144 k jsp=  0.25000  0.00000  0.00000 1 nev=   31
 bndfp: kpt   53 of  144 k jsp=  0.27500  0.00000  0.00000 1 nev=   31
 bndfp: kpt   54 of  144 k jsp=  0.30000  0.00000  0.00000 1 nev=   31
 bndfp: kpt   55 of  144 k jsp=  0.32500  0.00000  0.00000 1 nev=   31
 bndfp: kpt   56 of  144 k jsp=  0.35000  0.00000  0.00000 1 nev=   31
 bndfp: kpt   57 of  144 k jsp=  0.37500  0.00000  0.00000 1 nev=   31
 bndfp: kpt   58 of  144 k jsp=  0.40000  0.00000  0.00000 1 nev=   31
 bndfp: kpt   59 of  144 k jsp=  0.42500  0.00000  0.00000 1 nev=   31
 bndfp: kpt   60 of  144 k jsp=  0.45000  0.00000  0.00000 1 nev=   31
 bndfp: kpt   61 of  144 k jsp=  0.47500  0.00000  0.00000 1 nev=   31
 bndfp: kpt   62 of  144 k jsp=  0.50000  0.00000  0.00000 1 nev=   31
 bndfp: kpt   63 of  144 k jsp=  0.52500  0.00000  0.00000 1 nev=   31
 bndfp: kpt   64 of  144 k jsp=  0.55000  0.00000  0.00000 1 nev=   31
 bndfp: kpt   65 of  144 k jsp=  0.57500  0.00000  0.00000 1 nev=   31
 bndfp: kpt   66 of  144 k jsp=  0.60000  0.00000  0.00000 1 nev=   31
 bndfp: kpt   67 of  144 k jsp=  0.62500  0.00000  0.00000 1 nev=   31
 bndfp: kpt   68 of  144 k jsp=  0.65000  0.00000  0.00000 1 nev=   31
 bndfp: kpt   69 of  144 k jsp=  0.67500  0.00000  0.00000 1 nev=   31
 bndfp: kpt   70 of  144 k jsp=  0.70000  0.00000  0.00000 1 nev=   31
 bndfp: kpt   71 of  144 k jsp=  0.72500  0.00000  0.00000 1 nev=   31
 bndfp: kpt   72 of  144 k jsp=  0.75000  0.00000  0.00000 1 nev=   31
 bndfp: kpt   73 of  144 k jsp=  0.77500  0.00000  0.00000 1 nev=   31
 bndfp: kpt   74 of  144 k jsp=  0.80000  0.00000  0.00000 1 nev=   31
 bndfp: kpt   75 of  144 k jsp=  0.82500  0.00000  0.00000 1 nev=   31
 bndfp: kpt   76 of  144 k jsp=  0.85000  0.00000  0.00000 1 nev=   31
 bndfp: kpt   77 of  144 k jsp=  0.87500  0.00000  0.00000 1 nev=   31
 bndfp: kpt   78 of  144 k jsp=  0.90000  0.00000  0.00000 1 nev=   31
 bndfp: kpt   79 of  144 k jsp=  0.92500  0.00000  0.00000 1 nev=   31
 bndfp: kpt   80 of  144 k jsp=  0.95000  0.00000  0.00000 1 nev=   31
 bndfp: kpt   81 of  144 k jsp=  0.97500  0.00000  0.00000 1 nev=   31
 bndfp: kpt   82 of  144 k jsp=  1.00000  0.00000  0.00000 1 nev=   31
 bndfp: kpt   83 of  144 k jsp=  1.00000  0.00000  0.00000 1 nev=   31
 bndfp: kpt   84 of  144 k jsp=  1.00000  0.02500  0.00000 1 nev=   31
 bndfp: kpt   85 of  144 k jsp=  1.00000  0.05000  0.00000 1 nev=   31
 bndfp: kpt   86 of  144 k jsp=  1.00000  0.07500  0.00000 1 nev=   31
 bndfp: kpt   87 of  144 k jsp=  1.00000  0.10000  0.00000 1 nev=   31
 bndfp: kpt   88 of  144 k jsp=  1.00000  0.12500  0.00000 1 nev=   31
 bndfp: kpt   89 of  144 k jsp=  1.00000  0.15000  0.00000 1 nev=   31
 bndfp: kpt   90 of  144 k jsp=  1.00000  0.17500  0.00000 1 nev=   31
 bndfp: kpt   91 of  144 k jsp=  1.00000  0.20000  0.00000 1 nev=   31
 bndfp: kpt   92 of  144 k jsp=  1.00000  0.22500  0.00000 1 nev=   31
 bndfp: kpt   93 of  144 k jsp=  1.00000  0.25000  0.00000 1 nev=   31
 bndfp: kpt   94 of  144 k jsp=  1.00000  0.27500  0.00000 1 nev=   31
 bndfp: kpt   95 of  144 k jsp=  1.00000  0.30000  0.00000 1 nev=   31
 bndfp: kpt   96 of  144 k jsp=  1.00000  0.32500  0.00000 1 nev=   31
 bndfp: kpt   97 of  144 k jsp=  1.00000  0.35000  0.00000 1 nev=   31
 bndfp: kpt   98 of  144 k jsp=  1.00000  0.37500  0.00000 1 nev=   31
 bndfp: kpt   99 of  144 k jsp=  1.00000  0.40000  0.00000 1 nev=   31
 bndfp: kpt  100 of  144 k jsp=  1.00000  0.42500  0.00000 1 nev=   31
 bndfp: kpt  101 of  144 k jsp=  1.00000  0.45000  0.00000 1 nev=   31
 bndfp: kpt  102 of  144 k jsp=  1.00000  0.47500  0.00000 1 nev=   31
 bndfp: kpt  103 of  144 k jsp=  1.00000  0.50000  0.00000 1 nev=   31
 bndfp: kpt  104 of  144 k jsp=  1.00000  0.50000  0.00000 1 nev=   31
 bndfp: kpt  105 of  144 k jsp=  0.97500  0.48750  0.00000 1 nev=   31
 bndfp: kpt  106 of  144 k jsp=  0.95000  0.47500  0.00000 1 nev=   31
 bndfp: kpt  107 of  144 k jsp=  0.92500  0.46250  0.00000 1 nev=   31
 bndfp: kpt  108 of  144 k jsp=  0.90000  0.45000  0.00000 1 nev=   31
 bndfp: kpt  109 of  144 k jsp=  0.87500  0.43750  0.00000 1 nev=   31
 bndfp: kpt  110 of  144 k jsp=  0.85000  0.42500  0.00000 1 nev=   31
 bndfp: kpt  111 of  144 k jsp=  0.82500  0.41250  0.00000 1 nev=   31
 bndfp: kpt  112 of  144 k jsp=  0.80000  0.40000  0.00000 1 nev=   31
 bndfp: kpt  113 of  144 k jsp=  0.77500  0.38750  0.00000 1 nev=   31
 bndfp: kpt  114 of  144 k jsp=  0.75000  0.37500  0.00000 1 nev=   31
 bndfp: kpt  115 of  144 k jsp=  0.72500  0.36250  0.00000 1 nev=   31
 bndfp: kpt  116 of  144 k jsp=  0.70000  0.35000  0.00000 1 nev=   31
 bndfp: kpt  117 of  144 k jsp=  0.67500  0.33750  0.00000 1 nev=   31
 bndfp: kpt  118 of  144 k jsp=  0.65000  0.32500  0.00000 1 nev=   31
 bndfp: kpt  119 of  144 k jsp=  0.62500  0.31250  0.00000 1 nev=   31
 bndfp: kpt  120 of  144 k jsp=  0.60000  0.30000  0.00000 1 nev=   31
 bndfp: kpt  121 of  144 k jsp=  0.57500  0.28750  0.00000 1 nev=   31
 bndfp: kpt  122 of  144 k jsp=  0.55000  0.27500  0.00000 1 nev=   31
 bndfp: kpt  123 of  144 k jsp=  0.52500  0.26250  0.00000 1 nev=   31
 bndfp: kpt  124 of  144 k jsp=  0.50000  0.25000  0.00000 1 nev=   31
 bndfp: kpt  125 of  144 k jsp=  0.47500  0.23750  0.00000 1 nev=   31
 bndfp: kpt  126 of  144 k jsp=  0.45000  0.22500  0.00000 1 nev=   31
 bndfp: kpt  127 of  144 k jsp=  0.42500  0.21250  0.00000 1 nev=   31
 bndfp: kpt  128 of  144 k jsp=  0.40000  0.20000  0.00000 1 nev=   31
 bndfp: kpt  129 of  144 k jsp=  0.37500  0.18750  0.00000 1 nev=   31
 bndfp: kpt  130 of  144 k jsp=  0.35000  0.17500  0.00000 1 nev=   31
 bndfp: kpt  131 of  144 k jsp=  0.32500  0.16250  0.00000 1 nev=   31
 bndfp: kpt  132 of  144 k jsp=  0.30000  0.15000  0.00000 1 nev=   31
 bndfp: kpt  133 of  144 k jsp=  0.27500  0.13750  0.00000 1 nev=   31
 bndfp: kpt  134 of  144 k jsp=  0.25000  0.12500  0.00000 1 nev=   31
 bndfp: kpt  135 of  144 k jsp=  0.22500  0.11250  0.00000 1 nev=   31
 bndfp: kpt  136 of  144 k jsp=  0.20000  0.10000  0.00000 1 nev=   31
 bndfp: kpt  137 of  144 k jsp=  0.17500  0.08750  0.00000 1 nev=   31
 bndfp: kpt  138 of  144 k jsp=  0.15000  0.07500  0.00000 1 nev=   31
 bndfp: kpt  139 of  144 k jsp=  0.12500  0.06250  0.00000 1 nev=   31
 bndfp: kpt  140 of  144 k jsp=  0.10000  0.05000  0.00000 1 nev=   31
 bndfp: kpt  141 of  144 k jsp=  0.07500  0.03750  0.00000 1 nev=   31
 bndfp: kpt  142 of  144 k jsp=  0.05000  0.02500  0.00000 1 nev=   31
 bndfp: kpt  143 of  144 k jsp=  0.02500  0.01250  0.00000 1 nev=   31
 bndfp: kpt  144 of  144 k jsp=  0.00000  0.00000  0.00000 1 nev=   31
Exit 0 procid= 0 plot band mode done
