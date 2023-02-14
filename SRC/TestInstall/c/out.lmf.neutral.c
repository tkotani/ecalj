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
 cmdl for python=/home/takao/ecalj/SRC/TestInstall/bin/ctrl2ctrlp.py  --no-iactiv c -vzbak=0<ctrl.c >ctrlp.c
rval2: STRUC_NSPEC             requ n= 1 val= 1.00000000
rval2: STRUC_NBAS              requ n= 1 val= 1.00000000
rval2: HAM_NSPIN               defa n= 1 val= 2.00000000
rval2: IO_VERBOS               defa n= 1 val= 30.00000000
rval2: IO_TIM                  defa n= 1 val= 1.00000000
rval2: STRUC_ALAT              ---- n= 1 val= 7.93700526
rval2: STRUC_DALAT             ---- n= 0 val= 
rval2: STRUC_PLAT              requ n= 9 val= 1.00000000  1.00000000  0.00000000  1.00000000  0.00000000  1.00000000  0.00000000  1.00000000  1.00000000
rval2: OPTIONS_HF              defa n= 1 val= 0.00000000
rval2: HAM_REL                 defa n= 1 val= 1.00000000
rval2: HAM_SO                  defa n= 1 val= 0.00000000
rval2: HAM_SOCAXIS             defa n= 3 val= 0.00000000  0.00000000  1.00000000
rval2: HAM_GMAX                defa n= 1 val= 0.00000000
rval2: HAM_FTMESH              defa n= 3 val= 50.00000000  50.00000000  50.00000000
rval2: HAM_TOL                 defa n= 1 val= 0.00000100
rval2: HAM_FRZWF               defa n= 1 val= 0.00000000
rval2: HAM_XCFUN               defa n= 1 val= 2.00000000
rval2: HAM_FORCES              defa n= 1 val= 12.00000000
rval2: HAM_RDSIG               defa n= 1 val= 1.00000000
rval2: HAM_ScaledSigma         defa n= 1 val= 1.00000000
rval2: HAM_EWALD               defa n= 1 val= 0.00000000
rval2: HAM_OVEPS               defa n= 1 val= 0.00000010
rval2: HAM_PWMODE              defa n= 1 val= 0.00000000
rval2: HAM_PWEMAX              defa n= 1 val= 0.00000000
rval2: HAM_READP               defa n= 1 val= 0.00000000
rval2: HAM_V0FIX               defa n= 1 val= 0.00000000
rval2: HAM_PNUFIX              defa n= 1 val= 0.00000000
rval2: SPEC_Z@1                ---- n= 1 val= 6.00000000
rval2: SPEC_R@1                ---- n= 1 val= 3.00000000
rval2: SPEC_R/W@1              ---- n= 0 val= 
rval2: SPEC_R/A@1              ---- n= 0 val= 
rval2: SPEC_A@1                defa n= 1 val= 0.02000000
rval2: SPEC_NR@1               defa n= 1 val= 0.00000000
rval2: SPEC_RSMH@1             ---- n= 4 val= 1.30000000  1.10000000 -1.00000000 -1.00000000
rval2: SPEC_EH@1               requ n= 2 val= -0.70000000 -0.20000000
rval2: SPEC_RSMH2@1            ---- n= 2 val= 0.80000000  0.80000000
rval2: SPEC_EH2@1              requ n= 2 val= -1.50000000 -1.00000000
rval2: SPEC_LMX@1              defa n= 1 val= 999.00000000
rval2: SPEC_LMXA@1             defa n= 1 val= 3.00000000
rval2: SPEC_LMXL@1             defa n= 1 val= 3.00000000
rval2: SPEC_P@1                ---- n= 4 val= 2.90000000  2.85000000  3.18000000  4.12000000
rval2: SPEC_Q@1                ---- n= 0 val= 
rval2: SPEC_MMOM@1             ---- n= 2 val= 0.00000000  2.00000000
rval2: SPEC_NMCORE@1           defa n= 1 val= 0.00000000
rval2: SPEC_PZ@1               ---- n= 0 val= 
rval2: SPEC_LFOCA@1            defa n= 1 val= 0.00000000
rval2: SPEC_KMXA@1             defa n= 1 val= 3.00000000
rval2: SPEC_RSMA@1             defa n= 1 val= 1.20000000
rval2: SPEC_IDMOD@1            ---- n= 2 val= 0.00000000  1.00000000
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
rval2: BZ_NKABC                ---- n= 3 val= 4.00000000  4.00000000  4.00000000
rval2: BZ_BZJOB                ---- n= 1 val= 0.00000000
rval2: BZ_METAL                defa n= 1 val= 2.00000000
rval2: BZ_TETRA                defa n= 1 val= 1.00000000
rval2: BZ_N                    defa n= 1 val= 0.00000000
rval2: BZ_W                    defa n= 1 val= 0.00400000
rval2: BZ_ZBAK                 defa n= 1 val= 0.00000000
rval2: BZ_SAVDOS               defa n= 1 val= 0.00000000
rval2: BZ_NPTS                 defa n= 1 val= 200.00000000
rval2: BZ_DOSMAX               defa n= 1 val= 2.93992268
rval2: BZ_EFMAX                defa n= 1 val= 5.00000000
rval2: BZ_NEVMX                defa n= 1 val= 5.00000000
rval2: BZ_FSMOM                defa n= 1 val= -99999.00000000
rval2: BZ_FSMOMMETHOD          defa n= 1 val= 0.00000000
rval2: EWALD_AS                defa n= 1 val= 2.00000000
rval2: EWALD_TOL               defa n= 1 val= 0.00000001
rval2: EWALD_NKDMX             defa n= 1 val= 300.00000000
rval2: ITER_NIT                defa n= 1 val= 10.00000000
rval2: ITER_NRMIX              defa n= 1 val= 80.00000000
rval2: ITER_CONV               defa n= 1 val= 0.00001000
rval2: ITER_CONVC              defa n= 1 val= 0.00050000
rval2: ITER_UMIX               defa n= 1 val= 0.50000000
rval2: ITER_TOLU               defa n= 1 val= 0.00000000
rval2: DYN_MODE                defa n= 1 val= 0.00000000
rval2: DYN_NIT                 defa n= 1 val= 1.00000000
rval2: DYN_HESS                defa n= 1 val= 1.00000000
rval2: DYN_XTOL                defa n= 1 val= 0.00100000
rval2: DYN_GTOL                defa n= 1 val= 0.00000000
rval2: DYN_STEP                defa n= 1 val= 0.01500000
rval2: DYN_NKILL               defa n= 1 val= 0.00000000
mixing parameters: A/B nmix wt: 0 2 1.000000  1.000000  0.000000 beta elin wc killj=  0.500000 -1.000000 -1
 ===> for --jobgw, pwmode is switched to be  0

 ... Species  1
  bndfp (warning): no sigm file found ... LDA calculation only
pnuall: j isp pnu= 1 1 2.900000  2.850000  3.180000  4.120000
pnzall: j isp  pz= 1 1 0.000000  0.000000  0.000000  0.000000
pnuall: j isp pnu= 1 2 2.900000  2.850000  3.180000  4.120000
pnzall: j isp  pz= 1 2 0.000000  0.000000  0.000000  0.000000


mmm === MTO setting ===
mmm ispec lmxb lpz nkapii nkaphh=    1    1    0    2    2
mmm rsmh1    1  1.30  1.10
mmm   eh1    1 -0.70 -0.20
mmm rsmh2    1  0.80  0.80
mmm  eh2     1 -1.50 -1.00
mmm lh       1  1
xxx isp pz= 1 0.000000  0.000000  0.000000  0.000000
 goto freats

conf:------------------------------------------------------
conf:SPEC_ATOM= C : --- Table for atomic configuration ---
conf:  isp  l  int(P) int(P)z    Qval     Qcore   CoreConf
conf:    1  0       2  0         1.000    1.000 => 1,
conf:    1  1       2  0         2.000    0.000 => 
conf:    1  2       3  0         0.000    0.000 => 
conf:    1  3       4  0         0.000    0.000 => 
conf:    2  0       2  0         1.000    1.000 => 1,
conf:    2  1       2  0         0.000    0.000 => 
conf:    2  2       3  0         0.000    0.000 => 
conf:    2  3       4  0         0.000    0.000 => 
usedQ=     2.000     2.000     0.000     0.000
conf: Species  C        Z=  6.00 Qc=  2.000 R=  3.000000 Q=  0.000000 nsp= 2 mom=  2.000000
conf: rmt rmax a=  3.000000  19.671121  0.020000 nrmt nr= 369 463
 goto atomc xxx
 atomsc nmcore=           0

 end of atomsc xxxxx
 vsum=  -59.746958882741843                1
sumev= -2.876387 etot= -74.994908 eref=  0.000000 etot-eref= -74.994908

 Free-atom wavefunctions:
 valence:      eval       node at      max at       c.t.p.   rho(r>rmt)       pnu
   2s      -1.07383         0.379       1.210       1.959     0.041897       2.913  0
   2p      -0.46569         0.000       1.175       2.960     0.100835       2.887  0
   3d       0.06953         0.000      13.212      19.671*    0.998898       3.277  0
   4f       0.11071         0.000      14.510      19.671*    0.999980       4.134  0

 core:        ecore       node at      max at       c.t.p.   rho(r>rmt)
   1s     -19.88880         0.000       0.175       0.353     0.000000

 spin 2:
   2s      -0.87119         0.378       1.233       2.011     0.054560       2.908  0
   2p      -0.27759         0.000       1.240       3.400     0.153961       2.872  0
   3d       0.07469         0.000      13.625      19.671*    0.999525       3.243  0
   4f       0.11337         0.000      14.632      19.671*    0.999988       4.128  0

 core:        ecore       node at      max at       c.t.p.   rho(r>rmt)
   1s     -19.81912         0.000       0.175       0.353     0.000000
 tailsm: init

 tailsm: fit tails to 6 smoothed hankels, rmt= 3.00000, rsm= 1.50000
    q(fit):     0.243570    rms diff:   0.000004
    fit: r>rmt  0.243570   r<rmt  1.753709   qtot  1.997279
    rho: r>rmt  0.243570   r<rmt  2.756430   qtot  3.000000
 tailsm:  fit tails to        6 functions with

 rsm=  0.15000D+01 rms error=  0.37119D-05

 tailsm: spin 2 ...
    q(fit):     0.054561    rms diff:   0.000002
    fit: r>rmt  0.054561   r<rmt  0.609878   qtot  0.664439
    rho: r>rmt  0.054561   r<rmt  0.945439   qtot  1.000000
conf: Core rhoc(rmt)= 0.000000 spillout= 0.000000
 end of freats: spid nmcore=C                  0
OK! end of LMFA ======================
 m_lmfinit:program LMF
 m_lmfinit:program LMF
INFO: Ubuntu 20.04.4 LTS \n \l
INFO: GNU Fortran (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0
INFO: -O2 -g -fimplicit-none -finit-integer=NaN -finit-real=NaN -JOBJ.gfortran -IOBJ.gfortran
INFO: MATH: -lmkl_rt
INFO: git: commit 913e769c0a5a77a2254ce7ce7011c5bc7fb5168a
INFO:    : Author: Takao Kotani <takaokotani@gmail.com>
INFO:    : Date:   Mon Feb 13 19:43:59 2023 +0900
INFO: linked at Tue Feb 14 12:44:43 JST 2023
===START LMF with   ===
 m_lmfinit:program LMF
 mpisize=           4
 m_lmfinit:program LMF
 cmdl for python=/home/takao/ecalj/SRC/TestInstall/bin/ctrl2ctrlp.py  --no-iactiv c -vzbak=0<ctrl.c >ctrlp.c
rval2: STRUC_NSPEC             requ n= 1 val= 1.00000000
rval2: STRUC_NBAS              requ n= 1 val= 1.00000000
rval2: HAM_NSPIN               defa n= 1 val= 2.00000000
rval2: IO_VERBOS               defa n= 1 val= 30.00000000
rval2: IO_TIM                  defa n= 1 val= 1.00000000
rval2: STRUC_ALAT              ---- n= 1 val= 7.93700526
rval2: STRUC_DALAT             ---- n= 0 val= 
rval2: STRUC_PLAT              requ n= 9 val= 1.00000000  1.00000000  0.00000000  1.00000000  0.00000000  1.00000000  0.00000000  1.00000000  1.00000000
rval2: OPTIONS_HF              defa n= 1 val= 0.00000000
rval2: HAM_REL                 defa n= 1 val= 1.00000000
rval2: HAM_SO                  defa n= 1 val= 0.00000000
rval2: HAM_SOCAXIS             defa n= 3 val= 0.00000000  0.00000000  1.00000000
rval2: HAM_GMAX                defa n= 1 val= 0.00000000
rval2: HAM_FTMESH              defa n= 3 val= 50.00000000  50.00000000  50.00000000
rval2: HAM_TOL                 defa n= 1 val= 0.00000100
rval2: HAM_FRZWF               defa n= 1 val= 0.00000000
rval2: HAM_XCFUN               defa n= 1 val= 2.00000000
rval2: HAM_FORCES              defa n= 1 val= 12.00000000
rval2: HAM_RDSIG               defa n= 1 val= 1.00000000
rval2: HAM_ScaledSigma         defa n= 1 val= 1.00000000
rval2: HAM_EWALD               defa n= 1 val= 0.00000000
rval2: HAM_OVEPS               defa n= 1 val= 0.00000010
rval2: HAM_PWMODE              defa n= 1 val= 0.00000000
rval2: HAM_PWEMAX              defa n= 1 val= 0.00000000
rval2: HAM_READP               defa n= 1 val= 0.00000000
rval2: HAM_V0FIX               defa n= 1 val= 0.00000000
rval2: HAM_PNUFIX              defa n= 1 val= 0.00000000
rval2: SPEC_Z@1                ---- n= 1 val= 6.00000000
rval2: SPEC_R@1                ---- n= 1 val= 3.00000000
rval2: SPEC_R/W@1              ---- n= 0 val= 
rval2: SPEC_R/A@1              ---- n= 0 val= 
rval2: SPEC_A@1                defa n= 1 val= 0.02000000
rval2: SPEC_NR@1               defa n= 1 val= 0.00000000
rval2: SPEC_RSMH@1             ---- n= 4 val= 1.30000000  1.10000000 -1.00000000 -1.00000000
rval2: SPEC_EH@1               requ n= 2 val= -0.70000000 -0.20000000
rval2: SPEC_RSMH2@1            ---- n= 2 val= 0.80000000  0.80000000
rval2: SPEC_EH2@1              requ n= 2 val= -1.50000000 -1.00000000
rval2: SPEC_LMX@1              defa n= 1 val= 999.00000000
rval2: SPEC_LMXA@1             defa n= 1 val= 3.00000000
rval2: SPEC_LMXL@1             defa n= 1 val= 3.00000000
rval2: SPEC_P@1                ---- n= 4 val= 2.90000000  2.85000000  3.18000000  4.12000000
rval2: SPEC_Q@1                ---- n= 0 val= 
rval2: SPEC_MMOM@1             ---- n= 2 val= 0.00000000  2.00000000
rval2: SPEC_NMCORE@1           defa n= 1 val= 0.00000000
mixing parameters: A/B nmix wt: 0 2 1.000000  1.000000  0.000000 beta elin wc killj=  0.500000 -1.000000 -1
mixing parameters: A/B nmix wt: 0 2 1.000000  1.000000  0.000000 beta elin wc killj=  0.500000 -1.000000 -1
rval2: SPEC_PZ@1               ---- n= 0 val= 

 ... Species  1

 ... Species  1
mixing parameters: A/B nmix wt: 0 2 1.000000  1.000000  0.000000 beta elin wc killj=  0.500000 -1.000000 -1
rval2: SPEC_LFOCA@1            defa n= 1 val= 0.00000000

 ... Species  1
rval2: SPEC_KMXA@1             defa n= 1 val= 3.00000000
rval2: SPEC_RSMA@1             defa n= 1 val= 1.20000000
rval2: SPEC_IDMOD@1            ---- n= 2 val= 0.00000000  1.00000000
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
rval2: BZ_NKABC                ---- n= 3 val= 4.00000000  4.00000000  4.00000000
rval2: BZ_BZJOB                ---- n= 1 val= 0.00000000
rval2: BZ_METAL                defa n= 1 val= 2.00000000
rval2: BZ_TETRA                defa n= 1 val= 1.00000000
rval2: BZ_N                    defa n= 1 val= 0.00000000
rval2: BZ_W                    defa n= 1 val= 0.00400000
rval2: BZ_ZBAK                 defa n= 1 val= 0.00000000
rval2: BZ_SAVDOS               defa n= 1 val= 0.00000000
rval2: BZ_NPTS                 defa n= 1 val= 200.00000000
rval2: BZ_DOSMAX               defa n= 1 val= 2.93992268
rval2: BZ_EFMAX                defa n= 1 val= 5.00000000
rval2: BZ_NEVMX                defa n= 1 val= 5.00000000
rval2: BZ_FSMOM                defa n= 1 val= -99999.00000000
rval2: BZ_FSMOMMETHOD          defa n= 1 val= 0.00000000
rval2: EWALD_AS                defa n= 1 val= 2.00000000
rval2: EWALD_TOL               defa n= 1 val= 0.00000001
rval2: EWALD_NKDMX             defa n= 1 val= 300.00000000
rval2: ITER_NIT                defa n= 1 val= 10.00000000
rval2: ITER_NRMIX              defa n= 1 val= 80.00000000
rval2: ITER_CONV               defa n= 1 val= 0.00001000
rval2: ITER_CONVC              defa n= 1 val= 0.00050000
rval2: ITER_UMIX               defa n= 1 val= 0.50000000
rval2: ITER_TOLU               defa n= 1 val= 0.00000000
rval2: DYN_MODE                defa n= 1 val= 0.00000000
rval2: DYN_NIT                 defa n= 1 val= 1.00000000
rval2: DYN_HESS                defa n= 1 val= 1.00000000
rval2: DYN_XTOL                defa n= 1 val= 0.00100000
rval2: DYN_GTOL                defa n= 1 val= 0.00000000
rval2: DYN_STEP                defa n= 1 val= 0.01500000
rval2: DYN_NKILL               defa n= 1 val= 0.00000000
mixing parameters: A/B nmix wt: 0 2 1.000000  1.000000  0.000000 beta elin wc killj=  0.500000 -1.000000 -1
 ===> for --jobgw, pwmode is switched to be  0

 ... Species  1
  bndfp (warning): no sigm file found ... LDA calculation only
pnuall: j isp pnu= 1 1 2.900000  2.850000  3.180000  4.120000
pnzall: j isp  pz= 1 1 0.000000  0.000000  0.000000  0.000000
pnuall: j isp pnu= 1 2 2.900000  2.850000  3.180000  4.120000
pnzall: j isp  pz= 1 2 0.000000  0.000000  0.000000  0.000000
 imx=           2           2           2
 imx=           3           3           3
 imx=           2           2           2
 imx=           3           3           3
 imx=           2           2           2
 imx=           3           3           3


mmm === MTO setting ===
mmm ispec lmxb lpz nkapii nkaphh=    1    1    0    2    2
mmm rsmh1    1  1.30  1.10
mmm   eh1    1 -0.70 -0.20
mmm rsmh2    1  0.80  0.80
mmm  eh2     1 -1.50 -1.00
mmm lh       1  1

                Plat                                  Qlat
   1.000000   1.000000   0.000000        0.500000   0.500000  -0.500000
   1.000000   0.000000   1.000000        0.500000  -0.500000   0.500000
   0.000000   1.000000   1.000000       -0.500000   0.500000   0.500000
  Cell vol=   1000.000000
 imx=           2           2           2
 imx=           3           3           3

 LATTC:  as= 2.000   tol= 1.00E-08   alat= 7.93701   awald= 0.200
         r1=  3.459   nkd=  79      q1=  2.571   nkg= 137
SGROUP:  1 symmetry operations from 0
 SYMLAT: Bravais system is cubic       with 48 symmetry operations.
 SYMCRY: crystal invariant under  48 symmetry operations for tol=  0.000100
 ig  group op
   1  i*i
   2  i
   3  r3(-1,1,1)
   4  i*r3(-1,1,1)
   5  r3(1,-1,-1)
   6  i*r3(1,-1,-1)
   7  r3d
   8  i*r3d
   9  r3(-1,-1,-1)
  10  i*r3(-1,-1,-1)
  11  r2z
  12  mz
  13  r4z
  14  i*r4z
  15  r4(0,0,-1)
  16  i*r4(0,0,-1)
  17  r3(-1,-1,1)
  18  i*r3(-1,-1,1)
  19  r3(1,1,-1)
  20  i*r3(1,1,-1)
  21  r2(0,1,1)
  22  m(0,1,1)
  23  r2(1,0,-1)
  24  m(1,0,-1)
  25  r2y
  26  my
  27  r4y
  28  i*r4y
  29  r4(0,-1,0)
  30  i*r4(0,-1,0)
  31  r2(1,-1,0)
  32  m(1,-1,0)
  33  r2x
  34  mx
  35  r4(-1,0,0)
  36  i*r4(-1,0,0)
  37  r4x
  38  i*r4x
  39  r3(-1,1,-1)
  40  i*r3(-1,1,-1)
  41  r3(1,-1,1)
  42  i*r3(1,-1,1)
  43  r2(1,0,1)
  44  m(1,0,1)
  45  r2(0,1,-1)
  46  m(0,1,-1)
  47  r2(1,1,0)
  48  m(1,1,0)
GROUPG: the following are sufficient to generate the space group:
 Generators:trans(cart)  = i*r3(-1,1,1) r4z
 Generators::trans(frac) = i*r3(-1,1,1) r4z
MKSYM: found  48  space group operations
SPLCLS: ibas iclass ispec label(ispec)
 SPLCLS     1    1    1     C
 BZMESH:      8 irreducible QP from    4   4   4 shift=FFF
 TETIRR: sorting      384 tetrahedra ...
 SGVSYM: 1207 symmetry stars found for 45911 reciprocal lattice vectors

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol= 1.00E-06
 spec      l    rsm    eh     gmax    last term    cutoff
  C        0    1.30  -0.70   5.718    1.05E-06    3143 
  C        1    1.10  -0.20   7.226    1.06E-06    6375 
  C        0    0.80  -1.50   9.292    1.08E-06   13539 
  C        1    0.80  -1.00  10.038    1.00E-06   16961 
 m_qplistinit:start

 iors  : read rst restart file (binary mesh density)
 iors  : empty file ... nothing read

rdovfa: read and overlap free-atom densities (mesh density) ...
 rdovfa: expected C,       read C        with rmt=  3.0000  mesh   369  0.020
  ovlpfa: overlap smooth part of FA densities

 Free atom and overlapped crystal site charges:
   ib    true(FA)    smooth(FA)  true(OV)    smooth(OV)    local
    1    3.701869    2.363587    3.701843    2.363561    1.338282
 amom    1.810990    1.143831    1.810990    1.143831    0.667159

 Smooth charge on mesh:            2.661718    moment    1.332841
 Sum of local charges:             1.338282    moments   0.667159
 Total valence charge:             4.000000    moment    2.000000
 Sum of core charges:              2.000000    moment    0.000000
 Sum of nuclear charges:          -6.000000
 Homogeneous background:           0.000000
 Deviation from neutrality:       -0.000000

 Basis, after reading restart file
 site spec        pos (Cartesian coordinates)         pos (multiples of plat)
   1         C  0.000000   0.000000   0.000000    0.000000   0.000000   0.000000

--- BNDFP:  begin iteration 1 of 10
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.006607
   smooth rhoves      2.063910   charge     2.661718
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -1.109320 -1.522780 -0.218128
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -0.385581 -0.423770 -0.164568
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.90000   2.90000   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.18000   3.18000   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.12000   4.12000   0.00000   0.00000
     potential shift to crystal energy zero:    0.000003
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -3.931613       -10.430885       -14.362499
   rhoval*ves             -5.058553        -5.181775       -10.240327
   psnuc*ves               9.186374      -278.836573      -269.650199
   utot                    2.063910      -142.009174      -139.945263
   rho*exc                -1.494901        -8.107094        -9.601995
   rho*vxc                -1.946550       -10.686870       -12.633421
   valence chg             2.661718         1.338282         4.000000
   valence mag             1.332841         0.667159         2.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
   hom background     0.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.5622

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.430662;   4.000000 electrons
         Sum occ. bands:   -2.743234, incl. Bloechl correction: -0.000179
         Mag. moment:       2.000000
 bndfp:Generating TDOS: efermi= -0.430662  dos window emin emax=  -1.050710  2.509261


 m_bandcal_2nd: to fill eigenfunctions**2 up to Efermi
       contr. to mm extrapolated for r>rmt:   0.163680 est. true mm = 1.960270
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85523  sum tc=    31.38701  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78555  sum tc=    31.54499  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.686857    3.838105   -0.151248      1.796590    2.141140   -0.344551
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.743234 Vin*nin=     -14.362499 Ek=Eb-Vin*nin=      11.619265
 Ek(core)=      62.932005 Exc=      -9.601995 Ees=    -139.945263 Eharris=     -74.995989

 mkekin:
   nout*Vin = smpart,onsite,total=:     -7.491328     -6.855971    -14.347299
    E_B(band energy sum)=   -2.743234  E_B-nout*Vin=   11.604065

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.008341
   smooth rhoves      3.531040   charge     4.151248
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.415659 -3.336154 -0.237894
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -0.665883 -0.689339 -0.174538
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91338   2.90760   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16628   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10633   4.10242   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -9.130203        -5.227578       -14.357782
   rhoval*ves             -4.661420        -5.587089       -10.248509
   psnuc*ves              11.723500      -281.355190      -269.631690
   utot                    3.531040      -143.471140      -139.940099
   rho*exc                -3.081542        -6.510575        -9.592117
   rho*vxc                -4.025493        -8.594994       -12.620487
   valence chg             4.151248        -0.151248         4.000000
   valence mag             2.344550        -0.344551         2.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
   hom background     0.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       11.604065 Ekcore=        62.932006 Ektot    =       74.536071
 Exc=      -9.592117 Ees   =      -139.940099 EKohnSham=      -74.996145
 Magnetic moment=     2.000000
  smvxcm: all smrho_w is positive
  smvxcm: all smrho_w is positive
  smvxcm: all smrho_w is positive

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00   -0.00   -0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force =   0.000000 mRy/au (site 1 )
  Symmetrize forces ...
 mixrealsmooth= T
 wgtsmooth=   2.8284271247461901E-003
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 mixrho: sought 2 iter from file mixm ; read 0 RMS DQ= 8.23E-3
 AMIX: nmix=0 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 4.12D-03
 mixrho: add corrections to qcell smrho =  0.74319D-07  0.37159D-10

 iors  : write rst restart file (binary mesh density)

   it  1  of 10    ehf=     -74.995989   ehk=     -74.996145
h mmom= 2.0000 ehf(eV)=-1020.380425 ehk(eV)=-1020.382556 sev(eV)=-37.323889

--- BNDFP:  begin iteration 2 of 10
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.007474
   smooth rhoves      2.744964   charge     3.406483
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -1.724195 -2.376346 -0.229068
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -0.520836 -0.552450 -0.169961
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91338   2.90760   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16628   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10633   4.10242   0.00000   0.00000
     potential shift to crystal energy zero:    0.000004
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -6.342509        -8.017597       -14.360106
   rhoval*ves             -4.965008        -5.279444       -10.244452
   psnuc*ves              10.454937      -280.095881      -269.640945
   utot                    2.744964      -142.687663      -139.942698
   rho*exc                -2.245031        -7.351973        -9.597004
   rho*vxc                -2.928796        -9.698090       -12.626886
   valence chg             3.406483         0.593517         4.000000
   valence mag             1.838695         0.161304         2.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
   hom background     0.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.0679

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.431787;   4.000000 electrons
         Sum occ. bands:   -2.748502, incl. Bloechl correction: -0.000182
         Mag. moment:       2.000000
 bndfp:Generating TDOS: efermi= -0.431787  dos window emin emax=  -1.051996  2.508135

       contr. to mm extrapolated for r>rmt:   0.163407 est. true mm = 1.960158
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85779  sum tc=    31.38712  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78816  sum tc=    31.54494  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.687503    3.843381   -0.155877      1.796752    2.129366   -0.332614
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.748502 Vin*nin=     -14.360106 Ek=Eb-Vin*nin=      11.611604
 Ek(core)=      62.932005 Exc=      -9.597004 Ees=    -139.942698 Eharris=     -74.996093

 mkekin:
   nout*Vin = smpart,onsite,total=:     -8.378297     -5.986619    -14.364916
    E_B(band energy sum)=   -2.748502  E_B-nout*Vin=   11.616414

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.008327
   smooth rhoves      3.521831   charge     4.155877
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.413446 -3.333529 -0.237811
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -0.674101 -0.699705 -0.174472
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91349   2.90773   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16605   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10625   4.10242   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -9.152986        -5.213019       -14.366005
   rhoval*ves             -4.663530        -5.590989       -10.254520
   psnuc*ves              11.707192      -281.353960      -269.646768
   utot                    3.521831      -143.472474      -139.950644
   rho*exc                -3.087546        -6.506412        -9.593958
   rho*vxc                -4.033234        -8.589682       -12.622916
   valence chg             4.155877        -0.155877         4.000000
   valence mag             2.332614        -0.332614         2.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
   hom background     0.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       11.616414 Ekcore=        62.932059 Ektot    =       74.548472
 Exc=      -9.593958 Ees   =      -139.950644 EKohnSham=      -74.996130
 Magnetic moment=     2.000000
  smvxcm: all smrho_w is positive
  smvxcm: all smrho_w is positive
  smvxcm: all smrho_w is positive
 Maximum Harris force =   0.000000 mRy/au (site 1 )

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
 Maximum Harris force =   0.000000 mRy/au (site 1 )
   1   -0.00   -0.00   -0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force =   0.000000 mRy/au (site 1 )
  Symmetrize forces ...
 wgtsmooth=   2.8284271247461901E-003
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 mixrho: sought 2 iter from file mixm ; read 1 RMS DQ= 4.12E-3  last it= 8.23E-3
 AMIX: nmix=1 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 2.06D-03
   tj:-1.00050
 mixrho: add corrections to qcell smrho =  0.95182D-08  0.47591D-11

 iors  : write rst restart file (binary mesh density)

   it  2  of 10    ehf=     -74.996093   ehk=     -74.996130
 From last iter    ehf=     -74.995989   ehk=     -74.996145
 diffe(q)= -0.000104 (0.004120)    tol= 0.000010 (0.000500)   more=T
i mmom= 2.0000 ehf(eV)=-1020.381841 ehk(eV)=-1020.382343 sev(eV)=-37.395569

--- BNDFP:  begin iteration 3 of 10
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.008328
   smooth rhoves      3.522038   charge     4.156065
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.413626 -3.333780 -0.237813
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -0.674141 -0.699744 -0.174474
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91349   2.90773   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16605   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10625   4.10242   0.00000   0.00000
     potential shift to crystal energy zero:    0.000005
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -9.153736        -5.212271       -14.366007
   rhoval*ves             -4.663429        -5.591094       -10.254523
   psnuc*ves              11.707506      -281.354225      -269.646719
   utot                    3.522038      -143.472659      -139.950621
   rho*exc                -3.087767        -6.506189        -9.593956
   rho*vxc                -4.033523        -8.589389       -12.622912
   valence chg             4.156065        -0.156065         4.000000
   valence mag             2.332738        -0.332738         2.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
   hom background     0.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.0868

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.431877;   4.000000 electrons
         Sum occ. bands:   -2.749580, incl. Bloechl correction: -0.000186
         Mag. moment:       2.000000
 bndfp:Generating TDOS: efermi= -0.431877  dos window emin emax=  -1.052215  2.508045

       contr. to mm extrapolated for r>rmt:   0.163406 est. true mm = 1.959998
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85748  sum tc=    31.38660  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78781  sum tc=    31.54437  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.687633    3.847282   -0.159649      1.796592    2.119617   -0.323026
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.749580 Vin*nin=     -14.366007 Ek=Eb-Vin*nin=      11.616427
 Ek(core)=      62.932032 Exc=      -9.593956 Ees=    -139.950621 Eharris=     -74.996118

 mkekin:
   nout*Vin = smpart,onsite,total=:     -9.174978     -5.194807    -14.369786
    E_B(band energy sum)=   -2.749580  E_B-nout*Vin=   11.620206

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.008290
   smooth rhoves      3.515406   charge     4.159649
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.411330 -3.330903 -0.237828
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -0.680567 -0.707929 -0.174472
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91356   2.90781   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16589   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10619   4.10242   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -9.171543        -5.196916       -14.368459
   rhoval*ves             -4.665869        -5.590439       -10.256308
   psnuc*ves              11.696681      -281.346079      -269.649398
   utot                    3.515406      -143.468259      -139.952853
   rho*exc                -3.091897        -6.502552        -9.594448
   rho*vxc                -4.038831        -8.584729       -12.623561
   valence chg             4.159649        -0.159649         4.000000
   valence mag             2.323026        -0.323026         2.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
   hom background     0.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       11.620206 Ekcore=        62.930978 Ektot    =       74.551184
 Exc=      -9.594448 Ees   =      -139.952853 EKohnSham=      -74.996118
 Magnetic moment=     2.000000
  smvxcm: all smrho_w is positive
  smvxcm: all smrho_w is positive
  smvxcm: all smrho_w is positive
 Maximum Harris force =   0.000000 mRy/au (site 1 )

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00   -0.00   -0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force =   0.000000 mRy/au (site 1 )
  Symmetrize forces ...
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 wgtsmooth=   2.8284271247461901E-003
 mixrho: sought 2 iter from file mixm ; read 2 RMS DQ= 6.10E-5  last it= 4.12E-3
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 3.05D-05
   tj:-1.32030   0.66035
 mixrho: add corrections to qcell smrho =  0.62815D-07  0.31407D-10

 iors  : write rst restart file (binary mesh density)

   it  3  of 10    ehf=     -74.996118   ehk=     -74.996118
 From last iter    ehf=     -74.996093   ehk=     -74.996130
 diffe(q)= -0.000025 (0.000061)    tol= 0.000010 (0.000500)   more=T
i mmom= 2.0000 ehf(eV)=-1020.382181 ehk(eV)=-1020.382179 sev(eV)=-37.410234

--- BNDFP:  begin iteration 4 of 10
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.008296
   smooth rhoves      3.516503   charge     4.159013
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.411693 -3.331355 -0.237825
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -0.679468 -0.706529 -0.174472
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91356   2.90781   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16589   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10619   4.10242   0.00000   0.00000
     potential shift to crystal energy zero:    0.000005
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -9.168406        -5.199632       -14.368039
   rhoval*ves             -4.665470        -5.590532       -10.256002
   psnuc*ves              11.698476      -281.348080      -269.649604
   utot                    3.516503      -143.469306      -139.952803
   rho*exc                -3.091161        -6.503223        -9.594383
   rho*vxc                -4.037884        -8.585591       -12.623475
   valence chg             4.159013        -0.159013         4.000000
   valence mag             2.324659        -0.324659         2.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
   hom background     0.00000   deviation from neutrality:      0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.1067

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.431601;   4.000000 electrons
         Sum occ. bands:   -2.748478, incl. Bloechl correction: -0.000186
         Mag. moment:       2.000000
 bndfp:Generating TDOS: efermi= -0.431601  dos window emin emax=  -1.051927  2.508322

       contr. to mm extrapolated for r>rmt:   0.163546 est. true mm = 1.959960
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85678  sum tc=    31.38650  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78710  sum tc=    31.54427  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.687410    3.843876   -0.156466      1.796414    2.117499   -0.321085
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.748478 Vin*nin=     -14.368039 Ek=Eb-Vin*nin=      11.619561
 Ek(core)=      62.931505 Exc=      -9.594383 Ees=    -139.952803 Eharris=     -74.996120

 mkekin:
   nout*Vin = smpart,onsite,total=:     -9.160098     -5.206122    -14.366220
    E_B(band energy sum)=   -2.748478  E_B-nout*Vin=   11.617742

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.008266
   smooth rhoves      3.513827   charge     4.156466
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.407917 -3.326145 -0.237866
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -0.679707 -0.707080 -0.174495
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91355   2.90780   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16591   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10620   4.10242   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -9.157930        -5.208957       -14.366887
   rhoval*ves             -4.667240        -5.588005       -10.255246
   psnuc*ves              11.694894      -281.340877      -269.645984
   utot                    3.513827      -143.464441      -139.950615
   rho*exc                -3.087623        -6.506395        -9.594018
   rho*vxc                -4.033225        -8.589768       -12.622993
   valence chg             4.156466        -0.156466         4.000000
   valence mag             2.321085        -0.321085         2.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
   hom background     0.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       11.617742 Ekcore=        62.930771 Ektot    =       74.548513
 Exc=      -9.594018 Ees   =      -139.950615 EKohnSham=      -74.996120
 Magnetic moment=     2.000000
  smvxcm: all smrho_w is positive
  smvxcm: all smrho_w is positive
  smvxcm: all smrho_w is positive
 Maximum Harris force =   0.000000 mRy/au (site 1 )

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
   1    0.00    0.00    0.00     0.00    0.00    0.00    -0.00   -0.00   -0.00
 shift forces to make zero average correction:           -0.00   -0.00   -0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force =   0.000000 mRy/au (site 1 )
  Symmetrize forces ...
 wgtsmooth=   2.8284271247461901E-003
 mixrho: sought 2 iter from file mixm ; read 3 RMS DQ= 2.13E-5  last it= 6.10E-5
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 1.06D-05
   tj:-0.14043   0.00544
 mixrho: add corrections to qcell smrho =  0.41475D-07  0.20738D-10

 iors  : write rst restart file (binary mesh density)

   it  4  of 10    ehf=     -74.996120   ehk=     -74.996120
 From last iter    ehf=     -74.996118   ehk=     -74.996118
 diffe(q)= -0.000002 (0.000021)    tol= 0.000010 (0.000500)   more=F
c mmom= 2.0000 ehf(eV)=-1020.382207 ehk(eV)=-1020.382206 sev(eV)=-37.395238
Exit 0 procid= 0 OK! end of LMF ======================
