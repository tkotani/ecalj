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
 cmdl for python=/home/takao/ecalj/SRC/TestInstall/bin/ctrl2ctrlp.py  --no-iactiv c -vzbak=1<ctrl.c >ctrlp.c
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
rval2: BZ_ZBAK                 defa n= 1 val= 1.00000000
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
 mpisize=           4
 m_lmfinit:program LMF
 m_lmfinit:program LMF
 cmdl for python=/home/takao/ecalj/SRC/TestInstall/bin/ctrl2ctrlp.py  --no-iactiv c -vzbak=1<ctrl.c >ctrlp.c
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
mixing parameters: A/B nmix wt: 0 2 1.000000  1.000000  0.000000 beta elin wc killj=  0.500000 -1.000000 -1

 ... Species  1
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
mixing parameters: A/B nmix wt: 0 2 1.000000  1.000000  0.000000 beta elin wc killj=  0.500000 -1.000000 -1
rval2: SPEC_LFOCA@1            defa n= 1 val= 0.00000000
mixing parameters: A/B nmix wt: 0 2 1.000000  1.000000  0.000000 beta elin wc killj=  0.500000 -1.000000 -1


 ... Species  1
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
rval2: BZ_ZBAK                 defa n= 1 val= 1.00000000
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
 imx=           2           2           2
 imx=           3           3           3
 imx=           2           2           2
 imx=           3           3           3
 imx=           2           2           2
 imx=           3           3           3
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
 Uniform density added to neutralize background q=  1.000000

 Smooth charge on mesh:            1.661718    moment    1.332841
 Sum of local charges:             1.338282    moments   0.667159
 Total valence charge:             3.000000    moment    2.000000
 Sum of core charges:              2.000000    moment    0.000000
 Sum of nuclear charges:          -6.000000
 Homogeneous background:           1.000000
 Deviation from neutrality:       -0.000000

 Basis, after reading restart file
 site spec        pos (Cartesian coordinates)         pos (multiples of plat)
   1         C  0.000000   0.000000   0.000000    0.000000   0.000000   0.000000

--- BNDFP:  begin iteration 1 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.006607
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      2.063910   charge     2.661718
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      201254  -4.9989572695961402E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.49989572696961405E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -1.109298 -1.522749 -0.218032
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -0.385563 -0.423750 -0.164244
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       128
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       128
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       128
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       128
 vxcnsp (warning): negative rho: min val =  -4.61E-04
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.90000   2.90000   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.18000   3.18000   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.12000   4.12000   0.00000   0.00000
     potential shift to crystal energy zero:    0.000003
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -3.733857       -10.403599       -14.137456
   rhoval*ves             -5.058553        -5.181775       -10.240327
   psnuc*ves               9.186374      -278.836573      -269.650199
   utot                    2.063910      -142.009174      -139.945263
   rho*exc                -1.494861        -8.100688        -9.595549
   rho*vxc                -1.946499       -10.678430       -12.624928
   valence chg             1.661718         1.338282         3.000000
   valence mag             1.332841         0.667159         2.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.7637

 bzwts: --- Tetrahedron Integration ---
 ... only filled or empty bands encountered: ev= -0.824805 ec= -0.772717
 VBmax= -0.824805 CBmin= -0.772717 gap =  0.052087 Ry =   0.708689 eV
 BZINTS: Fermi energy:     -0.824805;   3.000000 electrons
         Sum occ. bands:   -3.277726, incl. Bloechl correction:  0.000000
         Mag. moment:      -1.000000
 bndfp:Generating TDOS: efermi= -0.824805  dos window emin emax=  -1.427957  2.115118


 m_bandcal_2nd: to fill eigenfunctions**2 up to Efermi
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85491  sum tc=    31.38764  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78517  sum tc=    31.54593  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.130156  386.510556 -384.380399     -0.144414 -327.846615  327.702201
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -3.277726 Vin*nin=     -14.137456 Ek=Eb-Vin*nin=      10.859730
 Ek(core)=      62.932005 Exc=      -9.595549 Ees=    -139.945263 Eharris=     -75.749078

 mkekin:
   nout*Vin = smpart,onsite,total=:  -1237.408125   1225.635237    -11.772888
    E_B(band energy sum)=   -3.277726  E_B-nout*Vin=    8.495162

 m_mkpot_energyterms
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.205342
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves     38.909184   charge   388.380399
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -179.107593 -116.221773 -0.227596
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -2291.203071 -3162.637447 -0.404877
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.89128   2.74142   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff          -2698.892631      2687.560586       -11.332044
   rhoval*ves             82.992199       -91.768732        -8.776534
   psnuc*ves              -5.173830      -258.804129      -263.977959
   utot                   38.909184      -175.286431      -136.377246
   rho*exc             -2470.310664      2462.083812        -8.226851
   rho*vxc             -3278.859220      3268.028480       -10.830740
   valence chg           387.380399      -384.380399         3.000000
   valence mag          -328.702201       327.702201        -1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=        8.495162 Ekcore=        62.933572 Ektot    =       71.428735
 Exc=      -8.226851 Ees   =      -136.377246 EKohnSham=      -73.175363
 Magnetic moment=    -1.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      201254  -4.9989874524264150E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.49989874525264153E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      201254  -4.9989270867658643E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.49989270868658646E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      201254  -4.9989572695961391E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.49989572696961394E-3
 Maximum Harris force =   0.000000 mRy/au (site 1 )

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
   1   -0.00   -0.00   -0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force =   0.000000 mRy/au (site 1 )
  Symmetrize forces ...
 mixrealsmooth= T
 wgtsmooth=   2.8284271247461901E-003
 mixrho: sought 2 iter from file mixm ; read 0 RMS DQ= 3.66E+0
 AMIX: nmix=0 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 1.83D+00
 mixrho: warning. negative smrho; isp number min=       1   93323 -0.24846D-03
 mixrho: warning. negative smrho; isp number min=       2   60589 -0.24856D-03
 mixrho: add corrections to qcell smrho =  0.30298D-07  0.15149D-10
 mixrho: warning. negative smrho; isp number min=       1   93323 -0.24846D-03
 mixrho: warning. negative smrho; isp number min=       1   93323 -0.24846D-03
 mixrho: warning. negative smrho; isp number min=       1   93323 -0.24846D-03
 mixrho: warning. negative smrho; isp number min=       2   60589 -0.24856D-03
 mixrho: warning. negative smrho; isp number min=       2   60589 -0.24856D-03
 mixrho: warning. negative smrho; isp number min=       2   60589 -0.24856D-03

 iors  : write rst restart file (binary mesh density)

   it  1  of 10    ehf=     -75.749078   ehk=     -73.175363
h mmom=-1.0000 ehf(eV)=-1030.626808 ehk(eV)=-995.609353 sev(eV)=-44.596081

--- BNDFP:  begin iteration 2 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.099368
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      8.892064   charge   195.521059
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      153912  -2.4856375869779305E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.24856375870779303E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -73.344149 -48.831148 -0.245284
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -913.875938 -1260.393217 -0.337806
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0        64
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0        64
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0        64
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0        64
 vxcnsp (warning): negative rho: min val =  -9.48E-05
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.89128   2.74142   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
     potential shift to crystal energy zero:    0.000151
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff          -1389.461763      1376.441884       -13.019879
   rhoval*ves             15.777856       -25.671329        -9.893473
   psnuc*ves               2.006272      -268.820350      -266.814079
   utot                    8.892064      -147.245840      -138.353776
   rho*exc              -987.220087       978.344548        -8.875539
   rho*vxc             -1309.224365      1297.544822       -11.679543
   valence chg           194.521059      -191.521059         3.000000
   valence mag          -163.684680       164.184680         0.500000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:      0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.4501

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.720401;   3.000000 electrons
         Sum occ. bands:   -3.283332, incl. Bloechl correction: -0.000018
         Mag. moment:       1.000000
 bndfp:Generating TDOS: efermi= -0.720401  dos window emin emax=  -1.349253  2.219522

       contr. to mm extrapolated for r>rmt:  -0.052767 est. true mm =-0.991699
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.56207  sum tc=    31.48865  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.52663  sum tc=    31.57660  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.897679    6.697088   -3.799409     -0.938932   -0.697979   -0.240952
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -3.283332 Vin*nin=     -13.019879 Ek=Eb-Vin*nin=       9.736547
 Ek(core)=      62.932636 Exc=      -8.875539 Ees=    -138.353776 Eharris=     -74.560131

 mkekin:
   nout*Vin = smpart,onsite,total=:    -32.819422     18.717570    -14.101852
    E_B(band energy sum)=   -3.283332  E_B-nout*Vin=   10.818520

 m_mkpot_energyterms
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.114016
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.285346   charge     7.799409
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -3.640144 -4.711126 -0.139094
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -4.068828 -5.379646 -0.180306
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92271   2.91873   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -28.173204        14.288457       -13.884747
   rhoval*ves             -4.065857        -6.547570       -10.613428
   psnuc*ves              14.636549      -283.137305      -268.500756
   utot                    5.285346      -144.842438      -139.557092
   rho*exc                -7.708971        -1.168246        -8.877217
   rho*vxc               -10.090771        -1.594247       -11.685018
   valence chg             6.799409        -3.799409         3.000000
   valence mag            -0.759048        -0.240952        -1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       10.818520 Ekcore=        63.065249 Ektot    =       73.883769
 Exc=      -8.877217 Ees   =      -139.557092 EKohnSham=      -74.550540
 Magnetic moment=    -1.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      153912  -2.4856526355463983E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.24856526356463980E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      153912  -2.4856225384094628E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.24856225385094625E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      153912  -2.4856375869779305E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.24856375870779303E-3

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
 Maximum Harris force =   0.000000 mRy/au (site 1 )
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
 mixrho: sought 2 iter from file mixm ; read 1 RMS DQ= 1.81E+0  last it= 3.66E+0
 AMIX: nmix=1 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 9.04D-01
   tj: 0.33092
 mixrho: warning. negative smrho; isp number min=       1   92435 -0.16534D-03
 mixrho: warning. negative smrho; isp number min=       1   92435 -0.16534D-03
 mixrho: add corrections to qcell smrho =  0.77475D-07  0.38737D-10
 mixrho: warning. negative smrho; isp number min=       1   92435 -0.16534D-03
 mixrho: warning. negative smrho; isp number min=       2   59461 -0.16440D-03
 mixrho: warning. negative smrho; isp number min=       2   59461 -0.16440D-03
 mixrho: warning. negative smrho; isp number min=       2   59461 -0.16440D-03
 mixrho: warning. negative smrho; isp number min=       1   92435 -0.16534D-03
 mixrho: warning. negative smrho; isp number min=       2   59461 -0.16440D-03

 iors  : write rst restart file (binary mesh density)

   it  2  of 10    ehf=     -74.560131   ehk=     -74.550540
 From last iter    ehf=     -75.749078   ehk=     -73.175363
 diffe(q)=  1.188947 (1.808465)    tol= 0.000010 (0.000500)   more=T
i mmom=-1.0000 ehf(eV)=-1014.450237 ehk(eV)=-1014.319734 sev(eV)=-44.672355

--- BNDFP:  begin iteration 3 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.104268
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      4.192906   charge   132.720806
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      151896  -1.6533819262874995E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.16533819263874996E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -46.537078 -32.008653 -0.223141
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -536.842477 -741.177660 -0.306197
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92271   2.91873   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
     potential shift to crystal energy zero:    0.000103
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff           -971.582823       958.305306       -13.277517
   rhoval*ves              2.154216       -12.296575       -10.142359
   psnuc*ves               6.231595      -273.651462      -267.419867
   utot                    4.192906      -142.974018      -138.781113
   rho*exc              -583.379555       574.536976        -8.842578
   rho*vxc              -773.186314       761.550341       -11.635973
   valence chg           131.720806      -128.720806         3.000000
   valence mag          -109.179668       109.177859        -0.001809
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.4438

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.646428;   3.000000 electrons
         Sum occ. bands:   -3.147175, incl. Bloechl correction: -0.000023
         Mag. moment:       1.000000
 bndfp:Generating TDOS: efermi= -0.646428  dos window emin emax=  -1.276999  2.293495

       contr. to mm extrapolated for r>rmt:   0.034243 est. true mm = 0.993974
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.45100  sum tc=    31.49394  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.44268  sum tc=    31.51704  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.907928    6.731021   -3.823093      0.959730    2.095335   -1.135605
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -3.147175 Vin*nin=     -13.277517 Ek=Eb-Vin*nin=      10.130341
 Ek(core)=      62.998923 Exc=      -8.842578 Ees=    -138.781113 Eharris=     -74.494427

 mkekin:
   nout*Vin = smpart,onsite,total=:    -31.710448     17.619456    -14.090992
    E_B(band energy sum)=   -3.147175  E_B-nout*Vin=   10.943817

 m_mkpot_energyterms
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.112609
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.512023   charge     7.823093
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.042408 -6.941980 -0.172180
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -2.752062 -3.265474 -0.139132
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92208   2.92034   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -28.397203        14.515305       -13.881898
   rhoval*ves             -3.855586        -6.740103       -10.595689
   psnuc*ves              14.879632      -283.304591      -268.424959
   utot                    5.512023      -145.022347      -139.510324
   rho*exc                -7.794470        -1.091589        -8.886059
   rho*vxc               -10.207453        -1.489242       -11.696695
   valence chg             6.823093        -3.823093         3.000000
   valence mag             2.135605        -1.135605         1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       10.943817 Ekcore=        63.010975 Ektot    =       73.954791
 Exc=      -8.886059 Ees   =      -139.510324 EKohnSham=      -74.441592
 Magnetic moment=     1.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      151896  -1.6533919099703814E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.16533919100703815E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      151896  -1.6533719426046177E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.16533719427046177E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      151896  -1.6533819262874995E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.16533819263874996E-3
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
 mixrho: sought 2 iter from file mixm ; read 2 RMS DQ= 1.21E+0  last it= 1.81E+0
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 6.03D-01
   tj:-0.07971   0.23828
 mixrho: warning. negative smrho; isp number min=       1   88535 -0.11834D-03
 mixrho: warning. negative smrho; isp number min=       1   88535 -0.11834D-03
 mixrho: warning. negative smrho; isp number min=       2   59749 -0.11860D-03
 mixrho: warning. negative smrho; isp number min=       2   59749 -0.11860D-03
 mixrho: warning. negative smrho; isp number min=       1   88535 -0.11834D-03
 mixrho: warning. negative smrho; isp number min=       2   59749 -0.11860D-03
 mixrho: add corrections to qcell smrho =  0.32690D-07  0.16345D-10
 mixrho: warning. negative smrho; isp number min=       1   88535 -0.11834D-03
 mixrho: warning. negative smrho; isp number min=       2   59749 -0.11860D-03

 iors  : write rst restart file (binary mesh density)

   it  3  of 10    ehf=     -74.494427   ehk=     -74.441592
 From last iter    ehf=     -74.560131   ehk=     -74.550540
 diffe(q)=  0.065705 (1.206786)    tol= 0.000010 (0.000500)   more=T
i mmom= 1.0000 ehf(eV)=-1013.556271 ehk(eV)=-1012.837409 sev(eV)=-42.819840

--- BNDFP:  begin iteration 4 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.106416
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      3.097188   charge    97.613905
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      148284  -1.1859896474127728E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.11859896475127728E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -34.193498 -24.702575 -0.215237
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -347.843026 -481.320288 -0.278762
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92208   2.92034   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
     potential shift to crystal energy zero:    0.000076
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff           -712.103453       698.646721       -13.456732
   rhoval*ves             -2.502218        -7.768902       -10.271120
   psnuc*ves               8.696595      -276.407709      -267.711114
   utot                    3.097188      -142.088305      -138.991117
   rho*exc              -382.036523       373.179239        -8.857284
   rho*vxc              -506.022863       494.366102       -11.656761
   valence chg            96.613905       -93.613905         3.000000
   valence mag           -77.483537        78.042557         0.559020
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.4375

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.651371;   3.000000 electrons
         Sum occ. bands:   -3.100739, incl. Bloechl correction: -0.000022
         Mag. moment:       1.000000
 bndfp:Generating TDOS: efermi= -0.651371  dos window emin emax=  -1.282501  2.288552

       contr. to mm extrapolated for r>rmt:   0.031971 est. true mm = 0.994234
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.40462  sum tc=    31.45830  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.37206  sum tc=    31.53274  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.907057    6.690772   -3.783715      0.962264    2.201512   -1.239248
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -3.100739 Vin*nin=     -13.456732 Ek=Eb-Vin*nin=      10.355992
 Ek(core)=      63.004963 Exc=      -8.857284 Ees=    -138.991117 Eharris=     -74.487446

 mkekin:
   nout*Vin = smpart,onsite,total=:    -31.889729     17.918650    -13.971079
    E_B(band energy sum)=   -3.100739  E_B-nout*Vin=   10.870340

 m_mkpot_energyterms
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.112434
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.573297   charge     7.783715
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.078880 -7.008142 -0.171751
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -2.642670 -3.104214 -0.140251
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92280   2.91966   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -28.183264        14.361431       -13.821833
   rhoval*ves             -3.810308        -6.733269       -10.543577
   psnuc*ves              14.956903      -283.264774      -268.307871
   utot                    5.573297      -144.999021      -139.425724
   rho*exc                -7.721550        -1.156400        -8.877950
   rho*vxc               -10.112356        -1.573660       -11.686016
   valence chg             6.783715        -3.783715         3.000000
   valence mag             2.239248        -1.239248         1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       10.870340 Ekcore=        62.991040 Ektot    =       73.861380
 Exc=      -8.877950 Ees   =      -139.425724 EKohnSham=      -74.442294
 Magnetic moment=     1.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      148284  -1.1859968212113441E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.11859968213113441E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      148284  -1.1859824736142016E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.11859824737142016E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      148284  -1.1859896474127728E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.11859896475127728E-3

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
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 wgtsmooth=   2.8284271247461901E-003
 mixrho: sought 2 iter from file mixm ; read 3 RMS DQ= 8.68E-1  last it= 1.21E+0
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 4.34D-01
   tj:-1.89092  -0.24126
 mixrho: warning. negative smrho; isp number min=       2   75363 -0.27386D-05
 mixrho: warning. negative smrho; isp number min=       2   75363 -0.27386D-05
 mixrho: warning. negative smrho; isp number min=       2   75363 -0.27386D-05
 mixrho: add corrections to qcell smrho =  0.21987D-07  0.10994D-10
 mixrho: warning. negative smrho; isp number min=       2   75363 -0.27386D-05

 iors  : write rst restart file (binary mesh density)

   it  4  of 10    ehf=     -74.487446   ehk=     -74.442294
 From last iter    ehf=     -74.494427   ehk=     -74.441592
 diffe(q)=  0.006981 (0.867844)    tol= 0.000010 (0.000500)   more=T
i mmom= 1.0000 ehf(eV)=-1013.461294 ehk(eV)=-1012.846965 sev(eV)=-42.188039

--- BNDFP:  begin iteration 5 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111951
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.681114   charge     7.656800
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       75363  -2.7385995448656113E-006
  smvxcm: enforce positive smrho_w. Add srshift= 0.27385995548656111E-5
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       864
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       864
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       864
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.533428 -7.691917 -0.182120
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -2.073162 -2.276635 -0.131597
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       864
 vxcnsp (warning): negative rho: min val =  -2.47E-02
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92280   2.91966   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
     potential shift to crystal energy zero:    0.000008
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -27.576036        13.721890       -13.854146
   rhoval*ves             -3.713846        -6.795382       -10.509228
   psnuc*ves              15.076074      -283.360795      -268.284721
   utot                    5.681114      -145.078088      -139.396974
   rho*exc                -7.606590        -1.318406        -8.924996
   rho*vxc                -9.968553        -1.781954       -11.750506
   valence chg             6.656800        -3.656800         3.000000
   valence mag             3.203559        -1.645424         1.558135
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.5024

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.641524;   3.000000 electrons
         Sum occ. bands:   -2.963194, incl. Bloechl correction: -0.000020
         Mag. moment:       1.000000
 bndfp:Generating TDOS: efermi= -0.641524  dos window emin emax=  -1.274184  2.298398

       contr. to mm extrapolated for r>rmt:   0.030009 est. true mm = 0.994602
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.29856  sum tc=    31.39771  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.22464  sum tc=    31.56334  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905558    6.881738   -3.976180      0.964593    2.105896   -1.141303
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.963194 Vin*nin=     -13.854146 Ek=Eb-Vin*nin=      10.890952
 Ek(core)=      62.997950 Exc=      -8.924996 Ees=    -139.396974 Eharris=     -74.433069

 mkekin:
   nout*Vin = smpart,onsite,total=:    -28.985230     15.325216    -13.660014
    E_B(band energy sum)=   -2.963194  E_B-nout*Vin=   10.696820

 m_mkpot_energyterms
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111738
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.815373   charge     7.976180
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.157014 -7.099561 -0.171472
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -2.827198 -3.356416 -0.142057
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92454   2.92055   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -29.355603        15.676303       -13.679300
   rhoval*ves             -3.623991        -6.796900       -10.420891
   psnuc*ves              15.254736      -283.317041      -268.062305
   utot                    5.815373      -145.056971      -139.241598
   rho*exc                -7.984212        -0.874749        -8.858962
   rho*vxc               -10.455978        -1.205004       -11.660981
   valence chg             6.976180        -3.976180         3.000000
   valence mag             2.141303        -1.141303         1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       10.696820 Ekcore=        62.961046 Ektot    =       73.657866
 Exc=      -8.858962 Ees   =      -139.241598 EKohnSham=      -74.442694
 Magnetic moment=     1.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       75363  -2.7385668432044742E-006
  smvxcm: enforce positive smrho_w. Add srshift= 0.27385668532044740E-5
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       75363  -2.7386322465267484E-006
  smvxcm: enforce positive smrho_w. Add srshift= 0.27386322565267482E-5
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       75363  -2.7385995448656113E-006
  smvxcm: enforce positive smrho_w. Add srshift= 0.27385995548656111E-5

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00   -0.00   -0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force =   0.000000 mRy/au (site 1 )
  Symmetrize forces ...
 wgtsmooth=   2.8284271247461901E-003
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 mixrho: sought 2 iter from file mixm ; read 3 RMS DQ= 4.06E-3  last it= 8.68E-1
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 2.03D-03
   tj:-4.14970   2.99047
 mixrho: warning. negative smrho; isp number min=       2   14753 -0.45008D-06
 mixrho: warning. negative smrho; isp number min=       2   14753 -0.45008D-06
 mixrho: warning. negative smrho; isp number min=       2   14753 -0.45008D-06
 mixrho: add corrections to qcell smrho =  0.42671D-07  0.21336D-10
 mixrho: warning. negative smrho; isp number min=       2   14753 -0.45008D-06

 iors  : write rst restart file (binary mesh density)

   it  5  of 10    ehf=     -74.433069   ehk=     -74.442694
 From last iter    ehf=     -74.487446   ehk=     -74.442294
 diffe(q)=  0.054377 (0.004062)    tol= 0.000010 (0.000500)   more=T
i mmom= 1.0000 ehf(eV)=-1012.721448 ehk(eV)=-1012.852404 sev(eV)=-40.316626

--- BNDFP:  begin iteration 6 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111698
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.794619   charge     8.339203
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       14753  -4.5008034849417714E-007
  smvxcm: enforce positive smrho_w. Add srshift= 0.45008035849417712E-6
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       160
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       160
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       160
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.260969 -7.178972 -0.172136
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -3.348900 -4.096873 -0.143025
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       160
 vxcnsp (warning): negative rho: min val =  -8.46E-04
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92454   2.92055   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
     potential shift to crystal energy zero:    0.000008
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -31.558867        17.872409       -13.686458
   rhoval*ves             -3.645169        -6.780394       -10.425562
   psnuc*ves              15.234407      -283.338713      -268.104306
   utot                    5.794619      -145.059553      -139.264934
   rho*exc                -8.609868        -0.252088        -8.861956
   rho*vxc               -11.275845        -0.389153       -11.664998
   valence chg             7.339203        -4.339203         3.000000
   valence mag             1.835248        -0.815646         1.019602
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:      0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.4340

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.623568;   3.000000 electrons
         Sum occ. bands:   -2.982078, incl. Bloechl correction: -0.000021
         Mag. moment:       1.000000
 bndfp:Generating TDOS: efermi= -0.623568  dos window emin emax=  -1.256305  2.316354

       contr. to mm extrapolated for r>rmt:   0.032902 est. true mm = 0.994404
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.31469  sum tc=    31.42826  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.26617  sum tc=    31.53657  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905891    6.938226   -4.032335      0.961502    1.954184   -0.992681
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.982078 Vin*nin=     -13.686458 Ek=Eb-Vin*nin=      10.704380
 Ek(core)=      62.979509 Exc=      -8.861956 Ees=    -139.264934 Eharris=     -74.443002

 mkekin:
   nout*Vin = smpart,onsite,total=:    -29.841405     16.125220    -13.716185
    E_B(band energy sum)=   -2.982078  E_B-nout*Vin=   10.734107

 m_mkpot_energyterms
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111800
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.797310   charge     8.032335
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.094275 -6.988416 -0.172086
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -2.976858 -3.580726 -0.141092
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92430   2.92159   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -29.674691        15.963307       -13.711384
   rhoval*ves             -3.637995        -6.811908       -10.449903
   psnuc*ves              15.232616      -283.342527      -268.109912
   utot                    5.797310      -145.077218      -139.279908
   rho*exc                -8.071133        -0.790792        -8.861925
   rho*vxc               -10.569142        -1.095715       -11.664857
   valence chg             7.032335        -4.032335         3.000000
   valence mag             1.992681        -0.992681         1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       10.734107 Ekcore=        62.964833 Ektot    =       73.698939
 Exc=      -8.861925 Ees   =      -139.279908 EKohnSham=      -74.442893
 Magnetic moment=     1.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       14753  -4.5007987375582294E-007
  smvxcm: enforce positive smrho_w. Add srshift= 0.45007988375582292E-6
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       14753  -4.5008082323253128E-007
  smvxcm: enforce positive smrho_w. Add srshift= 0.45008083323253126E-6
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       14753  -4.5008034849417708E-007
  smvxcm: enforce positive smrho_w. Add srshift= 0.45008035849417707E-6
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
 mixrho: sought 2 iter from file mixm ; read 3 RMS DQ= 2.55E-3  last it= 4.06E-3
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 1.27D-03
   tj:-0.12269  -0.00374
 mixrho: add corrections to qcell smrho =  0.43124D-07  0.21562D-10

 iors  : write rst restart file (binary mesh density)

   it  6  of 10    ehf=     -74.443002   ehk=     -74.442893
 From last iter    ehf=     -74.433069   ehk=     -74.442694
 diffe(q)= -0.009933 (0.002550)    tol= 0.000010 (0.000500)   more=T
i mmom= 1.0000 ehf(eV)=-1012.856592 ehk(eV)=-1012.855118 sev(eV)=-40.573560

--- BNDFP:  begin iteration 7 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111746
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.814663   charge     8.064426
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.113231 -7.012553 -0.171597
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -3.008614 -3.623047 -0.142374
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92430   2.92159   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
     potential shift to crystal energy zero:    0.000008
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -29.872962        16.181196       -13.691766
   rhoval*ves             -3.625296        -6.809167       -10.434462
   psnuc*ves              15.254621      -283.355837      -268.101216
   utot                    5.814663      -145.082502      -139.267839
   rho*exc                -8.121845        -0.737187        -8.859032
   rho*vxc               -10.635600        -1.025381       -11.660981
   valence chg             7.064426        -4.064426         3.000000
   valence mag             1.968927        -0.991299         0.977628
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.4276

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.620802;   3.000000 electrons
         Sum occ. bands:   -2.979959, incl. Bloechl correction: -0.000021
         Mag. moment:       1.000000
 bndfp:Generating TDOS: efermi= -0.620802  dos window emin emax=  -1.253561  2.319121

       contr. to mm extrapolated for r>rmt:   0.033360 est. true mm = 0.994357
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.31295  sum tc=    31.43015  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.26691  sum tc=    31.53372  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905712    6.941540   -4.035829      0.960997    1.935879   -0.974882
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.979959 Vin*nin=     -13.691766 Ek=Eb-Vin*nin=      10.711807
 Ek(core)=      62.972170 Exc=      -8.859032 Ees=    -139.267839 Eharris=     -74.442893

 mkekin:
   nout*Vin = smpart,onsite,total=:    -29.732591     16.021014    -13.711577
    E_B(band energy sum)=   -2.979959  E_B-nout*Vin=   10.731618

 m_mkpot_energyterms
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111807
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.798710   charge     8.035829
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.082815 -6.969600 -0.172272
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -2.991819 -3.604024 -0.141042
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92425   2.92170   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -29.694440        15.984641       -13.709799
   rhoval*ves             -3.637563        -6.811416       -10.448979
   psnuc*ves              15.234983      -283.340155      -268.105171
   utot                    5.798710      -145.075785      -139.277075
   rho*exc                -8.074635        -0.786691        -8.861326
   rho*vxc               -10.573624        -1.090438       -11.664062
   valence chg             7.035829        -4.035829         3.000000
   valence mag             1.974882        -0.974882         1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       10.731618 Ekcore=        62.963877 Ektot    =       73.695494
 Exc=      -8.861326 Ees   =      -139.277075 EKohnSham=      -74.442907
 Magnetic moment=     1.000000
  smvxcm: all smrho_w is positive
  smvxcm: all smrho_w is positive
  smvxcm: all smrho_w is positive

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
 Maximum Harris force =   0.000000 mRy/au (site 1 )
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
 mixrho: sought 2 iter from file mixm ; read 3 RMS DQ= 2.35E-4  last it= 2.55E-3
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 1.17D-04
   tj:-0.04684  -0.02066
 mixrho: add corrections to qcell smrho =  0.41386D-07  0.20693D-10

 iors  : write rst restart file (binary mesh density)

   it  7  of 10    ehf=     -74.442893   ehk=     -74.442907
 From last iter    ehf=     -74.443002   ehk=     -74.442893
 diffe(q)=  0.000108 (0.000235)    tol= 0.000010 (0.000500)   more=T
i mmom= 1.0000 ehf(eV)=-1012.855116 ehk(eV)=-1012.855300 sev(eV)=-40.544725

--- BNDFP:  begin iteration 8 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111776
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.808401   charge     8.048600
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.089615 -6.978463 -0.171907
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -3.004265 -3.620357 -0.142041
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92425   2.92170   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
     potential shift to crystal energy zero:    0.000008
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -29.772988        16.073318       -13.699669
   rhoval*ves             -3.630169        -6.811251       -10.441420
   psnuc*ves              15.246971      -283.349345      -268.102374
   utot                    5.808401      -145.080298      -139.271897
   rho*exc                -8.093880        -0.765711        -8.859590
   rho*vxc               -10.598820        -1.062901       -11.661721
   valence chg             7.048600        -4.048600         3.000000
   valence mag             1.960148        -0.978313         0.981835
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:      0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.4367

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.620580;   3.000000 electrons
         Sum occ. bands:   -2.979106, incl. Bloechl correction: -0.000021
         Mag. moment:       1.000000
 bndfp:Generating TDOS: efermi= -0.620580  dos window emin emax=  -1.253338  2.319343

       contr. to mm extrapolated for r>rmt:   0.033442 est. true mm = 0.994346
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.31209  sum tc=    31.42985  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.26616  sum tc=    31.53366  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905645    6.940659   -4.035014      0.960903    1.934902   -0.973999
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.979106 Vin*nin=     -13.699669 Ek=Eb-Vin*nin=      10.720563
 Ek(core)=      62.968023 Exc=      -8.859590 Ees=    -139.271897 Eharris=     -74.442901

 mkekin:
   nout*Vin = smpart,onsite,total=:    -29.708991     15.999862    -13.709128
    E_B(band energy sum)=   -2.979106  E_B-nout*Vin=   10.730022

 m_mkpot_energyterms
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111809
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.799464   charge     8.035014
  smvxcm: all smrho_w is positive
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.080959 -6.966873 -0.172322
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -2.991790 -3.604267 -0.141046
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92424   2.92171   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -29.689605        15.980996       -13.708609
   rhoval*ves             -3.637179        -6.810890       -10.448070
   psnuc*ves              15.236108      -283.338776      -268.102668
   utot                    5.799464      -145.074833      -139.275369
   rho*exc                -8.072749        -0.788320        -8.861069
   rho*vxc               -10.571140        -1.092583       -11.663723
   valence chg             7.035014        -4.035014         3.000000
   valence mag             1.973999        -0.973999         1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       10.730022 Ekcore=        62.963506 Ektot    =       73.693528
 Exc=      -8.861069 Ees   =      -139.275369 EKohnSham=      -74.442910
 Magnetic moment=     1.000000
  smvxcm: all smrho_w is positive
  smvxcm: all smrho_w is positive
  smvxcm: all smrho_w is positive

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
 Maximum Harris force =   0.000000 mRy/au (site 1 )
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
 mixrho: sought 2 iter from file mixm ; read 3 RMS DQ= 1.10E-4  last it= 2.35E-4
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 5.51D-05
   tj:-0.72282  -0.00411
 mixrho: add corrections to qcell smrho =  0.39452D-07  0.19726D-10

 iors  : write rst restart file (binary mesh density)

   it  8  of 10    ehf=     -74.442901   ehk=     -74.442910
 From last iter    ehf=     -74.442893   ehk=     -74.442907
 diffe(q)= -0.000008 (0.000110)    tol= 0.000010 (0.000500)   more=F
c mmom= 1.0000 ehf(eV)=-1012.855221 ehk(eV)=-1012.855345 sev(eV)=-40.533124
Exit 0 procid= 0 OK! end of LMF ======================
