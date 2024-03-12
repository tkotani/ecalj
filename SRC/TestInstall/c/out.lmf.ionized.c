INFO: Ubuntu 20.04.6 LTS \n \l
INFO: GNU Fortran (Ubuntu 9.4.0-1ubuntu1~20.04.2) 9.4.0
INFO: 
INFO: MATH: 
INFO: git: commit afd2faa4e42483b6706cdd3e22e28d8c17a2bb6e
INFO:    : Author: Takao Kotani <takaokotani@gmail.com>
INFO:    : Date:   Tue Mar 12 12:07:47 2024 +0900
INFO: linked at Tue Mar 12 14:41:02 JST 2024
=== START LFMA ===
 mpisize=           1
cmdl for python=/home/takao/ecalj/SRC/TestInstall/bin/ctrl2ctrlp.py  c -vzbak=1<ctrl.c >ctrlp.c
m_lmfinit: LMFA
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
rval2: HAM_FTMESH              defa n= 3 val= 25.00000000  25.00000000  25.00000000
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
rval2: HAM_READPSKIPF          defa n= 1 val= 1.00000000
rval2: HAM_V0FIX               defa n= 1 val= 0.00000000
rval2: HAM_PNUFIX              defa n= 1 val= 0.00000000
=== SPEC =1
rval2: SPEC_ATOM@1             val=  C
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
rval2: SPEC_C-HOLE@1           val= 
rval2: SPEC_C-HQ@1             defa n= 2 val= -1.00000000  0.00000000
rval2: SPEC_EREF1              defa n= 1 val= 0.00000000
=== SITE =1
rval2: SITE_ATOM@1             val=  C
rval2: SITE_POS@1              ---- n= 3 val= 0.00000000  0.00000000  0.00000000
rval2: SITE_RELAX@1            defa n= 3 val= 1.00000000  1.00000000  1.00000000
rval2: SITE_AF@1               defa n= 1 val= 0.00000000
rval2: STR_RMAXS               ---- n= 0 val= 
rval2: STR_RMAX                ---- n= 0 val= 
rval2: STR_MXNBR               defa n= 1 val= 0.00000000
rval2: BZ_NKABC                ---- n= 3 val= 2.00000000  2.00000000  2.00000000
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
rval2: SYMGRP                  val=   find
rval2: SYMGRPAF                val= 
rval2: EWALD_AS                defa n= 1 val= 2.00000000
rval2: EWALD_TOL               defa n= 1 val= 0.00000001
rval2: EWALD_NKDMX             defa n= 1 val= 300.00000000
rval2: ITER_NIT                defa n= 1 val= 10.00000000
rval2: ITER_NRMIX              defa n= 1 val= 80.00000000
rval2: ITER_MIX                val=  A2
rval2: ITER_CONV               defa n= 1 val= 0.00001000
rval2: ITER_CONVC              defa n= 1 val= 0.00050000
rval2: ITER_UMIX               defa n= 1 val= 0.50000000
rval2: ITER_TOLU               defa n= 1 val= 0.00000000
rval2: ITER_b                  defa n= 1 val= 0.50000000
rval2: ITER_wc                 defa n= 1 val= -1.00000000
rval2: ITER_w                  defa n= 2 val= 1.00000000  1.00000000
rval2: ITER_k                  defa n= 1 val= -1.00000000
rval2: DYN_MODE                defa n= 1 val= 0.00000000
rval2: DYN_NIT                 defa n= 1 val= 1.00000000
rval2: DYN_HESS                defa n= 1 val= 1.00000000
rval2: DYN_XTOL                defa n= 1 val= 0.00100000
rval2: DYN_GTOL                defa n= 1 val= 0.00000000
rval2: DYN_STEP                defa n= 1 val= 0.01500000
rval2: DYN_NKILL               defa n= 1 val= 0.00000000
mixing param: A/B nmix wt= 0 2 1.000000  1.000000 beta wc killj=  0.500000 -1.000000 -1
 ===> for --jobgw, pwmode is switched to be  0
  bndfp (warning): no sigm file found ... LDA calculation only

pnu === pnu setting ===
pnu   ibas isp   pnu(0:lmxa)  pz(0:lmxa)
pnu:  1 1 2.900  2.850  3.180  4.120   0.000  0.000  0.000  0.000 pnu:  1 2 2.900  2.850  3.180  4.120   0.000  0.000  0.000  0.000

mto === MTO setting ===
mto ispec lmxb lpz nkapii nkaphh=    1    1    0    2    2
mto rsmh1    1  1.30  1.10
mto   eh1    1 -0.70 -0.20
mto rsmh2    1  0.80  0.80
mto  eh2     1 -1.50 -1.00
mto lh       1  1
freats:

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
 vsum=  -59.746958847471802                1
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
INFO: Ubuntu 20.04.6 LTS \n \l
INFO: GNU Fortran (Ubuntu 9.4.0-1ubuntu1~20.04.2) 9.4.0
INFO: 
INFO: MATH: 
INFO: git: commit afd2faa4e42483b6706cdd3e22e28d8c17a2bb6e
INFO:    : Author: Takao Kotani <takaokotani@gmail.com>
INFO:    : Date:   Tue Mar 12 12:07:47 2024 +0900
INFO: linked at Tue Mar 12 14:41:02 JST 2024
===START LMF with   c -vzbak=1 ===
mpisize=4
cmdl for python=/home/takao/ecalj/SRC/TestInstall/bin/ctrl2ctrlp.py  c -vzbak=1<ctrl.c >ctrlp.c
m_lmfinit: LMF
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
rval2: HAM_FTMESH              defa n= 3 val= 25.00000000  25.00000000  25.00000000
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
rval2: HAM_READPSKIPF          defa n= 1 val= 1.00000000
rval2: HAM_V0FIX               defa n= 1 val= 0.00000000
rval2: HAM_PNUFIX              defa n= 1 val= 0.00000000
=== SPEC =1
rval2: SPEC_ATOM@1             val=  C
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
rval2: SPEC_C-HOLE@1           val= 
rval2: SPEC_C-HQ@1             defa n= 2 val= -1.00000000  0.00000000
rval2: SPEC_EREF1              defa n= 1 val= 0.00000000
=== SITE =1
rval2: SITE_ATOM@1             val=  C
rval2: SITE_POS@1              ---- n= 3 val= 0.00000000  0.00000000  0.00000000
rval2: SITE_RELAX@1            defa n= 3 val= 1.00000000  1.00000000  1.00000000
rval2: SITE_AF@1               defa n= 1 val= 0.00000000
rval2: STR_RMAXS               ---- n= 0 val= 
rval2: STR_RMAX                ---- n= 0 val= 
rval2: STR_MXNBR               defa n= 1 val= 0.00000000
rval2: BZ_NKABC                ---- n= 3 val= 2.00000000  2.00000000  2.00000000
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
rval2: SYMGRP                  val=   find
rval2: SYMGRPAF                val= 
rval2: EWALD_AS                defa n= 1 val= 2.00000000
rval2: EWALD_TOL               defa n= 1 val= 0.00000001
rval2: EWALD_NKDMX             defa n= 1 val= 300.00000000
rval2: ITER_NIT                defa n= 1 val= 10.00000000
rval2: ITER_NRMIX              defa n= 1 val= 80.00000000
rval2: ITER_MIX                val=  A2
rval2: ITER_CONV               defa n= 1 val= 0.00001000
rval2: ITER_CONVC              defa n= 1 val= 0.00050000
rval2: ITER_UMIX               defa n= 1 val= 0.50000000
rval2: ITER_TOLU               defa n= 1 val= 0.00000000
rval2: ITER_b                  defa n= 1 val= 0.50000000
rval2: ITER_wc                 defa n= 1 val= -1.00000000
rval2: ITER_w                  defa n= 2 val= 1.00000000  1.00000000
rval2: ITER_k                  defa n= 1 val= -1.00000000
rval2: DYN_MODE                defa n= 1 val= 0.00000000
rval2: DYN_NIT                 defa n= 1 val= 1.00000000
rval2: DYN_HESS                defa n= 1 val= 1.00000000
rval2: DYN_XTOL                defa n= 1 val= 0.00100000
rval2: DYN_GTOL                defa n= 1 val= 0.00000000
rval2: DYN_STEP                defa n= 1 val= 0.01500000
rval2: DYN_NKILL               defa n= 1 val= 0.00000000
mixing param: A/B nmix wt= 0 2 1.000000  1.000000 beta wc killj=  0.500000 -1.000000 -1
 ===> for --jobgw, pwmode is switched to be  0
  bndfp (warning): no sigm file found ... LDA calculation only

pnu === pnu setting ===
pnu   ibas isp   pnu(0:lmxa)  pz(0:lmxa)
pnu:  1 1 2.900  2.850  3.180  4.120   0.000  0.000  0.000  0.000 pnu:  1 2 2.900  2.850  3.180  4.120   0.000  0.000  0.000  0.000

mto === MTO setting ===
mto ispec lmxb lpz nkapii nkaphh=    1    1    0    2    2
mto rsmh1    1  1.30  1.10
mto   eh1    1 -0.70 -0.20
mto rsmh2    1  0.80  0.80
mto  eh2     1 -1.50 -1.00
mto lh       1  1

                Plat                                  Qlat
   1.000000   1.000000   0.000000        0.500000   0.500000  -0.500000
   1.000000   0.000000   1.000000        0.500000  -0.500000   0.500000
   0.000000   1.000000   1.000000       -0.500000   0.500000   0.500000
  Cell vol=   1000.000000

m_lattic_init:  as= 2.000  tol= 1.00E-08  alat= 7.93701   awald= 0.200
         r1=  3.459   nkd=  79      q1=  2.571   nkq= 137
SpaceGroupSym of Lattice: ========start========================== 
 SYMGRP = find
  Generators except find: 
 sgroup:  1 symmetry operations from 0 generators
 symlat: Bravais system is cubic        with 48 symmetry operations.
 symcry: crystal invariant under 48 following symmetry operations for tol=  0.000100
   Enlarging ngen= 1  ng nggen= 48 6
   Enlarging ngen= 2  ng nggen= 48 48
 groupg: the following are sufficient to generate the space group:
  Generators:  trans(cart)= i*r3(-1,1,1) r4z
  Generators:: trans(frac)= i*r3(-1,1,1) r4z
 gensym: ig group ops (:vector means translation in cartesian)
    1  e
    2  i*r3(-1,1,1)
    3  r3(1,-1,-1)
    4  i
    5  r3(-1,1,1)
    6  i*r3(1,-1,-1)
    7  r4z
    8  m(0,1,1)
    9  r4(0,-1,0)
   10  i*r4z
   11  r2(0,1,1)
   12  i*r4(0,-1,0)
   13  r2z
   14  i*r3(-1,-1,-1)
   15  r3(-1,-1,1)
   16  mz
   17  r3(-1,-1,-1)
   18  i*r3(-1,-1,1)
   19  r4(0,0,-1)
   20  i*r4(-1,0,0)
   21  r2(1,0,-1)
   22  i*r4(0,0,-1)
   23  r4(-1,0,0)
   24  m(1,0,-1)
   25  i*r3(1,1,-1)
   26  r3(-1,1,-1)
   27  mx
   28  r3(1,1,-1)
   29  i*r3(-1,1,-1)
   30  r2x
   31  i*r4y
   32  r2(1,-1,0)
   33  i*r4x
   34  r4y
   35  m(1,-1,0)
   36  r4x
   37  r3d
   38  my
   39  r3(1,-1,1)
   40  i*r3d
   41  r2y
   42  i*r3(1,-1,1)
   43  r2(0,1,-1)
   44  m(0,1,-1)
   45  m(1,1,0)
   46  r2(1,1,0)
   47  r2(1,0,1)
   48  m(1,0,1)
 gensym: site permutation table for group operations ...
  ib/ig:  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48
      1:  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1
 splcls:  ibas iclass ispec label(ispec)
            1     1     1     C
SpaceGroupSym of Lattice: ========end =========================== 

 BZMESH: ngrp nq  48    3 QP from    2   2   2 shift=FFF
 TETIRR: sorting       48 tetrahedra ...

sugcut:  make orbital-dependent reciprocal vector cutoffs for tol= 1.00E-06
 spec      l    rsm    eh     gmax    last term    cutoff
  C        0    1.30  -0.70   5.718    1.05E-06    3143 
  C        1    1.10  -0.20   7.226    2.65E-06    5817*
  C        0    0.80  -1.50   9.292    4.01E-04    5817*
  C        1    0.80  -1.00  10.038    2.81E-03    5817*
m_qplistinit:start

iors: read rst restart file (binary mesh density)
 iors  : empty file ... nothing read

rdovfa: read and overlap free-atom densities (mesh density) ...
 rdovfa: expected C,       read C        with rmt=  3.0000  mesh   369  0.020
  ovlpfa: overlap smooth part of FA densities
 rrrrrrdensity111222xxx  -1.9311462682337907E-004   1.9559763616178717E-004   1.9990094016687600E-006   11.773298852858490     
 rrrrrrdensity111222xxx  -1.9311462682337907E-004   1.9559763616178717E-004   1.9990094016687600E-006   11.773298852858490     
 rrrrrrdensity111222xxx  -1.9311462682337907E-004   1.9559763616178717E-004   1.9990094016687600E-006   11.773298852858490     

 Free atom and overlapped crystal site charges:
   ib    true(FA)    smooth(FA)  true(OV)    smooth(OV)    local
 rrrrrrdensity111222xxx  -1.9311462682337907E-004   1.9559763616178717E-004   1.9990094016687600E-006   11.773298852858490     
    1    3.701869    2.363587    3.701843    2.363561    1.338282
 amom    1.810990    1.143831    1.810990    1.143831   -0.667159
 Uniform density added to neutralize background q=  1.000000

 Smooth charge on mesh:            1.661718    moment    1.332841
 Sum of local charges:             1.338282    moments  -0.667159
 Total valence charge:             3.000000    moment    0.665681
 Sum of core charges:              2.000000    moment    0.000000
 Sum of nuclear charges:          -6.000000
 Homogeneous background:           1.000000
 Deviation from neutrality:       -0.000000

 Basis, after reading restart file
 site spec        pos (Cartesian coordinates)         pos (multiples of plat)
   1         C  0.000000   0.000000   0.000000    0.000000   0.000000   0.000000

--- BNDFP:  begin iteration 1 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode

 site class  ilm      vval      ves(rmax)
   1     1     1    0.023423    0.006607

 site class  ilm      vval      ves(rmax)
   1     1     1    0.023423    0.006607
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.006607
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159

 site class  ilm      vval      ves(rmax)
   1     1     1    0.023423    0.006607

 site class  ilm      vval      ves(rmax)
   1     1     1    0.023423    0.006607
   smooth rhoves      2.063910   charge     2.661718
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 25082 -0.49988145199391450E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -1.109295 -1.522745 -0.218019
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -0.385560 -0.423747 -0.164195
  locpot:
   site  1  z=  6.000  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -7.369404  rep=  -5.466799  q =   3.699868
 spin 2:           -5.086143        -3.998153        1.888878
  total:          -12.455547        -9.464952        5.588746
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -7.369404  rep=  -5.466799  q =   3.699868
 spin 2:           -5.086143        -3.998153        1.888878
  total:          -12.455547        -9.464952        5.588746
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -1.402026  rep=  -1.020601  q =   1.697148
 spin 2:           -0.375091        -0.343662        0.553317
  total:           -1.777117        -1.364263        2.250464

 local terms:     true           smooth         local
 rhoeps:        -9.464952      -1.364263      -8.100688
 rhomu:         -7.369404      -1.402026      -5.967378
 spin2:         -5.086143      -0.375091      -4.711052
 total:        -12.455547      -1.777117     -10.678430
 val*vef       -14.138029      -3.734429     -10.403600
 val chg:        3.588746       2.250464       1.338282
 val mmom:        0.667159  core mmom:  -0.000000
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.90000   2.90000   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.18000   3.18000   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.12000   4.12000   0.00000   0.00000
     potential shift to crystal energy zero:    0.000000
 mkpot:
   Energy terms(Ry):       smooth           local           total
   rhoval*veff             -3.733888       -10.403600       -14.137487
   Eestatic                 2.063910      -142.009174      -139.945263
   rho*exc                 -1.494856        -8.100688        -9.595544
   rho*vxc                 -1.946492       -10.678430       -12.624921
   valence chg              1.661718         1.338282         3.000000
   valence mag             1.332841         0.667159         2.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     3 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.0485
 bzwts: --- Tetrahedron Integration ---
 ... only filled or empty bands encountered: ev= -0.824735 ec= -0.768586
 VBmax= -0.824735 CBmin= -0.768586 gap =  0.056150 Ry =   0.763959 eV
 BZINTS: Fermi energy:     -0.824735;   3.000000 electrons
         Sum occ. bands:   -3.222313, incl. Bloechl correction:  0.000000
         Mag. moment:      -1.000000
 bndfp:Generating TDOS: efermi(eV)= -11.221180 DOSwindow emin emax(eV)=  -19.369951  28.778820


 m_bandcal_2nd: to fill eigenfunctions**2 up to Efermi
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85491  sum tc=    31.38764  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78517  sum tc=    31.54593  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.124682  364.161358 -362.036676     -0.138457 -313.300209  313.161752
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -3.222313 Vin*nin=     -14.137487 Ek=Eb-Vin*nin=      10.915175
 Ek(core)=      62.932005 Exc=      -9.595544 Ees=    -139.945263 Eharris=     -75.693628

 mkekin:
   nout*Vin = smpart,onsite,total=:  -1165.826635   1154.205309    -11.621326
    E_B(band energy sum)=   -3.222313  E_B-nout*Vin=    8.399014

 m_mkpot_energyterms
 Energy for background charge q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.736142   -0.207662

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.736142   -0.207662

 site class  ilm      vval      ves(rmax)

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.736142   -0.207662
   1     1     1   -0.736142   -0.207662
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.207662
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves     36.495506   charge   366.036673
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 16 -0.16773259976949416E-6
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -152.338266 -96.188305 -0.226907
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -2141.658237 -2948.665992 -0.408385
  locpot:
   site  1  z=  6.000  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -5.212237  rep=  -3.961025  q =   1.993112
 spin 2:           -5.289593        -4.012799        2.131570
  total:          -10.501830        -7.973824        4.124682
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 dfrce job=          12
 vxcnsp: loc rmu=  -5.212237  rep=  -3.961025  q =   1.993112
 spin 2:           -5.289593        -4.012799        2.131570
  total:          -10.501830        -7.973824        4.124682
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 dfrce job=          12
 dfrce job=          12
 vxcnsp: loc rmu= -96.186912  rep=-152.336386  q =  25.430575
 spin 2:         ***********      ***********      338.730784
  total:         ***********      ***********      364.161358

 local terms:     true           smooth         local
 rhoeps:        -7.973824   -2293.754901    2285.781077
 rhomu:         -5.212237     -96.186912      90.974675
 spin2:         -5.289593   -2948.353518    2943.063925
 total:        -10.501830   -3044.540430    3034.038600
 val*vef       -11.161621   -2516.898414    2505.736794
 val chg:        2.124682     364.161358    -362.036676
 val mmom:      313.161752  core mmom:  -0.000000
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.89276   2.76013   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.12000   4.12000   0.00000   0.00000
 mkpot:
   Energy terms(Ry):       smooth           local           total
   rhoval*veff          -2517.210594      2505.736794       -11.473801
   Eestatic                36.495506      -172.921564      -136.426058
   rho*exc              -2293.996503      2285.781077        -8.215426
   rho*vxc              -3044.854297      3034.038600       -10.815697
   valence chg            365.036673      -362.036676         2.999997
   valence mag          -314.161749       313.161752        -0.999997
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=        8.399014 Ekcore=        62.933572 Ektot    =       71.332586
 Exc=      -8.215426 Ees   =      -136.426058 EKohnSham=      -73.308898
 Magnetic moment=    -0.999997
 dfrce job=          12
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 25082 -0.49988446982455385E-3
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 25082 -0.49987843416327514E-3
 Maximum Harris force =   0.000000 mRy/au (site 1 )
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 25082 -0.49988145199391450E-3
 Maximum Harris force =   0.000000 mRy/au (site 1 )

 Harris correction to forces: no shift of atomic density
  ib         delta-n dVes             delta-n dVxc               total
 Maximum Harris force =   0.000000 mRy/au (site 1 )
   1    0.00    0.00    0.00     0.00    0.00    0.00    -0.00   -0.00   -0.00
 shift forces to make zero average correction:           -0.00   -0.00   -0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force =   0.000000 mRy/au (site 1 )
  Symmetrize forces ...
 wgtsmooth=   8.0000000000000002E-003
 mixrho: sought 2 iter from file mixm ; read 0 RMS DQ= 9.54E+0
 mmom         1.332841   -314.161749
 AMIX: nmix=0 mmix=8  nelts= 33302  beta=0.50000  tm= 5.00000  rmsdel= 4.77D+00
 mixrho: qcell,correction,qmx,summ =  0.17512D-05  0.87558D-09 -0.18035D+03  0.28648D+04
 add q=  0.000002 to preserve neutrality
 mixrho: warning. negative smrho; isp number min=   18946 -0.24843D-03

iors: write rst restart file (binary mesh density)

   it  1  of 10    ehf=     -75.693628   ehk=     -73.308898
h mmom=-1.0000 ehf(eV)=-1029.872362 ehk(eV)=-997.426200 sev(eV)=-43.842140

--- BNDFP:  begin iteration 2 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.356360   -0.100527

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.356360   -0.100527

 site class  ilm      vval      ves(rmax)

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.356360   -0.100527
   1     1     1   -0.356360   -0.100527
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.100527
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      8.326727   charge   184.349197
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 18946 -0.24843197177024069E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -62.642520 -40.704435 -0.244834
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -854.219865 -1175.218090 -0.340641
  locpot:
   site  1  z=  6.000  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.225159  rep=  -4.661838  q =   2.846490
 spin 2:           -5.189235        -4.009699        2.010224
  total:          -11.414393        -8.671537        4.856714
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.225159  rep=  -4.661838  q =   2.846490
 spin 2:           -5.189235        -4.009699        2.010224
  total:          -11.414393        -8.671537        4.856714
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu= -40.638711  rep= -62.581851  q =  13.563861
 spin 2:         ***********      -854.083105      169.642050
  total:         ***********      -916.664955      183.205911

 local terms:     true           smooth         local
 rhoeps:        -8.671537    -916.664955     907.993419
 rhomu:         -6.225159     -40.638711      34.413552
 spin2:         -5.189235   -1175.027305    1169.838070
 total:        -11.414393   -1215.666016    1204.251623
 val*vef       -12.881477   -1296.204636    1283.323159
 val chg:        2.856714     183.205911    -180.349197
 val mmom:      156.914456  core mmom:  -0.000000
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.89276   2.76013   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.12000   4.12000   0.00000   0.00000
     potential shift to crystal energy zero:   -0.000000
 mkpot:
   Energy terms(Ry):       smooth           local           total
   rhoval*veff          -1296.369828      1283.323159       -13.046669
   Eestatic                 8.326727      -146.683663      -138.356935
   rho*exc               -916.862386       907.993419        -8.868967
   rho*vxc              -1215.922525      1204.251623       -11.670902
   valence chg            183.349197      -180.349197         3.000000
   valence mag          -156.414454       156.914456         0.500001
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:      0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     3 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.0472
 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.722129;   3.000000 electrons
         Sum occ. bands:   -3.290097, incl. Bloechl correction: -0.000017
         Mag. moment:       1.000000
 bndfp:Generating TDOS: efermi(eV)= -9.825146 DOSwindow emin emax(eV)=  -18.393768  30.174854


 m_bandcal_2nd: to fill eigenfunctions**2 up to Efermi
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.56622  sum tc=    31.48945  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.53090  sum tc=    31.57762  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.912562    7.084344   -4.171782      0.965736    2.265020   -1.299284
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -3.290097 Vin*nin=     -13.046669 Ek=Eb-Vin*nin=       9.756572
 Ek(core)=      62.932634 Exc=      -8.868967 Ees=    -138.356935 Eharris=     -74.536697

 mkekin:
   nout*Vin = smpart,onsite,total=:    -31.293933     16.890988    -14.402945
    E_B(band energy sum)=   -3.290097  E_B-nout*Vin=   11.112848

 m_mkpot_energyterms
 Energy for background charge q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.399090   -0.112581

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.399090   -0.112581

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.399090   -0.112581

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.399090   -0.112581
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.112581
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.472958   charge     8.171781
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4632 -0.23073295866394718E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.590833 -7.712180 -0.273559
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -2.977431 -3.510969 -0.260394
  locpot:
   site  1  z=  6.000  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.491907  rep=  -4.848543  q =   2.939149
 spin 2:           -5.220226        -4.048897        1.973413
  total:          -11.712133        -8.897440        4.912562
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 dfrce job=          12
 dfrce job=          12
 vxcnsp: loc rmu=  -6.491907  rep=  -4.848543  q =   2.939149
 spin 2:           -5.220226        -4.048897        1.973413
  total:          -11.712133        -8.897440        4.912562
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 dfrce job=          12
 vxcnsp: loc rmu=  -7.626099  rep=  -5.525511  q =   4.674682
 spin 2:           -3.440125        -2.921705        2.409662
  total:          -11.066224        -8.447217        7.084344

 local terms:     true           smooth         local
 rhoeps:        -8.897440      -8.447217      -0.450223
 rhomu:         -6.491907      -7.626099       1.134192
 spin2:         -5.220226      -3.440125      -1.780101
 total:        -11.712133     -11.066224      -0.645909
 val*vef       -13.999074     -30.512129      16.513055
 val chg:        2.912562       7.084344      -4.171782
 val mmom:       -1.299284  core mmom:   0.000000
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92164   2.91793   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.12000   0.00000   0.00000
 mkpot:
   Energy terms(Ry):       smooth           local           total
   rhoval*veff            -30.538542        16.513055       -14.025487
   Eestatic                 5.472958      -145.180775      -139.707818
   rho*exc                 -8.568263        -0.450223        -9.018486
   rho*vxc                -11.223150        -0.645909       -11.869059
   valence chg              7.171781        -4.171782         3.000000
   valence mag             2.299284        -1.299284         1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       11.112848 Ekcore=        63.067070 Ektot    =       74.179919
 Exc=      -9.018486 Ees   =      -139.707818 EKohnSham=      -74.546385
 Magnetic moment=     1.000000
 dfrce job=          12
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 18946 -0.24843346217096729E-3
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 18946 -0.24843048136951408E-3
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 18946 -0.24843197177024069E-3

 Harris correction to forces: no shift of atomic density
  ib         delta-n dVes             delta-n dVxc               total
 Maximum Harris force =   0.000000 mRy/au (site 1 )
   1   -0.00   -0.00   -0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force =   0.000000 mRy/au (site 1 )
  Symmetrize forces ...
 wgtsmooth=   8.0000000000000002E-003
 mixrho: sought 2 iter from file mixm ; read 1 RMS DQ= 4.72E+0  last it= 9.54E+0
 mmom      -156.414454      2.299284
 AMIX: nmix=1 mmix=8  nelts= 33302  beta=0.50000  tm= 5.00000  rmsdel= 2.36D+00
   tj: 0.33113
 mixrho: qcell,correction,qmx,summ =  0.68918D-06  0.34459D-09 -0.12143D+03  0.19442D+04
 mixrho: warning. negative smrho; isp number min=   18672 -0.24244D-03

iors: write rst restart file (binary mesh density)

   it  2  of 10    ehf=     -74.536697   ehk=     -74.546385
 From last iter    ehf=     -75.693628   ehk=     -73.308898
 diffe(q)=  1.156931 (4.723447)    tol= 0.000010 (0.000500)   more=T
i mmom= 1.0000 ehf(eV)=-1014.131387 ehk(eV)=-1014.263211 sev(eV)=-44.764395

--- BNDFP:  begin iteration 3 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.370650   -0.104558

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.370650   -0.104558

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.370650   -0.104558
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.104558

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.370650   -0.104558
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      3.943208   charge   125.428930
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 18672 -0.24243900358088253E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -41.770798 -28.577939 -0.267592
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -500.214588 -689.713061 -0.331535
  locpot:
   site  1  z=  6.000  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.314190  rep=  -4.724076  q =   2.877479
 spin 2:           -5.199434        -4.022743        1.997913
  total:          -11.513624        -8.746819        4.875392
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.314190  rep=  -4.724076  q =   2.877479
 spin 2:           -5.199434        -4.022743        1.997913
  total:          -11.513624        -8.746819        4.875392
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu= -28.506864  rep= -41.708616  q =  10.590990
 spin 2:         -689.563070      -500.106566      113.713333
  total:         -718.069933      -541.815182      124.304322

 local terms:     true           smooth         local
 rhoeps:        -8.746819    -541.815182     533.068362
 rhomu:         -6.314190     -28.506864      22.192673
 spin2:         -5.199434    -689.563070     684.363636
 total:        -11.513624    -718.069933     706.556309
 val*vef       -13.262149    -909.650195     896.388047
 val chg:        2.875392     124.304322    -121.428930
 val mmom:      104.001909  core mmom:  -0.000000
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92164   2.91793   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.12000   0.00000   0.00000
     potential shift to crystal energy zero:    0.000000
 mkpot:
   Energy terms(Ry):       smooth           local           total
   rhoval*veff           -909.749873       896.388047       -13.361826
   Eestatic                 3.943208      -142.772702      -138.829494
   rho*exc               -541.985386       533.068362        -8.917024
   rho*vxc               -718.291000       706.556309       -11.734691
   valence chg            124.428930      -121.428930         3.000000
   valence mag          -103.334689       104.001909         0.667219
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:      0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     3 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.0477
 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.681199;   3.000000 electrons
         Sum occ. bands:   -3.164408, incl. Bloechl correction: -0.000026
         Mag. moment:       1.000000
 bndfp:Generating TDOS: efermi(eV)= -9.268256 DOSwindow emin emax(eV)=  -17.849743  30.731744


 m_bandcal_2nd: to fill eigenfunctions**2 up to Efermi
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.45276  sum tc=    31.45890  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.41248  sum tc=    31.55368  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.906616    6.733418   -3.826802      0.962634    2.380926   -1.418292
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -3.164408 Vin*nin=     -13.361826 Ek=Eb-Vin*nin=      10.197418
 Ek(core)=      62.999882 Exc=      -8.917024 Ees=    -138.829494 Eharris=     -74.549217

 mkekin:
   nout*Vin = smpart,onsite,total=:    -30.961609     16.861816    -14.099793
    E_B(band energy sum)=   -3.164408  E_B-nout*Vin=   10.935385

 m_mkpot_energyterms
 Energy for background charge q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.400329   -0.112931

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.400329   -0.112931

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.400329   -0.112931

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.400329   -0.112931
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.112931
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.462981   charge     7.826802
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4480 -0.22233047784492541E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.337220 -7.383865 -0.272356
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -2.619758 -3.037668 -0.258167
  locpot:
   site  1  z=  6.000  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.466217  rep=  -4.830275  q =   2.934625
 spin 2:           -5.210940        -4.040646        1.971991
  total:          -11.677156        -8.870921        4.906616
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 dfrce job=          12
 vxcnsp: loc rmu=  -6.466217  rep=  -4.830275  q =   2.934625
 spin 2:           -5.210940        -4.040646        1.971991
  total:          -11.677156        -8.870921        4.906616
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 dfrce job=          12
 dfrce job=          12
 vxcnsp: loc rmu=  -7.299308  rep=  -5.273152  q =   4.557172
 spin 2:           -2.969373        -2.565916        2.176246
  total:          -10.268681        -7.839068        6.733418

 local terms:     true           smooth         local
 rhoeps:        -8.870921      -7.839068      -1.031853
 rhomu:         -6.466217      -7.299308       0.833092
 spin2:         -5.210940      -2.969373      -2.241567
 total:        -11.677156     -10.268681      -1.408475
 val*vef       -13.861896     -28.466862      14.604966
 val chg:        2.906616       6.733418      -3.826802
 val mmom:       -1.418292  core mmom:   0.000000
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92113   2.91761   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.12000   0.00000   0.00000
 mkpot:
   Energy terms(Ry):       smooth           local           total
   rhoval*veff            -28.493784        14.604966       -13.888818
   Eestatic                 5.462981      -144.966025      -139.503044
   rho*exc                 -7.956977        -1.031853        -8.988830
   rho*vxc                -10.421533        -1.408475       -11.830008
   valence chg              6.826802        -3.826802         3.000000
   valence mag             2.418292        -1.418292         1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       10.935385 Ekcore=        63.012586 Ektot    =       73.947971
 Exc=      -8.988830 Ees   =      -139.503044 EKohnSham=      -74.543903
 Magnetic moment=     1.000000
 dfrce job=          12
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 18672 -0.24244022827420104E-3
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 18672 -0.24243777888756401E-3
 Maximum Harris force =   0.000000 mRy/au (site 1 )
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 18672 -0.24243900358088253E-3
 Maximum Harris force =   0.000000 mRy/au (site 1 )

 Harris correction to forces: no shift of atomic density
  ib         delta-n dVes             delta-n dVxc               total
 Maximum Harris force =   0.000000 mRy/au (site 1 )
   1   -0.00   -0.00   -0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force =   0.000000 mRy/au (site 1 )
  Symmetrize forces ...
 wgtsmooth=   8.0000000000000002E-003
 mixrho: sought 2 iter from file mixm ; read 2 RMS DQ= 3.15E+0  last it= 4.72E+0
 mmom      -103.334689      2.418292
 AMIX: nmix=2 mmix=8  nelts= 33302  beta=0.50000  tm= 5.00000  rmsdel= 1.57D+00
   tj:-1.61267   0.04822
 mixrho: qcell,correction,qmx,summ =  0.87950D-07  0.43975D-10 -0.20516D+02  0.36744D+03
 mixrho: warning. negative smrho; isp number min=   12868 -0.22163D-03

iors: write rst restart file (binary mesh density)

   it  3  of 10    ehf=     -74.549217   ehk=     -74.543903
 From last iter    ehf=     -74.536697   ehk=     -74.546385
 diffe(q)= -0.012520 (3.149723)    tol= 0.000010 (0.000500)   more=T
i mmom= 1.0000 ehf(eV)=-1014.301736 ehk(eV)=-1014.229434 sev(eV)=-43.054303

--- BNDFP:  begin iteration 4 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396607   -0.111881

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396607   -0.111881
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111881

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396607   -0.111881

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396607   -0.111881
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      4.410236   charge    24.516239
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 12868 -0.22163487493657004E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -10.300642 -10.074606 -0.273741
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -41.812882 -58.639907 -0.281658
  locpot:
   site  1  z=  6.000  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.433717  rep=  -4.807485  q =   2.924519
 spin 2:           -5.206128        -4.035123        1.975097
  total:          -11.639845        -8.842608        4.899616
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.433717  rep=  -4.807485  q =   2.924519
 spin 2:           -5.206128        -4.035123        1.975097
  total:          -11.639845        -8.842608        4.899616
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -9.992746  rep= -10.237175  q =   5.371000
 spin 2:          -58.558522       -41.750523       18.044856
  total:          -68.551268       -51.987698       23.415855

 local terms:     true           smooth         local
 rhoeps:        -8.842608     -51.987698      43.145090
 rhomu:         -6.433717      -9.992746       3.559028
 spin2:         -5.206128     -58.558522      53.352394
 total:        -11.639845     -68.551268      56.911423
 val*vef       -13.718270    -143.282481     129.564211
 val chg:        2.899616      23.415855     -20.516239
 val mmom:       13.623279  core mmom:  -0.000000
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92113   2.91761   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.12000   0.00000   0.00000
     potential shift to crystal energy zero:   -0.000000
 mkpot:
   Energy terms(Ry):       smooth           local           total
   rhoval*veff           -143.316562       129.564211       -13.752351
   Eestatic                 4.410236      -143.746746      -139.336510
   rho*exc                -52.113524        43.145090        -8.968434
   rho*vxc                -68.714513        56.911423       -11.803090
   valence chg             23.516239       -20.516239         3.000000
   valence mag           -12.670921        13.623279         0.952358
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     3 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.0466
 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.628021;   3.000000 electrons
         Sum occ. bands:   -2.999883, incl. Bloechl correction: -0.000038
         Mag. moment:       1.000000
 bndfp:Generating TDOS: efermi(eV)= -8.544729 DOSwindow emin emax(eV)=  -17.142708  31.455271


 m_bandcal_2nd: to fill eigenfunctions**2 up to Efermi
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.31172  sum tc=    31.42541  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.26560  sum tc=    31.53018  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.901913    7.014394   -4.112480      0.958221    2.176173   -1.217952
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.999883 Vin*nin=     -13.752351 Ek=Eb-Vin*nin=      10.752469
 Ek(core)=      63.006225 Exc=      -8.968434 Ees=    -139.336510 Eharris=     -74.546250

 mkekin:
   nout*Vin = smpart,onsite,total=:    -32.469693     18.753112    -13.716580
    E_B(band energy sum)=   -2.999883  E_B-nout*Vin=   10.716698

 m_mkpot_energyterms
 Energy for background charge q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.398512   -0.112418

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.398512   -0.112418

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.398512   -0.112418

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.398512   -0.112418
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.112418
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.718821   charge     8.112480
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4334 -0.21919044238501460E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.404904 -7.443963 -0.272924
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -2.953053 -3.501804 -0.258797
  locpot:
   site  1  z=  6.000  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.436514  rep=  -4.809174  q =   2.930067
 spin 2:           -5.200513        -4.031341        1.971846
  total:          -11.637027        -8.840516        4.901913
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 dfrce job=          12
 dfrce job=          12
 vxcnsp: loc rmu=  -6.436514  rep=  -4.809174  q =   2.930067
 spin 2:           -5.200513        -4.031341        1.971846
  total:          -11.637027        -8.840516        4.901913
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 dfrce job=          12
 vxcnsp: loc rmu=  -7.358907  rep=  -5.340473  q =   4.595283
 spin 2:           -3.434068        -2.899627        2.419110
  total:          -10.792976        -8.240100        7.014394

 local terms:     true           smooth         local
 rhoeps:        -8.840516      -8.240100      -0.600416
 rhomu:         -6.436514      -7.358907       0.922393
 spin2:         -5.200513      -3.434068      -1.766445
 total:        -11.637027     -10.792976      -0.844051
 val*vef       -13.688900     -30.207278      16.518377
 val chg:        2.901913       7.014394      -4.112480
 val mmom:       -1.217952  core mmom:   0.000000
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92047   2.91740   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.12000   0.00000   0.00000
 mkpot:
   Energy terms(Ry):       smooth           local           total
   rhoval*veff            -30.235040        16.518377       -13.716663
   Eestatic                 5.718821      -144.976366      -139.257545
   rho*exc                 -8.357958        -0.600416        -8.958373
   rho*vxc                -10.945767        -0.844051       -11.789818
   valence chg              7.112480        -4.112480         3.000000
   valence mag             2.217951        -1.217952         1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       10.716698 Ekcore=        62.955595 Ektot    =       73.672293
 Exc=      -8.958373 Ees   =      -139.257545 EKohnSham=      -74.543625
 Magnetic moment=     1.000000
 dfrce job=          12
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 12868 -0.22163561613449751E-3
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 12868 -0.22163413373864257E-3
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 12868 -0.22163487493657004E-3

 Harris correction to forces: no shift of atomic density
  ib         delta-n dVes             delta-n dVxc               total
 Maximum Harris force =   0.000000 mRy/au (site 1 )
   1    0.00   -0.00   -0.00     0.00    0.00    0.00    -0.00    0.00    0.00
 shift forces to make zero average correction:           -0.00    0.00    0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force =   0.000000 mRy/au (site 1 )
  Symmetrize forces ...
 wgtsmooth=   8.0000000000000002E-003
 mixrho: sought 2 iter from file mixm ; read 3 RMS DQ= 4.43E-1  last it= 3.15E+0
 mmom       -12.670921      2.217951
 AMIX: nmix=2 mmix=8  nelts= 33302  beta=0.50000  tm= 5.00000  rmsdel= 2.21D-01
   tj: 0.35488  -0.32783
 mixrho: qcell,correction,qmx,summ =  0.76501D-07  0.38251D-10 -0.39607D+01  0.10876D+03
 mixrho: warning. negative smrho; isp number min=    4616 -0.21839D-03

iors: write rst restart file (binary mesh density)

   it  4  of 10    ehf=     -74.546250   ehk=     -74.543625
 From last iter    ehf=     -74.549217   ehk=     -74.543903
 diffe(q)=  0.002967 (0.442703)    tol= 0.000010 (0.000500)   more=T
i mmom= 1.0000 ehf(eV)=-1014.261371 ehk(eV)=-1014.225652 sev(eV)=-40.815805

--- BNDFP:  begin iteration 5 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.399779   -0.112776

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.399779   -0.112776

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.399779   -0.112776
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.112776

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.399779   -0.112776
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.577045   charge     7.960683
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4616 -0.21839494107975000E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.319817 -7.336852 -0.272290
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -2.812680 -3.313630 -0.258092
  locpot:
   site  1  z=  6.000  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.444450  rep=  -4.814854  q =   2.931057
 spin 2:           -5.203781        -4.034141        1.971531
  total:          -11.648231        -8.848995        4.902588
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.444450  rep=  -4.814854  q =   2.931057
 spin 2:           -5.203781        -4.034141        1.971531
  total:          -11.648231        -8.848995        4.902588
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -7.252499  rep=  -5.255919  q =   4.546656
 spin 2:           -3.246217        -2.759508        2.316614
  total:          -10.498716        -8.015427        6.863270

 local terms:     true           smooth         local
 rhoeps:        -8.848995      -8.015427      -0.833568
 rhomu:         -6.444450      -7.252499       0.808049
 spin2:         -5.203781      -3.246217      -1.957564
 total:        -11.648231     -10.498716      -1.149515
 val*vef       -13.737470     -29.273801      15.536331
 val chg:        2.902588       6.863270      -3.960683
 val mmom:       -1.270516  core mmom:  -0.000000
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92047   2.91740   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.12000   0.00000   0.00000
     potential shift to crystal energy zero:    0.000000
 mkpot:
   Energy terms(Ry):       smooth           local           total
   rhoval*veff            -29.301329        15.536331       -13.764999
   Eestatic                 5.577045      -144.912131      -139.335086
   rho*exc                 -8.132497        -0.833568        -8.966065
   rho*vxc                -10.650483        -1.149515       -11.799997
   valence chg              6.960683        -3.960683         3.000000
   valence mag             2.270249        -1.270516         0.999733
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:      0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     3 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.0447
 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.623985;   3.000000 electrons
         Sum occ. bands:   -2.988141, incl. Bloechl correction: -0.000039
         Mag. moment:       1.000000
 bndfp:Generating TDOS: efermi(eV)= -8.489819 DOSwindow emin emax(eV)=  -17.090741  31.510181


 m_bandcal_2nd: to fill eigenfunctions**2 up to Efermi
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.30359  sum tc=    31.42531  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.25693  sum tc=    31.53090  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905396    7.789190   -4.883794      0.955893    1.666987   -0.711094
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.988141 Vin*nin=     -13.764999 Ek=Eb-Vin*nin=      10.776858
 Ek(core)=      62.980910 Exc=      -8.966065 Ees=    -139.335086 Eharris=     -74.543383

 mkekin:
   nout*Vin = smpart,onsite,total=:    -33.882242     20.151737    -13.730505
    E_B(band energy sum)=   -2.988141  E_B-nout*Vin=   10.742364

 m_mkpot_energyterms
 Energy for background charge q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.395678   -0.111619

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.395678   -0.111619

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.395678   -0.111619

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.395678   -0.111619
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111619
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.934917   charge     8.883793
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4366 -0.22343615376883833E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.739599 -7.797707 -0.274164
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -3.930593 -4.868218 -0.261263
  locpot:
   site  1  z=  6.000  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.436852  rep=  -4.809784  q =   2.930645
 spin 2:           -5.206415        -4.035509        1.974752
  total:          -11.643267        -8.845293        4.905396
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.436852  rep=  -4.809784  q =   2.930645
 spin 2:           -5.206415        -4.035509        1.974752
  total:          -11.643267        -8.845293        4.905396
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 dfrce job=          12
 dfrce job=          12
 dfrce job=          12
 vxcnsp: loc rmu=  -7.711238  rep=  -5.674040  q =   4.728088
 spin 2:           -4.799311        -3.876316        3.061102
  total:          -12.510549        -9.550355        7.789190

 local terms:     true           smooth         local
 rhoeps:        -8.845293      -9.550355       0.705062
 rhomu:         -6.436852      -7.711238       1.274386
 spin2:         -5.206415      -4.799311      -0.407105
 total:        -11.643267     -12.510549       0.867282
 val*vef       -13.706758     -34.974504      21.267745
 val chg:        2.905396       7.789190      -4.883794
 val mmom:       -0.711094  core mmom:   0.000000
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92022   2.91625   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.12000   0.00000   0.00000
 mkpot:
   Energy terms(Ry):       smooth           local           total
   rhoval*veff            -35.002147        21.267745       -13.734402
   Eestatic                 5.934917      -145.214392      -139.279475
   rho*exc                 -9.670193         0.705062        -8.965130
   rho*vxc                -12.665925         0.867282       -11.798643
   valence chg              7.883793        -4.883794         3.000000
   valence mag             1.711093        -0.711094         1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       10.742364 Ekcore=        62.956211 Ektot    =       73.698576
 Exc=      -8.965130 Ees   =      -139.279475 EKohnSham=      -74.546029
 Magnetic moment=     1.000000
 dfrce job=          12
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4616 -0.21839560261813624E-3
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4620 -0.21839427954136375E-3
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4616 -0.21839494107975000E-3
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )

 Harris correction to forces: no shift of atomic density
  ib         delta-n dVes             delta-n dVxc               total
 Maximum Harris force =   0.000000 mRy/au (site 1 )
   1    0.00   -0.00   -0.00     0.00    0.00    0.00    -0.00    0.00    0.00
 shift forces to make zero average correction:           -0.00    0.00    0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force =   0.000000 mRy/au (site 1 )
  Symmetrize forces ...
 wgtsmooth=   8.0000000000000002E-003
 mixrho: sought 2 iter from file mixm ; read 3 RMS DQ= 1.87E-2  last it= 4.43E-1
 mmom         2.270249      1.711093
 AMIX: nmix=2 mmix=8  nelts= 33302  beta=0.50000  tm= 5.00000  rmsdel= 9.37D-03
   tj: 2.26115  -0.32326
 mixrho: qcell,correction,qmx,summ =  0.83646D-07  0.41823D-10 -0.34518D+01  0.10081D+03
 mixrho: warning. negative smrho; isp number min=    4350 -0.21607D-03

iors: write rst restart file (binary mesh density)

   it  5  of 10    ehf=     -74.543383   ehk=     -74.546029
 From last iter    ehf=     -74.546250   ehk=     -74.543625
 diffe(q)=  0.002868 (0.018734)    tol= 0.000010 (0.000500)   more=T
i mmom= 1.0000 ehf(eV)=-1014.222355 ehk(eV)=-1014.258364 sev(eV)=-40.656046

--- BNDFP:  begin iteration 6 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.401303   -0.113206

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.401303   -0.113206
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.113206

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.401303   -0.113206

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.401303   -0.113206
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.472603   charge     7.451782
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4350 -0.21607083878936038E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.133888 -7.128991 -0.271564
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -2.189974 -2.463487 -0.256283
  locpot:
   site  1  z=  6.000  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.444128  rep=  -4.814430  q =   2.930823
 spin 2:           -5.200600        -4.031891        1.970070
  total:          -11.644728        -8.846321        4.900893
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.444128  rep=  -4.814430  q =   2.930823
 spin 2:           -5.200600        -4.031891        1.970070
  total:          -11.644728        -8.846321        4.900893
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -7.045393  rep=  -5.070639  q =   4.469879
 spin 2:           -2.396916        -2.137374        1.882796
  total:           -9.442309        -7.208013        6.352675

 local terms:     true           smooth         local
 rhoeps:        -8.846321      -7.208013      -1.638308
 rhomu:         -6.444128      -7.045393       0.601265
 spin2:         -5.200600      -2.396916      -2.803684
 total:        -11.644728      -9.442309      -2.202419
 val*vef       -13.731069     -26.283641      12.552572
 val chg:        2.900893       6.352675      -3.451782
 val mmom:       -1.626329  core mmom:  -0.000000
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92022   2.91625   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.12000   0.00000   0.00000
     potential shift to crystal energy zero:   -0.000000
 mkpot:
   Energy terms(Ry):       smooth           local           total
   rhoval*veff            -26.311148        12.552572       -13.758576
   Eestatic                 5.472603      -144.788009      -139.315406
   rho*exc                 -7.323861        -1.638308        -8.962169
   rho*vxc                 -9.592478        -2.202419       -11.794897
   valence chg              6.451782        -3.451782         3.000000
   valence mag             2.626380        -1.626329         1.000050
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:      0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     3 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.0468
 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.625216;   3.000000 electrons
         Sum occ. bands:   -2.991780, incl. Bloechl correction: -0.000038
         Mag. moment:       1.000000
 bndfp:Generating TDOS: efermi(eV)= -8.506565 DOSwindow emin emax(eV)=  -17.107953  31.493435


 m_bandcal_2nd: to fill eigenfunctions**2 up to Efermi
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.30652  sum tc=    31.42632  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.25965  sum tc=    31.53192  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.906227    7.935228   -5.029001      0.955634    1.551712   -0.596077
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.991780 Vin*nin=     -13.758576 Ek=Eb-Vin*nin=      10.766796
 Ek(core)=      62.968561 Exc=      -8.962169 Ees=    -139.315406 Eharris=     -74.542218

 mkekin:
   nout*Vin = smpart,onsite,total=:    -34.081228     20.334948    -13.746279
    E_B(band energy sum)=   -2.991780  E_B-nout*Vin=   10.754499

 m_mkpot_energyterms
 Energy for background charge q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.395194   -0.111482

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.395194   -0.111482

 site class  ilm      vval      ves(rmax)

   1     1     1   -0.395194   -0.111482
 site class  ilm      vval      ves(rmax)
   1     1     1   -0.395194   -0.111482
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111482
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.967575   charge     9.029001
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4382 -0.22473503303078676E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.792097 -7.843389 -0.274443
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -4.136092 -5.160958 -0.261821
  locpot:
   site  1  z=  6.000  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.438022  rep=  -4.810679  q =   2.930931
 spin 2:           -5.207805        -4.036564        1.975296
  total:          -11.645828        -8.847243        4.906227
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.438022  rep=  -4.810679  q =   2.930931
 spin 2:           -5.207805        -4.036564        1.975296
  total:          -11.645828        -8.847243        4.906227
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 dfrce job=          12
 dfrce job=          12
 dfrce job=          12
 vxcnsp: loc rmu=  -7.756541  rep=  -5.726236  q =   4.743470
 spin 2:           -5.091665        -4.081528        3.191758
  total:          -12.848206        -9.807764        7.935228

 local terms:     true           smooth         local
 rhoeps:        -8.847243      -9.807764       0.960521
 rhomu:         -6.438022      -7.756541       1.318519
 spin2:         -5.207805      -5.091665      -0.116141
 total:        -11.645828     -12.848206       1.202378
 val*vef       -13.715814     -35.889257      22.173443
 val chg:        2.906227       7.935228      -5.029001
 val mmom:       -0.596077  core mmom:  -0.000000
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92025   2.91602   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.12000   0.00000   0.00000
 mkpot:
   Energy terms(Ry):       smooth           local           total
   rhoval*veff            -35.916874        22.173443       -13.743430
   Eestatic                 5.967575      -145.259338      -139.291763
   rho*exc                 -9.928189         0.960521        -8.967668
   rho*vxc                -13.004347         1.202378       -11.801969
   valence chg              8.029001        -5.029001         3.000000
   valence mag             1.596077        -0.596077         1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       10.754499 Ekcore=        62.958240 Ektot    =       73.712739
 Exc=      -8.967668 Ees   =      -139.291763 EKohnSham=      -74.546691
 Magnetic moment=     1.000000
 dfrce job=          12
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4350 -0.21607149278640725E-3
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4350 -0.21607018479231351E-3
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4350 -0.21607083878936038E-3
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )

 Harris correction to forces: no shift of atomic density
  ib         delta-n dVes             delta-n dVxc               total
 Maximum Harris force =   0.000000 mRy/au (site 1 )
   1    0.00   -0.00   -0.00     0.00    0.00    0.00    -0.00    0.00    0.00
 shift forces to make zero average correction:           -0.00    0.00    0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force =   0.000000 mRy/au (site 1 )
  Symmetrize forces ...
 wgtsmooth=   8.0000000000000002E-003
 mixrho: sought 2 iter from file mixm ; read 3 RMS DQ= 3.34E-2  last it= 1.87E-2
 mmom         2.626380      1.596077
 AMIX: nmix=2 mmix=8  nelts= 33302  beta=0.50000  tm= 5.00000  rmsdel= 1.67D-02
   tj: 3.16515  -0.02781
 mixrho: qcell,correction,qmx,summ =  0.98896D-07  0.49448D-10 -0.45914D+01  0.11862D+03
 mixrho: warning. negative smrho; isp number min=    4362 -0.22203D-03

iors: write rst restart file (binary mesh density)

   it  6  of 10    ehf=     -74.542218   ehk=     -74.546691
 From last iter    ehf=     -74.543383   ehk=     -74.546029
 diffe(q)=  0.001164 (0.033392)    tol= 0.000010 (0.000500)   more=T
i mmom= 1.0000 ehf(eV)=-1014.206513 ehk(eV)=-1014.267374 sev(eV)=-40.705563

--- BNDFP:  begin iteration 7 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode

 site class  ilm      vval      ves(rmax)

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396621   -0.111885
   1     1     1   -0.396621   -0.111885

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396621   -0.111885
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111885

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396621   -0.111885
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.856572   charge     8.591400
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4362 -0.22202510935297488E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.635265 -7.703722 -0.273661
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -3.530805 -4.301217 -0.260256
  locpot:
   site  1  z=  6.000  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.439367  rep=  -4.811520  q =   2.930893
 spin 2:           -5.206443        -4.035688        1.974111
  total:          -11.645809        -8.847208        4.905004
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.439367  rep=  -4.811520  q =   2.930893
 spin 2:           -5.206443        -4.035688        1.974111
  total:          -11.645809        -8.847208        4.905004
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -7.617847  rep=  -5.570179  q =   4.693367
 spin 2:           -4.232870        -3.476940        2.803037
  total:          -11.850718        -9.047119        7.496405

 local terms:     true           smooth         local
 rhoeps:        -8.847208      -9.047119       0.199911
 rhomu:         -6.439367      -7.617847       1.178481
 spin2:         -5.206443      -4.232870      -0.973572
 total:        -11.645809     -11.850718       0.204908
 val*vef       -13.719949     -33.142856      19.422907
 val chg:        2.905004       7.496405      -4.591400
 val mmom:       -0.933548  core mmom:  -0.000000
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92025   2.91602   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.12000   0.00000   0.00000
     potential shift to crystal energy zero:   -0.000000
 mkpot:
   Energy terms(Ry):       smooth           local           total
   rhoval*veff            -33.170380        19.422907       -13.747473
   Eestatic                 5.856572      -145.157062      -139.300490
   rho*exc                 -9.166070         0.199911        -8.966159
   rho*vxc                -12.004939         0.204908       -11.800031
   valence chg              7.591400        -4.591400         3.000000
   valence mag             1.933733        -0.933548         1.000185
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:      0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     3 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.0488
 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.624327;   3.000000 electrons
         Sum occ. bands:   -2.989491, incl. Bloechl correction: -0.000039
         Mag. moment:       1.000000
 bndfp:Generating TDOS: efermi(eV)= -8.494470 DOSwindow emin emax(eV)=  -17.094198  31.505530


 m_bandcal_2nd: to fill eigenfunctions**2 up to Efermi
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.30792  sum tc=    31.42789  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.26182  sum tc=    31.53291  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.904970    7.684598   -4.779628      0.956107    1.742891   -0.786784
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.989491 Vin*nin=     -13.747473 Ek=Eb-Vin*nin=      10.757982
 Ek(core)=      62.963401 Exc=      -8.966159 Ees=    -139.300490 Eharris=     -74.545266

 mkekin:
   nout*Vin = smpart,onsite,total=:    -34.161427     20.431079    -13.730348
    E_B(band energy sum)=   -2.989491  E_B-nout*Vin=   10.740857

 m_mkpot_energyterms
 Energy for background charge q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396059   -0.111726

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396059   -0.111726

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396059   -0.111726
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111726

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396059   -0.111726
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.905432   charge     8.779628
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4360 -0.22231433098108142E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.698143 -7.757858 -0.273886
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -3.790276 -4.669678 -0.260793
  locpot:
   site  1  z=  6.000  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.436906  rep=  -4.809800  q =   2.930539
 spin 2:           -5.206097        -4.035286        1.974431
  total:          -11.643002        -8.845086        4.904970
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 dfrce job=          12
 dfrce job=          12
 vxcnsp: loc rmu=  -6.436906  rep=  -4.809800  q =   2.930539
 spin 2:           -5.206097        -4.035286        1.974431
  total:          -11.643002        -8.845086        4.904970
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 dfrce job=          12
 vxcnsp: loc rmu=  -7.671767  rep=  -5.632880  q =   4.713744
 spin 2:           -4.601137        -3.736274        2.970854
  total:          -12.272904        -9.369154        7.684598

 local terms:     true           smooth         local
 rhoeps:        -8.845086      -9.369154       0.524069
 rhomu:         -6.436906      -7.671767       1.234861
 spin2:         -5.206097      -4.601137      -0.604959
 total:        -11.643002     -12.272904       0.629902
 val*vef       -13.705938     -34.319427      20.613489
 val chg:        2.904970       7.684598      -4.779628
 val mmom:       -0.786784  core mmom:   0.000000
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92024   2.91652   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.12000   0.00000   0.00000
 mkpot:
   Energy terms(Ry):       smooth           local           total
   rhoval*veff            -34.347046        20.613489       -13.733557
   Eestatic                 5.905432      -145.188148      -139.282716
   rho*exc                 -9.488419         0.524069        -8.964350
   rho*vxc                -12.427536         0.629902       -11.797634
   valence chg              7.779628        -4.779628         3.000000
   valence mag             1.786783        -0.786784         1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       10.740857 Ekcore=        62.960792 Ektot    =       73.701649
 Exc=      -8.964350 Ees   =      -139.282716 EKohnSham=      -74.545418
 Magnetic moment=     1.000000
 dfrce job=          12
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4362 -0.22202578128545082E-3
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4362 -0.22202443742049894E-3
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4362 -0.22202510935297488E-3

 Harris correction to forces: no shift of atomic density
  ib         delta-n dVes             delta-n dVxc               total
 Maximum Harris force =   0.000000 mRy/au (site 1 )
   1   -0.00   -0.00    0.00     0.00    0.00    0.00     0.00    0.00   -0.00
 shift forces to make zero average correction:            0.00    0.00   -0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force =   0.000000 mRy/au (site 1 )
  Symmetrize forces ...
 wgtsmooth=   8.0000000000000002E-003
 mixrho: sought 2 iter from file mixm ; read 3 RMS DQ= 4.39E-3  last it= 3.34E-2
 mmom         1.933733      1.786783
 AMIX: nmix=2 mmix=8  nelts= 33302  beta=0.50000  tm= 5.00000  rmsdel= 2.19D-03
   tj:-0.41820   0.54164
 mixrho: qcell,correction,qmx,summ =  0.98708D-07  0.49354D-10 -0.47291D+01  0.12077D+03
 mixrho: warning. negative smrho; isp number min=    4358 -0.22223D-03

iors: write rst restart file (binary mesh density)

   it  7  of 10    ehf=     -74.545266   ehk=     -74.545418
 From last iter    ehf=     -74.542218   ehk=     -74.546691
 diffe(q)= -0.003048 (0.004386)    tol= 0.000010 (0.000500)   more=T
i mmom= 1.0000 ehf(eV)=-1014.247986 ehk(eV)=-1014.250048 sev(eV)=-40.674419

--- BNDFP:  begin iteration 8 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396294   -0.111792

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396294   -0.111792

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396294   -0.111792
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111792

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396294   -0.111792
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.881088   charge     8.729063
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4358 -0.22222922494696091E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.682272 -7.744668 -0.273795
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -3.722217 -4.572857 -0.260633
  locpot:
   site  1  z=  6.000  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.438172  rep=  -4.810692  q =   2.930722
 spin 2:           -5.206382        -4.035569        1.974323
  total:          -11.644554        -8.846260        4.905045
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.438172  rep=  -4.810692  q =   2.930722
 spin 2:           -5.206382        -4.035569        1.974323
  total:          -11.644554        -8.846260        4.905045
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -7.658676  rep=  -5.617086  q =   4.708224
 spin 2:           -4.504363        -3.668251        2.925884
  total:          -12.163039        -9.285337        7.634107

 local terms:     true           smooth         local
 rhoeps:        -8.846260      -9.285337       0.439076
 rhomu:         -6.438172      -7.658676       1.220504
 spin2:         -5.206382      -4.504363      -0.702019
 total:        -11.644554     -12.163039       0.518485
 val*vef       -13.713551     -33.999608      20.286057
 val chg:        2.905045       7.634107      -4.729063
 val mmom:       -0.825942  core mmom:  -0.000000
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92024   2.91652   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.12000   0.00000   0.00000
     potential shift to crystal energy zero:   -0.000000
 mkpot:
   Energy terms(Ry):       smooth           local           total
   rhoval*veff            -34.027188        20.286057       -13.741131
   Eestatic                 5.881088      -145.173378      -139.292290
   rho*exc                 -9.404489         0.439076        -8.965413
   rho*vxc                -12.317525         0.518485       -11.799040
   valence chg              7.729063        -4.729063         3.000000
   valence mag             1.825940        -0.825942         0.999998
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:      0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     3 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.0466
 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.624795;   3.000000 electrons
         Sum occ. bands:   -2.990889, incl. Bloechl correction: -0.000039
         Mag. moment:       1.000000
 bndfp:Generating TDOS: efermi(eV)= -8.500838 DOSwindow emin emax(eV)=  -17.100312  31.499162


 m_bandcal_2nd: to fill eigenfunctions**2 up to Efermi
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.30953  sum tc=    31.42843  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.26350  sum tc=    31.53333  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.904951    7.664125   -4.759174      0.956214    1.760649   -0.804434
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.990889 Vin*nin=     -13.741131 Ek=Eb-Vin*nin=      10.750242
 Ek(core)=      62.962096 Exc=      -8.965413 Ees=    -139.292290 Eharris=     -74.545365

 mkekin:
   nout*Vin = smpart,onsite,total=:    -34.173935     20.440805    -13.733130
    E_B(band energy sum)=   -2.990889  E_B-nout*Vin=   10.742241

 m_mkpot_energyterms
 Energy for background charge q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396128   -0.111746

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396128   -0.111746

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396128   -0.111746
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111746

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396128   -0.111746
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.898853   charge     8.759174
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4362 -0.22260040390660834E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.692484 -7.753850 -0.273932
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -3.761274 -4.628249 -0.260808
  locpot:
   site  1  z=  6.000  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.437191  rep=  -4.809995  q =   2.930583
 spin 2:           -5.206079        -4.035292        1.974368
  total:          -11.643269        -8.845287        4.904951
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 dfrce job=          12
 dfrce job=          12
 vxcnsp: loc rmu=  -6.437191  rep=  -4.809995  q =   2.930583
 spin 2:           -5.206079        -4.035292        1.974368
  total:          -11.643269        -8.845287        4.904951
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 dfrce job=          12
 vxcnsp: loc rmu=  -7.667669  rep=  -5.627154  q =   4.712387
 spin 2:           -4.559613        -3.707198        2.951738
  total:          -12.227282        -9.334352        7.664125

 local terms:     true           smooth         local
 rhoeps:        -8.845287      -9.334352       0.489065
 rhomu:         -6.437191      -7.667669       1.230478
 spin2:         -5.206079      -4.559613      -0.646465
 total:        -11.643269     -12.227282       0.584013
 val*vef       -13.707075     -34.191036      20.483960
 val chg:        2.904951       7.664125      -4.759174
 val mmom:       -0.804434  core mmom:   0.000000
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92025   2.91658   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.12000   0.00000   0.00000
 mkpot:
   Energy terms(Ry):       smooth           local           total
   rhoval*veff            -34.218670        20.483960       -13.734710
   Eestatic                 5.898853      -145.183724      -139.284871
   rho*exc                 -9.453759         0.489065        -8.964694
   rho*vxc                -12.382098         0.584013       -11.798085
   valence chg              7.759174        -4.759174         3.000000
   valence mag             1.804434        -0.804434         1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       10.742241 Ekcore=        62.961761 Ektot    =       73.704002
 Exc=      -8.964694 Ees   =      -139.284871 EKohnSham=      -74.545563
 Magnetic moment=     1.000000
 dfrce job=          12
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4358 -0.22222989777938548E-3
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4358 -0.22222855211453634E-3
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4358 -0.22222922494696091E-3

 Harris correction to forces: no shift of atomic density
  ib         delta-n dVes             delta-n dVxc               total
 Maximum Harris force =   0.000000 mRy/au (site 1 )
   1   -0.00   -0.00    0.00     0.00    0.00    0.00     0.00    0.00   -0.00
 shift forces to make zero average correction:            0.00    0.00   -0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force =   0.000000 mRy/au (site 1 )
  Symmetrize forces ...
 wgtsmooth=   8.0000000000000002E-003
 mixrho: sought 2 iter from file mixm ; read 3 RMS DQ= 6.38E-4  last it= 4.39E-3
 mmom         1.825940      1.804434
 AMIX: nmix=2 mmix=8  nelts= 33302  beta=0.50000  tm= 5.00000  rmsdel= 3.19D-04
   tj:-0.12580  -0.00501
 mixrho: qcell,correction,qmx,summ =  0.10800D-06  0.53999D-10 -0.47540D+01  0.12116D+03
 mixrho: warning. negative smrho; isp number min=    4360 -0.22246D-03

iors: write rst restart file (binary mesh density)

   it  8  of 10    ehf=     -74.545365   ehk=     -74.545563
 From last iter    ehf=     -74.545266   ehk=     -74.545418
 diffe(q)= -0.000098 (0.000638)    tol= 0.000010 (0.000500)   more=T
i mmom= 1.0000 ehf(eV)=-1014.249324 ehk(eV)=-1014.252019 sev(eV)=-40.693437

--- BNDFP:  begin iteration 9 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396184   -0.111762

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396184   -0.111762

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396184   -0.111762
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111762

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396184   -0.111762
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.891972   charge     8.754015
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4360 -0.22245567502426076E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.691224 -7.752905 -0.273879
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -3.755176 -4.619562 -0.260752
  locpot:
   site  1  z=  6.000  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.437607  rep=  -4.810292  q =   2.930643
 spin 2:           -5.206238        -4.035431        1.974364
  total:          -11.643845        -8.845723        4.905007
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.437607  rep=  -4.810292  q =   2.930643
 spin 2:           -5.206238        -4.035431        1.974364
  total:          -11.643845        -8.845723        4.905007
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -7.666802  rep=  -5.625952  q =   4.711674
 spin 2:           -4.550979        -3.701141        2.947348
  total:          -12.217781        -9.327093        7.659022

 local terms:     true           smooth         local
 rhoeps:        -8.845723      -9.327093       0.481369
 rhomu:         -6.437607      -7.666802       1.229195
 spin2:         -5.206238      -4.550979      -0.655259
 total:        -11.643845     -12.217781       0.573936
 val*vef       -13.709916     -34.157245      20.447329
 val chg:        2.905007       7.659022      -4.754015
 val mmom:       -0.808047  core mmom:  -0.000000
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92025   2.91658   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.12000   0.00000   0.00000
     potential shift to crystal energy zero:    0.000000
 mkpot:
   Energy terms(Ry):       smooth           local           total
   rhoval*veff            -34.184857        20.447329       -13.737528
   Eestatic                 5.891972      -145.180124      -139.288152
   rho*exc                 -9.446400         0.481369        -8.965030
   rho*vxc                -12.372467         0.573936       -11.798531
   valence chg              7.754015        -4.754015         3.000000
   valence mag             1.808034        -0.808047         0.999987
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:      0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     3 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.0454
 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.625086;   3.000000 electrons
         Sum occ. bands:   -2.991755, incl. Bloechl correction: -0.000039
         Mag. moment:       1.000000
 bndfp:Generating TDOS: efermi(eV)= -8.504794 DOSwindow emin emax(eV)=  -17.104150  31.495206


 m_bandcal_2nd: to fill eigenfunctions**2 up to Efermi
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.31041  sum tc=    31.42870  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.26441  sum tc=    31.53353  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.904989    7.662951   -4.757962      0.956250    1.763327   -0.807077
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.991755 Vin*nin=     -13.737528 Ek=Eb-Vin*nin=      10.745773
 Ek(core)=      62.961929 Exc=      -8.965030 Ees=    -139.288152 Eharris=     -74.545480

 mkekin:
   nout*Vin = smpart,onsite,total=:    -34.199951     20.464629    -13.735323
    E_B(band energy sum)=   -2.991755  E_B-nout*Vin=   10.743567

 m_mkpot_energyterms
 Energy for background charge q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396133   -0.111747

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396133   -0.111747

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396133   -0.111747
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111747

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396133   -0.111747
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.897677   charge     8.757962
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4362 -0.22266800526065766E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.693520 -7.755761 -0.273940
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -3.758639 -4.624251 -0.260815
  locpot:
   site  1  z=  6.000  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.437386  rep=  -4.810134  q =   2.930620
 spin 2:           -5.206145        -4.035351        1.974370
  total:          -11.643531        -8.845485        4.904989
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 dfrce job=          12
 vxcnsp: loc rmu=  -6.437386  rep=  -4.810134  q =   2.930620
 spin 2:           -5.206145        -4.035351        1.974370
  total:          -11.643531        -8.845485        4.904989
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 dfrce job=          12
 dfrce job=          12
 vxcnsp: loc rmu=  -7.669568  rep=  -5.628179  q =   4.713139
 spin 2:           -4.555596        -3.704548        2.949812
  total:          -12.225164        -9.332727        7.662951

 local terms:     true           smooth         local
 rhoeps:        -8.845485      -9.332727       0.487242
 rhomu:         -6.437386      -7.669568       1.232182
 spin2:         -5.206145      -4.555596      -0.650549
 total:        -11.643531     -12.225164       0.581633
 val*vef       -13.708118     -34.183242      20.475124
 val chg:        2.904989       7.662951      -4.757962
 val mmom:       -0.807077  core mmom:   0.000000
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92025   2.91659   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.12000   0.00000   0.00000
 mkpot:
   Energy terms(Ry):       smooth           local           total
   rhoval*veff            -34.210873        20.475124       -13.735750
   Eestatic                 5.897677      -145.184155      -139.286478
   rho*exc                 -9.452159         0.487242        -8.964917
   rho*vxc                -12.380012         0.581633       -11.798379
   valence chg              7.757962        -4.757962         3.000000
   valence mag             1.807077        -0.807077         1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       10.743567 Ekcore=        62.962234 Ektot    =       73.705801
 Exc=      -8.964917 Ees   =      -139.286478 EKohnSham=      -74.545594
 Magnetic moment=     1.000000
 dfrce job=          12
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4360 -0.22245634856014846E-3
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4358 -0.22245500148837305E-3
 Maximum Harris force =   0.000000 mRy/au (site 1 )
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4360 -0.22245567502426076E-3

 Harris correction to forces: no shift of atomic density
  ib         delta-n dVes             delta-n dVxc               total
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
   1   -0.00   -0.00    0.00     0.00    0.00    0.00     0.00    0.00   -0.00
 shift forces to make zero average correction:            0.00    0.00   -0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force =   0.000000 mRy/au (site 1 )
  Symmetrize forces ...
 wgtsmooth=   8.0000000000000002E-003
 mixrho: sought 2 iter from file mixm ; read 3 RMS DQ= 5.51E-5  last it= 6.38E-4
 mmom         1.808034      1.807077
 AMIX: nmix=2 mmix=8  nelts= 33302  beta=0.50000  tm= 5.00000  rmsdel= 2.76D-05
   tj:-0.46424   0.05211
 mixrho: qcell,correction,qmx,summ =  0.97195D-07  0.48597D-10 -0.47578D+01  0.12122D+03
 mixrho: warning. negative smrho; isp number min=    4362 -0.22261D-03

iors: write rst restart file (binary mesh density)

   it  9  of 10    ehf=     -74.545480   ehk=     -74.545594
 From last iter    ehf=     -74.545365   ehk=     -74.545563
 diffe(q)= -0.000116 (0.000055)    tol= 0.000010 (0.000500)   more=T
i mmom= 1.0000 ehf(eV)=-1014.250896 ehk(eV)=-1014.252441 sev(eV)=-40.705224

--- BNDFP:  begin iteration 10 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396144   -0.111750

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396144   -0.111750

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396144   -0.111750

 site class  ilm      vval      ves(rmax)
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111750
   1     1     1   -0.396144   -0.111750
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.896357   charge     8.757827
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4362 -0.22260961556898066E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.693341 -7.755463 -0.273923
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -3.758890 -4.624646 -0.260799
  locpot:
   site  1  z=  6.000  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.437439  rep=  -4.810173  q =   2.930626
 spin 2:           -5.206177        -4.035377        1.974371
  total:          -11.643616        -8.845550        4.904997
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.437439  rep=  -4.810173  q =   2.930626
 spin 2:           -5.206177        -4.035377        1.974371
  total:          -11.643616        -8.845550        4.904997
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -7.669294  rep=  -5.628020  q =   4.712921
 spin 2:           -4.556011        -3.704814        2.949903
  total:          -12.225305        -9.332834        7.662824

 local terms:     true           smooth         local
 rhoeps:        -8.845550      -9.332834       0.487284
 rhomu:         -6.437439      -7.669294       1.231856
 spin2:         -5.206177      -4.556011      -0.650166
 total:        -11.643616     -12.225305       0.581689
 val*vef       -13.708621     -34.182056      20.473436
 val chg:        2.904997       7.662824      -4.757827
 val mmom:       -0.806763  core mmom:  -0.000000
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92025   2.91659   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.12000   0.00000   0.00000
     potential shift to crystal energy zero:    0.000000
 mkpot:
   Energy terms(Ry):       smooth           local           total
   rhoval*veff            -34.209682        20.473436       -13.736247
   Eestatic                 5.896357      -145.183241      -139.286883
   rho*exc                 -9.452231         0.487284        -8.964947
   rho*vxc                -12.380109         0.581689       -11.798419
   valence chg              7.757827        -4.757827         3.000000
   valence mag             1.806759        -0.806763         0.999996
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:      0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     3 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.0281
 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.625188;   3.000000 electrons
         Sum occ. bands:   -2.992058, incl. Bloechl correction: -0.000039
         Mag. moment:       1.000000
 bndfp:Generating TDOS: efermi(eV)= -8.506182 DOSwindow emin emax(eV)=  -17.105496  31.493818


 m_bandcal_2nd: to fill eigenfunctions**2 up to Efermi
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.31070  sum tc=    31.42878  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.26470  sum tc=    31.53359  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905007    7.663791   -4.758784      0.956259    1.763479   -0.807220
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.992058 Vin*nin=     -13.736247 Ek=Eb-Vin*nin=      10.744188
 Ek(core)=      62.962081 Exc=      -8.964947 Ees=    -139.286883 Eharris=     -74.545561

 mkekin:
   nout*Vin = smpart,onsite,total=:    -34.213646     20.477515    -13.736131
    E_B(band energy sum)=   -2.992058  E_B-nout*Vin=   10.744073

 m_mkpot_energyterms
 Energy for background charge q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396130   -0.111746
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111746

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396130   -0.111746
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396130   -0.111746
   smooth rhoves      5.897570   charge     8.758784
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4362 -0.22268842642608042E-3

 site class  ilm      vval      ves(rmax)
   1     1     1   -0.396130   -0.111746
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.694469 -7.757057 -0.273942
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -3.759306 -4.625078 -0.260818
  locpot:
   site  1  z=  6.000  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.437454  rep=  -4.810183  q =   2.930633
 spin 2:           -5.206175        -4.035377        1.974374
  total:          -11.643630        -8.845560        4.905007
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 vxcnsp: loc rmu=  -6.437454  rep=  -4.810183  q =   2.930633
 spin 2:           -5.206175        -4.035377        1.974374
  total:          -11.643630        -8.845560        4.905007
 dfrce job=          12
 mesh:   nth,nph= -32   0   gives  32  angular points,   nrad= 369
 dfrce job=          12
 vxcnsp: loc rmu=  -7.670860  rep=  -5.629126  q =   4.713635
 spin 2:           -4.556418        -3.705210        2.950156
  total:          -12.227278        -9.334336        7.663791

 local terms:     true           smooth         local
 rhoeps:        -8.845560      -9.334336       0.488776
 rhomu:         -6.437454      -7.670860       1.233406
 spin2:         -5.206175      -4.556418      -0.649758
 total:        -11.643630     -12.227278       0.583648
 val*vef       -13.708512     -34.188327      20.479815
 val chg:        2.905007       7.663791      -4.758784
 val mmom:       -0.807220  core mmom:  -0.000000
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92025   2.91659   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.12000   0.00000   0.00000
 mkpot:
   Energy terms(Ry):       smooth           local           total
   rhoval*veff            -34.215957        20.479815       -13.736142
   Eestatic                 5.897570      -145.184623      -139.287052
   rho*exc                 -9.453774         0.488776        -8.964998
   rho*vxc                -12.382134         0.583648       -11.798486
   valence chg              7.758784        -4.758784         3.000000
   valence mag             1.807220        -0.807220         1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       10.744073 Ekcore=        62.962375 Ektot    =       73.706448
 Exc=      -8.964998 Ees   =      -139.287052 EKohnSham=      -74.545603
 Magnetic moment=     1.000000
 dfrce job=          12
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4362 -0.22261028955814654E-3
 dfrce job=          12
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4362 -0.22260894157981478E-3
smvxcm: smrho_w<minimumrho  number,min(smrho_w)= 4362 -0.22260961556898066E-3
 Maximum Harris force =   0.000000 mRy/au (site 1 )

 Harris correction to forces: no shift of atomic density
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00    0.00    0.00     0.00    0.00    0.00     0.00   -0.00   -0.00
 shift forces to make zero average correction:            0.00   -0.00   -0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force =   0.000000 mRy/au (site 1 )
  Symmetrize forces ...
 wgtsmooth=   8.0000000000000002E-003
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 mixrho: sought 2 iter from file mixm ; read 3 RMS DQ= 1.53E-5  last it= 5.51E-5
 mmom         1.806759      1.807220
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts= 33302  beta=0.50000  tm= 5.00000  rmsdel= 7.65D-06
   tj:-0.13385
 mixrho: qcell,correction,qmx,summ =  0.10277D-06  0.51385D-10 -0.47586D+01  0.12123D+03
 mixrho: warning. negative smrho; isp number min=    4362 -0.22266D-03

iors: write rst restart file (binary mesh density)

   it 10  of 10    ehf=     -74.545561   ehk=     -74.545603
 From last iter    ehf=     -74.545480   ehk=     -74.545594
 diffe(q)= -0.000080 (0.000015)    tol= 0.000010 (0.000500)   more=F
x mmom= 1.0000 ehf(eV)=-1014.251991 ehk(eV)=-1014.252564 sev(eV)=-40.709347
 OK! end of LMF ======================
