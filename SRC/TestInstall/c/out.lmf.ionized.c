INFO: Ubuntu 20.04.4 LTS \n \l
INFO: GNU Fortran (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0
INFO: -O2 -g -fimplicit-none -finit-integer=NaN -finit-real=NaN -JOBJ.gfortran -IOBJ.gfortran
INFO: MATH: -lmkl_rt
INFO: git: commit 895c0dab24c443e7d7fb788e9870fc04186617c7
INFO:    : Author: Takao Kotani <takaokotani@gmail.com>
INFO:    : Date:   Tue Mar 28 19:00:32 2023 +0900
INFO: linked at Wed Mar 29 10:09:58 JST 2023
=== START LFMA ===
 mpisize=           1
m_lmfinit: LMFA
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
rval2: HAM_V0FIX               defa n= 1 val= 0.00000000
rval2: HAM_PNUFIX              defa n= 1 val= 0.00000000
=== SPEC =1
 cccccccc SPEC_ATOM@1ch=### C###
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
 cccccccc SPEC_C-HOLE@1ch=######
rval2: SPEC_C-HQ@1             defa n= 2 val= -1.00000000  0.00000000
rval2: SPEC_EREF1              defa n= 1 val= 0.00000000
=== SITE =1
 cccccccc SITE_ATOM@1ch=### C###
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
 cccccccc SYMGRPch=###  find###
 cccccccc SYMGRPAFch=######
rval2: EWALD_AS                defa n= 1 val= 2.00000000
rval2: EWALD_TOL               defa n= 1 val= 0.00000001
rval2: EWALD_NKDMX             defa n= 1 val= 300.00000000
rval2: ITER_NIT                defa n= 1 val= 10.00000000
rval2: ITER_NRMIX              defa n= 1 val= 80.00000000
 cccccccc ITER_MIXch=### A2###
rval2: ITER_CONV               defa n= 1 val= 0.00001000
rval2: ITER_CONVC              defa n= 1 val= 0.00050000
rval2: ITER_UMIX               defa n= 1 val= 0.50000000
rval2: ITER_TOLU               defa n= 1 val= 0.00000000
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
pnu list       ibas isp  pnu(0:lmxa) 
pnu: j isp pnu= 1 1 2.900  2.850  3.180  4.120
pnz: j isp  pz= 1 1 0.000  0.000  0.000  0.000
pnu: j isp pnu= 1 2 2.900  2.850  3.180  4.120
pnz: j isp  pz= 1 2 0.000  0.000  0.000  0.000

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
           1 --no-iactiv
           2 c
           3 -vzbak=1
INFO: Ubuntu 20.04.4 LTS \n \l
INFO: GNU Fortran (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0
           1 --no-iactiv
           2 c
           3 -vzbak=1
           1 --no-iactiv
           2 c
           3 -vzbak=1
           1 --no-iactiv
           2 c
           3 -vzbak=1
INFO: -O2 -g -fimplicit-none -finit-integer=NaN -finit-real=NaN -JOBJ.gfortran -IOBJ.gfortran
INFO: MATH: -lmkl_rt
INFO: git: commit 895c0dab24c443e7d7fb788e9870fc04186617c7
INFO:    : Author: Takao Kotani <takaokotani@gmail.com>
INFO:    : Date:   Tue Mar 28 19:00:32 2023 +0900
INFO: linked at Wed Mar 29 10:09:58 JST 2023
===START LMF with   --no-iactiv c -vzbak=1 ===
mpisize=4
m_lmfinit: LMF
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
 cccccccc SPEC_ATOM@1ch=### C###
rval2: HAM_GMAX                defa n= 1 val= 0.00000000
rval2: HAM_FTMESH              defa n= 3 val= 25.00000000  25.00000000  25.00000000
rval2: HAM_TOL                 defa n= 1 val= 0.00000100
 cccccccc SPEC_ATOM@1ch=### C###
 cccccccc SPEC_ATOM@1ch=### C###
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
=== SPEC =1
 cccccccc SPEC_C-HOLE@1ch=######
 cccccccc SPEC_ATOM@1ch=### C###
rval2: SPEC_Z@1                ---- n= 1 val= 6.00000000
 cccccccc SITE_ATOM@1ch=### C###
rval2: SPEC_R@1                ---- n= 1 val= 3.00000000
rval2: SPEC_R/W@1              ---- n= 0 val= 
rval2: SPEC_R/A@1              ---- n= 0 val= 
rval2: SPEC_A@1                defa n= 1 val= 0.02000000
 cccccccc SPEC_C-HOLE@1ch=######
 cccccccc SPEC_C-HOLE@1ch=######
rval2: SPEC_NR@1               defa n= 1 val= 0.00000000
 cccccccc SITE_ATOM@1ch=### C###
 cccccccc SITE_ATOM@1ch=### C###
rval2: SPEC_RSMH@1             ---- n= 4 val= 1.30000000  1.10000000 -1.00000000 -1.00000000
rval2: SPEC_EH@1               requ n= 2 val= -0.70000000 -0.20000000
rval2: SPEC_RSMH2@1            ---- n= 2 val= 0.80000000  0.80000000
rval2: SPEC_EH2@1              requ n= 2 val= -1.50000000 -1.00000000
rval2: SPEC_LMX@1              defa n= 1 val= 999.00000000
 cccccccc SYMGRPch=###  find###
rval2: SPEC_LMXA@1             defa n= 1 val= 3.00000000
 cccccccc SYMGRPAFch=######
rval2: SPEC_LMXL@1             defa n= 1 val= 3.00000000
rval2: SPEC_P@1                ---- n= 4 val= 2.90000000  2.85000000  3.18000000  4.12000000
 cccccccc ITER_MIXch=### A2###
rval2: SPEC_Q@1                ---- n= 0 val= 
 cccccccc SYMGRPch=###  find###
rval2: SPEC_MMOM@1             ---- n= 2 val= 0.00000000  2.00000000
 cccccccc SYMGRPch=###  find###
 cccccccc SYMGRPAFch=######
rval2: SPEC_NMCORE@1           defa n= 1 val= 0.00000000
 cccccccc SYMGRPAFch=######
rval2: SPEC_PZ@1               ---- n= 0 val= 
rval2: SPEC_LFOCA@1            defa n= 1 val= 0.00000000
rval2: SPEC_KMXA@1             defa n= 1 val= 3.00000000
 cccccccc ITER_MIXch=### A2###
 cccccccc ITER_MIXch=### A2###
rval2: SPEC_RSMA@1             defa n= 1 val= 1.20000000
rval2: SPEC_IDMOD@1            ---- n= 2 val= 0.00000000  1.00000000
rval2: SPEC_FRZWF@1            defa n= 1 val= 0.00000000
rval2: SPEC_IDU@1              ---- n= 0 val= 
rval2: SPEC_UH@1               ---- n= 0 val= 
rval2: SPEC_JH@1               ---- n= 0 val= 
 cccccccc SPEC_C-HOLE@1ch=######
rval2: SPEC_C-HQ@1             defa n= 2 val= -1.00000000  0.00000000
rval2: SPEC_EREF1              defa n= 1 val= 0.00000000
=== SITE =1
 cccccccc SITE_ATOM@1ch=### C###
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
 cccccccc SYMGRPch=###  find###
 cccccccc SYMGRPAFch=######
rval2: EWALD_AS                defa n= 1 val= 2.00000000
rval2: EWALD_TOL               defa n= 1 val= 0.00000001
rval2: EWALD_NKDMX             defa n= 1 val= 300.00000000
rval2: ITER_NIT                defa n= 1 val= 10.00000000
rval2: ITER_NRMIX              defa n= 1 val= 80.00000000
 cccccccc ITER_MIXch=### A2###
rval2: ITER_CONV               defa n= 1 val= 0.00001000
rval2: ITER_CONVC              defa n= 1 val= 0.00050000
rval2: ITER_UMIX               defa n= 1 val= 0.50000000
rval2: ITER_TOLU               defa n= 1 val= 0.00000000
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
pnu list       ibas isp  pnu(0:lmxa) 
pnu: j isp pnu= 1 1 2.900  2.850  3.180  4.120
pnz: j isp  pz= 1 1 0.000  0.000  0.000  0.000
pnu: j isp pnu= 1 2 2.900  2.850  3.180  4.120
pnz: j isp  pz= 1 2 0.000  0.000  0.000  0.000

mto === MTO setting ===
mto ispec lmxb lpz nkapii nkaphh=    1    1    0    2    2
mto rsmh1    1  1.30  1.10
mto   eh1    1 -0.70 -0.20
mto rsmh2    1  0.80  0.80

                Plat                                  Qlat

                Plat                                  Qlat

                Plat                                  Qlat
mto  eh2     1 -1.50 -1.00
mto lh       1  1
   1.000000   1.000000   0.000000        0.500000   0.500000  -0.500000
   1.000000   0.000000   1.000000        0.500000  -0.500000   0.500000
   1.000000   1.000000   0.000000        0.500000   0.500000  -0.500000
   1.000000   0.000000   1.000000        0.500000  -0.500000   0.500000
   1.000000   1.000000   0.000000        0.500000   0.500000  -0.500000
   1.000000   0.000000   1.000000        0.500000  -0.500000   0.500000
   0.000000   1.000000   1.000000       -0.500000   0.500000   0.500000

                Plat                                  Qlat
   1.000000   1.000000   0.000000        0.500000   0.500000  -0.500000
   0.000000   1.000000   1.000000       -0.500000   0.500000   0.500000
   0.000000   1.000000   1.000000       -0.500000   0.500000   0.500000
   1.000000   0.000000   1.000000        0.500000  -0.500000   0.500000
   0.000000   1.000000   1.000000       -0.500000   0.500000   0.500000
  Cell vol=   1000.000000
  Cell vol=   1000.000000
  Cell vol=   1000.000000
  Cell vol=   1000.000000

LATTC:  as= 2.000   tol= 1.00E-08   alat= 7.93701   awald= 0.200

         r1=  3.459   nkd=  79      q1=  2.571   nkq= 137

LATTC:  as= 2.000   tol= 1.00E-08   alat= 7.93701   awald= 0.200
         r1=  3.459   nkd=  79      q1=  2.571   nkq= 137
LATTC:  as= 2.000   tol= 1.00E-08   alat= 7.93701   awald= 0.200
         r1=  3.459   nkd=  79      q1=  2.571   nkq= 137
SpaceGroupSym: ======================================= 

LATTC:  as= 2.000   tol= 1.00E-08   alat= 7.93701   awald= 0.200
  Generators except find=
         r1=  3.459   nkd=  79      q1=  2.571   nkq= 137
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
SpaceGroupSym: ========end of SYM section============= 

 BZMESH:      3 irreducible QP from    2   2   2 shift=FFF
 TETIRR: sorting       48 tetrahedra ...
 SGVSYM: 188 symmetry stars found for 5817 reciprocal lattice vectors

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol= 1.00E-06
 spec      l    rsm    eh     gmax    last term    cutoff
  C        0    1.30  -0.70   5.718    1.05E-06    3143 
  C        1    1.10  -0.20   7.226    2.65E-06    5817*
  C        0    0.80  -1.50   9.292    4.01E-04    5817*
  C        1    0.80  -1.00  10.038    2.81E-03    5817*
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
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       25082  -4.9988145201717389E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.49988145202717392E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -1.109295 -1.522745 -0.218019
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -0.385560 -0.423747 -0.164195
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       128
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       128
 vxcnsp (warning): negative rho: min val =  -4.61E-04
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       128
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       128
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.90000   2.90000   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.18000   3.18000   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.12000   4.12000   0.00000   0.00000
     potential shift to crystal energy zero:    0.000002
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -3.733888       -10.403600       -14.137487
   rhoval*ves             -5.058556        -5.181775       -10.240330
   psnuc*ves               9.186377      -278.836573      -269.650196
   utot                    2.063910      -142.009174      -139.945263
   rho*exc                -1.494856        -8.100688        -9.595544
   rho*vxc                -1.946492       -10.678430       -12.624921
   valence chg             1.661718         1.338282         3.000000
   valence mag             1.332841         0.667159         2.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     3 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.0184

 bzwts: --- Tetrahedron Integration ---
 ... only filled or empty bands encountered: ev= -0.824735 ec= -0.768586
 VBmax= -0.824735 CBmin= -0.768586 gap =  0.056150 Ry =   0.763959 eV
 BZINTS: Fermi energy:     -0.824735;   3.000000 electrons
         Sum occ. bands:   -3.222313, incl. Bloechl correction:  0.000000
         Mag. moment:      -1.000000
 bndfp:Generating TDOS: efermi= -0.824735  dos window emin emax=  -1.423654  2.115188


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
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.207662
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves     36.495506   charge   366.036673
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=          16  -1.6773260017262086E-007
  smvxcm: enforce positive smrho_w. Add srshift= 0.16773261017262087E-6
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -152.338263 -96.188303 -0.226907
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -2141.658238 -2948.665992 -0.408385
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.89276   2.76013   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.12000   4.12000   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff          -2517.210592      2505.736791       -11.473801
   rhoval*ves             77.961664       -86.892741        -8.931077
   psnuc*ves              -4.970652      -258.950387      -263.921039
   utot                   36.495506      -172.921564      -136.426058
   rho*exc             -2293.996501      2285.781076        -8.215426
   rho*vxc             -3044.854295      3034.038598       -10.815697
   valence chg           365.036673      -362.036676         2.999997
   valence mag          -314.161749       313.161752        -0.999997
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=        8.399014 Ekcore=        62.933572 Ektot    =       71.332586
 Exc=      -8.215426 Ees   =      -136.426058 EKohnSham=      -73.308898
 Magnetic moment=    -0.999997
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       25082  -4.9988446984781346E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.49988446985781349E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       25082  -4.9987843418653442E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.49987843419653445E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       25082  -4.9988145201717400E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.49988145202717403E-3

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
 mixrealsmooth= T
 wgtsmooth=   8.0000000000000002E-003
 mixrho: sought 2 iter from file mixm ; read 0 RMS DQ= 9.54E+0
 AMIX: nmix=0 mmix=8  nelts= 33302  beta=0.50000  tm= 5.00000  rmsdel= 4.77D+00
 mixrho: warning. negative smrho; isp number min=       1   11628 -0.24843D-03
 mixrho: warning. negative smrho; isp number min=       1   11628 -0.24843D-03
 mixrho: add corrections to qcell smrho =  0.17512D-05  0.87558D-09
 mixrho: warning. negative smrho; isp number min=       2    7318 -0.24382D-03
 mixrho: warning. negative smrho; isp number min=       2    7318 -0.24382D-03
 add q=  0.000002 to preserve neutrality
 mixrho: warning. negative smrho; isp number min=       1   11628 -0.24843D-03
 mixrho: warning. negative smrho; isp number min=       2    7318 -0.24382D-03

 iors  : write rst restart file (binary mesh density)
 mixrho: warning. negative smrho; isp number min=       1   11628 -0.24843D-03
 mixrho: warning. negative smrho; isp number min=       2    7318 -0.24382D-03

   it  1  of 10    ehf=     -75.693628   ehk=     -73.308898
h mmom=-1.0000 ehf(eV)=-1029.872362 ehk(eV)=-997.426200 sev(eV)=-43.842140

--- BNDFP:  begin iteration 2 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.100527
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      8.326727   charge   184.349197
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       18946  -2.4843197179100451E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.24843197180100449E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -62.642520 -40.704435 -0.244834
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -854.219866 -1175.218090 -0.340641
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0        64
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0        64
 vxcnsp (warning): negative rho: min val =  -1.96E-04
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0        64
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0        64
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.89276   2.76013   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.12000   4.12000   0.00000   0.00000
     potential shift to crystal energy zero:   -0.001086
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff          -1296.369827      1283.323158       -13.046669
   rhoval*ves             14.545593       -24.473846        -9.928254
   psnuc*ves               2.107862      -268.893479      -266.785617
   utot                    8.326727      -146.683663      -138.356935
   rho*exc              -916.862385       907.993418        -8.868967
   rho*vxc             -1215.922524      1204.251622       -11.670902
   valence chg           183.349197      -180.349197         3.000000
   valence mag          -156.414454       156.914456         0.500001
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:      0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     3 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.0351

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.722129;   3.000000 electrons
         Sum occ. bands:   -3.290097, incl. Bloechl correction: -0.000017
         Mag. moment:       1.000000
 bndfp:Generating TDOS: efermi= -0.722129  dos window emin emax=  -1.351906  2.217793

       contr. to mm extrapolated for r>rmt:  -0.054502 est. true mm =-0.989849
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.56622  sum tc=    31.48945  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.53090  sum tc=    31.57762  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.895099    6.797414   -3.902315     -0.935348   -0.530782   -0.404566
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -3.290097 Vin*nin=     -13.046669 Ek=Eb-Vin*nin=       9.756572
 Ek(core)=      62.932634 Exc=      -8.868967 Ees=    -138.356935 Eharris=     -74.536697

 mkekin:
   nout*Vin = smpart,onsite,total=:    -32.521180     18.413046    -14.108134
    E_B(band energy sum)=   -3.290097  E_B-nout*Vin=   10.818037

 m_mkpot_energyterms
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.114454
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.255246   charge     7.902315
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        3176  -1.3866407181040646E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.13866407182040646E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -3.869309 -5.044852 -0.237526
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -4.089619 -5.373249 -0.253032
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92164   2.91793   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.12000   4.10242   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -28.805532        14.908810       -13.896722
   rhoval*ves             -4.093578        -6.523657       -10.617235
   psnuc*ves              14.604071      -283.105631      -268.501560
   utot                    5.255246      -144.814644      -139.559398
   rho*exc                -7.958928        -0.976854        -8.935783
   rho*vxc               -10.418101        -1.342865       -11.760966
   valence chg             6.902315        -3.902315         3.000000
   valence mag            -0.595434        -0.404566        -1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       10.818037 Ekcore=        63.067070 Ektot    =       73.885108
 Exc=      -8.935783 Ees   =      -139.559398 EKohnSham=      -74.610072
 Magnetic moment=    -1.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       18946  -2.4843346219173123E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.24843346220173120E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       18946  -2.4843048139027780E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.24843048140027777E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       18946  -2.4843197179100451E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.24843197180100449E-3
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
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
 wgtsmooth=   8.0000000000000002E-003
 mixrho: sought 2 iter from file mixm ; read 1 RMS DQ= 4.72E+0  last it= 9.54E+0
 AMIX: nmix=1 mmix=8  nelts= 33302  beta=0.50000  tm= 5.00000  rmsdel= 2.36D+00
   tj: 0.33081
 mixrho: warning. negative smrho; isp number min=       1   11532 -0.16531D-03
 mixrho: warning. negative smrho; isp number min=       2    7194 -0.20842D-03
 mixrho: warning. negative smrho; isp number min=       1   11532 -0.16531D-03
 mixrho: warning. negative smrho; isp number min=       2    7194 -0.20842D-03
 mixrho: add corrections to qcell smrho =  0.66791D-06  0.33395D-09
 mixrho: warning. negative smrho; isp number min=       1   11532 -0.16531D-03
 mixrho: warning. negative smrho; isp number min=       2    7194 -0.20842D-03
 mixrho: warning. negative smrho; isp number min=       1   11532 -0.16531D-03
 mixrho: warning. negative smrho; isp number min=       2    7194 -0.20842D-03

 iors  : write rst restart file (binary mesh density)

   it  2  of 10    ehf=     -74.536697   ehk=     -74.610072
 From last iter    ehf=     -75.693628   ehk=     -73.308898
 diffe(q)=  1.156931 (4.716723)    tol= 0.000010 (0.000500)   more=T
i mmom=-1.0000 ehf(eV)=-1014.131387 ehk(eV)=-1015.129724 sev(eV)=-44.764395

--- BNDFP:  begin iteration 3 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.105187
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      3.978639   charge   125.311214
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       18726  -2.0841839072584644E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.20841839073584644E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -40.425980 -27.285906 -0.251258
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -501.816838 -691.373215 -0.325123
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92164   2.91793   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.12000   4.10242   0.00000   0.00000
     potential shift to crystal energy zero:   -0.000722
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff           -907.707856       894.423367       -13.284489
   rhoval*ves              1.668265       -11.824004       -10.155740
   psnuc*ves               6.289013      -273.690861      -267.401848
   utot                    3.978639      -142.757433      -138.778794
   rho*exc              -542.242818       533.380694        -8.862124
   rho*vxc              -718.659122       706.997789       -11.661332
   valence chg           124.311214      -121.311214         3.000000
   valence mag          -104.278420       104.276531        -0.001889
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     3 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.0368

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.648748;   3.000000 electrons
         Sum occ. bands:   -3.154593, incl. Bloechl correction: -0.000028
         Mag. moment:       1.000000
 bndfp:Generating TDOS: efermi= -0.648748  dos window emin emax=  -1.279720  2.291175

       contr. to mm extrapolated for r>rmt:   0.035229 est. true mm = 0.992227
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.45435  sum tc=    31.49432  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.44605  sum tc=    31.51771  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905570    6.818720   -3.913150      0.956998    2.101150   -1.144152
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -3.154593 Vin*nin=     -13.284489 Ek=Eb-Vin*nin=      10.129896
 Ek(core)=      62.999833 Exc=      -8.862124 Ees=    -138.778794 Eharris=     -74.511189

 mkekin:
   nout*Vin = smpart,onsite,total=:    -31.713112     17.620382    -14.092730
    E_B(band energy sum)=   -3.154593  E_B-nout*Vin=   10.938137

 m_mkpot_energyterms
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.112982
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.484754   charge     7.913150
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        4346  -2.1521390439155894E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.21521390440155895E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.190740 -7.142707 -0.271191
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -2.868278 -3.410607 -0.256685
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92084   2.91944   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -28.966404        15.071787       -13.894617
   rhoval*ves             -3.883349        -6.712846       -10.596195
   psnuc*ves              14.852858      -283.272224      -268.419366
   utot                    5.484754      -144.992535      -139.507781
   rho*exc                -8.059017        -0.923859        -8.982876
   rho*vxc               -10.553314        -1.268881       -11.822196
   valence chg             6.913150        -3.913150         3.000000
   valence mag             2.144152        -1.144152         1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       10.938137 Ekcore=        63.012036 Ektot    =       73.950173
 Exc=      -8.982876 Ees   =      -139.507781 EKohnSham=      -74.540484
 Magnetic moment=     1.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       18726  -2.0841952205776447E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.20841952206776448E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       18726  -2.0841725939392840E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.20841725940392840E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       18726  -2.0841839072584644E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.20841839073584644E-3

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
 wgtsmooth=   8.0000000000000002E-003
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 mixrho: sought 2 iter from file mixm ; read 2 RMS DQ= 3.15E+0  last it= 4.72E+0
 AMIX: nmix=2 mmix=8  nelts= 33302  beta=0.50000  tm= 5.00000  rmsdel= 1.57D+00
   tj:-0.28865   0.21227
 mixrho: warning. negative smrho; isp number min=       1   10766 -0.22163D-03
 mixrho: warning. negative smrho; isp number min=       2    7226 -0.10869D-03
 mixrho: warning. negative smrho; isp number min=       1   10766 -0.22163D-03
 mixrho: warning. negative smrho; isp number min=       1   10766 -0.22163D-03
 mixrho: warning. negative smrho; isp number min=       2    7226 -0.10869D-03
 mixrho: warning. negative smrho; isp number min=       2    7226 -0.10869D-03
 mixrho: add corrections to qcell smrho =  0.43723D-06  0.21862D-09
 mixrho: warning. negative smrho; isp number min=       1   10766 -0.22163D-03
 mixrho: warning. negative smrho; isp number min=       2    7226 -0.10869D-03

 iors  : write rst restart file (binary mesh density)

   it  3  of 10    ehf=     -74.511189   ehk=     -74.540484
 From last iter    ehf=     -74.536697   ehk=     -74.610072
 diffe(q)=  0.025508 (3.146534)    tol= 0.000010 (0.000500)   more=T
i mmom= 1.0000 ehf(eV)=-1013.784330 ehk(eV)=-1014.182919 sev(eV)=-42.920761

--- BNDFP:  begin iteration 4 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.107728
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      2.935135   charge    83.085025
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       17992  -2.2162581010208213E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.22162581011208213E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -27.476433 -20.011996 -0.269226
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -279.529657 -386.521480 -0.314507
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92084   2.91944   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
     potential shift to crystal energy zero:   -0.000458
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff           -596.178776       582.665338       -13.513438
   rhoval*ves             -3.543568        -6.765925       -10.309492
   psnuc*ves               9.413837      -277.188826      -267.774990
   utot                    2.935135      -141.977376      -139.042241
   rho*exc              -307.006090       298.081634        -8.924455
   rho*vxc              -406.533476       394.788817       -11.744659
   valence chg            82.085025       -79.085025         3.000000
   valence mag           -65.509200        66.224668         0.715468
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     3 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.0368

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.655639;   3.000000 electrons
         Sum occ. bands:   -3.096185, incl. Bloechl correction: -0.000031
         Mag. moment:       1.000000
 bndfp:Generating TDOS: efermi= -0.655639  dos window emin emax=  -1.286891  2.284283

       contr. to mm extrapolated for r>rmt:   0.032287 est. true mm = 0.992461
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.39470  sum tc=    31.44898  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.35578  sum tc=    31.53824  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.904071    6.762346   -3.858275      0.960173    2.289366   -1.329193
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -3.096185 Vin*nin=     -13.513438 Ek=Eb-Vin*nin=      10.417253
 Ek(core)=      63.005946 Exc=      -8.924455 Ees=    -139.042241 Eharris=     -74.543498

 mkekin:
   nout*Vin = smpart,onsite,total=:    -31.943576     18.006786    -13.936790
    E_B(band energy sum)=   -3.096185  E_B-nout*Vin=   10.840605

 m_mkpot_energyterms
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.112822
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.546507   charge     7.858275
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        4412  -2.2035236553194769E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.22035236554194770E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.277762 -7.289794 -0.272449
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -2.693345 -3.149323 -0.258123
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92099   2.91795   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -28.680232        14.864914       -13.815318
   rhoval*ves             -3.840655        -6.686841       -10.527496
   psnuc*ves              14.933669      -283.199765      -268.266096
   utot                    5.546507      -144.943303      -139.396796
   rho*exc                -7.971107        -1.003598        -8.974705
   rho*vxc               -10.439117        -1.372255       -11.811372
   valence chg             6.858275        -3.858275         3.000000
   valence mag             2.329193        -1.329193         1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       10.840605 Ekcore=        62.987216 Ektot    =       73.827821
 Exc=      -8.974705 Ees   =      -139.396796 EKohnSham=      -74.543679
 Magnetic moment=     1.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       17992  -2.2162681020872674E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.22162681021872674E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       17992  -2.2162480999543752E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.22162481000543752E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       17992  -2.2162581010208213E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.22162581011208213E-3
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
 wgtsmooth=   8.0000000000000002E-003
 mixrho: sought 2 iter from file mixm ; read 3 RMS DQ= 2.02E+0  last it= 3.15E+0
 AMIX: nmix=2 mmix=8  nelts= 33302  beta=0.50000  tm= 5.00000  rmsdel= 1.01D+00
   tj:-1.23656  -0.22945
 mixrho: warning. negative smrho; isp number min=       1    4434 -0.28124D-03
 mixrho: warning. negative smrho; isp number min=       2    6430 -0.37405D-04
 mixrho: warning. negative smrho; isp number min=       1    4434 -0.28124D-03
 mixrho: warning. negative smrho; isp number min=       2    6430 -0.37405D-04
 mixrho: warning. negative smrho; isp number min=       1    4434 -0.28124D-03
 mixrho: warning. negative smrho; isp number min=       2    6430 -0.37405D-04
 mixrho: add corrections to qcell smrho =  0.78314D-07  0.39157D-10
 mixrho: warning. negative smrho; isp number min=       1    4434 -0.28124D-03
 mixrho: warning. negative smrho; isp number min=       2    6430 -0.37405D-04

 iors  : write rst restart file (binary mesh density)

   it  4  of 10    ehf=     -74.543498   ehk=     -74.543679
 From last iter    ehf=     -74.511189   ehk=     -74.540484
 diffe(q)= -0.032309 (2.016665)    tol= 0.000010 (0.000500)   more=T
i mmom= 1.0000 ehf(eV)=-1014.223919 ehk(eV)=-1014.226390 sev(eV)=-42.126079

--- BNDFP:  begin iteration 5 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.112386
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.636692   charge     7.707627
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       10864  -2.8124465536371067E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.28124465537371064E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.708213 -7.936304 -0.287902
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -2.138429 -2.345841 -0.265343
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       864
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       864
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       864
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       864
 vxcnsp (warning): negative rho: min val =  -2.45E-02
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92099   2.91795   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
     potential shift to crystal energy zero:    0.000011
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -27.930264        14.066061       -13.864203
   rhoval*ves             -3.758332        -6.747615       -10.505947
   psnuc*ves              15.031717      -283.298267      -268.266551
   utot                    5.636692      -145.022941      -139.386249
   rho*exc                -7.846641        -1.207656        -9.054297
   rho*vxc               -10.282145        -1.635847       -11.917992
   valence chg             6.707627        -3.707627         3.000000
   valence mag             3.259179        -1.703751         1.555428
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     3 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.0362

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.648467;   3.000000 electrons
         Sum occ. bands:   -2.987727, incl. Bloechl correction: -0.000037
         Mag. moment:       1.000000
 bndfp:Generating TDOS: efermi= -0.648467  dos window emin emax=  -1.280479  2.291455

       contr. to mm extrapolated for r>rmt:   0.034188 est. true mm = 0.992510
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.30255  sum tc=    31.39856  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.22920  sum tc=    31.56335  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.908730    8.378369   -5.469639      0.958322    1.505213   -0.546891
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.987727 Vin*nin=     -13.864203 Ek=Eb-Vin*nin=      10.876476
 Ek(core)=      62.996539 Exc=      -9.054297 Ees=    -139.386249 Eharris=     -74.567531

 mkekin:
   nout*Vin = smpart,onsite,total=:    -36.740596     22.993281    -13.747315
    E_B(band energy sum)=   -2.987727  E_B-nout*Vin=   10.759588

 m_mkpot_energyterms
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.110967
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      6.093302   charge     9.469639
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        4496  -2.2786960695573248E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.22786960696573249E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -6.184644 -8.357874 -0.275146
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -4.552646 -5.709168 -0.263308
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91971   2.91256   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -38.764624        25.021545       -13.743078
   rhoval*ves             -3.398327        -7.062109       -10.460436
   psnuc*ves              15.584930      -283.718240      -268.133310
   utot                    6.093302      -145.390175      -139.296873
   rho*exc               -10.737290         1.764428        -8.972862
   rho*vxc               -14.067042         2.258242       -11.808800
   valence chg             8.469639        -5.469639         3.000000
   valence mag             1.546891        -0.546891         1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       10.759588 Ekcore=        62.961913 Ektot    =       73.721500
 Exc=      -8.972862 Ees   =      -139.296873 EKohnSham=      -74.548235
 Magnetic moment=     1.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       10864  -2.8124538968238619E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.28124538969238617E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       10864  -2.8124392104503514E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.28124392105503512E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       10864  -2.8124465536371067E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.28124465537371064E-3

 Harris correction to forces: screened shift in core+nuclear density  
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
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 mixrho: sought 2 iter from file mixm ; read 3 RMS DQ= 3.78E-2  last it= 2.02E+0
 AMIX: nmix=2 mmix=8  nelts= 33302  beta=0.50000  tm= 5.00000  rmsdel= 1.89D-02
   tj: 1.73447  -1.10728
 mixrho: warning. negative smrho; isp number min=       1    4456 -0.26757D-03
 mixrho: warning. negative smrho; isp number min=       2    5988 -0.27682D-04
 mixrho: warning. negative smrho; isp number min=       1    4456 -0.26757D-03
 mixrho: warning. negative smrho; isp number min=       2    5988 -0.27682D-04
 mixrho: warning. negative smrho; isp number min=       1    4456 -0.26757D-03
 mixrho: warning. negative smrho; isp number min=       2    5988 -0.27682D-04
 mixrho: add corrections to qcell smrho =  0.93169D-07  0.46584D-10
 mixrho: warning. negative smrho; isp number min=       1    4456 -0.26757D-03
 mixrho: warning. negative smrho; isp number min=       2    5988 -0.27682D-04

 iors  : write rst restart file (binary mesh density)

   it  5  of 10    ehf=     -74.567531   ehk=     -74.548235
 From last iter    ehf=     -74.543498   ehk=     -74.543679
 diffe(q)= -0.024034 (0.037843)    tol= 0.000010 (0.000500)   more=T
i mmom= 1.0000 ehf(eV)=-1014.550914 ehk(eV)=-1014.288376 sev(eV)=-40.650422

--- BNDFP:  begin iteration 6 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.112116
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.700974   charge     8.312583
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       10444  -2.6756678331467474E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.26756678332467472E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.824665 -8.036399 -0.284794
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -2.975211 -3.491720 -0.265747
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       736
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       736
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       736
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       736
 vxcnsp (warning): negative rho: min val =  -1.81E-02
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91971   2.91256   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
     potential shift to crystal energy zero:    0.000008
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -31.506595        17.681306       -13.825288
   rhoval*ves             -3.715190        -6.775215       -10.490406
   psnuc*ves              15.117138      -283.325355      -268.208217
   utot                    5.700974      -145.050285      -139.349311
   rho*exc                -8.799875        -0.230301        -9.030176
   rho*vxc               -11.528118        -0.357507       -11.885625
   valence chg             7.312583        -4.312583         3.000000
   valence mag             2.649782        -1.238315         1.411468
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:      0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     3 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.0369

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.642709;   3.000000 electrons
         Sum occ. bands:   -2.989395, incl. Bloechl correction: -0.000038
         Mag. moment:       1.000000
 bndfp:Generating TDOS: efermi= -0.642709  dos window emin emax=  -1.274698  2.297214

       contr. to mm extrapolated for r>rmt:   0.033894 est. true mm = 0.992425
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.30494  sum tc=    31.40618  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.23891  sum tc=    31.55676  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.907070    8.011035   -5.103965      0.958531    1.732259   -0.773728
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.989395 Vin*nin=     -13.825288 Ek=Eb-Vin*nin=      10.835894
 Ek(core)=      62.979231 Exc=      -9.030176 Ees=    -139.349311 Eharris=     -74.564363

 mkekin:
   nout*Vin = smpart,onsite,total=:    -35.450302     21.713293    -13.737009
    E_B(band energy sum)=   -2.989395  E_B-nout*Vin=   10.747614

 m_mkpot_energyterms
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111321
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.999088   charge     9.103965
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        4454  -2.2817193067546507E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.22817193068546508E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -6.004152 -8.165955 -0.275141
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -4.073889 -5.035914 -0.262707
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91986   2.91398   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -36.418498        22.682895       -13.735603
   rhoval*ves             -3.479127        -6.976607       -10.455734
   psnuc*ves              15.477304      -283.598396      -268.121092
   utot                    5.999088      -145.287502      -139.288413
   rho*exc               -10.078041         1.107423        -8.970618
   rho*vxc               -13.201869         1.396025       -11.805844
   valence chg             8.103965        -5.103965         3.000000
   valence mag             1.773728        -0.773728         1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       10.747614 Ekcore=        62.962936 Ektot    =       73.710551
 Exc=      -8.970618 Ees   =      -139.288413 EKohnSham=      -74.548480
 Magnetic moment=     1.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       10444  -2.6756750756716280E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.26756750757716278E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       10444  -2.6756605906218673E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.26756605907218671E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       10444  -2.6756678331467479E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.26756678332467477E-3

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
 Maximum Harris force =   0.000000 mRy/au (site 1 )
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
 mixrho: sought 2 iter from file mixm ; read 3 RMS DQ= 1.61E-2  last it= 3.78E-2
 AMIX: nmix=2 mmix=8  nelts= 33302  beta=0.50000  tm= 5.00000  rmsdel= 8.03D-03
   tj:-1.70415  -0.01087
 mixrho: warning. negative smrho; isp number min=       1    4432 -0.23676D-03
 mixrho: warning. negative smrho; isp number min=       2    2992 -0.48929D-05
 mixrho: warning. negative smrho; isp number min=       1    4432 -0.23676D-03
 mixrho: warning. negative smrho; isp number min=       2    2992 -0.48929D-05
 mixrho: warning. negative smrho; isp number min=       1    4432 -0.23676D-03
 mixrho: warning. negative smrho; isp number min=       2    2992 -0.48929D-05
 mixrho: add corrections to qcell smrho =  0.10804D-06  0.54019D-10
 mixrho: warning. negative smrho; isp number min=       1    4432 -0.23676D-03
 mixrho: warning. negative smrho; isp number min=       2    2992 -0.48929D-05

 iors  : write rst restart file (binary mesh density)

   it  6  of 10    ehf=     -74.564363   ehk=     -74.548480
 From last iter    ehf=     -74.567531   ehk=     -74.548235
 diffe(q)=  0.003168 (0.016068)    tol= 0.000010 (0.000500)   more=T
i mmom= 1.0000 ehf(eV)=-1014.507804 ehk(eV)=-1014.291713 sev(eV)=-40.673106

--- BNDFP:  begin iteration 7 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111805
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.859789   charge     8.512463
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        7424  -2.3676456146221525E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.23676456147221525E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.785539 -7.955865 -0.277204
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -3.276363 -3.913931 -0.262666
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       320
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       320
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       320
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       320
 vxcnsp (warning): negative rho: min val =  -3.77E-03
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91986   2.91398   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
     potential shift to crystal energy zero:    0.000008
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -32.720462        18.976557       -13.743905
   rhoval*ves             -3.592681        -6.863623       -10.456305
   psnuc*ves              15.312259      -283.445324      -268.133065
   utot                    5.859789      -145.154474      -139.294685
   rho*exc                -9.061901         0.081589        -8.980313
   rho*vxc               -11.869797         0.050947       -11.818850
   valence chg             7.512463        -4.512463         3.000000
   valence mag             2.253283        -1.166432         1.086852
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     3 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.0380

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.629537;   3.000000 electrons
         Sum occ. bands:   -2.992874, incl. Bloechl correction: -0.000038
         Mag. moment:       1.000000
 bndfp:Generating TDOS: efermi= -0.629537  dos window emin emax=  -1.261555  2.310386

       contr. to mm extrapolated for r>rmt:   0.035470 est. true mm = 0.992239
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.31007  sum tc=    31.42363  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.25944  sum tc=    31.53922  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905756    7.792428   -4.886672      0.956769    1.717987   -0.761218
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.992874 Vin*nin=     -13.743905 Ek=Eb-Vin*nin=      10.751031
 Ek(core)=      62.971082 Exc=      -8.980313 Ees=    -139.294685 Eharris=     -74.552884

 mkekin:
   nout*Vin = smpart,onsite,total=:    -34.693111     20.951272    -13.741839
    E_B(band energy sum)=   -2.992874  E_B-nout*Vin=   10.748964

 m_mkpot_energyterms
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111608
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.929878   charge     8.886672
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        4384  -2.2350042179391759E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.22350042180391760E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.783467 -7.867362 -0.274110
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -3.898255 -4.813946 -0.261254
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92012   2.91588   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -35.023131        21.284167       -13.738963
   rhoval*ves             -3.538416        -6.922645       -10.461061
   psnuc*ves              15.398171      -283.519397      -268.121226
   utot                    5.929878      -145.221021      -139.291143
   rho*exc                -9.681722         0.715031        -8.966691
   rho*vxc               -12.681308         0.880597       -11.800711
   valence chg             7.886672        -4.886672         3.000000
   valence mag             1.761218        -0.761218         1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       10.748964 Ekcore=        62.962847 Ektot    =       73.711811
 Exc=      -8.966691 Ees   =      -139.291143 EKohnSham=      -74.546023
 Magnetic moment=     1.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        7424  -2.3676525981708427E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.23676525982708427E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        7424  -2.3676386310734625E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.23676386311734625E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        7424  -2.3676456146221527E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.23676456147221528E-3

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00    0.00   -0.00     0.00    0.00    0.00     0.00   -0.00    0.00
 shift forces to make zero average correction:            0.00   -0.00    0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force =   0.000000 mRy/au (site 1 )
  Symmetrize forces ...
 wgtsmooth=   8.0000000000000002E-003
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 mixrho: sought 2 iter from file mixm ; read 3 RMS DQ= 1.02E-2  last it= 1.61E-2
 AMIX: nmix=2 mmix=8  nelts= 33302  beta=0.50000  tm= 5.00000  rmsdel= 5.08D-03
   tj: 0.85181  -0.53390
 mixrho: warning. negative smrho; isp number min=       1    4428 -0.23220D-03
 mixrho: warning. negative smrho; isp number min=       2    2698 -0.29761D-05
 mixrho: warning. negative smrho; isp number min=       1    4428 -0.23220D-03
 mixrho: warning. negative smrho; isp number min=       1    4428 -0.23220D-03
 mixrho: warning. negative smrho; isp number min=       2    2698 -0.29761D-05
 mixrho: warning. negative smrho; isp number min=       2    2698 -0.29761D-05
 mixrho: add corrections to qcell smrho =  0.10879D-06  0.54396D-10
 mixrho: warning. negative smrho; isp number min=       1    4428 -0.23220D-03
 mixrho: warning. negative smrho; isp number min=       2    2698 -0.29761D-05

 iors  : write rst restart file (binary mesh density)

   it  7  of 10    ehf=     -74.552884   ehk=     -74.546023
 From last iter    ehf=     -74.564363   ehk=     -74.548480
 diffe(q)=  0.011478 (0.010157)    tol= 0.000010 (0.000500)   more=T
i mmom= 1.0000 ehf(eV)=-1014.351633 ehk(eV)=-1014.258274 sev(eV)=-40.720451

--- BNDFP:  begin iteration 8 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111732
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.872950   charge     8.766212
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        7126  -2.3220061901975655E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.23220061902975655E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.817089 -7.947382 -0.276169
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -3.677727 -4.489437 -0.262509
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       256
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       256
 vxcnsp (warning): negative rho: min val =  -2.45E-03
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       256
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       256
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92012   2.91588   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
     potential shift to crystal energy zero:    0.000006
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -34.275940        20.534037       -13.741903
   rhoval*ves             -3.584950        -6.872873       -10.457823
   psnuc*ves              15.330850      -283.456896      -268.126046
   utot                    5.872950      -145.164885      -139.291934
   rho*exc                -9.494816         0.519341        -8.975476
   rho*vxc               -12.436819         0.624412       -11.812407
   valence chg             7.766212        -4.766212         3.000000
   valence mag             1.970141        -0.913546         1.056595
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:      0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     3 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.0376

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.627889;   3.000000 electrons
         Sum occ. bands:   -2.992036, incl. Bloechl correction: -0.000039
         Mag. moment:       1.000000
 bndfp:Generating TDOS: efermi= -0.627889  dos window emin emax=  -1.259913  2.312033

       contr. to mm extrapolated for r>rmt:   0.035508 est. true mm = 0.992208
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.30999  sum tc=    31.42542  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.26085  sum tc=    31.53731  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905317    7.711783   -4.806466      0.956700    1.761064   -0.804364
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.992036 Vin*nin=     -13.741903 Ek=Eb-Vin*nin=      10.749867
 Ek(core)=      62.966964 Exc=      -8.975476 Ees=    -139.291934 Eharris=     -74.550579

 mkekin:
   nout*Vin = smpart,onsite,total=:    -34.442273     20.705375    -13.736897
    E_B(band energy sum)=   -2.992036  E_B-nout*Vin=   10.744861

 m_mkpot_energyterms
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111690
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.910005   charge     8.806466
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        4374  -2.2280020389389943E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.22280020390389943E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.739061 -7.816617 -0.273953
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -3.800587 -4.678291 -0.260943
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92017   2.91618   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -34.518145        20.781973       -13.736172
   rhoval*ves             -3.555389        -6.903806       -10.459195
   psnuc*ves              15.375398      -283.491600      -268.116202
   utot                    5.910005      -145.197703      -139.287699
   rho*exc                -9.539648         0.574105        -8.965543
   rho*vxc               -12.494907         0.695701       -11.799207
   valence chg             7.806466        -4.806466         3.000000
   valence mag             1.804364        -0.804364         1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       10.744861 Ekcore=        62.962723 Ektot    =       73.707585
 Exc=      -8.965543 Ees   =      -139.287699 EKohnSham=      -74.545657
 Magnetic moment=     1.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        7126  -2.3220131045662606E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.23220131046662606E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        7126  -2.3219992758288706E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.23219992759288707E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        7126  -2.3220061901975658E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.23220061902975658E-3

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
   1   -0.00    0.00   -0.00     0.00    0.00    0.00     0.00   -0.00    0.00
 shift forces to make zero average correction:            0.00   -0.00    0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force =   0.000000 mRy/au (site 1 )
  Symmetrize forces ...
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 wgtsmooth=   8.0000000000000002E-003
 mixrho: sought 2 iter from file mixm ; read 3 RMS DQ= 2.05E-3  last it= 1.02E-2
 AMIX: nmix=2 mmix=8  nelts= 33302  beta=0.50000  tm= 5.00000  rmsdel= 1.03D-03
   tj:-0.35784   0.09766
 mixrho: warning. negative smrho; isp number min=       1    4400 -0.22855D-03
 mixrho: warning. negative smrho; isp number min=       2    2218 -0.18698D-05
 mixrho: warning. negative smrho; isp number min=       1    4400 -0.22855D-03
 mixrho: warning. negative smrho; isp number min=       1    4400 -0.22855D-03
 mixrho: warning. negative smrho; isp number min=       2    2218 -0.18698D-05
 mixrho: warning. negative smrho; isp number min=       2    2218 -0.18698D-05
 mixrho: add corrections to qcell smrho =  0.10486D-06  0.52430D-10
 mixrho: warning. negative smrho; isp number min=       1    4400 -0.22855D-03
 mixrho: warning. negative smrho; isp number min=       2    2218 -0.18698D-05

 iors  : write rst restart file (binary mesh density)

   it  8  of 10    ehf=     -74.550579   ehk=     -74.545657
 From last iter    ehf=     -74.552884   ehk=     -74.546023
 diffe(q)=  0.002305 (0.002052)    tol= 0.000010 (0.000500)   more=T
i mmom= 1.0000 ehf(eV)=-1014.320271 ehk(eV)=-1014.253301 sev(eV)=-40.709045

--- BNDFP:  begin iteration 9 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111713
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.886125   charge     8.809766
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        6618  -2.2854772813891183E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.22854772814891183E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.788879 -7.892558 -0.275355
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -3.772875 -4.631727 -0.261917
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       224
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       224
 vxcnsp (warning): negative rho: min val =  -1.72E-03
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       224
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       224
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92017   2.91618   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
     potential shift to crystal energy zero:    0.000006
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -34.541641        20.799776       -13.741865
   rhoval*ves             -3.574667        -6.885203       -10.459871
   psnuc*ves              15.346916      -283.470418      -268.123502
   utot                    5.886125      -145.177811      -139.291686
   rho*exc                -9.561753         0.589650        -8.972103
   rho*vxc               -12.524285         0.716348       -11.807936
   valence chg             7.809766        -4.809766         3.000000
   valence mag             1.876003        -0.835791         1.040212
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:      0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     3 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.0378

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.626936;   3.000000 electrons
         Sum occ. bands:   -2.991602, incl. Bloechl correction: -0.000039
         Mag. moment:       1.000000
 bndfp:Generating TDOS: efermi= -0.626936  dos window emin emax=  -1.258962  2.312986

       contr. to mm extrapolated for r>rmt:   0.035614 est. true mm = 0.992196
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.30989  sum tc=    31.42635  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.26173  sum tc=    31.53615  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905146    7.683444   -4.778299      0.956581    1.770044   -0.813463
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.991602 Vin*nin=     -13.741865 Ek=Eb-Vin*nin=      10.750263
 Ek(core)=      62.964843 Exc=      -8.972103 Ees=    -139.291686 Eharris=     -74.548683

 mkekin:
   nout*Vin = smpart,onsite,total=:    -34.334225     20.599218    -13.735007
    E_B(band energy sum)=   -2.991602  E_B-nout*Vin=   10.743404

 m_mkpot_energyterms
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111719
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.903400   charge     8.778299
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        4374  -2.2295978904267503E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.22295978905267504E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.718682 -7.791033 -0.273991
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -3.770805 -4.638033 -0.260926
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92020   2.91634   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -34.340645        20.605376       -13.735269
   rhoval*ves             -3.561064        -6.897589       -10.458653
   psnuc*ves              15.367864      -283.481912      -268.114048
   utot                    5.903400      -145.189750      -139.286350
   rho*exc                -9.489486         0.524184        -8.965303
   rho*vxc               -12.429065         0.630179       -11.798887
   valence chg             7.778299        -4.778299         3.000000
   valence mag             1.813463        -0.813463         1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       10.743404 Ekcore=        62.962500 Ektot    =       73.705905
 Exc=      -8.965303 Ees   =      -139.286350 EKohnSham=      -74.545748
 Magnetic moment=     1.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        6618  -2.2854841191692938E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.22854841192692938E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        6618  -2.2854704436089425E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.22854704437089425E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        6618  -2.2854772813891180E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.22854772814891180E-3

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
   1   -0.00   -0.00   -0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force =   0.000000 mRy/au (site 1 )
  Symmetrize forces ...
 wgtsmooth=   8.0000000000000002E-003
 mixrho: sought 2 iter from file mixm ; read 3 RMS DQ= 8.80E-4  last it= 2.05E-3
 AMIX: nmix=2 mmix=8  nelts= 33302  beta=0.50000  tm= 5.00000  rmsdel= 4.40D-04
   tj:-2.81612   0.57016
 mixrho: warning. negative smrho; isp number min=       1    4376 -0.22333D-03
 mixrho: warning. negative smrho; isp number min=       2     406 -0.33290D-06
 mixrho: warning. negative smrho; isp number min=       1    4376 -0.22333D-03
 mixrho: warning. negative smrho; isp number min=       2     406 -0.33290D-06
 mixrho: add corrections to qcell smrho =  0.10191D-06  0.50955D-10
 mixrho: warning. negative smrho; isp number min=       1    4376 -0.22333D-03
 mixrho: warning. negative smrho; isp number min=       1    4376 -0.22333D-03
 mixrho: warning. negative smrho; isp number min=       2     406 -0.33290D-06
 mixrho: warning. negative smrho; isp number min=       2     406 -0.33290D-06

 iors  : write rst restart file (binary mesh density)

   it  9  of 10    ehf=     -74.548683   ehk=     -74.545748
 From last iter    ehf=     -74.550579   ehk=     -74.545657
 diffe(q)=  0.001896 (0.000880)    tol= 0.000010 (0.000500)   more=T
i mmom= 1.0000 ehf(eV)=-1014.294471 ehk(eV)=-1014.254541 sev(eV)=-40.703145

--- BNDFP:  begin iteration 10 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111725
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.904058   charge     8.761837
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        4782  -2.2333151725975047E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.22333151726975048E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.702228 -7.768636 -0.274137
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -3.756775 -4.620385 -0.260874
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       128
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       128
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       128
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       128
 vxcnsp (warning): negative rho: min val =  -3.91E-04
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92020   2.91634   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
     potential shift to crystal energy zero:    0.000006
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -34.236482        20.497614       -13.738867
   rhoval*ves             -3.560048        -6.900991       -10.461039
   psnuc*ves              15.368164      -283.487774      -268.119610
   utot                    5.904058      -145.194382      -139.290324
   rho*exc                -9.459003         0.492557        -8.966446
   rho*vxc               -12.389021         0.588598       -11.800423
   valence chg             7.761837        -4.761837         3.000000
   valence mag             1.817655        -0.807320         1.010334
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     3 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.0375

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.625519;   3.000000 electrons
         Sum occ. bands:   -2.991714, incl. Bloechl correction: -0.000039
         Mag. moment:       1.000000
 bndfp:Generating TDOS: efermi= -0.625519  dos window emin emax=  -1.257541  2.314403

       contr. to mm extrapolated for r>rmt:   0.035866 est. true mm = 0.992183
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.31019  sum tc=    31.42791  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.26389  sum tc=    31.53424  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905016    7.665384   -4.760368      0.956317    1.766086   -0.809768
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.991714 Vin*nin=     -13.738867 Ek=Eb-Vin*nin=      10.747153
 Ek(core)=      62.963672 Exc=      -8.966446 Ees=    -139.290324 Eharris=     -74.545946

 mkekin:
   nout*Vin = smpart,onsite,total=:    -34.234261     20.499090    -13.735171
    E_B(band energy sum)=   -2.991714  E_B-nout*Vin=   10.743457

 m_mkpot_energyterms
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst= 0.111742
 cell interaction energy from homogeneous background (q=  1.000000 ) is   0.290159
   smooth rhoves      5.898644   charge     8.760368
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        4366  -2.2289711944693715E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.22289711945693716E-3
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -5.697996 -7.762329 -0.273988
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -3.758784 -4.623762 -0.260872
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92025   2.91655   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff            -34.226477        20.490867       -13.735610
   rhoval*ves             -3.565150        -6.894078       -10.459228
   psnuc*ves              15.362437      -283.475729      -268.113292
   utot                    5.898644      -145.184903      -139.286260
   rho*exc                -9.456780         0.491714        -8.965066
   rho*vxc               -12.386091         0.587518       -11.798573
   valence chg             7.760368        -4.760368         3.000000
   valence mag             1.809768        -0.809768         1.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
   hom background     1.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       10.743457 Ekcore=        62.962151 Ektot    =       73.705608
 Exc=      -8.965066 Ees   =      -139.286260 EKohnSham=      -74.545717
 Magnetic moment=     1.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        4782  -2.2333219118879916E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.22333219119879917E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        4782  -2.2333084333070179E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.22333084334070179E-3
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        4782  -2.2333151725975047E-004
  smvxcm: enforce positive smrho_w. Add srshift= 0.22333151726975048E-3

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
 Maximum Harris force =   0.000000 mRy/au (site 1 )
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
 mixrho: sought 2 iter from file mixm ; read 3 RMS DQ= 8.29E-5  last it= 8.80E-4
 AMIX: nmix=2 mmix=8  nelts= 33302  beta=0.50000  tm= 5.00000  rmsdel= 4.14D-05
   tj:-0.00689   0.00234
 mixrho: warning. negative smrho; isp number min=       1    4374 -0.22311D-03
 mixrho: warning. negative smrho; isp number min=       2     126 -0.14649D-06
 mixrho: warning. negative smrho; isp number min=       1    4374 -0.22311D-03
 mixrho: warning. negative smrho; isp number min=       2     126 -0.14649D-06
 mixrho: warning. negative smrho; isp number min=       1    4374 -0.22311D-03
 mixrho: warning. negative smrho; isp number min=       2     126 -0.14649D-06
 mixrho: add corrections to qcell smrho =  0.10451D-06  0.52253D-10
 mixrho: warning. negative smrho; isp number min=       1    4374 -0.22311D-03
 mixrho: warning. negative smrho; isp number min=       2     126 -0.14649D-06

 iors  : write rst restart file (binary mesh density)

   it 10  of 10    ehf=     -74.545946   ehk=     -74.545717
 From last iter    ehf=     -74.548683   ehk=     -74.545748
 diffe(q)=  0.002737 (0.000083)    tol= 0.000010 (0.000500)   more=F
x mmom= 1.0000 ehf(eV)=-1014.257234 ehk(eV)=-1014.254117 sev(eV)=-40.704665
Exit 0 procid= 0 OK! end of LMF ======================
