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
rval2: BZ_ZBAK                 defa n= 1 val= 0.00000000
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
           3 -vzbak=0
           1 --no-iactiv
           2 c
           3 -vzbak=0
INFO: Ubuntu 20.04.4 LTS \n \l
INFO: GNU Fortran (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0
INFO: -O2 -g -fimplicit-none -finit-integer=NaN -finit-real=NaN -JOBJ.gfortran -IOBJ.gfortran
INFO: MATH: -lmkl_rt
INFO: git: commit 895c0dab24c443e7d7fb788e9870fc04186617c7
           1 --no-iactiv
           2 c
           3 -vzbak=0
INFO:    : Author: Takao Kotani <takaokotani@gmail.com>
INFO:    : Date:   Tue Mar 28 19:00:32 2023 +0900
INFO: linked at Wed Mar 29 10:09:58 JST 2023
===START LMF with   --no-iactiv c -vzbak=0 ===
           1 --no-iactiv
           2 c
           3 -vzbak=0
mpisize=4
m_lmfinit: LMF
cmdl for python=/home/takao/ecalj/SRC/TestInstall/bin/ctrl2ctrlp.py  --no-iactiv c -vzbak=0<ctrl.c >ctrlp.c
rval2: STRUC_NSPEC             requ n= 1 val= 1.00000000
rval2: STRUC_NBAS              requ n= 1 val= 1.00000000
rval2: HAM_NSPIN               defa n= 1 val= 2.00000000
rval2: IO_VERBOS               defa n= 1 val= 30.00000000
rval2: IO_TIM                  defa n= 1 val= 1.00000000
rval2: STRUC_ALAT              ---- n= 1 val= 7.93700526
 cccccccc SPEC_ATOM@1ch=### C###
 cccccccc SPEC_ATOM@1ch=### C###
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
rval2: HAM_FRZWF               defa n= 1 val= 0.00000000
rval2: HAM_XCFUN               defa n= 1 val= 2.00000000
rval2: HAM_FORCES              defa n= 1 val= 12.00000000
rval2: HAM_RDSIG               defa n= 1 val= 1.00000000
rval2: HAM_ScaledSigma         defa n= 1 val= 1.00000000
rval2: HAM_EWALD               defa n= 1 val= 0.00000000
rval2: HAM_OVEPS               defa n= 1 val= 0.00000010
 cccccccc SPEC_C-HOLE@1ch=######
rval2: HAM_PWMODE              defa n= 1 val= 0.00000000
rval2: HAM_PWEMAX              defa n= 1 val= 0.00000000
rval2: HAM_READP               defa n= 1 val= 0.00000000
 cccccccc SITE_ATOM@1ch=### C###
rval2: HAM_V0FIX               defa n= 1 val= 0.00000000
rval2: HAM_PNUFIX              defa n= 1 val= 0.00000000
 cccccccc SPEC_C-HOLE@1ch=######
=== SPEC =1
 cccccccc SPEC_ATOM@1ch=### C###
rval2: SPEC_Z@1                ---- n= 1 val= 6.00000000
 cccccccc SITE_ATOM@1ch=### C###
 cccccccc SPEC_C-HOLE@1ch=######
rval2: SPEC_R@1                ---- n= 1 val= 3.00000000
rval2: SPEC_R/W@1              ---- n= 0 val= 
 cccccccc SITE_ATOM@1ch=### C###
rval2: SPEC_R/A@1              ---- n= 0 val= 
rval2: SPEC_A@1                defa n= 1 val= 0.02000000
rval2: SPEC_NR@1               defa n= 1 val= 0.00000000
rval2: SPEC_RSMH@1             ---- n= 4 val= 1.30000000  1.10000000 -1.00000000 -1.00000000
rval2: SPEC_EH@1               requ n= 2 val= -0.70000000 -0.20000000
rval2: SPEC_RSMH2@1            ---- n= 2 val= 0.80000000  0.80000000
rval2: SPEC_EH2@1              requ n= 2 val= -1.50000000 -1.00000000
rval2: SPEC_LMX@1              defa n= 1 val= 999.00000000
 cccccccc SYMGRPch=###  find###
rval2: SPEC_LMXA@1             defa n= 1 val= 3.00000000
 cccccccc SYMGRPAFch=######
 cccccccc SYMGRPch=###  find###
 cccccccc SYMGRPAFch=######
rval2: SPEC_LMXL@1             defa n= 1 val= 3.00000000
rval2: SPEC_P@1                ---- n= 4 val= 2.90000000  2.85000000  3.18000000  4.12000000
 cccccccc ITER_MIXch=### A2###
 cccccccc SYMGRPch=###  find###
rval2: SPEC_Q@1                ---- n= 0 val= 
rval2: SPEC_MMOM@1             ---- n= 2 val= 0.00000000  2.00000000
 cccccccc ITER_MIXch=### A2###
 cccccccc SYMGRPAFch=######
rval2: SPEC_NMCORE@1           defa n= 1 val= 0.00000000
rval2: SPEC_PZ@1               ---- n= 0 val= 
rval2: SPEC_LFOCA@1            defa n= 1 val= 0.00000000
rval2: SPEC_KMXA@1             defa n= 1 val= 3.00000000
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
rval2: BZ_ZBAK                 defa n= 1 val= 0.00000000
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

                Plat                                  Qlat
mto   eh1    1 -0.70 -0.20
mto rsmh2    1  0.80  0.80
mto  eh2     1 -1.50 -1.00
mto lh       1  1

                Plat                                  Qlat
   1.000000   1.000000   0.000000        0.500000   0.500000  -0.500000
   1.000000   0.000000   1.000000        0.500000  -0.500000   0.500000

                Plat                                  Qlat
   1.000000   1.000000   0.000000        0.500000   0.500000  -0.500000
   1.000000   0.000000   1.000000        0.500000  -0.500000   0.500000
   0.000000   1.000000   1.000000       -0.500000   0.500000   0.500000

                Plat                                  Qlat
   0.000000   1.000000   1.000000       -0.500000   0.500000   0.500000
   1.000000   1.000000   0.000000        0.500000   0.500000  -0.500000
   1.000000   0.000000   1.000000        0.500000  -0.500000   0.500000
   0.000000   1.000000   1.000000       -0.500000   0.500000   0.500000
   1.000000   1.000000   0.000000        0.500000   0.500000  -0.500000
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

LATTC:  as= 2.000   tol= 1.00E-08   alat= 7.93701   awald= 0.200
         r1=  3.459   nkd=  79      q1=  2.571   nkq= 137
SpaceGroupSym: ======================================= 
  Generators except find=
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
     potential shift to crystal energy zero:    0.000002
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -3.931616       -10.430886       -14.362502
   rhoval*ves             -5.058556        -5.181775       -10.240330
   psnuc*ves               9.186377      -278.836573      -269.650196
   utot                    2.063910      -142.009174      -139.945263
   rho*exc                -1.494901        -8.107094        -9.601995
   rho*vxc                -1.946550       -10.686870       -12.633421
   valence chg             2.661718         1.338282         4.000000
   valence mag             1.332841         0.667159         2.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
   hom background     0.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     3 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.0184

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.430305;   4.000000 electrons
         Sum occ. bands:   -2.743701, incl. Bloechl correction: -0.000404
         Mag. moment:       2.000000
 bndfp:Generating TDOS: efermi= -0.430305  dos window emin emax=  -1.050710  2.509617


 m_bandcal_2nd: to fill eigenfunctions**2 up to Efermi
       contr. to mm extrapolated for r>rmt:   0.164809 est. true mm = 1.960686
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85523  sum tc=    31.38701  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78555  sum tc=    31.54499  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.686108    3.832377   -0.146269      1.795877    2.125262   -0.329386
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.743701 Vin*nin=     -14.362502 Ek=Eb-Vin*nin=      11.618800
 Ek(core)=      62.932005 Exc=      -9.601995 Ees=    -139.945263 Eharris=     -74.996454

 mkekin:
   nout*Vin = smpart,onsite,total=:     -7.479092     -6.863479    -14.342571
    E_B(band energy sum)=   -2.743701  E_B-nout*Vin=   11.598869

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.008292
   smooth rhoves      3.523650   charge     4.146269
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=           2  -1.4980078303893662E-005
  smvxcm: enforce positive smrho_w. Add srshift= 0.14980078313893662E-4
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.406849 -3.324277 -0.245380
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -0.672848 -0.698584 -0.192706
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91337   2.90760   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10669   4.12000   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -9.111301        -5.245867       -14.357167
   rhoval*ves             -4.666655        -5.580793       -10.247448
   psnuc*ves              11.713954      -281.339304      -269.625350
   utot                    3.523650      -143.460049      -139.936399
   rho*exc                -3.079697        -6.517779        -9.597477
   rho*vxc                -4.022861        -8.604560       -12.627421
   valence chg             4.146269        -0.146269         4.000000
   valence mag             2.329386        -0.329386         2.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
   hom background     0.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       11.598869 Ekcore=        62.932006 Ektot    =       74.530875
 Exc=      -9.597477 Ees   =      -139.936399 EKohnSham=      -75.003001
 Magnetic moment=     2.000000
  smvxcm: all smrho_w is positive
  smvxcm: all smrho_w is positive
  smvxcm: all smrho_w is positive
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00    -0.00   -0.00   -0.00
 shift forces to make zero average correction:           -0.00   -0.00   -0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force =   0.000000 mRy/au (site 1 )
  Symmetrize forces ...
 mixrealsmooth= T
 wgtsmooth=   8.0000000000000002E-003
 mixrho: sought 2 iter from file mixm ; read 0 RMS DQ= 2.25E-2
 AMIX: nmix=0 mmix=8  nelts= 33302  beta=0.50000  tm= 5.00000  rmsdel= 1.12D-02
 mixrho: warning. negative smrho; isp number min=       1       2 -0.57849D-05
 mixrho: add corrections to qcell smrho =  0.78196D-07  0.39098D-10
 mixrho: warning. negative smrho; isp number min=       1       2 -0.57849D-05
 mixrho: warning. negative smrho; isp number min=       1       2 -0.57849D-05

 iors  : write rst restart file (binary mesh density)
 mixrho: warning. negative smrho; isp number min=       1       2 -0.57849D-05

   it  1  of 10    ehf=     -74.996454   ehk=     -75.003001
h mmom= 2.0000 ehf(eV)=-1020.386748 ehk(eV)=-1020.475827 sev(eV)=-37.330254

--- BNDFP:  begin iteration 2 of 10
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.007450
   smooth rhoves      2.741764   charge     3.403994
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=           2  -5.7848615283650336E-006
  smvxcm: enforce positive smrho_w. Add srshift= 0.57848615383650334E-5
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -1.719970 -2.370560 -0.232668
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -0.523729 -0.556409 -0.179645
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91337   2.90760   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10669   4.12000   0.00000   0.00000
     potential shift to crystal energy zero:    0.000003
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -6.334130        -8.025471       -14.359601
   rhoval*ves             -4.966638        -5.277288       -10.243926
   psnuc*ves              10.450165      -280.087939      -269.637773
   utot                    2.741764      -142.682613      -139.940850
   rho*exc                -2.243699        -7.355137        -9.598836
   rho*vxc                -2.926969        -9.702285       -12.629254
   valence chg             3.403994         0.596006         4.000000
   valence mag             1.831113         0.168887         2.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
   hom background     0.00000   deviation from neutrality:      0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     3 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.0338

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.431869;   4.000000 electrons
         Sum occ. bands:   -2.750893, incl. Bloechl correction: -0.000417
         Mag. moment:       2.000000
 bndfp:Generating TDOS: efermi= -0.431869  dos window emin emax=  -1.052404  2.508053

       contr. to mm extrapolated for r>rmt:   0.164886 est. true mm = 1.960320
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85853  sum tc=    31.38725  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78893  sum tc=    31.54500  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.686333    3.844200   -0.157867      1.795434    2.110399   -0.314965
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.750893 Vin*nin=     -14.359601 Ek=Eb-Vin*nin=      11.608708
 Ek(core)=      62.932005 Exc=      -9.598836 Ees=    -139.940850 Eharris=     -74.998973

 mkekin:
   nout*Vin = smpart,onsite,total=:     -8.391033     -5.971242    -14.362275
    E_B(band energy sum)=   -2.750893  E_B-nout*Vin=   11.611382

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.008205
   smooth rhoves      3.510925   charge     4.157867
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=           2  -1.2713946590290186E-005
  smvxcm: enforce positive smrho_w. Add srshift= 0.12713946600290186E-4
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.407747 -3.326114 -0.244399
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -0.685472 -0.714307 -0.190767
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91339   2.90758   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10659   4.12000   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -9.167601        -5.198194       -14.365795
   rhoval*ves             -4.671225        -5.583285       -10.254511
   psnuc*ves              11.693076      -281.333696      -269.640620
   utot                    3.510925      -143.458491      -139.947565
   rho*exc                -3.093219        -6.504829        -9.598048
   rho*vxc                -4.040421        -8.587774       -12.628195
   valence chg             4.157867        -0.157867         4.000000
   valence mag             2.314965        -0.314965         2.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
   hom background     0.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       11.611382 Ekcore=        62.932247 Ektot    =       74.543629
 Exc=      -9.598048 Ees   =      -139.947565 EKohnSham=      -75.001984
 Magnetic moment=     2.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=           2  -5.7848784936464802E-006
  smvxcm: enforce positive smrho_w. Add srshift= 0.57848785036464800E-5
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=           2  -5.7848445630835878E-006
  smvxcm: enforce positive smrho_w. Add srshift= 0.57848445730835876E-5
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=           2  -5.7848615283650344E-006
  smvxcm: enforce positive smrho_w. Add srshift= 0.57848615383650343E-5
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00   -0.00    0.00     0.00    0.00    0.00    -0.00    0.00   -0.00
 shift forces to make zero average correction:           -0.00    0.00   -0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force =   0.000000 mRy/au (site 1 )
  Symmetrize forces ...
 wgtsmooth=   8.0000000000000002E-003
 mixrho: sought 2 iter from file mixm ; read 1 RMS DQ= 1.14E-2  last it= 2.25E-2
 AMIX: nmix=1 mmix=8  nelts= 33302  beta=0.50000  tm= 5.00000  rmsdel= 5.68D-03
   tj:-1.01724
 mixrho: warning. negative smrho; isp number min=       1       2 -0.12774D-04
 mixrho: warning. negative smrho; isp number min=       1       2 -0.12774D-04
 mixrho: add corrections to qcell smrho =  0.12107D-07  0.60533D-11
 mixrho: warning. negative smrho; isp number min=       1       2 -0.12774D-04
 mixrho: warning. negative smrho; isp number min=       1       2 -0.12774D-04

 iors  : write rst restart file (binary mesh density)

   it  2  of 10    ehf=     -74.998973   ehk=     -75.001984
 From last iter    ehf=     -74.996454   ehk=     -75.003001
 diffe(q)= -0.002520 (0.011360)    tol= 0.000010 (0.000500)   more=T
i mmom= 2.0000 ehf(eV)=-1020.421031 ehk(eV)=-1020.462000 sev(eV)=-37.428105

--- BNDFP:  begin iteration 3 of 10
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.008212
   smooth rhoves      3.517994   charge     4.164365
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=           2  -1.2773626699674988E-005
  smvxcm: enforce positive smrho_w. Add srshift= 0.12773626709674988E-4
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.413953 -3.334741 -0.244490
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -0.686917 -0.715712 -0.190850
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91339   2.90758   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10659   4.12000   0.00000   0.00000
     potential shift to crystal energy zero:    0.000004
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -9.193617        -5.172234       -14.365851
   rhoval*ves             -4.667801        -5.586802       -10.254603
   psnuc*ves              11.703789      -281.344201      -269.640412
   utot                    3.517994      -143.465502      -139.947508
   rho*exc                -3.100870        -6.497165        -9.598035
   rho*vxc                -4.050453        -8.577724       -12.628178
   valence chg             4.164365        -0.164365         4.000000
   valence mag             2.319136        -0.319136         2.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
   hom background     0.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     3 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.0348

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.432519;   4.000000 electrons
         Sum occ. bands:   -2.754425, incl. Bloechl correction: -0.000428
         Mag. moment:       2.000000
 bndfp:Generating TDOS: efermi= -0.432519  dos window emin emax=  -1.053150  2.507404

       contr. to mm extrapolated for r>rmt:   0.165232 est. true mm = 1.959950
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85917  sum tc=    31.38685  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78959  sum tc=    31.54442  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.686125    3.852689   -0.166564      1.794718    2.094822   -0.300103
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.754425 Vin*nin=     -14.365851 Ek=Eb-Vin*nin=      11.611426
 Ek(core)=      62.932126 Exc=      -9.598035 Ees=    -139.947508 Eharris=     -75.001990

 mkekin:
   nout*Vin = smpart,onsite,total=:     -9.220430     -5.150257    -14.370687
    E_B(band energy sum)=   -2.754425  E_B-nout*Vin=   11.616262

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.008101
   smooth rhoves      3.500087   charge     4.166564
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=           2  -1.0629546775462664E-005
  smvxcm: enforce positive smrho_w. Add srshift= 0.10629546785462664E-4
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.405804 -3.323868 -0.243540
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -0.696625 -0.728465 -0.188889
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91338   2.90753   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10651   4.12000   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -9.210169        -5.159085       -14.369254
   rhoval*ves             -4.676214        -5.581523       -10.257737
   psnuc*ves              11.676389      -281.320814      -269.644425
   utot                    3.500087      -143.451168      -139.951081
   rho*exc                -3.102429        -6.495081        -9.597510
   rho*vxc                -4.052334        -8.575170       -12.627504
   valence chg             4.166564        -0.166564         4.000000
   valence mag             2.300103        -0.300103         2.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
   hom background     0.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       11.616262 Ekcore=        62.931269 Ektot    =       74.547531
 Exc=      -9.597510 Ees   =      -139.951081 EKohnSham=      -75.001060
 Magnetic moment=     2.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=           2  -1.2773664724368281E-005
  smvxcm: enforce positive smrho_w. Add srshift= 0.12773664734368281E-4
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=           2  -1.2773588674981695E-005
  smvxcm: enforce positive smrho_w. Add srshift= 0.12773588684981695E-4
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=           2  -1.2773626699674988E-005
  smvxcm: enforce positive smrho_w. Add srshift= 0.12773626709674988E-4

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
 Maximum Harris force =   0.000000 mRy/au (site 1 )
   1    0.00    0.00   -0.00     0.00    0.00    0.00    -0.00   -0.00    0.00
 shift forces to make zero average correction:           -0.00   -0.00    0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force =   0.000000 mRy/au (site 1 )
  Symmetrize forces ...
 wgtsmooth=   8.0000000000000002E-003
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 mixrho: sought 2 iter from file mixm ; read 2 RMS DQ= 3.11E-4  last it= 1.14E-2
 AMIX: nmix=2 mmix=8  nelts= 33302  beta=0.50000  tm= 5.00000  rmsdel= 1.55D-04
   tj:-1.59829   0.81014
 mixrho: warning. negative smrho; isp number min=       1       2 -0.10827D-04
 mixrho: warning. negative smrho; isp number min=       1       2 -0.10827D-04
 mixrho: warning. negative smrho; isp number min=       1       2 -0.10827D-04
 mixrho: add corrections to qcell smrho =  0.71256D-07  0.35628D-10
 mixrho: warning. negative smrho; isp number min=       1       2 -0.10827D-04

 iors  : write rst restart file (binary mesh density)

   it  3  of 10    ehf=     -75.001990   ehk=     -75.001060
 From last iter    ehf=     -74.998973   ehk=     -75.001984
 diffe(q)= -0.003017 (0.000311)    tol= 0.000010 (0.000500)   more=T
i mmom= 2.0000 ehf(eV)=-1020.462080 ehk(eV)=-1020.449421 sev(eV)=-37.476153

--- BNDFP:  begin iteration 4 of 10
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.008109
   smooth rhoves      3.498542   charge     4.163162
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=           2  -1.0827496532760066E-005
  smvxcm: enforce positive smrho_w. Add srshift= 0.10827496542760066E-4
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.403639 -3.320813 -0.243597
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -0.694888 -0.726421 -0.189062
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91338   2.90753   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10651   4.12000   0.00000   0.00000
     potential shift to crystal energy zero:    0.000004
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -9.195725        -5.173140       -14.368864
   rhoval*ves             -4.676983        -5.580375       -10.257358
   psnuc*ves              11.674067      -281.318700      -269.644633
   utot                    3.498542      -143.449538      -139.950996
   rho*exc                -3.098527        -6.499059        -9.597585
   rho*vxc                -4.047234        -8.580368       -12.627602
   valence chg             4.163162        -0.163162         4.000000
   valence mag             2.300085        -0.300085         2.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
   hom background     0.00000   deviation from neutrality:      0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     3 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.0347

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.432187;   4.000000 electrons
         Sum occ. bands:   -2.753129, incl. Bloechl correction: -0.000427
         Mag. moment:       2.000000
 bndfp:Generating TDOS: efermi= -0.432187  dos window emin emax=  -1.052820  2.507736

       contr. to mm extrapolated for r>rmt:   0.165263 est. true mm = 1.959979
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85838  sum tc=    31.38668  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78881  sum tc=    31.54422  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.686046    3.847590   -0.161543      1.794716    2.093783   -0.299067
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.753129 Vin*nin=     -14.368864 Ek=Eb-Vin*nin=      11.615735
 Ek(core)=      62.931697 Exc=      -9.597585 Ees=    -139.950996 Eharris=     -75.001149

 mkekin:
   nout*Vin = smpart,onsite,total=:     -9.187782     -5.179653    -14.367435
    E_B(band energy sum)=   -2.753129  E_B-nout*Vin=   11.614306

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.008097
   smooth rhoves      3.499314   charge     4.161543
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=           2  -1.0783313901629765E-005
  smvxcm: enforce positive smrho_w. Add srshift= 0.10783313911629765E-4
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.401639 -3.318024 -0.243606
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -0.694423 -0.725973 -0.189033
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91340   2.90757   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10652   4.12000   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -9.187738        -5.180152       -14.367890
   rhoval*ves             -4.676881        -5.579782       -10.256663
   psnuc*ves              11.675509      -281.316930      -269.641421
   utot                    3.499314      -143.448356      -139.949042
   rho*exc                -3.096062        -6.501235        -9.597297
   rho*vxc                -4.043997        -8.583224       -12.627221
   valence chg             4.161543        -0.161543         4.000000
   valence mag             2.299067        -0.299067         2.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
   hom background     0.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       11.614306 Ekcore=        62.930904 Ektot    =       74.545209
 Exc=      -9.597297 Ees   =      -139.949042 EKohnSham=      -75.001129
 Magnetic moment=     2.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=           2  -1.0827528671391469E-005
  smvxcm: enforce positive smrho_w. Add srshift= 0.10827528681391469E-4
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=           2  -1.0827464394128663E-005
  smvxcm: enforce positive smrho_w. Add srshift= 0.10827464404128663E-4
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=           2  -1.0827496532760066E-005
  smvxcm: enforce positive smrho_w. Add srshift= 0.10827496542760066E-4

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
 wgtsmooth=   8.0000000000000002E-003
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 mixrho: sought 2 iter from file mixm ; read 3 RMS DQ= 4.23E-5  last it= 3.11E-4
 AMIX: nmix=2 mmix=8  nelts= 33302  beta=0.50000  tm= 5.00000  rmsdel= 2.12D-05
   tj: 0.04255   0.00318
 mixrho: warning. negative smrho; isp number min=       1       2 -0.10839D-04
 mixrho: warning. negative smrho; isp number min=       1       2 -0.10839D-04
 mixrho: warning. negative smrho; isp number min=       1       2 -0.10839D-04
 mixrho: add corrections to qcell smrho =  0.45001D-07  0.22500D-10
 mixrho: warning. negative smrho; isp number min=       1       2 -0.10839D-04

 iors  : write rst restart file (binary mesh density)

   it  4  of 10    ehf=     -75.001149   ehk=     -75.001129
 From last iter    ehf=     -75.001990   ehk=     -75.001060
 diffe(q)=  0.000842 (0.000042)    tol= 0.000010 (0.000500)   more=T
i mmom= 2.0000 ehf(eV)=-1020.450628 ehk(eV)=-1020.450367 sev(eV)=-37.458524

--- BNDFP:  begin iteration 5 of 10
 m_mkpot_init: Making one-particle potential ...
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.008105
   smooth rhoves      3.498096   charge     4.161273
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=           2  -1.0838567783782281E-005
  smvxcm: enforce positive smrho_w. Add srshift= 0.10838567793782281E-4
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.401820 -3.318271 -0.243605
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -0.694236 -0.725724 -0.189073
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91340   2.90757   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10652   4.12000   0.00000   0.00000
     potential shift to crystal energy zero:    0.000004
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -9.187293        -5.181031       -14.368324
   rhoval*ves             -4.677331        -5.579620       -10.256950
   psnuc*ves              11.673522      -281.316474      -269.642952
   utot                    3.498096      -143.448047      -139.949951
   rho*exc                -3.096055        -6.501402        -9.597457
   rho*vxc                -4.043994        -8.583438       -12.627432
   valence chg             4.161273        -0.161273         4.000000
   valence mag             2.299284        -0.299284         2.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
   hom background     0.00000   deviation from neutrality:     -0.00000
 m_bandcal_init: start
 bndfp: kpt     1 of     3 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.0350

 bzwts: --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.432246;   4.000000 electrons
         Sum occ. bands:   -2.753369, incl. Bloechl correction: -0.000427
         Mag. moment:       2.000000
 bndfp:Generating TDOS: efermi= -0.432246  dos window emin emax=  -1.052882  2.507677

       contr. to mm extrapolated for r>rmt:   0.165235 est. true mm = 1.959986
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85858  sum tc=    31.38675  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78901  sum tc=    31.54428  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.686096    3.848339   -0.162244      1.794751    2.094071   -0.299320
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702
 Eb(band sum)=       -2.753369 Vin*nin=     -14.368324 Ek=Eb-Vin*nin=      11.614955
 Ek(core)=      62.931301 Exc=      -9.597457 Ees=    -139.949951 Eharris=     -75.001153

 mkekin:
   nout*Vin = smpart,onsite,total=:     -9.189597     -5.178669    -14.368266
    E_B(band energy sum)=   -2.753369  E_B-nout*Vin=   11.614897

 m_mkpot_energyterms
  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=-0.008102
   smooth rhoves      3.499579   charge     4.162244
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=           2  -1.0812364788763121E-005
  smvxcm: enforce positive smrho_w. Add srshift= 0.10812364798763121E-4
  smvxc2: smooth isp rhoeps rhomu vxcavg= 1 -2.402299 -3.318948 -0.243612
  smvxc2: smooth isp rhoeps rhomu vxcavg= 2 -0.694700 -0.726277 -0.189057
  locpot:
   site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91340   2.90757   0.00000   0.00000
 ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
 ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
 ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10652   4.12000   0.00000   0.00000
  mkpot:
   Energy terms:           smooth           local           total
   rhoval*veff             -9.190744        -5.177530       -14.368274
   rhoval*ves             -4.676617        -5.580309       -10.256926
   psnuc*ves              11.675775      -281.318146      -269.642372
   utot                    3.499579      -143.449228      -139.949649
   rho*exc                -3.096999        -6.500414        -9.597413
   rho*vxc                -4.045225        -8.582149       -12.627374
   valence chg             4.162244        -0.162244         4.000000
   valence mag             2.299320        -0.299320         2.000000
   core charge             2.000000         0.000000         2.000000
   Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
   hom background     0.00000   deviation from neutrality:     -0.00000

 m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout
 Ek=       11.614897 Ekcore=        62.931024 Ektot    =       74.545920
 Exc=      -9.597413 Ees   =      -139.949649 EKohnSham=      -75.001142
 Magnetic moment=     2.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=           2  -1.0838599955710031E-005
  smvxcm: enforce positive smrho_w. Add srshift= 0.10838599965710031E-4
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=           2  -1.0838535611854531E-005
  smvxcm: enforce positive smrho_w. Add srshift= 0.10838535621854531E-4
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=           2  -1.0838567783782281E-005
  smvxcm: enforce positive smrho_w. Add srshift= 0.10838567793782281E-4

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
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 Maximum Harris force =   0.000000 mRy/au (site 1 )
 mixrho: sought 2 iter from file mixm ; read 3 RMS DQ= 1.22E-5  last it= 4.23E-5
 AMIX: nmix=2 mmix=8  nelts= 33302  beta=0.50000  tm= 5.00000  rmsdel= 6.09D-06
   tj: 0.13451  -0.00057
 mixrho: warning. negative smrho; isp number min=       1       2 -0.10822D-04
 mixrho: warning. negative smrho; isp number min=       1       2 -0.10822D-04
 mixrho: warning. negative smrho; isp number min=       1       2 -0.10822D-04
 mixrho: add corrections to qcell smrho =  0.45046D-07  0.22523D-10
 mixrho: warning. negative smrho; isp number min=       1       2 -0.10822D-04

 iors  : write rst restart file (binary mesh density)

   it  5  of 10    ehf=     -75.001153   ehk=     -75.001142
 From last iter    ehf=     -75.001149   ehk=     -75.001129
 diffe(q)= -0.000004 (0.000012)    tol= 0.000010 (0.000500)   more=F
c mmom= 2.0000 ehf(eV)=-1020.450685 ehk(eV)=-1020.450531 sev(eV)=-37.461791
Exit 0 procid= 0 OK! end of LMF ======================
