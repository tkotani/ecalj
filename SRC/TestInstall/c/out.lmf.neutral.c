INFO: Ubuntu 20.04.4 LTS \n \l
INFO: GNU Fortran (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0
INFO: -O2 -g -fimplicit-none -ffixed-line-length-132 -JOBJ.gfortran -IOBJ.gfortran
INFO: MATH: -lmkl_rt
INFO: git: commit 6bfdb2ee2eaa61e6c5ffd2b51e721f524be30478
INFO:    : Author: Takao Kotani <takaokotani@gmail.com>
INFO:    : Date:   Sun Jun 19 17:13:24 2022 +0900
INFO: linked at Sun Jun 19 18:36:33 JST 2022
=== START LFMA ===
 mpisize=           1
  Token           Input   cast  (size,min,read,def)     result
 IO_VERBOS         opt    i4       1,  1,   1,  0       30
 >> level: 1  CPUsec=      0.00  enter m_lmfinit
 CONST             opt    chr      1,  0,   1           ef0=0 nk=4 rsma=.6 lmxa=3  nit=10 hf=f zbak=0
  Var       Name                 Val
   1        t                   1.0000    
   2        f                   0.0000    
   3        pi                  3.1416    
   4        zbak                0.0000    
   5        ef0                 0.0000    
   6        nk                  4.0000    
   7        rsma               0.60000    
   8        lmxa                3.0000    
   9        nit                 10.000    
  10        hf                  0.0000    
 STRUC_ALAT        reqd   r8       1,  1,   1, --       7.9370052598409968
 STRUC_NBAS        reqd   i4       1,  1,   1, --       1
 STRUC_PLAT        reqd   r8v      9,  9,   9, --       1  1  0  1  0  1  0  1  1
 ... found 1 species in SPEC category
  --- Program Options ---
 OPTIONS_HF        opt    lg       1,  1,   1,  0       F
 HAM_NSPIN         opt    i4       1,  1,   1,  0       2
 HAM_REL           opt    i4       1,  1,   1,  0       1
 HAM_XCFUN         opt    i4       1,  1,   1,  0       2
 ===> for --jobgw, pwmode is switched to be  0
 SYMGRP            opt    chr      1,  0,   1           find
  --- Parameters for species data ---
 ... Species  1
 SPEC_ATOM         reqd   chr      1,  0,   1           C
 SPEC_ATOM_Z       reqd   r8       1,  1,   1, --       6
 SPEC_ATOM_R       reqd   r8       1,  1,   1, --       3
 SPEC_ATOM_A       opt    r8       1,  1,   1,  0       0.02
 SPEC_ATOM_RSMH    opt    r8v     10,  1,   4,  6       1.3  1.1000000000000001 -1 -1  0  0  0  0  0  0
 SPEC_ATOM_EH      opt    r8v     10,  4,   4,  6       -0.7 -0.2  0  0  0  0  0  0  0  0
 SPEC_ATOM_RSMH2   opt    r8v     10,  1,   2,  8       0.8  0.8  0  0  0  0  0  0  0  0
 SPEC_ATOM_EH2     reqd   r8v     10,  2,   2, --       -1.5 -1
 SPEC_ATOM_LMXA    opt    i4       1,  1,   1,  0       3
 SPEC_ATOM_P       opt    r8v      4,  1,   4,  0       2.8999999999999999  2.8500000000000001  3.1800000000000002  4.1200000000000001
 SPEC_ATOM_MMOM    opt    r8v      4,  1,   2,  2       0  2  0  0
 SPEC_ATOM_IDMOD   opt    i4v      4,  1,   2,  2       0  1  0  0
 SPEC_ATOM_C-HOLE  opt    chr      1,  0,   *
 SPEC_ATOM_EREF    opt    r8       1,  1,   1,  0       -74.9949000000000012
  --- Parameters for site ---
 ... Site  1
 SITE_ATOM         reqd   chr      1,  0,   1           C
 SITE_ATOM_POS     reqd   r8v      3,  1,   3, --       0  0  0
  --- Parameters for Ewald sums ---
  --- Parameters for iterations ---


mmm === MTO setting ===
mmm ispec lmxb lpzex nkapii nkaphh=    1    3    0    2    2
mmm rsmh1    1  1.30  1.10
mmm   eh1    1 -0.70 -0.20
mmm rsmh2    1  0.80  0.80
mmm  eh2     1 -1.50 -1.00
mmm lh       1  1
 >>      0.03   exit  m_lmfinit       0.03
 goto freats

ttt: pnu qat=  1  0     2.900     2.000
ttt: pnu qat=  1  1     2.850     2.000
ttt: pnu qat=  1  2     3.180     0.000
ttt: pnu qat=  1  3     4.120     0.000
ttt: pnu qat=  2  0     2.900     0.000
ttt: pnu qat=  2  1     2.850     2.000
ttt: pnu qat=  2  2     3.180     0.000
ttt: pnu qat=  2  3     4.120     0.000

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
conf:   rmt=3.000000  rmax=19.671121  a=0.02  nr=369  nr(rmax)=463
 goto atomc xxx
 atomsc nmcore=           0

  iter     qint         drho          vh0          rho0          vsum     beta
    1    6.000000   5.461E+02       30.0000    0.2984E+02      -12.0633   0.30
   50    6.000000   3.935E-05       29.2312    0.1279E+03      -59.7470   0.30

 end of atomsc xxxxx
 vsum=  -59.746958882741843                1
sumev= -2.876387 etot= -74.994908 eref= -74.994900 etot-eref= -0.000008
 Optimise free-atom basis for species C        Rmt=  3.000000
 l  it    Rsm      Eh     stiffR   stiffE      Eval      Exact     Pnu    Ql
 0   3   1.500  -0.941       0.0     29.9   -1.07382  -1.07383    2.91   1.00
 1   8   1.500  -0.376       0.0    134.4   -0.46556  -0.46569    2.89   2.00
 eigenvalue sum:  exact  -2.00520    opt basis  -2.00495    error 0.00025
 Optimise free-atom basis for species C        Rmt=  3.000000
 l  it    Rsm      Eh     stiffR   stiffE      Eval      Exact     Pnu    Ql
 0   4   1.500  -0.769       0.0     50.2   -0.87118  -0.87119    2.91   1.00
 1   9   1.500  -0.211       0.0    359.8   -0.27748  -0.27759    2.87   0.00
 eigenvalue sum:  exact  -0.87119    opt basis  -0.87118    error 0.00001
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
 tailsm: end
conf: Core rhoc(rmt)= 0.000000 spillout= 0.000000
 end of freats: spid nmcore=C                  0
OK! end of LMFA ======================
INFO: Ubuntu 20.04.4 LTS \n \l
INFO: GNU Fortran (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0
INFO: -O2 -g -fimplicit-none -ffixed-line-length-132 -JOBJ.gfortran -IOBJ.gfortran
INFO: MATH: -lmkl_rt
INFO: git: commit 6bfdb2ee2eaa61e6c5ffd2b51e721f524be30478
INFO:    : Author: Takao Kotani <takaokotani@gmail.com>
INFO:    : Date:   Sun Jun 19 17:13:24 2022 +0900
INFO: linked at Sun Jun 19 18:36:33 JST 2022
===START LMF with   ===
 mpisize=           4
  Token           Input   cast  (size,min,read,def)     result
 IO_VERBOS         opt    i4       1,  1,   1,  0       30
 >> level: 1  CPUsec=      0.00  enter m_lmfinit
 CONST             opt    chr      1,  0,   1           ef0=0 nk=4 rsma=.6 lmxa=3  nit=10 hf=f zbak=0
  Var       Name                 Val
   1        t                   1.0000    
   2        f                   0.0000    
   3        pi                  3.1416    
   4        zbak                0.0000    
   5        ef0                 0.0000    
   6        nk                  4.0000    
   7        rsma               0.60000    
   8        lmxa                3.0000    
   9        nit                 10.000    
  10        hf                  0.0000    
 STRUC_ALAT        reqd   r8       1,  1,   1, --       7.9370052598409968
 STRUC_NBAS        reqd   i4       1,  1,   1, --       1
 STRUC_PLAT        reqd   r8v      9,  9,   9, --       1  1  0  1  0  1  0  1  1
 ... found 1 species in SPEC category
  --- Program Options ---
 OPTIONS_HF        opt    lg       1,  1,   1,  0       F
 HAM_NSPIN         opt    i4       1,  1,   1,  0       2
 HAM_REL           opt    i4       1,  1,   1,  0       1
 HAM_FTMESH        reqd   i4v      3,  1,   3, --       50  50  50
 HAM_TOL           opt    r8       1,  1,   1,  0       0.000001
 HAM_FORCES        opt    i4       1,  1,   1,  0       12
 HAM_XCFUN         opt    i4       1,  1,   1,  0       2
 ===> for --jobgw, pwmode is switched to be  0
 SYMGRP            opt    chr      1,  0,   1           find
 SYMGRPAF          opt    chr      1,  0,   *
  --- Parameters for species data ---
 ... Species  1
 SPEC_ATOM         reqd   chr      1,  0,   1           C
 SPEC_ATOM_Z       reqd   r8       1,  1,   1, --       6
 SPEC_ATOM_R       reqd   r8       1,  1,   1, --       3
 SPEC_ATOM_A       opt    r8       1,  1,   1,  0       0.02
 SPEC_ATOM_RSMH    reqd   r8v     10,  1,   4,  6       1.3  1.1000000000000001 -1 -1  0  0  0  0  0  0
 SPEC_ATOM_EH      reqd   r8v     10,  4,   4,  6       -0.7 -0.2  0  0  0  0  0  0  0  0
 SPEC_ATOM_RSMH2   opt    r8v     10,  1,   2,  8       0.8  0.8  0  0  0  0  0  0  0  0
 SPEC_ATOM_EH2     reqd   r8v     10,  2,   2, --       -1.5 -1
 SPEC_ATOM_LMXA    opt    i4       1,  1,   1,  0       3
 SPEC_ATOM_P       opt    r8v      4,  1,   4,  0       2.8999999999999999  2.8500000000000001  3.1800000000000002  4.1200000000000001
 SPEC_ATOM_MMOM    opt    r8v      4,  1,   2,  2       0  2  0  0
 SPEC_ATOM_IDMOD   opt    i4v      4,  1,   2,  2       0  1  0  0
 SPEC_ATOM_C-HOLE  opt    chr      1,  0,   *
 SPEC_ATOM_EREF    opt    r8       1,  1,   1,  0       -74.9949000000000012
  --- Parameters for site ---
 ... Site  1
 SITE_ATOM         reqd   chr      1,  0,   1           C
 SITE_ATOM_POS     reqd   r8v      3,  1,   3, --       0  0  0
  --- Parameters for Brillouin zone integration ---
 BZ_NKABC          reqd   i4v      3,  1,   3, --       4  4  4
 BZ_BZJOB          opt    i4v      3,  1,   1,  0       0
 BZ_METAL          opt    i4       1,  1,   1,  0       2
 BZ_TETRA          opt    lg       1,  1,   1,  0       T
 BZ_W              opt    r8       1,  1,   1,  0       0.004
 BZ_ZBAK           opt    r8       1,  1,   1,  0       0
 BZ_SAVDOS         opt    i4       1,  1,   1,  0       0
 BZ_NPTS           opt    i4       1,  1,   1,  0       200
 BZ_EFMAX          opt    r8       1,  1,   1,  0       5
 BZ_NEVMX          opt    i4       1,  1,   1,  0       5
  --- Parameters for Ewald sums ---
  --- Parameters for iterations ---
 ITER_NIT          opt    i4       1,  1,   1,  0       10
 ITER_MIX          opt    chr      1,  0,   1           A2,b=.5
 ITER_CONV         opt    r8       1,  1,   1,  0       0.00001
 ITER_CONVC        opt    r8       1,  1,   1,  0       0.0005
  --- Parameters for dynamics and statics ---
 DYN_MSTAT         opt    ---                           missing
  bndfp (warning): no sigm file found ... LDA calculation only


mmm === MTO setting ===
mmm ispec lmxb lpzex nkapii nkaphh=    1    3    0    2    2
mmm rsmh1    1  1.30  1.10
mmm   eh1    1 -0.70 -0.20
mmm rsmh2    1  0.80  0.80
mmm  eh2     1 -1.50 -1.00
mmm lh       1  1
 >>      0.03   exit  m_lmfinit       0.03
 >> level: 1  CPUsec=      0.03  enter m_lattic_init

                Plat                                  Qlat
   1.000000   1.000000   0.000000        0.500000   0.500000  -0.500000
   1.000000   0.000000   1.000000        0.500000  -0.500000   0.500000
   0.000000   1.000000   1.000000       -0.500000   0.500000   0.500000
  Cell vol=   1000.000000

 LATTC:  as= 2.000   tol= 1.00E-08   alat= 7.93701   awald= 0.200
         r1=  3.459   nkd=  79      q1=  2.571   nkg= 137
 >>      0.03   exit  m_lattic_ini    0.00
 >> level: 1  CPUsec=      0.03  enter m_mksym_init
SGROUP:  1 symmetry operations from 0
 SYMLAT: Bravais system is cubic       with 48 symmetry operations.
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-4
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
  21  r2(-0,1,1)                         
  22  m(-0,1,1)                          
  23  r2(1,0,-1)                         
  24  m(1,0,-1)                          
  25  r2y                                
  26  my                                 
  27  r4y                                
  28  i*r4y                              
  29  r4(0,-1,0)                         
  30  i*r4(0,-1,0)                       
  31  r2(1,-1,-0)                        
  32  m(1,-1,-0)                         
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
  45  r2(-0,1,-1)                        
  46  m(-0,1,-1)                         
  47  r2(1,1,-0)                         
  48  m(1,1,-0)                          
GROUPG: the following are sufficient to generate the space group:
 Generators:trans(cart)  = i*r3(-1,1,1) r4z
 Generators::trans(frac) = i*r3(-1,1,1) r4z
MKSYM: found  48  space group operations
SPLCLS: ibas iclass ispec label(ispec)
 SPLCLS     1    1    1     C
 >>      0.04   exit  m_mksym_init    0.01
 >> level: 1  CPUsec=      0.04  enter m_mkqp_init
 BZMESH:      8 irreducible QP from    4   4   4 shift=FFF
 TETIRR: sorting      384 tetrahedra ...
 >>      0.04   exit  m_mkqp_init     0.00
 >> level: 1  CPUsec=      0.04  enter m_supot_init
gvlst2: gmax=  13.994 a.u. created 45911 vectors of 125000 ( 36 %)
 SGVSYM: 1207 symmetry stars found for 45911 reciprocal lattice vectors
 >>      0.14   exit  m_supot_init    0.10

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol= 1.00E-06
 spec      l    rsm    eh     gmax    last term    cutoff
  C        0    1.30  -0.70   5.718    1.05E-06    3143 
  C        1    1.10  -0.20   7.226    1.06E-06    6375 
  C        0    0.80  -1.50   9.292    1.08E-06   13539 
  C        1    0.80  -1.00  10.038    1.00E-06   16961 
 >> level: 1  CPUsec=      0.14  enter m_suham_init
 >>      0.14   exit  m_suham_init    0.00
 >> level: 1  CPUsec=      0.14  enter m_qplist_init
 m_qplistinit:start
 >>      0.14   exit  m_qplist_ini    0.00
 >> level: 1  CPUsec=      0.14  enter m_qplist_qpsdivider
 >>      0.14   exit  m_qplist_qps    0.00
 >> level: 1  CPUsec=      0.14  enter m_igv2xall_init
 >>      0.14   exit  m_igv2xall_i    0.00
 >> level: 1  CPUsec=      0.14  enter m_hamindex_init

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv foca   rfoca

 species data:  augmentation                           density
 C        3.000  1.200    3    3         3  0.750  1.500    0   1.200
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv foca   rfoca
 C        3.000  1.200    3    3         3  0.750  1.500    0   1.200
 >>      0.14   exit  m_hamindex_i    0.00

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv foca   rfoca
 C        3.000  1.200    3    3         3  0.750  1.500    0   1.200
 >> level: 1  CPUsec=      0.14  enter lmfp

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv foca   rfoca
 C        3.000  1.200    3    3         3  0.750  1.500    0   1.200

 iors  : read rst restart file (binary mesh density)
 iors  : empty file ... nothing read

 rdovfa: read and overlap free-atom densities (mesh density) ...
 rdovfa: expected C,       read C        with rmt=  3.0000  mesh   369  0.020

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
 smves:: avg es pot at rmt= 0.006607  avg sphere pot= 0.000000  vconst=-0.006607

 smooth rhoves      2.063910   charge     2.661718
 smvxcm: all smrho_w is positive
smooth isp rhoeps rhomu vxcavg= 1 -1.109320 -1.522780 -0.218128
smooth isp rhoeps rhomu vxcavg= 2 -0.385581 -0.423770 -0.164568

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.90000   2.90000   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.18000   3.18000   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.12000   4.12000   0.00000   0.00000
 potential shift to crystal energy zero:    0.000003

 Energy terms:             smooth           local           total
   rhoval*vef             -3.931613       -10.430885       -14.362499
   rhoval*ves             -5.058553        -5.181775       -10.240327
   psnuc*ves               9.186373      -278.836573      -269.650199
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
 ... Done MPI k-loop: elapsed time=   0.5836

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.430662;   4.000000 electrons
         Sum occ. bands:   -2.743234, incl. Bloechl correction: -0.000179
         Mag. moment:       2.000000
Generating TDOS: efermi, and dos window=   -0.4307  -1.0507   2.5093
  mmmmm m_bandcal_2nd
       contr. to mm extrapolated for r>rmt:   0.163680 est. true mm = 1.960270
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85523  sum tc=    31.38701  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78555  sum tc=    31.54499  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.686857    3.838105   -0.151248      1.796590    2.141140   -0.344551
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 Harris energy:
 sumev=       -2.743234  val*vef=     -14.362499   sumtv=      11.619265
 sumec=      -39.640777  cor*vef=    -102.572782   ttcor=      62.932005
 rhoeps=      -9.601995     utot=    -139.945263    ehar=     -74.995989

 srhov:     -7.491327     -6.855971    -14.347299 sumev=   -2.743234   sumtv=   11.604065
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.008341  avg sphere pot=-0.006607  vconst=-0.008341

 smooth rhoves      3.531040   charge     4.151248
 smvxcm: all smrho_w is positive
smooth isp rhoeps rhomu vxcavg= 1 -2.415659 -3.336154 -0.237894
smooth isp rhoeps rhomu vxcavg= 2 -0.665883 -0.689339 -0.174538

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91338   2.90760   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16628   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10633   4.10242   0.00000   0.00000

 Energy terms:             smooth           local           total
   rhoval*vef             -9.130203        -5.227578       -14.357782
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

 Kohn-Sham energy:
 sumtv=       11.604065  sumtc=        62.932006   ekin=       74.536071
 rhoep=       -9.592117   utot=      -139.940099   ehks=      -74.996145
 mag. mom=     2.000000
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
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5
 mixrealsmooth= T
 wgtsmooth=   2.8284271247461901E-003
 mixrho:  sought 2 iter from file mixm; read 0.  RMS DQ=8.23e-3
 AMIX: nmix=0 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 4.12D-03
 mixrho: add corrections to qcell smrho =  0.78410D-07  0.39205D-10

 iors  : write rst restart file (binary mesh density)

   it  1  of 10    ehf=      -0.001089   ehk=      -0.001245
h zbak= 0 mmom= 2.000000 ehf(eV)=-0.014815 ehk(eV)=-0.016946 sev(eV)=-37.323889

 --- BNDFP:  begin iteration 2 of 10
 m_mkpot_init: Making one-particle potential ...
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.007474  avg sphere pot=-0.008341  vconst=-0.007474

 smooth rhoves      2.744964   charge     3.406483
 smvxcm: all smrho_w is positive
smooth isp rhoeps rhomu vxcavg= 1 -1.724195 -2.376346 -0.229068
smooth isp rhoeps rhomu vxcavg= 2 -0.520836 -0.552450 -0.169961

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91338   2.90760   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16628   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10633   4.10242   0.00000   0.00000
 potential shift to crystal energy zero:    0.000004

 Energy terms:             smooth           local           total
   rhoval*vef             -6.342509        -8.017597       -14.360106
   rhoval*ves             -4.965008        -5.279444       -10.244452
   psnuc*ves              10.454937      -280.095881      -269.640945
   utot                    2.744964      -142.687662      -139.942698
   rho*exc                -2.245031        -7.351973        -9.597004
   rho*vxc                -2.928796        -9.698090       -12.626886
   valence chg             3.406483         0.593517         4.000000
   valence mag             1.838696         0.161304         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.5012

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.431787;   4.000000 electrons
         Sum occ. bands:   -2.748502, incl. Bloechl correction: -0.000182
         Mag. moment:       2.000000
Generating TDOS: efermi, and dos window=   -0.4318  -1.0520   2.5081
       contr. to mm extrapolated for r>rmt:   0.163407 est. true mm = 1.960158
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85779  sum tc=    31.38712  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78816  sum tc=    31.54494  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.687503    3.843381   -0.155877      1.796752    2.129366   -0.332614
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 Harris energy:
 sumev=       -2.748502  val*vef=     -14.360106   sumtv=      11.611604
 sumec=      -39.645952  cor*vef=    -102.577957   ttcor=      62.932005
 rhoeps=      -9.597004     utot=    -139.942698    ehar=     -74.996093

 srhov:     -8.378297     -5.986619    -14.364916 sumev=   -2.748502   sumtv=   11.616414
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.008327  avg sphere pot=-0.007474  vconst=-0.008327

 smooth rhoves      3.521831   charge     4.155877
 smvxcm: all smrho_w is positive
smooth isp rhoeps rhomu vxcavg= 1 -2.413446 -3.333529 -0.237811
smooth isp rhoeps rhomu vxcavg= 2 -0.674101 -0.699705 -0.174472

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91349   2.90773   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16605   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10625   4.10242   0.00000   0.00000

 Energy terms:             smooth           local           total
   rhoval*vef             -9.152986        -5.213018       -14.366005
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

 Kohn-Sham energy:
 sumtv=       11.616414  sumtc=        62.932059   ekin=       74.548472
 rhoep=       -9.593958   utot=      -139.950644   ehks=      -74.996130
 mag. mom=     2.000000
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
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5
 wgtsmooth=   2.8284271247461901E-003
 mixrho:  sought 2 iter from file mixm; read 1.  RMS DQ=4.12e-3  last it=8.23e-3
 AMIX: nmix=1 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 2.06D-03
   tj:-1.00050
 mixrho: add corrections to qcell smrho =  0.54206D-08  0.27103D-11

 iors  : write rst restart file (binary mesh density)

   it  2  of 10    ehf=      -0.001193   ehk=      -0.001230
 From last iter    ehf=      -0.001089   ehk=      -0.001245
 diffe(q)= -0.000104 (0.004120)    tol= 0.000010 (0.000500)   more=T
i zbak= 0 mmom= 2.000000 ehf(eV)=-0.016231 ehk(eV)=-0.016733 sev(eV)=-37.395569

 --- BNDFP:  begin iteration 3 of 10
 m_mkpot_init: Making one-particle potential ...
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.008328  avg sphere pot=-0.008327  vconst=-0.008328

 smooth rhoves      3.522038   charge     4.156065
 smvxcm: all smrho_w is positive
smooth isp rhoeps rhomu vxcavg= 1 -2.413627 -3.333780 -0.237813
smooth isp rhoeps rhomu vxcavg= 2 -0.674141 -0.699744 -0.174474

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91349   2.90773   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16605   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10625   4.10242   0.00000   0.00000
 potential shift to crystal energy zero:    0.000005

 Energy terms:             smooth           local           total
   rhoval*vef             -9.153737        -5.212270       -14.366007
   rhoval*ves             -4.663429        -5.591094       -10.254523
   psnuc*ves              11.707506      -281.354225      -269.646719
   utot                    3.522038      -143.472659      -139.950621
   rho*exc                -3.087767        -6.506189        -9.593956
   rho*vxc                -4.033524        -8.589389       -12.622912
   valence chg             4.156065        -0.156065         4.000000
   valence mag             2.332738        -0.332738         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:      0.00000
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.4711

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.431877;   4.000000 electrons
         Sum occ. bands:   -2.749580, incl. Bloechl correction: -0.000186
         Mag. moment:       2.000000
Generating TDOS: efermi, and dos window=   -0.4319  -1.0522   2.5080
       contr. to mm extrapolated for r>rmt:   0.163406 est. true mm = 1.959998
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85748  sum tc=    31.38660  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78781  sum tc=    31.54437  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.687633    3.847282   -0.159649      1.796592    2.119618   -0.323026
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 Harris energy:
 sumev=       -2.749580  val*vef=     -14.366007   sumtv=      11.616427
 sumec=      -39.645287  cor*vef=    -102.577319   ttcor=      62.932032
 rhoeps=      -9.593956     utot=    -139.950621    ehar=     -74.996118

 srhov:     -9.174979     -5.194807    -14.369786 sumev=   -2.749580   sumtv=   11.620206
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.008290  avg sphere pot=-0.008328  vconst=-0.008290

 smooth rhoves      3.515406   charge     4.159649
 smvxcm: all smrho_w is positive
smooth isp rhoeps rhomu vxcavg= 1 -2.411330 -3.330903 -0.237828
smooth isp rhoeps rhomu vxcavg= 2 -0.680567 -0.707929 -0.174472

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91356   2.90781   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16589   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10619   4.10242   0.00000   0.00000

 Energy terms:             smooth           local           total
   rhoval*vef             -9.171543        -5.196916       -14.368459
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

 Kohn-Sham energy:
 sumtv=       11.620206  sumtc=        62.930978   ekin=       74.551184
 rhoep=       -9.594448   utot=      -139.952853   ehks=      -74.996118
 mag. mom=     2.000000
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
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5
 wgtsmooth=   2.8284271247461901E-003
 mixrho:  sought 2 iter from file mixm; read 2.  RMS DQ=6.10e-5  last it=4.12e-3
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 3.05D-05
   tj:-1.32030   0.66035
 mixrho: add corrections to qcell smrho =  0.65515D-07  0.32758D-10

 iors  : write rst restart file (binary mesh density)

   it  3  of 10    ehf=      -0.001218   ehk=      -0.001218
 From last iter    ehf=      -0.001193   ehk=      -0.001230
 diffe(q)= -0.000025 (0.000061)    tol= 0.000010 (0.000500)   more=T
i zbak= 0 mmom= 2.000000 ehf(eV)=-0.016570 ehk(eV)=-0.016569 sev(eV)=-37.410234

 --- BNDFP:  begin iteration 4 of 10
 m_mkpot_init: Making one-particle potential ...
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.008296  avg sphere pot=-0.008290  vconst=-0.008296

 smooth rhoves      3.516503   charge     4.159013
 smvxcm: all smrho_w is positive
smooth isp rhoeps rhomu vxcavg= 1 -2.411693 -3.331355 -0.237825
smooth isp rhoeps rhomu vxcavg= 2 -0.679468 -0.706529 -0.174472

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91356   2.90781   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16589   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10619   4.10242   0.00000   0.00000
 potential shift to crystal energy zero:    0.000005

 Energy terms:             smooth           local           total
   rhoval*vef             -9.168407        -5.199632       -14.368039
   rhoval*ves             -4.665470        -5.590532       -10.256002
   psnuc*ves              11.698476      -281.348080      -269.649604
   utot                    3.516503      -143.469306      -139.952803
   rho*exc                -3.091161        -6.503222        -9.594383
   rho*vxc                -4.037885        -8.585591       -12.623475
   valence chg             4.159013        -0.159013         4.000000
   valence mag             2.324659        -0.324659         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.4519

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.431601;   4.000000 electrons
         Sum occ. bands:   -2.748478, incl. Bloechl correction: -0.000186
         Mag. moment:       2.000000
Generating TDOS: efermi, and dos window=   -0.4316  -1.0519   2.5083
       contr. to mm extrapolated for r>rmt:   0.163546 est. true mm = 1.959960
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85678  sum tc=    31.38650  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78710  sum tc=    31.54427  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.687410    3.843876   -0.156467      1.796414    2.117499   -0.321085
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 Harris energy:
 sumev=       -2.748478  val*vef=     -14.368039   sumtv=      11.619561
 sumec=      -39.643884  cor*vef=    -102.575389   ttcor=      62.931505
 rhoeps=      -9.594383     utot=    -139.952803    ehar=     -74.996120

 srhov:     -9.160098     -5.206122    -14.366220 sumev=   -2.748478   sumtv=   11.617742
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.008266  avg sphere pot=-0.008296  vconst=-0.008266

 smooth rhoves      3.513827   charge     4.156466
 smvxcm: all smrho_w is positive
smooth isp rhoeps rhomu vxcavg= 1 -2.407917 -3.326145 -0.237866
smooth isp rhoeps rhomu vxcavg= 2 -0.679707 -0.707080 -0.174495

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91355   2.90780   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16591   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10620   4.10242   0.00000   0.00000

 Energy terms:             smooth           local           total
   rhoval*vef             -9.157930        -5.208957       -14.366887
   rhoval*ves             -4.667240        -5.588005       -10.255246
   psnuc*ves              11.694894      -281.340878      -269.645984
   utot                    3.513827      -143.464441      -139.950615
   rho*exc                -3.087623        -6.506395        -9.594018
   rho*vxc                -4.033225        -8.589768       -12.622993
   valence chg             4.156466        -0.156467         4.000000
   valence mag             2.321085        -0.321085         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       11.617742  sumtc=        62.930771   ekin=       74.548513
 rhoep=       -9.594018   utot=      -139.950615   ehks=      -74.996120
 mag. mom=     2.000000
 smvxcm: all smrho_w is positive
 smvxcm: all smrho_w is positive
 smvxcm: all smrho_w is positive

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00    -0.00   -0.00   -0.00
 shift forces to make zero average correction:           -0.00   -0.00   -0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5
 wgtsmooth=   2.8284271247461901E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=2.13e-5  last it=6.10e-5
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 1.06D-05
   tj:-0.14043   0.00544
 mixrho: add corrections to qcell smrho =  0.41473D-07  0.20737D-10

 iors  : write rst restart file (binary mesh density)

   it  4  of 10    ehf=      -0.001220   ehk=      -0.001220
 From last iter    ehf=      -0.001218   ehk=      -0.001218
 diffe(q)= -0.000002 (0.000021)    tol= 0.000010 (0.000500)   more=F
c zbak= 0 mmom= 2.000000 ehf(eV)=-0.016596 ehk(eV)=-0.016596 sev(eV)=-37.395238
 >>     10.26   exit  lmfp           10.12
CPU time:   10.259s     Sun Jun 19 18:38:09 2022   on process=0

  ==== procid=0 ====     calls      == cpu time ===   depth 1
  entry   xxxx  xxxx                per call  total  (depth is by TIM= in ctrl.*.)
      0      0      0        1      10.26      10.26   main
      0      0      0        1       0.03       0.03   |--m_lmfinit
      0      0      0        1       0.00       0.00   |--m_lattic_init
      0      0      0        1       0.01       0.01   |--m_mksym_init
      0      0      0        1       0.00       0.00   |--m_mkqp_init
      0      0      0        1       0.10       0.10   |--m_supot_init
      0      0      0        1       0.00       0.00   |--m_suham_init
      0      0      0        1       0.00       0.00   |--m_qplist_init
      0      0      0        1       0.00       0.00   |--m_qplist_qpsdivider
      0      0      0        1       0.00       0.00   |--m_igv2xall_init
      0      0      0        1       0.00       0.00   |--m_hamindex_init
      0      0      0        1      10.12      10.12   `--lmfp
Exit 0 procid= 0 OK! end of LMF ======================
