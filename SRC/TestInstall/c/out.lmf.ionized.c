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
   4        zbak                1.0000    
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
   4        zbak                1.0000    
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
 BZ_ZBAK           opt    r8       1,  1,   1,  0       1
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
 >>      0.05   exit  m_mksym_init    0.02
 >> level: 1  CPUsec=      0.05  enter m_mkqp_init
 BZMESH:      8 irreducible QP from    4   4   4 shift=FFF
 TETIRR: sorting      384 tetrahedra ...
 >>      0.05   exit  m_mkqp_init     0.00
 >> level: 1  CPUsec=      0.05  enter m_supot_init
gvlst2: gmax=  13.994 a.u. created 45911 vectors of 125000 ( 36 %)
 SGVSYM: 1207 symmetry stars found for 45911 reciprocal lattice vectors
 >>      0.17   exit  m_supot_init    0.12

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol= 1.00E-06
 spec      l    rsm    eh     gmax    last term    cutoff
  C        0    1.30  -0.70   5.718    1.05E-06    3143 
  C        1    1.10  -0.20   7.226    1.06E-06    6375 
  C        0    0.80  -1.50   9.292    1.08E-06   13539 
  C        1    0.80  -1.00  10.038    1.00E-06   16961 
 >> level: 1  CPUsec=      0.17  enter m_suham_init
 >>      0.17   exit  m_suham_init    0.00
 >> level: 1  CPUsec=      0.17  enter m_qplist_init
 m_qplistinit:start
 >>      0.17   exit  m_qplist_ini    0.00
 >> level: 1  CPUsec=      0.17  enter m_qplist_qpsdivider
 >>      0.17   exit  m_qplist_qps    0.00
 >> level: 1  CPUsec=      0.17  enter m_igv2xall_init
 >>      0.17   exit  m_igv2xall_i    0.00
 >> level: 1  CPUsec=      0.17  enter m_hamindex_init
 >>      0.17   exit  m_hamindex_i    0.00

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv foca   rfoca

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv foca   rfoca
 C        3.000  1.200    3    3         3  0.750  1.500    0   1.200
 C        3.000  1.200    3    3         3  0.750  1.500    0   1.200
 >> level: 1  CPUsec=      0.17  enter lmfp

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv foca   rfoca

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv foca   rfoca
 C        3.000  1.200    3    3         3  0.750  1.500    0   1.200
 C        3.000  1.200    3    3         3  0.750  1.500    0   1.200

 iors  : read rst restart file (binary mesh density)
 iors  : empty file ... nothing read

 rdovfa: read and overlap free-atom densities (mesh density) ...
 rdovfa: expected C,       read C        with rmt=  3.0000  mesh   369  0.020

 Free atom and overlapped crystal site charges:
   ib    true(FA)    smooth(FA)  true(OV)    smooth(OV)    local
    1    3.701869    2.363587    3.701843    2.363561    1.338282
 amom    1.810990    1.143831    1.810990    1.143831    0.667159
 Uniform density added to neutralize background, q=1.000000

 Smooth charge on mesh:            1.661718    moment    1.332841
 Sum of local charges:             1.338282    moments   0.667159
 Total valence charge:             3.000000    moment    2.000000
 Sum of core charges:              2.000000    moment    0.000000
 Sum of nuclear charges:          -6.000000
 Homogeneous background:           1.000000
 Deviation from neutrality:       -0.000000
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462

 Basis, after reading restart file
 site spec        pos (Cartesian coordinates)         pos (multiples of plat)
   1         C  0.000000   0.000000   0.000000    0.000000   0.000000   0.000000

 --- BNDFP:  begin iteration 1 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.006607  avg sphere pot= 0.000000  vconst=-0.006607
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      2.063910   charge     2.661718
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      201254  -4.9989572695401780E-004
 smvxcm: enforce positive smrho_w. Add srshift=   4.9989572696401783E-004
smooth isp rhoeps rhomu vxcavg= 1 -1.109298 -1.522749 -0.218032
smooth isp rhoeps rhomu vxcavg= 2 -0.385563 -0.423750 -0.164244

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
 potential shift to crystal energy zero:    0.000003

 Energy terms:             smooth           local           total
   rhoval*vef             -3.733857       -10.403599       -14.137456
   rhoval*ves             -5.058553        -5.181775       -10.240327
   psnuc*ves               9.186373      -278.836573      -269.650199
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
 ... Done MPI k-loop: elapsed time=   0.6024

 BZWTS : --- Tetrahedron Integration ---
 ... only filled or empty bands encountered: ev= -0.824805 ec= -0.772717
 VBmax= -0.824805 CBmin= -0.772717 gap =  0.052087 Ry =   0.708689 eV
 BZINTS: Fermi energy:     -0.824805;   3.000000 electrons
         Sum occ. bands:   -3.277726, incl. Bloechl correction:  0.000000
         Mag. moment:      -1.000000
Generating TDOS: efermi, and dos window=   -0.8248  -1.4280   2.1151
  mmmmm m_bandcal_2nd
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85491  sum tc=    31.38764  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78517  sum tc=    31.54593  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.130156  386.510559 -384.380403     -0.144414 -327.846612  327.702198
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 Harris energy:
 sumev=       -3.277726  val*vef=     -14.137456   sumtv=      10.859730
 sumec=      -39.640082  cor*vef=    -102.572087   ttcor=      62.932005
 rhoeps=      -9.595549     utot=    -139.945263    ehar=     -75.749078
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462

 srhov:  -1237.408097   1225.635209    -11.772888 sumev=   -3.277726   sumtv=    8.495162
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt=-0.205342  avg sphere pot=-0.006607  vconst= 0.205342
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves     38.909181   charge   388.380403
 smvxcm: all smrho_w is positive
smooth isp rhoeps rhomu vxcavg= 1 -179.107614 -116.221791 -0.227596
smooth isp rhoeps rhomu vxcavg= 2 -2291.203065 -3162.637449 -0.404877

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.89128   2.74142   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000

 Energy terms:             smooth           local           total
   rhoval*vef          -2698.892702      2687.560651       -11.332050
   rhoval*ves             82.992193       -91.768733        -8.776540
   psnuc*ves              -5.173830      -258.804129      -263.977959
   utot                   38.909181      -175.286431      -136.377249
   rho*exc             -2470.310679      2462.083828        -8.226851
   rho*vxc             -3278.859241      3268.028501       -10.830739
   valence chg           387.380403      -384.380403         3.000000
   valence mag          -328.702198       327.702198        -1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=        8.495162  sumtc=        62.933572   ekin=       71.428734
 rhoep=       -8.226851   utot=      -136.377249   ehks=      -73.175366
 mag. mom=    -1.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      201254  -4.9989874523704528E-004
 smvxcm: enforce positive smrho_w. Add srshift=   4.9989874524704531E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      201254  -4.9989270867099021E-004
 smvxcm: enforce positive smrho_w. Add srshift=   4.9989270868099024E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      201254  -4.9989572695401769E-004
 smvxcm: enforce positive smrho_w. Add srshift=   4.9989572696401772E-004

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
 mixrho:  sought 2 iter from file mixm; read 0.  RMS DQ=3.66e0
 mixrho: (warning) scr. and lin-mixed densities had 0 and 72123 negative points
 AMIX: nmix=0 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 1.83D+00
 mixrho: warning. negative smrho; isp number min=       1   93323 -0.24846D-03
 mixrho: warning. negative smrho; isp number min=       2   60589 -0.24856D-03
 mixrho: warning. negative smrho; isp number min=       1   93323 -0.24846D-03
 mixrho: warning. negative smrho; isp number min=       2   60589 -0.24856D-03
 mixrho: add corrections to qcell smrho =  0.37000D-07  0.18500D-10
 mixrho: warning. negative smrho; isp number min=       1   93323 -0.24846D-03
 mixrho: warning. negative smrho; isp number min=       2   60589 -0.24856D-03
 mixrho: warning. negative smrho; isp number min=       1   93323 -0.24846D-03
 mixrho: warning. negative smrho; isp number min=       2   60589 -0.24856D-03

 iors  : write rst restart file (binary mesh density)

   it  1  of 10    ehf=      -0.754178   ehk=       1.819534
h zbak= 1 mmom=-1.000000 ehf(eV)=-10.261199 ehk(eV)= 24.756220 sev(eV)=-44.596083

 --- BNDFP:  begin iteration 2 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt=-0.099368  avg sphere pot= 0.205342  vconst= 0.099368
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      8.892063   charge   195.521060
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      153912  -2.4856368766157421E-004
 smvxcm: enforce positive smrho_w. Add srshift=   2.4856368767157418E-004
smooth isp rhoeps rhomu vxcavg= 1 -73.344158 -48.831155 -0.245284
smooth isp rhoeps rhomu vxcavg= 2 -913.875935 -1260.393218 -0.337806

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0        64
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0        64
 vxcnsp (warning): negative rho: min val =  -9.48E-05
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0        64
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0        64
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.89128   2.74142   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
 potential shift to crystal energy zero:    0.000151

 Energy terms:             smooth           local           total
   rhoval*vef          -1389.461778      1376.441898       -13.019880
   rhoval*ves             15.777854       -25.671329        -9.893475
   psnuc*ves               2.006272      -268.820350      -266.814078
   utot                    8.892063      -147.245840      -138.353777
   rho*exc              -987.220093       978.344554        -8.875539
   rho*vxc             -1309.224373      1297.544830       -11.679543
   valence chg           194.521060      -191.521060         3.000000
   valence mag          -163.684679       164.184679         0.500000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.4575

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.720401;   3.000000 electrons
         Sum occ. bands:   -3.283332, incl. Bloechl correction: -0.000018
         Mag. moment:       1.000000
Generating TDOS: efermi, and dos window=   -0.7204  -1.3493   2.2195
       contr. to mm extrapolated for r>rmt:  -0.052767 est. true mm =-0.991699
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.56207  sum tc=    31.48865  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.52663  sum tc=    31.57660  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.897679    6.697076   -3.799397     -0.938931   -0.697976   -0.240956
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 Harris energy:
 sumev=       -3.283332  val*vef=     -13.019880   sumtv=       9.736548
 sumec=      -41.088703  cor*vef=    -104.021339   ttcor=      62.932636
 rhoeps=      -8.875539     utot=    -138.353777    ehar=     -74.560131
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462

 srhov:    -32.819351     18.717500    -14.101851 sumev=   -3.283332   sumtv=   10.818519
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt=-0.114016  avg sphere pot= 0.099368  vconst= 0.114016
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.285341   charge     7.799397
 smvxcm: all smrho_w is positive
smooth isp rhoeps rhomu vxcavg= 1 -3.640136 -4.711115 -0.139094
smooth isp rhoeps rhomu vxcavg= 2 -4.068815 -5.379629 -0.180306

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92271   2.91873   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000

 Energy terms:             smooth           local           total
   rhoval*vef            -28.173130        14.288384       -13.884746
   rhoval*ves             -4.065862        -6.547566       -10.613427
   psnuc*ves              14.636543      -283.137298      -268.500755
   utot                    5.285341      -144.842432      -139.557091
   rho*exc                -7.708951        -1.168266        -8.877217
   rho*vxc               -10.090744        -1.594273       -11.685018
   valence chg             6.799397        -3.799397         3.000000
   valence mag            -0.759044        -0.240956        -1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       10.818519  sumtc=        63.065249   ekin=       73.883768
 rhoep=       -8.877217   utot=      -139.557091   ehks=      -74.550540
 mag. mom=    -1.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      153912  -2.4856519251819498E-004
 smvxcm: enforce positive smrho_w. Add srshift=   2.4856519252819495E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      153912  -2.4856218280495344E-004
 smvxcm: enforce positive smrho_w. Add srshift=   2.4856218281495341E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      153912  -2.4856368766157421E-004
 smvxcm: enforce positive smrho_w. Add srshift=   2.4856368767157418E-004

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
 mixrho:  sought 2 iter from file mixm; read 1.  RMS DQ=1.81e0  last it=3.66e0
 mixrho: (warning) scr. and lin-mixed densities had 0 and 70587 negative points
 AMIX: nmix=1 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 9.04D-01
   tj: 0.33092
 mixrho: warning. negative smrho; isp number min=       1   92435 -0.16534D-03
 mixrho: warning. negative smrho; isp number min=       1   92435 -0.16534D-03
 mixrho: warning. negative smrho; isp number min=       2   59461 -0.16440D-03
 mixrho: warning. negative smrho; isp number min=       2   59461 -0.16440D-03
 mixrho: add corrections to qcell smrho =  0.79684D-07  0.39842D-10
 mixrho: warning. negative smrho; isp number min=       1   92435 -0.16534D-03
 mixrho: warning. negative smrho; isp number min=       2   59461 -0.16440D-03
 mixrho: warning. negative smrho; isp number min=       1   92435 -0.16534D-03
 mixrho: warning. negative smrho; isp number min=       2   59461 -0.16440D-03

 iors  : write rst restart file (binary mesh density)

   it  2  of 10    ehf=       0.434769   ehk=       0.444360
 From last iter    ehf=      -0.754178   ehk=       1.819534
 diffe(q)=  1.188948 (1.808465)    tol= 0.000010 (0.000500)   more=T
i zbak= 1 mmom=-1.000000 ehf(eV)= 5.915385 ehk(eV)= 6.045877 sev(eV)=-44.672355

 --- BNDFP:  begin iteration 3 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt=-0.104268  avg sphere pot= 0.114016  vconst= 0.104268
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      4.192906   charge   132.720804
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      151896  -1.6533819138315456E-004
 smvxcm: enforce positive smrho_w. Add srshift=   1.6533819139315456E-004
smooth isp rhoeps rhomu vxcavg= 1 -46.537077 -32.008652 -0.223141
smooth isp rhoeps rhomu vxcavg= 2 -536.842468 -741.177648 -0.306197

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92271   2.91873   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
 potential shift to crystal energy zero:    0.000103

 Energy terms:             smooth           local           total
   rhoval*vef           -971.582737       958.305220       -13.277517
   rhoval*ves              2.154218       -12.296578       -10.142359
   psnuc*ves               6.231593      -273.651459      -267.419866
   utot                    4.192906      -142.974018      -138.781113
   rho*exc              -583.379544       574.536966        -8.842578
   rho*vxc              -773.186300       761.550328       -11.635972
   valence chg           131.720804      -128.720804         3.000000
   valence mag          -109.179666       109.177857        -0.001809
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.3660

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.646428;   3.000000 electrons
         Sum occ. bands:   -3.147176, incl. Bloechl correction: -0.000023
         Mag. moment:       1.000000
Generating TDOS: efermi, and dos window=   -0.6464  -1.2770   2.2935
       contr. to mm extrapolated for r>rmt:   0.034243 est. true mm = 0.993974
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.45100  sum tc=    31.49394  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.44268  sum tc=    31.51704  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.907928    6.731013   -3.823085      0.959730    2.095333   -1.135603
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 Harris energy:
 sumev=       -3.147176  val*vef=     -13.277517   sumtv=      10.130341
 sumec=      -40.893677  cor*vef=    -103.892600   ttcor=      62.998923
 rhoeps=      -8.842578     utot=    -138.781113    ehar=     -74.494426
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462

 srhov:    -31.710401     17.619409    -14.090992 sumev=   -3.147176   sumtv=   10.943816
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt=-0.112609  avg sphere pot= 0.104268  vconst= 0.112609
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.512019   charge     7.823085
 smvxcm: all smrho_w is positive
smooth isp rhoeps rhomu vxcavg= 1 -5.042399 -6.941968 -0.172180
smooth isp rhoeps rhomu vxcavg= 2 -2.752057 -3.265467 -0.139132

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92208   2.92034   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000

 Energy terms:             smooth           local           total
   rhoval*vef            -28.397153        14.515256       -13.881897
   rhoval*ves             -3.855589        -6.740100       -10.595689
   psnuc*ves              14.879627      -283.304586      -268.424959
   utot                    5.512019      -145.022343      -139.510324
   rho*exc                -7.794456        -1.091603        -8.886059
   rho*vxc               -10.207435        -1.489260       -11.696695
   valence chg             6.823085        -3.823085         3.000000
   valence mag             2.135603        -1.135603         1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       10.943816  sumtc=        63.010975   ekin=       73.954791
 rhoep=       -8.886059   utot=      -139.510324   ehks=      -74.441592
 mag. mom=     1.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      151896  -1.6533918975128657E-004
 smvxcm: enforce positive smrho_w. Add srshift=   1.6533918976128657E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      151896  -1.6533719301502255E-004
 smvxcm: enforce positive smrho_w. Add srshift=   1.6533719302502255E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      151896  -1.6533819138315456E-004
 smvxcm: enforce positive smrho_w. Add srshift=   1.6533819139315456E-004

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
 mixrho:  sought 2 iter from file mixm; read 2.  RMS DQ=1.21e0  last it=1.81e0
 mixrho: (warning) scr. and lin-mixed densities had 0 and 69471 negative points
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 6.03D-01
   tj:-0.07971   0.23828
 mixrho: warning. negative smrho; isp number min=       1   88535 -0.11834D-03
 mixrho: warning. negative smrho; isp number min=       1   88535 -0.11834D-03
 mixrho: warning. negative smrho; isp number min=       1   88535 -0.11834D-03
 mixrho: warning. negative smrho; isp number min=       2   59749 -0.11860D-03
 mixrho: warning. negative smrho; isp number min=       2   59749 -0.11860D-03
 mixrho: warning. negative smrho; isp number min=       2   59749 -0.11860D-03
 mixrho: add corrections to qcell smrho =  0.33920D-07  0.16960D-10
 mixrho: warning. negative smrho; isp number min=       1   88535 -0.11834D-03
 mixrho: warning. negative smrho; isp number min=       2   59749 -0.11860D-03

 iors  : write rst restart file (binary mesh density)

   it  3  of 10    ehf=       0.500474   ehk=       0.553308
 From last iter    ehf=       0.434769   ehk=       0.444360
 diffe(q)=  0.065704 (1.206786)    tol= 0.000010 (0.000500)   more=T
i zbak= 1 mmom= 1.000000 ehf(eV)= 6.809344 ehk(eV)= 7.528202 sev(eV)=-42.819841

 --- BNDFP:  begin iteration 4 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt=-0.106416  avg sphere pot= 0.112609  vconst= 0.106416
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      3.097188   charge    97.613934
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      148284  -1.1859897196601064E-004
 smvxcm: enforce positive smrho_w. Add srshift=   1.1859897197601064E-004
smooth isp rhoeps rhomu vxcavg= 1 -34.193503 -24.702576 -0.215237
smooth isp rhoeps rhomu vxcavg= 2 -347.843181 -481.320499 -0.278762

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92208   2.92034   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
 potential shift to crystal energy zero:    0.000076

 Energy terms:             smooth           local           total
   rhoval*vef           -712.103608       698.646877       -13.456732
   rhoval*ves             -2.502213        -7.768907       -10.271120
   psnuc*ves               8.696590      -276.407703      -267.711114
   utot                    3.097188      -142.088305      -138.991117
   rho*exc              -382.036683       373.179399        -8.857284
   rho*vxc              -506.023075       494.366314       -11.656761
   valence chg            96.613934       -93.613934         3.000000
   valence mag           -77.483565        78.042585         0.559020
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.4359

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.651371;   3.000000 electrons
         Sum occ. bands:   -3.100739, incl. Bloechl correction: -0.000022
         Mag. moment:       1.000000
Generating TDOS: efermi, and dos window=   -0.6514  -1.2825   2.2886
       contr. to mm extrapolated for r>rmt:   0.031971 est. true mm = 0.994234
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.40462  sum tc=    31.45830  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.37206  sum tc=    31.53274  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.907057    6.690766   -3.783710      0.962264    2.201510   -1.239247
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 Harris energy:
 sumev=       -3.100739  val*vef=     -13.456732   sumtv=      10.355992
 sumec=      -40.776687  cor*vef=    -103.781650   ttcor=      63.004963
 rhoeps=      -8.857284     utot=    -138.991117    ehar=     -74.487446
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462

 srhov:    -31.889693     17.918613    -13.971079 sumev=   -3.100739   sumtv=   10.870340
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt=-0.112434  avg sphere pot= 0.106416  vconst= 0.112434
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.573294   charge     7.783709
 smvxcm: all smrho_w is positive
smooth isp rhoeps rhomu vxcavg= 1 -5.078874 -7.008134 -0.171751
smooth isp rhoeps rhomu vxcavg= 2 -2.642666 -3.104209 -0.140251

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92280   2.91966   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000

 Energy terms:             smooth           local           total
   rhoval*vef            -28.183228        14.361395       -13.821833
   rhoval*ves             -3.810310        -6.733266       -10.543577
   psnuc*ves              14.956899      -283.264770      -268.307871
   utot                    5.573294      -144.999018      -139.425724
   rho*exc                -7.721540        -1.156410        -8.877950
   rho*vxc               -10.112343        -1.573673       -11.686016
   valence chg             6.783709        -3.783710         3.000000
   valence mag             2.239246        -1.239247         1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       10.870340  sumtc=        62.991040   ekin=       73.861380
 rhoep=       -8.877950   utot=      -139.425724   ehks=      -74.442294
 mag. mom=     1.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      148284  -1.1859968934601603E-004
 smvxcm: enforce positive smrho_w. Add srshift=   1.1859968935601603E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      148284  -1.1859825458600526E-004
 smvxcm: enforce positive smrho_w. Add srshift=   1.1859825459600527E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      148284  -1.1859897196601065E-004
 smvxcm: enforce positive smrho_w. Add srshift=   1.1859897197601066E-004

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
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=8.68e-1  last it=1.21e0
 mixrho: (warning) scr. and lin-mixed densities had 0 and 68313 negative points
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 4.34D-01
   tj:-1.89092  -0.24126
 mixrho: warning. negative smrho; isp number min=       2   75363 -0.27386D-05
 mixrho: warning. negative smrho; isp number min=       2   75363 -0.27386D-05
 mixrho: add corrections to qcell smrho =  0.21886D-07  0.10943D-10
 mixrho: warning. negative smrho; isp number min=       2   75363 -0.27386D-05
 mixrho: warning. negative smrho; isp number min=       2   75363 -0.27386D-05

 iors  : write rst restart file (binary mesh density)

   it  4  of 10    ehf=       0.507454   ehk=       0.552606
 From last iter    ehf=       0.500474   ehk=       0.553308
 diffe(q)=  0.006980 (0.867844)    tol= 0.000010 (0.000500)   more=T
i zbak= 1 mmom= 1.000000 ehf(eV)= 6.904319 ehk(eV)= 7.518645 sev(eV)=-42.188041

 --- BNDFP:  begin iteration 5 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt=-0.111951  avg sphere pot= 0.112434  vconst= 0.111951
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.681112   charge     7.656798
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       75363  -2.7386131252040168E-006
 smvxcm: enforce positive smrho_w. Add srshift=   2.7386131352040166E-006
smooth isp rhoeps rhomu vxcavg= 1 -5.533423 -7.691911 -0.182120
smooth isp rhoeps rhomu vxcavg= 2 -2.073163 -2.276637 -0.131597

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       864
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       864
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       864
 vxcnsp (warning): negative rho: min val =  -2.47E-02
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       864
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92280   2.91966   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
 potential shift to crystal energy zero:    0.000008

 Energy terms:             smooth           local           total
   rhoval*vef            -27.576023        13.721878       -13.854146
   rhoval*ves             -3.713848        -6.795380       -10.509228
   psnuc*ves              15.076072      -283.360792      -268.284720
   utot                    5.681112      -145.078086      -139.396974
   rho*exc                -7.606586        -1.318410        -8.924996
   rho*vxc                -9.968548        -1.781958       -11.750506
   valence chg             6.656798        -3.656798         3.000000
   valence mag             3.203555        -1.645420         1.558135
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.4259

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.641524;   3.000000 electrons
         Sum occ. bands:   -2.963194, incl. Bloechl correction: -0.000020
         Mag. moment:       1.000000
Generating TDOS: efermi, and dos window=   -0.6415  -1.2742   2.2984
       contr. to mm extrapolated for r>rmt:   0.030009 est. true mm = 0.994602
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.29856  sum tc=    31.39771  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.22464  sum tc=    31.56334  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905558    6.881739   -3.976180      0.964593    2.105896   -1.141302
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 Harris energy:
 sumev=       -2.963194  val*vef=     -13.854146   sumtv=      10.890951
 sumec=      -40.523196  cor*vef=    -103.521146   ttcor=      62.997950
 rhoeps=      -8.924996     utot=    -139.396974    ehar=     -74.433069
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462

 srhov:    -28.985226     15.325213    -13.660014 sumev=   -2.963194   sumtv=   10.696820
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt=-0.111738  avg sphere pot= 0.111951  vconst= 0.111738
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.815372   charge     7.976180
 smvxcm: all smrho_w is positive
smooth isp rhoeps rhomu vxcavg= 1 -5.157014 -7.099561 -0.171472
smooth isp rhoeps rhomu vxcavg= 2 -2.827199 -3.356417 -0.142057

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92454   2.92055   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000

 Energy terms:             smooth           local           total
   rhoval*vef            -29.355603        15.676303       -13.679300
   rhoval*ves             -3.623992        -6.796900       -10.420891
   psnuc*ves              15.254736      -283.317041      -268.062305
   utot                    5.815372      -145.056970      -139.241598
   rho*exc                -7.984212        -0.874749        -8.858962
   rho*vxc               -10.455978        -1.205003       -11.660981
   valence chg             6.976180        -3.976180         3.000000
   valence mag             2.141302        -1.141302         1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       10.696820  sumtc=        62.961046   ekin=       73.657866
 rhoep=       -8.858962   utot=      -139.241598   ehks=      -74.442694
 mag. mom=     1.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       75363  -2.7385804235145875E-006
 smvxcm: enforce positive smrho_w. Add srshift=   2.7385804335145873E-006
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       75363  -2.7386458268934460E-006
 smvxcm: enforce positive smrho_w. Add srshift=   2.7386458368934459E-006
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       75363  -2.7386131252040168E-006
 smvxcm: enforce positive smrho_w. Add srshift=   2.7386131352040166E-006

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
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=4.06e-3  last it=8.68e-1
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 2.03D-03
   tj:-4.14964   2.99042
 mixrho: warning. negative smrho; isp number min=       2   14753 -0.45008D-06
 mixrho: warning. negative smrho; isp number min=       2   14753 -0.45008D-06
 mixrho: add corrections to qcell smrho =  0.43065D-07  0.21533D-10
 mixrho: warning. negative smrho; isp number min=       2   14753 -0.45008D-06
 mixrho: warning. negative smrho; isp number min=       2   14753 -0.45008D-06

 iors  : write rst restart file (binary mesh density)

   it  5  of 10    ehf=       0.561831   ehk=       0.552206
 From last iter    ehf=       0.507454   ehk=       0.552606
 diffe(q)=  0.054377 (0.004062)    tol= 0.000010 (0.000500)   more=T
i zbak= 1 mmom= 1.000000 ehf(eV)= 7.644162 ehk(eV)= 7.513206 sev(eV)=-40.316627

 --- BNDFP:  begin iteration 6 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt=-0.111698  avg sphere pot= 0.111738  vconst= 0.111698
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.794618   charge     8.339198
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       14753  -4.5008240755622230E-007
 smvxcm: enforce positive smrho_w. Add srshift=   4.5008241755622229E-007
smooth isp rhoeps rhomu vxcavg= 1 -5.260970 -7.178975 -0.172136
smooth isp rhoeps rhomu vxcavg= 2 -3.348890 -4.096859 -0.143024

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       160
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       160
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       160
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       160
 vxcnsp (warning): negative rho: min val =  -8.46E-04
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92454   2.92055   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
 potential shift to crystal energy zero:    0.000008

 Energy terms:             smooth           local           total
   rhoval*vef            -31.558834        17.872375       -13.686459
   rhoval*ves             -3.645169        -6.780394       -10.425563
   psnuc*ves              15.234405      -283.338712      -268.104307
   utot                    5.794618      -145.059553      -139.264935
   rho*exc                -8.609860        -0.252097        -8.861956
   rho*vxc               -11.275834        -0.389165       -11.664999
   valence chg             7.339198        -4.339198         3.000000
   valence mag             1.835257        -0.815652         1.019606
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.3700

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.623568;   3.000000 electrons
         Sum occ. bands:   -2.982078, incl. Bloechl correction: -0.000021
         Mag. moment:       1.000000
Generating TDOS: efermi, and dos window=   -0.6236  -1.2563   2.3164
       contr. to mm extrapolated for r>rmt:   0.032902 est. true mm = 0.994404
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.31469  sum tc=    31.42826  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.26617  sum tc=    31.53657  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905891    6.938226   -4.032335      0.961502    1.954185   -0.992682
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 Harris energy:
 sumev=       -2.982078  val*vef=     -13.686459   sumtv=      10.704381
 sumec=      -40.580863  cor*vef=    -103.560372   ttcor=      62.979509
 rhoeps=      -8.861956     utot=    -139.264935    ehar=     -74.443002
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462

 srhov:    -29.841399     16.125214    -13.716184 sumev=   -2.982078   sumtv=   10.734106
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt=-0.111800  avg sphere pot= 0.111698  vconst= 0.111800
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.797310   charge     8.032335
 smvxcm: all smrho_w is positive
smooth isp rhoeps rhomu vxcavg= 1 -5.094275 -6.988417 -0.172086
smooth isp rhoeps rhomu vxcavg= 2 -2.976857 -3.580724 -0.141092

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92430   2.92159   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000

 Energy terms:             smooth           local           total
   rhoval*vef            -29.674689        15.963306       -13.711383
   rhoval*ves             -3.637995        -6.811908       -10.449903
   psnuc*ves              15.232616      -283.342527      -268.109911
   utot                    5.797310      -145.077218      -139.279907
   rho*exc                -8.071132        -0.790793        -8.861925
   rho*vxc               -10.569142        -1.095715       -11.664857
   valence chg             7.032335        -4.032335         3.000000
   valence mag             1.992682        -0.992682         1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       10.734106  sumtc=        62.964833   ekin=       73.698939
 rhoep=       -8.861925   utot=      -139.279907   ehks=      -74.442893
 mag. mom=     1.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       14753  -4.5008193281310037E-007
 smvxcm: enforce positive smrho_w. Add srshift=   4.5008194281310035E-007
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       14753  -4.5008288229934429E-007
 smvxcm: enforce positive smrho_w. Add srshift=   4.5008289229934428E-007
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       14753  -4.5008240755622236E-007
 smvxcm: enforce positive smrho_w. Add srshift=   4.5008241755622234E-007

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
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=2.55e-3  last it=4.06e-3
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 1.27D-03
   tj:-0.12269  -0.00374
 mixrho: add corrections to qcell smrho =  0.43090D-07  0.21545D-10

 iors  : write rst restart file (binary mesh density)

   it  6  of 10    ehf=       0.551898   ehk=       0.552007
 From last iter    ehf=       0.561831   ehk=       0.552206
 diffe(q)= -0.009933 (0.002550)    tol= 0.000010 (0.000500)   more=T
i zbak= 1 mmom= 1.000000 ehf(eV)= 7.509019 ehk(eV)= 7.510493 sev(eV)=-40.573558

 --- BNDFP:  begin iteration 7 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt=-0.111746  avg sphere pot= 0.111800  vconst= 0.111746
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.814663   charge     8.064426
 smvxcm: all smrho_w is positive
smooth isp rhoeps rhomu vxcavg= 1 -5.113231 -7.012554 -0.171597
smooth isp rhoeps rhomu vxcavg= 2 -3.008614 -3.623046 -0.142374

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92430   2.92159   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
 potential shift to crystal energy zero:    0.000008

 Energy terms:             smooth           local           total
   rhoval*vef            -29.872962        16.181196       -13.691766
   rhoval*ves             -3.625296        -6.809166       -10.434462
   psnuc*ves              15.254621      -283.355837      -268.101216
   utot                    5.814663      -145.082502      -139.267839
   rho*exc                -8.121845        -0.737187        -8.859032
   rho*vxc               -10.635600        -1.025381       -11.660981
   valence chg             7.064426        -4.064426         3.000000
   valence mag             1.968928        -0.991299         0.977628
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.4412

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.620802;   3.000000 electrons
         Sum occ. bands:   -2.979959, incl. Bloechl correction: -0.000021
         Mag. moment:       1.000000
Generating TDOS: efermi, and dos window=   -0.6208  -1.2536   2.3191
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
       contr. to mm extrapolated for r>rmt:   0.033360 est. true mm = 0.994357
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.31295  sum tc=    31.43015  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.26691  sum tc=    31.53372  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905712    6.941540   -4.035829      0.960997    1.935879   -0.974882
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 Harris energy:
 sumev=       -2.979959  val*vef=     -13.691766   sumtv=      10.711807
 sumec=      -40.579867  cor*vef=    -103.552037   ttcor=      62.972170
 rhoeps=      -8.859032     utot=    -139.267839    ehar=     -74.442893

 srhov:    -29.732591     16.021014    -13.711577 sumev=   -2.979959   sumtv=   10.731618
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt=-0.111807  avg sphere pot= 0.111746  vconst= 0.111807
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.798710   charge     8.035829
 smvxcm: all smrho_w is positive
smooth isp rhoeps rhomu vxcavg= 1 -5.082816 -6.969600 -0.172272
smooth isp rhoeps rhomu vxcavg= 2 -2.991819 -3.604024 -0.141042

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92425   2.92170   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000

 Energy terms:             smooth           local           total
   rhoval*vef            -29.694440        15.984642       -13.709799
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

 Kohn-Sham energy:
 sumtv=       10.731618  sumtc=        62.963876   ekin=       73.695494
 rhoep=       -8.861326   utot=      -139.277075   ehks=      -74.442907
 mag. mom=     1.000000
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
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=2.35e-4  last it=2.55e-3
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 1.17D-04
   tj:-0.04684  -0.02066
 mixrho: add corrections to qcell smrho =  0.41372D-07  0.20686D-10

 iors  : write rst restart file (binary mesh density)

   it  7  of 10    ehf=       0.552007   ehk=       0.551993
 From last iter    ehf=       0.551898   ehk=       0.552007
 diffe(q)=  0.000108 (0.000235)    tol= 0.000010 (0.000500)   more=T
i zbak= 1 mmom= 1.000000 ehf(eV)= 7.510494 ehk(eV)= 7.510311 sev(eV)=-40.544725

 --- BNDFP:  begin iteration 8 of 10
 m_mkpot_init: Making one-particle potential ...
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt=-0.111776  avg sphere pot= 0.111807  vconst= 0.111776
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.808401   charge     8.048600
 smvxcm: all smrho_w is positive
smooth isp rhoeps rhomu vxcavg= 1 -5.089615 -6.978464 -0.171907
smooth isp rhoeps rhomu vxcavg= 2 -3.004265 -3.620356 -0.142041

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92425   2.92170   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000
 potential shift to crystal energy zero:    0.000008

 Energy terms:             smooth           local           total
   rhoval*vef            -29.772988        16.073319       -13.699669
   rhoval*ves             -3.630169        -6.811251       -10.441420
   psnuc*ves              15.246971      -283.349345      -268.102374
   utot                    5.808401      -145.080298      -139.271897
   rho*exc                -8.093880        -0.765710        -8.859590
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
 ... Done MPI k-loop: elapsed time=   1.3620

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.620580;   3.000000 electrons
         Sum occ. bands:   -2.979106, incl. Bloechl correction: -0.000021
         Mag. moment:       1.000000
Generating TDOS: efermi, and dos window=   -0.6206  -1.2533   2.3193
       contr. to mm extrapolated for r>rmt:   0.033442 est. true mm = 0.994346
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.31209  sum tc=    31.42985  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.26616  sum tc=    31.53366  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905645    6.940659   -4.035014      0.960903    1.934902   -0.973999
  Symmetrize density..
 Make new boundary conditions for phi,phidot..

 Harris energy:
 sumev=       -2.979106  val*vef=     -13.699669   sumtv=      10.720563
 sumec=      -40.578253  cor*vef=    -103.546277   ttcor=      62.968023
 rhoeps=      -8.859590     utot=    -139.271897    ehar=     -74.442901
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462

 srhov:    -29.708991     15.999863    -13.709128 sumev=   -2.979106   sumtv=   10.730022
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 Energy for background charge  q=  0.100000D+01 radius r= 6.2035049089939989 E=9/5*q*q/r= 0.29015855172296462
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt=-0.111809  avg sphere pot= 0.111776  vconst= 0.111809
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.799464   charge     8.035014
 smvxcm: all smrho_w is positive
smooth isp rhoeps rhomu vxcavg= 1 -5.080959 -6.966873 -0.172322
smooth isp rhoeps rhomu vxcavg= 2 -2.991790 -3.604267 -0.141046

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.92424   2.92171   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.14758   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10242   4.10242   0.00000   0.00000

 Energy terms:             smooth           local           total
   rhoval*vef            -29.689605        15.980997       -13.708609
   rhoval*ves             -3.637179        -6.810890       -10.448070
   psnuc*ves              15.236108      -283.338776      -268.102668
   utot                    5.799464      -145.074833      -139.275369
   rho*exc                -8.072749        -0.788320        -8.861069
   rho*vxc               -10.571140        -1.092582       -11.663723
   valence chg             7.035014        -4.035014         3.000000
   valence mag             1.973999        -0.973999         1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       10.730022  sumtc=        62.963506   ekin=       73.693528
 rhoep=       -8.861069   utot=      -139.275369   ehks=      -74.442910
 mag. mom=     1.000000
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
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.10e-4  last it=2.35e-4
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 5.51D-05
   tj:-0.72282  -0.00411
 mixrho: add corrections to qcell smrho =  0.39446D-07  0.19723D-10

 iors  : write rst restart file (binary mesh density)

   it  8  of 10    ehf=       0.551999   ehk=       0.551990
 From last iter    ehf=       0.552007   ehk=       0.551993
 diffe(q)= -0.000008 (0.000110)    tol= 0.000010 (0.000500)   more=F
c zbak= 1 mmom= 1.000000 ehf(eV)= 7.510389 ehk(eV)= 7.510266 sev(eV)=-40.533124
 >>     20.50   exit  lmfp           20.33
CPU time:   20.502s     Sun Jun 19 18:38:31 2022   on process=0

  ==== procid=0 ====     calls      == cpu time ===   depth 1
  entry   xxxx  xxxx                per call  total  (depth is by TIM= in ctrl.*.)
      0      0      0        1      20.50      20.50   main
      0      0      0        1       0.03       0.03   |--m_lmfinit
      0      0      0        1       0.00       0.00   |--m_lattic_init
      0      0      0        1       0.02       0.02   |--m_mksym_init
      0      0      0        1       0.00       0.00   |--m_mkqp_init
      0      0      0        1       0.12       0.12   |--m_supot_init
      0      0      0        1       0.00       0.00   |--m_suham_init
      0      0      0        1       0.00       0.00   |--m_qplist_init
      0      0      0        1       0.00       0.00   |--m_qplist_qpsdivider
      0      0      0        1       0.00       0.00   |--m_igv2xall_init
      0      0      0        1       0.00       0.00   |--m_hamindex_init
      0      0      0        1      20.33      20.33   `--lmfp
Exit 0 procid= 0 OK! end of LMF ======================
