===START LMFA   =====================================
 mpisize=           1

mmm === MTO setting ===
mmm ispec lmxb lpzex nkapii nkaphh=    1    3    0    1    1
mmm rsmh1    1  2.50  2.50  1.00  0.00
mmm   eh1    1 -0.01 -0.01 -0.01 -0.01
mmm rsmh2    1  0.00  0.00  0.00  0.00
mmm  eh2     1  0.00  0.00  0.00  0.00
mmm pz       1  5.50  5.50  4.50  0.00
mmm lh       1  2  2

                Plat                                  Qlat
   0.000000   0.500000   0.500000       -1.000000   1.000000   1.000000
   0.500000   0.000000   0.500000        1.000000  -1.000000   1.000000
   0.500000   0.500000   0.000000        1.000000   1.000000  -1.000000
  Cell vol= 78.538660

 LATTC:  as= 2.000   tol=1.00E-12   alat= 6.79800   awald= 0.467
         r1=  2.252   nkd= 201      q1=  6.871   nkg= 331

 SGROUP: 1 symmetry operations from 0 generators
 ADDBAS: basis is already complete --- no sites added
 SYMLAT: Bravais system is cubic with 48 symmetry operations.
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-5
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
  21  r2(1,1,-0)                         
  22  m(1,1,-0)                          
  23  r2(1,-0,-1)                        
  24  m(1,-0,-1)                         
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
  47  r2(-0,1,1)                         
  48  m(-0,1,1)                          
 GROUPG: the following are sufficient to generate the space group:
 Generator(cart): i*r3(1,1,-1) r4x
 Generator(frac): i*r3(1,1,-1) r4x
 MKSYM:  found 48 space group operations ... includes inversion
 goto freats

ttt: pnu qat=  1  0     4.650     1.000
ttt: pnu qat=  1  1     4.340     0.000
ttt: pnu qat=  1  2     3.870    10.000
ttt: pnu qat=  1  3     4.110     0.000

conf:------------------------------------------------------
conf:SPEC_ATOM= A : --- Table for atomic configuration ---
conf:  isp  l  int(P) int(P)z    Qval     Qcore   CoreConf
conf:    1  0       4  5         1.000    6.000 => 1,2,3,
conf:    1  1       4  5         0.000   12.000 => 2,3,
conf:    1  2       3  4        10.000    0.000 => 
conf:    1  3       4  0         0.000    0.000 => 
usedQ=     1.000     0.000    10.000     0.000

conf: Species A:  Z=29  Qc=18  R=2.311271  Q=0
conf:   rmt=2.311271  rmax=48.805862  a=0.025  nr=393  nr(rmax)=515
 goto atomc xxx
 atomsc nmcore=           0

  iter     qint         drho          vh0          rho0          vsum     beta
    1   29.000000   4.725E+03      145.0000    0.1442E+03      -58.2772   0.30
   51   29.000000   4.211E-05      274.8263    0.2631E+05     -130.7915   0.30

 end of atomsc xxxxx
 vsum=  -130.79144076069792                1

 sumev=-4.333254  etot=-3304.416258  eref=-3304.434500  diff= 0.018242

 Free-atom wavefunctions:
 valence:      eval       node at      max at       c.t.p.   rho(r>rmt)
   4s      -0.36411         0.890       2.256       3.582     0.643062
   5s      -0.00028         3.669      10.794      19.873     0.990448
   4p      -0.06295         0.975       3.484       7.414     0.901829
   5p       0.00796         6.760      30.414      48.806*    0.999240
   3d      -0.39691         0.000       0.600       3.429     0.056076
   4d       0.01308         1.868      33.290      48.806*    0.999995
   4f       0.01948         0.000      35.393      48.806*    1.000000

 core:        ecore       node at      max at       c.t.p.   rho(r>rmt)
   1s    -649.07634         0.000       0.034       0.069     0.000000
   2s     -77.91382         0.070       0.197       0.308     0.000000
   2p     -67.32532         0.000       0.158       0.335     0.000000
   3s      -8.39248         0.288       0.614       0.895     0.000141
   3p      -5.29682         0.260       0.619       1.078     0.000727

 Optimise free-atom basis for species A, rmt=2.311271
 l  it    Rsm      Eh     stiffR   stiffE      Eval      Exact     Pnu    Ql
 0  10   1.156  -0.108       0.0   6566.8   -0.34267  -0.36411    4.76   1.00
 1  11   1.156  -0.100       0.0  -3498.9    0.17087  -0.06295    4.56   0.00
 2   6   1.156  -0.509       0.0     20.5   -0.39483  -0.39691    3.89  10.00
 eigenvalue sum:  exact  -4.33325    opt basis  -4.29094    error 0.04231
 tailsm: init

 tailsm: fit tails to 6 smoothed hankels, rmt= 2.31127, rsm= 1.15564
  ---E:energies of smHankels. C:fitting coeeficient for core tail. ---
 E:    -1.00000    -2.00000    -4.00000    -6.00000    -9.00000    -15.0000
 C:    -0.07160    10.75053    -187.492    1222.023    -4717.79    21166.81
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
 tailsm: end
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
Sum of reference energies:                  -3304.434500000000
OK! end of LMFA ======================
===START LMF   =====================================
 mpisize=           4
  bndfp (warning): no sigm file found ... LDA calculation only

mmm === MTO setting ===
mmm ispec lmxb lpzex nkapii nkaphh=    1    3    0    1    1
mmm rsmh1    1  2.50  2.50  1.00  0.00
mmm   eh1    1 -0.01 -0.01 -0.01 -0.01
mmm rsmh2    1  0.00  0.00  0.00  0.00
mmm  eh2     1  0.00  0.00  0.00  0.00
mmm pz       1  5.50  5.50  4.50  0.00
mmm lh       1  2  2

                Plat                                  Qlat
   0.000000   0.500000   0.500000       -1.000000   1.000000   1.000000
   0.500000   0.000000   0.500000        1.000000  -1.000000   1.000000
   0.500000   0.500000   0.000000        1.000000   1.000000  -1.000000
  Cell vol= 78.538660

 LATTC:  as= 2.000   tol=1.00E-12   alat= 6.79800   awald= 0.467
         r1=  2.252   nkd= 201      q1=  6.871   nkg= 331

 SGROUP: 1 symmetry operations from 0 generators
 ADDBAS: basis is already complete --- no sites added
 SYMLAT: Bravais system is cubic with 48 symmetry operations.
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-5
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
  21  r2(1,1,-0)                         
  22  m(1,1,-0)                          
  23  r2(1,-0,-1)                        
  24  m(1,-0,-1)                         
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
  47  r2(-0,1,1)                         
  48  m(-0,1,1)                          
 GROUPG: the following are sufficient to generate the space group:
 Generator(cart): i*r3(1,1,-1) r4x
 Generator(frac): i*r3(1,1,-1) r4x
 MKSYM:  found 48 space group operations ... includes inversion
 
 BZMESH:     60 irreducible QP from    8   8   8 shift=TTT
 TETIRR: sorting 3072 tetrahedra ...
 264 inequivalent ones found
 >> level: 1  CPUsec=      0.06  enter lmfp
 gen_hamindex: not readin QGpsi.

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv  kmxv foca   rfoca
 A        2.311  0.925    3    4         3  0.578  1.156    15    1   0.925

 MSHSIZ: mesh has 11 x 11 x 11 divisions; length 0.437, 0.437, 0.437
         generated from gmax = 9.0 a.u. : 941 vectors of 1331 (70%)
 goto end of reading rst or atm           1
 goto end of reading rst or atm           1
 goto end of reading rst or atm           1

 GVLST2: gmax = 9.0 a.u. created 941 vectors of 1331 (70%)
         (input) mesh has 11 x 11 x 11 divisions; length 0.437, 0.437, 0.437
 SGVSYM: 41 symmetry stars found for 941 reciprocal lattice vectors

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol= 1.00E-06
 spec      l    rsm    eh     gmax    last term    cutoff
  A        0*   2.50  -0.01   2.974    2.30E-05      27 
  A        1    2.50  -0.01   3.093    1.29E-06      51 
  A        2    1.00  -0.01   8.508    1.16E-06     869 
 goto end of reading rst or atm           1

 iors  : read restart file (binary, mesh density) 
 iors  : empty file ... nothing read

 rdovfa: read and overlap free-atom densities (mesh density) ...
 rdovfa: expected A,       read A        with rmt=  2.3113  mesh   393  0.025

 ovlpfa: overlap smooth part of FA densities
 total smooth Q = 4.646654

 Free atom and overlapped crystal site charges:
   ib    true(FA)    smooth(FA)  true(OV)    smooth(OV)    local
 end of reading rst or atm
 end of reading rst or atm
 end of reading rst or atm
    1    9.796164    3.442818   10.275300    3.921954    6.353346

 Smooth charge on mesh:            4.646654
 Sum of local charges:             6.353346
 Total valence charge:            11.000000
 Sum of core charges:             18.000000
 Sum of nuclear charges:         -29.000000
 Homogeneous background:           0.000000
 Deviation from neutrality:       -0.000000
 end of reading rst or atm
 m_qplistinit:start
m_qplist_qspdivider: rank,(iqini,ispini),(iqend,ispend)=    1, (   16 1),  (   30 1)
m_qplist_qspdivider: rank,(iqini,ispini),(iqend,ispend)=    2, (   31 1),  (   45 1)
m_qplist_qspdivider: rank,(iqini,ispini),(iqend,ispend)=    3, (   46 1),  (   60 1)
m_qplist_qspdivider: rank,(iqini,ispini),(iqend,ispend)=    0, (    1 1),  (   15 1)

 Basis, after reading restart file
 site spec        pos (Cartesian coordinates)         pos (multiples of plat)
   1  A         0.000000   0.000000   0.000000    0.000000   0.000000   0.000000

 --- BNDFP:  begin iteration 1 of 12
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.554993  avg sphere pot= 0.633521  vconst=-0.554993

 smooth rhoves     11.022237   charge     4.646654
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -3.843801   rhomu =   -5.010455  avg vxc =   -0.851784 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000013

 Energy terms:             smooth           local           total
   rhoval*vef            -12.156987      -177.336819      -189.493805
   rhoval*ves            -46.689417      -115.324371      -162.013788
   psnuc*ves              68.733890    -12976.662453    -12907.928563
   utot                   11.022237     -6545.993412     -6534.971175
   rho*exc                -3.843801      -126.414296      -130.258096
   rho*vxc                -5.010455      -167.409313      -172.419769
   valence chg             4.646654         6.353346        11.000000
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
 ... Done MPI k-loop: elapsed time=   0.0404

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.144577;  11.000000 electrons
         Sum occ. bands:   -0.853241, incl. Bloechl correction: -0.006586
Generating TDOS: efermi, and dos window=    0.1446  -0.5000   1.6446
  mmmmm m_bandcal_2nd

 mkrout:  Qtrue      sm,loc       local
   1    9.927753    3.113494    6.814259

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  pnunew: ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.508835   -0.349850    4.650000    4.662521    4.500000    4.662521
 0     0      ---       0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.534756   -0.166636    4.340000    4.404039    4.250000    4.404039
 1     0      ---       0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.860630   -0.041120    3.870000    3.867798    3.147584    3.867798
 2     0      ---       0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.023532   -0.061245    4.110000    4.125360    4.102416    4.110000

 Harris energy:
 sumev=       -0.853241  val*vef=    -189.493805   sumtv=     188.640564
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -130.258096     utot=   -6534.971175    ehar=   -3304.832068

 srhov:     -6.360832   -168.222947   -174.583779 sumev=   -0.853241   sumtv=  173.730538
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.677273  avg sphere pot= 0.653661  vconst=-0.677273

 smooth rhoves     13.178917   charge     4.185741
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -3.054117   rhomu =   -3.974963  avg vxc =   -0.866699 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=F
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647

 Energy terms:             smooth           local           total
   rhoval*vef             -7.030294      -175.113437      -182.143731
   rhoval*ves            -50.182882      -106.215168      -156.398050
   psnuc*ves              76.540716    -12962.871091    -12886.330375
   utot                   13.178917     -6534.543129     -6521.364212
   rho*exc                -3.054117      -125.587140      -128.641257
   rho*vxc                -3.974963      -166.302313      -170.277276
   valence chg             4.185741         6.814259        11.000000
   core charge            18.000000        -0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Kohn-Sham energy:
 sumtv=      173.730538  sumtc=      3171.756639   ekin=     3345.487177
 rhoep=     -128.641257   utot=     -6521.364212   ehks=    -3304.518292
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrealsmooth= T
 wgtsmooth=   2.7410122234342145E-002
 mixrho:  sought 3 iter from file mixm; read 0.  RMS DQ=3.60e-2
 AMIX: nmix=0 mmix=8  nelts=  2405  beta=1.00000  tm= 5.00000  rmsdel=3.60D-02
 mixrho: add corrections to qcell smrho = -0.23316D-07 -0.29688D-09
 unscreened rms difference:  smooth  0.045849   local  0.019256
   screened rms difference:  smooth  0.045849   local  0.019256   tot  0.035987

 iors  : write restart file (binary, mesh density) 

   it  1  of 12    ehf=      -0.397568   ehk=      -0.083792
h nk=8 bigbas=0 ehf(eV)=-5.4092284 ehk(eV)=-1.1400635 sev(eV)=-11.6090277

 --- BNDFP:  begin iteration 2 of 12
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.633049  avg sphere pot= 0.653661  vconst=-0.633049

 smooth rhoves     12.286148   charge     4.185741
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -3.109207   rhomu =   -4.047536  avg vxc =   -0.858397 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000645

 Energy terms:             smooth           local           total
   rhoval*vef             -7.713664      -175.301697      -183.015361
   rhoval*ves            -48.957953      -107.937662      -156.895615
   psnuc*ves              73.530250    -12964.307473    -12890.777223
   utot                   12.286148     -6536.122568     -6523.836419
   rho*exc                -3.109207      -125.879109      -128.988316
   rho*vxc                -4.047536      -166.689181      -170.736716
   valence chg             4.185741         6.814259        11.000000
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
 ... Done MPI k-loop: elapsed time=   0.0599

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.191168;  11.000000 electrons
         Sum occ. bands:   -7.093537, incl. Bloechl correction: -0.013293
Generating TDOS: efermi, and dos window=   -0.1912  -0.5000   1.3088

 mkrout:  Qtrue      sm,loc       local
   1   10.391031    1.889094    8.501937

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  pnunew: ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.422469   -0.510146    4.662521    4.665863    4.500000    4.665863
 0     0      ---       0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.259043   -0.416122    4.404039    4.372230    4.250000    4.372230
 1     0      ---       0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    9.699171   -0.662756    3.867798    3.901921    3.147584    3.901921
 2     0      ---       0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.010349   -0.615446    4.110000    4.111973    4.102416    4.110000

 Harris energy:
 sumev=       -7.093537  val*vef=    -183.015361   sumtv=     175.921823
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -128.988316     utot=   -6523.836419    ehar=   -3305.146273

 srhov:     -4.200079   -222.070552   -226.270631 sumev=   -7.093537   sumtv=  219.177094
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.400483  avg sphere pot= 0.664691  vconst=-0.400483

 smooth rhoves      4.818396   charge     2.498063
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -1.661083   rhomu =   -2.159363  avg vxc =   -0.746182 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=F
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647

 Energy terms:             smooth           local           total
   rhoval*vef             -3.180685      -205.534191      -208.714876
   rhoval*ves            -36.498704      -142.152533      -178.651237
   psnuc*ves              46.135497    -12991.549490    -12945.413994
   utot                    4.818396     -6566.851012     -6562.032615
   rho*exc                -1.661083      -131.002602      -132.663685
   rho*vxc                -2.159363      -173.451554      -175.610917
   valence chg             2.498063         8.501937        11.000000
   core charge            18.000000        -0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=      219.177094  sumtc=      3171.756639   ekin=     3390.933733
 rhoep=     -132.663685   utot=     -6562.032615   ehks=    -3303.762566
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 wgtsmooth=   2.7410122234342145E-002
 mixrho:  sought 3 iter from file mixm; read 1.  RMS DQ=9.09e-2  last it=3.60e-2
 AMIX: nmix=1 mmix=8  nelts=  2405  beta=1.00000  tm= 5.00000  rmsdel=9.09D-02
   tj: 0.82101
 mixrho: add corrections to qcell smrho = -0.11103D-07 -0.14137D-09
 unscreened rms difference:  smooth  0.024246   local  0.064870
   screened rms difference:  smooth  0.024246   local  0.064870   tot  0.090874

 iors  : write restart file (binary, mesh density) 

   it  2  of 12    ehf=      -0.711773   ehk=       0.671934
 From last iter    ehf=      -0.397568   ehk=      -0.083792
 diffe(q)= -0.314205 (0.090874)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf(eV)=-9.6842426 ehk(eV)=9.1421935 sev(eV)=-96.5132507

 --- BNDFP:  begin iteration 3 of 12
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.604915  avg sphere pot= 0.655635  vconst=-0.604915

 smooth rhoves     10.955967   charge     3.883666
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -2.817548   rhomu =   -3.666817  avg vxc =   -0.843500 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000597

 Energy terms:             smooth           local           total
   rhoval*vef             -6.576698      -181.158821      -187.735520
   rhoval*ves            -47.654907      -113.403665      -161.058572
   psnuc*ves              69.566842    -12968.686245    -12899.119404
   utot                   10.955967     -6541.044955     -6530.088988
   rho*exc                -2.817548      -126.698751      -129.516299
   rho*vxc                -3.666817      -167.770394      -171.437211
   valence chg             3.883666         7.116334        11.000000
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
 ... Done MPI k-loop: elapsed time=   0.0592

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.121679;  11.000000 electrons
         Sum occ. bands:   -4.729081, incl. Bloechl correction: -0.011778
Generating TDOS: efermi, and dos window=   -0.1217  -0.5000   1.3783

 mkrout:  Qtrue      sm,loc       local
   1   10.274203    2.258616    8.015587

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  pnunew: ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.448182   -0.471870    4.665863    4.661731    4.500000    4.661731
 0     0      ---       0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.303498   -0.374406    4.372230    4.371002    4.250000    4.371002
 1     0      ---       0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    9.508935   -0.428368    3.901921    3.890746    3.147584    3.890746
 2     0      ---       0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.013587   -0.419417    4.110000    4.116636    4.102416    4.110000

 Harris energy:
 sumev=       -4.729081  val*vef=    -187.735520   sumtv=     183.006439
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.516299     utot=   -6530.088988    ehar=   -3304.842208

 srhov:     -4.839652   -203.034353   -207.874005 sumev=   -4.729081   sumtv=  203.144925
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.476735  avg sphere pot= 0.666104  vconst=-0.476735

 smooth rhoves      6.772582   charge     2.984414
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -2.045659   rhomu =   -2.660386  avg vxc =   -0.784527 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=F
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647

 Energy terms:             smooth           local           total
   rhoval*vef             -4.183503      -195.306394      -199.489896
   rhoval*ves            -41.181371      -129.701838      -170.883209
   psnuc*ves              54.726535    -12980.016447    -12925.289912
   utot                    6.772582     -6554.859143     -6548.086561
   rho*exc                -2.045659      -129.270665      -131.316325
   rho*vxc                -2.660386      -171.163634      -173.824020
   valence chg             2.984414         8.015587        11.000000
   core charge            18.000000        -0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Kohn-Sham energy:
 sumtv=      203.144925  sumtc=      3171.756639   ekin=     3374.901564
 rhoep=     -131.316325   utot=     -6548.086561   ehks=    -3304.501321
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 wgtsmooth=   2.7410122234342145E-002
 mixrho:  sought 3 iter from file mixm; read 2.  RMS DQ=4.40e-2  last it=9.09e-2
 AMIX: nmix=2 mmix=8  nelts=  2405  beta=1.00000  tm= 5.00000  rmsdel=4.40D-02
   tj:-1.03927  -0.10625
 mixrho: add corrections to qcell smrho = -0.30441D-06 -0.38759D-08
 unscreened rms difference:  smooth  0.012620   local  0.032800
   screened rms difference:  smooth  0.012620   local  0.032800   tot  0.044025

 iors  : write restart file (binary, mesh density) 

   it  3  of 12    ehf=      -0.407708   ehk=      -0.066821
 From last iter    ehf=      -0.711773   ehk=       0.671934
 diffe(q)=  0.304065 (0.044025)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf(eV)=-5.54719 ehk(eV)=-.9091513 sev(eV)=-64.3429236

 --- BNDFP:  begin iteration 4 of 12
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.545146  avg sphere pot= 0.668895  vconst=-0.545146

 smooth rhoves      8.671497   charge     3.362221
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -2.347351   rhomu =   -3.053460  avg vxc =   -0.813370 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000549

 Energy terms:             smooth           local           total
   rhoval*vef             -4.920304      -185.951724      -190.872028
   rhoval*ves            -44.678446      -118.834502      -163.512948
   psnuc*ves              62.021441    -12969.566049    -12907.544608
   utot                    8.671497     -6544.200275     -6535.528778
   rho*exc                -2.347351      -127.805101      -130.152452
   rho*vxc                -3.053460      -169.227031      -172.280492
   valence chg             3.362221         7.637779        11.000000
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
 ... Done MPI k-loop: elapsed time=   0.0608

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.066062;  11.000000 electrons
         Sum occ. bands:   -1.734814, incl. Bloechl correction: -0.007558
Generating TDOS: efermi, and dos window=    0.0661  -0.5000   1.5661

 mkrout:  Qtrue      sm,loc       local
   1   10.006010    2.935456    7.070553

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  pnunew: ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.497640   -0.390258    4.661731    4.662100    4.500000    4.662100
 0     0      ---       0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.475293   -0.227619    4.371002    4.395126    4.250000    4.395126
 1     0      ---       0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    9.011970   -0.128804    3.890746    3.872855    3.147584    3.872855
 2     0      ---       0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021106   -0.143522    4.110000    4.123593    4.102416    4.110000

 Harris energy:
 sumev=       -1.734814  val*vef=    -190.872028   sumtv=     189.137214
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -130.152452     utot=   -6535.528778    ehar=   -3304.787376

 srhov:     -5.820463   -174.784474   -180.604937 sumev=   -1.734814   sumtv=  178.870123
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.633468  avg sphere pot= 0.658285  vconst=-0.633468

 smooth rhoves     11.622224   charge     3.929447
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -2.832043   rhomu =   -3.685389  avg vxc =   -0.850506 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=F
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647

 Energy terms:             smooth           local           total
   rhoval*vef             -6.380268      -178.810105      -185.190373
   rhoval*ves            -48.593219      -110.325860      -158.919079
   psnuc*ves              71.837667    -12965.241024    -12893.403357
   utot                   11.622224     -6537.783442     -6526.161218
   rho*exc                -2.832043      -126.303874      -129.135916
   rho*vxc                -3.685389      -167.247614      -170.933003
   valence chg             3.929447         7.070553        11.000000
   core charge            18.000000        -0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Kohn-Sham energy:
 sumtv=      178.870123  sumtc=      3171.756639   ekin=     3350.626762
 rhoep=     -129.135916   utot=     -6526.161218   ehks=    -3304.670372
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 wgtsmooth=   2.7410122234342145E-002
 mixrho:  sought 3 iter from file mixm; read 3.  RMS DQ=2.43e-2  last it=4.40e-2
 AMIX: nmix=3 mmix=8  nelts=  2405  beta=1.00000  tm= 5.00000  rmsdel=2.43D-02
   tj: 0.76321  -0.24310  -0.00278
 mixrho: add corrections to qcell smrho = -0.15776D-06 -0.20086D-08
 unscreened rms difference:  smooth  0.007483   local  0.019601
   screened rms difference:  smooth  0.007483   local  0.019601   tot  0.024318

 iors  : write restart file (binary, mesh density) 

   it  4  of 12    ehf=      -0.352876   ehk=      -0.235872
 From last iter    ehf=      -0.407708   ehk=      -0.066821
 diffe(q)=  0.054831 (0.024318)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf(eV)=-4.8011663 ehk(eV)=-3.2092281 sev(eV)=-23.6035374

 --- BNDFP:  begin iteration 5 of 12
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.570017  avg sphere pot= 0.662708  vconst=-0.570017

 smooth rhoves      9.528920   charge     3.555451
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -2.515365   rhomu =   -3.272549  avg vxc =   -0.825511 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000560

 Energy terms:             smooth           local           total
   rhoval*vef             -5.478201      -184.923907      -190.402109
   rhoval*ves            -45.930153      -117.336551      -163.266704
   psnuc*ves              64.987992    -12970.137903    -12905.149910
   utot                    9.528920     -6543.737227     -6534.208307
   rho*exc                -2.515365      -127.430514      -129.945879
   rho*vxc                -3.272549      -168.734285      -172.006835
   valence chg             3.555451         7.444549        11.000000
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
 ... Done MPI k-loop: elapsed time=   0.0618

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.014429;  11.000000 electrons
         Sum occ. bands:   -2.765369, incl. Bloechl correction: -0.009079
Generating TDOS: efermi, and dos window=   -0.0144  -0.5000   1.4856

 mkrout:  Qtrue      sm,loc       local
   1   10.107144    2.683239    7.423905

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  pnunew: ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.474968   -0.427641    4.662100    4.659693    4.500000    4.659693
 0     0      ---       0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.409567   -0.282707    4.395126    4.386369    4.250000    4.386369
 1     0      ---       0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    9.204393   -0.232830    3.872855    3.878975    3.147584    3.878975
 2     0      ---       0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.018215   -0.237482    4.110000    4.121228    4.102416    4.110000

 Harris energy:
 sumev=       -2.765369  val*vef=    -190.402109   sumtv=     187.636739
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.945879     utot=   -6534.208307    ehar=   -3304.760807

 srhov:     -5.487174   -184.461023   -189.948197 sumev=   -2.765369   sumtv=  187.182828
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.574641  avg sphere pot= 0.662487  vconst=-0.574641

 smooth rhoves      9.656570   charge     3.576095
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -2.531012   rhomu =   -3.292924  avg vxc =   -0.827227 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=F
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647

 Energy terms:             smooth           local           total
   rhoval*vef             -5.513162      -184.684059      -190.197221
   rhoval*ves            -46.124076      -116.991394      -163.115470
   psnuc*ves              65.437216    -12969.929462    -12904.492246
   utot                    9.656570     -6543.460428     -6533.803858
   rho*exc                -2.531012      -127.365102      -129.896114
   rho*vxc                -3.292924      -168.647975      -171.940899
   valence chg             3.576095         7.423905        11.000000
   core charge            18.000000        -0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Kohn-Sham energy:
 sumtv=      187.182828  sumtc=      3171.756639   ekin=     3358.939467
 rhoep=     -129.896114   utot=     -6533.803858   ehks=    -3304.760505
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 wgtsmooth=   2.7410122234342145E-002
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=9.70e-4  last it=2.43e-2
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: nmix=2 mmix=8  nelts=  2405  beta=1.00000  tm= 5.00000  rmsdel=9.70D-04
   tj:-0.04331  -0.00109
 mixrho: add corrections to qcell smrho = -0.12110D-07 -0.15419D-09
 unscreened rms difference:  smooth  0.000347   local  0.000768
   screened rms difference:  smooth  0.000347   local  0.000768   tot  0.000970

 iors  : write restart file (binary, mesh density) 

   it  5  of 12    ehf=      -0.326307   ehk=      -0.326005
 From last iter    ehf=      -0.352876   ehk=      -0.235872
 diffe(q)=  0.026570 (0.000970)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf(eV)=-4.4396633 ehk(eV)=-4.435552 sev(eV)=-37.6250635

 --- BNDFP:  begin iteration 6 of 12
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.571791  avg sphere pot= 0.662665  vconst=-0.571791

 smooth rhoves      9.571667   charge     3.561435
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -2.519250   rhomu =   -3.277600  avg vxc =   -0.826128 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000558

 Energy terms:             smooth           local           total
   rhoval*vef             -5.484067      -184.917312      -190.401379
   rhoval*ves            -45.999700      -117.282483      -163.282183
   psnuc*ves              65.143033    -12970.133038    -12904.990005
   utot                    9.571667     -6543.707760     -6534.136094
   rho*exc                -2.519250      -127.411797      -129.931048
   rho*vxc                -3.277600      -168.709607      -171.987206
   valence chg             3.561435         7.438565        11.000000
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
 ... Done MPI k-loop: elapsed time=   0.0618

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.020515;  11.000000 electrons
         Sum occ. bands:   -2.851611, incl. Bloechl correction: -0.009203
Generating TDOS: efermi, and dos window=   -0.0205  -0.5000   1.4795

 mkrout:  Qtrue      sm,loc       local
   1   10.116130    2.663360    7.452770

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  pnunew: ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.474100   -0.430105    4.659693    4.659887    4.500000    4.659887
 0     0      ---       0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.402209   -0.289089    4.386369    4.385081    4.250000    4.385081
 1     0      ---       0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    9.221874   -0.241356    3.878975    3.879541    3.147584    3.879541
 2     0      ---       0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.017946   -0.246180    4.110000    4.121017    4.102416    4.110000

 Harris energy:
 sumev=       -2.851611  val*vef=    -190.401379   sumtv=     187.549768
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.931048     utot=   -6534.136094    ehar=   -3304.760734

 srhov:     -5.458720   -185.248220   -190.706941 sumev=   -2.851611   sumtv=  187.855330
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.569663  avg sphere pot= 0.662829  vconst=-0.569663

 smooth rhoves      9.501721   charge     3.547230
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -2.506989   rhomu =   -3.261612  avg vxc =   -0.825222 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=F
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647

 Energy terms:             smooth           local           total
   rhoval*vef             -5.446393      -185.134697      -190.581090
   rhoval*ves            -45.901266      -117.530591      -163.431857
   psnuc*ves              64.904709    -12970.299804    -12905.395095
   utot                    9.501721     -6543.915197     -6534.413476
   rho*exc                -2.506989      -127.452153      -129.959142
   rho*vxc                -3.261612      -168.762844      -172.024457
   valence chg             3.547230         7.452770        11.000000
   core charge            18.000000        -0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Kohn-Sham energy:
 sumtv=      187.855330  sumtc=      3171.756639   ekin=     3359.611969
 rhoep=     -129.959142   utot=     -6534.413476   ehks=    -3304.760649
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 wgtsmooth=   2.7410122234342145E-002
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=6.83e-4  last it=9.70e-4
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: nmix=2 mmix=8  nelts=  2405  beta=1.00000  tm= 5.00000  rmsdel=6.83D-04
   tj: 0.33583   0.00509
 mixrho: add corrections to qcell smrho = -0.34565D-07 -0.44010D-09
 unscreened rms difference:  smooth  0.000201   local  0.000518
   screened rms difference:  smooth  0.000201   local  0.000518   tot  0.000683

 iors  : write restart file (binary, mesh density) 

   it  6  of 12    ehf=      -0.326234   ehk=      -0.326149
 From last iter    ehf=      -0.326307   ehk=      -0.326005
 diffe(q)=  0.000072 (0.000683)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf(eV)=-4.4386771 ehk(eV)=-4.4375148 sev(eV)=-38.7984487

 --- BNDFP:  begin iteration 7 of 12
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.571487  avg sphere pot= 0.662691  vconst=-0.571487

 smooth rhoves      9.560764   charge     3.558869
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -2.516925   rhomu =   -3.274566  avg vxc =   -0.825985 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000558

 Energy terms:             smooth           local           total
   rhoval*vef             -5.475824      -184.951073      -190.426897
   rhoval*ves            -45.985079      -117.318257      -163.303336
   psnuc*ves              65.106608    -12970.154285    -12905.047678
   utot                    9.560764     -6543.736271     -6534.175507
   rho*exc                -2.516925      -127.418209      -129.935134
   rho*vxc                -3.274566      -168.718059      -171.992625
   valence chg             3.558869         7.441131        11.000000
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
 ... Done MPI k-loop: elapsed time=   0.0624

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.019274;  11.000000 electrons
         Sum occ. bands:   -2.833626, incl. Bloechl correction: -0.009175
Generating TDOS: efermi, and dos window=   -0.0193  -0.5000   1.4807

 mkrout:  Qtrue      sm,loc       local
   1   10.114510    2.667674    7.446837

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  pnunew: ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.474528   -0.429518    4.659887    4.659941    4.500000    4.659941
 0     0      ---       0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.403127   -0.288311    4.385081    4.385194    4.250000    4.385194
 1     0      ---       0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    9.218864   -0.239543    3.879541    3.879438    3.147584    3.879438
 2     0      ---       0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.017991   -0.244580    4.110000    4.121059    4.102416    4.110000

 Harris energy:
 sumev=       -2.833626  val*vef=    -190.426897   sumtv=     187.593271
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.935134     utot=   -6534.175507    ehar=   -3304.760730

 srhov:     -5.465520   -185.071917   -190.537437 sumev=   -2.833626   sumtv=  187.703811
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.570633  avg sphere pot= 0.662778  vconst=-0.570633

 smooth rhoves      9.532766   charge     3.553163
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -2.511985   rhomu =   -3.268125  avg vxc =   -0.825624 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=F
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647

 Energy terms:             smooth           local           total
   rhoval*vef             -5.460586      -185.029146      -190.489732
   rhoval*ves            -45.945769      -117.409172      -163.354940
   psnuc*ves              65.011301    -12970.207442    -12905.196141
   utot                    9.532766     -6543.808307     -6534.275540
   rho*exc                -2.511985      -127.433643      -129.945628
   rho*vxc                -3.268125      -168.738413      -172.006538
   valence chg             3.553163         7.446837        11.000000
   core charge            18.000000        -0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Kohn-Sham energy:
 sumtv=      187.703811  sumtc=      3171.756639   ekin=     3359.460450
 rhoep=     -129.945628   utot=     -6534.275540   ehks=    -3304.760718
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 wgtsmooth=   2.7410122234342145E-002
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=2.57e-4  last it=6.83e-4
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: nmix=2 mmix=8  nelts=  2405  beta=1.00000  tm= 5.00000  rmsdel=2.57D-04
   tj:-0.28729   0.11022
 mixrho: add corrections to qcell smrho = -0.43586D-07 -0.55496D-09
 unscreened rms difference:  smooth  0.000082   local  0.000202
   screened rms difference:  smooth  0.000082   local  0.000202   tot  0.000257

 iors  : write restart file (binary, mesh density) 

   it  7  of 12    ehf=      -0.326230   ehk=      -0.326218
 From last iter    ehf=      -0.326234   ehk=      -0.326149
 diffe(q)=  0.000004 (0.000257)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf(eV)=-4.4386265 ehk(eV)=-4.4384605 sev(eV)=-38.5537531

 --- BNDFP:  begin iteration 8 of 12
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.571307  avg sphere pot= 0.662732  vconst=-0.571307

 smooth rhoves      9.554408   charge     3.557396
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -2.515588   rhomu =   -3.272822  avg vxc =   -0.825903 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000558

 Energy terms:             smooth           local           total
   rhoval*vef             -5.471168      -184.960527      -190.431695
   rhoval*ves            -45.976525      -117.329935      -163.306459
   psnuc*ves              65.085341    -12970.151648    -12905.066307
   utot                    9.554408     -6543.740791     -6534.186383
   rho*exc                -2.515588      -127.421102      -129.936690
   rho*vxc                -3.272822      -168.721864      -171.994687
   valence chg             3.557396         7.442604        11.000000
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
 ... Done MPI k-loop: elapsed time=   0.0626

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018741;  11.000000 electrons
         Sum occ. bands:   -2.825995, incl. Bloechl correction: -0.009164
Generating TDOS: efermi, and dos window=   -0.0187  -0.5000   1.4813

 mkrout:  Qtrue      sm,loc       local
   1   10.113771    2.669475    7.444296

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  pnunew: ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.474661   -0.429287    4.659941    4.659944    4.500000    4.659944
 0     0      ---       0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.403634   -0.287883    4.385194    4.385268    4.250000    4.385268
 1     0      ---       0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    9.217464   -0.238781    3.879438    3.879391    3.147584    3.879391
 2     0      ---       0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.018012   -0.243869    4.110000    4.121076    4.102416    4.110000

 Harris energy:
 sumev=       -2.825995  val*vef=    -190.431695   sumtv=     187.605700
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.936690     utot=   -6534.186383    ehar=   -3304.760733

 srhov:     -5.468001   -184.999057   -190.467058 sumev=   -2.825995   sumtv=  187.641063
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.571059  avg sphere pot= 0.662754  vconst=-0.571059

 smooth rhoves      9.546210   charge     3.555704
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -2.514113   rhomu =   -3.270898  avg vxc =   -0.825798 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=F
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647

 Energy terms:             smooth           local           total
   rhoval*vef             -5.466570      -184.986147      -190.452717
   rhoval*ves            -45.965088      -117.358947      -163.324035
   psnuc*ves              65.057508    -12970.170512    -12905.113004
   utot                    9.546210     -6543.764730     -6534.218519
   rho*exc                -2.514113      -127.425803      -129.939915
   rho*vxc                -3.270898      -168.728065      -171.998963
   valence chg             3.555704         7.444296        11.000000
   core charge            18.000000        -0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Kohn-Sham energy:
 sumtv=      187.641063  sumtc=      3171.756639   ekin=     3359.397703
 rhoep=     -129.939915   utot=     -6534.218519   ehks=    -3304.760732
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 wgtsmooth=   2.7410122234342145E-002
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=7.98e-5  last it=2.57e-4
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=  2405  beta=1.00000  tm= 5.00000  rmsdel=7.98D-05
   tj:-0.44767
 mixrho: add corrections to qcell smrho = -0.45494D-07 -0.57926D-09
 unscreened rms difference:  smooth  0.000026   local  0.000061
   screened rms difference:  smooth  0.000026   local  0.000061   tot  0.000080

 iors  : write restart file (binary, mesh density) 

   it  8  of 12    ehf=      -0.326233   ehk=      -0.326232
 From last iter    ehf=      -0.326230   ehk=      -0.326218
 diffe(q)= -0.000002 (0.000080)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf(eV)=-4.4386604 ehk(eV)=-4.438643 sev(eV)=-38.4499223

 --- BNDFP:  begin iteration 9 of 12
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.571245  avg sphere pot= 0.662743  vconst=-0.571245

 smooth rhoves      9.552147   charge     3.556842
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -2.515072   rhomu =   -3.272149  avg vxc =   -0.825875 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000558

 Energy terms:             smooth           local           total
   rhoval*vef             -5.469316      -184.966836      -190.436152
   rhoval*ves            -45.973562      -117.336602      -163.310164
   psnuc*ves              65.077856    -12970.154088    -12905.076232
   utot                    9.552147     -6543.745345     -6534.193198
   rho*exc                -2.515072      -127.422323      -129.937396
   rho*vxc                -3.272149      -168.723473      -171.995623
   valence chg             3.556842         7.443158        11.000000
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
 ... Done MPI k-loop: elapsed time=   0.0633

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018528;  11.000000 electrons
         Sum occ. bands:   -2.822934, incl. Bloechl correction: -0.009159
Generating TDOS: efermi, and dos window=   -0.0185  -0.5000   1.4815

 mkrout:  Qtrue      sm,loc       local
   1   10.113462    2.670237    7.443225

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  pnunew: ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.474718   -0.429201    4.659944    4.659947    4.500000    4.659947
 0     0      ---       0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.403844   -0.287714    4.385268    4.385300    4.250000    4.385300
 1     0      ---       0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    9.216880   -0.238474    3.879391    3.879372    3.147584    3.879372
 2     0      ---       0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.018021   -0.243583    4.110000    4.121084    4.102416    4.110000

 Harris energy:
 sumev=       -2.822934  val*vef=    -190.436152   sumtv=     187.613219
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.937396     utot=   -6534.193198    ehar=   -3304.760735

 srhov:     -5.469141   -184.967426   -190.436566 sumev=   -2.822934   sumtv=  187.613633
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.571238  avg sphere pot= 0.662745  vconst=-0.571238

 smooth rhoves      9.551870   charge     3.556775
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -2.515010   rhomu =   -3.272068  avg vxc =   -0.825871 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=F
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647

 Energy terms:             smooth           local           total
   rhoval*vef             -5.469098      -184.967195      -190.436294
   rhoval*ves            -45.973204      -117.337041      -163.310245
   psnuc*ves              65.076944    -12970.153809    -12905.076865
   utot                    9.551870     -6543.745425     -6534.193555
   rho*exc                -2.515010      -127.422442      -129.937452
   rho*vxc                -3.272068      -168.723629      -171.995698
   valence chg             3.556775         7.443225        11.000000
   core charge            18.000000        -0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Kohn-Sham energy:
 sumtv=      187.613633  sumtc=      3171.756639   ekin=     3359.370272
 rhoep=     -129.937452   utot=     -6534.193555   ehks=    -3304.760735
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 wgtsmooth=   2.7410122234342145E-002
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=2.35e-6  last it=7.98e-5
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=  2405  beta=1.00000  tm= 5.00000  rmsdel=2.35D-06
   tj:-0.02441
 mixrho: add corrections to qcell smrho = -0.45484D-07 -0.57913D-09
 unscreened rms difference:  smooth  0.000001   local  0.000002
   screened rms difference:  smooth  0.000001   local  0.000002   tot  0.000002

 iors  : write restart file (binary, mesh density) 

   it  9  of 12    ehf=      -0.326235   ehk=      -0.326235
 From last iter    ehf=      -0.326233   ehk=      -0.326232
 diffe(q)= -0.000002 (0.000002)    tol= 0.000010 (0.000010)   more=F
c nk=8 bigbas=0 ehf(eV)=-4.4386899 ehk(eV)=-4.4386878 sev(eV)=-38.4082693
 >>      0.96   exit  lmfp            0.90
OK! end of LMF ======================
===START LMF   =====================================
 mpisize=           4
  bndfp (warning): no sigm file found ... LDA calculation only

mmm === MTO setting ===
mmm ispec lmxb lpzex nkapii nkaphh=    1    4    0    2    2
mmm rsmh1    1  2.50  2.50  1.00  0.00  0.00
mmm   eh1    1 -0.01 -0.01 -0.01 -0.01 -0.01
mmm rsmh2    1  1.30  0.00  1.00  1.30  0.00
mmm  eh2     1 -1.00 -1.00 -1.00 -0.01 -0.01
mmm pz       1  5.50  5.50  4.50  0.00  0.00
mmm lh       1  2  3  2

                Plat                                  Qlat
   0.000000   0.500000   0.500000       -1.000000   1.000000   1.000000
   0.500000   0.000000   0.500000        1.000000  -1.000000   1.000000
   0.500000   0.500000   0.000000        1.000000   1.000000  -1.000000
  Cell vol= 78.538660

 LATTC:  as= 2.000   tol=1.00E-12   alat= 6.79800   awald= 0.467
         r1=  2.252   nkd= 201      q1=  6.871   nkg= 331

 SGROUP: 1 symmetry operations from 0 generators
 ADDBAS: basis is already complete --- no sites added
 SYMLAT: Bravais system is cubic with 48 symmetry operations.
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-5
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
  21  r2(1,1,-0)                         
  22  m(1,1,-0)                          
  23  r2(1,-0,-1)                        
  24  m(1,-0,-1)                         
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
  47  r2(-0,1,1)                         
  48  m(-0,1,1)                          
 GROUPG: the following are sufficient to generate the space group:
 Generator(cart): i*r3(1,1,-1) r4x
 Generator(frac): i*r3(1,1,-1) r4x
 MKSYM:  found 48 space group operations ... includes inversion
 
 BZMESH:     60 irreducible QP from    8   8   8 shift=TTT
 TETIRR: sorting 3072 tetrahedra ...
 264 inequivalent ones found
 >> level: 1  CPUsec=      0.11  enter lmfp
 gen_hamindex: not readin QGpsi.

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv  kmxv foca   rfoca
 A        2.311  0.925    4    4         4  0.578  1.156    15    1   0.925

 MSHSIZ: mesh has 11 x 11 x 11 divisions; length 0.437, 0.437, 0.437
         generated from gmax = 9.0 a.u. : 941 vectors of 1331 (70%)
 goto end of reading rst or atm           1
 goto end of reading rst or atm           1
 goto end of reading rst or atm           1

 GVLST2: gmax = 9.0 a.u. created 941 vectors of 1331 (70%)
         (input) mesh has 11 x 11 x 11 divisions; length 0.437, 0.437, 0.437
 SGVSYM: 41 symmetry stars found for 941 reciprocal lattice vectors

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol= 1.00E-06
 spec      l    rsm    eh     gmax    last term    cutoff
  A        0*   2.50  -0.01   2.974    2.30E-05      27 
  A        1    2.50  -0.01   3.093    1.29E-06      51 
  A        2    1.00  -0.01   8.508    1.16E-06     869 
  A        0    1.30  -1.00   5.718    2.28E-06     259 
  A        2    1.00  -1.00   8.508    1.16E-06     869 
  A        3    1.30  -0.01   6.806    2.09E-06     411 
 goto end of reading rst or atm           1

 iors  : read restart file (binary, mesh density) 
         use from  restart file: ef window, positions, pnu 
         ignore in restart file: *
         site   1, species A       : augmentation lmax changed from 3 to 4
         site   1, species A       : inflate local density from nlm= 16 to 25
 end of reading rst or atm
 m_qplistinit:start
 end of reading rst or atm
 end of reading rst or atm
 end of reading rst or atm
m_qplist_qspdivider: rank,(iqini,ispini),(iqend,ispend)=    2, (   31 1),  (   45 1)
m_qplist_qspdivider: rank,(iqini,ispini),(iqend,ispend)=    3, (   46 1),  (   60 1)
m_qplist_qspdivider: rank,(iqini,ispini),(iqend,ispend)=    0, (    1 1),  (   15 1)
m_qplist_qspdivider: rank,(iqini,ispini),(iqend,ispend)=    1, (   16 1),  (   30 1)

 Basis, after reading restart file
 site spec        pos (Cartesian coordinates)         pos (multiples of plat)
   1  A         0.000000   0.000000   0.000000    0.000000   0.000000   0.000000

 --- BNDFP:  begin iteration 1 of 12
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.571242  avg sphere pot= 0.662745  vconst=-0.571242

 smooth rhoves      9.552001   charge     3.556801
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -2.515033   rhomu =   -3.272098  avg vxc =   -0.825873 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000558

 Energy terms:             smooth           local           total
   rhoval*vef             -5.469166      -184.966732      -190.435899
   rhoval*ves            -45.973388      -117.336522      -163.309910
   psnuc*ves              65.077390    -12970.153415    -12905.076025
   utot                    9.552001     -6543.744968     -6534.192967
   rho*exc                -2.515033      -127.422363      -129.937395
   rho*vxc                -3.272098      -168.723525      -171.995622
   valence chg             3.556801         7.443199        11.000000
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
 ... Done MPI k-loop: elapsed time=   0.1069

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018806;  11.000000 electrons
         Sum occ. bands:   -2.825415, incl. Bloechl correction: -0.009153
Generating TDOS: efermi, and dos window=   -0.0188  -0.5000   1.4812
  mmmmm m_bandcal_2nd

 mkrout:  Qtrue      sm,loc       local
   1   10.128183    2.833411    7.294772

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  pnunew: ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.477550   -0.430054    4.659947    4.659574    4.500000    4.659574
 0     0      ---       0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.409683   -0.289034    4.385300    4.384836    4.250000    4.384836
 1     0      ---       0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    9.215375   -0.238670    3.879372    3.879149    3.147584    3.879149
 2     0      ---       0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021180   -0.245040    4.110000    4.121040    4.102416    4.110000
 4     1    0.004395   -0.241226    5.100000    5.085239    5.077979    5.100000

 Harris energy:
 sumev=       -2.825415  val*vef=    -190.435899   sumtv=     187.610483
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.937395     utot=   -6534.192967    ehar=   -3304.763240

 srhov:     -5.929810   -184.355896   -190.285706 sumev=   -2.825415   sumtv=  187.460291
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.581228  avg sphere pot= 0.653085  vconst=-0.581228

 smooth rhoves     10.109125   charge     3.705228
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -2.657844   rhomu =   -3.458503  avg vxc =   -0.832964 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=F
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647

 Energy terms:             smooth           local           total
   rhoval*vef             -6.009813      -184.303400      -190.313213
   rhoval*ves            -46.615397      -116.581954      -163.197351
   psnuc*ves              66.833647    -12971.741471    -12904.907824
   utot                   10.109125     -6544.161712     -6534.052587
   rho*exc                -2.657844      -127.269667      -129.927512
   rho*vxc                -3.458503      -168.523981      -171.982484
   valence chg             3.705228         7.294772        11.000000
   core charge            18.000000        -0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Kohn-Sham energy:
 sumtv=      187.460291  sumtc=      3171.756639   ekin=     3359.216930
 rhoep=     -129.927512   utot=     -6534.052587   ehks=    -3304.763169
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrealsmooth= T
 wgtsmooth=   2.7410122234342145E-002
 mixrho:  sought 3 iter from file mixm; read 0.  RMS DQ=4.03e-3
 AMIX: nmix=0 mmix=8  nelts=  2567  beta=1.00000  tm= 5.00000  rmsdel=4.03D-03
 mixrho: add corrections to qcell smrho = -0.28596D-07 -0.36410D-09
 unscreened rms difference:  smooth  0.003489   local  0.004828
   screened rms difference:  smooth  0.003489   local  0.004828   tot  0.004029

 iors  : write restart file (binary, mesh density) 

   it  1  of 12    ehf=      -0.328740   ehk=      -0.328669
i nk=8 bigbas=1 pwmode=0 oveps=0 ehf(eV)=-4.4727703 ehk(eV)=-4.4717991 sev(eV)=-38.4420351

 --- BNDFP:  begin iteration 2 of 12
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.581163  avg sphere pot= 0.653085  vconst=-0.581163

 smooth rhoves     10.109074   charge     3.705228
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -2.657570   rhomu =   -3.458141  avg vxc =   -0.833027 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000127

 Energy terms:             smooth           local           total
   rhoval*vef             -6.009705      -184.305076      -190.314780
   rhoval*ves            -46.616015      -116.582921      -163.198936
   psnuc*ves              66.834162    -12971.742020    -12904.907858
   utot                   10.109074     -6544.162471     -6534.053397
   rho*exc                -2.657570      -127.269970      -129.927540
   rho*vxc                -3.458141      -168.524383      -171.982525
   valence chg             3.705228         7.294772        11.000000
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
 ... Done MPI k-loop: elapsed time=   0.1141

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.020713;  11.000000 electrons
         Sum occ. bands:   -2.853660, incl. Bloechl correction: -0.009204
Generating TDOS: efermi, and dos window=   -0.0207  -0.5000   1.4793
  mmmmm m_bandcal_2nd

 mkrout:  Qtrue      sm,loc       local
   1   10.131732    2.824840    7.306891

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  pnunew: ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.476763   -0.430575    4.659574    4.659654    4.500000    4.659654
 0     0      ---       0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.407335   -0.290439    4.384836    4.384698    4.250000    4.384698
 1     0      ---       0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    9.222235   -0.241551    3.879149    3.878964    3.147584    3.878964
 2     0      ---       0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021025   -0.247691    4.110000    4.120990    4.102416    4.110000
 4     1    0.004375   -0.244081    5.100000    5.085219    5.077979    5.100000

 Harris energy:
 sumev=       -2.853660  val*vef=    -190.314780   sumtv=     187.461120
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.927540     utot=   -6534.053397    ehar=   -3304.763177

 srhov:     -5.992884   -184.684129   -190.677013 sumev=   -2.853660   sumtv=  187.823353
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.579170  avg sphere pot= 0.653068  vconst=-0.579170

 smooth rhoves     10.043246   charge     3.693109
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -2.647612   rhomu =   -3.445163  avg vxc =   -0.832142 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=F
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647

 Energy terms:             smooth           local           total
   rhoval*vef             -5.979839      -184.558522      -190.538361
   rhoval*ves            -46.526457      -116.862469      -163.388925
   psnuc*ves              66.612948    -12971.992208    -12905.379260
   utot                   10.043246     -6544.427338     -6534.384093
   rho*exc                -2.647612      -127.311359      -129.958971
   rho*vxc                -3.445163      -168.579039      -172.024202
   valence chg             3.693109         7.306891        11.000000
   core charge            18.000000        -0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Kohn-Sham energy:
 sumtv=      187.823353  sumtc=      3171.756639   ekin=     3359.579993
 rhoep=     -129.958971   utot=     -6534.384093   ehks=    -3304.763071
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 wgtsmooth=   2.7410122234342145E-002
 mixrho:  sought 3 iter from file mixm; read 1.  RMS DQ=7.06e-4  last it=4.03e-3
 AMIX: nmix=1 mmix=8  nelts=  2567  beta=1.00000  tm= 5.00000  rmsdel=7.06D-04
   tj: 0.08237
 mixrho: add corrections to qcell smrho = -0.27092D-07 -0.34495D-09
 unscreened rms difference:  smooth  0.000177   local  0.000504
   screened rms difference:  smooth  0.000177   local  0.000504   tot  0.000706

 iors  : write restart file (binary, mesh density) 

   it  2  of 12    ehf=      -0.328677   ehk=      -0.328571
 From last iter    ehf=      -0.328740   ehk=      -0.328669
 diffe(q)=  0.000063 (0.000706)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=1 pwmode=0 oveps=0 ehf(eV)=-4.4719128 ehk(eV)=-4.4704729 sev(eV)=-38.8263262

 --- BNDFP:  begin iteration 3 of 12
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.579969  avg sphere pot= 0.653070  vconst=-0.579969

 smooth rhoves     10.060686   charge     3.694107
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -2.647618   rhomu =   -3.445160  avg vxc =   -0.832355 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000134

 Energy terms:             smooth           local           total
   rhoval*vef             -5.972883      -184.535218      -190.508101
   rhoval*ves            -46.555155      -116.812162      -163.367317
   psnuc*ves              66.676528    -12971.948724    -12905.272195
   utot                   10.060686     -6544.380443     -6534.319756
   rho*exc                -2.647618      -127.303347      -129.950964
   rho*vxc                -3.445160      -168.568437      -172.013597
   valence chg             3.694107         7.305893        11.000000
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
 ... Done MPI k-loop: elapsed time=   0.1143

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.014055;  11.000000 electrons
         Sum occ. bands:   -2.757338, incl. Bloechl correction: -0.009050
Generating TDOS: efermi, and dos window=   -0.0141  -0.5000   1.4859
  mmmmm m_bandcal_2nd

 mkrout:  Qtrue      sm,loc       local
   1   10.121519    2.844848    7.276671

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  pnunew: ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.478617   -0.427974    4.659654    4.659807    4.500000    4.659807
 0     0      ---       0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.413467   -0.285576    4.384698    4.385633    4.250000    4.385633
 1     0      ---       0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    9.203556   -0.231840    3.878964    3.878370    3.147584    3.878370
 2     0      ---       0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021441   -0.238755    4.110000    4.121233    4.102416    4.110000
 4     1    0.004439   -0.234481    5.100000    5.085324    5.077979    5.100000

 Harris energy:
 sumev=       -2.757338  val*vef=    -190.508101   sumtv=     187.750763
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.950964     utot=   -6534.319756    ehar=   -3304.763318

 srhov:     -6.020272   -183.674396   -189.694667 sumev=   -2.757338   sumtv=  186.937329
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.584657  avg sphere pot= 0.653080  vconst=-0.584657

 smooth rhoves     10.213952   charge     3.723329
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -2.672642   rhomu =   -3.477787  avg vxc =   -0.834273 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=F
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647

 Energy terms:             smooth           local           total
   rhoval*vef             -6.050087      -183.956396      -190.006483
   rhoval*ves            -46.758791      -116.183008      -162.941799
   psnuc*ves              67.186695    -12971.399062    -12904.212367
   utot                   10.213952     -6543.791035     -6533.577083
   rho*exc                -2.672642      -127.207021      -129.879663
   rho*vxc                -3.477787      -168.441259      -171.919046
   valence chg             3.723329         7.276671        11.000000
   core charge            18.000000        -0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Kohn-Sham energy:
 sumtv=      186.937329  sumtc=      3171.756639   ekin=     3358.693969
 rhoep=     -129.879663   utot=     -6533.577083   ehks=    -3304.762776
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 wgtsmooth=   2.7410122234342145E-002
 mixrho:  sought 3 iter from file mixm; read 2.  RMS DQ=1.64e-3  last it=7.06e-4
 AMIX: nmix=2 mmix=8  nelts=  2567  beta=1.00000  tm= 5.00000  rmsdel=1.64D-03
   tj: 0.69965  -0.00288
 mixrho: add corrections to qcell smrho = -0.27555D-07 -0.35084D-09
 unscreened rms difference:  smooth  0.000386   local  0.001200
   screened rms difference:  smooth  0.000386   local  0.001200   tot  0.001644

 iors  : write restart file (binary, mesh density) 

   it  3  of 12    ehf=      -0.328818   ehk=      -0.328276
 From last iter    ehf=      -0.328677   ehk=      -0.328571
 diffe(q)= -0.000141 (0.001644)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=1 pwmode=0 oveps=0 ehf(eV)=-4.4738332 ehk(eV)=-4.4664635 sev(eV)=-37.5157855

 --- BNDFP:  begin iteration 4 of 12
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.580840  avg sphere pot= 0.653072  vconst=-0.580840

 smooth rhoves     10.095156   charge     3.702237
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -2.655089   rhomu =   -3.454908  avg vxc =   -0.832803 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000135

 Energy terms:             smooth           local           total
   rhoval*vef             -6.000601      -184.377573      -190.378174
   rhoval*ves            -46.597890      -116.656710      -163.254600
   psnuc*ves              66.788201    -12971.812644    -12905.024443
   utot                   10.095156     -6544.234677     -6534.139521
   rho*exc                -2.655089      -127.279699      -129.934788
   rho*vxc                -3.454908      -168.537230      -171.992138
   valence chg             3.702237         7.297763        11.000000
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
 ... Done MPI k-loop: elapsed time=   0.1174

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018665;  11.000000 electrons
         Sum occ. bands:   -2.823670, incl. Bloechl correction: -0.009156
Generating TDOS: efermi, and dos window=   -0.0187  -0.5000   1.4813
  mmmmm m_bandcal_2nd

 mkrout:  Qtrue      sm,loc       local
   1   10.128469    2.831122    7.297347

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  pnunew: ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.477343   -0.429770    4.659807    4.659704    4.500000    4.659704
 0     0      ---       0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.409334   -0.288903    4.385633    4.385006    4.250000    4.385006
 1     0      ---       0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    9.216244   -0.238530    3.878370    3.878782    3.147584    3.878782
 2     0      ---       0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021154   -0.244903    4.110000    4.121067    4.102416    4.110000
 4     1    0.004395   -0.241090    5.100000    5.085253    5.077979    5.100000

 Harris energy:
 sumev=       -2.823670  val*vef=    -190.378174   sumtv=     187.554505
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.934788     utot=   -6534.139521    ehar=   -3304.763165

 srhov:     -6.001565   -184.368100   -190.369665 sumev=   -2.823670   sumtv=  187.545996
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.580908  avg sphere pot= 0.653065  vconst=-0.580908

 smooth rhoves     10.097097   charge     3.702653
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -2.655499   rhomu =   -3.455443  avg vxc =   -0.832819 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=F
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647

 Energy terms:             smooth           local           total
   rhoval*vef             -6.001920      -184.370822      -190.372742
   rhoval*ves            -46.600334      -116.649596      -163.249930
   psnuc*ves              66.794528    -12971.808066    -12905.013538
   utot                   10.097097     -6544.228831     -6534.131734
   rho*exc                -2.655499      -127.278566      -129.934066
   rho*vxc                -3.455443      -168.535736      -171.991180
   valence chg             3.702653         7.297347        11.000000
   core charge            18.000000        -0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Kohn-Sham energy:
 sumtv=      187.545996  sumtc=      3171.756639   ekin=     3359.302635
 rhoep=     -129.934066   utot=     -6534.131734   ehks=    -3304.763165
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 wgtsmooth=   2.7410122234342145E-002
 mixrho:  sought 3 iter from file mixm; read 3.  RMS DQ=1.87e-5  last it=1.64e-3
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=  2567  beta=1.00000  tm= 5.00000  rmsdel=1.87D-05
   tj:-0.01137
 mixrho: add corrections to qcell smrho = -0.24169D-07 -0.30774D-09
 unscreened rms difference:  smooth  0.000020   local  0.000015
   screened rms difference:  smooth  0.000020   local  0.000015   tot  0.000019

 iors  : write restart file (binary, mesh density) 

   it  4  of 12    ehf=      -0.328665   ehk=      -0.328665
 From last iter    ehf=      -0.328818   ehk=      -0.328276
 diffe(q)=  0.000153 (0.000019)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=1 pwmode=0 oveps=0 ehf(eV)=-4.4717536 ehk(eV)=-4.4717501 sev(eV)=-38.4182852

 --- BNDFP:  begin iteration 5 of 12
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.580867  avg sphere pot= 0.653065  vconst=-0.580867

 smooth rhoves     10.095876   charge     3.702418
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -2.655288   rhomu =   -3.455167  avg vxc =   -0.832805 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000135

 Energy terms:             smooth           local           total
   rhoval*vef             -6.001281      -184.375641      -190.376922
   rhoval*ves            -46.598717      -116.654784      -163.253501
   psnuc*ves              66.790469    -12971.812622    -12905.022153
   utot                   10.095876     -6544.233703     -6534.137827
   rho*exc                -2.655288      -127.279351      -129.934639
   rho*vxc                -3.455167      -168.536773      -171.991940
   valence chg             3.702418         7.297582        11.000000
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
 ... Done MPI k-loop: elapsed time=   0.1161

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018709;  11.000000 electrons
         Sum occ. bands:   -2.824260, incl. Bloechl correction: -0.009157
Generating TDOS: efermi, and dos window=   -0.0187  -0.5000   1.4813
  mmmmm m_bandcal_2nd

 mkrout:  Qtrue      sm,loc       local
   1   10.128519    2.830998    7.297521

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  pnunew: ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.477331   -0.429793    4.659704    4.659701    4.500000    4.659701
 0     0      ---       0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.409302   -0.288938    4.385006    4.385001    4.250000    4.385001
 1     0      ---       0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    9.216341   -0.238588    3.878782    3.878784    3.147584    3.878784
 2     0      ---       0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021151   -0.244957    4.110000    4.121066    4.102416    4.110000
 4     1    0.004394   -0.241147    5.100000    5.085252    5.077979    5.100000

 Harris energy:
 sumev=       -2.824260  val*vef=    -190.376922   sumtv=     187.552662
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.934639     utot=   -6534.137827    ehar=   -3304.763165

 srhov:     -6.001430   -184.374036   -190.375466 sumev=   -2.824260   sumtv=  187.551205
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.580879  avg sphere pot= 0.653065  vconst=-0.580879

 smooth rhoves     10.096160   charge     3.702479
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -2.655352   rhomu =   -3.455252  avg vxc =   -0.832807 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=F
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647

 Energy terms:             smooth           local           total
   rhoval*vef             -6.001482      -184.374477      -190.375960
   rhoval*ves            -46.599071      -116.653592      -163.252663
   psnuc*ves              66.791390    -12971.811706    -12905.020316
   utot                   10.096160     -6544.232649     -6534.136489
   rho*exc                -2.655352      -127.279168      -129.934520
   rho*vxc                -3.455252      -168.536530      -171.991782
   valence chg             3.702479         7.297521        11.000000
   core charge            18.000000        -0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Kohn-Sham energy:
 sumtv=      187.551205  sumtc=      3171.756639   ekin=     3359.307845
 rhoep=     -129.934520   utot=     -6534.136489   ehks=    -3304.763165
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 wgtsmooth=   2.7410122234342145E-002
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=3.14e-6  last it=1.87e-5
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=  2567  beta=1.00000  tm= 5.00000  rmsdel=3.14D-06
   tj:-0.19971
 mixrho: add corrections to qcell smrho = -0.26509D-07 -0.33752D-09
 unscreened rms difference:  smooth  0.000005   local  0.000002
   screened rms difference:  smooth  0.000005   local  0.000002   tot  0.000003

 iors  : write restart file (binary, mesh density) 

   it  5  of 12    ehf=      -0.328665   ehk=      -0.328665
 From last iter    ehf=      -0.328665   ehk=      -0.328665
 diffe(q)=  0.000000 (0.000003)    tol= 0.000010 (0.000010)   more=F
c nk=8 bigbas=1 pwmode=0 oveps=0 ehf(eV)=-4.4717483 ehk(eV)=-4.4717457 sev(eV)=-38.4263222
 >>      1.24   exit  lmfp            1.13
OK! end of LMF ======================
===START LMF   =====================================
 mpisize=           4
  bndfp (warning): no sigm file found ... LDA calculation only

mmm === MTO setting ===
mmm ispec lmxb lpzex nkapii nkaphh=    1    4    0    2    2
mmm rsmh1    1  2.50  2.50  1.00  0.00  0.00
mmm   eh1    1 -0.01 -0.01 -0.01 -0.01 -0.01
mmm rsmh2    1  1.30  0.00  1.00  1.30  0.00
mmm  eh2     1 -1.00 -1.00 -1.00 -0.01 -0.01
mmm pz       1  5.50  5.50  4.50  0.00  0.00
mmm lh       1  2  3  2

                Plat                                  Qlat
   0.000000   0.500000   0.500000       -1.000000   1.000000   1.000000
   0.500000   0.000000   0.500000        1.000000  -1.000000   1.000000
   0.500000   0.500000   0.000000        1.000000   1.000000  -1.000000
  Cell vol= 78.538660

 LATTC:  as= 2.000   tol=1.00E-12   alat= 6.79800   awald= 0.467
         r1=  2.252   nkd= 201      q1=  6.871   nkg= 331

 SGROUP: 1 symmetry operations from 0 generators
 ADDBAS: basis is already complete --- no sites added
 SYMLAT: Bravais system is cubic with 48 symmetry operations.
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-5
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
  21  r2(1,1,-0)                         
  22  m(1,1,-0)                          
  23  r2(1,-0,-1)                        
  24  m(1,-0,-1)                         
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
  47  r2(-0,1,1)                         
  48  m(-0,1,1)                          
 GROUPG: the following are sufficient to generate the space group:
 Generator(cart): i*r3(1,1,-1) r4x
 Generator(frac): i*r3(1,1,-1) r4x
 MKSYM:  found 48 space group operations ... includes inversion
 
 BZMESH:     60 irreducible QP from    8   8   8 shift=TTT
 TETIRR: sorting 3072 tetrahedra ...
 264 inequivalent ones found
 >> level: 1  CPUsec=      0.05  enter lmfp
 gen_hamindex: not readin QGpsi.

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv  kmxv foca   rfoca
 A        2.311  0.925    4    4         4  0.578  1.156    15    1   0.925

 MSHSIZ: mesh has 11 x 11 x 11 divisions; length 0.437, 0.437, 0.437
         generated from gmax = 9.0 a.u. : 941 vectors of 1331 (70%)
 goto end of reading rst or atm           1
 goto end of reading rst or atm           1
 goto end of reading rst or atm           1

 GVLST2: gmax = 9.0 a.u. created 941 vectors of 1331 (70%)
         (input) mesh has 11 x 11 x 11 divisions; length 0.437, 0.437, 0.437
 SGVSYM: 41 symmetry stars found for 941 reciprocal lattice vectors

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol= 1.00E-06
 spec      l    rsm    eh     gmax    last term    cutoff
  A        0*   2.50  -0.01   2.974    2.30E-05      27 
  A        1    2.50  -0.01   3.093    1.29E-06      51 
  A        2    1.00  -0.01   8.508    1.16E-06     869 
  A        0    1.30  -1.00   5.718    2.28E-06     259 
  A        2    1.00  -1.00   8.508    1.16E-06     869 
  A        3    1.30  -0.01   6.806    2.09E-06     411 
 goto end of reading rst or atm           1

 iors  : read restart file (binary, mesh density) 
         use from  restart file: ef window, positions, pnu 
         ignore in restart file: *
 end of reading rst or atm
 end of reading rst or atm
 m_qplistinit:start
  --- Readin syml file --- 
 end of reading rst or atm
 end of reading rst or atm
   41   0.5000   0.5000   0.5000    0.0000   0.0000   0.0000 L Gamma
   41   0.0000   0.0000   0.0000    1.0000   0.0000   0.0000 Gamma X
   21   1.0000   0.0000   0.0000    1.0000   0.5000   0.0000 X W
   41   1.0000   0.5000   0.0000    0.0000   0.0000   0.0000 W Gamma
nsyml nkp=    4  144
m_qplist_qspdivider: rank,(iqini,ispini),(iqend,ispend)=    3, (  109 1),  (  144 1)
m_qplist_qspdivider: rank,(iqini,ispini),(iqend,ispend)=    2, (   73 1),  (  108 1)
m_qplist_qspdivider: rank,(iqini,ispini),(iqend,ispend)=    1, (   37 1),  (   72 1)
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
m_qplist_qspdivider: rank,(iqini,ispini),(iqend,ispend)=    0, (    1 1),  (   36 1)

 Basis, after reading restart file
 site spec        pos (Cartesian coordinates)         pos (multiples of plat)
   1  A         0.000000   0.000000   0.000000    0.000000   0.000000   0.000000

 --- BNDFP:  begin iteration 1 of 1
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.580873  avg sphere pot= 0.653065  vconst=-0.580873

 smooth rhoves     10.095981   charge     3.702445
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -2.655321   rhomu =   -3.455211  avg vxc =   -0.832805 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000135

 Energy terms:             smooth           local           total
   rhoval*vef             -6.001387      -184.375221      -190.376608
   rhoval*ves            -46.598836      -116.654383      -163.253219
   psnuc*ves              66.790798    -12971.812428    -12905.021630
   utot                   10.095981     -6544.233406     -6534.137424
   rho*exc                -2.655321      -127.279285      -129.934607
   rho*vxc                -3.455211      -168.536686      -171.991897
   valence chg             3.702445         7.297555        11.000000
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
 ... Done MPI k-loop: elapsed time=   0.2732
  Writing bands to bands file ...
ikpoff=    2   41
ikpoff=    3   82
ikpoff=    4  103
ikpoff=    5  144
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
Exit 0 plot band mode done
 CPU time:    0.423s     Fri Oct  8 19:20:38 2021   on 

  ==== xxxxxxxxx ====     calls      == cpu time ===   depth 1
  entry   xxxx  xxxx                per call  total  (depth is by TIM= in ctrl.*.)
      0      0      0        1       0.42       0.42   main
      0      0    -10        0       0.00       0.00   `--lmfp
