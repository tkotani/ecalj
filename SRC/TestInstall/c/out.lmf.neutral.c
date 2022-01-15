=== START LFMA ===
 mpisize=           1
 HEADER sc C atom
idu uh jh=   0   0   0   0  0.00  0.00  0.00  0.00   0.00  0.00  0.00  0.00

mmm === MTO setting ===
mmm ispec lmxb lpzex nkapii nkaphh=    1    3    0    2    2
mmm rsmh1    1  1.30  1.10 -1.00 -1.00
mmm   eh1    1 -0.70 -0.20  0.00  0.00
mmm rsmh2    1  0.80  0.80  0.00  0.00
mmm  eh2     1 -1.50 -1.00  0.00  0.00
mmm pz       1  0.00  0.00  0.00  0.00
mmm lh       1  1  1
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

conf: Species C:  Z=6  Qc=2  R=3.000000  Q=0  mom=2
conf:   rmt=3.000000  rmax=19.671121  a=0.02  nr=369  nr(rmax)=463
 goto atomc xxx
 atomsc nmcore=           0

  iter     qint         drho          vh0          rho0          vsum     beta
    1    6.000000   5.461E+02       30.0000    0.2984E+02      -12.0633   0.30
   50    6.000000   3.935E-05       29.2312    0.1279E+03      -59.7470   0.30

 end of atomsc xxxxx
 vsum=  -59.746958882741843                1

 sumev=-2.876387  etot=-74.994908  eref=-74.994900  diff= -0.000008

 Optimise free-atom basis for species C, rmt=3
 l  it    Rsm      Eh     stiffR   stiffE      Eval      Exact     Pnu    Ql
 0   3   1.500  -0.941       0.0     29.9   -1.07382  -1.07383    2.91   1.00
 1   8   1.500  -0.376       0.0    134.4   -0.46556  -0.46569    2.89   2.00
 eigenvalue sum:  exact  -2.00520    opt basis  -2.00495    error 0.00025

 Optimise free-atom basis for species C, rmt=3
 l  it    Rsm      Eh     stiffR   stiffE      Eval      Exact     Pnu    Ql
 0   4   1.500  -0.769       0.0     50.2   -0.87118  -0.87119    2.91   1.00
 1   9   1.500  -0.211       0.0    359.8   -0.27748  -0.27759    2.87   0.00
 eigenvalue sum:  exact  -0.87119    opt basis  -0.87118    error 0.00001
 tailsm: init

 tailsm: fit tails to 6 smoothed hankels, rmt= 3.00000, rsm= 1.50000
    q(fit):     0.243570    rms diff:   0.000004
    fit: r>rmt  0.243570   r<rmt  1.753709   qtot  1.997279
    rho: r>rmt  0.243570   r<rmt  2.756430   qtot  3.000000

 tailsm: spin 2 ...
    q(fit):     0.054561    rms diff:   0.000002
    fit: r>rmt  0.054561   r<rmt  0.609878   qtot  0.664439
    rho: r>rmt  0.054561   r<rmt  0.945439   qtot  1.000000
 tailsm: end
conf: Core rhoc(rmt)= 0.000000 spillout= 0.000000
 end of freats: spid nmcore=C                  0
OK! end of LMFA ======================
===START LMF   ===
 mpisize=           4
 HEADER sc C atom
idu uh jh=   0   0   0   0  0.00  0.00  0.00  0.00   0.00  0.00  0.00  0.00
idu uh jh=   0   0   0   0  0.00  0.00  0.00  0.00   0.00  0.00  0.00  0.00
  bndfp (warning): no sigm file found ... LDA calculation only
idu uh jh=   0   0   0   0  0.00  0.00  0.00  0.00   0.00  0.00  0.00  0.00
idu uh jh=   0   0   0   0  0.00  0.00  0.00  0.00   0.00  0.00  0.00  0.00

mmm === MTO setting ===
mmm ispec lmxb lpzex nkapii nkaphh=    1    3    0    2    2
mmm rsmh1    1  1.30  1.10 -1.00 -1.00
mmm   eh1    1 -0.70 -0.20  0.00  0.00
mmm rsmh2    1  0.80  0.80  0.00  0.00
mmm  eh2     1 -1.50 -1.00  0.00  0.00
mmm pz       1  0.00  0.00  0.00  0.00
mmm lh       1  1  1

                Plat                                  Qlat
   1.000000   1.000000   0.000000        0.500000   0.500000  -0.500000
   1.000000   0.000000   1.000000        0.500000  -0.500000   0.500000
   0.000000   1.000000   1.000000       -0.500000   0.500000   0.500000
  Cell vol= 1000.000000

 LATTC:  as= 2.000   tol= 1.00E-08   alat= 7.93701   awald= 0.200
         r1=  3.459   nkd=  79      q1=  2.571   nkg= 137

 SGROUP: 1 symmetry operations from 0 generators
 ADDBAS: basis is already complete --- no sites added
 SYMLAT: Bravais system is cubic with 48 symmetry operations.
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-5
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
 Generator(cart): i*r3(-1,1,1) r4z
 Generator(frac): i*r3(-1,1,1) r4z
 MKSYM:  found 48 space group operations ... includes inversion
 
 BZMESH:      8 irreducible QP from    4   4   4 shift=FFF
 gen_hamindex: not readin QGpsi.

 GVLST2: gmax = 13.994 a.u. created 45911 vectors of 125000 (36%)
         (input) mesh has 50 x 50 x 50 divisions; length 0.224, 0.224, 0.224
 end of reading rst or atm
 SGVSYM: 1207 symmetry stars found for 45911 reciprocal lattice vectors

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol= 1.00E-06
 spec      l    rsm    eh     gmax    last term    cutoff
  C        0    1.30  -0.70   5.718    1.05E-06    3143 
  C        1    1.10  -0.20   7.226    1.06E-06    6375 
  C        0    0.80  -1.50   9.292    1.08E-06   13539 
  C        1    0.80  -1.00  10.038    1.00E-06   16961 
 end of reading rst or atm
 m_qplistinit:start
m_qplist_qspdivider: rank,(iqini,ispini),(iqend,ispend)=    2, (    5 1),  (    6 2)
m_qplist_qspdivider: rank,(iqini,ispini),(iqend,ispend)=    0, (    1 1),  (    2 2)
 goto end of reading rst or atm           0
 goto end of reading rst or atm           0

 iors  : read restart file (binary, mesh density) 
 iors  : empty file ... nothing read
 >> level: 1  CPUsec=      0.11  enter rdovfa

 rdovfa: read and overlap free-atom densities (mesh density) ...
 end of reading rst or atm
 end of reading rst or atm
m_qplist_qspdivider: rank,(iqini,ispini),(iqend,ispend)=    3, (    7 1),  (    8 2)
 goto end of reading rst or atm           0
m_qplist_qspdivider: rank,(iqini,ispini),(iqend,ispend)=    1, (    3 1),  (    4 2)
 goto end of reading rst or atm           0
 rdovfa: expected C,       read C        with rmt=  3.0000  mesh   369  0.020

 Free atom and overlapped crystal site charges:
   ib    true(FA)    smooth(FA)  true(OV)    smooth(OV)    local
 end of reading rst or atm
    1    3.701869    2.363587    3.701843    2.363561    1.338282
 amom    1.810990    1.143831    1.810990    1.143831    0.667159

 Smooth charge on mesh:            2.661718    moment    1.332841
 Sum of local charges:             1.338282    moments   0.667159
 Total valence charge:             4.000000    moment    2.000000
 Sum of core charges:              2.000000    moment    0.000000
 Sum of nuclear charges:          -6.000000
 Homogeneous background:           0.000000
 Deviation from neutrality:       -0.000000
 >>      0.16   exit  rdovfa          0.05
 end of reading rst or atm

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv  kmxv foca   rfoca
 C        3.000  1.200    3    3         3  0.750  1.500    15    0   1.200
 >> level: 1  CPUsec=      0.16  enter lmfp

 Basis, after reading restart file
 site spec        pos (Cartesian coordinates)         pos (multiples of plat)
   1  C         0.000000   0.000000   0.000000    0.000000   0.000000   0.000000

 --- BNDFP:  begin iteration 1 of 10

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1    0.377522    1.338282     2.00     6.00
 end of reading rst or atm
 end of reading rst or atm

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.006607  avg sphere pot= 0.019541  vconst=-0.006607
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      2.063910   charge     2.661718
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -1.494901 (  -1.109320,  -0.385581)
         rhomu =   -1.946550 (  -1.522780,  -0.423770)
       avg vxc =   -0.191348 (  -0.218128,  -0.164568)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.341410   -3.941483  -14.282893     -2.908175   -0.943828   -3.852003

 local terms:     true           smooth         local
 rhoeps:        -9.536045      -1.428951      -8.107094
 rhomu:         -7.422765      -1.450795      -5.971970
 spin2:         -5.125221      -0.410320      -4.714900
 total:        -12.547986      -1.861115     -10.686870
 val*vef       -14.282893      -6.925497      -7.357396
 val chg:        3.701843       2.363561       1.338282
 val mom:        1.810990       1.143831       0.667159    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.90000   2.90000   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.18000   3.18000   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.12000   4.12000   0.00000   0.00000
 potential shift to crystal energy zero:    0.000003

 potpus  spin 1 : pnu = 2.900000 2.850000 3.180000 4.120000

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000    -9.233051    7.589360    0.102212   -0.581578
 1      3.000000    1.000000    -5.887832    7.119808    0.144160   -0.533283
 2      3.000000    1.000000     4.727244   27.134669    0.465946   -0.095778
 3      3.000000    1.000000     7.577135   37.213692    0.543677   -0.062061

 potpus  spin 2 : pnu = 2.900000 2.850000 3.180000 4.120000

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000    -9.233051    7.437713    0.107908   -0.555888
 1      3.000000    1.000000    -5.887832    7.130054    0.153492   -0.500465
 2      3.000000    1.000000     4.727244   27.482688    0.468116   -0.093876
 3      3.000000    1.000000     7.577135   37.415024    0.544672   -0.061530

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
 ... Done MPI k-loop: elapsed time=   0.5750

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.430662;   4.000000 electrons
         Sum occ. bands:   -2.743234, incl. Bloechl correction: -0.000179
         Mag. moment:       2.000000
Generating TDOS: efermi, and dos window=   -0.4307  -0.9434   1.0693
  mmmmm m_bandcal_2nd


       contr. to mm extrapolated for r>rmt:   0.163680 est. true mm = 1.960270

 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85523  sum tc=    31.38701  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78555  sum tc=    31.54499  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.686857    3.838105   -0.151248      1.796590    2.141140   -0.344551

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  pnunew: ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.957935   -1.040174    2.900000    2.913380    2.500000    2.913380
 spn 2 0    0.945133   -0.837598    2.900000    2.907600    2.500000    2.907600
 1     1    1.783776   -0.432708    2.850000    2.886376    2.250000    2.850000
 spn 2 1    0.000000   -0.837769    2.850000    2.228737    2.250000    2.850000
 2     0    0.000011   -0.435910    3.180000    3.166283    3.147584    3.166283
 spn 2 0    0.000000   -0.837563    3.180000    3.124651    3.147584    3.147584
 3     0    0.000001   -0.441835    4.120000    4.106334    4.102416    4.106334

 spn 2 0    0.000000   -0.837445    4.120000    4.091601    4.102416    4.102416

 Harris energy:
 sumev=       -2.743234  val*vef=     -14.362499   sumtv=      11.619265
 sumec=      -39.640777  cor*vef=    -102.572782   ttcor=      62.932005
 rhoeps=      -9.601995     utot=    -139.945263    ehar=     -74.995989

 srhov:     -7.491327     -6.855971    -14.347299 sumev=   -2.743234   sumtv=   11.604065

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.042666   -0.151248     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.008341  avg sphere pot= 0.014058  vconst=-0.008341
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.531040   charge     4.151248
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -3.081542 (  -2.415659,  -0.665883)
         rhomu =   -4.025493 (  -3.336154,  -0.689339)
       avg vxc =   -0.206216 (  -0.237894,  -0.174538)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.331834   -3.946568  -14.278403     -7.320879   -1.729946   -9.050825

 local terms:     true           smooth         local
 rhoeps:        -9.525315      -3.014740      -6.510575
 rhomu:         -7.408851      -3.263040      -4.145811
 spin2:         -5.125107      -0.675925      -4.449183
 total:        -12.533959      -3.938965      -8.594994
 val*vef       -14.278403      -8.607534      -5.670869
 val chg:        3.686857       3.838105      -0.151248
 val mom:        1.796590       2.141140      -0.344551    core:   0.000000
 core chg:       2.000000       2.000000       0.000000
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
   1   -0.00   -0.00    0.00    -0.00    0.00   -0.00     0.00   -0.00    0.00
 shift forces to make zero average correction:            0.00   -0.00    0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 mixrealsmooth= T
 wgtsmooth=   2.8284271247461901E-003
 mixrho:  sought 2 iter from file mixm; read 0.  RMS DQ=8.24e-3
 AMIX: nmix=0 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 4.12D-03
 mixrho: add corrections to qcell smrho =  0.79714D-07  0.39857D-10
 unscreened rms difference:  smooth  0.012630   local  0.036279
   screened rms difference:  smooth  0.010354   local  0.036279   tot  0.008236
 mixrho: warning. negative smrho; isp number min=       2   37543 -0.55419D-05

 iors  : write restart file (binary, mesh density) 
 mixrho: warning. negative smrho; isp number min=       2   37543 -0.55419D-05
 mixrho: warning. negative smrho; isp number min=       2   37543 -0.55419D-05
 mixrho: warning. negative smrho; isp number min=       2   37543 -0.55419D-05

   it  1  of 10    ehf=      -0.001089   ehk=      -0.001245
h zbak=0 mmom=1.9999999 ehf(eV)=-.0148149 ehk(eV)=-.0169458 sev(eV)=-37.3238885

 --- BNDFP:  begin iteration 2 of 10

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1    0.167428    0.593517     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.008356  avg sphere pot= 0.016799  vconst=-0.008356
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      2.743369   charge     3.406483
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       37543  -5.5419248703921949E-006
 smvxcm: enforce positive smrho_w. Add srshift=   5.5419248803921948E-006
 smooth rhoeps =   -2.249370 (  -1.726578,  -0.522793)
         rhomu =   -2.934431 (  -2.379652,  -0.554779)
       avg vxc =   -0.201611 (  -0.230231,  -0.172991)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.336777   -3.945468  -14.282245     -4.949806   -1.314650   -6.264455

 local terms:     true           smooth         local
 rhoeps:        -9.532962      -2.180856      -7.352106
 rhomu:         -7.417505      -2.305487      -5.112018
 spin2:         -5.126437      -0.540189      -4.586248
 total:        -12.543942      -2.845676      -9.698265
 val*vef       -14.282245      -7.815003      -6.467242
 val chg:        3.697601       3.104084       0.593517
 val mom:        1.803790       1.642486       0.161304    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91338   2.90760   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16628   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10633   4.10242   0.00000   0.00000
 potential shift to crystal energy zero:    0.000003

 potpus  spin 1 : pnu = 2.913380 2.850000 3.166283 4.106334

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.750856    7.433928    0.093659   -0.587139
 1      3.000000    1.000000    -5.887832    7.121430    0.144054   -0.533608
 2      3.000000    1.000000     5.210637   28.053460    0.475959   -0.091976
 3      3.000000    1.000000     8.643874   39.422310    0.563199   -0.057686

 potpus  spin 2 : pnu = 2.907600 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.042852    7.339644    0.102864   -0.559273
 1      3.000000    1.000000    -5.887832    7.130575    0.153438   -0.500623
 2      3.000000    1.000000     6.000000   29.973164    0.494326   -0.084382
 3      3.000000    1.000000     9.000000   40.379855    0.570594   -0.055847

 Energy terms:             smooth           local           total
   rhoval*vef             -6.344798        -8.017788       -14.362585
   rhoval*ves             -4.963147        -5.280238       -10.243385
   psnuc*ves              10.449885      -280.096518      -269.646633
   utot                    2.743369      -142.688378      -139.945009
   rho*exc                -2.249370        -7.352106        -9.601476
   rho*vxc                -2.934431        -9.698265       -12.632696
   valence chg             3.406483         0.593517         4.000000
   valence mag             1.838696         0.161304         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:      0.00000
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.3513

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.431288;   4.000000 electrons
         Sum occ. bands:   -2.746491, incl. Bloechl correction: -0.000185
         Mag. moment:       2.000000
Generating TDOS: efermi, and dos window=   -0.4313  -0.9441   1.0687



       contr. to mm extrapolated for r>rmt:   0.164298 est. true mm = 1.959869

 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85641  sum tc=    31.38705  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78678  sum tc=    31.54486  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.686197    3.834521   -0.148324      1.795571    2.119216   -0.323645

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  pnunew: ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.958095   -1.040792    2.913380    2.913347    2.500000    2.913347
 spn 2 0    0.945312   -0.838791    2.907600    2.907508    2.500000    2.907508
 1     1    1.782777   -0.433430    2.850000    2.886204    2.250000    2.850000
 spn 2 1    0.000000   -0.838965    2.850000    2.228551    2.250000    2.850000
 2     0    0.000012   -0.436575    3.166283    3.166115    3.147584    3.166115
 spn 2 0    0.000000   -0.838754    3.147584    3.124605    3.147584    3.147584
 3     0    0.000001   -0.442404    4.106334    4.106278    4.102416    4.106278
 spn 2 0    0.000000   -0.838635    4.102416    4.091579    4.102416    4.102416

 Harris energy:
 sumev=       -2.746491  val*vef=     -14.362585   sumtv=      11.616094
 sumec=      -39.643186  cor*vef=    -102.575191   ttcor=      62.932005
 rhoeps=      -9.601476     utot=    -139.945009    ehar=     -74.998386

 srhov:     -8.353842     -5.999054    -14.352895 sumev=   -2.746491   sumtv=   11.606404

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.041841   -0.148324     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.008178  avg sphere pot= 0.014110  vconst=-0.008178
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.513418   charge     4.148324
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -3.076686 (  -2.402786,  -0.673900)
         rhomu =   -4.018948 (  -3.318772,  -0.700176)
       avg vxc =   -0.206362 (  -0.238091,  -0.174632)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.331797   -3.948938  -14.280736     -7.283370   -1.760793   -9.044163

 local terms:     true           smooth         local
 rhoeps:        -9.525009      -3.009751      -6.515258
 rhomu:         -7.408044      -3.245441      -4.162604
 spin2:         -5.125515      -0.686807      -4.438708
 total:        -12.533560      -3.932248      -8.601312
 val*vef       -14.280736      -8.610397      -5.670339
 val chg:        3.686197       3.834521      -0.148324
 val mom:        1.795571       2.119216      -0.323645    core:   0.000000
 core chg:       2.000000       2.000000       0.000000
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91335   2.90751   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16611   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10628   4.10242   0.00000   0.00000

 Energy terms:             smooth           local           total
   rhoval*vef             -9.123665        -5.236574       -14.360239
   rhoval*ves             -4.670952        -5.580284       -10.251236
   psnuc*ves              11.697788      -281.331558      -269.633771
   utot                    3.513418      -143.455921      -139.942503
   rho*exc                -3.076686        -6.515258        -9.591944
   rho*vxc                -4.018948        -8.601312       -12.620260
   valence chg             4.148324        -0.148324         4.000000
   valence mag             2.323645        -0.323645         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       11.606404  sumtc=        62.931909   ekin=       74.538313
 rhoep=       -9.591944   utot=      -139.942503   ehks=      -74.996134
 mag. mom=     2.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       37543  -5.5418886634155606E-006
 smvxcm: enforce positive smrho_w. Add srshift=   5.5418886734155604E-006
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       37543  -5.5419610773688293E-006
 smvxcm: enforce positive smrho_w. Add srshift=   5.5419610873688291E-006
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       37543  -5.5419248703921949E-006
 smvxcm: enforce positive smrho_w. Add srshift=   5.5419248803921948E-006

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00   -0.00    0.00    -0.00   -0.00    0.00     0.00    0.00   -0.00
 shift forces to make zero average correction:            0.00    0.00   -0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 wgtsmooth=   2.8284271247461901E-003
 mixrho:  sought 2 iter from file mixm; read 1.  RMS DQ=4.07e-3  last it=8.24e-3
 AMIX: nmix=1 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 2.03D-03
   tj:-0.97509
 mixrho: warning. negative smrho; isp number min=       2   44533 -0.94132D-05
 mixrho: warning. negative smrho; isp number min=       2   44533 -0.94132D-05
 mixrho: add corrections to qcell smrho =  0.71182D-08  0.35591D-11
 unscreened rms difference:  smooth  0.006257   local  0.017827
   screened rms difference:  smooth  0.005200   local  0.017827   tot  0.004069
 mixrho: warning. negative smrho; isp number min=       2   44533 -0.94132D-05
 mixrho: warning. negative smrho; isp number min=       2   44533 -0.94132D-05

 iors  : write restart file (binary, mesh density) 

   it  2  of 10    ehf=      -0.003486   ehk=      -0.001234
 From last iter    ehf=      -0.001089   ehk=      -0.001245
 diffe(q)= -0.002397 (0.004069)    tol= 0.000010 (0.000500)   more=T
i zbak=0 mmom=1.9999999 ehf(eV)=-.0474332 ehk(eV)=-.0167959 sev(eV)=-37.3682129

 --- BNDFP:  begin iteration 3 of 10

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.039235   -0.139084     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.009537  avg sphere pot= 0.014144  vconst=-0.009537
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.500589   charge     4.139084
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       44533  -9.4131530733723433E-006
 smvxcm: enforce positive smrho_w. Add srshift=   9.4131530833723431E-006
 smooth rhoeps =   -3.073142 (  -2.397938,  -0.675204)
         rhomu =   -4.014208 (  -3.312055,  -0.702153)
       avg vxc =   -0.210039 (  -0.240013,  -0.180065)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.332114   -3.951103  -14.283217     -7.251885   -1.756739   -9.008624

 local terms:     true           smooth         local
 rhoeps:        -9.528581      -3.002397      -6.526184
 rhomu:         -7.410743      -3.235795      -4.174949
 spin2:         -5.127466      -0.686773      -4.440693
 total:        -12.538209      -3.922567      -8.615642
 val*vef       -14.283217      -8.602686      -5.680531
 val chg:        3.691283       3.830367      -0.139084
 val mom:        1.795674       2.113278      -0.317605    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91335   2.90751   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16611   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10628   4.10242   0.00000   0.00000
 potential shift to crystal energy zero:    0.000004

 potpus  spin 1 : pnu = 2.913347 2.850000 3.166115 4.106278

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.746573    7.438818    0.093604   -0.587462
 1      3.000000    1.000000    -5.887832    7.123416    0.143913   -0.534049
 2      3.000000    1.000000     5.217017   28.057543    0.476030   -0.091971
 3      3.000000    1.000000     8.648776   39.426713    0.563255   -0.057682

 potpus  spin 2 : pnu = 2.907508 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.032189    7.343775    0.102866   -0.559473
 1      3.000000    1.000000    -5.887832    7.131448    0.153332   -0.500934
 2      3.000000    1.000000     6.000000   29.967250    0.494286   -0.084410
 3      3.000000    1.000000     9.000000   40.375635    0.570571   -0.055857

 Energy terms:             smooth           local           total
   rhoval*vef             -9.089594        -5.274594       -14.364188
   rhoval*ves             -4.673402        -5.576195       -10.249597
   psnuc*ves              11.674581      -281.317259      -269.642677
   utot                    3.500589      -143.446727      -139.946137
   rho*exc                -3.073142        -6.526184        -9.599326
   rho*vxc                -4.014208        -8.615642       -12.629850
   valence chg             4.139084        -0.139084         4.000000
   valence mag             2.317605        -0.317605         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.2884

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.432140;   4.000000 electrons
         Sum occ. bands:   -2.750740, incl. Bloechl correction: -0.000190
         Mag. moment:       2.000000
Generating TDOS: efermi, and dos window=   -0.4321  -0.9450   1.0679



       contr. to mm extrapolated for r>rmt:   0.164440 est. true mm = 1.959648

 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85808  sum tc=    31.38702  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78851  sum tc=    31.54457  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.686239    3.838551   -0.152311      1.795208    2.105299   -0.310092

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  pnunew: ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.958243   -1.041713    2.913347    2.913376    2.500000    2.913376
 spn 2 0    0.945515   -0.840287    2.907508    2.907513    2.500000    2.907513
 1     1    1.782468   -0.434345    2.850000    2.886169    2.250000    2.850000
 spn 2 1    0.000000   -0.840462    2.850000    2.228274    2.250000    2.850000
 2     0    0.000012   -0.437445    3.166115    3.165917    3.147584    3.165917
 spn 2 0    0.000000   -0.840249    3.147584    3.124528    3.147584    3.147584
 3     0    0.000001   -0.443220    4.106278    4.106210    4.102416    4.106210
 spn 2 0    0.000000   -0.840130    4.102416    4.091543    4.102416    4.102416

 Harris energy:
 sumev=       -2.750740  val*vef=     -14.364188   sumtv=      11.613447
 sumec=      -39.646594  cor*vef=    -102.578551   ttcor=      62.931957
 rhoeps=      -9.599326     utot=    -139.946137    ehar=     -75.000059

 srhov:     -9.138261     -5.225703    -14.363964 sumev=   -2.750740   sumtv=   11.613223

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.042966   -0.152311     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.008107  avg sphere pot= 0.014138  vconst=-0.008107
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.502608   charge     4.152311
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -3.081314 (  -2.398954,  -0.682360)
         rhomu =   -4.024875 (  -3.313853,  -0.711023)
       avg vxc =   -0.206394 (  -0.238142,  -0.174645)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.334363   -3.951298  -14.285660     -7.273716   -1.791902   -9.065618

 local terms:     true           smooth         local
 rhoeps:        -9.525841      -3.014416      -6.511425
 rhomu:         -7.408634      -3.240513      -4.168122
 spin2:         -5.126020      -0.697709      -4.428311
 total:        -12.534655      -3.938222      -8.596433
 val*vef       -14.285660      -8.620844      -5.664816
 val chg:        3.686239       3.838551      -0.152311
 val mom:        1.795208       2.105299      -0.310092    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91338   2.90751   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16592   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10621   4.10242   0.00000   0.00000

 Energy terms:             smooth           local           total
   rhoval*vef             -9.145068        -5.220043       -14.365110
   rhoval*ves             -4.675419        -5.579738       -10.255157
   psnuc*ves              11.680634      -281.321867      -269.641233
   utot                    3.502608      -143.450802      -139.948195
   rho*exc                -3.081314        -6.511425        -9.592739
   rho*vxc                -4.024875        -8.596433       -12.621309
   valence chg             4.152311        -0.152311         4.000000
   valence mag             2.310092        -0.310092         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       11.613223  sumtc=        62.931591   ekin=       74.544815
 rhoep=       -9.592739   utot=      -139.948195   ehks=      -74.996120
 mag. mom=     2.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       44533  -9.4131679060451310E-006
 smvxcm: enforce positive smrho_w. Add srshift=   9.4131679160451308E-006
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       44533  -9.4131382406995539E-006
 smvxcm: enforce positive smrho_w. Add srshift=   9.4131382506995537E-006
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       44533  -9.4131530733723416E-006
 smvxcm: enforce positive smrho_w. Add srshift=   9.4131530833723414E-006

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00   -0.00     0.00    0.00    0.00    -0.00   -0.00   -0.00
 shift forces to make zero average correction:           -0.00   -0.00   -0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 wgtsmooth=   2.8284271247461901E-003
 mixrho:  sought 2 iter from file mixm; read 2.  RMS DQ=9.75e-5  last it=4.07e-3
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 4.88D-05
   tj:-1.56628   0.76322
 mixrho: warning. negative smrho; isp number min=       2   32143 -0.70645D-05
 mixrho: warning. negative smrho; isp number min=       2   32143 -0.70645D-05
 mixrho: warning. negative smrho; isp number min=       2   32143 -0.70645D-05
 mixrho: add corrections to qcell smrho =  0.70228D-07  0.35114D-10
 unscreened rms difference:  smooth  0.000167   local  0.000344
   screened rms difference:  smooth  0.000159   local  0.000344   tot  0.000098
 mixrho: warning. negative smrho; isp number min=       2   32143 -0.70645D-05

 iors  : write restart file (binary, mesh density) 

   it  3  of 10    ehf=      -0.005159   ehk=      -0.001220
 From last iter    ehf=      -0.003486   ehk=      -0.001234
 diffe(q)= -0.001673 (0.000098)    tol= 0.000010 (0.000500)   more=T
i zbak=0 mmom=1.9999999 ehf(eV)=-.0701905 ehk(eV)=-.0165924 sev(eV)=-37.4260212

 --- BNDFP:  begin iteration 4 of 10

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.044675   -0.158370     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.008865  avg sphere pot= 0.014112  vconst=-0.008865
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.509471   charge     4.158370
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       32143  -7.0644852045434207E-006
 smvxcm: enforce positive smrho_w. Add srshift=   7.0644852145434206E-006
 smooth rhoeps =   -3.093103 (  -2.408078,  -0.685025)
         rhomu =   -4.040289 (  -3.326522,  -0.713767)
       avg vxc =   -0.210518 (  -0.240299,  -0.180738)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.334205   -3.952229  -14.286434     -7.296017   -1.793895   -9.089913

 local terms:     true           smooth         local
 rhoeps:        -9.527505      -3.023265      -6.504240
 rhomu:         -7.409839      -3.250905      -4.158934
 spin2:         -5.126980      -0.698917      -4.428063
 total:        -12.536819      -3.949822      -8.586997
 val*vef       -14.286434      -8.627105      -5.659329
 val chg:        3.688748       3.847118      -0.158370
 val mom:        1.795172       2.110815      -0.315643    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91338   2.90751   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16592   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10621   4.10242   0.00000   0.00000
 potential shift to crystal energy zero:    0.000005

 potpus  spin 1 : pnu = 2.913376 2.850000 3.165917 4.106210

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.750341    7.439663    0.093569   -0.587534
 1      3.000000    1.000000    -5.887832    7.124156    0.143885   -0.534124
 2      3.000000    1.000000     5.224542   28.069923    0.476171   -0.091925
 3      3.000000    1.000000     8.654789   39.437838    0.563357   -0.057662

 potpus  spin 2 : pnu = 2.907513 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.032834    7.345133    0.102840   -0.559549
 1      3.000000    1.000000    -5.887832    7.132234    0.153295   -0.501025
 2      3.000000    1.000000     6.000000   29.964365    0.494266   -0.084423
 3      3.000000    1.000000     9.000000   40.373505    0.570560   -0.055862

 Energy terms:             smooth           local           total
   rhoval*vef             -9.170448        -5.196522       -14.366970
   rhoval*ves             -4.670326        -5.583711       -10.254037
   psnuc*ves              11.689268      -281.334356      -269.645088
   utot                    3.509471      -143.459033      -139.949563
   rho*exc                -3.093103        -6.504240        -9.597343
   rho*vxc                -4.040289        -8.586997       -12.627286
   valence chg             4.158370        -0.158370         4.000000
   valence mag             2.315643        -0.315643         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:      0.00000
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.2738

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.432167;   4.000000 electrons
         Sum occ. bands:   -2.750912, incl. Bloechl correction: -0.000189
         Mag. moment:       2.000000
Generating TDOS: efermi, and dos window=   -0.4322  -0.9450   1.0678



       contr. to mm extrapolated for r>rmt:   0.164138 est. true mm = 1.959736

 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85791  sum tc=    31.38678  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78831  sum tc=    31.54436  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.686696    3.842464   -0.155769      1.795598    2.108660   -0.313062

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  pnunew: ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.958255   -1.041819    2.913376    2.913436    2.500000    2.913436
 spn 2 0    0.945548   -0.840316    2.907513    2.907610    2.500000    2.907610
 1     1    1.782880   -0.434363    2.850000    2.886277    2.250000    2.850000
 spn 2 1    0.000000   -0.840490    2.850000    2.228204    2.250000    2.850000
 2     0    0.000012   -0.437475    3.165917    3.165886    3.147584    3.165886
 spn 2 0    0.000000   -0.840278    3.147584    3.124504    3.147584    3.147584
 3     0    0.000001   -0.443287    4.106210    4.106197    4.102416    4.106197
 spn 2 0    0.000000   -0.840159    4.102416    4.091531    4.102416    4.102416

 Harris energy:
 sumev=       -2.750912  val*vef=     -14.366970   sumtv=      11.616058
 sumec=      -39.646218  cor*vef=    -102.577992   ttcor=      62.931774
 rhoeps=      -9.597343     utot=    -139.949563    ehar=     -74.999074

 srhov:     -9.162780     -5.204415    -14.367195 sumev=   -2.750912   sumtv=   11.616282

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.043942   -0.155769     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.008157  avg sphere pot= 0.014127  vconst=-0.008157
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.505496   charge     4.155769
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -3.086105 (  -2.403150,  -0.682955)
         rhomu =   -4.031168 (  -3.319686,  -0.711481)
       avg vxc =   -0.206322 (  -0.238050,  -0.174593)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.336028   -3.951401  -14.287429     -7.286224   -1.793443   -9.079667

 local terms:     true           smooth         local
 rhoeps:        -9.526586      -3.019319      -6.507267
 rhomu:         -7.409525      -3.246481      -4.163044
 spin2:         -5.126109      -0.698178      -4.427931
 total:        -12.535634      -3.944659      -8.590975
 val*vef       -14.287429      -8.624673      -5.662756
 val chg:        3.686696       3.842464      -0.155769
 val mom:        1.795598       2.108660      -0.313062    core:   0.000000
 core chg:       2.000000       2.000000       0.000000
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91344   2.90761   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16589   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10620   4.10242   0.00000   0.00000

 Energy terms:             smooth           local           total
   rhoval*vef             -9.158996        -5.207763       -14.366759
   rhoval*ves             -4.672835        -5.583148       -10.255984
   psnuc*ves              11.683828      -281.328165      -269.644337
   utot                    3.505496      -143.455657      -139.950161
   rho*exc                -3.086105        -6.507267        -9.593372
   rho*vxc                -4.031168        -8.590975       -12.622142
   valence chg             4.155769        -0.155769         4.000000
   valence mag             2.313062        -0.313062         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       11.616282  sumtc=        62.931133   ekin=       74.547415
 rhoep=       -9.593372   utot=      -139.950161   ehks=      -74.996118
 mag. mom=     2.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       32143  -7.0644851741994039E-006
 smvxcm: enforce positive smrho_w. Add srshift=   7.0644851841994037E-006
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       32143  -7.0644852348874384E-006
 smvxcm: enforce positive smrho_w. Add srshift=   7.0644852448874382E-006
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       32143  -7.0644852045434216E-006
 smvxcm: enforce positive smrho_w. Add srshift=   7.0644852145434214E-006

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00    0.00   -0.00    -0.00   -0.00    0.00     0.00    0.00   -0.00
 shift forces to make zero average correction:            0.00    0.00   -0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 wgtsmooth=   2.8284271247461901E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.59e-5  last it=9.75e-5
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 7.96D-06
   tj:-0.04258   0.00403
 mixrho: warning. negative smrho; isp number min=       2   25103 -0.54913D-05
 mixrho: warning. negative smrho; isp number min=       2   25103 -0.54913D-05
 mixrho: warning. negative smrho; isp number min=       2   25103 -0.54913D-05
 mixrho: add corrections to qcell smrho =  0.42708D-07  0.21354D-10
 unscreened rms difference:  smooth  0.000027   local  0.000069
   screened rms difference:  smooth  0.000018   local  0.000069   tot  0.000016
 mixrho: warning. negative smrho; isp number min=       2   25103 -0.54913D-05

 iors  : write restart file (binary, mesh density) 

   it  4  of 10    ehf=      -0.004174   ehk=      -0.001218
 From last iter    ehf=      -0.005159   ehk=      -0.001220
 diffe(q)=  0.000985 (0.000016)    tol= 0.000010 (0.000500)   more=T
i zbak=0 mmom=1.9999999 ehf(eV)=-.0567885 ehk(eV)=-.0165663 sev(eV)=-37.4283652

 --- BNDFP:  begin iteration 5 of 10

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.044013   -0.156023     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.008648  avg sphere pot= 0.014124  vconst=-0.008648
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.505875   charge     4.156023
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       25103  -5.4913383430082621E-006
 smvxcm: enforce positive smrho_w. Add srshift=   5.4913383530082619E-006
 smooth rhoeps =   -3.089507 (  -2.405124,  -0.684383)
         rhomu =   -4.035584 (  -3.322442,  -0.713142)
       avg vxc =   -0.209903 (  -0.239890,  -0.179916)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.335218   -3.952022  -14.287239     -7.287239   -1.793571   -9.080810

 local terms:     true           smooth         local
 rhoeps:        -9.527374      -3.020373      -6.507001
 rhomu:         -7.409939      -3.247392      -4.162547
 spin2:         -5.126717      -0.698639      -4.428078
 total:        -12.536655      -3.946030      -8.590625
 val*vef       -14.287239      -8.625075      -5.662165
 val chg:        3.688179       3.844202      -0.156023
 val mom:        1.795400       2.108834      -0.313434    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91344   2.90761   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16589   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10620   4.10242   0.00000   0.00000
 potential shift to crystal energy zero:    0.000005

 potpus  spin 1 : pnu = 2.913436 2.850000 3.165886 4.106197

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.758148    7.439204    0.093527   -0.587562
 1      3.000000    1.000000    -5.887832    7.124401    0.143881   -0.534127
 2      3.000000    1.000000     5.225725   28.071691    0.476191   -0.091918
 3      3.000000    1.000000     8.655967   39.439902    0.563376   -0.057658

 potpus  spin 2 : pnu = 2.907610 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.043911    7.344163    0.102771   -0.559595
 1      3.000000    1.000000    -5.887832    7.132505    0.153291   -0.501029
 2      3.000000    1.000000     6.000000   29.963732    0.494261   -0.084426
 3      3.000000    1.000000     9.000000   40.372987    0.570557   -0.055863

 Energy terms:             smooth           local           total
   rhoval*vef             -9.161037        -5.206430       -14.367467
   rhoval*ves             -4.672024        -5.582919       -10.254942
   psnuc*ves              11.683773      -281.329376      -269.645604
   utot                    3.505875      -143.456148      -139.950273
   rho*exc                -3.089507        -6.507001        -9.596508
   rho*vxc                -4.035584        -8.590625       -12.626209
   valence chg             4.156023        -0.156023         4.000000
   valence mag             2.313434        -0.313434         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:      0.00000
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.3579

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.432082;   4.000000 electrons
         Sum occ. bands:   -2.750550, incl. Bloechl correction: -0.000188
         Mag. moment:       2.000000
Generating TDOS: efermi, and dos window=   -0.4321  -0.9449   1.0679



       contr. to mm extrapolated for r>rmt:   0.163985 est. true mm = 1.959789

 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85773  sum tc=    31.38672  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78811  sum tc=    31.54433  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.686893    3.843313   -0.156420      1.795804    2.110854   -0.315049

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  pnunew: ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.958246   -1.041767    2.913436    2.913465    2.500000    2.913465
 spn 2 0    0.945544   -0.840192    2.907610    2.907657    2.500000    2.907657
 1     1    1.783090   -0.434270    2.850000    2.886334    2.250000    2.850000
 spn 2 1    0.000000   -0.840366    2.850000    2.228207    2.250000    2.850000
 2     0    0.000012   -0.437392    3.165886    3.165885    3.147584    3.165885
 spn 2 0    0.000000   -0.840155    3.147584    3.124501    3.147584    3.147584
 3     0    0.000001   -0.443222    4.106197    4.106195    4.102416    4.106195
 spn 2 0    0.000000   -0.840036    4.102416    4.091529    4.102416    4.102416

 Harris energy:
 sumev=       -2.750550  val*vef=     -14.367467   sumtv=      11.616917
 sumec=      -39.645847  cor*vef=    -102.577300   ttcor=      62.931454
 rhoeps=      -9.596508     utot=    -139.950273    ehar=     -74.998410

 srhov:     -9.160619     -5.206997    -14.367616 sumev=   -2.750550   sumtv=   11.617066

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.044125   -0.156420     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.008184  avg sphere pot= 0.014121  vconst=-0.008184
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.507372   charge     4.156420
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -3.087133 (  -2.404710,  -0.682423)
         rhomu =   -4.032530 (  -3.321823,  -0.710706)
       avg vxc =   -0.206286 (  -0.238004,  -0.174568)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.336604   -3.951192  -14.287796     -7.290427   -1.791279   -9.081706

 local terms:     true           smooth         local
 rhoeps:        -9.526853      -3.020393      -6.506460
 rhomu:         -7.409885      -3.248679      -4.161206
 spin2:         -5.126101      -0.697403      -4.428697
 total:        -12.535985      -3.946082      -8.589903
 val*vef       -14.287796      -8.624715      -5.663081
 val chg:        3.686893       3.843313      -0.156420
 val mom:        1.795804       2.110854      -0.315049    core:   0.000000
 core chg:       2.000000       2.000000       0.000000
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91347   2.90766   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16588   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10619   4.10242   0.00000   0.00000

 Energy terms:             smooth           local           total
   rhoval*vef             -9.160986        -5.206091       -14.367077
   rhoval*ves             -4.671466        -5.584548       -10.256015
   psnuc*ves              11.686210      -281.331471      -269.645261
   utot                    3.507372      -143.458010      -139.950638
   rho*exc                -3.087133        -6.506460        -9.593593
   rho*vxc                -4.032530        -8.589903       -12.622433
   valence chg             4.156420        -0.156420         4.000000
   valence mag             2.315049        -0.315049         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       11.617066  sumtc=        62.931047   ekin=       74.548113
 rhoep=       -9.593593   utot=      -139.950638   ehks=      -74.996118
 mag. mom=     2.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       25103  -5.4913289460763761E-006
 smvxcm: enforce positive smrho_w. Add srshift=   5.4913289560763760E-006
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       25103  -5.4913477399401480E-006
 smvxcm: enforce positive smrho_w. Add srshift=   5.4913477499401478E-006
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       25103  -5.4913383430082621E-006
 smvxcm: enforce positive smrho_w. Add srshift=   5.4913383530082619E-006

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00    0.00   -0.00     0.00    0.00   -0.00    -0.00   -0.00    0.00
 shift forces to make zero average correction:           -0.00   -0.00    0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 wgtsmooth=   2.8284271247461901E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.21e-5  last it=1.59e-5
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 6.04D-06
   tj: 0.36774   0.04566
 mixrho: warning. negative smrho; isp number min=       2   22091 -0.49233D-05
 mixrho: warning. negative smrho; isp number min=       2   22091 -0.49233D-05
 mixrho: warning. negative smrho; isp number min=       2   22091 -0.49233D-05
 mixrho: add corrections to qcell smrho =  0.42499D-07  0.21250D-10
 unscreened rms difference:  smooth  0.000022   local  0.000052
   screened rms difference:  smooth  0.000014   local  0.000052   tot  0.000012
 mixrho: warning. negative smrho; isp number min=       2   22091 -0.49233D-05

 iors  : write restart file (binary, mesh density) 

   it  5  of 10    ehf=      -0.003510   ehk=      -0.001218
 From last iter    ehf=      -0.004174   ehk=      -0.001218
 diffe(q)=  0.000664 (0.000012)    tol= 0.000010 (0.000500)   more=T
i zbak=0 mmom=1.9999999 ehf(eV)=-.0477578 ehk(eV)=-.0165678 sev(eV)=-37.4234285

 --- BNDFP:  begin iteration 6 of 10

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.044022   -0.156053     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.008603  avg sphere pot= 0.014122  vconst=-0.008603
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.506505   charge     4.156053
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       22091  -4.9232910646349373E-006
 smvxcm: enforce positive smrho_w. Add srshift=   4.9232910746349371E-006
 smooth rhoeps =   -3.089327 (  -2.405449,  -0.683878)
         rhomu =   -4.035361 (  -3.322877,  -0.712484)
       avg vxc =   -0.209572 (  -0.239693,  -0.179452)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.335519   -3.951842  -14.287361     -7.288444   -1.792202   -9.080646

 local terms:     true           smooth         local
 rhoeps:        -9.527378      -3.020457      -6.506922
 rhomu:         -7.410030      -3.248043      -4.161987
 spin2:         -5.126632      -0.698105      -4.428527
 total:        -12.536662      -3.946148      -8.590514
 val*vef       -14.287361      -8.624791      -5.662570
 val chg:        3.688090       3.844143      -0.156053
 val mom:        1.795515       2.109780      -0.314265    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91347   2.90766   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16588   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10619   4.10242   0.00000   0.00000
 potential shift to crystal energy zero:    0.000005

 potpus  spin 1 : pnu = 2.913465 2.850000 3.165885 4.106195

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.762011    7.438891    0.093507   -0.587572
 1      3.000000    1.000000    -5.887832    7.124464    0.143881   -0.534124
 2      3.000000    1.000000     5.225759   28.071655    0.476191   -0.091919
 3      3.000000    1.000000     8.656105   39.440097    0.563378   -0.057658

 potpus  spin 2 : pnu = 2.907657 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.049411    7.343568    0.102739   -0.559613
 1      3.000000    1.000000    -5.887832    7.132572    0.153291   -0.501024
 2      3.000000    1.000000     6.000000   29.963645    0.494261   -0.084427
 3      3.000000    1.000000     9.000000   40.372900    0.570557   -0.055863

 Energy terms:             smooth           local           total
   rhoval*vef             -9.160761        -5.206715       -14.367476
   rhoval*ves             -4.671609        -5.583425       -10.255034
   psnuc*ves              11.684618      -281.330122      -269.645504
   utot                    3.506505      -143.456773      -139.950269
   rho*exc                -3.089327        -6.506922        -9.596249
   rho*vxc                -4.035361        -8.590514       -12.625874
   valence chg             4.156053        -0.156053         4.000000
   valence mag             2.314264        -0.314265         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.2792

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.432042;   4.000000 electrons
         Sum occ. bands:   -2.750379, incl. Bloechl correction: -0.000188
         Mag. moment:       2.000000
Generating TDOS: efermi, and dos window=   -0.4320  -0.9449   1.0680


       contr. to mm extrapolated for r>rmt:   0.163935 est. true mm = 1.959808


 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85767  sum tc=    31.38671  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78805  sum tc=    31.54434  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.686954    3.843440   -0.156486      1.795873    2.111605   -0.315732

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  pnunew: ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.958241   -1.041738    2.913465    2.913475    2.500000    2.913475
 spn 2 0    0.945540   -0.840137    2.907657    2.907673    2.500000    2.907673
 1     1    1.783160   -0.434227    2.850000    2.886353    2.250000    2.850000
 spn 2 1    0.000000   -0.840310    2.850000    2.228213    2.250000    2.850000
 2     0    0.000012   -0.437353    3.165885    3.165886    3.147584    3.165886
 spn 2 0    0.000000   -0.840099    3.147584    3.124502    3.147584    3.147584
 3     0    0.000001   -0.443188    4.106195    4.106195    4.102416    4.106195
 spn 2 0    0.000000   -0.839981    4.102416    4.091530    4.102416    4.102416

 Harris energy:
 sumev=       -2.750379  val*vef=     -14.367476   sumtv=      11.617097
 sumec=      -39.645723  cor*vef=    -102.576974   ttcor=      62.931250
 rhoeps=      -9.596249     utot=    -139.950269    ehar=     -74.998171

 srhov:     -9.160503     -5.207117    -14.367620 sumev=   -2.750379   sumtv=   11.617241

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.044144   -0.156486     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.008193  avg sphere pot= 0.014119  vconst=-0.008193
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.508030   charge     4.156486
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -3.087278 (  -2.405127,  -0.682152)
         rhomu =   -4.032726 (  -3.322389,  -0.710337)
       avg vxc =   -0.206274 (  -0.237988,  -0.174559)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.336770   -3.951089  -14.287859     -7.291457   -1.790230   -9.081687

 local terms:     true           smooth         local
 rhoeps:        -9.526930      -3.020553      -6.506377
 rhomu:         -7.409994      -3.249264      -4.160730
 spin2:         -5.126092      -0.697033      -4.429059
 total:        -12.536086      -3.946297      -8.589789
 val*vef       -14.287859      -8.624471      -5.663387
 val chg:        3.686954       3.843440      -0.156486
 val mom:        1.795873       2.111605      -0.315732    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91347   2.90767   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16589   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10620   4.10242   0.00000   0.00000

 Energy terms:             smooth           local           total
   rhoval*vef             -9.160953        -5.206173       -14.367125
   rhoval*ves             -4.671012        -5.584971       -10.255983
   psnuc*ves              11.687072      -281.332596      -269.645524
   utot                    3.508030      -143.458783      -139.950753
   rho*exc                -3.087278        -6.506377        -9.593655
   rho*vxc                -4.032726        -8.589789       -12.622515
   valence chg             4.156486        -0.156486         4.000000
   valence mag             2.315732        -0.315732         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       11.617241  sumtc=        62.931050   ekin=       74.548291
 rhoep=       -9.593655   utot=      -139.950753   ehks=      -74.996118
 mag. mom=     2.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       22091  -4.9232782876094991E-006
 smvxcm: enforce positive smrho_w. Add srshift=   4.9232782976094990E-006
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       22091  -4.9233038416603755E-006
 smvxcm: enforce positive smrho_w. Add srshift=   4.9233038516603753E-006
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       22091  -4.9232910646349373E-006
 smvxcm: enforce positive smrho_w. Add srshift=   4.9232910746349371E-006

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00   -0.00     0.00   -0.00    0.00    -0.00   -0.00   -0.00
 shift forces to make zero average correction:           -0.00   -0.00   -0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 wgtsmooth=   2.8284271247461901E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.11e-5  last it=1.21e-5
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: Reducing nmix to  0: t_j exceeds tm: tj=-10.12736
 AMIX: nmix=0 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 5.55D-06
 mixrho: warning. negative smrho; isp number min=       2   16247 -0.38105D-05
 mixrho: warning. negative smrho; isp number min=       2   16247 -0.38105D-05
 mixrho: add corrections to qcell smrho =  0.42820D-07  0.21410D-10
 mixrho: warning. negative smrho; isp number min=       2   16247 -0.38105D-05
 unscreened rms difference:  smooth  0.000020   local  0.000048
   screened rms difference:  smooth  0.000013   local  0.000048   tot  0.000011
 mixrho: warning. negative smrho; isp number min=       2   16247 -0.38105D-05

 iors  : write restart file (binary, mesh density) 

   it  6  of 10    ehf=      -0.003271   ehk=      -0.001218
 From last iter    ehf=      -0.003510   ehk=      -0.001218
 diffe(q)=  0.000239 (0.000011)    tol= 0.000010 (0.000500)   more=T
i zbak=0 mmom=1.9999999 ehf(eV)=-.0445023 ehk(eV)=-.01657 sev(eV)=-37.4211109

 --- BNDFP:  begin iteration 7 of 10

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.044083   -0.156269     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.008491  avg sphere pot= 0.014121  vconst=-0.008491
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.507132   charge     4.156269
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       16247  -3.8105427866069666E-006
 smvxcm: enforce positive smrho_w. Add srshift=   3.8105427966069665E-006
 smooth rhoeps =   -3.089062 (  -2.405704,  -0.683358)
         rhomu =   -4.035029 (  -3.323222,  -0.711808)
       avg vxc =   -0.208997 (  -0.239352,  -0.178643)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.336159   -3.951582  -14.287741     -7.289924   -1.791314   -9.081238

 local terms:     true           smooth         local
 rhoeps:        -9.527348      -3.020697      -6.506652
 rhomu:         -7.410157      -3.248799      -4.161358
 spin2:         -5.126470      -0.697673      -4.428797
 total:        -12.536627      -3.946472      -8.590155
 val*vef       -14.287741      -8.624718      -5.663022
 val chg:        3.687811       3.844080      -0.156269
 val mom:        1.795694       2.110692      -0.314998    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91347   2.90767   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16589   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10620   4.10242   0.00000   0.00000
 potential shift to crystal energy zero:    0.000005

 potpus  spin 1 : pnu = 2.913475 2.850000 3.165886 4.106195

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.763281    7.438925    0.093500   -0.587575
 1      3.000000    1.000000    -5.887832    7.124610    0.143880   -0.534122
 2      3.000000    1.000000     5.225708   28.071268    0.476187   -0.091921
 3      3.000000    1.000000     8.656099   39.439835    0.563376   -0.057658

 potpus  spin 2 : pnu = 2.907673 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.051218    7.343512    0.102728   -0.559617
 1      3.000000    1.000000    -5.887832    7.132731    0.153291   -0.501020
 2      3.000000    1.000000     6.000000   29.963354    0.494258   -0.084428
 3      3.000000    1.000000     9.000000   40.372644    0.570555   -0.055864

 Energy terms:             smooth           local           total
   rhoval*vef             -9.161134        -5.206504       -14.367637
   rhoval*ves             -4.671178        -5.584241       -10.255419
   psnuc*ves              11.685442      -281.331407      -269.645965
   utot                    3.507132      -143.457824      -139.950692
   rho*exc                -3.089062        -6.506652        -9.595714
   rho*vxc                -4.035029        -8.590155       -12.625184
   valence chg             4.156269        -0.156269         4.000000
   valence mag             2.314998        -0.314998         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.2873

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.431974;   4.000000 electrons
         Sum occ. bands:   -2.750085, incl. Bloechl correction: -0.000188
         Mag. moment:       2.000000
Generating TDOS: efermi, and dos window=   -0.4320  -0.9448   1.0680



       contr. to mm extrapolated for r>rmt:   0.163832 est. true mm = 1.959845

 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85754  sum tc=    31.38667  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78791  sum tc=    31.54433  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.687084    3.843899   -0.156815      1.796013    2.113074   -0.317061

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  pnunew: ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.958233   -1.041691    2.913475    2.913495    2.500000    2.913495
 spn 2 0    0.945535   -0.840040    2.907673    2.907705    2.500000    2.907705
 1     1    1.783303   -0.434153    2.850000    2.886392    2.250000    2.850000
 spn 2 1    0.000000   -0.840212    2.850000    2.228222    2.250000    2.850000
 2     0    0.000012   -0.437286    3.165886    3.165887    3.147584    3.165887
 spn 2 0    0.000000   -0.840003    3.147584    3.124503    3.147584    3.147584
 3     0    0.000001   -0.443133    4.106195    4.106195    4.102416    4.106195
 spn 2 0    0.000000   -0.839885    4.102416    4.091530    4.102416    4.102416

 Harris energy:
 sumev=       -2.750085  val*vef=     -14.367637   sumtv=      11.617552
 sumec=      -39.645456  cor*vef=    -102.576606   ttcor=      62.931150
 rhoeps=      -9.595714     utot=    -139.950692    ehar=     -74.997704

 srhov:     -9.160838     -5.206938    -14.367775 sumev=   -2.750085   sumtv=   11.617690

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.044237   -0.156815     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.008211  avg sphere pot= 0.014115  vconst=-0.008211
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.509363   charge     4.156815
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -3.087817 (  -2.406083,  -0.681734)
         rhomu =   -4.033443 (  -3.323695,  -0.709748)
       avg vxc =   -0.206249 (  -0.237957,  -0.174542)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.337136   -3.950915  -14.288051     -7.293948   -1.788574   -9.082522

 local terms:     true           smooth         local
 rhoeps:        -9.527098      -3.021122      -6.505976
 rhomu:         -7.410228      -3.250610      -4.159618
 spin2:         -5.126079      -0.696443      -4.429637
 total:        -12.536307      -3.947053      -8.589255
 val*vef       -14.288051      -8.624277      -5.663774
 val chg:        3.687084       3.843899      -0.156815
 val mom:        1.796013       2.113074      -0.317061    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91349   2.90771   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16589   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10619   4.10242   0.00000   0.00000

 Energy terms:             smooth           local           total
   rhoval*vef             -9.161756        -5.205529       -14.367286
   rhoval*ves             -4.670069        -5.585894       -10.255963
   psnuc*ves              11.688795      -281.334863      -269.646068
   utot                    3.509363      -143.460379      -139.951016
   rho*exc                -3.087817        -6.505976        -9.593793
   rho*vxc                -4.033443        -8.589255       -12.622697
   valence chg             4.156815        -0.156815         4.000000
   valence mag             2.317061        -0.317061         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       11.617690  sumtc=        62.931001   ekin=       74.548691
 rhoep=       -9.593793   utot=      -139.951016   ehks=      -74.996118
 mag. mom=     2.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       16247  -3.8105233579884562E-006
 smvxcm: enforce positive smrho_w. Add srshift=   3.8105233679884560E-006
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       16247  -3.8105622152254766E-006
 smvxcm: enforce positive smrho_w. Add srshift=   3.8105622252254764E-006
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       16247  -3.8105427866069662E-006
 smvxcm: enforce positive smrho_w. Add srshift=   3.8105427966069660E-006

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00    0.00   -0.00     0.00    0.00    0.00    -0.00   -0.00   -0.00
 shift forces to make zero average correction:           -0.00   -0.00   -0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 wgtsmooth=   2.8284271247461901E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.34e-5  last it=1.11e-5
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 6.71D-06
   tj: 3.36415
 mixrho: warning. negative smrho; isp number min=       2   26111 -0.59587D-05
 mixrho: warning. negative smrho; isp number min=       2   26111 -0.59587D-05
 mixrho: warning. negative smrho; isp number min=       2   26111 -0.59587D-05
 mixrho: add corrections to qcell smrho =  0.43042D-07  0.21521D-10
 unscreened rms difference:  smooth  0.000023   local  0.000060
   screened rms difference:  smooth  0.000014   local  0.000060   tot  0.000013
 mixrho: warning. negative smrho; isp number min=       2   26111 -0.59587D-05

 iors  : write restart file (binary, mesh density) 

   it  7  of 10    ehf=      -0.002804   ehk=      -0.001218
 From last iter    ehf=      -0.003271   ehk=      -0.001218
 diffe(q)=  0.000467 (0.000013)    tol= 0.000010 (0.000500)   more=T
i zbak=0 mmom=1.9999999 ehf(eV)=-.0381454 ehk(eV)=-.0165717 sev(eV)=-37.4171115

 --- BNDFP:  begin iteration 8 of 10

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.043901   -0.155624     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.008672  avg sphere pot= 0.014128  vconst=-0.008672
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.504703   charge     4.155624
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       26111  -5.9586887088873766E-006
 smvxcm: enforce positive smrho_w. Add srshift=   5.9586887188873764E-006
 smooth rhoeps =   -3.089246 (  -2.404552,  -0.684694)
         rhomu =   -4.035233 (  -3.321662,  -0.713571)
       avg vxc =   -0.210219 (  -0.240070,  -0.180367)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.334984   -3.952192  -14.287176     -7.285208   -1.794401   -9.079609

 local terms:     true           smooth         local
 rhoeps:        -9.527328      -3.019896      -6.507432
 rhomu:         -7.409838      -3.246434      -4.163405
 spin2:         -5.126755      -0.698966      -4.427789
 total:        -12.536593      -3.945399      -8.591194
 val*vef       -14.287176      -8.625104      -5.662072
 val chg:        3.688218       3.843842      -0.155624
 val mom:        1.795317       2.107877      -0.312560    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91349   2.90771   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16589   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10619   4.10242   0.00000   0.00000
 potential shift to crystal energy zero:    0.000005

 potpus  spin 1 : pnu = 2.913495 2.850000 3.165887 4.106195

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.765884    7.438435    0.093487   -0.587588
 1      3.000000    1.000000    -5.887832    7.124357    0.143881   -0.534130
 2      3.000000    1.000000     5.225669   28.071636    0.476191   -0.091918
 3      3.000000    1.000000     8.656143   39.440331    0.563379   -0.057657

 potpus  spin 2 : pnu = 2.907705 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.054951    7.342856    0.102705   -0.559642
 1      3.000000    1.000000    -5.887832    7.132455    0.153289   -0.501035
 2      3.000000    1.000000     6.000000   29.963753    0.494262   -0.084426
 3      3.000000    1.000000     9.000000   40.373022    0.570558   -0.055863

 Energy terms:             smooth           local           total
   rhoval*vef             -9.159923        -5.207568       -14.367491
   rhoval*ves             -4.672692        -5.582227       -10.254918
   psnuc*ves              11.682097      -281.326773      -269.644676
   utot                    3.504703      -143.454500      -139.949797
   rho*exc                -3.089246        -6.507432        -9.596678
   rho*vxc                -4.035233        -8.591194       -12.626427
   valence chg             4.155624        -0.155624         4.000000
   valence mag             2.312560        -0.312560         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:      0.00000
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.2658

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.432116;   4.000000 electrons
         Sum occ. bands:   -2.750699, incl. Bloechl correction: -0.000188
         Mag. moment:       2.000000
Generating TDOS: efermi, and dos window=   -0.4321  -0.9450   1.0679
       contr. to mm extrapolated for r>rmt:   0.164024 est. true mm = 1.959774

 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85781  sum tc=    31.38675  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78820  sum tc=    31.54435  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.686847    3.843311   -0.156464      1.795750    2.110288   -0.314538


 Symmetrize density..


 Make new boundary conditions for phi,phidot..
  pnunew: ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.958251   -1.041795    2.913495    2.913458    2.500000    2.913458
 spn 2 0    0.945548   -0.840239    2.907705    2.907645    2.500000    2.907645
 1     1    1.783036   -0.434307    2.850000    2.886319    2.250000    2.850000
 spn 2 1    0.000000   -0.840413    2.850000    2.228199    2.250000    2.850000
 2     0    0.000012   -0.437426    3.165887    3.165883    3.147584    3.165883
 spn 2 0    0.000000   -0.840202    3.147584    3.124500    3.147584    3.147584
 3     0    0.000001   -0.443248    4.106195    4.106195    4.102416    4.106195
 spn 2 0    0.000000   -0.840083    4.102416    4.091529    4.102416    4.102416

 Harris energy:
 sumev=       -2.750699  val*vef=     -14.367491   sumtv=      11.616792
 sumec=      -39.646009  cor*vef=    -102.577084   ttcor=      62.931076
 rhoeps=      -9.596678     utot=    -139.949797    ehar=     -74.998608


 srhov:     -9.160603     -5.207076    -14.367679 sumev=   -2.750699   sumtv=   11.616979

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.044138   -0.156464     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.008177  avg sphere pot= 0.014123  vconst=-0.008177
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.506863   charge     4.156464
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -3.087143 (  -2.404465,  -0.682678)
         rhomu =   -4.032538 (  -3.321495,  -0.711043)
       avg vxc =   -0.206295 (  -0.238016,  -0.174574)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.336490   -3.951291  -14.287781     -7.289904   -1.792243   -9.082147

 local terms:     true           smooth         local
 rhoeps:        -9.526804      -3.020393      -6.506411
 rhomu:         -7.409807      -3.248336      -4.161471
 spin2:         -5.126113      -0.697741      -4.428373
 total:        -12.535920      -3.946076      -8.589844
 val*vef       -14.287781      -8.625055      -5.662725
 val chg:        3.686847       3.843311      -0.156464
 val mom:        1.795750       2.110288      -0.314538    core:   0.000000
 core chg:       2.000000       2.000000       0.000000
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91346   2.90765   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16588   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10619   4.10242   0.00000   0.00000

 Energy terms:             smooth           local           total
   rhoval*vef             -9.161438        -5.205635       -14.367072
   rhoval*ves             -4.671809        -5.584255       -10.256065
   psnuc*ves              11.685534      -281.330761      -269.645226
   utot                    3.506863      -143.457508      -139.950646
   rho*exc                -3.087143        -6.506411        -9.593554
   rho*vxc                -4.032538        -8.589844       -12.622381
   valence chg             4.156464        -0.156464         4.000000
   valence mag             2.314538        -0.314538         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       11.616979  sumtc=        62.931103   ekin=       74.548082
 rhoep=       -9.593554   utot=      -139.950646   ehks=      -74.996118
 mag. mom=     2.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       26111  -5.9586821111232740E-006
 smvxcm: enforce positive smrho_w. Add srshift=   5.9586821211232738E-006
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       26111  -5.9586953066514801E-006
 smvxcm: enforce positive smrho_w. Add srshift=   5.9586953166514799E-006
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=       26111  -5.9586887088873775E-006
 smvxcm: enforce positive smrho_w. Add srshift=   5.9586887188873773E-006

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00   -0.00    -0.00   -0.00   -0.00    -0.00   -0.00    0.00
 shift forces to make zero average correction:           -0.00   -0.00    0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 wgtsmooth=   2.8284271247461901E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.43e-5  last it=1.34e-5
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 7.13D-06
   tj: 1.74562
 mixrho: warning. negative smrho; isp number min=       2    5999 -0.16213D-05
 mixrho: warning. negative smrho; isp number min=       2    5999 -0.16213D-05
 mixrho: warning. negative smrho; isp number min=       2    5999 -0.16213D-05
 mixrho: add corrections to qcell smrho =  0.43120D-07  0.21560D-10
 unscreened rms difference:  smooth  0.000025   local  0.000064
   screened rms difference:  smooth  0.000016   local  0.000064   tot  0.000014
 mixrho: warning. negative smrho; isp number min=       2    5999 -0.16213D-05

 iors  : write restart file (binary, mesh density) 

   it  8  of 10    ehf=      -0.003708   ehk=      -0.001218
 From last iter    ehf=      -0.002804   ehk=      -0.001218
 diffe(q)= -0.000904 (0.000014)    tol= 0.000010 (0.000500)   more=T
i zbak=0 mmom=1.9999999 ehf(eV)=-.0504439 ehk(eV)=-.0165651 sev(eV)=-37.4254661

 --- BNDFP:  begin iteration 9 of 10

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.044264   -0.156914     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.008323  avg sphere pot= 0.014112  vconst=-0.008323
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.510055   charge     4.156914
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        5999  -1.6212690244183436E-006
 smvxcm: enforce positive smrho_w. Add srshift=   1.6212690344183436E-006
 smooth rhoeps =   -3.088891 (  -2.407074,  -0.681817)
         rhomu =   -4.034846 (  -3.325069,  -0.709777)
       avg vxc =   -0.207598 (  -0.238575,  -0.176620)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.337329   -3.950908  -14.288237     -7.295197   -1.787448   -9.082646

 local terms:     true           smooth         local
 rhoeps:        -9.527390      -3.021528      -6.505862
 rhomu:         -7.410506      -3.251471      -4.159035
 spin2:         -5.126182      -0.696118      -4.430064
 total:        -12.536688      -3.947589      -8.589099
 val*vef       -14.288237      -8.624085      -5.664152
 val chg:        3.687458       3.844372      -0.156914
 val mom:        1.796092       2.113971      -0.317879    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91346   2.90765   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16588   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10619   4.10242   0.00000   0.00000
 potential shift to crystal energy zero:    0.000005

 potpus  spin 1 : pnu = 2.913458 2.850000 3.165883 4.106195

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.761028    7.439359    0.093511   -0.587561
 1      3.000000    1.000000    -5.887832    7.124850    0.143880   -0.534112
 2      3.000000    1.000000     5.225824   28.071091    0.476186   -0.091922
 3      3.000000    1.000000     8.656118   39.439502    0.563374   -0.057659

 potpus  spin 2 : pnu = 2.907645 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.047993    7.344075    0.102749   -0.559591
 1      3.000000    1.000000    -5.887832    7.132995    0.153293   -0.501002
 2      3.000000    1.000000     6.000000   29.963023    0.494255   -0.084430
 3      3.000000    1.000000     9.000000   40.372313    0.570553   -0.055865

 Energy terms:             smooth           local           total
   rhoval*vef             -9.162110        -5.205592       -14.367703
   rhoval*ves             -4.669377        -5.586439       -10.255816
   psnuc*ves              11.689487      -281.336447      -269.646960
   utot                    3.510055      -143.461443      -139.951388
   rho*exc                -3.088891        -6.505862        -9.594753
   rho*vxc                -4.034846        -8.589099       -12.623945
   valence chg             4.156914        -0.156914         4.000000
   valence mag             2.317879        -0.317879         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.3209

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.431825;   4.000000 electrons
         Sum occ. bands:   -2.749440, incl. Bloechl correction: -0.000187
         Mag. moment:       2.000000
Generating TDOS: efermi, and dos window=   -0.4318  -0.9447   1.0682

       contr. to mm extrapolated for r>rmt:   0.163639 est. true mm = 1.959918

 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 

 sum q= 1.00  sum ec=   -19.85728  sum tc=    31.38661  rho(rmax) 0.00000

 sum q= 1.00  sum ec=   -19.78763  sum tc=    31.54431  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.687325    3.844630   -0.157304      1.796279    2.115835   -0.319557

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  pnunew: ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.958215   -1.041579    2.913458    2.913532    2.500000    2.913532
 spn 2 0    0.945523   -0.839831    2.907645    2.907767    2.500000    2.907767
 1     1    1.783575   -0.433990    2.850000    2.886465    2.250000    2.850000
 spn 2 1    0.000000   -0.840002    2.850000    2.228248    2.250000    2.850000
 2     0    0.000012   -0.437138    3.165883    3.165892    3.147584    3.165892
 spn 2 0    0.000000   -0.839794    3.147584    3.124506    3.147584    3.147584
 3     0    0.000001   -0.443013    4.106195    4.106195    4.102416    4.106195
 spn 2 0    0.000000   -0.839677    4.102416    4.091531    4.102416    4.102416

 Harris energy:
 sumev=       -2.749440  val*vef=     -14.367703   sumtv=      11.618262
 sumec=      -39.644909  cor*vef=    -102.575998   ttcor=      62.931089
 rhoeps=      -9.594753     utot=    -139.951388    ehar=     -74.996790

 srhov:     -9.161679     -5.206161    -14.367840 sumev=   -2.749440   sumtv=   11.618399

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.044375   -0.157304     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.008247  avg sphere pot= 0.014106  vconst=-0.008247
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.512043   charge     4.157304
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -3.088635 (  -2.407762,  -0.680873)
         rhomu =   -4.034534 (  -3.325982,  -0.708552)
       avg vxc =   -0.206203 (  -0.237896,  -0.174510)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.337782   -3.950534  -14.288316     -7.298196   -1.785205   -9.083401

 local terms:     true           smooth         local
 rhoeps:        -9.527398      -3.021995      -6.505402
 rhomu:         -7.410654      -3.252972      -4.157683
 spin2:         -5.126046      -0.695245      -4.430801
 total:        -12.536701      -3.948216      -8.588484
 val*vef       -14.288316      -8.623587      -5.664729
 val chg:        3.687325       3.844630      -0.157304
 val mom:        1.796279       2.115835      -0.319557    core:   0.000000
 core chg:       2.000000       2.000000       0.000000
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91353   2.90777   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16589   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10619   4.10242   0.00000   0.00000

 Energy terms:             smooth           local           total
   rhoval*vef             -9.162577        -5.204916       -14.367493
   rhoval*ves             -4.668230        -5.587621       -10.255852
   psnuc*ves              11.692316      -281.339275      -269.646959
   utot                    3.512043      -143.463448      -139.951405
   rho*exc                -3.088635        -6.505402        -9.594037
   rho*vxc                -4.034534        -8.588484       -12.623018
   valence chg             4.157304        -0.157304         4.000000
   valence mag             2.319557        -0.319557         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       11.618399  sumtc=        62.930925   ekin=       74.549324
 rhoep=       -9.594037   utot=      -139.951405   ehks=      -74.996118
 mag. mom=     2.000000
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        5999  -1.6212365242470748E-006
 smvxcm: enforce positive smrho_w. Add srshift=   1.6212365342470749E-006
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        5999  -1.6213015245896123E-006
 smvxcm: enforce positive smrho_w. Add srshift=   1.6213015345896123E-006
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=        5999  -1.6212690244183436E-006
 smvxcm: enforce positive smrho_w. Add srshift=   1.6212690344183436E-006

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00    0.00   -0.00     0.00    0.00    0.00     0.00   -0.00   -0.00
 shift forces to make zero average correction:            0.00   -0.00   -0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 wgtsmooth=   2.8284271247461901E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.06e-5  last it=1.43e-5
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 5.32D-06
   tj:-1.00299
 mixrho: add corrections to qcell smrho =  0.43829D-07  0.21915D-10
 unscreened rms difference:  smooth  0.000018   local  0.000047
   screened rms difference:  smooth  0.000012   local  0.000047   tot  0.000011

 iors  : write restart file (binary, mesh density) 

   it  9  of 10    ehf=      -0.001890   ehk=      -0.001218
 From last iter    ehf=      -0.003708   ehk=      -0.001218
 diffe(q)=  0.001818 (0.000011)    tol= 0.000010 (0.000500)   more=T
i zbak=0 mmom=1.9999999 ehf(eV)=-.0257112 ehk(eV)=-.0165766 sev(eV)=-37.4083338

 --- BNDFP:  begin iteration 10 of 10

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.044621   -0.158177     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.008068  avg sphere pot= 0.014093  vconst=-0.008068
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.516467   charge     4.158177
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -3.089732 (  -2.410558,  -0.679174)
         rhomu =   -4.036016 (  -3.329807,  -0.706209)
       avg vxc =   -0.207653 (  -0.238394,  -0.176911)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.339365   -3.949587  -14.288951     -7.305894   -1.779216   -9.085111

 local terms:     true           smooth         local
 rhoeps:        -9.527548      -3.023213      -6.504336
 rhomu:         -7.411212      -3.256943      -4.154268
 spin2:         -5.125694      -0.692909      -4.432786
 total:        -12.536906      -3.949852      -8.587054
 val*vef       -14.288951      -8.622506      -5.666445
 val chg:        3.686989       3.845166      -0.158177
 val mom:        1.796839       2.120741      -0.323903    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91353   2.90777   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16589   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10619   4.10242   0.00000   0.00000
 potential shift to crystal energy zero:    0.000005

 potpus  spin 1 : pnu = 2.913532 2.850000 3.165892 4.106195

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.770774    7.438745    0.093462   -0.587575
 1      3.000000    1.000000    -5.887832    7.125243    0.143882   -0.534089
 2      3.000000    1.000000     5.225476   28.069834    0.476173   -0.091928
 3      3.000000    1.000000     8.656140   39.438976    0.563371   -0.057661

 potpus  spin 2 : pnu = 2.907767 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.062063    7.342664    0.102670   -0.559613
 1      3.000000    1.000000    -5.887832    7.133430    0.153301   -0.500960
 2      3.000000    1.000000     6.000000   29.962633    0.494251   -0.084432
 3      3.000000    1.000000     9.000000   40.371862    0.570550   -0.055866

 Energy terms:             smooth           local           total
   rhoval*vef             -9.164078        -5.203842       -14.367920
   rhoval*ves             -4.665458        -5.590797       -10.256255
   psnuc*ves              11.698391      -281.347090      -269.648699
   utot                    3.516467      -143.468944      -139.952477
   rho*exc                -3.089732        -6.504336        -9.594068
   rho*vxc                -4.036016        -8.587054       -12.623070
   valence chg             4.158177        -0.158177         4.000000
   valence mag             2.323903        -0.323903         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.3438

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.431607;   4.000000 electrons
         Sum occ. bands:   -2.748498, incl. Bloechl correction: -0.000186
         Mag. moment:       2.000000
Generating TDOS: efermi, and dos window=   -0.4316  -0.9444   1.0684
       contr. to mm extrapolated for r>rmt:   0.163490 est. true mm = 1.959962

 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 

 sum q= 1.00  sum ec=   -19.85680  sum tc=    31.38651  rho(rmax) 0.00000

 sum q= 1.00  sum ec=   -19.78710  sum tc=    31.54430  rho(rmax) 0.00000

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.687484    3.845631   -0.158147      1.796472    2.118629   -0.322156

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  pnunew: ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.958206   -1.041400    2.913532    2.913565    2.500000    2.913565
 spn 2 0    0.945505   -0.839515    2.907767    2.907814    2.500000    2.907814
 1     1    1.783760   -0.433767    2.850000    2.886526    2.250000    2.850000
 spn 2 1    0.000000   -0.839686    2.850000    2.228276    2.250000    2.850000
 2     0    0.000012   -0.436924    3.165892    3.165896    3.147584    3.165896
 spn 2 0    0.000000   -0.839479    3.147584    3.124510    3.147584    3.147584
 3     0    0.000001   -0.442812    4.106195    4.106195    4.102416    4.106195
 spn 2 0    0.000000   -0.839361    4.102416    4.091532    4.102416    4.102416

 Harris energy:
 sumev=       -2.748498  val*vef=     -14.367920   sumtv=      11.619422
 sumec=      -39.643908  cor*vef=    -102.574915   ttcor=      62.931007
 rhoeps=      -9.594068     utot=    -139.952477    ehar=     -74.996117


 srhov:     -9.166285     -5.200607    -14.366893 sumev=   -2.748498   sumtv=   11.618394

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.044612   -0.158147     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode
 smves:: avg es pot at rmt= 0.008273  avg sphere pot= 0.014100  vconst=-0.008273
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.514340   charge     4.158147
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -3.089903 (  -2.409766,  -0.680136)
         rhomu =   -4.036217 (  -3.328722,  -0.707495)
       avg vxc =   -0.206175 (  -0.237858,  -0.174493)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.338053   -3.950070  -14.288123     -7.304047   -1.782199   -9.086246

 local terms:     true           smooth         local
 rhoeps:        -9.527543      -3.023301      -6.504242
 rhomu:         -7.410903      -3.255764      -4.155139
 spin2:         -5.125988      -0.694184      -4.431804
 total:        -12.536891      -3.949947      -8.586944
 val*vef       -14.288123      -8.623846      -5.664277
 val chg:        3.687484       3.845631      -0.158147
 val mom:        1.796472       2.118629      -0.322156    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
  ibas l=  1  0 pnu(1:nsp) pnz(1:nsp)=   2.91356   2.90781   0.00000   0.00000
  ibas l=  1  1 pnu(1:nsp) pnz(1:nsp)=   2.85000   2.85000   0.00000   0.00000
  ibas l=  1  2 pnu(1:nsp) pnz(1:nsp)=   3.16590   3.14758   0.00000   0.00000
  ibas l=  1  3 pnu(1:nsp) pnz(1:nsp)=   4.10619   4.10242   0.00000   0.00000

 Energy terms:             smooth           local           total
   rhoval*vef             -9.165384        -5.201878       -14.367262
   rhoval*ves             -4.666768        -5.588697       -10.255464
   psnuc*ves              11.695447      -281.342339      -269.646891
   utot                    3.514340      -143.465518      -139.951178
   rho*exc                -3.089903        -6.504242        -9.594145
   rho*vxc                -4.036217        -8.586944       -12.623161
   valence chg             4.158147        -0.158147         4.000000
   valence mag             2.322156        -0.322156         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       11.618394  sumtc=        62.930810   ekin=       74.549204
 rhoep=       -9.594145   utot=      -139.951178   ehks=      -74.996119
 mag. mom=     2.000000
 smvxcm: all smrho_w is positive
 smvxcm: all smrho_w is positive
 smvxcm: all smrho_w is positive

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00   -0.00   -0.00    -0.00   -0.00    0.00     0.00    0.00   -0.00
 shift forces to make zero average correction:            0.00    0.00   -0.00

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 wgtsmooth=   2.8284271247461901E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.51e-5  last it=1.06e-5
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel= 7.56D-06
   tj: 0.58922
 mixrho: add corrections to qcell smrho =  0.43000D-07  0.21500D-10
 unscreened rms difference:  smooth  0.000027   local  0.000061
   screened rms difference:  smooth  0.000019   local  0.000061   tot  0.000015

 iors  : write restart file (binary, mesh density) 

   it 10  of 10    ehf=      -0.001217   ehk=      -0.001219
 From last iter    ehf=      -0.001890   ehk=      -0.001218
 diffe(q)=  0.000673 (0.000015)    tol= 0.000010 (0.000500)   more=F
x zbak=0 mmom=1.9999999 ehf(eV)=-.0165534 ehk(eV)=-.0165812 sev(eV)=-37.3955195
 >>     25.09   exit  lmfp           24.92
OK! end of LMF ======================
