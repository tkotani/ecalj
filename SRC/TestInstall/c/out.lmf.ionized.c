===START LMFA   =====================================
 mpisize=           1
 HEADER sc C atom
 nkaphh(j) = nkapii(j) + lpzex(j)           2           2           0
  mpipid=           0
mmm === MTO setting ===
mmm ispec lmxb lpzex nkapii nkaphh=    1    3    0    2    2
mmm rsmh1     1  1.30  1.10 -1.00 -1.00
mmm   eh1     1 -0.70 -0.20  0.00  0.00
mmm rsmh2    1  0.80  0.80  0.00  0.00
mmm  eh2     1 -1.50 -1.00  0.00  0.00
mmm pz       1  0.00  0.00  0.00  0.00
mmm lh       1  1  1

                Plat                                  Qlat
   1.000000   1.000000   0.000000        0.500000   0.500000  -0.500000
   1.000000   0.000000   1.000000        0.500000  -0.500000   0.500000
   0.000000   1.000000   1.000000       -0.500000   0.500000   0.500000
  Cell vol= 1000.000000

 LATTC:  as= 2.000   tol=1.00E-08   alat= 7.93701   awald= 0.200
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
===START LMF   =====================================
 mpisize=           4
 HEADER sc C atom
 nkaphh(j) = nkapii(j) + lpzex(j)           2           2           0
 nkaphh(j) = nkapii(j) + lpzex(j)           2           2           0
 nkaphh(j) = nkapii(j) + lpzex(j)           2           2           0
 nkaphh(j) = nkapii(j) + lpzex(j)           2           2           0
  bndfp (warning): no sigm file found ... LDA calculation only
  mpipid=           2
  mpipid=           1
  mpipid=           3
  mpipid=           0
mmm === MTO setting ===
mmm ispec lmxb lpzex nkapii nkaphh=    1    3    0    2    2
mmm rsmh1     1  1.30  1.10 -1.00 -1.00
mmm   eh1     1 -0.70 -0.20  0.00  0.00
mmm rsmh2    1  0.80  0.80  0.00  0.00
mmm  eh2     1 -1.50 -1.00  0.00  0.00
mmm pz       1  0.00  0.00  0.00  0.00
mmm lh       1  1  1

                Plat                                  Qlat
   1.000000   1.000000   0.000000        0.500000   0.500000  -0.500000
   1.000000   0.000000   1.000000        0.500000  -0.500000   0.500000
   0.000000   1.000000   1.000000       -0.500000   0.500000   0.500000
  Cell vol= 1000.000000

 LATTC:  as= 2.000   tol=1.00E-08   alat= 7.93701   awald= 0.200
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
 
 BZMESH:  8 irreducible QP from 64 ( 4 4 4 )  shift= F F F
 >> level: 1  CPUsec=      0.04  enter lmfp

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv  kmxv foca   rfoca
 C        3.000  1.200    3    3         3  0.750  1.500    15    0   1.200

 GVLST2: gmax = 13.994 a.u. created 45911 vectors of 125000 (36%)
         (input) mesh has 50 x 50 x 50 divisions; length 0.224, 0.224, 0.224
 SGVSYM: 1207 symmetry stars found for 45911 reciprocal lattice vectors

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol= 1.00E-06
 spec      l    rsm    eh     gmax    last term    cutoff
  C        0    1.30  -0.70   5.718    1.05E-06    3143 
  C        1    1.10  -0.20   7.226    1.06E-06    6375 
  C        0    0.80  -1.50   9.292    1.08E-06   13539 
  C        1    0.80  -1.00  10.038    1.00E-06   16961 
 gen_hamindex: not readin QGpsi.

 iors  : read restart file (binary, mesh density) 
 iors  : empty file ... nothing read

 rdovfa: read and overlap free-atom densities (mesh density) ...
 rdovfa: expected C,       read C        with rmt=  3.0000  mesh   369  0.020
 Uniform density added to neutralize background, q=1.000000

 Smooth charge on mesh:            1.661718    moment    1.332841
 Sum of local charges:             1.338282    moments   0.667159
 Total valence charge:             3.000000    moment    2.000000
 Sum of core charges:              2.000000    moment    0.000000
 Sum of nuclear charges:          -6.000000
 Homogeneous background:           1.000000
 Deviation from neutrality:       -0.000000
 -------- qplist --------           0

 Basis, after reading restart file
 site spec        pos (Cartesian coordinates)         pos (multiples of plat)
   1  C         0.000000   0.000000   0.000000    0.000000   0.000000   0.000000
 --- BNDFP:  begin iteration 1 of 10
bndfp:start

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1    0.377522    1.338282     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 esmsmves
   50   50   50     45911     45911
   7.93701   7.93701   0.00000
   7.93701   0.00000   7.93701
   0.00000   7.93701   7.93701
effective screening medium method jesm=  0
 jtresm,tresm=  0     0.000
 z0esm,z1esm,z2esm=     3.969     0.000     0.000 sa0=     63.00

 jesm=0 return

 smves:: avg es pot at rmt= 0.006607  avg sphere pot= 0.019541  vconst=-0.006607
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      2.063910   charge     2.661718
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      201254  -4.9989572695401780E-004
 smvxcm: enforce positive smrho_w. Add srshift=   4.9989572696401783E-004
 smooth rhoeps =   -1.494861 (  -1.109298,  -0.385563)
         rhomu =   -1.946499 (  -1.522749,  -0.423750)
       avg vxc =   -0.191138 (  -0.218032,  -0.164244)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       128
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       128
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       128
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       128
 vxcnsp (warning): negative rho: min val =  -4.61E-04
 qqx nkapi nkape=           2           2
 qqx nkapi nkape=           2           2
 qqx nkapi nkape=           2           2

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.261804   -3.876233  -14.138038     -2.842621   -0.891813   -3.734435

 local terms:     true           smooth         local
 rhoeps:        -9.464952      -1.364263      -8.100688
 rhomu:         -7.369404      -1.402026      -5.967378
 spin2:         -5.086143      -0.375091      -4.711052
 total:        -12.455547      -1.777117     -10.678430
 val*vef       -14.138038      -6.807929      -7.330109
 val chg:        3.588746       2.250464       1.338282
 val mom:        1.810990       1.143831       0.667159    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 qqx nkapi nkape=           2           2
 potential shift to crystal energy zero:    0.000003

 potpus  spin 1 : pnu = 2.900000 2.850000 3.180000 4.120000

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000    -9.233051    7.622132    0.101666   -0.583561
 1      3.000000    1.000000    -5.887832    7.139046    0.143256   -0.535857
 2      3.000000    1.000000     4.727244   27.072461    0.465408   -0.096156
 3      3.000000    1.000000     7.577135   37.163438    0.543349   -0.062204

 potpus  spin 2 : pnu = 2.900000 2.850000 3.180000 4.120000

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000    -9.233051    7.498218    0.106861   -0.559305
 1      3.000000    1.000000    -5.887832    7.161467    0.151762   -0.504954
 2      3.000000    1.000000     4.727244   27.365019    0.467061   -0.094577
 3      3.000000    1.000000     7.577135   37.316589    0.544000   -0.061810

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
 ... Start MPI k-loop ...
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   0.6251

 BZWTS : --- Tetrahedron Integration ---
 ... only filled or empty bands encountered:  ev=-0.824805  ec=-0.772717
 VBmax = -0.824805  CBmin = -0.772717  gap = 0.052087 Ry = 0.70869 eV
 BZINTS: Fermi energy:     -0.824805;   3.000000 electrons
         Sum occ. bands:   -3.277726, incl. Bloechl correction:  0.000000
         Mag. moment:      -1.000000
Generating TDOS: efermi, and dos window=   -0.8248  -0.5000   0.6752
 bndfp: kpt    1 of    8 k isp=  0.0000  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    1 of    8 k isp=  0.0000  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    2 of    8 k isp=  0.1250  0.1250 -0.1250 1 ndimh nev=    8    5
 bndfp: kpt    2 of    8 k isp=  0.1250  0.1250 -0.1250 2 ndimh nev=    8    5
 bndfp: kpt    3 of    8 k isp=  0.2500  0.2500 -0.2500 1 ndimh nev=    8    5
 bndfp: kpt    3 of    8 k isp=  0.2500  0.2500 -0.2500 2 ndimh nev=    8    5
 bndfp: kpt    4 of    8 k isp=  0.2500  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    4 of    8 k isp=  0.2500  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    5 of    8 k isp=  0.3750  0.1250 -0.1250 1 ndimh nev=    8    5
 bndfp: kpt    5 of    8 k isp=  0.3750  0.1250 -0.1250 2 ndimh nev=    8    5
 bndfp: kpt    6 of    8 k isp=  0.5000  0.2500 -0.2500 1 ndimh nev=    8    5
 bndfp: kpt    6 of    8 k isp=  0.5000  0.2500 -0.2500 2 ndimh nev=    8    5
 bndfp: kpt    7 of    8 k isp=  0.5000  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    7 of    8 k isp=  0.5000  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    8 of    8 k isp=  0.5000  0.2500  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    8 of    8 k isp=  0.5000  0.2500  0.0000 2 ndimh nev=    8    5
  mmmmm m_bandcal_2nd
 mkrout: site(class) decomposed charge and magnetic moment. class->lmchk

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.130156  386.510557 -384.380401     -0.144414 -327.846611  327.702197
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85491  sum tc=    31.38764  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78517  sum tc=    31.54593  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.992767   -1.058527    2.900000    2.891276    2.500000    2.891276
 spn 2 0    1.134455   -0.954472    2.900000    2.741424    2.500000    2.741424
 1     1    0.000104   -1.058680    2.850000    2.201384    2.250000    2.850000
 spn 2 1    0.002815   -1.380732    2.850000    2.133553    2.250000    2.850000
 2     0    0.000000   -1.058192    3.180000    3.117620    3.147584    3.147584
 spn 2 0    0.000012   -1.345455    3.180000    3.101060    3.147584    3.147584
 3     0    0.000000   -1.058779    4.120000    4.088321    4.102416    4.102416
 spn 2 0    0.000004   -1.377084    4.120000    4.079830    4.102416    4.102416

 Harris energy:
 sumev=       -3.277726  val*vef=     -14.137456   sumtv=      10.859730
 sumec=      -39.640082  cor*vef=    -102.572087   ttcor=      62.932005
 rhoeps=      -9.595549     utot=    -139.945263    ehar=     -75.749078
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      201254  -4.9989874523704528E-004
 smvxcm: enforce positive smrho_w. Add srshift=   4.9989874524704531E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      201254  -4.9989270867099021E-004
 smvxcm: enforce positive smrho_w. Add srshift=   4.9989270868099024E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      201254  -4.9989572695401769E-004
 smvxcm: enforce positive smrho_w. Add srshift=   4.9989572696401772E-004

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00   -0.00    0.00    -0.00    0.00    0.00     0.00   -0.00   -0.00
 shift forces to make zero average correction:            0.00   -0.00   -0.00

 srhov:  -1237.408090   1225.635202    -11.772888 sumev=   -3.277726   sumtv=    8.495162

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1 -108.431709 -384.380401     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 jesm=0 return

 smves:: avg es pot at rmt=-0.205342  avg sphere pot= 0.030809  vconst= 0.205342
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves     38.909181   charge   388.380401
 smvxcm: all smrho_w is positive
 smooth rhoeps =-2470.310661 (-179.107611,***********)
         rhomu =-3278.859216 (-116.221789,***********)
       avg vxc =   -0.316237 (  -0.227596,  -0.404877)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -5.572123   -5.624673  -11.196796    -82.114234-2616.757716-2698.871951

 local terms:     true           smooth         local
 rhoeps:        -7.980819   -2470.064629    2462.083810
 rhomu:         -5.216693    -116.220437     111.003744
 spin2:         -5.294318   -3162.319051    3157.024733
 total:        -10.511011   -3278.539488    3268.028477
 val*vef       -11.196796   -3195.937112    3184.740316
 val chg:        2.130156     386.510557    -384.380401
 val mom:       -0.144414    -327.846611     327.702197    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef          -2698.892672      2687.560622       -11.332050
   rhoval*ves             82.992193       -91.768733        -8.776540
   psnuc*ves              -5.173830      -258.804129      -263.977959
   utot                   38.909181      -175.286431      -136.377249
   rho*exc             -2470.310661      2462.083810        -8.226851
   rho*vxc             -3278.859216      3268.028477       -10.830739
   valence chg           387.380401      -384.380401         3.000000
   valence mag          -328.702197       327.702197        -1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=        8.495162  sumtc=        62.933572   ekin=       71.428734
 rhoep=       -8.226851   utot=      -136.377249   ehks=      -73.175366
 mag. mom=    -1.000000

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrealsmooth= T
 wgtsmooth=   2.8284271247461901E-003
 mixrho:  sought 2 iter from file mixm; read 0.  RMS DQ=3.66e0
 mixrho: (warning) scr. and lin-mixed densities had 0 and 72123 negative points
 AMIX: nmix=0 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel=1.83D+00
 mixrho: add corrections to qcell smrho =  0.38444D-07  0.19222D-10
 mixrho: warning. negative smrho; isp number min=       1   90779 -0.67714D-03
 mixrho: warning. negative smrho; isp number min=       2   71667 -0.67036D-03
 unscreened rms difference:  smooth  6.407083   local 11.944385
   screened rms difference:  smooth  4.853815   local 11.944385   tot  3.656525
 mixrho: warning. negative smrho; isp number min=       1   90779 -0.67714D-03
 mixrho: warning. negative smrho; isp number min=       1   90779 -0.67714D-03
 mixrho: warning. negative smrho; isp number min=       1   90779 -0.67714D-03
 mixrho: warning. negative smrho; isp number min=       2   71667 -0.67036D-03
 mixrho: warning. negative smrho; isp number min=       2   71667 -0.67036D-03
 mixrho: warning. negative smrho; isp number min=       2   71667 -0.67036D-03

 iors  : write restart file (binary, mesh density) 

   it  1  of 10    ehf=      -0.754178   ehk=       1.819534
h zbak=1 mmom=-1.0000001 ehf=-.7541783 ehk=1.8195343
 --- BNDFP:  begin iteration 2 of 10
bndfp:start

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1  -54.027093 -191.521059     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 jesm=0 return

 smves:: avg es pot at rmt=-0.033200  avg sphere pot= 0.025175  vconst= 0.033200
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      8.694471   charge   195.521059
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      162446  -6.7713818817799780E-004
 smvxcm: enforce positive smrho_w. Add srshift=   6.7713818818799783E-004
 smooth rhoeps = -987.886461 ( -73.662017,-914.224444)
         rhomu =-1310.095948 ( -49.167920,***********)
       avg vxc =   -0.348692 (  -0.312582,  -0.384802)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 qqx nkapi nkape=           2           2
 qqx nkapi nkape=           2           2

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.282610   -4.877596  -13.160206    -55.594948-1311.702000-1367.296948

 local terms:     true           smooth         local
 rhoeps:        -8.864856    -987.434257     978.569401
 rhomu:         -6.359894     -48.963793      42.603899
 spin2:         -5.306301   -1260.544174    1255.237873
 total:        -11.666195   -1309.507967    1297.841772
 val*vef       -13.160206   -1293.661847    1280.501641
 val chg:        3.128514     194.649573    -191.521059
 val mom:        0.833288    -163.351390     164.184678    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 qqx nkapi nkape=           2           2
 potential shift to crystal energy zero:    0.000157

 potpus  spin 1 : pnu = 2.891276 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000    -8.438830    8.030573    0.095463   -0.636040
 1      3.000000    1.000000    -5.887832    7.215357    0.125265   -0.609243
 2      3.000000    1.000000     6.000000   28.948505    0.488585   -0.089186
 3      3.000000    1.000000     9.000000   39.802437    0.568028   -0.057151

 potpus  spin 2 : pnu = 2.741424 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000    -2.842555    8.735619    0.149589   -0.577371
 1      3.000000    1.000000    -5.887832    7.089103    0.131674   -0.585229
 2      3.000000    1.000000     6.000000   29.309516    0.491109   -0.087353
 3      3.000000    1.000000     9.000000   40.051424    0.569415   -0.056555
 qqx nkapi nkape=           2           2

 Energy terms:             smooth           local           total
   rhoval*vef          -1367.410357      1354.106742       -13.303615
   rhoval*ves             15.851667       -25.851346        -9.999679
   psnuc*ves               1.537274      -268.894255      -267.356981
   utot                    8.694471      -147.372801      -138.678330
   rho*exc              -987.886461       978.569401        -9.317060
   rho*vxc             -1310.095948      1297.841772       -12.254176
   valence chg           194.521059      -191.521059         3.000000
   valence mag          -163.684678       164.184678         0.500000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000
 ... Start MPI k-loop ...
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.4256

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.654075;   3.000000 electrons
         Sum occ. bands:   -3.083069, incl. Bloechl correction: -0.000035
         Mag. moment:       1.000000
Generating TDOS: efermi, and dos window=   -0.6541  -0.5000   0.8459
 bndfp: kpt    1 of    8 k isp=  0.0000  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    1 of    8 k isp=  0.0000  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    2 of    8 k isp=  0.1250  0.1250 -0.1250 1 ndimh nev=    8    5
 bndfp: kpt    2 of    8 k isp=  0.1250  0.1250 -0.1250 2 ndimh nev=    8    5
 bndfp: kpt    3 of    8 k isp=  0.2500  0.2500 -0.2500 1 ndimh nev=    8    5
 bndfp: kpt    3 of    8 k isp=  0.2500  0.2500 -0.2500 2 ndimh nev=    8    5
 bndfp: kpt    4 of    8 k isp=  0.2500  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    4 of    8 k isp=  0.2500  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    5 of    8 k isp=  0.3750  0.1250 -0.1250 1 ndimh nev=    8    5
 bndfp: kpt    5 of    8 k isp=  0.3750  0.1250 -0.1250 2 ndimh nev=    8    5
 bndfp: kpt    6 of    8 k isp=  0.5000  0.2500 -0.2500 1 ndimh nev=    8    5
 bndfp: kpt    6 of    8 k isp=  0.5000  0.2500 -0.2500 2 ndimh nev=    8    5
 bndfp: kpt    7 of    8 k isp=  0.5000  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    7 of    8 k isp=  0.5000  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    8 of    8 k isp=  0.5000  0.2500  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    8 of    8 k isp=  0.5000  0.2500  0.0000 2 ndimh nev=    8    5
 mkrout: site(class) decomposed charge and magnetic moment. class->lmchk

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.861901    5.816744   -2.954843     -0.908817    0.058785   -0.967603
       contr. to mm extrapolated for r>rmt:  -0.076136 est. true mm =-0.984953
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.43137  sum tc=    31.48010  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.39613  sum tc=    31.56723  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.976542   -1.269486    2.891276    2.915580    2.500000    2.915580
 spn 2 0    0.967266   -1.158972    2.741424    2.912383    2.500000    2.912383
 1     1    0.000000   -1.269535    2.850000    2.176030    2.250000    2.850000
 spn 2 1    0.918091   -0.550292    2.850000    2.887892    2.250000    2.850000
 2     0    0.000000   -1.269472    3.147584    3.109380    3.147584    3.147584
 spn 2 0    0.000002   -0.561226    3.147584    3.152233    3.147584    3.152233
 3     0    0.000000   -0.561226    4.102416    4.102416    4.102416    4.102416
 spn 2 0    0.000000   -0.570775    4.102416    4.101489    4.102416    4.102416

 Harris energy:
 sumev=       -3.083069  val*vef=     -13.303615   sumtv=      10.220546
 sumec=      -40.827502  cor*vef=    -103.760168   ttcor=      62.932666
 rhoeps=      -9.317060     utot=    -138.678330    ehar=     -74.842179
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      162446  -6.7714226801995122E-004
 smvxcm: enforce positive smrho_w. Add srshift=   6.7714226802995125E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      162446  -6.7713410833604426E-004
 smvxcm: enforce positive smrho_w. Add srshift=   6.7713410834604429E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      162446  -6.7713818817799769E-004
 smvxcm: enforce positive smrho_w. Add srshift=   6.7713818818799772E-004

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00   -0.00     0.00   -0.00    0.00    -0.00    0.00   -0.00
 shift forces to make zero average correction:           -0.00    0.00   -0.00

 srhov:    -26.126073     12.723637    -13.402436 sumev=   -3.083069   sumtv=   10.319367

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.833546   -2.954843     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 jesm=0 return

 smves:: avg es pot at rmt=-0.117243  avg sphere pot= 0.012540  vconst= 0.117243
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.079474   charge     6.954843
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -6.266271 (  -3.367044,  -2.899227)
         rhomu =   -8.197705 (  -4.492651,  -3.705053)
       avg vxc =   -0.172017 (  -0.148403,  -0.195632)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -5.020915   -8.467592  -13.488507    -13.049081  -10.244679  -23.293760

 local terms:     true           smooth         local
 rhoeps:        -8.774720      -6.242230      -2.532490
 rhomu:         -5.220061      -4.487780      -0.732280
 spin2:         -6.330468      -3.678793      -2.651674
 total:        -11.550528      -8.166574      -3.383955
 val*vef       -13.488507     -12.606435      -0.882071
 val chg:        2.861901       5.816744      -2.954843
 val mom:       -0.908817       0.058785      -0.967603    core:   0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef            -23.310811         9.805232       -13.505580
   rhoval*ves             -4.308560        -6.016766       -10.325327
   psnuc*ves              14.467508      -282.347306      -267.879798
   utot                    5.079474      -144.182036      -139.102562
   rho*exc                -6.266271        -2.532490        -8.798761
   rho*vxc                -8.197705        -3.383955       -11.581659
   valence chg             5.954843        -2.954843         3.000000
   valence mag            -0.032397        -0.967603        -1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       10.319367  sumtc=        63.047332   ekin=       73.366700
 rhoep=       -8.798761   utot=      -139.102562   ehks=      -74.534624
 mag. mom=    -1.000000

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 wgtsmooth=   2.8284271247461901E-003
 mixrho:  sought 2 iter from file mixm; read 1.  RMS DQ=1.81e0  last it=3.66e0
 mixrho: (warning) scr. and lin-mixed densities had 0 and 77427 negative points
 AMIX: nmix=1 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel=9.07D-01
   tj: 0.33171
 mixrho: warning. negative smrho; isp number min=       1   88727 -0.55344D-03
 mixrho: warning. negative smrho; isp number min=       2   71163 -0.55106D-03
 mixrho: warning. negative smrho; isp number min=       1   88727 -0.55344D-03
 mixrho: warning. negative smrho; isp number min=       1   88727 -0.55344D-03
 mixrho: warning. negative smrho; isp number min=       2   71163 -0.55106D-03
 mixrho: warning. negative smrho; isp number min=       2   71163 -0.55106D-03
 mixrho: add corrections to qcell smrho =  0.75812D-07  0.37906D-10
 unscreened rms difference:  smooth  3.178741   local  5.929540
   screened rms difference:  smooth  2.389284   local  5.929540   tot  1.814945
 mixrho: warning. negative smrho; isp number min=       1   88727 -0.55344D-03
 mixrho: warning. negative smrho; isp number min=       2   71163 -0.55106D-03

 iors  : write restart file (binary, mesh density) 

   it  2  of 10    ehf=       0.152721   ehk=       0.460276
 From last iter    ehf=      -0.754178   ehk=       1.819534
 diffe(q)=  0.906899 (1.814945)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=-.9999999 ehf=.1527211 ehk=.4602763
 --- BNDFP:  begin iteration 3 of 10
bndfp:start

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1  -36.252607 -128.512144     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 jesm=0 return

 smves:: avg es pot at rmt=-0.048228  avg sphere pot= 0.020952  vconst= 0.048228
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      4.074898   charge   132.512144
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      159890  -5.5344165688249567E-004
 smvxcm: enforce positive smrho_w. Add srshift=   5.5344165689249570E-004
 smooth rhoeps = -582.769383 ( -46.657436,-536.111947)
         rhomu = -772.367305 ( -32.187613,-740.179692)
       avg vxc =   -0.327379 (  -0.296346,  -0.358413)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 qqx nkapi nkape=           2           2

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -7.193507   -6.076154  -13.269661    -47.755199 -908.052239 -955.807438
 qqx nkapi nkape=           2           2

 local terms:     true           smooth         local
 rhoeps:        -8.839687    -582.419486     573.579799
 rhomu:         -5.979621     -32.025175      26.045554
 spin2:         -5.652773    -739.887438     734.234665
 total:        -11.632394    -771.912613     760.280220
 val*vef       -13.269661    -769.576328     756.306667
 val chg:        3.082514     131.594658    -128.512144
 val mom:        0.251168    -108.748301     108.999469    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 qqx nkapi nkape=           2           2
 potential shift to crystal energy zero:    0.000104
 qqx nkapi nkape=           2           2

 potpus  spin 1 : pnu = 2.915580 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -11.045169    7.782118    0.083218   -0.638250
 1      3.000000    1.000000    -5.887832    7.208279    0.126990   -0.601296
 2      3.000000    1.000000     6.000000   28.993545    0.488796   -0.088973
 3      3.000000    1.000000     9.000000   39.821569    0.568094   -0.057109

 potpus  spin 2 : pnu = 2.912383 2.850000 3.152233 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.622304    7.717818    0.086443   -0.630762
 1      3.000000    1.000000    -5.887832    7.143653    0.129413   -0.592961
 2      3.000000    1.000000     5.787123   28.723361    0.485559   -0.089790
 3      3.000000    1.000000     9.000000   39.946470    0.568836   -0.056804

 Energy terms:             smooth           local           total
   rhoval*vef           -955.917164       942.524353       -13.392811
   rhoval*ves              2.353685       -12.490430       -10.136745
   psnuc*ves               5.796112      -273.432818      -267.636707
   utot                    4.074898      -142.961624      -138.886726
   rho*exc              -582.769383       573.579799        -9.189584
   rho*vxc              -772.367305       760.280220       -12.087085
   valence chg           131.512144      -128.512144         3.000000
   valence mag          -109.000690       108.999469        -0.001221
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000
 ... Start MPI k-loop ...
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.4060

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.613981;   3.000000 electrons
         Sum occ. bands:   -3.039431, incl. Bloechl correction: -0.000038
         Mag. moment:       1.000000
Generating TDOS: efermi, and dos window=   -0.6140  -0.5000   0.8860
 bndfp: kpt    1 of    8 k isp=  0.0000  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    1 of    8 k isp=  0.0000  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    2 of    8 k isp=  0.1250  0.1250 -0.1250 1 ndimh nev=    8    5
 bndfp: kpt    2 of    8 k isp=  0.1250  0.1250 -0.1250 2 ndimh nev=    8    5
 bndfp: kpt    3 of    8 k isp=  0.2500  0.2500 -0.2500 1 ndimh nev=    8    5
 bndfp: kpt    3 of    8 k isp=  0.2500  0.2500 -0.2500 2 ndimh nev=    8    5
 bndfp: kpt    4 of    8 k isp=  0.2500  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    4 of    8 k isp=  0.2500  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    5 of    8 k isp=  0.3750  0.1250 -0.1250 1 ndimh nev=    8    5
 bndfp: kpt    5 of    8 k isp=  0.3750  0.1250 -0.1250 2 ndimh nev=    8    5
 bndfp: kpt    6 of    8 k isp=  0.5000  0.2500 -0.2500 1 ndimh nev=    8    5
 bndfp: kpt    6 of    8 k isp=  0.5000  0.2500 -0.2500 2 ndimh nev=    8    5
 bndfp: kpt    7 of    8 k isp=  0.5000  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    7 of    8 k isp=  0.5000  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    8 of    8 k isp=  0.5000  0.2500  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    8 of    8 k isp=  0.5000  0.2500  0.0000 2 ndimh nev=    8    5
 mkrout: site(class) decomposed charge and magnetic moment. class->lmchk

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.888660    6.200188   -3.311529      0.948290    2.207880   -1.259590
       contr. to mm extrapolated for r>rmt:   0.041076 est. true mm = 0.989366
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.38754  sum tc=    31.49414  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.37770  sum tc=    31.52008  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.975617   -1.230289    2.915580    2.915439    2.500000    2.915439
 spn 2 0    0.970185   -1.194599    2.912383    2.915020    2.500000    2.915020
 1     1    0.942857   -0.614541    2.850000    2.896732    2.250000    2.850000
 spn 2 1    0.000000   -1.194657    2.850000    2.183406    2.250000    2.850000
 2     0    0.000001   -0.621674    3.147584    3.147292    3.147584    3.147584
 spn 2 0    0.000000   -1.194591    3.152233    3.111784    3.147584    3.147584
 3     0    0.000000   -0.628372    4.102416    4.099594    4.102416    4.102416
 spn 2 0    0.000000   -0.628372    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -3.039431  val*vef=     -13.392811   sumtv=      10.353379
 sumec=      -40.765244  cor*vef=    -103.755216   ttcor=      62.989972
 rhoeps=      -9.189584     utot=    -138.886726    ehar=     -74.732958
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      159890  -5.5344500102559775E-004
 smvxcm: enforce positive smrho_w. Add srshift=   5.5344500103559778E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      159890  -5.5343831273939369E-004
 smvxcm: enforce positive smrho_w. Add srshift=   5.5343831274939372E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      159890  -5.5344165688249577E-004
 smvxcm: enforce positive smrho_w. Add srshift=   5.5344165689249580E-004

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00   -0.00     0.00   -0.00    0.00    -0.00    0.00    0.00
 shift forces to make zero average correction:           -0.00    0.00    0.00

 srhov:    -27.800982     14.114355    -13.686627 sumev=   -3.039431   sumtv=   10.647195

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.934165   -3.311529     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 jesm=0 return

 smves:: avg es pot at rmt=-0.114173  avg sphere pot= 0.011672  vconst= 0.114173
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.376394   charge     7.311528
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -6.929373 (  -4.671359,  -2.258013)
         rhomu =   -9.073254 (  -6.462564,  -2.610690)
       avg vxc =   -0.164114 (  -0.182948,  -0.145281)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.796926   -4.853016  -13.649942    -17.396479   -8.031223  -25.427702

 local terms:     true           smooth         local
 rhoeps:        -8.823007      -6.910652      -1.912356
 rhomu:         -6.411461      -6.444436       0.032975
 spin2:         -5.202622      -2.604577      -2.598045
 total:        -11.614083      -9.049013      -2.565070
 val*vef       -13.649942     -13.195625      -0.454317
 val chg:        2.888660       6.200188      -3.311529
 val mom:        0.948290       2.207880      -1.259590    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef            -25.441802        11.777735       -13.664067
   rhoval*ves             -4.022324        -6.408136       -10.430460
   psnuc*ves              14.775112      -282.866415      -268.091302
   utot                    5.376394      -144.637276      -139.260881
   rho*exc                -6.929373        -1.912356        -8.841728
   rho*vxc                -9.073254        -2.565070       -11.638324
   valence chg             6.311528        -3.311529         3.000000
   valence mag             2.259590        -1.259590         1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       10.647195  sumtc=        63.014225   ekin=       73.661420
 rhoep=       -8.841728   utot=      -139.260881   ehks=      -74.441190
 mag. mom=     1.000000

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 wgtsmooth=   2.8284271247461901E-003
 mixrho:  sought 2 iter from file mixm; read 2.  RMS DQ=1.21e0  last it=1.81e0
 mixrho: (warning) scr. and lin-mixed densities had 0 and 76683 negative points
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel=6.04D-01
   tj:-0.35718   0.20388
 mixrho: warning. negative smrho; isp number min=       1   84005 -0.46656D-03
 mixrho: warning. negative smrho; isp number min=       1   84005 -0.46656D-03
 mixrho: warning. negative smrho; isp number min=       2   72315 -0.46816D-03
 mixrho: warning. negative smrho; isp number min=       2   72315 -0.46816D-03
 mixrho: warning. negative smrho; isp number min=       1   84005 -0.46656D-03
 mixrho: warning. negative smrho; isp number min=       2   72315 -0.46816D-03
 mixrho: add corrections to qcell smrho =  0.23171D-07  0.11586D-10
 unscreened rms difference:  smooth  2.119091   local  3.943770
   screened rms difference:  smooth  1.588453   local  3.943770   tot  1.208543
 mixrho: warning. negative smrho; isp number min=       1   84005 -0.46656D-03
 mixrho: warning. negative smrho; isp number min=       2   72315 -0.46816D-03

 iors  : write restart file (binary, mesh density) 

   it  3  of 10    ehf=       0.261942   ehk=       0.553710
 From last iter    ehf=       0.152721   ehk=       0.460276
 diffe(q)=  0.109221 (1.208543)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999999 ehf=.2619417 ehk=.5537104
 --- BNDFP:  begin iteration 4 of 10
bndfp:start

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1  -22.661263  -80.332086     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 jesm=0 return

 smves:: avg es pot at rmt=-0.063724  avg sphere pot= 0.017210  vconst= 0.063724
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      2.855821   charge    84.332086
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      156320  -4.6816419028767562E-004
 smvxcm: enforce positive smrho_w. Add srshift=   4.6816419029767559E-004
 smooth rhoeps = -312.432544 ( -29.626751,-282.805793)
         rhomu = -413.692863 ( -21.973799,-391.719064)
       avg vxc =   -0.315696 (  -0.296423,  -0.334969)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 qqx nkapi nkape=           2           2
 qqx nkapi nkape=           2           2

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
 qqx nkapi nkape=           2           2
   1      -8.547284   -4.936333  -13.483617    -40.213897 -563.011227 -603.225123

 local terms:     true           smooth         local
 rhoeps:        -8.871941    -312.154399     303.282458
 rhomu:         -6.396950     -21.821959      15.425009
 spin2:         -5.279730    -391.509660     386.229929
 total:        -11.676680    -413.331618     401.654938
 val*vef       -13.483617    -416.930868     403.447251
 val chg:        3.042073      83.374159     -80.332086
 val mom:        0.875050     -65.578276      66.453325    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 qqx nkapi nkape=           2           2
 potential shift to crystal energy zero:    0.000066

 potpus  spin 1 : pnu = 2.915439 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -11.025930    7.832608    0.082668   -0.641435
 1      3.000000    1.000000    -5.887832    7.242956    0.125815   -0.605309
 2      3.000000    1.000000     6.000000   28.907326    0.488224   -0.089412
 3      3.000000    1.000000     9.000000   39.764142    0.567799   -0.057245

 potpus  spin 2 : pnu = 2.915020 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.968861    7.654644    0.086188   -0.623003
 1      3.000000    1.000000    -5.887832    7.147888    0.132026   -0.581039
 2      3.000000    1.000000     6.000000   29.211805    0.490222   -0.087880
 3      3.000000    1.000000     9.000000   39.959337    0.568829   -0.056781

 Energy terms:             smooth           local           total
   rhoval*vef           -603.314254       589.736211       -13.578043
   rhoval*ves             -3.564301        -6.725629       -10.289930
   psnuc*ves               9.275943      -277.188925      -267.912982
   utot                    2.855821      -141.957277      -139.101456
   rho*exc              -312.432544       303.282458        -9.150086
   rho*vxc              -413.692863       401.654938       -12.037925
   valence chg            83.332086       -80.332086         3.000000
   valence mag           -65.686142        66.453325         0.767184
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000
 ... Start MPI k-loop ...
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.4258

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.633928;   3.000000 electrons
         Sum occ. bands:   -3.022290, incl. Bloechl correction: -0.000033
         Mag. moment:       1.000000
Generating TDOS: efermi, and dos window=   -0.6339  -0.5000   0.8661
 bndfp: kpt    1 of    8 k isp=  0.0000  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    1 of    8 k isp=  0.0000  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    2 of    8 k isp=  0.1250  0.1250 -0.1250 1 ndimh nev=    8    5
 bndfp: kpt    2 of    8 k isp=  0.1250  0.1250 -0.1250 2 ndimh nev=    8    5
 bndfp: kpt    3 of    8 k isp=  0.2500  0.2500 -0.2500 1 ndimh nev=    8    5
 bndfp: kpt    3 of    8 k isp=  0.2500  0.2500 -0.2500 2 ndimh nev=    8    5
 bndfp: kpt    4 of    8 k isp=  0.2500  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    4 of    8 k isp=  0.2500  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    5 of    8 k isp=  0.3750  0.1250 -0.1250 1 ndimh nev=    8    5
 bndfp: kpt    5 of    8 k isp=  0.3750  0.1250 -0.1250 2 ndimh nev=    8    5
 bndfp: kpt    6 of    8 k isp=  0.5000  0.2500 -0.2500 1 ndimh nev=    8    5
 bndfp: kpt    6 of    8 k isp=  0.5000  0.2500 -0.2500 2 ndimh nev=    8    5
 bndfp: kpt    7 of    8 k isp=  0.5000  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    7 of    8 k isp=  0.5000  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    8 of    8 k isp=  0.5000  0.2500  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    8 of    8 k isp=  0.5000  0.2500  0.0000 2 ndimh nev=    8    5
 mkrout: site(class) decomposed charge and magnetic moment. class->lmchk

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.892764    6.390516   -3.497752      0.955260    2.443426   -1.488166
       contr. to mm extrapolated for r>rmt:   0.035531 est. true mm = 0.990791
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.35295  sum tc=    31.45001  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.31267  sum tc=    31.54186  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.976703   -1.251787    2.915439    2.916867    2.500000    2.916867
 spn 2 0    0.968752   -1.136104    2.915020    2.914175    2.500000    2.914175
 1     1    0.947308   -0.634396    2.850000    2.899398    2.250000    2.850000
 spn 2 1    0.000000   -1.136167    2.850000    2.187180    2.250000    2.850000
 2     0    0.000001   -0.641581    3.147584    3.145983    3.147584    3.147584
 spn 2 0    0.000000   -1.136095    3.147584    3.112895    3.147584    3.147584
 3     0    0.000000   -0.648371    4.102416    4.099066    4.102416    4.102416
 spn 2 0    0.000000   -0.648371    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -3.022290  val*vef=     -13.578043   sumtv=      10.555753
 sumec=      -40.665625  cor*vef=    -103.667736   ttcor=      63.002111
 rhoeps=      -9.150086     utot=    -139.101456    ehar=     -74.693679
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      156320  -4.6816702037621956E-004
 smvxcm: enforce positive smrho_w. Add srshift=   4.6816702038621953E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      156320  -4.6816136019913168E-004
 smvxcm: enforce positive smrho_w. Add srshift=   4.6816136020913165E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      156320  -4.6816419028767562E-004
 smvxcm: enforce positive smrho_w. Add srshift=   4.6816419029767559E-004

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00   -0.00     0.00   -0.00   -0.00    -0.00    0.00    0.00
 shift forces to make zero average correction:           -0.00    0.00    0.00

 srhov:    -29.460730     15.798496    -13.662234 sumev=   -3.022290   sumtv=   10.639944

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.986698   -3.497752     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 jesm=0 return

 smves:: avg es pot at rmt=-0.113509  avg sphere pot= 0.011411  vconst= 0.113509
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.489739   charge     7.497752
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -7.261686 (  -5.010142,  -2.251544)
         rhomu =   -9.511144 (  -6.950154,  -2.560990)
       avg vxc =   -0.162730 (  -0.180184,  -0.145276)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.858528   -4.779307  -13.637835    -18.663150   -7.905423  -26.568573

 local terms:     true           smooth         local
 rhoeps:        -8.826208      -7.243728      -1.582479
 rhomu:         -6.423276      -6.933317       0.510041
 spin2:         -5.195012      -2.554573      -2.640438
 total:        -11.618288      -9.487890      -2.130397
 val*vef       -13.637835     -13.535694      -0.102141
 val chg:        2.892764       6.390516      -3.497752
 val mom:        0.955260       2.443426      -1.488166    core:   0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef            -26.582289        12.930712       -13.651577
   rhoval*ves             -3.924784        -6.488460       -10.413245
   psnuc*ves              14.904262      -282.950231      -268.045969
   utot                    5.489739      -144.719346      -139.229607
   rho*exc                -7.261686        -1.582479        -8.844165
   rho*vxc                -9.511144        -2.130397       -11.641541
   valence chg             6.497752        -3.497752         3.000000
   valence mag             2.488166        -1.488166         1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       10.639944  sumtc=        62.991873   ekin=       73.631817
 rhoep=       -8.844165   utot=      -139.229607   ehks=      -74.441955
 mag. mom=     1.000000

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 wgtsmooth=   2.8284271247461901E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=7.43e-1  last it=1.21e0
 mixrho: (warning) scr. and lin-mixed densities had 0 and 75435 negative points
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel=3.71D-01
   tj:-1.26385  -0.14422
 mixrho: warning. negative smrho; isp number min=       1   64623 -0.25492D-03
 mixrho: warning. negative smrho; isp number min=       1   64623 -0.25492D-03
 mixrho: warning. negative smrho; isp number min=       1   64623 -0.25492D-03
 mixrho: warning. negative smrho; isp number min=       2   77697 -0.25984D-03
 mixrho: warning. negative smrho; isp number min=       2   77697 -0.25984D-03
 mixrho: warning. negative smrho; isp number min=       2   77697 -0.25984D-03
 mixrho: add corrections to qcell smrho =  0.29874D-07  0.14937D-10
 unscreened rms difference:  smooth  1.302397   local  2.425235
   screened rms difference:  smooth  0.975128   local  2.425235   tot  0.742992
 mixrho: warning. negative smrho; isp number min=       1   64623 -0.25492D-03
 mixrho: warning. negative smrho; isp number min=       2   77697 -0.25984D-03

 iors  : write restart file (binary, mesh density) 

   it  4  of 10    ehf=       0.301221   ehk=       0.552945
 From last iter    ehf=       0.261942   ehk=       0.553710
 diffe(q)=  0.039279 (0.742992)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999999 ehf=.301221 ehk=.5529447
 --- BNDFP:  begin iteration 5 of 10
bndfp:start

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.017793   -3.607981     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 jesm=0 return

 smves:: avg es pot at rmt=-0.094218  avg sphere pot= 0.011124  vconst= 0.094218
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.543774   charge     7.607981
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      142320  -2.5983803228458172E-004
 smvxcm: enforce positive smrho_w. Add srshift=   2.5983803229458170E-004
 smooth rhoeps =   -7.673044 (  -5.556509,  -2.116535)
         rhomu =  -10.053335 (  -7.726833,  -2.326502)
       avg vxc =   -0.262015 (  -0.274052,  -0.249978)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       768
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       768
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       768
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       768
 vxcnsp (warning): negative rho: min val =  -2.12E-02
 qqx nkapi nkape=           2           2
 qqx nkapi nkape=           2           2

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -9.905262   -3.871092  -13.776354    -20.089074   -7.191326  -27.280400

 local terms:     true           smooth         local
 qqx nkapi nkape=           2           2
 rhoeps:        -8.917033      -7.541778      -1.375254
 rhomu:         -6.797254      -7.625084       0.827830
 spin2:         -4.942290      -2.258092      -2.684198
 total:        -11.739544      -9.883176      -1.856368
 val*vef       -13.776354     -13.812027       0.035673
 val chg:        2.959255       6.567236      -3.607981
 val mom:        1.451232       3.083787      -1.632555    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 qqx nkapi nkape=           2           2
 potential shift to crystal energy zero:    0.000006

 potpus  spin 1 : pnu = 2.916867 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -11.224426    7.890294    0.081094   -0.645119
 1      3.000000    1.000000    -5.887832    7.299283    0.124502   -0.609077
 2      3.000000    1.000000     6.000000   28.778403    0.487298   -0.090089
 3      3.000000    1.000000     9.000000   39.670954    0.567292   -0.057471

 potpus  spin 2 : pnu = 2.914175 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.855570    7.684795    0.087707   -0.614957
 1      3.000000    1.000000    -5.887832    7.202121    0.133927   -0.570419
 2      3.000000    1.000000     6.000000   29.161710    0.489523   -0.088195
 3      3.000000    1.000000     9.000000   39.885408    0.568282   -0.056972

 Energy terms:             smooth           local           total
   rhoval*vef            -27.323926        13.504025       -13.819901
   rhoval*ves             -3.844184        -6.623527       -10.467711
   psnuc*ves              14.931733      -283.152102      -268.220370
   utot                    5.543774      -144.887815      -139.344040
   rho*exc                -7.673044        -1.375254        -9.048298
   rho*vxc               -10.053335        -1.856368       -11.909703
   valence chg             6.607981        -3.607981         3.000000
   valence mag             3.165203        -1.632555         1.532648
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000
 ... Start MPI k-loop ...
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.4352

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.644712;   3.000000 electrons
         Sum occ. bands:   -2.981350, incl. Bloechl correction: -0.000027
         Mag. moment:       1.000000
Generating TDOS: efermi, and dos window=   -0.6447  -0.5000   0.8553
 bndfp: kpt    1 of    8 k isp=  0.0000  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    1 of    8 k isp=  0.0000  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    2 of    8 k isp=  0.1250  0.1250 -0.1250 1 ndimh nev=    8    5
 bndfp: kpt    2 of    8 k isp=  0.1250  0.1250 -0.1250 2 ndimh nev=    8    5
 bndfp: kpt    3 of    8 k isp=  0.2500  0.2500 -0.2500 1 ndimh nev=    8    5
 bndfp: kpt    3 of    8 k isp=  0.2500  0.2500 -0.2500 2 ndimh nev=    8    5
 bndfp: kpt    4 of    8 k isp=  0.2500  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    4 of    8 k isp=  0.2500  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    5 of    8 k isp=  0.3750  0.1250 -0.1250 1 ndimh nev=    8    5
 bndfp: kpt    5 of    8 k isp=  0.3750  0.1250 -0.1250 2 ndimh nev=    8    5
 bndfp: kpt    6 of    8 k isp=  0.5000  0.2500 -0.2500 1 ndimh nev=    8    5
 bndfp: kpt    6 of    8 k isp=  0.5000  0.2500 -0.2500 2 ndimh nev=    8    5
 bndfp: kpt    7 of    8 k isp=  0.5000  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    7 of    8 k isp=  0.5000  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    8 of    8 k isp=  0.5000  0.2500  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    8 of    8 k isp=  0.5000  0.2500  0.0000 2 ndimh nev=    8    5
 mkrout: site(class) decomposed charge and magnetic moment. class->lmchk

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905198    7.988297   -5.083099      0.957405    1.560283   -0.602879
       contr. to mm extrapolated for r>rmt:   0.035312 est. true mm = 0.992716
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.30517  sum tc=    31.40600  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.23524  sum tc=    31.56264  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.978171   -1.265419    2.916867    2.919428    2.500000    2.919428
 spn 2 0    0.973896   -1.070867    2.914175    2.912894    2.500000    2.912894
 1     1    0.953130   -0.645062    2.850000    2.903795    2.250000    2.850000
 spn 2 1    0.000001   -1.070941    2.850000    2.188818    2.250000    2.850000
 2     0    0.000001   -0.652245    3.147584    3.144447    3.147584    3.147584
 spn 2 0    0.000000   -1.070839    3.147584    3.113245    3.147584    3.147584
 3     0    0.000000   -0.659144    4.102416    4.098415    4.102416    4.102416
 spn 2 0    0.000000   -0.659144    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -2.981350  val*vef=     -13.819901   sumtv=      10.838551
 sumec=      -40.540411  cor*vef=    -103.537370   ttcor=      62.996960
 rhoeps=      -9.048298     utot=    -139.344040    ehar=     -74.556828
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      142320  -2.5983959081804086E-004
 smvxcm: enforce positive smrho_w. Add srshift=   2.5983959082804083E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      142320  -2.5983647375112258E-004
 smvxcm: enforce positive smrho_w. Add srshift=   2.5983647376112256E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      142320  -2.5983803228458172E-004
 smvxcm: enforce positive smrho_w. Add srshift=   2.5983803229458170E-004

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00   -0.00   -0.00     0.00    0.00   -0.00    -0.00   -0.00    0.00
 shift forces to make zero average correction:           -0.00   -0.00    0.00

 srhov:    -34.522811     20.833928    -13.688883 sumev=   -2.981350   sumtv=   10.707533

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.433916   -5.083099     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 jesm=0 return

 smves:: avg es pot at rmt=-0.111291  avg sphere pot= 0.010055  vconst= 0.111291
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      6.006901   charge     9.083099
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -9.911449 (  -5.799753,  -4.111696)
         rhomu =  -12.984074 (  -7.860713,  -5.123362)
       avg vxc =   -0.159256 (  -0.175713,  -0.142799)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.909124   -4.769529  -13.678654    -21.190432  -15.068879  -36.259311

 local terms:     true           smooth         local
 rhoeps:        -8.842725      -9.895730       1.053006
 rhomu:         -6.438329      -7.845655       1.407326
 spin2:         -5.201584      -5.118066      -0.083518
 total:        -11.639913     -12.963721       1.323808
 val*vef       -13.678654     -16.566131       2.887477
 val chg:        2.905198       7.988297      -5.083099
 val mom:        0.957405       1.560283      -0.602879    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef            -36.271646        22.580611       -13.691035
   rhoval*ves             -3.483147        -6.950713       -10.433861
   psnuc*ves              15.496949      -283.583721      -268.086772
   utot                    6.006901      -145.267217      -139.260316
   rho*exc                -9.911449         1.053006        -8.858443
   rho*vxc               -12.984074         1.323808       -11.660267
   valence chg             8.083099        -5.083099         3.000000
   valence mag             1.602879        -0.602879         1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       10.707533  sumtc=        62.968642   ekin=       73.676176
 rhoep=       -8.858443   utot=      -139.260316   ehks=      -74.442584
 mag. mom=     1.000000

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 wgtsmooth=   2.8284271247461901E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.18e-2  last it=7.43e-1
 mixrho: (warning) scr. and lin-mixed densities had 0 and 66411 negative points
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel=5.89D-03
   tj:-1.02004   0.64031
 mixrho: warning. negative smrho; isp number min=       1   62895 -0.17369D-03
 mixrho: warning. negative smrho; isp number min=       1   62895 -0.17369D-03
 mixrho: warning. negative smrho; isp number min=       1   62895 -0.17369D-03
 mixrho: warning. negative smrho; isp number min=       2   74013 -0.17638D-03
 mixrho: warning. negative smrho; isp number min=       2   74013 -0.17638D-03
 mixrho: warning. negative smrho; isp number min=       2   74013 -0.17638D-03
 mixrho: add corrections to qcell smrho =  0.52642D-07  0.26321D-10
 unscreened rms difference:  smooth  0.020798   local  0.039624
   screened rms difference:  smooth  0.016291   local  0.039624   tot  0.011782
 mixrho: warning. negative smrho; isp number min=       1   62895 -0.17369D-03
 mixrho: warning. negative smrho; isp number min=       2   74013 -0.17638D-03

 iors  : write restart file (binary, mesh density) 

   it  5  of 10    ehf=       0.438072   ehk=       0.552316
 From last iter    ehf=       0.301221   ehk=       0.552945
 diffe(q)=  0.136851 (0.011782)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999999 ehf=.4380719 ehk=.5523162
 --- BNDFP:  begin iteration 6 of 10
bndfp:start

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.535968   -5.444863     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 jesm=0 return

 smves:: avg es pot at rmt=-0.099835  avg sphere pot= 0.010459  vconst= 0.099835
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.778362   charge     9.444863
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      136908  -1.7637642596807923E-004
 smvxcm: enforce positive smrho_w. Add srshift=   1.7637642597807923E-004
 smooth rhoeps =  -10.750459 (  -5.973674,  -4.776785)
         rhomu =  -14.085098 (  -8.003478,  -6.081621)
       avg vxc =   -0.247367 (  -0.255053,  -0.239681)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       416
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       416
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       416
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       416
 vxcnsp (warning): negative rho: min val =  -6.50E-03

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -9.221886   -4.497551  -13.719436    -21.092451  -17.372530  -38.464981

 local terms:     true           smooth         local
 rhoeps:        -8.876283     -10.659581       1.783298
 rhomu:         -6.558382      -7.935064       1.376682
 spin2:         -5.125987      -6.032172       0.906185
 total:        -11.684369     -13.967236       2.282867
 val*vef       -13.719436     -17.735107       4.015671
 val chg:        2.941746       8.386609      -5.444863
 val mom:        1.112151       1.294300      -0.182148    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 qqx nkapi nkape=           2           2
 qqx nkapi nkape=           2           2
 potential shift to crystal energy zero:    0.000009

 potpus  spin 1 : pnu = 2.919428 2.850000 3.147584 4.102416
 qqx nkapi nkape=           2           2
 qqx nkapi nkape=           2           2

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -11.597648    7.862673    0.079685   -0.644868
 1      3.000000    1.000000    -5.887832    7.301093    0.124800   -0.607538
 2      3.000000    1.000000     6.000000   28.779693    0.487266   -0.090090
 3      3.000000    1.000000     9.000000   39.666912    0.567255   -0.057482

 potpus  spin 2 : pnu = 2.912894 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.687853    7.747989    0.087560   -0.619478
 1      3.000000    1.000000    -5.887832    7.225210    0.132185   -0.576918
 2      3.000000    1.000000     6.000000   29.075368    0.488961   -0.088627
 3      3.000000    1.000000     9.000000   39.829383    0.567994   -0.057105

 Energy terms:             smooth           local           total
   rhoval*vef            -38.499683        24.745498       -13.754184
   rhoval*ves             -3.672182        -6.782641       -10.454823
   psnuc*ves              15.228905      -283.397084      -268.168178
   utot                    5.778362      -145.089863      -139.311501
   rho*exc               -10.750459         1.783298        -8.967161
   rho*vxc               -14.085098         2.282867       -11.802231
   valence chg             8.444863        -5.444863         3.000000
   valence mag             1.347798        -0.182148         1.165649
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000
 ... Start MPI k-loop ...
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.4191

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.628319;   3.000000 electrons
         Sum occ. bands:   -2.978450, incl. Bloechl correction: -0.000027
         Mag. moment:       1.000000
Generating TDOS: efermi, and dos window=   -0.6283  -0.5000   0.8717
 bndfp: kpt    1 of    8 k isp=  0.0000  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    1 of    8 k isp=  0.0000  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    2 of    8 k isp=  0.1250  0.1250 -0.1250 1 ndimh nev=    8    5
 bndfp: kpt    2 of    8 k isp=  0.1250  0.1250 -0.1250 2 ndimh nev=    8    5
 bndfp: kpt    3 of    8 k isp=  0.2500  0.2500 -0.2500 1 ndimh nev=    8    5
 bndfp: kpt    3 of    8 k isp=  0.2500  0.2500 -0.2500 2 ndimh nev=    8    5
 bndfp: kpt    4 of    8 k isp=  0.2500  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    4 of    8 k isp=  0.2500  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    5 of    8 k isp=  0.3750  0.1250 -0.1250 1 ndimh nev=    8    5
 bndfp: kpt    5 of    8 k isp=  0.3750  0.1250 -0.1250 2 ndimh nev=    8    5
 bndfp: kpt    6 of    8 k isp=  0.5000  0.2500 -0.2500 1 ndimh nev=    8    5
 bndfp: kpt    6 of    8 k isp=  0.5000  0.2500 -0.2500 2 ndimh nev=    8    5
 bndfp: kpt    7 of    8 k isp=  0.5000  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    7 of    8 k isp=  0.5000  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    8 of    8 k isp=  0.5000  0.2500  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    8 of    8 k isp=  0.5000  0.2500  0.0000 2 ndimh nev=    8    5
 mkrout: site(class) decomposed charge and magnetic moment. class->lmchk

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.902359    7.175377   -4.273017      0.957897    1.901923   -0.944026
       contr. to mm extrapolated for r>rmt:   0.034914 est. true mm = 0.992811
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.30538  sum tc=    31.42154  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.25170  sum tc=    31.54304  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.977139   -1.249444    2.919428    2.920768    2.500000    2.920768
 spn 2 0    0.972231   -1.100368    2.912894    2.916840    2.500000    2.916840
 1     1    0.952989   -0.628635    2.850000    2.905102    2.250000    2.850000
 spn 2 1    0.000000   -1.100416    2.850000    2.186412    2.250000    2.850000
 2     0    0.000001   -0.637874    3.147584    3.144636    3.147584    3.147584
 spn 2 0    0.000000   -1.100357    3.147584    3.112387    3.147584    3.147584
 3     0    0.000000   -0.647794    4.102416    4.098387    4.102416    4.102416
 spn 2 0    0.000000   -0.647794    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -2.978450  val*vef=     -13.754184   sumtv=      10.775734
 sumec=      -40.557071  cor*vef=    -103.539879   ttcor=      62.982808
 rhoeps=      -8.967161     utot=    -139.311501    ehar=     -74.520120
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      136908  -1.7637748588773781E-004
 smvxcm: enforce positive smrho_w. Add srshift=   1.7637748589773781E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      136908  -1.7637536604842065E-004
 smvxcm: enforce positive smrho_w. Add srshift=   1.7637536605842065E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      136908  -1.7637642596807923E-004
 smvxcm: enforce positive smrho_w. Add srshift=   1.7637642597807923E-004

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00   -0.00    0.00    -0.00    0.00   -0.00     0.00   -0.00    0.00
 shift forces to make zero average correction:            0.00   -0.00    0.00

 srhov:    -31.743743     18.077093    -13.666650 sumev=   -2.978450   sumtv=   10.688199

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.205396   -4.273017     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 jesm=0 return

 smves:: avg es pot at rmt=-0.112032  avg sphere pot= 0.010640  vconst= 0.112032
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1  -0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.797161   charge     8.273017
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -8.492232 (  -5.290879,  -3.201353)
         rhomu =  -11.121854 (  -7.246932,  -3.874921)
       avg vxc =   -0.159249 (  -0.175708,  -0.142790)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.886252   -4.779704  -13.665957    -19.516734  -11.651449  -31.168183

 local terms:     true           smooth         local
 rhoeps:        -8.838176      -8.475981      -0.362195
 rhomu:         -6.434382      -7.231531       0.797149
 spin2:         -5.199558      -3.869280      -1.330278
 total:        -11.633940     -11.100811      -0.533130
 val*vef       -13.665957     -14.877027       1.211070
 val chg:        2.902359       7.175377      -4.273017
 val mom:        0.957897       1.901923      -0.944026    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef            -31.180910        17.502190       -13.678720
   rhoval*ves             -3.655904        -6.770165       -10.426069
   psnuc*ves              15.250226      -283.306448      -268.056222
   utot                    5.797161      -145.038306      -139.241145
   rho*exc                -8.492232        -0.362195        -8.854427
   rho*vxc               -11.121854        -0.533130       -11.654983
   valence chg             7.273017        -4.273017         3.000000
   valence mag             1.944026        -0.944026         1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       10.688199  sumtc=        62.964581   ekin=       73.652780
 rhoep=       -8.854427   utot=      -139.241145   ehks=      -74.442792
 mag. mom=     1.000000

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 wgtsmooth=   2.8284271247461901E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=9.84e-3  last it=1.18e-2
 mixrho: (warning) scr. and lin-mixed densities had 0 and 62311 negative points
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel=4.92D-03
   tj: 0.36956  -0.00253
 mixrho: warning. negative smrho; isp number min=       1   60883 -0.15313D-03
 mixrho: warning. negative smrho; isp number min=       1   60883 -0.15313D-03
 mixrho: warning. negative smrho; isp number min=       1   60883 -0.15313D-03
 mixrho: warning. negative smrho; isp number min=       2   73677 -0.15570D-03
 mixrho: warning. negative smrho; isp number min=       2   73677 -0.15570D-03
 mixrho: warning. negative smrho; isp number min=       2   73677 -0.15570D-03
 mixrho: add corrections to qcell smrho =  0.46824D-07  0.23412D-10
 unscreened rms difference:  smooth  0.017045   local  0.032733
   screened rms difference:  smooth  0.014440   local  0.032733   tot  0.009843
 mixrho: warning. negative smrho; isp number min=       1   60883 -0.15313D-03
 mixrho: warning. negative smrho; isp number min=       2   73677 -0.15570D-03

 iors  : write restart file (binary, mesh density) 

   it  6  of 10    ehf=       0.474780   ehk=       0.552108
 From last iter    ehf=       0.438072   ehk=       0.552316
 diffe(q)=  0.036708 (0.009843)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999999 ehf=.4747801 ehk=.5521077
 --- BNDFP:  begin iteration 7 of 10
bndfp:start

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.290702   -4.575421     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 jesm=0 return

 smves:: avg es pot at rmt=-0.101534  avg sphere pot= 0.010555  vconst= 0.101534
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.779312   charge     8.575421
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      134560  -1.5569732187924198E-004
 smvxcm: enforce positive smrho_w. Add srshift=   1.5569732188924198E-004
 smooth rhoeps =   -9.132373 (  -5.624843,  -3.507529)
         rhomu =  -11.961470 (  -7.688355,  -4.273116)
       avg vxc =   -0.241137 (  -0.249717,  -0.232556)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       416
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       416
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       416
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       416
 vxcnsp (warning): negative rho: min val =  -5.87E-03
 qqx nkapi nkape=           2           2
 qqx nkapi nkape=           2           2
 qqx nkapi nkape=           2           2

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -9.188023   -4.524534  -13.712557    -20.386695  -12.643953  -33.030648

 local terms:     true           smooth         local
 rhoeps:        -8.871098      -9.051501       0.180403
 rhomu:         -6.545462      -7.626628       1.081167
 spin2:         -5.132063      -4.229964      -0.902100
 total:        -11.677525     -11.856592       0.179067
 val*vef       -13.712557     -15.618866       1.906309
 val chg:        2.936429       7.511850      -4.575421
 val mom:        1.097891       1.949585      -0.851693    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 qqx nkapi nkape=           2           2
 potential shift to crystal energy zero:    0.000008

 potpus  spin 1 : pnu = 2.920768 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -11.802420    7.852952    0.078807   -0.645576
 1      3.000000    1.000000    -5.887832    7.303637    0.124707   -0.607876
 2      3.000000    1.000000     6.000000   28.773121    0.487219   -0.090125
 3      3.000000    1.000000     9.000000   39.662175    0.567229   -0.057493

 potpus  spin 2 : pnu = 2.916840 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -11.220580    7.711897    0.084988   -0.621483
 1      3.000000    1.000000    -5.887832    7.228350    0.131975   -0.577697
 2      3.000000    1.000000     6.000000   29.064410    0.488887   -0.088683
 3      3.000000    1.000000     9.000000   39.822054    0.567955   -0.057122

 Energy terms:             smooth           local           total
   rhoval*vef            -33.063021        19.318055       -13.744966
   rhoval*ves             -3.663360        -6.787620       -10.450980
   psnuc*ves              15.221984      -283.361950      -268.139966
   utot                    5.779312      -145.074785      -139.295473
   rho*exc                -9.132373         0.180403        -8.951969
   rho*vxc               -11.961470         0.179067       -11.782403
   valence chg             7.575421        -4.575421         3.000000
   valence mag             2.002835        -0.851693         1.151142
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000
 ... Start MPI k-loop ...
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.4105

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.628465;   3.000000 electrons
         Sum occ. bands:   -2.981018, incl. Bloechl correction: -0.000026
         Mag. moment:       1.000000
Generating TDOS: efermi, and dos window=   -0.6285  -0.5000   0.8715
 bndfp: kpt    1 of    8 k isp=  0.0000  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    1 of    8 k isp=  0.0000  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    2 of    8 k isp=  0.1250  0.1250 -0.1250 1 ndimh nev=    8    5
 bndfp: kpt    2 of    8 k isp=  0.1250  0.1250 -0.1250 2 ndimh nev=    8    5
 bndfp: kpt    3 of    8 k isp=  0.2500  0.2500 -0.2500 1 ndimh nev=    8    5
 bndfp: kpt    3 of    8 k isp=  0.2500  0.2500 -0.2500 2 ndimh nev=    8    5
 bndfp: kpt    4 of    8 k isp=  0.2500  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    4 of    8 k isp=  0.2500  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    5 of    8 k isp=  0.3750  0.1250 -0.1250 1 ndimh nev=    8    5
 bndfp: kpt    5 of    8 k isp=  0.3750  0.1250 -0.1250 2 ndimh nev=    8    5
 bndfp: kpt    6 of    8 k isp=  0.5000  0.2500 -0.2500 1 ndimh nev=    8    5
 bndfp: kpt    6 of    8 k isp=  0.5000  0.2500 -0.2500 2 ndimh nev=    8    5
 bndfp: kpt    7 of    8 k isp=  0.5000  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    7 of    8 k isp=  0.5000  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    8 of    8 k isp=  0.5000  0.2500  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    8 of    8 k isp=  0.5000  0.2500  0.0000 2 ndimh nev=    8    5
 mkrout: site(class) decomposed charge and magnetic moment. class->lmchk

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.903008    7.198155   -4.295148      0.957980    1.852251   -0.894272
       contr. to mm extrapolated for r>rmt:   0.034984 est. true mm = 0.992964
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.30803  sum tc=    31.42285  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.25513  sum tc=    31.54289  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.977096   -1.249738    2.920768    2.921135    2.500000    2.921135
 spn 2 0    0.972514   -1.102509    2.916840    2.917233    2.500000    2.917233
 1     1    0.953397   -0.628770    2.850000    2.905587    2.250000    2.850000
 spn 2 1    0.000000   -1.102556    2.850000    2.186094    2.250000    2.850000
 2     0    0.000001   -0.638216    3.147584    3.144541    3.147584    3.147584
 spn 2 0    0.000000   -1.102499    3.147584    3.112277    3.147584    3.147584
 3     0    0.000000   -0.648622    4.102416    4.098334    4.102416    4.102416
 spn 2 0    0.000000   -0.648622    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -2.981018  val*vef=     -13.744966   sumtv=      10.763949
 sumec=      -40.563153  cor*vef=    -103.536847   ttcor=      62.973694
 rhoeps=      -8.951969     utot=    -139.295473    ehar=     -74.509799
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      134560  -1.5569825693425172E-004
 smvxcm: enforce positive smrho_w. Add srshift=   1.5569825694425172E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      134560  -1.5569638682423221E-004
 smvxcm: enforce positive smrho_w. Add srshift=   1.5569638683423222E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      134560  -1.5569732187924195E-004
 smvxcm: enforce positive smrho_w. Add srshift=   1.5569732188924196E-004

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00   -0.00   -0.00    -0.00    0.00    0.00     0.00   -0.00   -0.00
 shift forces to make zero average correction:            0.00   -0.00   -0.00

 srhov:    -31.428782     17.749164    -13.679618 sumev=   -2.981018   sumtv=   10.698600

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.211639   -4.295148     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 jesm=0 return

 smves:: avg es pot at rmt=-0.111975  avg sphere pot= 0.010617  vconst= 0.111975
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.804101   charge     8.295148
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -8.526520 (  -5.273943,  -3.252577)
         rhomu =  -11.166589 (  -7.214877,  -3.951713)
       avg vxc =   -0.158916 (  -0.175322,  -0.142509)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.888684   -4.785075  -13.673759    -19.433051  -11.862910  -31.295961

 local terms:     true           smooth         local
 rhoeps:        -8.839825      -8.510383      -0.329441
 rhomu:         -6.435546      -7.199559       0.764013
 spin2:         -5.200561      -3.946135      -1.254426
 total:        -11.636107     -11.145694      -0.490413
 val*vef       -13.673759     -14.913939       1.240180
 val chg:        2.903008       7.198155      -4.295148
 val mom:        0.957980       1.852251      -0.894272    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef            -31.308621        17.622166       -13.686455
   rhoval*ves             -3.648071        -6.784003       -10.432074
   psnuc*ves              15.256272      -283.326591      -268.070319
   utot                    5.804101      -145.055297      -139.251196
   rho*exc                -8.526520        -0.329441        -8.855961
   rho*vxc               -11.166589        -0.490413       -11.657002
   valence chg             7.295148        -4.295148         3.000000
   valence mag             1.894272        -0.894272         1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       10.698600  sumtc=        62.965740   ekin=       73.664340
 rhoep=       -8.855961   utot=      -139.251196   ehks=      -74.442817
 mag. mom=     1.000000

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 wgtsmooth=   2.8284271247461901E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=2.02e-3  last it=9.84e-3
 mixrho: (warning) scr. and lin-mixed densities had 0 and 60379 negative points
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel=1.01D-03
   tj: 0.23217   0.26824
 mixrho: warning. negative smrho; isp number min=       1   59413 -0.13963D-03
 mixrho: warning. negative smrho; isp number min=       2   73029 -0.14203D-03
 mixrho: warning. negative smrho; isp number min=       1   59413 -0.13963D-03
 mixrho: warning. negative smrho; isp number min=       1   59413 -0.13963D-03
 mixrho: warning. negative smrho; isp number min=       2   73029 -0.14203D-03
 mixrho: warning. negative smrho; isp number min=       2   73029 -0.14203D-03
 mixrho: add corrections to qcell smrho =  0.46123D-07  0.23062D-10
 unscreened rms difference:  smooth  0.003533   local  0.006632
   screened rms difference:  smooth  0.003479   local  0.006632   tot  0.002023
 mixrho: warning. negative smrho; isp number min=       1   59413 -0.13963D-03
 mixrho: warning. negative smrho; isp number min=       2   73029 -0.14203D-03

 iors  : write restart file (binary, mesh density) 

   it  7  of 10    ehf=       0.485101   ehk=       0.552083
 From last iter    ehf=       0.474780   ehk=       0.552108
 diffe(q)=  0.010320 (0.002023)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999999 ehf=.4851006 ehk=.5520825
 --- BNDFP:  begin iteration 8 of 10
bndfp:start

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.272126   -4.509570     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 jesm=0 return

 smves:: avg es pot at rmt=-0.102536  avg sphere pot= 0.010579  vconst= 0.102536
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.776001   charge     8.509570
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      132442  -1.4203276147674727E-004
 smvxcm: enforce positive smrho_w. Add srshift=   1.4203276148674727E-004
 smooth rhoeps =   -9.004061 (  -5.558439,  -3.445623)
         rhomu =  -11.793054 (  -7.600480,  -4.192574)
       avg vxc =   -0.237370 (  -0.245928,  -0.228812)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       352
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       352
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       352
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       352
 vxcnsp (warning): negative rho: min val =  -4.96E-03
 qqx nkapi nkape=           2           2
 qqx nkapi nkape=           2           2
 qqx nkapi nkape=           2           2

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -9.143312   -4.565135  -13.708447    -20.188961  -12.430913  -32.619873

 local terms:     true           smooth         local
 rhoeps:        -8.867396      -8.929527       0.062131
 rhomu:         -6.529632      -7.543353       1.013721
 spin2:         -5.142957      -4.153040      -0.989917
 total:        -11.672589     -11.696393       0.023804
 val*vef       -13.708447     -15.463215       1.754768
 val chg:        2.933297       7.442868      -4.509570
 val mom:        1.076894       1.943576      -0.866682    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 qqx nkapi nkape=           2           2
 potential shift to crystal energy zero:    0.000008

 potpus  spin 1 : pnu = 2.921135 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -11.859693    7.850777    0.078567   -0.645745
 1      3.000000    1.000000    -5.887832    7.304875    0.124684   -0.607932
 2      3.000000    1.000000     6.000000   28.770301    0.487196   -0.090140
 3      3.000000    1.000000     9.000000   39.659868    0.567216   -0.057499

 potpus  spin 2 : pnu = 2.917233 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -11.276385    7.712482    0.084653   -0.622091
 1      3.000000    1.000000    -5.887832    7.230709    0.131812   -0.578309
 2      3.000000    1.000000     6.000000   29.056071    0.488833   -0.088725
 3      3.000000    1.000000     9.000000   39.816618    0.567927   -0.057135

 Energy terms:             smooth           local           total
   rhoval*vef            -32.650748        18.911390       -13.739358
   rhoval*ves             -3.665964        -6.783696       -10.449661
   psnuc*ves              15.217966      -283.345904      -268.127937
   utot                    5.776001      -145.064800      -139.288799
   rho*exc                -9.004061         0.062131        -8.941930
   rho*vxc               -11.793054         0.023804       -11.769251
   valence chg             7.509570        -4.509570         3.000000
   valence mag             1.995105        -0.866682         1.128422
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000
 ... Start MPI k-loop ...
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.4384

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.627683;   3.000000 electrons
         Sum occ. bands:   -2.981584, incl. Bloechl correction: -0.000026
         Mag. moment:       1.000000
Generating TDOS: efermi, and dos window=   -0.6277  -0.5000   0.8723
 bndfp: kpt    1 of    8 k isp=  0.0000  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    1 of    8 k isp=  0.0000  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    2 of    8 k isp=  0.1250  0.1250 -0.1250 1 ndimh nev=    8    5
 bndfp: kpt    2 of    8 k isp=  0.1250  0.1250 -0.1250 2 ndimh nev=    8    5
 bndfp: kpt    3 of    8 k isp=  0.2500  0.2500 -0.2500 1 ndimh nev=    8    5
 bndfp: kpt    3 of    8 k isp=  0.2500  0.2500 -0.2500 2 ndimh nev=    8    5
 bndfp: kpt    4 of    8 k isp=  0.2500  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    4 of    8 k isp=  0.2500  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    5 of    8 k isp=  0.3750  0.1250 -0.1250 1 ndimh nev=    8    5
 bndfp: kpt    5 of    8 k isp=  0.3750  0.1250 -0.1250 2 ndimh nev=    8    5
 bndfp: kpt    6 of    8 k isp=  0.5000  0.2500 -0.2500 1 ndimh nev=    8    5
 bndfp: kpt    6 of    8 k isp=  0.5000  0.2500 -0.2500 2 ndimh nev=    8    5
 bndfp: kpt    7 of    8 k isp=  0.5000  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    7 of    8 k isp=  0.5000  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    8 of    8 k isp=  0.5000  0.2500  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    8 of    8 k isp=  0.5000  0.2500  0.0000 2 ndimh nev=    8    5
 mkrout: site(class) decomposed charge and magnetic moment. class->lmchk

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.903106    7.160358   -4.257251      0.958102    1.854958   -0.896856
       contr. to mm extrapolated for r>rmt:   0.034941 est. true mm = 0.993043
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.30897  sum tc=    31.42405  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.25709  sum tc=    31.54188  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.977019   -1.249048    2.921135    2.921381    2.500000    2.921381
 spn 2 0    0.972502   -1.104557    2.917233    2.917618    2.500000    2.917618
 1     1    0.953584   -0.627978    2.850000    2.905889    2.250000    2.850000
 spn 2 1    0.000000   -1.104602    2.850000    2.185870    2.250000    2.850000
 2     0    0.000001   -0.637701    3.147584    3.144507    3.147584    3.147584
 spn 2 0    0.000000   -1.104547    3.147584    3.112197    3.147584    3.147584
 3     0    0.000000   -0.648541    4.102416    4.098306    4.102416    4.102416
 spn 2 0    0.000000   -0.648541    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -2.981584  val*vef=     -13.739358   sumtv=      10.757774
 sumec=      -40.566058  cor*vef=    -103.535775   ttcor=      62.969717
 rhoeps=      -8.941930     utot=    -139.288799    ehar=     -74.503237
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      132442  -1.4203361427715382E-004
 smvxcm: enforce positive smrho_w. Add srshift=   1.4203361428715382E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      132442  -1.4203190867634070E-004
 smvxcm: enforce positive smrho_w. Add srshift=   1.4203190868634070E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      132442  -1.4203276147674724E-004
 smvxcm: enforce positive smrho_w. Add srshift=   1.4203276148674725E-004

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00   -0.00    -0.00    0.00    0.00     0.00   -0.00   -0.00
 shift forces to make zero average correction:            0.00   -0.00   -0.00

 srhov:    -31.184704     17.501907    -13.682796 sumev=   -2.981584   sumtv=   10.701212

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.200948   -4.257251     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 jesm=0 return

 smves:: avg es pot at rmt=-0.111992  avg sphere pot= 0.010641  vconst= 0.111992
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.795712   charge     8.257251
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -8.460527 (  -5.240382,  -3.220145)
         rhomu =  -11.079953 (  -7.170180,  -3.909773)
       avg vxc =   -0.158776 (  -0.175149,  -0.142403)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.888637   -4.787126  -13.675762    -19.312339  -11.746790  -31.059129

 local terms:     true           smooth         local
 rhoeps:        -8.840191      -8.444405      -0.395786
 rhomu:         -6.435836      -7.154882       0.719046
 spin2:         -5.200753      -3.904195      -1.296558
 total:        -11.636589     -11.059077      -0.577512
 val*vef       -13.675762     -14.833223       1.157461
 val chg:        2.903106       7.160358      -4.257251
 val mom:        0.958102       1.854958      -0.896856    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef            -31.071786        17.383331       -13.688455
   rhoval*ves             -3.653971        -6.779704       -10.433675
   psnuc*ves              15.245396      -283.319033      -268.073638
   utot                    5.795712      -145.049369      -139.253656
   rho*exc                -8.460527        -0.395786        -8.856313
   rho*vxc               -11.079953        -0.577512       -11.657465
   valence chg             7.257251        -4.257251         3.000000
   valence mag             1.896856        -0.896856         1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       10.701212  sumtc=        62.965926   ekin=       73.667138
 rhoep=       -8.856313   utot=      -139.253656   ehks=      -74.442831
 mag. mom=     1.000000

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 wgtsmooth=   2.8284271247461901E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.80e-3  last it=2.02e-3
 mixrho: (warning) scr. and lin-mixed densities had 0 and 58981 negative points
 AMIX: Reducing nmix to  1: t_j exceeds tm: tj=-11.12024   0.11003
 AMIX: Reducing nmix to  0: t_j exceeds tm: tj= -7.39079
 AMIX: nmix=0 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel=9.02D-04
 mixrho: warning. negative smrho; isp number min=       1   54997 -0.10499D-03
 mixrho: warning. negative smrho; isp number min=       2   71037 -0.10695D-03
 mixrho: warning. negative smrho; isp number min=       1   54997 -0.10499D-03
 mixrho: warning. negative smrho; isp number min=       1   54997 -0.10499D-03
 mixrho: warning. negative smrho; isp number min=       2   71037 -0.10695D-03
 mixrho: warning. negative smrho; isp number min=       2   71037 -0.10695D-03
 mixrho: add corrections to qcell smrho =  0.43967D-07  0.21984D-10
 unscreened rms difference:  smooth  0.003160   local  0.005910
   screened rms difference:  smooth  0.003122   local  0.005910   tot  0.001805
 mixrho: warning. negative smrho; isp number min=       1   54997 -0.10499D-03
 mixrho: warning. negative smrho; isp number min=       2   71037 -0.10695D-03

 iors  : write restart file (binary, mesh density) 

   it  8  of 10    ehf=       0.491663   ehk=       0.552069
 From last iter    ehf=       0.485101   ehk=       0.552083
 diffe(q)=  0.006562 (0.001805)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999999 ehf=.4916627 ehk=.5520688
 --- BNDFP:  begin iteration 9 of 10
bndfp:start

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.236537   -4.383411     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 jesm=0 return

 smves:: avg es pot at rmt=-0.105085  avg sphere pot= 0.010610  vconst= 0.105085
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.778636   charge     8.383411
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      126034  -1.0694823803267970E-004
 smvxcm: enforce positive smrho_w. Add srshift=   1.0694823804267970E-004
 smooth rhoeps =   -8.749905 (  -5.408753,  -3.341152)
         rhomu =  -11.459368 (  -7.397747,  -4.061621)
       avg vxc =   -0.226720 (  -0.235160,  -0.218281)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       256
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       256
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       256
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       256
 vxcnsp (warning): negative rho: min val =  -2.40E-03
 qqx nkapi nkape=           2           2

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
 qqx nkapi nkape=           2           2
   1      -9.017926   -4.679506  -13.697432    -19.745418  -12.087224  -31.832642

 local terms:     true           smooth         local
 rhoeps:        -8.858197      -8.691206      -0.166991
 rhomu:         -6.485859      -7.352186       0.866327
 spin2:         -5.174452      -4.031052      -1.143400
 total:        -11.660311     -11.383239      -0.277073
 val*vef       -13.697432     -15.152753       1.455321
 val chg:        2.925343       7.308753      -4.383411
 val mom:        1.017498       1.899267      -0.881769    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 qqx nkapi nkape=           2           2
 potential shift to crystal energy zero:    0.000008

 potpus  spin 1 : pnu = 2.921381 2.850000 3.147584 4.102416
 qqx nkapi nkape=           2           2

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -11.898360    7.852023    0.078389   -0.645898
 1      3.000000    1.000000    -5.887832    7.307983    0.124639   -0.608009
 2      3.000000    1.000000     6.000000   28.763451    0.487140   -0.090177
 3      3.000000    1.000000     9.000000   39.654100    0.567182   -0.057513

 potpus  spin 2 : pnu = 2.917618 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -11.331491    7.721264    0.084185   -0.623450
 1      3.000000    1.000000    -5.887832    7.237096    0.131383   -0.579913
 2      3.000000    1.000000     6.000000   29.033941    0.488688   -0.088836
 3      3.000000    1.000000     9.000000   39.802156    0.567852   -0.057169

 Energy terms:             smooth           local           total
   rhoval*vef            -31.859427        18.135174       -13.724253
   rhoval*ves             -3.663617        -6.782352       -10.445969
   psnuc*ves              15.220889      -283.333920      -268.113031
   utot                    5.778636      -145.058136      -139.279500
   rho*exc                -8.749905        -0.166991        -8.916896
   rho*vxc               -11.459368        -0.277073       -11.736440
   valence chg             7.383411        -4.383411         3.000000
   valence mag             1.945980        -0.881769         1.064211
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000
 ... Start MPI k-loop ...
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.4228

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.625250;   3.000000 electrons
         Sum occ. bands:   -2.982511, incl. Bloechl correction: -0.000025
         Mag. moment:       1.000000
Generating TDOS: efermi, and dos window=   -0.6252  -0.5000   0.8748
 bndfp: kpt    1 of    8 k isp=  0.0000  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    1 of    8 k isp=  0.0000  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    2 of    8 k isp=  0.1250  0.1250 -0.1250 1 ndimh nev=    8    5
 bndfp: kpt    2 of    8 k isp=  0.1250  0.1250 -0.1250 2 ndimh nev=    8    5
 bndfp: kpt    3 of    8 k isp=  0.2500  0.2500 -0.2500 1 ndimh nev=    8    5
 bndfp: kpt    3 of    8 k isp=  0.2500  0.2500 -0.2500 2 ndimh nev=    8    5
 bndfp: kpt    4 of    8 k isp=  0.2500  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    4 of    8 k isp=  0.2500  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    5 of    8 k isp=  0.3750  0.1250 -0.1250 1 ndimh nev=    8    5
 bndfp: kpt    5 of    8 k isp=  0.3750  0.1250 -0.1250 2 ndimh nev=    8    5
 bndfp: kpt    6 of    8 k isp=  0.5000  0.2500 -0.2500 1 ndimh nev=    8    5
 bndfp: kpt    6 of    8 k isp=  0.5000  0.2500 -0.2500 2 ndimh nev=    8    5
 bndfp: kpt    7 of    8 k isp=  0.5000  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    7 of    8 k isp=  0.5000  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    8 of    8 k isp=  0.5000  0.2500  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    8 of    8 k isp=  0.5000  0.2500  0.0000 2 ndimh nev=    8    5
 mkrout: site(class) decomposed charge and magnetic moment. class->lmchk

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.903390    7.070304   -4.166914      0.958416    1.863349   -0.904933
       contr. to mm extrapolated for r>rmt:   0.034835 est. true mm = 0.993251
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.31073  sum tc=    31.42692  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.26170  sum tc=    31.53849  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.976822   -1.246867    2.921381    2.922025    2.500000    2.922025
 spn 2 0    0.972487   -1.110120    2.917618    2.918638    2.500000    2.918638
 1     1    0.954081   -0.625522    2.850000    2.906687    2.250000    2.850000
 spn 2 1    0.000000   -1.110165    2.850000    2.185284    2.250000    2.850000
 2     0    0.000001   -0.636096    3.147584    3.144423    3.147584    3.147584
 spn 2 0    0.000000   -1.110113    3.147584    3.111985    3.147584    3.147584
 3     0    0.000000   -0.648169    4.102416    4.098234    4.102416    4.102416
 spn 2 0    0.000000   -0.648169    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -2.982511  val*vef=     -13.724253   sumtv=      10.741742
 sumec=      -40.572433  cor*vef=    -103.540254   ttcor=      62.967821
 rhoeps=      -8.916896     utot=    -139.279500    ehar=     -74.486833
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      126034  -1.0694887972841057E-004
 smvxcm: enforce positive smrho_w. Add srshift=   1.0694887973841057E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      126034  -1.0694759633694883E-004
 smvxcm: enforce positive smrho_w. Add srshift=   1.0694759634694884E-004
 smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=      126034  -1.0694823803267970E-004
 smvxcm: enforce positive smrho_w. Add srshift=   1.0694823804267970E-004

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00   -0.00     0.00    0.00    0.00    -0.00   -0.00   -0.00
 shift forces to make zero average correction:           -0.00   -0.00   -0.00

 srhov:    -30.631941     16.941647    -13.690293 sumev=   -2.982511   sumtv=   10.707782

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.175465   -4.166914     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 jesm=0 return

 smves:: avg es pot at rmt=-0.112025  avg sphere pot= 0.010693  vconst= 0.112025
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.777530   charge     8.166914
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -8.303613 (  -5.161554,  -3.142059)
         rhomu =  -10.873966 (  -7.065392,  -3.808574)
       avg vxc =   -0.158413 (  -0.174703,  -0.142122)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.888098   -4.792684  -13.680783    -19.028772  -11.467369  -30.496141

 local terms:     true           smooth         local
 rhoeps:        -8.841080      -8.287535      -0.553545
 rhomu:         -6.436521      -7.050150       0.613628
 spin2:         -5.201235      -3.802998      -1.398237
 total:        -11.637756     -10.853147      -0.784609
 val*vef       -13.680783     -14.639596       0.958813
 val chg:        2.903390       7.070304      -4.166914
 val mom:        0.958416       1.863349      -0.904933    core:   0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef            -30.508787        16.815324       -13.693463
   rhoval*ves             -3.666285        -6.771380       -10.437665
   psnuc*ves              15.221346      -283.301460      -268.080114
   utot                    5.777530      -145.036420      -139.258890
   rho*exc                -8.303613        -0.553545        -8.857157
   rho*vxc               -10.873966        -0.784609       -11.658575
   valence chg             7.166914        -4.166914         3.000000
   valence mag             1.904933        -0.904933         1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       10.707782  sumtc=        62.965405   ekin=       73.673187
 rhoep=       -8.857157   utot=      -139.258890   ehks=      -74.442861
 mag. mom=     1.000000

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 wgtsmooth=   2.8284271247461901E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.50e-3  last it=1.80e-3
 mixrho: (warning) scr. and lin-mixed densities had 0 and 53905 negative points
 AMIX: Reducing nmix to  1: t_j exceeds tm: tj=  6.06996  -6.02988
 AMIX: nmix=1 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel=7.49D-04
   tj:-3.80834
 mixrho: add corrections to qcell smrho =  0.45733D-07  0.22866D-10
 unscreened rms difference:  smooth  0.002619   local  0.004961
   screened rms difference:  smooth  0.002614   local  0.004961   tot  0.001498
 mixrho: all smrho are positive for isp=  1
 mixrho: all smrho are positive for isp=  2

 iors  : write restart file (binary, mesh density) 

   it  9  of 10    ehf=       0.508067   ehk=       0.552039
 From last iter    ehf=       0.491663   ehk=       0.552069
 diffe(q)=  0.016404 (0.001498)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999999 ehf=.5080672 ehk=.5520394
 --- BNDFP:  begin iteration 10 of 10
bndfp:start

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.089709   -3.862916     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 jesm=0 return

 smves:: avg es pot at rmt=-0.114161  avg sphere pot= 0.010810  vconst= 0.114161
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.751404   charge     7.862916
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -7.763231 (  -4.859755,  -2.903476)
         rhomu =  -10.164310 (  -6.658518,  -3.505791)
       avg vxc =   -0.174927 (  -0.185230,  -0.164625)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 qqx nkapi nkape=           2           2
 qqx nkapi nkape=           2           2

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.712621   -4.964243  -13.676864    -18.013313  -10.601449  -28.614762

 local terms:     true           smooth         local
 rhoeps:        -8.833421      -7.747393      -1.086028
 rhomu:         -6.378339      -6.644196       0.265857
 spin2:         -5.249097      -3.499597      -1.749500
 total:        -11.627435     -10.143793      -1.483643
 val*vef       -13.676864     -13.950343       0.273478
 val chg:        2.897318       6.760235      -3.862916
 val mom:        0.875455       1.812914      -0.937459    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 qqx nkapi nkape=           2           2
 potential shift to crystal energy zero:    0.000008

 potpus  spin 1 : pnu = 2.922025 2.850000 3.147584 4.102416
 qqx nkapi nkape=           2           2

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.000626    7.862262    0.077858   -0.646620
 1      3.000000    1.000000    -5.887832    7.321284    0.124385   -0.608637
 2      3.000000    1.000000     6.000000   28.733313    0.486902   -0.090341
 3      3.000000    1.000000     9.000000   39.629850    0.567042   -0.057573

 potpus  spin 2 : pnu = 2.918638 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -11.480095    7.749813    0.082906   -0.627242
 1      3.000000    1.000000    -5.887832    7.258185    0.130166   -0.584398
 2      3.000000    1.000000     6.000000   28.966368    0.488235   -0.089180
 3      3.000000    1.000000     9.000000   39.757033    0.567616   -0.057277

 Energy terms:             smooth           local           total
   rhoval*vef            -28.624833        14.937865       -13.686968
   rhoval*ves             -3.681976        -6.758341       -10.440317
   psnuc*ves              15.184783      -283.269525      -268.084742
   utot                    5.751404      -145.013933      -139.262530
   rho*exc                -7.763231        -1.086028        -8.849259
   rho*vxc               -10.164310        -1.483643       -11.647952
   valence chg             6.862916        -3.862916         3.000000
   valence mag             1.847295        -0.937459         0.909837
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000
 ... Start MPI k-loop ...
  m_bandcal_init: start
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     1 of     8 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 bndfp: kpt     2 of     8 k=  0.1250  0.1250 -0.1250 ndimh = nmto+napw =     8    8    0
 ... Done MPI k-loop: elapsed time=   1.4349

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.619097;   3.000000 electrons
         Sum occ. bands:   -2.984279, incl. Bloechl correction: -0.000021
         Mag. moment:       1.000000
Generating TDOS: efermi, and dos window=   -0.6191  -0.5000   0.8809
 bndfp: kpt    1 of    8 k isp=  0.0000  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    1 of    8 k isp=  0.0000  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    2 of    8 k isp=  0.1250  0.1250 -0.1250 1 ndimh nev=    8    5
 bndfp: kpt    2 of    8 k isp=  0.1250  0.1250 -0.1250 2 ndimh nev=    8    5
 bndfp: kpt    3 of    8 k isp=  0.2500  0.2500 -0.2500 1 ndimh nev=    8    5
 bndfp: kpt    3 of    8 k isp=  0.2500  0.2500 -0.2500 2 ndimh nev=    8    5
 bndfp: kpt    4 of    8 k isp=  0.2500  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    4 of    8 k isp=  0.2500  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    5 of    8 k isp=  0.3750  0.1250 -0.1250 1 ndimh nev=    8    5
 bndfp: kpt    5 of    8 k isp=  0.3750  0.1250 -0.1250 2 ndimh nev=    8    5
 bndfp: kpt    6 of    8 k isp=  0.5000  0.2500 -0.2500 1 ndimh nev=    8    5
 bndfp: kpt    6 of    8 k isp=  0.5000  0.2500 -0.2500 2 ndimh nev=    8    5
 bndfp: kpt    7 of    8 k isp=  0.5000  0.0000  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    7 of    8 k isp=  0.5000  0.0000  0.0000 2 ndimh nev=    8    5
 bndfp: kpt    8 of    8 k isp=  0.5000  0.2500  0.0000 1 ndimh nev=    8    5
 bndfp: kpt    8 of    8 k isp=  0.5000  0.2500  0.0000 2 ndimh nev=    8    5
 mkrout: site(class) decomposed charge and magnetic moment. class->lmchk

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905933    6.971795   -4.065863      0.960412    1.911738   -0.951326
       contr. to mm extrapolated for r>rmt:   0.033850 est. true mm = 0.994262
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.31527  sum tc=    31.43297  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.27324  sum tc=    31.52999  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..
  ebar: 
    without lo    : ebar = center of gravity of occupied states
    with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.
                    ebar for valence is at the center of gravity of occ. states.
    with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ
                    ebar for valence is at the Fermi energy.
  idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.976525   -1.241798    2.922025    2.924293    2.500000    2.924293
 spn 2 0    0.972760   -1.123200    2.918638    2.921974    2.500000    2.921974
 1     1    0.956647   -0.619276    2.850000    2.909734    2.250000    2.850000
 spn 2 1    0.000000   -1.123235    2.850000    2.183675    2.250000    2.850000
 2     0    0.000001   -0.633246    3.147584    3.144001    3.147584    3.147584
 spn 2 0    0.000000   -1.123198    3.147584    3.111382    3.147584    3.147584
 3     0    0.000000   -0.650135    4.102416    4.097916    4.102416    4.102416
 spn 2 0    0.000000   -0.650135    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -2.984279  val*vef=     -13.686968   sumtv=      10.702689
 sumec=      -40.588519  cor*vef=    -103.555130   ttcor=      62.966611
 rhoeps=      -8.849259     utot=    -139.262530    ehar=     -74.442488
 smvxcm: all smrho_w is positive
 smvxcm: all smrho_w is positive
 smvxcm: all smrho_w is positive

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00    0.00   -0.00    -0.00    0.00   -0.00     0.00   -0.00    0.00
 shift forces to make zero average correction:            0.00   -0.00    0.00

 srhov:    -29.715127     15.988708    -13.726419 sumev=   -2.984279   sumtv=   10.742140

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.146959   -4.065863     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000
 jesm=0 return

 smves:: avg es pot at rmt=-0.111827  avg sphere pot= 0.010673  vconst= 0.111827
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.793410   charge     8.065863
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -8.127514 (  -5.093362,  -3.034152)
         rhomu =  -10.642928 (  -6.979208,  -3.663720)
       avg vxc =   -0.156673 (  -0.172495,  -0.140852)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.896603   -4.809273  -13.705876    -18.776606  -11.085585  -29.862191

 local terms:     true           smooth         local
 rhoeps:        -8.846670      -8.111870      -0.734800
 rhomu:         -6.441760      -6.964449       0.522688
 spin2:         -5.203339      -3.658223      -1.545117
 total:        -11.645100     -10.622671      -1.022429
 val*vef       -13.705876     -14.383444       0.677568
 val chg:        2.905933       6.971795      -4.065863
 val mom:        0.960412       1.911738      -0.951326    core:   0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef            -29.874624        16.156282       -13.718342
   rhoval*ves             -3.641155        -6.815163       -10.456318
   psnuc*ves              15.227975      -283.343008      -268.115033
   utot                    5.793410      -145.079086      -139.285675
   rho*exc                -8.127514        -0.734800        -8.862314
   rho*vxc               -10.642928        -1.022429       -11.665357
   valence chg             7.065863        -4.065863         3.000000
   valence mag             1.951326        -0.951326         1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       10.742140  sumtc=        62.962964   ekin=       73.705103
 rhoep=       -8.862314   utot=      -139.285675   ehks=      -74.442886
 mag. mom=     1.000000

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 wgtsmooth=   2.8284271247461901E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.36e-3  last it=1.50e-3
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.50000  tm= 5.00000  rmsdel=6.78D-04
   tj: 0.03512   0.39677
 mixrho: warning. negative smrho; isp number min=       1   30191 -0.34589D-04
 mixrho: warning. negative smrho; isp number min=       2   58939 -0.35966D-04
 mixrho: warning. negative smrho; isp number min=       1   30191 -0.34589D-04
 mixrho: warning. negative smrho; isp number min=       1   30191 -0.34589D-04
 mixrho: warning. negative smrho; isp number min=       2   58939 -0.35966D-04
 mixrho: warning. negative smrho; isp number min=       2   58939 -0.35966D-04
 mixrho: add corrections to qcell smrho =  0.43134D-07  0.21567D-10
 unscreened rms difference:  smooth  0.002350   local  0.004682
   screened rms difference:  smooth  0.002268   local  0.004682   tot  0.001356
 mixrho: warning. negative smrho; isp number min=       1   30191 -0.34589D-04
 mixrho: warning. negative smrho; isp number min=       2   58939 -0.35966D-04

 iors  : write restart file (binary, mesh density) 

   it 10  of 10    ehf=       0.552412   ehk=       0.552014
 From last iter    ehf=       0.508067   ehk=       0.552039
 diffe(q)=  0.044345 (0.001356)    tol= 0.000010 (0.000500)   more=F
x zbak=1 mmom=.9999999 ehf=.5524119 ehk=.552014
 >>     26.14   exit  lmfp           26.10
OK! end of LMF ======================
