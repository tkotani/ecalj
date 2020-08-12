rdcmd:  lmfa --no-iactiv c -vzbak=0
HOST_INFORMATION platform: gfortran
HOST_INFORMATION compiler version: gcc version 4.3.4 (Ubuntu 4.3.4-10ubuntu1) 
HOST_INFORMATION FFLAGS (<=120): -O3 -fomit-frame-pointer -funroll-loops -ffast-math -ffixed-line-length-132 -DHASIARGC -DHASGETARG -DFDATE -DHASCPUTIME 
HOST_INFORMATION LIBLOC (<=120): /usr/lib64/libfftw3.a /usr/lib64/liblapack.so.3gf /usr/lib64/libblas.a
HOST_INFORMATION uname -a (<=120): Linux takaot61 2.6.32-25-generic #45-Ubuntu SMP Sat Oct 16 19:52:42 UTC 2010 x86_64 GNU/Linux
HOST_INFORMATION /etc/issue: Ubuntu 10.04.1 LTS \n \l
HOST_INFORMATION git branch: refs/heads/master
HOST_INFORMATION git commit: a175f89014e679a92f35536e8d0914a2588d1654
HOST_INFORMATION linked at: Thu Dec 23 00:35:06 JST 2010
 -----------------------  START LMFA     -----------------------
 ptenv() is called with EXT=c
 ptenv() not supported, but continue.
 HEADER sc C atom
 C        xxx            1           1
  mxcst switch =           1           0 F F F
  LMFA  vn 7.00(LMFA 7.0)  verb 30,40,60
 pot:      spin-pol, XC:BH
 end of rdctrl2 in imfav7
 lattic:

                Plat                                  Qlat
   1.000000   1.000000   0.000000        0.500000   0.500000  -0.500000
   1.000000   0.000000   1.000000        0.500000  -0.500000   0.500000
   0.000000   1.000000   1.000000       -0.500000   0.500000   0.500000
  Cell vol= 1000.000000

 LATTC:  as= 2.000   tol=1.00E-08   alat= 7.93701   awald= 0.200
         r1=  3.459   nkd=  79      q1=  2.571   nkg= 137
 goto mksym

 SGROUP: 1 symmetry operations from 0 generators
 SYMLAT: Bravais system is cubic with 48 symmetry operations.
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-5
 GROUPG: the following are sufficient to generate the space group:
         i*r3(-1,1,1) r4z
         i*r3(-1,1,1) r4z
 MKSYM:  found 48 space group operations ... includes inversion
 zzz nclass=           1
 end of mksym x
 goto defspc
 end of defspc
 goto freeat

conf:SPEC_ATOM= C : --- Table for atomic configuration ---
conf int(P)z = int(P) where P is replaced by PZ if it is semicore
conf:  isp  l  int(P) int(P)z    Qval    Qcore   CoreConf
conf:    1  0       2  2        1.000    1.000 => 1,
conf:    1  1       2  2        2.000    0.000 => 
conf:    1  2       3  3        0.000    0.000 => 
conf:    1  3       4  4        0.000    0.000 => 
conf:-----------------------------------------------------
conf:    2  0       2  2        1.000    1.000 => 1,
conf:    2  1       2  2        0.000    0.000 => 
conf:    2  2       3  3        0.000    0.000 => 
conf:    2  3       4  4        0.000    0.000 => 
conf:-----------------------------------------------------

 Species C:  Z=6  Qc=2  R=3.000000  Q=0  mom=2
 mesh:   rmt=3.000000  rmax=19.671121  a=0.02  nr=369  nr(rmax)=463
  Pl=  2.5     2.5     3.5     4.5     spn 2   2.5     2.5     3.5     4.5    
  Ql=  1.0     2.0     0.0     0.0     spn 2   1.0     0.0     0.0     0.0    

  iter     qint         drho          vh0          rho0          vsum     beta
    1    6.000000   5.461E+02       30.0000    0.2984E+02      -12.0633   0.30
   50    6.000000   3.935E-05       29.2312    0.1279E+03      -59.7470   0.30


 sumev=-2.876387  etot=-74.994908  eref=-74.994900  diff= -0.000008

 Optimise free-atom basis for species C, rmt=3
 l  it    Rsm      Eh     stiffR   stiffE      Eval      Exact     Pnu    Ql
 0  41   1.425  -0.888      15.8     35.1   -1.07382  -1.07383    2.91   1.00
 1  31   1.429  -0.336      65.9    157.3   -0.46562  -0.46569    2.89   2.00
 eigenvalue sum:  exact  -2.00520    opt basis  -2.00507    error 0.00013

 Optimise free-atom basis for species C, rmt=3
 l  it    Rsm      Eh     stiffR   stiffE      Eval      Exact     Pnu    Ql
 0  41   1.461  -0.748      17.8     52.8   -0.87118  -0.87119    2.91   1.00
 1  28   1.494  -0.209      80.9    335.8   -0.27747  -0.27759    2.87   0.00
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
 end of freats: spid=C       

  Write mtopara.* ...
 Exit 0 LMFA 
rdcmd:  lmf  --no-iactiv c -vzbak=0
HOST_INFORMATION platform: gfortran
HOST_INFORMATION compiler version: gcc version 4.3.4 (Ubuntu 4.3.4-10ubuntu1) 
HOST_INFORMATION FFLAGS (<=120): -O3 -fomit-frame-pointer -funroll-loops -ffast-math -ffixed-line-length-132 -DHASIARGC -DHASGETARG -DFDATE -DHASCPUTIME 
HOST_INFORMATION LIBLOC (<=120): /usr/lib64/libfftw3.a /usr/lib64/liblapack.so.3gf /usr/lib64/libblas.a
HOST_INFORMATION uname -a (<=120): Linux takaot61 2.6.32-25-generic #45-Ubuntu SMP Sat Oct 16 19:52:42 UTC 2010 x86_64 GNU/Linux
HOST_INFORMATION /etc/issue: Ubuntu 10.04.1 LTS \n \l
HOST_INFORMATION git branch: refs/heads/master
HOST_INFORMATION git commit: a175f89014e679a92f35536e8d0914a2588d1654
HOST_INFORMATION linked at: Thu Dec 23 00:35:06 JST 2010
 -----------------------  START LMF      -----------------------
 ptenv() is called with EXT=c
 ptenv() not supported, but continue.
 HEADER sc C atom
 C        xxx            1           1
  mxcst switch =           1           0 F F F
  LMF  vn 7.00(LMF 7.0)  verb 30,40,60
 special:  forces
 pot:      spin-pol, XC:BH
 bz:       metal(2), tetra, invit, fixed-spin-mom 
 goto setcg
 lattic:

                Plat                                  Qlat
   1.000000   1.000000   0.000000        0.500000   0.500000  -0.500000
   1.000000   0.000000   1.000000        0.500000  -0.500000   0.500000
   0.000000   1.000000   1.000000       -0.500000   0.500000   0.500000
  Cell vol= 1000.000000

 LATTC:  as= 2.000   tol=1.00E-08   alat= 7.93701   awald= 0.200
         r1=  3.459   nkd=  79      q1=  2.571   nkg= 137

 SGROUP: 1 symmetry operations from 0 generators
 SYMLAT: Bravais system is cubic with 48 symmetry operations.
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-5
 GROUPG: the following are sufficient to generate the space group:
         i*r3(-1,1,1) r4z
         i*r3(-1,1,1) r4z
 MKSYM:  found 48 space group operations ... includes inversion
 zzz nclass=           1
 
 lstar xxx=          -2
 BZMESH:  8 irreducible QP from 64 ( 4 4 4 )  shift= F F F
 lstar xxx=          -2

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv  kmxv foca   rfoca
 C        3.000  1.200    3    3         3  0.750  1.500    15    0   1.200

 gvlist: cutoff radius  13.994 gives  45911   recips of max 125000
 SGVSYM: 1207 symmetry stars found for 45911 reciprocal lattice vectors
 

 Makidx:  hamiltonian dimensions Low, Int, High, Negl: 8 0 24 0
 suham :  16 augmentation channels, 16 local potential channels  Maximum lmxa=3

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol= 1.00E-06
 spec      l    rsm    eh     gmax    last term    cutoff
  C        0    1.30  -0.70   5.718    1.05E-06    3143 
  C        1    1.10  -0.20   7.226    1.06E-06    6375 
  C        0    0.80  -1.50   9.292    1.08E-06   13539 
  C        1    0.80  -1.00  10.038    1.00E-06   16961 

 iors  : read restart file (binary, mesh density) 
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

 --- BNDFP:  begin iteration 1 of 10 ---
 ttt nevmx w=           5  4.00000000000000008E-003

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1    0.377522    1.338282     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt= 0.006607  avg sphere pot= 0.019541  vconst=-0.006607
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      2.063910   charge     2.661718
 smooth rhoeps =   -1.494901 (  -1.109320,  -0.385581)
         rhomu =   -1.946550 (  -1.522780,  -0.423770)
       avg vxc =   -0.191348 (  -0.218128,  -0.164568)
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -1.494901 (  -1.109320,  -0.385581)
         rhomu =   -1.946550 (  -1.522780,  -0.423770)
       avg vxc =   -0.191348 (  -0.218128,  -0.164568)

 locpot:
  i job kmax lfltwf=           0           1           3 T

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.341410   -3.941483  -14.282893     -2.908176   -0.943828   -3.852003

 local terms:     true           smooth         local
 rhoeps:        -9.536045      -1.428951      -8.107094
 rhomu:         -7.422765      -1.450795      -5.971970
 spin2:         -5.125221      -0.410320      -4.714901
 total:        -12.547986      -1.861115     -10.686870
 val*vef       -14.282893      -6.925497      -7.357396
 val chg:        3.701843       2.363561       1.338282
 val mom:        1.810990       1.143831       0.667159    core:  -0.000000
 core chg:       2.000000       2.000000      -0.000000
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
   core charge             2.000000        -0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000

 Incompatible or missing qp weights file ...
 Start first of two band passes ...
 end of suham2
 -------- qplist --------
    1   0.000   0.000   0.000
    2   0.125   0.125  -0.125
    3   0.250   0.250  -0.250
    4   0.250   0.000   0.000
    5   0.375   0.125  -0.125
    6   0.500   0.250  -0.250
    7   0.500   0.000   0.000
    8   0.500   0.250   0.000
 sigmamode= F
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.0407 -0.4290 -0.4290 -0.4290  0.1321  0.5284  0.5284  0.5284
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.8387 -0.2374 -0.2374 -0.2374  0.2088  0.6249  0.6249  0.6249
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.430662;   4.000000 electrons
         Sum occ. bands:   -2.743234, incl. Bloechl correction: -0.000179
         Mag. moment:       2.000000

 Saved qp weights ...
 Start second band pass ...
 -------- qplist --------
    1   0.000   0.000   0.000
    2   0.125   0.125  -0.125
    3   0.250   0.250  -0.250
    4   0.250   0.000   0.000
    5   0.375   0.125  -0.125
    6   0.500   0.250  -0.250
    7   0.500   0.000   0.000
    8   0.500   0.250   0.000
 sigmamode= F
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.0407 -0.4290 -0.4290 -0.4290  0.1321  0.5284  0.5284  0.5284
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.8387 -0.2374 -0.2374 -0.2374  0.2088  0.6249  0.6249  0.6249
 Est Ef = -0.431 < evl(4)=-0.429 ... using qval=4.0, revise to -0.4290
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.430662;   4.000000 electrons
         Sum occ. bands:   -2.743234, incl. Bloechl correction: -0.000179
         Mag. moment:       2.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.686857    3.838105   -0.151248      1.796590    2.141140   -0.344551
       contr. to mm extrapolated for r>rmt:   0.163680 est. true mm = 1.960270
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85523  sum tc=    31.38701  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78555  sum tc=    31.54499  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.957935   -1.038744    2.900000    2.914622    2.500000    2.914622
 spn 2 0    0.945133   -0.836887    2.900000    2.908183    2.500000    2.908183
 1     1    1.783776   -0.429110    2.850000    2.889411    2.250000    2.850000
 spn 2 1    0.000000   -1.101867    2.850000    2.167565    2.250000    2.850000
 2     0    0.000011   -0.813692    3.180000    3.132790    3.147584    3.147584
 spn 2 0    0.000000   -1.175532    3.180000    3.108512    3.147584    3.147584
 3     0    0.000001   -0.763293    4.120000    4.096170    4.102416    4.102416
 spn 2 0    0.000000   -1.128828    4.120000    4.085128    4.102416    4.102416

 Harris energy:
 sumev=       -2.743234  val*vef=     -14.362499   sumtv=      11.619265
 sumec=      -39.640777  cor*vef=    -102.572782   ttcor=      62.932005
 rhoeps=      -9.601995     utot=    -139.945263    ehar=     -74.995989
 smvxcm: all smrho_w is positive
 smvxcm: all smrho_w is positive
 smvxcm: all smrho_w is positive

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00   -0.00    0.00    -0.00    0.00   -0.00
 shift forces to make zero average correction:           -0.00    0.00   -0.00

 srhov:     -7.491327     -6.855971    -14.347299 sumev=   -2.743234   sumtv=   11.604065

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.042666   -0.151248     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt= 0.008341  avg sphere pot= 0.014058  vconst=-0.008341
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.531040   charge     4.151248
 smooth rhoeps =   -3.081542 (  -2.415659,  -0.665883)
         rhomu =   -4.025493 (  -3.336154,  -0.689339)
       avg vxc =   -0.206216 (  -0.237894,  -0.174538)
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -3.081542 (  -2.415659,  -0.665883)
         rhomu =   -4.025493 (  -3.336154,  -0.689339)
       avg vxc =   -0.206216 (  -0.237894,  -0.174538)

 locpot:
  i job kmax lfltwf=           0           0           3 F

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

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
 val mom:        1.796590       2.141140      -0.344551    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000

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
 mag. mom=     2.000000  (bands)        2.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
 mixrho: sum smrho  init = 0.249660D+03-0.119196D-25 0.249660D+03       0
 mixrho: sum smrnew new  = 0.405987D+03-0.228641D-16 0.405987D+03       0
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 mixrho: dqsum rmsuns=  0.14895D-02  0.89310D-02  0.15637D-18
 mixrealsmooth= T
 wgtsmooth=  2.82842712474619005E-003
 mixrho:  sought 2 iter from file mixm; read 0.  RMS DQ=8.24e-3
 AMIX: nmix=0 mmix=8  nelts=252052  beta=0.5  tm=5  rmsdel=4.12e-3
 mixrealsmooth= T
 smrho qcell: add correction to smrho=  8.19793182227357420E-008  4.09896591113678777E-011
 unscreened rms difference:  smooth  0.012630   local  0.036279
   screened rms difference:  smooth  0.010354   local  0.036279   tot  0.008236
 mixrho: all smrho is positive for isp=           1
 mixrho: warning. negative smrho; isp number min=           2       37543 -5.54192309110997225E-006

 iors  : write restart file (binary, mesh density) 

   it  1  of 10    ehf=      -0.001089   ehk=      -0.001245
h zbak=0 mmom=1.9999999 ehf=-.0010889 ehk=-.0012455

 --- BNDFP:  begin iteration 2 of 10 ---
 ttt nevmx w=           5  4.00000000000000008E-003

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1    0.167428    0.593517     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt= 0.008356  avg sphere pot= 0.016799  vconst=-0.008356
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      2.743369   charge     3.406483
 smvxcm (warning) mesh density negative at 37543 points:  rhomin=-5.54e-6
 smooth rhoeps =   -2.247142 (  -1.725344,  -0.521798)
         rhomu =   -2.931586 (  -2.377899,  -0.553688)
       avg vxc =   -0.191396 (  -0.225052,  -0.157740)
 smvxcm: negative smrho_w number,min(smrho_w)=       37543 -5.54192309110997225E-006
 smvxcm: enforce positive smrho_w. Add srshift=  5.54192310110997207E-006
 smooth rhoeps =   -2.249370 (  -1.726578,  -0.522793)
         rhomu =   -2.934431 (  -2.379652,  -0.554779)
       avg vxc =   -0.201611 (  -0.230231,  -0.172991)

 locpot:
  i job kmax lfltwf=           0           1           3 T

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

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
 potential shift to crystal energy zero:    0.000003

 potpus  spin 1 : pnu = 2.914622 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.915144    7.418308    0.092820   -0.587642
 1      3.000000    1.000000    -5.887832    7.121430    0.144054   -0.533608
 2      3.000000    1.000000     6.000000   29.630954    0.492412   -0.085937
 3      3.000000    1.000000     9.000000   40.185000    0.569710   -0.056283

 potpus  spin 2 : pnu = 2.908183 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.110313    7.331763    0.102466   -0.559526
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

 Read qp weights ...  ef=-0.430662
 end of suham2
 -------- qplist --------
    1   0.000   0.000   0.000
    2   0.125   0.125  -0.125
    3   0.250   0.250  -0.250
    4   0.250   0.000   0.000
    5   0.375   0.125  -0.125
    6   0.500   0.250  -0.250
    7   0.500   0.000   0.000
    8   0.500   0.250   0.000
 sigmamode= F
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.0413 -0.4295 -0.4295 -0.4295  0.1220  0.5243  0.5243  0.5243
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.8399 -0.2386 -0.2386 -0.2386  0.1936  0.6187  0.6187  0.6187
 Est Ef = -0.431 < evl(4)=-0.430 ... using qval=4.0, revise to -0.4295
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.431288;   4.000000 electrons
         Sum occ. bands:   -2.746491, incl. Bloechl correction: -0.000185
         Mag. moment:       2.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.686196    3.834421   -0.148225      1.795571    2.119153   -0.323582
       contr. to mm extrapolated for r>rmt:   0.164298 est. true mm = 1.959869
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85641  sum tc=    31.38705  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78678  sum tc=    31.54486  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.958094   -1.039238    2.914622    2.914698    2.500000    2.914698
 spn 2 0    0.945312   -0.837837    2.908183    2.908290    2.500000    2.908290
 1     1    1.782776   -0.429972    2.850000    2.889129    2.250000    2.850000
 spn 2 1    0.000000   -1.116964    2.850000    2.165383    2.250000    2.850000
 2     0    0.000012   -0.812935    3.147584    3.132784    3.147584    3.147584
 spn 2 0    0.000000   -1.179203    3.147584    3.108380    3.147584    3.147584
 3     0    0.000001   -0.762080    4.102416    4.096176    4.102416    4.102416
 spn 2 0    0.000000   -1.131551    4.102416    4.085079    4.102416    4.102416

 Harris energy:
 sumev=       -2.746491  val*vef=     -14.362585   sumtv=      11.616094
 sumec=      -39.643186  cor*vef=    -102.575191   ttcor=      62.932005
 rhoeps=      -9.601476     utot=    -139.945009    ehar=     -74.998386
 smvxcm: negative smrho_w number,min(smrho_w)=       37543 -5.54188688412260432E-006
 smvxcm: enforce positive smrho_w. Add srshift=  5.54188689412260414E-006
 smvxcm: negative smrho_w number,min(smrho_w)=       37543 -5.54195929809734018E-006
 smvxcm: enforce positive smrho_w. Add srshift=  5.54195930809734000E-006
 smvxcm: negative smrho_w number,min(smrho_w)=       37543 -5.54192309110997225E-006
 smvxcm: enforce positive smrho_w. Add srshift=  5.54192310110997207E-006

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00   -0.00   -0.00    -0.00   -0.00    0.00     0.00    0.00   -0.00
 shift forces to make zero average correction:            0.00    0.00   -0.00

 srhov:     -8.353510     -5.999383    -14.352893 sumev=   -2.746491   sumtv=   11.606401

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.041813   -0.148225     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt= 0.008178  avg sphere pot= 0.014111  vconst=-0.008178
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.513369   charge     4.148225
 smooth rhoeps =   -3.076559 (  -2.402682,  -0.673876)
         rhomu =   -4.018781 (  -3.318628,  -0.700153)
       avg vxc =   -0.206362 (  -0.238091,  -0.174632)
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -3.076559 (  -2.402682,  -0.673876)
         rhomu =   -4.018781 (  -3.318628,  -0.700153)
       avg vxc =   -0.206362 (  -0.238091,  -0.174632)

 locpot:
  i job kmax lfltwf=           0           0           3 F

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.331796   -3.948938  -14.280734     -7.283005   -1.760723   -9.043728

 local terms:     true           smooth         local
 rhoeps:        -9.525008      -3.009624      -6.515384
 rhomu:         -7.408043      -3.245297      -4.162746
 spin2:         -5.125515      -0.686784      -4.438731
 total:        -12.533558      -3.932081      -8.601477
 val*vef       -14.280734      -8.610255      -5.670479
 val chg:        3.686196       3.834421      -0.148225
 val mom:        1.795571       2.119153      -0.323582    core:   0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef             -9.123231        -5.237007       -14.360238
   rhoval*ves             -4.670978        -5.580258       -10.251236
   psnuc*ves              11.697716      -281.331484      -269.633768
   utot                    3.513369      -143.455871      -139.942502
   rho*exc                -3.076559        -6.515384        -9.591943
   rho*vxc                -4.018781        -8.601477       -12.620259
   valence chg             4.148225        -0.148225         4.000000
   valence mag             2.323582        -0.323582         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       11.606401  sumtc=        62.931909   ekin=       74.538310
 rhoep=       -9.591943   utot=      -139.942502   ehks=      -74.996135
 mag. mom=     2.000000  (bands)        2.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
 mixrho: sum smrho  init = 0.327824D+03 0.582685D-26 0.327824D+03     405
 mixrho: sum smrnew new  = 0.404488D+03-0.199976D-16 0.404488D+03       0
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 mixrho: dqsum rmsuns=  0.74174D-03  0.44238D-02  0.15172D-18
 mixrealsmooth= T
 wgtsmooth=  2.82842712474619005E-003
 mixrho:  sought 2 iter from file mixm; read 1.  RMS DQ=4.07e-3  last it=8.24e-3
 AMIX: nmix=1 mmix=8  nelts=252052  beta=0.5  tm=5  rmsdel=2.03e-3
   tj:-0.97479
 mixrealsmooth= T
 smrho qcell: add correction to smrho=  3.22408738595569844E-009  1.61204369297784949E-012
 unscreened rms difference:  smooth  0.006256   local  0.017824
   screened rms difference:  smooth  0.005199   local  0.017824   tot  0.004069
 mixrho: all smrho is positive for isp=           1
 mixrho: warning. negative smrho; isp number min=           2       44533 -9.41266562620377340E-006

 iors  : write restart file (binary, mesh density) 

   it  2  of 10    ehf=      -0.003486   ehk=      -0.001235
 From last iter    ehf=      -0.001089   ehk=      -0.001245
 diffe(q)= -0.002397 (0.004069)    tol= 0.000010 (0.000500)   more=T
i zbak=0 mmom=1.9999999 ehf=-.0034863 ehk=-.0012345

 --- BNDFP:  begin iteration 3 of 10 ---
 ttt nevmx w=           5  4.00000000000000008E-003

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.039176   -0.138874     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt= 0.009537  avg sphere pot= 0.014145  vconst=-0.009537
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.500418   charge     4.138874
 smvxcm (warning) mesh density negative at 44533 points:  rhomin=-9.41e-6
 smooth rhoeps =   -3.068949 (  -2.395540,  -0.673409)
         rhomu =   -4.008868 (  -3.308656,  -0.700212)
       avg vxc =   -0.196229 (  -0.232460,  -0.159997)
 smvxcm: negative smrho_w number,min(smrho_w)=       44533 -9.41266562620377340E-006
 smvxcm: enforce positive smrho_w. Add srshift=  9.41266563620377322E-006
 smooth rhoeps =   -3.072886 (  -2.397729,  -0.675156)
         rhomu =   -4.013872 (  -3.311765,  -0.702108)
       avg vxc =   -0.210038 (  -0.240011,  -0.180064)

 locpot:
  i job kmax lfltwf=           0           1           3 T

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.332114   -3.951102  -14.283216     -7.251151   -1.756599   -9.007750

 local terms:     true           smooth         local
 rhoeps:        -9.528581      -3.002141      -6.526440
 rhomu:         -7.410744      -3.235504      -4.175239
 spin2:         -5.127466      -0.686727      -4.440739
 total:        -12.538209      -3.922231      -8.615978
 val*vef       -14.283216      -8.602434      -5.680782
 val chg:        3.691284       3.830158      -0.138874
 val mom:        1.795674       2.113144      -0.317470    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000004

 potpus  spin 1 : pnu = 2.914698 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.925427    7.421852    0.092692   -0.588008
 1      3.000000    1.000000    -5.887832    7.123415    0.143913   -0.534049
 2      3.000000    1.000000     6.000000   29.622566    0.492354   -0.085978
 3      3.000000    1.000000     9.000000   40.179010    0.569678   -0.056297

 potpus  spin 2 : pnu = 2.908290 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.122767    7.333209    0.102332   -0.559812
 1      3.000000    1.000000    -5.887832    7.131447    0.153332   -0.500934
 2      3.000000    1.000000     6.000000   29.967250    0.494286   -0.084410
 3      3.000000    1.000000     9.000000   40.375636    0.570571   -0.055857

 Energy terms:             smooth           local           total
   rhoval*vef             -9.088720        -5.275466       -14.364187
   rhoval*ves             -4.673486        -5.576110       -10.249596
   psnuc*ves              11.674323      -281.316998      -269.642675
   utot                    3.500418      -143.446554      -139.946136
   rho*exc                -3.072886        -6.526440        -9.599326
   rho*vxc                -4.013872        -8.615978       -12.629850
   valence chg             4.138874        -0.138874         4.000000
   valence mag             2.317470        -0.317470         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000

 Read qp weights ...  ef=-0.431288
 end of suham2
 -------- qplist --------
    1   0.000   0.000   0.000
    2   0.125   0.125  -0.125
    3   0.250   0.250  -0.250
    4   0.250   0.000   0.000
    5   0.375   0.125  -0.125
    6   0.500   0.250  -0.250
    7   0.500   0.000   0.000
    8   0.500   0.250   0.000
 sigmamode= F
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.0423 -0.4304 -0.4304 -0.4304  0.1175  0.5221  0.5221  0.5221
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.8414 -0.2400 -0.2400 -0.2400  0.1863  0.6155  0.6155  0.6155
 Est Ef = -0.431 < evl(4)=-0.430 ... using qval=4.0, revise to -0.4304
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.432140;   4.000000 electrons
         Sum occ. bands:   -2.750741, incl. Bloechl correction: -0.000190
         Mag. moment:       2.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.686238    3.838429   -0.152191      1.795208    2.105238   -0.310031
       contr. to mm extrapolated for r>rmt:   0.164440 est. true mm = 1.959648
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85808  sum tc=    31.38702  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78851  sum tc=    31.54457  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.958242   -1.040109    2.914698    2.914772    2.500000    2.914772
 spn 2 0    0.945515   -0.839206    2.908290    2.908400    2.500000    2.908400
 1     1    1.782468   -0.431027    2.850000    2.888983    2.250000    2.850000
 spn 2 1    0.000000   -1.123875    2.850000    2.164446    2.250000    2.850000
 2     0    0.000012   -0.811887    3.147584    3.132792    3.147584    3.147584
 spn 2 0    0.000000   -1.179356    3.147584    3.108378    3.147584    3.147584
 3     0    0.000001   -0.760727    4.102416    4.096182    4.102416    4.102416
 spn 2 0    0.000000   -1.132304    4.102416    4.085064    4.102416    4.102416

 Harris energy:
 sumev=       -2.750741  val*vef=     -14.364187   sumtv=      11.613446
 sumec=      -39.646594  cor*vef=    -102.578551   ttcor=      62.931957
 rhoeps=      -9.599326     utot=    -139.946136    ehar=     -75.000059
 smvxcm: negative smrho_w number,min(smrho_w)=       44533 -9.41268045838445528E-006
 smvxcm: enforce positive smrho_w. Add srshift=  9.41268046838445511E-006
 smvxcm: negative smrho_w number,min(smrho_w)=       44533 -9.41265079402309321E-006
 smvxcm: enforce positive smrho_w. Add srshift=  9.41265080402309303E-006
 smvxcm: negative smrho_w number,min(smrho_w)=       44533 -9.41266562620377509E-006
 smvxcm: enforce positive smrho_w. Add srshift=  9.41266563620377492E-006

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00   -0.00   -0.00     0.00   -0.00   -0.00    -0.00    0.00    0.00
 shift forces to make zero average correction:           -0.00    0.00    0.00

 srhov:     -9.137640     -5.226323    -14.363963 sumev=   -2.750741   sumtv=   11.613222

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.042932   -0.152191     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt= 0.008107  avg sphere pot= 0.014138  vconst=-0.008107
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.502549   charge     4.152191
 smooth rhoeps =   -3.081163 (  -2.398838,  -0.682324)
         rhomu =   -4.024677 (  -3.313691,  -0.710986)
       avg vxc =   -0.206393 (  -0.238142,  -0.174645)
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -3.081163 (  -2.398838,  -0.682324)
         rhomu =   -4.024677 (  -3.313691,  -0.710986)
       avg vxc =   -0.206393 (  -0.238142,  -0.174645)

 locpot:
  i job kmax lfltwf=           0           0           3 F

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.334362   -3.951297  -14.285660     -7.273304   -1.791791   -9.065095

 local terms:     true           smooth         local
 rhoeps:        -9.525840      -3.014264      -6.511576
 rhomu:         -7.408634      -3.240351      -4.168283
 spin2:         -5.126020      -0.697672      -4.428348
 total:        -12.534654      -3.938023      -8.596631
 val*vef       -14.285660      -8.620675      -5.664984
 val chg:        3.686238       3.838429      -0.152191
 val mom:        1.795208       2.105238      -0.310031    core:   0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef             -9.144544        -5.220566       -14.365110
   rhoval*ves             -4.675449        -5.579708       -10.255157
   psnuc*ves              11.680546      -281.321777      -269.641231
   utot                    3.502549      -143.450743      -139.948194
   rho*exc                -3.081163        -6.511576        -9.592739
   rho*vxc                -4.024677        -8.596631       -12.621308
   valence chg             4.152191        -0.152191         4.000000
   valence mag             2.310031        -0.310031         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       11.613222  sumtc=        62.931591   ekin=       74.544813
 rhoep=       -9.592739   utot=      -139.948194   ehks=      -74.996120
 mag. mom=     2.000000  (bands)        2.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
 mixrho: sum smrho  init = 0.403522D+03-0.477894D-26 0.403522D+03     863
 mixrho: sum smrnew new  = 0.403889D+03 0.139170D-16 0.403889D+03       0
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 mixrho: dqsum rmsuns=  0.13317D-04  0.11828D-03  0.34570D-19
 mixrealsmooth= T
 wgtsmooth=  2.82842712474619005E-003
 mixrho:  sought 2 iter from file mixm; read 2.  RMS DQ=9.77e-5  last it=4.07e-3
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.5  tm=5  rmsdel=4.89e-5
   tj:-1.56876   0.76419
 mixrealsmooth= T
 smrho qcell: add correction to smrho=  7.14611994578895349E-008  3.57305997289447778E-011
 unscreened rms difference:  smooth  0.000167   local  0.000345
   screened rms difference:  smooth  0.000160   local  0.000345   tot  0.000098
 mixrho: all smrho is positive for isp=           1
 mixrho: warning. negative smrho; isp number min=           2       32143 -7.06353618157365483E-006

 iors  : write restart file (binary, mesh density) 

   it  3  of 10    ehf=      -0.005159   ehk=      -0.001220
 From last iter    ehf=      -0.003486   ehk=      -0.001235
 diffe(q)= -0.001672 (0.000098)    tol= 0.000010 (0.000500)   more=T
i zbak=0 mmom=1.9999999 ehf=-.0051587 ehk=-.0012195

 --- BNDFP:  begin iteration 4 of 10 ---
 ttt nevmx w=           5  4.00000000000000008E-003

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.044666   -0.158338     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt= 0.008865  avg sphere pot= 0.014112  vconst=-0.008865
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.509498   charge     4.158338
 smvxcm (warning) mesh density negative at 32143 points:  rhomin=-7.06e-6
 smooth rhoeps =   -3.090103 (  -2.406413,  -0.683690)
         rhomu =   -4.036438 (  -3.324142,  -0.712296)
       avg vxc =   -0.200143 (  -0.235499,  -0.164787)
 smvxcm: negative smrho_w number,min(smrho_w)=       32143 -7.06353618157365483E-006
 smvxcm: enforce positive smrho_w. Add srshift=  7.06353619157365465E-006
 smooth rhoeps =   -3.093054 (  -2.408040,  -0.685014)
         rhomu =   -4.040224 (  -3.326468,  -0.713756)
       avg vxc =   -0.210519 (  -0.240300,  -0.180739)

 locpot:
  i job kmax lfltwf=           0           1           3 T

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.334206   -3.952230  -14.286436     -7.295881   -1.793864   -9.089745

 local terms:     true           smooth         local
 rhoeps:        -9.527503      -3.023215      -6.504288
 rhomu:         -7.409837      -3.250851      -4.158986
 spin2:         -5.126980      -0.698906      -4.428074
 total:        -12.536817      -3.949757      -8.587060
 val*vef       -14.286436      -8.627029      -5.659407
 val chg:        3.688745       3.847083      -0.158338
 val mom:        1.795171       2.110797      -0.315627    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000005

 potpus  spin 1 : pnu = 2.914772 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.935323    7.422137    0.092626   -0.588098
 1      3.000000    1.000000    -5.887832    7.124156    0.143885   -0.534124
 2      3.000000    1.000000     6.000000   29.620271    0.492339   -0.085989
 3      3.000000    1.000000     9.000000   40.177358    0.569670   -0.056301

 potpus  spin 2 : pnu = 2.908400 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.135568    7.333166    0.102235   -0.559934
 1      3.000000    1.000000    -5.887832    7.132234    0.153295   -0.501025
 2      3.000000    1.000000     6.000000   29.964362    0.494266   -0.084423
 3      3.000000    1.000000     9.000000   40.373502    0.570560   -0.055862

 Energy terms:             smooth           local           total
   rhoval*vef             -9.170281        -5.196692       -14.366973
   rhoval*ves             -4.670316        -5.583726       -10.254042
   psnuc*ves              11.689311      -281.334400      -269.645089
   utot                    3.509498      -143.459063      -139.949565
   rho*exc                -3.093054        -6.504288        -9.597341
   rho*vxc                -4.040224        -8.587060       -12.627284
   valence chg             4.158338        -0.158338         4.000000
   valence mag             2.315626        -0.315627         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.43214
 end of suham2
 -------- qplist --------
    1   0.000   0.000   0.000
    2   0.125   0.125  -0.125
    3   0.250   0.250  -0.250
    4   0.250   0.000   0.000
    5   0.375   0.125  -0.125
    6   0.500   0.250  -0.250
    7   0.500   0.000   0.000
    8   0.500   0.250   0.000
 sigmamode= F
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.0424 -0.4304 -0.4304 -0.4304  0.1203  0.5235  0.5235  0.5235
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.8414 -0.2399 -0.2399 -0.2399  0.1916  0.6179  0.6179  0.6179
 Est Ef = -0.432 < evl(4)=-0.430 ... using qval=4.0, revise to -0.4304
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.432167;   4.000000 electrons
         Sum occ. bands:   -2.750914, incl. Bloechl correction: -0.000189
         Mag. moment:       2.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.686695    3.842349   -0.155654      1.795598    2.108601   -0.313003
       contr. to mm extrapolated for r>rmt:   0.164138 est. true mm = 1.959736
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85791  sum tc=    31.38678  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78831  sum tc=    31.54436  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.958254   -1.040287    2.914772    2.914770    2.500000    2.914770
 spn 2 0    0.945548   -0.839343    2.908400    2.908408    2.500000    2.908408
 1     1    1.782880   -0.431080    2.850000    2.889061    2.250000    2.850000
 spn 2 1    0.000000   -1.123162    2.850000    2.164513    2.250000    2.850000
 2     0    0.000012   -0.811700    3.147584    3.132787    3.147584    3.147584
 spn 2 0    0.000000   -1.178997    3.147584    3.108375    3.147584    3.147584
 3     0    0.000001   -0.760545    4.102416    4.096178    4.102416    4.102416
 spn 2 0    0.000000   -1.131767    4.102416    4.085065    4.102416    4.102416

 Harris energy:
 sumev=       -2.750914  val*vef=     -14.366973   sumtv=      11.616059
 sumec=      -39.646219  cor*vef=    -102.577993   ttcor=      62.931774
 rhoeps=      -9.597341     utot=    -139.949565    ehar=     -74.999073
 smvxcm: negative smrho_w number,min(smrho_w)=       32143 -7.06353614326287047E-006
 smvxcm: enforce positive smrho_w. Add srshift=  7.06353615326287029E-006
 smvxcm: negative smrho_w number,min(smrho_w)=       32143 -7.06353621988443919E-006
 smvxcm: enforce positive smrho_w. Add srshift=  7.06353622988443901E-006
 smvxcm: negative smrho_w number,min(smrho_w)=       32143 -7.06353618157365483E-006
 smvxcm: enforce positive smrho_w. Add srshift=  7.06353619157365465E-006

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00   -0.00   -0.00    -0.00   -0.00   -0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:     -9.162362     -5.204836    -14.367198 sumev=   -2.750914   sumtv=   11.616284

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.043909   -0.155654     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt= 0.008157  avg sphere pot= 0.014128  vconst=-0.008157
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.505437   charge     4.155654
 smooth rhoeps =   -3.085961 (  -2.403040,  -0.682922)
         rhomu =   -4.030979 (  -3.319532,  -0.711446)
       avg vxc =   -0.206322 (  -0.238050,  -0.174593)
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -3.085961 (  -2.403040,  -0.682922)
         rhomu =   -4.030979 (  -3.319532,  -0.711446)
       avg vxc =   -0.206322 (  -0.238050,  -0.174593)

 locpot:
  i job kmax lfltwf=           0           0           3 F

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.336029   -3.951402  -14.287431     -7.285830   -1.793339   -9.079169

 local terms:     true           smooth         local
 rhoeps:        -9.526586      -3.019175      -6.507411
 rhomu:         -7.409525      -3.246327      -4.163198
 spin2:         -5.126109      -0.698143      -4.427966
 total:        -12.535634      -3.944470      -8.591164
 val*vef       -14.287431      -8.624513      -5.662918
 val chg:        3.686695       3.842349      -0.155654
 val mom:        1.795598       2.108601      -0.313003    core:   0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef             -9.158498        -5.208263       -14.366761
   rhoval*ves             -4.672865        -5.583121       -10.255985
   psnuc*ves              11.683739      -281.328078      -269.644338
   utot                    3.505437      -143.455599      -139.950162
   rho*exc                -3.085961        -6.507411        -9.593372
   rho*vxc                -4.030979        -8.591164       -12.622142
   valence chg             4.155654        -0.155654         4.000000
   valence mag             2.313003        -0.313003         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       11.616284  sumtc=        62.931133   ekin=       74.547417
 rhoep=       -9.593372   utot=      -139.950162   ehks=      -74.996118
 mag. mom=     2.000000  (bands)        2.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
 mixrho: sum smrho  init = 0.404623D+03 0.124526D-25 0.404623D+03       0
 mixrho: sum smrnew new  = 0.404291D+03 0.189437D-16 0.404291D+03       0
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 mixrho: dqsum rmsuns= -0.26840D-05  0.19413D-04  0.10216D-18
 mixrealsmooth= T
 wgtsmooth=  2.82842712474619005E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.63e-5  last it=9.77e-5
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.5  tm=5  rmsdel=8.13e-6
   tj:-0.04190   0.00413
 mixrealsmooth= T
 smrho qcell: add correction to smrho=  4.16627998767182817E-008  2.08313999383591448E-011
 unscreened rms difference:  smooth  0.000027   local  0.000071
   screened rms difference:  smooth  0.000019   local  0.000071   tot  0.000016
 mixrho: all smrho is positive for isp=           1
 mixrho: warning. negative smrho; isp number min=           2       25103 -5.49172626856645950E-006

 iors  : write restart file (binary, mesh density) 

   it  4  of 10    ehf=      -0.004173   ehk=      -0.001218
 From last iter    ehf=      -0.005159   ehk=      -0.001220
 diffe(q)=  0.000985 (0.000016)    tol= 0.000010 (0.000500)   more=T
i zbak=0 mmom=1.9999999 ehf=-.0041735 ehk=-.0012176

 --- BNDFP:  begin iteration 5 of 10 ---
 ttt nevmx w=           5  4.00000000000000008E-003

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.043980   -0.155906     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt= 0.008648  avg sphere pot= 0.014124  vconst=-0.008648
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.505817   charge     4.155906
 smvxcm (warning) mesh density negative at 25103 points:  rhomin=-5.49e-6
 smooth rhoeps =   -3.087071 (  -2.403754,  -0.683318)
         rhomu =   -4.032444 (  -3.320480,  -0.711964)
       avg vxc =   -0.201659 (  -0.236329,  -0.166988)
 smvxcm: negative smrho_w number,min(smrho_w)=       25103 -5.49172626856645950E-006
 smvxcm: enforce positive smrho_w. Add srshift=  5.49172627856645933E-006
 smooth rhoeps =   -3.089360 (  -2.405011,  -0.684348)
         rhomu =   -4.035391 (  -3.322285,  -0.713106)
       avg vxc =   -0.209903 (  -0.239890,  -0.179916)

 locpot:
  i job kmax lfltwf=           0           1           3 T

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.335217   -3.952022  -14.287239     -7.286837   -1.793465   -9.080302

 local terms:     true           smooth         local
 rhoeps:        -9.527374      -3.020226      -6.507148
 rhomu:         -7.409938      -3.247234      -4.162704
 spin2:         -5.126717      -0.698603      -4.428113
 total:        -12.536654      -3.945837      -8.590817
 val*vef       -14.287239      -8.624911      -5.662328
 val chg:        3.688178       3.844084      -0.155906
 val mom:        1.795400       2.108774      -0.313374    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000005

 potpus  spin 1 : pnu = 2.914770 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.935134    7.422442    0.092625   -0.588101
 1      3.000000    1.000000    -5.887832    7.124401    0.143881   -0.534127
 2      3.000000    1.000000     6.000000   29.619726    0.492335   -0.085991
 3      3.000000    1.000000     9.000000   40.176920    0.569667   -0.056302

 potpus  spin 2 : pnu = 2.908408 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.136510    7.333384    0.102227   -0.559941
 1      3.000000    1.000000    -5.887832    7.132505    0.153291   -0.501029
 2      3.000000    1.000000     6.000000   29.963731    0.494261   -0.084426
 3      3.000000    1.000000     9.000000   40.372987    0.570557   -0.055863

 Energy terms:             smooth           local           total
   rhoval*vef             -9.160529        -5.206938       -14.367467
   rhoval*ves             -4.672053        -5.582890       -10.254943
   psnuc*ves              11.683687      -281.329289      -269.645602
   utot                    3.505817      -143.456089      -139.950273
   rho*exc                -3.089360        -6.507148        -9.596508
   rho*vxc                -4.035391        -8.590817       -12.626208
   valence chg             4.155906        -0.155906         4.000000
   valence mag             2.313374        -0.313374         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.432167
 end of suham2
 -------- qplist --------
    1   0.000   0.000   0.000
    2   0.125   0.125  -0.125
    3   0.250   0.250  -0.250
    4   0.250   0.000   0.000
    5   0.375   0.125  -0.125
    6   0.500   0.250  -0.250
    7   0.500   0.000   0.000
    8   0.500   0.250   0.000
 sigmamode= F
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.0423 -0.4303 -0.4303 -0.4303  0.1223  0.5244  0.5244  0.5244
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.8413 -0.2397 -0.2397 -0.2397  0.1952  0.6195  0.6195  0.6195
 Est Ef = -0.432 < evl(4)=-0.430 ... using qval=4.0, revise to -0.4303
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.432082;   4.000000 electrons
         Sum occ. bands:   -2.750551, incl. Bloechl correction: -0.000188
         Mag. moment:       2.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.686892    3.843208   -0.156316      1.795804    2.110796   -0.314992
       contr. to mm extrapolated for r>rmt:   0.163986 est. true mm = 1.959789
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85773  sum tc=    31.38672  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78812  sum tc=    31.54433  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.958245   -1.040279    2.914770    2.914761    2.500000    2.914761
 spn 2 0    0.945544   -0.839285    2.908408    2.908402    2.500000    2.908402
 1     1    1.783090   -0.430998    2.850000    2.889107    2.250000    2.850000
 spn 2 1    0.000000   -1.121472    2.850000    2.164740    2.250000    2.850000
 2     0    0.000012   -0.811647    3.147584    3.132784    3.147584    3.147584
 spn 2 0    0.000000   -1.178899    3.147584    3.108373    3.147584    3.147584
 3     0    0.000001   -0.760508    4.102416    4.096176    4.102416    4.102416
 spn 2 0    0.000000   -1.131440    4.102416    4.085068    4.102416    4.102416

 Harris energy:
 sumev=       -2.750551  val*vef=     -14.367467   sumtv=      11.616917
 sumec=      -39.645847  cor*vef=    -102.577301   ttcor=      62.931454
 rhoeps=      -9.596508     utot=    -139.950273    ehar=     -74.998410
 smvxcm: negative smrho_w number,min(smrho_w)=       25103 -5.49171687385401163E-006
 smvxcm: enforce positive smrho_w. Add srshift=  5.49171688385401146E-006
 smvxcm: negative smrho_w number,min(smrho_w)=       25103 -5.49173566327890652E-006
 smvxcm: enforce positive smrho_w. Add srshift=  5.49173567327890634E-006
 smvxcm: negative smrho_w number,min(smrho_w)=       25103 -5.49172626856645865E-006
 smvxcm: enforce positive smrho_w. Add srshift=  5.49172627856645848E-006

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00   -0.00   -0.00    -0.00   -0.00   -0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:     -9.160157     -5.207459    -14.367616 sumev=   -2.750551   sumtv=   11.617065

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.044096   -0.156316     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt= 0.008184  avg sphere pot= 0.014121  vconst=-0.008184
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.507319   charge     4.156316
 smooth rhoeps =   -3.087001 (  -2.404606,  -0.682394)
         rhomu =   -4.032356 (  -3.321679,  -0.710677)
       avg vxc =   -0.206286 (  -0.238004,  -0.174568)
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -3.087001 (  -2.404606,  -0.682394)
         rhomu =   -4.032356 (  -3.321679,  -0.710677)
       avg vxc =   -0.206286 (  -0.238004,  -0.174568)

 locpot:
  i job kmax lfltwf=           0           0           3 F

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.336603   -3.951192  -14.287796     -7.290060   -1.791191   -9.081251

 local terms:     true           smooth         local
 rhoeps:        -9.526852      -3.020261      -6.506592
 rhomu:         -7.409884      -3.248535      -4.161349
 spin2:         -5.126100      -0.697374      -4.428727
 total:        -12.535984      -3.945908      -8.590076
 val*vef       -14.287796      -8.624568      -5.663227
 val chg:        3.686892       3.843208      -0.156316
 val mom:        1.795804       2.110796      -0.314992    core:   0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef             -9.160532        -5.206545       -14.367077
   rhoval*ves             -4.671494        -5.584522       -10.256015
   psnuc*ves              11.686131      -281.331391      -269.645260
   utot                    3.507319      -143.457956      -139.950638
   rho*exc                -3.087001        -6.506592        -9.593592
   rho*vxc                -4.032356        -8.590076       -12.622432
   valence chg             4.156316        -0.156316         4.000000
   valence mag             2.314992        -0.314992         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       11.617065  sumtc=        62.931047   ekin=       74.548112
 rhoep=       -9.593592   utot=      -139.950638   ehks=      -74.996118
 mag. mom=     2.000000  (bands)        2.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
 mixrho: sum smrho  init = 0.404330D+03-0.714784D-26 0.404330D+03       0
 mixrho: sum smrnew new  = 0.404457D+03-0.509573D-17 0.404457D+03       0
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 mixrho: dqsum rmsuns=  0.40944D-06  0.15523D-04 -0.55686D-19
 mixrealsmooth= T
 wgtsmooth=  2.82842712474619005E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.21e-5  last it=1.63e-5
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.5  tm=5  rmsdel=6.04e-6
   tj: 0.36431   0.04617
 mixrealsmooth= T
 smrho qcell: add correction to smrho=  4.15057974101351590E-008  2.07528987050675853E-011
 unscreened rms difference:  smooth  0.000022   local  0.000052
   screened rms difference:  smooth  0.000014   local  0.000052   tot  0.000012
 mixrho: all smrho is positive for isp=           1
 mixrho: warning. negative smrho; isp number min=           2       22091 -4.92033476710347751E-006

 iors  : write restart file (binary, mesh density) 

   it  5  of 10    ehf=      -0.003510   ehk=      -0.001218
 From last iter    ehf=      -0.004173   ehk=      -0.001218
 diffe(q)=  0.000663 (0.000012)    tol= 0.000010 (0.000500)   more=T
i zbak=0 mmom=1.9999999 ehf=-.0035103 ehk=-.0012177

 --- BNDFP:  begin iteration 6 of 10 ---
 ttt nevmx w=           5  4.00000000000000008E-003

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.043991   -0.155945     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt= 0.008603  avg sphere pot= 0.014123  vconst=-0.008603
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.506456   charge     4.155945
 smvxcm (warning) mesh density negative at 22091 points:  rhomin=-4.92e-6
 smooth rhoeps =   -3.087141 (  -2.404220,  -0.682921)
         rhomu =   -4.032540 (  -3.321116,  -0.711425)
       avg vxc =   -0.202107 (  -0.236535,  -0.167679)
 smvxcm: negative smrho_w number,min(smrho_w)=       22091 -4.92033476710347751E-006
 smvxcm: enforce positive smrho_w. Add srshift=  4.92033477710347733E-006
 smooth rhoeps =   -3.089189 (  -2.405345,  -0.683845)
         rhomu =   -4.035180 (  -3.322731,  -0.712449)
       avg vxc =   -0.209571 (  -0.239692,  -0.179450)

 locpot:
  i job kmax lfltwf=           0           1           3 T

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.335521   -3.951841  -14.287362     -7.288074   -1.792100   -9.080174

 local terms:     true           smooth         local
 rhoeps:        -9.527378      -3.020320      -6.507058
 rhomu:         -7.410030      -3.247899      -4.162131
 spin2:         -5.126631      -0.698071      -4.428560
 total:        -12.536661      -3.945969      -8.590692
 val*vef       -14.287362      -8.624637      -5.662725
 val chg:        3.688088       3.844033      -0.155945
 val mom:        1.795515       2.109726      -0.314211    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000005

 potpus  spin 1 : pnu = 2.914761 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.933864    7.422618    0.092632   -0.588095
 1      3.000000    1.000000    -5.887832    7.124464    0.143881   -0.534124
 2      3.000000    1.000000     6.000000   29.619617    0.492334   -0.085992
 3      3.000000    1.000000     9.000000   40.176821    0.569667   -0.056302

 potpus  spin 2 : pnu = 2.908402 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.135839    7.333510    0.102231   -0.559936
 1      3.000000    1.000000    -5.887832    7.132572    0.153291   -0.501024
 2      3.000000    1.000000     6.000000   29.963644    0.494261   -0.084427
 3      3.000000    1.000000     9.000000   40.372899    0.570557   -0.055863

 Energy terms:             smooth           local           total
   rhoval*vef             -9.160288        -5.207188       -14.367477
   rhoval*ves             -4.671634        -5.583402       -10.255035
   psnuc*ves              11.684545      -281.330049      -269.645504
   utot                    3.506456      -143.456725      -139.950270
   rho*exc                -3.089189        -6.507058        -9.596247
   rho*vxc                -4.035180        -8.590692       -12.625872
   valence chg             4.155945        -0.155945         4.000000
   valence mag             2.314211        -0.314211         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000

 Read qp weights ...  ef=-0.432082
 end of suham2
 -------- qplist --------
    1   0.000   0.000   0.000
    2   0.125   0.125  -0.125
    3   0.250   0.250  -0.250
    4   0.250   0.000   0.000
    5   0.375   0.125  -0.125
    6   0.500   0.250  -0.250
    7   0.500   0.000   0.000
    8   0.500   0.250   0.000
 sigmamode= F
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.0423 -0.4303 -0.4303 -0.4303  0.1230  0.5247  0.5247  0.5247
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.8412 -0.2396 -0.2396 -0.2396  0.1966  0.6200  0.6200  0.6200
 Est Ef = -0.432 < evl(4)=-0.430 ... using qval=4.0, revise to -0.4303
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.432042;   4.000000 electrons
         Sum occ. bands:   -2.750379, incl. Bloechl correction: -0.000188
         Mag. moment:       2.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.686954    3.843342   -0.156388      1.795873    2.111553   -0.315680
       contr. to mm extrapolated for r>rmt:   0.163935 est. true mm = 1.959808
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85767  sum tc=    31.38671  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78805  sum tc=    31.54434  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.958241   -1.040266    2.914761    2.914756    2.500000    2.914756
 spn 2 0    0.945540   -0.839252    2.908402    2.908399    2.500000    2.908399
 1     1    1.783160   -0.430959    2.850000    2.889123    2.250000    2.850000
 spn 2 1    0.000000   -1.120613    2.850000    2.164860    2.250000    2.850000
 2     0    0.000012   -0.811634    3.147584    3.132783    3.147584    3.147584
 spn 2 0    0.000000   -1.178886    3.147584    3.108372    3.147584    3.147584
 3     0    0.000001   -0.760501    4.102416    4.096176    4.102416    4.102416
 spn 2 0    0.000000   -1.131327    4.102416    4.085069    4.102416    4.102416

 Harris energy:
 sumev=       -2.750379  val*vef=     -14.367477   sumtv=      11.617097
 sumec=      -39.645723  cor*vef=    -102.576974   ttcor=      62.931250
 rhoeps=      -9.596247     utot=    -139.950270    ehar=     -74.998170
 smvxcm: negative smrho_w number,min(smrho_w)=       22091 -4.92032197224494394E-006
 smvxcm: enforce positive smrho_w. Add srshift=  4.92032198224494376E-006
 smvxcm: negative smrho_w number,min(smrho_w)=       22091 -4.92034756196201108E-006
 smvxcm: enforce positive smrho_w. Add srshift=  4.92034757196201091E-006
 smvxcm: negative smrho_w number,min(smrho_w)=       22091 -4.92033476710347751E-006
 smvxcm: enforce positive smrho_w. Add srshift=  4.92033477710347733E-006

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00   -0.00   -0.00     0.00    0.00   -0.00    -0.00   -0.00    0.00
 shift forces to make zero average correction:           -0.00   -0.00    0.00

 srhov:     -9.160072     -5.207549    -14.367621 sumev=   -2.750379   sumtv=   11.617241

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.044116   -0.156388     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt= 0.008193  avg sphere pot= 0.014119  vconst=-0.008193
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.507985   charge     4.156388
 smooth rhoeps =   -3.087155 (  -2.405031,  -0.682124)
         rhomu =   -4.032565 (  -3.322256,  -0.710309)
       avg vxc =   -0.206274 (  -0.237988,  -0.174559)
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -3.087155 (  -2.405031,  -0.682124)
         rhomu =   -4.032565 (  -3.322256,  -0.710309)
       avg vxc =   -0.206274 (  -0.237988,  -0.174559)

 locpot:
  i job kmax lfltwf=           0           0           3 F

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.336770   -3.951089  -14.287859     -7.291115   -1.790146   -9.081261

 local terms:     true           smooth         local
 rhoeps:        -9.526930      -3.020430      -6.506500
 rhomu:         -7.409994      -3.249130      -4.160864
 spin2:         -5.126091      -0.697005      -4.429087
 total:        -12.536086      -3.946135      -8.589951
 val*vef       -14.287859      -8.624333      -5.663527
 val chg:        3.686954       3.843342      -0.156388
 val mom:        1.795873       2.111553      -0.315680    core:   0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef             -9.160527        -5.206599       -14.367126
   rhoval*ves             -4.671035        -5.584949       -10.255983
   psnuc*ves              11.687004      -281.332527      -269.645524
   utot                    3.507985      -143.458738      -139.950754
   rho*exc                -3.087155        -6.506500        -9.593655
   rho*vxc                -4.032565        -8.589951       -12.622515
   valence chg             4.156388        -0.156388         4.000000
   valence mag             2.315680        -0.315680         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       11.617241  sumtc=        62.931050   ekin=       74.548291
 rhoep=       -9.593655   utot=      -139.950754   ehks=      -74.996118
 mag. mom=     2.000000  (bands)        2.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
 mixrho: sum smrho  init = 0.404385D+03-0.968499D-27 0.404385D+03       0
 mixrho: sum smrnew new  = 0.404504D+03-0.127402D-16 0.404504D+03       0
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 mixrho: dqsum rmsuns=  0.44318D-06  0.14154D-04 -0.66892D-19
 mixrealsmooth= T
 wgtsmooth=  2.82842712474619005E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.11e-5  last it=1.21e-5
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: Reducing nmix to  0: t_j exceeds tm: tj=-10.08168
 AMIX: nmix=0 mmix=8  nelts=252052  beta=0.5  tm=5  rmsdel=5.54e-6
 mixrealsmooth= T
 smrho qcell: add correction to smrho=  4.18402003599105399E-008  2.09201001799552749E-011
 unscreened rms difference:  smooth  0.000020   local  0.000048
   screened rms difference:  smooth  0.000013   local  0.000048   tot  0.000011
 mixrho: all smrho is positive for isp=           1
 mixrho: warning. negative smrho; isp number min=           2       16247 -3.80813807009813665E-006

 iors  : write restart file (binary, mesh density) 

   it  6  of 10    ehf=      -0.003270   ehk=      -0.001218
 From last iter    ehf=      -0.003510   ehk=      -0.001218
 diffe(q)=  0.000241 (0.000011)    tol= 0.000010 (0.000500)   more=T
i zbak=0 mmom=1.9999999 ehf=-.0032696 ehk=-.0012179

 --- BNDFP:  begin iteration 7 of 10 ---
 ttt nevmx w=           5  4.00000000000000008E-003

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.044054   -0.156167     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt= 0.008490  avg sphere pot= 0.014121  vconst=-0.008490
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.507085   charge     4.156167
 smvxcm (warning) mesh density negative at 16247 points:  rhomin=-3.81e-6
 smooth rhoeps =   -3.087349 (  -2.404736,  -0.682612)
         rhomu =   -4.032814 (  -3.321833,  -0.710981)
       avg vxc =   -0.203155 (  -0.236966,  -0.169345)
 smvxcm: negative smrho_w number,min(smrho_w)=       16247 -3.80813807009813665E-006
 smvxcm: enforce positive smrho_w. Add srshift=  3.80813808009813647E-006
 smooth rhoeps =   -3.088931 (  -2.405604,  -0.683327)
         rhomu =   -4.034858 (  -3.323082,  -0.711776)
       avg vxc =   -0.208996 (  -0.239351,  -0.178641)

 locpot:
  i job kmax lfltwf=           0           1           3 T

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.336160   -3.951582  -14.287741     -7.289568   -1.791221   -9.080789

 local terms:     true           smooth         local
 rhoeps:        -9.527348      -3.020566      -6.506781
 rhomu:         -7.410157      -3.248660      -4.161497
 spin2:         -5.126469      -0.697641      -4.428828
 total:        -12.536626      -3.946301      -8.590324
 val*vef       -14.287741      -8.624572      -5.663169
 val chg:        3.687810       3.843976      -0.156167
 val mom:        1.795694       2.110640      -0.314946    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000005

 potpus  spin 1 : pnu = 2.914756 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.933283    7.422828    0.092634   -0.588093
 1      3.000000    1.000000    -5.887832    7.124610    0.143880   -0.534122
 2      3.000000    1.000000     6.000000   29.619323    0.492331   -0.085993
 3      3.000000    1.000000     9.000000   40.176573    0.569665   -0.056303

 potpus  spin 2 : pnu = 2.908399 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.135490    7.333705    0.102233   -0.559932
 1      3.000000    1.000000    -5.887832    7.132731    0.153291   -0.501020
 2      3.000000    1.000000     6.000000   29.963353    0.494258   -0.084428
 3      3.000000    1.000000     9.000000   40.372644    0.570555   -0.055864

 Energy terms:             smooth           local           total
   rhoval*vef             -9.160685        -5.206953       -14.367638
   rhoval*ves             -4.671202        -5.584218       -10.255420
   psnuc*ves              11.685372      -281.331336      -269.645964
   utot                    3.507085      -143.457777      -139.950692
   rho*exc                -3.088931        -6.506781        -9.595712
   rho*vxc                -4.034858        -8.590324       -12.625182
   valence chg             4.156167        -0.156167         4.000000
   valence mag             2.314945        -0.314946         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.432042
 end of suham2
 -------- qplist --------
    1   0.000   0.000   0.000
    2   0.125   0.125  -0.125
    3   0.250   0.250  -0.250
    4   0.250   0.000   0.000
    5   0.375   0.125  -0.125
    6   0.500   0.250  -0.250
    7   0.500   0.000   0.000
    8   0.500   0.250   0.000
 sigmamode= F
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.0422 -0.4302 -0.4302 -0.4302  0.1244  0.5254  0.5254  0.5254
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.8411 -0.2395 -0.2395 -0.2395  0.1992  0.6212  0.6212  0.6212
 Est Ef = -0.432 < evl(4)=-0.430 ... using qval=4.0, revise to -0.4302
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.431974;   4.000000 electrons
         Sum occ. bands:   -2.750086, incl. Bloechl correction: -0.000188
         Mag. moment:       2.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.687083    3.843805   -0.156722      1.796013    2.113022   -0.317009
       contr. to mm extrapolated for r>rmt:   0.163832 est. true mm = 1.959845
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85755  sum tc=    31.38667  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78791  sum tc=    31.54433  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.958233   -1.040250    2.914756    2.914749    2.500000    2.914749
 spn 2 0    0.945535   -0.839200    2.908399    2.908394    2.500000    2.908394
 1     1    1.783303   -0.430892    2.850000    2.889155    2.250000    2.850000
 spn 2 1    0.000000   -1.118897    2.850000    2.165099    2.250000    2.850000
 2     0    0.000012   -0.811602    3.147584    3.132781    3.147584    3.147584
 spn 2 0    0.000000   -1.178847    3.147584    3.108369    3.147584    3.147584
 3     0    0.000001   -0.760483    4.102416    4.096174    4.102416    4.102416
 spn 2 0    0.000000   -1.131104    4.102416    4.085072    4.102416    4.102416

 Harris energy:
 sumev=       -2.750086  val*vef=     -14.367638   sumtv=      11.617552
 sumec=      -39.645456  cor*vef=    -102.576606   ttcor=      62.931150
 rhoeps=      -9.595712     utot=    -139.950692    ehar=     -74.997703
 smvxcm: negative smrho_w number,min(smrho_w)=       16247 -3.80811862699468268E-006
 smvxcm: enforce positive smrho_w. Add srshift=  3.80811863699468250E-006
 smvxcm: negative smrho_w number,min(smrho_w)=       16247 -3.80815751320159062E-006
 smvxcm: enforce positive smrho_w. Add srshift=  3.80815752320159044E-006
 smvxcm: negative smrho_w number,min(smrho_w)=       16247 -3.80813807009813665E-006
 smvxcm: enforce positive smrho_w. Add srshift=  3.80813808009813647E-006

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00   -0.00   -0.00     0.00    0.00   -0.00    -0.00   -0.00    0.00
 shift forces to make zero average correction:           -0.00   -0.00    0.00

 srhov:     -9.160423     -5.207352    -14.367776 sumev=   -2.750086   sumtv=   11.617690

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.044210   -0.156722     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt= 0.008211  avg sphere pot= 0.014115  vconst=-0.008211
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.509319   charge     4.156721
 smooth rhoeps =   -3.087699 (  -2.405990,  -0.681709)
         rhomu =   -4.033287 (  -3.323565,  -0.709722)
       avg vxc =   -0.206249 (  -0.237956,  -0.174542)
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -3.087699 (  -2.405990,  -0.681709)
         rhomu =   -4.033287 (  -3.323565,  -0.709722)
       avg vxc =   -0.206249 (  -0.237956,  -0.174542)

 locpot:
  i job kmax lfltwf=           0           0           3 F

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.337136   -3.950915  -14.288051     -7.293617   -1.788496   -9.082113

 local terms:     true           smooth         local
 rhoeps:        -9.527098      -3.021003      -6.506095
 rhomu:         -7.410228      -3.250480      -4.159748
 spin2:         -5.126079      -0.696417      -4.429662
 total:        -12.536307      -3.946897      -8.589410
 val*vef       -14.288051      -8.624143      -5.663908
 val chg:        3.687083       3.843805      -0.156722
 val mom:        1.796013       2.113022      -0.317009    core:   0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef             -9.161347        -5.205939       -14.367286
   rhoval*ves             -4.670091        -5.585872       -10.255964
   psnuc*ves              11.688729      -281.334797      -269.646068
   utot                    3.509319      -143.460335      -139.951016
   rho*exc                -3.087699        -6.506095        -9.593793
   rho*vxc                -4.033287        -8.589410       -12.622697
   valence chg             4.156721        -0.156722         4.000000
   valence mag             2.317009        -0.317009         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       11.617690  sumtc=        62.931001   ekin=       74.548691
 rhoep=       -9.593793   utot=      -139.951016   ehks=      -74.996118
 mag. mom=     2.000000  (bands)        2.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
 mixrho: sum smrho  init = 0.404445D+03-0.140900D-26 0.404445D+03       0
 mixrho: sum smrnew new  = 0.404608D+03 0.155430D-16 0.404608D+03       0
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 mixrho: dqsum rmsuns=  0.55474D-06  0.16092D-04 -0.21977D-18
 mixrealsmooth= T
 wgtsmooth=  2.82842712474619005E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.34e-5  last it=1.11e-5
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=252052  beta=0.5  tm=5  rmsdel=6.7e-6
   tj: 3.35624
 mixrealsmooth= T
 smrho qcell: add correction to smrho=  4.20249813837259012E-008  2.10124906918629567E-011
 unscreened rms difference:  smooth  0.000023   local  0.000060
   screened rms difference:  smooth  0.000014   local  0.000060   tot  0.000013
 mixrho: all smrho is positive for isp=           1
 mixrho: warning. negative smrho; isp number min=           2       26087 -5.94796278093544935E-006

 iors  : write restart file (binary, mesh density) 

   it  7  of 10    ehf=      -0.002803   ehk=      -0.001218
 From last iter    ehf=      -0.003270   ehk=      -0.001218
 diffe(q)=  0.000467 (0.000013)    tol= 0.000010 (0.000500)   more=T
i zbak=0 mmom=1.9999999 ehf=-.0028026 ehk=-.001218

 --- BNDFP:  begin iteration 8 of 10 ---
 ttt nevmx w=           5  4.00000000000000008E-003

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.043869   -0.155513     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt= 0.008671  avg sphere pot= 0.014128  vconst=-0.008671
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.504660   charge     4.155513
 smvxcm (warning) mesh density negative at 26087 points:  rhomin=-5.95e-6
 smooth rhoeps =   -3.086620 (  -2.403083,  -0.683537)
         rhomu =   -4.031850 (  -3.319560,  -0.712290)
       avg vxc =   -0.201431 (  -0.236187,  -0.166674)
 smvxcm: negative smrho_w number,min(smrho_w)=       26087 -5.94796278093544935E-006
 smvxcm: enforce positive smrho_w. Add srshift=  5.94796279093544918E-006
 smooth rhoeps =   -3.089101 (  -2.404447,  -0.684653)
         rhomu =   -4.035043 (  -3.321516,  -0.713527)
       avg vxc =   -0.210213 (  -0.240067,  -0.180358)

 locpot:
  i job kmax lfltwf=           0           1           3 T

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.334989   -3.952190  -14.287179     -7.284839   -1.794280   -9.079119

 local terms:     true           smooth         local
 rhoeps:        -9.527328      -3.019756      -6.507572
 rhomu:         -7.409839      -3.246291      -4.163548
 spin2:         -5.126754      -0.698924      -4.427830
 total:        -12.536592      -3.945215      -8.591377
 val*vef       -14.287179      -8.624941      -5.662238
 val chg:        3.688215       3.843728      -0.155513
 val mom:        1.795319       2.107833      -0.312514    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000005

 potpus  spin 1 : pnu = 2.914749 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.932258    7.422684    0.092640   -0.588095
 1      3.000000    1.000000    -5.887832    7.124359    0.143881   -0.534130
 2      3.000000    1.000000     6.000000   29.619780    0.492335   -0.085991
 3      3.000000    1.000000     9.000000   40.176974    0.569668   -0.056302

 potpus  spin 2 : pnu = 2.908394 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.134965    7.333548    0.102235   -0.559941
 1      3.000000    1.000000    -5.887832    7.132456    0.153289   -0.501035
 2      3.000000    1.000000     6.000000   29.963750    0.494262   -0.084426
 3      3.000000    1.000000     9.000000   40.373020    0.570558   -0.055863

 Energy terms:             smooth           local           total
   rhoval*vef             -9.159432        -5.208061       -14.367493
   rhoval*ves             -4.672712        -5.582209       -10.254922
   psnuc*ves              11.682032      -281.326712      -269.644680
   utot                    3.504660      -143.454461      -139.949801
   rho*exc                -3.089101        -6.507572        -9.596673
   rho*vxc                -4.035043        -8.591377       -12.626420
   valence chg             4.155513        -0.155513         4.000000
   valence mag             2.312514        -0.312514         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000

 Read qp weights ...  ef=-0.431974
 end of suham2
 -------- qplist --------
    1   0.000   0.000   0.000
    2   0.125   0.125  -0.125
    3   0.250   0.250  -0.250
    4   0.250   0.000   0.000
    5   0.375   0.125  -0.125
    6   0.500   0.250  -0.250
    7   0.500   0.000   0.000
    8   0.500   0.250   0.000
 sigmamode= F
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.0424 -0.4303 -0.4303 -0.4303  0.1217  0.5242  0.5242  0.5242
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.8413 -0.2397 -0.2397 -0.2397  0.1941  0.6190  0.6190  0.6190
 Est Ef = -0.432 < evl(4)=-0.430 ... using qval=4.0, revise to -0.4303
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.432116;   4.000000 electrons
         Sum occ. bands:   -2.750697, incl. Bloechl correction: -0.000188
         Mag. moment:       2.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.686848    3.843218   -0.156370      1.795751    2.110246   -0.314495
       contr. to mm extrapolated for r>rmt:   0.164023 est. true mm = 1.959775
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85781  sum tc=    31.38675  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78820  sum tc=    31.54435  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.958250   -1.040293    2.914749    2.914765    2.500000    2.914765
 spn 2 0    0.945548   -0.839313    2.908394    2.908405    2.500000    2.908405
 1     1    1.783037   -0.431033    2.850000    2.889094    2.250000    2.850000
 spn 2 1    0.000000   -1.122191    2.850000    2.164638    2.250000    2.850000
 2     0    0.000012   -0.811655    3.147584    3.132785    3.147584    3.147584
 spn 2 0    0.000000   -1.178896    3.147584    3.108374    3.147584    3.147584
 3     0    0.000001   -0.760507    4.102416    4.096177    4.102416    4.102416
 spn 2 0    0.000000   -1.131526    4.102416    4.085066    4.102416    4.102416

 Harris energy:
 sumev=       -2.750697  val*vef=     -14.367493   sumtv=      11.616795
 sumec=      -39.646007  cor*vef=    -102.577082   ttcor=      62.931075
 rhoeps=      -9.596673     utot=    -139.949801    ehar=     -74.998603
 smvxcm: negative smrho_w number,min(smrho_w)=       26087 -5.94795611891675129E-006
 smvxcm: enforce positive smrho_w. Add srshift=  5.94795612891675111E-006
 smvxcm: negative smrho_w number,min(smrho_w)=       26087 -5.94796944295414742E-006
 smvxcm: enforce positive smrho_w. Add srshift=  5.94796945295414725E-006
 smvxcm: negative smrho_w number,min(smrho_w)=       26087 -5.94796278093544935E-006
 smvxcm: enforce positive smrho_w. Add srshift=  5.94796279093544918E-006

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00   -0.00   -0.00     0.00    0.00    0.00    -0.00   -0.00   -0.00
 shift forces to make zero average correction:           -0.00   -0.00   -0.00

 srhov:     -9.160177     -5.207502    -14.367679 sumev=   -2.750697   sumtv=   11.616982

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.044111   -0.156370     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt= 0.008177  avg sphere pot= 0.014123  vconst=-0.008177
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.506828   charge     4.156370
 smooth rhoeps =   -3.087025 (  -2.404377,  -0.682648)
         rhomu =   -4.032383 (  -3.321372,  -0.711011)
       avg vxc =   -0.206295 (  -0.238016,  -0.174574)
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -3.087025 (  -2.404377,  -0.682648)
         rhomu =   -4.032383 (  -3.321372,  -0.711011)
       avg vxc =   -0.206295 (  -0.238016,  -0.174574)

 locpot:
  i job kmax lfltwf=           0           0           3 F

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.336493   -3.951289  -14.287782     -7.289584   -1.792149   -9.081733

 local terms:     true           smooth         local
 rhoeps:        -9.526805      -3.020275      -6.506530
 rhomu:         -7.409808      -3.248213      -4.161595
 spin2:         -5.126113      -0.697709      -4.428404
 total:        -12.535921      -3.945922      -8.589999
 val*vef       -14.287782      -8.624916      -5.662865
 val chg:        3.686848       3.843218      -0.156370
 val mom:        1.795751       2.110246      -0.314495    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef             -9.161024        -5.206050       -14.367074
   rhoval*ves             -4.671825        -5.584240       -10.256065
   psnuc*ves              11.685480      -281.330710      -269.645230
   utot                    3.506828      -143.457475      -139.950647
   rho*exc                -3.087025        -6.506530        -9.593555
   rho*vxc                -4.032383        -8.589999       -12.622383
   valence chg             4.156370        -0.156370         4.000000
   valence mag             2.314495        -0.314495         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       11.616982  sumtc=        62.931102   ekin=       74.548084
 rhoep=       -9.593555   utot=      -139.950647   ehks=      -74.996118
 mag. mom=     2.000000  (bands)        2.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
 mixrho: sum smrho  init = 0.404252D+03-0.109228D-26 0.404252D+03       0
 mixrho: sum smrnew new  = 0.404429D+03 0.322300D-16 0.404429D+03       0
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 mixrho: dqsum rmsuns=  0.85718D-06  0.17718D-04 -0.23716D-18
 mixrealsmooth= T
 wgtsmooth=  2.82842712474619005E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.43e-5  last it=1.34e-5
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=252052  beta=0.5  tm=5  rmsdel=7.13e-6
   tj: 1.73918
 mixrealsmooth= T
 smrho qcell: add correction to smrho=  4.20897572905865047E-008  2.10448786452932578E-011
 unscreened rms difference:  smooth  0.000025   local  0.000064
   screened rms difference:  smooth  0.000016   local  0.000064   tot  0.000014
 mixrho: all smrho is positive for isp=           1
 mixrho: warning. negative smrho; isp number min=           2        6143 -1.63536824012925451E-006

 iors  : write restart file (binary, mesh density) 

   it  8  of 10    ehf=      -0.003703   ehk=      -0.001218
 From last iter    ehf=      -0.002803   ehk=      -0.001218
 diffe(q)= -0.000900 (0.000014)    tol= 0.000010 (0.000500)   more=T
i zbak=0 mmom=1.9999999 ehf=-.003703 ehk=-.0012175

 --- BNDFP:  begin iteration 9 of 10 ---
 ttt nevmx w=           5  4.00000000000000008E-003

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.044237   -0.156815     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt= 0.008324  avg sphere pot= 0.014112  vconst=-0.008324
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.509988   charge     4.156815
 smvxcm (warning) mesh density negative at 6143 points:  rhomin=-1.64e-6
 smooth rhoeps =   -3.088096 (  -2.406600,  -0.681496)
         rhomu =   -4.033814 (  -3.324389,  -0.709424)
       avg vxc =   -0.205042 (  -0.237593,  -0.172491)
 smvxcm: negative smrho_w number,min(smrho_w)=        6143 -1.63536824012925451E-006
 smvxcm: enforce positive smrho_w. Add srshift=  1.63536825012925454E-006
 smooth rhoeps =   -3.088773 (  -2.406970,  -0.681803)
         rhomu =   -4.034691 (  -3.324925,  -0.709766)
       avg vxc =   -0.207607 (  -0.238580,  -0.176634)

 locpot:
  i job kmax lfltwf=           0           1           3 T

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.337321   -3.950913  -14.288234     -7.294826   -1.787401   -9.082227

 local terms:     true           smooth         local
 rhoeps:        -9.527389      -3.021403      -6.505986
 rhomu:         -7.410503      -3.251321      -4.159182
 spin2:         -5.126183      -0.696104      -4.430079
 total:        -12.536687      -3.947425      -8.589262
 val*vef       -14.288234      -8.623957      -5.664277
 val chg:        3.687459       3.844275      -0.156815
 val mom:        1.796089       2.113894      -0.317805    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000005

 potpus  spin 1 : pnu = 2.914765 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.934398    7.422941    0.092628   -0.588089
 1      3.000000    1.000000    -5.887832    7.124849    0.143880   -0.534113
 2      3.000000    1.000000     6.000000   29.618909    0.492328   -0.085995
 3      3.000000    1.000000     9.000000   40.176201    0.569663   -0.056304

 potpus  spin 2 : pnu = 2.908405 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.136164    7.333812    0.102230   -0.559921
 1      3.000000    1.000000    -5.887832    7.132993    0.153293   -0.501002
 2      3.000000    1.000000     6.000000   29.963024    0.494255   -0.084430
 3      3.000000    1.000000     9.000000   40.372315    0.570553   -0.055865

 Energy terms:             smooth           local           total
   rhoval*vef             -9.161695        -5.206008       -14.367703
   rhoval*ves             -4.669414        -5.586401       -10.255814
   psnuc*ves              11.689390      -281.336341      -269.646951
   utot                    3.509988      -143.461371      -139.951383
   rho*exc                -3.088773        -6.505986        -9.594759
   rho*vxc                -4.034691        -8.589262       -12.623953
   valence chg             4.156815        -0.156815         4.000000
   valence mag             2.317805        -0.317805         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.432116
 end of suham2
 -------- qplist --------
    1   0.000   0.000   0.000
    2   0.125   0.125  -0.125
    3   0.250   0.250  -0.250
    4   0.250   0.000   0.000
    5   0.375   0.125  -0.125
    6   0.500   0.250  -0.250
    7   0.500   0.000   0.000
    8   0.500   0.250   0.000
 sigmamode= F
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.0421 -0.4301 -0.4301 -0.4301  0.1271  0.5266  0.5266  0.5266
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.8409 -0.2392 -0.2392 -0.2392  0.2045  0.6235  0.6235  0.6235
 Est Ef = -0.432 < evl(4)=-0.430 ... using qval=4.0, revise to -0.4301
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.431826;   4.000000 electrons
         Sum occ. bands:   -2.749445, incl. Bloechl correction: -0.000187
         Mag. moment:       2.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.687323    3.844532   -0.157208      1.796277    2.115764   -0.319487
       contr. to mm extrapolated for r>rmt:   0.163641 est. true mm = 1.959918
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85729  sum tc=    31.38661  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78763  sum tc=    31.54431  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.958214   -1.040201    2.914765    2.914732    2.500000    2.914732
 spn 2 0    0.945523   -0.839078    2.908405    2.908385    2.500000    2.908385
 1     1    1.783573   -0.430743    2.850000    2.889216    2.250000    2.850000
 spn 2 1    0.000000   -1.114342    2.850000    2.165752    2.250000    2.850000
 2     0    0.000012   -0.811553    3.147584    3.132777    3.147584    3.147584
 spn 2 0    0.000000   -1.178820    3.147584    3.108364    3.147584    3.147584
 3     0    0.000001   -0.760465    4.102416    4.096171    4.102416    4.102416
 spn 2 0    0.000000   -1.130672    4.102416    4.085077    4.102416    4.102416

 Harris energy:
 sumev=       -2.749445  val*vef=     -14.367703   sumtv=      11.618257
 sumec=      -39.644913  cor*vef=    -102.576002   ttcor=      62.931089
 rhoeps=      -9.594759     utot=    -139.951383    ehar=     -74.996796
 smvxcm: negative smrho_w number,min(smrho_w)=        6119 -1.63533582404296487E-006
 smvxcm: enforce positive smrho_w. Add srshift=  1.63533583404296491E-006
 smvxcm: negative smrho_w number,min(smrho_w)=        6143 -1.63540065621554436E-006
 smvxcm: enforce positive smrho_w. Add srshift=  1.63540066621554439E-006
 smvxcm: negative smrho_w number,min(smrho_w)=        6143 -1.63536824012925472E-006
 smvxcm: enforce positive smrho_w. Add srshift=  1.63536825012925476E-006

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00   -0.00   -0.00     0.00   -0.00   -0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:     -9.161268     -5.206571    -14.367840 sumev=   -2.749445   sumtv=   11.618394

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.044348   -0.157208     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt= 0.008247  avg sphere pot= 0.014107  vconst=-0.008247
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.511979   charge     4.157208
 smooth rhoeps =   -3.088512 (  -2.407658,  -0.680854)
         rhomu =   -4.034373 (  -3.325837,  -0.708535)
       avg vxc =   -0.206203 (  -0.237896,  -0.174510)
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -3.088512 (  -2.407658,  -0.680854)
         rhomu =   -4.034373 (  -3.325837,  -0.708535)
       avg vxc =   -0.206203 (  -0.237896,  -0.174510)

 locpot:
  i job kmax lfltwf=           0           0           3 F

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.337777   -3.950537  -14.288314     -7.297838   -1.785153   -9.082991

 local terms:     true           smooth         local
 rhoeps:        -9.527395      -3.021872      -6.505523
 rhomu:         -7.410651      -3.252826      -4.157825
 spin2:         -5.126046      -0.695228      -4.430818
 total:        -12.536697      -3.948054      -8.588643
 val*vef       -14.288314      -8.623461      -5.664853
 val chg:        3.687323       3.844532      -0.157208
 val mom:        1.796277       2.115764      -0.319487    core:   0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef             -9.162168        -5.205324       -14.367492
   rhoval*ves             -4.668266        -5.587587       -10.255853
   psnuc*ves              11.692224      -281.339176      -269.646952
   utot                    3.511979      -143.463382      -139.951403
   rho*exc                -3.088512        -6.505523        -9.594035
   rho*vxc                -4.034373        -8.588643       -12.623016
   valence chg             4.157208        -0.157208         4.000000
   valence mag             2.319487        -0.319487         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       11.618394  sumtc=        62.930925   ekin=       74.549319
 rhoep=       -9.594035   utot=      -139.951403   ehks=      -74.996118
 mag. mom=     2.000000  (bands)        2.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
 mixrho: sum smrho  init = 0.404664D+03 0.297840D-25 0.404664D+03       0
 mixrho: sum smrnew new  = 0.404793D+03 0.663365D-16 0.404793D+03       0
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 mixrho: dqsum rmsuns=  0.39286D-06  0.12514D-04 -0.41291D-18
 mixrealsmooth= T
 wgtsmooth=  2.82842712474619005E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.07e-5  last it=1.43e-5
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=252052  beta=0.5  tm=5  rmsdel=5.33e-6
   tj:-0.98157
 mixrealsmooth= T
 smrho qcell: add correction to smrho=  4.27498093036327020E-008  2.13749046518163549E-011
 unscreened rms difference:  smooth  0.000018   local  0.000047
   screened rms difference:  smooth  0.000012   local  0.000047   tot  0.000011
 mixrho: all smrho is positive for isp=           1
 mixrho: all smrho is positive for isp=           2

 iors  : write restart file (binary, mesh density) 

   it  9  of 10    ehf=      -0.001896   ehk=      -0.001218
 From last iter    ehf=      -0.003703   ehk=      -0.001218
 diffe(q)=  0.001807 (0.000011)    tol= 0.000010 (0.000500)   more=T
i zbak=0 mmom=1.9999999 ehf=-.0018956 ehk=-.0012184

 --- BNDFP:  begin iteration 10 of 10 ---
 ttt nevmx w=           5  4.00000000000000008E-003

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.044589   -0.158063     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt= 0.008075  avg sphere pot= 0.014094  vconst=-0.008075
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.516257   charge     4.158063
 smooth rhoeps =   -3.089589 (  -2.410372,  -0.679217)
         rhomu =   -4.035826 (  -3.329551,  -0.706276)
       avg vxc =   -0.207603 (  -0.238374,  -0.176832)
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -3.089589 (  -2.410372,  -0.679217)
         rhomu =   -4.035826 (  -3.329551,  -0.706276)
       avg vxc =   -0.207603 (  -0.238374,  -0.176832)

 locpot:
  i job kmax lfltwf=           0           1           3 T

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.339310   -3.949622  -14.288932     -7.305314   -1.779376   -9.084690

 local terms:     true           smooth         local
 rhoeps:        -9.527544      -3.023064      -6.504479
 rhomu:         -7.411193      -3.256682      -4.154511
 spin2:         -5.125707      -0.692974      -4.432733
 total:        -12.536900      -3.949656      -8.587244
 val*vef       -14.288932      -8.622432      -5.666500
 val chg:        3.687001       3.845064      -0.158063
 val mom:        1.796818       2.120511      -0.323693    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000005

 potpus  spin 1 : pnu = 2.914732 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.930053    7.423658    0.092651   -0.588060
 1      3.000000    1.000000    -5.887832    7.125232    0.143882   -0.534089
 2      3.000000    1.000000     6.000000   29.618318    0.492323   -0.085998
 3      3.000000    1.000000     9.000000   40.175642    0.569660   -0.056305

 potpus  spin 2 : pnu = 2.908385 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.133894    7.334302    0.102248   -0.559882
 1      3.000000    1.000000    -5.887832    7.133418    0.153301   -0.500961
 2      3.000000    1.000000     6.000000   29.962643    0.494251   -0.084432
 3      3.000000    1.000000     9.000000   40.371874    0.570550   -0.055866

 Energy terms:             smooth           local           total
   rhoval*vef             -9.163666        -5.204242       -14.367908
   rhoval*ves             -4.665582        -5.590661       -10.256243
   psnuc*ves              11.698097      -281.346745      -269.648648
   utot                    3.516257      -143.468703      -139.952445
   rho*exc                -3.089589        -6.504479        -9.594068
   rho*vxc                -4.035826        -8.587244       -12.623070
   valence chg             4.158063        -0.158063         4.000000
   valence mag             2.323693        -0.323693         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000

 Read qp weights ...  ef=-0.431826
 end of suham2
 -------- qplist --------
    1   0.000   0.000   0.000
    2   0.125   0.125  -0.125
    3   0.250   0.250  -0.250
    4   0.250   0.000   0.000
    5   0.375   0.125  -0.125
    6   0.500   0.250  -0.250
    7   0.500   0.000   0.000
    8   0.500   0.250   0.000
 sigmamode= F
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.0420 -0.4299 -0.4299 -0.4299  0.1293  0.5277  0.5277  0.5277
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.8406 -0.2388 -0.2388 -0.2388  0.2082  0.6252  0.6252  0.6252
 Est Ef = -0.432 < evl(4)=-0.430 ... using qval=4.0, revise to -0.4299
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0
  mode napw           0           0
  end of hambls mode=           0

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.431612;   4.000000 electrons
         Sum occ. bands:   -2.748519, incl. Bloechl correction: -0.000186
         Mag. moment:       2.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.687484    3.845529   -0.158045      1.796472    2.118540   -0.322068
       contr. to mm extrapolated for r>rmt:   0.163490 est. true mm = 1.959962
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85682  sum tc=    31.38651  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78712  sum tc=    31.54431  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.958205   -1.040071    2.914732    2.914725    2.500000    2.914725
 spn 2 0    0.945506   -0.838837    2.908385    2.908375    2.500000    2.908375
 1     1    1.783760   -0.430541    2.850000    2.889261    2.250000    2.850000
 spn 2 1    0.000000   -1.112748    2.850000    2.165955    2.250000    2.850000
 2     0    0.000012   -0.811417    3.147584    3.132774    3.147584    3.147584
 spn 2 0    0.000000   -1.178611    3.147584    3.108363    3.147584    3.147584
 3     0    0.000001   -0.760336    4.102416    4.096169    4.102416    4.102416
 spn 2 0    0.000000   -1.130314    4.102416    4.085079    4.102416    4.102416

 Harris energy:
 sumev=       -2.748519  val*vef=     -14.367908   sumtv=      11.619390
 sumec=      -39.643935  cor*vef=    -102.574942   ttcor=      62.931007
 rhoeps=      -9.594068     utot=    -139.952445    ehar=     -74.996117
 smvxcm: all smrho_w is positive
 smvxcm: all smrho_w is positive
 smvxcm: all smrho_w is positive

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00    0.00   -0.00     0.00   -0.00   -0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:     -9.165760     -5.201172    -14.366933 sumev=   -2.748519   sumtv=   11.618414

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.044584   -0.158045     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt= 0.008273  avg sphere pot= 0.014100  vconst=-0.008273
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |

 smooth rhoves      3.514289   charge     4.158045
 smooth rhoeps =   -3.089768 (  -2.409646,  -0.680123)
         rhomu =   -4.036040 (  -3.328554,  -0.707486)
       avg vxc =   -0.206175 (  -0.237858,  -0.174492)
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -3.089768 (  -2.409646,  -0.680123)
         rhomu =   -4.036040 (  -3.328554,  -0.707486)
       avg vxc =   -0.206175 (  -0.237858,  -0.174492)

 locpot:
  i job kmax lfltwf=           0           0           3 F

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1     -10.338060   -3.950078  -14.288138     -7.303621   -1.782171   -9.085792

 local terms:     true           smooth         local
 rhoeps:        -9.527545      -3.023166      -6.504379
 rhomu:         -7.410905      -3.255596      -4.155309
 spin2:         -5.125989      -0.694174      -4.431815
 total:        -12.536894      -3.949771      -8.587123
 val*vef       -14.288138      -8.623694      -5.664444
 val chg:        3.687484       3.845529      -0.158045
 val mom:        1.796472       2.118540      -0.322068    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef             -9.164930        -5.202346       -14.367277
   rhoval*ves             -4.666792        -5.588685       -10.255477
   psnuc*ves              11.695370      -281.342289      -269.646920
   utot                    3.514289      -143.465487      -139.951199
   rho*exc                -3.089768        -6.504379        -9.594147
   rho*vxc                -4.036040        -8.587123       -12.623164
   valence chg             4.158045        -0.158045         4.000000
   valence mag             2.322068        -0.322068         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       11.618414  sumtc=        62.930813   ekin=       74.549227
 rhoep=       -9.594147   utot=      -139.951199   ehks=      -74.996119
 mag. mom=     2.000000  (bands)        2.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
 mixrho: sum smrho  init = 0.405110D+03-0.273912D-26 0.405110D+03       0
 mixrho: sum smrnew new  = 0.405007D+03-0.321609D-16 0.405007D+03       0
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 mixrho: dqsum rmsuns= -0.17940D-07  0.17768D-04  0.34922D-18
 mixrealsmooth= T
 wgtsmooth=  2.82842712474619005E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.42e-5  last it=1.07e-5
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=252052  beta=0.5  tm=5  rmsdel=7.08e-6
   tj: 0.57244
 mixrealsmooth= T
 smrho qcell: add correction to smrho=  4.20770874531850581E-008  2.10385437265925338E-011
 unscreened rms difference:  smooth  0.000025   local  0.000057
   screened rms difference:  smooth  0.000018   local  0.000057   tot  0.000014
 mixrho: all smrho is positive for isp=           1
 mixrho: all smrho is positive for isp=           2

 iors  : write restart file (binary, mesh density) 

   it 10  of 10    ehf=      -0.001217   ehk=      -0.001219
 From last iter    ehf=      -0.001896   ehk=      -0.001218
 diffe(q)=  0.000679 (0.000014)    tol= 0.000010 (0.000500)   more=F
x zbak=0 mmom=1.9999999 ehf=-.0012168 ehk=-.0012187
 Exit 0 LMF 
