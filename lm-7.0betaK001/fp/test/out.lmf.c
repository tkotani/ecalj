rdcmd:  lmfa --no-iactiv c -vzbak=1
 -----------------------  START LMFA (80000K)  -----------------------
 HEADER sc C atom

 LMFA:     alat = 7.93701  nbas = 1  nspec = 1  vn 7.00(LMFA 7.0)  verb 30,40,|
 pot:      spin-pol, XC:BH

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

 tailsm: fit tails to 6 smoothed hankels, rmt= 3.00000, rsm= 1.50000
    q(fit):     0.243570    rms diff:   0.000004
    fit: r>rmt  0.243570   r<rmt  1.753709   qtot  1.997279
    rho: r>rmt  0.243570   r<rmt  2.756430   qtot  3.000000

 tailsm: spin 2 ...
    q(fit):     0.054561    rms diff:   0.000002
    fit: r>rmt  0.054561   r<rmt  0.609878   qtot  0.664439
    rho: r>rmt  0.054561   r<rmt  0.945439   qtot  1.000000
 
  Write mtopara.* ...
 Exit 0 LMFA 
 wkinfo:  used   101 K  workspace of 80000 K   in   0 K calls
rdcmd:  lmf  --no-iactiv c -vzbak=1
 -----------------------  START LMF (80000K)  -----------------------
 HEADER sc C atom

 LMF:      alat = 7.93701  nbas = 1  nspec = 1  vn 7.00(LMF 7.0)  verb 30,40,60
 special:  forces
 pot:      spin-pol, XC:BH
 bz:       metal(2), tetra, invit 

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
 
 BZMESH:  8 irreducible QP from 64 ( 4 4 4 )  shift= F F F

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
 Uniform density added to neutralize background, q=1.000000

 Smooth charge on mesh:            1.661718    moment    1.332840
 Sum of local charges:             1.338282    moments   0.667159
 Total valence charge:             3.000000    moment    2.000000
 Sum of core charges:              2.000000    moment    0.000000
 Sum of nuclear charges:          -6.000000
 Homogeneous background:           1.000000
 Deviation from neutrality:        0.000000

 --- BNDFP:  begin iteration 1 of 10 ---

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1    0.377522    1.338282     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt= 0.006607  avg sphere pot= 0.019541  vconst=-0.006607
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      2.063910   charge     2.661718
 smvxcm (warning) mesh density negative at 201254 points:  rhomin=-5e-4
 smooth rhoeps =   -1.385796 (  -1.043174,  -0.342622)
         rhomu =   -1.793806 (  -1.418944,  -0.374861)
       avg vxc =   -0.074422 (  -0.080699,  -0.068146)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T

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
 val mom:        1.810990       1.143831       0.667159    core:   0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000003

 potpus  spin 1 : pnu = 2.900000 2.850000 3.180000 4.120000
 l        enu         v           c          srdel        qpar        ppar
 0     -1.049913   -1.260952   -1.152144    0.197581      2.7872    0.898890
 1     -0.465402   -1.039466   -0.463120    0.173954     19.0465    0.926619
 2     -0.311939   -0.604280    1.529962    0.353093     17.1183    6.114474
 3     -0.096237   -0.537542    3.382071    0.400473     24.4391   10.634342

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000    -9.233051    7.622125    0.101666   -0.583560
 1      3.000000    1.000000    -5.887832    7.139035    0.143256   -0.535856
 2      3.000000    1.000000     4.727244   27.072458    0.465408   -0.096156
 3      3.000000    1.000000     7.577135   37.163437    0.543349   -0.062204

 potpus  spin 2 : pnu = 2.900000 2.850000 3.180000 4.120000
 l        enu         v           c          srdel        qpar        ppar
 0     -0.838165   -1.073432   -0.951603    0.208583      2.8002    0.932015
 1     -0.258802   -0.900692   -0.256240    0.184285     18.9762    0.977113
 2     -0.186358   -0.480549    1.677819    0.356098     17.0208    6.250867
 3      0.023022   -0.419226    3.516831    0.401776     24.3826   10.731834

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000    -9.233051    7.498210    0.106861   -0.559304
 1      3.000000    1.000000    -5.887832    7.161455    0.151762   -0.504954
 2      3.000000    1.000000     4.727244   27.365015    0.467061   -0.094577
 3      3.000000    1.000000     7.577135   37.316588    0.544000   -0.061810

 Energy terms:             smooth           local           total
   rhoval*vef             -3.772262       -10.403599       -14.175861
   rhoval*ves             -5.058553        -5.181775       -10.240327
   psnuc*ves               9.186373      -278.836573      -269.650199
   utot                    2.063910      -142.009174      -139.945263
   rho*exc                -1.385796        -8.100688        -9.486484
   rho*vxc                -1.793806       -10.678430       -12.472235
   valence chg             1.661718         1.338282         3.000000
   valence mag             1.332840         0.667159         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Incompatible or missing qp weights file ...
 Start first of two band passes ...
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.0284 -0.4081 -0.4081 -0.4081  0.3060  0.7027  0.7027  0.7027
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.8236 -0.2132 -0.2132 -0.2132  0.3259  0.7397  0.7397  0.7397

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.409798;   3.000000 electrons
         Sum occ. bands:   -2.261681, incl. Bloechl correction: -0.000048
         Mag. moment:       1.000000

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.0284 -0.4081 -0.4081 -0.4081  0.3060  0.7027  0.7027  0.7027
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.8236 -0.2132 -0.2132 -0.2132  0.3259  0.7397  0.7397  0.7397
 Est Ef = -0.410 < evl(3)=-0.408 ... using qval=3.0, revise to -0.4081

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.409798;   3.000000 electrons
         Sum occ. bands:   -2.261681, incl. Bloechl correction: -0.000048
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.842313    5.512615   -2.670302      0.935009    2.546033   -1.611024
       contr. to mm extrapolated for r>rmt:   0.050945 est. true mm = 0.985954

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.966101   -1.025877    2.900000    2.921600    2.500000    2.921600
 spn 2 0    0.953652   -0.821345    2.900000    2.914124    2.500000    2.914124
 1     1    0.922558   -0.406129    2.850000    2.903768    2.250000    2.850000
 spn 2 1    0.000000   -0.978105    2.850000    2.184937    2.250000    2.850000
 2     0    0.000002   -0.824147    3.180000    3.130921    3.147584    3.147584
 spn 2 0    0.000000   -1.162024    3.180000    3.107693    3.147584    3.147584
 3     0    0.000000   -0.793682    4.120000    4.094807    4.102416    4.102416
 spn 2 0    0.000000   -1.113934    4.120000    4.084684    4.102416    4.102416

 Harris energy:
 sumev=       -2.261681  val*vef=     -14.175861   sumtv=      11.914180
 sumec=      -39.640082  cor*vef=    -102.572087   ttcor=      62.932005
 rhoeps=      -9.486484     utot=    -139.945263    ehar=     -74.585563

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:    -13.118729      1.457883    -11.660846 sumev=   -2.261681   sumtv=    9.399165

 Kohn-Sham energy:
 sumtv=        9.399165  sumtc=        62.933572   ekin=       72.332737
 rhoep=       -8.678280   utot=      -138.059302   ehks=      -74.404845
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho:  sought 2 iter from file mixm; read 0.  RMS DQ=3.81e-2
 mixrho: (warning) scr. and lin-mixed densities had 81509 and 91835 negative poi
 AMIX: nmix=0 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=1.9e-2
 unscreened rms difference:  smooth  0.030368   local  0.096573
   screened rms difference:  smooth  0.030891   local  0.096573   tot  0.038099

 iors  : write restart file (binary, mesh density) 

   it  1  of 10    ehf=       0.409337   ehk=       0.590055
h zbak=1 mmom=.9999997 ehf=.4093368 ehk=.5900552

 --- BNDFP:  begin iteration 2 of 10 ---

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.187878   -0.666010     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.014474  avg sphere pot= 0.015037  vconst= 0.014474
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      3.778396   charge     4.666010
 smvxcm (warning) mesh density negative at 185118 points:  rhomin=-5.14e-4
 smooth rhoeps =   -3.420968 (  -2.507556,  -0.913412)
         rhomu =   -4.468228 (  -3.467846,  -1.000382)
       avg vxc =   -0.101336 (  -0.110205,  -0.092468)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -9.349535   -4.220007  -13.569542     -8.476395   -2.790222  -11.266617

 local terms:     true           smooth         local
 rhoeps:        -9.156633      -3.385059      -5.771575
 rhomu:         -6.881568      -3.431946      -3.449622
 spin2:         -5.168106      -0.991153      -4.176953
 total:        -12.049674      -4.423099      -7.626575
 val*vef       -13.569542      -9.208251      -4.361291
 val chg:        3.368708       4.034718      -0.666010
 val mom:        1.373000       1.844932      -0.471932    core:   0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000005

 potpus  spin 1 : pnu = 2.921600 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.148604   -1.360414   -1.262054    0.187585      2.7952    0.851964
 1     -0.582397   -1.083687   -0.580401    0.162680     19.0172    0.873038
 2     -0.625497   -0.625497    1.567371    0.365149     16.4460    6.754050
 3     -0.552379   -0.552379    3.496115    0.415130     23.4912   11.794126

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -11.932903    7.544922    0.082919   -0.619156
 1      3.000000    1.000000    -5.887832    7.148355    0.133970   -0.572585
 2      3.000000    1.000000     6.000000   29.270167    0.490327   -0.087640
 3      3.000000    1.000000     9.000000   39.976093    0.568742   -0.056759

 potpus  spin 2 : pnu = 2.914124 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -0.986629   -1.211627   -1.102730    0.197065      2.8041    0.882737
 1     -0.417464   -0.981526   -0.415233    0.171961     19.1506    0.913607
 2     -0.543112   -0.543112    1.675505    0.368193     16.3652    6.899349
 3     -0.476373   -0.476373    3.591354    0.416560     23.4410   11.904382

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.848787    7.461678    0.091932   -0.594062
 1      3.000000    1.000000    -5.887832    7.106278    0.141617   -0.543421
 2      3.000000    1.000000     6.000000   29.565557    0.491984   -0.086250
 3      3.000000    1.000000     9.000000   40.140837    0.569482   -0.056386

 Energy terms:             smooth           local           total
   rhoval*vef            -11.346630        -2.302928       -13.649558
   rhoval*ves             -4.805565        -5.246408       -10.051973
   psnuc*ves              12.362357      -280.744245      -268.381888
   utot                    3.778396      -142.995326      -139.216930
   rho*exc                -3.420968        -5.771575        -9.192542
   rho*vxc                -4.468228        -7.626575       -12.094803
   valence chg             3.666010        -0.666010         3.000000
   valence mag             1.971932        -0.471932         1.500000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.409798
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1504 -0.5311 -0.5311 -0.5311  0.2537  0.6335  0.6335  0.6335
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.9847 -0.3695 -0.3695 -0.3695  0.2859  0.6774  0.6774  0.6774

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.531737;   3.000000 electrons
         Sum occ. bands:   -2.666741, incl. Bloechl correction: -0.000036
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.872794    5.802109   -2.929314      0.946953    2.066495   -1.119542
       contr. to mm extrapolated for r>rmt:   0.043140 est. true mm = 0.990093

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.970177   -1.147871    2.921600    2.922271    2.500000    2.922271
 spn 2 0    0.962920   -0.982124    2.914124    2.918077    2.500000    2.918077
 1     1    0.939695   -0.525414    2.850000    2.908175    2.250000    2.850000
 spn 2 1    0.000000   -0.377045    2.850000    2.889339    2.250000    2.850000
 2     0    0.000001   -0.858896    3.147584    3.129995    3.147584    3.147584
 spn 2 0    0.000000   -1.249458    3.147584    3.106631    3.147584    3.147584
 3     0    0.000000   -0.833626    4.102416    4.094141    4.102416    4.102416
 spn 2 0    0.000000   -1.187697    4.102416    4.084347    4.102416    4.102416

 Harris energy:
 sumev=       -2.666741  val*vef=     -13.649558   sumtv=      10.982817
 sumec=      -40.255445  cor*vef=    -103.188197   ttcor=      62.932752
 rhoeps=      -9.192542     utot=    -139.216930    ehar=     -74.493904

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:    -18.875487      6.127728    -12.747759 sumev=   -2.666741   sumtv=   10.081017

 Kohn-Sham energy:
 sumtv=       10.081017  sumtc=        63.002032   ekin=       73.083050
 rhoep=       -8.772588   utot=      -138.744266   ehks=      -74.433804
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho:  sought 2 iter from file mixm; read 1.  RMS DQ=2.08e-2  last it=3.81e-2
 mixrho: (warning) scr. and lin-mixed densities had 79289 and 86897 negative poi
 AMIX: nmix=1 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=1.04e-2
   tj:-0.88407
 unscreened rms difference:  smooth  0.019296   local  0.050145
   screened rms difference:  smooth  0.019570   local  0.050145   tot  0.020774

 iors  : write restart file (binary, mesh density) 

   it  2  of 10    ehf=       0.500996   ehk=       0.561096
 From last iter    ehf=       0.409337   ehk=       0.590055
 diffe(q)=  0.091660 (0.020774)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999997 ehf=.5009963 ehk=.5610956

 --- BNDFP:  begin iteration 3 of 10 ---

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.789334   -2.798117     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.057324  avg sphere pot= 0.011121  vconst= 0.057324
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.615135   charge     6.798117
 smvxcm (warning) mesh density negative at 160958 points:  rhomin=-4.44e-4
 smooth rhoeps =   -6.100132 (  -4.104622,  -1.995510)
         rhomu =   -7.981001 (  -5.668033,  -2.312968)
       avg vxc =   -0.126766 (  -0.134349,  -0.119183)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.687584   -4.690419  -13.378004    -15.189240   -7.055214  -22.244454

 local terms:     true           smooth         local
 rhoeps:        -8.897480      -6.056710      -2.840771
 rhomu:         -6.470434      -5.632478      -0.837957
 spin2:         -5.239586      -2.292395      -2.947191
 total:        -11.710021      -7.924873      -3.785148
 val*vef       -13.378004     -11.738845      -1.639158
 val chg:        3.088280       5.886397      -2.798117
 val mom:        0.971650       2.053651      -1.082002    core:   0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000007

 potpus  spin 1 : pnu = 2.922271 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.241766   -1.432594   -1.345014    0.177637      2.7754    0.827102
 1     -0.671844   -1.112552   -0.670067    0.153507     18.7775    0.834664
 2     -0.627670   -0.627670    1.536944    0.361807     16.5355    6.601124
 3     -0.547444   -0.547444    3.479492    0.413524     23.5480   11.672343

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.040133    7.737471    0.078744   -0.642103
 1      3.000000    1.000000    -5.887832    7.226181    0.126410   -0.603226
 2      3.000000    1.000000     6.000000   28.951712    0.488481   -0.089192
 3      3.000000    1.000000     9.000000   39.791828    0.567911   -0.057183

 potpus  spin 2 : pnu = 2.918077 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.116110   -1.320451   -1.224118    0.185886      2.7879    0.850937
 1     -0.545258   -1.037947   -0.543293    0.161392     18.9905    0.867339
 2     -0.576377   -0.576377    1.611942    0.364615     16.4600    6.730826
 3     -0.502096   -0.502096    3.542630    0.414846     23.5015   11.772163

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -11.397896    7.615396    0.084739   -0.620665
 1      3.000000    1.000000    -5.887832    7.156893    0.132909   -0.576778
 2      3.000000    1.000000     6.000000   29.219724    0.490027   -0.087885
 3      3.000000    1.000000     9.000000   39.942403    0.568602   -0.056835

 Energy terms:             smooth           local           total
   rhoval*vef            -22.331729         8.866431       -13.465299
   rhoval*ves             -3.787815        -6.360784       -10.148598
   psnuc*ves              15.018086      -282.793938      -267.775852
   utot                    5.615135      -144.577361      -138.962225
   rho*exc                -6.100132        -2.840771        -8.940902
   rho*vxc                -7.981001        -3.785148       -11.766149
   valence chg             5.798117        -2.798117         3.000000
   valence mag             2.110985        -1.082002         1.028983
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.531737
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.2416 -0.6219 -0.6219 -0.6219  0.2455  0.6072  0.6072  0.6072
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1147 -0.4961 -0.4961 -0.4961  0.2685  0.6389  0.6389  0.6389

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.622095;   3.000000 electrons
         Sum occ. bands:   -2.978330, incl. Bloechl correction: -0.000020
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.897244    6.366148   -3.468904      0.957682    1.925889   -0.968207
       contr. to mm extrapolated for r>rmt:   0.035865 est. true mm = 0.993548

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.974594   -1.239612    2.922271    2.924404    2.500000    2.924404
 spn 2 0    0.969781   -1.112490    2.918077    2.921491    2.500000    2.921491
 1     1    0.952869   -0.616291    2.850000    2.912583    2.250000    2.850000
 spn 2 1    0.000000   -0.489999    2.850000    2.907449    2.250000    2.850000
 2     0    0.000001   -0.877005    3.147584    3.128869    3.147584    3.147584
 spn 2 0    0.000000   -1.305446    3.147584    3.105591    3.147584    3.147584
 3     0    0.000000   -0.850571    4.102416    4.093546    4.102416    4.102416
 spn 2 0    0.000000   -0.850571    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -2.978330  val*vef=     -13.465299   sumtv=      10.486969
 sumec=      -40.655878  cor*vef=    -103.623252   ttcor=      62.967374
 rhoeps=      -8.940902     utot=    -138.962225    ehar=     -74.448784

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:    -25.724852     12.112041    -13.612810 sumev=   -2.978330   sumtv=   10.634481

 Kohn-Sham energy:
 sumtv=       10.634481  sumtc=        63.010522   ekin=       73.645003
 rhoep=       -8.847377   utot=      -139.240236   ehks=      -74.442609
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho:  sought 2 iter from file mixm; read 2.  RMS DQ=7.81e-3  last it=2.08e-2
 mixrho: (warning) scr. and lin-mixed densities had 73335 and 78033 negative poi
 AMIX: nmix=2 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=3.91e-3
   tj:-1.62977   0.52986
 unscreened rms difference:  smooth  0.008208   local  0.015840
   screened rms difference:  smooth  0.008261   local  0.015840   tot  0.007815

 iors  : write restart file (binary, mesh density) 

   it  3  of 10    ehf=       0.546116   ehk=       0.552291
 From last iter    ehf=       0.500996   ehk=       0.561096
 diffe(q)=  0.045119 (0.007815)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999997 ehf=.5461156 ehk=.5522907

 --- BNDFP:  begin iteration 4 of 10 ---

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.129282   -4.003202     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.101126  avg sphere pot= 0.010036  vconst= 0.101126
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      6.121084   charge     8.003202
 smvxcm (warning) mesh density negative at 115474 points:  rhomin=-2.54e-4
 smooth rhoeps =   -7.937048 (  -5.042891,  -2.894157)
         rhomu =  -10.392239 (  -6.924380,  -3.467859)
       avg vxc =   -0.146685 (  -0.153231,  -0.140139)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.724151   -4.925595  -13.649747    -18.891420  -10.650914  -29.542334

 local terms:     true           smooth         local
 rhoeps:        -8.817440      -7.905752      -0.911688
 rhomu:         -6.365961      -6.900671       0.534710
 spin2:         -5.240462      -3.451063      -1.789399
 total:        -11.606423     -10.351734      -1.254689
 val*vef       -13.649747     -13.869937       0.220190
 val chg:        2.902644       6.905846      -4.003202
 val mom:        0.862716       1.968531      -1.105815    core:   0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000008

 potpus  spin 1 : pnu = 2.924404 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.256663   -1.442009   -1.358226    0.174028      2.7664    0.819842
 1     -0.685725   -1.103557   -0.684019    0.150398     18.5474    0.825603
 2     -0.604535   -0.604535    1.541491    0.359643     16.5913    6.511007
 3     -0.519349   -0.519349    3.491475    0.412339     23.5887   11.588598

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.393692    7.828873    0.076089   -0.649887
 1      3.000000    1.000000    -5.887832    7.303683    0.123844   -0.612108
 2      3.000000    1.000000     6.000000   28.757186    0.487200   -0.090191
 3      3.000000    1.000000     9.000000   39.661164    0.567263   -0.057492

 potpus  spin 2 : pnu = 2.921491 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.138829   -1.338586   -1.246349    0.182143      2.7802    0.842399
 1     -0.567641   -1.035066   -0.565755    0.158129     18.7687    0.857714
 2     -0.561715   -0.561715    1.607739    0.362427     16.5158    6.637677
 3     -0.482429   -0.482429    3.545866    0.413639     23.5428   11.685981

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -11.915744    7.690194    0.081089   -0.628990
 1      3.000000    1.000000    -5.887832    7.229096    0.130216   -0.585466
 2      3.000000    1.000000     6.000000   29.020977    0.488736   -0.088877
 3      3.000000    1.000000     9.000000   39.808681    0.567944   -0.057148

 Energy terms:             smooth           local           total
   rhoval*vef            -29.583389        15.892553       -13.690835
   rhoval*ves             -3.417660        -7.040020       -10.457680
   psnuc*ves              15.659828      -283.683140      -268.023311
   utot                    6.121084      -145.361580      -139.240496
   rho*exc                -7.937048        -0.911688        -8.848736
   rho*vxc               -10.392239        -1.254689       -11.646928
   valence chg             7.003202        -4.003202         3.000000
   valence mag             1.993735        -1.105815         0.887920
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.622095
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.2573 -0.6355 -0.6355 -0.6355  0.2667  0.6245  0.6245  0.6245
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1393 -0.5176 -0.5176 -0.5176  0.2843  0.6494  0.6494  0.6494

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.635665;   3.000000 electrons
         Sum occ. bands:   -3.032151, incl. Bloechl correction: -0.000019
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905383    6.701348   -3.795966      0.961295    1.927862   -0.966567
       contr. to mm extrapolated for r>rmt:   0.033138 est. true mm = 0.994433

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.976363   -1.255510    2.924404    2.925567    2.500000    2.925567
 spn 2 0    0.972043   -1.137380    2.921491    2.922876    2.500000    2.922876
 1     1    0.956976   -0.630539    2.850000    2.914260    2.250000    2.850000
 spn 2 1    0.000000   -0.502974    2.850000    2.917467    2.250000    2.850000
 2     0    0.000001   -0.862141    3.147584    3.128269    3.147584    3.147584
 spn 2 0    0.000000   -1.302766    3.147584    3.105023    3.147584    3.147584
 3     0    0.000000   -0.833223    4.102416    4.093250    4.102416    4.102416
 spn 2 0    0.000000   -0.833223    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -3.032151  val*vef=     -13.690835   sumtv=      10.658684
 sumec=      -40.648068  cor*vef=    -103.637018   ttcor=      62.988950
 rhoeps=      -8.848736     utot=    -139.240496    ehar=     -74.441598

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:    -28.966571     15.159185    -13.807385 sumev=   -3.032151   sumtv=   10.775234

 Kohn-Sham energy:
 sumtv=       10.775234  sumtc=        62.965499   ekin=       73.740733
 rhoep=       -8.865854   utot=      -139.317671   ehks=      -74.442792
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.75e-3  last it=7.81e-3
 mixrho: (warning) scr. and lin-mixed densities had 34189 and 53569 negative poi
 AMIX: nmix=2 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=8.77e-4
   tj:-0.32458   0.18866
 unscreened rms difference:  smooth  0.001097   local  0.005252
   screened rms difference:  smooth  0.001091   local  0.005252   tot  0.001754

 iors  : write restart file (binary, mesh density) 

   it  4  of 10    ehf=       0.553302   ehk=       0.552108
 From last iter    ehf=       0.546116   ehk=       0.552291
 diffe(q)=  0.007186 (0.001754)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999997 ehf=.5533021 ehk=.5521079

 --- BNDFP:  begin iteration 5 of 10 ---

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.058329   -3.751680     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.102813  avg sphere pot= 0.010764  vconst= 0.102813
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.747113   charge     7.751680
 smvxcm (warning) mesh density negative at 123196 points:  rhomin=-1.77e-4
 smooth rhoeps =   -7.594174 (  -4.833922,  -2.760252)
         rhomu =   -9.942615 (  -6.637359,  -3.305257)
       avg vxc =   -0.133395 (  -0.141408,  -0.125382)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.874507   -4.828913  -13.703420    -17.889614  -10.057369  -27.946983

 local terms:     true           smooth         local
 rhoeps:        -8.851187      -7.571479      -1.279709
 rhomu:         -6.434948      -6.618056       0.183108
 spin2:         -5.215892      -3.295195      -1.920698
 total:        -11.650841      -9.913251      -1.737589
 val*vef       -13.703420     -13.705142       0.001722
 val chg:        2.924597       6.676277      -3.751680
 val mom:        0.941696       1.936120      -0.994425    core:   0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000008

 potpus  spin 1 : pnu = 2.925567 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.241480   -1.430553   -1.345510    0.175244      2.7692    0.822792
 1     -0.671771   -1.095931   -0.670040    0.151515     18.5518    0.831028
 2     -0.600168   -0.600168    1.547022    0.359801     16.5857    6.520809
 3     -0.515159   -0.515159    3.495665    0.412350     23.5875   11.591996

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.594754    7.800803    0.075774   -0.647054
 1      3.000000    1.000000    -5.887832    7.302159    0.124763   -0.607668
 2      3.000000    1.000000     6.000000   28.776537    0.487250   -0.090105
 3      3.000000    1.000000     9.000000   39.665019    0.567249   -0.057486

 potpus  spin 2 : pnu = 2.922876 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.114123   -1.318272   -1.224567    0.183488      2.7832    0.845608
 1     -0.544635   -1.018999   -0.542720    0.159338     18.7594    0.863860
 2     -0.548373   -0.548373    1.622027    0.362562     16.5107    6.647134
 3     -0.469191   -0.469191    3.558807    0.413629     23.5422   11.688284

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.138526    7.660696    0.080668   -0.626106
 1      3.000000    1.000000    -5.887832    7.232167    0.131211   -0.580889
 2      3.000000    1.000000     6.000000   29.039020    0.488767   -0.088802
 3      3.000000    1.000000     9.000000   39.810585    0.567916   -0.057147

 Energy terms:             smooth           local           total
   rhoval*vef            -27.977473        14.243532       -13.733940
   rhoval*ves             -3.690235        -6.776790       -10.467024
   psnuc*ves              15.184461      -283.316349      -268.131888
   utot                    5.747113      -145.046569      -139.299456
   rho*exc                -7.594174        -1.279709        -8.873883
   rho*vxc                -9.942615        -1.737589       -11.680205
   valence chg             6.751680        -3.751680         3.000000
   valence mag             1.973229        -0.994425         0.978805
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.635665
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.2433 -0.6209 -0.6209 -0.6209  0.2814  0.6406  0.6406  0.6406
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1161 -0.4939 -0.4939 -0.4939  0.3055  0.6731  0.6731  0.6731
 Est Ef = -0.636 < evl(3)=-0.621 ... using qval=3.0, revise to -0.6209

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.621063;   3.000000 electrons
         Sum occ. bands:   -2.980443, incl. Bloechl correction: -0.000019
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.904129    6.774947   -3.870817      0.960777    1.957739   -0.996963
       contr. to mm extrapolated for r>rmt:   0.033623 est. true mm = 0.994399

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.976089   -1.241489    2.925567    2.925559    2.500000    2.925559
 spn 2 0    0.971676   -1.114162    2.922876    2.922839    2.500000    2.922839
 1     1    0.956363   -0.615771    2.850000    2.914258    2.250000    2.850000
 spn 2 1    0.000000   -0.478457    2.850000    2.917895    2.250000    2.850000
 2     0    0.000001   -0.857646    3.147584    3.128280    3.147584    3.147584
 spn 2 0    0.000000   -1.289060    3.147584    3.105038    3.147584    3.147584
 3     0    0.000000   -0.829260    4.102416    4.093243    4.102416    4.102416
 spn 2 0    0.000000   -0.829260    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -2.980443  val*vef=     -13.733940   sumtv=      10.753497
 sumec=      -40.576430  cor*vef=    -103.553655   ttcor=      62.977225
 rhoeps=      -8.873883     utot=    -139.299456    ehar=     -74.442616

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:    -28.601811     14.912850    -13.688961 sumev=   -2.980443   sumtv=   10.708518

 Kohn-Sham energy:
 sumtv=       10.708518  sumtc=        62.961444   ekin=       73.669962
 rhoep=       -8.857716   utot=      -139.255175   ehks=      -74.442929
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.36e-3  last it=1.75e-3
 mixrho: (warning) scr. and lin-mixed densities had 49969 and 58405 negative poi
 AMIX: nmix=2 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=6.78e-4
   tj: 0.13778  -0.15634
 unscreened rms difference:  smooth  0.001417   local  0.002824
   screened rms difference:  smooth  0.001424   local  0.002824   tot  0.001356

 iors  : write restart file (binary, mesh density) 

   it  5  of 10    ehf=       0.552284   ehk=       0.551971
 From last iter    ehf=       0.553302   ehk=       0.552108
 diffe(q)= -0.001018 (0.001356)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999997 ehf=.5522837 ehk=.5519706

 --- BNDFP:  begin iteration 6 of 10 ---

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.108456   -3.929374     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.110419  avg sphere pot= 0.010680  vconst= 0.110419
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.793050   charge     7.929374
 smvxcm (warning) mesh density negative at 62184 points:  rhomin=-1e-4
 smooth rhoeps =   -7.886451 (  -4.990338,  -2.896112)
         rhomu =  -10.326579 (  -6.847613,  -3.478966)
       avg vxc =   -0.153249 (  -0.162774,  -0.143724)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.892618   -4.828657  -13.721275    -18.462511  -10.574861  -29.037372

 local terms:     true           smooth         local
 rhoeps:        -8.842388      -7.867302      -0.975087
 rhomu:         -6.431071      -6.830851       0.399780
 spin2:         -5.208427      -3.470936      -1.737491
 total:        -11.639498     -10.301787      -1.337711
 val*vef       -13.721275     -14.069565       0.348290
 val chg:        2.900880       6.830254      -3.929374
 val mom:        0.943732       1.940407      -0.996675    core:   0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000008

 potpus  spin 1 : pnu = 2.925559 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.243261   -1.431076   -1.346671    0.174650      2.7671    0.822028
 1     -0.672996   -1.093503   -0.671275    0.151042     18.5076    0.829868
 2     -0.595242   -0.595242    1.548274    0.359381     16.5961    6.504359
 3     -0.509205   -0.509205    3.498123    0.412097     23.5958   11.575167

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.593275    7.821659    0.075591   -0.648006
 1      3.000000    1.000000    -5.887832    7.317366    0.124373   -0.608875
 2      3.000000    1.000000     6.000000   28.740566    0.486985   -0.090297
 3      3.000000    1.000000     9.000000   39.638260    0.567102   -0.057551

 potpus  spin 2 : pnu = 2.922839 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.115128   -1.317749   -1.224818    0.182804      2.7809    0.844747
 1     -0.544948   -1.014931   -0.543046    0.158798     18.7131    0.862477
 2     -0.542025   -0.542025    1.624363    0.362105     16.5218    6.628988
 3     -0.461722   -0.461722    3.562411    0.413350     23.5513   11.669657

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.132430    7.683375    0.080471   -0.627115
 1      3.000000    1.000000    -5.887832    7.247575    0.130765   -0.582189
 2      3.000000    1.000000     6.000000   28.999746    0.488479   -0.089006
 3      3.000000    1.000000     9.000000   39.781072    0.567754   -0.057218

 Energy terms:             smooth           local           total
   rhoval*vef            -29.055475        15.316065       -13.739410
   rhoval*ves             -3.650702        -6.830026       -10.480728
   psnuc*ves              15.236802      -283.368584      -268.131782
   utot                    5.793050      -145.099305      -139.306255
   rho*exc                -7.886451        -0.975087        -8.861538
   rho*vxc               -10.326579        -1.337711       -11.664290
   valence chg             6.929374        -3.929374         3.000000
   valence mag             1.975894        -0.996675         0.979219
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.621063
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.2448 -0.6219 -0.6219 -0.6219  0.2789  0.6436  0.6436  0.6436
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1167 -0.4939 -0.4939 -0.4939  0.3053  0.6779  0.6779  0.6779

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.622171;   3.000000 electrons
         Sum occ. bands:   -2.983623, incl. Bloechl correction: -0.000023
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905112    6.844412   -3.939301      0.960937    1.947651   -0.986713
       contr. to mm extrapolated for r>rmt:   0.033400 est. true mm = 0.994338

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.976391   -1.243036    2.925559    2.925784    2.500000    2.925784
 spn 2 0    0.972087   -1.114828    2.922839    2.923122    2.500000    2.923122
 1     1    0.956633   -0.617282    2.850000    2.914299    2.250000    2.850000
 spn 2 1    0.000000   -0.467949    2.850000    2.926867    2.250000    2.850000
 2     0    0.000001   -0.852931    3.147584    3.128250    3.147584    3.147584
 spn 2 0    0.000000   -1.285626    3.147584    3.104903    3.147584    3.147584
 3     0    0.000000   -0.824674    4.102416    4.093204    4.102416    4.102416
 spn 2 0    0.000000   -0.824674    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -2.983623  val*vef=     -13.739410   sumtv=      10.755787
 sumec=      -40.575865  cor*vef=    -103.545199   ttcor=      62.969334
 rhoeps=      -8.861538     utot=    -139.306255    ehar=     -74.442671

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:    -29.149976     15.439687    -13.710289 sumev=   -2.983623   sumtv=   10.726666

 Kohn-Sham energy:
 sumtv=       10.726666  sumtc=        62.959222   ekin=       73.685888
 rhoep=       -8.860241   utot=      -139.268552   ehks=      -74.442905
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=2.58e-4  last it=1.36e-3
 mixrho: (warning) scr. and lin-mixed densities had 24759 and 28937 negative poi
 AMIX: nmix=2 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=1.29e-4
   tj:-0.30421  -0.10845
 unscreened rms difference:  smooth  0.000268   local  0.000522
   screened rms difference:  smooth  0.000265   local  0.000522   tot  0.000258

 iors  : write restart file (binary, mesh density) 

   it  6  of 10    ehf=       0.552229   ehk=       0.551995
 From last iter    ehf=       0.552284   ehk=       0.551971
 diffe(q)= -0.000055 (0.000258)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999997 ehf=.552229 ehk=.5519946

 --- BNDFP:  begin iteration 7 of 10 ---

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.121482   -3.975551     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.113138  avg sphere pot= 0.010735  vconst= 0.113138
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.772272   charge     7.975551
 smvxcm (warning) mesh density negative at 36336 points:  rhomin=-3.9e-5
 smooth rhoeps =   -7.972604 (  -5.035372,  -2.937232)
         rhomu =  -10.439822 (  -6.907806,  -3.532016)
       avg vxc =   -0.162388 (  -0.172791,  -0.151985)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.907289   -4.805660  -13.712949    -18.597992  -10.717804  -29.315796

 local terms:     true           smooth         local
 rhoeps:        -8.842784      -7.955834      -0.886950
 rhomu:         -6.438527      -6.892472       0.453945
 spin2:         -5.201542      -3.525633      -1.675909
 total:        -11.640070     -10.418105      -1.221965
 val*vef       -13.712949     -14.202377       0.489428
 val chg:        2.897933       6.873484      -3.975551
 val mom:        0.957042       1.942694      -0.985652    core:   0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000008

 potpus  spin 1 : pnu = 2.925784 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.242707   -1.430572   -1.346258    0.174566      2.7668    0.821966
 1     -0.672522   -1.092326   -0.670803    0.150976     18.4929    0.829918
 2     -0.593472   -0.593472    1.548819    0.359245     16.5992    6.499615
 3     -0.507067   -0.507067    3.498913    0.412002     23.5988   11.569358

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.632836    7.825052    0.075419   -0.648113
 1      3.000000    1.000000    -5.887832    7.322467    0.124318   -0.608908
 2      3.000000    1.000000     6.000000   28.729975    0.486891   -0.090356
 3      3.000000    1.000000     9.000000   39.628778    0.567043   -0.057575

 potpus  spin 2 : pnu = 2.923122 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.112790   -1.315406   -1.222633    0.182664      2.7804    0.844566
 1     -0.542741   -1.011529   -0.540843    0.158676     18.6942    0.862368
 2     -0.537905   -0.537905    1.626879    0.361927     16.5258    6.622611
 3     -0.457125   -0.457125    3.565253    0.413226     23.5552   11.662016

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.178848    7.687921    0.080237   -0.627329
 1      3.000000    1.000000    -5.887832    7.253910    0.130664   -0.582356
 2      3.000000    1.000000     6.000000   28.985708    0.488357   -0.089083
 3      3.000000    1.000000     9.000000   39.768681    0.567676   -0.057249

 Energy terms:             smooth           local           total
   rhoval*vef            -29.328131        15.602814       -13.725316
   rhoval*ves             -3.661785        -6.804157       -10.465941
   psnuc*ves              15.206329      -283.325461      -268.119132
   utot                    5.772272      -145.064809      -139.292537
   rho*exc                -7.972604        -0.886950        -8.859555
   rho*vxc               -10.439822        -1.221965       -11.661786
   valence chg             6.975551        -3.975551         3.000000
   valence mag             1.980275        -0.985652         0.994623
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.622171
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.2443 -0.6212 -0.6212 -0.6212  0.2803  0.6482  0.6482  0.6482
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1145 -0.4913 -0.4913 -0.4913  0.3097  0.6857  0.6857  0.6857
 Est Ef = -0.622 < evl(3)=-0.621 ... using qval=3.0, revise to -0.6212

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.621507;   3.000000 electrons
         Sum occ. bands:   -2.980217, incl. Bloechl correction: -0.000022
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905588    6.912713   -4.007126      0.960964    1.942074   -0.981109
       contr. to mm extrapolated for r>rmt:   0.033368 est. true mm = 0.994332

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.976520   -1.242567    2.925784    2.925922    2.500000    2.925922
 spn 2 0    0.972312   -1.112565    2.923122    2.923334    2.500000    2.923334
 1     1    0.956755   -0.616770    2.850000    2.914375    2.250000    2.850000
 spn 2 1    0.000000   -0.460766    2.850000    2.930728    2.250000    2.850000
 2     0    0.000001   -0.851425    3.147584    3.128228    3.147584    3.147584
 spn 2 0    0.000000   -1.283993    3.147584    3.104798    3.147584    3.147584
 3     0    0.000000   -0.827170    4.102416    4.093083    4.102416    4.102416
 spn 2 0    0.000000   -0.827170    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -2.980217  val*vef=     -13.725316   sumtv=      10.745099
 sumec=      -40.574052  cor*vef=    -103.538330   ttcor=      62.964278
 rhoeps=      -8.859555     utot=    -139.292537    ehar=     -74.442714

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:    -29.494602     15.783893    -13.710709 sumev=   -2.980217   sumtv=   10.730492

 Kohn-Sham energy:
 sumtv=       10.730492  sumtc=        62.960882   ekin=       73.691374
 rhoep=       -8.861021   utot=      -139.273254   ehks=      -74.442900
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=3.28e-4  last it=2.58e-4
 mixrho: (warning) scr. and lin-mixed densities had 12413 and 16127 negative poi
 AMIX: nmix=2 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=1.64e-4
   tj: 0.16655  -0.28847
 unscreened rms difference:  smooth  0.000339   local  0.000708
   screened rms difference:  smooth  0.000335   local  0.000708   tot  0.000328

 iors  : write restart file (binary, mesh density) 

   it  7  of 10    ehf=       0.552186   ehk=       0.552000
 From last iter    ehf=       0.552229   ehk=       0.551995
 diffe(q)= -0.000043 (0.000328)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999997 ehf=.5521864 ehk=.552

 --- BNDFP:  begin iteration 8 of 10 ---

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.137913   -4.033795     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.114530  avg sphere pot= 0.010712  vconst= 0.114530
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1  0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.785249   charge     8.033795
 smvxcm (warning) mesh density negative at 1592 points:  rhomin=-2.54e-6
 smooth rhoeps =   -8.071874 (  -5.085983,  -2.985892)
         rhomu =  -10.570172 (  -6.975150,  -3.595023)
       avg vxc =   -0.174936 (  -0.186153,  -0.163719)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.907681   -4.799030  -13.706711    -18.775103  -10.898371  -29.673474

 local terms:     true           smooth         local
 rhoeps:        -8.841566      -8.055900      -0.785667
 rhomu:         -6.439172      -6.960309       0.521137
 spin2:         -5.199318      -3.589168      -1.610151
 total:        -11.638491     -10.549477      -1.089014
 val*vef       -13.706711     -14.324621       0.617910
 val chg:        2.895488       6.929284      -4.033795
 val mom:        0.960133       1.941347      -0.981213    core:   0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000008

 potpus  spin 1 : pnu = 2.925922 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.242959   -1.430778   -1.346556    0.174478      2.7666    0.821780
 1     -0.672837   -1.092015   -0.671120    0.150896     18.4850    0.829721
 2     -0.592714   -0.592714    1.548884    0.359167     16.6011    6.496726
 3     -0.506107   -0.506107    3.499161    0.411951     23.6004   11.566158

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.657386    7.827264    0.075298   -0.648306
 1      3.000000    1.000000    -5.887832    7.325206    0.124252   -0.609106
 2      3.000000    1.000000     6.000000   28.723625    0.486840   -0.090391
 3      3.000000    1.000000     9.000000   39.623622    0.567012   -0.057588

 potpus  spin 2 : pnu = 2.923334 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.112533   -1.315108   -1.222465    0.182542      2.7802    0.844266
 1     -0.542610   -1.010469   -0.540714    0.158560     18.6846    0.862044
 2     -0.536342   -0.536342    1.627550    0.361827     16.5281    6.618831
 3     -0.455304   -0.455304    3.566153    0.413161     23.5572   11.657867

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.213849    7.690231    0.080050   -0.627610
 1      3.000000    1.000000    -5.887832    7.257135    0.130568   -0.582641
 2      3.000000    1.000000     6.000000   28.977505    0.488291   -0.089127
 3      3.000000    1.000000     9.000000   39.762032    0.567636   -0.057266

 Energy terms:             smooth           local           total
   rhoval*vef            -29.683250        15.966729       -13.716521
   rhoval*ves             -3.649728        -6.808347       -10.458075
   psnuc*ves              15.220226      -283.329397      -268.109171
   utot                    5.785249      -145.068872      -139.283623
   rho*exc                -8.071874        -0.785667        -8.857541
   rho*vxc               -10.570172        -1.089014       -11.659186
   valence chg             7.033795        -4.033795         3.000000
   valence mag             1.979524        -0.981213         0.998310
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.621507
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.2446 -0.6214 -0.6214 -0.6214  0.2791  0.6481  0.6481  0.6481
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1143 -0.4911 -0.4911 -0.4911  0.3102  0.6875  0.6875  0.6875
 Est Ef = -0.622 < evl(3)=-0.621 ... using qval=3.0, revise to -0.6214

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.621760;   3.000000 electrons
         Sum occ. bands:   -2.980622, incl. Bloechl correction: -0.000022
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905810    6.943036   -4.037227      0.960913    1.933718   -0.972804
       contr. to mm extrapolated for r>rmt:   0.033383 est. true mm = 0.994296

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.976596   -1.242889    2.925922    2.925992    2.500000    2.925992
 spn 2 0    0.972448   -1.112409    2.923334    2.923450    2.500000    2.923450
 1     1    0.956765   -0.617155    2.850000    2.914362    2.250000    2.850000
 spn 2 1    0.000000   -0.459408    2.850000    2.931726    2.250000    2.850000
 2     0    0.000001   -0.850789    3.147584    3.128217    3.147584    3.147584
 spn 2 0    0.000000   -1.283915    3.147584    3.104737    3.147584    3.147584
 3     0    0.000000   -0.827129    4.102416    4.093058    4.102416    4.102416
 spn 2 0    0.000000   -0.827129    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -2.980622  val*vef=     -13.716521   sumtv=      10.735899
 sumec=      -40.576064  cor*vef=    -103.538644   ttcor=      62.962580
 rhoeps=      -8.857541     utot=    -139.283623    ehar=     -74.442685

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:    -29.705953     15.990573    -13.715380 sumev=   -2.980622   sumtv=   10.734758

 Kohn-Sham energy:
 sumtv=       10.734758  sumtc=        62.961819   ekin=       73.696577
 rhoep=       -8.861651   utot=      -139.277816   ehks=      -74.442890
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=9.00e-5  last it=3.28e-4
 AMIX: nmix=2 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=4.5e-5
   tj:-0.10896   0.11307
 unscreened rms difference:  smooth  0.000056   local  0.000205
   screened rms difference:  smooth  0.000043   local  0.000205   tot  0.000090

 iors  : write restart file (binary, mesh density) 

   it  8  of 10    ehf=       0.552215   ehk=       0.552010
 From last iter    ehf=       0.552186   ehk=       0.552000
 diffe(q)=  0.000029 (0.000090)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999997 ehf=.5522152 ehk=.55201

 --- BNDFP:  begin iteration 9 of 10 ---

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.136527   -4.028884     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.113527  avg sphere pot= 0.010706  vconst= 0.113527
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.786643   charge     8.028884
 smvxcm (warning) mesh density negative at 5838 points:  rhomin=-4.22e-6
 smooth rhoeps =   -8.063164 (  -5.078934,  -2.984230)
         rhomu =  -10.558673 (  -6.965084,  -3.593589)
       avg vxc =   -0.168977 (  -0.180720,  -0.157234)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.903274   -4.798810  -13.702084    -18.749204  -10.893279  -29.642483

 local terms:     true           smooth         local
 rhoeps:        -8.842521      -8.047266      -0.795254
 rhomu:         -6.439386      -6.950252       0.510866
 spin2:         -5.200323      -3.587829      -1.612494
 total:        -11.639709     -10.538081      -1.101628
 val*vef       -13.702084     -14.310298       0.608214
 val chg:        2.898947       6.927832      -4.028884
 val mom:        0.959763       1.937738      -0.977975    core:   0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000008

 potpus  spin 1 : pnu = 2.925992 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.242744   -1.430809   -1.346495    0.174563      2.7669    0.821836
 1     -0.672772   -1.092434   -0.671053    0.150957     18.4913    0.829854
 2     -0.593476   -0.593476    1.548638    0.359226     16.5996    6.499031
 3     -0.507014   -0.507014    3.498744    0.411987     23.5993   11.568522

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.669652    7.823638    0.075278   -0.648210
 1      3.000000    1.000000    -5.887832    7.323024    0.124302   -0.608959
 2      3.000000    1.000000     6.000000   28.728677    0.486877   -0.090364
 3      3.000000    1.000000     9.000000   39.627380    0.567032   -0.057579

 potpus  spin 2 : pnu = 2.923450 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.112407   -1.315317   -1.222560    0.182640      2.7807    0.844294
 1     -0.542704   -1.011126   -0.540806    0.158626     18.6914    0.862183
 2     -0.537308   -0.537308    1.627150    0.361891     16.5266    6.621364
 3     -0.456428   -0.456428    3.565579    0.413201     23.5559   11.660505

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.233148    7.685818    0.080000   -0.627534
 1      3.000000    1.000000    -5.887832    7.254865    0.130623   -0.582496
 2      3.000000    1.000000     6.000000   28.983006    0.488332   -0.089098
 3      3.000000    1.000000     9.000000   39.766220    0.567659   -0.057256

 Energy terms:             smooth           local           total
   rhoval*vef            -29.653300        15.940367       -13.712933
   rhoval*ves             -3.648946        -6.805268       -10.454214
   psnuc*ves              15.222233      -283.326015      -268.103783
   utot                    5.786643      -145.065642      -139.278998
   rho*exc                -8.063164        -0.795254        -8.858419
   rho*vxc               -10.558673        -1.101628       -11.660302
   valence chg             7.028884        -4.028884         3.000000
   valence mag             1.976252        -0.977975         0.998277
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.62176
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.2445 -0.6214 -0.6214 -0.6214  0.2798  0.6479  0.6479  0.6479
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1143 -0.4912 -0.4912 -0.4912  0.3113  0.6875  0.6875  0.6875
 Est Ef = -0.622 < evl(3)=-0.621 ... using qval=3.0, revise to -0.6214

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.621706;   3.000000 electrons
         Sum occ. bands:   -2.980506, incl. Bloechl correction: -0.000022
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905727    6.938691   -4.032964      0.960920    1.937475   -0.976555
       contr. to mm extrapolated for r>rmt:   0.033380 est. true mm = 0.994300

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.976557   -1.242773    2.925992    2.925963    2.500000    2.925963
 spn 2 0    0.972403   -1.112438    2.923450    2.923421    2.500000    2.923421
 1     1    0.956766   -0.617021    2.850000    2.914386    2.250000    2.850000
 spn 2 1    0.000000   -0.459910    2.850000    2.931375    2.250000    2.850000
 2     0    0.000001   -0.852122    3.147584    3.128183    3.147584    3.147584
 spn 2 0    0.000000   -1.284570    3.147584    3.104752    3.147584    3.147584
 3     0    0.000000   -0.829279    4.102416    4.093027    4.102416    4.102416
 spn 2 0    0.000000   -0.829279    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -2.980506  val*vef=     -13.712933   sumtv=      10.732427
 sumec=      -40.577170  cor*vef=    -103.539370   ttcor=      62.962200
 rhoeps=      -8.858419     utot=    -139.278998    ehar=     -74.442791

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:    -29.675325     15.961676    -13.713649 sumev=   -2.980506   sumtv=   10.733143

 Kohn-Sham energy:
 sumtv=       10.733143  sumtc=        62.962566   ekin=       73.695709
 rhoep=       -8.861454   utot=      -139.277151   ehks=      -74.442897
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=3.56e-5  last it=9.00e-5
 mixrho: (warning) scr. and lin-mixed densities had 0 and 459 negative points
 AMIX: nmix=2 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=1.78e-5
   tj:-0.19837  -0.08027
 unscreened rms difference:  smooth  0.000041   local  0.000090
   screened rms difference:  smooth  0.000032   local  0.000090   tot  0.000036

 iors  : write restart file (binary, mesh density) 

   it  9  of 10    ehf=       0.552109   ehk=       0.552003
 From last iter    ehf=       0.552215   ehk=       0.552010
 diffe(q)= -0.000106 (0.000036)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999997 ehf=.5521094 ehk=.5520032

 --- BNDFP:  begin iteration 10 of 10 ---

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.137742   -4.033191     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.112920  avg sphere pot= 0.010695  vconst= 0.112920
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.790553   charge     8.033191
 smooth rhoeps =   -8.070539 (  -5.082106,  -2.988433)
         rhomu =  -10.568319 (  -6.969112,  -3.599207)
       avg vxc =   -0.167516 (  -0.179952,  -0.155079)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.901192   -4.797213  -13.698405    -18.758926  -10.909115  -29.668041

 local terms:     true           smooth         local
 rhoeps:        -8.843470      -8.054791      -0.788679
 rhomu:         -6.440219      -6.954345       0.514126
 spin2:         -5.200714      -3.593577      -1.607137
 total:        -11.640933     -10.547923      -1.093010
 val*vef       -13.698405     -14.315223       0.616817
 val chg:        2.901531       6.934723      -4.033191
 val mom:        0.960413       1.937238      -0.976825    core:   0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000008

 potpus  spin 1 : pnu = 2.925963 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.242568   -1.430725   -1.346350    0.174620      2.7671    0.821938
 1     -0.672621   -1.092640   -0.670902    0.151005     18.4952    0.829983
 2     -0.593923   -0.593923    1.548516    0.359263     16.5986    6.500552
 3     -0.507549   -0.507549    3.498500    0.412008     23.5985   11.569990

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.664556    7.822010    0.075316   -0.648095
 1      3.000000    1.000000    -5.887832    7.321660    0.124342   -0.608827
 2      3.000000    1.000000     6.000000   28.731988    0.486900   -0.090347
 3      3.000000    1.000000     9.000000   39.629693    0.567044   -0.057573

 potpus  spin 2 : pnu = 2.923421 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.112225   -1.315234   -1.222412    0.182700      2.7808    0.844397
 1     -0.542556   -1.011369   -0.540657    0.158676     18.6953    0.862313
 2     -0.537756   -0.537756    1.627042    0.361930     16.5256    6.622958
 3     -0.456969   -0.456969    3.565349    0.413223     23.5552   11.662066

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.228324    7.684188    0.080040   -0.627423
 1      3.000000    1.000000    -5.887832    7.253534    0.130664   -0.582374
 2      3.000000    1.000000     6.000000   28.986449    0.488356   -0.089080
 3      3.000000    1.000000     9.000000   39.768674    0.567672   -0.057250

 Energy terms:             smooth           local           total
   rhoval*vef            -29.679358        15.969603       -13.709755
   rhoval*ves             -3.645334        -6.805185       -10.450519
   psnuc*ves              15.226440      -283.327524      -268.101084
   utot                    5.790553      -145.066355      -139.275801
   rho*exc                -8.070539        -0.788679        -8.859218
   rho*vxc               -10.568319        -1.093010       -11.661330
   valence chg             7.033191        -4.033191         3.000000
   valence mag             1.976106        -0.976825         0.999281
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.621706
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.2443 -0.6212 -0.6212 -0.6212  0.2803  0.6479  0.6479  0.6479
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1141 -0.4910 -0.4910 -0.4910  0.3127  0.6881  0.6881  0.6881
 Est Ef = -0.622 < evl(3)=-0.621 ... using qval=3.0, revise to -0.6212

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.621548;   3.000000 electrons
         Sum occ. bands:   -2.979968, incl. Bloechl correction: -0.000021
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905678    6.939139   -4.033461      0.960904    1.936963   -0.976059
       contr. to mm extrapolated for r>rmt:   0.033405 est. true mm = 0.994310

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.976536   -1.242580    2.925963    2.925952    2.500000    2.925952
 spn 2 0    0.972387   -1.112228    2.923421    2.923418    2.500000    2.923418
 1     1    0.956755   -0.616829    2.850000    2.914393    2.250000    2.850000
 spn 2 1    0.000000   -0.460297    2.850000    2.930940    2.250000    2.850000
 2     0    0.000001   -0.852551    3.147584    3.128186    3.147584    3.147584
 spn 2 0    0.000000   -1.285060    3.147584    3.104753    3.147584    3.147584
 3     0    0.000000   -0.829611    4.102416    4.093032    4.102416    4.102416
 spn 2 0    0.000000   -0.829611    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -2.979968  val*vef=     -13.709755   sumtv=      10.729787
 sumec=      -40.577486  cor*vef=    -103.539869   ttcor=      62.962383
 rhoeps=      -8.859218     utot=    -139.275801    ehar=     -74.442849

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:    -29.680051     15.968348    -13.711703 sumev=   -2.979968   sumtv=   10.731736

 Kohn-Sham energy:
 sumtv=       10.731736  sumtc=        62.963108   ekin=       73.694844
 rhoep=       -8.861285   utot=      -139.276460   ehks=      -74.442902
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=2.40e-5  last it=3.56e-5
 AMIX: nmix=2 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=1.2e-5
   tj: 0.46727  -0.15976
 unscreened rms difference:  smooth  0.000027   local  0.000045
   screened rms difference:  smooth  0.000024   local  0.000045   tot  0.000024

 iors  : write restart file (binary, mesh density) 

   it 10  of 10    ehf=       0.552051   ehk=       0.551998
 From last iter    ehf=       0.552109   ehk=       0.552003
 diffe(q)= -0.000059 (0.000024)    tol= 0.000010 (0.000500)   more=F
x zbak=1 mmom=.9999997 ehf=.5520506 ehk=.5519984
 Exit 0 LMF 
 wkinfo:  used 14885 K  workspace of 80000 K   in  68 K calls
