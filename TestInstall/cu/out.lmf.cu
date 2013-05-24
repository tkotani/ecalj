rdcmd:  lmfa --no-iactiv cu -vnk=8 -vbigbas=f
 -----------------------  START LMFA (80000K)  -----------------------
 ptenv() is called with EXT=cu
 ptenv() not supported, but continue.

 rdctrl: reset global max nl from 5 to 4

 LMFA:     alat = 6.798  nbas = 1  nspec = 1  vn 7.00(LMFA 7.0)  verb 31,30
 pot:      XC:BH

                Plat                                  Qlat
   0.000000   0.500000   0.500000       -1.000000   1.000000   1.000000
   0.500000   0.000000   0.500000        1.000000  -1.000000   1.000000
   0.500000   0.500000   0.000000        1.000000   1.000000  -1.000000
  Cell vol= 78.538660

 LATTC:  as= 2.000   tol=1.00E-12   alat= 6.79800   awald= 0.467
         r1=  2.252   nkd= 201      q1=  6.871   nkg= 331

 SGROUP: 1 symmetry operations from 0 generators
 SYMLAT: Bravais system is cubic with 48 symmetry operations.
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-5
 GROUPG: the following are sufficient to generate the space group:
         i*r3(1,1,-1) r4x
         i*r3(1,1,-1) r4x
 MKSYM:  found 48 space group operations ... includes inversion

conf:SPEC_ATOM= A : --- Table for atomic configuration ---
conf int(P)z = int(P) where P is replaced by PZ if it is semicore
conf:  isp  l  int(P) int(P)z    Qval    Qcore   CoreConf
conf:    1  0       4  4        1.000    6.000 => 1,2,3,
conf:    1  1       4  4        0.000   12.000 => 2,3,
conf:    1  2       3  3       10.000    0.000 => 
conf:    1  3       4  4        0.000    0.000 => 
conf:-----------------------------------------------------

 Species A:  Z=29  Qc=18  R=2.311271  Q=0
 mesh:   rmt=2.311271  rmax=48.805862  a=0.025  nr=393  nr(rmax)=515
  Pl=  4.5     4.5     3.5     4.5    
  Ql=  1.0     0.0     10.0    0.0    

  iter     qint         drho          vh0          rho0          vsum     beta
    1   29.000000   4.725E+03      145.0000    0.1442E+03      -58.2772   0.30
   51   29.000000   4.211E-05      274.8263    0.2631E+05     -130.7915   0.30


 sumev=-4.333254  etot=-3304.416258  eref=-3304.434500  diff= 0.018242

 Free-atom wavefunctions:
 valence:      eval       node at      max at       c.t.p.   rho(r>rmt)
   4s      -0.36411         0.890       2.256       3.582     0.643062
   4p      -0.06295         0.975       3.484       7.414     0.901829
   3d      -0.39691         0.000       0.600       3.429     0.056076
   4f       0.01948         0.000      35.393      48.806*    1.000000

 core:        ecore       node at      max at       c.t.p.   rho(r>rmt)
   1s    -649.07634         0.000       0.034       0.069     0.000000
   2s     -77.91382         0.070       0.197       0.308     0.000000
   2p     -67.32532         0.000       0.158       0.335     0.000000
   3s      -8.39248         0.288       0.614       0.895     0.000141
   3p      -5.29682         0.260       0.619       1.078     0.000727

 Optimise free-atom basis for species A, rmt=2.311271
 l  it    Rsm      Eh     stiffR   stiffE      Eval      Exact     Pnu    Ql
 0   9   2.311  -0.297     106.7    449.5   -0.36395  -0.36411    4.76   1.00
 ... rsm exceeded rmt*2/3 .. repeat with rsm=rmt
 0   5   1.541  -0.133     106.7   3896.8   -0.35292  -0.36411    4.76   1.00
 1  11   2.311  -0.100     147.4     45.6   -0.04993  -0.06295    4.56   0.00
 ... rsm exceeded rmt*2/3 .. repeat with rsm=rmt
 1   1   1.541  -0.100     147.4  -1923.8    0.05244  -0.06295    4.56   0.00
 2  27   0.962  -0.116     158.7    107.6   -0.39670  -0.39691    3.89  10.00
 eigenvalue sum:  exact  -4.33325    opt basis  -4.31995    error 0.01330

 tailsm: fit tails to 6 smoothed hankels, rmt= 2.31127, rsm= 1.15564
 E:    -1.00000    -2.00000    -4.00000    -6.00000    -9.00000    -15.0000
 C:    -0.07160    10.75053    -187.492    1222.023    -4717.79    21166.81
        r          rho         fit         diff
    2.311271    0.017797    0.017766    0.000031
    2.967767    0.005662    0.005658    0.000005
    3.810725    0.001517    0.001518   -0.000001
    4.893104    0.000305    0.000305    0.000000
    6.282906    0.000041    0.000041   -0.000001
    8.067448    0.000003    0.000003    0.000000
    q(fit):     1.203836    rms diff:   0.000016
    fit: r>rmt  1.203836   r<rmt  3.442816   qtot  4.646652
    rho: r>rmt  1.203836   r<rmt  9.796164   qtot 11.000000

 coretail: q=0.00392, rho(rmt)=0.00465.  Fit with Hankel e=-24.082  coeff=764.|
      r            rhoc          fit
    2.311271    0.02095279    0.02095279
    2.429779    0.01229068    0.01231367
    2.753317    0.00285262    0.00285190
    3.119934    0.00054243    0.00053465
    3.535366    0.00008235    0.00007888
    4.006112    0.00000969    0.00000887
    4.539536    0.00000085    0.00000073
    5.143985    0.00000005    0.00000004

  Not write mtopara.* when there is PZ.

 Sum of reference energies: -3304.4345
 Exit 0 LMFA 
 wkinfo:  used    94 K  workspace of 80000 K   in   0 K calls
rdcmd:  lmf  --no-iactiv cu -vnk=8 -vbigbas=f
 -----------------------  START LMF (80000K)  -----------------------
 ptenv() is called with EXT=cu
 ptenv() not supported, but continue.

 rdctrl: reset global max nl from 5 to 4

 LMF:      alat = 6.798  nbas = 1  nspec = 1  vn 7.00(LMF 7.0)  verb 31,30
 pot:      XC:BH
 bz:       metal(2), tetra, invit 

                Plat                                  Qlat
   0.000000   0.500000   0.500000       -1.000000   1.000000   1.000000
   0.500000   0.000000   0.500000        1.000000  -1.000000   1.000000
   0.500000   0.500000   0.000000        1.000000   1.000000  -1.000000
  Cell vol= 78.538660

 LATTC:  as= 2.000   tol=1.00E-12   alat= 6.79800   awald= 0.467
         r1=  2.252   nkd= 201      q1=  6.871   nkg= 331

 SGROUP: 1 symmetry operations from 0 generators
 SYMLAT: Bravais system is cubic with 48 symmetry operations.
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-5
 GROUPG: the following are sufficient to generate the space group:
         i*r3(1,1,-1) r4x
         i*r3(1,1,-1) r4x
 MKSYM:  found 48 space group operations ... includes inversion
 
 BZMESH:  60 irreducible QP from 512 ( 8 8 8 )  shift= T T T
 TETIRR: sorting 3072 tetrahedra ...
 264 inequivalent ones found

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv  kmxv foca   rfoca
 A        2.311  0.925    3    4         3  0.578  1.156    15    1   0.925

 MSHSIZ: mesh has 11 x 11 x 11 divisions; length 0.437, 0.437, 0.437
         generated from gmax = 9.0 a.u. : 941 vectors of 1331 (70%)

 gvlist: cutoff radius   9.000 gives    941   recips of max   1331
 SGVSYM: 41 symmetry stars found for 941 reciprocal lattice vectors
 

 Makidx:  hamiltonian dimensions Low, Int, High, Negl: 18 0 7 7
 kappa   Low   Int   High  L+I  L+I+H  Neglected
   1       9     0     7     9    16       0
   2       9     0     0     9     9       7
  all     18     0     7    18    25       7
 suham :  16 augmentation channels, 16 local potential channels  Maximum lmxa=3

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol= 1.00E-06
 spec      l    rsm    eh     gmax    last term    cutoff
  A        0*   2.50  -0.01   2.974    2.30E-05      27 
  A        1    2.50  -0.01   3.093    1.29E-06      51 
  A        2    1.00  -0.01   8.508    1.16E-06     869 

 iors  : read restart file (binary, mesh density) 
 iors  : empty file ... nothing read

 rdovfa: read and overlap free-atom densities (mesh density) ...
 rdovfa: expected A,       read A        with rmt=  2.3113  mesh   393  0.025

 ovlpfa: overlap smooth part of FA densities
 site   1  spec  1  pos  0.0000  0.0000  0.0000  Qsmooth 4.646654
 total smooth Q = 4.646654

 Free atom and overlapped crystal site charges:
   ib    true(FA)    smooth(FA)  true(OV)    smooth(OV)    local
    1    9.796164    3.442818   10.275300    3.921954    6.353346

 Smooth charge on mesh:            4.646654
 Sum of local charges:             6.353346
 Total valence charge:            11.000000
 Sum of core charges:             18.000000
 Sum of nuclear charges:         -29.000000
 Homogeneous background:           0.000000
 Deviation from neutrality:        0.000000

 --- BNDFP:  begin iteration 1 of 12 ---

 avg es pot at rmt= 0.554993  avg sphere pot= 0.633521  vconst=-0.554993

 smooth rhoves     11.022237   charge     4.646654
 smooth rhoeps =   -3.843801   rhomu =   -5.010456  avg vxc =   -0.851784 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000012

 Energy terms:             smooth           local           total
   rhoval*vef            -12.156988      -177.336818      -189.493806
   rhoval*ves            -46.689418      -115.324370      -162.013788
   psnuc*ves              68.733892    -12976.662455    -12907.928563
   utot                   11.022237     -6545.993412     -6534.971175
   rho*exc                -3.843801      -126.414296      -130.258097
   rho*vxc                -5.010456      -167.409313      -172.419769
   valence chg             4.646654         6.353346        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Incompatible or missing qp weights file ...

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.6627 -0.0591 -0.0536 -0.0536  0.0167  0.0167  1.7979  1.9471  1.9471
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.3283 -0.1257 -0.0727 -0.0255  0.0508  0.0920  0.6928  1.2801  1.4813
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4087 -0.1094 -0.0570 -0.0368  0.0286  0.0633  0.7989  1.4016  1.7830
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.3857 -0.1324 -0.0620  0.0030  0.0284  0.0660  0.9994  1.1878  1.4473
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.2227 -0.1325 -0.1100  0.0005  0.0467  0.2005  0.6490  0.8702  1.1936
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4328 -0.0995 -0.0481 -0.0481  0.0391  0.0391  0.7797  1.7311  1.7311

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.144577;  11.000000 electrons
         Sum occ. bands:   -0.853241, incl. Bloechl correction: -0.006586

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.6627 -0.0591 -0.0536 -0.0536  0.0167  0.0167  1.7979  1.9471  1.9471
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.3283 -0.1257 -0.0727 -0.0255  0.0508  0.0920  0.6928  1.2801  1.4813
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4087 -0.1094 -0.0570 -0.0368  0.0286  0.0633  0.7989  1.4016  1.7830
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.3857 -0.1324 -0.0620  0.0030  0.0284  0.0660  0.9994  1.1878  1.4473
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.2227 -0.1325 -0.1100  0.0005  0.0467  0.2005  0.6490  0.8702  1.1936
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4328 -0.0995 -0.0481 -0.0481  0.0391  0.0391  0.7797  1.7311  1.7311

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.144577;  11.000000 electrons
         Sum occ. bands:   -0.853241, incl. Bloechl correction: -0.006586

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1    9.927753    3.113493    6.814260

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.509191   -0.293352    4.650000    4.686057    4.500000    4.686057
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.531005   -0.280671    4.340000    4.364506    4.250000    4.364506
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.679694   -0.050852    3.870000    3.857413    3.147584    3.857413
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.023532   -0.620420    4.110000    4.109925    4.102416    4.110000

 Harris energy:
 sumev=       -0.853241  val*vef=    -189.493806   sumtv=     188.640564
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -130.258097     utot=   -6534.971175    ehar=   -3304.832068

 srhov:     -6.360828   -168.222982   -174.583810 sumev=   -0.853241   sumtv=  173.730568

 Kohn-Sham energy:
 sumtv=      173.730568  sumtc=      3171.756639   ekin=     3345.487208
 rhoep=     -128.641260   utot=     -6521.364242   ehks=    -3304.518293
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 0.  RMS DQ=3.93e-2
 charges:       old           new         screened      rms diff       lin mix
 smooth       4.646654      4.185740      4.185740      0.038395      4.185740
 site    1    6.353346      6.814260      6.814260      0.019256      6.814260
 AMIX: nmix=0 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=3.93e-2
 unscreened rms difference:  smooth  0.045849   local  0.019256   tot
   screened rms difference:  smooth  0.038395   local  0.019256   tot  0.039302

 iors  : write restart file (binary, mesh density) 

   it  1  of 12    ehf=      -0.397568   ehk=      -0.083793
h nk=8 bigbas=0 ehf=-.397568 ehk=-.0837934

 --- BNDFP:  begin iteration 2 of 12 ---

 avg es pot at rmt= 0.633049  avg sphere pot= 0.653661  vconst=-0.633049

 smooth rhoves     12.286140   charge     4.185740
 smooth rhoeps =   -3.109206   rhomu =   -4.047533  avg vxc =   -0.858397 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000645

 Energy terms:             smooth           local           total
   rhoval*vef             -7.713657      -175.301725      -183.015382
   rhoval*ves            -48.957947      -107.937688      -156.895635
   psnuc*ves              73.530227    -12964.307484    -12890.777257
   utot                   12.286140     -6536.122586     -6523.836446
   rho*exc                -3.109206      -125.879113      -128.988318
   rho*vxc                -4.047533      -166.689186      -170.736719
   valence chg             4.185740         6.814260        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=0.144577

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7748 -0.6825 -0.6802 -0.6802 -0.6430 -0.6430  1.6026  1.7558  1.7558
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.7179 -0.6847 -0.6654 -0.6541 -0.6178 -0.3403  0.4593  1.0915  1.2747
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.7184 -0.6784 -0.6695 -0.6492 -0.6265 -0.4571  0.5708  1.2078  1.5794
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.7146 -0.6988 -0.6593 -0.6431 -0.6231 -0.4289  0.7940  1.0031  1.2312
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.6997 -0.6980 -0.6768 -0.6419 -0.6231 -0.1362  0.4145  0.6649  0.9900
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.7269 -0.6758 -0.6758 -0.6287 -0.6287 -0.4938  0.5479  1.5355  1.5355

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.191164;  11.000000 electrons
         Sum occ. bands:   -7.093537, incl. Bloechl correction: -0.013294

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.391065    1.889160    8.501905

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.417713   -0.500906    4.686057    4.669917    4.500000    4.669917
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.253604   -0.580929    4.364506    4.321916    4.250000    4.321916
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.859413   -0.665867    3.857413    3.897507    3.147584    3.897507
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.010349   -0.855436    4.110000    4.106422    4.102416    4.110000

 Harris energy:
 sumev=       -7.093537  val*vef=    -183.015382   sumtv=     175.921845
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -128.988318     utot=   -6523.836446    ehar=   -3305.146280

 srhov:     -4.200311   -222.071601   -226.271912 sumev=   -7.093537   sumtv=  219.178375

 Kohn-Sham energy:
 sumtv=      219.178375  sumtc=      3171.756639   ekin=     3390.935015
 rhoep=     -132.663784   utot=     -6562.033736   ehks=    -3303.762505
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 1.  RMS DQ=9.93e-2  last it=3.93e-2
 mixrho: (warning) scr. and lin-mixed densities had 43 and 43 negative points
 charges:       old           new         screened      rms diff       lin mix
 smooth       4.185740      2.498095      2.498095      0.032750      2.498095
 site    1    6.814260      8.501905      8.501905      0.064870      8.501905
 AMIX: nmix=1 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=9.93e-2
   tj: 0.82102
 unscreened rms difference:  smooth  0.024244   local  0.064870   tot
   screened rms difference:  smooth  0.032750   local  0.064870   tot  0.099256

 iors  : write restart file (binary, mesh density) 

   it  2  of 12    ehf=      -0.711780   ehk=       0.671995
 From last iter    ehf=      -0.397568   ehk=      -0.083793
 diffe(q)= -0.314212 (0.099256)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf=-.7117796 ehk=.6719954

 --- BNDFP:  begin iteration 3 of 12 ---

 avg es pot at rmt= 0.604916  avg sphere pot= 0.655635  vconst=-0.604916

 smooth rhoves     10.956051   charge     3.883688
 smooth rhoeps =   -2.817571   rhomu =   -3.666846  avg vxc =   -0.843501 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000597

 Energy terms:             smooth           local           total
   rhoval*vef             -6.576790      -181.158615      -187.735405
   rhoval*ves            -47.654984      -113.403490      -161.058475
   psnuc*ves              69.567087    -12968.686260    -12899.119174
   utot                   10.956051     -6541.044875     -6530.088824
   rho*exc                -2.817571      -126.698711      -129.516282
   rho*vxc                -3.666846      -167.770342      -171.437189
   valence chg             3.883688         7.116312        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.191164

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7459 -0.4492 -0.4459 -0.4459 -0.3978 -0.3978  1.6618  1.8140  1.8140
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.5240 -0.4584 -0.4444 -0.4218 -0.3681 -0.2472  0.5246  1.1465  1.3362
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.5486 -0.4463 -0.4409 -0.4311 -0.3824 -0.3326  0.6354  1.2648  1.6406
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.5296 -0.4910 -0.4336 -0.4016 -0.3771 -0.3164  0.8534  1.0565  1.2958
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.4870 -0.4768 -0.4626 -0.4012 -0.3741 -0.0587  0.4793  0.7235  1.0493
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5673 -0.4409 -0.4409 -0.4023 -0.3806 -0.3806  0.6134  1.5939  1.5939

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.121671;  11.000000 electrons
         Sum occ. bands:   -4.729158, incl. Bloechl correction: -0.011779

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.274241    2.258615    8.015627

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.444930   -0.471542    4.669917    4.661876    4.500000    4.661876
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.298579   -0.559380    4.321916    4.315866    4.250000    4.315866
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.922377   -0.432983    3.897507    3.884967    3.147584    3.884967
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.013587   -0.809393    4.110000    4.107081    4.102416    4.110000

 Harris energy:
 sumev=       -4.729158  val*vef=    -187.735405   sumtv=     183.006247
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.516282     utot=   -6530.088824    ehar=   -3304.842220

 srhov:     -4.839716   -203.035058   -207.874774 sumev=   -4.729158   sumtv=  203.145616

 Kohn-Sham energy:
 sumtv=      203.145616  sumtc=      3171.756639   ekin=     3374.902255
 rhoep=     -131.316411   utot=     -6548.087137   ehks=    -3304.501293
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 2.  RMS DQ=4.81e-2  last it=9.93e-2
 mixrho: (warning) scr. and lin-mixed densities had 13 and 13 negative points
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.883688      2.984374      2.984374      0.016786      2.984374
 site    1    7.116312      8.015627      8.015627      0.032802      8.015627
 AMIX: nmix=2 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=4.81e-2
   tj:-1.03939  -0.10632
 add q= -0.000001 to preserve neutrality
 unscreened rms difference:  smooth  0.012621   local  0.032802   tot
   screened rms difference:  smooth  0.016786   local  0.032802   tot  0.048087

 iors  : write restart file (binary, mesh density) 

   it  3  of 12    ehf=      -0.407720   ehk=      -0.066793
 From last iter    ehf=      -0.711780   ehk=       0.671995
 diffe(q)=  0.304060 (0.048087)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf=-.4077199 ehk=-.0667929

 --- BNDFP:  begin iteration 4 of 12 ---

 avg es pot at rmt= 0.545124  avg sphere pot= 0.668904  vconst=-0.545124

 smooth rhoves      8.670801   charge     3.362076
 smooth rhoeps =   -2.347228   rhomu =   -3.053301  avg vxc =   -0.813360 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000549

 Energy terms:             smooth           local           total
   rhoval*vef             -4.919950      -185.951567      -190.871517
   rhoval*ves            -44.677352      -118.835010      -163.512362
   psnuc*ves              62.018954    -12969.563968    -12907.545014
   utot                    8.670801     -6544.199489     -6535.528688
   rho*exc                -2.347228      -127.805291      -130.152519
   rho*vxc                -3.053301      -169.227280      -172.280580
   valence chg             3.362076         7.637924        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.121671

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.6943 -0.1479 -0.1429 -0.1429 -0.0767 -0.0767  1.7567  1.9062  1.9062
 Est Ef = -0.122 < evl(5)=-0.077 ... using qval=11.0, revise to -0.0767
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.3699 -0.2013 -0.1586 -0.1153 -0.0430  0.0028  0.6426  1.2386  1.4371
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4446 -0.1874 -0.1450 -0.1264 -0.0644 -0.0306  0.7502  1.3595  1.7391
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4219 -0.2145 -0.1469 -0.0883 -0.0631 -0.0275  0.9552  1.1468  1.4017
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.2764 -0.2096 -0.1913 -0.0899 -0.0480  0.1255  0.5981  0.8253  1.1496
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4683 -0.1753 -0.1373 -0.1373 -0.0554 -0.0554  0.7303  1.6886  1.6886

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.066107;  11.000000 electrons
         Sum occ. bands:   -1.734343, incl. Bloechl correction: -0.007557

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.005947    2.935569    7.070378

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.497310   -0.352713    4.661876    4.677994    4.500000    4.677994
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.470940   -0.365528    4.315866    4.349210    4.250000    4.349210
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.758350   -0.137620    3.884967    3.863156    3.147584    3.863156
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021110   -0.677000    4.110000    4.109222    4.102416    4.110000

 Harris energy:
 sumev=       -1.734343  val*vef=    -190.871517   sumtv=     189.137174
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -130.152519     utot=   -6535.528688    ehar=   -3304.787394

 srhov:     -5.820552   -174.778370   -180.598922 sumev=   -1.734343   sumtv=  178.864579

 Kohn-Sham energy:
 sumtv=      178.864579  sumtc=      3171.756639   ekin=     3350.621218
 rhoep=     -129.135489   utot=     -6526.156007   ehks=    -3304.670278
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 3.  RMS DQ=2.66e-2  last it=4.81e-2
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.362076      3.929623      3.929623      0.009672      3.929623
 site    1    7.637924      7.070378      7.070378      0.019610      7.070378
 AMIX: nmix=3 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=2.66e-2
   tj: 0.76324  -0.24309  -0.00278
 unscreened rms difference:  smooth  0.007487   local  0.019610   tot
   screened rms difference:  smooth  0.009672   local  0.019610   tot  0.026573

 iors  : write restart file (binary, mesh density) 

   it  4  of 12    ehf=      -0.352894   ehk=      -0.235778
 From last iter    ehf=      -0.407720   ehk=      -0.066793
 diffe(q)=  0.054826 (0.026573)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf=-.3528937 ehk=-.2357778

 --- BNDFP:  begin iteration 5 of 12 ---

 avg es pot at rmt= 0.570014  avg sphere pot= 0.662710  vconst=-0.570014

 smooth rhoves      9.528839   charge     3.555443
 smooth rhoeps =   -2.515360   rhomu =   -3.272543  avg vxc =   -0.825511 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000560

 Energy terms:             smooth           local           total
   rhoval*vef             -5.478219      -184.922914      -190.401134
   rhoval*ves            -45.930034      -117.335734      -163.265768
   psnuc*ves              64.987712    -12970.136584    -12905.148872
   utot                    9.528839     -6543.736159     -6534.207320
   rho*exc                -2.515360      -127.430480      -129.945839
   rho*vxc                -3.272543      -168.734240      -172.006783
   valence chg             3.555443         7.444557        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=0.066107

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7153 -0.2526 -0.2483 -0.2483 -0.1888 -0.1888  1.7210  1.8717  1.8717
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4121 -0.2880 -0.2597 -0.2221 -0.1560 -0.0975  0.5961  1.2031  1.3985
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4748 -0.2770 -0.2491 -0.2327 -0.1755 -0.1419  0.7052  1.3231  1.7016
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4530 -0.3099 -0.2467 -0.1973 -0.1721 -0.1373  0.9157  1.1118  1.3611
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3389 -0.2998 -0.2860 -0.1980 -0.1618  0.0486  0.5511  0.7856  1.1108
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4977 -0.2604 -0.2428 -0.2428 -0.1689 -0.1689  0.6844  1.6526  1.6526

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.014427;  11.000000 electrons
         Sum occ. bands:   -2.765423, incl. Bloechl correction: -0.009080

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.107186    2.683244    7.423942

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.473421   -0.417406    4.677994    4.664139    4.500000    4.664139
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.405007   -0.445721    4.349210    4.334398    4.250000    4.334398
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.839017   -0.241044    3.863156    3.869539    3.147584    3.869539
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.018216   -0.727732    4.110000    4.108439    4.102416    4.110000

 Harris energy:
 sumev=       -2.765423  val*vef=    -190.401134   sumtv=     187.635711
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.945839     utot=   -6534.207320    ehar=   -3304.760809

 srhov:     -5.487236   -184.462988   -189.950223 sumev=   -2.765423   sumtv=  187.184801

 Kohn-Sham energy:
 sumtv=      187.184801  sumtc=      3171.756639   ekin=     3358.941440
 rhoep=     -129.896286   utot=     -6533.805664   ehks=    -3304.760509
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=1.06e-3  last it=2.66e-2
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.555443      3.576057      3.576057      0.000274      3.576057
 site    1    7.444557      7.423942      7.423942      0.000766      7.423942
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: nmix=2 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=1.06e-3
   tj:-0.04404  -0.00157
 unscreened rms difference:  smooth  0.000346   local  0.000766   tot
   screened rms difference:  smooth  0.000274   local  0.000766   tot  0.001055

 iors  : write restart file (binary, mesh density) 

   it  5  of 12    ehf=      -0.326309   ehk=      -0.326009
 From last iter    ehf=      -0.352894   ehk=      -0.235778
 diffe(q)=  0.026584 (0.001055)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf=-.3263093 ehk=-.3260093

 --- BNDFP:  begin iteration 6 of 12 ---

 avg es pot at rmt= 0.571782  avg sphere pot= 0.662667  vconst=-0.571782

 smooth rhoves      9.571474   charge     3.561414
 smooth rhoeps =   -2.519240   rhomu =   -3.277586  avg vxc =   -0.826125 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000558

 Energy terms:             smooth           local           total
   rhoval*vef             -5.484081      -184.917589      -190.401670
   rhoval*ves            -45.999373      -117.283030      -163.282403
   psnuc*ves              65.142320    -12970.133252    -12904.990932
   utot                    9.571474     -6543.708141     -6534.136667
   rho*exc                -2.519240      -127.411874      -129.931113
   rho*vxc                -3.277586      -168.709707      -171.987293
   valence chg             3.561414         7.438586        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.014427

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7173 -0.2613 -0.2570 -0.2570 -0.1979 -0.1979  1.7177  1.8684  1.8684
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4163 -0.2953 -0.2680 -0.2309 -0.1652 -0.1053  0.5922  1.1999  1.3950
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4778 -0.2844 -0.2576 -0.2414 -0.1846 -0.1508  0.7015  1.3199  1.6981
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4560 -0.3179 -0.2549 -0.2062 -0.1810 -0.1460  0.9123  1.1087  1.3576
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3447 -0.3074 -0.2939 -0.2069 -0.1711  0.0428  0.5472  0.7822  1.1074
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5006 -0.2674 -0.2515 -0.2515 -0.1782 -0.1782  0.6806  1.6493  1.6493

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.020483;  11.000000 electrons
         Sum occ. bands:   -2.851265, incl. Bloechl correction: -0.009203

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.116131    2.663451    7.452680

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.472518   -0.420342    4.664139    4.664126    4.500000    4.664126
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.397597   -0.454541    4.334398    4.332643    4.250000    4.332643
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.848578   -0.249221    3.869539    3.870439    3.147584    3.870439
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.017948   -0.733445    4.110000    4.108346    4.102416    4.110000

 Harris energy:
 sumev=       -2.851265  val*vef=    -190.401670   sumtv=     187.550405
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.931113     utot=   -6534.136667    ehar=   -3304.760736

 srhov:     -5.458915   -185.245576   -190.704491 sumev=   -2.851265   sumtv=  187.853226

 Kohn-Sham energy:
 sumtv=      187.853226  sumtc=      3171.756639   ekin=     3359.609865
 rhoep=     -129.958969   utot=     -6534.411548   ehks=    -3304.760652
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=7.40e-4  last it=1.06e-3
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.561414      3.547320      3.547320      0.000259      3.547320
 site    1    7.438586      7.452680      7.452680      0.000514      7.452680
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: nmix=2 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=7.4e-4
   tj: 0.32880   0.00545
 unscreened rms difference:  smooth  0.000200   local  0.000514   tot
   screened rms difference:  smooth  0.000259   local  0.000514   tot  0.000740

 iors  : write restart file (binary, mesh density) 

   it  6  of 12    ehf=      -0.326236   ehk=      -0.326152
 From last iter    ehf=      -0.326309   ehk=      -0.326009
 diffe(q)=  0.000073 (0.000740)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf=-.3262364 ehk=-.3261525

 --- BNDFP:  begin iteration 7 of 12 ---

 avg es pot at rmt= 0.571479  avg sphere pot= 0.662692  vconst=-0.571479

 smooth rhoves      9.560599   charge     3.558853
 smooth rhoeps =   -2.516916   rhomu =   -3.274555  avg vxc =   -0.825983 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000558

 Energy terms:             smooth           local           total
   rhoval*vef             -5.475845      -184.951259      -190.427105
   rhoval*ves            -45.984798      -117.318680      -163.303478
   psnuc*ves              65.105996    -12970.154482    -12905.048486
   utot                    9.560599     -6543.736581     -6534.175982
   rho*exc                -2.516916      -127.418278      -129.935195
   rho*vxc                -3.274555      -168.718150      -171.992705
   valence chg             3.558853         7.441147        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.020483

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7170 -0.2595 -0.2552 -0.2552 -0.1960 -0.1960  1.7182  1.8689  1.8689
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4155 -0.2938 -0.2663 -0.2290 -0.1633 -0.1037  0.5930  1.2005  1.3956
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4773 -0.2829 -0.2558 -0.2396 -0.1827 -0.1489  0.7022  1.3205  1.6987
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4555 -0.3162 -0.2532 -0.2043 -0.1791 -0.1442  0.9130  1.1092  1.3582
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3436 -0.3058 -0.2923 -0.2050 -0.1691  0.0440  0.5479  0.7828  1.1080
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5001 -0.2660 -0.2497 -0.2497 -0.1762 -0.1762  0.6814  1.6498  1.6498

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.019243;  11.000000 electrons
         Sum occ. bands:   -2.833300, incl. Bloechl correction: -0.009175

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.114511    2.667762    7.446749

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.472962   -0.419362    4.664126    4.664349    4.500000    4.664349
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.398510   -0.453484    4.332643    4.332816    4.250000    4.332816
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.847705   -0.247414    3.870439    3.870340    3.147584    3.870340
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.017993   -0.732733    4.110000    4.108357    4.102416    4.110000

 Harris energy:
 sumev=       -2.833300  val*vef=    -190.427105   sumtv=     187.593804
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.935195     utot=   -6534.175982    ehar=   -3304.760733

 srhov:     -5.465713   -185.069263   -190.534976 sumev=   -2.833300   sumtv=  187.701676

 Kohn-Sham energy:
 sumtv=      187.701676  sumtc=      3171.756639   ekin=     3359.458315
 rhoep=     -129.945456   utot=     -6534.273580   ehks=    -3304.760721
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=2.75e-4  last it=7.40e-4
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.558853      3.553251      3.553251      0.000103      3.553251
 site    1    7.441147      7.446749      7.446749      0.000198      7.446749
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: nmix=2 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=2.75e-4
   tj:-0.27367   0.11136
 unscreened rms difference:  smooth  0.000081   local  0.000198   tot
   screened rms difference:  smooth  0.000103   local  0.000198   tot  0.000275

 iors  : write restart file (binary, mesh density) 

   it  7  of 12    ehf=      -0.326233   ehk=      -0.326221
 From last iter    ehf=      -0.326236   ehk=      -0.326152
 diffe(q)=  0.000004 (0.000275)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf=-.3262327 ehk=-.3262211

 --- BNDFP:  begin iteration 8 of 12 ---

 avg es pot at rmt= 0.571303  avg sphere pot= 0.662733  vconst=-0.571303

 smooth rhoves      9.554396   charge     3.557414
 smooth rhoeps =   -2.515610   rhomu =   -3.272851  avg vxc =   -0.825903 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000558

 Energy terms:             smooth           local           total
   rhoval*vef             -5.471294      -184.960161      -190.431455
   rhoval*ves            -45.976455      -117.329765      -163.306220
   psnuc*ves              65.085247    -12970.151434    -12905.066187
   utot                    9.554396     -6543.740599     -6534.186203
   rho*exc                -2.515610      -127.421078      -129.936688
   rho*vxc                -3.272851      -168.721833      -171.994684
   valence chg             3.557414         7.442586        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.019243

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7168 -0.2587 -0.2544 -0.2544 -0.1952 -0.1952  1.7185  1.8692  1.8692
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4151 -0.2931 -0.2655 -0.2282 -0.1625 -0.1030  0.5933  1.2007  1.3959
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4770 -0.2822 -0.2551 -0.2388 -0.1819 -0.1481  0.7025  1.3207  1.6990
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4552 -0.3156 -0.2525 -0.2036 -0.1783 -0.1434  0.9132  1.1095  1.3585
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3431 -0.3052 -0.2916 -0.2042 -0.1683  0.0444  0.5483  0.7831  1.1083
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4998 -0.2654 -0.2489 -0.2489 -0.1754 -0.1754  0.6817  1.6501  1.6501

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018730;  11.000000 electrons
         Sum occ. bands:   -2.825938, incl. Bloechl correction: -0.009164

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.113797    2.669497    7.444300

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.473097   -0.419036    4.664349    4.664396    4.500000    4.664396
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.399003   -0.452916    4.332816    4.332919    4.250000    4.332919
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.847155   -0.246692    3.870340    3.870282    3.147584    3.870282
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.018013   -0.732361    4.110000    4.108363    4.102416    4.110000

 Harris energy:
 sumev=       -2.825938  val*vef=    -190.431455   sumtv=     187.605517
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.936688     utot=   -6534.186203    ehar=   -3304.760735

 srhov:     -5.468090   -184.999129   -190.467219 sumev=   -2.825938   sumtv=  187.641281

 Kohn-Sham energy:
 sumtv=      187.641281  sumtc=      3171.756639   ekin=     3359.397920
 rhoep=     -129.939954   utot=     -6534.218700   ehks=    -3304.760734
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=8.82e-5  last it=2.75e-4
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.557414      3.555700      3.555700      0.000033      3.555700
 site    1    7.442586      7.444300      7.444300      0.000062      7.444300
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=8.82e-5
   tj:-0.46979
 unscreened rms difference:  smooth  0.000026   local  0.000062   tot
   screened rms difference:  smooth  0.000033   local  0.000062   tot  0.000088

 iors  : write restart file (binary, mesh density) 

   it  8  of 12    ehf=      -0.326235   ehk=      -0.326234
 From last iter    ehf=      -0.326233   ehk=      -0.326221
 diffe(q)= -0.000002 (0.000088)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf=-.3262352 ehk=-.3262339

 --- BNDFP:  begin iteration 9 of 12 ---

 avg es pot at rmt= 0.571240  avg sphere pot= 0.662744  vconst=-0.571240

 smooth rhoves      9.552091   charge     3.556850
 smooth rhoeps =   -2.515086   rhomu =   -3.272167  avg vxc =   -0.825874 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000558

 Energy terms:             smooth           local           total
   rhoval*vef             -5.469413      -184.966553      -190.435966
   rhoval*ves            -45.973434      -117.336533      -163.309966
   psnuc*ves              65.077616    -12970.153887    -12905.076271
   utot                    9.552091     -6543.745210     -6534.193119
   rho*exc                -2.515086      -127.422320      -129.937406
   rho*vxc                -3.272167      -168.723469      -171.995637
   valence chg             3.556850         7.443150        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.01873

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7168 -0.2584 -0.2541 -0.2541 -0.1948 -0.1948  1.7186  1.8693  1.8693
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4150 -0.2929 -0.2652 -0.2279 -0.1621 -0.1027  0.5934  1.2008  1.3960
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4769 -0.2820 -0.2548 -0.2385 -0.1815 -0.1478  0.7026  1.3208  1.6991
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4551 -0.3153 -0.2522 -0.2032 -0.1780 -0.1431  0.9134  1.1096  1.3586
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3429 -0.3049 -0.2913 -0.2039 -0.1680  0.0447  0.5484  0.7832  1.1084
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4997 -0.2651 -0.2486 -0.2486 -0.1751 -0.1751  0.6818  1.6502  1.6502

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018513;  11.000000 electrons
         Sum occ. bands:   -2.822819, incl. Bloechl correction: -0.009160

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.113484    2.670273    7.443211

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.473158   -0.418903    4.664396    4.664418    4.500000    4.664418
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.399217   -0.452680    4.332919    4.332963    4.250000    4.332963
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.846904   -0.246385    3.870282    3.870257    3.147584    3.870257
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.018022   -0.732210    4.110000    4.108365    4.102416    4.110000

 Harris energy:
 sumev=       -2.822819  val*vef=    -190.435966   sumtv=     187.613148
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.937406     utot=   -6534.193119    ehar=   -3304.760737

 srhov:     -5.469249   -184.966960   -190.436209 sumev=   -2.822819   sumtv=  187.613391

 Kohn-Sham energy:
 sumtv=      187.613391  sumtc=      3171.756639   ekin=     3359.370030
 rhoep=     -129.937449   utot=     -6534.193318   ehks=    -3304.760737
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=2.34e-6  last it=8.82e-5
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.556850      3.556789      3.556789      0.000006      3.556789
 site    1    7.443150      7.443211      7.443211      0.000002      7.443211
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=2.34e-6
   tj:-0.01976
 unscreened rms difference:  smooth  0.000006   local  0.000002   tot
   screened rms difference:  smooth  0.000006   local  0.000002   tot  0.000002

 iors  : write restart file (binary, mesh density) 

   it  9  of 12    ehf=      -0.326237   ehk=      -0.326237
 From last iter    ehf=      -0.326235   ehk=      -0.326234
 diffe(q)= -0.000002 (0.000002)    tol= 0.000010 (0.000010)   more=F
c nk=8 bigbas=0 ehf=-.3262374 ehk=-.3262372
 Exit 0 LMF 
 wkinfo:  used   406 K  workspace of 80000 K   in  74 K calls
rdcmd:  rm mixm.cu
rdcmd:  lmf  --no-iactiv cu -vnk=8 -vbigbas=t -vpwmode=0 -voveps=0d-7
 -----------------------  START LMF (80000K)  -----------------------
 ptenv() is called with EXT=cu
 ptenv() not supported, but continue.

 LMF:      alat = 6.798  nbas = 1  nspec = 1  vn 7.00(LMF 7.0)  verb 31,30
 pot:      XC:BH
 bz:       metal(3), tetra, invit 

                Plat                                  Qlat
   0.000000   0.500000   0.500000       -1.000000   1.000000   1.000000
   0.500000   0.000000   0.500000        1.000000  -1.000000   1.000000
   0.500000   0.500000   0.000000        1.000000   1.000000  -1.000000
  Cell vol= 78.538660

 LATTC:  as= 2.000   tol=1.00E-12   alat= 6.79800   awald= 0.467
         r1=  2.252   nkd= 201      q1=  6.871   nkg= 331

 SGROUP: 1 symmetry operations from 0 generators
 SYMLAT: Bravais system is cubic with 48 symmetry operations.
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-5
 GROUPG: the following are sufficient to generate the space group:
         i*r3(1,1,-1) r4x
         i*r3(1,1,-1) r4x
 MKSYM:  found 48 space group operations ... includes inversion
 
 BZMESH:  60 irreducible QP from 512 ( 8 8 8 )  shift= T T T
 TETIRR: sorting 3072 tetrahedra ...
 264 inequivalent ones found

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv  kmxv foca   rfoca
 A        2.311  0.925    4    4         4  0.578  1.156    15    1   0.925

 MSHSIZ: mesh has 11 x 11 x 11 divisions; length 0.437, 0.437, 0.437
         generated from gmax = 9.0 a.u. : 941 vectors of 1331 (70%)

 gvlist: cutoff radius   9.000 gives    941   recips of max   1331
 SGVSYM: 41 symmetry stars found for 941 reciprocal lattice vectors
 

 Makidx:  hamiltonian dimensions Low, Int, High, Negl: 31 0 28 16
 kappa   Low   Int   High  L+I  L+I+H  Neglected
   1       9     0    16     9    25       0
   2      13     0    12    13    25       0
   3       9     0     0     9     9      16
  all     31     0    28    31    59      16
 suham :  25 augmentation channels, 25 local potential channels  Maximum lmxa=4

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol= 1.00E-06
 spec      l    rsm    eh     gmax    last term    cutoff
  A        0*   2.50  -0.01   2.974    2.30E-05      27 
  A        1    2.50  -0.01   3.093    1.29E-06      51 
  A        2    1.00  -0.01   8.508    1.16E-06     869 
  A        0    1.30  -1.00   5.718    2.28E-06     259 
  A        2    1.00  -1.00   8.508    1.16E-06     869 
  A        3    1.30  -0.01   6.806    2.09E-06     411 

 iors  : read restart file (binary, mesh density) 
         use from  restart file: ef window, positions, pnu 
         ignore in restart file: *
         site   1, species A       : augmentation lmax changed from 3 to 4
         site   1, species A       : inflate local density from nlm= 16 to 25

 --- BNDFP:  begin iteration 1 of 12 ---

 avg es pot at rmt= 0.571237  avg sphere pot= 0.662746  vconst=-0.571237

 smooth rhoves      9.551947   charge     3.556810
 smooth rhoeps =   -2.515047   rhomu =   -3.272117  avg vxc =   -0.825872 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000558

 Energy terms:             smooth           local           total
   rhoval*vef             -5.469267      -184.966399      -190.435666
   rhoval*ves            -45.973260      -117.336408      -163.309668
   psnuc*ves              65.077154    -12970.153157    -12905.076004
   utot                    9.551947     -6543.744783     -6534.192836
   rho*exc                -2.515047      -127.422356      -129.937403
   rho*vxc                -3.272117      -168.723515      -171.995632
   valence chg             3.556810         7.443190        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7176 -0.2587 -0.2544 -0.2544 -0.1946 -0.1946  1.5859  1.7373  1.8052
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4154 -0.2932 -0.2656 -0.2280 -0.1622 -0.1029  0.5909  1.1772  1.3092
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4775 -0.2825 -0.2549 -0.2386 -0.1817 -0.1479  0.6977  1.2862  1.5948
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4557 -0.3158 -0.2525 -0.2034 -0.1781 -0.1433  0.9069  1.0961  1.2465
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3433 -0.3051 -0.2917 -0.2038 -0.1683  0.0443  0.5466  0.7740  1.0777
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5003 -0.2656 -0.2487 -0.2487 -0.1751 -0.1751  0.6763  1.6276  1.6276

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018792;  11.000000 electrons
         Sum occ. bands:   -2.825304, incl. Bloechl correction: -0.009153

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7176 -0.2587 -0.2544 -0.2544 -0.1946 -0.1946  1.5859  1.7373  1.8052
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4154 -0.2932 -0.2656 -0.2280 -0.1622 -0.1029  0.5909  1.1772  1.3092
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4775 -0.2825 -0.2549 -0.2386 -0.1817 -0.1479  0.6977  1.2862  1.5948
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4557 -0.3158 -0.2525 -0.2034 -0.1781 -0.1433  0.9069  1.0961  1.2465
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3433 -0.3051 -0.2917 -0.2038 -0.1683  0.0443  0.5466  0.7740  1.0777
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5003 -0.2656 -0.2487 -0.2487 -0.1751 -0.1751  0.6763  1.6276  1.6276

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018792;  11.000000 electrons
         Sum occ. bands:   -2.825304, incl. Bloechl correction: -0.009153

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.128160    2.833312    7.294848

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.480584   -0.533465    4.664418    4.612337    4.500000    4.612337
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.407336   -0.469584    4.332963    4.328212    4.250000    4.328212
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.971687   -0.255973    3.870257    3.858616    3.147584    3.858616
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021180   -0.806341    4.110000    4.106719    4.102416    4.110000
 4     1    0.004395   -0.781220    5.100000    5.079490    5.077979    5.100000

 Harris energy:
 sumev=       -2.825304  val*vef=    -190.435666   sumtv=     187.610362
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.937403     utot=   -6534.192836    ehar=   -3304.763238

 srhov:     -5.929501   -184.355314   -190.284815 sumev=   -2.825304   sumtv=  187.459511

 Kohn-Sham energy:
 sumtv=      187.459511  sumtc=      3171.756639   ekin=     3359.216151
 rhoep=     -129.927456   utot=     -6534.051861   ehks=    -3304.763166
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 0.  RMS DQ=4.37e-3
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.556810      3.705152      3.705152      0.003469      3.705152
 site    1    7.443190      7.294848      7.294848      0.004825      7.294848
 AMIX: nmix=0 mmix=8  nelts=2178  beta=1  tm=5  rmsdel=4.37e-3
 unscreened rms difference:  smooth  0.003486   local  0.004825   tot
   screened rms difference:  smooth  0.003469   local  0.004825   tot  0.004371

 iors  : write restart file (binary, mesh density) 

   it  1  of 12    ehf=      -0.328738   ehk=      -0.328666
i nk=8 bigbas=1 pwmode=0 oveps=0 ehf=-.3287376 ehk=-.3286662

 --- BNDFP:  begin iteration 2 of 12 ---

 avg es pot at rmt= 0.581158  avg sphere pot= 0.653091  vconst=-0.581158

 smooth rhoves     10.108795   charge     3.705152
 smooth rhoeps =   -2.657496   rhomu =   -3.458045  avg vxc =   -0.833024 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000127

 Energy terms:             smooth           local           total
   rhoval*vef             -6.009415      -184.304839      -190.314253
   rhoval*ves            -46.615710      -116.582740      -163.198450
   psnuc*ves              66.833301    -12971.740403    -12904.907102
   utot                   10.108795     -6544.161572     -6534.052776
   rho*exc                -2.657496      -127.270006      -129.927501
   rho*vxc                -3.458045      -168.524429      -171.982474
   valence chg             3.705152         7.294848        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7177 -0.2622 -0.2578 -0.2578 -0.1970 -0.1970  1.5859  1.7372  1.8045
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4165 -0.2957 -0.2683 -0.2308 -0.1654 -0.1056  0.5898  1.1765  1.3087
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4780 -0.2853 -0.2576 -0.2415 -0.1847 -0.1509  0.6967  1.2855  1.5944
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4562 -0.3188 -0.2550 -0.2066 -0.1810 -0.1463  0.9061  1.0955  1.2460
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3450 -0.3076 -0.2945 -0.2062 -0.1720  0.0424  0.5456  0.7733  1.0771
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5008 -0.2682 -0.2516 -0.2516 -0.1782 -0.1782  0.6752  1.6269  1.6269

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.020714;  11.000000 electrons
         Sum occ. bands:   -2.853789, incl. Bloechl correction: -0.009204

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7177 -0.2622 -0.2578 -0.2578 -0.1970 -0.1970  1.5859  1.7372  1.8045
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4165 -0.2957 -0.2683 -0.2308 -0.1654 -0.1056  0.5898  1.1765  1.3087
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4780 -0.2853 -0.2576 -0.2415 -0.1847 -0.1509  0.6967  1.2855  1.5944
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4562 -0.3188 -0.2550 -0.2066 -0.1810 -0.1463  0.9061  1.0955  1.2460
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3450 -0.3076 -0.2945 -0.2062 -0.1720  0.0424  0.5456  0.7733  1.0771
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5008 -0.2682 -0.2516 -0.2516 -0.1782 -0.1782  0.6752  1.6269  1.6269

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.020714;  11.000000 electrons
         Sum occ. bands:   -2.853789, incl. Bloechl correction: -0.009204

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.131732    2.825132    7.306600

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.479387   -0.535597    4.612337    4.611672    4.500000    4.611672
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.404856   -0.472668    4.328212    4.327646    4.250000    4.327646
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.969521   -0.258624    3.858616    3.858674    3.147584    3.858674
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021024   -0.810002    4.110000    4.106662    4.102416    4.110000
 4     1    0.004375   -0.784068    5.100000    5.079474    5.077979    5.100000

 Harris energy:
 sumev=       -2.853789  val*vef=    -190.314253   sumtv=     187.460464
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.927501     utot=   -6534.052776    ehar=   -3304.763174

 srhov:     -5.994048   -184.684580   -190.678628 sumev=   -2.853789   sumtv=  187.824839

 Kohn-Sham energy:
 sumtv=      187.824839  sumtc=      3171.756639   ekin=     3359.581478
 rhoep=     -129.959113   utot=     -6534.385433   ehks=    -3304.763067
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 1.  RMS DQ=7.67e-4  last it=4.37e-3
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.705152      3.693400      3.693400      0.000208      3.693400
 site    1    7.294848      7.306600      7.306600      0.000500      7.306600
 AMIX: nmix=1 mmix=8  nelts=2178  beta=1  tm=5  rmsdel=7.67e-4
   tj: 0.08056
 unscreened rms difference:  smooth  0.000174   local  0.000500   tot
   screened rms difference:  smooth  0.000208   local  0.000500   tot  0.000767

 iors  : write restart file (binary, mesh density) 

   it  2  of 12    ehf=      -0.328674   ehk=      -0.328567
 From last iter    ehf=      -0.328738   ehk=      -0.328666
 diffe(q)=  0.000064 (0.000767)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=1 pwmode=0 oveps=0 ehf=-.328674 ehk=-.3285671

 --- BNDFP:  begin iteration 3 of 12 ---

 avg es pot at rmt= 0.579980  avg sphere pot= 0.653052  vconst=-0.579980

 smooth rhoves     10.061393   charge     3.694347
 smooth rhoeps =   -2.647850   rhomu =   -3.445462  avg vxc =   -0.832366 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000134

 Energy terms:             smooth           local           total
   rhoval*vef             -5.973948      -184.535193      -190.509140
   rhoval*ves            -46.555926      -116.812277      -163.368203
   psnuc*ves              66.678712    -12971.953009    -12905.274297
   utot                   10.061393     -6544.382643     -6534.321250
   rho*exc                -2.647850      -127.303260      -129.951109
   rho*vxc                -3.445462      -168.568326      -172.013789
   valence chg             3.694347         7.305653        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7163 -0.2523 -0.2479 -0.2479 -0.1865 -0.1865  1.5877  1.7387  1.8068
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4124 -0.2876 -0.2588 -0.2208 -0.1547 -0.0967  0.5937  1.1793  1.3115
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4754 -0.2769 -0.2478 -0.2315 -0.1743 -0.1406  0.7004  1.2883  1.5968
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4535 -0.3098 -0.2457 -0.1963 -0.1708 -0.1362  0.9093  1.0982  1.2490
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3388 -0.2991 -0.2856 -0.1960 -0.1612  0.0488  0.5494  0.7764  1.0801
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4983 -0.2604 -0.2416 -0.2416 -0.1675 -0.1675  0.6790  1.6297  1.6297

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.014006;  11.000000 electrons
         Sum occ. bands:   -2.756738, incl. Bloechl correction: -0.009049

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7163 -0.2523 -0.2479 -0.2479 -0.1865 -0.1865  1.5877  1.7387  1.8068
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4124 -0.2876 -0.2588 -0.2208 -0.1547 -0.0967  0.5937  1.1793  1.3115
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4754 -0.2769 -0.2478 -0.2315 -0.1743 -0.1406  0.7004  1.2883  1.5968
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4535 -0.3098 -0.2457 -0.1963 -0.1708 -0.1362  0.9093  1.0982  1.2490
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3388 -0.2991 -0.2856 -0.1960 -0.1612  0.0488  0.5494  0.7764  1.0801
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4983 -0.2604 -0.2416 -0.2416 -0.1675 -0.1675  0.6790  1.6297  1.6297

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.014006;  11.000000 electrons
         Sum occ. bands:   -2.756738, incl. Bloechl correction: -0.009049

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.121437    2.845312    7.276124

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.481245   -0.529907    4.611672    4.613320    4.500000    4.613320
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.411022   -0.464864    4.327646    4.329215    4.250000    4.329215
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.950064   -0.248938    3.858674    3.858073    3.147584    3.858073
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021443   -0.799134    4.110000    4.106882    4.102416    4.110000
 4     1    0.004440   -0.772503    5.100000    5.079581    5.077979    5.100000

 Harris energy:
 sumev=       -2.756738  val*vef=    -190.509140   sumtv=     187.752402
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.951109     utot=   -6534.321250    ehar=   -3304.763318

 srhov:     -6.021913   -183.666474   -189.688387 sumev=   -2.756738   sumtv=  186.931649

 Kohn-Sham energy:
 sumtv=      186.931649  sumtc=      3171.756639   ekin=     3358.688288
 rhoep=     -129.879167   utot=     -6533.571887   ehks=    -3304.762766
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 2.  RMS DQ=1.80e-3  last it=7.67e-4
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.694347      3.723876      3.723876      0.000529      3.723876
 site    1    7.305653      7.276124      7.276124      0.001211      7.276124
 AMIX: nmix=2 mmix=8  nelts=2178  beta=1  tm=5  rmsdel=1.8e-3
   tj: 0.70189  -0.00506
 unscreened rms difference:  smooth  0.000390   local  0.001211   tot
   screened rms difference:  smooth  0.000529   local  0.001211   tot  0.001801

 iors  : write restart file (binary, mesh density) 

   it  3  of 12    ehf=      -0.328818   ehk=      -0.328266
 From last iter    ehf=      -0.328674   ehk=      -0.328567
 diffe(q)= -0.000143 (0.001801)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=1 pwmode=0 oveps=0 ehf=-.3288175 ehk=-.3282662

 --- BNDFP:  begin iteration 4 of 12 ---

 avg es pot at rmt= 0.580860  avg sphere pot= 0.653051  vconst=-0.580860

 smooth rhoves     10.096266   charge     3.702580
 smooth rhoeps =   -2.655417   rhomu =   -3.455336  avg vxc =   -0.832820 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000135

 Energy terms:             smooth           local           total
   rhoval*vef             -6.002045      -184.376183      -190.378228
   rhoval*ves            -46.599143      -116.655488      -163.254631
   psnuc*ves              66.791674    -12971.816331    -12905.024657
   utot                   10.096266     -6544.235910     -6534.139644
   rho*exc                -2.655417      -127.279393      -129.934810
   rho*vxc                -3.455336      -168.536831      -171.992167
   valence chg             3.702580         7.297420        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7173 -0.2591 -0.2547 -0.2547 -0.1937 -0.1937  1.5865  1.7376  1.8052
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4152 -0.2932 -0.2653 -0.2277 -0.1621 -0.1029  0.5910  1.1774  1.3096
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4772 -0.2826 -0.2546 -0.2384 -0.1815 -0.1477  0.6978  1.2864  1.5951
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4553 -0.3160 -0.2521 -0.2034 -0.1778 -0.1431  0.9070  1.0963  1.2470
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3431 -0.3049 -0.2917 -0.2031 -0.1686  0.0444  0.5468  0.7742  1.0780
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5000 -0.2658 -0.2485 -0.2485 -0.1749 -0.1749  0.6764  1.6277  1.6277

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018651;  11.000000 electrons
         Sum occ. bands:   -2.823576, incl. Bloechl correction: -0.009156

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7173 -0.2591 -0.2547 -0.2547 -0.1937 -0.1937  1.5865  1.7376  1.8052
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4152 -0.2932 -0.2653 -0.2277 -0.1621 -0.1029  0.5910  1.1774  1.3096
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4772 -0.2826 -0.2546 -0.2384 -0.1815 -0.1477  0.6978  1.2864  1.5951
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4553 -0.3160 -0.2521 -0.2034 -0.1778 -0.1431  0.9070  1.0963  1.2470
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3431 -0.3049 -0.2917 -0.2031 -0.1686  0.0444  0.5468  0.7742  1.0780
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5000 -0.2658 -0.2485 -0.2485 -0.1749 -0.1749  0.6764  1.6277  1.6277

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018651;  11.000000 electrons
         Sum occ. bands:   -2.823576, incl. Bloechl correction: -0.009156

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.128443    2.831457    7.296986

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.479978   -0.533700    4.613320    4.612251    4.500000    4.612251
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.406875   -0.470206    4.329215    4.328151    4.250000    4.328151
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.963081   -0.255629    3.858073    3.858466    3.147584    3.858466
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021154   -0.806683    4.110000    4.106730    4.102416    4.110000
 4     1    0.004395   -0.780578    5.100000    5.079507    5.077979    5.100000

 Harris energy:
 sumev=       -2.823576  val*vef=    -190.378228   sumtv=     187.554652
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.934810     utot=   -6534.139644    ehar=   -3304.763162

 srhov:     -6.003000   -184.365805   -190.368804 sumev=   -2.823576   sumtv=  187.545228

 Kohn-Sham energy:
 sumtv=      187.545228  sumtc=      3171.756639   ekin=     3359.301868
 rhoep=     -129.934006   utot=     -6534.131024   ehks=    -3304.763162
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 3.  RMS DQ=2.21e-5  last it=1.80e-3
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.702580      3.703014      3.703014      0.000017      3.703014
 site    1    7.297420      7.296986      7.296986      0.000016      7.296986
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=2178  beta=1  tm=5  rmsdel=2.21e-5
   tj:-0.01228
 unscreened rms difference:  smooth  0.000021   local  0.000016   tot
   screened rms difference:  smooth  0.000017   local  0.000016   tot  0.000022

 iors  : write restart file (binary, mesh density) 

   it  4  of 12    ehf=      -0.328662   ehk=      -0.328662
 From last iter    ehf=      -0.328818   ehk=      -0.328266
 diffe(q)=  0.000155 (0.000022)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=1 pwmode=0 oveps=0 ehf=-.3286622 ehk=-.328662

 --- BNDFP:  begin iteration 5 of 12 ---

 avg es pot at rmt= 0.580888  avg sphere pot= 0.653045  vconst=-0.580888

 smooth rhoves     10.096995   charge     3.702758
 smooth rhoeps =   -2.655612   rhomu =   -3.455591  avg vxc =   -0.832823 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000135

 Energy terms:             smooth           local           total
   rhoval*vef             -6.002699      -184.374100      -190.376800
   rhoval*ves            -46.599988      -116.653393      -163.253381
   psnuc*ves              66.793979    -12971.815985    -12905.022006
   utot                   10.096995     -6544.234689     -6534.137694
   rho*exc                -2.655612      -127.279025      -129.934637
   rho*vxc                -3.455591      -168.536346      -171.991937
   valence chg             3.702758         7.297242        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7173 -0.2591 -0.2548 -0.2548 -0.1938 -0.1938  1.5865  1.7376  1.8051
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4153 -0.2932 -0.2654 -0.2278 -0.1621 -0.1029  0.5910  1.1774  1.3095
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4772 -0.2827 -0.2546 -0.2385 -0.1815 -0.1478  0.6978  1.2863  1.5951
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4554 -0.3160 -0.2522 -0.2035 -0.1779 -0.1432  0.9070  1.0963  1.2469
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3431 -0.3050 -0.2918 -0.2031 -0.1687  0.0443  0.5467  0.7742  1.0780
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5000 -0.2659 -0.2486 -0.2486 -0.1749 -0.1749  0.6763  1.6277  1.6277

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018702;  11.000000 electrons
         Sum occ. bands:   -2.824270, incl. Bloechl correction: -0.009157

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7173 -0.2591 -0.2548 -0.2548 -0.1938 -0.1938  1.5865  1.7376  1.8051
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4153 -0.2932 -0.2654 -0.2278 -0.1621 -0.1029  0.5910  1.1774  1.3095
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4772 -0.2827 -0.2546 -0.2385 -0.1815 -0.1478  0.6978  1.2863  1.5951
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4554 -0.3160 -0.2522 -0.2035 -0.1779 -0.1432  0.9070  1.0963  1.2469
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3431 -0.3050 -0.2918 -0.2031 -0.1687  0.0443  0.5467  0.7742  1.0780
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5000 -0.2659 -0.2486 -0.2486 -0.1749 -0.1749  0.6763  1.6277  1.6277

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018702;  11.000000 electrons
         Sum occ. bands:   -2.824270, incl. Bloechl correction: -0.009157

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.128504    2.831326    7.297178

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.479955   -0.533756    4.612251    4.612234    4.500000    4.612234
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.406826   -0.470263    4.328151    4.328141    4.250000    4.328141
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.963382   -0.255690    3.858466    3.858479    3.147584    3.858479
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021151   -0.806765    4.110000    4.106729    4.102416    4.110000
 4     1    0.004395   -0.780659    5.100000    5.079506    5.077979    5.100000

 Harris energy:
 sumev=       -2.824270  val*vef=    -190.376800   sumtv=     187.552529
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.934637     utot=   -6534.137694    ehar=   -3304.763162

 srhov:     -6.002879   -184.372691   -190.375570 sumev=   -2.824270   sumtv=  187.551300

 Kohn-Sham energy:
 sumtv=      187.551300  sumtc=      3171.756639   ekin=     3359.307939
 rhoep=     -129.934540   utot=     -6534.136561   ehks=    -3304.763162
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=3.06e-6  last it=2.21e-5
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.702758      3.702822      3.702822      0.000008      3.702822
 site    1    7.297242      7.297178      7.297178      0.000002      7.297178
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=2178  beta=1  tm=5  rmsdel=3.06e-6
   tj:-0.15999
 unscreened rms difference:  smooth  0.000009   local  0.000002   tot
   screened rms difference:  smooth  0.000008   local  0.000002   tot  0.000003

 iors  : write restart file (binary, mesh density) 

   it  5  of 12    ehf=      -0.328662   ehk=      -0.328662
 From last iter    ehf=      -0.328662   ehk=      -0.328662
 diffe(q)=  0.000000 (0.000003)    tol= 0.000010 (0.000010)   more=F
c nk=8 bigbas=1 pwmode=0 oveps=0 ehf=-.3286618 ehk=-.3286617
 Exit 0 LMF 
 wkinfo:  used   626 K  workspace of 80000 K   in 242 K calls
rdcmd:  lmf  --no-iactiv cu -vnk=8 -vbigbas=t -vpwmode=0 -voveps=0d-7 --band:fn=syml
 -----------------------  START LMF (80000K)  -----------------------
 ptenv() is called with EXT=cu
 ptenv() not supported, but continue.

 LMF:      alat = 6.798  nbas = 1  nspec = 1  vn 7.00(LMF 7.0)  verb 31,30
 pot:      XC:BH
 bz:       metal(3), tetra, invit 

                Plat                                  Qlat
   0.000000   0.500000   0.500000       -1.000000   1.000000   1.000000
   0.500000   0.000000   0.500000        1.000000  -1.000000   1.000000
   0.500000   0.500000   0.000000        1.000000   1.000000  -1.000000
  Cell vol= 78.538660

 LATTC:  as= 2.000   tol=1.00E-12   alat= 6.79800   awald= 0.467
         r1=  2.252   nkd= 201      q1=  6.871   nkg= 331

 SGROUP: 1 symmetry operations from 0 generators
 SYMLAT: Bravais system is cubic with 48 symmetry operations.
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-5
 GROUPG: the following are sufficient to generate the space group:
         i*r3(1,1,-1) r4x
         i*r3(1,1,-1) r4x
 MKSYM:  found 48 space group operations ... includes inversion
 
 BZMESH:  60 irreducible QP from 512 ( 8 8 8 )  shift= T T T
 TETIRR: sorting 3072 tetrahedra ...
 264 inequivalent ones found

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv  kmxv foca   rfoca
 A        2.311  0.925    4    4         4  0.578  1.156    15    1   0.925

 MSHSIZ: mesh has 11 x 11 x 11 divisions; length 0.437, 0.437, 0.437
         generated from gmax = 9.0 a.u. : 941 vectors of 1331 (70%)

 gvlist: cutoff radius   9.000 gives    941   recips of max   1331
 SGVSYM: 41 symmetry stars found for 941 reciprocal lattice vectors
 

 Makidx:  hamiltonian dimensions Low, Int, High, Negl: 31 0 28 16
 kappa   Low   Int   High  L+I  L+I+H  Neglected
   1       9     0    16     9    25       0
   2      13     0    12    13    25       0
   3       9     0     0     9     9      16
  all     31     0    28    31    59      16
 suham :  25 augmentation channels, 25 local potential channels  Maximum lmxa=4

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol= 1.00E-06
 spec      l    rsm    eh     gmax    last term    cutoff
  A        0*   2.50  -0.01   2.974    2.30E-05      27 
  A        1    2.50  -0.01   3.093    1.29E-06      51 
  A        2    1.00  -0.01   8.508    1.16E-06     869 
  A        0    1.30  -1.00   5.718    2.28E-06     259 
  A        2    1.00  -1.00   8.508    1.16E-06     869 
  A        3    1.30  -0.01   6.806    2.09E-06     411 

 iors  : read restart file (binary, mesh density) 
         use from  restart file: ef window, positions, pnu 
         ignore in restart file: *

 --- BNDFP:  begin iteration 1 of 12 ---

 avg es pot at rmt= 0.580894  avg sphere pot= 0.653044  vconst=-0.580894

 smooth rhoves     10.097116   charge     3.702792
 smooth rhoeps =   -2.655652   rhomu =   -3.455643  avg vxc =   -0.832822 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000135

 Energy terms:             smooth           local           total
   rhoval*vef             -6.002841      -184.373715      -190.376556
   rhoval*ves            -46.600123      -116.653034      -163.253156
   psnuc*ves              66.794355    -12971.816011    -12905.021656
   utot                   10.097116     -6544.234522     -6534.137406
   rho*exc                -2.655652      -127.278965      -129.934618
   rho*vxc                -3.455643      -168.536268      -171.991911
   valence chg             3.702792         7.297208        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read efermi from weights file : ef = -0.018702
 suqlst:  generate bands, mode 1

 suqlst:  nq= 41   q1= 0.5000 0.5000 0.5000   q2= 0.0000 0.0000 0.0000
 bndfp:  kpt 1 of 41, k=  -0.50000  -0.50000  -0.50000
 -0.4126 -0.2588 -0.2588 -0.1500 -0.1500 -0.1050  0.2487  1.5877  1.5973
 bndfp:  kpt 11 of 41, k=  0.37500  0.37500  0.37500
 -0.4449 -0.2536 -0.2536 -0.2248 -0.1624 -0.1624  0.4812  1.6094  1.6094
 bndfp:  kpt 21 of 41, k=  0.25000  0.25000  0.25000
 -0.5703 -0.2760 -0.2447 -0.2447 -0.1870 -0.1870  0.8928  1.6566  1.6566
 bndfp:  kpt 31 of 41, k=  0.12500  0.12500  0.12500
 -0.6860 -0.2652 -0.2497 -0.2497 -0.1953 -0.1953  1.3655  1.7127  1.7482
 bndfp:  kpt 41 of 41, k=  0.00000  0.00000  0.00000
 -0.7279 -0.2568 -0.2568 -0.2568 -0.1930 -0.1930  1.6567  1.8387  1.8387

 Read efermi from weights file : ef = -0.018702

 suqlst:  nq= 41   q1= 0.0000 0.0000 0.0000   q2= 1.0000 0.0000 0.0000
 bndfp:  kpt 1 of 41, k=  0.00000  0.00000  0.00000
 -0.7279 -0.2568 -0.2568 -0.2568 -0.1930 -0.1930  1.6567  1.8387  1.8387
 bndfp:  kpt 11 of 41, k=  0.25000  0.00000  0.00000
 -0.6722 -0.2732 -0.2433 -0.2433 -0.2075 -0.1870  1.5242  1.5425  1.5425
 bndfp:  kpt 21 of 41, k=  0.50000  0.00000  0.00000
 -0.5227 -0.3114 -0.2283 -0.2064 -0.2064 -0.1723  1.2328  1.2328  1.3539
 bndfp:  kpt 31 of 41, k=  0.75000  0.00000  0.00000
 -0.4048 -0.3480 -0.1614 -0.1614 -0.1572 -0.1230  0.8184  1.0147  1.0147
 bndfp:  kpt 41 of 41, k=  -1.00000  0.00000  0.00000
 -0.3961 -0.3627 -0.1509 -0.1394 -0.1394  0.0772  0.4995  0.9336  0.9336

 Read efermi from weights file : ef = -0.018702

 suqlst:  nq= 21   q1= 1.0000 0.0000 0.0000   q2= 1.0000 0.5000 0.0000
 bndfp:  kpt 1 of 21, k=  -1.00000  0.00000  0.00000
 -0.3961 -0.3627 -0.1509 -0.1394 -0.1394  0.0772  0.4995  0.9336  0.9336
 bndfp:  kpt 11 of 21, k=  -1.00000  0.25000  0.00000
 -0.3718 -0.3463 -0.2102 -0.1698 -0.1394  0.1977  0.5208  0.7107  0.7756
 bndfp:  kpt 21 of 21, k=  -1.00000  0.50000  0.00000
 -0.3401 -0.2953 -0.2953 -0.1964 -0.1394  0.4291  0.4291  0.5910  0.6156

 Read efermi from weights file : ef = -0.018702

 suqlst:  nq= 41   q1= 1.0000 0.5000 0.0000   q2= 0.0000 0.0000 0.0000
 bndfp:  kpt 1 of 41, k=  -1.00000  0.50000  0.00000
 -0.3401 -0.2953 -0.2953 -0.1964 -0.1394  0.4291  0.4291  0.5910  0.6156
 bndfp:  kpt 11 of 41, k=  0.75000  0.37500  0.00000
 -0.3461 -0.3172 -0.2822 -0.1886 -0.1617  0.0625  0.6342  0.7293  0.9129
 bndfp:  kpt 21 of 41, k=  0.50000  0.25000  0.00000
 -0.4790 -0.3053 -0.2569 -0.2057 -0.2018 -0.1412  0.9570  1.0059  1.4114
 bndfp:  kpt 31 of 41, k=  0.25000  0.12500  0.00000
 -0.6587 -0.2738 -0.2436 -0.2428 -0.2061 -0.1870  1.3635  1.4235  1.7178
 bndfp:  kpt 41 of 41, k=  0.00000  0.00000  0.00000
 -0.7279 -0.2568 -0.2568 -0.2568 -0.1930 -0.1930  1.6567  1.8387  1.8387

 Read efermi from weights file : ef = -0.018702
 Exit 0 bndfp 
 wkinfo:  used   555 K  workspace of 80000 K   in  53 K calls
