 -----------------------  START LMF (80000K)  -----------------------
 HEADER Cu

 rdctrl: reset global max nl from 5 to 6

 LMF:      alat = 6.798  nbas = 1  nspec = 1  vn 7.00(LMF 7.0)  verb 31,20
 special:  APW basis
 pot:      XC:BH
 bz:       metal(3), tetra, invit 

                Plat                                  Qlat
   0.000000   0.500000   0.500000       -1.000000   1.000000   1.000000
   0.500000   0.000000   0.500000        1.000000  -1.000000   1.000000
   0.500000   0.500000   0.000000        1.000000   1.000000  -1.000000
  Cell vol= 78.538660

 LATTC:  as= 2.000   tol=1.00E-08   alat= 6.79800   awald= 0.467
         r1=  1.959   nkd= 135      q1=  5.910   nkg= 181

 SGROUP: 1 symmetry operations from 0 generators
 SYMLAT: Bravais system is cubic with 48 symmetry operations.
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-5
 GROUPG: the following are sufficient to generate the space group:
         i*r3(1,1,-1) r4x
         i*r3(1,1,-1) r4x
 MKSYM:  found 48 space group operations ... includes inversion
 
 BZMESH:  29 irreducible QP from 512 ( 8 8 8 )  shift= F F F
 TETIRR: sorting 3072 tetrahedra ...
 76 inequivalent ones found

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv  kmxv foca   rfoca
 Cu       2.280  0.912    5    5         5  0.570  1.140    15    1   0.912

 gvlist: cutoff radius   9.803 gives   1243   recips of max   3375
 SGVSYM: 53 symmetry stars found for 1243 reciprocal lattice vectors
 

 Makidx:  hamiltonian dimensions Low, Int, High, Negl: 9 0 27 0
 suham :  36 augmentation channels, 36 local potential channels  Maximum lmxa=5

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol= 1.00E-06
 spec      l    rsm    eh     gmax    last term    cutoff
  Cu       0    1.50  -0.28   4.956    2.32E-06     169 
  Cu       1    1.50  -0.10   5.245    1.10E-06     181 
  Cu       2    0.95  -0.11   8.973    1.81E-06     941 

 suham:  q-dependent PW basis with  Emin = 0 < E < 5.
         Est. min,max PW dimension = 12,18.  Use npwpad = 3 => ndham = 30

 iors  : read restart file (binary, mesh density) 
 iors  : empty file ... nothing read

 rdovfa: read and overlap free-atom densities (mesh density) ...
 rdovfa: expected Cu,      read Cu       with rmt=  2.2800  mesh   393  0.025

 ovlpfa: overlap smooth part of FA densities
 site   1  spec  1  pos  0.0000  0.0000  0.0000  Qsmooth 4.729619
 total smooth Q = 4.729619

 Free atom and overlapped crystal site charges:
   ib    true(FA)    smooth(FA)  true(OV)    smooth(OV)    local
    1    9.758197    3.487816   10.209423    3.939042    6.270381

 Smooth charge on mesh:            4.729619
 Sum of local charges:             6.270381
 Total valence charge:            11.000000
 Sum of core charges:             18.000000
 Sum of nuclear charges:         -29.000000
 Homogeneous background:           0.000000
 Deviation from neutrality:        0.000000

 --- BNDFP:  begin iteration 1 of 20 ---

 avg es pot at rmt= 0.554330  avg sphere pot= 0.625056  vconst=-0.554330

 site  1  z= 29.0  rmt= 2.28000  nr=393   a=0.025  nlml=36  rg=0.570  Vfloat=T
 sm core charge = 0.295519 (sphere) + 0.00535 (spillout) = 0.300868
 potential shift to crystal energy zero:    0.000086


 subzi: tetrahedron integration of bands; tetrahedron integration of density

 Start first of two band passes ...
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.95D-06    2  0.57D-03    3  0.57D-03    4  0.57D-03    5  0.40D-01
    6  0.40D-01    7  0.51D-01    8  0.90D-01    9  0.90D-01   10  0.90D-01
   11  0.40D+00   12  0.40D+00   13  0.40D+00   14  0.47D+00   15  0.59D+00
   16  0.10D+05   17  0.10D+05   18  0.10D+05   19  0.10D+05   20  0.10D+05
   21  0.10D+05   22  0.10D+05   23  0.10D+05   24  0.10D+05
 eigenvalue=
 -0.6650 -0.0472 -0.0472 -0.0472  0.0270  0.0270  1.7100  1.9162  1.9162
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.4185 -0.0983 -0.0478 -0.0233  0.0207  0.0764  0.8953  1.2978  1.7913
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.2385 -0.1775  0.0760  0.0911  0.0911  0.1426  0.6393  1.0060  1.0060

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.151736;  11.000000 electrons
         Sum occ. bands:   -0.756621, incl. Bloechl correction: -0.006665

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.95D-06    2  0.57D-03    3  0.57D-03    4  0.57D-03    5  0.40D-01
    6  0.40D-01    7  0.51D-01    8  0.90D-01    9  0.90D-01   10  0.90D-01
   11  0.40D+00   12  0.40D+00   13  0.40D+00   14  0.47D+00   15  0.59D+00
   16  0.10D+05   17  0.10D+05   18  0.10D+05   19  0.10D+05   20  0.10D+05
   21  0.10D+05   22  0.10D+05   23  0.10D+05   24  0.10D+05
 eigenvalue=
 -0.6650 -0.0472 -0.0472 -0.0472  0.0270  0.0270  1.7100  1.9162  1.9162
 (warning) DOS window (-1,0) reset to (-1.1650,0.6517)
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.4185 -0.0983 -0.0478 -0.0233  0.0207  0.0764  0.8953  1.2978  1.7913
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.2385 -0.1775  0.0760  0.0911  0.0911  0.1426  0.6393  1.0060  1.0060

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.151736;  11.000000 electrons
         Sum occ. bands:   -0.756621, incl. Bloechl correction: -0.006665

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local
   1    9.841193    3.278284    6.562910

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:Cu      
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.488976   -0.340549    4.690000    4.656558    4.500000    4.656558
 1     0    0.510458   -0.159503    4.420000    4.397640    4.250000    4.397640
 2     0    8.810553   -0.031904    3.880000    3.867564    3.147584    3.867564
 3     0    0.024903   -0.052671    4.120000    4.125488    4.102416    4.125488
 4     0    0.005037   -0.038215    5.100000    5.087104    5.077979    5.087104
 5     0    0.001265   -0.016835    6.100000    6.067507    6.062833    6.067507

 Harris energy:
 sumev=       -0.756621  val*vef=    -189.393464   sumtv=     188.636843
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -130.258100     utot=   -6534.971246    ehar=   -3304.835864

 ekin=3344.586423  rho*v=-6649.079139  ehf=-3304.835864  ehks=-3304.492716
  
 mixrho:  sought 8 iter from file mixm; read 8.  RMS DQ=3.77e-2
 charges:       old           new         screened      rms diff       lin mix
 smooth       4.729619      4.437091      4.437091      0.045494      4.437091
 site    1    6.270381      6.562910      6.562910      0.014896      6.562910
 AMIX: condition of normal eqns >100000. Reducing nmix to 7
 AMIX: condition of normal eqns >100000. Reducing nmix to 6
 AMIX: condition of normal eqns >100000. Reducing nmix to 5
 AMIX: condition of normal eqns >100000. Reducing nmix to 4
 AMIX: condition of normal eqns >100000. Reducing nmix to 3
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=2678  beta=1  tm=5  rmsdel=3.77e-2
   tj: 1.00000
 unscreened rms difference:  smooth  0.045494   local  0.014896

 iors  : write restart file (binary, mesh density) 

   it  1  of 20    ehf=   -3304.835864   ehk=   -3304.492716
h ehf=-3304.8358636 ehk=-3304.4927165

 --- BNDFP:  begin iteration 2 of 20 ---

 avg es pot at rmt= 0.584471  avg sphere pot= 0.639061  vconst=-0.584471

 site  1  z= 29.0  rmt= 2.28000  nr=393   a=0.025  nlml=36  rg=0.570  Vfloat=T
 sm core charge = 0.295519 (sphere) + 0.00535 (spillout) = 0.300868
 potential shift to crystal energy zero:    0.000102


 subzi: tetrahedron integration of bands; tetrahedron integration of density

 Start first of two band passes ...
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.94D-06    2  0.51D-03    3  0.51D-03    4  0.51D-03    5  0.34D-01
    6  0.34D-01    7  0.45D-01    8  0.78D-01    9  0.78D-01   10  0.78D-01
   11  0.36D+00   12  0.36D+00   13  0.36D+00   14  0.42D+00   15  0.53D+00
   16  0.10D+05   17  0.10D+05   18  0.10D+05   19  0.10D+05   20  0.10D+05
   21  0.10D+05   22  0.10D+05   23  0.10D+05   24  0.10D+05
 eigenvalue=
 -0.7177 -0.2504 -0.2504 -0.2504 -0.1871 -0.1871  1.6698  1.8595  1.8595
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.4818 -0.2794 -0.2464 -0.2286 -0.1918 -0.1374  0.8081  1.2272  1.7467
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.3885 -0.3555 -0.1451 -0.1336 -0.1336  0.0874  0.5208  0.9468  0.9468

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.013372;  11.000000 electrons
         Sum occ. bands:   -2.748811, incl. Bloechl correction: -0.009406

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.94D-06    2  0.51D-03    3  0.51D-03    4  0.51D-03    5  0.34D-01
    6  0.34D-01    7  0.45D-01    8  0.78D-01    9  0.78D-01   10  0.78D-01
   11  0.36D+00   12  0.36D+00   13  0.36D+00   14  0.42D+00   15  0.53D+00
   16  0.10D+05   17  0.10D+05   18  0.10D+05   19  0.10D+05   20  0.10D+05
   21  0.10D+05   22  0.10D+05   23  0.10D+05   24  0.10D+05
 eigenvalue=
 -0.7177 -0.2504 -0.2504 -0.2504 -0.1871 -0.1871  1.6698  1.8595  1.8595
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.4818 -0.2794 -0.2464 -0.2286 -0.1918 -0.1374  0.8081  1.2272  1.7467
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.3885 -0.3555 -0.1451 -0.1336 -0.1336  0.0874  0.5208  0.9468  0.9468

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.013372;  11.000000 electrons
         Sum occ. bands:   -2.748811, incl. Bloechl correction: -0.009406

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local
   1   10.060610    2.897891    7.162718

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:Cu      
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.457684   -0.418743    4.656558    4.654853    4.500000    4.654853
 1     0    0.388448   -0.280459    4.397640    4.380403    4.250000    4.380403
 2     0    9.192066   -0.231959    3.867564    3.878939    3.147584    3.878939
 3     0    0.017427   -0.238670    4.125488    4.121226    4.102416    4.121226
 4     0    0.003909   -0.236466    5.087104    5.085272    5.077979    5.085272
 5     0    0.001076   -0.220914    6.067507    6.066539    6.062833    6.066539

 Harris energy:
 sumev=       -2.748811  val*vef=    -190.117846   sumtv=     187.369035
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.921166     utot=   -6533.968244    ehar=   -3304.763735

 ekin=3359.119828  rho*v=-6663.883561  ehf=-3304.763735  ehks=-3304.763733
  
 mixrho:  sought 8 iter from file mixm; read 8.  RMS DQ=1.73e-5  last it=3.77e-2
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.837641      3.837281      3.837281      0.000011      3.837281
 site    1    7.162359      7.162718      7.162718      0.000020      7.162718
 AMIX: condition of normal eqns >100000. Reducing nmix to 7
 AMIX: condition of normal eqns >100000. Reducing nmix to 6
 AMIX: condition of normal eqns >100000. Reducing nmix to 5
 AMIX: condition of normal eqns >100000. Reducing nmix to 4
 AMIX: condition of normal eqns >100000. Reducing nmix to 3
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=2678  beta=1  tm=5  rmsdel=1.73e-5
   tj:-0.00022
 unscreened rms difference:  smooth  0.000011   local  0.000020

 iors  : write restart file (binary, mesh density) 

   it  2  of 20    ehf=   -3304.763735   ehk=   -3304.763733
 From last iter    ehf=   -3304.835864   ehk=   -3304.492716
 diffe(q)=  0.072128 (0.000017)    tol= 0.000010 (0.000010)   more=T
i ehf=-3304.7637353 ehk=-3304.7637328

 --- BNDFP:  begin iteration 3 of 20 ---

 avg es pot at rmt= 0.584421  avg sphere pot= 0.639112  vconst=-0.584421

 site  1  z= 29.0  rmt= 2.28000  nr=393   a=0.025  nlml=36  rg=0.570  Vfloat=T
 sm core charge = 0.295519 (sphere) + 0.00535 (spillout) = 0.300868
 potential shift to crystal energy zero:    0.000102


 subzi: tetrahedron integration of bands; tetrahedron integration of density

 Start first of two band passes ...
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.94D-06    2  0.51D-03    3  0.51D-03    4  0.51D-03    5  0.35D-01
    6  0.35D-01    7  0.45D-01    8  0.79D-01    9  0.79D-01   10  0.79D-01
   11  0.36D+00   12  0.36D+00   13  0.36D+00   14  0.42D+00   15  0.53D+00
   16  0.10D+05   17  0.10D+05   18  0.10D+05   19  0.10D+05   20  0.10D+05
   21  0.10D+05   22  0.10D+05   23  0.10D+05   24  0.10D+05
 eigenvalue=
 -0.7177 -0.2504 -0.2504 -0.2504 -0.1872 -0.1872  1.6705  1.8605  1.8605
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.4818 -0.2794 -0.2465 -0.2287 -0.1919 -0.1375  0.8075  1.2266  1.7471
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.3886 -0.3555 -0.1452 -0.1336 -0.1336  0.0874  0.5202  0.9471  0.9471

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.013455;  11.000000 electrons
         Sum occ. bands:   -2.749399, incl. Bloechl correction: -0.009403

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.94D-06    2  0.51D-03    3  0.51D-03    4  0.51D-03    5  0.35D-01
    6  0.35D-01    7  0.45D-01    8  0.79D-01    9  0.79D-01   10  0.79D-01
   11  0.36D+00   12  0.36D+00   13  0.36D+00   14  0.42D+00   15  0.53D+00
   16  0.10D+05   17  0.10D+05   18  0.10D+05   19  0.10D+05   20  0.10D+05
   21  0.10D+05   22  0.10D+05   23  0.10D+05   24  0.10D+05
 eigenvalue=
 -0.7177 -0.2504 -0.2504 -0.2504 -0.1872 -0.1872  1.6705  1.8605  1.8605
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.4818 -0.2794 -0.2465 -0.2287 -0.1919 -0.1375  0.8075  1.2266  1.7471
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.3886 -0.3555 -0.1452 -0.1336 -0.1336  0.0874  0.5202  0.9471  0.9471

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.013455;  11.000000 electrons
         Sum occ. bands:   -2.749399, incl. Bloechl correction: -0.009403

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local
   1   10.060799    2.898227    7.162572

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:Cu      
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.457632   -0.418755    4.654853    4.654848    4.500000    4.654848
 1     0    0.388387   -0.280487    4.380403    4.380393    4.250000    4.380393
 2     0    9.192361   -0.232021    3.878939    3.878939    3.147584    3.878939
 3     0    0.017434   -0.238697    4.121226    4.121224    4.102416    4.121224
 4     0    0.003909   -0.236459    5.085272    5.085272    5.077979    5.085272
 5     0    0.001076   -0.220909    6.066539    6.066539    6.062833    6.066539

 Harris energy:
 sumev=       -2.749399  val*vef=    -190.115850   sumtv=     187.366451
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.920980     utot=   -6533.965861    ehar=   -3304.763750

 ekin=3359.131366  rho*v=-6663.895117  ehf=-3304.763750  ehks=-3304.763750
  
 mixrho:  sought 8 iter from file mixm; read 8.  RMS DQ=2.03e-5  last it=1.73e-5
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.837147      3.837428      3.837428      0.000010      3.837428
 site    1    7.162853      7.162572      7.162572      0.000021      7.162572
 AMIX: condition of normal eqns >100000. Reducing nmix to 7
 AMIX: condition of normal eqns >100000. Reducing nmix to 6
 AMIX: condition of normal eqns >100000. Reducing nmix to 5
 AMIX: condition of normal eqns >100000. Reducing nmix to 4
 AMIX: condition of normal eqns >100000. Reducing nmix to 3
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=2678  beta=1  tm=5  rmsdel=2.03e-5
   tj: 0.54081
 unscreened rms difference:  smooth  0.000010   local  0.000021

 iors  : write restart file (binary, mesh density) 

   it  3  of 20    ehf=   -3304.763750   ehk=   -3304.763750
 From last iter    ehf=   -3304.763735   ehk=   -3304.763733
 diffe(q)= -0.000015 (0.000020)    tol= 0.000010 (0.000010)   more=T
i ehf=-3304.7637504 ehk=-3304.7637503

 --- BNDFP:  begin iteration 4 of 20 ---

 avg es pot at rmt= 0.584441  avg sphere pot= 0.639088  vconst=-0.584441

 site  1  z= 29.0  rmt= 2.28000  nr=393   a=0.025  nlml=36  rg=0.570  Vfloat=T
 sm core charge = 0.295519 (sphere) + 0.00535 (spillout) = 0.300868
 potential shift to crystal energy zero:    0.000102


 subzi: tetrahedron integration of bands; tetrahedron integration of density

 Start first of two band passes ...
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.94D-06    2  0.51D-03    3  0.51D-03    4  0.51D-03    5  0.35D-01
    6  0.35D-01    7  0.45D-01    8  0.79D-01    9  0.79D-01   10  0.79D-01
   11  0.36D+00   12  0.36D+00   13  0.36D+00   14  0.42D+00   15  0.53D+00
   16  0.10D+05   17  0.10D+05   18  0.10D+05   19  0.10D+05   20  0.10D+05
   21  0.10D+05   22  0.10D+05   23  0.10D+05   24  0.10D+05
 eigenvalue=
 -0.7177 -0.2504 -0.2504 -0.2504 -0.1871 -0.1871  1.6705  1.8605  1.8605
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.4818 -0.2794 -0.2464 -0.2286 -0.1918 -0.1375  0.8075  1.2266  1.7471
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.3885 -0.3554 -0.1452 -0.1336 -0.1336  0.0874  0.5203  0.9471  0.9471

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.013421;  11.000000 electrons
         Sum occ. bands:   -2.748892, incl. Bloechl correction: -0.009402

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.94D-06    2  0.51D-03    3  0.51D-03    4  0.51D-03    5  0.35D-01
    6  0.35D-01    7  0.45D-01    8  0.79D-01    9  0.79D-01   10  0.79D-01
   11  0.36D+00   12  0.36D+00   13  0.36D+00   14  0.42D+00   15  0.53D+00
   16  0.10D+05   17  0.10D+05   18  0.10D+05   19  0.10D+05   20  0.10D+05
   21  0.10D+05   22  0.10D+05   23  0.10D+05   24  0.10D+05
 eigenvalue=
 -0.7177 -0.2504 -0.2504 -0.2504 -0.1871 -0.1871  1.6705  1.8605  1.8605
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.4818 -0.2794 -0.2464 -0.2286 -0.1918 -0.1375  0.8075  1.2266  1.7471
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.3885 -0.3554 -0.1452 -0.1336 -0.1336  0.0874  0.5203  0.9471  0.9471

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.013421;  11.000000 electrons
         Sum occ. bands:   -2.748892, incl. Bloechl correction: -0.009402

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local
   1   10.060738    2.898358    7.162380

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:Cu      
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.457643   -0.418743    4.654848    4.654850    4.500000    4.654850
 1     0    0.388420   -0.280463    4.380393    4.380399    4.250000    4.380399
 2     0    9.192253   -0.231969    3.878939    3.878936    3.147584    3.878936
 3     0    0.017436   -0.238650    4.121224    4.121226    4.102416    4.121226
 4     0    0.003910   -0.236408    5.085272    5.085272    5.077979    5.085272
 5     0    0.001076   -0.220857    6.066539    6.066540    6.062833    6.066540

 Harris energy:
 sumev=       -2.748892  val*vef=    -190.117371   sumtv=     187.368479
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.921120     utot=   -6533.967748    ehar=   -3304.763749

 ekin=3359.126145  rho*v=-6663.889894  ehf=-3304.763749  ehks=-3304.763749
  
 mixrho:  sought 8 iter from file mixm; read 8.  RMS DQ=7.93e-6  last it=2.03e-5
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.837349      3.837620      3.837620      0.000007      3.837620
 site    1    7.162651      7.162380      7.162380      0.000011      7.162380
 AMIX: condition of normal eqns >100000. Reducing nmix to 7
 AMIX: condition of normal eqns >100000. Reducing nmix to 6
 AMIX: condition of normal eqns >100000. Reducing nmix to 5
 AMIX: condition of normal eqns >100000. Reducing nmix to 4
 AMIX: condition of normal eqns >100000. Reducing nmix to 3
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=2678  beta=1  tm=5  rmsdel=7.93e-6
   tj:-0.28582
 unscreened rms difference:  smooth  0.000007   local  0.000011

 iors  : write restart file (binary, mesh density) 

   it  4  of 20    ehf=   -3304.763749   ehk=   -3304.763749
 From last iter    ehf=   -3304.763750   ehk=   -3304.763750
 diffe(q)=  0.000001 (0.000008)    tol= 0.000010 (0.000010)   more=F
c ehf=-3304.7637492 ehk=-3304.7637492
 Exit 0 LMF 
 wkinfo:  used  1003 K  workspace of 80000 K   in  72 K calls
