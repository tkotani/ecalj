HOST_INFORMATION platform: gfortran
HOST_INFORMATION compiler version: gcc バージョン 4.6.1 (Ubuntu/Linaro 4.6.1-9ubuntu3) 
HOST_INFORMATION FFLAGS (<=120): -O2 -fomit-frame-pointer -funroll-loops -ffast-math -ffixed-line-length-132 -DHASIARGC -DHASGETARG -DFDATE -DHASCPUTIME 
HOST_INFORMATION LIBLOC (<=120): /usr/lib/libfftw3.so.3 /usr/lib/liblapack.so.3gf /usr/lib/libblas.so.3gf
HOST_INFORMATION uname -a (<=120): Linux TT4 3.0.0-20-generic #34-Ubuntu SMP Tue May 1 17:24:39 UTC 2012 x86_64 x86_64 x86_64 GNU/Linux
HOST_INFORMATION /etc/issue: Ubuntu 11.10 \n \l
HOST_INFORMATION git branch: refs/heads/newaniso
HOST_INFORMATION git commit: ff4f0b57505a84e0eaabe44d640190810c88782c
HOST_INFORMATION linked at: Tue Jun 19 14:13:38 JST 2012
 -----------------------  START LMF      -----------------------
 ptenv() is called with EXT=cu
 ptenv() not supported, but continue.
 HEADER Cu
 Cu       xxx            1           1

 rdctrl: reset global max nl from 5 to 6
  mxcst switch =           1           0 F F F
  LMF  vn 7.00(LMF 7.0)  verb 31,20
 special:  APW basis
 bz:       metal(3), tetra, invit, fixed-spin-mom 
 goto setcg
 lattic:

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
 zzz nclass=           1
 
 lstar xxx=          -2
 BZMESH:  29 irreducible QP from 512 ( 8 8 8 )  shift= F F F
 lstar xxx=          -2
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
 rdovfa: expected Cu,      read Cu       with rmt=  2.2800  mesh   655  0.015

 ovlpfa: overlap smooth part of FA densities
 site   1  spec  1  pos  0.0000  0.0000  0.0000  Qsmooth 4.728719
 total smooth Q = 4.728719

 Free atom and overlapped crystal site charges:
   ib    true(FA)    smooth(FA)  true(OV)    smooth(OV)    local
    1    9.758197    3.486916   10.209422    3.938142    6.271281

 Smooth charge on mesh:            4.728719
 Sum of local charges:             6.271281
 Total valence charge:            11.000000
 Sum of core charges:             18.000000
 Sum of nuclear charges:         -29.000000
 Homogeneous background:           0.000000
 Deviation from neutrality:        0.000000

 --- BNDFP:  begin iteration 1 of 20 ---
 ttt nevmx w=           0  5.00000000000000010E-003

 avg es pot at rmt= 0.554288  avg sphere pot= 0.625098  vconst=-0.554288
 smvxcm: all smrho_w is positive
  i job kmax lfltwf(FRZWF see locpot.F)=           0           1           5 T

 site  1  z= 29.0  rmt= 2.28000  nr=655   a=0.015  nlml=36  rg=0.570  Vfloat=T
 sm core charge = 0.295521 (sphere) + 0.00535 (spillout) = 0.30087
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===
 potential shift to crystal energy zero:    0.000086


 subzi: tetrahedron integration of bands; tetrahedron integration of density

 Start first of two band passes ...
 end of suham2
 -------- qplist --------
    1   0.000   0.000   0.000
    2  -0.125   0.125   0.125
    3  -0.250   0.250   0.250
    4  -0.375   0.375   0.375
    5  -0.500   0.500   0.500
    6   0.000   0.000   0.250
    7  -0.125   0.125   0.375
    8  -0.250   0.250   0.500
    9  -0.375   0.375   0.625
   10  -0.500   0.500   0.750
   11  -0.625   0.625   0.875
   12  -0.750   0.750   1.000
   13   0.000   0.000   0.500
   14  -0.125   0.125   0.625
   15  -0.250   0.250   0.750
   16  -0.375   0.375   0.875
   17  -0.500   0.500   1.000
   18   0.000   0.000   0.750
   19  -0.125   0.125   0.875
   20  -0.250   0.250   1.000
   21   0.000   0.000   1.000
   22   0.000   0.250   0.500
   23  -0.125   0.375   0.625
   24  -0.250   0.500   0.750
   25   0.000   0.250   0.750
   26  -0.125   0.375   0.875
   27  -0.250   0.500   1.000
   28   0.000   0.250   1.000
   29   0.000   0.500   1.000
 sigmamode= F

  --- Hamiltonian index ---
  ib l  k offl(iorb)+1  offl(iorb)+2*l+1  trim(spec)
  ngrp=     -999999
  end of hambls mode=           0
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.95D-06    2  0.57D-03    3  0.57D-03    4  0.57D-03    5  0.40D-01
    6  0.40D-01    7  0.51D-01    8  0.90D-01    9  0.90D-01   10  0.90D-01
   11  0.40D+00   12  0.40D+00   13  0.40D+00   14  0.47D+00   15  0.59D+00
   ... skip larger eigenvalues ...
 eigenvalue=
 -0.6650 -0.0472 -0.0472 -0.0472  0.0270  0.0270  1.7100  1.9162  1.9162
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.4185 -0.0983 -0.0478 -0.0233  0.0207  0.0764  0.8953  1.2978  1.7913
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.2385 -0.1775  0.0760  0.0911  0.0911  0.1426  0.6393  1.0060  1.0060

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.151737;  11.000000 electrons
         Sum occ. bands:   -0.756612, incl. Bloechl correction: -0.006665

 Saved qp weights ...
 Start second band pass ...
 -------- qplist --------
    1   0.000   0.000   0.000
    2  -0.125   0.125   0.125
    3  -0.250   0.250   0.250
    4  -0.375   0.375   0.375
    5  -0.500   0.500   0.500
    6   0.000   0.000   0.250
    7  -0.125   0.125   0.375
    8  -0.250   0.250   0.500
    9  -0.375   0.375   0.625
   10  -0.500   0.500   0.750
   11  -0.625   0.625   0.875
   12  -0.750   0.750   1.000
   13   0.000   0.000   0.500
   14  -0.125   0.125   0.625
   15  -0.250   0.250   0.750
   16  -0.375   0.375   0.875
   17  -0.500   0.500   1.000
   18   0.000   0.000   0.750
   19  -0.125   0.125   0.875
   20  -0.250   0.250   1.000
   21   0.000   0.000   1.000
   22   0.000   0.250   0.500
   23  -0.125   0.375   0.625
   24  -0.250   0.500   0.750
   25   0.000   0.250   0.750
   26  -0.125   0.375   0.875
   27  -0.250   0.500   1.000
   28   0.000   0.250   1.000
   29   0.000   0.500   1.000
 sigmamode= F
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.95D-06    2  0.57D-03    3  0.57D-03    4  0.57D-03    5  0.40D-01
    6  0.40D-01    7  0.51D-01    8  0.90D-01    9  0.90D-01   10  0.90D-01
   11  0.40D+00   12  0.40D+00   13  0.40D+00   14  0.47D+00   15  0.59D+00
   ... skip larger eigenvalues ...
 eigenvalue=
 -0.6650 -0.0472 -0.0472 -0.0472  0.0270  0.0270  1.7100  1.9162  1.9162
 (warning) DOS window (-1,0) reset to (-1.1650,0.6517)
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.4185 -0.0983 -0.0478 -0.0233  0.0207  0.0764  0.8953  1.2978  1.7913
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.2385 -0.1775  0.0760  0.0911  0.0911  0.1426  0.6393  1.0060  1.0060

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.151737;  11.000000 electrons
         Sum occ. bands:   -0.756612, incl. Bloechl correction: -0.006665

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local
   1    9.841195    3.278269    6.562926

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:Cu      
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.488975   -0.340547    4.690000    4.656558    4.500000    4.656558
 1     0    0.510456   -0.159501    4.420000    4.397639    4.250000    4.397639
 2     0    8.810559   -0.031904    3.880000    3.867564    3.147584    3.867564
 3     0    0.024903   -0.052670    4.120000    4.125488    4.102416    4.125488
 4     0    0.005037   -0.038215    5.100000    5.087103    5.077979    5.087103
 5     0    0.001265   -0.016834    6.100000    6.067507    6.062833    6.067507

 Harris energy:
 sumev=       -0.756612  val*vef=    -189.393534   sumtv=     188.636922
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756520
 rhoeps=    -130.258096     utot=   -6534.971163    ehar=   -3304.835818

 avg es pot at rmt= 0.689901  avg sphere pot= 0.633467  vconst=-0.689901
 smvxcm: all smrho_w is positive
  i job kmax lfltwf(FRZWF see locpot.F)=           0           0           5 F

 site  1  z= 29.0  rmt= 2.28000  nr=655   a=0.015  nlml=36  rg=0.570  Vfloat=F
 sm core charge = 0.295521 (sphere) + 0.00535 (spillout) = 0.30087
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===


 ekin=3344.586478  rho*v=-6649.079153 ehf =-3304.835818  ehks =-3304.492675
 mixrho: sum smrho  init = 0.203205D+03-0.693212D-28 0.203205D+03       0
 mixrho: sum smrnew new  = 0.190672D+03 0.244466D-17 0.190672D+03       0
  
 mixrho: dqsum rmsuns= -0.37134D-02  0.45470D-01  0.13184D-18
 mixrealsmooth= T
 wgtsmooth=  1.72132593164774084E-002
 mixrho:  sought 8 iter from file mixm; read 0.  RMS DQ=3.44e-2
 charges:       old           new         screened      rms diff       lin mix
 smooth       4.728719      4.437074      4.437074      0.045470      4.437074
 site    1    6.271281      6.562926      6.562926      0.014872      6.562926
 AMIX: nmix=0 mmix=8  nelts=5333  beta=1  tm=5  rmsdel=3.44e-2
 mixrealsmooth= T
 smrho qcell: add correction to smrho= -5.57011645696547930E-008 -7.09219692408622063E-010
 mixrho: all smrho is positive for isp=           1

 iors  : write restart file (binary, mesh density) 

   it  1  of 20    ehf=   -3304.835818   ehk=   -3304.492675
h ehf=-3304.8358175 ehk=-3304.4926747

 --- BNDFP:  begin iteration 2 of 20 ---
 ttt nevmx w=           0  5.00000000000000010E-003

 avg es pot at rmt= 0.689901  avg sphere pot= 0.633467  vconst=-0.689901
 smvxcm: all smrho_w is positive
  i job kmax lfltwf(FRZWF see locpot.F)=           0           1           5 T

 site  1  z= 29.0  rmt= 2.28000  nr=655   a=0.015  nlml=36  rg=0.570  Vfloat=T
 sm core charge = 0.295521 (sphere) + 0.00535 (spillout) = 0.30087
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===
 potential shift to crystal energy zero:    0.000110


 subzi: tetrahedron integration of bands; tetrahedron integration of density

 Start first of two band passes ...
 end of suham2
 -------- qplist --------
    1   0.000   0.000   0.000
    2  -0.125   0.125   0.125
    3  -0.250   0.250   0.250
    4  -0.375   0.375   0.375
    5  -0.500   0.500   0.500
    6   0.000   0.000   0.250
    7  -0.125   0.125   0.375
    8  -0.250   0.250   0.500
    9  -0.375   0.375   0.625
   10  -0.500   0.500   0.750
   11  -0.625   0.625   0.875
   12  -0.750   0.750   1.000
   13   0.000   0.000   0.500
   14  -0.125   0.125   0.625
   15  -0.250   0.250   0.750
   16  -0.375   0.375   0.875
   17  -0.500   0.500   1.000
   18   0.000   0.000   0.750
   19  -0.125   0.125   0.875
   20  -0.250   0.250   1.000
   21   0.000   0.000   1.000
   22   0.000   0.250   0.500
   23  -0.125   0.375   0.625
   24  -0.250   0.500   0.750
   25   0.000   0.250   0.750
   26  -0.125   0.375   0.875
   27  -0.250   0.500   1.000
   28   0.000   0.250   1.000
   29   0.000   0.500   1.000
 sigmamode= F
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.93D-06    2  0.33D-03    3  0.33D-03    4  0.33D-03    5  0.21D-01
    6  0.21D-01    7  0.31D-01    8  0.45D-01    9  0.45D-01   10  0.45D-01
   11  0.24D+00   12  0.24D+00   13  0.24D+00   14  0.27D+00   15  0.35D+00
   ... skip larger eigenvalues ...
 eigenvalue=
 -0.9319 -0.9319 -0.9319 -0.8970 -0.8970 -0.8327  1.6065  1.7159  1.7159
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.9461 -0.9274 -0.9181 -0.8984 -0.8907 -0.5508  0.6217  1.0694  1.6576
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.9790 -0.9711 -0.8735 -0.8689 -0.8689 -0.0335  0.2689  0.8090  0.8090

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.252804;  11.000000 electrons
         Sum occ. bands:   -9.639146, incl. Bloechl correction: -0.013835

 Saved qp weights ...
 Start second band pass ...
 -------- qplist --------
    1   0.000   0.000   0.000
    2  -0.125   0.125   0.125
    3  -0.250   0.250   0.250
    4  -0.375   0.375   0.375
    5  -0.500   0.500   0.500
    6   0.000   0.000   0.250
    7  -0.125   0.125   0.375
    8  -0.250   0.250   0.500
    9  -0.375   0.375   0.625
   10  -0.500   0.500   0.750
   11  -0.625   0.625   0.875
   12  -0.750   0.750   1.000
   13   0.000   0.000   0.500
   14  -0.125   0.125   0.625
   15  -0.250   0.250   0.750
   16  -0.375   0.375   0.875
   17  -0.500   0.500   1.000
   18   0.000   0.000   0.750
   19  -0.125   0.125   0.875
   20  -0.250   0.250   1.000
   21   0.000   0.000   1.000
   22   0.000   0.250   0.500
   23  -0.125   0.375   0.625
   24  -0.250   0.500   0.750
   25   0.000   0.250   0.750
   26  -0.125   0.375   0.875
   27  -0.250   0.500   1.000
   28   0.000   0.250   1.000
   29   0.000   0.500   1.000
 sigmamode= F
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.93D-06    2  0.33D-03    3  0.33D-03    4  0.33D-03    5  0.21D-01
    6  0.21D-01    7  0.31D-01    8  0.45D-01    9  0.45D-01   10  0.45D-01
   11  0.24D+00   12  0.24D+00   13  0.24D+00   14  0.27D+00   15  0.35D+00
   ... skip larger eigenvalues ...
 eigenvalue=
 -0.9319 -0.9319 -0.9319 -0.8970 -0.8970 -0.8327  1.6065  1.7159  1.7159
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.9461 -0.9274 -0.9181 -0.8984 -0.8907 -0.5508  0.6217  1.0694  1.6576
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.9790 -0.9711 -0.8735 -0.8689 -0.8689 -0.0335  0.2689  0.8090  0.8090

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.252804;  11.000000 electrons
         Sum occ. bands:   -9.639146, incl. Bloechl correction: -0.013835

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local
   1   10.453936    1.990699    8.463236

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:Cu      
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.421972   -0.548056    4.656558    4.675732    4.500000    4.675732
 1     0    0.220875   -0.490003    4.397639    4.366703    4.250000    4.366703
 2     0    9.803395   -0.912551    3.867564    3.908754    3.147584    3.908754
 3     0    0.005603   -0.826635    4.125488    4.108379    4.102416    4.108379
 4     0    0.001555   -0.899999    5.087103    5.079159    5.077979    5.079159
 5     0    0.000535   -0.897556    6.067507    6.063251    6.062833    6.063251

 Harris energy:
 sumev=       -9.639146  val*vef=    -181.399258   sumtv=     171.760112
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756520
 rhoeps=    -128.563181     utot=   -6520.515972    ehar=   -3305.562522

 avg es pot at rmt= 0.372049  avg sphere pot= 0.636842  vconst=-0.372049
 smvxcm: all smrho_w is positive
  i job kmax lfltwf(FRZWF see locpot.F)=           0           0           5 F

 site  1  z= 29.0  rmt= 2.28000  nr=655   a=0.015  nlml=36  rg=0.570  Vfloat=F
 sm core charge = 0.295521 (sphere) + 0.00535 (spillout) = 0.30087
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===


 ekin=3400.623324  rho*v=-6703.686458 ehf =-3305.562522  ehks =-3303.063134
 mixrho: sum smrho  init = 0.190672D+03-0.545596D-27 0.190672D+03       0
 mixrho: sum smrnew new  = 0.109011D+03 0.364971D-17 0.109011D+03       0
  
 mixrho: dqsum rmsuns= -0.24196D-01  0.24682D-01 -0.20633D-19
 mixrealsmooth= T
 wgtsmooth=  1.72132593164774084E-002
 mixrho:  sought 8 iter from file mixm; read 1.  RMS DQ=1.01e-1  last it=3.44e-2
 charges:       old           new         screened      rms diff       lin mix
 smooth       4.437074      2.536763      2.536763      0.024682      2.536763
 site    1    6.562926      8.463236      8.463236      0.078768      8.463236
 AMIX: nmix=1 mmix=8  nelts=5333  beta=1  tm=5  rmsdel=1.01e-1
   tj: 0.81751
 mixrealsmooth= T
 smrho qcell: add correction to smrho= -1.98727940770027089E-008 -2.53031996933688123E-010
 mixrho: all smrho is positive for isp=           1

 iors  : write restart file (binary, mesh density) 

   it  2  of 20    ehf=   -3305.562522   ehk=   -3303.063134
 From last iter    ehf=   -3304.835818   ehk=   -3304.492675
 diffe(q)= -0.726704 (0.101381)    tol= 0.000010 (0.000010)   more=T
i ehf=-3305.5625216 ehk=-3303.063134

 --- BNDFP:  begin iteration 3 of 20 ---
 ttt nevmx w=           0  5.00000000000000010E-003

 avg es pot at rmt= 0.631894  avg sphere pot= 0.634083  vconst=-0.631894
 smvxcm: all smrho_w is positive
  i job kmax lfltwf(FRZWF see locpot.F)=           0           1           5 T

 site  1  z= 29.0  rmt= 2.28000  nr=655   a=0.015  nlml=36  rg=0.570  Vfloat=T
 sm core charge = 0.295521 (sphere) + 0.00535 (spillout) = 0.30087
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===
 potential shift to crystal energy zero:    0.000105


 subzi: tetrahedron integration of bands; tetrahedron integration of density

 Start first of two band passes ...
 end of suham2
 -------- qplist --------
    1   0.000   0.000   0.000
    2  -0.125   0.125   0.125
    3  -0.250   0.250   0.250
    4  -0.375   0.375   0.375
    5  -0.500   0.500   0.500
    6   0.000   0.000   0.250
    7  -0.125   0.125   0.375
    8  -0.250   0.250   0.500
    9  -0.375   0.375   0.625
   10  -0.500   0.500   0.750
   11  -0.625   0.625   0.875
   12  -0.750   0.750   1.000
   13   0.000   0.000   0.500
   14  -0.125   0.125   0.625
   15  -0.250   0.250   0.750
   16  -0.375   0.375   0.875
   17  -0.500   0.500   1.000
   18   0.000   0.000   0.750
   19  -0.125   0.125   0.875
   20  -0.250   0.250   1.000
   21   0.000   0.000   1.000
   22   0.000   0.250   0.500
   23  -0.125   0.375   0.625
   24  -0.250   0.500   0.750
   25   0.000   0.250   0.750
   26  -0.125   0.375   0.875
   27  -0.250   0.500   1.000
   28   0.000   0.250   1.000
   29   0.000   0.500   1.000
 sigmamode= F
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.94D-06    2  0.44D-03    3  0.44D-03    4  0.44D-03    5  0.29D-01
    6  0.29D-01    7  0.40D-01    8  0.67D-01    9  0.67D-01   10  0.67D-01
   11  0.31D+00   12  0.31D+00   13  0.31D+00   14  0.36D+00   15  0.46D+00
   ... skip larger eigenvalues ...
 eigenvalue=
 -0.7667 -0.5029 -0.5029 -0.5029 -0.4521 -0.4521  1.6438  1.8052  1.8052
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.5792 -0.5061 -0.4873 -0.4843 -0.4551 -0.3795  0.7220  1.1562  1.7106
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.5908 -0.5812 -0.4182 -0.4098 -0.4098  0.0359  0.4030  0.8900  0.8900

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.142459;  11.000000 electrons
         Sum occ. bands:   -5.274909, incl. Bloechl correction: -0.012364

 Saved qp weights ...
 Start second band pass ...
 -------- qplist --------
    1   0.000   0.000   0.000
    2  -0.125   0.125   0.125
    3  -0.250   0.250   0.250
    4  -0.375   0.375   0.375
    5  -0.500   0.500   0.500
    6   0.000   0.000   0.250
    7  -0.125   0.125   0.375
    8  -0.250   0.250   0.500
    9  -0.375   0.375   0.625
   10  -0.500   0.500   0.750
   11  -0.625   0.625   0.875
   12  -0.750   0.750   1.000
   13   0.000   0.000   0.500
   14  -0.125   0.125   0.625
   15  -0.250   0.250   0.750
   16  -0.375   0.375   0.875
   17  -0.500   0.500   1.000
   18   0.000   0.000   0.750
   19  -0.125   0.125   0.875
   20  -0.250   0.250   1.000
   21   0.000   0.000   1.000
   22   0.000   0.250   0.500
   23  -0.125   0.375   0.625
   24  -0.250   0.500   0.750
   25   0.000   0.250   0.750
   26  -0.125   0.375   0.875
   27  -0.250   0.500   1.000
   28   0.000   0.250   1.000
   29   0.000   0.500   1.000
 sigmamode= F
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.94D-06    2  0.44D-03    3  0.44D-03    4  0.44D-03    5  0.29D-01
    6  0.29D-01    7  0.40D-01    8  0.67D-01    9  0.67D-01   10  0.67D-01
   11  0.31D+00   12  0.31D+00   13  0.31D+00   14  0.36D+00   15  0.46D+00
   ... skip larger eigenvalues ...
 eigenvalue=
 -0.7667 -0.5029 -0.5029 -0.5029 -0.4521 -0.4521  1.6438  1.8052  1.8052
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.5792 -0.5061 -0.4873 -0.4843 -0.4551 -0.3795  0.7220  1.1562  1.7106
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.5908 -0.5812 -0.4182 -0.4098 -0.4098  0.0359  0.4030  0.8900  0.8900

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.142459;  11.000000 electrons
         Sum occ. bands:   -5.274909, incl. Bloechl correction: -0.012364

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local
   1   10.273018    2.478563    7.794455

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:Cu      
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.429288   -0.484088    4.675732    4.658674    4.500000    4.658674
 1     0    0.283944   -0.393195    4.366703    4.366592    4.250000    4.366592
 2     0    9.545093   -0.482273    3.908754    3.892305    3.147584    3.892305
 3     0    0.011079   -0.468597    4.108379    4.115929    4.102416    4.115929
 4     0    0.002780   -0.483678    5.079159    5.082929    5.077979    5.082929
 5     0    0.000835   -0.473242    6.063251    6.065297    6.062833    6.065297

 Harris energy:
 sumev=       -5.274909  val*vef=    -188.047626   sumtv=     182.772718
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756520
 rhoeps=    -129.432876     utot=   -6529.984448    ehar=   -3304.888087

 avg es pot at rmt= 0.477136  avg sphere pot= 0.641278  vconst=-0.477136
 smvxcm: all smrho_w is positive
  i job kmax lfltwf(FRZWF see locpot.F)=           0           0           5 F

 site  1  z= 29.0  rmt= 2.28000  nr=655   a=0.015  nlml=36  rg=0.570  Vfloat=F
 sm core charge = 0.295521 (sphere) + 0.00535 (spillout) = 0.30087
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===


 ekin=3376.827204  rho*v=-6681.250065 ehf =-3304.888087  ehks =-3304.422861
 mixrho: sum smrho  init = 0.175769D+03 0.802469D-27 0.175769D+03       0
 mixrho: sum smrnew new  = 0.137750D+03 0.339595D-17 0.137750D+03       0
  
 mixrho: dqsum rmsuns= -0.11265D-01  0.11425D-01 -0.12824D-19
 mixrealsmooth= T
 wgtsmooth=  1.72132593164774084E-002
 mixrho:  sought 8 iter from file mixm; read 2.  RMS DQ=4.17e-2  last it=1.01e-1
 charges:       old           new         screened      rms diff       lin mix
 smooth       4.090277      3.205545      3.205545      0.011425      3.205545
 site    1    6.909723      7.794455      7.794455      0.034207      7.794455
 AMIX: nmix=2 mmix=8  nelts=5333  beta=1  tm=5  rmsdel=4.17e-2
   tj:-0.73043  -0.05116
 mixrealsmooth= T
 smrho qcell: add correction to smrho= -2.40979860954837477E-007 -3.06829604342187205E-009
 mixrho: all smrho is positive for isp=           1

 iors  : write restart file (binary, mesh density) 

   it  3  of 20    ehf=   -3304.888087   ehk=   -3304.422861
 From last iter    ehf=   -3305.562522   ehk=   -3303.063134
 diffe(q)=  0.674435 (0.041697)    tol= 0.000010 (0.000010)   more=T
i ehf=-3304.8880869 ehk=-3304.4228609

 --- BNDFP:  begin iteration 4 of 20 ---
 ttt nevmx w=           0  5.00000000000000010E-003

 avg es pot at rmt= 0.543010  avg sphere pot= 0.644918  vconst=-0.543010
 smvxcm: all smrho_w is positive
  i job kmax lfltwf(FRZWF see locpot.F)=           0           1           5 T

 site  1  z= 29.0  rmt= 2.28000  nr=655   a=0.015  nlml=36  rg=0.570  Vfloat=T
 sm core charge = 0.295521 (sphere) + 0.00535 (spillout) = 0.30087
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===
 potential shift to crystal energy zero:    0.000100


 subzi: tetrahedron integration of bands; tetrahedron integration of density

 Start first of two band passes ...
 end of suham2
 -------- qplist --------
    1   0.000   0.000   0.000
    2  -0.125   0.125   0.125
    3  -0.250   0.250   0.250
    4  -0.375   0.375   0.375
    5  -0.500   0.500   0.500
    6   0.000   0.000   0.250
    7  -0.125   0.125   0.375
    8  -0.250   0.250   0.500
    9  -0.375   0.375   0.625
   10  -0.500   0.500   0.750
   11  -0.625   0.625   0.875
   12  -0.750   0.750   1.000
   13   0.000   0.000   0.500
   14  -0.125   0.125   0.625
   15  -0.250   0.250   0.750
   16  -0.375   0.375   0.875
   17  -0.500   0.500   1.000
   18   0.000   0.000   0.750
   19  -0.125   0.125   0.875
   20  -0.250   0.250   1.000
   21   0.000   0.000   1.000
   22   0.000   0.250   0.500
   23  -0.125   0.375   0.625
   24  -0.250   0.500   0.750
   25   0.000   0.250   0.750
   26  -0.125   0.375   0.875
   27  -0.250   0.500   1.000
   28   0.000   0.250   1.000
   29   0.000   0.500   1.000
 sigmamode= F
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.94D-06    2  0.56D-03    3  0.56D-03    4  0.56D-03    5  0.39D-01
    6  0.39D-01    7  0.50D-01    8  0.91D-01    9  0.91D-01   10  0.91D-01
   11  0.40D+00   12  0.40D+00   13  0.40D+00   14  0.46D+00   15  0.59D+00
   ... skip larger eigenvalues ...
 eigenvalue=
 -0.6761 -0.0678 -0.0678 -0.0678  0.0065  0.0065  1.6999  1.9075  1.9075
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.4300 -0.1174 -0.0677 -0.0436 -0.0000  0.0556  0.8818  1.2855  1.7822
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.2547 -0.1966  0.0552  0.0694  0.0694  0.1313  0.6241  0.9953  0.9953

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.133045;  11.000000 electrons
         Sum occ. bands:   -0.962555, incl. Bloechl correction: -0.006820

 Saved qp weights ...
 Start second band pass ...
 -------- qplist --------
    1   0.000   0.000   0.000
    2  -0.125   0.125   0.125
    3  -0.250   0.250   0.250
    4  -0.375   0.375   0.375
    5  -0.500   0.500   0.500
    6   0.000   0.000   0.250
    7  -0.125   0.125   0.375
    8  -0.250   0.250   0.500
    9  -0.375   0.375   0.625
   10  -0.500   0.500   0.750
   11  -0.625   0.625   0.875
   12  -0.750   0.750   1.000
   13   0.000   0.000   0.500
   14  -0.125   0.125   0.625
   15  -0.250   0.250   0.750
   16  -0.375   0.375   0.875
   17  -0.500   0.500   1.000
   18   0.000   0.000   0.750
   19  -0.125   0.125   0.875
   20  -0.250   0.250   1.000
   21   0.000   0.000   1.000
   22   0.000   0.250   0.500
   23  -0.125   0.375   0.625
   24  -0.250   0.500   0.750
   25   0.000   0.250   0.750
   26  -0.125   0.375   0.875
   27  -0.250   0.500   1.000
   28   0.000   0.250   1.000
   29   0.000   0.500   1.000
 sigmamode= F
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.94D-06    2  0.56D-03    3  0.56D-03    4  0.56D-03    5  0.39D-01
    6  0.39D-01    7  0.50D-01    8  0.91D-01    9  0.91D-01   10  0.91D-01
   11  0.40D+00   12  0.40D+00   13  0.40D+00   14  0.46D+00   15  0.59D+00
   ... skip larger eigenvalues ...
 eigenvalue=
 -0.6761 -0.0678 -0.0678 -0.0678  0.0065  0.0065  1.6999  1.9075  1.9075
 (warning) DOS window (-1,0) reset to (-1.1761,0.6330)
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.4300 -0.1174 -0.0677 -0.0436 -0.0000  0.0556  0.8818  1.2855  1.7822
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.2547 -0.1966  0.0552  0.0694  0.0694  0.1313  0.6241  0.9953  0.9953

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.133045;  11.000000 electrons
         Sum occ. bands:   -0.962555, incl. Bloechl correction: -0.006820

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local
   1    9.855952    3.261239    6.594714

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:Cu      
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.487308   -0.353357    4.658674    4.656541    4.500000    4.656541
 1     0    0.502065   -0.174962    4.366592    4.396391    4.250000    4.396391
 2     0    8.835915   -0.051917    3.892305    3.868334    3.147584    3.868334
 3     0    0.024398   -0.071741    4.115929    4.125224    4.102416    4.125224
 4     0    0.004985   -0.058075    5.082929    5.086988    5.077979    5.086988
 5     0    0.001282   -0.037146    6.065297    6.067445    6.062833    6.067445

 Harris energy:
 sumev=       -0.962555  val*vef=    -190.404774   sumtv=     189.442219
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756520
 rhoeps=    -130.256522     utot=   -6535.769395    ehar=   -3304.827179

 avg es pot at rmt= 0.683774  avg sphere pot= 0.633966  vconst=-0.683774
 smvxcm: all smrho_w is positive
  i job kmax lfltwf(FRZWF see locpot.F)=           0           0           5 F

 site  1  z= 29.0  rmt= 2.28000  nr=655   a=0.015  nlml=36  rg=0.570  Vfloat=F
 sm core charge = 0.295521 (sphere) + 0.00535 (spillout) = 0.30087
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===


 ekin=3345.170049  rho*v=-6649.688679 ehf =-3304.827179  ehks =-3304.518630
 mixrho: sum smrho  init = 0.156035D+03-0.586321D-27 0.156035D+03       0
 mixrho: sum smrnew new  = 0.189306D+03-0.159652D-16 0.189306D+03       0
  
 mixrho: dqsum rmsuns=  0.98581D-02  0.10080D-01 -0.24287D-19
 mixrealsmooth= T
 wgtsmooth=  1.72132593164774084E-002
 mixrho:  sought 8 iter from file mixm; read 3.  RMS DQ=3.25e-2  last it=4.17e-2
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.631043      4.405287      4.405287      0.010080      4.405287
 site    1    7.368957      6.594714      6.594714      0.028507      6.594714
 AMIX: nmix=3 mmix=8  nelts=5333  beta=1  tm=5  rmsdel=3.25e-2
   tj: 0.72976  -0.16331  -0.00000
 mixrealsmooth= T
 smrho qcell: add correction to smrho= -1.27398140925549797E-007 -1.62210738354781012E-009
 mixrho: all smrho is positive for isp=           1

 iors  : write restart file (binary, mesh density) 

   it  4  of 20    ehf=   -3304.827179   ehk=   -3304.518630
 From last iter    ehf=   -3304.888087   ehk=   -3304.422861
 diffe(q)=  0.060908 (0.032470)    tol= 0.000010 (0.000010)   more=T
i ehf=-3304.8271788 ehk=-3304.5186303

 --- BNDFP:  begin iteration 5 of 20 ---
 ttt nevmx w=           0  5.00000000000000010E-003

 avg es pot at rmt= 0.583885  avg sphere pot= 0.638832  vconst=-0.583885
 smvxcm: all smrho_w is positive
  i job kmax lfltwf(FRZWF see locpot.F)=           0           1           5 T

 site  1  z= 29.0  rmt= 2.28000  nr=655   a=0.015  nlml=36  rg=0.570  Vfloat=T
 sm core charge = 0.295521 (sphere) + 0.00535 (spillout) = 0.30087
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===
 potential shift to crystal energy zero:    0.000102


 subzi: tetrahedron integration of bands; tetrahedron integration of density

 Start first of two band passes ...
 end of suham2
 -------- qplist --------
    1   0.000   0.000   0.000
    2  -0.125   0.125   0.125
    3  -0.250   0.250   0.250
    4  -0.375   0.375   0.375
    5  -0.500   0.500   0.500
    6   0.000   0.000   0.250
    7  -0.125   0.125   0.375
    8  -0.250   0.250   0.500
    9  -0.375   0.375   0.625
   10  -0.500   0.500   0.750
   11  -0.625   0.625   0.875
   12  -0.750   0.750   1.000
   13   0.000   0.000   0.500
   14  -0.125   0.125   0.625
   15  -0.250   0.250   0.750
   16  -0.375   0.375   0.875
   17  -0.500   0.500   1.000
   18   0.000   0.000   0.750
   19  -0.125   0.125   0.875
   20  -0.250   0.250   1.000
   21   0.000   0.000   1.000
   22   0.000   0.250   0.500
   23  -0.125   0.375   0.625
   24  -0.250   0.500   0.750
   25   0.000   0.250   0.750
   26  -0.125   0.375   0.875
   27  -0.250   0.500   1.000
   28   0.000   0.250   1.000
   29   0.000   0.500   1.000
 sigmamode= F
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.94D-06    2  0.51D-03    3  0.51D-03    4  0.51D-03    5  0.35D-01
    6  0.35D-01    7  0.45D-01    8  0.78D-01    9  0.78D-01   10  0.78D-01
   11  0.36D+00   12  0.36D+00   13  0.36D+00   14  0.42D+00   15  0.53D+00
   ... skip larger eigenvalues ...
 eigenvalue=
 -0.7168 -0.2452 -0.2452 -0.2452 -0.1813 -0.1813  1.6704  1.8607  1.8607
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.4805 -0.2747 -0.2412 -0.2232 -0.1862 -0.1319  0.8099  1.2287  1.7474
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.3845 -0.3509 -0.1392 -0.1278 -0.1278  0.0884  0.5235  0.9479  0.9479

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.009792;  11.000000 electrons
         Sum occ. bands:   -2.696693, incl. Bloechl correction: -0.009350

 Saved qp weights ...
 Start second band pass ...
 -------- qplist --------
    1   0.000   0.000   0.000
    2  -0.125   0.125   0.125
    3  -0.250   0.250   0.250
    4  -0.375   0.375   0.375
    5  -0.500   0.500   0.500
    6   0.000   0.000   0.250
    7  -0.125   0.125   0.375
    8  -0.250   0.250   0.500
    9  -0.375   0.375   0.625
   10  -0.500   0.500   0.750
   11  -0.625   0.625   0.875
   12  -0.750   0.750   1.000
   13   0.000   0.000   0.500
   14  -0.125   0.125   0.625
   15  -0.250   0.250   0.750
   16  -0.375   0.375   0.875
   17  -0.500   0.500   1.000
   18   0.000   0.000   0.750
   19  -0.125   0.125   0.875
   20  -0.250   0.250   1.000
   21   0.000   0.000   1.000
   22   0.000   0.250   0.500
   23  -0.125   0.375   0.625
   24  -0.250   0.500   0.750
   25   0.000   0.250   0.750
   26  -0.125   0.375   0.875
   27  -0.250   0.500   1.000
   28   0.000   0.250   1.000
   29   0.000   0.500   1.000
 sigmamode= F
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.94D-06    2  0.51D-03    3  0.51D-03    4  0.51D-03    5  0.35D-01
    6  0.35D-01    7  0.45D-01    8  0.78D-01    9  0.78D-01   10  0.78D-01
   11  0.36D+00   12  0.36D+00   13  0.36D+00   14  0.42D+00   15  0.53D+00
   ... skip larger eigenvalues ...
 eigenvalue=
 -0.7168 -0.2452 -0.2452 -0.2452 -0.1813 -0.1813  1.6704  1.8607  1.8607
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.4805 -0.2747 -0.2412 -0.2232 -0.1862 -0.1319  0.8099  1.2287  1.7474
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.3845 -0.3509 -0.1392 -0.1278 -0.1278  0.0884  0.5235  0.9479  0.9479

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.009792;  11.000000 electrons
         Sum occ. bands:   -2.696693, incl. Bloechl correction: -0.009350

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local
   1   10.054917    2.908407    7.146510

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:Cu      
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.458594   -0.417172    4.656541    4.654884    4.500000    4.654884
 1     0    0.391469   -0.277780    4.396391    4.380823    4.250000    4.380823
 2     0    9.182226   -0.226726    3.868334    3.878624    3.147584    3.878624
 3     0    0.017608   -0.233836    4.125224    4.121347    4.102416    4.121347
 4     0    0.003938   -0.231284    5.086988    5.085324    5.077979    5.085324
 5     0    0.001082   -0.215596    6.067445    6.066567    6.062833    6.066567

 Harris energy:
 sumev=       -2.696693  val*vef=    -190.146711   sumtv=     187.450018
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756520
 rhoeps=    -129.931782     utot=   -6534.038420    ehar=   -3304.763664

 avg es pot at rmt= 0.587252  avg sphere pot= 0.639004  vconst=-0.587252
 smvxcm: all smrho_w is positive
  i job kmax lfltwf(FRZWF see locpot.F)=           0           0           5 F

 site  1  z= 29.0  rmt= 2.28000  nr=655   a=0.015  nlml=36  rg=0.570  Vfloat=F
 sm core charge = 0.295521 (sphere) + 0.00535 (spillout) = 0.30087
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===


 ekin=3358.698119  rho*v=-6663.461505 ehf =-3304.763664  ehks =-3304.763386
 mixrho: sum smrho  init = 0.164795D+03 0.525431D-27 0.164795D+03       0
 mixrho: sum smrnew new  = 0.165594D+03-0.165528D-16 0.165594D+03       0
  
 mixrho: dqsum rmsuns=  0.23659D-03  0.24293D-03 -0.27504D-19
 mixrealsmooth= T
 wgtsmooth=  1.72132593164774084E-002
 mixrho:  sought 8 iter from file mixm; read 4.  RMS DQ=1.00e-3  last it=3.25e-2
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.834908      3.853490      3.853490      0.000243      3.853490
 site    1    7.165092      7.146510      7.146510      0.000851      7.146510
 AMIX: nmix=4 mmix=8  nelts=5333  beta=1  tm=5  rmsdel=1e-3
   tj: 0.39679   0.67043  -0.14970   0.00002
 mixrealsmooth= T
 smrho qcell: add correction to smrho= -1.13463975104366455E-007 -1.44468946286300318E-009
 mixrho: all smrho is positive for isp=           1

 iors  : write restart file (binary, mesh density) 

   it  5  of 20    ehf=   -3304.763664   ehk=   -3304.763386
 From last iter    ehf=   -3304.827179   ehk=   -3304.518630
 diffe(q)=  0.063515 (0.001004)    tol= 0.000010 (0.000010)   more=T
i ehf=-3304.7636643 ehk=-3304.763386

 --- BNDFP:  begin iteration 6 of 20 ---
 ttt nevmx w=           0  5.00000000000000010E-003

 avg es pot at rmt= 0.583943  avg sphere pot= 0.638853  vconst=-0.583943
 smvxcm: all smrho_w is positive
  i job kmax lfltwf(FRZWF see locpot.F)=           0           1           5 T

 site  1  z= 29.0  rmt= 2.28000  nr=655   a=0.015  nlml=36  rg=0.570  Vfloat=T
 sm core charge = 0.295521 (sphere) + 0.00535 (spillout) = 0.30087
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===
 potential shift to crystal energy zero:    0.000102


 subzi: tetrahedron integration of bands; tetrahedron integration of density

 Start first of two band passes ...
 end of suham2
 -------- qplist --------
    1   0.000   0.000   0.000
    2  -0.125   0.125   0.125
    3  -0.250   0.250   0.250
    4  -0.375   0.375   0.375
    5  -0.500   0.500   0.500
    6   0.000   0.000   0.250
    7  -0.125   0.125   0.375
    8  -0.250   0.250   0.500
    9  -0.375   0.375   0.625
   10  -0.500   0.500   0.750
   11  -0.625   0.625   0.875
   12  -0.750   0.750   1.000
   13   0.000   0.000   0.500
   14  -0.125   0.125   0.625
   15  -0.250   0.250   0.750
   16  -0.375   0.375   0.875
   17  -0.500   0.500   1.000
   18   0.000   0.000   0.750
   19  -0.125   0.125   0.875
   20  -0.250   0.250   1.000
   21   0.000   0.000   1.000
   22   0.000   0.250   0.500
   23  -0.125   0.375   0.625
   24  -0.250   0.500   0.750
   25   0.000   0.250   0.750
   26  -0.125   0.375   0.875
   27  -0.250   0.500   1.000
   28   0.000   0.250   1.000
   29   0.000   0.500   1.000
 sigmamode= F
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.94D-06    2  0.51D-03    3  0.51D-03    4  0.51D-03    5  0.35D-01
    6  0.35D-01    7  0.45D-01    8  0.79D-01    9  0.79D-01   10  0.79D-01
   11  0.36D+00   12  0.36D+00   13  0.36D+00   14  0.42D+00   15  0.53D+00
   ... skip larger eigenvalues ...
 eigenvalue=
 -0.7169 -0.2456 -0.2456 -0.2456 -0.1819 -0.1819  1.6709  1.8614  1.8614
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.4806 -0.2752 -0.2417 -0.2237 -0.1868 -0.1324  0.8092  1.2279  1.7478
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.3849 -0.3513 -0.1398 -0.1283 -0.1283  0.0883  0.5227  0.9481  0.9481

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.010171;  11.000000 electrons
         Sum occ. bands:   -2.701581, incl. Bloechl correction: -0.009352

 Saved qp weights ...
 Start second band pass ...
 -------- qplist --------
    1   0.000   0.000   0.000
    2  -0.125   0.125   0.125
    3  -0.250   0.250   0.250
    4  -0.375   0.375   0.375
    5  -0.500   0.500   0.500
    6   0.000   0.000   0.250
    7  -0.125   0.125   0.375
    8  -0.250   0.250   0.500
    9  -0.375   0.375   0.625
   10  -0.500   0.500   0.750
   11  -0.625   0.625   0.875
   12  -0.750   0.750   1.000
   13   0.000   0.000   0.500
   14  -0.125   0.125   0.625
   15  -0.250   0.250   0.750
   16  -0.375   0.375   0.875
   17  -0.500   0.500   1.000
   18   0.000   0.000   0.750
   19  -0.125   0.125   0.875
   20  -0.250   0.250   1.000
   21   0.000   0.000   1.000
   22   0.000   0.250   0.500
   23  -0.125   0.375   0.625
   24  -0.250   0.500   0.750
   25   0.000   0.250   0.750
   26  -0.125   0.375   0.875
   27  -0.250   0.500   1.000
   28   0.000   0.250   1.000
   29   0.000   0.500   1.000
 sigmamode= F
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.94D-06    2  0.51D-03    3  0.51D-03    4  0.51D-03    5  0.35D-01
    6  0.35D-01    7  0.45D-01    8  0.79D-01    9  0.79D-01   10  0.79D-01
   11  0.36D+00   12  0.36D+00   13  0.36D+00   14  0.42D+00   15  0.53D+00
   ... skip larger eigenvalues ...
 eigenvalue=
 -0.7169 -0.2456 -0.2456 -0.2456 -0.1819 -0.1819  1.6709  1.8614  1.8614
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.4806 -0.2752 -0.2417 -0.2237 -0.1868 -0.1324  0.8092  1.2279  1.7478
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.3849 -0.3513 -0.1398 -0.1283 -0.1283  0.0883  0.5227  0.9481  0.9481

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.010171;  11.000000 electrons
         Sum occ. bands:   -2.701581, incl. Bloechl correction: -0.009352

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local
   1   10.055555    2.907931    7.147624

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:Cu      
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.458474   -0.417323    4.654884    4.654880    4.500000    4.654880
 1     0    0.391168   -0.278035    4.380823    4.380782    4.250000    4.380782
 2     0    9.183296   -0.227218    3.878624    3.878650    3.147584    3.878650
 3     0    0.017599   -0.234264    4.121347    4.121336    4.102416    4.121336
 4     0    0.003937   -0.231709    5.085324    5.085320    5.077979    5.085320
 5     0    0.001082   -0.216033    6.066567    6.066565    6.062833    6.066565

 Harris energy:
 sumev=       -2.701581  val*vef=    -190.144398   sumtv=     187.442817
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756520
 rhoeps=    -129.930789     utot=   -6534.032221    ehar=   -3304.763673

 avg es pot at rmt= 0.587023  avg sphere pot= 0.638966  vconst=-0.587023
 smvxcm: all smrho_w is positive
  i job kmax lfltwf(FRZWF see locpot.F)=           0           0           5 F

 site  1  z= 29.0  rmt= 2.28000  nr=655   a=0.015  nlml=36  rg=0.570  Vfloat=F
 sm core charge = 0.295521 (sphere) + 0.00535 (spillout) = 0.30087
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===


 ekin=3358.741597  rho*v=-6663.505045 ehf =-3304.763673  ehks =-3304.763449
 mixrho: sum smrho  init = 0.164806D+03-0.121425D-26 0.164806D+03       0
 mixrho: sum smrnew new  = 0.165546D+03-0.170383D-16 0.165546D+03       0
  
 mixrho: dqsum rmsuns=  0.21927D-03  0.22412D-03  0.16343D-19
 mixrealsmooth= T
 wgtsmooth=  1.72132593164774084E-002
 mixrho:  sought 8 iter from file mixm; read 5.  RMS DQ=9.06e-4  last it=1.00e-3
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.835155      3.852376      3.852376      0.000224      3.852376
 site    1    7.164845      7.147624      7.147624      0.000770      7.147624
 AMIX: condition of normal eqns >100000. Reducing nmix to 4
 AMIX: condition of normal eqns >100000. Reducing nmix to 3
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: nmix=2 mmix=8  nelts=5333  beta=1  tm=5  rmsdel=9.06e-4
   tj:-4.95700  -0.01351
 mixrealsmooth= T
 smrho qcell: add correction to smrho= -2.25884172699863939E-007 -2.87608894212328699E-009
 mixrho: all smrho is positive for isp=           1

 iors  : write restart file (binary, mesh density) 

   it  6  of 20    ehf=   -3304.763673   ehk=   -3304.763449
 From last iter    ehf=   -3304.763664   ehk=   -3304.763386
 diffe(q)= -0.000009 (0.000906)    tol= 0.000010 (0.000010)   more=T
i ehf=-3304.7636733 ehk=-3304.7634487

 --- BNDFP:  begin iteration 7 of 20 ---
 ttt nevmx w=           0  5.00000000000000010E-003

 avg es pot at rmt= 0.584578  avg sphere pot= 0.638845  vconst=-0.584578
 smvxcm: all smrho_w is positive
  i job kmax lfltwf(FRZWF see locpot.F)=           0           1           5 T

 site  1  z= 29.0  rmt= 2.28000  nr=655   a=0.015  nlml=36  rg=0.570  Vfloat=T
 sm core charge = 0.295521 (sphere) + 0.00535 (spillout) = 0.30087
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===
 potential shift to crystal energy zero:    0.000102


 subzi: tetrahedron integration of bands; tetrahedron integration of density

 Start first of two band passes ...
 end of suham2
 -------- qplist --------
    1   0.000   0.000   0.000
    2  -0.125   0.125   0.125
    3  -0.250   0.250   0.250
    4  -0.375   0.375   0.375
    5  -0.500   0.500   0.500
    6   0.000   0.000   0.250
    7  -0.125   0.125   0.375
    8  -0.250   0.250   0.500
    9  -0.375   0.375   0.625
   10  -0.500   0.500   0.750
   11  -0.625   0.625   0.875
   12  -0.750   0.750   1.000
   13   0.000   0.000   0.500
   14  -0.125   0.125   0.625
   15  -0.250   0.250   0.750
   16  -0.375   0.375   0.875
   17  -0.500   0.500   1.000
   18   0.000   0.000   0.750
   19  -0.125   0.125   0.875
   20  -0.250   0.250   1.000
   21   0.000   0.000   1.000
   22   0.000   0.250   0.500
   23  -0.125   0.375   0.625
   24  -0.250   0.500   0.750
   25   0.000   0.250   0.750
   26  -0.125   0.375   0.875
   27  -0.250   0.500   1.000
   28   0.000   0.250   1.000
   29   0.000   0.500   1.000
 sigmamode= F
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.94D-06    2  0.51D-03    3  0.51D-03    4  0.51D-03    5  0.35D-01
    6  0.35D-01    7  0.45D-01    8  0.79D-01    9  0.79D-01   10  0.79D-01
   11  0.36D+00   12  0.36D+00   13  0.36D+00   14  0.42D+00   15  0.53D+00
   ... skip larger eigenvalues ...
 eigenvalue=
 -0.7176 -0.2496 -0.2496 -0.2496 -0.1863 -0.1863  1.6706  1.8606  1.8606
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.4816 -0.2787 -0.2457 -0.2278 -0.1910 -0.1367  0.8078  1.2268  1.7472
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.3880 -0.3548 -0.1443 -0.1327 -0.1327  0.0876  0.5207  0.9472  0.9472

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.012904;  11.000000 electrons
         Sum occ. bands:   -2.741362, incl. Bloechl correction: -0.009395

 Saved qp weights ...
 Start second band pass ...
 -------- qplist --------
    1   0.000   0.000   0.000
    2  -0.125   0.125   0.125
    3  -0.250   0.250   0.250
    4  -0.375   0.375   0.375
    5  -0.500   0.500   0.500
    6   0.000   0.000   0.250
    7  -0.125   0.125   0.375
    8  -0.250   0.250   0.500
    9  -0.375   0.375   0.625
   10  -0.500   0.500   0.750
   11  -0.625   0.625   0.875
   12  -0.750   0.750   1.000
   13   0.000   0.000   0.500
   14  -0.125   0.125   0.625
   15  -0.250   0.250   0.750
   16  -0.375   0.375   0.875
   17  -0.500   0.500   1.000
   18   0.000   0.000   0.750
   19  -0.125   0.125   0.875
   20  -0.250   0.250   1.000
   21   0.000   0.000   1.000
   22   0.000   0.250   0.500
   23  -0.125   0.375   0.625
   24  -0.250   0.500   0.750
   25   0.000   0.250   0.750
   26  -0.125   0.375   0.875
   27  -0.250   0.500   1.000
   28   0.000   0.250   1.000
   29   0.000   0.500   1.000
 sigmamode= F
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.94D-06    2  0.51D-03    3  0.51D-03    4  0.51D-03    5  0.35D-01
    6  0.35D-01    7  0.45D-01    8  0.79D-01    9  0.79D-01   10  0.79D-01
   11  0.36D+00   12  0.36D+00   13  0.36D+00   14  0.42D+00   15  0.53D+00
   ... skip larger eigenvalues ...
 eigenvalue=
 -0.7176 -0.2496 -0.2496 -0.2496 -0.1863 -0.1863  1.6706  1.8606  1.8606
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.4816 -0.2787 -0.2457 -0.2278 -0.1910 -0.1367  0.8078  1.2268  1.7472
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.3880 -0.3548 -0.1443 -0.1327 -0.1327  0.0876  0.5207  0.9472  0.9472

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.012904;  11.000000 electrons
         Sum occ. bands:   -2.741362, incl. Bloechl correction: -0.009395

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local
   1   10.059917    2.899949    7.159968

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:Cu      
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.457782   -0.418517    4.654880    4.654867    4.500000    4.654867
 1     0    0.388860   -0.280075    4.380782    4.380472    4.250000    4.380472
 2     0    9.190824   -0.231214    3.878650    3.878892    3.147584    3.878892
 3     0    0.017461   -0.237952    4.121336    4.121245    4.102416    4.121245
 4     0    0.003914   -0.235658    5.085320    5.085281    5.077979    5.085281
 5     0    0.001077   -0.220088    6.066565    6.066544    6.062833    6.066544

 Harris energy:
 sumev=       -2.741362  val*vef=    -190.126130   sumtv=     187.384768
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756520
 rhoeps=    -129.922720     utot=   -6533.982258    ehar=   -3304.763692

 avg es pot at rmt= 0.584880  avg sphere pot= 0.639044  vconst=-0.584880
 smvxcm: all smrho_w is positive
  i job kmax lfltwf(FRZWF see locpot.F)=           0           0           5 F

 site  1  z= 29.0  rmt= 2.28000  nr=655   a=0.015  nlml=36  rg=0.570  Vfloat=F
 sm core charge = 0.295521 (sphere) + 0.00535 (spillout) = 0.30087
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===


 ekin=3359.063089  rho*v=-6663.826775 ehf =-3304.763692  ehks =-3304.763686
 mixrho: sum smrho  init = 0.164988D+03-0.918037D-28 0.164988D+03       0
 mixrho: sum smrnew new  = 0.165016D+03 0.140913D-16 0.165016D+03       0
  
 mixrho: dqsum rmsuns=  0.82062D-05  0.48422D-04 -0.51378D-19
 mixrealsmooth= T
 wgtsmooth=  1.72132593164774084E-002
 mixrho:  sought 8 iter from file mixm; read 6.  RMS DQ=1.47e-4  last it=9.06e-4
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.839387      3.840032      3.840032      0.000048      3.840032
 site    1    7.160613      7.159968      7.159968      0.000126      7.159968
 AMIX: condition of normal eqns >100000. Reducing nmix to 5
 AMIX: condition of normal eqns >100000. Reducing nmix to 4
 AMIX: condition of normal eqns >100000. Reducing nmix to 3
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=5333  beta=1  tm=5  rmsdel=1.47e-4
   tj:-0.17342
 mixrealsmooth= T
 smrho qcell: add correction to smrho= -3.15239923054377869E-009 -4.01381843612913807E-011
 mixrho: all smrho is positive for isp=           1

 iors  : write restart file (binary, mesh density) 

   it  7  of 20    ehf=   -3304.763692   ehk=   -3304.763686
 From last iter    ehf=   -3304.763673   ehk=   -3304.763449
 diffe(q)= -0.000018 (0.000147)    tol= 0.000010 (0.000010)   more=T
i ehf=-3304.7636917 ehk=-3304.7636859

 --- BNDFP:  begin iteration 8 of 20 ---
 ttt nevmx w=           0  5.00000000000000010E-003

 avg es pot at rmt= 0.584508  avg sphere pot= 0.639058  vconst=-0.584508
 smvxcm: all smrho_w is positive
  i job kmax lfltwf(FRZWF see locpot.F)=           0           1           5 T

 site  1  z= 29.0  rmt= 2.28000  nr=655   a=0.015  nlml=36  rg=0.570  Vfloat=T
 sm core charge = 0.295521 (sphere) + 0.00535 (spillout) = 0.30087
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===
 potential shift to crystal energy zero:    0.000102


 subzi: tetrahedron integration of bands; tetrahedron integration of density

 Start first of two band passes ...
 end of suham2
 -------- qplist --------
    1   0.000   0.000   0.000
    2  -0.125   0.125   0.125
    3  -0.250   0.250   0.250
    4  -0.375   0.375   0.375
    5  -0.500   0.500   0.500
    6   0.000   0.000   0.250
    7  -0.125   0.125   0.375
    8  -0.250   0.250   0.500
    9  -0.375   0.375   0.625
   10  -0.500   0.500   0.750
   11  -0.625   0.625   0.875
   12  -0.750   0.750   1.000
   13   0.000   0.000   0.500
   14  -0.125   0.125   0.625
   15  -0.250   0.250   0.750
   16  -0.375   0.375   0.875
   17  -0.500   0.500   1.000
   18   0.000   0.000   0.750
   19  -0.125   0.125   0.875
   20  -0.250   0.250   1.000
   21   0.000   0.000   1.000
   22   0.000   0.250   0.500
   23  -0.125   0.375   0.625
   24  -0.250   0.500   0.750
   25   0.000   0.250   0.750
   26  -0.125   0.375   0.875
   27  -0.250   0.500   1.000
   28   0.000   0.250   1.000
   29   0.000   0.500   1.000
 sigmamode= F
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.94D-06    2  0.51D-03    3  0.51D-03    4  0.51D-03    5  0.35D-01
    6  0.35D-01    7  0.45D-01    8  0.79D-01    9  0.79D-01   10  0.79D-01
   11  0.36D+00   12  0.36D+00   13  0.36D+00   14  0.42D+00   15  0.53D+00
   ... skip larger eigenvalues ...
 eigenvalue=
 -0.7178 -0.2506 -0.2506 -0.2506 -0.1874 -0.1874  1.6705  1.8604  1.8604
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.4819 -0.2796 -0.2467 -0.2289 -0.1921 -0.1377  0.8074  1.2266  1.7471
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.3887 -0.3557 -0.1454 -0.1339 -0.1339  0.0874  0.5201  0.9470  0.9470

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.013582;  11.000000 electrons
         Sum occ. bands:   -2.751237, incl. Bloechl correction: -0.009405

 Saved qp weights ...
 Start second band pass ...
 -------- qplist --------
    1   0.000   0.000   0.000
    2  -0.125   0.125   0.125
    3  -0.250   0.250   0.250
    4  -0.375   0.375   0.375
    5  -0.500   0.500   0.500
    6   0.000   0.000   0.250
    7  -0.125   0.125   0.375
    8  -0.250   0.250   0.500
    9  -0.375   0.375   0.625
   10  -0.500   0.500   0.750
   11  -0.625   0.625   0.875
   12  -0.750   0.750   1.000
   13   0.000   0.000   0.500
   14  -0.125   0.125   0.625
   15  -0.250   0.250   0.750
   16  -0.375   0.375   0.875
   17  -0.500   0.500   1.000
   18   0.000   0.000   0.750
   19  -0.125   0.125   0.875
   20  -0.250   0.250   1.000
   21   0.000   0.000   1.000
   22   0.000   0.250   0.500
   23  -0.125   0.375   0.625
   24  -0.250   0.500   0.750
   25   0.000   0.250   0.750
   26  -0.125   0.375   0.875
   27  -0.250   0.500   1.000
   28   0.000   0.250   1.000
   29   0.000   0.500   1.000
 sigmamode= F
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.94D-06    2  0.51D-03    3  0.51D-03    4  0.51D-03    5  0.35D-01
    6  0.35D-01    7  0.45D-01    8  0.79D-01    9  0.79D-01   10  0.79D-01
   11  0.36D+00   12  0.36D+00   13  0.36D+00   14  0.42D+00   15  0.53D+00
   ... skip larger eigenvalues ...
 eigenvalue=
 -0.7178 -0.2506 -0.2506 -0.2506 -0.1874 -0.1874  1.6705  1.8604  1.8604
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.4819 -0.2796 -0.2467 -0.2289 -0.1921 -0.1377  0.8074  1.2266  1.7471
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.3887 -0.3557 -0.1454 -0.1339 -0.1339  0.0874  0.5201  0.9470  0.9470

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.013582;  11.000000 electrons
         Sum occ. bands:   -2.751237, incl. Bloechl correction: -0.009405

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local
   1   10.060997    2.897854    7.163142

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:Cu      
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.457600   -0.418812    4.654867    4.654848    4.500000    4.654848
 1     0    0.388282   -0.280583    4.380472    4.380380    4.250000    4.380380
 2     0    9.192703   -0.232205    3.878892    3.878950    3.147584    3.878950
 3     0    0.017427   -0.238868    4.121245    4.121220    4.102416    4.121220
 4     0    0.003908   -0.236641    5.085281    5.085270    5.077979    5.085270
 5     0    0.001076   -0.221096    6.066544    6.066538    6.062833    6.066538

 Harris energy:
 sumev=       -2.751237  val*vef=    -190.113560   sumtv=     187.362323
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756520
 rhoeps=    -129.920581     utot=   -6533.961971    ehar=   -3304.763709

 avg es pot at rmt= 0.584338  avg sphere pot= 0.639065  vconst=-0.584338
 smvxcm: all smrho_w is positive
  i job kmax lfltwf(FRZWF see locpot.F)=           0           0           5 F

 site  1  z= 29.0  rmt= 2.28000  nr=655   a=0.015  nlml=36  rg=0.570  Vfloat=F
 sm core charge = 0.295521 (sphere) + 0.00535 (spillout) = 0.30087
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===


 ekin=3359.146492  rho*v=-6663.910200 ehf =-3304.763709  ehks =-3304.763708
 mixrho: sum smrho  init = 0.164924D+03 0.170591D-27 0.164924D+03       0
 mixrho: sum smrnew new  = 0.164879D+03 0.756802D-17 0.164879D+03       0
  
 mixrho: dqsum rmsuns= -0.13157D-04  0.13410D-04  0.10913D-19
 mixrealsmooth= T
 wgtsmooth=  1.72132593164774084E-002
 mixrho:  sought 8 iter from file mixm; read 7.  RMS DQ=5.08e-5  last it=1.47e-4
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.837891      3.836858      3.836858      0.000013      3.836858
 site    1    7.162109      7.163142      7.163142      0.000040      7.163142
 AMIX: condition of normal eqns >100000. Reducing nmix to 6
 AMIX: condition of normal eqns >100000. Reducing nmix to 5
 AMIX: condition of normal eqns >100000. Reducing nmix to 4
 AMIX: condition of normal eqns >100000. Reducing nmix to 3
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: nmix=2 mmix=8  nelts=5333  beta=1  tm=5  rmsdel=5.08e-5
   tj: 0.03869   0.04471
 mixrealsmooth= T
 smrho qcell: add correction to smrho= -4.10325728950056146E-009 -5.22450633701545939E-011
 mixrho: all smrho is positive for isp=           1

 iors  : write restart file (binary, mesh density) 

   it  8  of 20    ehf=   -3304.763709   ehk=   -3304.763708
 From last iter    ehf=   -3304.763692   ehk=   -3304.763686
 diffe(q)= -0.000017 (0.000051)    tol= 0.000010 (0.000010)   more=T
i ehf=-3304.7637088 ehk=-3304.7637081

 --- BNDFP:  begin iteration 9 of 20 ---
 ttt nevmx w=           0  5.00000000000000010E-003

 avg es pot at rmt= 0.584479  avg sphere pot= 0.639060  vconst=-0.584479
 smvxcm: all smrho_w is positive
  i job kmax lfltwf(FRZWF see locpot.F)=           0           1           5 T

 site  1  z= 29.0  rmt= 2.28000  nr=655   a=0.015  nlml=36  rg=0.570  Vfloat=T
 sm core charge = 0.295521 (sphere) + 0.00535 (spillout) = 0.30087
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===
 potential shift to crystal energy zero:    0.000102


 subzi: tetrahedron integration of bands; tetrahedron integration of density

 Start first of two band passes ...
 end of suham2
 -------- qplist --------
    1   0.000   0.000   0.000
    2  -0.125   0.125   0.125
    3  -0.250   0.250   0.250
    4  -0.375   0.375   0.375
    5  -0.500   0.500   0.500
    6   0.000   0.000   0.250
    7  -0.125   0.125   0.375
    8  -0.250   0.250   0.500
    9  -0.375   0.375   0.625
   10  -0.500   0.500   0.750
   11  -0.625   0.625   0.875
   12  -0.750   0.750   1.000
   13   0.000   0.000   0.500
   14  -0.125   0.125   0.625
   15  -0.250   0.250   0.750
   16  -0.375   0.375   0.875
   17  -0.500   0.500   1.000
   18   0.000   0.000   0.750
   19  -0.125   0.125   0.875
   20  -0.250   0.250   1.000
   21   0.000   0.000   1.000
   22   0.000   0.250   0.500
   23  -0.125   0.375   0.625
   24  -0.250   0.500   0.750
   25   0.000   0.250   0.750
   26  -0.125   0.375   0.875
   27  -0.250   0.500   1.000
   28   0.000   0.250   1.000
   29   0.000   0.500   1.000
 sigmamode= F
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.94D-06    2  0.51D-03    3  0.51D-03    4  0.51D-03    5  0.35D-01
    6  0.35D-01    7  0.45D-01    8  0.79D-01    9  0.79D-01   10  0.79D-01
   11  0.36D+00   12  0.36D+00   13  0.36D+00   14  0.42D+00   15  0.53D+00
   ... skip larger eigenvalues ...
 eigenvalue=
 -0.7177 -0.2504 -0.2504 -0.2504 -0.1871 -0.1871  1.6705  1.8605  1.8605
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.4818 -0.2794 -0.2465 -0.2287 -0.1919 -0.1375  0.8075  1.2266  1.7471
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.3886 -0.3555 -0.1452 -0.1336 -0.1336  0.0874  0.5202  0.9471  0.9471

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.013446;  11.000000 electrons
         Sum occ. bands:   -2.749249, incl. Bloechl correction: -0.009403

 Saved qp weights ...
 Start second band pass ...
 -------- qplist --------
    1   0.000   0.000   0.000
    2  -0.125   0.125   0.125
    3  -0.250   0.250   0.250
    4  -0.375   0.375   0.375
    5  -0.500   0.500   0.500
    6   0.000   0.000   0.250
    7  -0.125   0.125   0.375
    8  -0.250   0.250   0.500
    9  -0.375   0.375   0.625
   10  -0.500   0.500   0.750
   11  -0.625   0.625   0.875
   12  -0.750   0.750   1.000
   13   0.000   0.000   0.500
   14  -0.125   0.125   0.625
   15  -0.250   0.250   0.750
   16  -0.375   0.375   0.875
   17  -0.500   0.500   1.000
   18   0.000   0.000   0.750
   19  -0.125   0.125   0.875
   20  -0.250   0.250   1.000
   21   0.000   0.000   1.000
   22   0.000   0.250   0.500
   23  -0.125   0.375   0.625
   24  -0.250   0.500   0.750
   25   0.000   0.250   0.750
   26  -0.125   0.375   0.875
   27  -0.250   0.500   1.000
   28   0.000   0.250   1.000
   29   0.000   0.500   1.000
 sigmamode= F
 bndfp:  kpt 1 of 29, k=  0.00000  0.00000  0.00000   ndimh = 24
 zhev_tk: ovlmat=
    1  0.94D-06    2  0.51D-03    3  0.51D-03    4  0.51D-03    5  0.35D-01
    6  0.35D-01    7  0.45D-01    8  0.79D-01    9  0.79D-01   10  0.79D-01
   11  0.36D+00   12  0.36D+00   13  0.36D+00   14  0.42D+00   15  0.53D+00
   ... skip larger eigenvalues ...
 eigenvalue=
 -0.7177 -0.2504 -0.2504 -0.2504 -0.1871 -0.1871  1.6705  1.8605  1.8605
 bndfp:  kpt 11 of 29, k=  0.37500  -0.37500  -0.12500   ndimh = 25
 -0.4818 -0.2794 -0.2465 -0.2287 -0.1919 -0.1375  0.8075  1.2266  1.7471
 bndfp:  kpt 21 of 29, k=  0.00000  0.00000  -1.00000   ndimh = 23
 -0.3886 -0.3555 -0.1452 -0.1336 -0.1336  0.0874  0.5202  0.9471  0.9471

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.013446;  11.000000 electrons
         Sum occ. bands:   -2.749249, incl. Bloechl correction: -0.009403

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local
   1   10.060772    2.898298    7.162474

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:Cu      
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.457637   -0.418756    4.654848    4.654851    4.500000    4.654851
 1     0    0.388401   -0.280483    4.380380    4.380397    4.250000    4.380397
 2     0    9.192313   -0.232005    3.878950    3.878938    3.147584    3.878938
 3     0    0.017434   -0.238683    4.121220    4.121225    4.102416    4.121225
 4     0    0.003910   -0.236443    5.085270    5.085272    5.077979    5.085272
 5     0    0.001076   -0.220892    6.066538    6.066539    6.062833    6.066539

 Harris energy:
 sumev=       -2.749249  val*vef=    -190.117889   sumtv=     187.368641
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756520
 rhoeps=    -129.921067     utot=   -6533.967798    ehar=   -3304.763705

 avg es pot at rmt= 0.584452  avg sphere pot= 0.639062  vconst=-0.584452
 smvxcm: all smrho_w is positive
  i job kmax lfltwf(FRZWF see locpot.F)=           0           0           5 F

 site  1  z= 29.0  rmt= 2.28000  nr=655   a=0.015  nlml=36  rg=0.570  Vfloat=F
 sm core charge = 0.295521 (sphere) + 0.00535 (spillout) = 0.30087
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===


 ekin=3359.128252  rho*v=-6663.891956 ehf =-3304.763705  ehks =-3304.763705
 mixrho: sum smrho  init = 0.164914D+03-0.199779D-27 0.164914D+03       0
 mixrho: sum smrnew new  = 0.164908D+03 0.181775D-17 0.164908D+03       0
  
 mixrho: dqsum rmsuns= -0.18874D-05  0.19484D-05  0.41811D-19
 mixrealsmooth= T
 wgtsmooth=  1.72132593164774084E-002
 mixrho:  sought 8 iter from file mixm; read 8.  RMS DQ=6.36e-6  last it=5.08e-5
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.837674      3.837526      3.837526      0.000002      3.837526
 site    1    7.162326      7.162474      7.162474      0.000006      7.162474
 AMIX: condition of normal eqns >100000. Reducing nmix to 7
 AMIX: condition of normal eqns >100000. Reducing nmix to 6
 AMIX: condition of normal eqns >100000. Reducing nmix to 5
 AMIX: condition of normal eqns >100000. Reducing nmix to 4
 AMIX: condition of normal eqns >100000. Reducing nmix to 3
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=5333  beta=1  tm=5  rmsdel=6.36e-6
   tj:-0.13742
 mixrealsmooth= T
 smrho qcell: add correction to smrho= -4.48256542995295604E-009 -5.70746356919923363E-011
 mixrho: all smrho is positive for isp=           1

 iors  : write restart file (binary, mesh density) 

   it  9  of 20    ehf=   -3304.763705   ehk=   -3304.763705
 From last iter    ehf=   -3304.763709   ehk=   -3304.763708
 diffe(q)=  0.000004 (0.000006)    tol= 0.000010 (0.000010)   more=F
c ehf=-3304.7637048 ehk=-3304.7637047
 Exit 0 LMF 
