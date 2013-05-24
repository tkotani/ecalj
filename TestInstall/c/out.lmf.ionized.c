rdcmd:  lmfa --no-iactiv c -vzbak=1
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
rdcmd:  lmf  --no-iactiv c -vzbak=1
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
 Uniform density added to neutralize background, q=1.000000

 Smooth charge on mesh:            1.661718    moment    1.332841
 Sum of local charges:             1.338282    moments   0.667159
 Total valence charge:             3.000000    moment    2.000000
 Sum of core charges:              2.000000    moment    0.000000
 Sum of nuclear charges:          -6.000000
 Homogeneous background:           1.000000
 Deviation from neutrality:       -0.000000

 --- BNDFP:  begin iteration 1 of 10 ---
 ttt nevmx w=           5  4.00000000000000008E-003

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
         rhomu =   -1.806539 (  -1.431678,  -0.374861)
       avg vxc =   -0.082065 (  -0.095983,  -0.068146)
 smvxcm: negative smrho_w number,min(smrho_w)=      201254 -4.99895726979550110E-004
 smvxcm: enforce positive smrho_w. Add srshift=  4.99895726989550140E-004
 smooth rhoeps =   -1.494861 (  -1.109298,  -0.385563)
         rhomu =   -1.946499 (  -1.522749,  -0.423750)
       avg vxc =   -0.191138 (  -0.218032,  -0.164244)

 locpot:
  i job kmax lfltwf=           0           1           3 T

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 === rho1 valence true density ===
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       128
 vxcnsp (warning): negative rho: min val =  -4.61E-04
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

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
 core chg:       2.000000       2.000000      -0.000000
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
   core charge             2.000000        -0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

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
 -1.0534 -0.6779 -0.4297 -0.4297 -0.4297  0.2977  0.2977  0.2977
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1365 -0.8248 -0.2365 -0.2365 -0.2365  0.2428  0.2428  0.2428
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
 ... only filled or empty bands encountered:  ev=-0.824805  ec=-0.772717
 VBmax = -0.824805  CBmin = -0.772717  gap = 0.052087 Ry = 0.70839 eV
 BZINTS: Fermi energy:     -0.824805;   3.000000 electrons
         Sum occ. bands:   -3.277726, incl. Bloechl correction:  0.000000
         Mag. moment:      -1.000000

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
 -1.0534 -0.6779 -0.4297 -0.4297 -0.4297  0.2977  0.2977  0.2977
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1365 -0.8248 -0.2365 -0.2365 -0.2365  0.2428  0.2428  0.2428
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
 ... only filled or empty bands encountered:  ev=-0.824805  ec=-0.772717
 VBmax = -0.824805  CBmin = -0.772717  gap = 0.052087 Ry = 0.70839 eV
 BZINTS: Fermi energy:     -0.824805;   3.000000 electrons
         Sum occ. bands:   -3.277726, incl. Bloechl correction:  0.000000
         Mag. moment:      -1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.130156  386.510555 -384.380398     -0.144414 -327.846613  327.702199
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -19.85491  sum tc=    31.38764  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -19.78517  sum tc=    31.54593  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.992767   -0.956642    2.900000    2.967335    2.500000    2.967335
 spn 2 0    1.134455   -0.773670    2.900000    2.946917    2.500000    2.946917
 1     1    0.000104   -0.505501    2.850000    2.802066    2.250000    2.850000
 spn 2 1    0.002815   -0.125968    2.850000    2.943095    2.250000    2.850000
 2     0    0.000000   -1.059853    3.180000    3.117538    3.147584    3.147584
 spn 2 0    0.000012   -1.006449    3.180000    3.114300    3.147584    3.147584
 3     0    0.000000   -0.979574    4.120000    4.090136    4.102416    4.102416
 spn 2 0    0.000004   -0.926739    4.120000    4.088649    4.102416    4.102416

 Harris energy:
 sumev=       -3.277726  val*vef=     -14.137456   sumtv=      10.859730
 sumec=      -39.640082  cor*vef=    -102.572087   ttcor=      62.932005
 rhoeps=      -9.595549     utot=    -139.945263    ehar=     -75.749078
 smvxcm: negative smrho_w number,min(smrho_w)=      201254 -4.99898745262577810E-004
 smvxcm: enforce positive smrho_w. Add srshift=  4.99898745272577840E-004
 smvxcm: negative smrho_w number,min(smrho_w)=      201254 -4.99892708696522410E-004
 smvxcm: enforce positive smrho_w. Add srshift=  4.99892708706522440E-004
 smvxcm: negative smrho_w number,min(smrho_w)=      201254 -4.99895726979550110E-004
 smvxcm: enforce positive smrho_w. Add srshift=  4.99895726989550140E-004

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00    -0.00   -0.00    0.00     0.00    0.00   -0.00
 shift forces to make zero average correction:            0.00    0.00   -0.00

 srhov:  -1237.408080   1225.635192    -11.772888 sumev=   -3.277726   sumtv=    8.495162

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1 -108.431708 -384.380398     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.205342  avg sphere pot= 0.030809  vconst= 0.205342
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves     38.909181   charge   388.380398
 smooth rhoeps =-2470.310652 (-179.107596,***********)
         rhomu =-3278.859204 (-116.221776,***********)
       avg vxc =   -0.316237 (  -0.227596,  -0.404877)
 smvxcm: all smrho_w is positive
 smooth rhoeps =-2470.310652 (-179.107596,***********)
         rhomu =-3278.859204 (-116.221776,***********)
       avg vxc =   -0.316237 (  -0.227596,  -0.404877)

 locpot:
  i job kmax lfltwf=           0           0           3 F

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -5.572123   -5.624673  -11.196796    -82.114225-2616.757709-2698.871934

 local terms:     true           smooth         local
 rhoeps:        -7.980819   -2470.064620    2462.083801
 rhomu:         -5.216693    -116.220424     111.003732
 spin2:         -5.294318   -3162.319051    3157.024733
 total:        -10.511011   -3278.539476    3268.028465
 val*vef       -11.196796   -3195.937100    3184.740304
 val chg:        2.130156     386.510555    -384.380398
 val mom:       -0.144414    -327.846613     327.702199    core:   0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef          -2698.892656      2687.560606       -11.332050
   rhoval*ves             82.992192       -91.768732        -8.776540
   psnuc*ves              -5.173830      -258.804129      -263.977959
   utot                   38.909181      -175.286431      -136.377249
   rho*exc             -2470.310652      2462.083801        -8.226851
   rho*vxc             -3278.859204      3268.028465       -10.830739
   valence chg           387.380398      -384.380398         3.000000
   valence mag          -328.702199       327.702199        -1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=        8.495162  sumtc=        62.933572   ekin=       71.428734
 rhoep=       -8.226851   utot=      -136.377249   ehks=      -73.175366
 mag. mom=    -1.000000  (bands)       -1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
 mixrho: sum smrho  init = 0.187160D+03 0.111757D-25 0.262834D+03   99349
 mixrho: sum smrnew new  = 0.366739D+04 0.200179D-15 0.366739D+04       0
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho: dqsum rmsuns=  0.38572D+00  0.45305D+01 -0.14805D-16
 mixrealsmooth= T
 wgtsmooth=  2.82842712474619005E-003
 mixrho:  sought 2 iter from file mixm; read 0.  RMS DQ=3.66e0
 mixrho: (warning) scr. and lin-mixed densities had 0 and 72123 negative points
 AMIX: nmix=0 mmix=8  nelts=252052  beta=0.5  tm=5  rmsdel=1.83e0
 mixrealsmooth= T
 smrho qcell: add correction to smrho=  4.13685654621076537E-008  2.06842827310538307E-011
 unscreened rms difference:  smooth  6.407083   local 11.944385
   screened rms difference:  smooth  4.853815   local 11.944385   tot  3.656525
 mixrho: warning. negative smrho; isp number min=           1       90779 -6.77138186736459175E-004
 mixrho: warning. negative smrho; isp number min=           2       71667 -6.70355045056268586E-004

 iors  : write restart file (binary, mesh density) 

   it  1  of 10    ehf=      -0.754178   ehk=       1.819534
h zbak=1 mmom=-1.0000001 ehf=-.7541783 ehk=1.8195343

 --- BNDFP:  begin iteration 2 of 10 ---
 ttt nevmx w=           5  4.00000000000000008E-003

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1  -54.027093 -191.521058     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.033200  avg sphere pot= 0.025175  vconst= 0.033200
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      8.694471   charge   195.521058
 smvxcm (warning) mesh density negative at 162446 points:  rhomin=-6.77e-4
 smooth rhoeps = -987.546423 ( -73.508422,-914.038001)
         rhomu =-1309.654603 ( -48.994865,***********)
       avg vxc =   -0.192850 (  -0.156479,  -0.229221)
 smvxcm: negative smrho_w number,min(smrho_w)=      162446 -6.77138186736338178E-004
 smvxcm: enforce positive smrho_w. Add srshift=  6.77138186746338209E-004
 smooth rhoeps = -987.886458 ( -73.662011,-914.224446)
         rhomu =-1310.095943 ( -49.167915,***********)
       avg vxc =   -0.348692 (  -0.312582,  -0.384802)

 locpot:
  i job kmax lfltwf=           0           1           3 T

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.282610   -4.877596  -13.160206    -55.594942-1311.701998-1367.296941

 local terms:     true           smooth         local
 rhoeps:        -8.864856    -987.434253     978.569397
 rhomu:         -6.359894     -48.963788      42.603894
 spin2:         -5.306301   -1260.544174    1255.237873
 total:        -11.666195   -1309.507962    1297.841767
 val*vef       -13.160206   -1293.661842    1280.501636
 val chg:        3.128514     194.649572    -191.521058
 val mom:        0.833288    -163.351391     164.184679    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000157

 potpus  spin 1 : pnu = 2.967335 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -29.131051    7.148885    0.041667   -0.661517
 1      3.000000    1.000000    -5.887832    7.215357    0.125265   -0.609243
 2      3.000000    1.000000     6.000000   28.948505    0.488585   -0.089186
 3      3.000000    1.000000     9.000000   39.802438    0.568028   -0.057151

 potpus  spin 2 : pnu = 2.946917 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -17.822414    7.206386    0.062673   -0.637493
 1      3.000000    1.000000    -5.887832    7.089103    0.131674   -0.585229
 2      3.000000    1.000000     6.000000   29.309516    0.491109   -0.087353
 3      3.000000    1.000000     9.000000   40.051424    0.569415   -0.056555

 Energy terms:             smooth           local           total
   rhoval*vef          -1367.410350      1354.106735       -13.303615
   rhoval*ves             15.851667       -25.851346        -9.999679
   psnuc*ves               1.537274      -268.894255      -267.356981
   utot                    8.694471      -147.372801      -138.678330
   rho*exc              -987.886458       978.569397        -9.317060
   rho*vxc             -1310.095943      1297.841767       -12.254176
   valence chg           194.521058      -191.521058         3.000000
   valence mag          -163.684679       164.184679         0.500000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Read qp weights ...  ef=-0.824805
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
 -1.2698 -0.6530 -0.6530 -0.6530 -0.0541  0.4096  0.4096  0.4096
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1592 -0.5475 -0.5475 -0.5475  0.0306  0.4324  0.4324  0.4324
 Est Ef = -0.825 < evl(3)=-0.653 ... using qval=3.0, revise to -0.6530
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
 BZINTS: Fermi energy:     -0.654075;   3.000000 electrons
         Sum occ. bands:   -3.082995, incl. Bloechl correction: -0.000035
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.861574    5.615180   -2.753606     -0.908714    0.296359   -1.205073
       contr. to mm extrapolated for r>rmt:  -0.076189 est. true mm =-0.984902
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.43137  sum tc=    31.48010  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.39613  sum tc=    31.56723  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.976430   -1.258901    2.967335    2.926427    2.500000    2.926427
 spn 2 0    0.967059   -1.152994    2.946917    2.918327    2.500000    2.918327
 1     1    0.000000   -1.151375    2.850000    2.207623    2.250000    2.850000
 spn 2 1    0.918082   -0.542737    2.850000    2.895343    2.250000    2.850000
 2     0    0.000000   -1.421803    3.147584    3.103507    3.147584    3.147584
 spn 2 0    0.000002   -0.819844    3.147584    3.131716    3.147584    3.147584
 3     0    0.000000   -0.819844    4.102416    4.102416    4.102416    4.102416
 spn 2 0    0.000000   -0.777422    4.102416    4.095389    4.102416    4.102416

 Harris energy:
 sumev=       -3.082995  val*vef=     -13.303615   sumtv=      10.220620
 sumec=      -40.827502  cor*vef=    -103.760168   ttcor=      62.932666
 rhoeps=      -9.317060     utot=    -138.678330    ehar=     -74.842105
 smvxcm: negative smrho_w number,min(smrho_w)=      162446 -6.77142266578283146E-004
 smvxcm: enforce positive smrho_w. Add srshift=  6.77142266588283176E-004
 smvxcm: negative smrho_w number,min(smrho_w)=      162446 -6.77134106894393102E-004
 smvxcm: enforce positive smrho_w. Add srshift=  6.77134106904393133E-004
 smvxcm: negative smrho_w number,min(smrho_w)=      162446 -6.77138186736338070E-004
 smvxcm: enforce positive smrho_w. Add srshift=  6.77138186746338100E-004

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00    -0.00   -0.00   -0.00
 shift forces to make zero average correction:           -0.00   -0.00   -0.00

 srhov:    -24.709830     11.296011    -13.413818 sumev=   -3.082995   sumtv=   10.330824

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.776778   -2.753606     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.117576  avg sphere pot= 0.012802  vconst= 0.117576
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      4.984091   charge     6.753606
 smooth rhoeps =   -5.969489 (  -3.356446,  -2.613043)
         rhomu =   -7.809206 (  -4.517619,  -3.291587)
       avg vxc =   -0.171717 (  -0.148481,  -0.194953)
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -5.969489 (  -3.356446,  -2.613043)
         rhomu =   -7.809206 (  -4.517619,  -3.291587)
       avg vxc =   -0.171717 (  -0.148481,  -0.194953)

 locpot:
  i job kmax lfltwf=           0           0           3 F

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -5.015239   -8.483005  -13.498244    -13.056465   -9.100329  -22.156794

 local terms:     true           smooth         local
 rhoeps:        -8.775882      -5.945400      -2.830482
 rhomu:         -5.219355      -4.512720      -0.706635
 spin2:         -6.332715      -3.265293      -3.067422
 total:        -11.552070      -7.778013      -3.774057
 val*vef       -13.498244     -12.283640      -1.214604
 val chg:        2.861574       5.615180      -2.753606
 val mom:       -0.908714       0.296359      -1.205073    core:   0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef            -22.173849         8.658530       -13.515319
   rhoval*ves             -4.373934        -5.959643       -10.333577
   psnuc*ves              14.342115      -282.234549      -267.892434
   utot                    4.984091      -144.097096      -139.113006
   rho*exc                -5.969489        -2.830482        -8.799971
   rho*vxc                -7.809206        -3.774057       -11.583263
   valence chg             5.753606        -2.753606         3.000000
   valence mag             0.205073        -1.205073        -1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Kohn-Sham energy:
 sumtv=       10.330824  sumtc=        63.047332   ekin=       73.378156
 rhoep=       -8.799971   utot=      -139.113006   ehks=      -74.534821
 mag. mom=     1.000000  (bands)       -1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
 mixrho: sum smrho  init = 0.192727D+04-0.111133D-25 0.200901D+04   79215
 mixrho: sum smrnew new  = 0.372417D+03 0.198109D-16 0.372417D+03       0
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho: dqsum rmsuns= -0.18877D+00  0.22499D+01  0.26890D-18
 mixrealsmooth= T
 wgtsmooth=  2.82842712474619005E-003
 mixrho:  sought 2 iter from file mixm; read 1.  RMS DQ=1.82e0  last it=3.66e0
 mixrho: (warning) scr. and lin-mixed densities had 0 and 77403 negative points
 AMIX: nmix=1 mmix=8  nelts=252052  beta=0.5  tm=5  rmsdel=9.08e-1
   tj: 0.33193
 mixrealsmooth= T
 smrho qcell: add correction to smrho=  9.49614786804886535E-009  4.74807393402443408E-012
 unscreened rms difference:  smooth  3.181809   local  5.935998
   screened rms difference:  smooth  2.391301   local  5.935998   tot  1.816759
 mixrho: warning. negative smrho; isp number min=           1       88727 -5.53466594462584806E-004
 mixrho: warning. negative smrho; isp number min=           2       71163 -5.51088494175321249E-004

 iors  : write restart file (binary, mesh density) 

   it  2  of 10    ehf=       0.152795   ehk=       0.460079
 From last iter    ehf=      -0.754178   ehk=       1.819534
 diffe(q)=  0.906974 (1.816759)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=-.9999999 ehf=.1527953 ehk=.4600793

 --- BNDFP:  begin iteration 3 of 10 ---
 ttt nevmx w=           5  4.00000000000000008E-003

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1  -36.239533 -128.465799     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.048341  avg sphere pot= 0.021041  vconst= 0.048341
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      4.073958   charge   132.465799
 smvxcm (warning) mesh density negative at 159890 points:  rhomin=-5.53e-4
 smooth rhoeps = -582.251555 ( -46.559922,-535.691633)
         rhomu = -771.687195 ( -32.075626,-739.611569)
       avg vxc =   -0.180195 (  -0.148305,  -0.212086)
 smvxcm: negative smrho_w number,min(smrho_w)=      159890 -5.53466594463211800E-004
 smvxcm: enforce positive smrho_w. Add srshift=  5.53466594473211830E-004
 smooth rhoeps = -582.513852 ( -46.678786,-535.835067)
         rhomu = -772.027721 ( -32.212230,-739.815491)
       avg vxc =   -0.327366 (  -0.296363,  -0.358370)

 locpot:
  i job kmax lfltwf=           0           1           3 T

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -7.191790   -6.080942  -13.272732    -47.679690 -906.292752 -953.972443

 local terms:     true           smooth         local
 rhoeps:        -8.839982    -582.163925     573.323943
 rhomu:         -5.979477     -32.049753      26.070276
 spin2:         -5.653307    -739.523237     733.869929
 total:        -11.632785    -771.572990     759.940206
 val*vef       -13.272732    -769.195520     755.922788
 val chg:        3.082380     131.548178    -128.465799
 val mom:        0.251396    -108.687034     108.938430    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000104

 potpus  spin 1 : pnu = 2.926427 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.747401    7.666583    0.076309   -0.641935
 1      3.000000    1.000000    -5.887832    7.208214    0.127012   -0.601194
 2      3.000000    1.000000     6.000000   28.994092    0.488798   -0.088970
 3      3.000000    1.000000     9.000000   39.821781    0.568094   -0.057109

 potpus  spin 2 : pnu = 2.918327 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -11.434439    7.655653    0.082770   -0.632869
 1      3.000000    1.000000    -5.887832    7.143641    0.129413   -0.592962
 2      3.000000    1.000000     6.000000   29.160511    0.490066   -0.088102
 3      3.000000    1.000000     9.000000   39.946646    0.568838   -0.056804

 Energy terms:             smooth           local           total
   rhoval*vef           -954.082206       940.686296       -13.395911
   rhoval*ves              2.394991       -12.534448       -10.139457
   psnuc*ves               5.752924      -273.393662      -267.640738
   utot                    4.073958      -142.964055      -138.890097
   rho*exc              -582.513852       573.323943        -9.189909
   rho*vxc              -772.027721       759.940206       -12.087515
   valence chg           131.465799      -128.465799         3.000000
   valence mag          -108.939485       108.938430        -0.001055
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Read qp weights ...  ef=-0.654075
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
 -1.2303 -0.6126 -0.6126 -0.6126 -0.0274  0.4405  0.4405  0.4405
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1947 -0.5801 -0.5801 -0.5801  0.0537  0.4542  0.4542  0.4542
 Est Ef = -0.654 < evl(3)=-0.613 ... using qval=3.0, revise to -0.6126
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
 BZINTS: Fermi energy:     -0.613752;   3.000000 electrons
         Sum occ. bands:   -3.038813, incl. Bloechl correction: -0.000038
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.888571    6.188544   -3.299972      0.948219    2.201976   -1.253757
       contr. to mm extrapolated for r>rmt:   0.041135 est. true mm = 0.989354
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.38687  sum tc=    31.49394  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.37713  sum tc=    31.51993  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.975582   -1.219783    2.926427    2.925780    2.500000    2.925780
 spn 2 0    0.970176   -1.188875    2.918327    2.920587    2.500000    2.920587
 1     1    0.942812   -0.604234    2.850000    2.906655    2.250000    2.850000
 spn 2 1    0.000000   -0.877471    2.850000    2.356994    2.250000    2.850000
 2     0    0.000001   -0.840709    3.147584    3.130628    3.147584    3.147584
 spn 2 0    0.000000   -1.337307    3.147584    3.105963    3.147584    3.147584
 3     0    0.000000   -0.802047    4.102416    4.094620    4.102416    4.102416
 spn 2 0    0.000000   -0.802047    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -3.038813  val*vef=     -13.395911   sumtv=      10.357098
 sumec=      -40.763993  cor*vef=    -103.753966   ttcor=      62.989973
 rhoeps=      -9.189909     utot=    -138.890097    ehar=     -74.732937
 smvxcm: negative smrho_w number,min(smrho_w)=      159890 -5.53469938754726524E-004
 smvxcm: enforce positive smrho_w. Add srshift=  5.53469938764726554E-004
 smvxcm: negative smrho_w number,min(smrho_w)=      159890 -5.53463250171697075E-004
 smvxcm: enforce positive smrho_w. Add srshift=  5.53463250181697106E-004
 smvxcm: negative smrho_w number,min(smrho_w)=      159890 -5.53466594463211800E-004
 smvxcm: enforce positive smrho_w. Add srshift=  5.53466594473211830E-004

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00   -0.00     0.00    0.00    0.00    -0.00   -0.00   -0.00
 shift forces to make zero average correction:           -0.00   -0.00   -0.00

 srhov:    -27.698739     14.014049    -13.684690 sumev=   -3.038813   sumtv=   10.645877

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.930905   -3.299972     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.114187  avg sphere pot= 0.011682  vconst= 0.114187
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1  -0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.373447   charge     7.299972
 smooth rhoeps =   -6.909452 (  -4.656689,  -2.252763)
         rhomu =   -9.047070 (  -6.441934,  -2.605136)
       avg vxc =   -0.164128 (  -0.182967,  -0.145289)
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -6.909452 (  -4.656689,  -2.252763)
         rhomu =   -9.047070 (  -6.441934,  -2.605136)
       avg vxc =   -0.164128 (  -0.182967,  -0.145289)

 locpot:
  i job kmax lfltwf=           0           0           3 F

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.795789   -4.853153  -13.648942    -17.343249   -8.015060  -25.358309

 local terms:     true           smooth         local
 rhoeps:        -8.822773      -6.890714      -1.932059
 rhomu:         -6.411171      -6.423786       0.012615
 spin2:         -5.202603      -2.599021      -2.603582
 total:        -11.613774      -9.022807      -2.590967
 val*vef       -13.648942     -13.171821      -0.477121
 val chg:        2.888571       6.188544      -3.299972
 val mom:        0.948219       2.201976      -1.253757    core:   0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef            -25.372420        11.709342       -13.663079
   rhoval*ves             -4.024701        -6.405003       -10.429704
   psnuc*ves              14.771594      -282.860706      -268.089112
   utot                    5.373447      -144.632855      -139.259408
   rho*exc                -6.909452        -1.932059        -8.841511
   rho*vxc                -9.047070        -2.590967       -11.638037
   valence chg             6.299972        -3.299972         3.000000
   valence mag             2.253757        -1.253757         1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       10.645877  sumtc=        63.013864   ekin=       73.659741
 rhoep=       -8.841511   utot=      -139.259408   ehks=      -74.441178
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
 mixrho: sum smrho  init = 0.140789D+04-0.442965D-26 0.147081D+04   78291
 mixrho: sum smrnew new  = 0.534608D+03 0.163544D-16 0.534608D+03       0
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho: dqsum rmsuns= -0.12517D+00  0.14980D+01 -0.81564D-19
 mixrealsmooth= T
 wgtsmooth=  2.82842712474619005E-003
 mixrho:  sought 2 iter from file mixm; read 2.  RMS DQ=1.21e0  last it=1.82e0
 mixrho: (warning) scr. and lin-mixed densities had 0 and 76683 negative points
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.5  tm=5  rmsdel=6.04e-1
   tj:-0.43709   0.19366
 mixrealsmooth= T
 smrho qcell: add correction to smrho=  5.66376598953866051E-008  2.83188299476933103E-011
 unscreened rms difference:  smooth  2.118468   local  3.942372
   screened rms difference:  smooth  1.588171   local  3.942372   tot  1.208167
 mixrho: warning. negative smrho; isp number min=           1       83501 -4.59777418595691664E-004
 mixrho: warning. negative smrho; isp number min=           2       72387 -4.61695730056266014E-004

 iors  : write restart file (binary, mesh density) 

   it  3  of 10    ehf=       0.261963   ehk=       0.553722
 From last iter    ehf=       0.152795   ehk=       0.460079
 diffe(q)=  0.109168 (1.208167)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999999 ehf=.2619633 ehk=.5537221

 --- BNDFP:  begin iteration 4 of 10 ---
 ttt nevmx w=           5  4.00000000000000008E-003

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1  -21.595334  -76.553466     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.064933  avg sphere pot= 0.016920  vconst= 0.064933
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      2.852902   charge    80.553466
 smvxcm (warning) mesh density negative at 155888 points:  rhomin=-4.62e-4
 smooth rhoeps = -292.906145 ( -28.243534,-264.662611)
         rhomu = -387.801776 ( -21.073830,-366.727946)
       avg vxc =   -0.167014 (  -0.147042,  -0.186987)
 smvxcm: negative smrho_w number,min(smrho_w)=      155888 -4.61695730056266014E-004
 smvxcm: enforce positive smrho_w. Add srshift=  4.61695730066265990E-004
 smooth rhoeps = -293.118360 ( -28.343981,-264.774379)
         rhomu = -388.077522 ( -21.195727,-366.881795)
       avg vxc =   -0.314592 (  -0.296349,  -0.332836)

 locpot:
  i job kmax lfltwf=           0           1           3 T

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.652613   -4.848656  -13.501269    -39.414430 -534.915539 -574.329970

 local terms:     true           smooth         local
 rhoeps:        -8.875442    -292.845505     283.970064
 rhomu:         -6.429888     -21.044621      14.614732
 spin2:         -5.251602    -366.678545     361.426943
 total:        -11.681490    -387.723165     376.041675
 val*vef       -13.501269    -391.601418     378.100149
 val chg:        3.038906      79.592371     -76.553466
 val mom:        0.923677     -62.203564      63.127241    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000063

 potpus  spin 1 : pnu = 2.925780 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.632147    7.729098    0.076121   -0.645187
 1      3.000000    1.000000    -5.887832    7.245653    0.125729   -0.605594
 2      3.000000    1.000000     6.000000   28.900756    0.488179   -0.089446
 3      3.000000    1.000000     9.000000   39.759690    0.567776   -0.057256

 potpus  spin 2 : pnu = 2.920587 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -11.774389    7.590379    0.082707   -0.624366
 1      3.000000    1.000000    -5.887832    7.148537    0.132226   -0.580130
 2      3.000000    1.000000     6.000000   29.215438    0.490231   -0.087864
 3      3.000000    1.000000     9.000000   39.960040    0.568827   -0.056780

 Energy terms:             smooth           local           total
   rhoval*vef           -574.417405       560.823886       -13.593520
   rhoval*ves             -3.841706        -6.459369       -10.301074
   psnuc*ves               9.547510      -277.479048      -267.931539
   utot                    2.852902      -141.969208      -139.116306
   rho*exc              -293.118360       283.970064        -9.148296
   rho*vxc              -388.077522       376.041675       -12.035847
   valence chg            79.553466       -76.553466         3.000000
   valence mag           -62.300081        63.127241         0.827160
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Read qp weights ...  ef=-0.613752
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
 -1.2536 -0.6345 -0.6345 -0.6345 -0.0095  0.4576  0.4576  0.4576
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1316 -0.5160 -0.5160 -0.5160  0.0859  0.4965  0.4965  0.4965
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
 BZINTS: Fermi energy:     -0.635410;   3.000000 electrons
         Sum occ. bands:   -3.020563, incl. Bloechl correction: -0.000033
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.893040    6.398979   -3.505939      0.955754    2.455268   -1.499514
       contr. to mm extrapolated for r>rmt:   0.035142 est. true mm = 0.990895
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.35043  sum tc=    31.44673  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.30785  sum tc=    31.54390  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.976770   -1.243730    2.925780    2.926739    2.500000    2.926739
 spn 2 0    0.968643   -1.125228    2.920587    2.920028    2.500000    2.920028
 1     1    0.947626   -0.626319    2.850000    2.909012    2.250000    2.850000
 spn 2 1    0.000000   -1.048144    2.850000    2.210961    2.250000    2.850000
 2     0    0.000001   -0.853015    3.147584    3.130138    3.147584    3.147584
 spn 2 0    0.000000   -1.295773    3.147584    3.106170    3.147584    3.147584
 3     0    0.000000   -0.816657    4.102416    4.094301    4.102416    4.102416
 spn 2 0    0.000000   -0.816657    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -3.020563  val*vef=     -13.593520   sumtv=      10.572957
 sumec=      -40.658273  cor*vef=    -103.660203   ttcor=      63.001929
 rhoeps=      -9.148296     utot=    -139.116306    ehar=     -74.689717
 smvxcm: negative smrho_w number,min(smrho_w)=      155888 -4.61698520025606422E-004
 smvxcm: enforce positive smrho_w. Add srshift=  4.61698520035606398E-004
 smvxcm: negative smrho_w number,min(smrho_w)=      155888 -4.61692940086925660E-004
 smvxcm: enforce positive smrho_w. Add srshift=  4.61692940096925636E-004
 smvxcm: negative smrho_w number,min(smrho_w)=      155888 -4.61695730056266068E-004
 smvxcm: enforce positive smrho_w. Add srshift=  4.61695730066266044E-004

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00   -0.00   -0.00     0.00    0.00    0.00    -0.00   -0.00   -0.00
 shift forces to make zero average correction:           -0.00   -0.00   -0.00

 srhov:    -29.533336     15.873644    -13.659692 sumev=   -3.020563   sumtv=   10.639129

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.989007   -3.505939     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.113467  avg sphere pot= 0.011396  vconst= 0.113467
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.496561   charge     7.505939
 smooth rhoeps =   -7.276413 (  -5.026310,  -2.250103)
         rhomu =   -9.530565 (  -6.973448,  -2.557117)
       avg vxc =   -0.162625 (  -0.179970,  -0.145279)
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -7.276413 (  -5.026310,  -2.250103)
         rhomu =   -9.530565 (  -6.973448,  -2.557117)
       avg vxc =   -0.162625 (  -0.179970,  -0.145279)

 locpot:
  i job kmax lfltwf=           0           0           3 F

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.863018   -4.773680  -13.636698    -18.723924   -7.895084  -26.619008

 local terms:     true           smooth         local
 rhoeps:        -8.826422      -7.258506      -1.567916
 rhomu:         -6.424113      -6.956701       0.532588
 spin2:         -5.194457      -2.550677      -2.643780
 total:        -11.618570      -9.507377      -2.111192
 val*vef       -13.636698     -13.549006      -0.087692
 val chg:        2.893040       6.398979      -3.505939
 val mom:        0.955754       2.455268      -1.499514    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef            -26.632698        12.982283       -13.650415
   rhoval*ves             -3.918691        -6.493087       -10.411778
   psnuc*ves              14.911813      -282.954851      -268.043039
   utot                    5.496561      -144.723969      -139.227409
   rho*exc                -7.276413        -1.567916        -8.844329
   rho*vxc                -9.530565        -2.111192       -11.641757
   valence chg             6.505939        -3.505939         3.000000
   valence mag             2.499514        -1.499514         1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       10.639129  sumtc=        62.990624   ekin=       73.629753
 rhoep=       -8.844329   utot=      -139.227409   ehks=      -74.441985
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
 mixrho: sum smrho  init = 0.107834D+04-0.518992D-26 0.111966D+04   77259
 mixrho: sum smrnew new  = 0.562841D+03 0.433058D-16 0.562841D+03       0
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho: dqsum rmsuns= -0.73048D-01  0.87571D+00 -0.99988D-19
 mixrealsmooth= T
 wgtsmooth=  2.82842712474619005E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=7.07e-1  last it=1.21e0
 mixrho: (warning) scr. and lin-mixed densities had 0 and 75219 negative points
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.5  tm=5  rmsdel=3.53e-1
   tj:-1.03982  -0.16655
 mixrealsmooth= T
 smrho qcell: add correction to smrho=  3.44219586168037495E-008  1.72109793084018781E-011
 unscreened rms difference:  smooth  1.238438   local  2.306315
   screened rms difference:  smooth  0.927102   local  2.306315   tot  0.706533
 mixrho: warning. negative smrho; isp number min=           1       64623 -2.55808646163801286E-004
 mixrho: warning. negative smrho; isp number min=           2       77793 -2.60775372558712812E-004

 iors  : write restart file (binary, mesh density) 

   it  4  of 10    ehf=       0.305183   ehk=       0.552915
 From last iter    ehf=       0.261963   ehk=       0.553722
 diffe(q)=  0.043220 (0.706533)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999999 ehf=.3051833 ehk=.5529153

 --- BNDFP:  begin iteration 5 of 10 ---
 ttt nevmx w=           5  4.00000000000000008E-003

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.025626   -3.635750     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.094065  avg sphere pot= 0.011063  vconst= 0.094065
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.568279   charge     7.635750
 smvxcm (warning) mesh density negative at 142416 points:  rhomin=-2.61e-4
 smooth rhoeps =   -7.599117 (  -5.479070,  -2.120047)
         rhomu =   -9.957336 (  -7.619010,  -2.338326)
       avg vxc =   -0.126902 (  -0.142144,  -0.111661)
 smvxcm: negative smrho_w number,min(smrho_w)=      142416 -2.60775372558712812E-004
 smvxcm: enforce positive smrho_w. Add srshift=  2.60775372568712788E-004
 smooth rhoeps =   -7.704688 (  -5.535724,  -2.168964)
         rhomu =  -10.094191 (  -7.694592,  -2.399600)
       avg vxc =   -0.262286 (  -0.274243,  -0.250329)

 locpot:
  i job kmax lfltwf=           0           1           3 T

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 === rho1 valence true density ===
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       768
 vxcnsp (warning): negative rho: min val =  -2.14E-02
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -9.915885   -3.859139  -13.775024    -20.027358   -7.404885  -27.432243

 local terms:     true           smooth         local
 rhoeps:        -8.917356      -7.572903      -1.344453
 rhomu:         -6.800559      -7.592436       0.791876
 spin2:         -4.939435      -2.330921      -2.608514
 total:        -11.739995      -9.923357      -1.816638
 val*vef       -13.775024     -13.832957       0.057932
 val chg:        2.959541       6.595291      -3.635750
 val mom:        1.455950       3.027648      -1.571698    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000006

 potpus  spin 1 : pnu = 2.926739 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.803646    7.789231    0.074880   -0.648502
 1      3.000000    1.000000    -5.887832    7.299331    0.124481   -0.609178
 2      3.000000    1.000000     6.000000   28.777953    0.487297   -0.090091
 3      3.000000    1.000000     9.000000   39.670830    0.567291   -0.057471

 potpus  spin 2 : pnu = 2.920028 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -11.688573    7.618104    0.083945   -0.617015
 1      3.000000    1.000000    -5.887832    7.201795    0.133956   -0.570310
 2      3.000000    1.000000     6.000000   29.163079    0.489530   -0.088189
 3      3.000000    1.000000     9.000000   39.886121    0.568284   -0.056971

 Energy terms:             smooth           local           total
   rhoval*vef            -27.475882        13.657197       -13.818685
   rhoval*ves             -3.825157        -6.640807       -10.465964
   psnuc*ves              14.961715      -283.177942      -268.216226
   utot                    5.568279      -144.909374      -139.341095
   rho*exc                -7.704688        -1.344453        -9.049141
   rho*vxc               -10.094191        -1.816638       -11.910829
   valence chg             6.635750        -3.635750         3.000000
   valence mag             3.109671        -1.571698         1.537973
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.63541
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
 -1.2660 -0.6444 -0.6444 -0.6444  0.0386  0.5074  0.5074  0.5074
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.0709 -0.4509 -0.4509 -0.4509 -0.0192  0.5215  0.5215  0.5215
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
 BZINTS: Fermi energy:     -0.645140;   3.000000 electrons
         Sum occ. bands:   -2.981818, incl. Bloechl correction: -0.000026
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905115    7.955530   -5.050414      0.957569    1.587216   -0.629646
       contr. to mm extrapolated for r>rmt:   0.035149 est. true mm = 0.992719
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.30575  sum tc=    31.40592  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.23559  sum tc=    31.56320  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.978180   -1.257329    2.926739    2.928028    2.500000    2.928028
 spn 2 0    0.973773   -1.056402    2.920028    2.926084    2.500000    2.926084
 1     1    0.953162   -0.637063    2.850000    2.912034    2.250000    2.850000
 spn 2 1    0.000001   -0.843830    2.850000    2.284166    2.250000    2.850000
 2     0    0.000001   -0.857687    3.147584    3.129311    3.147584    3.147584
 spn 2 0    0.000000   -1.282113    3.147584    3.104622    3.147584    3.147584
 3     0    0.000000   -0.825422    4.102416    4.093772    4.102416    4.102416
 spn 2 0    0.000000   -0.825422    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -2.981818  val*vef=     -13.818685   sumtv=      10.836867
 sumec=      -40.541341  cor*vef=    -103.537588   ttcor=      62.996247
 rhoeps=      -9.049141     utot=    -139.341095    ehar=     -74.557122
 smvxcm: negative smrho_w number,min(smrho_w)=      142416 -2.60776936634210890E-004
 smvxcm: enforce positive smrho_w. Add srshift=  2.60776936644210866E-004
 smvxcm: negative smrho_w number,min(smrho_w)=      142416 -2.60773808483214680E-004
 smvxcm: enforce positive smrho_w. Add srshift=  2.60773808493214656E-004
 smvxcm: negative smrho_w number,min(smrho_w)=      142416 -2.60775372558712758E-004
 smvxcm: enforce positive smrho_w. Add srshift=  2.60775372568712734E-004

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00   -0.00   -0.00     0.00    0.00    0.00    -0.00   -0.00   -0.00
 shift forces to make zero average correction:           -0.00   -0.00   -0.00

 srhov:    -34.424322     20.735212    -13.689110 sumev=   -2.981818   sumtv=   10.707292

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.424696   -5.050414     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.111319  avg sphere pot= 0.010079  vconst= 0.111319
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.998008   charge     9.050414
 smooth rhoeps =   -9.853879 (  -5.789124,  -4.064755)
         rhomu =  -12.908559 (  -7.852340,  -5.056219)
       avg vxc =   -0.159260 (  -0.175687,  -0.142833)
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -9.853879 (  -5.789124,  -4.064755)
         rhomu =  -12.908559 (  -7.852340,  -5.056219)
       avg vxc =   -0.159260 (  -0.175687,  -0.142833)

 locpot:
  i job kmax lfltwf=           0           0           3 F

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.910617   -4.767885  -13.678502    -21.164166  -14.887139  -36.051305

 local terms:     true           smooth         local
 rhoeps:        -8.842706      -9.838148       0.995442
 rhomu:         -6.438564      -7.837290       1.398726
 spin2:         -5.201327      -5.050898      -0.150429
 total:        -11.639891     -12.888188       1.248297
 val*vef       -13.678502     -16.497963       2.819461
 val chg:        2.905115       7.955530      -5.050414
 val mom:        0.957569       1.587216      -0.629646    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef            -36.063650        22.372758       -13.690892
   rhoval*ves             -3.490474        -6.943279       -10.433753
   psnuc*ves              15.486490      -283.573862      -268.087372
   utot                    5.998008      -145.258571      -139.260563
   rho*exc                -9.853879         0.995442        -8.858438
   rho*vxc               -12.908559         1.248297       -11.660262
   valence chg             8.050414        -5.050414         3.000000
   valence mag             1.629646        -0.629646         1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       10.707292  sumtc=        62.969119   ekin=       73.676410
 rhoep=       -8.858438   utot=      -139.260563   ehks=      -74.442590
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
 mixrho: sum smrho  init = 0.609089D+03 0.788506D-26 0.619550D+03   70407
 mixrho: sum smrnew new  = 0.605004D+03 0.113154D-15 0.605004D+03       0
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho: dqsum rmsuns=  0.14147D-02  0.13993D-01  0.18048D-17
 mixrealsmooth= T
 wgtsmooth=  2.82842712474619005E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.12e-2  last it=7.07e-1
 mixrho: (warning) scr. and lin-mixed densities had 0 and 66411 negative points
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.5  tm=5  rmsdel=5.6e-3
   tj:-1.06135   0.63363
 mixrealsmooth= T
 smrho qcell: add correction to smrho=  4.54608972688674839E-008  2.27304486344337477E-011
 unscreened rms difference:  smooth  0.019789   local  0.037510
   screened rms difference:  smooth  0.015692   local  0.037510   tot  0.011193
 mixrho: warning. negative smrho; isp number min=           1       62847 -1.71937969314202299E-004
 mixrho: warning. negative smrho; isp number min=           2       74013 -1.74578364513986210E-004

 iors  : write restart file (binary, mesh density) 

   it  5  of 10    ehf=       0.437778   ehk=       0.552310
 From last iter    ehf=       0.305183   ehk=       0.552915
 diffe(q)=  0.132595 (0.011193)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999999 ehf=.4377782 ehk=.5523101

 --- BNDFP:  begin iteration 6 of 10 ---
 ttt nevmx w=           5  4.00000000000000008E-003

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.540350   -5.460398     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.099953  avg sphere pot= 0.010432  vconst= 0.099953
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.789593   charge     9.460398
 smvxcm (warning) mesh density negative at 136860 points:  rhomin=-1.75e-4
 smooth rhoeps =  -10.706392 (  -5.932004,  -4.774388)
         rhomu =  -14.028247 (  -7.943146,  -6.085101)
       avg vxc =   -0.129311 (  -0.139430,  -0.119191)
 smvxcm: negative smrho_w number,min(smrho_w)=      136860 -1.74578364513986210E-004
 smvxcm: enforce positive smrho_w. Add srshift=  1.74578364523986213E-004
 smooth rhoeps =  -10.774885 (  -5.968141,  -4.806744)
         rhomu =  -14.117115 (  -7.990764,  -6.126351)
       avg vxc =   -0.246980 (  -0.254589,  -0.239372)

 locpot:
  i job kmax lfltwf=           0           1           3 T

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 === rho1 valence true density ===
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       416
 vxcnsp (warning): negative rho: min val =  -6.21E-03
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -9.209690   -4.508389  -13.718079    -21.075384  -17.491307  -38.566691

 local terms:     true           smooth         local
 rhoeps:        -8.875510     -10.684863       1.809353
 rhomu:         -6.553891      -7.923027       1.369135
 spin2:         -5.129437      -6.077335       0.947898
 total:        -11.683328     -14.000362       2.317034
 val*vef       -13.718079     -17.759405       4.041326
 val chg:        2.941335       8.401733      -5.460398
 val mom:        1.105605       1.265370      -0.159765    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000009

 potpus  spin 1 : pnu = 2.928028 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -13.041167    7.772499    0.074170   -0.647768
 1      3.000000    1.000000    -5.887832    7.301201    0.124800   -0.607534
 2      3.000000    1.000000     6.000000   28.779474    0.487264   -0.090091
 3      3.000000    1.000000     9.000000   39.666703    0.567253   -0.057482

 potpus  spin 2 : pnu = 2.926084 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.685995    7.600883    0.078953   -0.624326
 1      3.000000    1.000000    -5.887832    7.225677    0.132156   -0.577026
 2      3.000000    1.000000     6.000000   29.073793    0.488950   -0.088635
 3      3.000000    1.000000     9.000000   39.828296    0.567988   -0.057107

 Energy terms:             smooth           local           total
   rhoval*vef            -38.601174        24.848565       -13.752609
   rhoval*ves             -3.663103        -6.791136       -10.454239
   psnuc*ves              15.242288      -283.409175      -268.166887
   utot                    5.789593      -145.100155      -139.310563
   rho*exc               -10.774885         1.809353        -8.965532
   rho*vxc               -14.117115         2.317034       -11.800081
   valence chg             8.460398        -5.460398         3.000000
   valence mag             1.318376        -0.159765         1.158611
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.64514
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
 -1.2494 -0.6274 -0.6274 -0.6274  0.1165  0.5499  0.5499  0.5499
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1012 -0.4803 -0.4803 -0.4803  0.1278  0.5785  0.5785  0.5785
 Est Ef = -0.645 < evl(3)=-0.627 ... using qval=3.0, revise to -0.6274
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
 BZINTS: Fermi energy:     -0.628040;   3.000000 electrons
         Sum occ. bands:   -2.978472, incl. Bloechl correction: -0.000026
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.902331    7.159012   -4.256681      0.957908    1.905463   -0.947555
       contr. to mm extrapolated for r>rmt:   0.034912 est. true mm = 0.992820
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.30547  sum tc=    31.42184  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.25208  sum tc=    31.54271  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.977118   -1.243197    2.928028    2.926823    2.500000    2.926823
 spn 2 0    0.972211   -1.093704    2.926084    2.923769    2.500000    2.923769
 1     1    0.953001   -0.621182    2.850000    2.912080    2.250000    2.850000
 spn 2 1    0.000000   -1.530493    2.850000    2.128768    2.250000    2.850000
 2     0    0.000001   -0.846859    3.147584    3.129196    3.147584    3.147584
 spn 2 0    0.000000   -1.288057    3.147584    3.104781    3.147584    3.147584
 3     0    0.000000   -0.816210    4.102416    4.093679    4.102416    4.102416
 spn 2 0    0.000000   -0.816210    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -2.978472  val*vef=     -13.752609   sumtv=      10.774137
 sumec=      -40.557546  cor*vef=    -103.540234   ttcor=      62.982689
 rhoeps=      -8.965532     utot=    -139.310563    ehar=     -74.519269
 smvxcm: negative smrho_w number,min(smrho_w)=      136860 -1.74579413670924457E-004
 smvxcm: enforce positive smrho_w. Add srshift=  1.74579413680924460E-004
 smvxcm: negative smrho_w number,min(smrho_w)=      136860 -1.74577315357047963E-004
 smvxcm: enforce positive smrho_w. Add srshift=  1.74577315367047966E-004
 smvxcm: negative smrho_w number,min(smrho_w)=      136860 -1.74578364513986210E-004
 smvxcm: enforce positive smrho_w. Add srshift=  1.74578364523986213E-004

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00    0.00   -0.00    -0.00    0.00    0.00     0.00   -0.00   -0.00
 shift forces to make zero average correction:            0.00   -0.00   -0.00

 srhov:    -31.678031     18.011071    -13.666961 sumev=   -2.978472   sumtv=   10.688489

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.200788   -4.256681     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.112048  avg sphere pot= 0.010653  vconst= 0.112048
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1  -0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.792165   charge     8.256681
 smooth rhoeps =   -8.464319 (  -5.278468,  -3.185850)
         rhomu =  -11.085227 (  -7.230896,  -3.854331)
       avg vxc =   -0.159233 (  -0.175691,  -0.142776)
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -8.464319 (  -5.278468,  -3.185850)
         rhomu =  -11.085227 (  -7.230896,  -3.854331)
       avg vxc =   -0.159233 (  -0.175691,  -0.142776)

 locpot:
  i job kmax lfltwf=           0           0           3 F

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.886079   -4.780137  -13.666216    -19.473005  -11.593770  -31.066775

 local terms:     true           smooth         local
 rhoeps:        -8.838184      -8.448062      -0.390122
 rhomu:         -6.434379      -7.215492       0.781113
 spin2:         -5.199573      -3.848686      -1.350887
 total:        -11.633952     -11.064177      -0.569774
 val*vef       -13.666216     -14.844347       1.178132
 val chg:        2.902331       7.159012      -4.256681
 val mom:        0.957908       1.905463      -0.947555    core:   0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef            -31.079507        17.400524       -13.678983
   rhoval*ves             -3.659839        -6.766476       -10.426315
   psnuc*ves              15.244170      -283.300639      -268.056470
   utot                    5.792165      -145.033557      -139.241392
   rho*exc                -8.464319        -0.390122        -8.854440
   rho*vxc               -11.085227        -0.569774       -11.655001
   valence chg             7.256681        -4.256681         3.000000
   valence mag             1.947555        -0.947555         1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       10.688489  sumtc=        62.964550   ekin=       73.653038
 rhoep=       -8.854440   utot=      -139.241392   ehks=      -74.442794
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
 mixrho: sum smrho  init = 0.611173D+03-0.696051D-26 0.617682D+03   68019
 mixrho: sum smrnew new  = 0.575265D+03 0.136251D-15 0.575265D+03       0
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho: dqsum rmsuns= -0.12037D-02  0.12340D-01  0.69184D-19
 mixrealsmooth= T
 wgtsmooth=  2.82842712474619005E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.01e-2  last it=1.12e-2
 mixrho: (warning) scr. and lin-mixed densities had 0 and 62215 negative points
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.5  tm=5  rmsdel=5.04e-3
   tj: 0.37484  -0.00304
 mixrealsmooth= T
 smrho qcell: add correction to smrho=  3.76850124439442880E-008  1.88425062219721487E-011
 unscreened rms difference:  smooth  0.017452   local  0.033580
   screened rms difference:  smooth  0.014761   local  0.033580   tot  0.010081
 mixrho: warning. negative smrho; isp number min=           1       60787 -1.52763684605454291E-004
 mixrho: warning. negative smrho; isp number min=           2       73725 -1.55327535150920890E-004

 iors  : write restart file (binary, mesh density) 

   it  6  of 10    ehf=       0.475631   ehk=       0.552106
 From last iter    ehf=       0.437778   ehk=       0.552310
 diffe(q)=  0.037853 (0.010081)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999999 ehf=.475631 ehk=.5521059

 --- BNDFP:  begin iteration 7 of 10 ---
 ttt nevmx w=           5  4.00000000000000008E-003

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.285859   -4.558252     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.101555  avg sphere pot= 0.010542  vconst= 0.101555
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.785453   charge     8.558252
 smvxcm (warning) mesh density negative at 134512 points:  rhomin=-1.55e-4
 smooth rhoeps =   -9.039132 (  -5.578527,  -3.460605)
         rhomu =  -11.839963 (  -7.627922,  -4.212041)
       avg vxc =   -0.128167 (  -0.139445,  -0.116889)
 smvxcm: negative smrho_w number,min(smrho_w)=      134512 -1.55327535150920890E-004
 smvxcm: enforce positive smrho_w. Add srshift=  1.55327535160920894E-004
 smooth rhoeps =   -9.099005 (  -5.610296,  -3.488709)
         rhomu =  -11.917639 (  -7.669888,  -4.247751)
       avg vxc =   -0.241035 (  -0.249616,  -0.232455)

 locpot:
  i job kmax lfltwf=           0           1           3 T

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 === rho1 valence true density ===
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       416
 vxcnsp (warning): negative rho: min val =  -5.86E-03
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -9.188253   -4.523822  -13.712074    -20.348104  -12.576537  -32.924641

 local terms:     true           smooth         local
 rhoeps:        -8.870969      -9.018307       0.147338
 rhomu:         -6.545324      -7.608281       1.062957
 spin2:         -5.132032      -4.204706      -0.927326
 total:        -11.677356     -11.812987       0.135630
 val*vef       -13.712074     -15.570119       1.858044
 val chg:        2.936324       7.494575      -4.558252
 val mom:        1.097635       1.951831      -0.854196    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000008

 potpus  spin 1 : pnu = 2.926823 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.818816    7.789622    0.074924   -0.647636
 1      3.000000    1.000000    -5.887832    7.303713    0.124701   -0.607902
 2      3.000000    1.000000     6.000000   28.772869    0.487217   -0.090126
 3      3.000000    1.000000     9.000000   39.662035    0.567228   -0.057494

 potpus  spin 2 : pnu = 2.923769 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.286341    7.633889    0.080450   -0.623986
 1      3.000000    1.000000    -5.887832    7.228401    0.131975   -0.577696
 2      3.000000    1.000000     6.000000   29.064314    0.488886   -0.088683
 3      3.000000    1.000000     9.000000   39.821933    0.567954   -0.057122

 Energy terms:             smooth           local           total
   rhoval*vef            -32.956966        19.212530       -13.744436
   rhoval*ves             -3.658165        -6.792423       -10.450588
   psnuc*ves              15.229072      -283.368189      -268.139118
   utot                    5.785453      -145.080306      -139.294853
   rho*exc                -9.099005         0.147338        -8.951667
   rho*vxc               -11.917639         0.135630       -11.782009
   valence chg             7.558252        -4.558252         3.000000
   valence mag             2.005106        -0.854196         1.150910
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.62804
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
 -1.2500 -0.6279 -0.6279 -0.6279  0.1293  0.5576  0.5576  0.5576
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1028 -0.4816 -0.4816 -0.4816  0.1309  0.5838  0.5838  0.5838
 Est Ef = -0.628 < evl(3)=-0.628 ... using qval=3.0, revise to -0.6279
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
 BZINTS: Fermi energy:     -0.628506;   3.000000 electrons
         Sum occ. bands:   -2.981149, incl. Bloechl correction: -0.000026
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.903011    7.193201   -4.290190      0.957999    1.853134   -0.895136
       contr. to mm extrapolated for r>rmt:   0.034970 est. true mm = 0.992969
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.30815  sum tc=    31.42288  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.25525  sum tc=    31.54291  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.977094   -1.244206    2.926823    2.926759    2.500000    2.926759
 spn 2 0    0.972506   -1.095408    2.923769    2.924023    2.500000    2.924023
 1     1    0.953411   -0.621877    2.850000    2.912304    2.250000    2.850000
 spn 2 1    0.000000   -1.536323    2.850000    2.128312    2.250000    2.850000
 2     0    0.000001   -0.847461    3.147584    3.129121    3.147584    3.147584
 spn 2 0    0.000000   -1.291166    3.147584    3.104646    3.147584    3.147584
 3     0    0.000000   -0.817359    4.102416    4.093630    4.102416    4.102416
 spn 2 0    0.000000   -0.817359    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -2.981149  val*vef=     -13.744436   sumtv=      10.763287
 sumec=      -40.563397  cor*vef=    -103.537017   ttcor=      62.973619
 rhoeps=      -8.951667     utot=    -139.294853    ehar=     -74.509614
 smvxcm: negative smrho_w number,min(smrho_w)=      134512 -1.55328467967095102E-004
 smvxcm: enforce positive smrho_w. Add srshift=  1.55328467977095105E-004
 smvxcm: negative smrho_w number,min(smrho_w)=      134512 -1.55326602334746679E-004
 smvxcm: enforce positive smrho_w. Add srshift=  1.55326602344746682E-004
 smvxcm: negative smrho_w number,min(smrho_w)=      134512 -1.55327535150920890E-004
 smvxcm: enforce positive smrho_w. Add srshift=  1.55327535160920894E-004

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00    0.00   -0.00    -0.00   -0.00    0.00     0.00    0.00   -0.00
 shift forces to make zero average correction:            0.00    0.00   -0.00

 srhov:    -31.401746     17.721664    -13.680081 sumev=   -2.981149   sumtv=   10.698932

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.210240   -4.290190     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.111980  avg sphere pot= 0.010621  vconst= 0.111980
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.802412   charge     8.290190
 smooth rhoeps =   -8.518079 (  -5.270107,  -3.247972)
         rhomu =  -11.155512 (  -7.209904,  -3.945609)
       avg vxc =   -0.158906 (  -0.175308,  -0.142504)
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -8.518079 (  -5.270107,  -3.247972)
         rhomu =  -11.155512 (  -7.209904,  -3.945609)
       avg vxc =   -0.158906 (  -0.175308,  -0.142504)

 locpot:
  i job kmax lfltwf=           0           0           3 F

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.888958   -4.785059  -13.674017    -19.419364  -11.845669  -31.265033

 local terms:     true           smooth         local
 rhoeps:        -8.839872      -8.501943      -0.337929
 rhomu:         -6.435609      -7.194588       0.758979
 spin2:         -5.200560      -3.940030      -1.260530
 total:        -11.636169     -11.134618      -0.501551
 val*vef       -13.674017     -14.904158       1.230141
 val chg:        2.903011       7.193201      -4.290190
 val mom:        0.957999       1.853134      -0.895136    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef            -31.277693        17.590980       -13.686713
   rhoval*ves             -3.649362        -6.782918       -10.432280
   psnuc*ves              15.254187      -283.324985      -268.070798
   utot                    5.802412      -145.053952      -139.251539
   rho*exc                -8.518079        -0.337929        -8.856008
   rho*vxc               -11.155512        -0.501551       -11.657063
   valence chg             7.290190        -4.290190         3.000000
   valence mag             1.895136        -0.895136         1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       10.698932  sumtc=        62.965796   ekin=       73.664729
 rhoep=       -8.856008   utot=      -139.251539   ehks=      -74.442818
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
 mixrho: sum smrho  init = 0.597710D+03 0.119919D-25 0.603236D+03   66699
 mixrho: sum smrnew new  = 0.574083D+03-0.112316D-15 0.574083D+03       0
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho: dqsum rmsuns= -0.26806D-03  0.23641D-02  0.13434D-17
 mixrealsmooth= T
 wgtsmooth=  2.82842712474619005E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.92e-3  last it=1.01e-2
 mixrho: (warning) scr. and lin-mixed densities had 0 and 60283 negative points
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.5  tm=5  rmsdel=9.59e-4
   tj: 0.21818   0.27205
 mixrealsmooth= T
 smrho qcell: add correction to smrho=  3.58767850983099379E-008  1.79383925491549729E-011
 unscreened rms difference:  smooth  0.003343   local  0.006298
   screened rms difference:  smooth  0.003294   local  0.006298   tot  0.001917
 mixrho: warning. negative smrho; isp number min=           1       59293 -1.39449783317469414E-004
 mixrho: warning. negative smrho; isp number min=           2       73029 -1.41855471483889070E-004

 iors  : write restart file (binary, mesh density) 

   it  7  of 10    ehf=       0.485286   ehk=       0.552082
 From last iter    ehf=       0.475631   ehk=       0.552106
 diffe(q)=  0.009655 (0.001917)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999999 ehf=.4852865 ehk=.5520821

 --- BNDFP:  begin iteration 8 of 10 ---
 ttt nevmx w=           5  4.00000000000000008E-003

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.268554   -4.496906     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.102544  avg sphere pot= 0.010570  vconst= 0.102544
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.780045   charge     8.496906
 smvxcm (warning) mesh density negative at 132322 points:  rhomin=-1.42e-4
 smooth rhoeps =   -8.925420 (  -5.517621,  -3.407799)
         rhomu =  -11.690620 (  -7.546681,  -4.143939)
       avg vxc =   -0.128181 (  -0.139424,  -0.116938)
 smvxcm: negative smrho_w number,min(smrho_w)=      132322 -1.41855471483889070E-004
 smvxcm: enforce positive smrho_w. Add srshift=  1.41855471493889073E-004
 smooth rhoeps =   -8.979549 (  -5.546355,  -3.433193)
         rhomu =  -11.760848 (  -7.584643,  -4.176205)
       avg vxc =   -0.237316 (  -0.245875,  -0.228756)

 locpot:
  i job kmax lfltwf=           0           1           3 T

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 === rho1 valence true density ===
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       352
 vxcnsp (warning): negative rho: min val =  -4.98E-03
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -9.144899   -4.563241  -13.708139    -20.154145  -12.387327  -32.541472

 local terms:     true           smooth         local
 rhoeps:        -8.867357      -8.905095       0.037737
 rhomu:         -6.529996      -7.527565       0.997569
 spin2:         -5.142546      -4.136726      -1.005820
 total:        -11.672542     -11.664291      -0.008251
 val*vef       -13.708139     -15.427729       1.719590
 val chg:        2.933238       7.430145      -4.496906
 val mom:        1.077358       1.943472      -0.866114    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000008

 potpus  spin 1 : pnu = 2.926759 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.807262    7.791889    0.074955   -0.647661
 1      3.000000    1.000000    -5.887832    7.304937    0.124678   -0.607958
 2      3.000000    1.000000     6.000000   28.770079    0.487195   -0.090141
 3      3.000000    1.000000     9.000000   39.659756    0.567215   -0.057499

 potpus  spin 2 : pnu = 2.924023 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.328990    7.635887    0.080200   -0.624531
 1      3.000000    1.000000    -5.887832    7.230701    0.131815   -0.578295
 2      3.000000    1.000000     6.000000   29.056164    0.488833   -0.088724
 3      3.000000    1.000000     9.000000   39.816631    0.567927   -0.057135

 Energy terms:             smooth           local           total
   rhoval*vef            -32.572326        18.833296       -13.739029
   rhoval*ves             -3.662557        -6.786802       -10.449359
   psnuc*ves              15.222648      -283.350015      -268.127367
   utot                    5.780045      -145.068408      -139.288363
   rho*exc                -8.979549         0.037737        -8.941811
   rho*vxc               -11.760848        -0.008251       -11.769099
   valence chg             7.496906        -4.496906         3.000000
   valence mag             1.995060        -0.866114         1.128946
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.628506
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
 -1.2493 -0.6272 -0.6272 -0.6272  0.1403  0.5638  0.5638  0.5638
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1048 -0.4835 -0.4835 -0.4835  0.1435  0.5904  0.5904  0.5904
 Est Ef = -0.629 < evl(3)=-0.627 ... using qval=3.0, revise to -0.6272
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
 BZINTS: Fermi energy:     -0.627752;   3.000000 electrons
         Sum occ. bands:   -2.981706, incl. Bloechl correction: -0.000026
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.903108    7.156136   -4.253028      0.958122    1.856204   -0.898082
       contr. to mm extrapolated for r>rmt:   0.034924 est. true mm = 0.993046
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.30907  sum tc=    31.42404  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.25716  sum tc=    31.54193  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.977018   -1.243883    2.926759    2.926658    2.500000    2.926658
 spn 2 0    0.972493   -1.097868    2.924023    2.923951    2.500000    2.923951
 1     1    0.953596   -0.621298    2.850000    2.912417    2.250000    2.850000
 spn 2 1    0.000000   -1.538046    2.850000    2.128225    2.250000    2.850000
 2     0    0.000001   -0.847324    3.147584    3.129072    3.147584    3.147584
 spn 2 0    0.000000   -1.291318    3.147584    3.104647    3.147584    3.147584
 3     0    0.000000   -0.817648    4.102416    4.093596    4.102416    4.102416
 spn 2 0    0.000000   -0.817648    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -2.981706  val*vef=     -13.739029   sumtv=      10.757323
 sumec=      -40.566231  cor*vef=    -103.535938   ttcor=      62.969708
 rhoeps=      -8.941811     utot=    -139.288363    ehar=     -74.503143
 smvxcm: negative smrho_w number,min(smrho_w)=      132322 -1.41856323199467913E-004
 smvxcm: enforce positive smrho_w. Add srshift=  1.41856323209467916E-004
 smvxcm: negative smrho_w number,min(smrho_w)=      132322 -1.41854619768310255E-004
 smvxcm: enforce positive smrho_w. Add srshift=  1.41854619778310258E-004
 smvxcm: negative smrho_w number,min(smrho_w)=      132322 -1.41855471483889097E-004
 smvxcm: enforce positive smrho_w. Add srshift=  1.41855471493889101E-004

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00   -0.00   -0.00    -0.00   -0.00    0.00     0.00    0.00   -0.00
 shift forces to make zero average correction:            0.00    0.00   -0.00

 srhov:    -31.161159     17.477989    -13.683170 sumev=   -2.981706   sumtv=   10.701464

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.199757   -4.253028     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.111996  avg sphere pot= 0.010644  vconst= 0.111996
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1  -0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.794281   charge     8.253028
 smooth rhoeps =   -8.453399 (  -5.237493,  -3.215906)
         rhomu =  -11.070602 (  -7.166543,  -3.904059)
       avg vxc =   -0.158769 (  -0.175137,  -0.142401)
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -8.453399 (  -5.237493,  -3.215906)
         rhomu =  -11.070602 (  -7.166543,  -3.904059)
       avg vxc =   -0.158769 (  -0.175137,  -0.142401)

 locpot:
  i job kmax lfltwf=           0           0           3 F

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.888923   -4.787033  -13.675956    -19.302226  -11.730681  -31.032907

 local terms:     true           smooth         local
 rhoeps:        -8.840229      -8.437278      -0.402951
 rhomu:         -6.435895      -7.151248       0.715352
 spin2:         -5.200743      -3.898479      -1.302264
 total:        -11.636639     -11.049727      -0.586912
 val*vef       -13.675956     -14.824974       1.149018
 val chg:        2.903108       7.156136      -4.253028
 val mom:        0.958122       1.856204      -0.898082    core:   0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef            -31.045564        17.356915       -13.688649
   rhoval*ves             -3.655069        -6.778759       -10.433828
   psnuc*ves              15.243630      -283.317644      -268.074013
   utot                    5.794281      -145.048201      -139.253920
   rho*exc                -8.453399        -0.402951        -8.856350
   rho*vxc               -11.070602        -0.586912       -11.657514
   valence chg             7.253028        -4.253028         3.000000
   valence mag             1.898082        -0.898082         1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       10.701464  sumtc=        62.965975   ekin=       73.667439
 rhoep=       -8.856350   utot=      -139.253920   ehks=      -74.442831
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
 mixrho: sum smrho  init = 0.593248D+03 0.297322D-26 0.598177D+03   65955
 mixrho: sum smrnew new  = 0.571944D+03-0.162830D-15 0.571944D+03       0
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho: dqsum rmsuns= -0.24388D-03  0.21428D-02  0.97698D-18
 mixrealsmooth= T
 wgtsmooth=  2.82842712474619005E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.73e-3  last it=1.92e-3
 mixrho: (warning) scr. and lin-mixed densities had 0 and 58981 negative points
 AMIX: Reducing nmix to  1: t_j exceeds tm: tj=-11.03506   0.05512
 AMIX: Reducing nmix to  0: t_j exceeds tm: tj= -8.76362
 AMIX: nmix=0 mmix=8  nelts=252052  beta=0.5  tm=5  rmsdel=8.67e-4
 mixrealsmooth= T
 smrho qcell: add correction to smrho=  3.68831223340748693E-008  1.84415611670374402E-011
 unscreened rms difference:  smooth  0.003030   local  0.005691
   screened rms difference:  smooth  0.002994   local  0.005691   tot  0.001734
 mixrho: warning. negative smrho; isp number min=           1       54937 -1.04846429104787109E-004
 mixrho: warning. negative smrho; isp number min=           2       71037 -1.06804420186464401E-004

 iors  : write restart file (binary, mesh density) 

   it  8  of 10    ehf=       0.491757   ehk=       0.552069
 From last iter    ehf=       0.485286   ehk=       0.552082
 diffe(q)=  0.006471 (0.001734)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999999 ehf=.491757 ehk=.5520687

 --- BNDFP:  begin iteration 9 of 10 ---
 ttt nevmx w=           5  4.00000000000000008E-003

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.234155   -4.374967     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.105096  avg sphere pot= 0.010607  vconst= 0.105096
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.779965   charge     8.374967
 smvxcm (warning) mesh density negative at 125974 points:  rhomin=-1.07e-4
 smooth rhoeps =   -8.694477 (  -5.380217,  -3.314260)
         rhomu =  -11.387191 (  -7.360189,  -4.027002)
       avg vxc =   -0.128661 (  -0.139816,  -0.117506)
 smvxcm: negative smrho_w number,min(smrho_w)=      125974 -1.06804420186464401E-004
 smvxcm: enforce positive smrho_w. Add srshift=  1.06804420196464405E-004
 smooth rhoeps =   -8.734117 (  -5.401279,  -3.332838)
         rhomu =  -11.438632 (  -7.388020,  -4.050612)
       avg vxc =   -0.226669 (  -0.235111,  -0.218228)

 locpot:
  i job kmax lfltwf=           0           1           3 T

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 === rho1 valence true density ===
 vxcns2 (warning): nr*np=     11808  negative density # of points=         0       256
 vxcnsp (warning): negative rho: min val =  -2.41E-03
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -9.018854   -4.678500  -13.697354    -19.722959  -12.057394  -31.780353

 local terms:     true           smooth         local
 rhoeps:        -8.858181      -8.675480      -0.182701
 rhomu:         -6.486060      -7.342501       0.856441
 spin2:         -5.174232      -4.020083      -1.154149
 total:        -11.660293     -11.362584      -0.297708
 val*vef       -13.697354     -15.130960       1.433605
 val chg:        2.925294       7.300261      -4.374967
 val mom:        1.017740       1.899838      -0.882098    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000008

 potpus  spin 1 : pnu = 2.926658 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.789073    7.796714    0.075001   -0.647685
 1      3.000000    1.000000    -5.887832    7.308025    0.124635   -0.608022
 2      3.000000    1.000000     6.000000   28.763317    0.487139   -0.090178
 3      3.000000    1.000000     9.000000   39.654028    0.567181   -0.057514

 potpus  spin 2 : pnu = 2.923951 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.316938    7.650179    0.080039   -0.625719
 1      3.000000    1.000000    -5.887832    7.237099    0.131385   -0.579904
 2      3.000000    1.000000     6.000000   29.033982    0.488687   -0.088836
 3      3.000000    1.000000     9.000000   39.802155    0.567852   -0.057169

 Energy terms:             smooth           local           total
   rhoval*vef            -31.807119        18.082962       -13.724157
   rhoval*ves             -3.662453        -6.783429       -10.445882
   psnuc*ves              15.222383      -283.335275      -268.112892
   utot                    5.779965      -145.059352      -139.279387
   rho*exc                -8.734117        -0.182701        -8.916818
   rho*vxc               -11.438632        -0.297708       -11.736340
   valence chg             7.374967        -4.374967         3.000000
   valence mag             1.946571        -0.882098         1.064473
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.627752
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
 -1.2471 -0.6247 -0.6247 -0.6247  0.1703  0.5812  0.5812  0.5812
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1104 -0.4886 -0.4886 -0.4886  0.1784  0.6087  0.6087  0.6087
 Est Ef = -0.628 < evl(3)=-0.625 ... using qval=3.0, revise to -0.6247
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
 BZINTS: Fermi energy:     -0.625280;   3.000000 electrons
         Sum occ. bands:   -2.982556, incl. Bloechl correction: -0.000025
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.903388    7.066817   -4.163429      0.958428    1.864225   -0.905797
       contr. to mm extrapolated for r>rmt:   0.034826 est. true mm = 0.993254
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.31077  sum tc=    31.42691  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.26171  sum tc=    31.53851  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.976820   -1.242551    2.926658    2.926401    2.500000    2.926401
 spn 2 0    0.972480   -1.104694    2.923951    2.923783    2.500000    2.923783
 1     1    0.954088   -0.619284    2.850000    2.912726    2.250000    2.850000
 spn 2 1    0.000000   -1.132091    2.850000    2.180339    2.250000    2.850000
 2     0    0.000001   -0.846934    3.147584    3.128927    3.147584    3.147584
 spn 2 0    0.000000   -1.291774    3.147584    3.104655    3.147584    3.147584
 3     0    0.000000   -0.818571    4.102416    4.093498    4.102416    4.102416
 spn 2 0    0.000000   -0.818571    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -2.982556  val*vef=     -13.724157   sumtv=      10.741601
 sumec=      -40.572479  cor*vef=    -103.540320   ttcor=      62.967841
 rhoeps=      -8.916818     utot=    -139.279387    ehar=     -74.486763
 smvxcm: negative smrho_w number,min(smrho_w)=      125974 -1.06805061007503253E-004
 smvxcm: enforce positive smrho_w. Add srshift=  1.06805061017503257E-004
 smvxcm: negative smrho_w number,min(smrho_w)=      125974 -1.06803779365425563E-004
 smvxcm: enforce positive smrho_w. Add srshift=  1.06803779375425566E-004
 smvxcm: negative smrho_w number,min(smrho_w)=      125974 -1.06804420186464415E-004
 smvxcm: enforce positive smrho_w. Add srshift=  1.06804420196464418E-004

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1   -0.00   -0.00   -0.00    -0.00   -0.00   -0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:    -30.610878     16.920415    -13.690463 sumev=   -2.982556   sumtv=   10.707907

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.174482   -4.163429     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.112028  avg sphere pot= 0.010696  vconst= 0.112028
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.776380   charge     8.163429
 smooth rhoeps =   -8.297733 (  -5.159056,  -3.138677)
         rhomu =  -10.866253 (  -7.062205,  -3.804048)
       avg vxc =   -0.158408 (  -0.174696,  -0.142121)
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -8.297733 (  -5.159056,  -3.138677)
         rhomu =  -10.866253 (  -7.062205,  -3.804048)
       avg vxc =   -0.158408 (  -0.174696,  -0.142121)

 locpot:
  i job kmax lfltwf=           0           0           3 F

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.888247   -4.792633  -13.680880    -19.019991  -11.454600  -30.474591

 local terms:     true           smooth         local
 rhoeps:        -8.841097      -8.281655      -0.559441
 rhomu:         -6.436551      -7.046964       0.610413
 spin2:         -5.201227      -3.798470      -1.402757
 total:        -11.637779     -10.845434      -0.792345
 val*vef       -13.680880     -14.632777       0.951897
 val chg:        2.903388       7.066817      -4.163429
 val mom:        0.958428       1.864225      -0.905797    core:   0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef            -30.487237        16.793677       -13.693560
   rhoval*ves             -3.667176        -6.770567       -10.437743
   psnuc*ves              15.219935      -283.300210      -268.080274
   utot                    5.776380      -145.035388      -139.259008
   rho*exc                -8.297733        -0.559441        -8.857174
   rho*vxc               -10.866253        -0.792345       -11.658598
   valence chg             7.163429        -4.163429         3.000000
   valence mag             1.905797        -0.905797         1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       10.707907  sumtc=        62.965415   ekin=       73.673322
 rhoep=       -8.857174   utot=      -139.259008   ehks=      -74.442861
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
 mixrho: sum smrho  init = 0.582596D+03-0.492091D-26 0.585986D+03   62887
 mixrho: sum smrnew new  = 0.566827D+03 0.532446D-16 0.566827D+03       0
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho: dqsum rmsuns= -0.21154D-03  0.18008D-02 -0.31394D-18
 mixrealsmooth= T
 wgtsmooth=  2.82842712474619005E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.46e-3  last it=1.73e-3
 mixrho: (warning) scr. and lin-mixed densities had 0 and 53905 negative points
 AMIX: Reducing nmix to  1: t_j exceeds tm: tj= 13.85692 -11.27722
 AMIX: nmix=1 mmix=8  nelts=252052  beta=0.5  tm=5  rmsdel=7.29e-4
   tj:-3.92182
 mixrealsmooth= T
 smrho qcell: add correction to smrho=  4.51779920140893410E-008  2.25889960070446757E-011
 unscreened rms difference:  smooth  0.002547   local  0.004836
   screened rms difference:  smooth  0.002542   local  0.004836   tot  0.001458
 mixrho: all smrho is positive for isp=           1
 mixrho: all smrho is positive for isp=           2

 iors  : write restart file (binary, mesh density) 

   it  9  of 10    ehf=       0.508137   ehk=       0.552039
 From last iter    ehf=       0.491757   ehk=       0.552069
 diffe(q)=  0.016380 (0.001458)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999999 ehf=.5081366 ehk=.5520395

 --- BNDFP:  begin iteration 10 of 10 ---
 ttt nevmx w=           5  4.00000000000000008E-003

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.087303   -3.854390     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.114386  avg sphere pot= 0.010826  vconst= 0.114386
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.746080   charge     7.854390
 smooth rhoeps =   -7.749294 (  -4.852003,  -2.897291)
         rhomu =  -10.146029 (  -6.648028,  -3.498001)
       avg vxc =   -0.175816 (  -0.185981,  -0.165652)
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -7.749294 (  -4.852003,  -2.897291)
         rhomu =  -10.146029 (  -6.648028,  -3.498001)
       avg vxc =   -0.175816 (  -0.185981,  -0.165652)

 locpot:
  i job kmax lfltwf=           0           1           3 T

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.704411   -4.972297  -13.676709    -17.983955  -10.578687  -28.562642

 local terms:     true           smooth         local
 rhoeps:        -8.832884      -7.733447      -1.099437
 rhomu:         -6.375556      -6.633738       0.258182
 spin2:         -5.251166      -3.491761      -1.759405
 total:        -11.626722     -10.125499      -1.501223
 val*vef       -13.676709     -13.936760       0.260051
 val chg:        2.896660       6.751051      -3.854390
 val mom:        0.871778       1.812198      -0.940420    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000
 potential shift to crystal energy zero:    0.000008

 potpus  spin 1 : pnu = 2.926401 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.742744    7.816758    0.075050   -0.648088
 1      3.000000    1.000000    -5.887832    7.321585    0.124382   -0.608636
 2      3.000000    1.000000     6.000000   28.732686    0.486897   -0.090345
 3      3.000000    1.000000     9.000000   39.629302    0.567039   -0.057574

 potpus  spin 2 : pnu = 2.923783 2.850000 3.147584 4.102416

 l,E        a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.288712    7.693896    0.079540   -0.629154
 1      3.000000    1.000000    -5.887832    7.258719    0.130135   -0.584511
 2      3.000000    1.000000     6.000000   28.964677    0.488224   -0.089189
 3      3.000000    1.000000     9.000000   39.755929    0.567611   -0.057280

 Energy terms:             smooth           local           total
   rhoval*vef            -28.572486        14.885901       -13.686585
   rhoval*ves             -3.686222        -6.754226       -10.440447
   psnuc*ves              15.178381      -283.263021      -268.084640
   utot                    5.746080      -145.008623      -139.262544
   rho*exc                -7.749294        -1.099437        -8.848732
   rho*vxc               -10.146029        -1.501223       -11.647252
   valence chg             6.854390        -3.854390         3.000000
   valence mag             1.846230        -0.940420         0.905811
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.62528
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
 -1.2418 -0.6186 -0.6186 -0.6186  0.2812  0.6499  0.6499  0.6499
  mode napw           0           0
  end of hambls mode=           0
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1237 -0.5005 -0.5005 -0.5005  0.3094  0.6854  0.6854  0.6854
 Est Ef = -0.625 < evl(3)=-0.619 ... using qval=3.0, revise to -0.6186
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
 BZINTS: Fermi energy:     -0.618955;   3.000000 electrons
         Sum occ. bands:   -2.984378, incl. Bloechl correction: -0.000021
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905960    6.973334   -4.067374      0.960394    1.911038   -0.950644
       contr. to mm extrapolated for r>rmt:   0.033868 est. true mm = 0.994262
 getcor:  qcore=  2.00  qsc=  2.00  konf = 2  2  3  4 
 sum q= 1.00  sum ec=   -20.31528  sum tc=    31.43310  rho(rmax) 0.00000
 sum q= 1.00  sum ec=   -20.27345  sum tc=    31.52976  rho(rmax) 0.00000

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.976529   -1.240010    2.926401    2.925958    2.500000    2.925958
 spn 2 0    0.972783   -1.121814    2.923783    2.923663    2.500000    2.923663
 1     1    0.956647   -0.614284    2.850000    2.914348    2.250000    2.850000
 spn 2 1    0.000000   -0.467222    2.850000    2.933368    2.250000    2.850000
 2     0    0.000001   -0.850938    3.147584    3.128152    3.147584    3.147584
 spn 2 0    0.000000   -1.291608    3.147584    3.104626    3.147584    3.147584
 3     0    0.000000   -0.829272    4.102416    4.092985    4.102416    4.102416
 spn 2 0    0.000000   -0.829272    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -2.984378  val*vef=     -13.686585   sumtv=      10.702208
 sumec=      -40.588739  cor*vef=    -103.555366   ttcor=      62.966626
 rhoeps=      -8.848732     utot=    -139.262544    ehar=     -74.442441
 smvxcm: all smrho_w is positive
 smvxcm: all smrho_w is positive
 smvxcm: all smrho_w is positive

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00   -0.00   -0.00     0.00    0.00   -0.00    -0.00   -0.00    0.00
 shift forces to make zero average correction:           -0.00   -0.00    0.00

 srhov:    -29.712303     15.985091    -13.727211 sumev=   -2.984378   sumtv=   10.742834

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.147385   -4.067374     2.00     6.00

 after vesgcomp: forces are:
   1    0.000000    0.000000    0.000000

 avg es pot at rmt=-0.111828  avg sphere pot= 0.010674  vconst= 0.111828
 average electrostatic potential at MT boundaries after shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159

 smooth rhoves      5.792884   charge     8.067374
 smooth rhoeps =   -8.130327 (  -5.094313,  -3.036014)
         rhomu =  -10.646619 (  -6.980340,  -3.666279)
       avg vxc =   -0.156663 (  -0.172492,  -0.140835)
 smvxcm: all smrho_w is positive
 smooth rhoeps =   -8.130327 (  -5.094313,  -3.036014)
         rhomu =  -10.646619 (  -6.980340,  -3.666279)
       avg vxc =   -0.156663 (  -0.172492,  -0.140835)

 locpot:
  i job kmax lfltwf=           0           0           3 F

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=F
 === rho1 valence true density ===
 === rho2 valence counter density ===
 === rhol1 valence+core density ===
 === rho2 ->valence+smooth core density ===

 ilm                   rho*vtrue                              rho*vsm
             spin1       spin2       tot           spin1       spin2       tot
   1      -8.896464   -4.809956  -13.706420    -18.778853  -11.092398  -29.871251

 local terms:     true           smooth         local
 rhoeps:        -8.846752      -8.114689      -0.732063
 rhomu:         -6.441776      -6.965582       0.523806
 spin2:         -5.203431      -3.660787      -1.542644
 total:        -11.645207     -10.626369      -1.018838
 val*vef       -13.706420     -14.387478       0.681058
 val chg:        2.905960       6.973334      -4.067374
 val mom:        0.960394       1.911038      -0.950644    core:  -0.000000
 core chg:       2.000000       2.000000       0.000000

 Energy terms:             smooth           local           total
   rhoval*vef            -29.883680        16.164798       -13.718882
   rhoval*ves             -3.641492        -6.815268       -10.456760
   psnuc*ves              15.227259      -283.342877      -268.115618
   utot                    5.792884      -145.079072      -139.286189
   rho*exc                -8.130327        -0.732063        -8.862390
   rho*vxc               -10.646619        -1.018838       -11.665458
   valence chg             7.067374        -4.067374         3.000000
   valence mag             1.950644        -0.950644         1.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:     -0.00000

 Kohn-Sham energy:
 sumtv=       10.742834  sumtc=        62.962862   ekin=       73.705695
 rhoep=       -8.862390   utot=      -139.286189   ehks=      -74.442884
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
 mixrho: sum smrho  init = 0.543789D+03-0.672898D-27 0.543789D+03       0
 mixrho: sum smrnew new  = 0.563626D+03-0.177545D-15 0.563626D+03       0
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho: dqsum rmsuns=  0.21298D-03  0.17297D-02 -0.29075D-17
 mixrealsmooth= T
 wgtsmooth=  2.82842712474619005E-003
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.41e-3  last it=1.46e-3
 AMIX: nmix=2 mmix=8  nelts=252052  beta=0.5  tm=5  rmsdel=7.07e-4
   tj: 0.07092   0.38421
 mixrealsmooth= T
 smrho qcell: add correction to smrho=  3.66357824077567784E-008  1.83178912038783933E-011
 unscreened rms difference:  smooth  0.002446   local  0.004882
   screened rms difference:  smooth  0.002369   local  0.004882   tot  0.001414
 mixrho: warning. negative smrho; isp number min=           1       31067 -3.51962521678338483E-005
 mixrho: warning. negative smrho; isp number min=           2       59275 -3.65796172025638816E-005

 iors  : write restart file (binary, mesh density) 

   it 10  of 10    ehf=       0.552459   ehk=       0.552016
 From last iter    ehf=       0.508137   ehk=       0.552039
 diffe(q)=  0.044322 (0.001414)    tol= 0.000010 (0.000500)   more=F
x zbak=1 mmom=.9999999 ehf=.5524587 ehk=.5520161
 Exit 0 LMF 
