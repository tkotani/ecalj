 -----------------------  START LMFA (80000K)  -----------------------
 HEADER Cu

 rdctrl: reset global max nl from 5 to 6

 LMFA:     alat = 6.798  nbas = 1  nspec = 1  vn 7.00(LMFA 7.0)  verb 31,20
 pot:      XC:BH

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
 
conf:SPEC_ATOM= Cu : --- Table for atomic configuration ---
conf int(P)z = int(P) where P is replaced by PZ if it is semicore
conf:  isp  l  int(P) int(P)z    Qval    Qcore   CoreConf
conf:    1  0       4  4        1.000    6.000 => 1,2,3,
conf:    1  1       4  4        0.000   12.000 => 2,3,
conf:    1  2       3  3       10.000    0.000 => 
conf:    1  3       4  4        0.000    0.000 => 
conf:    1  4       5  5        0.000    0.000 => 
conf:    1  5       6  6        0.000    0.000 => 
conf:-----------------------------------------------------

 Species Cu:  Z=29  Qc=18  R=2.280000  Q=0
 mesh:   rmt=2.280000  rmax=48.145529  a=0.025  nr=393  nr(rmax)=515
  Pl=  4.5     4.5     3.5     4.5     5.5     6.5    
  Ql=  1.0     0.0     10.0    0.0     0.0     0.0    

  iter     qint         drho          vh0          rho0          vsum     beta
    1   29.000000   4.725E+03      145.0000    0.1442E+03      -58.2780   0.30
   51   29.000000   4.214E-05      274.8263    0.2633E+05     -130.7925   0.30


 sumev=-4.333255  etot=-3304.416258  eref=0.000000

 Free-atom wavefunctions:
 valence:      eval       node at      max at       c.t.p.   rho(r>rmt)
   4s      -0.36411         0.890       2.256       3.581     0.655341
   4p      -0.06295         0.975       3.484       7.413     0.906557
   3d      -0.39691         0.000       0.600       3.429     0.058644
   4f       0.02001         0.000      34.923      48.146*    1.000000
   5g       0.02754         0.000      36.368      48.146*    1.000000
   65       0.03607         0.000      37.464      48.146*    1.000000

 core:        ecore       node at      max at       c.t.p.   rho(r>rmt)
   1s    -649.07634         0.000       0.034       0.069     0.000000
   2s     -77.91382         0.070       0.197       0.308     0.000000
   2p     -67.32532         0.000       0.158       0.335     0.000000
   3s      -8.39248         0.288       0.614       0.895     0.000167
   3p      -5.29682         0.260       0.619       1.078     0.000836

 Optimise free-atom basis for species Cu, rmt=2.28
 l  it    Rsm      Eh     stiffR   stiffE      Eval      Exact     Pnu    Ql
 0   9   2.280  -0.284     108.7    530.9   -0.36392  -0.36411    4.75   1.00
 ... rsm exceeded rmt*2/3 .. repeat with rsm=rmt
 0   5   1.520  -0.129     108.7   4188.6   -0.35161  -0.36411    4.75   1.00
 1  11   2.280  -0.100     159.1      6.2   -0.04874  -0.06295    4.55   0.00
 ... rsm exceeded rmt*2/3 .. repeat with rsm=rmt
 1   1   1.520  -0.100     159.1  -2032.6    0.05991  -0.06295    4.55   0.00
 2  27   0.951  -0.112     175.3    114.6   -0.39668  -0.39691    3.89  10.00
 l=4  increase Pnu=   5.091  to    5.100
 l=5  increase Pnu=   6.069  to    6.100
 eigenvalue sum:  exact  -4.33326    opt basis  -4.31841    error 0.01484

 tailsm: fit tails to 6 smoothed hankels, rmt= 2.28000, rsm= 1.14000
 E:    -1.00000    -2.00000    -4.00000    -6.00000    -9.00000    -15.0000
 C:    -0.07287    10.70262    -181.908    1157.772    -4329.13    18357.95
        r          rho         fit         diff
    2.280000    0.018887    0.018854    0.000034
    2.927614    0.006048    0.006042    0.000005
    3.759167    0.001640    0.001641   -0.000001
    4.826901    0.000336    0.000336    0.000000
    6.197900    0.000046    0.000047   -0.000001
    7.958297    0.000004    0.000003    0.000000
    q(fit):     1.241803    rms diff:   0.000017
    fit: r>rmt  1.241803   r<rmt  3.487816   qtot  4.729619
    rho: r>rmt  1.241803   r<rmt  9.758197   qtot 11.000000

 coretail: q=0.00464, rho(rmt)=0.00535.  Fit with Hankel e=-24.113  coeff=770.|
      r            rhoc          fit
    2.280000    0.02411513    0.02411513
    2.396905    0.01425235    0.01427892
    2.716066    0.00337589    0.00337554
    3.077722    0.00065680    0.00064769
    3.487533    0.00010234    0.00009811
    3.951910    0.00001240    0.00001137
    4.478117    0.00000113    0.00000097
    5.074388    0.00000007    0.00000006
 
  Write mtopara.* ...

 Sum of reference energies: 0
 Exit 0 LMFA 
 wkinfo:  used    94 K  workspace of 80000 K   in   0 K calls
