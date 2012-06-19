HOST_INFORMATION platform: gfortran
HOST_INFORMATION compiler version: gcc バージョン 4.6.1 (Ubuntu/Linaro 4.6.1-9ubuntu3) 
HOST_INFORMATION FFLAGS (<=120): -O2 -fomit-frame-pointer -funroll-loops -ffast-math -ffixed-line-length-132 -DHASIARGC -DHASGETARG -DFDATE -DHASCPUTIME 
HOST_INFORMATION LIBLOC (<=120): /usr/lib/libfftw3.so.3 /usr/lib/liblapack.so.3gf /usr/lib/libblas.so.3gf
HOST_INFORMATION uname -a (<=120): Linux TT4 3.0.0-20-generic #34-Ubuntu SMP Tue May 1 17:24:39 UTC 2012 x86_64 x86_64 x86_64 GNU/Linux
HOST_INFORMATION /etc/issue: Ubuntu 11.10 \n \l
HOST_INFORMATION git branch: refs/heads/newaniso
HOST_INFORMATION git commit: ff4f0b57505a84e0eaabe44d640190810c88782c
HOST_INFORMATION linked at: Tue Jun 19 14:13:38 JST 2012
 -----------------------  START LMFA     -----------------------
 ptenv() is called with EXT=cu
 ptenv() not supported, but continue.
 HEADER Cu
 Cu       xxx            1           1

 rdctrl: reset global max nl from 5 to 6
  mxcst switch =           1           0 F F F
  LMFA  vn 7.00(LMFA 7.0)  verb 31,20
 end of rdctrl2 in imfav7
 lattic:

                Plat                                  Qlat
   0.000000   0.500000   0.500000       -1.000000   1.000000   1.000000
   0.500000   0.000000   0.500000        1.000000  -1.000000   1.000000
   0.500000   0.500000   0.000000        1.000000   1.000000  -1.000000
  Cell vol= 78.538660

 LATTC:  as= 2.000   tol=1.00E-08   alat= 6.79800   awald= 0.467
         r1=  1.959   nkd= 135      q1=  5.910   nkg= 181
 goto mksym

 SGROUP: 1 symmetry operations from 0 generators
 SYMLAT: Bravais system is cubic with 48 symmetry operations.
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-5
 GROUPG: the following are sufficient to generate the space group:
         i*r3(1,1,-1) r4x
         i*r3(1,1,-1) r4x
 MKSYM:  found 48 space group operations ... includes inversion
 zzz nclass=           1
 end of mksym x
 goto defspc
 end of defspc
 goto freeat

ttt: pnu qat=  1  0     4.690     1.000
ttt: pnu qat=  1  1     4.420     0.000
ttt: pnu qat=  1  2     3.880    10.000
ttt: pnu qat=  1  3     4.120     0.000
ttt: pnu qat=  1  4     5.100     0.000
ttt: pnu qat=  1  5     6.100     0.000
 NOTE: when we have two valence: P and PZ, We assume eigen(PZ) is deeper than eigen(P).
=== Charge for l     : Qtot=Qv=     0   1.000   1.000
=== Charge for l     : Qtot=Qv=     1   0.000   0.000
=== Charge for l     : Qtot=Qv=     2  10.000  10.000
=== Charge for l     : Qtot=Qv=     3   0.000   0.000
=== Charge for l     : Qtot=Qv=     4   0.000   0.000
=== Charge for l     : Qtot=Qv=     5   0.000   0.000

conf:SPEC_ATOM= Cu : --- Table for atomic configuration ---
conf When int(P)z .ne. int(P), Qval: Q for MTOcore(PZ)+MTO(P)
conf:  isp  l  int(P) int(P)z    Qval     Qcore   CoreConf
conf:    1  0       4  4         1.000    6.000 => 1,2,3,
conf:    1  1       4  4         0.000   12.000 => 2,3,
conf:    1  2       3  3        10.000    0.000 => 
conf:    1  3       4  4         0.000    0.000 => 
conf:    1  4       5  5         0.000    0.000 => 
conf:    1  5       6  6         0.000    0.000 => 
conf:-----------------------------------------------------

 Species Cu:  Z=29  Qc=18  R=2.280000  Q=0
 mesh:   rmt=2.280000  rmax=48.629375  a=0.015  nr=655  nr(rmax)=859

  iter     qint         drho          vh0          rho0          vsum     beta
 NOTE: rhocor: core density is spin-independent now for any MMOM, june2012.
       We use spin-avaraged potential to calculate core density rhoc.
       Thus diff of core eigen is pot. due to valence electron. I/O:iofa.F
    1   29.000000   7.875E+03      145.0000    0.1442E+03      -58.3040   0.30
   52   29.000000   4.902E-05      274.8263    0.2697E+05     -130.8188   0.30


 sumev=-4.333259  etot=-3304.416215  eref=0.000000

 Free-atom wavefunctions:
 valence:      eval       node at      max at       c.t.p.   rho(r>rmt)
   4s      -0.36411         0.890       2.256       3.581     0.655341
   4p      -0.06295         0.975       3.484       7.412     0.906557
   3d      -0.39692         0.000       0.600       3.428     0.058644
   4f       0.01964         0.000      35.282      48.629*    1.000000
   5g       0.02702         0.000      36.739      48.629*    1.000000
   65       0.03539         0.000      37.851      48.629*    1.000000

 core:        ecore       node at      max at       c.t.p.   rho(r>rmt)
   1s    -649.07635         0.000       0.034       0.069     0.000000
   2s     -77.91381         0.070       0.197       0.308     0.000000
   2p     -67.32533         0.000       0.158       0.335     0.000000
   3s      -8.39248         0.288       0.614       0.895     0.000167
   3p      -5.29682         0.260       0.619       1.078     0.000836

 Optimise free-atom basis for species Cu, rmt=2.28
 l  it    Rsm      Eh     stiffR   stiffE      Eval      Exact     Pnu    Ql
 0  10   1.140  -0.104       0.0   6532.8   -0.34064  -0.36411    4.75   1.00
 1  11   1.140  -0.100       0.0  -3648.1    0.18334  -0.06295    4.55   0.00
 2   7   1.140  -0.497       0.0     21.3   -0.39473  -0.39692    3.89  10.00
 l=4  increase Pnu=   5.091  to    5.100
 l=5  increase Pnu=   6.069  to    6.100
 eigenvalue sum:  exact  -4.33326    opt basis  -4.28791    error 0.04535
 tailsm: init

 tailsm: fit tails to 6 smoothed hankels, rmt= 2.28000, rsm= 1.14000
  ---E:energies of smHankels. C:fitting coeeficient for core tail. ---
 E:    -1.00000    -2.00000    -4.00000    -6.00000    -9.00000    -15.0000
 C:    -0.07224    10.68573    -181.432    1154.646    -4317.83    18314.70
        r          rho         fit         diff
    2.280000    0.018887    0.018853    0.000034
    2.649002    0.009672    0.009671    0.000002
    3.077722    0.004737    0.004737    0.000000
    3.575823    0.002169    0.002173   -0.000004
    4.154534    0.000906    0.000904    0.000002
    4.826900    0.000336    0.000336   -0.000000
    5.608078    0.000108    0.000110   -0.000002
    6.515677    0.000029    0.000029   -0.000000
    7.570157    0.000006    0.000006    0.000000
    8.795288    0.000001    0.000001    0.000000
    q(fit):     1.241803    rms diff:   0.000017
    fit: r>rmt  1.241803   r<rmt  3.486916   qtot  4.728719
    rho: r>rmt  1.241803   r<rmt  9.758197   qtot 11.000000
 tailsm: end

 coretail: q=0.00464, rho(rmt)=0.00535.  Fit with Hankel e=-24.113  coeff=770.4
      r            rhoc          fit
    2.280000    0.02411508    0.02411508
    2.457586    0.01084310    0.01086777
    2.649002    0.00457091    0.00457610
    2.855327    0.00179819    0.00179084
    3.077722    0.00065680    0.00064766
    3.317437    0.00022148    0.00021513
    3.575823    0.00006853    0.00006520
    3.854333    0.00001932    0.00001790
    4.154534    0.00000493    0.00000442
    4.478116    0.00000113    0.00000097
    4.826900    0.00000023    0.00000019
    5.202849    0.00000004    0.00000003
 end of freats: spid=Cu      

  Write mtopara.* ...

 Sum of reference energies: 0
 Exit 0 LMFA 
