Cstruc_makeinit.py , make *_init, *_mpibc1.  2010-06-24 13:00:30.897363
Cstruc_makeinit.py , make *_init, *_mpibc1.  2010-06-18 11:22:12.488426
C this is translated from lstra.F in lmv7-beta
C Kino April2009.
      module m_struc_def

      type s_ri
      integer  ::   t
      integer(8)  ::   i
      real(8)  ::   r
      end type

Cr --- explanation of  ubz ---
Cr   off offe  name    purpose
Cr    2    1   def      Uncertainty in Fermi level
Cr    3    2   dosw     Energy window over which DOS accumulated
Cr    5    3   ef       Fermi level
Cr    6    4   efmax    eigenvalues above efmax are not needed
Cr    7    5   fsmom    fixed-spin moment (fixed spin moment method)
Cr    8    6   lcond    conductivity, and directions
Cr   12    7   lio      1  read qp
Cr                      2  write qp
Cr                      4  write dos
Cr                      8  density matrix
Cr                      16 save evcs
Cr                      32 save dos weights
Cr                      NB: 4,32 out of date ; see ldos
Cr   13    8   lmet     0 insulator   1 3-point scheme  2 save evcs
Cr                      3 two passes  4 wt from prior iter
Cr                      (numbers, not powers of 2)
Cr   14    9   lmull    not used
Cr   15   10   lopt     2 BZMAP option
Cr   16   11   lshft    shift mesh in each of three directions
Cr   19   12   n        Polynomial order for Methfessel-Paxton sampling
Cr   20   13   ndos     No. DOS points (sampling integration, and lmdos)
Cr   21   14   nevmx    eigenvectors above nevmx are not needed
Cr   22   15   nkabc    number of divisions in qp mesh
Cr   25   16   nkp      Number of qp
Cr   26   17   ntet     number of tetrahedra
Cr   27   18   odos     offset to total dos
Cr   28   19   oidtet   offset to tetrahedron corners
Cr   29   20   oipq     offset to ipq array made in bzmesh
Cr   30   21   opdos    offset to partial dos
Cr   31   22   oqp      offset to qp array
Cr   32   23   ostar    offset to array containing info on star of k
Cr   33   24   owtkp    offset to array containing qp weights
Cr   34   25   range    number of FWHM for sampling integration
Cr   35   26   semsh    nz, modec, emin,emax, ecc,eps delta
Cr   45   27   size     size of this structure
Cr   46   28   stnr     parameters for Stoner model (not used)
Cr   47   29   w        Line broadening for sampling integration
Cr   48   30   zval     No. electrons to accumulate in BZ integration
C ----------------------------------------------------------------
C       data  /2,3,5,6,7,8,12,13,14,15,16,19,20,21,22,25,26,27,28,        29,30,31,32,33,34,35,45,46,47,48,49/
C       data /4,4,4,4,4,4,18*2,4,4,-1,4,4,4/
C  1   2  4        def  1
C  2   3  4       dosw  2
C  3   5  4         ef  1
C  4   6  4      efmax  1
C  5   7  4      fsmom  1
C  6   8  4      lcond  4
C  7  12  2        lio  1
C  8  13  2       lmet  1
C  9  14  2      lmull  1
C 10  15  2       lopt  1
C 11  16  2      lshft  3
C 12  19  2          n  1
C 13  20  2       ndos  1
C 14  21  2      nevmx  1
C 15  22  2      nkabc  3
C 16  25  2        nkp  1
C 17  26  2       ntet  1
C 18  27  2       odos  1 possible_structure
C 19  28  2     oidtet  1 possible_structure
C 20  29  2       oipq  1 possible_structure
C 21  30  2      opdos  1 possible_structure
C 22  31  2        oqp  1 possible_structure
C 23  32  2      ostar  1 possible_structure
C 24  33  2      owtkp  1 possible_structure
C 25  34  4      range  1
C 26  35  4      semsh 10
C 27  45 -1       size  1
C 28  46  4       stnr  1
C 29  47  4          w  1
C 30  48  4       zval  1
      type s_bz
      real(8)  ::   size
      real(8)   ::   def
      real(8)   ::   dosw(2)
      real(8)   ::   ef
      real(8)   ::   efmax
      real(8)   ::   fsmom
      real(8)   ::   lcond(4)
      integer(8)   ::   lio
      integer(8)   ::   lmet
      integer(8)   ::   lmull
      integer(8)   ::   lopt
      integer(8)   ::   lshft(3)
      integer(8)   ::   n
      integer(8)   ::   ndos
      integer(8)   ::   nevmx
      integer(8)   ::   nkabc(3)
      integer(8)   ::   nkp
      integer(8)   ::   ntet
CCHSTRUC       integer(8)   ::   odos ! not used because of undef
      integer(8)   ::   oidtet
CCHSTRUC       integer(8)   ::   oipq
       integer , pointer ::  iv_p_oipq (:)
CCHSTRUC       integer(8)   ::   opdos ! not used because of undef
CCHSTRUC       integer(8)   ::   oqp
       real(8) , pointer ::  rv_p_oqp (:)
      integer(8)   ::   ostar
CCHSTRUC       integer(8)   ::   owtkp
       real(8) , pointer ::  rv_p_owtkp (:)
      real(8)   ::   range
      real(8)   ::   semsh(10)
      real(8)   ::   stnr
      real(8)   ::   w
      real(8)   ::   zval
      end type s_bz

Cr --- explanation of  uarray ---
Cr    off lsti  name    purpose
Cr    2    1   ncl     (O-N) number of clusters
Cr    3    2   nclasp  (pgf) number of classes including PL -1,npl
Cr    4    3   nofgl   (pgf) number of left  PL for off-diagonal GF
Cr    5    4   nofgr   (pgf) number of right PL for off-diagonal GF
Cr    6    5   npadl   (pgf) number of basis atoms in PL -2,-1,0
Cr    7    6   npadr   (pgf) number of basis atoms in PL npl-1,npl,npl+1
Cr    8    7   npl     (pgf) number of principal layers
Cr    9    8   oclabl  offset to class label array (r8tos8 format)
Cr   10    9   ocllst  (O-N) list of sites belonging to each cluster
Cr   11   10   oclp    (O-N) cluster-related dimensioning parameters
Cr   12   11   oclssl  (O-N) site ib belongs to cluster clssl(ib)
Cr   13   12   ogroup
Cr   14   13   ogrp2
Cr   15   14   ohave   Switches containing what info available, by class
Cr                     1 P,Q  2 pp  4 sop  8 vintra  16 pmpol  32 gradrm
Cr   16   15   oics    Index to which species each class belongs
Cr   17   16   oifrlx  (dynamics)
Cr   18   17   oipc    index to which class each site belongs
Cr   19   18   oipcp   index to which class each site belongs
Cr   20   19   oips    index to which species each site belongs
Cr   21   20   omxcst  mixing constraint
Cr   22   21   onrc    number of sites in each class
Cr   23   22   onrcp   number of sites in each class
Cr   24   23   opgfsl
Cr   25   24   opgfvl
Cr   26   25   opgord
Cr   27   26   opgplp
Cr   28   27   opos
Cr   29   28   ormax   sphere radius, by class
Cr   30   29   osordn  structure for order-N information
Cr   31   30   ostni   Stoner
Cr   32   31   ozos    Stoner DOS
Cr   33   32   size    size of this structure
C ----------------------------------------------------------------
C       data  /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,        22,23,24,25,26,27,28,29,30,31,32,33,34/
C       data / 31*2,-1/
C  1   2  2        ncl  1
C  2   3  2     nclasp  1
C  3   4  2      nofgl  1
C  4   5  2      nofgr  1
C  5   6  2      npadl  1
C  6   7  2      npadr  1
C  7   8  2        npl  1
C  8   9  2     oclabl  1 possible_structure
C  9  10  2     ocllst  1 possible_structure
C 10  11  2       oclp  1 possible_structure
C 11  12  2     oclssl  1 possible_structure
C 12  13  2     ogroup  1 possible_structure
C 13  14  2      ogrp2  1 possible_structure
C 14  15  2      ohave  1 possible_structure
C 15  16  2       oics  1 possible_structure
C 16  17  2     oifrlx  1 possible_structure
C 17  18  2       oipc  1 possible_structure
C 18  19  2      oipcp  1 possible_structure
C 19  20  2       oips  1 possible_structure
C 20  21  2     omxcst  1 possible_structure
C 21  22  2       onrc  1 possible_structure
C 22  23  2      onrcp  1 possible_structure
C 23  24  2     opgfsl  1 possible_structure
C 24  25  2     opgfvl  1 possible_structure
C 25  26  2     opgord  1 possible_structure
C 26  27  2     opgplp  1 possible_structure
C 27  28  2       opos  1 possible_structure
C 28  29  2      ormax  1 possible_structure
C 29  30  2     osordn  1 possible_structure
C 30  31  2      ostni  1 possible_structure
C 31  32  2       ozos  1 possible_structure
C 32  33 -1       size  1
      type s_array
      real(8)  ::   size
      integer(8)   ::   ncl
      integer(8)   ::   nclasp
      integer(8)   ::   nofgl
      integer(8)   ::   nofgr
      integer(8)   ::   npadl
      integer(8)   ::   npadr
      integer(8)   ::   npl
CCHSTRUC       integer(8)   ::   oclabl
       real(8) , pointer ::  rv_p_oclabl (:)
CCHSTRUC       integer(8)   ::   ocllst ! not used because of undef
CCHSTRUC       integer(8)   ::   oclp ! not used because of undef
CCHSTRUC       integer(8)   ::   oclssl ! not used because of undef
      integer(8)   ::   ogroup
CCHSTRUC       integer(8)   ::   ogrp2 ! not used because of undef
CCHSTRUC       integer(8)   ::   ohave
       integer , pointer ::  iv_p_ohave (:)
CCHSTRUC       integer(8)   ::   oics
       integer , pointer ::  iv_p_oics (:)
CCHSTRUC       integer(8)   ::   oifrlx ! not used because of undef
      integer(8)   ::   oipc
      integer(8)   ::   oipcp
      integer(8)   ::   oips
CCHSTRUC       integer(8)   ::   omxcst ! not used because of undef
      integer(8)   ::   onrc
      integer(8)   ::   onrcp
CCHSTRUC       integer(8)   ::   opgfsl ! not used because of undef
CCHSTRUC       integer(8)   ::   opgfvl ! not used because of undef
CCHSTRUC       integer(8)   ::   opgord ! not used because of undef
CCHSTRUC       integer(8)   ::   opgplp ! not used because of undef
      integer(8)   ::   opos
      integer(8)   ::   ormax
CCHSTRUC       integer(8)   ::   osordn ! not used because of undef
CCHSTRUC       integer(8)   ::   ostni ! not used because of undef
CCHSTRUC       integer(8)   ::   ozos ! not used because of undef
      end type s_array

Cr --- explanation of  uctrl ---
Cr   off offe  name    purpose
Cr    2    1   defm    Lines for deformation when relaxing plat
Cr    8    2   elin    ccor linearization energy, 2C hamiltonian
Cr    9    3   lasa    1 Make V from P,Q  2 Make pp  4 ccor  8 free atm
Cr                     16 map  32 nonspherical mpol moms 64 MT corr
Cr                     128 interpretation of sphere Q2; see newrho.f
Cr                     256 how ASA Q1,Q2 are accumulated; see makwts.f
Cr                     512 (spin pol) alpha repsn = (gamma(1) + gamma(nsp))/2
Cr   10    4   lbas    1 Hamiltonian has no screening transformation
Cr                     2 Hamiltonian is nfp style
Cr                     16 freeze phi,phidot for all species
Cr   11    5   lcd     1 freeze core
Cr                     2 non-self-consistent Harris
Cr                     4 represent full potential density on a uniform mesh
Cr                     8 represent full potential density via TCF
Cr                    16 unused
Cr                    32 unused
Cr                    64 (molecules) XC potential by FFT
Cr   12    6   lcgf    Green's function
Cr                    >0 Green's function calculation
Cr   13    7   ldos    1 make dos
Cr                     2 generate weights for partial dos
Cr                     4 generate weights for m-decompos'n of pdos
Cr                     8 BZ map
Cr   14    8   lfp     switches for full-potential
Cr                     1 program uses full potential
Cr   15    9   lfrce   How forces are calculated
Cr                     0 do not calculate forces
Cr                     1 shift in FA density
Cr                     2 shift in core+nucleus
Cr                    10 added for shift to be screened
Cr   16   10   lgen3   switches for third generation LMTO
Cr                     1 3rd generation LMTO
Cr                     2 make structure constants by Ewald summation
Cr   17   11   lham    1 (ASA) 2-center
Cr                     1 (molecules) two-panel
Cr                     2 (ASA) 2-c + pert. corr
Cr                     4 (ASA) auto-down-fold
Cr                     8 (ASA) change rep interactively
Cr                    16 (ASA) suppress d hybridization
Cr                    32 (ASA) preserve ortho. evecs
Cr                    64 (ASA) save evecs to disk
Cr                   128 (ASA) gamma-rep
Cr                   256       use true spherical harmonics
Cr                   512       not used
Cr                  1024       not used
Cr   18   12   lmet    1 metal     2 tetrahedron
Cr                     4 (GF) V-shift1 is zero  8 (GF) V-shift2 is zero
Cr   19   13   lncol   1 noncollinear magnetism
Cr                     2 spin spirals
Cr                     4 spin-orbit coupling
Cr                     8 External magnetic field
Cr                    16 mag. forces
Cr                    32 spin-orbit coupling, LzSz only
Cr                    64 spin-orbit coupling, LzSz + (L.S-LzSz) pert
Cr   20   14   loptc   1 generate Im(eps) (optics package)
Cr                     2 generate Im(eps) w/ on-the-fly sampling integration
Cr                    -1 generate joint DOS
Cr                     Add 10 for SHG
Cr   21   15   lordn   1 Embedded GF  4 Vanderbilt  8 Embedded cluster
Cr   22   16   lpgf    layer GF.  First entry:
Cr                     1: Generate diagonal GF and output density
Cr                     2: Find Fermi left- and right-bulk Ef
Cr                     3: Trace out band structure for spec'd energy mesh
Cr                        In this case make Im(z)=0.
Cr                     4: Ditto for right bulk.
Cr                     5: Calculate current through structure
Cr                        Second entry:
Cr                     0  Make GF using standard approach
Cr                     1  Use LU decomposition
Cr   24   17   lqp     1 do not add inversion 2 inverse iteration
Cr   25   18   lrel    1 scalar relativistic
Cr                     2 fully relativistic
Cr   26   19   lrs     switches concerning restart mode.
Cr                     1 Read from restart file
Cr                     2 Read from restart file, ascii mode
Cr                     4 Read from restart file, invoke smshft
Cr                     8 Write new density to restart file
Cr                    16 Write new density to restart file, ascii format
Cr                    32 read site positions from input file
Cr                    64 read starting fermi level from input file
Cr                   128 read starting pnu level from input file
Cr                   256 rotate local density after reading
Cr   27   20   lscr    0 do nothing
Cr                     1 Make P0(0)
Cr                     2 Screen output q and ves
Cr                     3 Screen output ves only
Cr                     4 Use model response to screen output q
Cr                       Add 1 to combine mode 1 with another mode
Cr                       Add 10*k to compute intra-site contribution to
Cr                       vbare each kth iteration
Cr                       Add 100*k to compute response function only
Cr                       each kth iteration
Cr   28   21   lstonr  second digit for graphical output
Cr   31   22   lstr    Note: no longer used; see str->lshow,str->lequiv
Cr                     1 print strux   2 find equiv strux
Cr   32   23   lsx     1 Calculate screened exchange sigma
Cr                    10's digit nonzero make vintra to include in sigma
Cr   33   24   ltb     1 overlap        2 crystal-field     4 ovlp+CF
Cr                     8 add ebarLL    16 forces           32 fij
Cr                    64 not used     128 pressure        256 evdisc
Cr                   512 pair pot    1024 TrH & local E  2048 local rho
Cr                  2^12 Hubbard U   2^13 No Madelung    2^14 wgt avg U
Cr                  2^15 L>0 estat   2^16 disc read incr 2^17 gamma-pt
Cr   34   25   lves    1 take ves as input
Cr   35   26   lxcf    parameter defining XC functional
Cr                     1s digit:
Cr                     1 for Ceperly-Alder
Cr                     2 for Barth-Hedin (ASW fit)
Cr                     3 for PW91
Cr                     4 for PBE
Cr                     10s digit for GGAs
Cr                     0 for LSDA
Cr                     1 for LMH
Cr                     2 for PW91
Cr                     3 for PBE
Cr                     4 for PBE with Becke exchange
Cr   36   27   maxit   max. no.  iterations in self-consistency cycle
Cr   37   28   mdprm   arg 1: 1 new dynamics  2  restart dynamics
Cr                            4 relax with conjugate gradients
Cr                            5 relax with variable metric
Cr                            6 relax with Broyden
Cr                     arg 2: statics: switch
Cr                            1 read hessian matrix
Cr                            dynamics:
Cr                              number of iterations between printouts.
Cr                     arg 3: (stat) relaxation x-tolerance
Cr                            (dyn)  temperature
Cr                     arg 4: (stat) relaxation g-tolerance
Cr                            (dyn)  time step
Cr                     arg 5: (stat) step length
Cr                            (dyn)  relaxation time
Cr                     arg 6: (stat) Remove hessian after this many steps
Cr                            (dyn)  --
Cr   43   29   modep.. which dimensions are periodic
Cr   46   30   nbas    size of basis
Cr   47   31   nbasp   size of padded basis (layer geometry)
Cr   48   32   nclass  size of class
Cr   49   33   nesabc  mesh for determining empty spheres
Cr   52   34   nitmv   max number of mol-dynamics iterations
Cr   53   35   nl      1 + maximum lmxa
Cr   54   36   nmap    number of maps (ASA)
Cr   55   37   npl     number of principal layers (PGF)
Cr   56   38   nsite   number of sites
Cr   57   39   nspec   number of species
Cr   58   40   nspin   number of spins
Cr   59   41   nvario  number of variables to output
Cr   60   42   omax1   sphere overlap constraints, type 1
Cr   63   43   omax2   sphere overlap constraints, type 2
Cr   66   44   pfloat  stop after:  1 show  2 atom  4 bands
Cr   67   45   quit    stop after:  1 show  2 atom  4 bands
Cr   68   46   rmaxes  upper limit to ES radius when finding new empty spheres
Cr   69   47   rmines  lower limit to ES radius when finding new empty spheres
Cr   70   48   sclwsr  scale wsr until reaching this fractional vol
Cr                     10s digit used for assymetric treatment of ES:
Cr                       0 ES and other sites are treated symmetrically
Cr                       1 all sites with z>0 are resized first; then
Cr                         all sites are resized.
Cr                       2 all sites with z>0 are resized first; then
Cr                         the ES sites only are resized
Cr   71   49   sdmod   spin dynamics mode
Cr   72   50   sdprm   spin dynamics parameters
Cr   77   51   sdxsi   Bulgac and Kusnezov global deamons
Cr   81   52   size    size of this structure
Cr   82   53   smalit  parameters for small iterations
Cr   84   54   tol     1 q- tolerance for self-consistency
Cr                     2 e- tolerance for self-consistency
Cr   87   55   wsrmax  constraint on size of largest WS sphere
Cr   88   56   zbak    background charge and MT correction parameters
Cr                     zbak(1) = uniform background charge included
Cr                               in electrostatics but not in xc pot.
Cr                     zbak(2) = charge used to make MT correction
Cr                               (ASA only)
C ----------------------------------------------------------------
C       data  /2,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,24,25,        26,27,28,31,32,33,34,35,36,37,43,46,47,48,49,52,53,54,55,56,57,        58,59,60,63,66,67,68,69,70,71,72,77,81,82,84,87,88,90/
C       data / 4,4,25*2,4,13*2,4,4,2,2,4,4,4,2,4,4,-1,2,4,4,4/
C  1   2  4       defm  6
C  2   8  4       elin  1
C  3   9  2       lasa  1
C  4  10  2       lbas  1
C  5  11  2        lcd  1
C  6  12  2       lcgf  1
C  7  13  2       ldos  1
C  8  14  2        lfp  1
C  9  15  2      lfrce  1
C 10  16  2      lgen3  1
C 11  17  2       lham  1
C 12  18  2       lmet  1
C 13  19  2      lncol  1
C 14  20  2      loptc  1
C 15  21  2      lordn  1
C 16  22  2       lpgf  2
C 17  24  2        lqp  1
C 18  25  2       lrel  1
C 19  26  2        lrs  1
C 20  27  2       lscr  1
C 21  28  2     lstonr  3
C 22  31  2       lstr  1
C 23  32  2        lsx  1
C 24  33  2        ltb  1
C 25  34  2       lves  1
C 26  35  2       lxcf  1
C 27  36  2      maxit  1
C 28  37  4      mdprm  6
C 29  43  2    modep..  3
C 30  46  2       nbas  1
C 31  47  2      nbasp  1
C 32  48  2     nclass  1
C 33  49  2     nesabc  3
C 34  52  2      nitmv  1
C 35  53  2         nl  1
C 36  54  2       nmap  1
C 37  55  2        npl  1
C 38  56  2      nsite  1
C 39  57  2      nspec  1
C 40  58  2      nspin  1
C 41  59  2     nvario  1
C 42  60  4      omax1  3 possible_structure
C 43  63  4      omax2  3 possible_structure
C 44  66  2     pfloat  1
C 45  67  2       quit  1
C 46  68  4     rmaxes  1
C 47  69  4     rmines  1
C 48  70  4     sclwsr  1
C 49  71  2      sdmod  1
C 50  72  4      sdprm  5
C 51  77  4      sdxsi  4
C 52  81 -1       size  1
C 53  82  2     smalit  2
C 54  84  4        tol  3
C 55  87  4     wsrmax  1
C 56  88  4       zbak  2
      type s_ctrl
      real(8)  ::   size
      real(8)   ::   defm(6)
      real(8)   ::   elin
      integer(8)   ::   lasa
      integer(8)   ::   lbas
      integer(8)   ::   lcd
      integer(8)   ::   lcgf
      integer(8)   ::   ldos
      integer(8)   ::   lfp
      integer(8)   ::   lfrce
      integer(8)   ::   lgen3
      integer(8)   ::   lham
      integer(8)   ::   lmet
      integer(8)   ::   lncol
      integer(8)   ::   loptc
      integer(8)   ::   lordn
      integer(8)   ::   lpgf(2)
      integer(8)   ::   lqp
      integer(8)   ::   lrel
      integer(8)   ::   lrs
      integer(8)   ::   lscr
      integer(8)   ::   lstonr(3)
      integer(8)   ::   lstr
      integer(8)   ::   lsx
      integer(8)   ::   ltb
      integer(8)   ::   lves
      integer(8)   ::   lxcf
      integer(8)   ::   maxit
      real(8)   ::   mdprm(6)
      integer(8)   ::   modep(3)
      integer(8)   ::   nbas
      integer(8)   ::   nbasp
      integer(8)   ::   nclass
      integer(8)   ::   nesabc(3)
      integer(8)   ::   nitmv
      integer(8)   ::   nl
      integer(8)   ::   nmap
      integer(8)   ::   npl
      integer(8)   ::   nsite
      integer(8)   ::   nspec
      integer(8)   ::   nspin
      integer(8)   ::   nvario
      real(8)   ::   omax1(3)
      real(8)   ::   omax2(3)
      integer(8)   ::   pfloat
      integer(8)   ::   quit
      real(8)   ::   rmaxes
      real(8)   ::   rmines
      real(8)   ::   sclwsr
      integer(8)   ::   sdmod
      real(8)   ::   sdprm(5)
      real(8)   ::   sdxsi(4)
      integer(8)   ::   smalit(2)
      real(8)   ::   tol(3)
      real(8)   ::   wsrmax
      real(8)   ::   zbak(2)
      end type s_ctrl

Cr --- explanation of  ugw ---
Cr   off offe  name    purpose
Cr    2    1   delre   spacing on real axis for energy contour
Cr                     delre(1) = dw(0)
Cr                     delre(2) = omegac (coff to parabolic term)
Cr    4    2   deltaw  width in finite diff for sigma energy derivative
Cr    5    3   deltax  Im. E (broadening delta) in computation of x0
Cr    7    4   ecuts   energy cutoff for calculating sigma_nn'
Cr    7    5   gcutb   G-vector cutoff for basis
Cr    8    6   gcutx   G-vector cutoff for response function
Cr    9    7   gsmear  smearing in pole of Green's function
Cr   10    8   lgw     specifies mode for GW computation
Cr                     1s digit
Cr                     0 standard mode
Cr                     1 Faleev's fast energy integration
Cr   11    9   lmxpb   default l-cutoff for product basis
Cr   12   10   mksig   Mode for generating self-energy sigma matrix
Cr                     0 do not create sigma matrix
Cr                     1 Make sigma_nn'(E_f)
Cr                     2 Make sigma_nn'(E_n/2+E_n'/2)
Cr                     3 Make sigma_nn'(E_n)/2 + sigma_nn'(E_n')/2
Cr                     lsig>0 implies 1s digit lgw>0
Cr   13   11   nband   number of bands in perturbation calculation
Cr   14   12   nime    # energies along Im axis for energy contour
Cr   15   13   nkabc   the number of k-points for GW calculation
Cr   18   14   pb1     default product basis 1
Cr   19   15   pb2     default product basis 2
Cr   20   16   pbtol   product-basis overlap cutoff
Cr   21   17   qoffp   qp offset parameter
Cr   22   18   rdsig   parameters for how to read,interpolate SE matrix
Cr   23   19   size    size of this structure
C -----------------------------------------------------------------------
C       data  /2,4,5,6,7,8,9,10,11,12,13,14,15,18,19,20,21,22,23,24/
C       data  /4,4,4,4,4,4,4,2,2,2,2,2,2,1,1,4,4,2,-1/
C  1   2  4      delre  2
C  2   4  4     deltaw  1
C  3   5  4     deltax  1
C  4   6  4      ecuts  1
C  5   7  4      gcutb  1
C  6   8  4      gcutx  1
C  7   9  4     gsmear  1
C  8  10  2        lgw  1
C  9  11  2      lmxpb  1
C 10  12  2      mksig  1
C 11  13  2      nband  1
C 12  14  2       nime  1
C 13  15  2      nkabc  3
C 14  18  1        pb1  1
C 15  19  1        pb2  1
C 16  20  4      pbtol  1
C 17  21  4      qoffp  1
C 18  22  2      rdsig  1
C 19  23 -1       size  1
      type s_gw
      real(8)  ::   size
      real(8)   ::   delre(2)
      real(8)   ::   deltaw
      real(8)   ::   deltax
      real(8)   ::   ecuts
      real(8)   ::   gcutb
      real(8)   ::   gcutx
      real(8)   ::   gsmear
      integer(8)   ::   lgw
      integer(8)   ::   lmxpb
      integer(8)   ::   mksig
      integer(8)   ::   nband
      integer(8)   ::   nime
      integer(8)   ::   nkabc(3)
      real(8)   ::   pb1 !string
      real(8)   ::   pb2 !string
      real(8)   ::   pbtol
      real(8)   ::   qoffp
      integer(8)   ::   rdsig
      end type s_gw

Cr --- explanation of  uham ---
Cr   off offe  name    purpose
Cr    2    1   alfsi  stability factor alfsi*(S+)*S (molecules)
Cr    3    2   amgm   magnetization
Cr    4    3   bandw  Maximum size of off-diagonal, packed storage
Cr    5    4   dabc   Spacing for real-space mesh (molecules)
Cr    8    5   ehf    Harris-Foulkes energy
Cr    9    6   ehk    Hohnberg-Kohn energy
Cr   10    7   elind  Lindhard screening parameter
Cr   11    8   eterms terms making up the total energy.   For FP:
Cr                     1  =  ehf
Cr                     2  =  eks
Cr                     3  =  utot
Cr                     4  =  rhoves int (rhoval x Ves) (true-sm)
Cr                     5  =  cpnves int (core+nucleus x Ves) (true-sm)
Cr                     6  =  rhoeps int (rho * exc)
Cr                     7  =  rhomu  int (rho * vxc)
Cr                     8  =  sumec
Cr                     9  =  sumtc
Cr                     10 =  xcore  rhoc * total potential
Cr                     11 =  valvef int (rhov * vef) (true-sm)
Cr                     12 =  sumt0  sum of foca core energies
Cr                     13 =  sumev  band structure energy
Cr                     13 =  dq1
Cr                     14 =  dq2
Cr                     15 =  amom   system magnetic moment
Cr                     16 =  sumev  sum of single-particle evals
Cr                     17 =  rinvxt input density * external pot
Cr                     18 =  rouvxt output density * external pot
Co                     19 =  rhosig trace of self-energy over occ states
Cr   31    9   hord   order of polynomial approximation to ham.
Cr                    In 2nd gen LMTO, relevant only in GF context.
Cr   32   10   kmto   envelope kinetic energies making up chord LMTO
Cr   38   11   lasa   not used; see ctrl lasa
Cr   39   12   ldham  vector describing hamiltonian dimensions:
Cr                    1: ldim   = dimension of lmto basis
Cr                    2: lidim  = ldim + size of downfolding block
Cr                    3: lidhim = lidim + size of higher block
Cr                    4: nspc   = number of coupled spins
Cr                              = 1 unless noncollinear magnetism
Cr                    5: ldimc  = ldim * nspc
Cr                    6: lidimc = lidim * nspc
Cr                    7: lihdimc = lihdim * nspc
Cr                    8: nspx   = number of separate spin channels
Cr                              = nsp if nspc is not 2; else 1
Cr   55   13   lgen3   shouldn't be used; just a copy of ctrl->lgen3
Cr   56   14   lham    1 (ASA) 2-center
Cr                     1 (molecules) two-panel
Cr                     2 (ASA) 2-c + pert. corr
Cr                     4 (ASA) auto-down-fold
Cr                     8 (ASA) change rep interactively
Cr                    16 (ASA) suppress d hybridization
Cr                    32 (ASA) preserve ortho. evecs
Cr                    64 (ASA) save evecs to disk
Cr                   128 (ASA) gamma-rep
Cr                   256       use true spherical harmonics
Cr                   512       not used
Cr                  1024       not used
Cr                       NB: ham->lham and ctrl->lham should be duplicates
Cr   57   15   lmaxu   dimensioning parameter for LDA+U potential
Cr   58   16   lmxax   largest augmentation lmax in basis
Cr   59   17   lncol   1 noncollinear magnetism
Cr                     2 spin spirals
Cr                     4 spin-orbit coupling
Cr                     8 External magnetic field
Cr                    16 mag. forces
Cr                    32 spin-orbit coupling, LzSz only
Cr                    64 spin-orbit coupling, LzSz + (L.S-LzSz) pert
Cr                    NB: ham->lncol and ctrl->lncol should be duplicates
Cr   60   18   lsig   parameters concerning self-energy
Cr                     1 read sigma, Assume r.s. sigma is real
Cr                     2 read sigma
Cr   61   19   ltb    not used; see ctrl->ltb
Cr   62   20   lxcf   parameter defining local XC functional.
Cr                    Not used: see ctrl->lxcf
Cr   63   21   nbf    number of channels for magnetic field
Cr   64   22   ndham  Largest dimension of lmto+PW basis
Cr   65   23   ndhrs  dimension of site block for r.s. h or sigma
Cr   66   24   ndofH  leading dimension to ooffH
Cr   67   25   neula  number of channels for euler angles
Cr                    0 -> No euler angles defined
Cr   68   26   nkaph  number of repetitions of one l-quant.n. in basis
Cr   69   27   nlibu  Number of LDA+U blocks
Cr   70   28   nmto   Polynomial order to NMTO
Cr   71   29   npwmin Estimate for minimum number of PWs to be calculated
Cr                    (PW basis is q-dependent; max size not known a priori)
Cr   72   30   npwpad Padding to be added to estimated max basis dimension
Cr                    and subtracted from min basis dimension
Cr                    (PW basis is q-dependent; max size not known a priori)
Cr   73   31   nqsig  Number of k-points for which self-energy known
Cr                    (Used for interpolating between qpoints)
Cr   74   32   obdots mag. field in local coordinates, dotted with
Cr                    Paulie matrices (pointer in w array)
Cr   75   33   oeula  Euler angles (pointer in w array)
Cr   76   34   ohrs   r.s. h or sigma (pointer in w array)
Cr   77   35   oiaxs  neighbor table (pointer in w array)
Cr   78   36   oindxo orbital permutation table (pointer in w array)
Cr   79   37   olmxa  vector of augmentation l-cutoffs (pointer in w array)
Cr   80   38   omagf  magnetic field (pointer in w array)
Cr   81   39   ontabs table of number-of-pairs per site (pointer in w array)
Cr   82   40   ooffH  table of offsets to hamiltonian matrix (pointer in w array)
Cr                    (pointer in w array)
Cr                    (Used for interpolating between qpoints)
Cr   83   41   oqsig  list of qp at which sigma can be computed
Cr   84   42   oveps  When diagonalizing hamiltonian, discard part of hibert space
Cr                    corresponding to evals of overlap < oveps
Cr   85   43   pmax   global minimum allowed values for pnu
Cr   95   44   pmin   global minimum allowed values for pnu
Cr  105   45   pwemax High Energy cutoff for PW part of basis
Cr  106   46   pwemin Low Energy cutoff for PW part of basis
Cr  107   47   pwmode Controls PW part of basis
Cr                    0 => no PW part of basis
Cr                    1 => include PWs in basis
Cr                    2 => include only PWs in basis
Cr  108   48   qpoff  qp offset when generating transformation of self-energy
Cr  111   49   qss    spin spiral
Cr  115   50   rsrnge cutoff in connecting vector length for r.s. sigma
Cr  116   51   rsstol tolerance in Bloch sum error for r.s. sigma
Cr  117   52   seref  Sum of reference energies
Cr  118   53   sigp   parameters for approximating self-energy sigma
Cr                    arg 1: mode : specifies how to set its diagonal part
Cr                           for states above the high-energy cutoff
Cr                           0 constrain sigii to be > asig+bsig*e
Cr                           1 constrain sigii to be = asig+bsig*e
Cr                           2 constrain sigii to be > asig and < bsig
Cr                           3 constraint same as mode 1.
Cr                             Mode 3 differs in that the least-squares fit to
Cr                             sigii (for informational purposes only, to help
Cr                             estimate asig and bsig) is done for states between
Cr                             efit and nmax or emax
Cr                    arg 2: nmin : sigma for states 1..nmin are approximated by sigii
Cr                    arg 3: emin : (used only if nmin<0)
Cr                                : sigma for levels e<emin are approximated by sigii
Cr                    arg 4: nmax : sigma for levels i>nmax are approximated by
Cr                                  sigii AND constrained according to mode
Cr                    arg 5: emax : (used only if nmax<=0)
Cr                                : sigma for levels e<emax are approximated by
Cr                                  sigii AND constrained according to mode
Cr                    arg 6: asig : constraint used to approximate
Cr                                  sigii = asig + E * bsig  or
Cr                                  asig < sigii < bsig
Cr                    arg 7: bsig : constraint used to approximate
Cr                                  sigii = asig + E * bsig  or
Cr                                  asig < sigii < bsig
Cr                    arg 8: efit : (mode 3) energy minimium (not used?)
Cr                                  for fitting asig and bsig
Cr                    arg 9: ? lwrite:T, write sigii to file (not used?)
Cr                    10 not used now
Cr  128   54   size   size of this structure
Cr  129   55   thrpv  3 PV for cell
Cr  130   56   udiag  diagonal-only LDA+U
C ----------------------------------------------------------------
C       data  /2,3,4,5,8,9,10,11,31,32,38,39,55,56,57,58,59,60,61,        62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,        83,84,85,95,105,106,107,108,111,115,116,117,118,128,129,130,131/
C       data /4,4,2,4,4,4,4,4,2,4,31*2,5*4,2,6*4,-1,4,2/
C  1   2  4      alfsi  1
C  2   3  4       amgm  1
C  3   4  2      bandw  1
C  4   5  4       dabc  3
C  5   8  4        ehf  1
C  6   9  4        ehk  1
C  7  10  4      elind  1
C  8  11  4     eterms 20
C  9  31  2       hord  1
C 10  32  4       kmto  6
C 11  38  2       lasa  1
C 12  39  2      ldham 16
C 13  55  2      lgen3  1
C 14  56  2       lham  1
C 15  57  2      lmaxu  1
C 16  58  2      lmxax  1
C 17  59  2      lncol  1
C 18  60  2       lsig  1
C 19  61  2        ltb  1
C 20  62  2       lxcf  1
C 21  63  2        nbf  1
C 22  64  2      ndham  1
C 23  65  2      ndhrs  1
C 24  66  2      ndofH  1
C 25  67  2      neula  1
C 26  68  2      nkaph  1
C 27  69  2      nlibu  1
C 28  70  2       nmto  1
C 29  71  2     npwmin  1
C 30  72  2     npwpad  1
C 31  73  2      nqsig  1
C 32  74  2     obdots  1 possible_structure
C 33  75  2      oeula  1 possible_structure
C 34  76  2       ohrs  1 possible_structure
C 35  77  2      oiaxs  1 possible_structure
C 36  78  2     oindxo  1 possible_structure
C 37  79  2      olmxa  1 possible_structure
C 38  80  2      omagf  1 possible_structure
C 39  81  2     ontabs  1 possible_structure
C 40  82  2      ooffH  1 possible_structure
C 41  83  2      oqsig  1 possible_structure
C 42  84  4      oveps  1 possible_structure
C 43  85  4       pmax 10
C 44  95  4       pmin 10
C 45 105  4     pwemax  1
C 46 106  4     pwemin  1
C 47 107  2     pwmode  1
C 48 108  4      qpoff  3
C 49 111  4        qss  4
C 50 115  4     rsrnge  1
C 51 116  4     rsstol  1
C 52 117  4      seref  1
C 53 118  4       sigp 10
C 54 128 -1       size  1
C 55 129  4      thrpv  1
C 56 130  2      udiag  1
      type s_ham
      real(8)  ::   size
      real(8)   ::   alfsi
      real(8)   ::   amgm
      integer(8)   ::   bandw
      real(8)   ::   dabc(3)
      real(8)   ::   ehf
      real(8)   ::   ehk
      real(8)   ::   elind
      real(8)   ::   eterms(20)
      integer(8)   ::   hord
      real(8)   ::   kmto(6)
      integer(8)   ::   lasa
      integer(8)   ::   ldham(16)
      integer(8)   ::   lgen3
      integer(8)   ::   lham
      integer(8)   ::   lmaxu
      integer(8)   ::   lmxax
      integer(8)   ::   lncol
      integer(8)   ::   lsig
      integer(8)   ::   ltb
      integer(8)   ::   lxcf
      integer(8)   ::   nbf
      integer(8)   ::   ndham
      integer(8)   ::   ndhrs
      integer(8)   ::   ndofH
      integer(8)   ::   neula
      integer(8)   ::   nkaph
      integer(8)   ::   nlibu
      integer(8)   ::   nmto
      integer(8)   ::   npwmin
      integer(8)   ::   npwpad
      integer(8)   ::   nqsig
CCHSTRUC       integer(8)   ::   obdots ! not used because of undef
CCHSTRUC       integer(8)   ::   oeula ! not used because of undef
      integer(8)   ::   ohrs
      integer(8)   ::   oiaxs
CCHSTRUC       integer(8)   ::   oindxo
       integer , pointer ::  iv_p_oindxo (:)
CCHSTRUC       integer(8)   ::   olmxa ! not used because of undef
CCHSTRUC       integer(8)   ::   omagf ! not used because of undef
      integer(8)   ::   ontabs
CCHSTRUC       integer(8)   ::   ooffH
       integer , pointer ::  iv_p_ooffH (:)
CCHSTRUC       integer(8)   ::   oqsig
       real(8) , pointer ::  rv_p_oqsig (:)
      real(8)   ::   oveps
      real(8)   ::   pmax(10)
      real(8)   ::   pmin(10)
      real(8)   ::   pwemax
      real(8)   ::   pwemin
      integer(8)   ::   pwmode
      real(8)   ::   qpoff(3)
      real(8)   ::   qss(4)
      real(8)   ::   rsrnge
      real(8)   ::   rsstol
      real(8)   ::   seref
      real(8)   ::   sigp(10)
      real(8)   ::   thrpv
      integer(8)   ::   udiag
      end type s_ham

Cr --- explanation of  ulat ---
Cr   off offe  name    purpose
Cr    2    1    alat    lattice parameter, in a.u.
Cr    3    2    as      dimensionless Ewald smoothing parameter
Cr    4    3    avw     average MT radius
Cr    5    4    awald   Ewald smoothing parameter
Cr    6    5    dist    deformation parameters (cf lattdf.f)
Cr   15    6    gam     lattice shear parms: gam, gx,gy,gz
Cr   19    7    gmax    cutoff gmax for Fourier transforms
Cr   20    8    ldist   switch specifying what kind of dist (cf lattdf.f)
Cr   21    9    nabc    no. divisions for F.T. mesh
Cr   24   10    ng      no. G vectors
Cr   25   11    nkd     no. direct latt. vecs. for Ewald sum
Cr   26   12    nkdmx   dimensioning for arrays holding latt. vecs
Cr   27   13    nkq     no. reciprocal latt. vecs. for Ewald sum
Cr   28   14    nkqmx   dimensioning for arrays holding latt. vecs
Cr   29   15    npgrp   Number of point symmetry group operations
Cr   30   16    nsgrp   Number of space symmetry group operations
Cr   31   17    oag     offset to symmetry group translations
Cr   32   18    obgv    phase factor sum for symmetrization of mesh rho
Cr   33   19    ocg     offset to Clebsch Gordan coeffs
Cr   34   20    ocy     offset to Ylm normalization constants
Cr   35   21    odlv    offset to direct lattice vector
Cr   36   22    ogv     offset to list of F.T. G vectors
Cr   37   23    oidxcg  offset to Clebsch Gordan indxcg
Cr   38   24    oips0   pointer to first vec in star (for symm mesh rho)
Cr   39   25    oistab  offset to site permutations table for group ops
Cr   40   26    ojcg    offset to Clebsch Gordan jcg
Cr   41   27    okv     offset to indices in list of F.T. G vectors
Cr   42   28    opos    offset to site positions
Cr   43   29    oqlv    offset to Ewald reciprocal lattice vectors
Cr   44   30    osymgr  offset to symmetry group rotation matrices
Cr   45   31    plat..  lattice vectors, units of alat
Cr   54   32    plat0.. lattice vectors before distortion
Cr   63   33    plat2.. secondary lattice vecs used in various contexts
Cr   72   34    plate.. order-N
Cr   81   35    platl.. pgf
Cr   90   36    platr.. pgf
Cr   99   37    qlat..  reciprocal lattice vectors, units 2pi/a
Cr  108   38    rpad..  truncate Ewald to rpad*rmax when lattice vector
Cr                      list has to be padded in order to include at
Cr                      least one lattice vector
Cr  109   39    size    size of this structure
Cr  110   40    slat    superlattice vectors
Cr  119   41    tol     Ewald tolerance
Cr  120   42    tolft   FT mesh tolerance
Cr  121   43    vol     cell volume
C ----------------------------------------------------------------
C       data  /2,3,4,5,6,15,19,20,21,24,25,26,27,28,29,30,31,32,33,        34,35,36,37,38,39,40,41,42,43,44,45,54,63,72,81,90,99,108,109,        110,119,120,121,122/
C       data / 4,4,4,4,4,4,4,23*2,4,4,4,4,4,4,4,4,-1,4,4,4,4/
C  1   2  4       alat  1
C  2   3  4         as  1
C  3   4  4        avw  1
C  4   5  4      awald  1
C  5   6  4       dist  9
C  6  15  4        gam  4
C  7  19  4       gmax  1
C  8  20  2      ldist  1
C  9  21  2       nabc  3
C 10  24  2         ng  1
C 11  25  2        nkd  1
C 12  26  2      nkdmx  1
C 13  27  2        nkq  1
C 14  28  2      nkqmx  1
C 15  29  2      npgrp  1
C 16  30  2      nsgrp  1
C 17  31  2        oag  1 possible_structure
C 18  32  2       obgv  1 possible_structure
C 19  33  2        ocg  1 possible_structure
C 20  34  2        ocy  1 possible_structure
C 21  35  2       odlv  1 possible_structure
C 22  36  2        ogv  1 possible_structure
C 23  37  2     oidxcg  1 possible_structure
C 24  38  2      oips0  1 possible_structure
C 25  39  2     oistab  1 possible_structure
C 26  40  2       ojcg  1 possible_structure
C 27  41  2        okv  1 possible_structure
C 28  42  2       opos  1 possible_structure
C 29  43  2       oqlv  1 possible_structure
C 30  44  2     osymgr  1 possible_structure
C 31  45  4     plat..  9
C 32  54  4    plat0..  9
C 33  63  4    plat2..  9
C 34  72  4    plate..  9
C 35  81  4    platl..  9
C 36  90  4    platr..  9
C 37  99  4     qlat..  9
C 38 108  4     rpad..  1
C 39 109 -1       size  1
C 40 110  4       slat  9
C 41 119  4        tol  1
C 42 120  4      tolft  1
C 43 121  4        vol  1
      type s_lat
      real(8)  ::   size
      real(8)   ::   alat
      real(8)   ::   as
      real(8)   ::   avw
      real(8)   ::   awald
      real(8)   ::   dist(9)
      real(8)   ::   gam(4)
      real(8)   ::   gmax
      integer(8)   ::   ldist
      integer(8)   ::   nabc(3)
      integer(8)   ::   ng
      integer(8)   ::   nkd
      integer(8)   ::   nkdmx
      integer(8)   ::   nkq
      integer(8)   ::   nkqmx
      integer(8)   ::   npgrp
      integer(8)   ::   nsgrp
      integer(8)   ::   oag
      integer(8)   ::   obgv
      integer(8)   ::   ocg
      integer(8)   ::   ocy
      integer(8)   ::   odlv
      integer(8)   ::   ogv
      integer(8)   ::   oidxcg
CCHSTRUC       integer(8)   ::   oips0
       integer , pointer ::  iv_p_oips0 (:)
CCHSTRUC       integer(8)   ::   oistab
       integer , pointer ::  iv_p_oistab (:)
      integer(8)   ::   ojcg
      integer(8)   ::   okv
      integer(8)   ::   opos
CCHSTRUC       integer(8)   ::   oqlv
       real(8) , pointer ::  rv_p_oqlv (:)
      integer(8)   ::   osymgr
      real(8)   ::   plat(9)
      real(8)   ::   plat0(9)
      real(8)   ::   plat2(9)
      real(8)   ::   plate(9)
      real(8)   ::   platl(9)
      real(8)   ::   platr(9)
      real(8)   ::   qlat(9)
      real(8)   ::   rpad
      real(8)   ::   slat(9)
      real(8)   ::   tol
      real(8)   ::   tolft
      real(8)   ::   vol
      end type s_lat

Cr --- explanation of  umix ---
Cr   off offe  name    purpose
Cr    2    1   b       mixing beta
Cr    3    2   bl      previous mixing beta
Cr    4    3   bv      extra potential mixing
C     5    4   elind   Lindhard energy for model screening
Cr    6    5   fn      mixing file name
Cr    7    6   kill    kill the mixing file after k iterations
Cr    8    7   lxpot   decouple potential and the charge
Cr                     1: mix vin and v(qmix)
Cr                     2: mix vin and v(qout)
Cr    9    8   mmix    maximum number to mix
Cr   10    9   mode    1 Anderson 2 Broyden
Cr   11   10   model   previous mixing mode
Cr   12   11   n       Number of iterations for this species
Cr   13   12   nitu    max number of LDA+U itreations
Cr   14   13   nmix    actual number mixed
Cr   15   14   nsave   number of iterations to save on disk
Cr   16   15   r..     expression for rmscst
Cr   19   16   rms1    1st rms error
Cr   20   17   rms2    2nd rms error
Cr   21   18   size    total size of this struc
Cr   22   19   tj..    Anderson t's
Cr   32   20   tolu    tolerance for LDA+U
Cr   33   21   umix    mixing parameter for LDA+U
Cr   34   22   w..     Linear mixing weights
Cr   37   23   wc      Broyden weight
C ----------------------------------------------------------------
C       data  /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,19,20,21,22,32,        33,34,37,38/
C       data  /4,4,4,4,1,2,2,2,2,2,2,2,2,2,1,4,4,-1,4,4,4,4,4/
C  1   2  4          b  1
C  2   3  4         bl  1
C  3   4  4         bv  1
C  4   5  4      elind  1
C  5   6  1         fn  1
C  6   7  2       kill  1
C  7   8  2      lxpot  1
C  8   9  2       mmix  1
C  9  10  2       mode  1
C 10  11  2      model  1
C 11  12  2          n  1
C 12  13  2       nitu  1
C 13  14  2       nmix  1
C 14  15  2      nsave  1
C 15  16  1        r..  3
C 16  19  4       rms1  1
C 17  20  4       rms2  1
C 18  21 -1       size  1
C 19  22  4       tj.. 10
C 20  32  4       tolu  1
C 21  33  4       umix  1
C 22  34  4        w..  3
C 23  37  4         wc  1
      type s_mix
      real(8)  ::   size
      real(8)   ::   b
      real(8)   ::   bl
      real(8)   ::   bv
      real(8)   ::   elind
      real(8)   ::   fn !string
      integer(8)   ::   kill
      integer(8)   ::   lxpot
      integer(8)   ::   mmix
      integer(8)   ::   mode
      integer(8)   ::   model
      integer(8)   ::   n
      integer(8)   ::   nitu
      integer(8)   ::   nmix
      integer(8)   ::   nsave
      real(8)   ::   r(3) !string
      real(8)   ::   rms1
      real(8)   ::   rms2
      real(8)   ::   tj(10)
      real(8)   ::   tolu
      real(8)   ::   umix
      real(8)   ::   w(3)
      real(8)   ::   wc
      end type s_mix

Cr --- explanation of  umove ---
Cr   off offe  name    purpose
Cr    2    1   ct      coefficients for global deamons thermostat
Cr    5    2   gyro    gyromagnetic ratio (magnetic dynamics)
Cr    6    3   kt      temperature, atomic units
Cr    7    4   modt    mode for thermostat
Cr   10    5   nmodt   number of thermostat modes
Cr   11    6   prmint  Parameters for numerical integration
Cr                     For Bulirsch-Stoer integration:
Cr                     1:   mode: (1 for BS-integration)
Cr                     2:   rst:  1 for start new integration, 0 to continue
Cr                     3:   ts0:  minimum time step size
Cr                     4:   tol:  tolerance in integration errors
Cr                     5:   mx:   order of rational function extrapolation
Cr                     6:   mi:   number of midpoint rules
Cr                     7-17 nseq: sequence of no. midpoint divisions
Cr                     18:        offset to bs workspace
Cr   31    7   size
Cr   32    8   tnow    duration of simulation to this point, in units of ts
Cr   33    9   ts      suggested time step size
Cr   34   10   tsequ   initial equilibration time before stats accumulated
Cr   35   11   tstot   duration of total simulation, in units of ts
C ----------------------------------------------------------------
C       data  /2,5,6,7,10,11,31,32,33,34,35,36/
C       data / 4,4,4,2,2,4,-1,4,4,4,4 /
C  1   2  4         ct  3
C  2   5  4       gyro  1
C  3   6  4         kt  1
C  4   7  2       modt  3
C  5  10  2      nmodt  1
C  6  11  4     prmint 20
C  7  31 -1       size  1
C  8  32  4       tnow  1
C  9  33  4         ts  1
C 10  34  4      tsequ  1
C 11  35  4      tstot  1
      type s_move
      real(8)  ::   size
      real(8)   ::   ct(3)
      real(8)   ::   gyro
      real(8)   ::   kt
      integer(8)   ::   modt(3)
      integer(8)   ::   nmodt
      real(8)   ::   prmint(20)
      real(8)   ::   tnow
      real(8)   ::   ts
      real(8)   ::   tsequ
      real(8)   ::   tstot
      end type s_move

Cr --- explanation of  uoptic ---
Cr   off offe  name    purpose
Cr    2    1   axes   (abc) axes in xi^(abc) for nonlinear optics
Cr   20    2   cll    core level l quantum number
Cr   21    3   cln    core level n quantum number
Cr   22    4   cls    1 or 2 for EELS or XANES
Cr   23    5   clsite site id for core level
Cr   24    6   esciss Energy shift (scissors operator)
Cr   25    7   mode   1 calculate band-to-band contribution
Cr   26    8   nchi2  number of (abc) for which to calculate xi^(abc)
Cr   27    9   ne     number of energy points
Cr   28   10   ocrng  range of occupied bands
Cr   30   11   size   size of this structure
Cr   31   12   unrng  range of unoccupied bands
Cr   33   13   window energy window for epsilon
Cr   35   ----------- points past end of structure
C       data  /2,20,21,22,23,24,25,26,27,28,30,31,33,35/
C       data / 2,2,2,2,2,4,2,2,2,2,-1,2,4/
C  1   2  2       axes 18
C  2  20  2        cll  1
C  3  21  2        cln  1
C  4  22  2        cls  1
C  5  23  2     clsite  1
C  6  24  4     esciss  1
C  7  25  2       mode  1
C  8  26  2      nchi2  1
C  9  27  2         ne  1
C 10  28  2      ocrng  2 possible_structure
C 11  30 -1       size  1
C 12  31  2      unrng  2
C 13  33  4     window  2
      type s_optic
      real(8)  ::   size
      integer(8)   ::   axes(18)
      integer(8)   ::   cll
      integer(8)   ::   cln
      integer(8)   ::   cls
      integer(8)   ::   clsite
      real(8)   ::   esciss
      integer(8)   ::   mode
      integer(8)   ::   nchi2
      integer(8)   ::   ne
      integer(8)   ::   ocrng(2)
      integer(8)   ::   unrng(2)
      real(8)   ::   window(2)
      end type s_optic

Cr --- explanation of  uordn ---
Cr   off offe  name    purpose
Cr    2    1   efre    effective energy at which free electrons
Cr                     start propagating
Cr    3    2   mode    order-N mode
Cr    4    3   ncl     Number of clusters
Cr    5    4   ndofh   leading dimension to offcH
Cr    6    5   oclp
Cr    7    6   oclssl
Cr    8    7   oiaxg   neighbor table for Green's function
Cr    9    8   omapgv  maps GF pairs into iax table
Cr   10    9   ooffch  offset to table of Hamiltonian offsets
Cr   11   10   rmaxg   ranges of Green's function:
Cr                     (1) = range for direct inversion
Cr                     (2) = range for iterative inversion
Cr                     (3) = range for free-electron GF
Cr   14   11   size    size of this structure
C ----------------------------------------------------------------
C       data  /2,3,4,5,6,7,8,9,10,11,14,15/
C       data / 4,2,2,2,2,2,2,2,2,4,-1/
C  1   2  4       efre  1
C  2   3  2       mode  1
C  3   4  2        ncl  1
C  4   5  2      ndofh  1
C  5   6  2       oclp  1 possible_structure
C  6   7  2     oclssl  1 possible_structure
C  7   8  2      oiaxg  1 possible_structure
C  8   9  2     omapgv  1 possible_structure
C  9  10  2     ooffch  1 possible_structure
C 10  11  4      rmaxg  3
C 11  14 -1       size  1
      type s_ordn
      real(8)  ::   size
      real(8)   ::   efre
      integer(8)   ::   mode
      integer(8)   ::   ncl
      integer(8)   ::   ndofh
CCHSTRUC       integer(8)   ::   oclp ! not used because of undef
CCHSTRUC       integer(8)   ::   oclssl ! not used because of undef
CCHSTRUC       integer(8)   ::   oiaxg ! not used because of undef
CCHSTRUC       integer(8)   ::   omapgv ! not used because of undef
CCHSTRUC       integer(8)   ::   ooffch ! not used because of undef
      real(8)   ::   rmaxg(3)
      end type s_ordn

Cr --- explanation of  upot ---
Cr   off offe  name    purpose
Cr    2    1   bfield  global magnetic field direction and amplitude (4)
Cr    6    2   nlma    number of augmentation channels in system
Cr    7    3   nlml    (FP) total number of site density channels
Cr    8    4   nrhos   (FP) total number of site density channels
Cr    9    5   oaamom  offset to array containing local ASA mag.mom
Cr   10    6   obxc    offset to array containing local ASA Bxc dir.
Cr   11    7   odddpf  3rd derivative of potential functions (mkptfp)
Cr   12    8   oddpf   2nd derivative of potential functions (mkptfp)
Cr   13    9   oddpfr   2nd derivative of potential functions (mkptfp),
Cr                     fully relativistic case
Cr   14   10   odel    (tbe) potential shifts.
Cr   15   11   odpf    energy derivative of potential functions (mkptfp)
Cr   16   12   odpfr   energy derivative of potential functions (mkptfp),
Cr                     fully relativistic case
Cr   17   13   ofes    electrostatic contribution to forces
Cr   18   14   ogma    gamma-alpha in potential functions
Cr   19   15   ogmar   (gamma-alpha)*P^alpha/P^gamma fully relativistic
Cr   20   16   ogrrme  radial matrix elements of gradient.
Cr   21   17   ohab    not used now
Cr   22   18   oivso   (tbe) spin orbit
Cr   23   19   omad    Madelung matrix
Cr   24   20   oorhat  pointers to atom densities
Cr   25   21   opalp   ASA potential functions, alpha repsn (mkptfp)
Cr   26   22   opapg   p^alpha/p^gamma fully relativistic
Cr   27   23   opdel   ASA Delta, in hamiltonian order (for LDA+U)
Cr   28   24   opf     ASA potential functions (mkptfp)
Cr   29   25   opfnc   noncollinear ASA potential functions
Cr   30   26   opfr    ASA potential functions (mkptfp),
Cr                     fully relativistic case
Cr   31   27   opmpol  multipole integrals of phi,phidot
Cr   32   28   opnu    P-nu (Methfessel's log derivative function)
Cr   33   29   opp     potential parameters
Cr   34   30   oppi    not used now
Cr   35   31   oppn    NMTO generation potential parameters
Cr   36   32   opprel  Dirac potential parameters
Cr   37   33   opti    inverse potential functions
Cr   38   34   oqc     Sphere core charge
Cr   39   35   oqmom   offset to multipole moments
Cr   40   36   oqnu    Energy moments of the charge density
Cr   41   37   oqpp    Multipole moments of the nonspherical charge
Cr   42   38   oqt     Sphere total charge
Cr   43   39   orhos   spin-density matrix
Cr   44   40   orhrmx  sphere electron density at rmax
Cr   45   41   osab    not used now
Cr   46   42   osgw    structure containing GW parameters
Cr   47   43   osmpot  offset to smoothed potential
Cr   48   44   osmrho  offset to smoothed density
Cr   49   45   osop    spin-orbit parameters
Cr   50   46   osoptc  structure containing optical matrix elements
Cr   51   47   osrout  offset to smoothed output density
Cr   52   48   otau    not used now
Cr   53   49   ovab    not used now
Cr   54   50   ovdif   difference ves(q,rmax) - vactual(rmax)
Cr   55   51   oves    electrostatic potential at rmax
Cr   56   52   ovintr  intra-atomic W for screening, monopole approx
Cr   57   53   ovrmax  total potential at rmax
Cr   58   54   ovshf   ASA constant potential shifts
Cr   59   55   size    size of this structure
Cr   60   56   vconst  Constant estat potential shifts, used where
Cr                     the Fermi level is specified and the potential
Cr                     adjusts to it.
Cr                     vconst(1) = potential shift
Cr                     vconst(2) = potential shift of L end region (PGF)
Cr                     vconst(3) = potential shift of R end region (PGF)
Cr   63   57   vmtz    ASA muffin-tin zero
Cr   64   58   vmtz0   Input ASA muffin-tin zero
C -----------------------------------------------------------------
C       data  /2,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,        24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,        45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,63,64,65/
C       data / 4,53*2,-1,3*4/
C  1   2  4     bfield  4
C  2   6  2       nlma  1
C  3   7  2       nlml  1
C  4   8  2      nrhos  1
C  5   9  2     oaamom  1 possible_structure
C  6  10  2       obxc  1 possible_structure
C  7  11  2     odddpf  1 possible_structure
C  8  12  2      oddpf  1 possible_structure
C  9  13  2     oddpfr  1 possible_structure
C 10  14  2       odel  1 possible_structure
C 11  15  2       odpf  1 possible_structure
C 12  16  2      odpfr  1 possible_structure
C 13  17  2       ofes  1 possible_structure
C 14  18  2       ogma  1 possible_structure
C 15  19  2      ogmar  1 possible_structure
C 16  20  2     ogrrme  1 possible_structure
C 17  21  2       ohab  1 possible_structure
C 18  22  2      oivso  1 possible_structure
C 19  23  2       omad  1 possible_structure
C 20  24  2     oorhat  1 possible_structure
C 21  25  2      opalp  1 possible_structure
C 22  26  2      opapg  1 possible_structure
C 23  27  2      opdel  1 possible_structure
C 24  28  2        opf  1 possible_structure
C 25  29  2      opfnc  1 possible_structure
C 26  30  2       opfr  1 possible_structure
C 27  31  2     opmpol  1 possible_structure
C 28  32  2       opnu  1 possible_structure
C 29  33  2        opp  1 possible_structure
C 30  34  2       oppi  1 possible_structure
C 31  35  2       oppn  1 possible_structure
C 32  36  2     opprel  1 possible_structure
C 33  37  2       opti  1 possible_structure
C 34  38  2        oqc  1 possible_structure
C 35  39  2      oqmom  1 possible_structure
C 36  40  2       oqnu  1 possible_structure
C 37  41  2       oqpp  1 possible_structure
C 38  42  2        oqt  1 possible_structure
C 39  43  2      orhos  1 possible_structure
C 40  44  2     orhrmx  1 possible_structure
C 41  45  2       osab  1 possible_structure
C 42  46  2       osgw  1 possible_structure
C 43  47  2     osmpot  1 possible_structure
C 44  48  2     osmrho  1 possible_structure
C 45  49  2       osop  1 possible_structure
C 46  50  2     osoptc  1 possible_structure
C 47  51  2     osrout  1 possible_structure
C 48  52  2       otau  1 possible_structure
C 49  53  2       ovab  1 possible_structure
C 50  54  2      ovdif  1 possible_structure
C 51  55  2       oves  1 possible_structure
C 52  56  2     ovintr  1 possible_structure
C 53  57  2     ovrmax  1 possible_structure
C 54  58  2      ovshf  1 possible_structure
C 55  59 -1       size  1
C 56  60  4     vconst  3
C 57  63  4       vmtz  1
C 58  64  4      vmtz0  1
      type s_pot
      real(8)  ::   size
      real(8)   ::   bfield(4)
      integer(8)   ::   nlma
      integer(8)   ::   nlml
      integer(8)   ::   nrhos
CCHSTRUC       integer(8)   ::   oaamom ! not used because of undef
CCHSTRUC       integer(8)   ::   obxc
       real(8) , pointer ::  rv_p_obxc (:)
CCHSTRUC       integer(8)   ::   odddpf ! not used because of undef
CCHSTRUC       integer(8)   ::   oddpf ! not used because of undef
CCHSTRUC       integer(8)   ::   oddpfr ! not used because of undef
CCHSTRUC       integer(8)   ::   odel ! not used because of undef
CCHSTRUC       integer(8)   ::   odpf ! not used because of undef
CCHSTRUC       integer(8)   ::   odpfr ! not used because of undef
CCHSTRUC       integer(8)   ::   ofes ! not used because of undef
CCHSTRUC       integer(8)   ::   ogma ! not used because of undef
CCHSTRUC       integer(8)   ::   ogmar ! not used because of undef
CCHSTRUC       integer(8)   ::   ogrrme
       real(8) , pointer ::  rv_p_ogrrme (:)
CCHSTRUC       integer(8)   ::   ohab ! not used because of undef
CCHSTRUC       integer(8)   ::   oivso ! not used because of undef
CCHSTRUC       integer(8)   ::   omad
       real(8) , pointer ::  rv_p_omad (:)
CCHSTRUC       integer(8)   ::   oorhat
       integer , pointer ::  iv_p_oorhat (:)
CCHSTRUC       integer(8)   ::   opalp ! not used because of undef
CCHSTRUC       integer(8)   ::   opapg ! not used because of undef
CCHSTRUC       integer(8)   ::   opdel ! not used because of undef
CCHSTRUC       integer(8)   ::   opf ! not used because of undef
CCHSTRUC       integer(8)   ::   opfnc ! not used because of undef
CCHSTRUC       integer(8)   ::   opfr ! not used because of undef
CCHSTRUC       integer(8)   ::   opmpol
       real(8) , pointer ::  rv_p_opmpol (:)
CCHSTRUC       integer(8)   ::   opnu
       real(8) , pointer ::  rv_p_opnu (:)
CCHSTRUC       integer(8)   ::   opp
       real(8) , pointer ::  rv_p_opp (:)
CCHSTRUC       integer(8)   ::   oppi ! not used because of undef
CCHSTRUC       integer(8)   ::   oppn ! not used because of undef
      integer(8)   ::   opprel
CCHSTRUC       integer(8)   ::   opti ! not used because of undef
      integer(8)   ::   oqc
CCHSTRUC       integer(8)   ::   oqmom ! not used because of undef
      integer(8)   ::   oqnu
      integer(8)   ::   oqpp
      integer(8)   ::   oqt
      integer(8)   ::   orhos
      integer(8)   ::   orhrmx
CCHSTRUC       integer(8)   ::   osab ! not used because of undef
CCHSTRUC       integer(8)   ::   osgw ! not used because of undef
CCHSTRUC       integer(8)   ::   osmpot
       complex(8) , pointer ::  zv_p_osmpot (:)
CCHSTRUC       integer(8)   ::   osmrho
       complex(8) , pointer ::  zv_p_osmrho (:)
      integer(8)   ::   osop
CCHSTRUC       integer(8)   ::   osoptc ! not used because of undef
CCHSTRUC       integer(8)   ::   osrout ! not used because of undef
CCHSTRUC       integer(8)   ::   otau ! not used because of undef
CCHSTRUC       integer(8)   ::   ovab ! not used because of undef
CCHSTRUC       integer(8)   ::   ovdif
       real(8) , pointer ::  rv_p_ovdif (:)
CCHSTRUC       integer(8)   ::   oves
       real(8) , pointer ::  rv_p_oves (:)
CCHSTRUC       integer(8)   ::   ovintr
       real(8) , pointer ::  rv_p_ovintr (:)
CCHSTRUC       integer(8)   ::   ovrmax
       real(8) , pointer ::  rv_p_ovrmax (:)
CCHSTRUC       integer(8)   ::   ovshf ! not used because of undef
      real(8)   ::   vconst(3)
      real(8)   ::   vmtz
      real(8)   ::   vmtz0
      end type s_pot

Cr --- explanation of  uspec ---
Cr   off offe  name    purpose
Cr   20     1  a      a for mesh
Cr  111..   2  alpha  screening alphas
Cr   24     3  beta   mixing beta
Cr  211..   4  chfa   coefficients to fit of free-atom density tails
Cr   46     5  colxbs color of atom (for graphics)
Cr   17     6  coreh  core hole channel (char: eg '1s')
Cr   18     7  coreq  coreq(1) = charge in core hole channel
Cr                    coreq(2) = moment in core hole channel (nsp=2 only)
Cr   26     8  ctail  coefficients to fit of free-atom core tail by unsm. Hankel
Cr   37     9  dv     constant shift in potential
Cr   53    10  eh3    sm Hankel energy for high local orbitals
Cr  341    11  ehvl   l-dependent sm Hankel energies for val-lap basis
Cr  291    12  enu    linearization energies
Cr   22    13  eref   reference energy
Cr   27    14  etail  energy to fit of free-atom core tail
Cr   23    15  etf    Thomas-Fermi screening length
Cr  171..  16  exi    Hankel energies for fit to c.d.
Cr                    For free atoms, fit to free-atom density tails.
Cr   39    17  group  group
Cr   40    18  grp2   grp2
Cr  161..  19  hcr    kinky wave hard sphere radii, atomic units
Cr  121..  20  idmod  idmol(l) controls how linearization energy is
Cr                    determined for an orbital of angular momentum l
Cr                    0 float to center of gravity of occupied part of band
Cr                    1 freeze to spec'd logarithmic derivative
Cr                    2 freeze to spec'd linearization energy
Cr                    3 float to natural center of band (C parameter)
Cr  321    21  idu    identifies l-channels with Hubbard U (LDA+U)
Cr  131..  22  idxdn  idxdn
Cr  333..  23  iq1    (tbfit) switches to freeze TB site energies
Cr  337..  24  ivso   (tbfit) switches to freeze TB spin-orbit params
Cr  329    25  jh     LDA+U J parameters for each l-channel
Cr   12    26  kmxh   k cutoffs for head augmentation expansion
Cr   14    27  kmxt   k cutoffs for tail augmentation expansion
Cr   56    28  kmxv   k-cutoff for 1-center projection of free-atom rho
Cr   25    29  lfoca  switch specifying treatment of core density
Cr   11    30  lmxa   l cutoff for augmentation expansion
Cr   31    31  lmxb   highest l in basis
Cr   32    32  lmxf   fitting l
Cr   13    33  lmxl   cutoff for local density
Cr   43    34  lmxpb  l-cutoff for product basis
Cr   49    35  lxi    l-cutoffs for interstital product basis
Cr    3    36  mass   atomic mass
Cr   38    37  mxcst  species-specific constraints of various types:
Cr                    bit function
Cr                    1   suppresses c.d. mixing for this species (ASA)
Cr                    2   Exclude this species when auto-resizing sphere radii
Cr                    4   Freeze augmentation w.f. for this species (FP)
Cr   28    38  name   species name
Cr    6    39  naug   number of k,L for augmentation
Cr  181    40  ngcut  orbital-dependent G cutoffs (for nfp basis)
Cr    7    41  norb   number of orbitals in basis
Cr   16    42  norp   number of parameters needed to define orbital
Cr   21    43  nr     nr for mesh
Cr    5    44  ntorb  number of orbital types in basis
Cr   41    45  nxi    Number of energies in fit of free-atom density tails
Cr  231..  46  orbp   (lmf) parameters defining orbital shape
Cr   29    47  orhoc  pointer to core density
Cr   71..  48  p      log derivative parameter
Cr   44    49  pb1    product basis for linear response
Cr   45    50  pb2    second product basis for linear response
Cr  301... 51  pz     log derivative parameter, semicore states
Cr   91..  52  q      starting q's (charges)
Cr   42    53  qc     core charge
Cr  181..  54  qpol   (tbe) 10 polarisability parameters
Cr   50    55  radxbs radius of atom (for graphics)
Cr   57    56  rcfa   renormalization radius of free atom density, and width
Cr   51    57  rcut   range over which true, smoothed Hankels differ---may
Cr                    be used in the calculation of XC potential
Cr   36    58  rfoca  smoothing radius for frozen core overlap approx
Cr   10    59  rg     rsm for gaussians to fix multipole moments
Cr   35    60  rham   range of hamiltonian
Cr   52    61  rint   range of interstitial product basis
Cr    4    62  rmt    augmentation radius
Cr   54    63  rs3    Lower bound to rsm for local orbital
Cr    9    64  rsma   rsm for augmentation expansion
Cr    8    65  rsmfa  rsm to fit free atom density
Cr   15    66  rsmv   rsm for one-center projection of free-atom density
Cr                    In TCF, smoothing radius to calculate exchange
Cr                    (artificially smoothed density)
Cr    1    67  size   leading dimension of sspec
Cr   30    68  stc    core kinetic energy
Cr   33    69  stni   (tbe) Stoner I
Cr   34    70  stnm   Stoner mmax
Cr  325    71  uh     Hubbard U
Cr   55    72  vmtz   Asymptotic potential for fitting functions at rmt
Cr  341    73  vso    (tbe) spin-orbit parameters
Cr    2    74  z      atomic number
C ----------------------------------------------------------------
C       data  /20,111,24,211,46,17,18,26,37,53,341,291,22,27,23,171,        39,40,161,121,321,131,333,337,329,12,14,56,25,11,31,32,13,43,49,        3,38,28,6,181,7,16,21,5,41,231,29,71,44,45,301,91,42,181,50,57,        51,36,10,35,52,4,54,9,8,15,1,30,33,34,325,55,341,2/
C       data / 5*4,1,10*4,2,2,4,5*2,4,10*2,4,2,1,7*2,4,2,                  4,1,1,16*4,-1,7*4/
C  1   1 -1       size  1
C  2   2  4          z  1
C  3   3  4       mass  1
C  4   4  4        rmt  1
C  5   5  2      ntorb  1
C  6   6  2       naug  1
C  7   7  2       norb  1
C  8   8  4      rsmfa  1
C  9   9  4       rsma  1
C 10  10  4         rg  1
C 11  11  2       lmxa  1
C 12  12  2       kmxh  1
C 13  13  2       lmxl  1
C 14  14  2       kmxt  1
C 15  15  4       rsmv  1
C 16  16  2       norp  1
C 17  17  1      coreh  1
C 18  18  4      coreq  2
C 19  20  4          a  1
C 20  21  2         nr  1
C 21  22  4       eref  1
C 22  23  4        etf  1
C 23  24  4       beta  1
C 24  25  2      lfoca  1
C 25  26  4      ctail  1
C 26  27  4      etail  1
C 27  28  1       name  1
C 28  29  2      orhoc  1 possible_structure
C 29  30  4        stc  1
C 30  31  2       lmxb  1
C 31  32  2       lmxf  1
C 32  33  4       stni  1
C 33  34  4       stnm  1
C 34  35  4       rham  1
C 35  36  4      rfoca  1
C 36  37  4         dv  1
C 37  38  2      mxcst  1
C 38  39  2      group  1
C 39  40  2       grp2  1
C 40  41  2        nxi  1
C 41  42  4         qc  1
C 42  43  2      lmxpb  1
C 43  44  1        pb1  1
C 44  45  1        pb2  1
C 45  46  4     colxbs  3
C 46  49  2        lxi  1
C 47  50  4     radxbs  1
C 48  51  4       rcut  1
C 49  52  4       rint  1
C 50  53  4        eh3  1
C 51  54  4        rs3  1
C 52  55  4       vmtz  1
C 53  56  2       kmxv  1
C 54  57  4       rcfa  2
C 55  71  4          p 20
C 56  91  4          q 20
C 57 111  4      alpha 10
C 58 121  2      idmod 10
C 59 131  2      idxdn 30
C 60 161  4        hcr 10
C 61 171  4        exi 10
C 62 181  2      ngcut 30
C 63 181  4       qpol 10
C 64 211  4       chfa 20
C 65 231  4       orbp 60 possible_structure
C 66 291  4        enu 10
C 67 301  4         pz 20
C 68 321  2        idu  4
C 69 325  4         uh  4
C 70 329  4         jh  4
C 71 333  2        iq1  4
C 72 337  2       ivso  4
C 73 341  4       ehvl 10
C 74 341  4        vso  4
      type s_spec
      real(8)  ::   size
      real(8)   ::   z
      real(8)   ::   mass
      real(8)   ::   rmt
      integer(8)   ::   ntorb
      integer(8)   ::   naug
      integer(8)   ::   norb
      real(8)   ::   rsmfa
      real(8)   ::   rsma
      real(8)   ::   rg
      integer(8)   ::   lmxa
      integer(8)   ::   kmxh
      integer(8)   ::   lmxl
      integer(8)   ::   kmxt
      real(8)   ::   rsmv
      integer(8)   ::   norp
      real(8)   ::   coreh !string
      real(8)   ::   coreq(2)
      real(8)   ::   a
      integer(8)   ::   nr
      real(8)   ::   eref
      real(8)   ::   etf
      real(8)   ::   beta
      integer(8)   ::   lfoca
      real(8)   ::   ctail
      real(8)   ::   etail
      real(8)   ::   name !string
      integer(8)   ::   orhoc
      real(8)   ::   stc
      integer(8)   ::   lmxb
      integer(8)   ::   lmxf
      real(8)   ::   stni
      real(8)   ::   stnm
      real(8)   ::   rham
      real(8)   ::   rfoca
      real(8)   ::   dv
      integer(8)   ::   mxcst
      integer(8)   ::   group
      integer(8)   ::   grp2
      integer(8)   ::   nxi
      real(8)   ::   qc
      integer(8)   ::   lmxpb
      real(8)   ::   pb1 !string
      real(8)   ::   pb2 !string
      real(8)   ::   colxbs(3)
      integer(8)   ::   lxi
      real(8)   ::   radxbs
      real(8)   ::   rcut
      real(8)   ::   rint
      real(8)   ::   eh3
      real(8)   ::   rs3
      real(8)   ::   vmtz
      integer(8)   ::   kmxv
      real(8)   ::   rcfa(2)
      real(8)   ::   p(20)
      real(8)   ::   q(20)
      real(8)   ::   alpha(10)
      integer(8)   ::   idmod(10)
      integer(8)   ::   idxdn(30)
      real(8)   ::   hcr(10)
      real(8)   ::   exi(10)
      integer(8)   ::   ngcut(30)
      real(8)   ::   qpol(10)
      real(8)   ::   chfa(20)
      real(8)   ::   orbp(60)
      real(8)   ::   enu(10)
      real(8)   ::   pz(20)
      integer(8)   ::   idu(4)
      real(8)   ::   uh(4)
      real(8)   ::   jh(4)
      integer(8)   ::   iq1(4)
      integer(8)   ::   ivso(4)
      real(8)   ::   ehvl(10)
      real(8)   ::   vso(4)
      end type s_spec

Cr --- explanation of  usite ---
Cr   off offe  name    purpose
Cr    2    1   bfield  site magnetic field (not used)
Cr    5    2   clabel  string labelling class
Cr    6    3   class   class index
Cr    7    4   cli     (RGF) cluster index
Cr    8    5   delta   (tbe) vector of on-site energy shifts
Cr   14    6   dpole   point charge dipole (molecules)
Cr   17    7   eula    Euler angles (noncollinear magnetism)
Cr   20    8   force   force (dynamics)
Cr   23    9   mpole   point charge monopole (molecules)
Cr   24   10   ndelta  size of delta
Cr   25   11   norb    dimension contributed by this site to hamilt.
Cr   26   12   offh    hamiltonian offset
Cr   27   13   orho    pointer to MT density
Cr   28   14   ov0     pointer to potential that defines wave functions
Cr   29   15   ov1     pointer to spherical part of MT potential
Cr   30   16   pl      (PGF) principal layer index
Cr   31   17   plv     (PGF) principal layer potential index
Cr   32   18   pnu     log derivative parameter
Cr   52   19   pos     true coordinates of atom
Cr   55   20   pos0    coordinates of atom before distortion
Cr   58   21   pz      log derivative parameter, semicore states
Cr   78   22   relax   (dynamics) flags which coordinates to relax
Cr   81   23   sid     (RGF) Site id, used to merge coincident clusters
Cr   82   24   size    size of this structure
Cr   83   25   spec    species index
Cr   84   26   vel     velocity
Cr   87   27   vshft   constant potential shift for this site
C ----------------------------------------------------------------
C       data  /2,5,6,7,8,14,17,20,23,24,25,26,27,28,29,30,31,32,52,        55,58,78,81,82,83,84,87,88/
C       data / 4,1,2,2,5*4,8*2,4*4,2,2,-1,2,4,4/
C  1   2  4     bfield  3
C  2   5  1     clabel  1
C  3   6  2      class  1
C  4   7  2        cli  1
C  5   8  4      delta  6
C  6  14  4      dpole  3
C  7  17  4       eula  3
C  8  20  4      force  3
C  9  23  4      mpole  1
C 10  24  2     ndelta  1
C 11  25  2       norb  1
C 12  26  2       offh  1 possible_structure
C 13  27  2       orho  1 possible_structure
C 14  28  2        ov0  1 possible_structure
C 15  29  2        ov1  1 possible_structure
C 16  30  2         pl  1
C 17  31  2        plv  1
C 18  32  4        pnu 20
C 19  52  4        pos  3
C 20  55  4       pos0  3
C 21  58  4         pz 20
C 22  78  2      relax  3
C 23  81  2        sid  1
C 24  82 -1       size  1
C 25  83  2       spec  1
C 26  84  4        vel  3
C 27  87  4      vshft  1
      type s_site
      real(8)  ::   size
      real(8)   ::   bfield(3)
      real(8)   ::   clabel !string
      integer(8)   ::   class
      integer(8)   ::   cli
      real(8)   ::   delta(6)
      real(8)   ::   dpole(3)
      real(8)   ::   eula(3)
      real(8)   ::   force(3)
      real(8)   ::   mpole
      integer(8)   ::   ndelta
      integer(8)   ::   norb
CCHSTRUC       integer(8)   ::   offh ! not used because of undef
CCHSTRUC       integer(8)   ::   orho ! not used because of undef
      integer(8)   ::   ov0
      integer(8)   ::   ov1
      integer(8)   ::   pl
      integer(8)   ::   plv
      real(8)   ::   pnu(20)
      real(8)   ::   pos(3)
      real(8)   ::   pos0(3)
      real(8)   ::   pz(20)
      integer(8)   ::   relax(3)
      integer(8)   ::   sid
      integer(8)   ::   spec
      real(8)   ::   vel(3)
      real(8)   ::   vshft
      end type s_site

Cr --- explanation of  ustr ---
Cr   off offe  name    purpose
Cr    2    1   adec    log spacing in distance for TCF fit
Cr    3    2   amode   mode for calculating screening alpha's
Cr                    1s digit
Cr                     0: alpha=alp0, (which is input)
Cr                     1: alpha=gamma, calculated from pot pars.
Cr                        NO LONGER SUPPORTED
Cr                     2: alpha calculated for hard spheres.
Cr                        NO LONGER SUPPORTED; see mktral
Cr                    10s digit 0: adot=Andersen's localized
Cr                     0: adot=Andersen's localized
Cr                     1: adot so Kanpur notes Eq 3.87 is zero
Cr                     2: adot calculated for hard spheres.
Cr                        NO LONGER SUPPORTED; see mktral
Cr    4    3   drwats  padding beyond cluster for Watson sphere
Cr    5    4   iinv    Parameters for iterative inversion of tb-strx
Cr                     arg 1: number of sites for cutoff
Cr                     arg 2: number of iterations
Cr                     arg 3: allowed tolerance
Cr                     arg 4-5: not used now.
Cr   10    5   ivl     Type of val-lap functions
Cr   11    6   kaps    Kinetic energies at which strux evaluated (NMTO)
Cr   17    7   lequiv  look for equivalent connecting vectors
Cr   18    8   lmaxw   lmax for Watson sphere
Cr   19    9   loka    conventions for Hankel and Bessel functions
Cr                     0 Methfessel's conventions
Cr                     1 Andersen's 2nd gen LMTO conventions
Cr   20   10   lshow   show strux
Cr   21   11   mxnbr   maximum number of neighbors
Cr   22   12   nalf    polynomial order for d-dependence of TCF fit coffs
Cr   23   13   nbisi   mesh parameters for TCF
Cr   26   14   ncupl   reduced polynomial order in TCF
Cr   27   15   ndust   no. iterations to improve on TCF
Cr   28   16   nitab   ?? same as nttab
Cr   29   17   nkaps   number of kaps
Cr   30   18   nttab   number of pairs in this table
Cr   31   19   oadot   offset to screening alpha's, e. derivative
Cr   32   20   oalp    offset to screening alpha's
Cr   33   21   oiax    offset to neighbor table
Cr   34   22   ontab   offset to number-of-pairs table
Cr   35   23   os      offset to strux
Cr   36   24   osdot   offset to strux energy derivative
Cr   37   25   rfit    ?? (not used) default hard core radius
Cr   38   26   rmax    range of strux (tbe: range of hamiltonian)
Cr   39   27   rmaxg   range of val-lap strux
Cr   40   28   size    size of this structure
Cr   41   29   skmsh   parameters for 3rd gen LMTO energy interpolation
Cr                     not used?
Cr   47   30   wztcf   weighting of z-points in TCF
C ----------------------------------------------------------------
C       data  /2,3,4,5,10,11,17,18,19,20,21,22,23,26,27,28,29,30,31,        32,33,34,35,36,37,38,39,40,41,47,48/
C       data / 4,2,4,4,2,4,18*2,4,4,4,-1,4,4/
C  1   2  4       adec  1
C  2   3  2      amode  1
C  3   4  4     drwats  1
C  4   5  4       iinv  5
C  5  10  2        ivl  1
C  6  11  4       kaps  6
C  7  17  2     lequiv  1
C  8  18  2      lmaxw  1
C  9  19  2       loka  1
C 10  20  2      lshow  1
C 11  21  2      mxnbr  1
C 12  22  2       nalf  1
C 13  23  2      nbisi  3
C 14  26  2      ncupl  1
C 15  27  2      ndust  1
C 16  28  2      nitab  1
C 17  29  2      nkaps  1
C 18  30  2      nttab  1
C 19  31  2      oadot  1 possible_structure
C 20  32  2       oalp  1 possible_structure
C 21  33  2       oiax  1 possible_structure
C 22  34  2      ontab  1 possible_structure
C 23  35  2         os  1 possible_structure
C 24  36  2      osdot  1 possible_structure
C 25  37  4       rfit  1
C 26  38  4       rmax  1
C 27  39  4      rmaxg  1
C 28  40 -1       size  1
C 29  41  4      skmsh  6
C 30  47  4      wztcf  1
      type s_str
      real(8)  ::   size
      real(8)   ::   adec
      integer(8)   ::   amode
      real(8)   ::   drwats
      real(8)   ::   iinv(5)
      integer(8)   ::   ivl
      real(8)   ::   kaps(6)
      integer(8)   ::   lequiv
      integer(8)   ::   lmaxw
      integer(8)   ::   loka
      integer(8)   ::   lshow
      integer(8)   ::   mxnbr
      integer(8)   ::   nalf
      integer(8)   ::   nbisi(3)
      integer(8)   ::   ncupl
      integer(8)   ::   ndust
      integer(8)   ::   nitab
      integer(8)   ::   nkaps
      integer(8)   ::   nttab
CCHSTRUC       integer(8)   ::   oadot ! not used because of undef
CCHSTRUC       integer(8)   ::   oalp ! not used because of undef
CCHSTRUC       integer(8)   ::   oiax ! not used because of undef
CCHSTRUC       integer(8)   ::   ontab ! not used because of undef
CCHSTRUC       integer(8)   ::   os ! not used because of undef
CCHSTRUC       integer(8)   ::   osdot ! not used because of undef
      real(8)   ::   rfit
      real(8)   ::   rmax
      real(8)   ::   rmaxg
      real(8)   ::   skmsh(6)
      real(8)   ::   wztcf
      end type s_str

Cr --- explanation of  utb ---
Cr   off offe  name    purpose
Cr    2    1   alam   (tbfit) lambda variable, Levenburg-Marquardt (Numerical Recipes)
Cr    3    2   alsc   (tbfit) scale factor for Levenburg-Marquardt lambda
Cr    4    3   ebfit  (tbfit) energy window for fitting eigenvalues
Cr                    lower and upper limit of eigenvalues to be fit
Cr    6    4   fmode  (tbfit) mode for fitting TB parmeters
Cr                    =0 fit parameters used in tbham, =1 not implemented yet
Cr    7    5   nbfit  (tbfit) min/max band numbers to fit
Cr                    lower and upper band indices to be fit for all k-points
Cr    9    6   rmfit  (tbfit) real space range to fit TB parameters (fmode=1)
Cr   11    7   size   size of this structure
C ----------------------------------------------------------------
C       data  /2,3,4,6,7,9,11,12/
C       data  /4,4,4,2,2,4,-1/
C  1   2  4       alam  1
C  2   3  4       alsc  1
C  3   4  4      ebfit  2
C  4   6  2      fmode  1
C  5   7  2      nbfit  2
C  6   9  4      rmfit  2
C  7  11 -1       size  1
      type s_tb
      real(8)  ::   size
      real(8)   ::   alam
      real(8)   ::   alsc
      real(8)   ::   ebfit(2)
      integer(8)   ::   fmode
      integer(8)   ::   nbfit(2)
      real(8)   ::   rmfit(2)
      end type s_tb
      end module m_struc_def



c------------------------------------------------------------------
      module m_struc_func

c      interface struc_eval_io
c      module procedure struc_eval_io_r8, struc_eval_io_i8
c     .     , struc_eval_io_r8v, struc_eval_io_i8v
c      end interface struc_eval_io

      contains

      integer function uarray_size()
#if MPI|MPIK|MPIADDRESS
      use m_struc_def, only: s_array
      implicit none
      type(s_array):: v(2)
      integer:: iadd(2), ierr,imod,n
      call mpi_get_address(v(1),iadd(1),ierr)
      call mpi_get_address(v(2),iadd(2),ierr)
      n=iadd(2)-iadd(1)
      imod=mod(n,8)
      n=n/8
      if (imod.ne.0) n=n+1
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      uarray_size=n
C        write(*,*) 'uarray_size=',n
#else  /* MPIADDRESS*/
      integer:: n
      n=-9999 !n= 32
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      uarray_size=n
#endif  /* MPIADDRESS*/
      end function uarray_size

      integer function ubz_size()
#if MPI|MPIK|MPIADDRESS
      use m_struc_def, only: s_bz
      implicit none
      type(s_bz):: v(2)
      integer:: iadd(2), ierr,imod,n
      call mpi_get_address(v(1),iadd(1),ierr)
      call mpi_get_address(v(2),iadd(2),ierr)
      n=iadd(2)-iadd(1)
      imod=mod(n,8)
      n=n/8
      if (imod.ne.0) n=n+1
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      ubz_size=n
C        write(*,*) 'ubz_size=',n
#else  /* MPIADDRESS*/
      integer:: n
      n=-9999 !n= 47
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      ubz_size=n
#endif  /* MPIADDRESS*/
      end function ubz_size

      integer function uctrl_size()
#if MPI|MPIK|MPIADDRESS
      use m_struc_def, only: s_ctrl
      implicit none
      type(s_ctrl):: v(2)
      integer:: iadd(2), ierr,imod,n
      call mpi_get_address(v(1),iadd(1),ierr)
      call mpi_get_address(v(2),iadd(2),ierr)
      n=iadd(2)-iadd(1)
      imod=mod(n,8)
      n=n/8
      if (imod.ne.0) n=n+1
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      uctrl_size=n
C        write(*,*) 'uctrl_size=',n
#else  /* MPIADDRESS*/
      integer:: n
      n=-9999 !n= 88
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      uctrl_size=n
#endif  /* MPIADDRESS*/
      end function uctrl_size

      integer function ugw_size()
#if MPI|MPIK|MPIADDRESS
      use m_struc_def, only: s_gw
      implicit none
      type(s_gw):: v(2)
      integer:: iadd(2), ierr,imod,n
      call mpi_get_address(v(1),iadd(1),ierr)
      call mpi_get_address(v(2),iadd(2),ierr)
      n=iadd(2)-iadd(1)
      imod=mod(n,8)
      n=n/8
      if (imod.ne.0) n=n+1
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      ugw_size=n
C        write(*,*) 'ugw_size=',n
#else  /* MPIADDRESS*/
      integer:: n
      n=-9999 !n= 22
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      ugw_size=n
#endif  /* MPIADDRESS*/
      end function ugw_size

      integer function uham_size()
#if MPI|MPIK|MPIADDRESS
      use m_struc_def, only: s_ham
      implicit none
      type(s_ham):: v(2)
      integer:: iadd(2), ierr,imod,n
      call mpi_get_address(v(1),iadd(1),ierr)
      call mpi_get_address(v(2),iadd(2),ierr)
      n=iadd(2)-iadd(1)
      imod=mod(n,8)
      n=n/8
      if (imod.ne.0) n=n+1
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      uham_size=n
C        write(*,*) 'uham_size=',n
#else  /* MPIADDRESS*/
      integer:: n
      n=-9999 !n= 129
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      uham_size=n
#endif  /* MPIADDRESS*/
      end function uham_size

      integer function ulat_size()
#if MPI|MPIK|MPIADDRESS
      use m_struc_def, only: s_lat
      implicit none
      type(s_lat):: v(2)
      integer:: iadd(2), ierr,imod,n
      call mpi_get_address(v(1),iadd(1),ierr)
      call mpi_get_address(v(2),iadd(2),ierr)
      n=iadd(2)-iadd(1)
      imod=mod(n,8)
      n=n/8
      if (imod.ne.0) n=n+1
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      ulat_size=n
C        write(*,*) 'ulat_size=',n
#else  /* MPIADDRESS*/
      integer:: n
      n=-9999 !n= 120
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      ulat_size=n
#endif  /* MPIADDRESS*/
      end function ulat_size

      integer function umix_size()
#if MPI|MPIK|MPIADDRESS
      use m_struc_def, only: s_mix
      implicit none
      type(s_mix):: v(2)
      integer:: iadd(2), ierr,imod,n
      call mpi_get_address(v(1),iadd(1),ierr)
      call mpi_get_address(v(2),iadd(2),ierr)
      n=iadd(2)-iadd(1)
      imod=mod(n,8)
      n=n/8
      if (imod.ne.0) n=n+1
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      umix_size=n
C        write(*,*) 'umix_size=',n
#else  /* MPIADDRESS*/
      integer:: n
      n=-9999 !n= 36
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      umix_size=n
#endif  /* MPIADDRESS*/
      end function umix_size

      integer function umove_size()
#if MPI|MPIK|MPIADDRESS
      use m_struc_def, only: s_move
      implicit none
      type(s_move):: v(2)
      integer:: iadd(2), ierr,imod,n
      call mpi_get_address(v(1),iadd(1),ierr)
      call mpi_get_address(v(2),iadd(2),ierr)
      n=iadd(2)-iadd(1)
      imod=mod(n,8)
      n=n/8
      if (imod.ne.0) n=n+1
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      umove_size=n
C        write(*,*) 'umove_size=',n
#else  /* MPIADDRESS*/
      integer:: n
      n=-9999 !n= 34
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      umove_size=n
#endif  /* MPIADDRESS*/
      end function umove_size

      integer function uoptic_size()
#if MPI|MPIK|MPIADDRESS
      use m_struc_def, only: s_optic
      implicit none
      type(s_optic):: v(2)
      integer:: iadd(2), ierr,imod,n
      call mpi_get_address(v(1),iadd(1),ierr)
      call mpi_get_address(v(2),iadd(2),ierr)
      n=iadd(2)-iadd(1)
      imod=mod(n,8)
      n=n/8
      if (imod.ne.0) n=n+1
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      uoptic_size=n
C        write(*,*) 'uoptic_size=',n
#else  /* MPIADDRESS*/
      integer:: n
      n=-9999 !n= 33
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      uoptic_size=n
#endif  /* MPIADDRESS*/
      end function uoptic_size

      integer function uordn_size()
#if MPI|MPIK|MPIADDRESS
      use m_struc_def, only: s_ordn
      implicit none
      type(s_ordn):: v(2)
      integer:: iadd(2), ierr,imod,n
      call mpi_get_address(v(1),iadd(1),ierr)
      call mpi_get_address(v(2),iadd(2),ierr)
      n=iadd(2)-iadd(1)
      imod=mod(n,8)
      n=n/8
      if (imod.ne.0) n=n+1
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      uordn_size=n
C        write(*,*) 'uordn_size=',n
#else  /* MPIADDRESS*/
      integer:: n
      n=-9999 !n= 13
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      uordn_size=n
#endif  /* MPIADDRESS*/
      end function uordn_size

      integer function upot_size()
#if MPI|MPIK|MPIADDRESS
      use m_struc_def, only: s_pot
      implicit none
      type(s_pot):: v(2)
      integer:: iadd(2), ierr,imod,n
      call mpi_get_address(v(1),iadd(1),ierr)
      call mpi_get_address(v(2),iadd(2),ierr)
      n=iadd(2)-iadd(1)
      imod=mod(n,8)
      n=n/8
      if (imod.ne.0) n=n+1
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      upot_size=n
C        write(*,*) 'upot_size=',n
#else  /* MPIADDRESS*/
      integer:: n
      n=-9999 !n= 63
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      upot_size=n
#endif  /* MPIADDRESS*/
      end function upot_size

      integer function usite_size()
#if MPI|MPIK|MPIADDRESS
      use m_struc_def, only: s_site
      implicit none
      type(s_site):: v(2)
      integer:: iadd(2), ierr,imod,n
      call mpi_get_address(v(1),iadd(1),ierr)
      call mpi_get_address(v(2),iadd(2),ierr)
      n=iadd(2)-iadd(1)
      imod=mod(n,8)
      n=n/8
      if (imod.ne.0) n=n+1
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      usite_size=n
C        write(*,*) 'usite_size=',n
#else  /* MPIADDRESS*/
      integer:: n
      n=-9999 !n= 86
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      usite_size=n
#endif  /* MPIADDRESS*/
      end function usite_size

      integer function uspec_size()
#if MPI|MPIK|MPIADDRESS
      use m_struc_def, only: s_spec
      implicit none
      type(s_spec):: v(2)
      integer:: iadd(2), ierr,imod,n
      call mpi_get_address(v(1),iadd(1),ierr)
      call mpi_get_address(v(2),iadd(2),ierr)
      n=iadd(2)-iadd(1)
      imod=mod(n,8)
      n=n/8
      if (imod.ne.0) n=n+1
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      uspec_size=n
C        write(*,*) 'uspec_size=',n
#else  /* MPIADDRESS*/
      integer:: n
      n=-9999 !n= 352
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      uspec_size=n
#endif  /* MPIADDRESS*/
      end function uspec_size

      integer function ustr_size()
#if MPI|MPIK|MPIADDRESS
      use m_struc_def, only: s_str
      implicit none
      type(s_str):: v(2)
      integer:: iadd(2), ierr,imod,n
      call mpi_get_address(v(1),iadd(1),ierr)
      call mpi_get_address(v(2),iadd(2),ierr)
      n=iadd(2)-iadd(1)
      imod=mod(n,8)
      n=n/8
      if (imod.ne.0) n=n+1
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      ustr_size=n
C        write(*,*) 'ustr_size=',n
#else  /* MPIADDRESS*/
      integer:: n
      n=-9999 !n= 46
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      ustr_size=n
#endif  /* MPIADDRESS*/
      end function ustr_size

      integer function utb_size()
#if MPI|MPIK|MPIADDRESS
      use m_struc_def, only: s_tb
      implicit none
      type(s_tb):: v(2)
      integer:: iadd(2), ierr,imod,n
      call mpi_get_address(v(1),iadd(1),ierr)
      call mpi_get_address(v(2),iadd(2),ierr)
      n=iadd(2)-iadd(1)
      imod=mod(n,8)
      n=n/8
      if (imod.ne.0) n=n+1
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      utb_size=n
C        write(*,*) 'utb_size=',n
#else  /* MPIADDRESS*/
      integer:: n
      n=-9999 !n= 10
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
      if (mod(n,2).eq.1) n=n+1
#endif
      utb_size=n
#endif  /* MPIADDRESS*/
      end function utb_size



       subroutine mpibc1_s_ri(struc,mlog,funnam,label)
       use m_struc_def, only: s_ri
       implicit none
#if MPI|MPIK
      include 'mpif.h'
      integer numprocs, ierr
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      character*256 strn
      logical lgunit
      integer procid,master,i_data_size
      integer:: n=0, cast=4
C use i_data_size in order to cast an address to integer
#endif
C ... Passed parameters
      logical mlog
       type(s_ri):: struc
      character funnam*(*), label*(*)
#if MPI|MPIK
       master=0
       call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
       call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
       call mpi_bcast(struc%t, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%i, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%r, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)

      if (mlog) then
        call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
        call strcop(shortname(procid),name,10,'.',ierr)
        namelen(procid) = ierr-1
        call gettime(datim)
        strn = ' '//funnam//' '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//' bcast '//label//
     .  ' (%i %?#n==2#int##%?#n==4#d.p.##%?#n==6#d.c.##)'
        call awrit6(strn,' ',-256,lgunit(3),procid,numprocs,n,cast,cast,
     .  cast)
      endif
#endif
      end subroutine mpibc1_s_ri

       subroutine mpibc1_s_bz(struc,mlog,funnam,label)
       use m_struc_def, only: s_bz
       implicit none
#if MPI|MPIK
      include 'mpif.h'
      integer numprocs, ierr
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      character*256 strn
      logical lgunit
      integer procid,master,i_data_size
      integer:: n=0, cast=4
C use i_data_size in order to cast an address to integer
#endif
C ... Passed parameters
      logical mlog
       type(s_bz):: struc
      character funnam*(*), label*(*)
#if MPI|MPIK
       master=0
       call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
       call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
       call mpi_bcast(struc%size, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%def, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%dosw) ! dosw(2)
       call mpi_bcast(struc%dosw, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ef, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%efmax, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%fsmom, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%lcond) ! lcond(4)
       call mpi_bcast(struc%lcond, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lio, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lmet, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lmull, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lopt, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%lshft) ! lshft(3)
       call mpi_bcast(struc%lshft, i_data_size,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%n, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ndos, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nevmx, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%nkabc) ! nkabc(3)
       call mpi_bcast(struc%nkabc, i_data_size,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nkp, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ntet, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%oidtet, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
Cstruc_makeinit.py skip iv_p_oipq(:) because it is a pointer.
Cstruc_makeinit.py skip rv_p_oqp(:) because it is a pointer.
       call mpi_bcast(struc%ostar, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
Cstruc_makeinit.py skip rv_p_owtkp(:) because it is a pointer.
       call mpi_bcast(struc%range, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%semsh) ! semsh(10)
       call mpi_bcast(struc%semsh, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%stnr, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%w, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%zval, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)

      if (mlog) then
        call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
        call strcop(shortname(procid),name,10,'.',ierr)
        namelen(procid) = ierr-1
        call gettime(datim)
        strn = ' '//funnam//' '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//' bcast '//label//
     .  ' (%i %?#n==2#int##%?#n==4#d.p.##%?#n==6#d.c.##)'
        call awrit6(strn,' ',-256,lgunit(3),procid,numprocs,n,cast,cast,
     .  cast)
      endif
#endif
      end subroutine mpibc1_s_bz

       subroutine mpibc1_s_array(struc,mlog,funnam,label)
       use m_struc_def, only: s_array
       implicit none
#if MPI|MPIK
      include 'mpif.h'
      integer numprocs, ierr
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      character*256 strn
      logical lgunit
      integer procid,master,i_data_size
      integer:: n=0, cast=4
C use i_data_size in order to cast an address to integer
#endif
C ... Passed parameters
      logical mlog
       type(s_array):: struc
      character funnam*(*), label*(*)
#if MPI|MPIK
       master=0
       call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
       call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
       call mpi_bcast(struc%size, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ncl, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nclasp, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nofgl, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nofgr, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%npadl, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%npadr, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%npl, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
Cstruc_makeinit.py skip rv_p_oclabl(:) because it is a pointer.
       call mpi_bcast(struc%ogroup, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
Cstruc_makeinit.py skip iv_p_ohave(:) because it is a pointer.
Cstruc_makeinit.py skip iv_p_oics(:) because it is a pointer.
       call mpi_bcast(struc%oipc, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%oipcp, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%oips, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%onrc, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%onrcp, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%opos, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ormax, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)

      if (mlog) then
        call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
        call strcop(shortname(procid),name,10,'.',ierr)
        namelen(procid) = ierr-1
        call gettime(datim)
        strn = ' '//funnam//' '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//' bcast '//label//
     .  ' (%i %?#n==2#int##%?#n==4#d.p.##%?#n==6#d.c.##)'
        call awrit6(strn,' ',-256,lgunit(3),procid,numprocs,n,cast,cast,
     .  cast)
      endif
#endif
      end subroutine mpibc1_s_array

       subroutine mpibc1_s_ctrl(struc,mlog,funnam,label)
       use m_struc_def, only: s_ctrl
       implicit none
#if MPI|MPIK
      include 'mpif.h'
      integer numprocs, ierr
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      character*256 strn
      logical lgunit
      integer procid,master,i_data_size
      integer:: n=0, cast=4
C use i_data_size in order to cast an address to integer
#endif
C ... Passed parameters
      logical mlog
       type(s_ctrl):: struc
      character funnam*(*), label*(*)
#if MPI|MPIK
       master=0
       call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
       call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
       call mpi_bcast(struc%size, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%defm) ! defm(6)
       call mpi_bcast(struc%defm, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%elin, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lasa, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lbas, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lcd, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lcgf, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ldos, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lfp, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lfrce, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lgen3, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lham, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lmet, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lncol, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%loptc, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lordn, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%lpgf) ! lpgf(2)
       call mpi_bcast(struc%lpgf, i_data_size,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lqp, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lrel, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lrs, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lscr, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%lstonr) ! lstonr(3)
       call mpi_bcast(struc%lstonr, i_data_size,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lstr, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lsx, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ltb, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lves, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lxcf, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%maxit, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%mdprm) ! mdprm(6)
       call mpi_bcast(struc%mdprm, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%modep) ! modep(3)
       call mpi_bcast(struc%modep, i_data_size,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nbas, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nbasp, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nclass, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%nesabc) ! nesabc(3)
       call mpi_bcast(struc%nesabc, i_data_size,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nitmv, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nl, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nmap, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%npl, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nsite, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nspec, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nspin, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nvario, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%omax1) ! omax1(3)
       call mpi_bcast(struc%omax1, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%omax2) ! omax2(3)
       call mpi_bcast(struc%omax2, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%pfloat, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%quit, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%rmaxes, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%rmines, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%sclwsr, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%sdmod, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%sdprm) ! sdprm(5)
       call mpi_bcast(struc%sdprm, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%sdxsi) ! sdxsi(4)
       call mpi_bcast(struc%sdxsi, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%smalit) ! smalit(2)
       call mpi_bcast(struc%smalit, i_data_size,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%tol) ! tol(3)
       call mpi_bcast(struc%tol, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%wsrmax, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%zbak) ! zbak(2)
       call mpi_bcast(struc%zbak, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)

      if (mlog) then
        call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
        call strcop(shortname(procid),name,10,'.',ierr)
        namelen(procid) = ierr-1
        call gettime(datim)
        strn = ' '//funnam//' '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//' bcast '//label//
     .  ' (%i %?#n==2#int##%?#n==4#d.p.##%?#n==6#d.c.##)'
        call awrit6(strn,' ',-256,lgunit(3),procid,numprocs,n,cast,cast,
     .  cast)
      endif
#endif
      end subroutine mpibc1_s_ctrl

       subroutine mpibc1_s_gw(struc,mlog,funnam,label)
       use m_struc_def, only: s_gw
       implicit none
#if MPI|MPIK
      include 'mpif.h'
      integer numprocs, ierr
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      character*256 strn
      logical lgunit
      integer procid,master,i_data_size
      integer:: n=0, cast=4
C use i_data_size in order to cast an address to integer
#endif
C ... Passed parameters
      logical mlog
       type(s_gw):: struc
      character funnam*(*), label*(*)
#if MPI|MPIK
       master=0
       call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
       call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
       call mpi_bcast(struc%size, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%delre) ! delre(2)
       call mpi_bcast(struc%delre, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%deltaw, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%deltax, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ecuts, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%gcutb, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%gcutx, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%gsmear, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lgw, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lmxpb, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%mksig, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nband, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nime, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%nkabc) ! nkabc(3)
       call mpi_bcast(struc%nkabc, i_data_size,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%pb1, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%pb2, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%pbtol, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%qoffp, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%rdsig, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)

      if (mlog) then
        call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
        call strcop(shortname(procid),name,10,'.',ierr)
        namelen(procid) = ierr-1
        call gettime(datim)
        strn = ' '//funnam//' '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//' bcast '//label//
     .  ' (%i %?#n==2#int##%?#n==4#d.p.##%?#n==6#d.c.##)'
        call awrit6(strn,' ',-256,lgunit(3),procid,numprocs,n,cast,cast,
     .  cast)
      endif
#endif
      end subroutine mpibc1_s_gw

       subroutine mpibc1_s_ham(struc,mlog,funnam,label)
       use m_struc_def, only: s_ham
       implicit none
#if MPI|MPIK
      include 'mpif.h'
      integer numprocs, ierr
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      character*256 strn
      logical lgunit
      integer procid,master,i_data_size
      integer:: n=0, cast=4
C use i_data_size in order to cast an address to integer
#endif
C ... Passed parameters
      logical mlog
       type(s_ham):: struc
      character funnam*(*), label*(*)
#if MPI|MPIK
       master=0
       call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
       call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
       call mpi_bcast(struc%size, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%alfsi, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%amgm, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%bandw, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%dabc) ! dabc(3)
       call mpi_bcast(struc%dabc, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ehf, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ehk, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%elind, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%eterms) ! eterms(20)
       call mpi_bcast(struc%eterms, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%hord, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%kmto) ! kmto(6)
       call mpi_bcast(struc%kmto, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lasa, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%ldham) ! ldham(16)
       call mpi_bcast(struc%ldham, i_data_size,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lgen3, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lham, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lmaxu, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lmxax, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lncol, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lsig, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ltb, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lxcf, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nbf, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ndham, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ndhrs, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ndofH, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%neula, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nkaph, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nlibu, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nmto, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%npwmin, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%npwpad, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nqsig, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ohrs, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%oiaxs, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
Cstruc_makeinit.py skip iv_p_oindxo(:) because it is a pointer.
       call mpi_bcast(struc%ontabs, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
Cstruc_makeinit.py skip iv_p_ooffH(:) because it is a pointer.
Cstruc_makeinit.py skip rv_p_oqsig(:) because it is a pointer.
       call mpi_bcast(struc%oveps, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%pmax) ! pmax(10)
       call mpi_bcast(struc%pmax, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%pmin) ! pmin(10)
       call mpi_bcast(struc%pmin, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%pwemax, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%pwemin, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%pwmode, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%qpoff) ! qpoff(3)
       call mpi_bcast(struc%qpoff, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%qss) ! qss(4)
       call mpi_bcast(struc%qss, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%rsrnge, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%rsstol, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%seref, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%sigp) ! sigp(10)
       call mpi_bcast(struc%sigp, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%thrpv, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%udiag, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)

      if (mlog) then
        call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
        call strcop(shortname(procid),name,10,'.',ierr)
        namelen(procid) = ierr-1
        call gettime(datim)
        strn = ' '//funnam//' '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//' bcast '//label//
     .  ' (%i %?#n==2#int##%?#n==4#d.p.##%?#n==6#d.c.##)'
        call awrit6(strn,' ',-256,lgunit(3),procid,numprocs,n,cast,cast,
     .  cast)
      endif
#endif
      end subroutine mpibc1_s_ham

       subroutine mpibc1_s_lat(struc,mlog,funnam,label)
       use m_struc_def, only: s_lat
       implicit none
#if MPI|MPIK
      include 'mpif.h'
      integer numprocs, ierr
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      character*256 strn
      logical lgunit
      integer procid,master,i_data_size
      integer:: n=0, cast=4
C use i_data_size in order to cast an address to integer
#endif
C ... Passed parameters
      logical mlog
       type(s_lat):: struc
      character funnam*(*), label*(*)
#if MPI|MPIK
       master=0
       call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
       call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
       call mpi_bcast(struc%size, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%alat, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%as, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%avw, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%awald, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%dist) ! dist(9)
       call mpi_bcast(struc%dist, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%gam) ! gam(4)
       call mpi_bcast(struc%gam, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%gmax, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ldist, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%nabc) ! nabc(3)
       call mpi_bcast(struc%nabc, i_data_size,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ng, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nkd, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nkdmx, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nkq, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nkqmx, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%npgrp, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nsgrp, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%oag, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%obgv, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ocg, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ocy, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%odlv, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ogv, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%oidxcg, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
Cstruc_makeinit.py skip iv_p_oips0(:) because it is a pointer.
Cstruc_makeinit.py skip iv_p_oistab(:) because it is a pointer.
       call mpi_bcast(struc%ojcg, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%okv, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%opos, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
Cstruc_makeinit.py skip rv_p_oqlv(:) because it is a pointer.
       call mpi_bcast(struc%osymgr, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%plat) ! plat(9)
       call mpi_bcast(struc%plat, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%plat0) ! plat0(9)
       call mpi_bcast(struc%plat0, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%plat2) ! plat2(9)
       call mpi_bcast(struc%plat2, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%plate) ! plate(9)
       call mpi_bcast(struc%plate, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%platl) ! platl(9)
       call mpi_bcast(struc%platl, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%platr) ! platr(9)
       call mpi_bcast(struc%platr, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%qlat) ! qlat(9)
       call mpi_bcast(struc%qlat, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%rpad, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%slat) ! slat(9)
       call mpi_bcast(struc%slat, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%tol, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%tolft, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%vol, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)

      if (mlog) then
        call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
        call strcop(shortname(procid),name,10,'.',ierr)
        namelen(procid) = ierr-1
        call gettime(datim)
        strn = ' '//funnam//' '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//' bcast '//label//
     .  ' (%i %?#n==2#int##%?#n==4#d.p.##%?#n==6#d.c.##)'
        call awrit6(strn,' ',-256,lgunit(3),procid,numprocs,n,cast,cast,
     .  cast)
      endif
#endif
      end subroutine mpibc1_s_lat

       subroutine mpibc1_s_mix(struc,mlog,funnam,label)
       use m_struc_def, only: s_mix
       implicit none
#if MPI|MPIK
      include 'mpif.h'
      integer numprocs, ierr
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      character*256 strn
      logical lgunit
      integer procid,master,i_data_size
      integer:: n=0, cast=4
C use i_data_size in order to cast an address to integer
#endif
C ... Passed parameters
      logical mlog
       type(s_mix):: struc
      character funnam*(*), label*(*)
#if MPI|MPIK
       master=0
       call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
       call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
       call mpi_bcast(struc%size, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%b, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%bl, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%bv, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%elind, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%fn, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%kill, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lxpot, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%mmix, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%mode, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%model, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%n, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nitu, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nmix, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nsave, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%r) ! r(3)
       call mpi_bcast(struc%r, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%rms1, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%rms2, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%tj) ! tj(10)
       call mpi_bcast(struc%tj, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%tolu, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%umix, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%w) ! w(3)
       call mpi_bcast(struc%w, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%wc, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)

      if (mlog) then
        call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
        call strcop(shortname(procid),name,10,'.',ierr)
        namelen(procid) = ierr-1
        call gettime(datim)
        strn = ' '//funnam//' '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//' bcast '//label//
     .  ' (%i %?#n==2#int##%?#n==4#d.p.##%?#n==6#d.c.##)'
        call awrit6(strn,' ',-256,lgunit(3),procid,numprocs,n,cast,cast,
     .  cast)
      endif
#endif
      end subroutine mpibc1_s_mix

       subroutine mpibc1_s_move(struc,mlog,funnam,label)
       use m_struc_def, only: s_move
       implicit none
#if MPI|MPIK
      include 'mpif.h'
      integer numprocs, ierr
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      character*256 strn
      logical lgunit
      integer procid,master,i_data_size
      integer:: n=0, cast=4
C use i_data_size in order to cast an address to integer
#endif
C ... Passed parameters
      logical mlog
       type(s_move):: struc
      character funnam*(*), label*(*)
#if MPI|MPIK
       master=0
       call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
       call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
       call mpi_bcast(struc%size, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%ct) ! ct(3)
       call mpi_bcast(struc%ct, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%gyro, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%kt, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%modt) ! modt(3)
       call mpi_bcast(struc%modt, i_data_size,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nmodt, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%prmint) ! prmint(20)
       call mpi_bcast(struc%prmint, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%tnow, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ts, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%tsequ, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%tstot, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)

      if (mlog) then
        call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
        call strcop(shortname(procid),name,10,'.',ierr)
        namelen(procid) = ierr-1
        call gettime(datim)
        strn = ' '//funnam//' '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//' bcast '//label//
     .  ' (%i %?#n==2#int##%?#n==4#d.p.##%?#n==6#d.c.##)'
        call awrit6(strn,' ',-256,lgunit(3),procid,numprocs,n,cast,cast,
     .  cast)
      endif
#endif
      end subroutine mpibc1_s_move

       subroutine mpibc1_s_optic(struc,mlog,funnam,label)
       use m_struc_def, only: s_optic
       implicit none
#if MPI|MPIK
      include 'mpif.h'
      integer numprocs, ierr
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      character*256 strn
      logical lgunit
      integer procid,master,i_data_size
      integer:: n=0, cast=4
C use i_data_size in order to cast an address to integer
#endif
C ... Passed parameters
      logical mlog
       type(s_optic):: struc
      character funnam*(*), label*(*)
#if MPI|MPIK
       master=0
       call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
       call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
       call mpi_bcast(struc%size, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%axes) ! axes(18)
       call mpi_bcast(struc%axes, i_data_size,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%cll, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%cln, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%cls, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%clsite, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%esciss, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%mode, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nchi2, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ne, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%ocrng) ! ocrng(2)
       call mpi_bcast(struc%ocrng, i_data_size,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%unrng) ! unrng(2)
       call mpi_bcast(struc%unrng, i_data_size,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%window) ! window(2)
       call mpi_bcast(struc%window, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)

      if (mlog) then
        call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
        call strcop(shortname(procid),name,10,'.',ierr)
        namelen(procid) = ierr-1
        call gettime(datim)
        strn = ' '//funnam//' '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//' bcast '//label//
     .  ' (%i %?#n==2#int##%?#n==4#d.p.##%?#n==6#d.c.##)'
        call awrit6(strn,' ',-256,lgunit(3),procid,numprocs,n,cast,cast,
     .  cast)
      endif
#endif
      end subroutine mpibc1_s_optic

       subroutine mpibc1_s_ordn(struc,mlog,funnam,label)
       use m_struc_def, only: s_ordn
       implicit none
#if MPI|MPIK
      include 'mpif.h'
      integer numprocs, ierr
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      character*256 strn
      logical lgunit
      integer procid,master,i_data_size
      integer:: n=0, cast=4
C use i_data_size in order to cast an address to integer
#endif
C ... Passed parameters
      logical mlog
       type(s_ordn):: struc
      character funnam*(*), label*(*)
#if MPI|MPIK
       master=0
       call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
       call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
       call mpi_bcast(struc%size, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%efre, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%mode, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ncl, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ndofh, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%rmaxg) ! rmaxg(3)
       call mpi_bcast(struc%rmaxg, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)

      if (mlog) then
        call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
        call strcop(shortname(procid),name,10,'.',ierr)
        namelen(procid) = ierr-1
        call gettime(datim)
        strn = ' '//funnam//' '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//' bcast '//label//
     .  ' (%i %?#n==2#int##%?#n==4#d.p.##%?#n==6#d.c.##)'
        call awrit6(strn,' ',-256,lgunit(3),procid,numprocs,n,cast,cast,
     .  cast)
      endif
#endif
      end subroutine mpibc1_s_ordn

       subroutine mpibc1_s_pot(struc,mlog,funnam,label)
       use m_struc_def, only: s_pot
       implicit none
#if MPI|MPIK
      include 'mpif.h'
      integer numprocs, ierr
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      character*256 strn
      logical lgunit
      integer procid,master,i_data_size
      integer:: n=0, cast=4
C use i_data_size in order to cast an address to integer
#endif
C ... Passed parameters
      logical mlog
       type(s_pot):: struc
      character funnam*(*), label*(*)
#if MPI|MPIK
       master=0
       call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
       call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
       call mpi_bcast(struc%size, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%bfield) ! bfield(4)
       call mpi_bcast(struc%bfield, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nlma, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nlml, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nrhos, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
Cstruc_makeinit.py skip rv_p_obxc(:) because it is a pointer.
Cstruc_makeinit.py skip rv_p_ogrrme(:) because it is a pointer.
Cstruc_makeinit.py skip rv_p_omad(:) because it is a pointer.
Cstruc_makeinit.py skip iv_p_oorhat(:) because it is a pointer.
Cstruc_makeinit.py skip rv_p_opmpol(:) because it is a pointer.
Cstruc_makeinit.py skip rv_p_opnu(:) because it is a pointer.
Cstruc_makeinit.py skip rv_p_opp(:) because it is a pointer.
       call mpi_bcast(struc%opprel, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%oqc, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%oqnu, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%oqpp, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%oqt, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%orhos, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%orhrmx, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
Cstruc_makeinit.py skip zv_p_osmpot(:) because it is a pointer.
Cstruc_makeinit.py skip zv_p_osmrho(:) because it is a pointer.
       call mpi_bcast(struc%osop, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
Cstruc_makeinit.py skip rv_p_ovdif(:) because it is a pointer.
Cstruc_makeinit.py skip rv_p_oves(:) because it is a pointer.
Cstruc_makeinit.py skip rv_p_ovintr(:) because it is a pointer.
Cstruc_makeinit.py skip rv_p_ovrmax(:) because it is a pointer.
       i_data_size=size(struc%vconst) ! vconst(3)
       call mpi_bcast(struc%vconst, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%vmtz, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%vmtz0, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)

      if (mlog) then
        call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
        call strcop(shortname(procid),name,10,'.',ierr)
        namelen(procid) = ierr-1
        call gettime(datim)
        strn = ' '//funnam//' '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//' bcast '//label//
     .  ' (%i %?#n==2#int##%?#n==4#d.p.##%?#n==6#d.c.##)'
        call awrit6(strn,' ',-256,lgunit(3),procid,numprocs,n,cast,cast,
     .  cast)
      endif
#endif
      end subroutine mpibc1_s_pot

       subroutine mpibc1_s_spec(struc,mlog,funnam,label)
       use m_struc_def, only: s_spec
       implicit none
#if MPI|MPIK
      include 'mpif.h'
      integer numprocs, ierr
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      character*256 strn
      logical lgunit
      integer procid,master,i_data_size
      integer:: n=0, cast=4
C use i_data_size in order to cast an address to integer
#endif
C ... Passed parameters
      logical mlog
       type(s_spec):: struc
      character funnam*(*), label*(*)
#if MPI|MPIK
       master=0
       call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
       call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
       call mpi_bcast(struc%size, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%z, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%mass, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%rmt, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ntorb, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%naug, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%norb, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%rsmfa, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%rsma, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%rg, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lmxa, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%kmxh, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lmxl, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%kmxt, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%rsmv, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%norp, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%coreh, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%coreq) ! coreq(2)
       call mpi_bcast(struc%coreq, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%a, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nr, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%eref, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%etf, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%beta, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lfoca, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ctail, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%etail, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%name, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%orhoc, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%stc, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lmxb, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lmxf, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%stni, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%stnm, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%rham, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%rfoca, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%dv, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%mxcst, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%group, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%grp2, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nxi, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%qc, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lmxpb, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%pb1, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%pb2, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%colxbs) ! colxbs(3)
       call mpi_bcast(struc%colxbs, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lxi, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%radxbs, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%rcut, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%rint, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%eh3, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%rs3, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%vmtz, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%kmxv, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%rcfa) ! rcfa(2)
       call mpi_bcast(struc%rcfa, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%p) ! p(20)
       call mpi_bcast(struc%p, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%q) ! q(20)
       call mpi_bcast(struc%q, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%alpha) ! alpha(10)
       call mpi_bcast(struc%alpha, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%idmod) ! idmod(10)
       call mpi_bcast(struc%idmod, i_data_size,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%idxdn) ! idxdn(30)
       call mpi_bcast(struc%idxdn, i_data_size,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%hcr) ! hcr(10)
       call mpi_bcast(struc%hcr, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%exi) ! exi(10)
       call mpi_bcast(struc%exi, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%ngcut) ! ngcut(30)
       call mpi_bcast(struc%ngcut, i_data_size,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%qpol) ! qpol(10)
       call mpi_bcast(struc%qpol, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%chfa) ! chfa(20)
       call mpi_bcast(struc%chfa, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%orbp) ! orbp(60)
       call mpi_bcast(struc%orbp, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%enu) ! enu(10)
       call mpi_bcast(struc%enu, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%pz) ! pz(20)
       call mpi_bcast(struc%pz, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%idu) ! idu(4)
       call mpi_bcast(struc%idu, i_data_size,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%uh) ! uh(4)
       call mpi_bcast(struc%uh, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%jh) ! jh(4)
       call mpi_bcast(struc%jh, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%iq1) ! iq1(4)
       call mpi_bcast(struc%iq1, i_data_size,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%ivso) ! ivso(4)
       call mpi_bcast(struc%ivso, i_data_size,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%ehvl) ! ehvl(10)
       call mpi_bcast(struc%ehvl, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%vso) ! vso(4)
       call mpi_bcast(struc%vso, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)

      if (mlog) then
        call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
        call strcop(shortname(procid),name,10,'.',ierr)
        namelen(procid) = ierr-1
        call gettime(datim)
        strn = ' '//funnam//' '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//' bcast '//label//
     .  ' (%i %?#n==2#int##%?#n==4#d.p.##%?#n==6#d.c.##)'
        call awrit6(strn,' ',-256,lgunit(3),procid,numprocs,n,cast,cast,
     .  cast)
      endif
#endif
      end subroutine mpibc1_s_spec

       subroutine mpibc1_s_site(struc,mlog,funnam,label)
       use m_struc_def, only: s_site
       implicit none
#if MPI|MPIK
      include 'mpif.h'
      integer numprocs, ierr
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      character*256 strn
      logical lgunit
      integer procid,master,i_data_size
      integer:: n=0, cast=4
C use i_data_size in order to cast an address to integer
#endif
C ... Passed parameters
      logical mlog
       type(s_site):: struc
      character funnam*(*), label*(*)
#if MPI|MPIK
       master=0
       call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
       call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
       call mpi_bcast(struc%size, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%bfield) ! bfield(3)
       call mpi_bcast(struc%bfield, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%clabel, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%class, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%cli, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%delta) ! delta(6)
       call mpi_bcast(struc%delta, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%dpole) ! dpole(3)
       call mpi_bcast(struc%dpole, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%eula) ! eula(3)
       call mpi_bcast(struc%eula, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%force) ! force(3)
       call mpi_bcast(struc%force, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%mpole, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ndelta, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%norb, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ov0, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ov1, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%pl, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%plv, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%pnu) ! pnu(20)
       call mpi_bcast(struc%pnu, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%pos) ! pos(3)
       call mpi_bcast(struc%pos, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%pos0) ! pos0(3)
       call mpi_bcast(struc%pos0, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%pz) ! pz(20)
       call mpi_bcast(struc%pz, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%relax) ! relax(3)
       call mpi_bcast(struc%relax, i_data_size,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%sid, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%spec, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%vel) ! vel(3)
       call mpi_bcast(struc%vel, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%vshft, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)

      if (mlog) then
        call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
        call strcop(shortname(procid),name,10,'.',ierr)
        namelen(procid) = ierr-1
        call gettime(datim)
        strn = ' '//funnam//' '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//' bcast '//label//
     .  ' (%i %?#n==2#int##%?#n==4#d.p.##%?#n==6#d.c.##)'
        call awrit6(strn,' ',-256,lgunit(3),procid,numprocs,n,cast,cast,
     .  cast)
      endif
#endif
      end subroutine mpibc1_s_site

       subroutine mpibc1_s_str(struc,mlog,funnam,label)
       use m_struc_def, only: s_str
       implicit none
#if MPI|MPIK
      include 'mpif.h'
      integer numprocs, ierr
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      character*256 strn
      logical lgunit
      integer procid,master,i_data_size
      integer:: n=0, cast=4
C use i_data_size in order to cast an address to integer
#endif
C ... Passed parameters
      logical mlog
       type(s_str):: struc
      character funnam*(*), label*(*)
#if MPI|MPIK
       master=0
       call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
       call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
       call mpi_bcast(struc%size, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%adec, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%amode, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%drwats, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%iinv) ! iinv(5)
       call mpi_bcast(struc%iinv, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ivl, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%kaps) ! kaps(6)
       call mpi_bcast(struc%kaps, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lequiv, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lmaxw, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%loka, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%lshow, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%mxnbr, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nalf, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%nbisi) ! nbisi(3)
       call mpi_bcast(struc%nbisi, i_data_size,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ncupl, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%ndust, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nitab, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nkaps, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%nttab, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%rfit, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%rmax, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%rmaxg, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%skmsh) ! skmsh(6)
       call mpi_bcast(struc%skmsh, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%wztcf, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)

      if (mlog) then
        call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
        call strcop(shortname(procid),name,10,'.',ierr)
        namelen(procid) = ierr-1
        call gettime(datim)
        strn = ' '//funnam//' '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//' bcast '//label//
     .  ' (%i %?#n==2#int##%?#n==4#d.p.##%?#n==6#d.c.##)'
        call awrit6(strn,' ',-256,lgunit(3),procid,numprocs,n,cast,cast,
     .  cast)
      endif
#endif
      end subroutine mpibc1_s_str

       subroutine mpibc1_s_tb(struc,mlog,funnam,label)
       use m_struc_def, only: s_tb
       implicit none
#if MPI|MPIK
      include 'mpif.h'
      integer numprocs, ierr
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      character*256 strn
      logical lgunit
      integer procid,master,i_data_size
      integer:: n=0, cast=4
C use i_data_size in order to cast an address to integer
#endif
C ... Passed parameters
      logical mlog
       type(s_tb):: struc
      character funnam*(*), label*(*)
#if MPI|MPIK
       master=0
       call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
       call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
       call mpi_bcast(struc%size, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%alam, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%alsc, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%ebfit) ! ebfit(2)
       call mpi_bcast(struc%ebfit, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
       call mpi_bcast(struc%fmode, 1,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%nbfit) ! nbfit(2)
       call mpi_bcast(struc%nbfit, i_data_size,MPI_INTEGER8
     , , master, MPI_COMM_WORLD,ierr)
       i_data_size=size(struc%rmfit) ! rmfit(2)
       call mpi_bcast(struc%rmfit, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)

      if (mlog) then
        call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
        call strcop(shortname(procid),name,10,'.',ierr)
        namelen(procid) = ierr-1
        call gettime(datim)
        strn = ' '//funnam//' '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//' bcast '//label//
     .  ' (%i %?#n==2#int##%?#n==4#d.p.##%?#n==6#d.c.##)'
        call awrit6(strn,' ',-256,lgunit(3),procid,numprocs,n,cast,cast,
     .  cast)
      endif
#endif
      end subroutine mpibc1_s_tb

      end module m_struc_func
c------------------------------------------------------------------




c------------------------------------------------------------------


C$$$        subroutine uarray(struc,offe,nli,off,cast,nelt)
C$$$        use m_struc_def
C$$$        implicit none
C$$$        type(s_array):: struc
C$$$        integer:: offe(1),nli,is,off,cast,nelt
C$$$        integer::val
C$$$        if (offe(1).eq.-1) then
C$$$          call uarray_init(struc)
C$$$        else
C$$$          write(*,*) 'uarray offe=',offe(1),' not suppoted'
C$$$          stop
C$$$        endif
C$$$        end subroutine uarray



C$$$        subroutine ubz(struc,offe,nli,off,cast,nelt)
C$$$        use m_struc_def
C$$$        implicit none
C$$$        type(s_bz):: struc
C$$$        integer:: offe(1),nli,is,off,cast,nelt
C$$$        integer::val
C$$$        if (offe(1).eq.-1) then
C$$$          call ubz_init(struc)
C$$$        else
C$$$          write(*,*) 'ubz offe=',offe(1),' not suppoted'
C$$$          stop
C$$$        endif
C$$$        end subroutine ubz



C$$$        subroutine uctrl(struc,offe,nli,off,cast,nelt)
C$$$        use m_struc_def
C$$$        implicit none
C$$$        type(s_ctrl):: struc
C$$$        integer:: offe(1),nli,is,off,cast,nelt
C$$$        integer::val
C$$$        if (offe(1).eq.-1) then
C$$$          call uctrl_init(struc)
C$$$        else
C$$$          write(*,*) 'uctrl offe=',offe(1),' not suppoted'
C$$$          stop
C$$$        endif
C$$$        end subroutine uctrl

C$$$        subroutine ugw_init(struc)
C$$$        use m_struc_def
C$$$         use m_struc_func
C$$$        implicit none
C$$$        type(s_gw):: struc
C$$$        struc%size=ugw_size()
C$$$        struc%delre(:)=0.0d0
C$$$        struc%deltaw=0.0d0
C$$$        struc%deltax=0.0d0
C$$$        struc%ecuts=0.0d0
C$$$        struc%gcutb=0.0d0
C$$$        struc%gcutx=0.0d0
C$$$        struc%gsmear=0.0d0
C$$$        struc%lgw=0
C$$$        struc%lmxpb=0
C$$$        struc%mksig=0
C$$$        struc%nband=0
C$$$        struc%nime=0
C$$$        struc%nkabc(:)=0
C$$$        struc%pb1=0.0d0
C$$$        struc%pb2=0.0d0
C$$$        struc%pbtol=0.0d0
C$$$        struc%qoffp=0.0d0
C$$$        struc%rdsig=0
C$$$        end subroutine ugw_init

C$$$        subroutine ugw_show(struc)
C$$$        use m_struc_def
C$$$        implicit none
C$$$        type(s_gw):: struc
C$$$        write(*,*) struc
C$$$        end subroutine ugw_show

C$$$        subroutine ugw(struc,offe,nli,off,cast,nelt)
C$$$        use m_struc_def
C$$$        implicit none
C$$$        type(s_gw):: struc
C$$$        integer:: offe(1),nli,is,off,cast,nelt
C$$$        integer::val
C$$$        if (offe(1).eq.-1) then
C$$$          call ugw_init(struc)
C$$$        else
C$$$          write(*,*) 'ugw offe=',offe(1),' not suppoted'
C$$$          stop
C$$$        endif
C$$$        end subroutine ugw



C$$$        subroutine uham(struc,offe,nli,off,cast,nelt)
C$$$        use m_struc_def
C$$$        implicit none
C$$$        type(s_ham):: struc
C$$$        integer:: offe(1),nli,is,off,cast,nelt
C$$$        integer::val
C$$$        if (offe(1).eq.-1) then
C$$$          call uham_init(struc)
C$$$        else
C$$$          write(*,*) 'uham offe=',offe(1),' not suppoted'
C$$$          stop
C$$$        endif
C$$$        end subroutine uham



C$$$        subroutine ulat(struc,offe,nli,off,cast,nelt)
C$$$        use m_struc_def
C$$$        implicit none
C$$$        type(s_lat):: struc
C$$$        integer:: offe(1),nli,is,off,cast,nelt
C$$$        integer::val
C$$$        if (offe(1).eq.-1) then
C$$$          call ulat_init(struc)
C$$$        else
C$$$          write(*,*) 'ulat offe=',offe(1),' not suppoted'
C$$$          stop
C$$$        endif
C$$$        end subroutine ulat



C$$$        subroutine umix(struc,offe,nli,off,cast,nelt)
C$$$        use m_struc_def
C$$$        implicit none
C$$$        type(s_mix):: struc
C$$$        integer:: offe(1),nli,is,off,cast,nelt
C$$$        integer::val
C$$$        if (offe(1).eq.-1) then
C$$$          call umix_init(struc)
C$$$        else
C$$$          write(*,*) 'umix offe=',offe(1),' not suppoted'
C$$$          stop
C$$$        endif
C$$$        end subroutine umix



C$$$        subroutine umove(struc,offe,nli,off,cast,nelt)
C$$$        use m_struc_def
C$$$        implicit none
C$$$        type(s_move):: struc
C$$$        integer:: offe(1),nli,is,off,cast,nelt
C$$$        integer::val
C$$$        if (offe(1).eq.-1) then
C$$$          call umove_init(struc)
C$$$        else
C$$$          write(*,*) 'umove offe=',offe(1),' not suppoted'
C$$$          stop
C$$$        endif
C$$$        end subroutine umove



C$$$        subroutine uoptic(struc,offe,nli,off,cast,nelt)
C$$$        use m_struc_def
C$$$        implicit none
C$$$        type(s_optic):: struc
C$$$        integer:: offe(1),nli,is,off,cast,nelt
C$$$        integer::val
C$$$        if (offe(1).eq.-1) then
C$$$          call uoptic_init(struc)
C$$$        else
C$$$          write(*,*) 'uoptic offe=',offe(1),' not suppoted'
C$$$          stop
C$$$        endif
C$$$        end subroutine uoptic



C$$$        subroutine uordn(struc,offe,nli,off,cast,nelt)
C$$$        use m_struc_def
C$$$        implicit none
C$$$        type(s_ordn):: struc
C$$$        integer:: offe(1),nli,is,off,cast,nelt
C$$$        integer::val
C$$$        if (offe(1).eq.-1) then
C$$$          call uordn_init(struc)
C$$$        else
C$$$          write(*,*) 'uordn offe=',offe(1),' not suppoted'
C$$$          stop
C$$$        endif
C$$$        end subroutine uordn



C$$$        subroutine upot(struc,offe,nli,off,cast,nelt)
C$$$        use m_struc_def
C$$$        implicit none
C$$$        type(s_pot):: struc
C$$$        integer:: offe(1),nli,is,off,cast,nelt
C$$$        integer::val
C$$$        if (offe(1).eq.-1) then
C$$$          call upot_init(struc)
C$$$        else
C$$$          write(*,*) 'upot offe=',offe(1),' not suppoted'
C$$$          stop
C$$$        endif
C$$$        end subroutine upot



C$$$        subroutine usite(struc,offe,nli,is,off,cast,nelt)
C$$$        use m_struc_def
C$$$        implicit none
C$$$        type(s_site):: struc
C$$$        integer:: offe(1),nli,is,off,cast,nelt
C$$$        integer::val
C$$$        if (offe(1).eq.-1) then
C$$$          call usite_init(struc)
C$$$        else
C$$$          write(*,*) 'usite offe=',offe(1),' not suppoted'
C$$$          stop
C$$$        endif
C$$$        end subroutine usite



C$$$        subroutine uspec(struc,offe,nli,is,off,cast,nelt)
C$$$        use m_struc_def
C$$$        implicit none
C$$$        type(s_spec):: struc
C$$$        integer:: offe(1),nli,is,off,cast,nelt
C$$$        integer::val
C$$$        if (offe(1).eq.-1) then
C$$$          call uspec_init(struc)
C$$$        else
C$$$          write(*,*) 'uspec offe=',offe(1),' not suppoted'
C$$$          stop
C$$$        endif
C$$$        end subroutine uspec



C$$$        subroutine ustr(struc,offe,nli,off,cast,nelt)
C$$$        use m_struc_def
C$$$        implicit none
C$$$        type(s_str):: struc
C$$$        integer:: offe(1),nli,is,off,cast,nelt
C$$$        integer::val
C$$$        if (offe(1).eq.-1) then
C$$$          call ustr_init(struc)
C$$$        else
C$$$          write(*,*) 'ustr offe=',offe(1),' not suppoted'
C$$$          stop
C$$$        endif
C$$$        end subroutine ustr



C$$$        subroutine utb(struc,offe,nli,off,cast,nelt)
C$$$        use m_struc_def
C$$$        implicit none
C$$$        type(s_tb):: struc
C$$$        integer:: offe(1),nli,is,off,cast,nelt
C$$$        integer::val
C$$$        if (offe(1).eq.-1) then
C$$$          call utb_init(struc)
C$$$        else
C$$$          write(*,*) 'utb offe=',offe(1),' not suppoted'
C$$$          stop
C$$$        endif
C$$$        end subroutine utb


















       subroutine uri_init (struc)
       use m_struc_def
       use m_struc_func
       implicit none
       type(s_ri):: struc
       struc%t=0
       struc%i=0
       struc%r=0.0d0
       end subroutine uri_init

       subroutine ubz_init (struc)
       use m_struc_def
       use m_struc_func
       implicit none
       type(s_bz):: struc
       struc%size= ubz_size()
       struc%def=0.0d0
       struc%dosw=0.0d0
       struc%ef=0.0d0
       struc%efmax=0.0d0
       struc%fsmom=0.0d0
       struc%lcond=0.0d0
       struc%lio=0
       struc%lmet=0
       struc%lmull=0
       struc%lopt=0
       struc%lshft=0
       struc%n=0
       struc%ndos=0
       struc%nevmx=0
       struc%nkabc=0
       struc%nkp=0
       struc%ntet=0
       struc%oidtet=0
       struc%ostar=0
       struc%range=0.0d0
       struc%semsh=0.0d0
       struc%stnr=0.0d0
       struc%w=0.0d0
       struc%zval=0.0d0
       end subroutine ubz_init

       subroutine uarray_init (struc)
       use m_struc_def
       use m_struc_func
       implicit none
       type(s_array):: struc
       struc%size= uarray_size()
       struc%ncl=0
       struc%nclasp=0
       struc%nofgl=0
       struc%nofgr=0
       struc%npadl=0
       struc%npadr=0
       struc%npl=0
       struc%ogroup=0
       struc%oipc=0
       struc%oipcp=0
       struc%oips=0
       struc%onrc=0
       struc%onrcp=0
       struc%opos=0
       struc%ormax=0
       end subroutine uarray_init

       subroutine uctrl_init (struc)
       use m_struc_def
       use m_struc_func
       implicit none
       type(s_ctrl):: struc
       struc%size= uctrl_size()
       struc%defm=0.0d0
       struc%elin=0.0d0
       struc%lasa=0
       struc%lbas=0
       struc%lcd=0
       struc%lcgf=0
       struc%ldos=0
       struc%lfp=0
       struc%lfrce=0
       struc%lgen3=0
       struc%lham=0
       struc%lmet=0
       struc%lncol=0
       struc%loptc=0
       struc%lordn=0
       struc%lpgf=0
       struc%lqp=0
       struc%lrel=0
       struc%lrs=0
       struc%lscr=0
       struc%lstonr=0
       struc%lstr=0
       struc%lsx=0
       struc%ltb=0
       struc%lves=0
       struc%lxcf=0
       struc%maxit=0
       struc%mdprm=0.0d0
       struc%modep=0
       struc%nbas=0
       struc%nbasp=0
       struc%nclass=0
       struc%nesabc=0
       struc%nitmv=0
       struc%nl=0
       struc%nmap=0
       struc%npl=0
       struc%nsite=0
       struc%nspec=0
       struc%nspin=0
       struc%nvario=0
       struc%omax1=0.0d0
       struc%omax2=0.0d0
       struc%pfloat=0
       struc%quit=0
       struc%rmaxes=0.0d0
       struc%rmines=0.0d0
       struc%sclwsr=0.0d0
       struc%sdmod=0
       struc%sdprm=0.0d0
       struc%sdxsi=0.0d0
       struc%smalit=0
       struc%tol=0.0d0
       struc%wsrmax=0.0d0
       struc%zbak=0.0d0
       end subroutine uctrl_init

       subroutine ugw_init (struc)
       use m_struc_def
       use m_struc_func
       implicit none
       type(s_gw):: struc
       struc%size= ugw_size()
       struc%delre=0.0d0
       struc%deltaw=0.0d0
       struc%deltax=0.0d0
       struc%ecuts=0.0d0
       struc%gcutb=0.0d0
       struc%gcutx=0.0d0
       struc%gsmear=0.0d0
       struc%lgw=0
       struc%lmxpb=0
       struc%mksig=0
       struc%nband=0
       struc%nime=0
       struc%nkabc=0
       struc%pb1=0.0d0
       struc%pb2=0.0d0
       struc%pbtol=0.0d0
       struc%qoffp=0.0d0
       struc%rdsig=0
       end subroutine ugw_init

       subroutine uham_init (struc)
       use m_struc_def
       use m_struc_func
       implicit none
       type(s_ham):: struc
       struc%size= uham_size()
       struc%alfsi=0.0d0
       struc%amgm=0.0d0
       struc%bandw=0
       struc%dabc=0.0d0
       struc%ehf=0.0d0
       struc%ehk=0.0d0
       struc%elind=0.0d0
       struc%eterms=0.0d0
       struc%hord=0
       struc%kmto=0.0d0
       struc%lasa=0
       struc%ldham=0
       struc%lgen3=0
       struc%lham=0
       struc%lmaxu=0
       struc%lmxax=0
       struc%lncol=0
       struc%lsig=0
       struc%ltb=0
       struc%lxcf=0
       struc%nbf=0
       struc%ndham=0
       struc%ndhrs=0
       struc%ndofH=0
       struc%neula=0
       struc%nkaph=0
       struc%nlibu=0
       struc%nmto=0
       struc%npwmin=0
       struc%npwpad=0
       struc%nqsig=0
       struc%ohrs=0
       struc%oiaxs=0
       struc%ontabs=0
       struc%oveps=0.0d0
       struc%pmax=0.0d0
       struc%pmin=0.0d0
       struc%pwemax=0.0d0
       struc%pwemin=0.0d0
       struc%pwmode=0
       struc%qpoff=0.0d0
       struc%qss=0.0d0
       struc%rsrnge=0.0d0
       struc%rsstol=0.0d0
       struc%seref=0.0d0
       struc%sigp=0.0d0
       struc%thrpv=0.0d0
       struc%udiag=0
       end subroutine uham_init

       subroutine ulat_init (struc)
       use m_struc_def
       use m_struc_func
       implicit none
       type(s_lat):: struc
       struc%size= ulat_size()
       struc%alat=0.0d0
       struc%as=0.0d0
       struc%avw=0.0d0
       struc%awald=0.0d0
       struc%dist=0.0d0
       struc%gam=0.0d0
       struc%gmax=0.0d0
       struc%ldist=0
       struc%nabc=0
       struc%ng=0
       struc%nkd=0
       struc%nkdmx=0
       struc%nkq=0
       struc%nkqmx=0
       struc%npgrp=0
       struc%nsgrp=0
       struc%oag=0
       struc%obgv=0
       struc%ocg=0
       struc%ocy=0
       struc%odlv=0
       struc%ogv=0
       struc%oidxcg=0
       struc%ojcg=0
       struc%okv=0
       struc%opos=0
       struc%osymgr=0
       struc%plat=0.0d0
       struc%plat0=0.0d0
       struc%plat2=0.0d0
       struc%plate=0.0d0
       struc%platl=0.0d0
       struc%platr=0.0d0
       struc%qlat=0.0d0
       struc%rpad=0.0d0
       struc%slat=0.0d0
       struc%tol=0.0d0
       struc%tolft=0.0d0
       struc%vol=0.0d0
       end subroutine ulat_init

       subroutine umix_init (struc)
       use m_struc_def
       use m_struc_func
       implicit none
       type(s_mix):: struc
       struc%size= umix_size()
       struc%b=0.0d0
       struc%bl=0.0d0
       struc%bv=0.0d0
       struc%elind=0.0d0
       struc%fn=0.0d0
       struc%kill=0
       struc%lxpot=0
       struc%mmix=0
       struc%mode=0
       struc%model=0
       struc%n=0
       struc%nitu=0
       struc%nmix=0
       struc%nsave=0
       struc%r=0.0d0
       struc%rms1=0.0d0
       struc%rms2=0.0d0
       struc%tj=0.0d0
       struc%tolu=0.0d0
       struc%umix=0.0d0
       struc%w=0.0d0
       struc%wc=0.0d0
       end subroutine umix_init

       subroutine umove_init (struc)
       use m_struc_def
       use m_struc_func
       implicit none
       type(s_move):: struc
       struc%size= umove_size()
       struc%ct=0.0d0
       struc%gyro=0.0d0
       struc%kt=0.0d0
       struc%modt=0
       struc%nmodt=0
       struc%prmint=0.0d0
       struc%tnow=0.0d0
       struc%ts=0.0d0
       struc%tsequ=0.0d0
       struc%tstot=0.0d0
       end subroutine umove_init

       subroutine uoptic_init (struc)
       use m_struc_def
       use m_struc_func
       implicit none
       type(s_optic):: struc
       struc%size= uoptic_size()
       struc%axes=0
       struc%cll=0
       struc%cln=0
       struc%cls=0
       struc%clsite=0
       struc%esciss=0.0d0
       struc%mode=0
       struc%nchi2=0
       struc%ne=0
       struc%ocrng=0
       struc%unrng=0
       struc%window=0.0d0
       end subroutine uoptic_init

       subroutine uordn_init (struc)
       use m_struc_def
       use m_struc_func
       implicit none
       type(s_ordn):: struc
       struc%size= uordn_size()
       struc%efre=0.0d0
       struc%mode=0
       struc%ncl=0
       struc%ndofh=0
       struc%rmaxg=0.0d0
       end subroutine uordn_init

       subroutine upot_init (struc)
       use m_struc_def
       use m_struc_func
       implicit none
       type(s_pot):: struc
       struc%size= upot_size()
       struc%bfield=0.0d0
       struc%nlma=0
       struc%nlml=0
       struc%nrhos=0
       struc%opprel=0
       struc%oqc=0
       struc%oqnu=0
       struc%oqpp=0
       struc%oqt=0
       struc%orhos=0
       struc%orhrmx=0
       struc%osop=0
       struc%vconst=0.0d0
       struc%vmtz=0.0d0
       struc%vmtz0=0.0d0
       end subroutine upot_init

       subroutine uspec_init (struc)
       use m_struc_def
       use m_struc_func
       implicit none
       type(s_spec):: struc
       struc%size= uspec_size()
       struc%z=0.0d0
       struc%mass=0.0d0
       struc%rmt=0.0d0
       struc%ntorb=0
       struc%naug=0
       struc%norb=0
       struc%rsmfa=0.0d0
       struc%rsma=0.0d0
       struc%rg=0.0d0
       struc%lmxa=0
       struc%kmxh=0
       struc%lmxl=0
       struc%kmxt=0
       struc%rsmv=0.0d0
       struc%norp=0
       struc%coreh=0.0d0
       struc%coreq=0.0d0
       struc%a=0.0d0
       struc%nr=0
       struc%eref=0.0d0
       struc%etf=0.0d0
       struc%beta=0.0d0
       struc%lfoca=0
       struc%ctail=0.0d0
       struc%etail=0.0d0
       struc%name=0.0d0
       struc%orhoc=0
       struc%stc=0.0d0
       struc%lmxb=0
       struc%lmxf=0
       struc%stni=0.0d0
       struc%stnm=0.0d0
       struc%rham=0.0d0
       struc%rfoca=0.0d0
       struc%dv=0.0d0
       struc%mxcst=0
       struc%group=0
       struc%grp2=0
       struc%nxi=0
       struc%qc=0.0d0
       struc%lmxpb=0
       struc%pb1=0.0d0
       struc%pb2=0.0d0
       struc%colxbs=0.0d0
       struc%lxi=0
       struc%radxbs=0.0d0
       struc%rcut=0.0d0
       struc%rint=0.0d0
       struc%eh3=0.0d0
       struc%rs3=0.0d0
       struc%vmtz=0.0d0
       struc%kmxv=0
       struc%rcfa=0.0d0
       struc%p=0.0d0
       struc%q=0.0d0
       struc%alpha=0.0d0
       struc%idmod=0
       struc%idxdn=0
       struc%hcr=0.0d0
       struc%exi=0.0d0
       struc%ngcut=0
       struc%qpol=0.0d0
       struc%chfa=0.0d0
       struc%orbp=0.0d0
       struc%enu=0.0d0
       struc%pz=0.0d0
       struc%idu=0
       struc%uh=0.0d0
       struc%jh=0.0d0
       struc%iq1=0
       struc%ivso=0
       struc%ehvl=0.0d0
       struc%vso=0.0d0
       end subroutine uspec_init

       subroutine usite_init (struc)
       use m_struc_def
       use m_struc_func
       implicit none
       type(s_site):: struc
       struc%size= usite_size()
       struc%bfield=0.0d0
       struc%clabel=0.0d0
       struc%class=0
       struc%cli=0
       struc%delta=0.0d0
       struc%dpole=0.0d0
       struc%eula=0.0d0
       struc%force=0.0d0
       struc%mpole=0.0d0
       struc%ndelta=0
       struc%norb=0
       struc%ov0=0
       struc%ov1=0
       struc%pl=0
       struc%plv=0
       struc%pnu=0.0d0
       struc%pos=0.0d0
       struc%pos0=0.0d0
       struc%pz=0.0d0
       struc%relax=0
       struc%sid=0
       struc%spec=0
       struc%vel=0.0d0
       struc%vshft=0.0d0
       end subroutine usite_init

       subroutine ustr_init (struc)
       use m_struc_def
       use m_struc_func
       implicit none
       type(s_str):: struc
       struc%size= ustr_size()
       struc%adec=0.0d0
       struc%amode=0
       struc%drwats=0.0d0
       struc%iinv=0.0d0
       struc%ivl=0
       struc%kaps=0.0d0
       struc%lequiv=0
       struc%lmaxw=0
       struc%loka=0
       struc%lshow=0
       struc%mxnbr=0
       struc%nalf=0
       struc%nbisi=0
       struc%ncupl=0
       struc%ndust=0
       struc%nitab=0
       struc%nkaps=0
       struc%nttab=0
       struc%rfit=0.0d0
       struc%rmax=0.0d0
       struc%rmaxg=0.0d0
       struc%skmsh=0.0d0
       struc%wztcf=0.0d0
       end subroutine ustr_init

       subroutine utb_init (struc)
       use m_struc_def
       use m_struc_func
       implicit none
       type(s_tb):: struc
       struc%size= utb_size()
       struc%alam=0.0d0
       struc%alsc=0.0d0
       struc%ebfit=0.0d0
       struc%fmode=0
       struc%nbfit=0
       struc%rmfit=0.0d0
       end subroutine utb_init

