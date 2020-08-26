!! explanations are old documents.
!! A problem in these structures are 'environmental data and iteration data are mixed up'!
!! Check I/O for bndfp.F
      module m_struc_def
!! --- explanation of  ubz ---
Cr name    purpose
Cr def      Uncertainty in Fermi level
Cr dosw     Energy window over which DOS accumulated
Cr ef       Fermi level
Cr efmax    eigenvalues above efmax are not needed
Cr fsmom    fixed-spin moment (fixed spin moment method)
Cr lcond    conductivity, and directions
cCr lio      1  read qp
cCr          2  write qp
cCr          4  write dos
cCr          8  density matrix
cCr          16 save evcs
cCr          32 save dos weights
cCr          NB: 4,32 out of date ; see ldos
Cr lmet     0 insulator   1 3-point scheme  2 save evcs
Cr          3 two passes  4 wt from prior iter
Cr          (numbers, not powers of 2)
Cr lmull    not used
Cr lopt     2 BZMAP option
Cr lshft    shift mesh in each of three directions
Cr n        Polynomial order for Methfessel-Paxton sampling
Cr ndos     No. DOS points (sampling integration, and lmdos)
Cr nevmx    eigenvectors above nevmx are not needed
Cr nkabc    number of divisions in qp mesh
Cr nkp      Number of qp
Cr ntet     number of tetrahedra
Cr odos     offset to total dos
Cr oidtet   offset to tetrahedron corners
Cr oipq     offset to ipq array made in bzmesh
Cr opdos    offset to partial dos
Cr oqp      offset to qp array
Cr ostar    offset to array containing info on star of k
Cr owtkp    offset to array containing qp weights
Cr range    number of FWHM for sampling integration
Cr semsh    nz, modec, emin,emax, ecc,eps delta
Cr size     size of this structure
Cr stnr     parameters for Stoner model (not used)
C     r w        Line broadening for sampling integration
      
c$$$      type s_bz
c$$$      real(8)  ::   size
c$$$      real(8)   ::   def
c$$$      real(8)   ::   dosw(2)
c$$$      real(8):: dosmax !mar2013
c$$$      real(8)   ::   ef
c$$$      real(8)   ::   efmax
c$$$      real(8)   ::   fsmom
c$$$      integer   ::   fsmommethod !Dec 2010 takao
c$$$      real(8)   ::   lcond(4)
c$$$c      integer   ::   lio
c$$$      integer   ::   lmet
c$$$      integer   ::   lmull
c$$$c      integer   ::   lopt
c$$$      integer   ::   lshft(3)
c$$$      integer   ::   n
c$$$      integer   ::   ndos
c$$$      integer   ::   nevmx
c$$$      integer   ::   nkabc(3)
c$$$      integer   ::   nkp
c$$$      integer   ::   ntet
c$$$c      integer , allocatable ::  iv_a_oidtet (:)
c$$$c      integer , allocatable ::  iv_a_oipq (:)
c$$$c      real(8) , pointer ::  rv_p_oqp (:) =>NULL()
c$$$c      integer , allocatable ::  iv_a_ostar (:)
c$$$c      real(8) , allocatable ::  rv_a_owtkp (:)
c$$$      real(8)   ::   range
c$$$      real(8)   ::   w
c$$$      end type s_bz

!!--- explanation of  uarray ---
Cr    name    purpose
Cr   nclasp  (pgf) number of classes including PL -1,npl
Cr   oclabl  offset to class label array (r8tos8 format)
Cr   oics    Index to which species each class belongs
Cr   oipc    index to which class each site belongs
Cr   oipcp   index to which class each site belongs
Cr   oips    index to which species each site belongs
Cr   onrc    number of sites in each class
Cr   onrcp   number of sites in each class
Cr   ormax   sphere radius, by class
Cr   size    size of this structure
      type s_array
      real(8)  ::    size
      integer  ::   nclasp
      end type s_array

!! --- explanation of  uctrl ---
Cr   off offe  name    purpose
Cr    2    1   defm    Lines for deformation when relaxing plat
Cr    8    2   elin    ccor linearization energy, 2C hamiltonian
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
Cr   24   17   noinv    1 do not add inversion !noinv BZ_NOINV
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
Cr                     103 for PBE-GGA
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
      type s_ctrl
      real(8)  ::   size
      real(8)   ::   defm(6)
      real(8)   ::   elin
      integer   ::   lbas
c      integer   ::   lcd
c      integer   ::   lcgf
      integer   ::   ldos
      integer   ::   lfp
      integer   ::   lfrce
      integer   ::   lgen3
      integer   ::   lham
      integer   ::   lmet
c      integer   ::   lncol
c      integer   ::   loptc
c      integer   ::   lordn
c      integer   ::   lpgf(2)
      integer   ::   noinv!,lqp2
      integer   ::   lrel
      integer   ::   lrs
      integer   ::   lscr
c      integer   ::   lstonr(3)
      integer   ::   lstr
      integer   ::   lsx
c      integer   ::   ltb
      integer   ::   lves
      integer   ::   lxcf
      integer   ::   maxit
      real(8)   ::   mdprm(6)
      integer   ::   modep(3)
      integer   ::   nbas
      integer   ::   nbasp
      integer   ::   nclass
      integer   ::   nesabc(3)
      integer   ::   nitmv
      integer   ::   nl
      integer   ::   nmap
c      integer   ::   npl
      integer   ::   nsite
      integer   ::   nspec
      integer   ::   nspin
      integer   ::   nvario
      real(8)   ::   omax1(3)
      real(8)   ::   omax2(3)
      integer   ::   pfloat
      integer   ::   quit
      real(8)   ::   rmaxes
      real(8)   ::   rmines
      real(8)   ::   sclwsr
      integer   ::   sdmod
      real(8)   ::   sdprm(5)
      real(8)   ::   sdxsi(4)
      integer   ::   smalit(2)
      real(8)   ::   tol(3)
      real(8)   ::   wsrmax
      real(8)   ::   zbak(2)
      end type s_ctrl

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
      type s_ham
      real(8)  ::   size
      real(8)   ::   alfsi
      real(8)   ::   amgm
      integer   ::   bandw
      real(8)   ::   dabc(3)
      real(8)   ::   ehf
      real(8)   ::   ehk
      real(8)   ::   elind
      real(8)   ::   eterms(20)
      integer   ::   hord
c      real(8)   ::   kmto(6)
c      integer   ::   lasa
      integer   ::   ldham(16)
      integer   ::   lgen3
      integer   ::   lham
c      integer   ::   lmaxu
      integer   ::   lmxax
c      integer   ::   lncol
      integer   ::   lsig
c      integer   ::   ltb
      integer   ::   lxcf
      integer   ::   nbf
      integer   ::   ndham
      integer   ::   ndhrs
c      integer   ::   ndofH
      integer   ::   neula
      integer   ::   nkaph
      integer   ::   nlibu
c      integer   ::   nmto
      integer   ::   npwmin
      integer   ::   npwpad
      integer   ::   nqsig
      real(8)   ::   oveps,delta_stabilize !takao sep2010
      ! H-> H+ delta_stabilize*O^-1
      real(8)   ::   pmax(10)
      real(8)   ::   pmin(10)
      real(8)   ::   pwemax
      real(8)   ::   pwemin
      integer   ::   pwmode
      real(8)   ::   qpoff(3)
      real(8)   ::   qss(4)
      real(8)   ::   rsrnge
      real(8)   ::   rsstol
      real(8)   ::   seref
      real(8)   ::   sigp(10)
      real(8)   ::   thrpv
      integer   ::   udiag
      real(8):: scaledsigma
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
      type s_lat
      real(8)  ::   size
      real(8)   ::   alat
      real(8)   ::   as
      real(8)   ::   avw
      real(8)   ::   awald
      real(8)   ::   dist(9)
      real(8)   ::   gam(4)
      real(8)   ::   gmax
      integer   ::   ldist
      integer   ::   nabc(3)
      integer   ::   ng
      integer   ::   nkd
      integer   ::   nkdmx
      integer   ::   nkq
      integer   ::   nkqmx
      integer   ::   npgrp
      integer   ::   nsgrp
      real(8)   ::   plat(3,3)
      real(8)   ::   plat0(3,3)
      real(8)   ::   plat2(3,3)
      real(8)   ::   plate(3,3)
      real(8)   ::   platl(3,3)
      real(8)   ::   platr(3,3)
      real(8)   ::   qlat(3,3)
      real(8)   ::   rpad
      real(8)   ::   slat(3,3)
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
      type s_mix
      real(8)  ::   size
      real(8)   ::   b
      real(8)   ::   bl
      real(8)   ::   bv
      real(8)   ::   elind
      character(8)   ::   fn !string
      integer   ::   kill
      integer   ::   lxpot
      integer   ::   mmix
      integer   ::   mode
      integer   ::   model
      integer   ::   n
      integer   ::   nitu
      integer   ::   nmix
      integer   ::   nsave
      character(8)   ::   r(3) !string
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
c$$$      type s_move
c$$$      real(8)  ::   size
c$$$      real(8)   ::   ct(3)
c$$$      real(8)   ::   gyro
c$$$      real(8)   ::   kt
c$$$      integer   ::   modt(3)
c$$$      integer   ::   nmodt
c$$$      real(8)   ::   prmint(20)
c$$$      real(8)   ::   tnow
c$$$      real(8)   ::   ts
c$$$      real(8)   ::   tsequ
c$$$      real(8)   ::   tstot
c$$$      end type s_move

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

c$$$      type s_optic
c$$$      real(8)  ::   size
c$$$      integer   ::   axes(18)
c$$$      integer   ::   cll
c$$$      integer   ::   cln
c$$$      integer   ::   cls
c$$$      integer   ::   clsite
c$$$      real(8)   ::   esciss
c$$$      integer   ::   mode
c$$$      integer   ::   nchi2
c$$$      integer   ::   ne
c$$$      integer   ::   ocrng(2)
c$$$      integer   ::   unrng(2)
c$$$      real(8)   ::   window(2)
c$$$      end type s_optic

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
c$$$      type s_ordn
c$$$      real(8)  ::   size
c$$$      real(8)   ::   efre
c$$$      integer   ::   mode
c$$$c      integer   ::   ncl
c$$$      integer   ::   ndofh
c$$$      real(8)   ::   rmaxg(3)
c$$$      end type s_ordn

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
      type s_rv1
      real(8),allocatable:: v(:)
      end type s_rv1

      type s_pot
      real(8)  ::   size
      real(8)   ::   bfield(4)
      integer   ::   nlma
      integer   ::   nlml
      integer   ::   nrhos
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
      type s_spec
      real(8), allocatable :: rv_a_orhoc(:)
      real(8)  ::   size
      real(8)   ::   z
      real(8)   ::   mass
      real(8)   ::   rmt
      integer   ::   ntorb
      integer   ::   naug
      integer   ::   norb
      real(8)   ::   rsmfa
      real(8)   ::   rsma
      real(8)   ::   rg
      integer   ::   lmxa
      integer   ::   kmxh
      integer   ::   lmxl
      integer   ::   kmxt
      real(8)   ::   rsmv
      integer   ::   norp
      character(8)   ::   coreh !string
      real(8)   ::   coreq(2)
      real(8)   ::   a
      integer   ::   nr
      real(8)   ::   eref
      real(8)   ::   etf
      real(8)   ::   beta
      integer   ::   lfoca
      real(8)   ::   ctail
      real(8)   ::   etail
      character(8)   ::   name !string
      real(8)   ::   stc
      integer   ::   lmxb
      integer   ::   lmxf
      real(8)   ::   stni
      real(8)   ::   stnm
      real(8)   ::   rham
      real(8)   ::   rfoca
      real(8)   ::   dv
      integer   ::   mxcst
      integer   ::   group
      integer   ::   grp2
      integer   ::   nxi
      real(8)   ::   qc
      integer   ::   lmxpb
      character(8)   ::   pb1 !string
      character(8)   ::   pb2 !string
      real(8)   ::   colxbs(3)
      integer   ::   lxi
      real(8)   ::   radxbs
      real(8)   ::   rcut
      real(8)   ::   rint
      real(8)   ::   eh3
      real(8)   ::   rs3
      real(8)   ::   vmtz
      integer   ::   kmxv
      real(8)   ::   rcfa(2)
      real(8)   ::   p(20)
      real(8)   ::   q(20)
c      real(8)   ::   alpha(10)
      integer   ::   idmod(10)
      integer   ::   idxdn(30)
      real(8)   ::   hcr(10)
      real(8)   ::   exi(10)
      integer   ::   ngcut(30)
      real(8)   ::   qpol(10)
      real(8)   ::   chfa(20)
      real(8)   ::   orbp(60)
      real(8)   ::   enu(10)
      real(8)   ::   pz(20)
      integer   ::   idu(4)
      real(8)   ::   uh(4)
      real(8)   ::   jh(4)
      integer   ::   iq1(4)
      integer   ::   ivso(4)
      real(8)   ::   ehvl(10)
      real(8)   ::   vso(4)
      integer :: nmcore !jun2012takao
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
      type s_site
      integer  :: iantiferro !may2015
      real(8)  ::   size
      real(8)   ::   bfield(3)
      character(8)   ::   clabel !string
      integer   ::   class
      integer   ::   cli
      real(8)   ::   delta(6)
      real(8)   ::   dpole(3)
      real(8)   ::   eula(3)
      real(8)   ::   force(3)
      real(8)   ::   mpole
      integer   ::   ndelta
      integer   ::   norb
      real(8) , allocatable ::  rv_a_ov0 (:)
      real(8) , allocatable ::  rv_a_ov1 (:)
      integer   ::   pl
      integer   ::   plv
      real(8)   ::   pnu(20)
      real(8)   ::   pos(3)
      real(8)   ::   pos0(3)
      real(8)   ::   pz(20)
      integer   ::   relax(3)
      integer   ::   sid
      integer   ::   spec
      real(8)   ::   vel(3)
      real(8)   ::   vshft
      end type s_site
      end module m_struc_def




C------------------------------------------------------------------
      Module m_struc_func
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

c$$$      integer function ubz_size()
c$$$#if MPI|MPIK|MPIADDRESS
c$$$      use m_struc_def, only: s_bz
c$$$      implicit none
c$$$      type(s_bz):: v(2)
c$$$      integer:: iadd(2), ierr,imod,n
c$$$      call mpi_get_address(v(1),iadd(1),ierr)
c$$$      call mpi_get_address(v(2),iadd(2),ierr)
c$$$      n=iadd(2)-iadd(1)
c$$$      imod=mod(n,8)
c$$$      n=n/8
c$$$      if (imod.ne.0) n=n+1
c$$$#ifdef ALIGNMENT128
c$$$      if (mod(n,2).eq.1) n=n+1
c$$$#endif
c$$$      ubz_size=n
c$$$C        write(*,*) 'ubz_size=',n
c$$$#else  /* MPIADDRESS*/
c$$$      integer:: n
c$$$      n=-9999 !n= 47
c$$$C  +2 is margin
c$$$C         n=n+2
c$$$#ifdef ALIGNMENT128
c$$$      if (mod(n,2).eq.1) n=n+1
c$$$#endif
c$$$      ubz_size=n
c$$$#endif  /* MPIADDRESS*/
c$$$      end function ubz_size

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

c$$$      integer function ugw_size()
c$$$#if MPI|MPIK|MPIADDRESS
c$$$      use m_struc_def, only: s_gw
c$$$      implicit none
c$$$      type(s_gw):: v(2)
c$$$      integer:: iadd(2), ierr,imod,n
c$$$      call mpi_get_address(v(1),iadd(1),ierr)
c$$$      call mpi_get_address(v(2),iadd(2),ierr)
c$$$      n=iadd(2)-iadd(1)
c$$$      imod=mod(n,8)
c$$$      n=n/8
c$$$      if (imod.ne.0) n=n+1
c$$$#ifdef ALIGNMENT128
c$$$      if (mod(n,2).eq.1) n=n+1
c$$$#endif
c$$$      ugw_size=n
c$$$C        write(*,*) 'ugw_size=',n
c$$$#else  /* MPIADDRESS*/
c$$$      integer:: n
c$$$      n=-9999 !n= 22
c$$$C  +2 is margin
c$$$C         n=n+2
c$$$#ifdef ALIGNMENT128
c$$$      if (mod(n,2).eq.1) n=n+1
c$$$#endif
c$$$      ugw_size=n
c$$$#endif  /* MPIADDRESS*/
c$$$      end function ugw_size

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

c$$$      subroutine mpibc1_s_bz(struc,mlog,funnam,label)
c$$$      use m_struc_def, only: s_bz
c$$$      implicit none
c$$$#if MPI|MPIK
c$$$      include 'mpif.h'
c$$$      integer numprocs, ierr
c$$$      integer MAX_PROCS
c$$$      parameter (MAX_PROCS = 100)
c$$$      integer resultlen
c$$$      character*(MPI_MAX_PROCESSOR_NAME) name
c$$$      character*10 shortname(0:MAX_PROCS-1)
c$$$      character*26 datim
c$$$      integer namelen(0:MAX_PROCS-1)
c$$$      character*256 strn
c$$$      logical lgunit
c$$$      integer procid,master,i_data_size
c$$$      integer:: n=0, cast=4
c$$$C use i_data_size in order to cast an address to integer
c$$$#endif
c$$$C ... Passed parameters
c$$$      logical mlog
c$$$      type(s_bz):: struc
c$$$      character funnam*(*), label*(*)
c$$$#if MPI|MPIK
c$$$      master=0
c$$$      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
c$$$      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
c$$$      call mpi_bcast(struc%size, 1,MPI_REAL8
c$$$     , , master, MPI_COMM_WORLD,ierr)
c$$$      call mpi_bcast(struc%def, 1,MPI_REAL8
c$$$     , , master, MPI_COMM_WORLD,ierr)
c$$$      i_data_size=size(struc%dosw) ! dosw(2)
c$$$      call mpi_bcast(struc%dosw, i_data_size,MPI_REAL8
c$$$     , , master, MPI_COMM_WORLD,ierr)
c$$$      call mpi_bcast(struc%ef, 1,MPI_REAL8
c$$$     , , master, MPI_COMM_WORLD,ierr)
c$$$      call mpi_bcast(struc%efmax, 1,MPI_REAL8
c$$$     , , master, MPI_COMM_WORLD,ierr)
c$$$      call mpi_bcast(struc%fsmom, 1,MPI_REAL8
c$$$     , , master, MPI_COMM_WORLD,ierr)
c$$$      i_data_size=size(struc%lcond) ! lcond(4)
c$$$      call mpi_bcast(struc%lcond, i_data_size,MPI_REAL8
c$$$     , , master, MPI_COMM_WORLD,ierr)
c$$$c      call mpi_bcast(struc%lio, 1,MPI_INTEGER
c$$$c     , , master, MPI_COMM_WORLD,ierr)
c$$$      call mpi_bcast(struc%lmet, 1,MPI_INTEGER
c$$$     , , master, MPI_COMM_WORLD,ierr)
c$$$      call mpi_bcast(struc%lmull, 1,MPI_INTEGER
c$$$     , , master, MPI_COMM_WORLD,ierr)
c$$$c      call mpi_bcast(struc%lopt, 1,MPI_INTEGER
c$$$c     , , master, MPI_COMM_WORLD,ierr)
c$$$      i_data_size=size(struc%lshft) ! lshft(3)
c$$$      call mpi_bcast(struc%lshft, i_data_size,MPI_INTEGER
c$$$     , , master, MPI_COMM_WORLD,ierr)
c$$$      call mpi_bcast(struc%n, 1,MPI_INTEGER
c$$$     , , master, MPI_COMM_WORLD,ierr)
c$$$      call mpi_bcast(struc%ndos, 1,MPI_INTEGER
c$$$     , , master, MPI_COMM_WORLD,ierr)
c$$$      call mpi_bcast(struc%nevmx, 1,MPI_INTEGER
c$$$     , , master, MPI_COMM_WORLD,ierr)
c$$$      i_data_size=size(struc%nkabc) ! nkabc(3)
c$$$      call mpi_bcast(struc%nkabc, i_data_size,MPI_INTEGER
c$$$     , , master, MPI_COMM_WORLD,ierr)
c$$$      call mpi_bcast(struc%nkp, 1,MPI_INTEGER
c$$$     , , master, MPI_COMM_WORLD,ierr)
c$$$      call mpi_bcast(struc%ntet, 1,MPI_INTEGER
c$$$     , , master, MPI_COMM_WORLD,ierr)
c$$$      call mpi_bcast(struc%range, 1,MPI_REAL8
c$$$     , , master, MPI_COMM_WORLD,ierr)
c$$$c      i_data_size=size(struc%semsh) ! semsh(10)
c$$$c      call mpi_bcast(struc%semsh, i_data_size,MPI_REAL8
c$$$c     , , master, MPI_COMM_WORLD,ierr)
c$$$c      call mpi_bcast(struc%stnr, 1,MPI_REAL8
c$$$c     , , master, MPI_COMM_WORLD,ierr)
c$$$      call mpi_bcast(struc%w, 1,MPI_REAL8
c$$$     , , master, MPI_COMM_WORLD,ierr)
c$$$c      call mpi_bcast(struc%zval, 1,MPI_REAL8
c$$$c     , , master, MPI_COMM_WORLD,ierr)
c$$$
c$$$      if (mlog) then
c$$$        call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
c$$$        call strcop(shortname(procid),name,10,'.',ierr)
c$$$        namelen(procid) = ierr-1
c$$$        call gettime(datim)
c$$$        strn = ' '//funnam//' '//datim//' Process %i of %i on '
c$$$     .  //shortname(procid)(1:namelen(procid))//' bcast '//label//
c$$$     .  ' (%i %?#n==2#int##%?#n==4#d.p.##%?#n==6#d.c.##)'
c$$$        call awrit6(strn,' ',-256,lgunit(3),procid,numprocs,n,cast,cast,
c$$$     .  cast)
c$$$      endif
c$$$#endif
c$$$      end subroutine mpibc1_s_bz


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
c      call mpi_bcast(struc%lasa, 1,MPI_INTEGER
c     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lbas, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
c      call mpi_bcast(struc%lcd, 1,MPI_INTEGER
c     , , master, MPI_COMM_WORLD,ierr)
c      call mpi_bcast(struc%lcgf, 1,MPI_INTEGER
c     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%ldos, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lfp, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lfrce, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lgen3, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lham, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lmet, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
c      call mpi_bcast(struc%lncol, 1,MPI_INTEGER
c     , , master, MPI_COMM_WORLD,ierr)
c      call mpi_bcast(struc%loptc, 1,MPI_INTEGER
c     , , master, MPI_COMM_WORLD,ierr)
c      call mpi_bcast(struc%lordn, 1,MPI_INTEGER
c     , , master, MPI_COMM_WORLD,ierr)
c      i_data_size=size(struc%lpgf) ! lpgf(2)
c      call mpi_bcast(struc%lpgf, i_data_size,MPI_INTEGER
c     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%noinv, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
c      call mpi_bcast(struc%lqp2, 1,MPI_INTEGER
c     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lrel, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lrs, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lscr, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
c      i_data_size=size(struc%lstonr) ! lstonr(3)
c      call mpi_bcast(struc%lstonr, i_data_size,MPI_INTEGER
c     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lstr, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lsx, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
c      call mpi_bcast(struc%ltb, 1,MPI_INTEGER
c     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lves, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lxcf, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%maxit, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      i_data_size=size(struc%mdprm) ! mdprm(6)
      call mpi_bcast(struc%mdprm, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      i_data_size=size(struc%modep) ! modep(3)
      call mpi_bcast(struc%modep, i_data_size,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nbas, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nbasp, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nclass, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      i_data_size=size(struc%nesabc) ! nesabc(3)
      call mpi_bcast(struc%nesabc, i_data_size,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nitmv, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nl, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nmap, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
c      call mpi_bcast(struc%npl, 1,MPI_INTEGER
c     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nsite, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nspec, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nspin, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nvario, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      i_data_size=size(struc%omax1) ! omax1(3)
      call mpi_bcast(struc%omax1, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      i_data_size=size(struc%omax2) ! omax2(3)
      call mpi_bcast(struc%omax2, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%pfloat, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%quit, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%rmaxes, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%rmines, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%sclwsr, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%sdmod, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      i_data_size=size(struc%sdprm) ! sdprm(5)
      call mpi_bcast(struc%sdprm, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      i_data_size=size(struc%sdxsi) ! sdxsi(4)
      call mpi_bcast(struc%sdxsi, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      i_data_size=size(struc%smalit) ! smalit(2)
      call mpi_bcast(struc%smalit, i_data_size,MPI_INTEGER
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
      call mpi_bcast(struc%bandw, 1,MPI_INTEGER
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
      call mpi_bcast(struc%hord, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
c      i_data_size=size(struc%kmto) ! kmto(6)
c      call mpi_bcast(struc%kmto, i_data_size,MPI_REAL8
c     , , master, MPI_COMM_WORLD,ierr)
c      call mpi_bcast(struc%lasa, 1,MPI_INTEGER
c     , , master, MPI_COMM_WORLD,ierr)
      i_data_size=size(struc%ldham) ! ldham(16)
      call mpi_bcast(struc%ldham, i_data_size,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lgen3, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lham, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
c      call mpi_bcast(struc%lmaxu, 1,MPI_INTEGER
c     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lmxax, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
c      call mpi_bcast(struc%lncol, 1,MPI_INTEGER
c     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lsig, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
c      call mpi_bcast(struc%ltb, 1,MPI_INTEGER
c     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lxcf, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nbf, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%ndham, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%ndhrs, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
c      call mpi_bcast(struc%ndofH, 1,MPI_INTEGER
c     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%neula, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nkaph, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nlibu, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
c      call mpi_bcast(struc%nmto, 1,MPI_INTEGER
c     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%npwmin, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%npwpad, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nqsig, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%oveps, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%delta_stabilize, 1,MPI_REAL8 !takao sep2010
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
      call mpi_bcast(struc%pwmode, 1,MPI_INTEGER
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
      call mpi_bcast(struc%udiag, 1,MPI_INTEGER
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
      call mpi_bcast(struc%ldist, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      i_data_size=size(struc%nabc) ! nabc(3)
      call mpi_bcast(struc%nabc, i_data_size,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%ng, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nkd, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nkdmx, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nkq, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nkqmx, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%npgrp, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nsgrp, 1,MPI_INTEGER
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
      call mpi_bcast(struc%kill, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lxpot, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%mmix, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%mode, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%model, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%n, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nitu, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nmix, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nsave, 1,MPI_INTEGER
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
      call mpi_bcast(struc%nlma, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nlml, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nrhos, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
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
      call mpi_bcast(struc%ntorb, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%naug, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%norb, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%rsmfa, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%rsma, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%rg, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lmxa, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%kmxh, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lmxl, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%kmxt, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%rsmv, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%norp, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%coreh, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      i_data_size=size(struc%coreq) ! coreq(2)
      call mpi_bcast(struc%coreq, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%a, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nr, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%eref, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%etf, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%beta, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lfoca, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%ctail, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%etail, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%name, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
cki       call mpi_bcast(struc%orhoc, 1,MPI_INTEGER
cki     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%stc, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lmxb, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lmxf, 1,MPI_INTEGER
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
      call mpi_bcast(struc%mxcst, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%group, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%grp2, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%nxi, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%qc, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lmxpb, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%pb1, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%pb2, 1,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      i_data_size=size(struc%colxbs) ! colxbs(3)
      call mpi_bcast(struc%colxbs, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%lxi, 1,MPI_INTEGER
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
      call mpi_bcast(struc%kmxv, 1,MPI_INTEGER
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
c      i_data_size=size(struc%alpha) ! alpha(10)
c      call mpi_bcast(struc%alpha, i_data_size,MPI_REAL8
c     , , master, MPI_COMM_WORLD,ierr)
      i_data_size=size(struc%idmod) ! idmod(10)
      call mpi_bcast(struc%idmod, i_data_size,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      i_data_size=size(struc%idxdn) ! idxdn(30)
      call mpi_bcast(struc%idxdn, i_data_size,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      i_data_size=size(struc%hcr) ! hcr(10)
      call mpi_bcast(struc%hcr, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      i_data_size=size(struc%exi) ! exi(10)
      call mpi_bcast(struc%exi, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      i_data_size=size(struc%ngcut) ! ngcut(30)
      call mpi_bcast(struc%ngcut, i_data_size,MPI_INTEGER
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
      call mpi_bcast(struc%idu, i_data_size,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      i_data_size=size(struc%uh) ! uh(4)
      call mpi_bcast(struc%uh, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      i_data_size=size(struc%jh) ! jh(4)
      call mpi_bcast(struc%jh, i_data_size,MPI_REAL8
     , , master, MPI_COMM_WORLD,ierr)
      i_data_size=size(struc%iq1) ! iq1(4)
      call mpi_bcast(struc%iq1, i_data_size,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      i_data_size=size(struc%ivso) ! ivso(4)
      call mpi_bcast(struc%ivso, i_data_size,MPI_INTEGER
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
      call mpi_bcast(struc%class, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%cli, 1,MPI_INTEGER
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
      call mpi_bcast(struc%ndelta, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%norb, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
Cstruc_makeinit.py skip rv_p_ov0(:) because it is a pointer.
Cstruc_makeinit.py skip rv_p_ov1(:) because it is a pointer.
      call mpi_bcast(struc%pl, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%plv, 1,MPI_INTEGER
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
      call mpi_bcast(struc%relax, i_data_size,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%sid, 1,MPI_INTEGER
     , , master, MPI_COMM_WORLD,ierr)
      call mpi_bcast(struc%spec, 1,MPI_INTEGER
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
      end module m_struc_func
c------------------------------------------------------------------

c$$$      subroutine ubz_init (struc)
c$$$      use m_struc_def
c$$$      use m_struc_func
c$$$      implicit none
c$$$      type(s_bz):: struc
c$$$      struc%size= ubz_size()
c$$$      struc%def=0.0d0
c$$$      struc%dosw=0.0d0
c$$$      struc%ef=0.0d0
c$$$      struc%efmax=0.0d0
c$$$      struc%fsmom=0.0d0
c$$$      struc%lcond=0.0d0
c$$$      struc%lmet=0
c$$$      struc%lmull=0
c$$$      struc%lshft=0
c$$$      struc%n=0
c$$$      struc%ndos=0
c$$$      struc%nevmx=0
c$$$      struc%nkabc=0
c$$$      struc%nkp=0
c$$$      struc%ntet=0
c$$$      struc%range=0.0d0
c$$$      struc%w=0.0d0
c$$$      end subroutine ubz_init

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
      struc%fn=''!0.0d0
      struc%kill=0
      struc%lxpot=0
      struc%mmix=0
      struc%mode=0
      struc%model=0
      struc%n=0
      struc%nitu=0
      struc%nmix=0
      struc%nsave=0
      struc%r='' !0.0d0
      struc%rms1=0.0d0
      struc%rms2=0.0d0
      struc%tj=0.0d0
      struc%tolu=0.0d0
      struc%umix=0.0d0
      struc%w=0.0d0
      struc%wc=0.0d0
      end subroutine umix_init

