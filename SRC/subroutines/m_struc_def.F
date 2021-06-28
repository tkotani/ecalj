!! explanations are old, some are obsolate 
      module m_struc_def
      public s_rv1, s_nv2, s_spec, s_site
      private
      
      type s_rv1
      real(8),allocatable:: v(:)
      end type s_rv1
      
      type s_nv2
      integer,allocatable:: nv2(:,:)
      end type s_nv2
      
!! --- explanation of s_spec --- (some are obsolate)
Cr  name    purpose
Cr  a      a for mesh
Cr  alpha  screening alphas
Cr  beta   mixing beta
Cr  chfa   coefficients to fit of free-atom density tails
Cr  colxbs color of atom (for graphics)
Cr  coreh  core hole channel (char: eg '1s')
Cr  coreq  coreq(1) = charge in core hole channel
Cr         coreq(2) = moment in core hole channel (nsp=2 only)
Cr  ctail  coefficients to fit of free-atom core tail by unsm. Hankel
Cr  dv     constant shift in potential
Cr  eh3    sm Hankel energy for high local orbitals
Cr  ehvl   l-dependent sm Hankel energies for val-lap basis
Cr  enu    linearization energies
Cr  eref   reference energy
Cr  etail  energy to fit of free-atom core tail
Cr  etf    Thomas-Fermi screening length
Cr  exi    Hankel energies for fit to c.d.
Cr         For free atoms, fit to free-atom density tails.
Cr  group  group
Cr  grp2   grp2
Cr  hcr    kinky wave hard sphere radii, atomic units
Cr  idmod  idmol(l) controls how linearization energy is
Cr         determined for an orbital of angular momentum l
Cr         0 float to center of gravity of occupied part of band
Cr         1 freeze to spec'd logarithmic derivative
Cr         2 freeze to spec'd linearization energy
Cr         3 float to natural center of band (C parameter)
Cr  idu    identifies l-channels with Hubbard U (LDA+U)
Cr  idxdn  idxdn
Cr  jh     LDA+U J parameters for each l-channel
Cr  kmxh   k cutoffs for head augmentation expansion
Cr  kmxt   k cutoffs for tail augmentation expansion
Cr  kmxv   k-cutoff for 1-center projection of free-atom rho
Cr  lfoca  switch specifying treatment of core density
Cr  lmxa   l cutoff for augmentation expansion
Cr  lmxb   highest l in basis
Cr  lmxf   fitting l
Cr  lmxl   cutoff for local density
Cr  lmxpb  l-cutoff for product basis
Cr  lxi    l-cutoffs for interstital product basis
Cr  mass   atomic mass
Cr  mxcst  species-specific constraints of various types:
Cr         bit function
Cr         1   suppresses c.d. mixing for this species (ASA)
Cr         2   Exclude this species when auto-resizing sphere radii
Cr         4   Freeze augmentation w.f. for this species (FP)
Cr  name   species name
Cr  naug   number of k,L for augmentation
Cr  ngcut  orbital-dependent G cutoffs (for nfp basis)
Cr  norb   number of orbitals in basis
Cr  norp   number of parameters needed to define orbital
Cr  nr     nr for mesh
Cr  ntorb  number of orbital types in basis
Cr  nxi    Number of energies in fit of free-atom density tails
Cr  orbp   (lmf) parameters defining orbital shape
Cr  orhoc  pointer to core density
Cr  p      log derivative parameter
Cr  pb1    product basis for linear response
Cr  pb2    second product basis for linear response
Cr  pz     log derivative parameter, semicore states
Cr  q      starting q's (charges)
Cr  qc     core charge
Cr  qpol   (tbe) 10 polarisability parameters
Cr  radxbs radius of atom (for graphics)
Cr  rcfa   renormalization radius of free atom density, and width
Cr  rcut   range over which true, smoothed Hankels differ---may
Cr         be used in the calculation of XC potential
Cr  rfoca  smoothing radius for frozen core overlap approx
Cr  rg     rsm for gaussians to fix multipole moments
Cr  rham   range of hamiltonian
Cr  rint   range of interstitial product basis
Cr  rmt    augmentation radius
Cr  rs3    Lower bound to rsm for local orbital
Cr  rsma   rsm for augmentation expansion
Cr  rsmfa  rsm to fit free atom density
Cr  rsmv   rsm for one-center projection of free-atom density
Cr         In TCF, smoothing radius to calculate exchange
Cr         (artificially smoothed density)
Cr  size   leading dimension of sspec
Cr  stc    core kinetic energy
Cr  stni   (tbe) Stoner I
Cr  stnm   Stoner mmax
Cr  uh     Hubbard U
Cr  vmtz   Asymptotic potential for fitting functions at rmt
Cr  vso    (tbe) spin-orbit parameters
Cr  z      atomic number
      type s_spec
      real(8), allocatable :: rv_a_orhoc(:)
c      real(8)  ::   size
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
      integer   ::   idmod(10)
      integer   ::   idxdn(30)
      real(8)   ::   hcr(10)
      real(8)   ::   exi(10)
      integer   ::   ngcut(30)
      real(8)   ::   chfa(20)
      real(8)   ::   orbp(60)
      real(8)   ::   enu(10)
      real(8)   ::   pz(20)
      integer   ::   idu(4)
      real(8)   ::   uh(4)
      real(8)   ::   jh(4)
      real(8)   ::   ehvl(10)
      integer :: nmcore !jun2012takao
      end type s_spec

      
!! --- explanation of s_site --- (some are obsolate)
Cr   bfield  site magnetic field (not used)
Cr   clabel  string labelling class
Cr   class   class index
Cr   eula    Euler angles (noncollinear magnetism)
Cr   force   force (dynamics)
Cr   norb    dimension contributed by this site to hamilt.
Cr   offh    hamiltonian offset
Cr   orho    pointer to MT density
Cr   ov0     pointer to potential that defines wave functions
Cr   ov1     pointer to spherical part of MT potential
Cr   pnu     log derivative parameter
Cr   pos     true coordinates of atom
Cr   pos0    coordinates of atom before distortion
Cr   pz      log derivative parameter, semicore states
Cr   relax   (dynamics) flags which coordinates to relax
Cr   size    size of this structure
Cr   spec    species index
Cr   vel     velocity
Cr   vshft   constant potential shift for this site
      type s_site
      integer  :: iantiferro !may2015
c      real(8)  ::   size
      real(8)   ::   bfield(3)
      character(8)   ::   clabel !string
      integer   ::   class
      real(8)   ::   eula(3)
      real(8)   ::   force(3)
      integer   ::   norb
      real(8) , allocatable ::  rv_a_ov0 (:)
      real(8) , allocatable ::  rv_a_ov1 (:)
      real(8)   ::   pnu(20)
      real(8)   ::   pos(3)
      real(8)   ::   pos0(3)
      real(8)   ::   pz(20)
      integer   ::   relax(3)
      integer   ::   spec
      real(8)   ::   vel(3)
      real(8)   ::   vshft
      end type s_site
      end module m_struc_def
