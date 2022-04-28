!! explanations are old, some are obsolate 
      module m_struc_def
      public s_rv1, s_nv2, s_cv1, s_cv2,s_cv3,s_cv4,s_spec, s_site,s_sblock
      private
      
      type s_rv1
      real(8),allocatable:: v(:)
      end type s_rv1

      type s_cv1
      complex(8),allocatable:: cv(:)
      end type s_cv1
      
      type s_sblock
      complex(8),allocatable:: sdiag(:,:) !(1,1) and (2,2) spin block as sdiag(:,isp)
      complex(8),allocatable:: soffd(:,:) !(1,2) and (2,1) spin block as soffd(:,isp)
      end type s_sblock
      
      type s_nv2
      integer,allocatable:: nv2(:,:)
      end type s_nv2

      type s_cv2
      complex(8),allocatable:: cv2(:,:)
      end type s_cv2
      
      type s_cv3
      complex(8),allocatable:: cv3(:,:,:)
      end type s_cv3
      
      type s_cv4
      complex(8),allocatable:: cv4(:,:,:,:)
      end type s_cv4
      
!! --- explanation of s_spec --- (some are obsolate)
Cr  name    purpose
Cr  a      a for mesh
Cr  alpha  screening alphas
Cr  chfa   coefficients to fit of free-atom density tails
Cr  coreh  core hole channel (char: eg '1s')
Cr  coreq  coreq(1) = charge in core hole channel
Cr         coreq(2) = moment in core hole channel (nsp=2 only)
Cr  ctail  coefficients to fit of free-atom core tail by unsm. Hankel
Cr  eh3    sm Hankel energy for high local orbitals

Cr  eref   reference energy
Cr  etail  energy to fit of free-atom core tail
Cr  exi    Hankel energies for fit to c.d.
Cr         For free atoms, fit to free-atom density tails.
Cr  idmod  idmol(l) controls how linearization energy is
Cr         determined for an orbital of angular momentum l
Cr         0 float to center of gravity of occupied part of band
Cr         1 freeze to spec'd logarithmic derivative
Cr         2 freeze to spec'd linearization energy
Cr         3 float to natural center of band (C parameter)
Cr  idu    identifies l-channels with Hubbard U (LDA+U)
Cr  idxdn  idxdn
Cr  jh     LDA+U J parameters for each l-channel
Cr  kmxt   k cutoffs for tail augmentation expansion
Cr  kmxv   k-cutoff for 1-center projection of free-atom rho
Cr  lfoca  switch specifying treatment of core density
Cr  lmxa   l cutoff for augmentation expansion
Cr  lmxb   highest l in basis
Cr  lmxl   cutoff for local density
Cr  name   species name
Cr  ngcut  orbital-dependent G cutoffs (for nfp basis)
Cr  nr     nr for mesh
Cr  nxi    Number of energies in fit of free-atom density tails
Cr  orhoc  pointer to core density
Cr  p      log derivative parameter
Cr  pz     log derivative parameter, semicore states
Cr  q      starting q's (charges)
Cr  qc     core charge
Cr  rcfa   renormalization radius of free atom density, and width
Cr  rfoca  smoothing radius for frozen core overlap approx
Cr  rg     rsm for gaussians to fix multipole moments
Cr  rmt    augmentation radius
Cr  rs3    Lower bound to rsm for local orbital
Cr  rsma   rsm for augmentation expansion
Cr  rsmfa  rsm to fit free atom density
Cr  rsmv   rsm for one-center projection of free-atom density
Cr         In TCF, smoothing radius to calculate exchange
Cr         (artificially smoothed density)
      type s_spec
      real(8), allocatable :: rv_a_orhoc(:)
      character(8)   ::   name !string
      real(8)   ::   z !atomic number
      real(8)   ::   mass
      real(8)   ::   rmt
      real(8)   ::   rsmfa
      real(8)   ::   rsma
      real(8)   ::   rg
      integer   ::   lmxa
      integer   ::   lmxl
      integer   ::   kmxt
      real(8)   ::   rsmv
      character(8)   ::   coreh !string
      real(8)   ::   coreq(2)
      real(8)   ::   a
      integer   ::   nr
      real(8)   ::   eref
      integer   ::   lfoca
      real(8)   ::   ctail
      real(8)   ::   etail
      real(8)   ::   stc !core kinetic energy
      integer   ::   lmxb
      real(8)   ::   rfoca
      integer   ::   nxi
      real(8)   ::   qc
      real(8)   ::   eh3
      real(8)   ::   rs3
      real(8)   ::   vmtz!Asymptotic potential for fitting functions at rmt
      integer   ::   kmxv
      real(8)   ::   rcfa(2)
      real(8)   ::   p(20)
      real(8)   ::   q(20)
      integer   ::   idmod(10)
      integer   ::   idxdn(30)
      real(8)   ::   exi(10)
      integer   ::   ngcut(30)
      real(8)   ::   chfa(20)
      real(8)   ::   pz(20)
      integer   ::   idu(4)
      real(8)   ::   uh(4)!Hubbard U
      real(8)   ::   jh(4)
      integer :: nmcore !jun2012takao
      end type s_spec
!      
      type s_site
      integer  :: iantiferro !may2015
      character(8)   ::   clabel !labelling class
      integer   ::   class !class index
      real(8)   ::   force(3) !force   force (dynamics)
      real(8) , allocatable ::  rv_a_ov0 (:)!pointer to potential that defines wave functions
      real(8) , allocatable ::  rv_a_ov1 (:)!pointer to spherical part of MT potential
      real(8)   ::   pnu(20)!log derivative parameter
      real(8)   ::   pos(3) !true coordinates of atom
      real(8)   ::   pos0(3) !coordinates of atom before distortion
      real(8)   ::   pz(20) !log derivative parameter, semicore states
      integer   ::   relax(3) !(dynamics) flags which coordinates to relax
      integer   ::   spec !species index
      real(8)   ::   vel(3) !velocity
      real(8)   ::   vshft !constant potential shift for this site
      end type s_site
      end module m_struc_def
