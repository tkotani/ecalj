!! explanations are old, some are obsolate
module m_struc_def
  integer,parameter::  n0=10,nkap0=3
  !      use m_lmfinit,only:n0,nkap0

  public s_rv1, s_nv2, s_cv1, s_cv2,s_cv3,s_cv4, s_spec, s_site, s_sblock
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


  !! --- explanation of s_spec ---
  !r  idmod  idmol(l) controls how linearization energy is
  !r         determined for an orbital of angular momentum l
  !r         0 float to center of gravity of occupied part of band
  !r         1 freeze to spec'd logarithmic derivative
  !r         2 freeze to spec'd linearization energy
  !r         3 float to natural center of band (C parameter)
  ! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

  type s_spec               !I think all of them are fixed during iteration cycle of lmf-MPIK
     real(8), allocatable :: rv_a_orhoc(:) !pointer to core density
     real(8)   ::   z !atomic number
     real(8)   ::   rmt !augmentation radius
     real(8)   ::   rsmfa !rsm to fit free atom density
     real(8)   ::   rsma !rsm for augmentation expansion
     real(8)   ::   rg !rsm for gaussians to fix multipole moments
     integer ::     lmxa !  l cutoff for augmentation expansion
     integer ::     lmxl !cutoff for local density
     integer ::     kmxt !  k cutoffs for tail augmentation expansion
     real(8) ::   rsmv !  rsmv  =rmt*.5d0 in defspc from m_lmfinit. smoothing radius of gaussian
     character(8)   ::   coreh !core hole channel (char: eg '1s')
     real(8)   :: coreq(2)!coreq(1)=charge in core hole channel;coreq(2):moment in core hole channel (nsp=2 only)
     real(8)   ::   a !a for mesh
     integer   ::   nr !nr for mesh
     real(8)   ::   eref!reference energy
     integer   ::   lfoca !switch specifying treatment of core density
     real(8)   ::   ctail !coefficients to fit of free-atom core tail by unsm. Hankel
     real(8)   ::   etail !energy to fit of free-atom core tail
     real(8)   ::   stc !core kinetic energy
     integer   ::   lmxb !  highest l in basis
     real(8)   ::   rfoca ! smoothing radius for frozen core overlap approx
     integer   ::   nxi !Number of energies in fit of free-atom density tails
     real(8)   ::   qc !core charge
     real(8)   ::   eh3!sm Hankel energy for high local orbitals
     real(8)   ::   rs3 !Lower bound to rsm for local orbital
     integer   ::   kmxv !  k-cutoff for 1-center projection of free-atom rho
     real(8)   ::   rcfa(2) !renormalization radius of free atom density, and width
     real(8)   ::   q(n0,2)  !starting q's (charges)
     real(8)   ::   exi(n0) !Hankel energies for fit to c.d.;  For free atoms, fit to free-atom density tails.
     real(8)   ::   chfa(n0,2) ! coefficients to fit of free-atom density tails
!     integer   ::   idu(4) !identifies l-channels with Hubbard U (LDA+U)
!     real(8)   ::   uh(4) !Hubbard U
!     real(8)   ::   jh(4) !LDA+U J parameters for each l-channel
     integer   :: nmcore !jun2012takao
     real(8)   ::  p(n0) !log derivative for spec taken from ctrl file !shown by >lmfa si |grep conf
     real(8)   ::   pz(n0) !log derivative for spec taken from ctrl file
!     character(8)  ::   name !species name
!     real(8)   ::   vmtz!Asymptotic potential for fitting functions at rmt
!     integer   ::   ngcut(n0,nkap0) !orbital-dependent G cutoffs (for nfp basis)-->m_sugcut
!     integer   ::   idmod(n0) !see m_lmfinit.F
  end type s_spec

  type s_site
     integer :: iantiferro  !fixed during lmf-MPIK.
     integer ::   spec     !fixed. species index
     integer ::   class    !fixed. class index
     integer ::   relax(3) !fixed.(dynamics) flags which coordinates to relax

     real(8) ::   pos(3)  !fixed during do 1000 in lmfp.F  Coordinates of atom
     real(8) ::   pos0(3) ! atomic pos in previous loop of do 2000 in lmfp.F smshit (for MD)

     real(8) ::   force(3)  ! Force
     real(8) , allocatable ::  rv_a_ov0(:)! pointer to potential that defines wave functions
     real(8) , allocatable ::  rv_a_ov1(:)! pointer to spherical part of MT potential
     real(8) ::   pnu(n0,2) ! log derivative parameter
     real(8) ::   pz(n0,2)  ! log derivative parameter for LO
  end type s_site
end module m_struc_def
