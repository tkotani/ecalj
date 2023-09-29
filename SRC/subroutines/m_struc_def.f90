!>definision of structures. pointer array
module m_struc_def 
  integer,parameter::  n0=10,nkap0=3
  public s_rv1, s_nv2, s_cv1, s_cv2,s_cv3,s_cv4, s_sblock, s_spec,s_rv5,s_rv4,s_cv5
  private
  type s_rv1
     real(8),allocatable:: v(:)
  end type s_rv1
  type s_rv4
     real(8),allocatable:: v(:,:,:,:)
  end type s_rv4
  type s_rv5
     real(8),allocatable:: v(:,:,:,:,:)
  end type s_rv5
  type s_cv5
     complex(8),allocatable:: cv(:,:,:,:,:)
  end type s_cv5
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
  type s_spec  
     real(8),allocatable :: rv_a_orhoc(:) !pointer to core density
     real(8):: rsmfa !rsm to fit free atom density
     real(8):: ctail !coefficients to fit of free-atom core tail by unsm. Hankel
     real(8):: etail !energy to fit of free-atom core tail
     real(8):: stc   !core kinetic energy
     real(8):: qc    !core charge
     integer:: nxi    ! Number of energies in fit of free-atom density tails
     real(8):: exi(n0)! Hankel energies for fit to c.d.; fit to free-atom density tails.
     real(8):: chfa(n0,2) ! coefficients to fit of free-atom density tails
  end type s_spec
end module m_struc_def
!     real(8):: z     !atomic number
!     real(8):: rg    !rsm for gaussians to fix multipole moments
!     integer:: lmxa  !  l cutoff for augmentation expansion
!     integer:: lmxl  !  l cutoff for local density and potential
!     integer:: lmxb  !  highest l in basis
!     integer:: kmxt  !  k cutoffs for tail augmentation expansion of P_kl
!     real(8):: rsmv  !  rsmv  =rmt*.5d0 in defspc from m_lmfinit. smoothing radius of gaussian
!     real(8):: a  !a for mesh
!     integer:: lfoca ! switch specifying treatment of core density
!     real(8):: rfoca ! smoothing radius for frozen core overlap approx
! followings are given at lmfa. And read by rdovfa->iofa, I think.     
