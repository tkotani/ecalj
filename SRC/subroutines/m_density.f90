module m_density
  use m_struc_def,only: s_rv1
  type(s_rv1), allocatable :: v1pot(:),v0pot(:)
  type(s_rv1), allocatable ::  orhoat(:,:) ! atom densities  i/o
  complex(8) , allocatable ::  osmrho(:) ! smooth density  i/o
  real(8),allocatable,target ::   pnuall(:,:,:) ! log derivative parameter
  real(8),allocatable,target ::   pnzall(:,:,:)  ! log derivative parameter for LO
end module m_density
!o   orhoat(1:3,ibas)%v : atomic density in 3-component form (true rho, smoothed rho, core rho)
!o   smrho :smoothed interstitial density
!o         :* for smrho = smoothed mesh density, smrho is complex and
!o         :  smrho = smrho(k1,k2,k3), k1,k2,k3 are in smpot

