module m_density !all are changing during iteration. Not protected!!!!
  use m_struc_def,only: s_rv1
  type(s_rv1), allocatable :: v0pot(:) !v0pot(ib)%v MT potential that defines wave functions
  type(s_rv1), allocatable :: v1pot(:) !v1pot(ib)%v  MT potential of spherical part
  type(s_rv1), allocatable ::  orhoat(:,:) 
!  orhoat(1:3,ibas)%v : atomic density in 3-component form (true rho, smoothed rho, core rho)
  complex(8) , allocatable ::  osmrho(:,:) ! smoothed interstitial density, smrho(k1,k2,k3,nsp)
  real(8),allocatable,target ::   pnuall(:,:,:) ! log derivative parameter
  real(8),allocatable,target ::   pnzall(:,:,:)  ! log derivative parameter for LocalOrbital
end module m_density
