module m_density
  use m_struc_def,only: s_rv1
  type(s_rv1), allocatable ::  orhoat(:,:) ! atom densities  i/o
  complex(8) , allocatable ::  osmrho(:) ! smooth density  i/o
end module m_density
!o   orhoat(1:3,ibas)%v : atomic density in 3-component form (true rho, smoothed rho, core rho)
!o   smrho :smoothed interstitial density
!o         :* for smrho = smoothed mesh density, smrho is complex and
!o         :  smrho = smrho(k1,k2,k3), k1,k2,k3 are in smpot

