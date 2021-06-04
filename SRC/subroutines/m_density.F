      module m_density
      use m_struc_def,only: s_rv1
      type(s_rv1), allocatable ::  orhoat(:,:) ! atom densities  i/o
      complex(8) , allocatable ::  osmrho(:) ! smooth density  i/o
      real(8):: eh,eks,utot,valves,cpnves,rhoexc,rhovxc,sumec,
     &         sumtc,xcore,valvef,amom,sumt0,sumev !rhosig,
      end
Co   orhoat(1:3,ibas)%v : atomic density in 3-component form (true rho, smoothed rho, core rho)
Co   smrho :smoothed interstitial density
Co         :* for smrho = smoothed mesh density, smrho is complex and
Co         :  smrho = smrho(k1,k2,k3), k1,k2,k3 are in smpot
