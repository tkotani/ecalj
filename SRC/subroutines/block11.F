!!   For now, SO by site, not by class
      if (lwtkb .eq. 1 .and. lso .ne. 0) then
         allocate(ips_iv(nbas))
         do i_spackv=1,nbas
            ips_iv(i_spackv)=ssite ( i_spackv ) %spec 
         enddo
!! orbital moment
         call iorbtm ( sspec , ips_iv , nl , nl , nbas , nsp , orbtm_rv  )
         if (allocated(ips_iv)) deallocate(ips_iv)
      endif
      if (ipr.gt.0) write (stdl,"('nv nkp',i5,'  ebot',f9.4,'   qval',f10.4,'  qbg',f8.4,'  esmear',f8.4)")
     & nkp,ebot,qval,qbg,bz_w
!! --- Interpolate density to Fermi energy ---
      sev = sumev(1,1)
      if (lmet .eq. 4) then
         call mshn3p ( nbas , ssite , sspec , lmet , lrout , lfrce , qval 
     .     - qbg , ef0 , def , sumqv , sumev , n1 , n2 , n3 , k1 , k2 , 
     .     k3 , srout_zv , sv_p_oqkkl , frc , lrep )
C   ... Store val q & magnetic moment in sumqv(1) and sumqv(2)
         sumqv(2,1) = sumqv(1,1) - sumqv(1,2)
         sumqv(1,1) = sumqv(1,1) + sumqv(1,2)
C   ... Eigenvalue sum including entropy term
         sev = sumev(1,1) + sumev(2,1)
C   ... Remake sev,ef linearly interpolating tabulated sampling DOS
         sev00 = sev
         ef00  = ef0
         if (ldos .ne. 0) then
            call efldos ( qval , nsp , emin , emax , ndos , dos_rv , eferm, sev1 )
            sev   = sev1
            ef0  = eferm
         endif
         if (ipr .gt. 30 .and. ldos .ne. 0)
     .     write(stdo,"(' ipol:  sev=',f12.6,'   ef=',f12.6:/
     .     ' dos:   sev=',f12.6,'   ef=',f12.6/
     .     ' use:   sev=',f12.6,'   ef=',f12.6)") sev00,ef00,sev1,eferm,sev,ef0
         if (ipr.gt.0) write (stdl,"('nf EF:',3f9.5,'    EB:',3f12.5)") ef00,eferm,ef0,sev00,sev1,sev
         bz_dosw = dosw
         bz_ef=ef0
         bz_def=def
         if (lrep .eq. 1) then
            ef0 = -1
           write(stdo,*)'Input Fermi energy was too far off, repeat band pass'
            call mpi_barrier(MPI_comm_world,ierr)  
            goto 99
         endif
      endif
