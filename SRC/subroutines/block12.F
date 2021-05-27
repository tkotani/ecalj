      goto99=.false.
      call mpi_barrier(MPI_comm_world,ierr)
      if (master_mpi) ierr=unlink('MagField') !delete
!! --- BZ integration for fermi level, band sum and qp weights ---
      if(lmet .ge. 0 .and. (lmet .ne. 4 .or. ltet)) then
         dosrng = 8
         if (bz_n .lt. 0) dosrng = 16
         if( debug) print *, 'tttttttttt 99999999 6 call bzwtsf'
         if(bz_fsmommethod == 1) then !takao dec2010
!     ! vnow june22013 !vnow !(in Ry) contains magnetic field
!     ! For eigenvalus, add  -vnow/2 for isp=1, and +vnow/2 for isp=2.
            call bzwtsf2 ( ndham , ndham , nsp , nspc , nkabc ( 1 ) , nkabc 
     .           ( 2 ) , nkabc ( 3 ) , nkp , ntet , iv_a_oidtet , qval - qbg , 
     .           fsmom , lmet.ne.0 , ltet , bz_n , ndos , bz_w 
     .            , dosrng , rv_a_owtkp , evlall ,  lswtk , rv_a_oswtk 
     .           , eferm , sev , rv_a_owtkb , sumqv ( 1 , 2 ) , lwtkb ,lfill,vnow)
         else
!     ! vnow june22013
            call bzwtsf ( ndham , ndham , nsp , nspc , nkabc ( 1 ) , nkabc 
     .           ( 2 ) , nkabc ( 3 ) , nkp , ntet , iv_a_oidtet , qval - qbg , 
     .           fsmom , lmet.ne.0 , ltet , bz_n , ndos , bz_w 
     .           , dosrng , rv_a_owtkp , evlall , lswtk , rv_a_oswtk 
     .           , eferm , sev , rv_a_owtkb , sumqv ( 1 , 2 ) , lwtkb ,lfill, vnow)
         endif
!     ! june2013 magfield is added
         if(fsmom/=NULLR.and.master_mpi) then
            open(newunit=ifimag,file='MagField',status='unknown')
            write(ifimag,"(d23.16,' !(in Ry) -vnow/2 for isp=1, +vnow/2 for isp=2')")vnow
            close(ifimag)
         endif
!     !         Store val charge & magnetic moment in sumqv(1..2)
         if (lmet .ne. 4) then
            sumqv(1,1) = sumqv(1,2)
            sumqv(2,1) = sumqv(2,2)
         endif
         if (lmet .ne. 4) then
            ef0 = eferm
            bz_ef=ef0
         endif
         if (lmet .gt. 0) then
            if (master_mpi) then
               ifi=ifile_handle()
               open(ifi,file='wkp.'//trim(sname),form='unformatted') 
               i = iobzwt ( 0 , ndhamx , nkp , nspx , eferm , rv_a_owtkb , -ifi )
               close(ifi) !call fclr('wkp',ifi)
            endif
         endif
         goto99=.false.
         if (lwtkb .eq. -1 .and. lrout .gt. 0) then
            call info0(20,0,0,' Start second band pass ...')
            lwtkb = 1
c            if (nspc .eq. 2) lswtk = 1 !commnet out nov2015  right?
            goto99=.true.
         endif
         if (lwtkb .eq. 2 .and. lrout .gt. 0) then
            call info0(20,0,0,' New pass with constrained weights ...')
            goto99=.true.
         endif
         call mpi_barrier(MPI_comm_world,ierr)  
         if(goto99) goto 99
      endif
!! == goto99 ix99 loop over =============================
