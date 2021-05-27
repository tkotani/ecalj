!! ===   Repeat loop for printout and goto 99 ===
!!   jsp=isp in the collinear case; jsp=1 in the noncollinear
!!     Thus jsp should be used in place of isp
!     !     isp serves as a flag for the noncollinear case
!! block6      
      goto99=.false.
      do  iq = 1, nkp
         qp=qplist(:,iq)
         do isp = 1, nspx
            jsp = isp
            if(iprint()>20.and.ix99==1) then
               write(*,'(" bndfp: kpt",i5," of",i5, " k isp=",3f8.4,i2," ndimh nev=",2i5)')
     &           iq,nkp,qp,jsp,ndimhx_(iq),nevls(iq,jsp)
            endif   
            if(iprint()>=35.and.ix99==1) then
               write(stdo,"(9f8.4)") (evlall(i,jsp,iq), i=1,nevls(iq,jsp))
            endif
            evl(1:ndhamx,jsp) = evlall(1:ndhamx,jsp,iq)
            if (mod(iq,10) .ne. 1) call pshpr(iprint()-6)
            if(debug) call info5(30,0,0,' bndfp 299loop:  kpt %i of %i, k=%3:2,5;5d', iq,nkp,qp,0,0)
            if(debug) print *,'eeeee11111111 ',procid,numprocs
            if(numprocs==1) then
               pidorigin = master
            else
               pidorigin=-999
               do idx=0,numprocs-1
                  if (kpproc(idx)<=iq .and. iq<= kpproc(idx+1)-1) then
                     pidorigin=idx
                     exit
                  endif
               enddo   
            endif
            nev_iq    = nev_(iq) 
            if (mod(iq,10) .ne. 1) call poppr
            ebot = dmin1(ebot,evl(1,jsp))
            i = max(1,nint(qval-qbg)/(3-nspc))
            evtop = max(evtop,evl(i,jsp))
            ecbot = min(ecbot,evl(i+1,jsp))
            if (lmet .eq. 0 .and. iq .eq. 1 .and. jsp .eq. 1) ef0 = evtop
            if(debug) print *,'eeeee44444444444 plbnd=',plbnd
            
            if (plbnd .eq. 0) then
               if (ipr.ge.10 .and. iq.eq.1 .and. ipr.gt.0)
     .           write (stdl,"('fp evl',8f8.4)") (evl(i,jsp),i=1,nev_iq)
               if (lwtkb .ne. -1) then ! .and. .not. lwndow) then
                  if (iq .eq. 1 .and. jsp .eq. nsp ) then !
!! Energy-Window reset. dosw=(emin,emax)
                     ef00 = ef0
                     call fixef0(qval-qbg,jsp,1,nev_iq,ndhamx,evl,dosw,ef0)
                     if(master_mpi) then
                        open(newunit=ifid,file='ewindow.'//trim(sname)) 
                        write(ifid,"(3d23.15)") dosw(1:2),ef0
                        close(ifid) !call fclose(ifid)
                     endif
                     if (jsp .eq. 2 .and. ef00 .ne. ef0 .and.
     .                    lwtkb .eq. 0 .and. lmet .gt. 0 .and. lrout .ne. 0) then
                        if (master_mpi) call info0(10,1,1,
     .                       ' ... Fermi level reset in second spin'//
     .                       ' channel ... restart band pass')
                        goto99=.true.
                        exit    !this is the case of make co test at ecalj/TestInstall/
                     endif
                  endif

C     Check for cases when nevmx is too small : i=2 => fatal error
                  i = 0
                  if (nevmx.ge.0 .and. lmet .ne. 0) then
                     dum = evl(max(nev_iq,1),jsp)
                     if (.not. ltet .and. ef0+5*bz_w .gt. dum) i=2
                     if (lmet.eq.4 .and. ef0+def+5*bz_w .gt.dum)i=2
                  endif
                  if(i .eq. 2) then
                     write(stdo,"(a,f13.5,a, f13.5)")
     &                 'evl(nev='//trim(i2char(nev_iq))//')=',
     &                 evl(max(nev_iq,1),jsp),' but ef0=',ef0
                     call rx('bndfp:... restart with larger nevmx: bndfp')
                  endif
               endif
            endif
         enddo                  ! end second loop over isp
      enddo                     !end second loop over iq
      call mpi_barrier(MPI_comm_world,ierr) !takao added I think 'safer' by this barrier. jan2015
      if(goto99) goto 99
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
