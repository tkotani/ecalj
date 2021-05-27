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
c            if(iprint()>20.and.ix99==1) then
c               write(*,'(" bndfp: kpt",i5," of",i5, " k isp=",3f8.4,i2," ndimh nev=",2i5)')
c     &           iq,nkp,qp,jsp,ndimhx_(iq),nevls(iq,jsp)
c            endif   
c            if(iprint()>=35.and.ix99==1) then
c               write(stdo,"(9f8.4)") (evlall(i,jsp,iq), i=1,nevls(iq,jsp))
c            endif
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
