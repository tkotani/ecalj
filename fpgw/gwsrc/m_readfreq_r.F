      module m_readfreq_r
      use m_genallcf_v3,only: niw
!     ! Readin WV.d (dimension file)
!     direct access files WVR and WVI which include W-V.
      integer,protected::nprecx,mrecl,nblochpmx,nwp,niwt, nqnum, nw_i,nw
      real(8),allocatable,protected:: freq_r(:)
      real(8),allocatable::freqx(:),freqw(:),wwx(:)

      contains
      subroutine readfreq_r()
      open(newunit=ifwd,file='WV.d')
      read(ifwd,*) nprecx,mrecl,nblochpmx,nwp,niwt, nqnum, nw_i
      write(6,"(' Readin WV.d =', 10i8)") nprecx,mrecl,nblochpmx,nwp,niwt, nqnum, nw_i
      close(ifwd)
      nw = nwp-1 
c      if(niwt /= niw) call rx( 'hsfp0_sc: wrong niw')
!     ! Energy mesh along real axis. Read 'freq_r' 
!     !  NOTE nw_i=-nw for non-timereversal case.
!     !       nw_i=0 for time-reversal case.
!     !  NOTE: We assume freq_r(i) == -freq_r(-i) in this code. feb2006
!     !  NOTE: this program assumes freq_r(iw)=freq_r(-iw). freq_r(iw <0) is redundant.
c     write(6,'("    niw nw dw   =",2i6,f13.6)') niw,nw,dw
      open(newunit=ififr,file='freq_r')
      read(ififr,*)nwxx
      if(nwxx/= nw+1) call rx( ' freq_r nw /=nw')
      allocate(freq_r(nw_i:nw)) !freq_r(0)=0d0
      do iw= nw_i,nw  
         read(ififr,*) freq_r(iw)
      enddo
      close(ififr)
      if(nw_i/=0) then
         if(nw/= -nw_i)        call rx( "sxcf_fal3_scz: nw/=-nw_i")
         if(freq_r(0)/=0d0)    call rx( "sxcf_fal3_scz: freq_r(0)/=0")
         if( sum(abs( freq_r(1:nw)+freq_r(-1:-nw:-1)))/=0) 
     &        call rx( "sxcf_fal3_scz: freq_r /= -freq_r")
      endif
!! Generate frequencies x between (0,1) and w=(1-x)/x for Gaussian integral along im axis
      if(niw/=0) then
        allocate(freqx(niw),freqw(niw),wwx(niw)) !,expa(niw))
        call freq01x  (niw,     !ua,
     o   freqx,freqw,wwx)       !,expa)
      endif
      end subroutine
      end module
