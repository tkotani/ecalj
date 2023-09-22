!>read freq_r ! Readin WV.d (dimension file)  direct access files WVR and WVI which include W-V.
module m_readfreq_r
  use m_genallcf_v3,only: niw
  use m_readgwinput,only: ua_
  implicit none
  integer:: mrecl
  integer,protected::nprecx,nblochpmx,nwp,niwt, nqnum, nw_i,nw,npm
  real(8),allocatable,protected:: freq_r(:)
  real(8),protected,allocatable::freqx(:),freqw(:),wwx(:),expa_(:)
contains
  subroutine readfreq_r()
    integer:: ififr,ifwd,ix,nwxx,iw
    open(newunit=ifwd,file='WV.d')
    read(ifwd,*) nprecx,mrecl,nblochpmx,nwp,niwt, nqnum, nw_i
    write(6,"(' Readin WV.d =', 10i8)") nprecx,mrecl,nblochpmx,nwp,niwt, nqnum, nw_i
    close(ifwd)
    nw = nwp-1
    !      if(niwt /= niw) call rx( 'hsfp0_sc: wrong niw')
    !     ! Energy mesh along real axis. Read 'freq_r'
    !     !  NOTE nw_i=-nw for non-timereversal case.
    !     !       nw_i=0 for time-reversal case.
    !     !  NOTE: We assume freq_r(i) == -freq_r(-i) in this code. feb2006
    !     !  NOTE: this program assumes freq_r(iw)=freq_r(-iw). freq_r(iw <0) is redundant.
    !     write(6,'("    niw nw dw   =",2i6,f13.6)') niw,nw,dw
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
       if( sum(abs( freq_r(1:nw)+freq_r(-1:-nw:-1)))/=0) &
            call rx( "sxcf_fal3_scz: freq_r /= -freq_r")
    endif
    !! Generate frequencies x between (0,1) and w=(1-x)/x for Gaussian integral along im axis
    if(niw/=0) then
       allocate(freqx(niw),freqw(niw),wwx(niw),expa_(niw))
       call gauss   (niw,0d0,1d0,freqx,wwx) !!w = 1/(1+x) !    generate gaussian points
       write(6,"(' --- readfreq:  ix    x    freqw(a.u.)---')")
       do ix = 1,niw
          freqw(ix)= (1.d0 - freqx(ix)) / freqx(ix)
          write(6,"('            ',i4,2f9.4)") ix,freqx(ix),freqw(ix)
          expa_(ix)= dexp(-ua_**2*freqw(ix)*freqw(ix))
       end do
    endif
    npm = 1                   ! npm=1    Timeveversal case
    if(nw_i/=0) npm = 2       ! npm=2 No TimeReversal case. Need negative energy part of W(omega)
  end subroutine readfreq_r
end module m_readfreq_r
