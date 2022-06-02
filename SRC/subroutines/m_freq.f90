! Frequency mesh generator
!! - OUTPUT
!!   - fhris :histgram bins to accumlate im part
!!   - freq_r: omega along real axis
!!   - freq_i: omega along imag axis
!!   - wiw: integration weight along im axis
!!   - npm: npm=1 means only positive omega;npm=2 means positive and negative omega.
!! - NOTE: change of frequency mesh defined here may destroy consistency or not. Need check
module m_freq
  use m_read_bzdata,only: &
       dq_,qbz,nqbz
  use m_qbze,only: &
       nqbze,nqibze,qbze,qibze
  use m_genallcf_v3,only: niw_in=>niw,ecore,nctot,nspin
  use m_readhbe,only: nband
  implicit none
  !!--------------
  public:: getfreq2, getfreq3, getfreq
  real(8),allocatable,protected,public:: frhis(:),freq_r(:),freq_i(:),wiw(:), freqx(:),wx(:)
  integer,protected,public:: nwhis, npm, nw_i, nw, niw
  !!
  private
  real(8),private:: emin,emax,omg2max

contains

  !> Get data set for m_freq. All arguments are input.
  subroutine getfreq3(lqall,epsmode,realomega,imagomega,ua,iprint)!,tetra
    intent(in)::        lqall,epsmode,realomega,imagomega,ua,iprint
    real(8):: wemax
    integer:: iq
    !      logical,optional:: npmtwo !! Added Aug2017 for hmagnon
    !      integer:: niw !,nw_input
    logical:: realomega,imagomega,iprint,epsmode,lqall
    real(8):: omg2max,ua
    !! We get frhis,freq_r,freq_i, nwhis,nw,npm,wiw  by getfreq
    call findemaxmin(nband,qbze,nqbze,nspin, emax,emin)
    if (nctot > 0) Emin=minval(ecore(:,1:nspin))
    omg2max = (Emax-Emin)*.5d0+.2d0
    ! (in Hartree) covers all relevant omega, +.2 for margin
    if(iprint) write(6,"(' emin emax omega2max=',3f13.5)") emin, emax, omg2max
    wemax=1d10!dummy
    if( .NOT. epsmode) call getwemax(lqall,wemax) !wemax is to determine nw !real axis divisions
    niw=niw_in
    if( .NOT. imagomega) niw=1  !dummy
    call getfreq(epsmode,realomega,imagomega,omg2max,wemax,niw,ua)!tetra,
  end subroutine getfreq3
  !----------------
  subroutine getfreq2(epsmode,realomega,imagomega,ua,iprint,npmtwo) !,tetra
    intent(in)::        epsmode,realomega,imagomega,ua,iprint,npmtwo
    real(8)::emax2,emin2,wemax
    real(8),allocatable:: qbz2(:,:)
    integer:: iq
    logical,optional:: npmtwo !! Added Aug2017 for hmagnon
    !      integer:: niw !,nw_input
    logical:: realomega,imagomega,iprint,epsmode,qbzreg
    real(8):: omg2max,ua
    call Findemaxmin(nband,qbze,nqbze,nspin, emax,emin)
    if( .NOT. qbzreg()) then
       allocate(qbz2(3,nqbz))
       do iq=1,nqbz
          qbz2(:,iq)=qbz(:,iq)+dq_
       enddo
       call Findemaxmin(nband,qbz2,nqbz,nspin ,emax2,emin2)
       emax=max(emax,emax2)
       emin=min(emin,emin2)
       deallocate(qbz2)
    endif
    if(nctot > 0) Emin=minval(ecore(:,1:nspin))
    omg2max = (Emax-Emin)*.5d0+.2d0
    ! (in Hartree) covers all relevant omega, +.2 for margin
    call Getwemax(.true.,wemax) !wemax is to determine nw !real axis divisions
    niw=niw_in
    call Getfreq(.false.,realomega,imagomega,omg2max,wemax,niw,ua)!,npmtwo,tetra,
    if(iprint) write(6,"(' wemax=  ',f13.4)") wemax
    if(iprint) write(6,"(' emin emax omega2max=',3f13.5)") emin, emax, omg2max
  end subroutine getfreq2
  !----------------
  subroutine Getfreq(epsmode,realomega,imagomega,omg2max,wemax,niw,ua,npmtwo) !,tetra
    use m_keyvalue,only: getkeyvalue
    intent(in)::       epsmode,realomega,imagomega,omg2max,wemax,niw,ua,npmtwo
    integer:: niw !,nw_input
    logical:: realomega,imagomega,epsmode
    real(8):: omg2max,ua
    real(8),allocatable:: expa(:)
    logical:: timereversal,onceww
    integer:: nw2,iw,ihis
    real(8)::omg_c,dw,omg2,wemax
    real(8), allocatable :: freqr2(:)  ,frhis_tmp(:)
    real(8)::  pi = 4d0*datan(1d0), aa,bb,ratio,oratio,daa
    integer::nee,noo,ifif,ifile_handle
    logical,optional:: npmtwo !! Added Aug2017 for hmagnon
    logical:: npm2
    logical,save:: done=.false.
    if(done) call rx('gerfreq is already done') !sanity check
    done =.true.
    nw=-99999 !for sanity check
    !      nw = nw_input
    !! Histogram bin divisions
    !! We first accumulate Imaginary parts.
    !! Then it is K-K transformed to obtain real part.

    !      call getkeyvalue("GWinput","dw",dw )
    !      call getkeyvalue("GWinput","omg_c",omg_c )
    !      write(6,"('dw, omg_c= ',2f13.5)") dw, omg_c
    call getkeyvalue("GWinput","HistBin_ratio",oratio, default=1.03d0)
    call getkeyvalue("GWinput","HistBin_dw",dw, default=1d-5) !a.u.
    aa = oratio-1d0
    bb = dw/aa
    iw = 0d0
    do
       iw=iw+1
       if( bb*( exp(aa*(iw-1)) - 1d0 ) >omg2max+1d-6) exit
    enddo
    nwhis = iw+2 !+2 for margin. Necessary?
    allocate(frhis(1:nwhis+1))
    do iw = 1,nwhis+1
       frhis(iw) = bb*( exp(aa*(iw-1)) - 1d0 )
    enddo
    write(6,"('dw, omg_ratio, nwhis= ',d9.2,f13.5,i6)") dw, aa,nwhis

    !! Determine nw. Is this correct?
    do iw=3,nwhis
       omg2 = (frhis(iw-2)+frhis(iw-1))/2d0
       if (omg2 > wemax/2d0 ) then !>dw*(nw_input-3)) then !omg is in unit of Hartree
          nw=iw
          exit
       endif
    enddo
    !! document need to be fixed...
    !      nw=nw2-1      ! nw+1 is how many points of real omega we use
    ! for dressed coulomb line W(iw=0:nw) iw=0 corresponds omg=0
    ! maximum nw=nw2-1 because nwhis=nw2-1
    !! document need to be fixed...
    ! w is chosen from condition that frhis_m(nw-3)<dw*(nw_input-3) <frhis_m(nw-2).
    ! ere frhis_m(iw)= (freqr2(iw)+freqr2(iw+1))/2d0
    ! w was constructed such that omg=dw*(nw-2)> all relevant frequensies needed
    ! for correlation Coulomb Wc(omg),
    ! and one more point omg=dw*(nw-1) needed for extrapolation.
    ! Now, frhis_m(nw-1)> all relevent frequensies for Wc(omg)
    ! and one more point omg=frhis_m(nw) needed for extropolation
    !! Determine freq_r
    if(epsmode) then
       nw  = nwhis-1
    endif
    allocate(freq_r(0:nw))
    freq_r(0)=0d0
    do iw=1,nw
       freq_r(iw)=(frhis(iw)+frhis(iw+1))/2d0
    enddo
    !! Timereversal=F is implimented only for tetra=T and sergeyv=T
    !! nw_i and npm
    npm=1
    nw_i=0
    npm2=.false.
    if(present(npmtwo)) then
       if(npmtwo) npm2= .TRUE. 
    endif
    if( .NOT. timereversal() .OR. npm2)  then
       write(6,"('TimeReversal off mode')")
       npm=2
       nw_i=-nw
       !         if(.not.tetra)   call rx( ' tetra=T for timereversal=off')
    endif
    write(6,*)'Timereversal=',Timereversal()
    !! Determine freq_i  : gaussian frequencies x between (0,1) and w=(1-x)/x
    if (imagomega) then
       write(6,*)' freqimg: niw =',niw
       allocate( freq_i(niw) ,freqx(niw),wx(niw),expa(niw) )
       call freq01 (niw,ua, &
            freqx,freq_i,wx,expa)
       allocate(wiw(niw))
       do iw=1,niw
          wiw(iw)=wx(iw)/(2d0*pi*freqx(iw)*freqx(iw))
       enddo
       deallocate(wx,expa) !freqx,
    endif
    !! Plot frhis
    if(onceww(1)) then
       write(6,*)' we set frhis nwhis noo-->nee=',nwhis,noo,nee
       write(6,*)' --- Frequency bins to accumulate Im part  (a.u.) are ---- '
       do ihis= 1, nwhis !min(10,nwhis)
          write(6,"(' ihis Init  End=', i5,2f18.11)") ihis,frhis(ihis),frhis(ihis+1)
       enddo
    endif
  end subroutine getfreq
end module m_freq
