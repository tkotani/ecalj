!> Frequency mesh generator
module m_freq 
!! - OUTPUT
!!   - fhris :histgram bins to accumlate im part
!!   - freq_r: omega along real axis
!!   - freq_i: omega along imag axis
!!   - wiw: integration weight along im axis
!!   - npm: npm=1 means only positive omega;npm=2 means positive and negative omega.
!! - NOTE: change of frequency mesh defined here may destroy consistency or not. Need check
  use m_read_bzdata,only: dq_,qbz,nqbz
  use m_qbze,only:  nqbze,nqibze,qbze,qibze
  use m_genallcf_v3,only: niw_in=>niw,ecore,nctot,nspin,nband
  !use m_readhbe,only: nband
  implicit none
  public:: getfreq2, getfreq3, getfreq,freq01
  real(8),allocatable,protected,public::frhis(:),freq_r(:),freq_i(:),wiw(:),freqx(:),wx(:)
  integer,protected,public:: nwhis, npm, nw_i, nw, niw
  private
  real(8),private:: emin,emax,omg2max
contains
  subroutine getfreq3(lqall,epsmode,realomega,imagomega,ua,iprint) ! Get data set for m_freq. All arguments are input.
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
    logical:: realomega,imagomega,iprint,epsmode,qbzreg
    real(8):: omg2max,ua
    call findemaxmin(nband,qbze,nqbze,nspin, emax,emin)
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
    integer::ifif
    logical,optional:: npmtwo !! Added Aug2017 for hmagnon
    logical:: npm2
    logical,save:: done=.false.
    if(done) call rx('gerfreq is already done') !sanity check
    done =.true.
    nw=-99999 !for sanity check
    !! Histogram bin divisions
    !! We first accumulate Imaginary parts.
    !! Then it is K-K transformed to obtain real part.
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
       if (omg2 > wemax/2d0 ) then ! dw*(nw_input-3)) then !omg is in unit of Hartree
          nw=iw
          exit
       endif
    enddo
    !! document need to be fixed...
    ! nw+1 is how many points of real omega we use
    ! Screened Coulomb W(iw=0:nw) iw=0 corresponds omg=0
    ! maximum nw=nw2-1 because nwhis=nw2-1
    ! W is chosen from condition that frhis_m(nw-3)<dw*(nw_input-3) <frhis_m(nw-2).
    ! ere frhis_m(iw)= (freqr2(iw)+freqr2(iw+1))/2d0
    ! w was constructed such that omg=dw*(nw-2)> all relevant frequensies needed
    ! for correlation Coulomb Wc(omg),
    ! and one more point omg=dw*(nw-1) needed for extrapolation.
    ! Now, frhis_m(nw-1)> all relevent frequensies for Wc(omg)
    ! and one more point omg=frhis_m(nw) needed for extropolation
    !! Determine freq_r
    if(epsmode) nw  = nwhis-1
    allocate(freq_r(0:nw))
    freq_r(0)=0d0
    freq_r(1:nw)=(frhis(1:nw)+frhis(2:nw+1))/2d0
    !! Timereversal=F is implimented only for tetra=T and sergeyv=T
    !! nw_i and npm
    npm=1
    nw_i=0
    npm2=.false.
    if(present(npmtwo))then
       if(npmtwo) npm2= .TRUE.
    endif   
    if( .NOT. timereversal() .OR. npm2)  then
       write(6,"('TimeReversal off mode')")
       npm=2
       nw_i=-nw
       !  if(.not.tetra)   call rx( ' tetra=T for timereversal=off')
    endif
    write(6,*)'Timereversal=',Timereversal()
    !! Determine freq_i  : gaussian frequencies x between (0,1) and w=(1-x)/x
    if (imagomega) then
       write(6,*)' freqimg: niw =',niw
       allocate( freq_i(niw) ,freqx(niw),wx(niw),expa(niw),wiw(niw))
       call freq01 (niw,ua, freqx,freq_i,wx,expa)
       wiw=wx/(2d0*pi*freqx**2)
       deallocate(wx,expa) 
    endif
    if(onceww(1)) then !plot frhis
       write(6,*)' we set frhis nwhis =',nwhis 
       write(6,*)' --- Frequency bins to accumulate Im part  (a.u.) are ---- '
       write(6,"(' ihis Init  End=', i5,2f18.11)") (ihis,frhis(ihis),frhis(ihis+1),ihis=1,nwhis)
    endif
  end subroutine getfreq
  subroutine freq01 (nx,ua,  freqx,freqw,wx,expa)!Generates a gaussian point x between (0,1) and w = (1-x)/x 
    ! and the weights in x
    ! also generates expa = exp(-ua^2 w^2)
    ! nx    = no. gaussian points
    ! ua   = s.o.
    !o freqx = gaussian points
    !o freqw = (1-x)/x
    !o wx    = weights of gaussian points x
    ! expa  = s.o.
    ! originally 92.02.27 Ferdi.Aryasetiawan
    implicit real*8 (a-h,o-z)
    integer:: nx,ix
    real(8):: freqx(nx),freqw(nx),wx(nx),expa(nx)
    call gauss   (nx,0.d0,1.d0,freqx,wx)! generate gaussian points
    ua2    = ua*ua
    freqw  = (1d0 - freqx) / freqx
    expa   = exp(-ua2*freqw**2)
    !   write(6,"(' --- freq01x:  ix    x    freqw(a.u.)---')")
    !   do      ix = 1,nx
    !      freqw(ix)  = (1.d0 - freqx(ix)) / freqx(ix)
    !      write(6,"('            ',i4,2f9.4)") ix,freqx(ix),freqw(ix)
    !   end do
  end subroutine freq01
end module m_freq

subroutine freq01x (nx,    freqx,freqw,wx)
  implicit real*8 (a-h,o-z)
  real(8):: freqx(nx),freqw(nx),wx(nx),expa(nx)
  integer:: nx,ix
  call gauss (nx,0.d0,1.d0,freqx,wx)
  freqw  = (1d0 - freqx) / freqx
end subroutine freq01x

subroutine getwemax(lqall,wemax)!> In order to get |e_ip-ef| on real space integration ! too complicated ---> need to fix in future
  use m_read_bzdata,only: read_bzdata, nqibz,qibz,ginv,qbz,nqbz,wibz
  use m_genallcf_v3,only: nspin, konf,z,nl,natom,iclass,esmr,deltaw,nband,laf !,anfcond
  use m_keyvalue,only:getkeyvalue
  use m_readeigen,only: readeval !init* is alreaday called.
  use m_ReadEfermi,only: ef !ef is set at main routine
!  use m_anf,only:
  implicit none
  logical,intent(in):: lqall
  real(8),intent(out):: wemax
  integer:: ntq,i,nq,ifqpnt,ret,ipx1,itx1,ipx2,itx2
  integer,allocatable:: itq(:)
  real(8),allocatable:: q(:,:),eqt(:) !,q0p(:,:)
  logical ::  legas = .false., tetra,lqallxxx
  real(8):: omegav,omegac,eee,efnew,emaxv,eminc,ffac,we,valn
  integer:: ifief,ib,ip,iq,iqall,it,k,is,nspinmx
  lqallxxx=.true.
!  call anfcond()
  if( .NOT. lqall) then
     call getkeyvalue("GWinput","<QPNT>",unit=ifqpnt,status=ret)
     call readx   (ifqpnt,10)
!     iaf=-1
     read (ifqpnt,*) iqall!,iaf
     if (iqall == 1) then
        lqallxxx = .true.
     else
        lqallxxx = .false.
        call readx   (ifqpnt,100)
        read (ifqpnt,*) ntq
        allocate( itq(ntq) )
        read (ifqpnt,*) (itq(i),i=1,ntq)
        call readx   (ifqpnt,100)
        read (ifqpnt,*) nq
        allocate(q(3,nq))
        do       k = 1,nq
           read (ifqpnt,*) i,q(1,k),q(2,k),q(3,k)
           write(6,'(i3,3f13.6)') i,q(1,k),q(2,k),q(3,k)
        enddo
     endif
     close(ifqpnt) !Walter fixed this Jul2008
  endif
  if(lqallxxx) then
     ntq = nband
     allocate (itq(ntq))
     do i = 1, ntq
        itq(i) = i
     enddo
     nq  = nqibz
     allocate(q(3,nq))
     q = qibz(1:3,1:nqibz) !call dcopy   (3*nqibz,qibz,1,q,1)
  endif
  nspinmx = nspin
  if (laf) nspinmx =1
  !! for 1shot GW deltaw id for d\Sigma/d_omega
  allocate(eqt(nband))
  omegav =  0d0 !1d99 fixed oct.2003 takao
  omegac =  0d0 !-1d99
  do is = 1,nspinmx
     do ip = 1,nq
        eqt = readeval(q(1,ip),is)
        do it=1,ntq
           eee = eqt(itq(it)) - 2d0*deltaw  !2.d0*(-1d0-shtw)*deltaw
           !            write(6,*)' is ip it eee=',eee,eqt(itq(it))
           if(eee>=1d20-1d10) cycle !takao jun2009
           if( eee <ef .AND. eee< omegav )  then
              omegav = eee
              ipx1 = ip
              itx1 = it
           endif
           eee = eqt(itq(it)) + 2d0*deltaw  !+ 2.d0*(1d0-shtw)*deltaw
           !            write(6,*)' is ip it eee=',eee,eqt(itq(it))
           if( eee >ef .AND. eee> omegac )  then
              omegac = eee
              ipx2 = ip
              itx2 = it
           endif
        enddo
     enddo
  enddo
  !! for Gaussian smearing   !      if(GaussSmear()) then
  ffac=10d0 !This is OK?
  !      else   !        ffac=0.5d0  !      endif
  emaxv =  0d0 !-1d99 fixed oct.2003 takao
  eminc =  0d0 !1d99
  do is = 1, nspinmx
     do iq = 1, nqbz
        eqt= readeval(qbz(1,iq),is)
        do ib=1,nband
           eee = eqt(ib)
           if( eee <ef+ffac*esmr .AND. eee>emaxv) emaxv = eee
           if( eee >ef-ffac*esmr .AND. eee<eminc) eminc = eee
        enddo
     enddo
  enddo
  deallocate(eqt)
  we  = max(abs(emaxv - ef), abs(omegav-ef),abs(omegac- ef) , abs(eminc-ef) )
  wemax= we+ffac*esmr
end subroutine getwemax
subroutine findemaxmin(nband,qbz,nqbz,nspin, emax,emin)
  use m_readeigen, only: readeval
  implicit none
  integer :: nband,nqbz,nspin,isp,kx,i 
  real(8)::emax,emin,qbz(3,nqbz),eee
  real(8),allocatable:: ekxxx(:,:,:)
  allocate( ekxxx(nband,nqbz,nspin))
  do isp =1, nspin
     do kx = 1, nqbz
        ekxxx(1:nband,kx,isp) = readeval(qbz(:,kx), isp)
     enddo
  enddo
  Emax = maxval(ekxxx,mask=ekxxx<1d9) !not eee<1d9 corresponds to 1d20 for padding in lmf2gw.F and sugw.Fago
  Emin = minval(ekxxx)
  deallocate(ekxxx)
end subroutine findemaxmin
