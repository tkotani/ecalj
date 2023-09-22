!> Get effective W(q=0,omega) for GW.
module m_w0w0i
  use m_ll,only:ll
  use m_llw,only: llw,llwI,wmuk,ngbq0,w4pmode
  use m_read_bzdata,only: q0i,nq0i, ixyz
  use m_genallcf_v3,only:  tpioa
  !! See the Friedlich's paper.
  implicit none

  public:: W0w0i
  complex(8),allocatable,protected,public :: w0(:),w0i(:),llmat(:,:)

  private
contains
  !----------------------------------------------------------
  subroutine finalizew4p()
    !      use m_w0w0i,only: llmat
    integer:: i,igb,ifw4p
    real(8):: qv(3,3)
    complex(8),allocatable:: wmu(:,:)
    !! Finalize w4phonon
    !! wmuk(ix)= matmul(wmu,qv) ==> wmu= matmul(wmuk,qvinv)
    do i=1,3
       qv(:,i)= tpioa*q0i(:,ixyz(i))
       qv(:,i)= qv(:,i)/sqrt(sum(qv(:,i)**2))
    enddo
    call matinv(3,qv)
    allocate( wmu(2:ngbq0,3) )
    do igb=2,ngbq0
       wmu(igb,:) =matmul(wmuk(igb,:),qv)
    enddo
    open(newunit=ifw4p,file='W4PHONON.HeadWing',form='unformatted')
    write(ifw4p) llmat(1:3,1:3),ngbq0 !for q~0
    write(ifw4p) wmu(2:ngbq0,1:3) !for q~0
    close(ifw4p)
    deallocate(wmu)
  end subroutine finalizew4p
  !-------------------------------------------------------------
  subroutine w0w0i(nw_i,nw,nq0i,niw,q0i) !llw,llwI,
    !! Get w0 and w0i (diagonal element at Gamma point) for given llw and llwi
    !! Outputs w0,w0i,llmat. See use m_w0w0i at the begining of this routine.
    use m_read_bzdata,only:  wbz,lxklm,dmlx,epinvq0i,epinv,wklm
    intent(in)::     nw_i,nw,nq0i,niw,q0i !llw,llwI,
    real(8)   :: q0i(1:3,1:nq0i)
    integer:: nw_i,nw,nq0i,nq0ix,niw,ifidmlx,i,ifw0w0i,ixc,nlxklm
    logical:: readw0w0itest
    complex(8):: llmat_dummy(3,3)
    write(6,*) 'w0w0i:'
    write(6,*)' ==== newaniso mode W(0) divergent part ==== '
    !! == W(0) divergent part ==  starting from llw(iw,iq0),llwI(iw,iq0)
    !! === <e|L|e> (eq.36 in Friedrich paper) is expanded in YL -->stored in llwyl. ===
    allocate(w0(nw_i:nw),w0i(niw),llmat(3,3))
    write(6,*)' goto getw0 nq0i epinvq0i=',nq0i,epinvq0i,wklm
    !! wbz(1) is the weight for q=0 = 1/(n1*n2*n3)
    !! llmat is added. July2016. llw is calculated at iw=0, when nw_i<=0
    call getw0(llw, nw_i,nw,nq0i,dmlx,epinvq0i,wklm,wbz(1), lxklm,  q0i,epinv,w0,  llmat)
    call getw0(llwI,1,niw  ,nq0i,dmlx,epinvq0i,wklm,wbz(1), lxklm,  q0i,epinv,w0i, llmat_dummy)
    !$$$!! test mode
    !$$$      if(ixc/=1011) then
    !$$$      open(newunit=ifw0w0i,file='W0W0I',form='unformatted')
    !$$$      write(ifw0w0i) nw_i,nw,niw,nq0i
    !$$$      write(ifw0w0i) llw(nw_i:nw,1:nq0i)
    !$$$      write(ifw0w0i) llwI(1:niw,1:nq0i)
    !$$$      write(ifw0w0i) w0(nw_i:nw)
    !$$$      write(ifw0w0i) w0i(1:niw)
    !$$$      close(ifw0w0i)
    !$$$      endif
    do i=nw_i,nw
       write(6,"('w0 =',i4,2f13.4)")i,w0(i)
    enddo
    do i=1,niw
       write(6,"('w0i=',i4,2f13.4)")i,w0i(i)
    enddo
    !! modivy files WVR and WVI
    call ModifyWV0()
    !!
    if(w4pmode) call FinalizeW4p() !W for phonon mode finalized.
  end subroutine w0w0i
  !-------------------------------------------------
  subroutine modifyWV0()
    use m_qbze,only: &
         nqbze,nqibze,qbze,qibze
    use m_rdpp,only: nblochpmx,mrecl
    use m_freq,only: niw ,nw,nw_i
    integer:: ifrcwx,iq,ircw,iw,nini,nend,mreclx
    real(8)::q(3)
    complex(8),allocatable:: zw(:,:)
    character(10):: i2char
    mreclx=mrecl
    !! Read WVR and WVI at Gamma point, and give correct W(0) (averaged in the Gamma cell, where
    !! Gamma cell) is the micro cell of BZ including Gamma point).
    !! === w0,w0i are stored to zw for q=0 ===
    !! === w_ks*wk are stored to zw for iq >nqibz ===
    ! We assume iq=1 is for rank=0
    allocate( zw(nblochpmx,nblochpmx) )
    iq = 1             !iq=1 only 4pi/k**2 /eps part only ! iq = iqxini,iqxend
    q = qibze(:,iq)
    do ircw=1,2
       if (ircw==1) then
          nini=nw_i
          nend=nw
          open(newunit=ifrcwx,  file='WVR.'//i2char(iq), form='unformatted', &
               status='old',access='direct',recl=mreclx)
       elseif(ircw==2) then;  nini=1;      nend=niw;
          open(newunit=ifrcwx,  file='WVI.'//i2char(iq), form='unformatted', &
               status='old',access='direct',recl=mreclx)
       endif
       do iw=nini,nend
          read(ifrcwx, rec= iw-nini+1 ) zw !(1:ngb,1:ngb)
          if( iq==1 ) then
             if(ircw==1) zw(1,1) = w0(iw)
             if(ircw==2) zw(1,1) = w0i(iw)
          endif
          write(ifrcwx,rec=iw-nini+1) zw !(1:ngb,1:ngb)
       enddo
       close(ifrcwx)
    enddo
  end subroutine modifyWV0
  
  subroutine getw0(llw,ii,ie,nq0i,dmlx,epinvq0i,wklm,wbz,lmxax,q0i,epinv, w0,llmat)
    !! == Obtain effective screened Coulomnb interaction W-V(q,omega) at q=0 ==
    !! Roughly speaking, we average out W-V in the gamma cell (=microcell including Gamma point).
    !! Output
    !!     w0(ii:ie), others are inputs.
    !! INPUT
    !! ii: start index of iw
    !! iw: end index of iw
    implicit none
    integer,intent(in)::  ii,ie,nq0i,lmxax
    real(8),intent(in)::  epinvq0i(nq0i,nq0i),dmlx(nq0i,9), wklm((lmxax+1)**2),wbz,epinv(3,3,nq0i),q0i(3,nq0i)
    complex(8),intent(in):: llw(ii:ie,nq0i)
    complex(8),intent(out):: w0(ii:ie),llmat(3,3)
    integer:: lm,lm1,lm2,nlklm,iw,lm1x,lm2x
    complex(8):: llwyl(9,ii:ie),llw_invr(ii:ie,1:nq0i)
    real(8):: epinvq0i_m1(nq0i,nq0i)
    real(8),allocatable::  cg(:,:,:)
    complex(8)::   cgllw((lmxax+1)**2,(lmxax+1)**2),sss,sss2
    real(8):: r2s,emat(3,3),rrr(3),qnorm2(nq0i),ylm,rrr2(3)
    integer:: nlmm,lx,iq0i,iq0x, i,lxx
    real(8),allocatable:: cy(:),yl(:),yl2(:)
    real(8),parameter:: pi  = 4d0*datan(1d0),fpi = 4d0*pi
    write(6,"(' getw0 start:',4i3,d13.6)") nq0i,lmxax,ii,ie,sum(abs(llw(ii:ie,1:nq0i)))
    nlklm=(lmxax+1)**2
    epinvq0i_m1=epinvq0i
    call matinv(nq0i,epinvq0i_m1)
    llwyl=0d0
    do iw= ii,ie
       llw_invr(iw,1:nq0i) = matmul (epinvq0i_m1, llw(iw,1:nq0i))
       ! llw_invr: representation in a linear combination of invariat tensor of epsilon.
       lm=1
       llwyl(lm,iw) = sum(dmlx(1:nq0i,lm)*llw_invr(iw,1:nq0i))
       do lm=5,9
          llwyl(lm,iw)= sum(dmlx(1:nq0i,lm)*llw_invr(iw,1:nq0i))
       enddo
    enddo
    !! omega=0 3x3 matrix, only when ii<= 0 <= nw (only when incluging iw=0 for omega=0)
    if(ii<=0) then
       iw=0
       llmat=0d0
       do iq0i=1,nq0i
          llmat(1:3,1:3) = llmat(1:3,1:3) + llw_invr(iw,iq0i)*epinv(1:3,1:3,iq0i)
       enddo
    endif
    lxx=max(lmxax,2)
    allocate( cg((lxx+1)**2,(lxx+1)**2,(2*lxx+1)**2))
    call rotcg(lxx,(/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/),1,cg) !no rotation

    !$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !$$$      lx=2
    !$$$      allocate(cy((lx+1)**2),yl((lx+1)**2),yl2((lx+1)**2))
    !$$$      call sylmnc(cy,lx)
    !$$$      do iw= ii,ie
    !$$$        do iq0x=1,nq0i
    !$$$        sss=0d0
    !$$$        do iq0i=1,nq0i
    !$$$            sss= sss  + llw_invr(iw,iq0i)*sum(q0i(:,iq0x)*matmul(epinv(:,:,iq0i),q0i(:,iq0x)))
    !$$$     &                 /sum(q0i(:,iq0x)**2)
    !$$$        enddo
    !$$$          write(*,"(' ttt: epinv expansion=',2i4,2f10.5,2x,2d13.5)") iq0x,iw,sss,sss-llw(iw,iq0x)
    !$$$        enddo
    !$$$        write(*,*)
    !$$$      enddo
    !$$$      stop '---- ttt test1: reproduce llw at q0i -----------------------'

    !$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !$$$!! === test for one r vector as for <ehat| epinv|ehat> = \sum_lm dmlx(iq0i,lm) *Y_lm(ehat) ===
    !$$$!! ===== generate YL for a test vector rrr (rrr is ehat above).=====
    !$$$      lx=2
    !$$$      allocate(cy((lx+1)**2),yl((lx+1)**2),yl2((lx+1)**2))
    !$$$      call sylmnc(cy,lx)
    !$$$      do i=1,3
    !$$$         if(i==1) rrr =(/1d0,1d0,0d0/)
    !$$$         if(i==2) rrr =(/1d0,-1d0,0d0/)
    !$$$         if(i==3) rrr =(/0d0,0d0,1d0/)
    !$$$         rrr = rrr/sqrt(sum(rrr**2))
    !$$$c     write(*,"(' testttt: r=',3f10.5)") rrr
    !$$$         call sylm(rrr,yl,lx,r2s) !spherical factor Y( q+G )
    !$$$         do iw= ii,ie
    !$$$            sss=0d0
    !$$$            do iq0i=1,nq0i
    !$$$               sss = sss + llw_invr(iw,iq0i)*sum(rrr*matmul(epinv(:,:,iq0i),rrr))
    !$$$            enddo
    !$$$c            if(abs(llwyl(1,iw)*cy(1)*yl(1)+sum(llwyl(5:9,iw)*cy(5:9)*yl(5:9))-sss)>1d-12) then
    !$$$            write(*,"(' ttt: epinv expansion=',i3,2f10.5,2x,4d13.5)") iw,sss,
    !$$$     &      llwyl(1,iw)*cy(1)*yl(1)+sum(llwyl(5:9,iw)*cy(5:9)*yl(5:9))-sss, llw(iw,i)-sss
    !$$$c            endif
    !$$$         enddo
    !$$$      enddo
    !$$$      stop '---- ttt testxxx: reproduce llw at q0i -----------------------'

    !$$$      lx=lmxax
    !$$$      if(allocated(cy)) deallocate(cy)
    !$$$      if(allocated(yl)) deallocate(yl)
    !$$$      allocate(cy((lx+1)**2),yl((lx+1)**2))
    !$$$      call sylmnc(cy,lx)
    !$$$      rrr=(/-0.36d0,0.20d0,0.4d0/)
    !$$$      rrr=rrr/sqrt(sum(rrr**2))
    !$$$      call sylm(rrr,yl,lx,r2s) !spherical factor Y( q+G )
    !$$$      write(*,"(  ' test input: q=',3f10.5)") rrr

    w0=0d0
    do iw =ii,ie
       lm2x=0
       cgllw=0d0
       do lm2=1,nlklm
          if(mod(ll(lm2),2)==1) cycle  !only even l
          lm2x=lm2x+1
          lm1x=0
          do lm1=1,nlklm
             if(mod(ll(lm1),2)==1) cycle !only even l
             lm1x=lm1x+1
             cgllw(lm1x,lm2x) = cg(lm2,1,lm1)*llwyl(1,iw) + sum(cg(lm2,5:9,lm1)*llwyl(5:9,iw))
          enddo
       enddo
       !! === inversion of 1 = \sum_L1 llwyl(L1)Y_L1 * \sum_L2 K_L2 Y_L2===
       !! Both sides are multipled by Y_L and spherically integraled.
       !! warn: Klm should be not confused with wklm
       !! Klm is defeined in Christoph's paper PRB 125102-8
       call matcinv(lm2x,cgllw(1:lm2x,1:lm2x))
       cgllw(1:lm2x,1:lm2x)= sqrt(fpi)*cgllw(1:lm2x,1:lm2x)
       ! qrt(fpi) comes from \int d\Omega Y_0 (right hand side of inversion eq).

       !! === spherical integral of Klm ===
       lm1x=0
       do lm1=1, nlklm
          if(mod(ll(lm1),2)==1) cycle
          lm1x=lm1x+1
          w0(iw)= w0(iw)+ wklm(lm1)* fpi* cgllw(lm1x,1)/wbz      ! Klm=cgllw(:,1)
          if(lm1==1) w0(iw)= w0(iw) - wklm(lm1)*fpi**1.5d0/wbz   ! fpi**1.5 --> subtract Coulomb
          ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !          w0(iw)= wklm(1)*fpi**1.5*(1d0/llw(iw,1)-1d0)/wbz !test case of llw at q0i(:,iq0i=1).
       enddo
       !$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       !$$$        lx=lmxax
       !$$$        write(6,*)' lx=',lx
       !$$$        if(allocated(cy)) deallocate(cy)
       !$$$        if(allocated(yl)) deallocate(yl)
       !$$$        allocate(cy((lx+1)**2),yl((lx+1)**2))
       !$$$        call sylmnc(cy,lx)
       !$$$        rrr=(/1d0,1d0,0d0/)
       !$$$c        rrr=(/1d0,-1d0,0d0/)
       !$$$c        rrr=(/0d0,0d0,1d0/)
       !$$$        rrr=rrr/sqrt(sum(rrr**2))
       !$$$        call sylm(rrr,yl,lx,r2s) !spherical factor Y( q+G )
       !$$$c        write(*,"(  ' qqq test input: q=',3f10.5)") rrr
       !$$$        sss=0d0
       !$$$        lm2x=0
       !$$$        do lm2=1,nlklm
       !$$$           if(mod(ll(lm2),2)==1) cycle !only even l
       !$$$           lm2x=lm2x+1
       !$$$           sss= sss + cgllw(lm2x,1) *cy(lm2)*yl(lm2)
       !$$$        enddo
       !$$$        write(*,"(' qqq ep test:',i3,2f10.5,2x,2f10.5,2x,2f10.5)") iw,sss,
       !$$$     &       1d0/llw(iw,1)
       !$$$c     &      1d0/(llwyl(1,iw)*cy(1)*yl(1) + sum(llwyl(5:9,iw)*cy(5:9)*yl(5:9)))
       !$$$      stop ' --- test3 ttt vvvvvvvvvvvvvvvvv test end vvvvvvvvvvvvvvvvv ---'
    enddo
    !$$$      do iw=ii,ie
    !$$$        write(6,"('w0 & w0(1:nq0i)=',2d11.3,10(1x,2d11.3))")
    !$$$     &  w0(iw), (wklm(1)*fpi**1.5*(1d0/llw(iw,iq0i)-1d0)/wbz,iq0i=1,nq0i)
    !$$$  enddo
    deallocate(cg)
  end subroutine getw0
end module m_w0w0i
