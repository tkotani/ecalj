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
  integer:: lm,lm1,lm2,nlklm,iw,ll,lm1x,lm2x
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
