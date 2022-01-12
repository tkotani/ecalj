      subroutine getw0(llw,ii,ie,nq0i,dmlx,epinvq0i,wklm,wbz,lmxax,q0i,epinv, w0,llmat)
!! == Obtain effective screened Coulomnb interaction W-V at q=0 ==
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
!!
      integer:: lm,lm1,lm2,nlklm,iw,ll,lm1x,lm2x
      complex(8):: llwyl(9,ii:ie),llw_invr(ii:ie,1:nq0i)
      real(8):: epinvq0i_m1(nq0i,nq0i)
      real(8)::  cg((lmxax+1)**2,(lmxax+1)**2,(2*lmxax+1)**2)
      complex(8)::   cgllw((lmxax+1)**2,(lmxax+1)**2),sss,sss2
      real(8):: fpi,pi,r2s,emat(3,3),rrr(3),qnorm2(nq0i),ylm,rrr2(3)
      integer:: nlmm,lx,iq0i,iq0x, i
      real(8),allocatable:: cy(:),yl(:),yl2(:)
!!
      pi  = 4d0*datan(1d0)
      fpi = 4d0*pi
      nlklm=(lmxax+1)**2
      epinvq0i_m1=epinvq0i
c      write(6,"(' getw0 start1:',4i3,d13.6)") nq0i,lmxax,ii,ie
      write(6,"(' getw0 start:',4i3,d13.6)") nq0i,lmxax,ii,ie,sum(abs(llw(ii:ie,1:nq0i)))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      nlm=(lmxax+1)**2
c      do lm=1,nlm
c        if(abs(wklm(lm))>1d-10) write(6,*)lm,ll(lm),wklm(lm)
c      enddo
c      stop 'test end vvvvvvvvvvvvvvvvvvvvvv'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
c   write(*,"('qqq iw llw',i3,3(2d12.4,x),3x,20(2d12.4,2x))")iw,llw(iw,1:nq0i),llwyl(1,iw),llwyl(5:9,iw)
      enddo

!! omega=0 3x3 matrix, only when ii<= 0 <= nw (only when incluging iw=0 for omega=0)
      if(ii<=0) then
         iw=0
         llmat=0d0
         do iq0i=1,nq0i
            llmat(1:3,1:3) = llmat(1:3,1:3) + llw_invr(iw,iq0i)*epinv(1:3,1:3,iq0i)
         enddo
      endif

      call rotcg(lmxax,(/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/),1,cg) !no rotation

c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$      lx=2
c$$$      allocate(cy((lx+1)**2),yl((lx+1)**2),yl2((lx+1)**2))
c$$$      call sylmnc(cy,lx)
c$$$      do iw= ii,ie
c$$$        do iq0x=1,nq0i
c$$$        sss=0d0
c$$$        do iq0i=1,nq0i
c$$$            sss= sss  + llw_invr(iw,iq0i)*sum(q0i(:,iq0x)*matmul(epinv(:,:,iq0i),q0i(:,iq0x)))
c$$$     &                 /sum(q0i(:,iq0x)**2)
c$$$        enddo
c$$$          write(*,"(' ttt: epinv expansion=',2i4,2f10.5,2x,2d13.5)") iq0x,iw,sss,sss-llw(iw,iq0x)
c$$$        enddo
c$$$        write(*,*)
c$$$      enddo
c$$$      stop '---- ttt test1: reproduce llw at q0i -----------------------'

c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$!! === test for one r vector as for <ehat| epinv|ehat> = \sum_lm dmlx(iq0i,lm) *Y_lm(ehat) ===
c$$$!! ===== generate YL for a test vector rrr (rrr is ehat above).=====
c$$$      lx=2
c$$$      allocate(cy((lx+1)**2),yl((lx+1)**2),yl2((lx+1)**2))
c$$$      call sylmnc(cy,lx)
c$$$      do i=1,3
c$$$         if(i==1) rrr =(/1d0,1d0,0d0/)
c$$$         if(i==2) rrr =(/1d0,-1d0,0d0/)
c$$$         if(i==3) rrr =(/0d0,0d0,1d0/)
c$$$         rrr = rrr/sqrt(sum(rrr**2))
c$$$c     write(*,"(' testttt: r=',3f10.5)") rrr
c$$$         call sylm(rrr,yl,lx,r2s) !spherical factor Y( q+G )
c$$$         do iw= ii,ie
c$$$            sss=0d0
c$$$            do iq0i=1,nq0i
c$$$               sss = sss + llw_invr(iw,iq0i)*sum(rrr*matmul(epinv(:,:,iq0i),rrr))
c$$$            enddo
c$$$c            if(abs(llwyl(1,iw)*cy(1)*yl(1)+sum(llwyl(5:9,iw)*cy(5:9)*yl(5:9))-sss)>1d-12) then
c$$$            write(*,"(' ttt: epinv expansion=',i3,2f10.5,2x,4d13.5)") iw,sss,
c$$$     &      llwyl(1,iw)*cy(1)*yl(1)+sum(llwyl(5:9,iw)*cy(5:9)*yl(5:9))-sss, llw(iw,i)-sss
c$$$c            endif
c$$$         enddo
c$$$      enddo
c$$$      stop '---- ttt testxxx: reproduce llw at q0i -----------------------'

c$$$      lx=lmxax
c$$$      if(allocated(cy)) deallocate(cy)
c$$$      if(allocated(yl)) deallocate(yl)
c$$$      allocate(cy((lx+1)**2),yl((lx+1)**2))
c$$$      call sylmnc(cy,lx)
c$$$      rrr=(/-0.36d0,0.20d0,0.4d0/)
c$$$      rrr=rrr/sqrt(sum(rrr**2))
c$$$      call sylm(rrr,yl,lx,r2s) !spherical factor Y( q+G )
c$$$      write(*,"(  ' test input: q=',3f10.5)") rrr

!! ---
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
        !sqrt(fpi) comes from \int d\Omega Y_0 (right hand side of inversion eq).

!! === spherical integral of Klm ===
        lm1x=0
        do lm1=1, nlklm
          if(mod(ll(lm1),2)==1) cycle
          lm1x=lm1x+1
          w0(iw)= w0(iw)+ wklm(lm1)* fpi* cgllw(lm1x,1)/wbz      ! Klm=cgllw(:,1)
          if(lm1==1) w0(iw)= w0(iw) - wklm(lm1)*fpi**1.5d0/wbz   ! fpi**1.5 --> subtract Coulomb
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c          w0(iw)= wklm(1)*fpi**1.5*(1d0/llw(iw,1)-1d0)/wbz !test case of llw at q0i(:,iq0i=1).
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        enddo

c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$        lx=lmxax
c$$$        write(6,*)' lx=',lx
c$$$        if(allocated(cy)) deallocate(cy)
c$$$        if(allocated(yl)) deallocate(yl)
c$$$        allocate(cy((lx+1)**2),yl((lx+1)**2))
c$$$        call sylmnc(cy,lx)
c$$$        rrr=(/1d0,1d0,0d0/)
c$$$c        rrr=(/1d0,-1d0,0d0/)
c$$$c        rrr=(/0d0,0d0,1d0/)
c$$$        rrr=rrr/sqrt(sum(rrr**2))
c$$$        call sylm(rrr,yl,lx,r2s) !spherical factor Y( q+G )
c$$$c        write(*,"(  ' qqq test input: q=',3f10.5)") rrr
c$$$        sss=0d0
c$$$        lm2x=0
c$$$        do lm2=1,nlklm
c$$$           if(mod(ll(lm2),2)==1) cycle !only even l
c$$$           lm2x=lm2x+1
c$$$           sss= sss + cgllw(lm2x,1) *cy(lm2)*yl(lm2)
c$$$        enddo
c$$$        write(*,"(' qqq ep test:',i3,2f10.5,2x,2f10.5,2x,2f10.5)") iw,sss,
c$$$     &       1d0/llw(iw,1)
c$$$c     &      1d0/(llwyl(1,iw)*cy(1)*yl(1) + sum(llwyl(5:9,iw)*cy(5:9)*yl(5:9)))
c$$$      stop ' --- test3 ttt vvvvvvvvvvvvvvvvv test end vvvvvvvvvvvvvvvvv ---'
      enddo
c$$$      do iw=ii,ie
c$$$        write(6,"('w0 & w0(1:nq0i)=',2d11.3,10(1x,2d11.3))")
c$$$     &  w0(iw), (wklm(1)*fpi**1.5*(1d0/llw(iw,iq0i)-1d0)/wbz,iq0i=1,nq0i)
c$$$      enddo
      end
