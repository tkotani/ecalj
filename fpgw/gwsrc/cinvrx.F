      subroutine cinvrx(a,b,n,mx,m)
C- Given a complex matrix a, this returns its inverse matrix.
C----------------------------------------------------------------------
Cr     Given a complex matrix a, this returns its inverse matrix.
Cr     After calling of this program a is replaced by the lu
Cr     decomposition. No pivoting is employed.
Cr    In this version the operations running over the first index of
Cr    the matrix come up to the outermost loop, suitable for a big
Cr    matrix with rather small size regardint the first index.
Cr    coded by H.Akai., Feb. 4, 1992, Osaka
C----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 a(mx,n,n),b(mx,n,n)
      do 40 k=1,m
        do 10 i=1,n-1
        do 10 j=i+1,n
          a(k,j,i)=-a(k,j,i)/a(k,i,i)
          do 30 l=1,i-1
            a(k,j,l)=a(k,j,l)+a(k,j,i)*a(k,i,l)
   30     continue
        do 10 l=i+1,n
          a(k,j,l)=a(k,j,l)+a(k,j,i)*a(k,i,l)
   10   continue
      do 40 i=n,1,-1
        do 60 j=1,i-1
          b(k,i,j)=a(k,i,j)
        do 60 l=i+1,n
          b(k,i,j)=b(k,i,j)-a(k,i,l)*b(k,l,j)
   60   continue
        do 70 j=i+1,n
          b(k,i,j)=(0d0,0d0)
        do 70 l=i+1,n
          b(k,i,j)=b(k,i,j)-a(k,i,l)*b(k,l,j)
   70   continue
        b(k,i,i)=(1d0,0d0)
        do 100 l=i+1,n
          b(k,i,i)=b(k,i,i)-a(k,i,l)*b(k,l,i)
  100   continue
      do 40 j=1,n
        b(k,i,j)=b(k,i,j)/a(k,i,i)
   40 continue
      return
      end

c------------------------------------------------------------------
c$$$c------------------------------------------------------------------
c$$$      subroutine TESTrQGpsi
c$$$c- readin DATA4GW and so on and write some files for GW.
c$$$      implicit none
c$$$      integer(4) :: ifiqg,ngrp,ifi,iq,ig,nqnum,ngpmx
c$$$      real(8) :: QpGcut_psi
c$$$      real(8),allocatable :: q(:,:)
c$$$      integer(4),allocatable :: ngp(:),ngvecp(:,:,:)
c$$$
c$$$
c$$$c- Readin QGpsi
c$$$c True q is given by
c$$$c    True_q(1:3)     = 2*pi/alat * q(1:3)
c$$$c True G is given by
c$$$c    True_G(1:3,igp) = 2*pi/alat * matmul(qlat * ngvec(1:3,igp)) ,igp=1,ngp
c$$$c
c$$$c------------------------------------------------------------------
c$$$      ifiqg =101
c$$$      open(ifiqg, file='QGpsi',form='unformatted')
c$$$      read(ifiqg) nqnum, ngpmx, QpGcut_psi
c$$$      write(6,*)nqnum, ngpmx, QpGcut_psi
c$$$      allocate( q(3,nqnum),ngp(nqnum),ngvecp(3,ngpmx,nqnum) )
c$$$      do  iq    = 1, nqnum
c$$$        read (ifiqg) q(1:3,iq), ngp(iq)
c$$$        read (ifiqg) ngvecp(1:3,1:ngp(iq),iq)
c$$$        write(6,"('q=',3f10.4,' ngp=',i5)" ) q(1:3,iq), ngp(iq)
c$$$        do ig=1,ngp(iq)
c$$$          write(6,"(' ngvec=',3i5)")  ngvecp(1:3,ig,iq)
c$$$        enddo
c$$$      enddo
c$$$      close(ifiqg)
c$$$      end
