      integer(4) function iqagree(qb,qqr,nxx,nx,nq,plat,qlat)
      integer(4) :: nxx,nx(nq),ix,iq
      real(8) :: qb(3),qbo(3),qqr(3,nxx,nq),qqro(3)
     & ,plat(3,3),qlat(3,3),aa1(3)
      iqagree = 0
      do iq = 1,nq
        do ix = 1,nx(iq)
          do i=1,3
            aa1(i) = sum(plat(:,i)*(qqr(:,ix,iq)-qb))
          enddo
          if(sum(abs( aa1-nint(aa1) ))<1d-6) then
            iqagree = iq
            return
          endif
ccccccccccc
c       if(qb(1)==0.5d0.and.qb(2)==0.5d0.and.qb(3)==0.5d0) then
c       print *
c      write(6,"(' iq ix =',2i3)")iq,ix
c       write(6,"('  qb =',3f10.4)")qb
c       write(6,"('  qqr=',3f10.4)")qqr(:,ix,iq)
c       write(6,"('  qbo=',3f10.4)")qbo
c       write(6,"(' qqro=',3f10.4)")qqro
c       endif
ccccccccccc
        enddo
      enddo
      end

c-------------------------------------
      real(8) function polinta2(xx, x,y,n, xmin, xmax, d0min, d0max)
c Interpolation considering the zero derivative at xmin (and/or) at xmax.
c  d0min=T means that the derivative of the interpolated function at xmin is zero.
c  d0max=T means that the derivative of the interpolated function at xmax is zero.
      implicit none
      integer(4):: n ,imin,imax, id1,id2, nnn,ni,nx,i
      real(8) :: xx, x(n), y(n),xin(3*n),yin(3*n),xmin,xmax,
     &           polinta,eps=1d-13
      logical :: d0min, d0max
      if(d0min) then
        if(abs(x(1)-xmin)<eps) then
          ni=n-1
        else
          ni=n
        endif
        do i=1,ni
          xin(i) = xmin -(x(n+1-i)-xmin)
          yin(i) = y(n+1-i)
        enddo
      else
        ni = 0
      endif

      xin(ni+1:ni+n)= x(1:n)
      yin(ni+1:ni+n)= y(1:n)

      if(d0max) then
        if(abs(x(n)-xmax)<eps) then
          nx=n-1
        else
          nx=n
        endif
        do i=1,nx
          xin(i+ni+n) = xmax + (xmax - x(1+nx-i))
          yin(i+ni+n) = y(1+nx-i)
        enddo
      else
        nx = 0
      endif
      nnn = ni+ n + nx
      polinta2 = polinta(xx,xin,yin,nnn)
      end











