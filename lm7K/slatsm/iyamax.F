      integer function iyamax(n,zxr,zxi,incx)
C     implicit none
c
c     finds the index of element having max. absolute value.
c     Adapted from blas, but avoids separates real, imaginary parts
c     jack dongarra, 1/15/85.
c     modified 3/93 to return if incx .le. 0.
c
      double precision zxr(1),zxi(1)
      double precision smax
      integer i,incx,ix,n
c
      iyamax = 0
      if( n.lt.1 .or. incx.le.0 )return
      iyamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = dabs(zxr(1)) + dabs(zxi(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(zxr(ix))+dabs(zxi(ix)).le.smax) go to 5
         iyamax = i
         smax = dabs(zxr(ix))+dabs(zxi(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smax = dabs(zxr(1)) + dabs(zxi(1))
      do 30 i = 2,n
         if(dabs(zxr(i))+dabs(zxi(i)).le.smax) go to 30
         iyamax = i
         smax = dabs(zxr(i))+dabs(zxi(i))
   30 continue
      return
      end

