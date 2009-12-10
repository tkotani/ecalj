      subroutine  yscal(n,zar,zai,zxr,zxi,incx)
c
c     scales a complex vector by a constant, using real arithmetic
c     jack dongarra, 3/11/78.
c
      double precision zar,zai,zxr(1),zxi(1),tmp
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
C     if(incx.lt.0)ix = (-n+1)*incx + 1
      do 10 i = 1,n
        tmp      = zar*zxr(ix) - zai*zxi(ix)
        zxi(ix) = zar*zxi(ix) + zai*zxr(ix)
        zxr(ix) = tmp
        ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        tmp    = zar*zxr(i) - zai*zxi(i)
        zxi(i) = zar*zxi(i) + zai*zxr(i)
        zxr(i) = tmp
   30 continue
      return
      end

