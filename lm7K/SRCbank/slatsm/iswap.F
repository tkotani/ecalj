      subroutine iswap (n,dx,incx,dy,incy)
c
c     interchanges two vectors. Adapted from:
c     jack dongarra, linpack, 3/11/78.
c
      integer dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,n
c
      ix = 1
      iy = 1
      if (incx.lt.0) ix = (1-n)*incx + 1
      if (incy.lt.0) iy = (1-n)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      end

