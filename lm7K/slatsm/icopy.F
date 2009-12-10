
      subroutine icopy(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.  Adapted from:
c     jack dongarra, linpack, 3/11/78.
c
      integer dx(1),dy(1)
      integer i,incx,incy,ix,iy,n
      ix = 1
      iy = 1
      if (incx .lt. 0) ix = (1-n)*incx + 1
      if (incy .lt. 0) iy = (1-n)*incy + 1
      do  10  i = 1, n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      end

