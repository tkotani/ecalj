      subroutine daxpby(n,da,dx,incx,db,dy,incy)
C- Constant times a vector plus a constant times vector.
C ----------------------------------------------------------------
Ci Inputs
Ci   n,da,dx,incx,db,dy,incy
Co Outputs
Co   y replaced by da*dx + db*dy
Cr Remarks
Cr   Extension of daxpy, in which y is replaced by da*dx + dy
C ----------------------------------------------------------------
C      implicit none
C Passed parameters
      double precision dx(1),dy(1),da,db
      integer incx,incy,n
C Local parameters
      integer i,ix,iy
c
      ix = 1
      iy = 1
      if (incx .lt. 0) ix = (1-n)*incx + 1
      if (incy .lt. 0) iy = (1-n)*incy + 1
      do  10  i = 1, n
        dy(iy) = db*dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      end

