      subroutine yyxcpy(n,ar,ai,xr,xi,incx,yr,yi,incy,sw)
C- Complex daxpy, using real arithmetic and complex conjugate of x
C     implicit none
C Passed parameters
      logical sw
      integer n,incx,incy
      double precision ar,ai,xr(1),xi(1),yr(1),yi(1)
C Local variables
      integer ix,iy,i

      if (n .le. 0) return

      call daxpy(n, ar,xr,incx,yr,incy)
      if (sw) then
        call daxpy(n, ai,xi,incx,yr,incy)
        call daxpy(n,-ar,xi,incx,yi,incy)
        call daxpy(n, ai,xr,incx,yi,incy)
      endif
      end

