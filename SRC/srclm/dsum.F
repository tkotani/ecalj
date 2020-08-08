      double precision function dsum(n,dx,incx)
c
c     takes the sum of the values.  Adapted from:
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1)
      integer i,incx,n,nincx
c
      dsum = 0d0
      if (n .le. 0) return

      nincx = n*incx
      do  10  i = 1, nincx,incx
        dsum = dsum + dx(i)
   10 continue
      end


