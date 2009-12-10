      double complex function zsum(n,dx,incx)
c
c     takes the sum of the values.  Adapted from:
c     jack dongarra, linpack, 3/11/78.
c
      double complex dx(1)
      integer i,incx,n,nincx
c
      zsum = 0d0
      if (n .le. 0) return

      nincx = n*incx
      do  10  i = 1, nincx,incx
        zsum = zsum + dx(i)
   10 continue
      end


