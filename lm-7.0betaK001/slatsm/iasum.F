      integer function iasum(n,ix,incx)
c
c     takes the sum of the absolute values.  Adapted from:
c     jack dongarra, linpack, 3/11/78.
c
      integer ix(1)
      integer i,incx,n,nincx
c
      iasum = 0
      if (n .le. 0) return

      nincx = n*incx
      do  10  i = 1, nincx,incx
        iasum = iasum + iabs(ix(i))
   10 continue
      end


