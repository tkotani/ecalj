      integer function isum(n,idx,incx)
c
c     takes the sum of the values.  Adapted from:
c     jack dongarra, linpack, 3/11/78.
c
      integer idx(1)
      integer i,incx,n,nincx
c
      isum = 0
      if (n .le. 0) return
      nincx = n*incx
      do  10  i = 1, nincx,incx
   10 isum = isum + idx(i)
      end


