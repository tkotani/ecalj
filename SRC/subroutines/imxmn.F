      subroutine imxmn(n,idx,incx,imin,imax)
C- Finds maximum and minimum value of an integer array
Cu   23 Apr 02 Bug fix
C     implicit none
      integer n,idx(n),incx,imin,imax
      integer i,ix

      imin = idx(1)
      imax = idx(1)
      ix = 1
      do  10  i = 2, n
        ix = ix + incx
        imin = min(idx(ix),imin)
        imax = max(idx(ix),imax)
   10 continue
      end

