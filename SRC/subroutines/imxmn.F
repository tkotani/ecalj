c$$$      subroutine imxmn(n,idx,incx,imin,imax)
c$$$C- Finds maximum and minimum value of an integer array
c$$$Cu   23 Apr 02 Bug fix
c$$$C     implicit none
c$$$      integer n,idx(n),incx,imin,imax
c$$$      integer i,ix
c$$$
c$$$      imin = idx(1)
c$$$      imax = idx(1)
c$$$      ix = 1
c$$$      do  10  i = 2, n
c$$$        ix = ix + incx
c$$$        imin = min(idx(ix),imin)
c$$$        imax = max(idx(ix),imax)
c$$$   10 continue
c$$$      end
c$$$
