c--------------------------------------------------------------------
c routines to generate indices and dimensions
c--------------------------------------------------------------------
      integer function idxk (q,qibz,nqibz)
      implicit none
      integer(4) :: nqibz,ic
      real(8) :: q(3),qibz(3,nqibz)
      do ic = 1,nqibz
        if(sum(abs(qibz(:,ic)-q))< 1d-5) goto 1111
      end do
      stop ' idxk: no match '
 1111 continue
      idxk = ic
      return
      end

c--------------------------------------------------------------------
c routines to generate indices and dimensions
c cyclic search 8/18/2005
c--------------------------------------------------------------------
      integer function idxk_cyclic (q,qibz,nqibz)
      implicit none
      integer(4) :: nqibz,ic,icc
      integer(4),save:: is = 1
      real(8) :: q(3),qibz(3,nqibz)
      do ic = is, nqibz+is - 1
        icc=ic
        if(icc>nqibz ) icc= ic-nqibz
        if(sum(abs(qibz(:,icc)-q))< 1d-6) goto 1111
      end do
      stop ' idxk: no match '
 1111 continue
      idxk_cyclic = icc
      is   = icc
      return
      end
