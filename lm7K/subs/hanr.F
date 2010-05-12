      subroutine hanr(rsq,lmin,lmax,nrx,nr,e,xi)
C- Vector of unsmoothed hankel functions for l=0...lmax, negative e.
C ---------------------------------------------------------------
Ci Inputs
Ci   rsq,nr  vector of points r**2, and number of points.
Ci   e       energy of Hankel.
Ci   lmin:lmax generate xi from lmin:lmax.  lmin must be -1 or 0.
Ci   nrx     leading dimension of xi.
Co Outputs
Co   xi      radial Hankel functions/r**l for: xi(1..nr, lmin:lmax)
Co           Solid hankel function is hl(ilm) = xi(l)*cy(ilm)*yl(ilm)
Co           Energy derivative is    hlp(ilm) = x(l-1)/2*cy(ilm)*yl(ilm)
C ---------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nrx,nr,lmin,lmax
      double precision rsq(nr),e,xi(nrx,lmin:lmax)
C ... Local parameters
      double precision sre,r2,r,h0,xx,akap
      integer l,ir

      if (lmin .ne. 0 .and. lmin .ne. -1)
     .  call rx('hanr: input lmin must be -1 or 0')
      if (lmax .lt. lmin .or. nr .le. 0) return
      akap = dsqrt(-e)

C --- Make xi(lmin), xi(lmin+1) ---
C ... xi(-1) for lmax=-1 only
      if (lmin .eq. -1 .and. lmax .eq. -1) then
        do  1  ir = 1, nr
          r = dsqrt(rsq(ir))
          sre = akap*r
          h0 = dexp(-sre)/r
          xi(ir,-1) = -h0*sre/e
    1   continue
      elseif (lmin .eq. -1) then
        do  3  ir = 1, nr
          r = dsqrt(rsq(ir))
          sre = akap*r
          h0 = dexp(-sre)/r
          xi(ir,0) = h0
          xi(ir,-1) = -h0*sre/e
    3   continue
C ... xi(0) for lmax=0 only
      elseif (lmax .eq. 0) then
        do  5  ir = 1, nr
          r2 = rsq(ir)
          r = dsqrt(r2)
          xi(ir,0) = dexp(-akap*r)/r
    5   continue
      else
        do  10  ir = 1, nr
          r = dsqrt(rsq(ir))
          sre = akap*r
          h0 = dexp(-sre)/r
          xi(ir,0) = h0
          xi(ir,1) = h0*(1d0+sre)/rsq(ir)
   10   continue
      endif

C --- xi(*,lmin+2:lmax) by upward recursion ---
      do  30  l = lmin+2, lmax
        xx = 2*l-1
        do  32  ir = 1, nr
          xi(ir,l) = (xx*xi(ir,l-1) - e*xi(ir,l-2))/rsq(ir)
   32   continue
   30 continue

      end

