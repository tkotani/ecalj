subroutine hanr(rsq,lmin,lmax,nrx,nr,e,xi)
  !- Vector of unsmoothed hankel functions for l=0...lmax, negative e.
  ! ---------------------------------------------------------------
  !i Inputs
  !i   rsq,nr  vector of points r**2, and number of points.
  !i   e       energy of Hankel.
  !i   lmin:lmax generate xi from lmin:lmax.  lmin must be -1 or 0.
  !i   nrx     leading dimension of xi.
  !o Outputs
  !o   xi      radial Hankel functions/r**l for: xi(1..nr, lmin:lmax)
  !o           Solid hankel function is hl(ilm) = xi(l)*cy(ilm)*yl(ilm)
  !o           Energy derivative is    hlp(ilm) = x(l-1)/2*cy(ilm)*yl(ilm)
  ! ---------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: nrx,nr,lmin,lmax
  double precision :: rsq(nr),e,xi(nrx,lmin:lmax)
  ! ... Local parameters
  double precision :: sre,r2,r,h0,xx,akap
  integer :: l,ir

  if (lmin /= 0 .AND. lmin /= -1) &
       call rx('hanr: input lmin must be -1 or 0')
  if (lmax < lmin .OR. nr <= 0) return
  akap = dsqrt(-e)

  ! --- Make xi(lmin), xi(lmin+1) ---
  ! ... xi(-1) for lmax=-1 only
  if (lmin == -1 .AND. lmax == -1) then
     do  1  ir = 1, nr
        r = dsqrt(rsq(ir))
        sre = akap*r
        h0 = dexp(-sre)/r
        xi(ir,-1) = -h0*sre/e
1    enddo
  elseif (lmin == -1) then
     do  3  ir = 1, nr
        r = dsqrt(rsq(ir))
        sre = akap*r
        h0 = dexp(-sre)/r
        xi(ir,0) = h0
        xi(ir,-1) = -h0*sre/e
3    enddo
     ! ... xi(0) for lmax=0 only
  elseif (lmax == 0) then
     do  5  ir = 1, nr
        r2 = rsq(ir)
        r = dsqrt(r2)
        xi(ir,0) = dexp(-akap*r)/r
5    enddo
  else
     do  10  ir = 1, nr
        r = dsqrt(rsq(ir))
        sre = akap*r
        h0 = dexp(-sre)/r
        xi(ir,0) = h0
        xi(ir,1) = h0*(1d0+sre)/rsq(ir)
10   enddo
  endif

  ! --- xi(*,lmin+2:lmax) by upward recursion ---
  do  30  l = lmin+2, lmax
     xx = 2*l-1
     do  32  ir = 1, nr
        xi(ir,l) = (xx*xi(ir,l-1) - e*xi(ir,l-2))/rsq(ir)
32   enddo
30 enddo

end subroutine hanr

