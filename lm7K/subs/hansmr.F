      subroutine hansmr(r,e,a,xi,lmax)
C- Smoothed hankel functions for l=0...lmax, negative e.
C ---------------------------------------------------------------
Co  Outputs: xi(0:lmax)
Co
Cr  xi is the radial part divided by r**l.
Cr  A solid smoothed hankel is from xi as:
Cr  hl(ilm) = xi(l)*cy(ilm)*yl(ilm)

Cr  See J. Math. Phys. 39, 3393 (1998).
Cr  xi(l)= 2/sqrt(pi) * 2^l int_a^inf u^2l dexp(-r^2*u^2+kap2/4u^2) du
Cu Updates
Cu   04 Jul 08 (S. Lozovoi) Bug fix, lmax
Cu   11 May 07 (S. Lozovoi) small bug fixes
C ---------------------------------------------------------------
C     implicit none
      integer lmax,l,n,nmax
      double precision r,e,a,xi(0:lmax)
      double precision a0(0:40),a2,add,akap,al,cc,derfc,dudr,ema2r2,fac,
     .  gl,r2n,radd,rfac,rhs,rlim,srpi,sum,ta,ta2l,tol,u,uminus,uplus,w
      parameter (nmax=1000, tol=1d-20, srpi=1.77245385090551602729817d0)

      if (e .gt. 0d0) call rx('hansmr: e gt 0')
      if (lmax .gt. 40) call rx('hansmr: lmax gt 40')
C ... old tolerance
c     rlim = 5d0/(a+1d0)
      rlim = 1.5d0/a
      akap = dsqrt(dabs(e))
      ta = a+a
      a2 = a*a
      cc = 4d0*a2*a*dexp(e/(ta*ta))/srpi

c ---- return zero if exponential part very small -------
c      if (akap*r.gt.80d0) then
c        do  27  l = 0, lmax
c  27    xi(l)=0d0
c        return
c        endif

C  --- call bessl if exp(-a*a*r*r) is very small ---
      if (a*r .gt. 10d0) then
        call bessl(e*r*r,lmax,a0,xi)
        rfac = r
        do  28  l = 0, lmax
          rfac = rfac*(1d0/(r*r))
          xi(l) = rfac*xi(l)
   28   continue
        return
      endif

C --- Power series for small r ---
      if (r .gt. rlim) goto 90
      a0(0) = cc/(ta*a) - akap*derfc(akap/ta)
      rhs = cc
      fac = 1d0
      al = a0(0)
      do  10  l = 1, lmax
        al = -(e*al+rhs)/(2*l*(2*l+1))
        rhs = -rhs*a2/l
        fac = -2d0*fac*l
        a0(l) = fac*al
   10 continue
      ta2l = 1d0
      do  20  l = 0, lmax
        rhs = cc*ta2l
        sum = a0(l)
        add = sum
        r2n = 1d0
        do  21  n = 1, nmax
          add = -(e*add+rhs)/( 2*n*(2*n+(l+l+1)) )
          r2n = r2n*(r*r)
          radd = add*r2n
          sum = sum+radd
          if (dabs(radd) .lt. tol) goto 22
          rhs = -rhs*a2/n
   21   continue
        print *, 'hansmr (warning): power series did not converge'
   22   continue
        xi(l) = sum
        ta2l = ta2l*(2d0*a2)
   20 continue
      return

C --- Big r: make xi0,xi1 explicitly; the higher l by recursion ---
   90 ema2r2 = dexp(-a*a*r*r)
      uminus = derfc(akap/ta-r*a)*dexp(-akap*r)
      uplus = derfc(akap/ta+r*a)*dexp(+akap*r)
      u = .5d0*(uminus-uplus)
      w = -.5d0*(uminus+uplus)
      dudr = akap*w + ta*dexp(e/(ta*ta))*ema2r2/srpi
      xi(0) = u/r
      if (lmax .ge. 1) then
        xi(1) = (u/r - dudr)/(r*r)
        gl = cc*ema2r2
        if (lmax .ge. 2) then
          do  l = 2, lmax
            xi(l) = ((2*l-1)*xi(l-1) -e*xi(l-2) - gl)/(r*r)
            gl = 2d0*a2*gl
          enddo
        endif
      endif

      end

