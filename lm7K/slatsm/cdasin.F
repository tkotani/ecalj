      double complex function cdasin (zinp)
c Adapted from ...  w. fullerton, c3, los alamos scientific lab.
c ref -- l. l. pennisi, elements of complex variables,  holt, rinehart
c and winston, 1963.  page 126.
c
C     implicit none
      double complex zinp, z, z2, sqzp1, ci, cdlog, cdsqrt
      integer nterms,i
      double precision twoi,pi2,pi,rmin,d1mach,dlog,r
c     external cdlog, cdsqrt, d1mach
      data pi2 /1.5707963267 9489661923d0/
      data pi /3.1415926535 8979324d0/
      data ci /(0d0,1d0)/
      data nterms, rmin / 0, 0d0 /
c
      if (nterms .ne. 0) go to 10
c nterms = alog(eps)/alog(rmax)  where rmax = 0.1
      nterms = -0.4343d0*dlog(d1mach(3)) + 1d0
      rmin = dsqrt (6d0*d1mach(3))
c
 10   z = zinp
      r = cdabs (z)
      if (r .gt. 0.1d0) go to 30
c
      cdasin = z
      if (r .lt. rmin) return
c
      cdasin = (0d0, 0d0)
      z2 = z*z
      do 20 i=1,nterms
        twoi = 2*(nterms-i) + 1
        cdasin = 1d0/twoi + twoi*cdasin*z2/(twoi+1d0)
 20   continue
      cdasin = z*cdasin
      return
c
 30   if (dble(zinp) .lt. 0d0) z = -zinp
c
      sqzp1 = cdsqrt (z+1d0)
      if (dimag(sqzp1) .lt. 0d0) sqzp1 = -sqzp1
      cdasin = pi2 - ci * cdlog (z + sqzp1*cdsqrt(z-1d0))
c
      if (dble(cdasin) .gt. pi2) cdasin = pi - cdasin
      if (dble(cdasin) .le. (-pi2)) cdasin = -pi - cdasin
      if (dble(zinp) .lt. 0d0) cdasin = -cdasin
c
      end
      double complex function cdacos (z)
c Adapted from ...  w. fullerton, c3, los alamos scientific lab.
C     implicit none
      double complex z, cdasin
      double precision pi2
      external cdasin
      data pi2 /1.5707963267 9489661923d0/

      cdacos = pi2 - cdasin (z)
      end
C ... test dcacos
C      implicit none
C      double complex z,fz,cdacos
C
C      z = (1d0,2d0)
C      fz = cdacos(z)
C
C      print *, z, fz
C      z = fz
C      fz = cdcos(z)
C      print *, z, fz
C      end

