      double complex function dcsqrt (z)
c april 1977 version.  w. fullerton, c3, los alamos scientific lab.
C     implicit none
#if FUJITSU
      complex*16 z
#else
      double complex z
#endif

      double precision cdabs,dsqrt,xtmp,ytmp,x,y,r
c
      x = dble  (z)
      y = dimag (z)
      r = cdabs (z)
c
      if (r.eq.0.d0) dcsqrt = (0.0d0, 0.0d0)
      if (r.eq.0.d0) return
c
      xtmp = dsqrt ((r+dabs(x))*0.5d0)
      ytmp = y*0.5d0/xtmp
c
      if (x.ge.0.d0) dcsqrt = dcmplx (xtmp, ytmp)
      if (y.eq.0.d0) y = 1.d0
      if (x.lt.0.d0) dcsqrt = dcmplx (dabs(ytmp), dsign(xtmp,y))
c
      return
      end

