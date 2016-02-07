      subroutine dlmtor(r,a,lmax,f,fbar)
c  radial part of damped lmtos f and fbar, l=0 to lmax
      implicit real*8 (a-h,p-z), integer(o)
      dimension f(0:lmax),fbar(0:lmax)
      fbsrpi=0.564189835d0
      pi=4.d0*datan(1.d0)
      fbsqpi=1.d0/dsqrt(pi)
      z=a*r


      emz2=dexp(-z*z)
      erfc0=derfc(z)
      erfc1=-z*erfc0+fbsrpi*emz2
      erfc2=-0.5d0*z*erfc1+0.25d0*erfc0
      f(0)=erfc0/r
      fbar(0)=-erfc2/(a*a*r)
      ta2r=2.d0*a*a*r
      g=2.d0*a*emz2*fbsrpi/r
      flm2=fbsrpi*emz2/z-erfc0
      do 10 l=1,lmax
        f(l)=((l+l-1)/r)*f(l-1)+g
        fbar(l)=((l+l-1)/r)*fbar(l-1)-flm2
        flm2=f(l-1)
        g=g*ta2r
   10 continue
      return
      end
