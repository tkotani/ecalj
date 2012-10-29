      subroutine sylmnc(c,lmx)
c  normalization constants for the spherical harmonics
c  use together with sylm (from asw package)
      implicit real*8 (a-h,p-z), integer(o)
      dimension c(*) ! MIZUHO-IR
      tpi=8.d0*datan(1.d0)
      fpi=2.d0*tpi
      y0=1.d0/dsqrt(fpi)
      n=(lmx+1)**2
      c(1)=y0
      do 2 l=1,lmx
        lp1=l+1
        tlp1=l+lp1
        lav=l*lp1+1
        c(lav)=dsqrt(tlp1/fpi)
      do 2 m=1,l
        n2=lp1-m
        n1=n2+1
        n3=l+m
        fn2=dfloat(n2)
        do 1 i=n1,n3
          fn2=fn2*dfloat(i)
    1   continue
        i1=lav+m
        i2=lav-m
        c(i1)=dsqrt(tlp1/(fn2*tpi))
        c(i2)=c(i1)
    2 continue
      return
      end
