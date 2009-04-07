      subroutine sylm(r,yl,lmx,r2s)
C- Generate unnormalized spherical harmonic polynomials
C ----------------------------------------------------------------
Ci Inputs
Ci   r     :vector with 3 components
Ci   lmax  :maximum l for which ylm is calculated
Co Outputs:
Co   ylm   :unnormalized spherical harmonic polynomials
Co   r2s   :length of dr**2
Cr Remarks:
Cr   polar axis along 001 axis. (adapted from ASW programs)
Cr   use together with sylmnc:
Cr   The true spherical harmonic is:  Y_L = ylm(lm,r) / r^l
Cr   The first 9 values are for ylm/cy:
Cr     l=0:                 1
Cr     l=1:        y        z       x
Cr     l=2:  6xy  3yz  3zz/2-rr/2  3xz  3(xx-yy)
Cr   Factors cy are (F = 4*pi)
Cr     l=0:              sqrt(1/4F)
Cr     l=1:              sqrt(3/F)   sqrt(3/F)  sqrt(3/F)
Cr     l=2:  sqrt(5/12F) sqrt(5/3F)  sqrt(5/F)  sqrt(5/3F) sqrt(5/12F)  
Cu Updates
Cu   18 Sep 06 Set tolerance to smaller number for gen. Ylm
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer lmx
      double precision r2s
      double precision r(3),yl(*)
C Local parameters
      integer i,l,lav,lavml,lavmm,lavpl,lavpm,lm1,lmm,lp1,m,mp1,n,nt
      double precision r2,st,x,y,z,z2
      double precision c(15),s(15),p(15,15)
      equivalence (x,c(2)),(y,s(2)),(z,p(2,1))
      data c(1),s(1),p(1,1),p(2,2) /1.d0,0.d0,1.d0,1.d0/

      n = (lmx+1)**2
      yl(1) = 1.d0
      x = r(1)
      y = r(2)
      z = r(3)
      st = x*x + y*y
      z2 = z*z
      r2 = st+z2
      r2s = r2
      if (n .lt. 2) return
      if (r2 .gt. 1d-28) goto 1
      do  6  i = 2, n
    6 yl(i) = 0.d0
      return
    1 yl(2) = y
      yl(3) = z
      yl(4) = x
      nt = 1
      do  5  l = 2, lmx
        lp1 = l+1
        lm1 = l-1
        lav = l*lp1 + 1
        p(lp1,1) = ((l+lm1)*z*p(l,1) - lm1*r2*p(lm1,1)) / l
        yl(lav) = p(lp1,1)
        nt = nt+2
        p(lp1,lp1) = p(l,l)*nt
        c(lp1) = x*c(l) - y*s(l)
        s(lp1) = x*s(l) + y*c(l)
        lavpl = lav+l
        yl(lavpl) = p(lp1,lp1)*c(lp1)
        lavml = lav-l
        yl(lavml) = p(lp1,lp1)*s(lp1)
        if (st .gt. z2) goto 3
        do  2  lmm = 1, lm1
          m = l-lmm
          lavpm = lav+m
          lavmm = lav-m
          mp1 = m+1
          p(lp1,mp1) = (r2*(l+m)*p(l,mp1) - st*p(lp1,mp1+1)) / (z*(l-m))
          yl(lavpm) = p(lp1,mp1)*c(mp1)
          yl(lavmm) = p(lp1,mp1)*s(mp1)
    2   continue
        goto 5
    3   do  4  m = 1, lm1
          mp1 = m+1
          lavpm = lav+m
          lavmm = lav-m
          p(lp1,mp1) = ((lm1+m)*r2*p(l,m) - (lp1-m)*z*p(lp1,m))/st
          yl(lavpm) = p(lp1,mp1)*c(mp1)
          yl(lavmm) = p(lp1,mp1)*s(mp1)
    4   continue
    5 continue
      end

