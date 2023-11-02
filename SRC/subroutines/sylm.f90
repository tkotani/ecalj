subroutine sylm(r,yl,lmx,r2s)   !- Generate unnormalized spherical harmonic polynomials
  !i   r     :vector with 3 components
  !i   lmax  :maximum l for which ylm is calculated
  !o Outputs:
  !o   ylm   :unnormalized spherical harmonic polynomials
  !o   r2s   :length of dr**2
  !r Remarks:
  !r   polar axis along 001 axis. (adapted from ASW programs)
  !r   use together with sylmnc:
  !r   The true spherical harmonic is:  Y_L = ylm(lm,r) / r^l
  !r   The first 9 values are for ylm/cy:
  !r     l=0:                 1
  !r     l=1:        y        z       x
  !r     l=2:  6xy  3yz  3zz/2-rr/2  3xz  3(xx-yy)
  !r   Factors cy are (F = 4*pi)
  !r     l=0:              sqrt(1/4F)
  !r     l=1:              sqrt(3/F)   sqrt(3/F)  sqrt(3/F)
  !r     l=2:  sqrt(5/12F) sqrt(5/3F)  sqrt(5/F)  sqrt(5/3F) sqrt(5/12F)
   implicit none
  integer :: lmx
  double precision :: r2s
  double precision :: r(3),yl(*)
  integer :: i,l,lav,lavml,lavmm,lavpl,lavpm,lm1,lmm,lp1,m,mp1,n,nt
  double precision :: r2,st,x,y,z,z2
  double precision :: c(15),s(15),p(15,15)
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
  if (n < 2) return
  if (r2 > 1d-28) goto 1
  do  6  i = 2, n
     yl(i) = 0.d0
6 enddo
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
     if (st > z2) goto 3
     do  2  lmm = 1, lm1
        m = l-lmm
        lavpm = lav+m
        lavmm = lav-m
        mp1 = m+1
        p(lp1,mp1) = (r2*(l+m)*p(l,mp1) - st*p(lp1,mp1+1)) / (z*(l-m))
        yl(lavpm) = p(lp1,mp1)*c(mp1)
        yl(lavmm) = p(lp1,mp1)*s(mp1)
2    enddo
     goto 5
3    do  4  m = 1, lm1
        mp1 = m+1
        lavpm = lav+m
        lavmm = lav-m
        p(lp1,mp1) = ((lm1+m)*r2*p(l,m) - (lp1-m)*z*p(lp1,m))/st
        yl(lavpm) = p(lp1,mp1)*c(mp1)
        yl(lavmm) = p(lp1,mp1)*s(mp1)
4    enddo
5 enddo
end subroutine sylm
subroutine sylmnc(c,lmx)   !- Normalization constants for the spherical harmonics
  ! ----------------------------------------------------------------
  !i Inputs
  !i   lmx
  !o Outputs
  !o   C
  !r Remarks
  !r   use together with sylm (from ASW package)
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  integer :: lmx
  double precision :: c(1)
  ! Local parameters
  integer :: i,i1,i2,l,lav,lp1,m,n1,n2,n3
  double precision :: fn2,fpi,tlp1,tpi,y0
  tpi = 8.d0*datan(1.d0)
  fpi = 2.d0*tpi
  y0 = 1.d0/dsqrt(fpi)
  c(1) = y0
  do  21  l = 1, lmx
     lp1 = l+1
     tlp1 = l+lp1
     lav = l*lp1 + 1
     c(lav) = dsqrt(tlp1/fpi)
     do  2  m = 1, l
        n2 = lp1-m
        n1 = n2+1
        n3 = l+m
        fn2 = n2
        do  1  i = n1, n3
           fn2 = fn2*i
1       enddo
        i1 = lav+m
        i2 = lav-m
        c(i1) = dsqrt(tlp1/(fn2*tpi))
        c(i2) = c(i1)
2    enddo
21 enddo
end subroutine sylmnc


