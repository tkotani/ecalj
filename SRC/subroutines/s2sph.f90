subroutine s2sph(opt,nla,nlb,s,nds1,nds2,ndss1,ndss2,sph)
  !- Rotate a structure matrix from real to spherical harmonics
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   opt   :1s digit governs storage of the imaginary part of sph
  !i          0: real and imaginary parts separated
  !i          1: hamiltonian returned as complex*16
  !i          2: real and imaginary parts by columns
  !i         10s digit
  !i           0 if s is real, 1 if s is complex
  !i   nla   :number of l's to rotate in the 1st (augmentation) dimension
  !i   nlb   :number of l's to rotate in the 2st (basis) dimension
  !i   s     :real-space structure constant matrix
  !i   nds1  :leading dimension of s
  !i   nds2  :second dimension of s; not used unless s is complex
  !i   ndss1 :leading dimension of sph
  !i   ndss2 :second dimension of sph
  !o Outputs
  !o   sph   :u s u+, where u rotates s from real to spherical harmonics
  !o          See cb2sph for the explicit construction of u.
  !r Remarks
  !r   This routine rotates s to sph through a rotation u, s.t.
  !r      sph = u s u+, where u does the following rotation.
  !r   Let R_lm = real harmonic; Y_lm = spherical harmonic.  Then
  !r     Y_l,m = 1/sqrt(2)      (R_l,-m + i * R_l,m),  m<0
  !r     Y_l,m = R_l,m,                                m=0
  !r     Y_l,m = (-1)^m/sqrt(2) (R_l,-m - i * R_l,m),  m>0
  !r
  !r   In particular
  !r    l  m        R_l,m                 Y_l,m
  !r    1  -1   sqrt(3/4/pi)*y            sqrt(3/4/pi/2)(x + i*y)
  !r    1   1   sqrt(3/4/pi)*x            sqrt(3/4/pi/2)(-x + i*y)
  !r    2  -2   2*sqrt(15/16/pi)*x*y      sqrt(15/32/pi)*(x + i*y)^2
  !r    2   2   sqrt(15/16/pi)*(x*x-y*y)
  !r
  !r   The Y_l,m here are the same as those in Jackson, except that
  !r       Y_l,m = Y_l,-m (Jackson)
  !r
  !r   Rotation of the strux S.  S are one-center expansion coefficients
  !r   such as the expansion of function H around another site:
  !r      H_L(r-dr) = (h R_L) = sum_L' S_LL' J_L' = S^R_LL'(dr) (j R_L')
  !r   If H and J are to be rotated from real harmonics R_L to spherical
  !r   harmonics Y_L, we have
  !r     (h R_L) =  S^R_LL' (j R_L')
  !r     (h Y_L) =  S^Y_LL' (j Y_L')
  !r   Then
  !r     (h u R_L) =  S^Y_LL' (j u R_L')
  !r   So
  !r     (h R_L) =  u^-1 S^Y_LL' u (j R_L') -> u^-1 S^Y_LL' u = S^R_LL'
  !r   Therefore
  !r     S^Y_LL' = u S^R_LL' u^-1 ->  S^Y = u S^R u+
  !r
  !r   Other remarks.
  !r   s2sph performs the same function as cb2sph, but is much faster.
  !b Bugs
  !b   Not checked for the case nla ne nlb
  !u Updates
  !u   10 Oct 03 Extended to s being complex (s in kcplx mode 0)
  ! ----------------------------------------------------------------------
  implicit none
  integer :: opt,nla,nlb,nds1,nds2,ndss1,ndss2
  double precision :: s(nds1,nds2,2),sph(ndss1,ndss2,2)
  logical :: lscplx
  integer :: l,lma,lmb,nl2,lc,m,ndwk,ladd
  double precision :: sr2,sm2,spmr,spmi,sppr,sppi
  parameter (ndwk=64)
  double precision :: wk(ndwk,ndwk,2)
  m = max(nla*nla,nlb*nlb)
  if (m > ndwk) call rxi('s2sph: increase ndwk to ge',m)
  sr2 = dsqrt(0.5d0)
  lscplx = mod(opt/10,10).ne.0
  ladd = 0
  if (mod(opt/100,10) /= 0) ladd = 1
  ! ... Make us = u*s, where u is defined as in cb2sph
  nl2 = nlb*nlb
  do  10  l = 0, nla-1
     do  14  lmb = 1, nl2
        lc = (l+1)**2-l
        wk(lc,lmb,1) = s(lc,lmb,1)
        wk(lc,lmb,2) = 0
        if (lscplx) wk(lc,lmb,2) = s(lc,lmb,2)
        sm2 = -sr2
        do  12  m = 1, l
           wk(lc-m,lmb,1) =  sr2*s(lc+m,lmb,1)
           wk(lc-m,lmb,2) =  sr2*s(lc-m,lmb,1)
           wk(lc+m,lmb,1) =  sm2*s(lc+m,lmb,1)
           wk(lc+m,lmb,2) = -sm2*s(lc-m,lmb,1)
           if (lscplx) then
              wk(lc-m,lmb,1) = wk(lc-m,lmb,1) - sr2*s(lc+m,lmb,2)
              wk(lc-m,lmb,2) = wk(lc-m,lmb,2) + sr2*s(lc-m,lmb,2)
              wk(lc+m,lmb,1) = wk(lc+m,lmb,1) + sm2*s(lc+m,lmb,2)
              wk(lc+m,lmb,2) = wk(lc+m,lmb,2) - sm2*s(lc-m,lmb,2)
           endif
           sm2 = -sm2
12      enddo
14   enddo
10 enddo
  ! ... Overwrite wk with u*s*u+, where u is defined as in cb2wk
  nl2 = nla*nla
  do  20  l = 0, nlb-1
     do  24  lma = 1, nl2
        lc = (l+1)**2-l
        wk(lma,lc,1) = wk(lma,lc,1)
        wk(lma,lc,2) = wk(lma,lc,2)
        sm2 = -sr2
        do  22  m = 1, l
           spmr =  sr2*(wk(lma,lc+m,1) + wk(lma,lc-m,2))
           spmi = -sr2*(wk(lma,lc-m,1) - wk(lma,lc+m,2))
           sppr =  sm2*(wk(lma,lc+m,1) - wk(lma,lc-m,2))
           sppi =  sm2*(wk(lma,lc-m,1) + wk(lma,lc+m,2))
           wk(lma,lc-m,1) = spmr
           wk(lma,lc-m,2) = spmi
           wk(lma,lc+m,1) = sppr
           wk(lma,lc+m,2) = sppi
           sm2 = -sm2
22      enddo
24   enddo
20 enddo
  if (mod(opt,10) == 0) then
     call ymscop0(nla*nla,nlb*nlb,ndwk,ndss1,  wk,ndwk*ndwk,sph,ndss1*ndss2)
  elseif (mod(opt,10) == 1) then
     call rx('s2sph not ready for complex*16 format')
  else
     call ymscop0(nla*nla,nlb*nlb,ndwk,ndss1*2,wk,ndwk*ndwk,sph,ndss1)
  endif
end subroutine s2sph

subroutine ymscop0(nlma,nlmb,ndas,ndad, src,offsi,dest,offdi)
  implicit none
  integer:: nlma,nlmb,ndas,ndad,offsi,offdi,ia,ib
  real(8):: src(ndas,1),dest(ndad,1)
  if (nlma <= 0 .OR. nlmb <= 0) return
  do  ib = 1, nlmb
     do  ia = 1, nlma
        dest(ia, ib) = src(ia, ib)
        dest(ia+offdi,ib) = src(ia+offsi,ib)
     enddo
  enddo
end subroutine ymscop0
