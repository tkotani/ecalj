      subroutine s2sph(opt,nla,nlb,s,nds1,nds2,ndss1,ndss2,sph)
C- Rotate a structure matrix from real to spherical harmonics
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :1s digit governs storage of the imaginary part of sph
Ci          0: real and imaginary parts separated
Ci          1: hamiltonian returned as complex*16
Ci          2: real and imaginary parts by columns
Ci         10s digit
Ci           0 if s is real, 1 if s is complex
Ci   nla   :number of l's to rotate in the 1st (augmentation) dimension
Ci   nlb   :number of l's to rotate in the 2st (basis) dimension
Ci   s     :real-space structure constant matrix
Ci   nds1  :leading dimension of s
Ci   nds2  :second dimension of s; not used unless s is complex
Ci   ndss1 :leading dimension of sph
Ci   ndss2 :second dimension of sph
Co Outputs
Co   sph   :u s u+, where u rotates s from real to spherical harmonics
Co          See cb2sph for the explicit construction of u.
Cr Remarks
Cr   This routine rotates s to sph through a rotation u, s.t.
Cr      sph = u s u+, where u does the following rotation.
Cr   Let R_lm = real harmonic; Y_lm = spherical harmonic.  Then
Cr     Y_l,m = 1/sqrt(2)      (R_l,-m + i * R_l,m),  m<0
Cr     Y_l,m = R_l,m,                                m=0
Cr     Y_l,m = (-1)^m/sqrt(2) (R_l,-m - i * R_l,m),  m>0
Cr
Cr   In particular
Cr    l  m        R_l,m                 Y_l,m
Cr    1  -1   sqrt(3/4/pi)*y            sqrt(3/4/pi/2)(x + i*y)
Cr    1   1   sqrt(3/4/pi)*x            sqrt(3/4/pi/2)(-x + i*y)
Cr    2  -2   2*sqrt(15/16/pi)*x*y      sqrt(15/32/pi)*(x + i*y)^2
Cr    2   2   sqrt(15/16/pi)*(x*x-y*y)
Cr
Cr   The Y_l,m here are the same as those in Jackson, except that
Cr       Y_l,m = Y_l,-m (Jackson)
Cr
Cr   Rotation of the strux S.  S are one-center expansion coefficients
Cr   such as the expansion of function H around another site:
Cr      H_L(r-dr) = (h R_L) = sum_L' S_LL' J_L' = S^R_LL'(dr) (j R_L')
Cr   If H and J are to be rotated from real harmonics R_L to spherical
Cr   harmonics Y_L, we have
Cr     (h R_L) =  S^R_LL' (j R_L')
Cr     (h Y_L) =  S^Y_LL' (j Y_L')
Cr   Then
Cr     (h u R_L) =  S^Y_LL' (j u R_L')
Cr   So
Cr     (h R_L) =  u^-1 S^Y_LL' u (j R_L') -> u^-1 S^Y_LL' u = S^R_LL'
Cr   Therefore
Cr     S^Y_LL' = u S^R_LL' u^-1 ->  S^Y = u S^R u+
Cr
Cr   Other remarks.
Cr   s2sph performs the same function as cb2sph, but is much faster.
Cb Bugs
Cb   Not checked for the case nla ne nlb
Cu Updates
Cu   10 Oct 03 Extended to s being complex (s in kcplx mode 0)
C ----------------------------------------------------------------------
      implicit none
      integer opt,nla,nlb,nds1,nds2,ndss1,ndss2
      double precision s(nds1,nds2,2),sph(ndss1,ndss2,2)
      logical lscplx
      integer l,lma,lmb,nl2,lc,m,ndwk,ladd
      double precision sr2,sm2,spmr,spmi,sppr,sppi
      parameter (ndwk=64)
      double precision wk(ndwk,ndwk,2)
      m = max(nla*nla,nlb*nlb)
      if (m .gt. ndwk) call rxi('s2sph: increase ndwk to ge',m)
      sr2 = dsqrt(0.5d0)
      lscplx = mod(opt/10,10).ne.0
      ladd = 0
      if (mod(opt/100,10) .ne. 0) ladd = 1
C ... Make us = u*s, where u is defined as in cb2sph
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
   12     continue
   14   continue
   10 continue
C ... Overwrite wk with u*s*u+, where u is defined as in cb2wk
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
   22     continue
   24   continue
   20 continue
      if (mod(opt,10) .eq. 0) then
        call ymscop0(nla*nla,nlb*nlb,ndwk,ndss1,  wk,ndwk*ndwk,sph,ndss1*ndss2)
      elseif (mod(opt,10) .eq. 1) then
        call rx('s2sph not ready for complex*16 format')
      else
        call ymscop0(nla*nla,nlb*nlb,ndwk,ndss1*2,wk,ndwk*ndwk,sph,ndss1)
      endif
      end
      
      subroutine ymscop0(nlma,nlmb,ndas,ndad, src,offsi,dest,offdi)
      implicit none
      integer:: nlma,nlmb,ndas,ndad,offsi,offdi,ia,ib
      real(8):: src(ndas,1),dest(ndad,1)
       if (nlma .le. 0 .or. nlmb .le. 0) return
       do  ib = 1, nlmb
       do  ia = 1, nlma
           dest(ia, ib) = src(ia, ib)
           dest(ia+offdi,ib) = src(ia+offsi,ib)
       enddo
       enddo
      end
