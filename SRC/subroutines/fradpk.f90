      subroutine fradpk(kmax,rsma,lmxa,nr,rofi,fp,xp,vp,dp)
C- Set up radial pkl functions, grad of same, vals and slopes.
C ----------------------------------------------------------------------
Ci Inputs
Ci   kmax  :polynomial cutoff in augmentation expansion
Ci   rsma  :augmentation smoothing radius
Ci   lmxa  :augmentation L-cutoff
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Co Outputs
Co   fp    :(Pkl) * r
Co   xp    :(Pkl') * r
Co   vp    :Pkl at rofi(nr) = fp/rofi(nr)
Co   dp    :Pkl' at rofi(nr)
Cr Remarks
Cr   Radial functions are multiplied by r, but vp and dp are not.
Cu Updates
Cu  16 May 00 Adapted from nfp frapkl.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer kmax,lmxa,nr
      double precision rofi(nr),fp(nr,0:lmxa,0:kmax),rsma,
     .xp(nr,0:lmxa,0:kmax),vp(0:lmxa,0:kmax),dp(0:lmxa,0:kmax)
C ... Local parameters
      integer k0,l0,l,k,i
      parameter ( k0=30, l0=7 )
      double precision pkl(0:k0,0:l0),rmt,fac,xx,r,rlp1,rr,rr2

      if (kmax+1.gt.k0) call rxi('fradpk: need k0 .ge.',kmax+1)
      if (lmxa+1.gt.l0) call rxi('fradpk: need l0 .ge.',lmxa+1)
      if (lmxa .eq. -1) return

C ... Values and slopes of the pkl
      rmt = rofi(nr)
      call radpkl(rmt,rsma,kmax,lmxa+1,k0,pkl)
      do  l = 0, lmxa
        do  k = 0, kmax
          fac = (2*k+2*l+3d0)/(2*l+3d0)
          xx = pkl(k,l)
          dp(l,k) = 2*rmt**(l+1)/rsma**2*(xx-pkl(k,l+1)*fac*rsma)
     .    + l*rmt**(l-1)*xx
          vp(l,k) = rmt**l*pkl(k,l)
        enddo
      enddo

C ... Set up the functions pkl(r) and derivative on radial mesh
      do  i = 1, nr
        r = rofi(i)
        call radpkl(r,rsma,kmax,lmxa+1,k0,pkl)
        do  l = 0, lmxa
          rlp1 = 0d0
          if ( r .gt. 0d0 ) rlp1 = r**(l+1)
          do  k = 0, kmax
            fp(i,l,k) = pkl(k,l)*rlp1
            fac = (2*k+2*l+3d0)/(2*l+3d0)
            xx = pkl(k,l)
            rr = 0d0
            if ( r .gt. 0d0 ) rr = r**l
            rr2 = rr * r * r
c           xp(i,l,k) = 2*r**(l+2)/rsma**2*(xx-pkl(k,l+1)*fac*rsma)
c    .         + l*r**l*xx
            xp(i,l,k) = 2*rr2/rsma**2*(xx-pkl(k,l+1)*fac*rsma)
     .      + l*rr*xx
          enddo
        enddo
      enddo

      end

