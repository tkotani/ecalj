      subroutine gklft(v,rsm,e,tau,alat,kmax,nlm,k0,cy,gkl)
C- Returns one Fourier component of gaussians Gkl centered at tau.
C ----------------------------------------------------------------------
Ci Inputs
Ci   v     :reciprocal vector for which FT calc., units of 2*pi/alat
Ci   rsm   :smoothing radius
Ci   e     :gaussian scaled by exp(e*rsm**2/4)
Ci   tau   :center of gaussian, in units of alat
Ci   alat  :length scale of v, tau
Ci   kmax  :polynomial cutoff for gaussian
Ci   nlm   :L-cutoff for gaussian
Ci   k0    :leading dimension of gkl
Ci   cy    :Normalization constants for spherical harmonics
Co Outputs
Co   gkl   :FT of gaussian centered at tau
Cr Remarks
Cu Updates
Cu   23 Apr 00 Adapted from nfp gkl_ft.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer k0,kmax,nlm
      double precision alat,e,rsm,v(3),cy(1),tau(3)
      double complex gkl(0:k0,nlm)
C ... Local parameters
      integer nlm0,ilm,k,l,ll,lmax,m
      parameter ( nlm0=121 )
      double precision aa,fac,gam,pi,scalp,v2,yl(nlm0),vv(3)
      double complex phase,qfac

      if (nlm .eq. 0) return
      if (nlm .gt. nlm0) call rxi('gklft: increase nlm0, need',nlm)
      if (kmax .gt. k0)  call rxi('gklft: increase k0, need',kmax)
      pi = 4d0*datan(1d0)
      lmax = ll(nlm)
      gam = 0.25d0*rsm**2
      vv(1) = v(1)*2*pi/alat
      vv(2) = v(2)*2*pi/alat
      vv(3) = v(3)*2*pi/alat
      call sylm(vv,yl,lmax,v2)
      aa = dexp(gam*(e-v2))
      scalp = -2*pi*(tau(1)*v(1)+tau(2)*v(2)+tau(3)*v(3))
      phase = dcmplx(dcos(scalp),dsin(scalp))

      qfac = dcmplx(0d0,1d0)
      ilm = 0
      do  l = 0, lmax
        qfac = qfac*dcmplx(0d0,-1d0)
        do m = -l,l
          ilm = ilm+1
          gkl(0,ilm) = phase*aa*qfac*cy(ilm)*yl(ilm)
        enddo
      enddo

      fac = 1d0
      do  k = 1, kmax
        fac = -v2*fac
        do  ilm = 1, nlm
          gkl(k,ilm) = fac*gkl(0,ilm)
        enddo
      enddo

      end

