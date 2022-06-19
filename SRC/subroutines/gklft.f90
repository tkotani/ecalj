subroutine gklft(v,rsm,e,tau,alat,kmax,nlm,k0,cy,gkl)
  !- Returns one Fourier component of gaussians Gkl centered at tau.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   v     :reciprocal vector for which FT calc., units of 2*pi/alat
  !i   rsm   :smoothing radius
  !i   e     :gaussian scaled by exp(e*rsm**2/4)
  !i   tau   :center of gaussian, in units of alat
  !i   alat  :length scale of v, tau
  !i   kmax  :polynomial cutoff for gaussian
  !i   nlm   :L-cutoff for gaussian
  !i   k0    :leading dimension of gkl
  !i   cy    :Normalization constants for spherical harmonics
  !o Outputs
  !o   gkl   :FT of gaussian centered at tau
  !r Remarks
  !u Updates
  !u   23 Apr 00 Adapted from nfp gkl_ft.f
  ! ----------------------------------------------------------------------
  implicit none
  integer:: k0,nlm, ilm,k,l,ll,lmax,m,kmax
  real(8):: alat,e,rsm,v(3),cy(1),tau(3), aa,fac,gam,scalp,v2,yl(nlm),vv(3)
  complex(8):: phaseaa,qfac,img=(0d0,1d0), gkl(0:k0,nlm)
  real(8),parameter::pi = 4d0*datan(1d0)
  if (nlm == 0) return
  gam = 0.25d0*rsm**2
  call sylm(v(:)*2d0*pi/alat,yl, ll(nlm), v2)
  gkl(0,:) = exp(gam*(e-v2)) * exp(-2d0*img*pi*sum(tau*v)) &
       * [((-img)**ll(ilm)*cy(ilm)*yl(ilm),  ilm=1,nlm)]
  do k = 1, kmax
     gkl(k,:) = (-v2)**k *gkl(0,:)
  enddo
end subroutine gklft

