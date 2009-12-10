      subroutine makdla(kbak,l,hcr,slo,val,rmax,phia,dla)
C- Gets logarithmic derivative of wave function at hard core radius
C ----------------------------------------------------------------------
Ci Inputs:
Ci   kbak  :kinetic energy for back extrapolation
Ci   l     :angular momentum
Ci   hcr   :hard-core screening radius a, in atomic units
Ci   slo   :slope of g at rmax, where g = rmax * radial w.f.
Ci   val   :value of g at rmax, where g = rmax * radial w.f.
Ci   rmax  :radius at which val,slo are evaluated, in atomic units
Co Outputs:
Co   phia  :value of radial w.f. at a
Co         :If a=rmax, phia = val/rmax
Co   dla   :logarithmic derivative at a, i.e. a/phi*(dphi/dr)_a
Co         :If a=rmax, dla = log deriv of phi = rmax*slo/val - 1
Cr Remarks:
Cr This was adapted from the Stuttgart third-generation LMTO package.
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer l
      double precision kbak,hcr,slo,val,rmax,phia,dla
C Local variables:
      integer nlmax
      parameter(nlmax=20)
      double precision er2,fi(0:nlmax+1),gi(0:nlmax+1),wn,wj,dlr,dj,dn,
     .                 sigma,phi,rdphia,rdfi,rdgi
C External calls:
      external bessl2

      er2 = kbak*rmax*rmax
      call bessl2(er2,0,l+1,fi(0),gi(0))
C     phi,dlr are value and logarithmic derivative r/phi dphi/dr at rmax
C     free dlr=dj
      phi =  val/rmax
      dlr =  rmax*slo/val - 1
C     dj,dn are the logarithmic derivatives of Bessel and Hankels
      dj  = l-fi(l+1)/fi(l)/(l+l+1)*er2
      dn  = l-gi(l+1)/gi(l)*(l+l+1)
C     wj,wn are the amounts of Bessel and Hankel making up phi0
      wj  = (dlr-dj)*fi(l)*phi
      wn  = (dlr-dn)*gi(l)*phi
      sigma = hcr/rmax
      er2 = kbak*hcr**2
      call bessl2(er2,0,l+1,fi(0),gi(0))
      rdgi = l*gi(l) - gi(l+1)*(l+l+1)
      rdfi = l*fi(l) - fi(l+1)/(l+l+1)*er2

      phia   = 2d0*(wn*fi(l)*sigma**l - wj*gi(l)*sigma**(-l-1))
      rdphia = 2d0*(wn*rdfi*sigma**l  - wj*rdgi*sigma**(-l-1))
      dla    = rdphia/phia
      end

