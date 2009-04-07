      subroutine makptf(avw,kbak,l,hcr,slo,val,rmax,ptf,beta)
C- Makes potential function and downfolding beta
C ----------------------------------------------------------------------
Ci Inputs:
Ci   kbak  :kinetic energy of the envelope function
Ci   l     :angular momentum
Ci   hcr   :hard core sphere radius
Ci   slo   :slope of g at rmax
Ci   val   :value of g at rmax
Ci   rmax  :MT radius Wigner-Seitz sphere radius, in atomic units
Co Outputs:
Co   ptf   :potential function P^0 = W{phi,K}/W{phi,J}
Co         :In 2nd gen LMTO, where kbak=0, using Andersen 1980's conv.
Co         :P^0 -> 2(2l+1) (avw/rmax)^(2l+1) (D+l+1) / (D-l)
Co         :and (P^alpha)^-1 = (P^0)^-1 - alpha
Co   beta  :finds the specific alpha that matches to w.f. at hcr,
Co          optimal for downfolding.
Cr Remarks:
Cr   The potential function is 1/beta, where beta is obtainable
Cr   from the back extrapolated function phi^0, which has the form
Cr
Cr     phi^0 = const * (j^0(kr) - beta n^0(kr))
Cr
Cr   See Lecture Notes in Physics 535, after Eq. 34.
Cr   This was adapted from the Stuttgart third-generation LMTO package.
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer l
      double precision beta,avw,kbak,ptf,hcr,slo,val,rmax
C Local variables:
      integer nlmax
      parameter(nlmax=20)
      double precision er2,fi(0:nlmax+1),gi(0:nlmax+1),wn,wj,dphi
C External calls:
      external  bessl2

      er2 = kbak*rmax*rmax
      call bessl2(er2,0,l+1,fi(0),gi(0))
      dphi = rmax*slo / val - 1.d0
C ... wj = (D{phi}-D{J}), wn = (D{phi}-D{K})
      wn = (dphi-l) + gi(l+1)/gi(l)*(l+l+1)
      wj = (dphi-l) +er2*fi(l+1)/fi(l)/(l+l+1)
      ptf = (avw/rmax)**(l+l+1)*gi(l)/fi(l)*wn/wj
      if (hcr .ne. 0) then
        er2 = kbak*hcr*hcr
        call bessl2(er2,0,l,fi(0),gi(0))
        beta = fi(l)/gi(l)*(hcr/avw)**(l+l+1)
      else
        beta = fi(l)/gi(l)
      endif
      end

