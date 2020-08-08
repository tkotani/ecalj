      subroutine jxpos(rsmg,ej,kmax,lmax,k0,c)
C-  Coefficients to expand Bessel function into a sum of P-kl's
C ----------------------------------------------------------------------
Ci Inputs
Ci   rsmg  :smoothing radius of gaussian
Ci   ej    :energy of Bessel function
Ci   kmax  :polynomial cutoff
Ci   lmax  :maximum l for a given site
Ci   k0    :leading dimension of coefficient array c
Co Outputs
Co   c     :structure constants c(k,M,L); see Remarks
Cr Remarks
Cr   Expansion is:  J_l = sum(k) c(k,l) * P_kl
Cr   Polynomials P_kL are centered at the same site.
Cr   Radius rsmg determines the range of validity for the expansion.
Cr   Straightforward evaluation would be:
Cr      xx = dexp(-gamma*ej)
Cr      c(k,l) = xx*(-ej)**k/(((4*a*a)**k)*(a**l)*fac_k*dfac_l)
Cr   See also: hxpos.f
Cu Updates
Cu   02 May 02 Adapted from nfp jxpos.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer kmax,lmax,k0
      double precision rsmg,c(0:k0,0:lmax),ej
C ... Local parameters
      integer l,k
      double precision gamma,a,c0l

      gamma = rsmg*rsmg/4d0
      a = 1d0/rsmg
      c0l = dexp(-gamma*ej)
      do  l = 0, lmax
        c(0,l) = c0l
        do  k = 1, kmax
          c(k,l) = c(k-1,l)*(-ej)/(4*a*a*k)
        enddo
        c0l = c0l/(a*(2*l+3))
      enddo

c|      write(stdo,*) 'Coeffs:'
c|      do  l = 0, lmax
c|        write(stdo,200) l,(c(k,l),k=0,kmax)
c|  200   format(' l=',i2,8f13.7)
c|      enddo


      end

