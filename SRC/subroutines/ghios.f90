      subroutine ghios(rsmg,rsmh,eg,eh,nlmh,kmax,ndim,s)
C- Integrals between 3D Gaussians and smooth Hankels at same site.
C ----------------------------------------------------------------------
Ci Inputs
Ci   rsmg  :smoothing radius of gaussian
Ci   rsmh  :vector of l-dependent smoothing radii of smoothed hankel
Ci         :rsmh must be specified for 1..ll(nlmh)+1
Ci   eg    :gkL scaled by exp(eg*rsm**2/4)
Ci   eh    :vector of l-dependent energies of smoothed Hankel
Ci         :eh must be specified for 1..ll(nlmh)+1
Ci   nlmh  :L-cutoff for smoothed Hankel functions and P_kL
Ci   kmax  :polynomial cutoff
Ci   ndim  :nl*nl*nbas
Co Outputs
Co   s     :real-space structure constants c(k,M,L); see Remarks
Co         :s(ilm,*) for rsmh(ll(ilm)) <= 0 are UNTOUCHED
Cr Remarks
Cr   s(L,M,k) contains integral of G_L^*(r) (lap)^k H_M(r)
Cr
Cr   Only diagonal elements ilm=jlm are nonzero and returned.
Cr   Equivalent to ghiml called for pg=ph.
Cu Updates
Cu   18 May 00 Made rsmh,eh l-dependent
Cu   24 Apr 00 Adapted from nfp ghi_os.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer kmax,ndim,nlmh
      double precision eg,eh(*),rsmg,rsmh(*),s(ndim,0:kmax)
C ... Local parameters
      integer ktop0,ilm,k,ktop,l,ll,lmaxh,l1,l2,ilm1,ilm2
      parameter( ktop0=50 )
      double precision fac,fpi,gamg,gamh,rsmx,y0,h0k(0:ktop0),rsm,e

      if (nlmh .eq. 0) return
      fpi = 16d0*datan(1d0)
      y0 = 1d0/dsqrt(fpi)
      lmaxh = ll(nlmh)

C --- Loop over sequences of l with a common rsm,e ---
      l2 = -1
      do  20  l1  = 0, lmaxh
        if (l1 .le. l2) goto 20
        call gtbsl2(l1,lmaxh,eh,rsmh,l2)
        rsm  = rsmh(l1+1)
        e    = eh(l1+1)
        if (rsm .le. 0 .or. e .gt. 0) goto 20
        ilm1 = l1**2+1
        ilm2 = (l2+1)**2
        gamh = 0.25d0*rsm*rsm
        gamg = 0.25d0*rsmg*rsmg
        rsmx = 2d0*dsqrt(gamg+gamh)

C   ... Make hankels for l=0 and k=0..kmax
        ktop = l2+kmax
        if(ktop .gt. ktop0) call rxi('ghios: increase ktop0, need',ktop)
        call hklos(rsmx,e,ktop,h0k)

C   ... Evaluate what is left of Clebsch-Gordan sum
        fac = y0*dexp(gamg*(eg-e))
        do  ilm = ilm1, ilm2
          l = ll(ilm)
          do  k = 0, kmax
            s(ilm,k) = fac * (-1)**l * h0k(k+l)
          enddo
        enddo

   20 continue

      end

