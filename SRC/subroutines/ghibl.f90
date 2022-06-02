! FCPP#define F90 1
subroutine ghibl(pg,ph,q,rsmg,rsmh,eg,eh,nlmg,nlmh,kmax, &
     ndim1,ndim2,cg,indxcg,jcg,cy,s)!,slat
  !- Block of integrals between smooth hankels and gaussians with some power
  !  of the laplace operator.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ph    :Function is centered at ph; see Remarks
  !i   pg    :Function it expansed is at pg; see Remarks
  !i   q     :wave number for Bloch sum
  !i   rsmg  :smoothing radius of gaussian
  !i   rsmh  :vector of l-dependent smoothing radii of smoothed hankel
  !i         :rsmh must be specified for 1..ll(nlmh)+1
  !i   eg    :gkL scaled by exp(e*rsm**2/4)
  !i   eh    :vector of l-dependent energies of smoothed Hankel
  !i         :eh must be specified for 1..ll(nlmh)+1
  !i   nlmg  :L-cutoff for P_kL expansion
  !i   nlmh  :L-cutoff for smoothed Hankel functions
  !i   kmax  :polynomial cutoff
  !i   ndim1 :leading dimension of s
  !i   ndim2 :second dimension of s
  !i   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
  !i   indxcg:index for Clebsch Gordon coefficients
  !i   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
  !i   cy    :Normalization constants for spherical harmonics
  !i   slat  :struct containing information about the lattice
  !o Outputs
  !o   s     :integrals of gaussian and Hankels; see Remarks
  !r Remarks
  !r   s(L,M,k) contains integral of G_L^*(r-pg) (lap)^k H_M(r-ph)
  !r   See J. Math. Phys. {\bf 39},3393 (1998), Eq. 8.4
  !u Updates
  !u   18 May 00 Made rsmh,eh l-dependent
  !u   24 Apr 00 Adapted from nfp ghi_bl.f
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: nlmg,nlmh,kmax,ndim1,ndim2,jcg(1),indxcg(1)
  double precision :: rsmg,rsmh(1),eg,eh(1)
  double precision :: ph(3),pg(3),cg(1),cy(1),q(3)!,slat(1)
  double complex s(ndim1,ndim2,0:*)
  ! ... Local parameters
  integer :: nlm0,ktop0,icg,icg1,icg2,ii,ilg,ilh,ilm,indx,ip,jlm,k, &
       ktop,lg,lh,ll,lm,lmaxg,lmaxh,lmaxx,m,nlmx,l1,l2,ilm1,ilm2
  ! FCPP#if F90
  complex(8),allocatable:: hkl(:,:)
  ! FCPP#else
  ! FCPP      parameter (nlm0=144, ktop0=21)
  ! FCPP      double complex hkl(0:ktop0,nlm0)
  ! FCPP#endif
  double precision :: ee,fac,gamg,gamh,rsmx,dr(3),e,rsm

  if (nlmh == 0 .OR. nlmg == 0) return

  ! ... rsmh- and eh- independent setup
  do  1  m = 1, 3
     dr(m) = pg(m)-ph(m)
1 enddo
  lmaxh = ll(nlmh)
  lmaxg = ll(nlmg)
  lmaxx = lmaxg+lmaxh
  nlmx = (lmaxx+1)**2
  ktop = max0(lmaxg,lmaxh)+kmax
  ktop0 = ktop
  nlm0  = nlmx
  allocate( hkl(0:ktop0,nlm0))
  do    k = 0, kmax
     do    jlm = 1, nlmh
        do    ilm = 1, nlmg
           s(ilm,jlm,k) = dcmplx(0d0,0d0)
        enddo
     enddo
  enddo

  ! --- Loop over sequences of l with a common rsm,e ---
  l2 = -1
  do  20  l1 = 0, lmaxh
     if (l1 <= l2) goto 20
     call gtbsl2(l1,lmaxh,eh,rsmh,l2)
     rsm  = rsmh(l1+1)
     e    = eh(l1+1)
     if (rsm <= 0 .OR. e > 0) goto 20
     ilm1 = l1**2+1
     ilm2 = (l2+1)**2
     lmaxx= lmaxg+l2
     nlmx = (lmaxx+1)**2
     gamh = 0.25d0*rsm*rsm
     gamg = 0.25d0*rsmg*rsmg
     rsmx = 2d0*dsqrt(gamg+gamh)
     ktop = max0(lmaxg,l2)+kmax
     call hklbl(dr,rsmx,e,q,ktop,nlmx,ktop0,cy,hkl)!,slat

     !   ... Combine with Clebsch-Gordan coefficients
     ee = dexp(gamg*(eg-e))
     do  1111  ilg = 1, nlmg
        lg = ll(ilg)
        do  111  ilh = ilm1, ilm2
           lh = ll(ilh)
           ii = max0(ilg,ilh)
           indx = (ii*(ii-1))/2 + min0(ilg,ilh)
           icg1 = indxcg(indx)
           icg2 = indxcg(indx+1)-1
           do  11  icg = icg1, icg2
              ilm = jcg(icg)
              lm = ll(ilm)
              k = (lg+lh-lm)/2
              fac = ee*(-1d0)**lg*cg(icg)
              do  12  ip = 0, kmax
                 s(ilg,ilh,ip) = s(ilg,ilh,ip) + fac*hkl(k+ip,ilm)
12            enddo
11         enddo
111     enddo
1111 enddo
20 enddo
end subroutine ghibl

