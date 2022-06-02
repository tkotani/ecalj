! FCPP#define F90 1
subroutine ghigbl(pg,ph,q,rsmg,rsmh,eg,eh,nlmg,nlmh,kmax, &
     ndim1,ndim2,k0,cg,indxcg,jcg,cy,s,ds)!,slat
  !- Block of integrals between smooth hankels and gaussians with some
  !  power of the laplace operator, and gradients.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   pg    :Function it expansed is at pg; see Remarks
  !i   ph    :Function is centered at ph; see Remarks
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
  !i   ndim1 :leading dimension of s,ds
  !i   ndim2 :second dimension of s,ds
  !i   k0    :third dimenson of s,ds
  !i   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
  !i   indxcg:index for Clebsch Gordon coefficients
  !i   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
  !i   cy    :Normalization constants for spherical harmonics
  !i   slat  :struct containing information about the lattice
  !o Outputs
  !o   s     :integrals of gaussian and Hankels; see Remarks
  !o   ds    :gradients of s; see Remarks
  !r Remarks
  !r   s(L,M,k) contains integral of G_L^*(r-pg) (lap)^k H_M(r-ph)
  !r  ds((L,M,k) contains grad s wrt ph; take negative for grad wrt pg.
  !u Updates
  !u   25 May 00 Made rsmh,eh l-dependent
  !u   25 May 00 Adapted from nfp ghig_bl.f
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: k0,kmax,ndim1,ndim2,nlmg,nlmh,jcg(1),indxcg(1)
  double precision :: rsmg,rsmh(1),eg,eh(1)
  double precision :: ph(3),pg(3),cg(1),cy(1),q(3)!,slat(1)
  double complex s(ndim1,ndim2,0:k0),ds(ndim1,ndim2,0:k0,3)
  ! ... Local parameters
  integer :: icg,icg1,icg2,ii,ilg,ilh,ilm,ilm1,ilm2,indx,ip,jlm,k,ktop, &
       ktop0,l1,l2,lg,lh,ll,lm,lmaxg,lmaxh,lmaxx,m,nlm0,nlmx
  double precision :: ee,fac,gamg,gamh,rsmx,dr(3),e,rsm
  complex(8),allocatable:: hkl(:,:),dhkl(:,:,:)
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
  ktop0 = ktop+1
  nlm0  = (lmaxx+2)**2
  allocate( hkl(0:ktop0,nlm0),dhkl(0:ktop0,nlm0,3))
  do    k = 0, kmax
     do    jlm = 1, nlmh
        do    ilm = 1, nlmg
           s(ilm,jlm,k) = dcmplx(0d0,0d0)
           ds(ilm,jlm,k,1) = dcmplx(0d0,0d0)
           ds(ilm,jlm,k,2) = dcmplx(0d0,0d0)
           ds(ilm,jlm,k,3) = dcmplx(0d0,0d0)
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
     call hklgbl(dr,rsmx,e,q,ktop,nlmx,ktop0,nlm0,cy,hkl,dhkl)!,slat

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
                 s(ilg,ilh,ip)    = s(ilg,ilh,ip)    + fac*hkl(k+ip,ilm)
                 ds(ilg,ilh,ip,1) = ds(ilg,ilh,ip,1) + fac*dhkl(k+ip,ilm,1)
                 ds(ilg,ilh,ip,2) = ds(ilg,ilh,ip,2) + fac*dhkl(k+ip,ilm,2)
                 ds(ilg,ilh,ip,3) = ds(ilg,ilh,ip,3) + fac*dhkl(k+ip,ilm,3)
12            enddo
11         enddo
111     enddo
1111 enddo
20 enddo
  deallocate(hkl,dhkl)
end subroutine ghigbl

