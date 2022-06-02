! FCPP#define F90 1
subroutine gfigbl(pg,ph,rsmg,rsmh,nlmg,nlmh,kmax,ndim1,ndim2,kdim, &
     cg,indxcg,jcg,cy,s,ds)!,slat
  !- Integrals between smooth hankels(eh=0,q=0) and gaussians
  !  with some power of the laplace operator, and their gradients.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   pg    :center of gaussians
  !i   ph    :center of smooth Hankels
  !i   rsmg  :smoothing radius of gaussians
  !i   rsmh  :smoothing radius of Hankels
  !i   nlmg  :L-max for gaussians
  !i   nlmh  :L-max for Hankels
  !i   kmax  :cutoff in power of Laplace operator
  !i   ndim1 :leading dimensions of s,ds
  !i   ndim2 :second dimensions of s,ds
  !i   kdim  :dimensions s,ds
  !i   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
  !i   indxcg:index for Clebsch Gordon coefficients
  !i   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
  !i   cy    :Normalization constants for spherical harmonics
  !i   vol   :cell volume
  !i   slat  :struct containing information about the lattice
  !o Outputs
  !o   s     :integrals between Hankels and gaussians; see Remarks
  !o   ds    :gradient of s; see Remarks
  !r Remarks
  !r   s(L,M,k) contains integral of G_L^*(r-pg) (laplace)^k F_M(r-ph)
  !r            s,ds are generated for L=1..nlmg and M=1..nlmh
  !r   ds is gradient of s wrt pg; use -ds for grad wrt ph.
  !u Updates
  !u   28 Apr 97 adapted to handle case q=0 and e1 or e2 zero.
  !u   22 Apr 00 Adapted from nfp hhig_bl.f
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: jcg(*),indxcg(*),nlmg,nlmh,kmax,ndim1,ndim2,kdim
  double precision :: ph(3),pg(3),cg(1),cy(1),rsmg,rsmh !,slat(1)
  double complex s(ndim1,ndim2,0:kdim),ds(ndim1,ndim2,0:kdim,3)
  ! ... Local parameters
  integer :: nlm0,ktop0,m,lmaxh,ll,lmaxg,lmaxx,nlmx,nlmxp1,ktop,ktopp1, &
       ilm,kz,kx1,kx2,ky1,ky2,k,jlm,ilg,lg,ilh,lh,ii,indx,icg1,icg2, &
       icg,lm,ip,nlmxx
  double precision :: dr(3),gamh,gamg,rsmx,cz,cx1,cx2,cy1,cy2,fac
  ! FCPP#if F90
  complex(8),allocatable:: hkl(:,:),ghkl(:,:,:)
  ! FCPP#else
  ! FCPP      parameter (nlm0=196, ktop0=7)
  ! FCPP      double complex hkl(0:ktop0,nlm0),ghkl(0:ktop0,nlm0,3)
  ! FCPP#endif

  if (nlmh == 0 .OR. nlmg == 0) return
  gamh = 0.25d0*rsmh*rsmh
  gamg = 0.25d0*rsmg*rsmg
  rsmx = 2d0*dsqrt(gamg+gamh)
  do  1  m = 1, 3
     dr(m) = pg(m)-ph(m)
1 enddo
  lmaxh = ll(nlmh)
  lmaxg = ll(nlmg)
  lmaxx = lmaxg+lmaxh
  nlmx = (lmaxx+1)**2
  nlmxp1 = (lmaxx+2)**2
  ktop = max0(lmaxg,lmaxh) + kmax
  ktopp1 = ktop+1

  !     Memory allocation
  ! FCPP#if F90
  nlm0 = nlmxp1
  ktop0 = ktopp1
  allocate(hkl(0:ktop0,nlm0),ghkl(0:ktop0,nlm0,3))
  ! FCPP#endif
  if (nlmxp1 > nlm0)  call rxi('gfigbl: need nlm0 ge',nlmxp1)
  if (ktopp1 > ktop0) call rxi('gfigbl: need ktop0 ge',ktopp1)

  ! ... Make functions for connecting vector
  call fklbl(dr,rsmx,ktopp1,nlmxp1,ktop0,cy, hkl) !,slat

  ! ... Make gradients using CGs for p functions
  call dpzero(ghkl, 2*(ktop0+1)*nlm0*3)
  nlmxx = lmaxx*lmaxx
  do  ilm = 1, nlmx
     call scglp1(ilm,kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2)
     do  k = 0, ktop
        ghkl(k,ilm,1) = ghkl(k,ilm,1)-cx1*hkl(k,kx1)-cx2*hkl(k,kx2)
        ghkl(k,ilm,2) = ghkl(k,ilm,2)-cy1*hkl(k,ky1)-cy2*hkl(k,ky2)
        ghkl(k,ilm,3) = ghkl(k,ilm,3)-cz*hkl(k,kz)
        if (ilm <= nlmxx) then
           ghkl(k,kx1,1) = ghkl(k,kx1,1) - cx1*hkl(k+1,ilm)
           ghkl(k,kx2,1) = ghkl(k,kx2,1) - cx2*hkl(k+1,ilm)
           ghkl(k,ky1,2) = ghkl(k,ky1,2) - cy1*hkl(k+1,ilm)
           ghkl(k,ky2,2) = ghkl(k,ky2,2) - cy2*hkl(k+1,ilm)
           ghkl(k,kz,3)  = ghkl(k,kz,3)  - cz *hkl(k+1,ilm)
        endif
     enddo
  enddo

  do    k = 0, kmax
     do    jlm = 1, nlmh
        do    ilm = 1, nlmg
           s(ilm,jlm,k)    = dcmplx(0d0,0d0)
           ds(ilm,jlm,k,1) = dcmplx(0d0,0d0)
           ds(ilm,jlm,k,2) = dcmplx(0d0,0d0)
           ds(ilm,jlm,k,3) = dcmplx(0d0,0d0)
        enddo
     enddo
  enddo

  ! --- Combine with Clebsch-Gordan coefficients ---
  do  1111  ilg = 1, nlmg
     lg = ll(ilg)
     do  111  ilh = 1, nlmh
        lh = ll(ilh)
        ii = max0(ilg,ilh)
        indx = (ii*(ii-1))/2 + min0(ilg,ilh)
        icg1 = indxcg(indx)
        icg2 = indxcg(indx+1)-1
        do  11  icg = icg1, icg2
           ilm = jcg(icg)
           lm = ll(ilm)
           k = (lg+lh-lm)/2
           fac = (-1d0)**lg*cg(icg)
           do  12  ip = 0, kmax
              s(ilg,ilh,ip) = s(ilg,ilh,ip) + fac*hkl(k+ip,ilm)
              ds(ilg,ilh,ip,1) = ds(ilg,ilh,ip,1)+fac*ghkl(k+ip,ilm,1)
              ds(ilg,ilh,ip,2) = ds(ilg,ilh,ip,2)+fac*ghkl(k+ip,ilm,2)
              ds(ilg,ilh,ip,3) = ds(ilg,ilh,ip,3)+fac*ghkl(k+ip,ilm,3)
12         enddo
11      enddo
111  enddo
1111 enddo
  deallocate(hkl,ghkl)
end subroutine gfigbl

