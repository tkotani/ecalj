subroutine hhigbl(mode,p1,p2,q,rsm1,rsm2,e1,e2,nlm1,nlm2,kmax, &
     ndim1,ndim2,k0,cg,indxcg,jcg,cy,s,ds) !,slat
  use m_struc_def
  !- Integrals between smooth hankels with k-th power of Laplace operator
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :1's digit
  !i         :0 rsm1,rsm2,e1,e2 are scalars
  !i         :1 rsm1,rsm2,e1,e2 are l-dependent vectors
  !i         :10's digit
  !i         :1: do not calculate strux for any (l1,l2) pairs for
  !i         :   if rsm(l1) or rsm(l2) is zero.
  !i   p1    :first center
  !i   p2    :second center
  !i   q     :wave number for Bloch sum
  !i   rsm1  :smoothing radii of Hankels at p1 (l-dependent)
  !i   rsm2  :smoothing radii of Hankels at p2 (l-dependent)
  !i   e1    :energies of smooth Hankels at p1 (l-dependent)
  !i   e2    :energies of smooth Hankels at p2 (l-dependent)
  !i   nlm1  :L-cutoff for functions at p1
  !i   nlm2  :L-cutoff for functions at p2
  !i   kmax  :cutoff in power of Laplace operator
  !i   ndim1 :leading dimensions of s,ds
  !i   ndim2 :second dimensions of s,ds
  !i   k0    :dimensions s,ds
  !i   cg    :Clebsch Gordon coefficients (scg.f)
  !i   indxcg:index for Clebsch Gordon coefficients
  !i   jcg   :L q.n. for the C.G. coefficients (scg.f)
  !i   cy    :Normalization constants for spherical harmonics
  !i   slat  :struct containing information about the lattice
  !o Outputs
  !o   s     :integrals between smooth Bloch Hankels; see Remarks
  !o   ds    :gradient of s; see Remarks
  !r Remarks
  !r  s(L,M,k) contains integral of H_L^*(r-p1) (laplace)^k H_M(r-p2)
  !r  Gradient is wrt p1; use -ds for grad wrt p2.
  !u Updates
  !u   8 Jun 00  Added 10s digit mode.
  !u   23 May 00 Made rsm1,e1,rsm2,e2 l-dependent
  !u   22 Apr 00 Adapted from nfp hhig_bl.f
  !u   28 Apr 97 adapted to handle case q=0 and e1 or e2 zero.
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: mode,jcg(1),indxcg(1),nlm1,nlm2,kmax,ndim1,ndim2,k0
  real(8):: p1(3) , p2(3) , q(3) , cg(1) , cy(1) , rsm1(0:*) , &
       rsm2(0:*) , e1(0:*) , e2(0:*)
  !      type(s_lat)::slat

  double complex s(ndim1,ndim2,0:k0),ds(ndim1,ndim2,0:k0,3)
  ! ... Local parameters
  integer :: ilm,jlm,k,l1,l1t,l2,l2t,ll,lm11,lm12,lm21,lm22,lmx1,lmx2,m
  double precision :: dr(3)

  if (nlm1 == 0 .OR. nlm2 == 0) return

  do  1  m = 1, 3
     dr(m) = p1(m)-p2(m)
1 enddo
  lmx1 = ll(nlm1)
  lmx2 = ll(nlm2)

  do    k = 0, kmax
     do    jlm = 1, nlm2
        do    ilm = 1, nlm1
           s(ilm,jlm,k)    = dcmplx(0d0,0d0)
           ds(ilm,jlm,k,1) = dcmplx(0d0,0d0)
           ds(ilm,jlm,k,2) = dcmplx(0d0,0d0)
           ds(ilm,jlm,k,3) = dcmplx(0d0,0d0)
        enddo
     enddo
  enddo

  l1t = -1
  do  20  l1 = 0, lmx1
     if (l1 <= l1t) goto 20
     if (mod(mode,10) == 0) then
        l1t = lmx1
     else
        call gtbsl2(l1,lmx1,e1,rsm1,l1t)
     endif
     !       l1t = l1

     l2t = -1
     do  22  l2 = 0, lmx2
        if (l2 <= l2t) goto 22
        if (mod(mode,10) == 0) then
           l2t = lmx2
        else
           call gtbsl2(l2,lmx2,e2,rsm2,l2t)
        endif
        !         l2t = l2

        lm11 = l1**2+1
        lm12 = (l1t+1)**2
        lm21 = l2**2+1
        lm22 = (l2t+1)**2
        if (mode/10 == 1 .AND. rsm1(l1)*rsm2(l2) == 0) goto 22
        call phhigb(dr,q,rsm1(l1),rsm2(l2),e1(l1),e2(l2),lm11,lm12, &
             lm21,lm22,kmax,ndim1,ndim2,k0,cg,indxcg,jcg,cy,s,ds) !slat,
22   enddo
20 enddo
end subroutine hhigbl

subroutine phhigb(dr,q,rsm1,rsm2,e1,e2,mlm1,nlm1,mlm2,nlm2, &
     kmax,ndim1,ndim2,k0,cg,indxcg,jcg,cy,s,ds)!slat,
  use m_lattic,only:lat_vol
  use m_struc_def  !Cgetarg
  !- Integrals between smooth hankels with k-th power of Laplace operator
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   dr    :p1-p2
  !i   q     :wave number for Bloch sum
  !i   rsm1  :smoothing radius of Hankels at p1
  !i   rsm2  :smoothing radius of Hankels at p2
  !i   e1    :energy of smooth Hankels at p1
  !i   e2    :energy of smooth Hankels at p2
  !i   mlm1  :lower L-cutoff for functions at p1
  !i   nlm1  :upper L-cutoff for functions at p1
  !i   mlm2  :lower L-cutoff for functions at p2
  !i   nlm2  :upper L-cutoff for functions at p2
  !i   kmax  :cutoff in power of Laplace operator
  !i   ndim1 :leading dimensions of s,ds
  !i   ndim2 :second dimensions of s,ds
  !i   k0    :dimensions s,ds
  !i   cg    :Clebsch Gordon coefficients (scg.f)
  !i   indxcg:index for Clebsch Gordon coefficients
  !i   jcg   :L q.n. for the C.G. coefficients (scg.f)
  !i   cy    :Normalization constants for spherical harmonics
  !i   slat  :struct containing information about the lattice
  !o Outputs
  !o   s     :integrals between smooth Bloch Hankels; see Remarks
  !o   ds    :gradient of s; see Remarks
  !r Remarks
  !r  s(L,M,k) contains integral of H_L^*(r-p1) (laplace)^k H_M(r-p2)
  !r  Gradient is wrt p1; use -ds for grad wrt p2.
  !u Updates
  !u   22 Apr 00 Adapted from nfp hhig_bl.f
  !u   28 Apr 97 adapted to handle case q=0 and e1 or e2 zero.
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: jcg(1),indxcg(1),nlm1,nlm2,kmax,ndim1,ndim2,k0,mlm1,mlm2
  real(8):: q(3) , cg(1) , cy(1) , dr(3) , rsm1 , rsm2 , e1 , e2
  !      type(s_lat)::slat

  double complex s(ndim1,ndim2,0:k0),ds(ndim1,ndim2,0:k0,3)
  ! ... Local parameters
  integer :: nlm0,ktop0,lmax1,ll,lmax2,lmaxx,nlmx,nlmxp1,ktop,ktopp1, &
       k,ilm,kz,kx1,kx2,ky1,ky2,ilm1,l1,ilm2,l2,ii,indx,icg1,icg2, &
       icg,lm,ip
  double precision :: gam1,fpi,gam2,gamx,rsmx,qq,fac1,fac2,e,cz,cx1, &
       cx2,cy1,cy2,fac,add,vol,dgets
  parameter (nlm0=100, ktop0=10)
  double complex hkl1(0:ktop0,nlm0),hkl2(0:ktop0,nlm0), &
       ghkl(0:ktop0,nlm0,3),hsm(nlm0),hsmp(nlm0)

  vol = lat_vol

  fpi = 16d0*datan(1.d0)
  gam1 = 0.25d0*rsm1*rsm1
  gam2 = 0.25d0*rsm2*rsm2
  gamx = gam1+gam2
  rsmx = 2d0*dsqrt(gamx)
  qq = dsqrt(q(1)*q(1)+q(2)*q(2)+q(3)*q(3))
  lmax1 = ll(nlm1)
  lmax2 = ll(nlm2)
  lmaxx = lmax1+lmax2
  nlmx = (lmaxx+1)**2
  nlmxp1 = (lmaxx+2)**2
  ktop = max0(lmax1,lmax2)+kmax
  ktopp1 = ktop+1
  if (nlmxp1 > nlm0)  call rxi('hhigbl: need nlm0 ge',nlmxp1)
  if (ktopp1 > ktop0) call rxi('hhigbl: need ktop0 ge',ktopp1)

  ! --- Set up functions for connecting vector p2-p1 ---
  if (dabs(e1-e2) > 1d-5) then
     fac1 = dexp(gam2*(e2-e1))/(e1-e2)
     fac2 = dexp(gam1*(e1-e2))/(e2-e1)
     if (dabs(e1) > 1d-6) then
        call hklbl(dr,rsmx,e1,q,ktopp1,nlmxp1,ktop0,cy,hkl1) !slat,
     else
        if (qq > 1d-6) call rx('hhigbl: e1=0 only allowed if q=0')
        call fklbl(dr,rsmx,ktopp1,nlmxp1,ktop0,cy,hkl1) !slat,
     endif
     if (dabs(e2) > 1d-6) then
        call hklbl(dr,rsmx,e2,q,ktopp1,nlmxp1,ktop0,cy,hkl2) !slat,
     else
        if (qq > 1d-6) call rx('hhigbl: e2=0 only allowed if q=0')
        call fklbl(dr,rsmx,ktopp1,nlmxp1,ktop0,cy,hkl2) !slat,
     endif
     do  4  ilm = 1, nlmxp1
        do  5  k = 0, ktopp1
           hkl1(k,ilm) = fac1*hkl1(k,ilm) + fac2*hkl2(k,ilm)
5       enddo
4    enddo
  else
     e = .5d0*(e1+e2)
     if (qq < 1d-6 .AND. dabs(e) < 1d-6) &
          call rx('hhigbl: case q=0 and e1=e2=0 not available')
     call hklbl(dr,rsmx,e,q,ktopp1,nlmxp1,ktop0,cy,hkl2) !slat,
     call hsmbl(dr,rsmx,e,q,lmaxx+1,cy,hsm,hsmp) !slat,
     do  6  ilm = 1, nlmxp1
        hkl1(0,ilm) = hsmp(ilm) - gamx*hsm(ilm)
        do  7  k = 1, ktopp1
           hkl1(k,ilm) = -e*hkl1(k-1,ilm) - hkl2(k-1,ilm)
7       enddo
6    enddo
  endif

  ! ... Make gradients using CGs for p functions
  call dpzero(ghkl, 2*(ktop0+1)*nlm0*3)
  do  ilm = 1, nlmx
     call scglp1(ilm,kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2)
     do  k = 0, ktop
        ghkl(k,ilm,1) = ghkl(k,ilm,1)-cx1*hkl1(k,kx1)-cx2*hkl1(k,kx2)
        ghkl(k,ilm,2) = ghkl(k,ilm,2)-cy1*hkl1(k,ky1)-cy2*hkl1(k,ky2)
        ghkl(k,ilm,3) = ghkl(k,ilm,3)-cz*hkl1(k,kz)
        if (ilm <= lmaxx*lmaxx) then
           ghkl(k,kx1,1) = ghkl(k,kx1,1) - cx1*hkl1(k+1,ilm)
           ghkl(k,kx2,1) = ghkl(k,kx2,1) - cx2*hkl1(k+1,ilm)
           ghkl(k,ky1,2) = ghkl(k,ky1,2) - cy1*hkl1(k+1,ilm)
           ghkl(k,ky2,2) = ghkl(k,ky2,2) - cy2*hkl1(k+1,ilm)
           ghkl(k,kz,3)  = ghkl(k,kz,3)  - cz *hkl1(k+1,ilm)
        endif
     enddo
  enddo

  ! ... Combine with Clebsch-Gordan coefficients
  do  1111  ilm1 = mlm1, nlm1
     l1 = ll(ilm1)
     do  111  ilm2 = mlm2, nlm2
        l2 = ll(ilm2)
        ii = max0(ilm1,ilm2)
        indx = (ii*(ii-1))/2 + min0(ilm1,ilm2)
        icg1 = indxcg(indx)
        icg2 = indxcg(indx+1)-1
        do  11  icg = icg1, icg2
           ilm = jcg(icg)
           lm = ll(ilm)
           k = (l1+l2-lm)/2
           fac = fpi*(-1d0)**l1*cg(icg)
           do  12  ip = 0, kmax
              s(ilm1,ilm2,ip) = s(ilm1,ilm2,ip) + fac*hkl1(k+ip,ilm)
              ds(ilm1,ilm2,ip,1) = ds(ilm1,ilm2,ip,1) + fac*ghkl(k+ip,ilm,1)
              ds(ilm1,ilm2,ip,2) = ds(ilm1,ilm2,ip,2) + fac*ghkl(k+ip,ilm,2)
              ds(ilm1,ilm2,ip,3) = ds(ilm1,ilm2,ip,3) + fac*ghkl(k+ip,ilm,3)
12         enddo
11      enddo
111  enddo
1111 enddo

  ! ... Extra term when e1 or e2 is zero (but not both)
  if (mlm1 == 1 .AND. mlm2 == 1) then
     add = 0d0
     if (dabs(e2) < 1d-6) add = fpi*dexp(gam1*e1)/(vol*e1*e1)
     if (dabs(e1) < 1d-6) add = fpi*dexp(gam2*e2)/(vol*e2*e2)
     s(1,1,0) = s(1,1,0) + add
  endif

end subroutine phhigb

