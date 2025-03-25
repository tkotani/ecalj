!>Bloch sum of smooth Hankel and Gaussians. Expansion and Integrals.
module m_smhankel   !2023-4-28 memo: TK added qshortn(q). This allows q is not need to be within BZ.
  use m_ll,only:ll
  use m_factorial,only: factorial_init,factorial2,factorial
  use m_lmfinit,only: cg=>rv_a_ocg,indxcg=>iv_a_oidxcg,jcg=>iv_a_ojcg,cy=>rv_a_ocy,icgi,icge
  use m_lattic,only: vol=>lat_vol,plat=>lat_plat,qlat=>lat_qlat, awald=>lat_awald
  use m_lattic,only: dlv=>rv_a_odlv, nkd=>lat_nkd, qlv=>rv_a_oqlv,nkq=>lat_nkq
  ! JMP39: Bott, E., M. Methfessel, W. Krabs, and P. C. Schmidt.
  ! “Nonsingular Hankel Functions as a New Basis for Electronic Structure Calculations.”
  ! Journal of Mathematical Physics 39, no. 6 (June 1, 1998): 3393–3425.
  ! https://doi.org/doi:10.1063/1.532437.
  implicit none
  public hxpbl,hxpgbl,  hhibl,hhigbl,hhugbl, hgugbl, ggugbl, hxpos ! *bl means blochsum. *g* means gradient. xp means expansion. hhig means integral of h*h
         !xp=expansion !hh integral         !hg integral    !ggintegral   !expansion 
  private
contains
  subroutine hhugbl(mode,p1,p2,rsm1,rsm2,e1,e2,nlm1,nlm2,ndim1,ndim2, s,ds) ! Estatic energy integrals between Bloch Hankels, and gradients.
    !i   mode  :0 rsm1,rsm2,e1,e2 are scalars
    !i         :1 rsm1,rsm2,e1,e2 are l-dependent
    !i   p1    :first center
    !i   p2    :second center
    !i   rsm1  :smoothing radius of Hankels at p1 (l-dependent)
    !i   rsm2  :smoothing radius of Hankels at p2 (l-dependent)
    !i   e1    :energy  of Hankels at p1 (l-dependent)
    !i   e2    :energy  of Hankels at p2 (l-dependent)
    !i   nlm1  :L-max for  Hankels at p1
    !i   nlm2  :L-max for  Hankels at p2
    !i   ndim1 :leading dimensions of s,ds
    !i   ndim2 :second dimensions of s,ds
    !r Remarks
    !r   Gradient is wrt p1; use -ds for grad wrt p2.
    implicit none
    integer :: mode,nlm1,nlm2,ndim1,ndim2, kmax,kdim,i2,i1 ,l 
    real(8) :: add,q(3),gam1,gam2,xx1,xx2,zer(0:nlm1),bet1(nlm1),fac(nlm1), rsm1(0:*),rsm2(0:*),e1(0:*),e2(0:*),p1(3),p2(3)
    complex(8):: s(ndim1,ndim2),ds(ndim1,ndim2,3), wk(ndim1,ndim2),dwk(ndim1,ndim2,3)
    real(8),parameter:: pi = 4d0*datan(1d0), fpi = 4d0*pi
    q=0d0
    zer=0d0
    kmax = 0
    kdim = 0
    call hhigbl(mode,p1,p2,q,rsm1,rsm2,e1,e2, nlm1,nlm2,kmax,ndim1,ndim2,kdim,   s, ds)
    call hhigbl(mode,p1,p2,q,rsm1,rsm2,e1,zer,nlm1,nlm2,kmax,ndim1,ndim2,kdim,  wk,dwk)
    do  i2 = 1, nlm2
       l = merge(0,ll(i2),mode==0)
       s(:,i2)    = 8d0*pi/e2(l)*( s(:,i2)   - exp(e2(l)*rsm2(l)*rsm2(l)/4d0)* wk(:,i2)  )
       ds(:,i2,:) = 8d0*pi/e2(l)*(ds(:,i2,:) - exp(e2(l)*rsm2(l)*rsm2(l)/4d0)*dwk(:,i2,:))
    enddo
    gam1 = .25d0*rsm1(0)**2
    gam2 = .25d0*rsm2(0)**2
    s(1,1) = s(1,1) -2*fpi*dexp(gam1*e1(0))/(vol*e1(0))* fpi*dexp(gam2*e2(0))/(vol*e2(0)) *vol/e2(0) ! ... Extra term for l1=l2=0
  end subroutine hhugbl
  subroutine hgugbl(p1,p2,rsm1,rsm2,e1,nlm1,nlm2,ndim1,ndim2, s,ds)! Estatic energy integrals between Bloch Hankels and gaussians, and grad
    !i Inputs
    !i   p1    :first center
    !i   p2    :second center
    !i   rsm1  :smoothing radius of Hankels at p1
    !i   rsm2  :smoothing radius of gaussians at p2
    !i   e1    :energy  of Hankels at p1
    !i   nlm1  :L-max for  Hankels at p1
    !i   nlm2  :L-max for  gaussians at p2
    !i   ndim1 :leading dimensions of s,ds
    !i   ndim2 :second dimensions of s,ds
    !o Outputs
    !o   s     :integrals between Bloch Hankels and gaussians
    !o   ds    :gradient of s; see Remarks
    !r Remarks
    !r   Gradient is wrt p1; use -ds for grad wrt p2.
    implicit none
    integer :: nlm1,nlm2,ndim1,ndim2,l,kmax,kdim,ilm2,ilm1
    real(8):: rsm1,rsm2,p1(3),p2(3),e1
    complex(8):: s(ndim1,ndim2),ds(ndim1,ndim2,3)
    real(8) :: q(3),e2
    q=[0d0,0d0,0d0]
    kmax = 0
    e2 = 0d0
    call hhigbl(0,p1,p2,q,[(rsm1,l=0,ll(nlm1))],[(rsm2,l=0,ll(nlm2))],[(e1,l=0,ll(nlm1))],[(e2,l=0,ll(nlm2))],&
         nlm1,nlm2 ,kmax,ndim1,ndim2, k0=0,s=s,ds=ds)
    s(1:nlm1,1:nlm2)   = 2d0* s(1:nlm1,1:nlm2)
    ds(1:nlm1,1:nlm2,:)= 2d0*ds(1:nlm1,1:nlm2,:)
  end subroutine hgugbl
  subroutine hhigbl(mode,p1,p2,q,rsm1,rsm2,e1,e2,nlm1,nlm2,kmax,ndim1,ndim2,k0, s,ds)!Integrals between smooth hankels with k-th power of Laplace operator
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
    !o Outputs
    !o   s     :integrals between smooth Bloch Hankels; see Remarks
    !o   ds    :gradient of s; see Remarks
    !r Remarks
    !r  s(L,M,k) contains integral of H_L^*(r-p1) (laplace)^k H_M(r-p2)
    !r  Gradient is wrt p1; use -ds for grad wrt p2.
    implicit none
    integer :: mode,nlm1,nlm2,kmax,ndim1,ndim2,k0,ilm,jlm,k,l1,l1t,l2,l2t,lm11,lm12,lm21,lm22,lmx1,lmx2,m
    real(8):: p1(3),p2(3),q(3),rsm1(0:*) ,rsm2(0:*),e1(0:*),e2(0:*),dr(3)
    complex(8):: s(ndim1,ndim2,0:k0),ds(ndim1,ndim2,0:k0,3)
    if (nlm1 == 0 .OR. nlm2 == 0) return
    dr=p1-p2
    lmx1 = ll(nlm1)
    lmx2 = ll(nlm2)
    s=0d0
    ds=0d0
    l1t = -1
    do  l1 = 0, lmx1
       if (l1 <= l1t) cycle
       if (mod(mode,10) == 0) then
          l1t = lmx1
       else
          l1t = gtbsl2(l1,lmx1,e1,rsm1)
       endif
       l2t = -1
       do l2 = 0, lmx2
          if (l2 <= l2t) cycle
          if (mod(mode,10) == 0) then
             l2t = lmx2
          else
             l2t = gtbsl2(l2,lmx2,e2,rsm2)
          endif
          lm11 = l1**2+1
          lm12 = (l1t+1)**2
          lm21 = l2**2+1
          lm22 = (l2t+1)**2
          if (mode/10 == 1 .AND. rsm1(l1)*rsm2(l2) == 0) cycle
          call phhigb(dr,q,rsm1(l1),rsm2(l2),e1(l1),e2(l2),lm11,lm12, lm21,lm22,kmax,ndim1,ndim2,k0, s,ds) !cg,indxcg,jcg,cy,
       enddo
    enddo
  end subroutine hhigbl
  subroutine phhigb(dr,q,rsm1,rsm2,e1,e2,mlm1,nlm1,mlm2,nlm2,kmax,ndim1,ndim2,k0, s,ds)!Integrals between smooth hankels with k-th power of Laplace operator
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
    !o Outputs
    !o   s     :integrals between smooth Bloch Hankels; see Remarks
    !o   ds    :gradient of s; see Remarks
    !r Remarks
    !r  s(L,M,k) contains integral of H_L^*(r-p1) (laplace)^k H_M(r-p2)
    !r  Gradient is wrt p1; use -ds for grad wrt p2.
    implicit none
    integer :: nlm1,nlm2,kmax,ndim1,ndim2,k0,mlm1,mlm2 !jcg(1),indxcg(1),
    real(8):: q(3), dr(3),rsm1,rsm2,e1,e2 !cg(1),cy(1) 
    real(8) :: gam1,fpi,gam2,gamx,rsmx,qq,fac1,fac2,e,cz,cx1, cx2,cy1,cy2,fac,add
    complex(8):: s(ndim1,ndim2,0:k0),ds(ndim1,ndim2,0:k0,3)
    integer :: lmax1,lmax2,lmaxx,nlmx,nlmxp1,ktop,ktopp1, &
         k,ilm,kz,kx1,kx2,ky1,ky2,ilm1,l1,ilm2,l2,ii,indx,icg,ip
    fpi = 16d0*datan(1.d0)
    gam1 = .25d0*rsm1*rsm1
    gam2 = .25d0*rsm2*rsm2
    gamx = gam1+gam2
    rsmx = 2d0*dsqrt(gamx)
    qq = sum(q**2)**.5 
    lmax1 = ll(nlm1)
    lmax2 = ll(nlm2)
    lmaxx = lmax1+lmax2
    nlmx = (lmaxx+1)**2
    nlmxp1 = (lmaxx+2)**2
    ktop   = max0(lmax1,lmax2)+kmax !Note difference of kmax and ktop
    ktopp1 = ktop+1
    hklblock: block
      complex(8):: hkl1(0:ktopp1,nlmxp1),hkl2(0:ktopp1,nlmxp1), ghkl(0:ktopp1,nlmxp1,3),hsm(nlmxp1),hsmp(nlmxp1)
      if (dabs(e1-e2) > 1d-5) then !Set up functions for connecting vector p2-p1 ---
         if(qq > 1d-6.and.dabs(e1) <= 1d-6) call rx('hhigbl: e1=0 only allowed if q=0')
         if(qq > 1d-6.and.dabs(e2) <= 1d-6) call rx('hhigbl: e2=0 only allowed if q=0')
         if(dabs(e1) > 1d-6) then; call hklbl(dr,rsmx,e1,q,ktopp1,nlmxp1,ktopp1,hkl1) 
         else;                     call fklbl(dr,rsmx,ktopp1,nlmxp1,ktopp1,hkl1) ;   endif
         if(dabs(e2) > 1d-6) then; call hklbl(dr,rsmx,e2,q,ktopp1,nlmxp1,ktopp1,hkl2) 
         else;                     call fklbl(dr,rsmx,ktopp1,nlmxp1,ktopp1,hkl2) ;   endif
         hkl1(0:ktopp1,1:nlmxp1)=exp(gam2*(e2-e1))/(e1-e2)*hkl1(0:ktopp1,1:nlmxp1)+exp(gam1*(e1-e2))/(e2-e1)*hkl2(0:ktopp1,1:nlmxp1)
      else
         e = .5d0*(e1+e2)
         if(qq<1d-6 .AND. dabs(e) < 1d-6) call rx('hhigbl: case q=0 and e1=e2=0 not available')
         call hklbl(dr,rsmx,e,q,ktopp1,nlmxp1,ktopp1,hkl2) 
         call hsmbl(dr,rsmx,e,q,lmaxx+1,hsm,hsmp) 
         hkl1(0,:) = hsmp(:) - gamx*hsm(:)
         do k = 1, ktopp1
            hkl1(k,:) = -e*hkl1(k-1,:) - hkl2(k-1,:)
         enddo
      endif
      call ropylg2(lmaxx**2,ktop,nlmx,ktopp1,nlmxp1,hkl1, ghkl) !gradiend of hkl
      do  1111  ilm1 = mlm1, nlm1 ! ... Combine with Clebsch-Gordan coefficients
         l1 = ll(ilm1)
         do  111  ilm2 = mlm2, nlm2
            l2 = ll(ilm2)
            do  11  icg = icgi(ilm1,ilm2),icge(ilm1,ilm2) 
               ilm = jcg(icg)
               k = (l1+l2-ll(ilm))/2
               do ip = 0, kmax
                  s(ilm1,ilm2,ip)    = s(ilm1,ilm2,ip)    + fpi*(-1d0)**l1*cg(icg)*hkl1(k+ip,ilm)
                  ds(ilm1,ilm2,ip,:) = ds(ilm1,ilm2,ip,:) + fpi*(-1d0)**l1*cg(icg)*ghkl(k+ip,ilm,:)
               enddo
11          enddo
111      enddo
1111  enddo
    endblock hklblock
    if (mlm1 == 1 .AND. mlm2 == 1) then !Extra term when e1 or e2 is zero (but not both)
       add = 0d0
       if (dabs(e2) < 1d-6) add = fpi*dexp(gam1*e1)/(vol*e1*e1)
       if (dabs(e1) < 1d-6) add = fpi*dexp(gam2*e2)/(vol*e2*e2)
       s(1,1,0) = s(1,1,0) + add
    endif
  end subroutine phhigb
  subroutine hhibl(p1,p2,q,rsm1,rsm2,e1,e2,nlm1,nlm2,kmax, ndim1,ndim2, s) ! Integrals between smooth Hankels with k-th power of Laplace operator.
    !   do not calculate strux for any (l1,l2) pairs if rsm(l1) or rsm(l2) is zero.
    !i Inputs
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
    !i   ndim1 :leading dimension of s
    !i   ndim2 :second dimension of s
    !i   cg    :Clebsch Gordon coefficients (scg.f)
    !i   indxcg:index for Clebsch Gordon coefficients
    !i   jcg   :L q.n. for the C.G. coefficients (scg.f)
    !i   cy    :Normalization constants for spherical harmonics
    !o Outputs
    !o   s     :integrals; see Remarks
    !r Remarks
    !r   s(L,M) contains integral of H_L^*(r-p1) (laplace)^k H_M(r-p2)
    !r   Row L corresponds to p1 and col M corresponds to p2.
    !r   Strux s(L,M) is computed for dr=p1-p2
    !r   See JMP 39, 3393, Section 10
    !u Updates
    !u   19 May 00 Made rsm1,e1,rsm2,e2 l-dependent.  Elements for which
    !u             rsm1 =0 or rsm2 = 0 are not computed.
    implicit none
    integer :: nlm1,nlm2,kmax,ndim1,ndim2 !jcg(1),indxcg(1),
    real(8) :: p1(3),p2(3),q(3),rsm1(0:*),rsm2(0:*),e1(0:*),e2(0:*) !cg(1),cy(1),
    complex(8):: s(ndim1,ndim2,0:kmax)
    integer :: m,lmx1,lmx2,l1,l2,k,jlm,ilm,lm11,lm21,lm12,lm22,l1t,l2t
    real(8) :: dr(3)
    if (nlm1 == 0 .OR. nlm2 == 0) return
    dr = p1-p2
    lmx1 = ll(nlm1)
    lmx2 = ll(nlm2)
    s=0d0
    l1t = -1
    do 20  l1 = 0, lmx1
       if (l1 <= l1t) cycle
       l1t = gtbsl2(l1,lmx1,e1,rsm1)
       l2t = -1
       do l2 = 0, lmx2
          if (l2 <= l2t) cycle
          l2t= gtbsl2(l2,lmx2,e2,rsm2)
          lm11 = l1**2+1
          lm12 = (l1t+1)**2
          lm21 = l2**2+1
          lm22 = (l2t+1)**2
          if (rsm1(l1)*rsm2(l2) == 0) cycle
          call phhibl(dr,q,rsm1(l1),rsm2(l2),e1(l1),e2(l2),lm11,lm12,lm21,lm22,kmax,ndim1,ndim2,s) !,cg,indxcg,jcg,cy
       enddo
20  enddo
  end subroutine hhibl
  subroutine phhibl(dr,q,rsm1,rsm2,e1,e2,mlm1,nlm1,mlm2,nlm2,kmax, ndim1,ndim2, s)! Integrals between smooth Hankels with k-th power of Laplace operator.
    !i Inputs
    !i   dr    :p1-p2
    !i   q     :wave number for Bloch sum
    !i   rsm1  :smoothing radius of Hankels at p1
    !i   rsm2  :smoothing radius of Hankels at p2
    !i   e1    :energy of smooth Hankels at p1
    !i   e2    :energy of smooth Hankels at p2
    !i   nlm1  :L-cutoff for functions at p1
    !i   nlm2  :L-cutoff for functions at p2
    !i   kmax  :cutoff in power of Laplace operator
    !i   ndim1 :leading dimension of s
    !i   ndim2 :second dimension of s
    !i   cg    :Clebsch Gordon coefficients (scg.f)
    !i   indxcg:index for Clebsch Gordon coefficients
    !i   jcg   :L q.n. for the C.G. coefficients (scg.f)
    !i   cy    :Normalization constants for spherical harmonics
    !o Outputs
    !o   s     :integrals; see Remarks
    !r Remarks
    !r   s(L,M) contains integral of H_L^*(r-p1) (laplace)^k H_M(r-p2)
    implicit none
    integer :: mlm1,nlm1,mlm2,nlm2,kmax,ndim1,ndim2 !jcg(1),indxcg(1),
    real(8) :: dr(3),q(3),rsm1, rsm2,e1,e2 !,cg(1),cy(1)
    complex(8):: s(ndim1,ndim2,0:kmax)
    integer :: icg,ii,ilm,ilm1,ilm2,indx,ip,k, ktop,l1,l2,lm,lmax1,lmax2,lmaxx,nlmx
    real(8) :: fpi,e,fac,fac1,fac2,gam1,gam2,gamx,rsmx
    fpi = 16d0*datan(1.d0)
    gam1 = 0.25d0*rsm1*rsm1
    gam2 = 0.25d0*rsm2*rsm2
    gamx = gam1+gam2
    rsmx = 2d0*dsqrt(gamx)
    lmax1 = ll(nlm1)
    lmax2 = ll(nlm2)
    lmaxx = lmax1+lmax2
    nlmx = (lmaxx+1)**2
    ktop = max0(lmax1,lmax2)+kmax
    block
      complex(8):: hkl1(0:ktop,nlmx),hkl2(0:ktop,nlmx), hsm(nlmx),hsmp(nlmx)
      if (dabs(e1-e2) > 1d-5) then !Set up functions for connecting vector p1-p2
         call hklbl(dr,rsmx,e1,q,ktop,nlmx,ktop, hkl1) 
         call hklbl(dr,rsmx,e2,q,ktop,nlmx,ktop, hkl2) 
         hkl1(0:ktop,1:nlmx) = dexp(gam2*(e2-e1))/(e1-e2) *hkl1(0:ktop,1:nlmx) + dexp(gam1*(e1-e2))/(e2-e1) *hkl2(0:ktop,1:nlmx)
      else
         e = .5d0*(e1+e2)
         call hklbl(dr,rsmx,e,q,ktop,nlmx,ktop, hkl2) 
         call hsmbl(dr,rsmx,e,q,lmaxx, hsm,hsmp) 
         do  ilm = 1, nlmx
            hkl1(0,ilm) = hsmp(ilm) - gamx*hsm(ilm)
            do   k = 1, ktop
               hkl1(k,ilm) = -e*hkl1(k-1,ilm) - hkl2(k-1,ilm)
            enddo
         enddo
      endif
      do  1111  ilm1 = mlm1, nlm1 !Combine with Clebsch-Gordan coefficients
         l1 = ll(ilm1)
         do  111  ilm2 = mlm2, nlm2
            l2 = ll(ilm2)
            do icg = icgi(ilm1,ilm2),icge(ilm1,ilm2)
               ilm = jcg(icg)
               s(ilm1,ilm2,0:kmax) = s(ilm1,ilm2,0:kmax) + fpi*(-1d0)**l1*cg(icg)* [(hkl1((l1+l2-ll(ilm))/2+ip,ilm),ip=0,kmax)]
            enddo
111      enddo
1111  enddo
    endblock
  end subroutine phhibl
  subroutine gklbl(p,rsm,e,q,kmax,nlm,k0, gkl) ! Bloch-sums of k,L-dependent gaussians
    use m_lmfinit,only:alat=> lat_alat
    use m_shortn3_plat,only: shortn3_plat,nout,nlatout
    !i Inputs
    !i   p     :Function is centered at p
    !i   rsm   :smoothing radius
    !i   e     :gkL scaled by exp(e*rsm**2/4)
    !i   q     :wave number for Bloch sum (units of 2*pi/alat)
    !i   kmax  :polynomial cutoff
    !i   nlm   :L-cutoff for gkl
    !i   k0    :leading dimension of gkl
    !i   cy    :Normalization constants for spherical harmonics
    !o Outputs
    !o   gkl   :Bloch-summed Gaussians
    implicit none
    integer :: k0,kmax,nlm
    real(8):: e,rsm,q(3),p(3) ,pi,sp,p1(3),rwald,ppin(3) !, cy(1)
    complex(8):: gkl(0:k0,nlm),cfac,phase,img=(0d0,1d0)
    integer:: ilm,k,lmax,owk,l,m
    if (nlm == 0) return
    pi = 4d0*datan(1d0)
    ppin=matmul(transpose(qlat),p) 
    call shortn3_plat(ppin) 
    p1= matmul(plat,ppin+nlatout(:,1))
    phase = exp(img*2*pi*sum(q*(p-p1)))
    lmax = ll(nlm)
    rwald = 1d0/awald
    ! ... If the smoothing radius is larger than the Ewald parameter, make the summation in reciprocal space, else in real space
    if (rsm >= rwald) then
       call gklblq ( p1,rsm,q,kmax,nlm,k0,alat, vol,gkl )
    else
       call gklbld ( p1,rsm,q,kmax,nlm,k0,alat, gkl )
    endif
    cfac = exp(0.25d0*e*rsm*rsm)*phase
    do ilm=1,nlm
       l=ll(ilm)
       gkl(0:kmax,ilm) = cfac*cy(ilm)*(-img)**l*gkl(0:kmax,ilm)
    enddo
  end subroutine gklbl
  subroutine gklbld(p,rsm,q,kmax,nlm,k0,alat,gkl) ! Evaluate gkl in real space
    !i   p     :Function is centered at p
    !i   rsm   :smoothing radius
    !i   q     :wave number for Bloch sum
    !i   kmax  :polynomial cutoff
    !i   nlm   :L-cutoff for gkl
    !i   k0    :leading dimension of gkl
    !i   alat  :length scale of lattice and basis vectors, a.u.
    !i   dlv   :direct-space lattice vectors
    !i   nkd   :number of dlv
    !    wk    :work array of dimension (kmax+1)(lmax+1), lmax=ll(nlm)
    !o Outputs
    !o   gkl   :Bloch-summed Gaussians
    implicit none
    integer :: k0,kmax,nlm !nkd,
    real(8) :: yl(nlm),r(3),qdotr,r1,r2, p(3),q(3),alat,rsm,v(3)
    complex(8):: gkl(0:k0,nlm),cfac,img=(0d0,1d0)
    integer :: ilm,ir,k,l,lmax,m,nm
    real(8),parameter:: tpi = 8d0*datan(1d0)
    real(8),allocatable:: wk(:,:)
    lmax = ll(nlm)
    allocate(wk(0:kmax,0:lmax))
    gkl = 0d0
    do  ir = 1, nkd
       v=alat*(p(:)-dlv(:,ir))
       call sylm(v, yl,lmax,r2)
       call radgkl(r2**.5d0,rsm,kmax,lmax,kmax,wk)
       cfac = exp(img*tpi*sum(q(:)*dlv(:,ir)))
       do ilm=1,nlm
          l=ll(ilm)
          gkl(0:kmax,ilm) = gkl(0:kmax,ilm) + yl(ilm)*cfac*wk(0:kmax,l)*img**l
       enddo
    enddo
  end subroutine gklbld
  subroutine gklblq(p,rsm,q,kmax,nlm,k0,alat,vol, gkl) ! Evaluate gkl in reciprocal space
    use m_qplist,only:qshortn
    !i Inputs
    !i   p     :Function is centered at p
    !i   rsm   :smoothing radius
    !i   q     :wave number for Bloch sum
    !i   kmax  :polynomial cutoff
    !i   nlm   :L-cutoff for gkl
    !i   k0    :leading dimension of gkl
    !i   alat  :length scale of lattice and basis vectors, a.u.
    !i   nkq   :number of qlv
    !i   vol   :cell volume
    !o Outputs
    !o   gkl   :Bloch-summed Gaussians
    implicit none
    integer :: k0,kmax,nlm !,nkq
    real(8) :: alat,rsm,vol,q(3),p(3) 
    complex(8):: gkl(0:k0,nlm)
    integer :: ilm,ir,k,lmax
    real(8) :: a,gamma,r2,scalp,tpiba,vfac,r(3),yl(nlm)
    complex(8):: eiphi,add,add0
    complex(8):: img=(0d0,1d0)
    real(8),parameter:: tpi = 8d0*datan(1d0)
    real(8) :: qshortn_q(3)
    a = 1d0/rsm
    gamma = .25d0/(a*a)
    vfac = 1d0/vol
    tpiba = tpi/alat
    lmax = ll(nlm)
    gkl=0d0
    !MO moved qshortn outside of ir loop 2024-11-08
    qshortn_q = qshortn(q)
    do ir = 1,nkq
       ! r = tpiba*(qshortn(q)+qlv(:,ir)) !r = tpiba*(q(:)+qlv(:,ir))
       r = tpiba*(qshortn_q+qlv(:,ir)) !r = tpiba*(q(:)+qlv(:,ir))
       eiphi = exp(img*alat*sum(r*p))
       call sylm(r,yl,lmax,r2)
       do ilm = 1, nlm
          gkl(0:kmax,ilm)=gkl(0:kmax,ilm)+ eiphi*vfac*exp(-gamma*r2)*[((-r2)**k,k=0,kmax)]*yl(ilm)
       enddo
    enddo
  end subroutine gklblq
  subroutine hklbl(p,rsm,e,q,kmax,nlm,k0, hkl) !Bloch-sums of k,L-dependent smooth Hankel functions.
    use m_lmfinit,only: alat=>lat_alat
    use m_shortn3_plat,only: shortn3_plat,nout,nlatout
    use m_hsmq,only: hsmq
    use m_qplist,only:qshortn
    use m_ftox
    !i Inputs
    !i   p     :Function is centered at p
    !i   rsm   :smoothing radius
    !i   e     :energy of smoothed Hankel
    !i   q     :wave number for Bloch sum
    !i   kmax  :polynomial cutoff
    !i   nlm   :L-cutoff for hkl
    !i   k0    :leading dimension of hkl
    !i   cy    :Normalization constants for spherical harmonics
    !o Outputs
    !o   hkl   :Bloch-summed smoothed Hankels
    !r Remarks
    !r   H_kL = laplace^k H_L
    !r   Uses the recursion relation H_k+1,L = -e*H_kL - 4*pi*G_kL
    implicit none
    integer :: k0,kmax,nlm
    real(8):: e,rsm,q(3),p(3) 
    complex(8):: hkl(0:k0,nlm),gkl(0:k0,nlm), hsm(nlm),hsmp(nlm),phase
    integer:: ilm,job,k,lmax,nrx, owk,oyl
    real(8),parameter:: pi = 4d0*datan(1d0),fpi = 4d0*pi
    real(8) :: sp,p1(3),ppin(3)
    real(8),allocatable:: wk(:),yl(:)
    complex(8):: img=(0d0,1d0)
    if (nlm == 0) return
    lmax = ll(nlm)
    ppin=matmul(transpose(qlat),p) 
    call shortn3_plat(ppin) 
    p1= matmul(plat,ppin+nlatout(:,1)) !ppin (fractional) is shortened to be  p1
    sp = 2*pi*sum(q*(p-p1))
    phase = exp(img*sp) 
    nrx = max(nkd,nkq)
    allocate(yl(nrx*(lmax+1)**2))
    call hsmq ( 1,0,[ll(nlm)],[e],[rsm],0000,qshortn(q),p1,nrx,nlm,yl, awald,alat,qlv,nkq,dlv,nkd,vol, hsm,hsmp )
    if (rsm > 1d0/awald) then
       call gklbl(p1,rsm,e,q,kmax-1,nlm,k0, gkl) 
    else
       call gklq(lmax,rsm,q,p1,e, kmax - 1,k0,alat,nrx,yl,gkl )
    endif
    deallocate(yl)
    ! --- H_kL by recursion ---
    do    ilm = 1, nlm
       hkl(0,ilm) = hsm(ilm)
       do k = 1, kmax
          hkl(k,ilm) = -e*hkl(k-1,ilm) - fpi*gkl(k-1,ilm)
       enddo
    enddo
    if (sp /= 0) hkl(0:kmax,:) = phase*hkl(0:kmax,:)! ... Put in phase to undo shortening
!    write(6,ftox)'sssssssumhkl=',sum(hkl)
  end subroutine hklbl
  subroutine fklbl(p,rsm,kmax,nlm,k0, fkl) !Bloch sum of smooth Hankels for e=0 and q=(0,0,0).
    use m_lmfinit,only: alat=>lat_alat,tol=>lat_tol
    use m_shortn3_plat,only: shortn3_plat,nout,nlatout
    use m_hsmq,only: hsmqe0
    use m_qplist,only:qshortn
    !i Inputs
    !i   p     :Function is centered at p
    !i   rsm   :smoothing radius
    !i   kmax  :polynomial cutoff
    !i   nlm   :L-cutoff for gkl
    !i   k0    :leading dimension of gkl
    !i   cy    :Normalization constants for spherical harmonics
    !o Outputs
    !o   fkl   :Bloch-summed Hankels for q=0 and e=0
    !r Remarks
    !r   For (k=0,l=0) f equals the limit of hklbl minus the avg value.
    !r   For all other cases f is the limit of hklbl as e goes to zero.
    implicit none
    integer :: kmax,nlm,k0
    real(8):: p(3),rsm,ppin(3)
    integer::lmax,nrx,owk,oyl,job, ilm,k
    real(8) :: q(3),p1(3),fpi,y0,e
    complex(8):: fsm(nlm),fkl(0:k0,nlm),gkl(0:k0,nlm)
    real(8),allocatable:: wk(:),yl(:)
    if (nlm == 0) return
    fpi = 16d0*datan(1d0)
    y0 = 1d0/dsqrt(fpi)
    lmax = ll(nlm)
    e = 0d0
    q = 0d0
    ppin=matmul(transpose(qlat),p) 
    call shortn3_plat(ppin) 
    p1= matmul(plat,ppin+nlatout(:,1)) !p1 is shortened p for qlat modulo
    nrx = max(nkd,nkq)
    allocate( yl(nrx*(lmax+1)**2))
    call hsmqe0 ( lmax,rsm,0,qshortn(q),p1,nrx,nlm,yl, awald,alat,qlv,nkq,dlv,nkd,vol, fsm  )
    if (rsm > 1d0/awald) then; call gklbl(p1,rsm,e,q,kmax-1,nlm,k0, gkl) 
    else;                      call gklq(lmax,rsm,q,p1,e,kmax-1,k0,alat, nrx,yl,gkl )
    endif
    deallocate(yl)
    do ilm = 1, nlm ! ... Upward recursion in k: mainly sets fkl = -4*pi * g(k-1,l)
       fkl(0,ilm) = fsm(ilm)
       do    k = 1, kmax
          fkl(k,ilm) = -fpi*gkl(k-1,ilm) 
       enddo
    enddo
    fkl(1,1) = fkl(1,1) + fpi*y0/vol !! ... Add extra term to F(k=1,l=0)
  end subroutine fklbl
  subroutine gfigbl(pg,ph,rsmg,rsmh,nlmg,nlmh,kmax,ndim1,ndim2,kdim, s,ds) ! Integrals between smooth hankels(eh=0,q=0) and gaussians 
    !  with some power of the laplace operator, and their gradients.
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
    !o Outputs
    !o   s     :integrals between Hankels and gaussians; see Remarks
    !o   ds    :gradient of s; see Remarks
    !r Remarks
    !r   s(L,M,k) contains integral of G_L^*(r-pg) (laplace)^k F_M(r-ph)
    !r            s,ds are generated for L=1..nlmg and M=1..nlmh
    !r   ds is gradient of s wrt pg; use -ds for grad wrt ph.
    implicit none
    integer :: nlmg,nlmh,kmax,ndim1,ndim2,kdim !jcg(*),indxcg(*),
    real(8) :: ph(3),pg(3),rsmg,rsmh  !,cg(1),cy(1)
    complex(8):: s(ndim1,ndim2,0:kdim),ds(ndim1,ndim2,0:kdim,3)
    integer ::m,lmaxh,lmaxg,lmaxx,nlmx,nlmxp1,ktop,ktopp1, &
         ilm,kz,kx1,kx2,ky1,ky2,k,jlm,ilg,lg,ilh,lh,ii,indx, icg,lm,ip,nlmxx
    real(8) :: dr(3),gamh,gamg,rsmx,cz,cx1,cx2,cy1,cy2,fac
    complex(8),allocatable:: hkl(:,:),ghkl(:,:,:)
    if (nlmh == 0 .OR. nlmg == 0) return
    gamh = 0.25d0*rsmh*rsmh
    gamg = 0.25d0*rsmg*rsmg
    rsmx = 2d0*dsqrt(gamg+gamh)
    dr=pg-ph
    lmaxh = ll(nlmh)
    lmaxg = ll(nlmg)
    lmaxx = lmaxg+lmaxh
    nlmx = (lmaxx+1)**2
    nlmxp1 = (lmaxx+2)**2
    ktop = max0(lmaxg,lmaxh) + kmax
    ktopp1 = ktop+1
    allocate(hkl(0:ktopp1,nlmxp1),ghkl(0:ktopp1,nlmxp1,3))
    call fklbl(dr,rsmx,ktopp1,nlmxp1,ktopp1, hkl)  ! hkl for connecting vector dr
    call ropylg2(lmaxx*lmaxx,ktop,nlmx,ktopp1,nlmxp1,hkl, ghkl) !gradiend of hkl
    s=0d0
    ds=0d0
    do  1111  ilg = 1, nlmg !Combine with Clebsch-Gordan coefficients ---
       lg = ll(ilg)
       do  111  ilh = 1, nlmh
          lh = ll(ilh)
          do 11  icg = icgi(ilg,ilh),icge(ilg,ilh)
             ilm = jcg(icg)
             lm = ll(ilm)
             k = (lg+lh-lm)/2
             fac = (-1d0)**lg*cg(icg)
             do  12  ip = 0, kmax
                s(ilg,ilh,ip)    = s(ilg,ilh,ip)   + fac*hkl(k+ip,ilm)
                ds(ilg,ilh,ip,:) = ds(ilg,ilh,ip,:)+ fac*ghkl(k+ip,ilm,:)
12           enddo
11        enddo
111    enddo
1111 enddo
    deallocate(hkl,ghkl)
  end subroutine gfigbl
  subroutine ggugbl(p1,p2,rsm1,rsm2,nlm1,nlm2,ndim1,ndim2, s,ds) ! Estatic energy integrals between Bloch gaussians, and gradients.
    !i Inputs
    !i   p1    :first center
    !i   p2    :second center
    !i   rsm1  :smoothing radius of Gaussians at p1
    !i   rsm2  :smoothing radius of Gaussians at p2
    !i   e1    :energy  of Gaussians at p1
    !i   e2    :energy  of Gaussians at p2
    !i   nlm1  :L-max for  Gaussians at p1
    !i   nlm2  :L-max for  Gaussians at p2
    !i   ndim1 :leading dimensions of s,ds
    !i   ndim2 :second dimensions of s,ds
    !o Outputs
    !o   s     :integrals between Bloch Gaussians
    !o   ds    :gradient of s; see Remarks
    !r Remarks
    !r   Gradient is wrt p1; use -ds for grad wrt p2.
    implicit none
    integer :: nlm1,nlm2,ndim1,ndim2
    real(8):: rsm1,rsm2,p1(3),p2(3)
    complex(8):: s(ndim1,ndim2),ds(ndim1,ndim2,3)
    integer:: kmax,kdim,ilm2,ilm1
    kmax = 0
    kdim = 0
    call gfigbl(p1,p2,rsm1,rsm2,nlm1,nlm2,kmax,ndim1,ndim2,kdim, s,ds)  !See Eq.(11.48)
    do  ilm2 = 1, nlm2
       do  ilm1 = 1, nlm1
          s(ilm1,ilm2)    = 2d0*s(ilm1,ilm2)
          ds(ilm1,ilm2,:) = 2d0*ds(ilm1,ilm2,:)
       enddo
    enddo
  end subroutine ggugbl
  subroutine hklgbl(p,rsm,e,q,kmax,nlm,k0,nlm0,hkl,ghkl) ! Bloch-sums of k,L-dependent smooth hankel functions and gradients
    !i Inputs
    !i   p     :Function is centered at p
    !i   rsm   :smoothing radius
    !i   e     :energy of smoothed Hankel
    !i   q     :wave number for Bloch sum
    !i   kmax  :polynomial cutoff
    !i   nlm   :L-cutoff for hkl
    !i   k0    :leading dimension of hkl
    !i   cy    :Normalization constants for spherical harmonics
    !i   slat  :struct containing information about the lattice
    !o Outputs
    !o   hkl   :Bloch-summed smoothed Hankels
    !o   ghkl  :gradients of hkl
    !r Remarks
    !r   H_kL = laplace^k H_L
    !r   Uses the recursion relation H_k+1,L = -e*H_kL - 4*pi*G_kL
    !r   H_kL are made to kmax+1, lmax+1 in order to assemble grads.
    implicit none
    integer :: k0,kmax,nlm,nlm0,lmax,nlm1
    real(8) :: q(3),p(3),e,rsm
    complex(8):: hkl(0:k0,nlm0),ghkl(0:k0,nlm0,3)
    if (nlm == 0) return
    lmax = ll(nlm)
    nlm1 = (lmax+2)**2
    if (nlm1 > nlm0) call rxi('hklgbl: need nlm0 ge ',nlm1)
    if (kmax+1 > k0) call rxi('hklgbl: need k0 ge ',kmax+1)
    call hklbl(p,rsm,e,q,kmax+1,nlm1,k0, hkl) ! ... Make Hkl's up to one higher in l and k
    call ropylg2(lmax**2,kmax,nlm,k0,nlm0,hkl, ghkl) !gradiend of hkl
  end subroutine hklgbl
  subroutine ghibl(pg,ph,q,rsmg,rsmh,eg,eh,nlmg,nlmh,kmax,ndim1,ndim2, s) ! Block of integrals between smooth hankels and gaussians with some power of the laplace operator.
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
    implicit none
    integer :: nlmg,nlmh,kmax,ndim1,ndim2 !,jcg(1),indxcg(1)
    real(8) :: rsmg,rsmh(*),eg,eh(*), ph(3),pg(3),q(3) !,cg(1),cy(1)
    complex(8):: s(ndim1,ndim2,0:kmax)
    integer :: icg,ii,ilg,ilh,ilm,indx,ip,jlm,k,ktop,lg,lh,lm,lmaxg,lmaxh,lmaxx,m,nlmx,l1,l2,ilm1,ilm2
    complex(8),allocatable:: hkl(:,:)
    real(8) :: ee,fac,gamg,gamh,rsmx,dr(3),e,rsm
    if (nlmh == 0 .OR. nlmg == 0) return
    dr=pg-ph
    lmaxh = ll(nlmh) ! ... rsmh- and eh- independent setup
    lmaxg = ll(nlmg)
    lmaxx = lmaxg+lmaxh
    nlmx = (lmaxx+1)**2
    ktop = max0(lmaxg,lmaxh)+kmax
    allocate( hkl(0:ktop,nlmx))
    s(1:nlmg,1:nlmh,0:kmax) = 0d0
    l2 = -1
    do  20  l1 = 0, lmaxh ! Loop over sequences of l with a common rsm,e ---
       if (l1 <= l2) cycle
       l2=gtbsl2(l1,lmaxh,eh,rsmh)
       rsm  = rsmh(l1+1)
       e    = eh(l1+1)
       if (rsm <= 0 .OR. e > 0) cycle
       ilm1 = l1**2+1
       ilm2 = (l2+1)**2
       lmaxx= lmaxg+l2
       nlmx = (lmaxx+1)**2
       gamh = 0.25d0*rsm*rsm
       gamg = 0.25d0*rsmg*rsmg
       rsmx = 2d0*dsqrt(gamg+gamh)
       ktop = max0(lmaxg,l2)+kmax
       call hklbl(dr,rsmx,e,q,ktop,nlmx,ktop, hkl)
       ee = dexp(gamg*(eg-e)) ! Combine with Clebsch-Gordan coefficients
       do  1111   ilg = 1, nlmg
          lg = ll(ilg)
          do 111  ilh = ilm1, ilm2
             lh = ll(ilh)
             do icg = icgi(ilg,ilh),icge(ilg,ilh)
                ilm = jcg(icg)
                lm = ll(ilm)
                k = (lg+lh-lm)/2
                s(ilg,ilh,:) = s(ilg,ilh,:) + ee*(-1d0)**lg*cg(icg)*hkl(k:k+kmax,ilm)
             enddo
111       enddo
1111   enddo
20  enddo
  end subroutine ghibl
  subroutine ghigbl(pg,ph,q,rsmg,rsmh,eg,eh,nlmg,nlmh,kmax,ndim1,ndim2,k0, s,ds) ! Block of integrals between smooth hankels and gaussians with some power of the laplace operator, and gradients.
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
    !o Outputs
    !o   s     :integrals of gaussian and Hankels; see Remarks
    !o   ds    :gradients of s; see Remarks
    !r Remarks
    !r   s(L,M,k) contains integral of G_L^*(r-pg) (lap)^k H_M(r-ph)
    !r  ds((L,M,k) contains grad s wrt ph; take negative for grad wrt pg.
    !u Updates
    !u   25 May 00 Made rsmh,eh l-dependent
    implicit none
    integer :: k0,kmax,ndim1,ndim2,nlmg,nlmh,icg,ii,ilg,ilh,ilm,ilm1,ilm2,indx,ip,jlm,k,ktop,ktopp1,l1,l2,lg,lh,lm,&
         lmaxg,lmaxh,lmaxx,m,nlm1,nlmx
    real(8) :: rsmg,rsmh(*),eg,eh(*),ph(3),pg(3),q(3),ee,fac,gamg,gamh,rsmx,dr(3),e,rsm
    complex(8):: s(ndim1,ndim2,0:k0),ds(ndim1,ndim2,0:k0,3)
    complex(8),allocatable:: hkl(:,:),dhkl(:,:,:)
    if (nlmh == 0 .OR. nlmg == 0) return
    ! ... rsmh- and eh- independent setup
    dr= pg-ph
    lmaxh = ll(nlmh)
    lmaxg = ll(nlmg)
    lmaxx = lmaxg+lmaxh
    nlmx = (lmaxx+1)**2
    nlm1 = (lmaxx+2)**2
    ktop   = max0(lmaxg,lmaxh)+kmax
    ktopp1 = ktop+1
    allocate( hkl(0:ktopp1,nlm1),dhkl(0:ktopp1,nlm1,3))
    s=0d0
    ds=0d0
    l2 = -1
    do  20  l1 = 0, lmaxh ! --- Loop over sequences of l with a common rsm,e ---
       if (l1 <= l2) goto 20
       l2=gtbsl2(l1,lmaxh,eh,rsmh)
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
       call hklgbl(dr,rsmx,e,q, ktop,nlmx,ktopp1,nlm1, hkl,dhkl)
       ee = dexp(gamg*(eg-e))
       do  1111  ilg = 1, nlmg !   ... Combine with Clebsch-Gordan coefficients
          lg = ll(ilg)
          do  111  ilh = ilm1, ilm2
             lh = ll(ilh)
             do  11  icg = icgi(ilg,ilh),icge(ilg,ilh) 
                ilm = jcg(icg)
                lm = ll(ilm)
!                write(6,*)'xxxxxxxxvvvvvvvvv',ilg,ilh,lg,lh,ilm,nlm1
                k = (lg+lh-lm)/2
                fac = ee*(-1d0)**lg*cg(icg)
                do  ip = 0, kmax
                   s(ilg,ilh,ip)    = s(ilg,ilh,ip)    + fac*hkl(k+ip,ilm)
                   ds(ilg,ilh,ip,:) = ds(ilg,ilh,ip,:) + fac*dhkl(k+ip,ilm,:)
                enddo
11           enddo
111       enddo
1111   enddo
20  enddo
    deallocate(hkl,dhkl)
  end subroutine ghigbl
  subroutine hxpbl(ph,pg,q,rsmh,rsmg,eh,kmax,nlmh,nlmg,k0,ndim, c) ! Coefficients to expand smooth bloch hankels centered at ph into a sum of polynomials P_kL centered at pg.
    !i Inputs
    !i   ph    :Function is centered at ph; see Remarks
    !i   pg    :Function it expansed is at pg; see Remarks
    !i   q     :wave number for Bloch sum
    !i   rsmh  :vector of l-dependent smoothing radii of smoothed hankel
    !i         :rsmh must be specified for 1..ll(nlmh)+1
    !i   rsmg  :smoothing radius of gaussian
    !i   eh    :vector of l-dependent energies of smoothed Hankel
    !i         :eh must be specified for 1..ll(nlmh)+1
    !i   kmax  :polynomial cutoff
    !i   nlmh  :L-cutoff for smoothed Hankel functions being expanded
    !i   nlmg  :L-cutoff for P_kL expansion
    !i   k0    :leading dimension of coefficient array c
    !i   ndim  :second dimension of coefficient array c
    !i   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
    !i   indxcg:index for Clebsch Gordon coefficients
    !i   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
    !i   cy    :Normalization constants for spherical harmonics
    !i   slat  :struct containing information about the lattice
    !o Outputs
    !o   c     :cofficients c(k,M,L); see Remarks
    !o         :here k=0..kmax, M=1..nlmg, L=1..nlmh
    !r Remarks
    !r   Expansion is:  H_L(r-ph) = sum(kM) c(k,M,L) * P_kM(r-pg)
    !r   See J. Math. Phys. {\bf 39},3393 (1998), Eq. 12.17
    !r
    !r   As rsmg -> 0, expansion turns into a Taylor series of H_L.
    !r   As rsmg increases the error is spread out over a progressively
    !r   larger radius, thus reducing the error for larger r while
    !r   reducing accuracy at small r.
    !b Bugs
    !b   Doesn't properly handle case rsmh<0
    implicit none
    integer :: k0,kmax,ndim,nlmg,nlmh!,jcg(1),indxcg(1)
    real(8) :: eh(*),rsmg,rsmh(*),ph(3),pg(3),q(3) !,cg(1),cy(1)
    complex(8):: c(0:k0,ndim,nlmh)
    integer :: ndim1,ndim2,ilmg,ilmh,k,l,lmaxg,m,nm
    real(8) :: a,dfact,eg,fac,factk,fpi
    complex(8),allocatable:: s(:,:,:)
    if (nlmg == 0 .OR. nlmh == 0) return
    fpi = 16d0*datan(1d0)
    allocate(s(nlmg,nlmh,0:kmax))
    eg = 0d0
    call ghibl(pg,ph,q,rsmg,rsmh,eg,eh,nlmg,nlmh,kmax,nlmg,nlmh, s) !Integrals of gaussians and smoothed Hankels
    a = 1d0/rsmg
    lmaxg = ll(nlmg)
    call factorial_init(kmax,2*lmaxg+1)
    ilmg = 0
    do l = 0, lmaxg 
       nm = 2*l+1
       do m = 1, nm
          ilmg = ilmg+1
          do   k = 0, kmax ! ... Scale to get coefficients of the P_kL
             c(k,ilmg,:) = s(ilmg,:,k)*fpi / ((4*a*a)**k * a**l *factorial(k)*factorial2(2*l+1))
          enddo
       enddo
    enddo
    deallocate(s)
  end subroutine hxpbl
  subroutine hxpgbl(ph,pg,q,rsmh,rsmg,eh,kmax,nlmh,nlmg,k0,ndimh, ndimg, c,dc) ! Coefficients to expand smooth bloch hankels and grads centered at ph into a sum of polynomials P_kL centered at pg.
    !i Inputs
    !i   ph    :Function is centered at ph; see Remarks
    !i   pg    :Function it expansed is at pg; see Remarks
    !i   q     :wave number for Bloch sum
    !i   rsmh  :smoothing radius of smoothed Hankel (l-dependent)
    !i   rsmg  :smoothing radius of gaussian
    !i   eh    :energies of smoothed Hankel (l-dependent)
    !i   kmax  :polynomial cutoff
    !i   nlmh  :L-cutoff for smoothed Hankel functions being expanded
    !i   nlmg  :L-cutoff for P_kL expansion
    !i   k0    :leading dimension of coefficient arrays c,dc
    !i   ndimg :second dimension of coefficient arrays c,dc
    !i   ndimh :htird dimension of coefficient arrays c,dc
    !i   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
    !i   indxcg:index for Clebsch Gordon coefficients
    !i   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
    !i   cy    :Normalization constants for spherical harmonics
    !i   slat  :struct containing information about the lattice
    !o Outputs
    !o   c     :cofficients c(k,M,L); see Remarks
    !o   dc    :cofficients dc(k,M,L,1..3) (gradients of c); see Remarks
    !o         :here k=0..kmax, M=1..nlmg, L=1..nlmh
    !r Remarks
    !r   Expansion is:  H_L(r-ph) = sum(kM)  c(k,M,L) * P_kM(r-pg)
    !r             grad H_L(r-ph) = sum(kM) dc(k,M,L) * P_kM(r-pg)
    !r   As rsmg -> 0, expansion turns into a Taylor series of H_L.
    !r   As rsmg increases the error is spread out over a progressively
    !r   larger radius, thus reducing the error for larger r while
    !r   reducing accuracy at small r.
    !r
    !r   Grads are wrt to ph; take negative for grad wrt to pg.
    implicit none
    integer :: k0,kmax,ndimg,ndimh,nlmg,nlmh
    real(8) :: eh(*),rsmg,rsmh(*),ph(3),pg(3),q(3),a,dfact,eg,fac,factk,fpi
    complex(8):: c(0:k0,ndimg,ndimh),dc(0:k0,ndimg,ndimh,3)
    integer :: ndim1,ndim2,ilmg,ilmh,k,l,lmaxg,m,nm
    complex(8),allocatable:: s(:,:,:),ds(:,:,:,:)
    if (nlmg == 0 .OR. nlmh == 0) return
    fpi = 16d0*datan(1d0)
    ndim1 = nlmg
    ndim2 = nlmh
    allocate(s(ndim1,ndim2,0:kmax),ds(ndim1,ndim2,0:kmax,3))
    eg = 0d0
    call ghigbl(pg,ph,q,rsmg,rsmh,eg,eh,nlmg,nlmh,kmax,nlmg,nlmh,kmax,s,ds) ! ... Integrals of Hankels with Gaussians
    a = 1d0/rsmg ! ... Scale to get coefficients of the PkL
    lmaxg = ll(nlmg)
    call factorial_init(kmax,2*lmaxg+1)
    do ilmg = 1,(lmaxg+1)**2
       l=ll(ilmg)
       do  k = 0, kmax
          fac = fpi/ ((4*a*a)**k * a**l *factorial(k)*factorial2(2*l+1))
          c (k,ilmg,1:nlmh)     =   s(ilmg,1:nlmh,k)*fac
          dc(k,ilmg,1:nlmh,1:3) = -ds(ilmg,1:nlmh,k,1:3)*fac
       enddo
    enddo
    deallocate(s,ds)
  end subroutine hxpgbl
  subroutine gklq(lmax,rsm,q,p,e,kmax,k0,alat,nrx,yl, gkl)!Bloch sum of k,L-dependent gaussians (vectorizes)
    use m_ropyln,only: ropyln
    use m_shortn3_plat,only: shortn3_plat,nout,nlatout
    !i Inputs:
    !i  lmax   :l-cutoff for gkl
    !i   rsm   :smoothing radius
    !i   q     :wave number for Bloch sum (units of 2*pi/alat)
    !i   p     :connecting vector (units of alat)
    !i   e     :G_kL scaled by exp(e*rsm**2/4)
    !i   alat  :length scale of lattice and basis vectors, a.u.
    !i  dlv,nkd:direct lattice vectors, and number
    !i   nrx   :leading dimension of wk,yl
    !i   yl    :dimensioned nrx*(lmax+1)**2
    !i   nrx   :dimensions work arrays yl, and must be >= max(nkq,nkd)
    !i   k0    :leading dimension of gkl
    !o   gkl: G_kL * exp(e*rsm**2/4) generated for (0:kmax,0:lmax)
    implicit none
    integer :: k0,kmax,lmax,nrx,job, ilm,ir,k,l,m,nlm,job0,job1,li,le,lx
    real(8) :: alat,rsm,p(3),q(3),yl(nrx,(lmax+1)**2) !,dlv(3,nkd)
    complex(8) :: gkl(0:k0,(lmax+1)**2),img=(0d0,1d0),phase,sss(nkd)
    complex(8):: wkc1(nkd),wkc2(nkd),wkz(nkd)       !    wkz = wk(1:nkd,3)+img*wk(1:nkd,4) !phase
    real(8) :: qdotr,ta2,a2,g0fac,xx1,xx2,x1,sp,pp(3),wk1(nkd,0:lmax),wk2(nkd,0:lmax),wkfac(nkd),r2(nkd),e
    real(8),parameter:: pi  = 4*datan(1d0),tpi = 2*pi,y0  = 1/dsqrt(4*pi)
    if (kmax < 0 .OR. lmax < 0 .OR. rsm == 0d0) return
    nlm = (lmax+1)**2
    a2  = 1/rsm**2
    ta2 = 2*a2
    gkl=0d0
    !wkz =[(exp(img*tpi*sum(q*dlv(:,ir))),  ir=1,nkd)]
    !r2  = alat**2*[(sum((p-dlv(:,ir))**2),ir=1,nkd)] ! wk(1:nkd,1) !wk(:,1) is length**2 !wk(:,1)= alat**2*(p-dlv)**2 
    do ir=1,nkd
       wkz(ir)=exp(img*tpi*sum(q*dlv(:,ir)))
       r2(ir)=alat**2*sum((p-dlv(:,ir))**2)
    enddo   
    wkfac(1:nkd) = y0*dexp(-r2*a2) 
    kloop2:   do 301  k = 0, kmax, 2 ! --- Outer loop over k (in blocks of 2), and over l ---
       lloop: do  30  l = 0, lmax
          g0fac = 1/rsm*ta2**(l+1)/pi * dexp(e*rsm*rsm/4)
          if (k == 0) then !Make radial part of the G_kl(1..nkd) for k= 0, 1
             wk1(:,l) = [(g0fac*wkfac(ir),ir=1,nkd)]
             wk2(:,l) = [((ta2*r2(ir)-3-2*l)* ta2,ir=1,nkd)]*wk1(:,l)
          else ! Make radial part of the G_kl(1..nkd) for k, k+1 from k-1, k-2 and cos(q.dlv) * G_kl and sin(q.dlv) * G_kl
             associate( x => 2*(k-1)*(2*k + 2*l-1), y => 4*k + 2*l-1, x2 => 2*k*(2*(k+1) + 2*l-1), y2 => 4*(k+1) + 2*l-1)
               do  ir = 1, nkd
                  xx1 = ta2*((ta2*r2(ir)-y)*wk2(ir,l)  - x*ta2*wk1(ir,l))
                  xx2 = ta2*((ta2*r2(ir)-y2)*xx1       -x2*ta2*wk2(ir,l))
                  wk1(ir,l) = xx1
                  wk2(ir,l) = xx2
               enddo
             endassociate
          endif
          li=l*l+1
          le=l*l+2*l+1 !For each point, add G_kl Y_L exp(i q.dlv) into Bloch G_kL
          !  gkl(k,li:le) =             gkl(k,li:le)   +[(sum([(wkz(ir)*wk1(ir,l)*yl(ir,ilm),ir=nkd,1,-1)]), ilm=li,le)] !bugfix for nvfortran24.1
          !  if(k<kmax) gkl(k+1,li:le)= gkl(k+1,li:le) +[(sum([(wkz(ir)*wk2(ir,l)*yl(ir,ilm),ir=nkd,1,-1)]), ilm=li,le)]
          do lx=li,le
             sss=[(wkz(ir)*wk1(ir,l)*yl(ir,lx),ir=nkd,1,-1)]
             gkl(k,lx) =             gkl(k,lx)  + sum(sss)
             sss=[(wkz(ir)*wk2(ir,l)*yl(ir,lx),ir=nkd,1,-1)]
             if(k<kmax) gkl(k+1,lx)= gkl(k+1,lx)+ sum(sss)
          enddo
30     enddo lloop
301 enddo kloop2
  end subroutine gklq
  subroutine hxpos(rsmh,rsmg,eh,kmax,nlmh,k0, c) ! Coefficients to expand smooth hankels at (0,0,0) into P_kl's.
    !i Inputs
    !i   rsmh  :vector of l-dependent smoothing radii of smoothed hankel
    !i         :rsmh must be specified for 1..ll(nlmh)+1
    !i   rsmg  :smoothing radius of gaussian
    !i   eh    :vector of l-dependent energies of smoothed Hankel
    !i         :eh must be specified for 1..ll(nlmh)+1
    !i   kmax  :polynomial cutoff
    !i   nlmh  :L-cutoff for smoothed Hankel functions and P_kL
    !i   k0    :leading dimension of coefficient array c
    !o Outputs
    !o   c     :structure constants c(k,M,L); see Remarks
    !o         :c(*,ilm) for rsmh(ll(ilm))<=0 are SET TO ZERO
    !r Remarks
    !r   Expansion is:  H_L = sum(k) c(k,L) * P_kL
    !r   As rsmg -> 0, expansion turns into a Taylor series of H_L.
    !r   As rsmg increases the error is spread out over a progressively
    !r   larger radius, thus reducing the error for larger r while
    !r   reducing accuracy at small r.
    !r
    !r    Only diagonal elements ilmh=ilmg are nonzero and returned.
    !r    Routine is equivalent to hxpml with ph=pg.
    integer :: k0,kmax,nlmh,ik,i,ilm1,ilm2
    real(8) :: eh(*),rsmg,rsmh(*),c(0:k0,nlmh)
    integer :: ndim,ilm,k,l,lmax,m,nm,ktop,l1,l2,lmaxh
    real(8) :: a,dfact,eg,fac,factk,sig,rsm
    real(8):: s(nlmh,0:kmax),e,gamg,gamh,rsmx,gam,asm,akap,hh,gg
    real(8),parameter:: fpi = 16d0*datan(1d0), pi = 4d0*datan(1d0),y0 = 1d0/dsqrt(4*pi)
    if (nlmh == 0) return
    eg = 0d0
    lmaxh = ll(nlmh)
    l2 = -1
    do  20  l1  = 0, lmaxh ! --- Loop over sequences of l with a common rsm,e ---
       if (l1 <= l2) cycle
       l2=gtbsl2(l1,lmaxh,eh,rsmh)
       rsm  = rsmh(l1+1)
       e    = eh(l1+1)
       if (rsm <= 0 .OR. e > 0) cycle
       ilm1 = l1**2+1
       ilm2 = (l2+1)**2
       gamh = 0.25d0*rsm*rsm
       gamg = 0.25d0*rsmg*rsmg
       rsmx = 2d0*dsqrt(gamg+gamh)
       ktop = l2+kmax
       block 
         real(8) :: h0k(0:ktop) !   ... Make hankels for l=0 and k=0..kmax
         !call hklos(rsmx,e,ktop,h0k)
         gam = rsmx*rsmx/4d0
         asm = 0.5d0/dsqrt(gam)
         akap = sqrt(abs(e))
         hh = 4d0*asm**3*exp(gam*e)/sqrt(pi)/(2*asm*asm)-akap*erfc(akap/(2*asm))
         h0k(0) = hh*y0
         gg = exp(gam*e) * (asm*asm/pi)**1.5d0
         do k = 1,ktop ! ... Upward recursion for laplace**k h0
            hh = -e*hh - 4*pi*gg
            h0k(k) = hh*y0 ! k,L-dependent smooth hankel functions at (0,0,0)
            !r   Only the functions for l=0 are generated (remaining are zero). Equivalent to hkl_ml for p=(0,0,0) and lmax=0.
            gg = -2*asm*asm*(2*k+1)*gg
         enddo
         do  ilm = ilm1, ilm2 ! .. Evaluate what is left of Clebsch-Gordan sum
            l = ll(ilm) 
            s(ilm,0:kmax) = y0*dexp(gamg*(eg-e)) * (-1)**l * h0k(l+0:l+kmax) ! Make integrals with gaussians G_kL
         enddo
       endblock
20  enddo
!    call ghios(rsmg,rsmh,eg,eh,nlmh,kmax,nlmh,s)
    a = 1d0/rsmg    ! ... Scale integrals to get coefficients of the P_kL
    lmax = ll(nlmh)
    call factorial_init(kmax,2*lmax+1)
    c=0d0
    do  l = 0, lmax
       if(rsmh(l+1) <= 0) cycle
       do i=l**2+1,l**2+2*l+1
          c(0:kmax,i) = [(s(i,k)*fpi/( (4*a*a)**k * a**l* factorial(k)*factorial2(2*l+1)), k=0,kmax)]
       enddo   
    enddo
  end subroutine hxpos
  subroutine hsmbl(p,rsm,e,q,lmax, hsm,hsmp)  !Bloch-sum of smooth Hankel functions and energy derivatives
    use m_lmfinit,only: alat=>lat_alat,tol=>lat_tol
    use m_shortn3_plat,only: shortn3_plat,nout,nlatout
    !i Inputs
    !i   p     :Function is centered at p
    !i   rsm   :smoothing radius
    !i   e     :energy of smoothed Hankel
    !i   q     :wave number for Bloch sum
    !i   lmax  :l-cutoff for hsm
    !i   cy    :Normalization constants for spherical harmonics
    !i   slat  :struct containing information about the lattice
    !o Outputs
    !o   hsm   :Bloch-summed smoothed Hankels
    !o   hsmp  :Energy derivatives of hsm
    !r Remarks
    !r  Hankel functions evaluated by Ewald summation.
    !r  p in units of alat, qlv in units of 2 pi / alat
    !r  See also hsmq for a vectorized version.
    implicit none
    integer :: lmax
    real(8):: rsm,e,q(3),p(3),ppin(3)
    complex(8):: hsm(*),hsmp(*)
    integer:: ilm,l,m
    real(8) :: p1(3),sp,rwald,asm
    real(8),parameter:: pi = 4d0*datan(1d0)
    complex(8):: cfac,phase,img=(0d0,1d0)
    ppin=matmul(transpose(qlat),p) 
    call shortn3_plat(ppin) 
    p1= matmul(plat,ppin+nlatout(:,1)) ! p is shortened to be p1
    phase = exp(img*2*pi*sum(q*(p-p1)))
    rwald = 1d0/awald
    asm = 1d0/rsm
    if (rsm < rwald) then
       call hsmblq(p1,     e,q, awald, lmax, alat, vol, hsm,hsmp )
       call hsmbld(p1, rsm,e,q, awald, lmax, alat,      hsm,hsmp )
    else
       call hsmblq(p1,     e,q, asm,   lmax, alat, vol, hsm,hsmp )
    endif
    ! ... Multiply by phase to undo shortening
    do ilm=1,(lmax+1)**2
       l=ll(ilm)
       hsm(ilm)  = phase*cy(ilm)*hsm(ilm)*(-img)**l
       hsmp(ilm) = phase*cy(ilm)*hsmp(ilm)*(-img)**l
    enddo
  end subroutine hsmbl
  subroutine hsmblq(p,e,q,a,lmax,alat,vol, dl,dlp) ! k-space part of smooth hankel bloch sum
    use m_qplist,only:qshortn
    !i   p     :Function is centered at p (units of alat)
    !i   e     :energy of smoothed Hankel
    !i   q     :wave number for Bloch sum
    !i   a     :Ewald smoothing parameter
    !i   lmax  :l-cutoff for dl,dlp
    !i   alat  :length scale of lattice and basis vectors, a.u.
    !i   qlv   :reciprocal lattice vectors, units of 2pi/alat
    !i   nkq   :number of r.l.v.
    !i   vol   :cell volume
    !o Outputs
    !i   dl    :k-summed smoothed Bloch hankels
    !i   dlp   :energy derivative of dl
    implicit none
    integer :: lmax !nkq,
    real(8) :: a,alat,e,vol, q(3),p(3)
    complex(8):: dl((lmax+1)**2),dlp((lmax+1)**2),img=(0d0,1d0)
    integer ::nlm,ilm,ir
    real(8) :: r(3),tpi,gamma,fpibv,tpiba,scalp,r2,den0,den1
    real(8) :: yl((lmax+1)**2)
    complex(8):: eiphi
    real(8) :: qshortn_q(3)
    tpi = 8d0*datan(1d0)
    gamma = .25d0/(a*a)
    fpibv = 2d0*tpi/vol
    tpiba = tpi/alat
    nlm = (lmax+1)**2
    dl=0d0
    dlp=0d0
    !MO moved qshortn outside of ir loop 2024-11-08
    qshortn_q = qshortn(q)
    do ir = 1, nkq
       ! r(:) = tpiba*(qshortn(q)+qlv(:,ir))
       r(:) = tpiba*(qshortn_q+qlv(:,ir))
       eiphi = exp(img*alat*sum(r*p))
       call sylm(r,yl,lmax,r2)
       den0 = dexp(gamma*(e-r2))/(r2-e)
       dl(1:nlm)  = dl(1:nlm)  + yl(1:nlm)*eiphi*den0
       dlp(1:nlm) = dlp(1:nlm) + yl(1:nlm)*eiphi*den0/(r2-e)
    enddo
    dl(1:nlm)  = fpibv*dl(1:nlm)
    dlp(1:nlm) = fpibv*dlp(1:nlm) + gamma*dl(1:nlm)
  end subroutine hsmblq
  subroutine hsmbld(p,rsm,e,q,a,lmax,alat, dl,dlp) !Adds real space part of reduced structure constants (ewald).
    implicit none
    integer :: ilm,ir,l,lmax,m,nm
    real(8) :: q(3),p(3)
    complex(8):: dl(*),dlp(*)
    real(8) :: a,a2,akap,alat,asm,asm2,cc,ccsm,derfc,e,emkr,gl,qdotr,r1,r2,rsm,srpi,ta,ta2,tasm,tasm2,umins,uplus,tpi,kap
    real(8) :: yl((lmax+1)**2),chi1(-1:10),chi2(-1:10),r(3)
    complex(8):: cfac,zikap,expikr,zuplus,zerfc,img=(0d0,1d0)
    real(8),parameter :: srfmax=16d0, fmax=srfmax*srfmax
    if (e>0d0) call rx('EH >0 not supported') !lpos removed
    tpi = 8.d0*datan(1.d0)
    srpi = dsqrt(tpi/2.d0)
    akap = dsqrt(-e)
    ta = 2d0*a
    a2 = a*a
    ta2 = 2d0*a2
    cc = 4d0*a2*a*dexp(e/(ta*ta))/srpi
    asm = 1d0/rsm
    tasm = 2d0*asm
    asm2 = asm*asm
    tasm2 = 2d0*asm2
    ccsm = 4d0*asm2*asm*dexp(e/(tasm*tasm))/srpi
    do  20  ir = 1, nkd
       r = alat*(p-dlv(:,ir))
       call sylm(r,yl,lmax,r2)
       r1 = dsqrt(r2)
       if (r1 < 1d-6) then !Make the xi's from -1 to lmax ---
          chi1(-1:lmax) = [derfc(akap/ta)/akap,   ta*dexp(e/(2d0*ta2))/srpi -akap*derfc(akap/ta),      (0d0,l=1,lmax)]
          chi2(-1:lmax) = [derfc(akap/tasm)/akap, tasm*dexp(e/(2d0*tasm2))/srpi -akap*derfc(akap/tasm),(0d0,l=1,lmax)]
       else
          !         If large, these are -log(uplus),-log(umins); then chi->0
          !r        If both (akap*rsm/2 -/+ r/rsm) >> 1, we have
          !r        -log u(+/-) -> (akap*rsm/2 -/+ r/rsm)^2 -/+ akap*r   =  (akap*rsm/2)^2 + (r/rsm)^2 >> 1
          !r         u(+/-)     -> exp[-(akap*rsm/2)^2 - (r/rsm)^2] -> 0
          !r        Also if akap*r >> 1,   chi < dexp(-akap*r1) -> 0
          emkr = dexp(-akap*r1)
          if (.5d0*akap/a+r1*a > srfmax .AND. .5d0*akap/a-r1*a > srfmax .OR. akap*r1 > fmax) then
             chi1(-1:lmax) = 0
          else
             uplus = derfc(.5d0*akap/a+r1*a)/emkr
             umins = derfc(.5d0*akap/a-r1*a)*emkr
             chi1(-1:0) = [(umins+uplus)/(2.d0*akap), 0.5d0*(umins-uplus)/r1]
             gl = cc*dexp(-a2*r2)/ta2
             do l = 1, lmax !resursion
                chi1(l) = ((2*l-1)*chi1(l-1)-e*chi1(l-2)-ta2**(l-1)*gl)/r2
             enddo
          endif
          !        chi2 is complex; so is chi1, but the imaginary part    is the same, so the difference is real
          if(.5d0*akap/asm+r1*asm > srfmax .AND..5d0*akap/asm-r1*asm > srfmax .OR. akap*r1 > fmax)then
             chi2(-1:lmax) = 0
          else
             uplus = derfc(.5d0*akap/asm+r1*asm)/emkr
             umins = derfc(.5d0*akap/asm-r1*asm)*emkr
             chi2(-1:0) = [(umins+uplus)/(2d0*akap), 0.5d0*(umins-uplus)/r1]
             gl = ccsm*dexp(-asm2*r2)/tasm2
             do l = 1, lmax
                chi2(l) = ((2*l-1)*chi2(l-1)-e*chi2(l-2)-tasm2**(l-1)*gl)/r2
             enddo
          endif
       endif
       cfac = exp(img*tpi*(q(1)*dlv(1,ir)+q(2)*dlv(2,ir)+q(3)*dlv(3,ir)))
       do ilm = 1,(lmax+1)**2
          l=ll(ilm)
          dl(ilm)  = dl(ilm)  + yl(ilm)*(chi2(l)-chi1(l))*cfac  *img**l
          dlp(ilm) = dlp(ilm) + yl(ilm)*0.5d0*(chi2(l-1)-chi1(l-1))*cfac*img**l
       enddo
20  enddo
  end subroutine hsmbld
  subroutine ropylg2(lmax2,kmax,nlm,kmax0,nlm0,hkl, ghkl) !ghkl are derivatives of hkl wrt (x,y,z)
    use m_scg,only:scglp1
    implicit none !hkl(k,ilm) \propto r^k Y_ilm
    integer :: kmax,nlm,kmax0,nlm0,ilm,k,kx1,kx2,ky1,ky2,kz,lmax2,m,nlm1
    complex(8):: hkl(0:kmax0,nlm0),ghkl(0:kmax0,nlm0,3)
    real(8) :: cx1,cx2,cy1,cy2,cz
    ghkl=0d0 
    do  ilm = 1, nlm
       call scglp1(ilm,kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2) !see ropylg.f90
       do  k = 0, kmax
          ghkl(k,ilm,1:3)=[ghkl(k,ilm,1)-cx1*hkl(k,kx1)-cx2*hkl(k,kx2), &
                           ghkl(k,ilm,2)-cy1*hkl(k,ky1)-cy2*hkl(k,ky2), ghkl(k,ilm,3)-cz*hkl(k,kz)]
          if (ilm <= lmax2) then
             ghkl(k,kx1,1) = ghkl(k,kx1,1) - cx1*hkl(k+1,ilm)
             ghkl(k,kx2,1) = ghkl(k,kx2,1) - cx2*hkl(k+1,ilm)
             ghkl(k,ky1,2) = ghkl(k,ky1,2) - cy1*hkl(k+1,ilm)
             ghkl(k,ky2,2) = ghkl(k,ky2,2) - cy2*hkl(k+1,ilm)
             ghkl(k,kz,3) =  ghkl(k,kz,3)  - cz *hkl(k+1,ilm)
          endif
       enddo
    enddo
  end subroutine ropylg2
  pure function gtbsl2(l1,lmxh,eh,rsmh) result(l2) ! Returns the highest l with rsm,e common to those of a given l
    use m_nvfortran,only:findloc
    implicit none
    intent(in)::       l1,lmxh,eh,rsmh
    integer :: l1,l2,lmxh
    real(8) :: rsmh(0:lmxh),eh(0:lmxh), e,rsm
    logical:: lll(l1:lmxh-1)
    lll=[(rsmh(l2+1)/=rsmh(l1).or.eh(l2+1)/=eh(l1), l2=l1,lmxh-1)]
    l2 = findloc(lll,value=.true.,dim=1) + l1-1
    if(l2==l1-1) l2=lmxh
  end function gtbsl2
end module m_smhankel
