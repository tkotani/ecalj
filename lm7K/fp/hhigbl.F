Cgetarg...info...           structure ['phhigb', 'lat', 'slat']
Cgetarg...info...           use_to_add ['phhigb', 'm_struc_def']
      subroutine hhigbl(mode,p1,p2,q,rsm1,rsm2,e1,e2,nlm1,nlm2,kmax,
     .  ndim1,ndim2,k0,cg,indxcg,jcg,cy,slat,s,ds)
      use m_struc_def
C- Integrals between smooth hankels with k-th power of Laplace operator
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1's digit
Ci         :0 rsm1,rsm2,e1,e2 are scalars
Ci         :1 rsm1,rsm2,e1,e2 are l-dependent vectors
Ci         :10's digit
Ci         :1: do not calculate strux for any (l1,l2) pairs for
Ci         :   if rsm(l1) or rsm(l2) is zero.
Ci   p1    :first center
Ci   p2    :second center
Ci   q     :wave number for Bloch sum
Ci   rsm1  :smoothing radii of Hankels at p1 (l-dependent)
Ci   rsm2  :smoothing radii of Hankels at p2 (l-dependent)
Ci   e1    :energies of smooth Hankels at p1 (l-dependent)
Ci   e2    :energies of smooth Hankels at p2 (l-dependent)
Ci   nlm1  :L-cutoff for functions at p1
Ci   nlm2  :L-cutoff for functions at p2
Ci   kmax  :cutoff in power of Laplace operator
Ci   ndim1 :leading dimensions of s,ds
Ci   ndim2 :second dimensions of s,ds
Ci   k0    :dimensions s,ds
Ci   cg    :Clebsch Gordon coefficients (scg.f)
Ci   indxcg:index for Clebsch Gordon coefficients
Ci   jcg   :L q.n. for the C.G. coefficients (scg.f)
Ci   cy    :Normalization constants for spherical harmonics
Ci   slat  :struct containing information about the lattice
Co Outputs
Co   s     :integrals between smooth Bloch Hankels; see Remarks
Co   ds    :gradient of s; see Remarks
Cr Remarks
Cr  s(L,M,k) contains integral of H_L^*(r-p1) (laplace)^k H_M(r-p2)
Cr  Gradient is wrt p1; use -ds for grad wrt p2.
Cu Updates
Cu   8 Jun 00  Added 10s digit mode.
Cu   23 May 00 Made rsm1,e1,rsm2,e2 l-dependent
Cu   22 Apr 00 Adapted from nfp hhig_bl.f
Cu   28 Apr 97 adapted to handle case q=0 and e1 or e2 zero.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,jcg(1),indxcg(1),nlm1,nlm2,kmax,ndim1,ndim2,k0
Cgetarg       double precision p1(3),p2(3),q(3),cg(1),cy(1),slat(1),
Cgetarg      .  rsm1(0:*),rsm2(0:*),e1(0:*),e2(0:*)
       real(8):: p1(3) , p2(3) , q(3) , cg(1) , cy(1) , rsm1(0:*) , 
     . rsm2(0:*) , e1(0:*) , e2(0:*) 
       type(s_lat)::slat

      double complex s(ndim1,ndim2,0:k0),ds(ndim1,ndim2,0:k0,3)
C ... Local parameters
      integer ilm,jlm,k,l1,l1t,l2,l2t,ll,lm11,lm12,lm21,lm22,lmx1,lmx2,m
      double precision dr(3)

      if (nlm1 .eq. 0 .or. nlm2 .eq. 0) return

      do  1  m = 1, 3
    1 dr(m) = p1(m)-p2(m)
      lmx1 = ll(nlm1)
      lmx2 = ll(nlm2)

      do  3  k = 0, kmax
      do  3  jlm = 1, nlm2
      do  3  ilm = 1, nlm1
      s(ilm,jlm,k)    = dcmplx(0d0,0d0)
      ds(ilm,jlm,k,1) = dcmplx(0d0,0d0)
      ds(ilm,jlm,k,2) = dcmplx(0d0,0d0)
    3 ds(ilm,jlm,k,3) = dcmplx(0d0,0d0)

      l1t = -1
      do  20  l1 = 0, lmx1
        if (l1 .le. l1t) goto 20
        if (mod(mode,10) .eq. 0) then
          l1t = lmx1
        else
          call gtbsl2(l1,lmx1,e1,rsm1,l1t)
        endif
C       l1t = l1

        l2t = -1
        do  22  l2 = 0, lmx2
          if (l2 .le. l2t) goto 22
          if (mod(mode,10) .eq. 0) then
            l2t = lmx2
          else
            call gtbsl2(l2,lmx2,e2,rsm2,l2t)
          endif
C         l2t = l2

          lm11 = l1**2+1
          lm12 = (l1t+1)**2
          lm21 = l2**2+1
          lm22 = (l2t+1)**2
          if (mode/10.eq.1 .and. rsm1(l1)*rsm2(l2) .eq. 0) goto 22
          call phhigb(dr,q,rsm1(l1),rsm2(l2),e1(l1),e2(l2),lm11,lm12,
     .      lm21,lm22,kmax,ndim1,ndim2,k0,cg,indxcg,jcg,cy,slat,s,ds)
   22   continue
   20 continue
Cgetarg       end
       end subroutine hhigbl 

      subroutine phhigb(dr,q,rsm1,rsm2,e1,e2,mlm1,nlm1,mlm2,nlm2,
     .  kmax,ndim1,ndim2,k0,cg,indxcg,jcg,cy,slat,s,ds)
       
       use m_struc_def  !Cgetarg

C- Integrals between smooth hankels with k-th power of Laplace operator
C ----------------------------------------------------------------------
Ci Inputs
Ci   dr    :p1-p2
Ci   q     :wave number for Bloch sum
Ci   rsm1  :smoothing radius of Hankels at p1
Ci   rsm2  :smoothing radius of Hankels at p2
Ci   e1    :energy of smooth Hankels at p1
Ci   e2    :energy of smooth Hankels at p2
Ci   mlm1  :lower L-cutoff for functions at p1
Ci   nlm1  :upper L-cutoff for functions at p1
Ci   mlm2  :lower L-cutoff for functions at p2
Ci   nlm2  :upper L-cutoff for functions at p2
Ci   kmax  :cutoff in power of Laplace operator
Ci   ndim1 :leading dimensions of s,ds
Ci   ndim2 :second dimensions of s,ds
Ci   k0    :dimensions s,ds
Ci   cg    :Clebsch Gordon coefficients (scg.f)
Ci   indxcg:index for Clebsch Gordon coefficients
Ci   jcg   :L q.n. for the C.G. coefficients (scg.f)
Ci   cy    :Normalization constants for spherical harmonics
Ci   slat  :struct containing information about the lattice
Co Outputs
Co   s     :integrals between smooth Bloch Hankels; see Remarks
Co   ds    :gradient of s; see Remarks
Cr Remarks
Cr  s(L,M,k) contains integral of H_L^*(r-p1) (laplace)^k H_M(r-p2)
Cr  Gradient is wrt p1; use -ds for grad wrt p2.
Cu Updates
Cu   22 Apr 00 Adapted from nfp hhig_bl.f
Cu   28 Apr 97 adapted to handle case q=0 and e1 or e2 zero.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer jcg(1),indxcg(1),nlm1,nlm2,kmax,ndim1,ndim2,k0,mlm1,mlm2
Cgetarg       double precision q(3),cg(1),cy(1),slat(1),dr(3),rsm1,rsm2,e1,e2
       real(8):: q(3) , cg(1) , cy(1) , dr(3) , rsm1 , rsm2 , e1 , e2 
       type(s_lat)::slat

      double complex s(ndim1,ndim2,0:k0),ds(ndim1,ndim2,0:k0,3)
C ... Local parameters
      integer nlm0,ktop0,lmax1,ll,lmax2,lmaxx,nlmx,nlmxp1,ktop,ktopp1,
     .  k,ilm,kz,kx1,kx2,ky1,ky2,ilm1,l1,ilm2,l2,ii,indx,icg1,icg2,
     .  icg,lm,ip
      double precision gam1,fpi,gam2,gamx,rsmx,qq,fac1,fac2,e,cz,cx1,
     .  cx2,cy1,cy2,fac,add,vol,dgets
      parameter (nlm0=100, ktop0=10)
      double complex hkl1(0:ktop0,nlm0),hkl2(0:ktop0,nlm0),
     .  ghkl(0:ktop0,nlm0,3),hsm(nlm0),hsmp(nlm0)

Cgetarg       vol = dgets('lat vol',slat)
       vol = slat%vol 

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
      if (nlmxp1 .gt. nlm0)  call rxi('hhigbl: need nlm0 ge',nlmxp1)
      if (ktopp1 .gt. ktop0) call rxi('hhigbl: need ktop0 ge',ktopp1)

C --- Set up functions for connecting vector p2-p1 ---
      if (dabs(e1-e2) .gt. 1d-5) then
        fac1 = dexp(gam2*(e2-e1))/(e1-e2)
        fac2 = dexp(gam1*(e1-e2))/(e2-e1)
        if (dabs(e1) .gt. 1d-6) then
          call hklbl(dr,rsmx,e1,q,ktopp1,nlmxp1,ktop0,cy,slat,hkl1)
        else
          if (qq .gt. 1d-6) call rx('hhigbl: e1=0 only allowed if q=0')
          call fklbl(dr,rsmx,ktopp1,nlmxp1,ktop0,cy,slat,hkl1)
        endif
        if (dabs(e2) .gt. 1d-6) then
          call hklbl(dr,rsmx,e2,q,ktopp1,nlmxp1,ktop0,cy,slat,hkl2)
        else
          if (qq .gt. 1d-6) call rx('hhigbl: e2=0 only allowed if q=0')
          call fklbl(dr,rsmx,ktopp1,nlmxp1,ktop0,cy,slat,hkl2)
        endif
        do  4  ilm = 1, nlmxp1
          do  5  k = 0, ktopp1
            hkl1(k,ilm) = fac1*hkl1(k,ilm) + fac2*hkl2(k,ilm)
    5     continue
    4   continue
      else
        e = .5d0*(e1+e2)
        if (qq.lt.1d-6 .and. dabs(e) .lt. 1d-6)
     .    call rx('hhigbl: case q=0 and e1=e2=0 not available')
        call hklbl(dr,rsmx,e,q,ktopp1,nlmxp1,ktop0,cy,slat,hkl2)
        call hsmbl(dr,rsmx,e,q,lmaxx+1,cy,slat,hsm,hsmp)
        do  6  ilm = 1, nlmxp1
          hkl1(0,ilm) = hsmp(ilm) - gamx*hsm(ilm)
          do  7  k = 1, ktopp1
            hkl1(k,ilm) = -e*hkl1(k-1,ilm) - hkl2(k-1,ilm)
    7     continue
    6   continue
      endif

C ... Make gradients using CGs for p functions
      call dpzero(ghkl, 2*(ktop0+1)*nlm0*3)
      do  ilm = 1, nlmx
        call scglp1(ilm,kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2)
        do  k = 0, ktop
          ghkl(k,ilm,1) = ghkl(k,ilm,1)-cx1*hkl1(k,kx1)-cx2*hkl1(k,kx2)
          ghkl(k,ilm,2) = ghkl(k,ilm,2)-cy1*hkl1(k,ky1)-cy2*hkl1(k,ky2)
          ghkl(k,ilm,3) = ghkl(k,ilm,3)-cz*hkl1(k,kz)
          if (ilm .le. lmaxx*lmaxx) then
            ghkl(k,kx1,1) = ghkl(k,kx1,1) - cx1*hkl1(k+1,ilm)
            ghkl(k,kx2,1) = ghkl(k,kx2,1) - cx2*hkl1(k+1,ilm)
            ghkl(k,ky1,2) = ghkl(k,ky1,2) - cy1*hkl1(k+1,ilm)
            ghkl(k,ky2,2) = ghkl(k,ky2,2) - cy2*hkl1(k+1,ilm)
            ghkl(k,kz,3)  = ghkl(k,kz,3)  - cz *hkl1(k+1,ilm)
          endif
        enddo
      enddo

C ... Combine with Clebsch-Gordan coefficients
      do  11  ilm1 = mlm1, nlm1
      l1 = ll(ilm1)
      do  11  ilm2 = mlm2, nlm2
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
  12  ds(ilm1,ilm2,ip,3) = ds(ilm1,ilm2,ip,3) + fac*ghkl(k+ip,ilm,3)
  11  continue

C ... Extra term when e1 or e2 is zero (but not both)
      if (mlm1 .eq. 1 .and. mlm2 .eq. 1) then
        add = 0d0
        if (dabs(e2) .lt. 1d-6) add = fpi*dexp(gam1*e1)/(vol*e1*e1)
        if (dabs(e1) .lt. 1d-6) add = fpi*dexp(gam2*e2)/(vol*e2*e2)
        s(1,1,0) = s(1,1,0) + add
      endif

Cgetarg       end
       end subroutine phhigb 


