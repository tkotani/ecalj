CSFCPP#define F90 1
      subroutine hklbl(p,rsm,e,q,kmax,nlm,k0,cy,hkl) !slat,
      use m_lattic,only: rv_a_oqlv,rv_a_odlv
      use m_struc_def  !Cgetarg
      use m_lmfinit,only: lat_awald,lat_tol,lat_vol,lat_nkd, lat_nkq,
     &   lat_alat,lat_plat,lat_qlat
C- Bloch-sums of k,L-dependent smooth Hankel functions.
C ----------------------------------------------------------------------
Ci Inputs
Ci   p     :Function is centered at p
Ci   rsm   :smoothing radius
Ci   e     :energy of smoothed Hankel
Ci   q     :wave number for Bloch sum
Ci   kmax  :polynomial cutoff
Ci   nlm   :L-cutoff for hkl
Ci   k0    :leading dimension of hkl
Ci   cy    :Normalization constants for spherical harmonics
Ci   slat  :struct containing information about the lattice
Co Outputs
Co   hkl   :Bloch-summed smoothed Hankels
Cr Remarks
Cr   H_kL = laplace^k H_L
Cr   Uses the recursion relation H_k+1,L = -e*H_kL - 4*pi*G_kL
Cu Updates
Cu   24 Apr 00 Adapted from nfp hkl_bl.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer k0,kmax,nlm
      real(8):: e , rsm , q(3) , p(3) , cy(1)
c      type(s_lat)::slat
      double complex hkl(0:k0,nlm)
C ... Local parameters
      integer:: nlm0 , ilm , job , k , ll , lmax , nkd , nkq , nrx 
     .  , owk , oyl
      parameter (nlm0=144)
      double precision alat,awald,fpi,pi,sp,tol,vol,plat(3,3),
     .qlat(3,3),p1(3)
      double complex hsm(nlm0),hsmp(nlm0),phase,gklsav,gklnew
      real(8),allocatable:: wk(:),yl(:)
      double precision faca
      parameter (faca=1d0)

      if (nlm .eq. 0) return
      if (nlm .gt. nlm0) call rxi('increase nlm0 in hklbl need',nlm)
      pi = 4d0*datan(1d0)
      fpi = 4d0*pi
      lmax = ll(nlm)
C ... Shorten p, multiply by corresponding phase later
      alat=lat_alat
      plat=lat_plat
      qlat=lat_qlat
      call shorbz(p,p1,plat,qlat)
      sp = 2*pi*(q(1)*(p(1)-p1(1))+q(2)*(p(2)-p1(2))+q(3)*(p(3)-p1(3)))
      phase = dcmplx(dcos(sp),dsin(sp))
      awald=lat_awald
      tol=lat_tol
      vol=lat_vol
      nkd=lat_nkd
      nkq=lat_nkq
      nrx = max(nkd,nkq)
      allocate(wk(nrx*(2*lmax+10)),yl(nrx*(lmax+1)**2))
      call hsmq ( 1 , 0 , ll ( nlm ) , e , rsm , 0000 , q , p1 , nrx
     .  , nlm , wk , yl , awald , alat , rv_a_oqlv , nkq , rv_a_odlv 
     .  , nkd , vol , hsm , hsmp )
      if (rsm .gt. faca/awald) then
        call gklbl(p1,rsm,e,q,kmax-1,nlm,k0,cy, hkl) !slat,
      else
        job = 2
        call gklq ( lmax , rsm , q , p1 , e , kmax - 1 , k0 , alat ,
     .   rv_a_odlv , nkd , nrx , yl , wk , job , hkl )
      endif
      deallocate(wk,yl)
C --- H_kL by recursion ---
      do  2  ilm = 1, nlm
        gklsav = hkl(0,ilm)
        hkl(0,ilm) = hsm(ilm)
      do  2  k = 1, kmax
        gklnew = hkl(k,ilm)
        hkl(k,ilm) = -e*hkl(k-1,ilm) - fpi*gklsav
        gklsav = gklnew
    2 continue
C ... Put in phase to undo shortening
      if (sp .ne. 0) then
        do  20  ilm = 1,nlm
          do  22  k = 0,kmax
            hkl(k,ilm) = phase*hkl(k,ilm)
   22     continue
   20   continue
      endif
      end subroutine hklbl


