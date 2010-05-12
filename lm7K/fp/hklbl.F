Cgetarg...info...           structure ['hklbl', 'lat', 'slat']
Cgetarg...info...           use_to_add ['hklbl', 'm_struc_def']
#define F90 1
      subroutine hklbl(p,rsm,e,q,kmax,nlm,k0,cy,slat,hkl)

      use m_struc_def  !Cgetarg

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
C     implicit none
C ... Passed parameters
      integer k0,kmax,nlm
Cgetarg       double precision e,rsm,slat(1),q(3),p(3),cy(1)
      real(8):: e , rsm , q(3) , p(3) , cy(1)
      type(s_lat)::slat

      double complex hkl(0:k0,nlm)
C ... Local parameters
      integer nlm0,ilm,job,k,ll,lmax,nkd,nkq,nrx,odlv,oqlv,owk,oyl
      parameter (nlm0=144)
      double precision alat,awald,fpi,pi,sp,tol,vol,plat(3,3),
     .  qlat(3,3),p1(3)
      double complex hsm(nlm0),hsmp(nlm0),phase,gklsav,gklnew
C     double complex hsmx(nlm0),hsmpx(nlm0)
C ... Heap
      integer w(1)
      common /w/ w

#if F90
      real(8),allocatable:: wk(:),yl(:)
#endif
      double precision faca
      parameter (faca=1d0)

      if (nlm .eq. 0) return
      if (nlm .gt. nlm0) call rxi('increase nlm0 in hklbl need',nlm)

      pi = 4d0*datan(1d0)
      fpi = 4d0*pi
      lmax = ll(nlm)

C ... Shorten p, multiply by corresponding phase later
Cgetarg...info...           real(8) :: alat , real(8) :: slat%alat 1
Cgetarg...info...           real(8) :: plat(3,3) , real(8) :: slat%plat 9
Cgetarg...info...           real(8) :: qlat(3,3) , real(8) :: slat%qlat 9
Cgetarg       call upack('lat alat plat qlat',slat,alat,plat,qlat,0,0)

      alat=slat%alat
      i_copy_size=size(slat%plat)
      call dcopy(i_copy_size,slat%plat,1,plat,1)
      i_copy_size=size(slat%qlat)
      call dcopy(i_copy_size,slat%qlat,1,qlat,1)

      call shorbz(p,p1,plat,qlat)
      sp = 2*pi*(q(1)*(p(1)-p1(1))+q(2)*(p(2)-p1(2))+q(3)*(p(3)-p1(3)))
      phase = dcmplx(dcos(sp),dsin(sp))

C ... Use standard routines
C     call gklbl(p1,rsm,e,q,kmax-1,nlm,k0,cy,slat, hkl)
C     call hsmbl(p1,rsm,e,q,lmax,cy,slat, hsm,hsmp)

C ... Alternatively, use vectorized equivalents (about 2x faster)
Cgetarg...info...           real(8) :: awald , real(8) :: slat%awald 1
Cgetarg...info...           real(8) :: tol , real(8) :: slat%tol 1
Cgetarg...info...           real(8) :: vol , real(8) :: slat%vol 1
Cgetarg       call upack('lat awald tol vol',slat,awald,tol,vol,0,0)

      awald=slat%awald
      tol=slat%tol
      vol=slat%vol

Cgetarg...info...           integer :: nkd , integer(8) :: slat%nkd 1
Cgetarg...info...           integer :: nkq , integer(8) :: slat%nkq 1
Cgetarg...info...           integer :: odlv , integer(8) :: slat%odlv 1
Cgetarg...info...           integer :: oqlv , integer(8) :: slat%oqlv 1
Cgetarg       call upack('lat nkd nkq odlv oqlv',slat,nkd,nkq,odlv,oqlv,0)

      nkd=slat%nkd
      nkq=slat%nkq
      odlv=slat%odlv
      oqlv=slat%oqlv

      nrx = max(nkd,nkq)
#if F90
c      print *,'hklbl:aaa'
      allocate(wk(nrx*(2*lmax+10)),yl(nrx*(lmax+1)**2))
c      print *,'hklbl:bbb'
      call hsmq(1,0,ll(nlm),e,rsm,0000,q,p1,nrx,nlm,wk,yl,
     .  awald,alat,w(oqlv),nkq,w(odlv),nkd,vol,hsm,hsmp)
      if (rsm .gt. faca/awald) then
        call gklbl(p1,rsm,e,q,kmax-1,nlm,k0,cy,slat, hkl)
      else
        job = 2
        call gklq(lmax,rsm,q,p1,e,kmax-1,k0,alat,w(odlv),nkd,nrx,
     .    yl,wk,job,hkl)
      endif
      deallocate(wk,yl)
#else
      call defrr(owk, nrx*(2*lmax+10))
      call defrr(oyl, nrx*(lmax+1)**2)
      call hsmq(1,0,ll(nlm),e,rsm,0000,q,p1,nrx,nlm,w(owk),w(oyl),
     .  awald,alat,w(oqlv),nkq,w(odlv),nkd,vol,hsm,hsmp)
      if (rsm .gt. faca/awald) then
        call gklbl(p1,rsm,e,q,kmax-1,nlm,k0,cy,slat, hkl)
      else
        job = 2
        call gklq(lmax,rsm,q,p1,e,kmax-1,k0,alat,w(odlv),nkd,nrx,
     .    w(oyl),w(owk),job,hkl)
      endif
      call rlse(owk)
#endif

C  ... for debugging
C      do  10  ilm = 1, nlm
C        if (cdabs(hsm(ilm)-hsmx(ilm)) .gt. 1d-10) then
C          if (cdabs(hsm(ilm)/hsmx(ilm)-1) .gt. 1d-8) then
C     .      .and. cdabs(hsm(ilm)) .gt. .1d0
C     .      .and. cdabs(hsm(ilm)) .lt. 1.d0) then
C          print *, 'error at ilm=',ilm
C          print *, hsm(ilm)
C          print *, hsmx(ilm)
C          print *, '---'
C          stop
C        endif
C        endif
C   10 continue

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

Cgetarg       end
      end subroutine hklbl


