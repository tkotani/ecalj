#define F90 1
      subroutine fklbl(p,rsm,kmax,nlm,k0,cy,slat,fkl)

      use m_struc_def  !Cgetarg

C- Bloch sum of smooth Hankels for e=0 and q=(0,0,0).
C ----------------------------------------------------------------------
Ci Inputs
Ci   p     :Function is centered at p
Ci   rsm   :smoothing radius
Ci   kmax  :polynomial cutoff
Ci   nlm   :L-cutoff for gkl
Ci   k0    :leading dimension of gkl
Ci   cy    :Normalization constants for spherical harmonics
Ci   slat  :struct containing information about the lattice
Co Outputs
Co   fkl   :Bloch-summed Hankels for q=0 and e=0
Cr Remarks
Cr   For (k=0,l=0) f equals the limit of hklbl minus the avg value.
Cr   For all other cases f is the limit of hklbl as e goes to zero.
Cu Updates
Cu   24 Apr 00 Adapted from nfp fkl_bl.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer kmax,nlm,k0
      real(8):: p(3) , cy(1) , rsm
      type(s_lat)::slat

      double complex fkl(0:k0,nlm)
C ... Local parameters
      integer nlm0,lmax,ll,nkd,nkq,odlv,oqlv,nrx,owk,oyl,job,ilm,k
      parameter ( nlm0=196 )
      double precision q(3),alat,plat(9),qlat(9),p1(3)
      double precision faca,fpi,y0,e,vol,awald,tol
      double complex fsm(nlm0),gklsav,gklnew
      parameter (faca=1d0)
#if F90
      real(8),allocatable:: wk(:),yl(:)
#endif
C ... Heap
      integer w(1)
      common /w/ w

      if (nlm .eq. 0) return
      fpi = 16d0*datan(1d0)
      y0 = 1d0/dsqrt(fpi)
      if (nlm .gt. nlm0) call rx('fklbl: increase nlm0')
      lmax = ll(nlm)
      e = 0d0
      q(1) = 0d0
      q(2) = 0d0
      q(3) = 0d0

C ... Use standard routines
C     call gklbl(p,rsm,e,q,kmax-1,nlm,k0,cy,slat, fkl)
C     call fsmbl(p,rsm,lmax,cy,slat, fsm)

C ... Alternatively, use vectorized equivalents (about 2x faster)

      alat=slat%alat
      i_copy_size=size(slat%plat)
      call dcopy(i_copy_size,slat%plat,1,plat,1)
      i_copy_size=size(slat%qlat)
      call dcopy(i_copy_size,slat%qlat,1,qlat,1)


      awald=slat%awald
      tol=slat%tol
      vol=slat%vol


      nkd=slat%nkd
      nkq=slat%nkq
      odlv=slat%odlv
      oqlv=slat%oqlv

      call shorbz(p,p1,plat,qlat)
      nrx = max(nkd,nkq)
#if F90
      allocate(wk(nrx*(2*lmax+10)),yl(nrx*(lmax+1)**2))
      call hsmqe0(lmax,rsm,0,q,p1,nrx,nlm,wk,yl,
     .  awald,alat,w(oqlv),nkq,w(odlv),nkd,vol,fsm)
      if (rsm .gt. faca/awald) then
        call gklbl(p1,rsm,e,q,kmax-1,nlm,k0,cy,slat, fkl)
      else
        job = 2
        call gklq(lmax,rsm,q,p1,e,kmax-1,k0,alat,w(odlv),nkd,nrx,
     .    yl,wk,job,fkl)
      endif
      deallocate(wk,yl)
#else
      call defrr(owk, nrx*(2*lmax+10))
      call defrr(oyl, nrx*(lmax+1)**2)
      call hsmqe0(lmax,rsm,0,q,p1,nrx,nlm,w(owk),w(oyl),
     .  awald,alat,w(oqlv),nkq,w(odlv),nkd,vol,fsm)
      if (rsm .gt. faca/awald) then
        call gklbl(p1,rsm,e,q,kmax-1,nlm,k0,cy,slat, fkl)
      else
        job = 2
        call gklq(lmax,rsm,q,p1,e,kmax-1,k0,alat,w(odlv),nkd,nrx,
     .    w(oyl),w(owk),job,fkl)
      endif
      call rlse(owk)
#endif

C ... Upward recursion in k: mainly sets fkl = -4*pi * g(k-1,l)
      do  2  ilm = 1, nlm
        gklsav = fkl(0,ilm)
        fkl(0,ilm) = fsm(ilm)
      do  2  k = 1, kmax
        gklnew = fkl(k,ilm)
        fkl(k,ilm) = -fpi*gklsav
        gklsav = gklnew
  2   continue

C ... Add extra term to F(k=1,l=0)
      fkl(1,1) = fkl(1,1) + fpi*y0/vol

      end subroutine fklbl


