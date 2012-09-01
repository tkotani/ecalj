CSFCPP#define F90 1
      subroutine zhevx(n,lh,h,s,lov,lx,nmx,emx,nev,wk,linv,e,lz,z)
C- Eigenvalues and/or some eigenvectors of a Hermitian matrix
C ----------------------------------------------------------------
Ci Inputs:
Ci   n:    order of h and s
Ci   lh:   leading dimension of h and s
Ci   h:    hermitian matrix, dimensioned h(n,n)
Ci   s:    hermitian overlap matrix, (used only if lov is true)
Ci   nmx:  maximum number of eigenvectors to be found
Ci   emx:  eigenvalue limit for eigenvectors to be found
Ci         (not used if LAPACK zhegv is invoked)
Ci   wk:   work array of length at least 11n
Ci         NB: If LAPACK version is used, and eigenvectors are sought
Ci         wk should be dimensioned (n*nmx*2)
Ci   lov:  0 no overlap matrix
Ci         1 overlap matrix, return evecs of nonorthogonal H
Ci   lx:   if T, calls routines to exploit unit stride lengths (risc)
Ci         Not used if LAPACK zhegv is invoked.
Ci   linv: if T, using inverse iteration
Ci         Not used if LAPACK zhegv is invoked.
Ci   lz:   leading dimension of z
Co Outputs:
Co   e:    eigenvalues
Co   nev:  number of eigenvectors found
Co   z:    eigenvectors (1..nev)  (declared as z(n,*)
Co   s:    has been decomposed into and LL+ decomposition.
Co         You can call zhev2 to scale a vector by L
Cr Remarks:
Cr   z must be at least of dimension z(n,n), even though nev<n.
Cr   h and s are destroyed on exit.
Cr   Aborts on exit
Cp Procedures used:
Cp   (lapack)  zhegv
Cp   (eispack) htribk, htridx, imtql2, tqlrat
Cu Updates
Cu   24 Feb 07 Bug fix when nmx=0
Cu   17 May 03 Adapted from zhev, intended to supersede zhev.
Cu   14 Aug 02 Added zheev when lov is F; new zhev2.
Cu   21 Jan 02 Added code to invoke LAPACK zhegv in place of diagno
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      logical linv,lx
      integer lov,lh,lz
      integer n,nev,nmx
      double precision h(*),s(*),e(n),wk(*),z(2),emx

C$$$#if LAPACK
C Local parameters
      integer ier,lwork
      character jobz

      call tcn('zhev')
      if (nmx .le. 0) then
C       lh = n
        jobz = 'N'
        lwork = 4*n
        if (lov .gt. 0) then
          call zhegv(1,jobz,'U',n,h,lh,s,lh,e,wk(1+3*n),lwork,wk(1),ier)
          call rxx(ier.ne.0,'zhev: zhegv cannot find all evals')
        else
          call zheev(jobz,'U',n,h,lh,e,wk(1+3*n),lwork,wk(1),ier)
          call rxx(ier.ne.0,'zhev: zheev cannot find all evals')
        endif
        nev = 0
      else
        jobz = 'V'
        lwork = n*nmx
        if (lov .gt. 0) then
          call zhegv(1,jobz,'U',n,h,lh,s,lh,e,z,lwork,wk(1),ier)
          call rxx(ier.ne.0,'zhev: zhegv cannot find all evals')
        else
          call zheev(jobz,'U',n,h,lh,e,z,lwork,wk(1),ier)
          call rxx(ier.ne.0,'zhev: zheev cannot find all evals')
        endif
        call zmcpy('N',h,lh,1,z,lz,1,n,min(n,nmx))
C       call zprm('evecs',2,z,lz,n,nmx)
        nev = min(n,nmx)
      endif
C$$$#else
C$$$C Local parameters
C$$$      integer i,j,k,m,mdim,k1,k2,k3,k4,k5,k6,k7,k8,k9,n2,ier

C$$$      call tcn('zhev')
C$$$      n2 = n**2

C$$$C --- Take care of n=1 case ---
C$$$      if (n .eq. 1) then
C$$$        e(1) = h(1)
C$$$        if (lov .gt. 0) e(1) = h(1)/s(1)
C$$$        if (1 .gt. nev .or. e(1) .gt. emx) return
C$$$        z(1) = 1
C$$$        z(2) = 0
C$$$        return
C$$$      endif

C$$$C --- Separate real and imaginary parts ---
C$$$      mdim = 2*lh
C$$$      call ztoy(h,lh,n,n,0)
C$$$      if (lov .gt. 0) call ztoy(s,lh,n,n,0)

C$$$C --- Debugging: eigenvalues of overlap ---
C$$$C      call htridx(mdim,n,s(1),s(lh+1),e,wk,wk(3*n+1),wk(n+1))
C$$$C      do  11  j = 1, n
C$$$C   11 wk(j) = wk(j)**2
C$$$C      call tqlrat(n,e,wk,ier)
C$$$C      write(6,600) e
C$$$C  600 format(' evl='/(1x,1p,5e14.6))
C$$$C      call rx('eigenvalues of overlap')

C$$$C --- H <- S^-1/2  H  S^-1/2 ---
C$$$      if (lov .gt. 0) then
C$$$        call yyhchd(mdim,n,s,s(lh+1),wk,lx,.true.,ier)
C$$$        call rxx(ier.ne.0,'ZHEV: error in yyhchd')
C$$$        if (lx) then
C$$$          if (lz .ne. lh) call rx('zhev not ready for lz ne lh')
C$$$          call yyhrdx(mdim,n,h,h(lh+1),s,s(lh+1),z,z(lz+1))
C$$$        else
C$$$          call yyhred(mdim,n,h,h(lh+1),s,s(lh+1),.true.)
C$$$        endif
C$$$      endif

C$$$C --- Transform to tridiagonal matrix ---
C$$$      if (linv) then
C$$$        k1 = 1
C$$$        k2 = k1 + 3*n
C$$$        k3 = k2 + n
C$$$        k4 = k3 + n
C$$$        k5 = k4 + n
C$$$        k6 = k5 + n
C$$$        k7 = k6 + n
C$$$        k8 = k7 + n
C$$$        k9 = k8 + n
C$$$#ifndef GENERIC
C$$$        call htridx(mdim,n,h(1),h(lh+1),wk(k1),wk(k2),wk(k3),wk(n+1))
C$$$#else
C$$$        call htridi(mdim,n,h(1),h(lh+1),wk(k1),wk(k2),wk(k3),wk(n+1))
C$$$#endif
C$$$      else
C$$$#ifndef GENERIC
C$$$        call htridx(mdim,n,h(1),h(lh+1),e,wk,wk(3*n+1),wk(n+1))
C$$$#else
C$$$        call htridi(mdim,n,h(1),h(lh+1),e,wk,wk(3*n+1),wk(n+1))
C$$$#endif
C$$$      endif

C$$$C --- Eigenvalues only ---
C$$$      if (nmx .le. 0) then
C$$$        do  12  j = 1, n
C$$$   12   wk(j) = wk(j)**2
C$$$        call tqlrat(n,e,wk,ier)
C$$$        call rxx(ier.ne.0,'ZHEV: tqlrat cannot find all evals')
C$$$        nev = 0
C$$$        goto 100

C$$$C --- Eigenvalues and eigenvectors ---
C$$$      else if (linv) then
C$$$        call imtqlv(n,wk(k1),wk(k2),wk(k3),e,wk(k9),ier,wk(k4))
C$$$        call rxx(ier.ne.0,'ZHEV: imtqlv cannot find all evals')
C$$$C   --- Determine number of eigenvectors to be calculated ---
C$$$        nev = 1
C$$$        do  14  j = 2, n
C$$$          if (j .le. nmx .and. e(j-1) .le. emx) nev = j
C$$$   14   continue
C$$$        call tinvit(2*lz,n,wk(k1),wk(k2),wk(k3),nev,e,wk(k9),z,ier,
C$$$     .  wk(k4),wk(k5),wk(k6),wk(k7),wk(k8))
C$$$        call rxx(ier.ne.0,'ZHEV: tinvit cannot find all evecs')
C$$$      else
C$$$        do  17  j = 1, n
C$$$          k = (j-1)*2*lz
C$$$          m = k+n
C$$$          do  16  i = k+1, m
C$$$   16     z(i) = 0d0
C$$$          z(k+j) = 1d0
C$$$   17   continue
C$$$        call imtql2(2*lz,n,e,wk,z,ier)
C$$$C       call prmx('eval',e,n,n,1)
C$$$        call rxx(ier.ne.0,'ZHEV: imtql2 cannot find all evecs')

C$$$C   --- Determine number of eigenvectors to be calculated ---
C$$$        nev = 1
C$$$        do  15  j = 2, n
C$$$          if (j .le. nmx .and. e(j-1) .le. emx) nev = j
C$$$   15   continue
C$$$      endif
C$$$C     call prmx('eval',e,n,n,1)

C$$$      if (nev .gt. 0) then
C$$$#ifndef GENERIC
C$$$        if (lz .ne. lh) call rx('zhev not ready for lz ne lh')
C$$$        call htribx(mdim,n,h(1),h(lh+1),wk(n+1),nev,z(1),z(lz+1))
C$$$#else
C$$$        call htribk(mdim,n,h(1),h(lh+1),wk(n+1),nev,z(1),z(lz+1))
C$$$#endif

C$$$C --- Get the eigenvectors of H - E O ---
C$$$        if (lov .gt. 0) then
C$$$          if (lx) then
C$$$            call ymcpy(z,2*lz,1,lz,h,mdim,1,lh,n,n)
C$$$            call yympy(s,s(lh+1),mdim,1,h,h(lh+1),mdim,1,z,z(lz+1),2*lz,
C$$$     .      1,n,nev,n)
C$$$          else
C$$$            call yyhbak(mdim,n,s,s(lh+1),nev,z,z(lz+1),.true.)
C$$$          endif
C$$$        endif

C$$$C   --- Convert eigenvectors to double complex storage ---
C$$$        call ztoy(z,lz,n,nev,1)
C$$$C       call zprm('evecs',2,z,lz,n,nmx)
C$$$      endif

C$$$#endif

  100 call tcx('zhev')
      end

      subroutine zhev(n,h,s,lov,lx,nmx,emx,nev,wk,linv,ltime,e,z)
C- Eigenvalues and/or some eigenvectors of a Hermitian matrix
C ----------------------------------------------------------------
Ci Inputs:
Ci   n:    dimension of h
Ci   h,n:  hermitian matrix, dimensioned h(n,n)
Ci   s:    hermitian overlap matrix, (used only if lov is true)
Ci   nmx:  maximum number of eigenvectors to be found
Ci   emx:  eigenvalue limit for eigenvectors to be found
Ci   wk:   work array of length at least 11n
Ci   lov:  if T, non-orthogonal
Ci   lx:   if T, calls routines to exploit unit stride lengths (risc)
Ci         Not used if LAPACK zhegv is invoked.
Ci   linv: if T, using inverse iteration
Ci         Not used if LAPACK zhegv is invoked.
Co Outputs:
Co   e:    eigenvalues
Co   nev:  number of eigenvectors found
Co   z:    eigenvectors (1..nev)  (declared as z(n,*)
Co   s:    has been decomposed into and LL+ decomposition.
Co         You can call zhev2 to scale a vector by L
Cr Remarks:
Cr   z must be at least of dimension z(n,n), even though nev<n.
Cr   h and s are destroyed on exit.
Cr   Aborts on exit
Cp Procedures used:
Cp   (lapack)  zhegv
Cp   (eispack) htribk, htridx, imtql2, tqlrat
Cu Updates
Cu   14 Aug 02 Added zheev when lov is F; new zhev2.
Cu   21 Jan 02 Added code to invoke LAPACK zhegv in place of diagno
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      logical linv,lx
      integer n,nev,nmx,ltime
      double precision h(*),s(*),e(n),wk(*),z(*),emx
      logical lov

C$$$#if LAPACK
C Local parameters
      integer ier,lwork
      character jobz

      call tcn('zhev')
      if (nmx .le. 0) then
        jobz = 'N'
        lwork = 4*n
        if (lov) then
          call zhegv(1,jobz,'U',n,h,n,s,n,e,wk(1+3*n),lwork,wk(1),ier)
        else
          call zheev(jobz,'U',n,h,n,e,wk(1+3*n),lwork,wk(1),ier)
        endif
        nev = 0
      else
        jobz = 'V'
        lwork = n*min(n,nmx)
        if (lov) then
          call zhegv(1,jobz,'U',n,h,n,s,n,e,z,lwork,wk(1),ier)
        else
          call zheev(jobz,'U',n,h,n,e,z,lwork,wk(1),ier)
        endif
        call zcopy(n*min(n,nmx),h,1,z,1)
        nev = min(n,nmx)
      endif
      call rxx(ier.ne.0,'zhev: zhegv cannot find all evals')
C$$$#else
C$$$C Local parameters
C$$$      integer i,j,k,m,mdim,k1,k2,k3,k4,k5,k6,k7,k8,k9,n2,ier

C$$$C     call prm('(5f20.15)',s,n,n)
C$$$      call tcn('zhev')
C$$$      n2 = n**2

C$$$C --- Take care of n=1 case ---
C$$$      if (n .eq. 1) then
C$$$        e(1) = h(1)
C$$$        if (lov) e(1) = h(1)/s(1)
C$$$        if (1 .gt. nev .or. e(1) .gt. emx) return
C$$$        z(1) = 1
C$$$        z(2) = 0
C$$$        return
C$$$      endif

C$$$CC --- Use diagno if lx ---
C$$$C      if (lx) then
C$$$C        call cplx2r(n2,0,h,z)
C$$$C        if (lov) call cplx2r(n2,0,s,z)
C$$$C        i = 0
C$$$C        if (lov) i = 1
C$$$C        call diagno(n,h,s,wk,lx,i,linv,nmx,emx,nev,z,e)
C$$$C        call cplx2r(n2,1,z,h)
C$$$C        goto 40
C$$$C      endif

C$$$C --- Separate real and imaginary parts ---
C$$$      mdim = 2*n
C$$$      call ztoy(h,n,n,n,0)
C$$$      if (lov) call ztoy(s,n,n,n,0)
C$$$C      mdim = 2*n
C$$$C      do  10  j = 1, n
C$$$C        k = (j-1)*mdim + 1
C$$$C        l = k + n
C$$$C        call dcopy(n,h(k+1),2,wk,1)
C$$$C        call dcopy(n,h(k),2,h(k),1)
C$$$C        call dcopy(n,wk,1,h(l),1)
C$$$C        if (lov) then
C$$$C          call dcopy(n,s(k+1),2,wk,1)
C$$$C          call dcopy(n,s(k),2,s(k),1)
C$$$C          call dcopy(n,wk,1,s(l),1)
C$$$C        endif
C$$$C   10 continue


C$$$C --- Debugging: eigenvalues of overlap ---
C$$$C      call htridx(mdim,n,s(1),s(n+1),e,wk,wk(3*n+1),wk(n+1))
C$$$C      do  11  j = 1, n
C$$$C   11 wk(j) = wk(j)**2
C$$$C      call tqlrat(n,e,wk,ier)
C$$$C      write(6,600) e
C$$$C  600 format(' evl='/(1x,1p,5e14.6))
C$$$C      call rx('eigenvalues of overlap')

C$$$C --- H <- S^-1/2  H  S^-1/2 ---
C$$$      if (lov) then
C$$$        call yyhchd(mdim,n,s,s(n+1),wk,lx,.true.,ier)
C$$$        call rxx(ier.ne.0,'ZHEV: error in yyhchd')
C$$$        if (lx) then
C$$$          call yyhrdx(mdim,n,h,h(n+1),s,s(n+1),z,z(n+1))
C$$$        else
C$$$          call yyhred(mdim,n,h,h(n+1),s,s(n+1),.true.)
C$$$        endif
C$$$      endif

C$$$C --- Transform to tridiagonal matrix ---
C$$$      if (linv) then
C$$$        k1 = 1
C$$$        k2 = k1 + 3*n
C$$$        k3 = k2 + n
C$$$        k4 = k3 + n
C$$$        k5 = k4 + n
C$$$        k6 = k5 + n
C$$$        k7 = k6 + n
C$$$        k8 = k7 + n
C$$$        k9 = k8 + n
C$$$#ifndef GENERIC
C$$$        call htridx(mdim,n,h(1),h(n+1),wk(k1),wk(k2),wk(k3),wk(n+1))
C$$$#else
C$$$        call htridi(mdim,n,h(1),h(n+1),wk(k1),wk(k2),wk(k3),wk(n+1))
C$$$#endif
C$$$      else
C$$$#ifndef GENERIC
C$$$        call htridx(mdim,n,h(1),h(n+1),e,wk,wk(3*n+1),wk(n+1))
C$$$#else
C$$$        call htridi(mdim,n,h(1),h(n+1),e,wk,wk(3*n+1),wk(n+1))
C$$$#endif
C$$$      endif

C$$$C --- Eigenvalues only ---
C$$$      if (nmx .le. 0) then
C$$$        do  12  j = 1, n
C$$$   12   wk(j) = wk(j)**2
C$$$        call tqlrat(n,e,wk,ier)
C$$$        call rxx(ier.ne.0,'ZHEV: tqlrat cannot find all evals')
C$$$        nev = 0
C$$$        goto 100

C$$$C --- Eigenvalues and eigenvectors ---
C$$$      else if (linv) then
C$$$        call imtqlv(n,wk(k1),wk(k2),wk(k3),e,wk(k9),ier,wk(k4))
C$$$        call rxx(ier.ne.0,'ZHEV: imtqlv cannot find all evals')
C$$$C   --- Determine number of eigenvectors to be calculated ---
C$$$        nev = 1
C$$$        do  14  j = 2, n
C$$$          if (j .le. nmx .and. e(j-1) .le. emx) nev = j
C$$$   14   continue
C$$$        call tinvit(mdim,n,wk(k1),wk(k2),wk(k3),nev,e,wk(k9),z,ier,
C$$$     .  wk(k4),wk(k5),wk(k6),wk(k7),wk(k8))
C$$$        call rxx(ier.ne.0,'ZHEV: tinvit cannot find all evecs')
C$$$      else
C$$$        do  17  j = 1, n
C$$$          k = (j-1)*mdim
C$$$          m = k+n
C$$$          do  16  i = k+1, m
C$$$   16     z(i) = 0d0
C$$$          z(k+j) = 1d0
C$$$   17   continue
C$$$        call imtql2(mdim,n,e,wk,z,ier)
C$$$        call rxx(ier.ne.0,'ZHEV: imtql2 cannot find all evecs')

C$$$C   --- Determine number of eigenvectors to be calculated ---
C$$$        nev = 1
C$$$        do  15  j = 2, n
C$$$          if (j .le. nmx .and. e(j-1) .le. emx) nev = j
C$$$   15   continue
C$$$      endif

C$$$      if (nev .gt. 0) then
C$$$#ifndef GENERIC
C$$$        call htribx(mdim,n,h(1),h(n+1),wk(n+1),nev,z(1),z(n+1))
C$$$#else
C$$$        call htribk(mdim,n,h(1),h(n+1),wk(n+1),nev,z(1),z(n+1))
C$$$#endif

C$$$C --- Get the eigenvectors of H - E O ---
C$$$        if (lov) then
C$$$          if (lx) then
C$$$            call dcopy(n2*2,z,1,h,1)
C$$$            call yympy(s,s(n+1),mdim,1,h,h(n+1),mdim,1,z,z(n+1),mdim,1,
C$$$     .      n,nev,n)
C$$$          else
C$$$            call yyhbak(mdim,n,s,s(n+1),nev,z,z(n+1),.true.)
C$$$          endif
C$$$        endif

C$$$C   --- Convert eigenvectors to double complex storage ---
C$$$        call ztoy(z,n,n,nev,1)
C$$$C       call zprm('evecs',2,z,lz,n,min(nev,nmx))
C$$$      endif

C$$$C  40 if (nev .gt. 0 .and. ltime .ge. 0) print 337, nev,nmx,n,emx
C$$$C  337 format(' nev, nevmx, ndim=',3i4,' emx=',f10.5)
C$$$#endif

  100 call tcx('zhev')
      end

      subroutine zhev2(mode,ldhs,n,m,lx,h,s)
C- Scales a hermitian matrix h by Cholesky decomposition of overlap s
C ----------------------------------------------------------------
Ci Inputs:
Ci   mode :For all modes, see Remarks.
Ci        :0 do nothing to h
Ci        :2  replace h to orthogonal form L^-1 h L+^-1 or U+^-1 h U^-1
Ci        :   h is assumed to be hermitian
Ci        :4  replace h to orthogonal by L h L+ or U+ h U
Ci        :   Modes 2 and 4 perform inverse functions of each other
Ci        :10 replace h to U^-1 h
Ci        :20 replace h to U h
Ci        :   Modes 10 and 20 perform inverse functions of each other
Ci        :
Ci        :Adding 1 to any of the following modes causes zhev2 to
Ci        :Cholesky-decompose hermitian matrix s.
Ci        :h is not used in the decomposition of s.
Ci        :
Ci   n    :dimension of h and s
Ci   m    :second dimension of h (mode>=10)
Ci   lx   :if T, calls custom routines for Cholesky decomp'sn
Ci        :that exploit unit stride lengths
Ci        :(not used in LAPACK implementation)
Co Inputs/Outputs:
Cio  s    :Overlap matrix.  s should be in Cholesky-decomposed form as
Cio       :computed by zhev or zhev2; or else 1 should be added to mode
Cio       :to decompose s.  s is unchanged apart from possible decomps'n
Cio  h    :On input:
Cio       :  h is a hermitian matrix of dimension (n,n) for mode<10
Cio       :  h is a matrix of dimension (n,m) for mode>=10
Cio       :On output, h is scaled depending on mode (or mode-1 if s is C.D.)
Cio       :  mode           h is transformed into:
Cio       :   0             h is unchanged
Cio       :   2             L^-1 h L+^-1 or U+^-1 h U^-1
Cio       :   4             L h L+       or U+ h U
Cio       :  10             U^-1 h
Cr Remarks:
Cr   zhev2 performs one of several operations that deal with overlap
Cr   matrices in a generalized eigenvalue problem.  A hermitian overlap
Cr   matrix S can be decomposed into L L+ (LAPACK implementation: U+ U)
Cr   Then the eigenvalues e of (H - e S) are the eigenvalues of
Cr     H' = L^-1 H L+^-1 (H' = U+^-1 H U^-1)
Cr   and eigenvectors Z of (H - e S) are related to vectors Z' of H' as
Cr     Z = L+^-1 Z' (Z = U^-1 Z')
Cu Updates
Cu   13 Aug 02 first created
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer mode,ldhs,n,m
      double precision h(*),s(*)
      logical lx
C Local parameters
      integer ier,itype,lmode
      double complex one
      parameter (one=(1.0d+0,0.0d+0))

      call tcn('zhev2')

      lmode = mode

C$$$#ifndef LAPACK
C$$$      call rx('zhev2:  only implemented in LAPACK mode, sorry')
C$$$#endif

C --- Cholesky decompose s ---
      if (mod(lmode,2) .ne. 0) then
        call zpotrf('U',n,s,ldhs,ier)
        if (ier .ne. 0) then
          call rxi('zhev2: failed to c.d. s: ier = ',ier)
        endif
        lmode = mode-1
      endif

C --- do nothing to h ---
      if (lmode .eq. 0) then
        continue

C --- h <- U+^-1 h U^-1 ---
      elseif (lmode .eq. 2) then
        itype = 1
        call zhegst(itype,'U',n,h,ldhs,s,ldhs,ier)
        call rxx(ier.ne.0,'zhev: zhegv cannot find all evals')
        call z2herm('U',ldhs,n,h)

C --- h <- U+ h U ---
      elseif (lmode .eq. 4) then
C       Overwrite h with U+ h
        call ztrmm('L','U','C','N',n,n,one,s,ldhs,h,ldhs)
C       Overwrite U+ h with U+ h U
        call ztrmm('R','U','N','N',n,n,one,s,ldhs,h,ldhs)
C       print '(a)', '#U+ h U'
C       call zprm('(5f12.6)',h,n,n)

C --- h <- U^-1 h ---
      elseif (lmode .eq. 10) then
        call ztrsm('L','U','N','N',n,m,one,s,ldhs,h,ldhs)

C --- h <- U h ---
      elseif (lmode .eq. 20) then
        call ztrmm('L','U','N','N',n,m,one,s,ldhs,h,ldhs)

      else
        call rxi('zhev2: bad mode',mode)
      endif


      call tcx('zhev2')
      end

      subroutine zhev3(mode,n,h,s,z)
C- Scales a hermitian matrix h by evec decomposition of overlap s
C ----------------------------------------------------------------
Ci Inputs:
Ci   mode :For all modes, see Remarks.
Ci        :0 do nothing to h
Ci        :2  replace h to orthogonal form L^-1 h L+^-1 or U+^-1 h U^-1
Ci        :   h is assumed to be hermitian
Ci        :4  replace h to orthogonal by L h L+ or U+ h U
Ci        :   Modes 2 and 4 perform inverse functions of each other
Ci        :10 replace h to U^-1 h
Ci        :20 replace h to U h
Ci        :   Modes 10 and 20 perform inverse functions of each other
Ci        :
Ci        :Adding 1 to any of the following modes causes zhev3 to
Ci        :Cholesky-decompose hermitian matrix s.
Ci        :h is not used in the decomposition of s.
Ci        :
Ci   n    :dimension of h and s
Ci   m    :second dimension of h (mode>=10)
Ci   lx   :if T, calls custom routines for Cholesky decomp'sn
Ci        :that exploit unit stride lengths
Ci        :(not used in LAPACK implementation)
Co Inputs/Outputs:
Cio  s    :Overlap matrix.  s should be in Cholesky-decomposed form as
Cio       :computed by zhev or zhev3; or else 1 should be added to mode
Cio       :to decompose s.  s is unchanged apart from possible decomps'n
Cio  h    :On input:
Cio       :  h is a hermitian matrix of dimension (n,n) for mode<10
Cio       :  h is a matrix of dimension (n,m) for mode>=10
Cio       :On output, h is scaled depending on mode (or mode-1 if s is C.D.)
Cio       :  mode           h is transformed into:
Cio       :   0             h is unchanged
Cio       :   2             L^-1 h L+^-1 or U+^-1 h U^-1
Cio       :   4             L h L+       or U+ h U
Cio       :  10             U^-1 h
Cr Remarks:
Cr   zhev3 performs one of several operations that deal with overlap
Cr   matrices in a generalized eigenvalue problem.  A hermitian overlap
Cr   matrix S can be decomposed into L L+ (LAPACK implementation: U+ U)
Cr   Then the eigenvalues e of (H - e S) are the eigenvalues of
Cr     H' = L^-1 H L+^-1 (H' = U+^-1 H U^-1)
Cr   and eigenvectors Z of (H - e S) are related to vectors Z' of H' as
Cr     Z = L+^-1 Z' (Z = U^-1 Z')
Cu Updates
Cu   13 Aug 02 first created
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer mode,n
      double complex h(n,n),s(n,n),z(n,n)
C Local parameters
CSFCPP#if F90
      real(8),allocatable:: e(:),wk(:)
CSFCPP#else
CSFCPP      double precision e(n),wk(11*n)
CSFCPP#endif
      logical lx
      integer i,lmode
      double complex zer,one
      parameter (zer=(0d0,0d0),one=(1d0,0d0))
      parameter(lx=.true.)

      call tcn('zhev3')

C$$$#ifndef LAPACK
C$$$      call rx('zhev3:  only implemented in LAPACK mode, sorry')
C$$$#endif

      lmode = mode

C --- Eigenvectors and eigenvalues of s ---
      if (mod(mode,2) .ne. 0) then
CSFCPP#if F90
        allocate(e(n),wk(11*n))
CSFCPP#endif
        call zhev(n,s,s,.false.,lx,n,9d9,i,wk,.false.,0,e,z)
C      if (lmode .eq. 2) then
C      do  i = 1, n
C        do  j = 1, n
C          z(i,j) = z(i,j) / sqrt(e(j))
C        enddo
C      enddo
C      endif
C      if (lmode .eq. 4) then
C      do  i = 1, n
C        do  j = 1, n
C          z(i,j) = sqrt(e(i)) * z(i,j)
C        enddo
C      enddo
C      endif
CSFCPP#if F90
        deallocate(e,wk)
CSFCPP#endif
        lmode = mode-1
      endif

C --- do nothing to h ---
      if (lmode .eq. 0) then
        continue

C --- h <- z+ h z ---
      elseif (lmode .eq. 2) then
        call zgemm('N','N',n,n,n,one,h,n,z,n,zer,s,n)
        call zgemm('C','N',n,n,n,one,z,n,s,n,zer,h,n)
C       call zprm('z+ h z',2,h,n,n,n)

C --- h <- z h z+ ---
      elseif (lmode .eq. 4) then
        call zgemm('N','C',n,n,n,one,h,n,z,n,zer,s,n)
        call zgemm('N','N',n,n,n,one,z,n,s,n,zer,h,n)
C       call zprm('z h z+',2,h,n,n,n)

      else
        call rxi('zhev3: bad mode',mode)

      endif

      call tcx('zhev3')
      end

      subroutine zhevo(n,lh,h,s,lx,nmx,emx,epsovl,nevl,nev,e,eo,lz,z)
C- Eigenvalues and/or some eigenvectors of a Hermitian matrix in
C- reduced hilbert space, eliminating projections onto small overlap
C ----------------------------------------------------------------
Ci Inputs:
Ci   n:    order of h and s
Ci   lh:   leading dimension of h and s
Ci   h:    hermitian matrix, dimensioned h(n,n)
Ci   s:    hermitian overlap matrix
Ci   lx:   if T, calls routines to exploit unit stride lengths (risc)
Ci         Not used if diagonalization through LAPACK
Ci   nmx:  maximum number of eigenvectors to be found
Ci   emx:  eigenvalue limit for eigenvectors to be found
Ci         (not used if LAPACK zhegv is invoked)
Ci  epsovl:hilbert space belonging to evals of s < epsovl is discarded
Ci   lz:   leading dimension of z
Co Outputs:
Co   nevl: reduced dimension of hilbert space, and number of evals found
Co   nev:  number of eigenvectors found
Co   e:    eigenvalues
Co   eo:   eigenvalues of overlap matrix
Co   z:    eigenvectors (1..nev)
Cr Remarks:
Cr   h and s are destroyed on exit.
Cp Procedures used:
Cp   (lapack) zheev
Cp   (blas)   zgemm
Cb Bugs
Cb   Probably possible to improve memory usage.
Cb   Replace call to zheev with call to zhevx
Cu Updates
Cu   8  Jul 08 (T. Kotani) first implementation, adapted from zhev.
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      logical lx
      integer lh,lz
      integer n,nevl,nev,nmx
      double complex h(lh,n),s(lh,n),z(lz,n)
      double precision e(n),eo(n),emx,epsovl
CSFCPP#if F90
C Local parameters
      integer ier,lwork,i,ni
      character jobz
      complex(8),allocatable:: wk(:),hh(:,:)
      real(8),allocatable:: rwork(:)
      double complex zer,one
      parameter (zer=(0d0,0d0),one=(1d0,0d0))
C external
      integer iprint
      call tcn('zhevo')

C ... Eigenvalues of s
      jobz = 'V'
      lwork = n*n
      allocate(wk(lwork),rwork(max(1,3*n-2)))
      call zheev(jobz,'U',n,s,lh,eo,wk,lwork,rwork,ier)
      deallocate(wk,rwork)
C ... Find 1st eigenvalue ni of S exceeeding threshold
      do  i = 1, n
        if (eo(i) .gt. epsovl) then
          ni = i
          exit
        endif
      enddo

C ... Map s to reduced Hilbert space (dimension = nevl)
C     For now, use temporary zz.  Better: use s in place of zz.
      nevl = n-ni+1
C ... Construct projection matrix
      allocate(hh(n,nevl))
      do  i = ni, n
        hh(:,i-ni+1) = s(:,i)/sqrt(eo(i))
      enddo
      do  i = 1, n-ni+1
        s(:,i) = hh(:,i)
      enddo
      deallocate(hh)

C ... Hamiltonian  <zz | H | zz>
      allocate(hh(nevl,n))
C     call zprm3('h',2,h,lh,n,n)
C     call zprm3('zz',2,s,lh,n,nevl)
      call zgemm('C','N',nevl,n,n,one,s,lh,h,lh,zer,hh,nevl)
C     call zprm3('zz+*h',2,hh,nevl,nevl,n)
      call zgemm('N','N',nevl,nevl,n,one,hh,nevl,s,lh,zer,h,lh)
C     call zprm3('zz+ * h * zz',2,h,lh,nevl,nevl)

ctakao
      if(iprint()>41) then
        write(6,"(a)")' zhevo: Overlap matrix'
        write(6,"(5(i5,d10.2))") (i,eo(i),i=1,ni-1)
        write(6,*) '   --- '
        write(6,"(5(i5,d10.2))") (i,eo(i),i=ni,n)
      endif
c      write(6,"(5(a,3i5)") nevl,nmx,nev

C ... Diagonalize
      if (nmx .le. 0) then
        jobz = 'N'
        lwork = 4*nevl
        allocate(wk(lwork),rwork(max(1,3*nevl-2)))
        call zheev(jobz,'U',nevl,h,lh,e,wk,lwork,rwork,ier)
        nev = 0
      else
        jobz = 'V'
        nev = min(nmx,nevl)
        lwork = nevl*4 ! nevl*nev ---> bug when Hydrogen case with nev=1
c        print *,'nnnn lwork=',lwork,nevl,nev
        allocate(wk(lwork),rwork(max(1,3*nevl-2)))
        call zheev(jobz,'U',nevl,h,lh,e,wk,lwork,rwork,ier)
C       Rotate back to nonorthogonal evecs including overlap matrix
        call zgemm('N','N',n,nev,nevl,one,s,lh,h,lh,zer,z,lz)
C       call zprm3('z',2,z,lz,n,nev)
      endif
      deallocate(wk,rwork,hh)

      call tcx('zhevo')

CSFCPP#else
CSFCPP      call rx('zhevo only implemented with F90')
CSFCPP#endif

      end





ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C we assume LAPACK and F90
      subroutine zhev_tk(nlmto,n,h,s,nmx, epsovl,ncut,
     &  nevl,nev,e,eo,lz,z,ipr,delta_stabilize)
C- Eigenvalues and/or some eigenvectors of a Hermitian matrix (throw away poor linear dependency)
C ----------------------------------------------------------------
Ci Inputs:
Ci   nlmto: not used (dimension of MTO)
Ci   n:    dimension of h
Ci   h,n:  hermitian matrix, dimensioned h(n,n)
Ci   s:    hermitian overlap matrix, (used only if lov is true)
Ci   nmx:  maximum number of eigenvectors to be found
Ci   emx:  eigenvalue limit for eigenvectors to be found
Ci   wk:   work array of length at least 11n
Ci   epsovl: threshold to throw away linear dependence of the overlap matrix s.
Ci   ipr :print switch
Co Outputs:
Co   e:    eigenvalues
Co   nev:  number of eigenvectors found
Co   z:    eigenvectors (1..nev)  (declared as z(n,*)
Co   s:    has been decomposed into and LL+ decomposition.
Co         You can call zhev2 to scale a vector by L
Cr Remarks:
Cr   z must be at least of dimension z(n,n), even though nev<n.
Cr   h and s are destroyed on exit.
Cr   Aborts on exit
Cp Procedures used:
Cp   (lapack)  zhegv
Cr
Cr This routine throw away basis set which are poorly linear dependent on others.
Cr We first find the Hilbert space which is linearly independent.
Cr Then solve <zz|H|zz> of eigenvalue problem, where |zz> spans the Hilbert space.
Cr It eigenfucntion is moved back to the original Hamiltonian's eigenfunction.
Cr So this can be used just as a replacement of zhev.
Cr
Cu Updates
Cu   8  Jul 08 kotani renewed it. Linear-dependency removal
Cu   This is originally from Mark'z zhev.
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer n,nev,nmx,ltime,ngv,lz,nevl,ncut
      logical ipr 
      double precision e(n),wk(1),emx,eee,fac
      complex(8):: h(n,n),s(n,n),z(lz,n)
      complex(8),allocatable:: omat(:,:),imat(:,:),wk11(:),
     &zz(:,:),hh(:,:),omat2(:,:),pmat(:,:),omat2i(:,:),oo(:,:),
     &hhm(:,:),zzz(:),omatinv(:,:),oxx(:,:)
#ifndef USE_TRANSPOSE
      complex(8),allocatable:: h_t(:,:)
#endif
      integer:: i,j,nlmto,ni,nm,ik,ik2,k
      real(8):: omax_mto, omax_g,ex ,  ddd2,ddd,epsovl,epsx
      real(8),allocatable:: scale(:),rwork(:)
      real(8):: eo(n) ,delta_stabilize

      integer :: epsovl_mode=1

C Local parameters
      integer ier,lwork,ix
      character*1 jobz
c
      call tcn('zhev_tk')

C --- mode 1. Remove poor linear dependency. All basis are treated on the same footing.
C ... rescale omat
c      allocate(omat2(n,n))
c      omat2=s
      allocate(scale(n))
      do i=1,n
        scale(i)= 1d0/sqrt(s(i,i))  
      enddo

c This scaling is to keep the Hilbert space of MTO.
      scale(1:nlmto)=scale(1:nlmto)*100d0
      do i=1,n
        do j=1,n
c        omat(i,j)=scale(i)*scale(j)*omat(i,j) !this is based on  1d0/scale*|basis>
          s(i,j)= scale(i)*scale(j)*s(i,j)
          h(i,j)= scale(i)*scale(j)*h(i,j) !this is to isolate the degenerated space.
        enddo
      enddo

!! ... eigenvalue of s (omat)
      allocate(omat(n,n))
      omat = s !reserved
      jobz = 'V'
      lwork = n*n
      allocate(wk11(lwork),rwork(max(1,3*n-2)))
      call zheev(jobz,'U',n,omat,n,eo,wk11,lwork,rwork,ier)
      deallocate(wk11,rwork)
      if(ipr) then
        write(6,*)'zhev_tk: ovlmat='
        do ix=0,n,5
          if(eo(ix)>0.1d0) then
           write(6,*) '  ... skip larger eigenvalues ...'
           exit
          endif
        write(6,"(5(i5,d10.2))") (i,eo(i),i=ix+1,ix+5)
        enddo  
      endif
c      write(6,*)' zheev: diagonalization of overlap matrix info=',ier



!! === new mode to avoid linear depedency by replacing h and o ===
!! h= h + 99*2**(-o/epsovl)
!! o= o +    2**(-o/epsovl) ! note overlap matrix o is by s in this code.
!!                          ! oox -> 2**(-o/expovl) 
c      print *,' ier=',ier
c      print *,'epsovl=',epsovl
c      if(epsovl_mode==1) then
c      if(epsovl_mode==1.and.epsovl/=0d0) then
cc       print *,' zhev:epsovl_mode=1 epsovl=',epsovl
c      if(epsovl>1.or.epsovl<0) call rx('HAM_OVEPS: NOT in 0<OVEPS<1')
      if(delta_stabilize>0d0) then
       allocate(oxx(n,n))
       oxx=0d0
       do k=1,n
c         fac = 2**(-eo(k)/epsovl)
         fac = 1d0/eo(k)
       do j=1,n
       do i=1,n
         oxx(i,j) = oxx(i,j) + omat(i,k)*fac*dconjg(omat(j,k)) 
         !this is eigenfunction for original problem.
       enddo
       enddo
       enddo

       h = h + delta_stabilize *oxx
#ifndef USE_TRANSPOSE
       allocate(h_t(n,n))
       do j=1,n
       do i=1,n
c          H^T 
         h_t(i,j)=h(j,i)
       enddo
       enddo
       do j=1,n
       do i=1,n
        h(i,j)= (h(i,j) + dconjg(h_t(i,j)))*0.5d0 
            ! make sure h is hermitian
       enddo
       enddo
       deallocate(h_t)
#else
       h = (h + transpose(dconjg(h)))/2d0 !make sure h is hermitian
#endif
c       s= s + oxx 
       !This procedure push up null vector(pooly-linear dependent basis) to very high eigenvalues.
       ! delta_stabilie can be somehow unclear; if you have strange eigenvalue at low energy,
       ! enlarge this.

c       omat = s !reserved
c       jobz = 'V'
c       lwork = n*n
c       allocate(wk11(lwork),rwork(max(1,3*n-2)))
c       call zheev(jobz,'U',n,omat,n,eo,wk11,lwork,rwork,ier)
c       deallocate(wk11,rwork)
c       if(ipr) then
c        write(6,*)'zhev_tk: ovlmat (after modified)='
c        do ix=0,n,5
c          if(eo(ix)>0.1d0) then
c           write(6,*) '  ... skip larger eigenvalues ...'
c           exit
c          endif
c        write(6,"(5(i5,d10.2))") (i,eo(i),i=ix+1,ix+5)
c        enddo  
c       endif
       deallocate(oxx)
C$$$      if(stabilize_zhev) then
C$$$       allocate(omatinv(n,n))!,oox(n,n))
C$$$       omatinv=0d0
C$$$       do j=1,n
C$$$       do i=1,n
C$$$       do k=1,n
C$$$         omatinv(i,j) = omatinv(i,j) + omat(i,k)*1d0/eo(k)*dconjg(omat(j,k)) !this is eigenfunction for original problem.
C$$$       enddo
C$$$       enddo
C$$$       enddo
C$$$       delta_stabilize = 1d-10
C$$$       h= h + delta_stabilize*omatinv 
C$$$       !This procedure push up null vector(pooly-linear dependent basis) to very high eigenvalues.
C$$$       ! delta_stabilie can be somehow unclear; if you have strange eigenvalue at low energy,
C$$$       ! enlarge this.
C$$$       deallocate(omatinv)
C$$$      endif
      endif




ccbug in gfortran?
cc      do j=1,n
cc      do i=1,n
cc        omatinv(i,j) = sum( omat(i,:)*1d0/eo(:)*dconjg(omat(j,:)))
cc      enddo
cc      enddo
ccccccccccc check omatinv
c      oox = matmul(s,omatinv)
c      do j=1,n
c      do i=1,n
c        epsx=oox(i,j)
c        if(i==j) write(6,*) i,j,epsx,oox(i,j)
c        if(i/=j.and.epsx>1d-5) write(6,*) i,j,epsx !,oox(i,j)
c      enddo
c      enddo
c      stop 'xxxxxxxxxxxxxxxxxxxxxxxxxx'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



ccccccccccccccccccccccccc
c     print *,' ------- omat*s*omat ------- '
c      do ix=1,n
c      do j=1,n
c        write(6,"('vvv:  ',i5,2f18.8)")j,omat(j,ix)
c      enddo
c      write(6,"(5(i5,d10.2))") ix,
c     & dreal(sum(dconjg(omat(:,ix))*matmul(s,omat(:,ix))))
c      enddo
c      stop 
ccccccccccccccccccccccc

      if(ipr) write(6,*)'eigenvalue='
! We solve  (z* h z -  epsilon z* o z) a =0
      if(ncut==0) then
c        if(epsovl_mode==1 )then
c          ni=1
c        else
          do ix= 1,n
           if(eo(ix)>epsovl) then
           ni = ix ! take i = ni...n
          exit
          endif
          enddo
c        endif
      else
        write(6,"(a,i3)")' zhev_tk: ncut=',ncut
        ni = ncut+1
      endif

C ... We solve solution in the Hilbert space spanned by zz (dimension = nm)
      nm = n-ni+1            ! this is the dimension.
      nevl=nm
      allocate(zz(n,nm))     ! zz is the projection matrix  
      do ix=ni,n
        zz(:,ix-ni+1) = omat(:,ix)/sqrt(eo(ix))
      enddo

C ... Hamiltonian  <zz|H|zz>
      allocate(hh(nm,nm),hhm(nm,n))
c This failed in ifort when hm >600 or so.-->maybe need ulimit -s unlimited.
c      hh = matmul(dconjg(transpose(zz)),matmul(h,zz))
c In anyway, blas will be better.
      call zgemm('C','N',nm,n,n,(1d0,0d0),zz,n,h,n,(0d0,0d0),hhm,nm)
      call zgemm('N','N',nm,nm,n,(1d0,0d0),hhm,nm,zz,n,(0d0,0d0),hh,nm)
c      write(6,*)'sss=',sum(abs(
c     &       matmul(dconjg(transpose(zz)),matmul(h,zz))-hh))
      deallocate(hhm)


      e=99.d0 !initialization
!! === diagonalize ===
      if (nmx .le. 0) then
        jobz = 'N'
        lwork = 4*nm
        allocate(wk11(lwork),rwork(max(1,3*nm-2)))
        call tcn('zheeveonly')
        call zheev(jobz,'U',nm,hh,nm,e,wk11,lwork,rwork,ier)
        call tcx('zheeveonly')
        deallocate(wk11, rwork)
        nev = 0
      else
        jobz = 'V'
        nev = min(nmx,nevl)
        lwork = nevl*4 ! nevl*nev ---> bug when Hydrogen case with nev=1
        allocate(wk11(lwork),rwork(max(1,3*nevl-2)))
        call tcn('zheevz')
        call zheev(jobz,'U',nm,hh,nm,e,wk11,lwork,rwork,ier)
        call tcx('zheevz')
        deallocate(wk11,rwork)
        z=1d99
        do i=1,min(nmx,nm)
          do j=1,n
            z(j,i) = scale(j) *sum(zz(j,:)*hh(:,i))  !this is eigenfunction for original problem.
          enddo
        enddo
        nev = min(nmx,nm)
      endif
      deallocate(omat,scale,zz,hh)
      call rxx(ier.ne.0,'zhev: zhegv cannot find all evals 111')
      goto 999 !exit


C$$$C --- mode 2. rescpect MTO --- a branch. not good.
C$$$! I tested it; it seems not good. ---------------------------------
C$$$ 100  continue
C$$$C ... PW part ! ... project out mto part
C$$$      allocate(omat2(nlmto,nlmto),omat2i(nlmto,nlmto))
C$$$      omat2  = s(1:nlmto,1:nlmto) !omat2=<chi_i|chi_j>
C$$$      omat2i = omat2
C$$$      call matcinv(nlmto,omat2i)  !omat2i is inverse of <chi_i|chi_j>
C$$$      ngv = n - nlmto
C$$$      deallocate(omat2)
C$$$      allocate( omat(ngv,ngv) )
C$$$c      omat=0d0
C$$$c      do i=1,ngv
C$$$c         omat(i,i) = s(nlmto+i,nlmto+i)
C$$$c        print *,'qqq',i,omat(i,i)
C$$$c      enddo                
C$$$! <PW'|PW'> matrix. PW' = |PW> - |chi_i> O^(-1)_ij <chi_j|PW>
C$$$      omat = s(nlmto+1:n,nlmto+1:n)  ! n=nlmto+ngv
C$$$     .-  matmul( dconjg(transpose( s(1:nlmto,nlmto+1:n))) ,
C$$$     .matmul( omat2i, s(1:nlmto,nlmto+1:n))   )
C$$$C ... eigenvalue of <PW'|PW'>
C$$$      jobz = 'V'
C$$$      lwork = ngv*ngv
C$$$      allocate(wk11(lwork),rwork(max(1,3*ngv-2)))
C$$$      call zheev(jobz,'U',ngv,omat,ngv,eo,wk11,lwork,rwork,ier)
C$$$      deallocate(wk11,rwork)
C$$$      if(ipr) write(6,*)'zhev_tk: ovlmat='
C$$$      if(ipr) write(6,"(5(i5,d10.2))") (i,eo(i),i=1,ngv)
C$$$! skip low  eigenvalue
C$$$      do ix= 1,ngv
C$$$        if(eo(ix)>epsovl) then
C$$$          ni = ix ! we use i = ni...n. Skip 1...ni-1
C$$$          exit
C$$$        endif
C$$$      enddo
C$$$      nm = ngv-ni+1  +nlmto      ! this is the dimension.
C$$$C ... We solve solution in the Hilbert space spanned by zz (dimension = nm)
C$$$      print *,' nlmto ngv nm=',nlmto,ngv,nm
C$$$      allocate(zz(n,nm))
C$$$      zz=0d0
C$$$      do ix=1,nlmto
C$$$        zz(ix,ix) = 1d0
C$$$      enddo
C$$$      do ix= 1, nm-nlmto
C$$$c        zz(ix+nlmto,ix+nlmto) = omat(ix,ix) /sqrt(eo(ix))
C$$$c        print *,' ppp omat',ix,zz(ix,ix)
C$$$        zz(nlmto+1:n,ix+nlmto) = omat(:,ix) /sqrt(eo(ix))
C$$$        zz(1:nlmto,  ix+nlmto) = 
C$$$     &  - matmul(omat2i, matmul(s(1:nlmto,nlmto+1:n),omat(:,ix))) 
C$$$     &  /sqrt(eo(ix))
C$$$      enddo
C$$$C ... Hamiltonian  <zz|H|zz>  <zz|Ozz>
C$$$      allocate(hh(nm,nm),oo(nm,nm))
C$$$      hh = matmul(dconjg(transpose(zz)),matmul(h,zz))
C$$$      oo = matmul(dconjg(transpose(zz)),matmul(s,zz))
C$$$      e = 99.d0 !initialization
C$$$C ... diagonalize
C$$$      if (nmx .le. 0) then
C$$$        jobz = 'N'
C$$$        lwork = 4*nm
C$$$        allocate(wk11(lwork),rwork(max(1,3*nm-2)))
C$$$        call zhegv(1,jobz,'U',nm,hh,nm,oo,nm,e,wk11,lwork,rwork,ier)
C$$$        nev = 0
C$$$      else
C$$$        jobz = 'V'
C$$$        lwork = nm*nmx
C$$$        allocate(wk11(lwork),rwork(max(1,3*nm-2)))
C$$$        call zhegv(1,jobz,'U',nm,hh,nm,oo,nm,e,wk11,lwork,rwork,ier)
C$$$        z=1d99
C$$$        do i=1,min(nmx,nm)
C$$$          do j=1,n
C$$$            z(j,i) = sum(zz(j,:)*hh(:,i))  !this is eigenfunction for original problem.
C$$$          enddo
C$$$        enddo
C$$$        nev = min(nmx,nm)
C$$$      endif
C$$$ccccccccccccccccccccccccc
C$$$c      do i=1,nm
C$$$c        write(6,"(a,i5,12d13.5)") 'zhev mode2: eee=',i,e(i)
C$$$c      enddo
C$$$c      stop ' --- vvv:zhev --- '
C$$$ccccccccccccccccccccccccc
C$$$      deallocate(wk11,rwork)
C$$$      deallocate(omat,zz,oo,omat2i)
C$$$      call rxx(ier.ne.0,'zhev: zhegv cannot find all evals 222')
C$$$c      stop 'xxxxxxxxxxxxxxxxxxx'

C$$$! exit 
 999  continue
      call tcx('zhev_tk')
      end


!!-----------------------------------
      subroutine zhev_tk2(nlmto,n,h,s,nmx,nev, e,z,ipr,  savez,getz,ifig)
!!== Eigenvalues and/or some eigenvectors of a Hermitian matrix (weighted for first nlmto basis).==
!! ----------------------------------------------------------------
!! Inputs:
!!   nlmto:dimension of MTO space of 1:nlmto i respected when diagonalization.
!!     n:    dimension of h
!!   h,n:  hermitian matrix, dimensioned h(n,n)
!!   s:    hermitian overlap matrix, 
!!   nmx:  requested number of eigenvectors to be found (and eigenvalues). If nmx>n, nmx is taken to be n.
!!         if nmx=0, nev=n (see NOTE below).
!!   ipr :print switch
!!   ifig: eigenfunctions should be saved, or get.
!!   savez,getz
!! Outputs:
!!   e:    eigenvalues
!!   nev:  number of eigenvectors (=nmx) or (=n if nmx=0)
!!   z:    eigenvectors (1..nev)  (declared as z(n,*)
!!   h and s are destroyed on exit.
!!   july2012takao
!! Saved data:
!!      When savez=T and getz=F, eigenfunctions are saved to a file with ifig
!!      When savez=F and getz=T, eigenfunctions are returned. Read from afile with ifig
!! NOTE: this can be called in the loop of ikp,isp loop. Then data are appended to a file ifig.
!!   If nmx==0, all eigenvalues are returned but without eigenfunctions.
!! -----------------------------------------------------------------------
      implicit none
      integer n,nev,nmx,ltime,ngv,ncut
      logical ipr 
      complex(8) :: h(n,n),s(n,n),z(n,*)
      complex(8),allocatable:: work(:)
      integer:: i,j,nlmto,ni,ik,ik2,k
      real(8):: epsovl,epsx,e(n),emx,eee,fac,vldummy,vudummy,abstol
      real(8),allocatable:: rwork(:) 
      integer:: ier,lwork,ix,ifi,fopna,ifig,nmx0,ifail(n),nevx
      integer,allocatable:: iwork(:)
      character*1:: jobz
      logical ::savez,getz,nexist
      integer,save:: lworksave=0
!! ---
      if(getz) then !temporary fix to acceralate
        read(ifig) nev,nmx
        read(ifig) e(1:nev)
        read(ifig) z(1:n,1:nmx)
c        if(min(nmx,n)>=nev) goto 901
        return
      endif  
 901  continue
      call tcn('zhev_tk2')
c!! --- calculation
c      allocate(scale(n))
c      scale=1d0
c      do i=1,n
c        scale(i)= 1d0 !/sqrt(s(i,i))  
c      enddo
!! ==== This scaling is to keep the Hilbert space of MTO. ====
c      scale(1:nlmto)=scale(1:nlmto)*100d0
c      do i=1,n
c        do j=1,n
c          s(i,j)= scale(i)*scale(j)*s(i,j)
c          h(i,j)= scale(i)*scale(j)*h(i,j) !this is to isolate the degenerated space.
c        enddo
c      enddo
!! === diagonalize ===
      if(nmx==0) then
         jobz='N'
         nev=n
      else
        jobz = 'V'
        nev=nmx
      endif

c      lwork = max(1,2*n-1)
c      allocate(work(lwork),rwork(max(1,3*n-2)))
c      call zhegv(1,jobz,'U',n,h,n,s,n,e,work,lwork,rwork,ier)
c      call rxx(ier.ne.0,'zhegv: zhegv cannot find all evals 111')
c      deallocate(work,rwork)
c      do i=1,nmx
c          z(:,i) = h(:,i)  !this is eigenfunction for original problem.
c       enddo
c      z(1:n,1:nmx) = h(1:n,1:nmx)
c      print *,'nev n=',nev,n

      abstol=1d-10 ! OK?
      lwork=max(1,2*n,lworksave) !OK? efficient?
      allocate(work(lwork),rwork(7*n),iwork(5*n))
      call zhegvx(1,jobz,'I','U',n,h,n,s,n,vldummy,vudummy,1,nev,abstol,nevx,e,z,n,work,lwork,rwork,iwork,ifail,ier)
      lworksave= WORK(1)  !this is optimum lwork
c      print *,'nev nevx n=',nev,nevx,n
      call rxx(nev/=nevx,'zhev_tk2: nev /=nevx something wrong. ')
      call rxx(ier.ne.0, 'zhev_tk2: zhegvx cannot find all eigen.')
      deallocate(work,iwork,rwork)


c!! ==== scale back to the eigenfuntion in original problem ====
      if(savez) then
        write(ifig) nev,nmx
        write(ifig) e(1:nev)
        write(ifig) z(1:n,1:nmx)
      endif  
      call tcx('zhev_tk2')
      end

