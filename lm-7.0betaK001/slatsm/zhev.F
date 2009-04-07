#define F90 1
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
      double precision h(1),s(1),e(n),wk(1),z(2),emx

#if LAPACK
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
#else
C Local parameters
      integer i,j,k,m,mdim,k1,k2,k3,k4,k5,k6,k7,k8,k9,n2,ier

      call tcn('zhev')
      n2 = n**2

C --- Take care of n=1 case ---
      if (n .eq. 1) then
         e(1) = h(1)
         if (lov .gt. 0) e(1) = h(1)/s(1)
         if (1 .gt. nev .or. e(1) .gt. emx) return
         z(1) = 1
         z(2) = 0
         return
      endif

C --- Separate real and imaginary parts ---
      mdim = 2*lh
      call ztoy(h,lh,n,n,0)
      if (lov .gt. 0) call ztoy(s,lh,n,n,0)

C --- Debugging: eigenvalues of overlap ---
C      call htridx(mdim,n,s(1),s(lh+1),e,wk,wk(3*n+1),wk(n+1))
C      do  11  j = 1, n
C   11 wk(j) = wk(j)**2
C      call tqlrat(n,e,wk,ier)
C      write(6,600) e
C  600 format(' evl='/(1x,1p,5e14.6))
C      call rx('eigenvalues of overlap')

C --- H <- S^-1/2  H  S^-1/2 ---
      if (lov .gt. 0) then
        call yyhchd(mdim,n,s,s(lh+1),wk,lx,.true.,ier)
        call rxx(ier.ne.0,'ZHEV: error in yyhchd')
        if (lx) then
          if (lz .ne. lh) call rx('zhev not ready for lz ne lh')
          call yyhrdx(mdim,n,h,h(lh+1),s,s(lh+1),z,z(lz+1))
        else
          call yyhred(mdim,n,h,h(lh+1),s,s(lh+1),.true.)
        endif
      endif

C --- Transform to tridiagonal matrix ---
      if (linv) then
        k1 = 1
        k2 = k1 + 3*n
        k3 = k2 + n
        k4 = k3 + n
        k5 = k4 + n
        k6 = k5 + n
        k7 = k6 + n
        k8 = k7 + n
        k9 = k8 + n
#ifndef GENERIC
        call htridx(mdim,n,h(1),h(lh+1),wk(k1),wk(k2),wk(k3),wk(n+1))
#else
        call htridi(mdim,n,h(1),h(lh+1),wk(k1),wk(k2),wk(k3),wk(n+1))
#endif
      else
#ifndef GENERIC
        call htridx(mdim,n,h(1),h(lh+1),e,wk,wk(3*n+1),wk(n+1))
#else
        call htridi(mdim,n,h(1),h(lh+1),e,wk,wk(3*n+1),wk(n+1))
#endif
      endif

C --- Eigenvalues only ---
      if (nmx .le. 0) then
        do  12  j = 1, n
   12   wk(j) = wk(j)**2
        call tqlrat(n,e,wk,ier)
        call rxx(ier.ne.0,'ZHEV: tqlrat cannot find all evals')
        nev = 0
        goto 100

C --- Eigenvalues and eigenvectors ---
      else if (linv) then
        call imtqlv(n,wk(k1),wk(k2),wk(k3),e,wk(k9),ier,wk(k4))
        call rxx(ier.ne.0,'ZHEV: imtqlv cannot find all evals')
C   --- Determine number of eigenvectors to be calculated ---
        nev = 1
        do  14  j = 2, n
          if (j .le. nmx .and. e(j-1) .le. emx) nev = j
   14   continue
        call tinvit(2*lz,n,wk(k1),wk(k2),wk(k3),nev,e,wk(k9),z,ier,
     .    wk(k4),wk(k5),wk(k6),wk(k7),wk(k8))
        call rxx(ier.ne.0,'ZHEV: tinvit cannot find all evecs')
      else
        do  17  j = 1, n
          k = (j-1)*2*lz
          m = k+n
          do  16  i = k+1, m
   16     z(i) = 0d0
          z(k+j) = 1d0
   17   continue
        call imtql2(2*lz,n,e,wk,z,ier)
C       call prmx('eval',e,n,n,1)
        call rxx(ier.ne.0,'ZHEV: imtql2 cannot find all evecs')

C   --- Determine number of eigenvectors to be calculated ---
        nev = 1
        do  15  j = 2, n
          if (j .le. nmx .and. e(j-1) .le. emx) nev = j
   15   continue
      endif
C     call prmx('eval',e,n,n,1)

      if (nev .gt. 0) then
#ifndef GENERIC
        if (lz .ne. lh) call rx('zhev not ready for lz ne lh')
        call htribx(mdim,n,h(1),h(lh+1),wk(n+1),nev,z(1),z(lz+1))
#else
        call htribk(mdim,n,h(1),h(lh+1),wk(n+1),nev,z(1),z(lz+1))
#endif

C --- Get the eigenvectors of H - E O ---
        if (lov .gt. 0) then
          if (lx) then
            call ymcpy(z,2*lz,1,lz,h,mdim,1,lh,n,n)
            call yympy(s,s(lh+1),mdim,1,h,h(lh+1),mdim,1,z,z(lz+1),2*lz,
     .        1,n,nev,n)
          else
            call yyhbak(mdim,n,s,s(lh+1),nev,z,z(lz+1),.true.)
          endif
        endif

C   --- Convert eigenvectors to double complex storage ---
        call ztoy(z,lz,n,nev,1)
C       call zprm('evecs',2,z,lz,n,nmx)
      endif

#endif

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
      double precision h(1),s(1),e(n),wk(1),z(2),emx
      logical lov

#if LAPACK
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
#else
C Local parameters
      integer i,j,k,m,mdim,k1,k2,k3,k4,k5,k6,k7,k8,k9,n2,ier

C     call prm('(5f20.15)',s,n,n)
      call tcn('zhev')
      n2 = n**2

C --- Take care of n=1 case ---
      if (n .eq. 1) then
         e(1) = h(1)
         if (lov) e(1) = h(1)/s(1)
         if (1 .gt. nev .or. e(1) .gt. emx) return
         z(1) = 1
         z(2) = 0
         return
      endif

CC --- Use diagno if lx ---
C      if (lx) then
C        call cplx2r(n2,0,h,z)
C        if (lov) call cplx2r(n2,0,s,z)
C        i = 0
C        if (lov) i = 1
C        call diagno(n,h,s,wk,lx,i,linv,nmx,emx,nev,z,e)
C        call cplx2r(n2,1,z,h)
C        goto 40
C      endif

C --- Separate real and imaginary parts ---
      mdim = 2*n
      call ztoy(h,n,n,n,0)
      if (lov) call ztoy(s,n,n,n,0)
C      mdim = 2*n
C      do  10  j = 1, n
C        k = (j-1)*mdim + 1
C        l = k + n
C        call dcopy(n,h(k+1),2,wk,1)
C        call dcopy(n,h(k),2,h(k),1)
C        call dcopy(n,wk,1,h(l),1)
C        if (lov) then
C          call dcopy(n,s(k+1),2,wk,1)
C          call dcopy(n,s(k),2,s(k),1)
C          call dcopy(n,wk,1,s(l),1)
C        endif
C   10 continue


C --- Debugging: eigenvalues of overlap ---
C      call htridx(mdim,n,s(1),s(n+1),e,wk,wk(3*n+1),wk(n+1))
C      do  11  j = 1, n
C   11 wk(j) = wk(j)**2
C      call tqlrat(n,e,wk,ier)
C      write(6,600) e
C  600 format(' evl='/(1x,1p,5e14.6))
C      call rx('eigenvalues of overlap')

C --- H <- S^-1/2  H  S^-1/2 ---
      if (lov) then
        call yyhchd(mdim,n,s,s(n+1),wk,lx,.true.,ier)
        call rxx(ier.ne.0,'ZHEV: error in yyhchd')
        if (lx) then
          call yyhrdx(mdim,n,h,h(n+1),s,s(n+1),z,z(n+1))
        else
          call yyhred(mdim,n,h,h(n+1),s,s(n+1),.true.)
        endif
      endif

C --- Transform to tridiagonal matrix ---
      if (linv) then
        k1 = 1
        k2 = k1 + 3*n
        k3 = k2 + n
        k4 = k3 + n
        k5 = k4 + n
        k6 = k5 + n
        k7 = k6 + n
        k8 = k7 + n
        k9 = k8 + n
#ifndef GENERIC
        call htridx(mdim,n,h(1),h(n+1),wk(k1),wk(k2),wk(k3),wk(n+1))
#else
        call htridi(mdim,n,h(1),h(n+1),wk(k1),wk(k2),wk(k3),wk(n+1))
#endif
      else
#ifndef GENERIC
        call htridx(mdim,n,h(1),h(n+1),e,wk,wk(3*n+1),wk(n+1))
#else
        call htridi(mdim,n,h(1),h(n+1),e,wk,wk(3*n+1),wk(n+1))
#endif
      endif

C --- Eigenvalues only ---
      if (nmx .le. 0) then
        do  12  j = 1, n
   12   wk(j) = wk(j)**2
        call tqlrat(n,e,wk,ier)
        call rxx(ier.ne.0,'ZHEV: tqlrat cannot find all evals')
        nev = 0
        goto 100

C --- Eigenvalues and eigenvectors ---
      else if (linv) then
        call imtqlv(n,wk(k1),wk(k2),wk(k3),e,wk(k9),ier,wk(k4))
        call rxx(ier.ne.0,'ZHEV: imtqlv cannot find all evals')
C   --- Determine number of eigenvectors to be calculated ---
        nev = 1
        do  14  j = 2, n
          if (j .le. nmx .and. e(j-1) .le. emx) nev = j
   14   continue
        call tinvit(mdim,n,wk(k1),wk(k2),wk(k3),nev,e,wk(k9),z,ier,
     .    wk(k4),wk(k5),wk(k6),wk(k7),wk(k8))
        call rxx(ier.ne.0,'ZHEV: tinvit cannot find all evecs')
      else
        do  17  j = 1, n
          k = (j-1)*mdim
          m = k+n
          do  16  i = k+1, m
   16     z(i) = 0d0
          z(k+j) = 1d0
   17   continue
        call imtql2(mdim,n,e,wk,z,ier)
        call rxx(ier.ne.0,'ZHEV: imtql2 cannot find all evecs')

C   --- Determine number of eigenvectors to be calculated ---
        nev = 1
        do  15  j = 2, n
          if (j .le. nmx .and. e(j-1) .le. emx) nev = j
   15   continue
      endif

      if (nev .gt. 0) then
#ifndef GENERIC
        call htribx(mdim,n,h(1),h(n+1),wk(n+1),nev,z(1),z(n+1))
#else
        call htribk(mdim,n,h(1),h(n+1),wk(n+1),nev,z(1),z(n+1))
#endif

C --- Get the eigenvectors of H - E O ---
        if (lov) then
          if (lx) then
            call dcopy(n2*2,z,1,h,1)
            call yympy(s,s(n+1),mdim,1,h,h(n+1),mdim,1,z,z(n+1),mdim,1,
     .        n,nev,n)
          else
            call yyhbak(mdim,n,s,s(n+1),nev,z,z(n+1),.true.)
          endif
        endif

C   --- Convert eigenvectors to double complex storage ---
        call ztoy(z,n,n,nev,1)
C       call zprm('evecs',2,z,lz,n,min(nev,nmx))
      endif

C  40 if (nev .gt. 0 .and. ltime .ge. 0) print 337, nev,nmx,n,emx
C  337 format(' nev, nevmx, ndim=',3i4,' emx=',f10.5)
#endif

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
      double precision h(1),s(1)
      logical lx
C Local parameters
      integer ier,itype,lmode
      double complex one
      parameter (one=(1.0d+0,0.0d+0))

      call tcn('zhev2')

      lmode = mode

#ifndef LAPACK
      call rx('zhev2:  only implemented in LAPACK mode, sorry')
#endif

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
#if F90
       real(8),allocatable:: e(:),wk(:)
#else
      double precision e(n),wk(11*n)
#endif
      logical lx
      integer i,lmode
      double complex zer,one
      parameter (zer=(0d0,0d0),one=(1d0,0d0))
      parameter(lx=.true.)

      call tcn('zhev3')

#ifndef LAPACK
      call rx('zhev3:  only implemented in LAPACK mode, sorry')
#endif

      lmode = mode

C --- Eigenvectors and eigenvalues of s ---
      if (mod(mode,2) .ne. 0) then
#if F90
      allocate(e(n),wk(11*n))
#endif
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
#if F90
      deallocate(e,wk)
#endif
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
#if F90
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

#else
      call rx('zhevo only implemented with F90')
#endif

      end





ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C we assume LAPACK and F90
      subroutine zhev_tk(nlmto,n,h,s,nmx, epsovl,ncut,
     &   nevl,nev,e,eo,lz,z,ipr)
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
      double precision e(n),wk(1),emx
      complex(8):: h(n,n),s(n,n),z(lz,n)
      complex(8),allocatable:: omat(:,:),imat(:,:),wk11(:),
     & zz(:,:),hh(:,:),omat2(:,:),pmat(:,:),omat2i(:,:),oo(:,:),
     & hhm(:,:)
      integer:: i,j,nlmto,ni,nm,ik,ik2
      real(8):: omax_mto, omax_g,ex ,  ddd2,ddd,epsovl
      real(8),allocatable:: scale(:),rwork(:)
      real(8):: eo(1)
C Local parameters
      integer ier,lwork,ix
      character*1 jobz
c
      call tcn('zhev_tk')
ccccccccccccccccccccccccccccc
c      do i=1,n
c      do j=1,n
c        if(i<5) print *,' 444 sss=',j,i,abs(s(j,i))
c      enddo
c      enddo
ccccccccccccccccccccccccccc


! if mode 2. For mode 1, comment 'goto 100'
c      goto 100


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
      allocate(omat(n,n))
      omat = s !reserved

c$$$cccccccccccccccccccccccccccc
c$$$      do i=1,n
c$$$         ddd=0d0
c$$$         do j=1,n
c$$$          if(i<=5.and.j<=5) 
c$$$     &    write(6,"('mmm: sss=',2i5,f18.8)") j,i,abs(s(j,i))
c$$$          if(j==i) then
c$$$c            write(6,"(a,2i4,f9.3,i4,f9.3)")'mmm: diag ',i,i,abs(s(j,i))
c$$$            cycle
c$$$          endif
c$$$          if(abs(omat(j,i) )>ddd) then
c$$$              ddd=abs(omat(j,i))
c$$$              ik=j
c$$$          endif
c$$$         enddo
c$$$         ddd2=0d0
c$$$         do j=1,n
c$$$          if(j==i) cycle
c$$$          if(ddd>abs(omat(j,i)).and.abs(omat(j,i))>ddd2) then
c$$$              ddd2=abs(omat(j,i))
c$$$              ik2=j
c$$$          endif
c$$$         enddo
c$$$         write(6,"(a,2i4,f19.13,i4,f19.13)")'mmm: max '
c$$$     &   , i, ik,ddd,ik2,ddd2
c$$$      enddo

cccccccccccccccccccccccccccc
c      do j=1,n
c      print *,'qqq =',j,omat(1,j)
c      enddo
c      do i=2,n
c      do j=2,n
c         if (i/=j) omat(i,j)=0d0
c         if (i==j) omat(i,j)=0.00d0
c      enddo
c      enddo
c      do j=150,n
c         omat(1,j)=0d0
c         omat(j,1)=0d0
c      enddo
c      print *,'xxx sum abs=',sum(abs(omat(2:n,2:n)))
c      print *,'xxx sum abs=',sum(abs(omat(2:n,1)))
c      print *,'xxx sum abs=',sum(abs(omat(1,2:n)))
c      print *,'xxx sum abs=',abs(omat(1,1))
cc      s=omat
c      omat(2:nlmto,nlmto+1:n)=0d0
c      omat(nlmto+1:n,1)=0d0
c      stop 
ccccccccccccccccccccccccccc



C ... eigenvalue of s (omat)
      jobz = 'V'
      lwork = n*n
      allocate(wk11(lwork),rwork(max(1,3*n-2)))
      call zheev('V','U',n,omat,n,eo,wk11,lwork,rwork,ier)
      deallocate(wk11,rwork)
      if(ipr) write(6,*)'zhev_tk: ovlmat='
      if(ipr) write(6,"(5(i5,d10.2))") (i,eo(i),i=1,n)
c      write(6,*)' zheev: diagonalization of overlap matrix info=',ier

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
        do ix= 1,n
        if(eo(ix)>epsovl) then
           ni = ix ! take i = ni...n
           exit
        endif
        enddo
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

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      do ix=1,nm
c        write(6,"('hh dia',i5,2f13.4)") ix,hh(ix,ix)
c      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

C ... diagonalize
      if (nmx .le. 0) then
        jobz = 'N'
        lwork = 4*nm
        allocate(wk11(lwork),rwork(max(1,3*nm-2)))
        call zheev(jobz,'U',nm,hh,nm,e,wk11,lwork,rwork,ier)
        deallocate(wk11, rwork)
        nev = 0
      else
        jobz = 'V'
        nev = min(nmx,nevl)
        lwork = nevl*4 ! nevl*nev ---> bug when Hydrogen case with nev=1
        allocate(wk11(lwork),rwork(max(1,3*nevl-2)))
        call zheev(jobz,'U',nm,hh,nm,e,wk11,lwork,rwork,ier)
        deallocate(wk11,rwork)
        z=1d99
        do i=1,min(nmx,nm)
        do j=1,n
           z(j,i) = scale(j) *sum(zz(j,:)*hh(:,i))  !this is eigenfunction for original problem.
        enddo
        enddo
        nev = min(nmx,nm)
      endif
c      nev=nm
c      deallocate(wk11,rwork) ! move to the line after calling zheev 
      deallocate(omat,scale,zz,hh)

c      do i=1,nev
c      do j=1,nev
c        if(abs(sum( dconjg(z(:,i))*matmul(omat2,z(:,j)) ))>1d-10) then
c        write(6,"(' matooo ',2i4,2d13.6)")  i,j,
c     &   sum( dconjg(z(:,i))*matmul(omat2,z(:,j)) )
c        endif
c      enddo
c      enddo

cccccccccccccccccccccccc
c      do i=1,nm
c        write(6,"(a,i5,12d13.5)") 'eee=',i,e(i)
c      enddo
c      stop ' --- vvv:zhev --- '
ccccccccccccccccccccccccc
      call rxx(ier.ne.0,'zhev: zhegv cannot find all evals 111')
      goto 999 !exit




C --- mode 2. rescpect MTO --- a branch. not good.
! I tested it; it seems not good. ---------------------------------
 100  continue
C ... PW part ! ... project out mto part
      allocate(omat2(nlmto,nlmto),omat2i(nlmto,nlmto))
      omat2  = s(1:nlmto,1:nlmto) !omat2=<chi_i|chi_j>
      omat2i = omat2
      call matcinv(nlmto,omat2i)  !omat2i is inverse of <chi_i|chi_j>
      ngv = n - nlmto
      deallocate(omat2)
      allocate( omat(ngv,ngv) )
c      omat=0d0
c      do i=1,ngv
c         omat(i,i) = s(nlmto+i,nlmto+i)
c        print *,'qqq',i,omat(i,i)
c      enddo                
! <PW'|PW'> matrix. PW' = |PW> - |chi_i> O^(-1)_ij <chi_j|PW>
       omat = s(nlmto+1:n,nlmto+1:n)  ! n=nlmto+ngv
     .      -  matmul( dconjg(transpose( s(1:nlmto,nlmto+1:n))) ,
     .                   matmul( omat2i, s(1:nlmto,nlmto+1:n))   )
C ... eigenvalue of <PW'|PW'>
      jobz = 'V'
      lwork = ngv*ngv
      allocate(wk11(lwork),rwork(max(1,3*ngv-2)))
      call zheev(jobz,'U',ngv,omat,ngv,eo,wk11,lwork,rwork,ier)
      deallocate(wk11,rwork)
      if(ipr) write(6,*)'zhev_tk: ovlmat='
      if(ipr) write(6,"(5(i5,d10.2))") (i,eo(i),i=1,ngv)
! skip low  eigenvalue
      do ix= 1,ngv
        if(eo(ix)>epsovl) then
           ni = ix ! we use i = ni...n. Skip 1...ni-1
           exit
        endif
      enddo
      nm = ngv-ni+1  +nlmto      ! this is the dimension.
C ... We solve solution in the Hilbert space spanned by zz (dimension = nm)
      print *,' nlmto ngv nm=',nlmto,ngv,nm
      allocate(zz(n,nm))
      zz=0d0
      do ix=1,nlmto
         zz(ix,ix) = 1d0
      enddo
      do ix= 1, nm-nlmto
c        zz(ix+nlmto,ix+nlmto) = omat(ix,ix) /sqrt(eo(ix))
c        print *,' ppp omat',ix,zz(ix,ix)
        zz(nlmto+1:n,ix+nlmto) = omat(:,ix) /sqrt(eo(ix))
        zz(1:nlmto,  ix+nlmto) = 
     &    - matmul(omat2i, matmul(s(1:nlmto,nlmto+1:n),omat(:,ix))) 
     &    /sqrt(eo(ix))
      enddo
C ... Hamiltonian  <zz|H|zz>  <zz|Ozz>
      allocate(hh(nm,nm),oo(nm,nm))
      hh = matmul(dconjg(transpose(zz)),matmul(h,zz))
      oo = matmul(dconjg(transpose(zz)),matmul(s,zz))
      e = 99.d0 !initialization
C ... diagonalize
      if (nmx .le. 0) then
        jobz = 'N'
        lwork = 4*nm
        allocate(wk11(lwork),rwork(max(1,3*nm-2)))
        call zhegv(1,jobz,'U',nm,hh,nm,oo,nm,e,wk11,lwork,rwork,ier)
        nev = 0
      else
        jobz = 'V'
        lwork = nm*nmx
        allocate(wk11(lwork),rwork(max(1,3*nm-2)))
        call zhegv(1,jobz,'U',nm,hh,nm,oo,nm,e,wk11,lwork,rwork,ier)
        z=1d99
        do i=1,min(nmx,nm)
        do j=1,n
           z(j,i) = sum(zz(j,:)*hh(:,i))  !this is eigenfunction for original problem.
        enddo
        enddo
        nev = min(nmx,nm)
      endif
ccccccccccccccccccccccccc
c      do i=1,nm
c        write(6,"(a,i5,12d13.5)") 'zhev mode2: eee=',i,e(i)
c      enddo
c      stop ' --- vvv:zhev --- '
ccccccccccccccccccccccccc
      deallocate(wk11,rwork)
      deallocate(omat,zz,oo,omat2i)
      call rxx(ier.ne.0,'zhev: zhegv cannot find all evals 222')
c      stop 'xxxxxxxxxxxxxxxxxxx'

! exit 
 999  continue
      call tcx('zhev_tk')
      end


      subroutine matcinv(n,a)
C --- a inverse is returned.
      implicit none
      integer :: n, info, ipiv(n)
      complex(8):: a(n,n)
      complex(8),allocatable:: work(:)
      call zgetrf(n,n,a,n,ipiv,info)
      if(info/=0) then
        print *,' matcinv: zegtrf info=',info
        stop    ' matcinv: zegtrf '
      endif
      allocate(work(n*n))
      call zgetri(n,a,n,ipiv,work,n*n,info)
      deallocate(work)
      if(info/=0) then
        print *,'matcinv: zegtri info=',info
        stop    'matcinv: zegtri '
      endif
      end


