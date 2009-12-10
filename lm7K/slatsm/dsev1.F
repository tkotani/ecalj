      subroutine dsev1(n,h,o,wk,ipr,lx,lov,linv,nmx,emx,nev,z,e)
C- Diagonalize secular equation with overlap 
C----------------------------------------------------------------------
Ci Inputs
Ci    n:    dimension of hamiltonian
Ci    h,o:  hamiltonian, overlap matrices
Ci    wk:   work array length at least (5,11)*ndim (linv = (0,>0))
Ci    ipr:  verbosity
Ci    nmx:  maximum number of eigenvectors to be found
Ci    emx:  eigenvalue limit for eigenvectors to be found
Ci    lov:  true: overlap matrix, false, no overlap matrix
Ci   linv:  1 for inverse iteration
Ci    lx:   true, call x version for overlap handling
Co Outputs
Co    z:    eigenvectors; e, eigenvalues
Co    nev:  number of eigenvectors found
Cu Updates
Cu   07 Apr 07 Bug fix: returns true nev, not nmx in tql2 branch
Cr Remarks
Cr    h,o,z are dimensioned (n,n)
Cr    h,o are OVERWRITTEN in this routine
C----------------------------------------------------------------------
C     implicit none
C Passed parameters
      logical lov,lx
      integer n,ipr,nmx,nev,linv
      double precision h(n,n),o(n,n),z(n,n),e(n),wk(n,11),emx
C Local variables
      integer ierr,j,iprint

      call tcn('dsev1')
      nev = 0

C --- Eigenvalues of O^-1/2  H  O^-1/2 ---
      if (lov) then
        call dschd(n,n,o,wk,lx,ierr)
        call rxx(ierr.ne.0,'DSEV1: error in dschd')
        if (lx) then
          call dsredx(n,n,h,o,z)
        else
          call dsred(n,n,h,o)
        endif
      endif

      if (linv .eq. 1 .and. nmx .gt. 0) then
        call dtridx(n,n,h,wk,wk(1,4),wk(1,5),wk(1,2))
      endif

      if (nmx .le. 0) then
        call dtridx(n,n,h,e,wk,wk(1,4),wk(1,2))
        do  12  j = 1, n
   12   wk(j,1) = wk(j,1)**2
        call tqlrat(n,e,wk,ierr)
        call rxx(ierr.ne.0,'DSEV1: tqlrat cannot find all evals') 
        goto 100
      else if (linv .eq. 1) then
        call imtqlv(n,wk,wk(1,4),wk(1,5),e,wk(1,11),ierr,wk(1,6))
        call rxx(ierr.ne.0,'DSEV1: imtqlv cannot find all evals') 
C   ... Determine number of eigenvectors to be calculated
        nev = 1
        do  14  j = 2, n
   14   if (j .le. nmx .and. e(j-1) .le. emx) nev = j
        call tinvit(n,n,wk(1,1),wk(1,4),wk(1,5),nev,e,wk(1,11),z,
     .    ierr,wk(1,6),wk(1,7),wk(1,8),wk(1,9),wk(1,10))
        call rxx(ierr.ne.0,'DSEV1: tinvit cannot find all evecs') 
        call dtribx(n,n,h,wk(1,2),nev,z)
      else
        call tred2(n,n,h,e,wk(1,2),z)
        call tql2(n,n,e,wk(1,2),z,ierr)
        call rxx(ierr.ne.0,'DSEV1: tql2 cannot find all evecs') 
        nev = n
      endif

C --- Get the eigenvectors of H - E O ---
      if (.not. lov .or. nmx .le. 0) goto 100
      if (lx) then
        call dcopy(n*n,z,1,h,1)
        call dmpy(o,n,1,h,n,1,z,n,1,n,nev,n)
      else
        call dsbak(n,n,o,nev,z)
      endif

C --- Exit ---
  100 continue
      if(iprint() .ge. 60) print 600, e
  600 format(' evl='/(1x,8f10.5))
      call tcx('dsev1')
      end

