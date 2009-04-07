      subroutine diagno(n,h,o,wk,lx,lov,linv,nmx,emx,nev,z,e)
C- Diagonalize secular equation with overlap 
C----------------------------------------------------------------------
Ci Inputs
Ci    n:    dimension of hamiltonian
Ci    h,o:  hamiltonian, overlap matrices
Ci    wk:   work array length at least (3,11)*ndim (linv =0,>0))
Ci          (Special mods to htridi: need 5*ndim for linv=0)
Ci    nmx:  maximum number of eigenvectors to be found
Ci    emx:  eigenvalue limit for eigenvectors to be found
Ci    lov:  0 no overlap matrix
Ci          1 overlap matrix, return evecs of nonorthogonal H
Ci          2 overlap matrix, return evecs of orthogonal H
Ci          3 make orthogonal H, but do not diagonalize
Ci    lx:   true, call x version for overlap handling
Co Outputs
Co    z:    eigenvectors; e, eigenvalues
Co    nev:  number of eigenvectors found
Cr Remarks
Cr    h,o,z are dimensioned (n,n,2) (real followed by imag. parts)
Cr    h,o are OVERWRITTEN in this routine
C----------------------------------------------------------------------
C     implicit none
C Passed parameters
      logical lx
      integer lov,n,nmx,nev,linv
      double precision h(n,n,2),o(n,n,2),z(n,n,2),
     .                 e(n),wk(n,11),emx
C Local variables
      integer ierr,j,n2

      call tcn('diagno')
      nev = 0

C     call yprm('diagno: h',12,h,n*n,n,n,n)
C     call yprm('diagno: o',12,o,n*n,n,n,n)
C     call ztoyy(h,n,n,n,n,0,1)
C     call ztoyy(o,n,n,n,n,0,1)

C --- Eigenvalues of O^-1/2  H  O^-1/2 ---
      if (lov .gt. 0) then
        call yyhchd(n,n,o,o(1,1,2),wk,lx,.true.,ierr)
        call rxx(ierr.ne.0,'DIAGNO: error in yyhchd')
        if (lx) then
          call yyhrdx(n,n,h,h(1,1,2),o,o(1,1,2),z,z(1,1,2))
        else
          call yyhred(n,n,h,h(1,1,2),o,o(1,1,2),.true.)
        endif
      endif
      if (lov .eq. 3) goto 100

      if (linv .eq. 1 .and. nmx .gt. 0) then
#ifndef GENERIC
        call htridx(n,n,h,h(1,1,2),wk,wk(1,4),wk(1,5),wk(1,2))
#else
        call htridi(n,n,h,h(1,1,2),wk,wk(1,4),wk(1,5),wk(1,2))
#endif
      else
#ifndef GENERIC
        call htridx(n,n,h,h(1,1,2),e,wk,wk(1,4),wk(1,2))
#else
        call htridi(n,n,h,h(1,1,2),e,wk,wk,wk(1,2))
#endif
      endif

      if (nmx .le. 0) then
        do  12  j = 1, n
   12   wk(j,1) = wk(j,1)**2
        call tqlrat(n,e,wk,ierr)
        call rxx(ierr.ne.0,'DIAGNO: tqlrat cannot find all evals') 
        goto 100
      else if (linv .eq. 1) then
       call imtqlv(n,wk,wk(1,4),wk(1,5),e,wk(1,11),ierr,wk(1,6))
        call rxx(ierr.ne.0,'diagno: imtqlv cannot find all evals') 
C   --- Determine number of eigenvectors to be calculated ---
        nev = 1
        do  14  j = 2, n
          if (j .le. nmx .and. e(j-1) .le. emx) nev = j
   14   continue
        call tinvit(n,n,wk(1,1),wk(1,4),wk(1,5),nev,e,wk(1,11),z,
     .    ierr,wk(1,6),wk(1,7),wk(1,8),wk(1,9),wk(1,10))
        call rxx(ierr.ne.0,'DIAGNO: tinvit cannot find all evecs') 
      else
C --- Make unit eigenvector matrix ---
        call zinit(z,n**2)
C       call dvcpy(1d0,0,z,n+1,n)
        call dcopy(n,1d0,0,z,n+1)
        call imtql2(n,n,e,wk,z,ierr)
        call rxx(ierr.ne.0,'DIAGNO: imtql2 cannot find all evecs') 

C --- Determine number of eigenvectors to be calculated ---
        nev = 1
        do  15  j = 2, n
          if (j .le. nmx .and. e(j-1) .le. emx) nev = j
   15   continue
      endif

      if (nev .gt. 0) then
#ifndef GENERIC
        call htribx(n,n,h,h(1,1,2),wk(1,2),nev,z,z(1,1,2))
#else
        call htribk(n,n,h,h(1,1,2),wk(1,2),nev,z,z(1,1,2))
#endif

C --- Get the eigenvectors of H - E O ---
        if (lov .eq. 1) then
          if (lx) then
            n2 = n**2
            call dcopy(n2*2,z,1,h,1)
            call yympy(o,o(1,1,2),n,1,h,h(1,1,2),n,1,z,z(1,1,2),n,1,
     .        n,nev,n)
          else
            call yyhbak(n,n,o,o(1,1,2),nev,z,z(1,1,2),.true.)
          endif
        endif
      endif

  100 call tcx('diagno')
      end

