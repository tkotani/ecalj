      subroutine mksbet(ccor,ndim,xsi,xsidot,s,sdot,wk)
C- Transform structure constants to new representation
C ----------------------------------------------------------------
Ci Inputs
Ci   ccor, logical switch. If T combined correction turned on
Ci   ndim; xsi is beta - alpha; xsidot is w^-2(beta^dot - alpha^dot)
Ci   s, sdot, Bloch transformed structure constants and w^-2*S^dot
Ci   wk, work array 2*ndim**2
Co Outputs
Co   s and sdot are transformed;
Co   xsi becomes xsi^-1; xsidot becomes w^-2(xsi^-1)^dot
Cr Remarks
Cr   S and Sdot are input in representation alpha, and output
Cr   in representation beta. If ccor is F sdot, xsidot and wk
Cr   may point to any address.
Cm Memory
Cm   One double precision vector length ndim taken from the heap.
Cm   The routine needs one ndim X ndim complex work array.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      logical ccor
      integer ndim
      double precision s(ndim,ndim*2),sdot(ndim,ndim*2),
     .                 xsi(ndim),xsidot(ndim),wk(ndim,ndim*2)
C Local parameters
      integer n2,ierr
      integer LEFT,RIGHT
      parameter (LEFT=0,RIGHT=1)
C     integer fopna

C     call yprm('S',12,s,ndim*ndim,ndim,ndim,ndim)
      n2 = ndim**2

C      ierr = fopna('strx',-1,4)
C      rewind ierr
CC      print *, 'gamma-alpha, binary file strx'
CC      call ywrm(1,' ',1,ierr,' ',xsi,0,ndim,ndim,1)
C      print *, 'salp, binary file strx'
C      call ywrm(1,' ',2,ierr,' ',s,ndim*ndim,ndim,ndim,ndim)
C      call fclose(ierr)
C      pause

C --- Make xsi^-1 - S ---
      call invalp(ndim,1,xsi)
      call dscal(2*n2,-1d0,s,1)
      call daxpy(ndim,1d0,xsi,1,s,ndim+1)

C --- Invert xsi^-1 - S ---
C      call defi(owk,ndim)
C      call defdr(owk2,ndim)
C      call defdr(owk3,ndim)
C      call yyhifa(s,s(n2+1,1),ndim,ndim,w(owk),ierr)
C      if (ierr .ne. 0)
C     .  call fexit(-1,119, 'mksbet: (xsi^-1 - S) singular',0)
C      call yyhidi(s,s(n2+1,1),ndim,ndim,w(owk),0,0,w(owk2),w(owk3),1)
C      call rlse(owk)
C      do  1  j = 1, ndim
C      do  1  i = j, ndim
C        s(i,j) = s(j,i)
C        s(n2+i,j) = - s(n2+j,i)
C    1 continue
      call yqinv('h',s,n2,ndim,2,ndim,wk,ndim,ierr)
      if (ierr .ne. 0) call rx('mksbet: (xsi^-1 - S) singular')
C     call yprm('(xsi^-1-S)^-1',12,s,ndim*ndim,ndim,ndim,ndim)

      if (ccor) then
C --- define, M := (xsi^-1 - S)^-1; x := (xsi^-1)^dot; y := xsi^-1 .
C --- make x = - y * xsi^dot * y ---
       call scalpv(xsidot,ndim,1,1,xsi,LEFT)
       call scalpv(xsidot,ndim,1,1,xsi,RIGHT)
       call dscal(ndim,-1d0,xsidot,1)

C --- S^beta^dot -->  - y M (x - S^alpha^dot) M y ---
       call daxpy(ndim,-1d0,xsidot,1,sdot,ndim+1)
       call zmpy(sdot,ndim,1,n2,s,ndim,1,n2,wk,ndim,1,n2,ndim,ndim,ndim)
       call zmpy(s,ndim,1,n2,wk,ndim,1,n2,sdot,ndim,1,n2,ndim,ndim,ndim)
       call scalpv(sdot,ndim,ndim,1,xsi,LEFT)
       call scalpv(sdot(n2+1,1),ndim,ndim,1,xsi,LEFT)
       call scalpv(sdot,ndim,ndim,1,xsi,RIGHT)
       call scalpv(sdot(n2+1,1),ndim,ndim,1,xsi,RIGHT)

C --- S^beta^dot --> S^beta^dot + y M x ---
       call dcopy(2*n2,s,1,wk,1)
       call scalpv(wk,ndim,ndim,1,xsi,LEFT)
       call scalpv(wk(n2+1,1),ndim,ndim,1,xsi,LEFT)
       call scalpv(wk,ndim,ndim,1,xsidot,RIGHT)
       call scalpv(wk(n2+1,1),ndim,ndim,1,xsidot,RIGHT)
       call daxpy(2*n2,1d0,wk,1,sdot,1)

C --- S^beta^dot --> S^beta^dot + x M y ---
       call dcopy(2*n2,s,1,wk,1)
       call scalpv(wk,ndim,ndim,1,xsidot,LEFT)
       call scalpv(wk(n2+1,1),ndim,ndim,1,xsidot,LEFT)
       call scalpv(wk,ndim,ndim,1,xsi,RIGHT)
       call scalpv(wk(n2+1,1),ndim,ndim,1,xsi,RIGHT)
       call daxpy(2*n2,1d0,wk,1,sdot,1)

C --- S^beta^dot --> S^beta^dot - x ---
       call daxpy(ndim,-1d0,xsidot,1,sdot,ndim+1)
      endif

C --- Make xsi^-1 (xsi^-1 - S)^-1 xsi^-1 - xsi^-1 ---
      call scalpv(s,ndim,ndim,1,xsi,LEFT)
      call scalpv(s(n2+1,1),ndim,ndim,1,xsi,LEFT)
      call scalpv(s,ndim,ndim,1,xsi,RIGHT)
      call scalpv(s(n2+1,1),ndim,ndim,1,xsi,RIGHT)
      call daxpy(ndim,-1d0,xsi,1,s,ndim+1)

      end

