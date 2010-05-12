      subroutine prtev(t,n,evl,nmx,emx,nev)
C- Printout the eigenvalues
C ----------------------------------------------------------------------
Ci Inputs
Ci   evl   :eigenvalues
Ci   t     :eigenvectors
Ci   n     :dimension of hamiltonian,overlap
Ci   nmx:  :maximum number of eigenvectors to be found
Ci   emx:  :eigenvalue limit for eigenvectors to be found
Ci   nev   :actual number of eigenvectors generated
Cu Updates
Cu   28 Aug 04 prints out |evec| at high verbosities
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer n,nmx,nev
      double precision emx,evl(n)
      double complex t(n,n)
C ... Local parameters
      integer i,ipr,j,stdo,lgunit

      stdo = lgunit(1)
      call getpr(ipr)
C     ipr = 100

      if (ipr .ge. 30 .and. ipr .lt. 100) then
        j = min(9,n)
        if (ipr .ge. 35) j = n
        write(stdo,103) (evl(i), i=1,j)
  103   format(9f8.4)
#ifndef MPI
c #if .not. (MPI|MPIK)
        if (ipr .ge. 36 .and. nev .gt. 0) call awrit5(
     .    ' nev, nmx, dim=  %i  %i  %i  ev(nev) = %1;5d  emx '//
     .    '= %1;5g',' ',80,stdo,nev,nmx,n,evl(nev),emx)
        call ftflsh(stdo)
#endif
      endif

      if (ipr .ge. 100) then
        do  30  i = 1, min(n,nmx)
          write(stdo,863) i,evl(i),(cdabs(t(j,i)),j=1,n)
   30   continue
  863   format(' i=',i5,'   evl=',f12.6,'  abs(z)='/(8f9.4))
        call ftflsh(stdo)
      endif

      end

