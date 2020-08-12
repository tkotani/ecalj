      subroutine prtev(t,n,evl,nmx,emx,nev)
C- Printout the eigenvalues
C ----------------------------------------------------------------------
Ci Inputs
Ci   evl   :eigenvalues
Ci   t     :eigenvectors
Ci   n     :dimension of hamiltonian,overlap
Ci   nmx:  : number of eigenvectors
c??? Ci   emx:  :eigenvalue limit for eigenvectors to be found
Ci   nev   : number of eigenvalues
Cu Updates
Cu   28 Aug 04 prints out |evec| at high verbosities
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer n,nmx,nev
      double precision emx,evl(nev)
      double complex t(n,nmx)
C ... Local parameters
      integer i,ipr,j,stdo,lgunit
      stdo = lgunit(1)
      call getpr(ipr)
C     ipr = 100
      if (ipr .ge. 30 .and. ipr .lt. 100) then
        j = min(9,nmx)
        if (ipr .ge. 35) j = nev
        write(stdo,103) (evl(i), i=1,nev)
  103   format(9f8.4)
        write(6,"(' prtev: nev nmx ndim evl(nev)=',3i5,f12.5)") nev,nmx,n,evl(nev)
      endif
      if (ipr .ge. 100) then
        do i = 1, nmx
          write(stdo,863) i,evl(i),(cdabs(t(j,i)),j=1,n)
        enddo
  863   format(' i=',i5,'   evl=',f12.6,'  abs(z)='/(8f9.4))
      endif
      end

