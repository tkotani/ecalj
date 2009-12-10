      subroutine diagcv(s,h,t,n,evl,nmx,emx,nev)
C- Diagonalizes and returns the nev lowest eigenstates
C ----------------------------------------------------------------------
Ci Inputs
Ci   s     :overlap
Ci   h     :hamiltonian
Ci   n     :dimension of hamiltonian,overlap
Ci   nmx:  :maximum number of eigenvectors to be found
Ci   emx:  :eigenvalue limit for eigenvectors to be found
Ci   nev   :actual number of eigenvectors generated
Co Outputs
Co   evl   :eigenvalues
Co   t     :eigenvectors
Cr Remarks
Cr   z must be at least of dimension z(n,n), even though nev<n.
Cr   h and s are destroyed on exit.
Cu Updates
Cu   30 May 00 Adapted from nfp diagcv.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer n,nmx,nev
      double precision emx,evl(n)
      double complex s(n,n),h(n,n),t(n,n)
      logical lx,lnv
C ... Local parameters
      integer i,ipr,j,oww,stdo,lgunit
      integer w(1)
      common /w/ w

      stdo = lgunit(1)
      call getpr(ipr)
      lx = .true.
      lnv = .false.
      call defrr(oww,  11*n)
      call zhev(n,h,s,.true.,lx,nmx,emx,nev,w(oww),lnv,-1,evl,t)
      call rlse(oww)

      if (ipr .ge. 30 .and. ipr .lt. 100) then
        j = min(9,n)
        if (ipr .ge. 35) j = n
        write(stdo,103) (evl(i), i=1,j)
  103   format(9f8.4)
        if (ipr .ge. 36 .and. nev .gt. 0) call awrit5(
     .    ' nev, nmx, dim=  %i  %i  %i  ev(nev) = %1;5d  emx '//
     .    '= %1;5g',' ',80,stdo,nev,nmx,n,evl(nev),emx)
      endif

      if (ipr .ge. 100) then
        do  30  i = 1, n
   30   write(stdo,863) i,evl(i),(t(j,i),j=1,n)
  863   format(' i=',i5,'   evl=',f12.6/(8f9.4))
      endif
      call ftflsh(stdo)

      end

