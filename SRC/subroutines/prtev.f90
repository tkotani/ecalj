subroutine prtev(t,n,evl,nmx,nev)
  use m_lgunit,only:stdo
  !- Printout the eigenvalues
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   evl   :eigenvalues
  !i   t     :eigenvectors
  !i   n     :dimension of hamiltonian,overlap
  !i   nmx:  : number of eigenvectors
  !i   nev   : number of eigenvalues
  !u Updates
  !u   28 Aug 04 prints out |evec| at high verbosities
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: n,nmx,nev
  double precision :: evl(nev)
  double complex t(n,nmx)
  ! ... Local parameters
  integer :: i,ipr,j
  !      stdo = lgunit(1)
  call getpr(ipr)
  ipr = 30 !100
  if (ipr >= 30 .AND. ipr < 100) then
     j = min(9,nmx)
     if (ipr >= 35) j = nev
     write(stdo,103) (evl(i), i=1,nev)
103  format(9f8.4)
     write(6,"(' prtev: nev nmx ndim evl(nev)=',3i5,f12.5)") nev,nmx,n,evl(nev)
  endif
  if (ipr >= 100) then
     do i = 1, nmx
        write(stdo,863) i,evl(i),(cdabs(t(j,i)),j=1,n)
     enddo
863  format(' i=',i5,'   evl=',f12.6,'  abs(z)='/(8f9.4))
  endif
end subroutine prtev

