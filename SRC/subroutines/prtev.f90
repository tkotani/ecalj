subroutine prtev(t,n,evl,nmx,nev)!- Printout the eigenvalues
  use m_lgunit,only:stdo
  implicit none
  integer :: n,nmx,nev
  double precision :: evl(nev)
  double complex t(n,nmx)
  integer :: i,ipr,j
  write(stdo,"(9f8.4)") (evl(i), i=1,nev)
  write(stdo,"(' prtev: nev nmx ndim evl(nev)=',3i5,f12.5)") nev,nmx,n,evl(nev)
end subroutine prtev

