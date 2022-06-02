real(8) function dot3(n,a,b,c)
  !- Inner product of three functions
  !     implicit none
  integer :: n,i
  double precision :: a(n),b(n),c(n),xx
  xx = 0d0
  do    i = 1, n
     xx = xx + a(i)*b(i)*c(i)
  enddo
  dot3 = xx
END function dot3

