real(8) function dsum(n,dx,incx)

  !     takes the sum of the values.  Adapted from:
  !     jack dongarra, linpack, 3/11/78.

  double precision :: dx(1)
  integer :: i,incx,n,nincx

  dsum = 0d0
  if (n <= 0) return

  nincx = n*incx
  do  10  i = 1, nincx,incx
     dsum = dsum + dx(i)
10 enddo
END function dsum


