integer function mxint(n,ivec)
  !- Return the largest integer in ivec
  !     implicit none
  integer :: n,ivec(n)
  integer :: imax,i

  mxint = 0
  if (n <= 0) return
  imax = ivec(1)
  do   i = 1, n
     imax = max(imax,ivec(i))
  enddo
  mxint = imax
end function mxint

