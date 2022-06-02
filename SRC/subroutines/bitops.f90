integer function getdig(n,i,base)
  !- Extracts one digit from an integer
  ! ----------------------------------------------------------------
  !i Inputs
  !i   n,i,base
  !o Outputs
  !o   getdig = ith digit from n, base "base"; eg 4=getdig(12345,1,10)
  ! ----------------------------------------------------------------
  implicit none
  integer :: n,i,base
  getdig = mod(n/base**i,base)
end function getdig
