module NaNum
  !      use,intrinsic :: iso_fortran_env
  !c      use,intrinsic :: ieee_arithmetic
  !c      real(8),parameter :: NaN = transfer(-1_int64,0.0_real64) !not good for module in ifort
  real(8),parameter :: NaN = -9999999 !transfer(-1_int64,0.0_real64)
end module NaNum

!      program test1
!      Use NaNum,only: NaN
!      real(8):: a(8)=NaN,c=NaN
!      write(6,*) a
!      write(6,*) 'ccc=',c
!      end
