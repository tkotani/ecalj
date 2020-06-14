      module NaNum
      use,intrinsic :: iso_fortran_env
      use,intrinsic :: ieee_arithmetic
      real(8),parameter :: NaN = transfer(-1_int64,0.0_real64)
      end module
      
c      program test1
c      Use NaNum,only: NaN
c      real(8):: a(8)=NaN,c=NaN
c      write(6,*) a
c      write(6,*) 'ccc=',c
c      end
