      module NaNum
      use,intrinsic :: iso_fortran_env
      use,intrinsic :: ieee_arithmetic
      real(8),parameter :: NaN = transfer(-1_int64,0.0_real64)
      end module
      
      program test1
      Use NaNum,only: NaN
      real(8):: a(8)=NaN,c=NaN
      write(6,*) a
      write(6,*) 'ccc=',c
      end
