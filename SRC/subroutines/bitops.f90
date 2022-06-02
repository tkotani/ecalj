      integer function getdig(n,i,base)
C- Extracts one digit from an integer
C ----------------------------------------------------------------
Ci Inputs
Ci   n,i,base
Co Outputs
Co   getdig = ith digit from n, base "base"; eg 4=getdig(12345,1,10)
C ----------------------------------------------------------------
      implicit none
      integer n,i,base
      getdig = mod(n/base**i,base)
      end
