      subroutine dpdot(a,b,n,sum)
C     implicit none

      integer n,i
      double precision a(1),b(1),sum
#if APOLLO | HP
      double precision vec_$ddot
      sum = vec_$ddot(a,b,n)
#elif CRAY
      sum = sdot(n,a,1,b,1)
#elif BLAS
      sum = ddot(n,a,1,b,1)
#else
      sum = 0d0
      do  10  i = 1, n
   10 sum = sum + a(i)*b(i)
#endif
      end

