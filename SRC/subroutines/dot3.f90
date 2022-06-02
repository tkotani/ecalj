      double precision function dot3(n,a,b,c)
C- Inner product of three functions
C     implicit none
      integer n,i
      double precision a(n),b(n),c(n),xx

      xx = 0d0
      do    i = 1, n
       xx = xx + a(i)*b(i)*c(i)
      enddo
      dot3 = xx
      end

