      subroutine setfac(n,fac)
C- set up array of factorials.
C     implicit none
      integer n,i
      double precision fac(0:n)

      fac(0) = 1d0
      do    i = 1, n
         fac(i) = i*fac(i-1)
      enddo   
      end
      subroutine stdfac(n,df)
C- Set up array of double factorials.
c  for odd numbers,  makes 1*3*5*..*n
c  for even numbers, makes 2*4*6*..*n
C     implicit none
      integer n,i
      double precision df(0:n)

      df(0) = 1d0
      df(1) = 1d0
      do  i = 2, n
         df(i) = i*df(i-2)
      enddo   
      end


