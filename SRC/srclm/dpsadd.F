      subroutine dpsadd(adest,asrc,nel,n1,n2,fac)
C- shift and add. nel=number of elements, n1,n2= start in asrc,adest
C     implicit none
      integer n1,n2,i,iadd,ntop,nel
      double precision asrc(1),adest(1),fac
c No CRAY takao 2010jan BLAS is assumed
c #if CRAY
c      call scopy(nel,asrc(n2),1,adest(n1),1)
c      call saxpy(nel,fac,asrc(n2),1,adest(n1),1)
c #elif BLAS
c#if BLAS
      call daxpy(nel,fac,asrc(n2),1,adest(n1),1)
c#else
c      iadd = n1-n2
c      ntop = n2+nel-1
c      do  10  i = n2, ntop
c   10 adest(i+iadd) = adest(i+iadd) + fac*asrc(i)
c#endif
      end
C      program test
C      implicit none
C      double precision from(10), to(10)
C      integer i
C
C      do  10  i = 1, 10
C      to(i) = 100+i
C   10 from(i) = i
C
C      call dpsadd(to,from,4,3,5,2d0)
C      print 333, to
C  333 format(10f7.2)
C      end

