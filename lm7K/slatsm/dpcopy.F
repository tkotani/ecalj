      subroutine dpcopy(afrom,ato,n1,n2,fac)
C- Copy and scale a portion of a vector
C     implicit none
      integer n1,n2,i
      double precision afrom(1),ato(1),fac
      if (fac .ne. 1d0) goto 100
#if CRAY
      call scopy(n2-n1+1,afrom(n1),1,ato(n1),1)
#elif APOLLO | HP
        call vec_$dcopy(afrom(n1),ato(n1),n2-n1+1)
#elif BLAS
      call dcopy(n2-n1+1,afrom(n1),1,ato(n1),1)
#else
        do  10 i = n1, n2
   10   ato(i) = afrom(i)
#endif
        return

C --- fac not unity ---
  100 continue
      do  110  i = n1, n2
  110 ato(i) = fac*afrom(i)
      end
C      program test
C      implicit none
C      double precision from(10), to(10)
C      integer i
C
C      do  10  i = 1, 10
C      to(i) = 0
C   10 from(i) = i
C
C      call dpcopy(from,to,3,7,2d0)
C      print 333, to
C      call dpcopy(from,to,3,7,1d0)
C      print 333, to
C  333 format(10f7.3)
C      end

