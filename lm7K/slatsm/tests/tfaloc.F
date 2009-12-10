C this tests whether a compiler can handle pointers, and routine faloc.
      subroutine fmain
      implicit none
      integer nwords,i
      parameter (nwords=20)
      pointer (iarr, arr), (iarr2, arr2)
      double precision arr(20),arr2(20)
      pointer (ptr,v), (ptr2,v2)
      character a*12, v*12, z*1, v2*12
      data a /'abcdefghijkl'/

C      ptr = loc(a)
C      ptr = ptr + 4
C      ptr2 = malloc(12)
C      v2 = a
C      z = v(1:1)
C      print *, z
C      z = v2(5:5)
C      print *, z

      call faloc(iarr, 4, nwords)
      do  10  i = 1, 20
   10 arr(i) = i
      print 333, arr
  333 format(10f8.4)
      iarr2 = iarr
      print *
      print *, arr2(2), arr2(20)

      end
