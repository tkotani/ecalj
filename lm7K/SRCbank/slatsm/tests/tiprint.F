      subroutine fmain
      print *, 'test iprint, pshprt, popprt ...'
      print 333, iprint()
      call pshprt(1)
      call pshprt(2)
      call pshprt(3)
      call pshprt(4)
      call pshprt(5)
      call pshprt(6)
      call pshprt(7)
      print 333, iprint()
  333 format(7i4)
      print 333, (iprt(i), i=0,6)
      do 10 i = 1,6
      call popprt
      print 333, iprint()
   10 continue
      print 333, (iprt(i), i=0,5)
      print *, 'test setpr ...'
      call setpr(11)
      print 333, iprint(), (iprt(i), i=0,5)
      do 12 j = 1,6
      call popprt
      print 333, iprint(), (iprt(i), i=0,5)
   12 continue
      print *, 'test iprt, sprt ...'
      print 333, (iprt(i), i=0,5)
      call sprt(3,13)
      print 333, (iprt(i), i=0,5)
      call sprt(4,14)
      print 333, (iprt(i), i=0,5)
      call sprt(0,10)
      print 333, (iprt(i), i=0,5)
      call sprt(9,19)
      print 333, (iprt(i), i=0,5)
      call sprt(-1,99)
      print 333, (iprt(i), i=0,5)

      end
