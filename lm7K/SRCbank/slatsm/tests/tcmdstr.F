      subroutine fmain
C     implicit none
      character*(25) a
      integer iarg,nargf,nxargc
      logical cmdstr,lx
      character*50 s

      print *, 'to agree with test outfile, invoke with: a.out one 2 3'
      print *, 'start nargf=',nargf(),'  nxargc=',nxargc()
      s = 'xfirst xsecond ''this is xthird'''
      call acmdop(s,len(s),0)
      print *, 'added nargf=',nargf(),'  nxargc=',nxargc()

      do  10  iarg = 1, nxargc()


        a = '1234567890123456789012345678901234567890'
        lx = cmdstr(iarg,a)
        print 333, iarg, lx, a
  333   format(i4,l2,1x,a)
   10 continue

      end
