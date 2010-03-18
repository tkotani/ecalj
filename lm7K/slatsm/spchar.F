ctakao this routine looks buggy (sch.. is char*2). But it worked. So I keep this as it is. mar2009 
      subroutine spchar(i,ch)
C- Returns special characters
Ci Inputs
Ci   i=1  return backward slash in ch
Ci   i=2  return tab in ch
Ci   i=2  return newline in ch
Ci Outputs
Co   ch
C     implicit none
      integer i
      character*(1) ch
      character*1 sch*1,tch*1,nch*1

c #if CRAY
c      data sch /'\'/ tch /'	'/  there was nch /'^M'/
c       --> but nch=ctrl+M ---> this caused compilar error in machines.
c #else
      data sch /'\\'/ tch /'\t'/  nch /'\n'/
c #endif

      if (i .eq. 1) ch = sch
      if (i .eq. 2) ch = tch
      if (i .eq. 3) ch = nch

      end
C$$$#if TEST
C$$$      character ch*1
C$$$      call spchar(1,ch)
C$$$      print 333, 'backslash',ch
C$$$      call spchar(2,ch)
C$$$      print 333, 'tab',ch
C$$$      call spchar(3,ch)
C$$$      print 333, 'newline',ch
C$$$  333 format(a,' :',a1,':')
C$$$      end
C$$$#endif
