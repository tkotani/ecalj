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

#if CRAY
      data sch /'\'/ tch /'	'/  nch /''/
#else
      data sch /'\\'/ tch /'\t'/  nch /'\n'/
#endif

      if (i .eq. 1) ch = sch
      if (i .eq. 2) ch = tch
      if (i .eq. 3) ch = nch

      end
#if TEST
      character ch*1
      call spchar(1,ch)
      print 333, 'backslash',ch
      call spchar(2,ch)
      print 333, 'tab',ch
      call spchar(3,ch)
      print 333, 'newline',ch
  333 format(a,' :',a1,':')
      end
#endif
