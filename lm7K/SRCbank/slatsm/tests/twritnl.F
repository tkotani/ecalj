      subroutine fmain
      implicit none
      character*80 strn
      integer i1,i2,i1mach,mxlen,pad1,pad2
      strn = 'Now is the time for all good men to come to '//
     .  'the aid of their country'
      call skpblb(strn,len(strn),i2)
      i2 = i2+1
      i1 = 1
      mxlen = 29
      pad1 = 1
      pad2 = 5
      call writnl(strn,i1,i2,mxlen,' ',pad1,pad2,i1mach(2))

      print *, strn(i1:i2)
      end
