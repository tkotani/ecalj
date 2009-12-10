      subroutine fmain
C     implicit none
      character*80 strn
      integer nw,ifblk,iw
      integer i2,j1,j2

C ... test cases at end of word
      j1 = 4
      call nwordg('bz n',0,'; ',1,j1,j2)
      print '(2i4)', j1,j2
      j1 = 5
      call nwordg('bz n',0,'; ',1,j1,j2)
      print '(2i4)', j1,j2

      strn = ' Now ; is;the:time; for;all:good:Men:to;come;to;'//
     .  'the;aid:of:their:country'
      ifblk = 1
      call skpblb(strn,len(strn),i2)
      i2 = i2+1
      print '(1x,a)', strn(1:i2)

      call wrdsg(strn(1:i2),ifblk,':;',nw)
      print '(1x,a,i4)', '------------- nw=',nw
      do  10  iw = 1, nw
        call wordg(strn(1:i2),ifblk,':;',iw,j1,j2)
        print '(3i3,''|'',a,''|  term='',a)', iw,j1,j2,
     .  strn(j1:j2),strn(j2+1:j2+1)
   10 continue

      call wordg(strn(1:i2),ifblk,':;',nw+1,j1,j2)
      print '(2i4)', j1,j2

      call wrdsg(strn(1:i2),10+ifblk,'A-Za-z',nw)
      print '(1x,a,i4)', '------------- nw=',nw
      do  20  iw = 1, nw
        call wordg(strn(1:i2),10+ifblk,'A-Za-z',iw,j1,j2)
        print '(3i3,''|'',a,''|  term='',a)', iw,j1,j2,
     .    strn(j1:j2),strn(j2+1:j2+1)
   20 continue

      call wordg(strn(1:i2),ifblk,'A-Za-z',nw+1,j1,j2)
      print '(2i4)', j1,j2

      strn = 'abc565'
      call wordg(strn,100,'56',1,j1,j2)
      print '(2i4,a,a)', j1,j2, ' str = ',strn(j1:j2)
      call wordg(strn,100,'56',2,j1,j2)
      print '(2i4,a,a)', j1,j2, ' str = '
      call wordg(strn,100,'56',3,j1,j2)
      print '(2i4,a,a)', j1,j2, ' str = '
      call wordg(strn,100,'56',4,j1,j2)
      print '(2i4,a,a)', j1,j2, ' str = '


C ... this needs to be examined
C      call wordg(strn,110,'56',1,j1,j2)
C      print *, j1,j2, ' str = ',strn(j1:j2)
C      call wordg(strn,110,'56',2,j1,j2)
C      print *, j1,j2, ' str = '
C      call wordg(strn,110,'56',3,j1,j2)
C      print *, j1,j2, ' str = '
C      call wordg(strn,1103,'56',4,j1,j2)
C      print *, j1,j2, ' str = '

      end
