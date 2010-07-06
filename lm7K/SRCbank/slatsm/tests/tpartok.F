C tests partok
      subroutine fmain
      implicit none
      character*80 s
      integer partok,recoff,reclen,catbeg,catsiz,subsiz,i0
      double precision a,b,c(4)
      s = 'TEST a=1+1 b=pi/2 ,c:pi/2,1-1/2;((5+1)*5)/2 pi*(2-1)'
      recoff = 0
      reclen = len(s)
      catbeg = 0
      catsiz = len(s)
      subsiz = len(s)
      call partk0(recoff,reclen,1,-1,catbeg,catsiz,-1,1,.false.)
      call dcopy(4,-99d0,0,c,1)

      i0 = partok(s,'a=','=',a,' ',1,4,0,0)
      print 100, i0,a
      i0 = partok(s,'b=','=',b,' ',1,4,0,0)
      print 100, i0,b
      i0 = partok(s,'c:',': ,;',c,' ',4,4,0,0)
      print 100, i0,c
      call stlibv(12,11,.false.)
      i0 = partok(s,'c:',': ,;',c,' ',4,4,0,0)
      print 100, i0,c

  100 format(i4,4f18.12)
      end
