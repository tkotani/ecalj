      subroutine fmain
      implicit none
      character*90 strn
      integer it(5),i,ip,lstrn,i1mach,parg
      double precision wt(5)

      lstrn = len(strn)
      strn = ' B30,n=4,w=2,1,wa=9,fn=mxm,wc=11,b=1; A,k=3,bv=.11,.22'
      ip = 1
      i = parg(',w=',4,strn,ip,lstrn,',; ',2,2,it,wt)
      call awrit3('i=%i  it=%2:i  wt=%2:1d',' ',80,i1mach(2),i,it,wt)
      call dcopy(4,-99d0,0,wt,1)
      ip = 1
      i = parg(',w=',4,strn,ip,lstrn,',; ',2,3,it,wt)
      call awrit3('i=%i  it=%2:i  wt=%2:1d',' ',80,i1mach(2),i,it,wt)
      call dcopy(4,-99d0,0,wt,1)
      ip = 1
      i = parg('w=',4,strn,ip,lstrn,',; ',2,3,it,wt)
      call awrit3('i=%i  it=%2:i  wt=%2:1d',' ',80,i1mach(2),i,it,wt)
      end

