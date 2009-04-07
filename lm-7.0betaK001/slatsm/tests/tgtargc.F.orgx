      subroutine fmain
      implicit none
      integer n,i,nargc
      character strn*15

      print *, 'hello, world'
      n = nargc()
      print *, 'nargc returned ', n, ' arguments'
      do  i = 0, n-1
C       call getarg(i,strn)
C       print *, i, '(getarg) ', strn
        call gtargc(i,strn)
        print *, i, '(gtargc) ', strn
      enddo

      call cexit(-1,1)
C     call rx0('done')

      end
