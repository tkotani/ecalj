      subroutine fmain
      integer nargf
      character*40 strn

      nargs = nargf()
      print *, 'nargf found', nargs, ' arguments'
      do  10  iarg = 0, nargs-1
        call getarf(iarg,strn)
        print *, iarg, ' ', strn
   10 continue
      end
