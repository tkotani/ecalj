c$$$      subroutine getarf(iarg,strn)
c$$$      integer iarg
c$$$      character*(*) strn
c$$$      call getarg(iarg,strn)
c$$$      end

      integer function nargf()
      integer iargc
      nargf = iargc() + 1
      end

      subroutine ftime(datim)!fortran-callable date and time
      character datim*(*)
      call fdate(datim)
      datim=datim(1:24) !takao. If this is not, write(6,*) gives CR at the ene of datim*26.
      end

c$$$#if TESTF
c$$$      program test
c$$$      integer i
c$$$      character strn*50
c$$$      i = nargf()
c$$$      print *, 'no command line arguments + 1 = ', i
c$$$      call getarf(1,strn)
c$$$      print *, '1st command line argument = ', strn
c$$$      end
c$$$#endif
c$$$#if TESTC
c$$$      subroutine fmain
c$$$      integer i
c$$$      character strn*50
c$$$      i = nargf()
c$$$      print *, 'no command line arguments + 1 = ', i
c$$$      call getarf(1,strn)
c$$$      print *, '1st command line argument = ', strn
c$$$      end
c$$$#endif
c$$$
