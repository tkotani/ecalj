!$$$      subroutine getarf(iarg,strn)
!$$$      integer iarg
!$$$      character*(*) strn
!$$$      call getarg(iarg,strn)
!$$$      end

integer function nargf()
  integer :: iargc
  nargf = iargc() + 1
end function nargf

subroutine ftime(datim)!fortran-callable date and time
  character datim*(*)
  call fdate(datim)
  datim=datim(1:24) !takao. If this is not, write(6,*) gives CR at the ene of datim*26.
end subroutine ftime

!$$$#if TESTF
!$$$      program test
!$$$      integer i
!$$$      character strn*50
!$$$      i = nargf()
!$$$      print *, 'no command line arguments + 1 = ', i
!$$$      call getarf(1,strn)
!$$$      print *, '1st command line argument = ', strn
!$$$      end
!$$$#endif
!$$$#if TESTC
!$$$      subroutine fmain
!$$$      integer i
!$$$      character strn*50
!$$$      i = nargf()
!$$$      print *, 'no command line arguments + 1 = ', i
!$$$      call getarf(1,strn)
!$$$      print *, '1st command line argument = ', strn
!$$$      end
!$$$#endif
!$$$
