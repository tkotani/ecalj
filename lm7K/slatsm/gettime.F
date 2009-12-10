      subroutine gettime(datim)
C     implicit none
      character datim*(*)

      datim = ' '
#if SGI
      call fdate(datim)
#else
      call ftime(datim)
#endif
      end

