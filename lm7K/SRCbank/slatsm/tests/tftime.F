C Test ftime in fsubs.c
      subroutine fmain
C      implicit none
      integer lgunit
      character*26 datim

      datim = ' '
      call ftime(datim)
      call awrit0(datim,' ',26,lgunit(1))
      end
