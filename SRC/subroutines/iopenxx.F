      integer function iopenxx (nam)
      character*(*) nam
      integer(4),save :: ix !=10
      ix=ifile_handle()
c      print *, ' iopenxx ix=',ix
      open(ix,file=nam)
      iopenxx=ix
      end
