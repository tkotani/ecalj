      integer function iopenxx (nam)
      character*(*) nam
      integer(4),save :: ix !=10
c      ix=ifile_handle()
c      print *, ' iopenxx ix=',ix
      open(newunit=ix,file=nam)
      iopenxx=ix
      end
