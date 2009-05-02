      integer function iopenxx (nam)
      character*(*) nam
      integer(4),save :: ix=10
      ix=ix+1
c      print *, ' iopenxx ix=',ix
      open(ix,file=nam)
      iopenxx=ix
      end
