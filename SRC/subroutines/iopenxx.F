      integer function iopenxx (nam)
      character*(*) nam
      integer(4),save :: ix !=10
      open(newunit=ix,file=nam)
      iopenxx=ix
      end
