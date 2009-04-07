      subroutine fmain
      implicit none
      integer i1,i2,bitand,bitor,bitlow,i

      i1 = (1+2+8+32+128+512)*512 + (1+2+16+32+256)
      i2 = (1+2+8+32+128+256+512)*256 + (1+2+8+16)
      
      i = bitand(i1,i2)
      print *, i1,i2,i
      i = bitor(i1,i2)
      print *, i1,i2,i
      i1 = 764*2
      i = bitlow(i1)
      print *, i1,i
      print *, bitlow(0)
      end
