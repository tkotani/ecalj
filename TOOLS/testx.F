c need interface for gfortan--> See ecalj .
      program aaa
      call abc(b)
      end
      subroutine abc(b,ax)
      real(8),optional::ax
      print *,'pxxx=',present(ax)
      end
