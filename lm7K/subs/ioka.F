      integer function ioka()
C     implicit none
C- returns 0, or 1000 if OKA is defined
#if OKA
      ioka = 1000
#else
      ioka = 0
#endif
      end

