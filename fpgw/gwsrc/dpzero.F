      subroutine dpzero(a,n)
      real*8 a(n)
      do 1 i=1,n
        a(i)=0.d0
    1 continue
      return
      end
