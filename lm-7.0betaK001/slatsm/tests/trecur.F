      subroutine fmain
      implicit none
      integer n
      n = 10
      call trecur(n)
      end
      subroutine trecur(nn)
c     implicit none
      integer  nn
      integer n,m

      n = nn
      m = nn/2
      print *, 'enter trecur,nn,n,m=',nn,n,m
      if (nn .gt. 1) call trecur(nn/2)
      print *, ' exit trecur,nn,n,m=',nn,n,m
      end
