      subroutine icopy(n,ii,io)
c 92.03.15
c copy an integer array of size n ii into io
      implicit real*8 (a-h,o-z)
      dimension ii(n),io(n)
      do       i = 1,n
        io(i)    = ii(i)
      end do
      end






















