      subroutine fmain
c      program test
      implicit none
      integer i,n,m
      double precision z(1000),w(1000),f,fx,x,x0,x1
      fx(x) = (m+1)*x**m/2

    5 print *, 'enter n (no of points) and m (integrate (m+1)x**m/2)'
      read(*,*) n,m
      call mklegw(n,z,w,0)
      f = 0
      do  10  i = 1, n
   10 f = f + w(i)*fx(z(i))
      print *, 'integral I and 1-I are ', f, 1-f

      print *, 'enter x0,x1 to integrate between those boundaries:'
      read(*,*) x0,x1

      call gausq(n,x0,x1,z,w,0,000)
      f = 0
      do  20  i = 1, n
   20 f = f + w(i)*fx(z(i))
      print 333, f, f - (x1**(m+1) - x0**(m+1))/2
  333 format('I, abs error:', f25.15,g20.10)

      end

