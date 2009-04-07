      integer i,j,n,ier
      complex*16 a(4,4),z(4,4),wk(32),c(4,4),scale
      double precision e(4)
      
      call sa(a)
      call eigch(a,4,11,e,z,4,wk,ier)

      print *, 'eigenvalues: (see imsl)'
      print 333, e
      print *, 'eigenvectors:'
      print 333, z
  333 format(4f14.5)

      call sa(a)
      n = 4
      call zmpy(a,2*n,2,1,z,2*n,2,1,c,2*n,2,1,n,n,n)

      print *, 'H Z  -  E Z'
      do  10  i = 1, 4
      do  10  j = 1, 4
   10 print 333, c(i,j) - e(j)*z(i,j)

      end

      subroutine sa(a)
      complex*16 a(4,4)
      a(1,1) = 3
      a(2,2) = 3
      a(3,3) = 1
      a(4,4) = 1
      a(2,1) = 1
      a(3,1) = 0
      a(4,1) = (0d0,-2d0)
      a(3,2) = (0d0,+2d0)
      a(4,2) = 0
      a(4,3) = 1
      a(1,2) = 1
      a(1,3) = 0
      a(1,4) = (0d0,+2d0)
      a(2,3) = (0d0,-2d0)
      a(2,4) = 0
      a(3,4) = 1

      end

