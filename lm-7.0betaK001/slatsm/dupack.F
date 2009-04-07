      subroutine dupack(a,na)
C- Unpacks a double precision matrix
C ----------------------------------------------------------------
Ci Inputs
Ci   a,na
Co Outputs
Co   a is packed so that a_ij(i<j) is found by i + j(j+1)/2,
Co   counting from zero; i+j(j-1)/2, counting from 1 and a_ij(i>j) is
Co   obtained from a_ji
Cr Remarks
Cr   
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer na
      double precision a(0:na-1,0:na-1)
C Local parameters
      integer ii,j,i

C Copy backwards into the upper triangle first
      ii = (na*(na+1))/2
      do  100 j = na-1, 0, -1
      do  100 i = j, 0, -1
      ii = ii-1
      a(i,j) = a(ii,0)
  100 continue

C Copy to lower triangle
      do  200  i = 0, na-1
      do  200  j = 0, i-1
  200 a(i,j) = a(j,i)

c      print *, 'exit dupack: i,j='
c      print 333, a
      end

