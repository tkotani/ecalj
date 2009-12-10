      subroutine dpack(a,na)
C- Packs a double precision matrix
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

      ii = 0
      do  100 j = 0, na-1
      do  100 i = 0, j
      a(ii,0) =  a(i,j)
      ii = ii+1
  100 continue
      end

