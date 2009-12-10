      integer function itrpck(i,j,n)
C- Store dho in packed upper triangular form
C ----------------------------------------------------------------------
Ci Inputs:
Ci  i,j,n 
Co Outputs:
Co  an index
Cr Remarks
Cr  To store a symmetic n x n matrix with zero diagonal elements
Cr  an array of length n*(n-1)/2 is needed. To pack or to access the
Cr  (i,j) element of this array the pointer to the packed array is
Cr        n(i-1) - i(i-1)/2 + j - i (j > i)
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer i,j,n

      if (i .ge. j) then
        call rx(' ITRPCK: i>=j')
      else
        itrpck = (i-1)*n - i*(i-1)/2 + j - i
      endif
      end
