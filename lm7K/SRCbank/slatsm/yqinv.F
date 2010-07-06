      subroutine yqinv(cs,a,ofai,lda,nlev,n,w,ldw,ierr)
C- Inversion of a complex matrix using Strassen's algorithm
C ----------------------------------------------------------------
Ci Inputs:
Ci   cs:   :if 'h', a is assumed hermitian
Ci   a     :matrix to be inverted
Ci   ofai  :number of elements separating real, imaginary parts of a
Ci   lda   :leading dimension of a
Ci   n     :rank of the matrix to be inverted
Ci   nlev  :the maximum number of recursion levels allowed.
Ci          To avoid roundoff errors, nlev=2 is suggested.
Ci   w     :double precision work array of dimension ldw*(n+1)
Ci   ldw   :leading dimension of w
Co Outputs:
Co   a     :is overwritten by inverse of input a
Co   ierr  :returned nonzero if matrix was not fully inverted.
Cb Bugs
Cb   The algorithm fails if a22 is singular, even if a is not.
Cb   Similarly, if smaller subblocks are singular, yyqinv may fail
Cb   when called recursively.
Cr Remarks:
Cr   This is a front end for yyqinv, which which is passed the
Cr   real and imaginary parts of arrays.
C ----------------------------------------------------------------
C     implicit none
      character*(*) cs
      integer n,lda,ldw,ierr,nlev,ofai
      double precision a(lda),w(ldw,1)

      call yyqinv(cs,a,a(1+ofai),lda,nlev,n,w,ldw,ierr)

      end

