      subroutine yqinvb(cs,a,ofai,lda,n,nb,w,ldw,w2,b,ofbi,ldb,ierr)
C- Solution of a x = b by vectorizable multiplications and inversions
C ----------------------------------------------------------------
Ci Inputs:
Ci   cs:   : a string containing any of the following characters.
Ci          't'  solve b = x a instead of a x = b
Ci          'h'  a is assumed hermitian.
Ci          'b'  Assume partial inverse for a is already performed.
Ci               ar,ai must be preserved between successive calls.
Ci          '1'  ignored by yyqnvb, for compatibility with ysbnvb
Ci          '2'  ignored by yyqnvb, for compatibility with ysbnvb
Ci          '4'  Do multiplications using standard four real operations
Ci               (slower, but avoids additoins and subtractions that
Ci                can reduce machine precision)
Ci   a     :lhs of equation a x = b
Ci   ofai  :number of elements separating real, imaginary parts of a
Ci   lda   :leading dimension of ar,ai
Ci   n     :solve a x = b for matrix a(1:n,1:n)
Ci   w,ldw :a double precision work array of dimension ldw*(n+1)
Ci   w2    :a double precision work array of dimension nb*(n+1)*2
Ci          NB: w and w2 may use the same address space
Ci   b     :right hand side of equation a x = b
Ci   ofbi  :number of elements separating real, imaginary parts of b
Ci   ldb   :leading dimension of br,bi
Ci   nb    :the number of columns (rows, if cs contains 't') to solve.
Co Outputs:
Co   a     :is OVERWRITTEN, into a partially decomposed form
Co   ierr  :is returned nonzero if matrix was not successfully inverted
Co   br,bi :is OVERWRITTEN with a^-1 b (b a^-1 in the transpose case)
Cb Bugs
Cb   yyqnvb fails if a22 is singular, even if a is not.
Cb   Similarly, yyqinv may fail to invert a22 even if it is not singular.
Cr Remarks:
Cr   This is a front end for yyqnvb, which which is passed the
Cr   real and imaginary parts of arrays.
C ----------------------------------------------------------------
C     implicit none
      integer n,lda,ldw,ldb,ierr,nb,ofai,ofbi
      character cs*(*)
      double precision a(lda),b(ldb),w(ldw,*),w2(nb*(n+1))

      call yyqnvb(cs,a,a(1+ofai),lda,n,nb,w,ldw,w2,b,b(1+ofbi),ldb,ierr)
      end

