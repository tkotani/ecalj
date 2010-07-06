      subroutine ygemm(transa,transb,m,n,k,alpha,a,ofai,lda,b,ofbi,
     .ldb,beta,c,ofci,ldc)
C- Analog of zgemm, using real arithmetic
C ----------------------------------------------------------------
Ci Inputs:
Ci   transa: specifies the form of op( A ) to be used in multiplication
Ci           'N' or 'n',  op( A ) = A.
Ci           'T' or 't',  op( A ) = A', with A' = transpose A
Ci           'C' or 'c',  op( A ) = conjg( A' ).
Ci   transb: specifies the form of op( B ) to be used in multiplication
Ci           'N' or 'n',  op( B ) = B.
Ci           'T' or 't',  op( B ) = B', with B' = transpose B
Ci           'C' or 'c',  op( B ) = conjg( B' ).
Ci   m     :number of rows of the matrices  op( A )  and  C
Ci   n     :number of columns of the matrices  op( B )  and  C
Ci   k     :The length of the inner product
Ci   alpha :scaling factor scaling multiplication; see Outputs
Ci   a     :first matrix in multiplication
Ci   ofai  :offset to imaginary part of a
Ci   lda   :leading dimension of a
Ci   b     :second matrix in multiplication
Ci   ofbi  :offset to imaginary part of b
Ci   ldb   :leading dimension of b
Ci   beta  :amount of c to add into result; see Outputs
Ci   c     :to be added to matrix product; see Outputs
Ci          When beta is zero, c need not be set on input
Ci   ofci  :offset to imaginary part of c
Ci   ldc   :leading dimension of c
Co Outputs:
Co   c     :overwritten by alpha*op( A )*op( B ) + beta*C
C Remarks:
Cr   This is a front end for yygemm, which which is passed the
Cr   real and imaginary parts of arrays.
Cr
Cr   Unless to tell it not to, yygemm may perform the multiplication
Cr   in three real multiplications.  In doing so, it temporarily
Cr   overwrites ar,ai,br,bi with linear combinations eg ar+ai.
Cr   This means that the input arrays are NOT restored exactly intact
Cr   and some precision may be lost.
Cr
Cr   You can make yygemm do the multiplication in the normal 4 real
Cr   multiplications; this leaves untouched ar,ai,br,bi.  To accomplish
Cr   this, make transa two characters long, and put '4' in the second
Cr   character.  NB  if a and b point to the same address,
Cr   you MUST put '4' the second character of transa.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      character*(*)      transa, transb
      integer            m, n, k, lda, ldb, ldc, ofai, ofbi, ofci
      double precision   alpha, beta
      double precision   a(lda,1), b(ldb,1), c(ldc,1)

      call yygemm(transa,transb,m,n,k,alpha,a,a(1+ofai,1),lda,
     .b,b(1+ofbi,1),ldb,beta,c,c(1+ofci,1),ldc)

      end

