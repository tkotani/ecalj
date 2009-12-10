      subroutine yspbmm(mode,i1,i2,j1,j2,nc,a,ofai,lda,ija,offs,
     .  x,ofxi,ldx,b,ofbi,ldb)
C- Complex sparse-block-matrix dense-matrix multiply
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci          0  multiply contents into b
Ci          1  add product  a x to existing contents of b
Ci          2  add product -a x to existing contents of b
Ci          10s digit
Ci          0  a is real
Ci          1  only the diagonal parts of a are complex
Ci          2  a is complex
Ci   i1,i2 :calculate a(i,j)*x(j) for row subblocks i = i1..i2
Ci   j1,j2 :calculate a(i,j)*x(j) for col subblocks j = j1..j2
Ci   nc    :number of columns in x and the result matrix b
Ci   a     :sparse matrix, stored in block form by rows.
Ci          a consists of a vector of matrix subblocks:
Ci          a(*,*,i) = matrix subblock i
Ci   ofai  :number of elements separating real, imaginary parts of a
Ci   lda   :leading dimension of a
Ci   ija   :column index packed array pointer data to array a
Ci         ija(1,*) follows essentially the same conventions
Ci         as for scalar packed arrays (see da2spr)
Ci         except that indices now refer to matrix subblocks.
Ci         ija(1,1)= n+2, where n = max(number of rows, number of cols)
Ci         ija(1,i), i = 1,..n+1 = points to first entry in a for row i
Ci         ija(1,i), i = n+2... column index element a(i).  Thus
Co                   for row i, k ranges from ija(i) <= k < ija(i+1) and
Co                   sum_j a_ij x_j -> sum_k a_(ija(2,k)) x_(ija(1,k))
Ci         ija(2,*)  pointers to the matrix subblocks blocks in a:
Ci         ija(2,i), i=1..n  pointers to blocks on the diagonal of a
Ci         ija(2,i), i=n+2.. pointers to elements of a, grouped by rows
Ci   offs  :offsets to first entries in matrix subblocks
Ci          Thus the dimension of row i = offs(i+1) - offs(i)
Ci          If a consists of scalar subblocks, offs(i) = i-1.
Ci   ofbi  :number of elements separating real, imaginary parts of b
Ci   x     :dense matrix, and second operand
Ci   ofxi  :number of elements separating real, imaginary parts of x
Ci   ldx   :leading dimension of x
Co Outputs
Co   b     :result matrix
Co   ldb   :leading dimension of b
Cr Remarks
Cr   This routine multiplies a sparse matrix whose elements
Cr   are matrix subblocks, by a dense matrix.
Cu Updates
Cb Bugs
Cb   Never checked for mode=2
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer i1,i2,j1,j2,nc,lda,ldb,ldx,ija(2,*),offs(i2),mode
      integer ofai,ofbi,ofxi
      double precision a(lda),x(ldx),b(ldb)

      call yysbmm(mode,i1,i2,j1,j2,nc,a,a(1+ofai),lda,ija,offs,
     .  x,x(1+ofxi),ldx,b,b(1+ofbi),ldb)
      end

