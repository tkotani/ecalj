      subroutine yspr2a(mode,i1,i2,j1,j2,ap,ofapi,ldap,ija,offs,a,ofai,
     .  lda)
C- Unpacks a subblock of a complex sparse matrix into standard form
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit not used now
Ci          0  a is real
Ci          1  only the diagonal blocks of a are complex
Ci         >1  a is complex
Ci   i1,i2 :generate a(i,j) for row subblocks i = i1..i2
Ci   j1,j2 :generate a(i,j) for col subblocks j = j1..j2
Ci   ap    :sparse matrix, stored in block form by rows.
Ci          ap consists of a vector of matrix subblocks:
Ci          apr(*,*,i) = real part of matrix subblock i
Ci          api(*,*,i) = imaginary part of matrix subblock i
Ci   ofapi :number of elements separating real, imaginary parts of ap
Ci   ldap  :leading dimension of ap
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
Ci   ofai  :number of elements separating real, imaginary parts of a
Ci   lda   :leading dimension of ar,ai
Co Outputs
Co   a     :dense (normal) form of matrix subblock a(i1..i2,j1..j2)
Cr Remarks
Cr   The storage convention for packed form of matrices follows that of
Cr   Numerical Recipes (which see for further description), except that
Cr   here entries are themselves matrix subblocks, instead of scalars.
Cr   The matrix subblocks are storage in apr,api.
Cr   offs holds the information about the matrix subblock sizes.
Cr   ija(2,*) hold pointers the appropriate subblock in apr,api
C ----------------------------------------------------------------------
C     implicit none
      integer mode,i1,i2,j1,j2,ldap,lda,ija(2,*),offs(i2),ofapi,ofai
      double precision ap(ldap),a(lda)

      call yysp2a(mode,i1,i2,j1,j2,ap,ap(1+ofapi),ldap,ija,offs,
     .  a,a(1+ofai),lda)

      end

