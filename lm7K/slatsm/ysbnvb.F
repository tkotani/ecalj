      subroutine ysbnvb(cs,mode,ap,ofapi,ldap,ija,offs,nlev,
     .  nn2,nb,w,ldw,w2,a,ofai,lda,b,ofbi,ldb,ierr)
C- Solve a x = b, with a = complex block sparse matrix
C ----------------------------------------------------------------
Ci Inputs:
Ci   cs:  : a string containing any of the following characters.
Ci          't'  solve b = x a instead of a x = b
Ci          'h'  a is assumed hermitian.
Ci          'b'  Assume partial inverse for a is already performed
Ci               and solve a x = b for a new b
Ci               ar,ai must be preserved between successive calls.
Ci          '1'  Only return b1 (solution of ax=b for lower partition)
Ci               Useful if only rhs for this lower block is needed.
Ci          '2'  Skip final matrix multiplication to generate x2
Ci               (see Remarks).  Useful if only some portions of
Ci               ax=b for upper partition are needed.
Ci          '4'  Do multiplications using standard four real operations
Ci               (slower, but avoids additions and subtractions that
Ci                can reduce machine precision)
Ci   mode   :1s digit
Ci           1  only the diagonal parts subblocks of matrix are complex
Ci              (only case tested so far)
Ci           2  full matrix is complex
Ci   ap     :left hand side of equation a x = b, in packed form
Ci   ofapi  :number of elements separating real, imaginary parts of ap
Ci   ldap   :leading dimensions of apr,api
Ci   a      :a work array holding decomposition of a
Ci   ofai   :number of elements separating real, imaginary parts of a
Ci   lda    :leading dimension of ar,ai
Ci   ija    :column index packed array pointer data; see Remarks
Ci   offs   :offsets to first entries in matrix subblocks. NB: offs(1)=0
Ci           Subblock dimension of row(or col) i = offs(i+1) - offs(i)
Ci   nn1,nn2:range of subblocks which comprise the matrix to be inverted
Cr           the matrix subblock to be inverted consist of the
Cr           rows and columns  offs(nn1)+1...offs(nn2+1)
Ci   nlev   :the maximum number of recursion levels allowed in the
Ci           matrix inversion steps; see Remarks
Ci   w,ldw  :complex work array of dimension ldw*nd, nd= rank of a
Ci   w2     :a complex work array of dimension nb*nd
Ci           NB: w and w2 may use the same address space
Ci   b     :right hand side of equation a x = b
Ci   ofbi  :number of elements separating real, imaginary parts of b
Ci   ldb    :leading dimension of br,bi
Co Outputs:
Co   ar,ai  :partial matrix inverse; see Remarks
Co   ierr   :returned nonzero if matrix could not be inverted.
Co   b      :is OVERWRITTEN with a^-1 b (b a^-1 in the transpose case)
Co           NB: ysbnvb may return before completing generation of
Co           a^-1 b, depending on input cs; see Remarks.
Co   w2     :b2 - a21 x1; see Remarks
Cb Bugs
Cb   ysbnvb has not been tested for the hermitian case,
Cb   nor handle it efficiently.
Cr Remarks:
Cr   This is a front end for yysbib, which which is passed the
Cr   real and imaginary parts of arrays, which see for description
C ----------------------------------------------------------------
C     implicit none
      integer mode,nn2,nb,lda,ldb,ldw,ierr,nlev,ldap,ofapi,ofai,ofbi,
     .  ija(2,*),offs(nn2)
      double precision a(lda),b(ldb),w(ldw,*),ap(ldap),w2(*)
      character cs*(*)

      call yysbib(cs,mode,ap,ap(1+ofapi),ldap,ija,offs,nlev,
     .  nn2,nb,w,ldw,w2,a,a(1+ofai),lda,b,b(1+ofbi),ldb,ierr)
      end

