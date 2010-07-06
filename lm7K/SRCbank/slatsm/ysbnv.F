      subroutine ysbnv(mode,ap,ofapi,ldap,ija,offs,nlev,
     .nn1,nn2,w,ldw,a,ofai,lda,ierr)
C- Inversion of a complex block sparse matrix using Strassen's algorithm
C ----------------------------------------------------------------
Ci Inputs:
Ci   mode   :1s digit
Ci           1 only the diagonal parts subblocks of matrix are complex
Ci             (only case tested so far)
Ci           2 full matrix is complex
Ci           3 matrix is hermitian (not tested)
Ci          10s digit fixes how array is partitioned
Ci             It is passed as first argument to psybnv
Ci           0 partitions array into approximately equal numbers
Ci             of nonzero subblocks
Ci           1 partitions array into approximately equal dimensions
Ci   nnmin  :fewest subblocks for which ysbnv calls itself recursively.
Ci   ap     :array to invert, in packed form
Ci   ofapi  :number of elements separating real, imaginary parts of ap
Ci   ldap   :leading dimensions of apr,api
Ci   ofai   :number of elements separating real, imaginary parts of a
Ci   lda    :leading dimension of ar,ai
Ci   ija    :column index packed array pointer data; see Remarks
Ci   offs   :offsets to first entries in matrix subblocks. NB: offs(1)=0
Ci           Subblock dimension of row(or col) i = offs(i+1) - offs(i)
Ci   nn1,nn2:range of subblocks which comprise the matrix to be inverted
Ci           the matrix subblock to be inverted consist of the
Ci           rows and columns  offs(nn1)+1...offs(nn2+1)
Ci   nlev   :the maximum number of recursion levels allowed.
Ci           To avoid roundoff errors, nlev=2 is suggested.
Ci   w,ldw  :complex work array of dimension ldw*n
Co Outputs:
Co   a      :inverse of matrix, stored in a(i:j,i:j) with
Co           i = offs(nn1)+1 and j = offs(nn2+1).
Co   ierr is returned nonzero if matrix could not be inverted.
Cr Remarks:
Cr  *ysbnv uses the Strassen algorithm to invert a subblock of a matrix
Cr   stored in block packed form.  The inverse is not assumed to be
Cr   sparse, and is returned in conventional form in ar,ai.  Arrays
Cr   ija,offs,apr,api contain data for the matrix; see yysp2a for a
Cr   description of block matrix storage conventions, the use of these
Cr   arrays, and how a matrix subblock may be unpacked from them into
Cr   conventional form.
Cr
Cr  *The matrix to be inverted comprises the rows (and columns)
Cr   offs(nn1)+1...offs(nn2+1) of a.  ysqnv partitions these rows and
Cr   columns into four subblocks a11,a21,a12,a22.  (See below for how
Cr   the 1 and 2 partitions are apportioned.)  Let c be the inverse,
Cr   with subblocks c11,c21,c12,c22.  Then inversion proceeds by:
Cr      (1) invert   a22 (see below for recursive inversion)
Cr      (2) generate (c11)^-1 = (a11 - a12 a22^-1 a21)
Cr      (3) invert   (c11)^-1
Cr      (4) generate c21 = -(a22^-1 a21) c11 is
Cr      (5) generate c12 = -c11 a12 a22^-1 is generated
Cr      (6) generate c22 = a22^-1 + (a22^-1 a21 c11) (a12 a22^-1)
Cr   These steps require two inversions, three sparse and three normal
Cr   multiplications.
Cr
Cr  *Partitioning into subblocks.  Partitions 1 and 2 consist of
Cr   offs(nn1)+1...offs(nm1+1) and offs(nm1+1)+1...offs(nn2+1).
Cr   nm1 is set in psybnv according to the 10s digit of mode.
Cr
Cr  *Inversion of the subblocks.  It is more efficient to allow
Cr   inversion of a22 to proceed recursively, if your compiler allows
Cr   it.  Recursion proceeds provided nlev>0, the number of subblocks
Cr   exceeds nnmin and the dimension of the matrix to be inverted
Cr   exceeds nmin.  (c11)^-1 is inverted calling yyqinv, which can
Cr   also proceed recursively.
Cb Bugs
Cb   ysbnv has not been tested for the hermitian case, 
Cb   nor handle it efficiently.
C ----------------------------------------------------------------
C     implicit none
      integer mode,nn1,nn2,lda,ldw,ierr,nlev,ldap,ija(2,*),offs(nn2)
      integer ofapi,ofai
      double precision a(lda),w(ldw,*),ap(ldap)

      call yysbnv(mode,ap,ap(1+ofapi),ldap,ija,offs,nlev,
     .nn1,nn2,w,ldw,a,a(1+ofai),lda,ierr)

      end

