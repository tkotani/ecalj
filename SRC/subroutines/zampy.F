#define BLAS3 1
      subroutine zampy(a,nca,nra,nia,b,ncb,nrb,nib,c,ncc,nrc,nic,n,m,l)
C- complex matrix multiplication: matrix incremented, not overwritten 
C ----------------------------------------------------------------
Ci Inputs:
Ci   a,nca,nra is the left matrix and respectively the spacing
Ci      between column elements and row elements and between the
Ci      real and imaginary parts of a (in real words)
Ci   b,ncb,nrb is the right matrix and respectively the spacing
Ci      between column elements and row elements and between the
Ci      real and imaginary parts of b (in real words)
Ci   c,ncc,nrc is the product matrix and respectively the number of
Ci      between column elements and row elements and between the
Ci      real and imaginary parts of c (in real words)
Ci   n,m: the number of rows and columns, respectively, to calculate
Ci   l:   length of vector for matrix multiply
Co Outputs:
Co   product matrix stored added into c
Cr Remarks:
Cr   This is a general-purpose matrix multiplication routine,
Cr   multiplying a subblock of matrix a by a subblock of matrix b.
Cr   Normally matrix nc{a,b,c} is the row dimension of matrix {a,b,c}
Cr   and nr{a,b,c} is 1.  Reverse nr and nc for a transposed matrix.
Cr   Arrays are locally one-dimensional so as to optimize inner loop,
Cr   which is executed n*m*l times.  No attempt is made to optimize
Cr   the outer loops, executed n*m times.
Cr     Examples: product of complex matrix c = a*b  (arrays b,c
Cr     dimensioned complex*16; a real*8 with imag following real)
Cr     call zampy(a,n,1,ndim**2,b,2*n,2,1,c,2*n,2,1,n,n,n)
Cr     To generate c = a*b^T, use:
Cr     call zampy(a,2*n,2,1,b,2,2*n,1,c,2*n,2,1,n,n,n)
Cr   This version uses vectorizable BLAS-style daxpy loops.
C ----------------------------------------------------------------
C     implicit none
      integer nca,nra,nia,ncb,nrb,nib,ncc,nrc,nic,n,m,l
      double precision a(0:*), b(0:*), c(0:*)
      integer i,j,k,nrci,nccj,ncbj,nrcicj
      double precision ar,ai,xx
#if BLAS3
      integer lda,ldb
      character*1 ta,tb
#endif

#if BLAS3
C --- Test for real and imaginary parts separated ---
      if (nra .eq. 1) then
        lda = nca
        ta = 'n'
      elseif (nca .eq. 1) then
        lda = nra
        ta = 't'
      else
        lda = -1
      endif
      if (nrb .eq. 1) then
        ldb = ncb
        tb = 'n'
      elseif (ncb .eq. 1) then
        ldb = nrb
        tb = 't'
      else
        ldb = -1
      endif
      if (min(lda,ldb) .lt. 0 .or. nrc .ne. 1) goto 10
#if TIMING
      call dmpytm(0)
#endif
      xx = 1
#if PARALLEL
      call pp_$dgemm(ta,tb,n,m,l,xx,a,lda,b,ldb,xx,c,ncc)
      call pp_$dgemm(ta,tb,n,m,l,-xx,a(nia),lda,b(nib),ldb,xx,c,ncc)
      call pp_$dgemm(ta,tb,n,m,l,xx,a,lda,b(nib),ldb,xx,c(nic),ncc)
      call pp_$dgemm(ta,tb,n,m,l,xx,a(nia),lda,b,ldb,xx,c(nic),ncc)
#else
      call dgemm(ta,tb,n,m,l,xx,a,lda,b,ldb,xx,c,ncc)
      call dgemm(ta,tb,n,m,l,-xx,a(nia),lda,b(nib),ldb,xx,c,ncc)
      call dgemm(ta,tb,n,m,l,xx,a,lda,b(nib),ldb,xx,c(nic),ncc)
      call dgemm(ta,tb,n,m,l,xx,a(nia),lda,b,ldb,xx,c(nic),ncc)
#endif
#if TIMING
      call dmpytm(1)
#endif
      return

C --- Test for complex*16 matrices ---
   10 continue
      call rx('zampy not set for c*16 yet')
      goto 15
   15 continue
#endif

C --- Do multiplication ---
      do  20  k = l-1, 0, -1
        do  20  i = n-1, 0, -1
          nrci = nrc*i
#if BLAS
          call yyaxpy(m,a(nra*i+nca*k),a(nia+nra*i+nca*k),
     .    b(nrb*k),b(nib+nrb*k),ncb,c(nrci),c(nic+nrci),ncc,.true.)
#else
          ar = a(      nra*i + nca*k)
          ai = a(nia + nra*i + nca*k)
          nccj = -ncc
          ncbj = -ncb + nrb*k
          do  22  j = m-1, 0, -1 
            nccj = nccj + ncc
            ncbj = ncbj + ncb
            nrcicj = nrci + nccj
            c(nrcicj)     = c(nrcicj)     + ar*b(ncbj) - ai*b(nib+ncbj)
            c(nic+nrcicj) = c(nic+nrcicj) + ar*b(nib+ncbj) + ai*b(ncbj)
   22     continue
#endif
   20 continue
      end

