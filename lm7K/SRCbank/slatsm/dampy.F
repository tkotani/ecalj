#define BLAS3 1
      subroutine dampy(a,nca,nra,b,ncb,nrb,c,ncc,nrc,n,m,l)
C- matrix multiplication, adding result into destination
C ----------------------------------------------------------------
Ci Inputs:
Ci   a,nca,nra is the left matrix and respectively the spacing
Ci      between elements in adjacent columns and rows.
Ci   b,ncb,nrb is the right matrix and respectively the spacing
Ci      between elements in adjacent columns and rows.
Ci   c,ncc,nrc is the product matrix and respectively the spacing
Ci      between elements in adjacent columns and rows.
Ci   n,m: the number of rows and columns, respectively, to calculate
Ci   l:   length of vector for matrix multiply
Co Outputs:
Co   product matrix stored in c
Cr Remarks:
Cr   Adapted from dmpy.
Cr   This is a general-purpose matrix multiplication routine,
Cr   multiplying a subblock of matrix a by a subblock of matrix b.
Cr   Normally matrix nc{a,b,c} is the row dimension of matrix {a,b,c}
Cr   and nr{a,b,c} is 1.  Reverse nr and nc for a transposed matrix.
Cr   Arrays are locally one-dimensional so as to optimize inner loop,
Cr   which is executed n*m*l times.  No attempt is made to optimize
Cr   the outer loops, executed n*m times.
Cr     Examples: product of (n,l) subblock of a into (l,m) subblock of b
Cr   call dmpy(a,nrowa,1,b,nrowb,1,c,nrowc,1,n,m,l)
Cr     nrowa, nrowb, and nrowc are the leading dimensions of a, b and c.
Cr     To generate the tranpose of that product, use:
Cr   call dmpy(a,nrowa,1,b,nrowb,1,c,1,nrowc,n,m,l)
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nca,nra,ncb,nrb,ncc,nrc,n,m,l
      double precision a(0:*), b(0:*), c(0:*)
C Local parameters
      double precision ar
      integer i,j,k,nccj,ncbj
#if APOLLO | HP
      double precision vec_$ddot_i
#endif
#if BLAS3
      integer lda,ldb
      character*1 transa,transb
#endif

#if TIMING
      call dmpytm(0)
#endif

C      if (min(nra,nca,nrb,ncb,nrc,ncc,n,m,l) .lt. 0)
C     .  call rx('dmpy: bad integer argument')

#if BLAS3
      if (nra .eq. 1) then
        lda = nca
        transa = 'n'
      elseif (nca .eq. 1) then
        lda = nra
        transa = 't'
      else
        lda = -1
      endif
      if (nrb .eq. 1) then
        ldb = ncb
        transb = 'n'
      elseif (ncb .eq. 1) then
        ldb = nrb
        transb = 't'
      else
        ldb = -1
      endif
      if (min(lda,ldb) .lt. 0 .or. nrc .ne. 1) goto 10
#if PARALLEL
      call pp_$dgemm(transa,transb,n,m,l,1d0,a,lda,b,ldb,1d0,c,ncc)
#else
      call dgemm(transa,transb,n,m,l,1d0,a,lda,b,ldb,1d0,c,ncc)
#endif
#if TIMING
      call dmpytm(1)
#endif
      return
   10 continue
#endif

#if APOLLO | HP
      do  200  i = n-1, 0, -1
        do  200  j = m-1, 0, -1
          c(i*nrc+j*ncc) = c(i*nrc+j*ncc)
     .    + vec_$ddot_i(a(nra*i),nca,b(ncb*j),nrb,l)
  200 continue
#else

C --- Do multiplication ---
      do  20  k = l-1, 0, -1
        do  20  i = n-1, 0, -1
          ar = a(nra*i + nca*k)
          if (ar .eq. 0) goto 20
#if BLAS
          call daxpy(m,ar,b(nrb*k),ncb,c(nrc*i),ncc)
#else
          nccj = -ncc + nrc*i
          ncbj = -ncb + nrb*k
          do  25  j = m-1, 0, -1
            nccj = nccj + ncc
            ncbj = ncbj + ncb
            c(nccj) = c(nccj) + ar*b(ncbj)
   25     continue
#endif
   20 continue
#endif
#if TIMING
      call dmpytm(1)
#endif
      end

