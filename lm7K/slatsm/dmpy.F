#define BLAS3 1
      subroutine dmpy(a,nca,nra,b,ncb,nrb,c,ncc,nrc,n,m,l)
C- matrix multiplication
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
      integer i,j,k,nccj,nrci,nrcicj,ncbj
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
      if (min(lda,ldb) .lt. 0 .or. nrc .ne. 1) goto 11
#if PARALLEL
      call pp_$dgemm(transa,transb,n,m,l,1d0,a,lda,b,ldb,0d0,c,ncc)
#else
      call dgemm(transa,transb,n,m,l,1d0,a,lda,b,ldb,0d0,c,ncc)
#endif
#if TIMING
      call dmpytm(1)
#endif
      return
   11 continue
#endif

#if CRAY
      call mxma(a,nra,nca,b,nrb,ncb,c,nrc,ncc,n,l,m)
#elif APOLLO | HP
      do  200  i = n-1, 0, -1
        do  200  j = m-1, 0, -1
        c(i*nrc+j*ncc) = vec_$ddot_i(a(nra*i),nca,b(ncb*j),nrb,l)
C The equivalent in-line code ...
c        double precision sum
c        sum = 0
c        nakpi = nra*i
c        nbjpk = ncb*j
c        do  210  k = l-1, 0, -1
c          sum = sum + a(nakpi)*b(nbjpk)
c          nakpi = nakpi + nca
c          nbjpk = nbjpk + nrb
c  210   continue
  200 continue
#else

C --- Initialize array to zero ---
      do  10  i = n-1, 0, -1
        nrci = nrc*i
        nccj = -ncc
        do  10  j = m-1, 0, -1
        nccj = nccj + ncc
        nrcicj = nrci + nccj
        c(nrcicj) = 0
   10 continue

C --- Do multiplication ---
      do  20  k = l-1, 0, -1
        do  20  i = n-1, 0, -1
        ar = a(nra*i + nca*k)
        if (ar .eq. 0) goto 20
#if BLAS
        call daxpy(m,ar,b(nrb*k),ncb,c(nrc*i),ncc)
#else
        nrci = nrc*i
        nccj = -ncc
        ncbj = -ncb + nrb*k
        do  15  j = m-1, 0, -1
          nccj = nccj + ncc
          ncbj = ncbj + ncb
          nrcicj = nrci + nccj
          c(nrcicj) = c(nrcicj) + ar*b(ncbj)
   15   continue
#endif
   20 continue
#endif
#if TIMING
      call dmpytm(1)
#endif
      end
#if TIMING
      subroutine dmpytm(i)
      implicit none
      integer i
      integer i1mach,iprint
      double precision xx,nettim
      save nettim
      data nettim /0d0/

      if (iprint() .lt. 52) return
#if DEBUG
      if (i .eq. 0) then
        call cpudel(6,'Enter dmpytm',xx)
      else
        call cpudel(-1,'Enter dmpytm',xx)
      endif
#else
      call cpudel(-1,' ',xx)
#endif
      if (i .eq. 0) return
      nettim = nettim+xx
      print 333, xx, nettim
  333 format(22x,'Done dmpy:  time:',g10.3,' dmpy t:',g10.3)
      end
#endif

