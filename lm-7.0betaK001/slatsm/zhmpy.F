       subroutine zhmpy(a,nca,nra,nia,b,ncb,nrb,nib,c,ncc,lc,l)
C- complex matrix multiplication: result assumed hermitian
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
Ci   l:   length of vector for matrix multiply
Co Outputs:
Co   product matrix stored added into c(1:ncc,1:ncc)
Cr Remarks:
Cr   Adapted from zampy; calculates only for the lower triangular part
Cr   This version uses vectorizable BLAS-style daxpy loops.
C ----------------------------------------------------------------
C     implicit none
      logical lc
      integer nca,nra,nia,ncb,nrb,nib,ncc,nic,l
      double precision a(0:*), b(0:*), c(0:*)
      integer i,j,k,nccj,ncbj,nrcicj
      double precision ar,ai

      nic = ncc**2

C --- Do multiplication ---
      do  20  k = l-1, 0, -1
        do  20  i = ncc-1, 0, -1
#if BLAS
        call yyaxpy(i+1,a(nra*i+nca*k),a(nia+nra*i+nca*k),
     .    b(nrb*k),b(nib+nrb*k),ncb,c(i),c(nic+i),ncc,.true.)
#else
        ar = a(      nra*i + nca*k)
        ai = a(nia + nra*i + nca*k)
        nccj = -ncc
        ncbj = -ncb + nrb*k
        do  20  j = i, 0, -1 
        nccj = nccj + ncc
        ncbj = ncbj + ncb
        nrcicj = i + nccj
        c(nrcicj)     = c(nrcicj)     + ar*b(ncbj) - ai*b(nib+ncbj)
        c(nic+nrcicj) = c(nic+nrcicj) + ar*b(nib+ncbj) + ai*b(ncbj)
#endif
   20 continue

C --- Copy c into c+ ---
      if (.not. lc) return
      do  40  i = ncc-1, 0, -1
        do  40  j = ncc-1, i+1, -1 
        c(i+ncc*j) = c(j+ncc*i)
        c(nic+i+ncc*j) = -c(nic+j+ncc*i)
   40 continue

      end

