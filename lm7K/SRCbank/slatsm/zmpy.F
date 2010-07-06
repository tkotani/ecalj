      subroutine zmpy(a,nca,nra,nia,b,ncb,nrb,nib,c,ncc,nrc,nic,n,m,l)
C- General complex matrix multiplication
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
Co   product matrix stored in c
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
Cr     call zmpy(a,n,1,ndim**2,b,2*n,2,1,c,2*n,2,1,n,n,n)
Cr     To generate c = a*b^T, use:
Cr     call zmpy(a,2*n,2,1,b,2,2*n,1,c,2*n,2,1,n,n,n)
Cr   This version uses vectorizable BLAS-style daxpy loops.
C ----------------------------------------------------------------
C     implicit none
      integer nca,nra,nia,ncb,nrb,nib,ncc,nrc,nic,n,m,l
      double precision a(0:*), b(0:*), c(0:*)
      integer i,j,ii
      double precision ar,xx

C --- ci = ai *  bi ---
      call dmpy(a(nia),nca,nra,b(nib),ncb,nrb,c(nic),ncc,nrc,n,m,l)

C --- cr = ar *  br ---
      call dmpy(a,nca,nra,b,ncb,nrb,c,ncc,nrc,n,m,l)

C --- cr <- ar*br - ai*bi,  ci <- -ar*br - ai*bi  ---
      if (m*n .eq. nic) then
        do  12  ii = 0, m*n-1
          ar = c(ii)
          c(ii)     =  ar - c(nic+ii)
          c(nic+ii) = -ar - c(nic+ii)
   12   continue
      else
        do  10  j = 0, m-1
          do  10  i = 0, n-1
            ii = nrc*i+ncc*j
            ar = c(ii)
            c(ii)     =  ar - c(nic+ii)
            c(nic+ii) = -ar - c(nic+ii)
   10   continue
      endif

C --- ci += (ar+ai)(br+bi) ---
      xx = 1
      call dmadd(a,nca,nra,xx,a(nia),nca,nra,xx,a,nca,nra,n,l)
      call dmadd(b,ncb,nrb,xx,b(nib),ncb,nrb,xx,b,ncb,nrb,l,m)
      call dampy(a,nca,nra,b,ncb,nrb,c(nic),ncc,nrc,n,m,l)
      call dmadd(a,nca,nra,xx,a(nia),nca,nra,-xx,a,nca,nra,n,l)
      call dmadd(b,ncb,nrb,xx,b(nib),ncb,nrb,-xx,b,ncb,nrb,l,m)
      end

