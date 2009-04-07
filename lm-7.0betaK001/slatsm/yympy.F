      subroutine yympy(ar,ai,nca,nra,br,bi,ncb,nrb,cr,ci,ncc,nrc,n,m,l)
C- complex matrix multiplication
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
Cr   This is a general-purpose complex matrix multiplication routine,
Cr   multiplying a subblock of matrix a by a subblock of matrix b,
Cr   where the real and imaginary parts are separated.
Cr  
Cr   Normally matrix nc{a,b,c} is the row dimension of matrix {a,b,c}
Cr   and nr{a,b,c} is 1.  Reverse nr and nc for a transposed matrix.
Cr
Cr   Example: product of complex matrix c = a*b  (arrays b,c
Cr   dimensioned complex*16; a real*8 with imag following real)
Cr     call yympy(a,n,1,ndim**2,b,2*n,2,1,c,2*n,2,1,n,n,n)
Cr   Example: to generate c = a*b^T, use:
Cr     call yympy(a,2*n,2,1,b,2,2*n,1,c,2*n,2,1,n,n,n)
C ----------------------------------------------------------------
C     implicit none
      integer nca,nra,ncb,nrb,ncc,nrc,n,m,l
      double precision ar(0:*),ai(0:*), br(0:*),bi(0:*), cr(0:*),ci(0:*)
      integer i,j,ii
      double precision xx

C --- ci = ai *  bi ---
      call dmpy(ai,nca,nra,bi,ncb,nrb,ci,ncc,nrc,n,m,l)

C --- cr = ar *  br ---
      call dmpy(ar,nca,nra,br,ncb,nrb,cr,ncc,nrc,n,m,l)

C --- cr <- ar*br - ai*bi,  ci <- -ar*br - ai*bi  ---
      do  10  j = 0, m-1
        do  10  i = 0, n-1
          ii = nrc*i+ncc*j
          xx = cr(ii)
          cr(ii) =  xx - ci(ii)
          ci(ii) = -xx - ci(ii)
   10 continue

C --- ci += (ar+ai)(br+bi) ---
      xx = 1
      call dmadd(ar,nca,nra,xx,ai,nca,nra,xx,ar,nca,nra,n,l)
      call dmadd(br,ncb,nrb,xx,bi,ncb,nrb,xx,br,ncb,nrb,l,m)
      call dampy(ar,nca,nra,br,ncb,nrb,ci,ncc,nrc,n,m,l)
      call dmadd(ar,nca,nra,xx,ai,nca,nra,-xx,ar,nca,nra,n,l)
      call dmadd(br,ncb,nrb,xx,bi,ncb,nrb,-xx,br,ncb,nrb,l,m)
      end

