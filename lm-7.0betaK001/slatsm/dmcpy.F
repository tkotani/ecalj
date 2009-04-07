      subroutine dmcpy(a,nca,nra,b,ncb,nrb,n,m)
C- General matrix copy
C ----------------------------------------------------------------
Ci Inputs:
Ci   a,nca,nra is the source matrix and respectively the number of
Ci      elements separating columns and rows.
Ci   b,ncb,nrb is the dest. matrix and respectively the number of
Ci      elements separating columns and rows.
Ci   n,m: the number of columns and rows, respectively, to calculate
Co Outputs:
Co   result matrix stored in c
Cr Remarks:
Cr   This is a general-purpose matrix copy routine,
Cr   copying a subblock of matrix a to a subblock of matrix b.
Cr   Normally matrix nc{a,b} is the row dimension of matrix {a,b}
Cr   and nr{a,b} is 1.  Reverse nr and nc for a transposed matrix.
Cr   Arrays are locally one-dimensional so as to optimize inner loop.
Cr
Cr   Example: Set 3-by-2 block of matrix c to constant z
Cr     call dmcpy(z,0,0,c,nc,1,3,2)
Cr   Here scalar z is represented by an array of 0 dimension
Cr   Example: copy nbas-by-3 'wk' to its transpose 'bas', and the reverse:
Cr     call dmcpy(wk,nbas,1,bas,1,3,nbas,3)
Cr     call dmcpy(bas,1,3,wk,nbas,1,nbas,3)
C ----------------------------------------------------------------
C     implicit none
      integer nca,nra,ncb,nrb,n,m
      double precision a(0:*), b(0:*)
      integer i,j,ia,ib
      if (nra .eq. 1 .and. nrb .eq. 1) then
        do  10  j = 0, m-1
        ia = j*nca
        ib = j*ncb
        do  10  i = 0, n-1
          b(i+ib) = a(i+ia)
   10   continue
        return
      endif
      do  20  i = n-1, 0, -1
        ia = i*nra+m*nca
        ib = i*nrb+m*ncb
        do  20  j = m-1, 0, -1 
        ia = ia-nca
        ib = ib-ncb
        b(ib) = a(ia)
   20 continue
      end
C      subroutine dvcpy(a,nca,b,ncb,n)
CC- General vector copy (obsolete: should use blas dcopy)
CC ----------------------------------------------------------------
CCi Inputs:
CCi   a,nca is the source vector and the number of
CCi      elements separating each element in the vector
CCi   b,ncb,nrb is the destination vector and the number of
CCi      elements separating each element in the vector
CCi   n: the number of elements to calculate
CCo Outputs:
CCo   result matrix stored in b
CCr Remarks:
CCr   This is a general-purpose vectore copy routine
CCr   Example: Set all elements of 3-by-2 matrix c to -1.d0
CCr     call dvcpy(-1.d0,0,c,1,3*2)
CCr   Example: Set block (n,m) of array a(p,m)=0, and a(i,i)=1, i=1,m
CCr     call dmcpy(0.d0,0,0,a,p,1,n,m)
CCr     call dvcpy(1.d0,0,a,p+1,m)
CC ----------------------------------------------------------------
C      integer nca,ncb,n
C      double precision a(0:*), b(0:*)
C      integer i,ia,ib
C
C      ia = n*nca
C      ib = n*ncb
C      do  200  i = n-1, 0, -1
C        ia = ia-nca
C        ib = ib-ncb
C        b(ib) = a(ia)
C  200 continue
C
C      return
C      end

