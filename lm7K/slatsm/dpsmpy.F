      subroutine dpsmpy(a,iofa,jofa,b,ldb,c,ldc,nr,nc,l)
C- subblock matrix multiplication, (packed) (normal) -> (normal)
C ----------------------------------------------------------------
Ci Inputs:
Ci   a is the left matrix (packed)
Ci   b,ldb is the right matrix and the row dimension of b
Ci   c,ldc is the product matrix and the row dimension of c
Ci   nr,nc: the number of rows and columns, respectively, to calculate
Ci   l:   length of vector for matrix inner product
Co Outputs:
Co   product matrix stored in c
Cr Remarks:
Cr   Examples: product of (nr,l) subblock of a into (l,nc) subblock of b
Cr   call dmpy(a,nrowa,1,b,nrowb,1,c,nrowc,1,nr,nc,l)
Cr     nrowa, nrowb, and nrowc are the leading dimensions of a, b and c.
Cr     To generate the tranpose of that product, use:
Cr   call dmpy(a,nrowa,1,b,nrowb,1,c,1,nrowc,nr,nc,l)
C ----------------------------------------------------------------
C Passed Parameters
      integer ldb,ldc,nr,nc,l,iofa,jofa
      double precision a(1), b(ldb,1), c(ldc,1)
C Local parameters
      double precision sum
      integer i,j,k,offa

      do  20  i = 1, nr
        ia = i+iofa-1
        do  20  j = 1, nc
          sum = 0
          offa = (ia*(ia-1))/2 + jofa-1
          ktop = max(ia-jofa,0)
          do  21  k = 1, ktop
C         print *, i,j,k, k+offa, a(k+offa),b(k,j), '*'
            sum = sum + a(k+offa)*b(k,j)
   21     continue
          offa = ((ia-1)*ia)/2 + ia
          do  22  k = ktop+1, l
C         print *, i,j,k, offa, a(offa),b(k,j)
            sum = sum + a(offa)*b(k,j)
            offa = offa + k+jofa-1
   22     continue
          c(i,j) = sum
   20 continue
      end

