      subroutine dpmpy(a,b,nscb,nsrb,c,nscc,nsrc,nr,nc,l)
C- matrix multiplication, (packed) (normal) -> (normal)
C ----------------------------------------------------------------
Ci Inputs:
Ci   a is the left matrix (packed)
Ci   b,nscb,nsrb is the right matrix and respectively the spacing
Ci      between column elements and row elements.
Ci   c,nscc,nsrc is the product matrix and respectively the number of
Ci      elements separating columns and rows.
Ci   nr,nc: the number of rows and columns, respectively, to calculate
Ci   l:   length of vector for matrix multiply
Co Outputs:
Co   product matrix stored in c
Cr Remarks:
Cr   This is a general-purpose matrix multiplication routine,
Cr   multiplying a subblock of matrix a by a subblock of matrix b.
Cr   Normally matrix nc{a,b,c} is the row dimension of matrix {a,b,c}
Cr   and nr{a,b,c} is 1.  Reverse nr and nc for a transposed matrix.
Cr   Arrays are locally one-dimensional so as to optimize inner loop,
Cr   which is executed nr*nc*l times.  No attempt is made to optimize
Cr   the outer loops, executed nr*nc times.
Cr     Examples: product of (nr,l) subblock of a into (l,nc) subblock of b
Cr   call dmpy(a,nrowa,1,b,nrowb,1,c,nrowc,1,nr,nc,l)
Cr     nrowa, nrowb, and nrowc are the leading dimensions of a, b and c.
Cr     To generate the tranpose of that product, use:
Cr   call dmpy(a,nrowa,1,b,nrowb,1,c,1,nrowc,nr,nc,l)
C ----------------------------------------------------------------
C Passed Parameters
      integer nscb,nsrb,nscc,nsrc,nr,nc,l
      double precision a(0:*), b(0:*), c(0:*)
C Local parameters
      double precision sum
      integer i,j,k,offa,offb

      do  20  i = 0, nr-1
        do  20  j = 0, nc-1
          sum = 0
          offa = (i*(i+1))/2
          offb = nscb*j
          do  21  k = 0, i-1
            sum = sum + a(offa)*b(offb)
            offa = offa + 1
            offb = offb + nsrb
   21     continue
          do  22  k = i, l-1
            sum = sum + a(offa)*b(offb)
            offa = offa + k+1
            offb = offb + nsrb
   22     continue
          c(i*nsrc+j*nscc) = sum
   20 continue
      end

