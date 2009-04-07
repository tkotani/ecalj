      subroutine zmcpy(transa,a,nca,nra,b,ncb,nrb,n,m)
C- General matrix copy, complex matrices
C ----------------------------------------------------------------
Ci Inputs:
Ci   transa:if 'c' or 'C', copy complex conjugate of a into b
Ci   a     :source matrix
Ci   nca   :number of elements separating columns of a
Ci   nra   :number of elements separating rows of a
Ci   ncb   :number of elements separating columns of b
Ci   nrb   :number of elements separating rows of b
Ci   n     :number of rows to copy
Ci   m     :number of columns to copy
Co Outputs:
Co   b     :elements (1..n , 1..m) copied into b from a
Cr Remarks:
Cr   This is a general-purpose matrix copy routine,
Cr   copying a subblock of matrix a to a subblock of matrix b.
Cr
Cr   Normally matrix nc{a,b} is the row dimension of matrix {a,b}
Cr   and nr{a,b} is 1.  Reverse nr and nc for a transposed matrix.
Cr   Arrays are locally one-dimensional so as to optimize inner loop.
Cr
Cr   Example: Set 3-by-2 block of matrix c to constant z
Cr     call zmcpy('n',z,0,0,c,nc,1,3,2)
Cr   Here scalar z is passed an array of 0 dimension
Cr
Cr   Example: copy nbas-by-3 'wk' to its transpose 'bas', and the reverse:
Cr     call zmcpy('n',wk,nbas,1,bas,1,3,nbas,3)
Cr     call zmcpy('n',bas,1,3,wk,nbas,1,nbas,3)
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character*1 transa
      integer nca,nra,ncb,nrb,n,m
      double complex a(0:*), b(0:*)
C ... Local parameters
      logical lcc
      integer i,j,ia,ib

      lcc = transa .eq. 'c' .or. transa .eq. 'C'

C ... Fast copy (unit strides)
      if (nra .eq. 1 .and. nrb .eq. 1) then
        do  j = 0, m-1
          ia = j*nca
          ib = j*ncb
          if (lcc) then
            do  i = 0, n-1
              b(i+ib) = dconjg(a(i+ia))
            enddo
          else
            do  i = 0, n-1
              b(i+ib) = a(i+ia)
            enddo
          endif
        enddo
        return
      endif

C ... General copy
      do  i = n-1, 0, -1
        ia = i*nra+m*nca
        ib = i*nrb+m*ncb
        if (lcc) then
          do  j = m-1, 0, -1
            ia = ia-nca
            ib = ib-ncb
            b(ib) = dconjg(a(ia))
          enddo
        else
          do  j = m-1, 0, -1
            ia = ia-nca
            ib = ib-ncb
            b(ib) = a(ia)
          enddo
        endif
      enddo
      end

