      subroutine dmadd(a,nca,nra,scalea,b,ncb,nrb,scaleb,c,ncc,nrc,n,m)
C- General matrix addition
C ----------------------------------------------------------------
Ci Inputs:
Ci   a,nca,nra is the left matrix and respectively the number of
Ci      elements separating columns and rows.
Ci   b,ncb,nrb is the right matrix and respectively the number of
Ci      elements separating columns and rows.
Ci   c,ncc,nrc is the result matrix and respectively the number of
Ci      elements separating columns and rows.
Ci   n,m: the number of rows and columns, respectively, to calculate
Co Outputs:
Co   result matrix stored in c
Cr Remarks:
Cr   This is a general-purpose matrix linear combination routine,
Cr   adding a subblock of matrix a to a subblock of matrix b.
Cr   Normally matrix nc{a,b,c} is the row dimension of matrix {a,b,c}
Cr   and nr{a,b,c} is 1.  Reverse nr and nc for a transposed matrix.
Cr   Arrays are locally one-dimensional so as to optimize inner loop.
Cr
Cr   Destination matrix c can coincide with either a or b, provided that
Cr   the transpose of the coincident matrix is not taken.
Cr   Example: Add 3-by-2 block of (transpose of a - .5*b) into c
Cr     call dmadd(a,1,na,1.d0,b,nb,0,-.5d0,c,nc,1,3,2)
Cr   This version checks for unit stride lengths in (a,b,c)
C ----------------------------------------------------------------
C     implicit none
      integer nca,nra,ncb,nrb,ncc,nrc,n,m
      double precision a(0:*), b(0:*), c(0:*), scalea, scaleb
      integer i,j,ia,ib,ic

C --- Case unit stride length for each column ---
      if (nra .eq. 1 .and. nrb .eq. 1 .and. nrc .eq. 1) then
        if (scalea .eq. 0) then
          do  10  j = 0, m-1
            ia = nca * j
            ib = ncb * j
            ic = ncc * j
            do  10  i = 0, n-1
              c(i+ic) = b(i+ib)*scaleb
   10     continue
          return
        elseif (scaleb .eq. 0) then
          do  20  j = 0, m-1
            ia = nca * j
            ib = ncb * j
            ic = ncc * j
            do  20  i = 0, n-1
              c(i+ic) = a(i+ia)*scalea
   20     continue
          return
        elseif (scalea .eq. 1 .and. scaleb .eq. 1) then
          do  30  j = 0, m-1
            ia = nca * j
            ib = ncb * j
            ic = ncc * j
            do  30  i = 0, n-1
              c(i+ic) = a(i+ia) + b(i+ib)
   30     continue
          return
        elseif (scalea .eq. 1 .and. scaleb .eq. -1) then
          do  40  j = 0, m-1
            ia = nca * j
            ib = ncb * j
            ic = ncc * j
            do  40  i = 0, n-1
              c(i+ic) = a(i+ia) - b(i+ib)
   40     continue
          return
        else
          do  50  j = 0, m-1
            ia = nca * j
            ib = ncb * j
            ic = ncc * j
            do  50  i = 0, n-1
              c(i+ic) = a(i+ia)*scalea + b(i+ib)*scaleb
   50     continue
          return
        endif
      endif

C --- General case ---
      if (scalea .ne. 0 .and. scaleb .ne. 0) then
        do  210  i = n-1, 0, -1
          ia = i*nra+m*nca
          ib = i*nrb+m*ncb
          ic = i*nrc+m*ncc
          do  210  j = m-1, 0, -1 
            ia = ia-nca
            ib = ib-ncb
            ic = ic-ncc
            c(ic) = a(ia)*scalea + b(ib)*scaleb
 210    continue
      elseif (scalea .ne. 0) then
        do  220  i = n-1, 0, -1
          ia = i*nra+m*nca
          ic = i*nrc+m*ncc
          do  220  j = m-1, 0, -1 
            ia = ia-nca
            ic = ic-ncc
            c(ic) = a(ia)*scalea
  220   continue
      elseif (scaleb .ne. 0) then
        do  230  i = n-1, 0, -1
          ib = i*nrb+m*ncb
          ic = i*nrc+m*ncc
          do  230  j = m-1, 0, -1 
            ib = ib-ncb
            ic = ic-ncc
            c(ic) = b(ib)*scaleb
  230   continue
      endif
      end

