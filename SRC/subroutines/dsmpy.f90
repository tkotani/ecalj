      subroutine dsmpy(n,k,a,lda,b,ldb,beta,c,ldc)
C- Matrix multiplication: result assumed symmetric
C ----------------------------------------------------------------------
Ci Inputs
Ci   n     :number of rows of a and c and columns of b and c
Ci   k     :length of inner product
Ci   a,lda :first matrix in a*b product, and leading dimension
Ci   b,ldb :second matrix in a*b product, and leading dimension
Ci   beta  :add beta * c to result
Ci   c,ldc :resultant matrix, and leading dimension:
Ci          beta * (input c) is added to resultant c
Co Outputs
Co   c     :result matrix, a*b + beta*(input c)
Cr Remarks
Cr   It is assumed that c(i,j) = c(j,i).
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lda, ldb, ldc, n, k
      double precision a(lda,*), b(ldb,*), c(ldc,*), beta
C ... Local parameters
      integer i,j,kk,ncol,ic,i1,i2,mcol,ibar
      parameter (ncol=48)

C --- For each ncol columns, do ---
      do  10  ic = -ncol/2, n, ncol
        mcol = min(n-ic+1,ncol)
        ibar = ic+(mcol-1)/2

C  ...  Matrix subblock (1..ibar , ic..mcol)
        if (ic .gt. 0) then
          call dgemm('N','N',ibar,mcol,k,1d0,a,lda,
     .    b(1,ic),ldb,beta,c(1,ic),ldc)
        endif

C  ...  Matrix triangle (i1,i2)
        i1 = max(ic,ibar,0)+1
        i2 = ic+mcol-1

        if (beta .eq. 0d0) then
          do    j = i1, i2
            do  i = j, i2
               c(i,j) = 0
            enddo
            enddo
        elseif (beta .ne. 1d0) then
          do    j = i1, i2
            do  i = j, i2
               c(i,j) = beta*c(i,j)
            enddo
            enddo
        endif

        do  12  j = i1, i2
          do  14  kk = 1, k
            do  16  i = j, i2
              c(i,j) = c(i,j) + a(i,kk) * b(kk,j)
   16       continue
   14     continue
   12   continue

        do    j = i1, i2
          do  i = j, i2
             c(j,i) = c(i,j)
          enddo
          enddo

   10 continue

C...  Fill in lower half of triangle
      do    j = 1, n
        do  i = 1, j-1
           c(j,i) = c(i,j)
        enddo
      enddo
      end

