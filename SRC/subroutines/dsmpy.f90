subroutine dsmpy(n,k,a,lda,b,ldb,beta,c,ldc)
  !- Matrix multiplication: result assumed symmetric
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   n     :number of rows of a and c and columns of b and c
  !i   k     :length of inner product
  !i   a,lda :first matrix in a*b product, and leading dimension
  !i   b,ldb :second matrix in a*b product, and leading dimension
  !i   beta  :add beta * c to result
  !i   c,ldc :resultant matrix, and leading dimension:
  !i          beta * (input c) is added to resultant c
  !o Outputs
  !o   c     :result matrix, a*b + beta*(input c)
  !r Remarks
  !r   It is assumed that c(i,j) = c(j,i).
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: lda, ldb, ldc, n, k
  double precision :: a(lda,*), b(ldb,*), c(ldc,*), beta
  ! ... Local parameters
  integer :: i,j,kk,ncol,ic,i1,i2,mcol,ibar
  parameter (ncol=48)

  ! --- For each ncol columns, do ---
  do  10  ic = -ncol/2, n, ncol
     mcol = min(n-ic+1,ncol)
     ibar = ic+(mcol-1)/2

     !  ...  Matrix subblock (1..ibar , ic..mcol)
     if (ic > 0) then
        call dgemm('N','N',ibar,mcol,k,1d0,a,lda, &
             b(1,ic),ldb,beta,c(1,ic),ldc)
     endif

     !  ...  Matrix triangle (i1,i2)
     i1 = max(ic,ibar,0)+1
     i2 = ic+mcol-1

     if (beta == 0d0) then
        do    j = i1, i2
           do  i = j, i2
              c(i,j) = 0
           enddo
        enddo
     elseif (beta /= 1d0) then
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
16         enddo
14      enddo
12   enddo

     do    j = i1, i2
        do  i = j, i2
           c(j,i) = c(i,j)
        enddo
     enddo

10 enddo

  !...  Fill in lower half of triangle
  do    j = 1, n
     do  i = 1, j-1
        c(j,i) = c(i,j)
     enddo
  enddo
end subroutine dsmpy

