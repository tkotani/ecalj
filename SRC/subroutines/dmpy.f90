subroutine dmpy(a,nca,nra,b,ncb,nrb,c,ncc,nrc,n,m,l)
  !- matrix multiplication
  ! ----------------------------------------------------------------
  !i Inputs:
  !i   a,nca,nra is the left matrix and respectively the spacing
  !i      between elements in adjacent columns and rows.
  !i   b,ncb,nrb is the right matrix and respectively the spacing
  !i      between elements in adjacent columns and rows.
  !i   c,ncc,nrc is the product matrix and respectively the spacing
  !i      between elements in adjacent columns and rows.
  !i   n,m: the number of rows and columns, respectively, to calculate
  !i   l:   length of vector for matrix multiply
  !o Outputs:
  !o   product matrix stored in c
  !r Remarks:
  !r   This is a general-purpose matrix multiplication routine,
  !r   multiplying a subblock of matrix a by a subblock of matrix b.
  !r   Normally matrix nc{a,b,c} is the row dimension of matrix {a,b,c}
  !r   and nr{a,b,c} is 1.  Reverse nr and nc for a transposed matrix.
  !r   Arrays are locally one-dimensional so as to optimize inner loop,
  !r   which is executed n*m*l times.  No attempt is made to optimize
  !r   the outer loops, executed n*m times.
  !r     Examples: product of (n,l) subblock of a into (l,m) subblock of b
  !r   call dmpy(a,nrowa,1,b,nrowb,1,c,nrowc,1,n,m,l)
  !r     nrowa, nrowb, and nrowc are the leading dimensions of a, b and c.
  !r     To generate the tranpose of that product, use:
  !r   call dmpy(a,nrowa,1,b,nrowb,1,c,1,nrowc,n,m,l)
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  integer :: nca,nra,ncb,nrb,ncc,nrc,n,m,l
  double precision :: a(0:*), b(0:*), c(0:*)
  ! Local parameters
  double precision :: ar
  integer :: i,j,k,nccj,nrci,nrcicj,ncbj
  integer :: lda,ldb
  character(1) :: transa,transb

  if (nra == 1) then
     lda = nca
     transa = 'n'
  elseif (nca == 1) then
     lda = nra
     transa = 't'
  else
     lda = -1
  endif
  if (nrb == 1) then
     ldb = ncb
     transb = 'n'
  elseif (ncb == 1) then
     ldb = nrb
     transb = 't'
  else
     ldb = -1
  endif
  if (min(lda,ldb) < 0 .OR. nrc /= 1) goto 11
  ! if PARALLEL
  !      call pp_$dgemm(transa,transb,n,m,l,1d0,a,lda,b,ldb,0d0,c,ncc)
  ! else
  call dgemm(transa,transb,n,m,l,1d0,a,lda,b,ldb,0d0,c,ncc)
  ! endif
  return
11 continue
  ! endif
  ! --- Initialize array to zero ---
  do  10  i = n-1, 0, -1
     nrci = nrc*i
     nccj = -ncc
     do   j = m-1, 0, -1
        nccj = nccj + ncc
        nrcicj = nrci + nccj
        c(nrcicj) = 0
     enddo
10 enddo

  ! --- Do multiplication ---
  do  201  k = l-1, 0, -1
     do  20  i = n-1, 0, -1
        ar = a(nra*i + nca*k)
        if (ar == 0) goto 20
        call daxpy(m,ar,b(nrb*k),ncb,c(nrc*i),ncc)
20   enddo
201 enddo
end subroutine dmpy
