subroutine dsred(nm,n,hr,ar)
  !- Reduction of nonorthogonal symmetric matrix to orthogonal form
  ! ----------------------------------------------------------------
  !i Inputs
  !i   h,nm: hermitian matrix, declared as h(nm,*).  (Lower triangle only)
  !i   a: nonorthogonality matrix, Cholesky-decomposed by dschd into L(L+)
  !i   n:  order of h and a
  !o Outputs
  !o   H replaced by H'' = L^-1 H (L+)^-1
  !r Remarks
  !r   Makes h'ij  = (hij  - sum_k<i lik h'kj)/lii
  !r         h''ij = (h'ij - sum_k<j h''ik (l*)jk)/ljj
  !r   This version uses vectorizable BLAS-style daxpy loops.
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  integer :: n,nm
  double precision :: hr(nm,n),ar(nm,n)
  ! Local parameters
  integer :: i,j,k

  ! --- Make h' ---
  do  10  i = 1, n
     do  k = 1, i-1
        call daxpy(n,-ar(i,k),hr(k,1),nm,hr(i,1),nm)
     enddo
     call dscal(n,1/ar(i,i),hr(i,1),nm)
10 enddo

  ! --- Make h'' (lower triangle only) ---
  do  30  j = 1, n
     do    k = 1, j-1
        call daxpy(n-j+1,-ar(j,k),hr(j,k),1,hr(j,j),1)
     enddo
     call dscal(n-j+1,1/ar(j,j),hr(j,j),1)

     ! --- Copy lower triangle into upper ---
     do  i = j+1, n
        hr(j,i) =  hr(i,j)
     enddo
30 enddo

  !      print 337, hr,hi
  !      pause
  !  337 format(9f10.6)
end subroutine dsred

