      subroutine dsred(nm,n,hr,ar)
C- Reduction of nonorthogonal symmetric matrix to orthogonal form
C ----------------------------------------------------------------
Ci Inputs
Ci   h,nm: hermitian matrix, declared as h(nm,*).  (Lower triangle only)
Ci   a: nonorthogonality matrix, Cholesky-decomposed by dschd into L(L+)
Ci   n:  order of h and a
Co Outputs
Co   H replaced by H'' = L^-1 H (L+)^-1
Cr Remarks
Cr   Makes h'ij  = (hij  - sum_k<i lik h'kj)/lii
Cr         h''ij = (h'ij - sum_k<j h''ik (l*)jk)/ljj
Cr   This version uses vectorizable BLAS-style daxpy loops.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters 
      integer n,nm
      double precision hr(nm,n),ar(nm,n)
C Local parameters 
      integer i,j,k

C --- Make h' ---
      do  10  i = 1, n
        do  20  k = 1, i-1
   20   call daxpy(n,-ar(i,k),hr(k,1),nm,hr(i,1),nm)
        call dscal(n,1/ar(i,i),hr(i,1),nm)
   10 continue

C --- Make h'' (lower triangle only) ---
      do  30  j = 1, n
        do  40  k = 1, j-1
   40   call daxpy(n-j+1,-ar(j,k),hr(j,k),1,hr(j,j),1)
        call dscal(n-j+1,1/ar(j,j),hr(j,j),1)

C --- Copy lower triangle into upper ---
        do  50  i = j+1, n
   50   hr(j,i) =  hr(i,j)
   30 continue

c      print 337, hr,hi
c      pause
c  337 format(9f10.6)
      end

