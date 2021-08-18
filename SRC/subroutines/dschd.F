      subroutine dschd(nm,n,ar,wk,swI,ierr)
C- Cholesky decomposition of hermitian matrix
C ----------------------------------------------------------------
Ci Inputs
Ci   a,nm: hermitian matrix, declared as a(nm,*).  (Lower triangle only)
Ci   n:  order of a.
Ci   swI:true if to return L^-1
Ci   wk: real work array of dimension n (uneeded if swI is false)
Co Outputs
Co   A replaced by L or L^-1 if swI true
Co   ierr:nonzero if matrix not positive definite.
Cr Remarks
Cr   Makes ljj = (ajj - sum_k<j ljk (l+)jk)^1/2
Cr         lij = (aij - sum_k<j lik (l+)jk)/ljj for i>j
Cr   The strict upper triangle is unused.
Cr   This version uses BLAS-style daxpy loops with unit stride.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters 
      logical swI
      integer ierr,n,nm
      double precision ar(nm,n),wk(n)
C Local parameters 
      integer i,j,k
      double precision xr

      call tcn('dschd')

C --- Cholesky decomposition of a into L(L+) (lower triangle only) ---
      do  10  j = 1, n
        do  20  k = 1, j-1
          xr = ar(j,k)
          do    i = j, n
             ar(i,j) = ar(i,j) - xr*ar(i,k)
          enddo   
   20   continue
        ierr = j
        if (ar(j,j) .le. 0) return
        ar(j,j) = dsqrt(ar(j,j))
        call dscal(n-j,1/ar(j,j),ar(j+1,j),1)
   10 continue
      ierr = 0
      if (.not. swI) return

C --- Inversion of L (lower triangle only) ---
      do    j = 1, n
         ar(j,j) = 1/ar(j,j)
      enddo   

      do  40  j = 2, n
        call dcopy(n-j+1,ar(j,j-1),1,wk,1)
        call dpzero(ar(j,j-1),n-j+1)
C This loop parallelizable ...
        do  50  k = 1, j-1
          xr = ar(j-1,k)
          do  i = j, n
             ar(i,k) = ar(i,k) - xr*wk(i-j+1)
          enddo   
          ar(j,k) = ar(j,k)*ar(j,j)
   50   continue
   40 continue

      call tcx('dschd')
      end

