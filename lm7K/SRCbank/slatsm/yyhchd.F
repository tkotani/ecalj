C#define DDOT 1
      subroutine yyhchd(nm,n,ar,ai,wk,swI,sw,ierr)
C- Cholesky decomposition of hermitian matrix
C ----------------------------------------------------------------
Ci Inputs
Ci   a,nm: hermitian matrix, declared as a(nm,*).  (Lower triangle only)
Ci   n:  order of a.
Ci   sw: false if matrix is real (works only for BLAS version)
Ci   swI:true if to return L^-1
Ci   wk: real work array of precision 2*n (uneeded if swI is false)
Co Outputs
Co   A replaced by L or L^-1 if swI true
Co   ierr:nonzero if matrix not positive definite.
Cr Remarks
Cr   Makes ljj = (ajj - sum_k<j ljk (l+)jk)^1/2
Cr         lij = (aij - sum_k<j lik (l+)jk)/ljj for i>j
Cr   default version: uses BLAS-style daxpy loops with unit stride.
Cr                    The strict upper triangle is unused.
Cr   DDOT version: uses ddot-style loops with unit stride.
Cr                 l+ is accumulated in the uppper triangle.
Cr                 This version parallelizes well.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters 
      logical swI,sw
      integer ierr,n,nm
      double precision ar(nm,n),ai(nm,n),wk(n,2)
C Local parameters 
      integer i,j,k
      double precision xr,xi

      call tcn('yyhchd')

C$$$#ifndef DDOT
C$$$C --- Cholesky decomposition of a into L(L+) (lower triangle only) ---
C$$$      do  10  j = 1, n
C$$$        do  20  k = 1, j-1
C$$$#if BLAS
C$$$          call yyaxpy(n-j+1,-ar(j,k),ai(j,k),ar(j,k),ai(j,k),1,
C$$$     .    ar(j,j),ai(j,j),1,sw)
C$$$#else
C$$$          xr = ar(j,k)
C$$$          xi = ai(j,k)
C$$$          do  23  i = j, n
C$$$            ar(i,j) = ar(i,j) - xr*ar(i,k) - xi*ai(i,k)
C$$$            ai(i,j) = ai(i,j) - xr*ai(i,k) + xi*ar(i,k)
C$$$   23     continue
C$$$#endif
C$$$   20   continue
C$$$        ierr = j
C$$$        if (ar(j,j) .le. 0) return
C$$$        ar(j,j) = dsqrt(ar(j,j))
C$$$        call dscal(n-j,1/ar(j,j),ar(j+1,j),1)
C$$$        if (sw) call dscal(n-j,1/ar(j,j),ai(j+1,j),1)
C$$$   10 continue
C$$$#else
C --- Cholesky decomposition of a into L(L+) (upper triangle only) ---
      do  10  j = 1, n
        do  23  i = j, n
          xr = 0
          xi = 0
          do  20  k = 1, j-1
            xr = xr - ar(k,j)*ar(k,i) - ai(k,j)*ai(k,i)
            xi = xi - ar(k,j)*ai(k,i) + ai(k,j)*ar(k,i)
   20     continue
          ar(j,i) = ar(j,i) + xr
          ai(j,i) = ai(j,i) + xi
   23   continue
        ierr = j
        if (ar(j,j) .le. 0) return
        ar(j,j) = dsqrt(ar(j,j))
        call dscal(n-j,1/ar(j,j),ar(j,j+1),nm)
        if (sw) call dscal(n-j,1/ar(j,j),ai(j,j+1),nm)
   10 continue

C --- Swap to lower triangle ---
      do  25  j = 1, n
        do  25  i = j+1, n
          ar(i,j) =  ar(j,i)
          ai(i,j) = -ai(j,i)
   25 continue
C$$$#endif
      ierr = 0
      if (.not. swI) return

C --- Inversion of L (lower triangle only) ---
      do  30  j = 1, n
   30 ar(j,j) = 1/ar(j,j)

      do  40  j = 2, n
        call dcopy(n-j+1,ar(j,j-1),1,wk,1)
        call dcopy(n-j+1,ai(j,j-1),1,wk(1,2),1)
        call dpzero(ar(j,j-1),n-j+1)
        call dpzero(ai(j,j-1),n-j+1)
        do  50  k = 1, j-1
C$$$#if BLAS
          call yyaxpy(n-j+1,-ar(j-1,k),-ai(j-1,k),wk,wk(1,2),1,
     .    ar(j,k),ai(j,k),1,sw)
C$$$#else
C$$$          xr = ar(j-1,k)
C$$$          xi = ai(j-1,k)
C$$$          do  55  i = j, n
C$$$            ar(i,k) = ar(i,k) - xr*wk(i-j+1,1) + xi*wk(i-j+1,2)
C$$$            ai(i,k) = ai(i,k) - xr*wk(i-j+1,2) - xi*wk(i-j+1,1)
C$$$   55     continue
C$$$#endif
          ar(j,k) = ar(j,k)*ar(j,j)
          ai(j,k) = ai(j,k)*ar(j,j)
   50   continue
   40 continue

      call tcx('yyhchd')
      end

