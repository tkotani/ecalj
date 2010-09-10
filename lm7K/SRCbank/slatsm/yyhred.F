      subroutine yyhred(nm,n,hr,hi,ar,ai,sw)
C- Reduction of nonorthogonal hermitian matrix to orthogonal form
C ----------------------------------------------------------------
Ci Inputs
Ci   h,nm: hermitian matrix, declared as h(nm,*).  (Lower triangle only)
Ci   a: nonorthogonality matrix Cholesky-decomposed by yyhchd into L(L+)
Ci   n:  order of h and a
Ci   sw: false if h and a matrix are real.
Co Outputs
Co   H replaced by H'' = L^-1 H (L+)^-1
Cr Remarks
Cr   Makes h'ij  = (hij  - sum_k<i lik h'kj)/lii
Cr         h''ij = (h'ij - sum_k<j h''ik (l*)jk)/ljj
Cr   This version uses vectorizable BLAS-style daxpy loops.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters 
      logical sw
      integer n,nm
      double precision hr(nm,n),hi(nm,n),ar(nm,n),ai(nm,n)
C Local parameters 
      integer i,j,k

      call tcn('yyhred')

C --- Make h' ---
      do  10  i = 1, n
C       do  20  k = 1, i-1
C         call yyaxpy(n,-ar(i,k),-ai(i,k),hr(k,1),hi(k,1),nm,
C    .      hr(i,1),hi(i,1),nm,sw)
        do  20  j = 1, n
          do  20  k = 1, i-1
            hr(i,j) = hr(i,j) - ar(i,k)*hr(k,j) + ai(i,k)*hi(k,j) 
            hi(i,j) = hi(i,j) - ai(i,k)*hr(k,j) - ar(i,k)*hi(k,j) 
   20   continue
        call dscal(n,1/ar(i,i),hr(i,1),nm)
        if (sw) call dscal(n,1/ar(i,i),hi(i,1),nm)
   10 continue

C --- Make h'' (lower triangle only) ---
      do  30  j = 1, n
C        do  40  k = 1, j-1
C          call yyaxpy(n-j+1,-ar(j,k),ai(j,k),hr(j,k),hi(j,k),1,
C     .      hr(j,j),hi(j,j),1,sw)
        do  40  i = j, n
          do  40  k = 1, j-1
            hr(i,j) = hr(i,j) - ar(j,k)*hr(i,k) - ai(j,k)*hi(i,k)
            hi(i,j) = hi(i,j) - ar(j,k)*hi(i,k) + ai(j,k)*hr(i,k)
   40   continue
        call dscal(n-j+1,1/ar(j,j),hr(j,j),1)
        if (sw) call dscal(n-j+1,1/ar(j,j),hi(j,j),1)

C --- Copy lower triangle into upper ---
        do  50  i = j+1, n
          hr(j,i) =  hr(i,j)
          if (sw) hi(j,i) = -hi(i,j)
   50   continue
   30 continue

c      print 337, hr,hi
c      pause
c  337 format(9f10.6)

      call tcx('yyhred')
      end

