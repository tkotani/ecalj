      subroutine svbksb(nm,m,n,nrhs,w,u,v,b,x,wk)
C- Solves linear system A x = b, where A has been decomposed by svd
C ----------------------------------------------------------------
Ci Inputs
Ci   nm:  row dimension of arrays u and v, b and x
Ci    m:  number of rows of u.
Ci    n:  number of columns of u and the order of v.
Ci nhrs:  number of right-hand side vectors b
Ci u,w,v: the singular decomposition of A: A = u w v+
Ci        with w the eigenvalues of A
Ci   wk:  work array of dimension n
Ci    b:  right hand side of system of equations A x = b
Ci        b is left intact, unless b and x occupy the same address space
Co Outputs
Co    x:  solution to linear system A x = b
Co        x may occupy the same address space as b
Cr Remarks
Cr   Call this routine after eg calling svd:
Cr     call svd(nm,m,n,u,w,.true.,u,.true.,v,ierr,wk)
Cr   or for symmetric matrices, the eigenvalue routine, eg:
Cr     call rs(nm,n,a,w,1,v,wk,wk(1,2),ierr)
Cr     call svbksb(nm,n,n,m,w,v,v,b,x,wk)
Cr
Cr   Any zero eigenvalues w are here projected out of the Hilbert space
Cr   when solving for x. After calling svd, the caller may independently
Cr   zero out any eigenvalues in w (likely candidates are small elements
Cr   in w); this has the effect of excluding the portion of Hilbert
Cr   corresponding to that eigenvlaue in the solution for x
Cr   Hilbert space in the solution to x.
Cb Bugs
Cb   The loops in this routine have not been optimized.
C ----------------------------------------------------------------
C     implicit none
      integer m,nm,n,nrhs
      double precision u(nm,n),w(n),v(nm,n),b(nm,nrhs),x(nm,nrhs),wk(n)
      integer i,j,jj,irhs
      double precision s

      do  10  irhs = 1, nrhs

C --- wk = w^-1 u+ b ---
      do  12  j = 1, n
        s = 0d0
        if (w(j) .ne. 0d0) then
          do  11  i = 1, m
   11     s = s + u(i,j)*b(i,irhs)
          s = s/w(j)
        endif
        wk(j) = s
   12 continue

C --- x = v (w^-1 u+ b) ---
      do  14  j = 1, n
        s = 0d0
        do  13  jj = 1, n
   13   s = s + v(j,jj)*wk(jj)
        x(j,irhs) = s
   14 continue

   10 continue

      end

