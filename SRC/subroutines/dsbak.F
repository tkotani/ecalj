      subroutine dsbak(nm,n,ar,m,zr)
C- Back-transforms eigenvectors to nonorthogonal representation
C ----------------------------------------------------------------
Ci Inputs
Ci   z,nm: eigenvectors, declared as z(nm,*)
Ci   n: order of a and z
Ci   a: nonorthogonality matrix, Cholesky-decomposed by dschd into L(L+)
Ci   m: number of eigenvectors to be back transformed.
Co Outputs
Co   z transformed eigenvectors
Cr Remarks
Cr   Nonorthogonal eigenvectors are given by z <- (L+)^-1 z
Cr   This version uses vectorizable BLAS-style daxpy loops.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters 
      integer m,n,nm
      double precision ar(nm,n),zr(nm,m)
C Local parameters 
      integer nmi,k

      do  10  nmi = n, 1, -1
        do    k = n, nmi+1, -1
           call daxpy(m,-ar(k,nmi),zr(k,1),nm,zr(nmi,1),nm)
        enddo   
        call dscal(m,1/ar(nmi,nmi),zr(nmi,1),nm)
   10 continue

      end

