      subroutine yyhbak(nm,n,ar,ai,m,zr,zi,sw)
C- Back-transforms eigenvectors to nonorthogonal representation
C ----------------------------------------------------------------
Ci Inputs
Ci   z,nm: eigenvectors, declared as z(nm,*)
Ci   n: order of a and z
Ci   a: nonorthogonality matrix, Cholesky-decomposed by yyhchd into L(L+)
Ci   m: number of eigenvectors to be back transformed.
Ci   sw:false if a is real 
Co Outputs
Co   z transformed eigenvectors
Cr Remarks
Cr   Nonorthogonal eigenvectors are given by z <- (L+)^-1 z
Cr   This version uses vectorizable BLAS-style daxpy loops.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters 
      logical sw
      integer m,n,nm
      double precision ar(nm,n),ai(nm,n),zr(nm,m),zi(nm,m)
C Local parameters 
      integer nmi,j,k

      call tcn('yyhbak')

      do  10  nmi = n, 1, -1
C       do  20  k = n, nmi+1, -1
C         call yyaxpy(m,-ar(k,nmi),ai(k,nmi),zr(k,1),zi(k,1),nm,
C    .      zr(nmi,1),zi(nmi,1),nm,sw)
        do  20  j = 1, m
          do  20  k = n, nmi+1, -1
            zr(nmi,j) = zr(nmi,j) - ar(k,nmi)*zr(k,j) - ai(k,nmi)*zi(k,j)
            zi(nmi,j) = zi(nmi,j) - ar(k,nmi)*zi(k,j) + ai(k,nmi)*zr(k,j)
   20   continue
        call dscal(m,1/ar(nmi,nmi),zr(nmi,1),nm)
        call dscal(m,1/ar(nmi,nmi),zi(nmi,1),nm)
   10 continue

C     call yywrm(0,' ',12,6,'(5f16.10)',zr,zi,n,n,n)

      call tcx('yyhbak')
      end

