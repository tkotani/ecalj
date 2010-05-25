      subroutine hbak(nm,n,br,bi,dl,m,zr,zi)
C- Back-transforms eigenvectors to nonorthogonal representation
C ----------------------------------------------------------------
Ci Inputs
Ci   nm: true row dimension of arrays b and z
Ci   n:  order of the matrix system.
Ci   b:  Cholesky-decomposed nonorthogonality matrix, as decomposed
Ci       by reduch (in its strict lower triangle).
Ci   dl: contains further information about the transformation.
Ci   m:  number of eigenvectors to be back transformed.
Ci   z:  eigenvectors to be back transformed (first m columns)
Co Outputs
Co   z transformed eigenvectors
Cr Remarks
Cr   Adapted from rebakh, eispack
C ----------------------------------------------------------------
C Passed parameters 
      integer m,n,nm
      double precision br(nm,n),bi(nm,n),dl(n),zr(nm,m),zi(nm,m)
C Local parameters 
      integer i,j,k
      double precision xr,xi

      do 100 j = 1, m
        do 100 i = n, 1, -1
          xr = zr(i,j)
          xi = zi(i,j)
          do 60 k = i+1, n
            xr = xr - (br(k,i)*zr(k,j) + bi(k,i)*zi(k,j))
            xi = xi - (br(k,i)*zi(k,j) - bi(k,i)*zr(k,j))
   60     continue
          zr(i,j) = xr/dl(i)
          zi(i,j) = xi/dl(i)
  100 continue
      end

