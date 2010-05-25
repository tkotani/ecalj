      subroutine hchd(nm,n,br,bi,dl,ierr)
C- Cholesky decomposition of hermitian matrix
C ----------------------------------------------------------------
Ci Inputs
Ci   nm: true row dimension of arrays b and z, as declared by caller
Ci   n:  order of the matrix system.  If b has already been Cholesky
Ci       decomposed, n should be prefixed with a minus sign.
Ci   b:  nonorthogonality matrix (only full upper triangle is used)
Ci       If n is negative, both b and dl must be altered
Co Outputs
Co   b:  strict lower triangle of Cholesky factor l of b
Co       The strict upper triangle unaltered.
Co   dl: diagonal elements of l.
Co   ierr:nonzero if decomposition finds matrix not pos. def.
Cr Remarks
Cr   Adapted from eispack reduc for hermitian matrices
C ----------------------------------------------------------------
C Passed parameters 
      integer ierr,n,nm
      double precision br(nm,n),bi(nm,n),dl(n)
C Local parameters 
      integer i,im1,j,k
      double precision xr,xi,y

c     .......... form l in the arrays b and dl ..........
      do 80 i = 1, n
        im1 = i - 1
        do 80 j = i, n
          xr =  br(i,j)
          xi = -bi(i,j)
          do 20 k = 1, im1
c             x = x - b*(i,k) * b(j,k)
            xr = xr - (br(i,k)*br(j,k) + bi(i,k)*bi(j,k))
            xi = xi - (br(i,k)*bi(j,k) - bi(i,k)*br(j,k))
   20     continue
          if (j .ne. i) go to 60
          ierr = i
          if (xr .le. 0.0d0) return
          y = dsqrt(xr)
          dl(i) = y
          go to 80
   60     br(j,i) = xr / y
          bi(j,i) = xi / y
   80 continue
      ierr = 0
      end

