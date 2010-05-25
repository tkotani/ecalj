      subroutine hcinv(nm,n,br,bi,dl,lopt)
C- Obtain inverse from Cholesky decomposed matrix
C ----------------------------------------------------------------
Ci Inputs
Ci   nm: true row dimension of arrays b and z, as declared by caller
Ci   n:  order of the matrix system.
Ci   b:  Cholesky-decomposed nonorthogonality matrix
Ci       (strict lower triangle; see hchd)
Ci   dl: diagonal elements of l (see hchd).
Ci   lopt: false, only upper half of b^-1 is generated, and lower
Ci         half continues to hold the c.d. of b
Ci         true,  lower half of b also filled in
Co Outputs
Co   b:  Inverse of Cholesky matrix
Cr Remarks
Cr   Adapted from eispack reduc for hermitian matrices
C ----------------------------------------------------------------
C Passed parameters 
      integer n,nm
      logical lopt
      double precision br(nm,n),bi(nm,n),dl(n)
C Local parameters 
      integer i,im1,j,k
      double precision xr,xi,y
#if BLAS
      double precision ddot
#endif

C --- form inv(l) and store transpose in full upper triangle of b ---
C  (l^-1)_ij  =  (del_ij - sum_k<i  l_ik * (l^-1)_kj) / l_ii
      do 200 j = 1, n
        do 200 i = j, n
          y = dl(i)
#if BLAS
          if (j .eq. i) then
            xr = 1d0
            xi = 0
          else
            xr = - ddot(i-j,br(i,j),nm,br(j,j),nm)
     .      + ddot(i-j,bi(i,j),nm,bi(j,j),nm)
            xi = - ddot(i-j,br(i,j),nm,bi(j,j),nm)
     .      - ddot(i-j,bi(i,j),nm,br(j,j),nm)
          endif
#else
          xr = 0
          xi = 0
          if (j .eq. i) xr = 1d0
          im1 = i-1
          do 160 k = j, im1
            xr = xr - (br(i,k)*br(j,k) - bi(i,k)*bi(j,k))
            xi = xi - (br(i,k)*bi(j,k) + bi(i,k)*br(j,k))
  160     continue
#endif
          br(j,i) = xr / y
          bi(j,i) = xi / y
  200 continue
c
c     .......... pre-multiply by inv(l) and overwrite ..........
      do 300 i = 1, n
        do 300 j = i, n
C           x = x - l^dag_ik * l_kj; here i <= j <= k
#if BLAS
          br(i,j) = ddot(n+1-j,br(i,j),nm,br(j,j),nm) +
     .    ddot(n+1-j,bi(i,j),nm,bi(j,j),nm)
          bi(i,j) = ddot(n+1-j,br(i,j),nm,bi(j,j),nm) -
     .    ddot(n+1-j,bi(i,j),nm,br(j,j),nm)
#else
          xr = 0
          xi = 0
          do 260 k = j, n
            xr = xr + (br(i,k)*br(j,k) + bi(i,k)*bi(j,k))
            xi = xi + (br(i,k)*bi(j,k) - bi(i,k)*br(j,k))
  260     continue
          br(i,j) = xr
          bi(i,j) = xi
#endif
  300 continue

      if (.not. lopt) return
      do  400  i = 1, n
        do  400  j = i+1, n
          br(j,i) =  br(i,j)
          bi(j,i) = -bi(i,j)
  400 continue
      end

