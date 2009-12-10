      subroutine hred(nm,n,ar,ai,br,bi,dl)
C- Reduction of nonorthogonal hermitian matrix to orthogonal form
C ----------------------------------------------------------------
Ci Inputs
Ci   nm: true row dimension of arrays b and z, as declared by caller
Ci   n:  order of the matrix system.
Ci   a:  hamiltonian matrix (only full upper triangle is used)
Ci   b:  Cholesky-decomposed nonorthogonality matrix
Ci       (strict lower triangle; see hchd)
Ci   dl: diagonal elements of l (see hchd).
Co Outputs
Co   a:  reduced hermitian matrix:  a <-  ldag^-1 a l^-1
Co       The strict upper triangle unaltered.
Cr Remarks
Cr   Adapted from eispack reduc for hermitian matrices
C ----------------------------------------------------------------
C Passed parameters 
      integer n,nm
      double precision ar(nm,n),ai(nm,n),br(nm,n),bi(nm,n),dl(n)
C Local parameters 
      integer i,im1,j,j1,k
      double precision xr,xi,y

c     .......... form the transpose of the upper triangle of inv(l)*a
c                in the lower triangle of the array a ..........
      do 200 i = 1, n
         im1 = i - 1
         y = dl(i)
c
C  (l^-1 * a)_ij  =  (a_ij - sum_k<i  l_ik * (l^-1)_kj) / l_ii
         do 200 j = i, n
            xr = ar(i,j)
            xi = ai(i,j)
            do 160 k = 1, im1
              xr = xr - (br(i,k)*ar(j,k) - bi(i,k)*ai(j,k))
              xi = xi - (br(i,k)*ai(j,k) + bi(i,k)*ar(j,k))
  160       continue
            ar(j,i) = xr / y
            ai(j,i) = xi / y
  200 continue
c     .......... pre-multiply by inv(l) and overwrite ..........
      do 300 j = 1, n
         j1 = j - 1
c
         do 300 i = j, n
            xr =  ar(i,j)
            xi = -ai(i,j)
            im1 = i - 1
            do 220 k = j, im1
c             x = x - a(k,j) * b(i,k)
              xr = xr - (ar(k,j)*br(i,k) - ai(k,j)*bi(i,k))
              xi = xi - (ar(k,j)*bi(i,k) + ai(k,j)*br(i,k))
  220       continue
            do 260 k = 1, j1
c             x = x - a*(j,k) * b(i,k)
              xr = xr - (ar(j,k)*br(i,k) + ai(j,k)*bi(i,k))
              xi = xi - (ar(j,k)*bi(i,k) - ai(j,k)*br(i,k))
  260       continue
            ar(i,j) = xr / dl(i)
            ai(i,j) = xi / dl(i)
  300 continue
      end

