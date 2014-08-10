C#define NEWVERS
C#ifdef NEWVERS
      subroutine diagno(ndim,h,o,wk,iwk,z,eb)

C- Diagonalize secular equation with overlap
C----------------------------------------------------------------------
Ci Inputs
Ci    ndim: dimension of problem
Ci    h,o:  hamiltonian, overlap matrices
Ci    wk, work array length at least 11*ndim; iwk, work array
Co Outputs
Co    z:    eigenvectors; eb, eigenvalues
Cr Remarks
Cr    h,o,z are dimensioned (ndim,ndim,2) (real followed by imag. parts)
Cr    h,o are OVERWRITTEN in this routine
C----------------------------------------------------------------------
C Passed parameters
      integer ndim, iwk(ndim)
      double precision h(ndim,ndim,2),o(ndim,ndim,2),z(ndim,ndim,2),
     .         eb(ndim),wk(ndim,11)

C Local variables:
      integer ierr

C --- Make unit eigenvector matrix ---
      call zinit(z,ndim**2)
      call dvcpy(1.d0,0,z,ndim+1,ndim)

C --- Find the eigenvalues and vectors of O^-1/2  H  O^-1/2 ---
      call bchd(ndim,o,o(1,1,2),wk,ierr)
      if (ierr .ne. 0) stop 'DIAGNO: error in bchd'
      call bred(ndim,h,h(1,1,2),o,o(1,1,2),z,z(1,1,2),wk)
      call btridi(ndim,h,h(1,1,2),wk(1,2),wk(1,3),wk(1,4),wk(1,5))
      call imtqlv(ndim,wk(1,2),wk(1,3),wk(1,4),eb,iwk,ierr,wk(1,7))
      if (ierr .ne. 0) stop 'DIAGNO: error in imtqlv'
      call binvit(ndim,wk(1,2),wk(1,3),wk(1,4),ndim,eb,iwk,
     .           z,ierr,wk(1,7),wk(1,8),wk(1,9),wk(1,10),wk(1,11))
      if (ierr .ne. 0) stop 'DIAGNO: error in binvit'
      call btribk(ndim,h,h(1,1,2),wk(1,5),ndim,z,z(1,1,2))
C --- Get the eigenvectors of H - E O ---
      call bbak(ndim,o,o(1,1,2),wk,ndim,z,z(1,1,2))
      end
      subroutine bchd(n,br,bi,dl,ierr)
C- Cholesky decomposition of hermetian matrix
C ----------------------------------------------------------------
Ci Inputs
Ci   n:  order of the matrix system.  If b has already been Cholesky
Ci       decomposed, n should be prefixed with a minus sign.
Ci   b:  nonorthogonality matrix (only full upper triangle is used)
Ci       If n is negative, both b and dl must be altered
Co Outputs
Co   b:  strict lower triangle of Cholesky factor l of b
Co       The strict upper triangle is the transposed of the lower one
Co   dl: diagonal elements of l.
Co   ierr:nonzero if decomposition finds matrix not pos. def.
Cr Remarks
Cr   Adapted from hchd which is adapted from eispack reduc for
Cr   hermetian matrices. bchd is written for square matrices.
Cr   At the beginning the strict lower triangle is set to the
Cr   transposed of the strict upper one. Then it is possible to
Cr   change the indices in particular for the innermost loop so
Cr   that the first array index equals the loop variable.
C ----------------------------------------------------------------
C Passed parameters
      integer ierr,n
      double precision br(n,n),bi(n,n),dl(n)
C Local parameters
      integer i,im1,j,k
      double precision xr,xi,y

      do 40 j = 1, n
        do 40 i = 1, j
          br(j,i) = br(i,j)
          bi(j,i) = bi(i,j)
   40 continue

c     .......... form l in the arrays b and dl ..........
      do 80 i = 1, n
         im1 = i - 1
         do 80 j = i, n
            xr =  br(j,i)
            xi = -bi(j,i)
            do 20 k = 1, im1
c             x = x - b*(i,k) * b(j,k)
              xr = xr - (br(k,i)*br(k,j) + bi(k,i)*bi(k,j))
              xi = xi - (br(k,i)*bi(k,j) - bi(k,i)*br(k,j))
   20       continue
            if (j .ne. i) go to 60
            ierr = i
            if (xr .le. 0.0d0) return
            y = dsqrt(xr)
            dl(i) = y
            go to 80
   60       br(j,i) = xr / y
            bi(j,i) = xi / y
            br(i,j) = br(j,i)
            bi(i,j) = bi(j,i)
   80 continue
      ierr = 0
      end
      subroutine bred(n,ar,ai,br,bi,atr,ati,dl)
C- Reduction of nonorthogonal hermetian matrix to orthogonal form
C ----------------------------------------------------------------
Ci Inputs
Ci   n:  order of the matrix system.
Ci   a:  hamiltonian matrix(only full upper triangle is used)
Ci   b:  Cholesky-decomposed nonorthogonality matrix
Ci       (strict upper triangle; see hchd)
Ci   at: work array
Ci   dl: diagonal elements of l (see hchd).
Co Outputs
Co   a:  reduced hermetian matrix:  a <-  ldag^-1 a l^-1
Co       The strict upper triangle
Cr Remarks
Cr   Adapted from hred which is adapted from eispack reduc for
Cr   hermetian matrices. bred is written for square matrices.
Cr   At the beginning the strict upper triangle of at is set to the
Cr   transposed of the strict lower one of a. Then it is possible to
Cr   change the indices in particular for the innermost loop so
Cr   that the first array index equals the loop variable.
C ----------------------------------------------------------------
C Passed parameters
      integer n
      double precision ar(n,n),ai(n,n),br(n,n),bi(n,n),dl(n),
     .                 atr(n,n),ati(n,n)
C Local parameters
      integer i,im1,j,j1,k
      double precision xr,xi,y
c
      do 40 j = 1, n
        do 40 i = 1, j
          atr(i,j) = ar(i,j)
          ati(i,j) = -ai(i,j)
   40 continue
c
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
              xr = xr - (br(k,i)*atr(k,j) - bi(k,i)*ati(k,j))
              xi = xi - (br(k,i)*ati(k,j) + bi(k,i)*atr(k,j))
  160       continue
            ar(j,i) = xr / y
            ai(j,i) = xi / y
            atr(i,j) = ar(j,i)
            ati(i,j) = ai(j,i)
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
              xr = xr - (atr(k,j)*br(k,i) - ati(k,j)*bi(k,i))
              xi = xi - (atr(k,j)*bi(k,i) + ati(k,j)*br(k,i))
  220       continue
            do 260 k = 1, j1
c             x = x - a*(j,k) * b(i,k)
              xr = xr - (ar(k,j)*br(k,i) + ai(k,j)*bi(k,i))
              xi = xi - (ar(k,j)*bi(k,i) - ai(k,j)*br(k,i))
  260       continue
            ar(j,i) = xr / dl(i)
            ai(j,i) = xi / dl(i)
            atr(i,j) = ar(j,i)
            ati(i,j) = ai(j,i)
  300 continue
      end
      subroutine bbak(n,br,bi,dl,m,zr,zi)
C- Back-transforms eigenvectors to nonorthogonal representation
C ----------------------------------------------------------------
Ci Inputs
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
      integer m,n
      double precision br(n,n),bi(n,n),dl(n),zr(n,m),zi(n,m)
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
   60       continue
            zr(i,j) = xr/dl(i)
            zi(i,j) = xi/dl(i)
  100 continue
      end
C#elseC
C      subroutine diagno(ndim,h,o,wk,iwk,z,eb)
C
CC- Diagonalize secular equation with overlap
CC----------------------------------------------------------------------
CCi Inputs
CCi    ndim: dimension of problem
CCi    h,o:  hamiltonian, overlap matrices
CCi    wk, work array length at least 5*ndim
CCo Outputs
CCo    z:    eigenvectors; eb, eigenvalues
CCr Remarks
CCr    h,o,z are dimensioned (ndim,ndim,2) (real followed by imag. parts
CCr    h,o are OVERWRITTEN in this routine
CC----------------------------------------------------------------------
CC Passed parameters
C      integer ndim, iwk(ndim)
C      double precision h(ndim,ndim,2),o(ndim,ndim,2),z(ndim,ndim,2),
C     .         eb(ndim),wk(ndim,5)
C
CC Local variables:
C      integer ierr
C
CC --- Make unit eigenvector matrix ---
C      call zinit(z,ndim**2)
C      call dvcpy(1.d0,0,z,ndim+1,ndim)
C
CC --- Find the eigenvalues and vectors of O^-1/2  H  O^-1/2 ---
C      call hchd(ndim,ndim,o,o(1,1,2),wk,ierr)
C      if (ierr .ne. 0) stop 'DIAGNO: error in hchd'
C      call hred(ndim,ndim,h,h(1,1,2),o,o(1,1,2),wk)
C      call htridi(ndim,ndim,h,h(1,1,2),eb,wk(1,2),wk(1,3),wk(1,4))
C      call imtql2(ndim,ndim,eb,wk(1,2),z,ierr)
C      if (ierr .ne. 0) stop 'DIAGNO: error in imtql2'
C      call htribk(ndim,ndim,h,h(1,1,2),wk(1,4),ndim,z,z(1,1,2))
C
CC --- Get the eigenvectors of H - E O ---
C      call hbak(ndim,ndim,o,o(1,1,2),wk,ndim,z,z(1,1,2))
C      end
C      subroutine hchd(nm,n,br,bi,dl,ierr)
CC- Cholesky decomposition of hermetian matrix
CC ----------------------------------------------------------------
CCi Inputs
CCi   nm: true row dimension of arrays b and z, as declared by caller
CCi   n:  order of the matrix system.  If b has already been Cholesky
CCi       decomposed, n should be prefixed with a minus sign.
CCi   b:  nonorthogonality matrix (only full upper triangle is used)
CCi       If n is negative, both b and dl must be altered
CCo Outputs
CCo   b:  strict lower triangle of Cholesky factor l of b
CCo       The strict upper triangle unaltered.
CCo   dl: diagonal elements of l.
CCo   ierr:nonzero if decomposition finds matrix not pos. def.
CCr Remarks
CCr   Adapted from eispack reduc for hermetian matrices
CC ----------------------------------------------------------------
CC Passed parameters
C      integer ierr,n,nm
C      double precision br(nm,n),bi(nm,n),dl(n)
CC Local parameters
C      integer i,im1,j,k
C      double precision xr,xi,y
C
Cc     .......... form l in the arrays b and dl ..........
C      do 80 i = 1, n
C         im1 = i - 1
C         do 80 j = i, n
C            xr =  br(i,j)
C            xi = -bi(i,j)
C            do 20 k = 1, im1
Cc             x = x - b*(i,k) * b(j,k)
C              xr = xr - (br(i,k)*br(j,k) + bi(i,k)*bi(j,k))
C              xi = xi - (br(i,k)*bi(j,k) - bi(i,k)*br(j,k))
C   20       continue
C            if (j .ne. i) go to 60
C            ierr = i
C            if (xr .le. 0.0d0) return
C            y = dsqrt(xr)
C            dl(i) = y
C            go to 80
C   60       br(j,i) = xr / y
C            bi(j,i) = xi / y
C   80 continue
C      ierr = 0
C      end
C      subroutine hred(nm,n,ar,ai,br,bi,dl)
CC- Reduction of nonorthogonal hermetian matrix to orthogonal form
CC ----------------------------------------------------------------
CCi Inputs
CCi   nm: true row dimension of arrays b and z, as declared by caller
CCi   n:  order of the matrix system.
CCi   a:  hamiltonian matrix (only full upper triangle is used)
CCi   b:  Cholesky-decomposed nonorthogonality matrix
CCi       (strict lower triangle; see hchd)
CCi   dl: diagonal elements of l (see hchd).
CCo Outputs
CCo   a:  reduced hermetian matrix:  a <-  ldag^-1 a l^-1
CCo       The strict upper triangle unaltered.
CCr Remarks
CCr   Adapted from eispack reduc for hermetian matrices
CC ----------------------------------------------------------------
CC Passed parameters
C      integer n,nm
C      double precision ar(nm,n),ai(nm,n),br(nm,n),bi(nm,n),dl(n)
CC Local parameters
C      integer i,im1,j,j1,k
C      double precision xr,xi,y
C
Cc     .......... form the transpose of the upper triangle of inv(l)*a
Cc                in the lower triangle of the array a ..........
C      do 200 i = 1, n
C         im1 = i - 1
C         y = dl(i)
Cc
CC  (l^-1 * a)_ij  =  (a_ij - sum_k<i  l_ik * (l^-1)_kj) / l_ii
C         do 200 j = i, n
C            xr = ar(i,j)
C            xi = ai(i,j)
C            do 160 k = 1, im1
C              xr = xr - (br(i,k)*ar(j,k) - bi(i,k)*ai(j,k))
C              xi = xi - (br(i,k)*ai(j,k) + bi(i,k)*ar(j,k))
C  160       continue
C            ar(j,i) = xr / y
C            ai(j,i) = xi / y
C  200 continue
Cc     .......... pre-multiply by inv(l) and overwrite ..........
C      do 300 j = 1, n
C         j1 = j - 1
Cc
C         do 300 i = j, n
C            xr =  ar(i,j)
C            xi = -ai(i,j)
C            im1 = i - 1
C            do 220 k = j, im1
Cc             x = x - a(k,j) * b(i,k)
C              xr = xr - (ar(k,j)*br(i,k) - ai(k,j)*bi(i,k))
C              xi = xi - (ar(k,j)*bi(i,k) + ai(k,j)*br(i,k))
C  220       continue
C            do 260 k = 1, j1
Cc             x = x - a*(j,k) * b(i,k)
C              xr = xr - (ar(j,k)*br(i,k) + ai(j,k)*bi(i,k))
C              xi = xi - (ar(j,k)*bi(i,k) - ai(j,k)*br(i,k))
C  260       continue
C            ar(i,j) = xr / dl(i)
C            ai(i,j) = xi / dl(i)
C  300 continue
C      end
C      subroutine hbak(nm,n,br,bi,dl,m,zr,zi)
CC- Back-transforms eigenvectors to nonorthogonal representation
CC ----------------------------------------------------------------
CCi Inputs
CCi   nm: true row dimension of arrays b and z
CCi   n:  order of the matrix system.
CCi   b:  Cholesky-decomposed nonorthogonality matrix, as decomposed
CCi       by reduch (in its strict lower triangle).
CCi   dl: contains further information about the transformation.
CCi   m:  number of eigenvectors to be back transformed.
CCi   z:  eigenvectors to be back transformed (first m columns)
CCo Outputs
CCo   z transformed eigenvectors
CCr Remarks
CCr   Adapted from rebakh, eispack
CC ----------------------------------------------------------------
CC Passed parameters
C      integer m,n,nm
C      double precision br(nm,n),bi(nm,n),dl(n),zr(nm,m),zi(nm,m)
CC Local parameters
C      integer i,j,k
C      double precision xr,xi
C
C      do 100 j = 1, m
C         do 100 i = n, 1, -1
C            xr = zr(i,j)
C            xi = zi(i,j)
C            do 60 k = i+1, n
C              xr = xr - (br(k,i)*zr(k,j) + bi(k,i)*zi(k,j))
C              xi = xi - (br(k,i)*zi(k,j) - bi(k,i)*zr(k,j))
C   60       continue
C            zr(i,j) = xr/dl(i)
C            zi(i,j) = xi/dl(i)
C  100 continue
C      end
C#endif
      subroutine hcinv(nm,n,br,bi,dl,lopt)
C- Obtain inverse from Cholesky decomposed matrix
C ----------------------------------------------------------------
Ci Inputs
Ci   nm: true row dimension of arrays b and z, as declared by caller
Ci   n:  order of the matrix system.
Ci   b:  Cholesky-decomposed nonorthogonality matrix
Ci       (strict lower triangle; see hchd)
Ci   dl: diagonal elements of l (see hchd).
Ci   lopt: false, only upper half of b\-1 is generated, and lower
Ci         half continues to hold the c.d. of b
Ci         true,  lower half of b also filled in
Co Outputs
Co   b:  Inverse of Cholesky matrix
Cr Remarks
Cr   Adapted from eispack reduc for hermetian matrices
C ----------------------------------------------------------------
C Passed parameters
      integer n,nm
      logical lopt
      double precision br(nm,n),bi(nm,n),dl(n)
C Local parameters
      integer i,im1,j,jm1,k
      double precision xr,xi,y
C#ifdefC BLAS
C      double precision ddot
C#endif

C --- form inv(l) and store transpose in full upper triangle of b ---
C  (l\-1)_ij  =  (del_ij - sum_k<i  l_ik * (l\-1)_kj) / l_ii
      do 200 j = 1, n
        do 200 i = j, n
          y = dl(i)
C#ifdefC BLAS
C          if (j .eq. i) then
C            xr = 1d0
C            xi = 0
C          else
C            xr = - ddot(i-j,br(i,j),nm,br(j,j),nm)
C     .           + ddot(i-j,bi(i,j),nm,bi(j,j),nm)
C            xi = - ddot(i-j,br(i,j),nm,bi(j,j),nm)
C     .           - ddot(i-j,bi(i,j),nm,br(j,j),nm)
C          endif
C#else
          xr = 0
          xi = 0
          if (j .eq. i) xr = 1d0
          im1 = i-1
          do 160 k = j, im1
            xr = xr - (br(i,k)*br(j,k) - bi(i,k)*bi(j,k))
            xi = xi - (br(i,k)*bi(j,k) + bi(i,k)*br(j,k))
  160     continue
C#endif
          br(j,i) = xr / y
          bi(j,i) = xi / y
  200 continue
c
c     .......... pre-multiply by inv(l) and overwrite ..........
      do 300 i = 1, n
         do 300 j = i, n
C           x = x - l\dag_ik * l_kj; here i <= j <= k
C#ifdefC BLAS
C           br(i,j) = ddot(n+1-j,br(i,j),nm,br(j,j),nm) +
C     .               ddot(n+1-j,bi(i,j),nm,bi(j,j),nm)
C           bi(i,j) = ddot(n+1-j,br(i,j),nm,bi(j,j),nm) -
C     .               ddot(n+1-j,bi(i,j),nm,br(j,j),nm)
C#else
            xr = 0
            xi = 0
            do 260 k = j, n
              xr = xr + (br(i,k)*br(j,k) + bi(i,k)*bi(j,k))
              xi = xi + (br(i,k)*bi(j,k) - bi(i,k)*br(j,k))
  260       continue
            br(i,j) = xr
            bi(i,j) = xi
C#endif
  300 continue

      if (.not. lopt) return
      do  400  i = 1, n
       do  400  j = i+1, n
         br(j,i) =  br(i,j)
         bi(j,i) = -bi(i,j)
  400 continue
      end

      subroutine yhifa(ar,ai,lda,n,kpvt,info)
C- Factor an hermetian matrix with symmetric pivoting
      integer lda,n,kpvt(1),info
      double precision ar(lda,1),ai(lda,1)
c
c     to solve  a*x = b , follow yhifa by yhisl.
c     to compute  inverse(a)*c , follow yhifa by yhisl.
c     to compute  determinant(a) , follow yhifa by yhidi.
c     to compute  inertia(a) , follow yhifa by yhidi.
c     to compute  inverse(a) , follow yhifa by yhidi.
c
c     on entry
c
c        ar,ai   the hermitian matrix to be factored.
c                only the diagonal and upper triangle are used.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       a block diagonal matrix and the multipliers which
c                were used to obtain it.
c                the factorization can be written  a = u*d*ctrans(u)
c                where  u  is a product of permutation and unit
c                upper triangular matrices , ctrans(u) is the
c                conjugate transpose of  u , and  d  is block diagonal
c                with 1 by 1 and 2 by 2 blocks.
c
c        kpvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if the k-th pivot block is singular. this is
c                     not an error condition for this subroutine,
c                     but it does indicate that yhisl or yhidi may
c                     divide by zero if called.
c
c     Adapted from linpack's version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     internal variables
c
      double precision mulk(2),mulkm1(2),denom(2),ak(2),akm1(2),
     .  bk(2),bkm1(2),t(2)

      double precision absakk,alpha,colmax,rowmax,dcabs1
      integer imax,imaxp1,j,jj,jmax,k,km1,km2,kstep,iyamax
      logical swap

c
C --- Initialize ---
c     alpha is used in choosing pivot block size.
      alpha = (1.0d0 + dsqrt(17.0d0))/8.0d0
      info = 0

C --- Main loop on k, which goes from n to 1. ---
      k = n
   10 continue

C ---    Leave the loop if k=0 or k=1. ---
        if (k .eq. 0) go to 200
        if (k .gt. 1) go to 20
          kpvt(1) = 1
          if (dcabs1(ar,ai) .eq. 0.0d0) info = 1
          go to 200
   20   continue
c
c     This section of code determines the kind of
c     elimination to be performed.  when it is completed,
c     kstep will be set to the size of the pivot block, and
c     swap will be set to .true. if an interchange is
c     required.

        km1 = k - 1
        absakk = dcabs1(ar(k,k),ai(k,k))

C ---  Determine the largest off-diagonal element in column k ---
        imax = iyamax(k-1,ar(1,k),ai(1,k),1)
        colmax = dcabs1(ar(imax,k),ai(imax,k))
        if (absakk .lt. alpha*colmax) go to 30
          kstep = 1
          swap = .false.
          go to 90
   30   continue

C ---  Determine the largest off-diagonal element in row imax. ---
        rowmax = 0.0d0
        imaxp1 = imax + 1
        do 40 j = imaxp1, k
          rowmax = dmax1(rowmax,dcabs1(ar(imax,j),ai(imax,j)))
   40   continue
        if (imax .eq. 1) go to 50
          jmax = iyamax(imax-1,ar(1,imax),ai(1,imax),1)
          rowmax=dmax1(rowmax,dcabs1(ar(jmax,imax),ai(jmax,imax)))
   50   continue
        if (dcabs1(ar(imax,imax),ai(imax,imax)) .lt. alpha*rowmax)
     .    go to 60
          kstep = 1
          swap = .true.
          go to 80
   60   continue
        if (absakk .lt. alpha*colmax*(colmax/rowmax)) go to 70
          kstep = 1
          swap = .false.
          go to 80
   70   continue
        kstep = 2
        swap = imax .ne. km1
   80   continue
   90   continue
        if (dmax1(absakk,colmax) .ne. 0.0d0) go to 100

C ---   Column k is zero.  set info and iterate the loop. ---
        kpvt(k) = k
        info = k
        go to 190
  100   continue
        if (kstep .eq. 2) go to 140

C ---      1 x 1 pivot block. ---
          if (.not.swap) go to 120

C ---       Perform an interchange. ---
            call dswap(imax,ar(1,imax),1,ar(1,k),1)
            call dswap(imax,ai(1,imax),1,ai(1,k),1)
            do 110 jj = imax, k
              j = k + imax - jj
              t(1) = ar(j,k)
              t(2) = -ai(j,k)
              ar(j,k) =  ar(imax,j)
              ai(j,k) = -ai(imax,j)
              ar(imax,j) = t(1)
              ai(imax,j) = t(2)
  110       continue
  120       continue

C ---       Perform the elimination. ---
            do 130 jj = 1, km1
              j = k - jj
C     complex divide ... mulk = -a(j,k)/a(k,k) ...
              call cdiv(-ar(j,k),-ai(j,k),
     .          ar(k,k),ai(k,k),mulk(1),mulk(2))
c     call zaxpy(j,mulk*,ar(1,k),1,ar(1,j),1)
              call daxpy(j,mulk(1),ar(1,k),1,ar(1,j),1)
              call daxpy(j,mulk(2),ai(1,k),1,ar(1,j),1)
              call daxpy(j,mulk(1),ai(1,k),1,ai(1,j),1)
              call daxpy(j,-mulk(2),ar(1,k),1,ai(1,j),1)
              ai(j,j) = 0
              ar(j,k) = mulk(1)
              ai(j,k) = mulk(2)
  130       continue

C --- Set the pivot array. ---
            kpvt(k) = k
            if (swap) kpvt(k) = imax
            go to 190
  140     continue

C --- 2 x 2 pivot block. ---
          if (.not.swap) go to 160
c
C ---      Perform an interchange. ---
            call dswap(imax,ar(1,imax),1,ar(1,k-1),1)
            call dswap(imax,ai(1,imax),1,ai(1,k-1),1)
            do 150 jj = imax, km1
              j = km1 + imax - jj
              t(1) = ar(j,k-1)
              t(2) = -ai(j,k-1)
              ar(j,k-1) =  ar(imax,j)
              ai(j,k-1) = -ai(imax,j)
              ar(imax,j) = t(1)
              ai(imax,j) = t(2)
  150       continue
            t(1) = ar(k-1,k)
            t(2) = ai(k-1,k)
            ar(k-1,k) = ar(imax,k)
            ai(k-1,k) = ai(imax,k)
            ar(imax,k) = t(1)
            ai(imax,k) = t(2)
  160     continue

C ---       Perform the elimination. ---
          km2 = k - 2
          if (km2 .eq. 0) go to 180
c ... ak = a(k,k)/a(k-1,k),  akm1 = a(k-1,k-1)/dconjg(a(k-1,k)) ...
            call cdiv(ar(k,k),ai(k,k),ar(k-1,k),ai(k-1,k),ak(1),ak(2))
            call cdiv(ar(k-1,k-1),ai(k-1,k-1),ar(k-1,k),-ai(k-1,k),
     .        akm1(1),akm1(2))
c ... denom = 1.0d0 - ak*akm1 ...
            denom(1) = 1 - ak(1)*akm1(1) + ak(2)*akm1(2)
            denom(2) =   - ak(1)*akm1(2) - ak(2)*akm1(1)
            do 170 jj = 1, km2
              j = km1 - jj
c ... bk = a(j,k)/a(k-1,k),  bkm1 = a(j,k-1)/dconjg(a(k-1,k))
              call cdiv(ar(j,k),ai(j,k),
     .          ar(k-1,k),ai(k-1,k),bk(1),bk(2))
              call cdiv(ar(j,k-1),ai(j,k-1),
     .          ar(k-1,k),-ai(k-1,k),bkm1(1),bkm1(2))
c ... mulk = (akm1*bk - bkm1)/denom
              call cpy(akm1(1),akm1(2),bk(1),bk(2),mulk(1),mulk(2))
              mulk(1) = mulk(1) - bkm1(1)
              mulk(2) = mulk(2) - bkm1(2)
              call cdiv(mulk(1),mulk(2),denom(1),denom(2),
     .          mulk(1),mulk(2))
c ... mulkm1 = (ak*bkm1 - bk)/denom
              call cpy(ak(1),ak(2),bkm1(1),bkm1(2),mulkm1(1),mulkm1(2))
              mulkm1(1) = mulkm1(1) - bk(1)
              mulkm1(2) = mulkm1(2) - bk(2)
              call cdiv(mulkm1(1),mulkm1(2),denom(1),denom(2),
     .          mulkm1(1),mulkm1(2))
c     call zaxpy(j,mulk*,ar(1,k),1,ar(1,j),1)
              call daxpy(j,mulk(1),ar(1,k),1,ar(1,j),1)
              call daxpy(j,mulk(2),ai(1,k),1,ar(1,j),1)
              call daxpy(j,mulk(1),ai(1,k),1,ai(1,j),1)
              call daxpy(j,-mulk(2),ar(1,k),1,ai(1,j),1)
c     call zaxpy(j,mulkm1*,ar(1,k-1),1,ar(1,j),1)
              call daxpy(j,mulkm1(1),ar(1,k-1),1,ar(1,j),1)
              call daxpy(j,mulkm1(2),ai(1,k-1),1,ar(1,j),1)
              call daxpy(j,mulkm1(1),ai(1,k-1),1,ai(1,j),1)
              call daxpy(j,-mulkm1(2),ar(1,k-1),1,ai(1,j),1)
              ar(j,k) = mulk(1)
              ai(j,k) = mulk(2)
              ar(j,k-1) = mulkm1(1)
              ai(j,k-1) = mulkm1(2)
              ai(j,j) = 0
  170       continue
  180     continue

C ---    Set the pivot array. ---
          kpvt(k) = 1 - k
          if (swap) kpvt(k) = -imax
          kpvt(k-1) = kpvt(k)
  190   continue
        k = k - kstep
        go to 10
  200 continue
      end
      subroutine yhidi(ar,ai,lda,n,kpvt,det,inert,work,job)
C- Compute determinant, inertia and inverse using factors from yhifa
      integer lda,n,job
      double precision ar(lda,lda),ai(lda,lda),work(2,1)
      double precision det(2)
      integer kpvt(1),inert(3)
c
c     on entry
c
c        ar,ai   output from yhifa
c
c        lda     integer
c                the leading dimension of the array a.
c
c        n       integer
c                the order of the matrix a.
c
c        kpvt    integer(n)
c                the pivot vector from yhifa.
c
c        work    complex*16(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                job has the decimal expansion  abc  where
c                   if  c .ne. 0, the inverse is computed,
c                   if  b .ne. 0, the determinant is computed,
c                   if  a .ne. 0, the inertia is computed.
c
c                for example, job = 111  gives all three.
c
c     on return
c
c        variables not requested by job are not used.
c
c        a      contains the upper triangle of the inverse of
c               the original matrix.  the strict lower triangle
c               is never referenced.
c
c        det    double precision(2)
c               determinant of original matrix.
c               determinant = det(1) * 10.0**det(2)
c               with 1.0 .le. dabs(det(1)) .lt. 10.0
c               or det(1) = 0.0.
c
c        inert  integer(3)
c               the inertia of the original matrix.
c               inert(1)  =  number of positive eigenvalues.
c               inert(2)  =  number of negative eigenvalues.
c               inert(3)  =  number of zero eigenvalues.
c
c     error condition
c
c        a division by zero may occur if the inverse is requested
c        and  yhico  has set rcond .eq. 0.0
c        or  yhifa  has set  info .ne. 0 .
c
c     Adapted from linpack's version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab
c
c     internal variables.
c
      double precision ten,d,t,ak,akp1,akkp1(2),temp,cdabs2,ddot
      integer j,jb,k,km1,ks,kstep
      logical noinv,nodet,noert
c
      noinv = mod(job,10) .eq. 0
      nodet = mod(job,100)/10 .eq. 0
      noert = mod(job,1000)/100 .eq. 0
c
      if (nodet .and. noert) go to 140
        if (noert) go to 10
          inert(1) = 0
          inert(2) = 0
          inert(3) = 0
   10   continue
        if (nodet) go to 20
          det(1) = 1.0d0
          det(2) = 0.0d0
          ten = 10.0d0
   20   continue
        t = 0.0d0
        do 130 k = 1, n
          d = ar(k,k)

C ---       Check if 1 by 1 ---
          if (kpvt(k) .gt. 0) go to 50

c     2 by 2 block
c     use det (d  s)  =  (d/t * c - t) * t  ,  t = cdabs(s)
c                        (s  c)
c     to avoid underflow/overflow troubles.
c     take two passes through scaling.  use  t  for flag.
c
            if (t .ne. 0.0d0) go to 30
              t = cdabs2(ar(k,k+1),ai(k,k+1))
              d = (d/t)*ar(k+1,k+1) - t
              go to 40
   30       continue
            d = t
            t = 0.0d0
   40       continue
   50     continue
c
          if (noert) go to 60
            if (d .gt. 0.0d0) inert(1) = inert(1) + 1
            if (d .lt. 0.0d0) inert(2) = inert(2) + 1
            if (d .eq. 0.0d0) inert(3) = inert(3) + 1
   60     continue
c
          if (nodet) go to 120
            det(1) = d*det(1)
            if (det(1) .eq. 0.0d0) go to 110
   70       if (dabs(det(1)) .ge. 1.0d0) go to 80
              det(1) = ten*det(1)
              det(2) = det(2) - 1.0d0
              go to 70
   80       continue
   90       if (dabs(det(1)) .lt. ten) go to 100
              det(1) = det(1)/ten
              det(2) = det(2) + 1.0d0
              go to 90
  100       continue
  110       continue
  120       continue
  130    continue
  140 continue
c
C --- Compute inverse(a) ---
      if (noinv) go to 270
      k = 1
  150 if (k .gt. n) go to 260
        km1 = k - 1
        if (kpvt(k) .lt. 0) go to 180
c
C --- 1 by 1 ---
        ar(k,k) = 1/ar(k,k)
        ai(k,k) = 0
        if (km1 .lt. 1) go to 170
          call dcopy(km1,ar(1,k),1,work,2)
          call dcopy(km1,ai(1,k),1,work(2,1),2)
          do 160 j = 1, km1
            ar(j,k) = ddot(j,ar(1,j),1,work,2)
     .        +ddot(j,ai(1,j),1,work(2,1),2)
            ai(j,k) = ddot(j,ar(1,j),1,work(2,1),2)
     .        -ddot(j,ai(1,j),1,work,2)
c     call zaxpy(j-1,work(1,j),ar(1,j),1,ar(1,k),1)
            call daxpy(j-1,work(1,j),ar(1,j),1,ar(1,k),1)
            call daxpy(j-1,-work(2,j),ai(1,j),1,ar(1,k),1)
            call daxpy(j-1,work(1,j),ai(1,j),1,ai(1,k),1)
            call daxpy(j-1,work(2,j),ar(1,j),1,ai(1,k),1)
  160     continue
          ar(k,k) = ar(k,k)
     .      +ddot(km1,work,2,ar(1,k),1)
     .      +ddot(km1,work(2,1),2,ai(1,k),1)
          ai(k,k) = 0
  170   continue
        kstep = 1
        go to 220
  180   continue

C  ---2 by 2 ---
        t = cdabs2(ar(k,k+1),ai(k,k+1))
        ak = ar(k,k)/t
        akp1 = ar(k+1,k+1)/t
c ... akkp1 = a(k,k+1)/t
        akkp1(1) = ar(k,k+1)/t
        akkp1(2) = ai(k,k+1)/t
        d = t*(ak*akp1 - 1.0d0)
        ar(k,k) = akp1/d
        ai(k,k) = 0
        ar(k+1,k+1) = ak/d
        ai(k+1,k+1) = 0
        ar(k,k+1) = -akkp1(1)/d
        ai(k,k+1) = -akkp1(2)/d
        if (km1 .lt. 1) go to 210
          call dcopy(km1,ar(1,k+1),1,work,2)
          call dcopy(km1,ai(1,k+1),1,work(2,1),2)
          do 190 j = 1, km1
            ar(j,k+1) =
     .        ddot(j,ar(1,j),1,work,2) +
     .        ddot(j,ai(1,j),1,work(2,1),2)
            ai(j,k+1) =
     .        ddot(j,ar(1,j),1,work(2,1),2) -
     .        ddot(j,ai(1,j),1,work,2)
c     call zaxpy(j-1,work(1,j),ar(1,j),1,ar(1,k+1),1)
            call daxpy(j-1,work(1,j),ar(1,j),1,ar(1,k+1),1)
            call daxpy(j-1,-work(2,j),ai(1,j),1,ar(1,k+1),1)
            call daxpy(j-1,work(1,j),ai(1,j),1,ai(1,k+1),1)
            call daxpy(j-1,work(2,j),ar(1,j),1,ai(1,k+1),1)
  190     continue
          ar(k+1,k+1) = ar(k+1,k+1)
     .      +ddot(km1,work,2,ar(1,k+1),1)
     .      +ddot(km1,work(2,1),2,ai(1,k+1),1)
          ai(k+1,k+1) = 0
          ar(k,k+1) = ar(k,k+1)
     .      +ddot(km1,ar(1,k),1,ar(1,k+1),1)
     .      +ddot(km1,ai(1,k),1,ai(1,k+1),1)
          ai(k,k+1) = ai(k,k+1)
     .      +ddot(km1,ar(1,k),1,ai(1,k+1),1)
     .      -ddot(km1,ai(1,k),1,ar(1,k+1),1)
          call dcopy(km1,ar(1,k),1,work,2)
          call dcopy(km1,ai(1,k),1,work(2,1),2)
          do 200 j = 1, km1
            ar(j,k) =
     .        ddot(j,ar(1,j),1,work,2)
     .        +ddot(j,ai(1,j),1,work(2,1),2)
            ai(j,k) =
     .        ddot(j,ar(1,j),1,work(2,1),2)
     .        -ddot(j,ai(1,j),1,work,2)
c     call zaxpy(j-1,work(1,j),ar(1,j),1,ar(1,k),1)
            call daxpy(j-1,work(1,j),ar(1,j),1,ar(1,k),1)
            call daxpy(j-1,-work(2,j),ai(1,j),1,ar(1,k),1)
            call daxpy(j-1,work(1,j),ai(1,j),1,ai(1,k),1)
            call daxpy(j-1,work(2,j),ar(1,j),1,ai(1,k),1)
  200     continue
          ar(k,k) = ar(k,k)
     .      +ddot(km1,work,2,ar(1,k),1)
     .      +ddot(km1,work(2,1),2,ai(1,k),1)
          ai(k,k) = 0
  210   continue
        kstep = 2
  220   continue

C --- Swap ---
        ks = iabs(kpvt(k))
        if (ks .eq. k) go to 250
          call dswap(ks,ar(1,ks),1,ar(1,k),1)
          call dswap(ks,ai(1,ks),1,ai(1,k),1)
          do 230 jb = ks, k
            j = k + ks - jb
            temp = ar(j,k)
            ar(j,k) = ar(ks,j)
            ar(ks,j) = temp
            temp = -ai(j,k)
            ai(j,k) = -ai(ks,j)
            ai(ks,j) = temp
  230     continue
          if (kstep .eq. 1) go to 240
            temp = ar(ks,k+1)
            ar(ks,k+1) = ar(k,k+1)
            ar(k,k+1) = temp
            temp = ai(ks,k+1)
            ai(ks,k+1) = ai(k,k+1)
            ai(k,k+1) = temp
  240     continue
  250     continue
          k = k + kstep
          go to 150
  260   continue
  270 continue
      end
      integer function iyamax(n,zxr,zxi,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, 1/15/85.
c
      double precision zxr(1),zxi(1)
      double precision smax
      double precision cdabs2
c
      iyamax = 1
      if(n.le.1)return
      if(incx.eq.1)go to 20
c
c        code for increments not equal to 1
c
      ix = 1
      smax = cdabs2(zxr,zxi)
      ix = ix + incx
      do 10 i = 2,n
         if(cdabs2(zxr(ix),zxi(ix)).le.smax) go to 5
         iyamax = i
         smax = cdabs2(zxr(ix),zxi(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increments equal to 1
c
   20 smax = dcabs1(zxr,zxi)
      do 30 i = 2,n
         if(cdabs2(zxr(i),zxi(i)).le.smax) go to 30
         iyamax = i
         smax = cdabs2(zxr(i),zxi(i))
   30 continue
      return
      end
      double precision function dcabs1(zr,zi)
      double precision zr,zi
      dcabs1 = dabs(zr) + dabs(zi)
      end
      double precision function cdabs2(x1,x2)
      double precision x1,x2
      cdabs2 = dsqrt(x1**2 + x2**2)
      end
      subroutine cpy(tr,ti,dr,di,t1,t2)
C- complex multiply (t1,t2) = (tr,ti) * (dr,di)
      double precision tr,ti,dr,di,t1,t2
      double precision tmp
      tmp = tr*dr - ti*di
      t2  = tr*di + ti*dr
      t1 = tmp
      end
      subroutine cdiv(tr,ti,dr,di,t1,t2)
C- complex divide (t1,t2) = (tr,ti) / (dr,di)
Cr Remarks
Cr   Adapted from eispack.
Cr   It is permissible for (t1,t2) to occupy the same address space
Cr   as either (tr,ti) or (dr,di)
      double precision tr,ti,dr,di,t1,t2
      double precision rr,d,tmp

      if (dabs(di) .gt. dabs(dr)) then
        rr = dr / di
        d = dr * rr + di
        tmp = (tr * rr + ti) / d
        t2 = (ti * rr - tr) / d
        t1 = tmp
      else
        rr = di / dr
        d = dr + di * rr
        tmp = (tr + ti * rr) / d
        t2 = (ti - tr * rr) / d
        t1 = tmp
      endif
      end

C#ifndef SLATEC
C#ifdef NEWVERS
      subroutine btridi(n,ar,ai,d,e,e2,tau)
c
      integer i,j,k,l,n,ii,jp1
      double precision ar(n,n),ai(n,n),d(n),e(n),e2(n),tau(2,n)
      double precision f,g,h,fi,gi,hh,si,scale,pythag
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure tred1, num. math. 11, 181-195(1968)
c     by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a complex hermitian matrix
c     to a real symmetric tridiagonal matrix using
c     unitary similarity transformations.
c
c     on input
c
c        n is the order of the matrix.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex hermitian input matrix.
c          only the upper triangle of the matrix need be supplied.
c
c     on output
c
c        ar and ai contain information about the unitary trans-
c          formations used in the reduction in their full upper
c          triangles.  their strict lower triangles and the
c          diagonal of ar are unaltered.
c
c        d contains the diagonal elements of the the tridiagonal matrix
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c
c        tau contains further information about the transformations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version is adapted from the eispack htridi. In contrast
c     only the upper triangle of a is used for input and output
c
c     -----------------------------------------------------------------
c
      tau(1,n) = 1.0d0
      tau(2,n) = 0.0d0
c
      do 100 i = 1, n
  100 d(i) = ar(i,i)
c     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 1) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(ar(k,i)) + dabs(ai(k,i))
c
         if (scale .ne. 0.0d0) go to 140
         tau(1,l) = 1.0d0
         tau(2,l) = 0.0d0
  130    e(i) = 0.0d0
         e2(i) = 0.0d0
         go to 290
c
  140    do 150 k = 1, l
            ar(k,i) = ar(k,i) / scale
            ai(k,i) = ai(k,i) / scale
            h = h + ar(k,i) * ar(k,i) + ai(k,i) * ai(k,i)
  150    continue
c
         e2(i) = scale * scale * h
         g = dsqrt(h)
         e(i) = scale * g
         f = pythag(ar(l,i),ai(l,i))
c     .......... form next diagonal element of matrix t ......... .
         if (f .eq. 0.0d0) go to 160
         tau(1,l) = (ai(l,i) * tau(2,i) - ar(l,i) * tau(1,i)) / f
         si = (ar(l,i) * tau(2,i) + ai(l,i) * tau(1,i)) / f
         h = h + f * g
         g = 1.0d0 + g / f
         ar(l,i) = g * ar(l,i)
         ai(l,i) = g * ai(l,i)
         if (l .eq. 1) go to 270
         go to 170
  160    tau(1,l) = -tau(1,i)
         si = tau(2,i)
         ar(l,i) = g
  170    f = 0.0d0
c
         do 240 j = 1, l
            g = 0.0d0
            gi = 0.0d0
c     .......... form element of a*u ..........
            do 180 k = 1, j
               g = g + ar(k,j) * ar(k,i) + ai(k,j) * ai(k,i)
               gi = gi - ar(k,j) * ai(k,i) + ai(k,j) * ar(k,i)
  180       continue
c
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + ar(j,k) * ar(k,i) - ai(j,k) * ai(k,i)
               gi = gi - ar(j,k) * ai(k,i) - ai(j,k) * ar(k,i)
  200       continue
c     .......... form element of p ..........
  220       e(j) = g / h
            tau(2,j) = gi / h
            f = f + e(j) * ar(j,i) - tau(2,j) * ai(j,i)
  240    continue
c
         hh = f / (h + h)
c     .......... form reduced a ..........
         do 260 j = 1, l
            f = ar(j,i)
            g = e(j) - hh * f
            e(j) = g
            fi = -ai(j,i)
            gi = tau(2,j) - hh * fi
            tau(2,j) = -gi
c
            do 260 k = 1, j
               ar(k,j) = ar(k,j) - f * e(k) - g * ar(k,i)
     x                           + fi * tau(2,k) + gi * ai(k,i)
               ai(k,j) = ai(k,j) - f * tau(2,k) - g * ai(k,i)
     x                           - fi * e(k) - gi * ar(k,i)
  260    continue
c
  270    do 280 k = 1, l
            ar(k,i) = scale * ar(k,i)
            ai(k,i) = scale * ai(k,i)
  280    continue
c
         tau(2,l) = -si
  290    hh = d(i)
         d(i) = ar(i,i)
         ar(i,i) = hh
         ai(i,i) = scale * dsqrt(h)
  300 continue
c
      return
      end
      SUBROUTINE IMTQLV(N,D,E,E2,W,IND,IERR,RV1)
C
      INTEGER I,J,K,L,M,N,II,MML,TAG,IERR
      DOUBLE PRECISION D(N),E(N),E2(N),W(N),RV1(N)
      DOUBLE PRECISION B,C,F,G,P,R,S,TST1,TST2,PYTHAG
      INTEGER IND(N)
C
C     THIS ROUTINE IS A VARIANT OF  IMTQL1  WHICH IS A TRANSLATION OF
C     ALGOL PROCEDURE IMTQL1, NUM. MATH. 12, 377-383(1968) BY MARTIN AND
C     WILKINSON, AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
C
C     THIS ROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC TRIDIAGONAL
C     MATRIX BY THE IMPLICIT QL METHOD AND ASSOCIATES WITH THEM
C     THEIR CORRESPONDING SUBMATRIX INDICES.
C
C     ON INPUT
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2(1) IS ARBITRARY.
C
C     ON OUTPUT
C
C        D AND E ARE UNALTERED.
C
C        ELEMENTS OF E2, CORRESPONDING TO ELEMENTS OF E REGARDED
C          AS NEGLIGIBLE, HAVE BEEN REPLACED BY ZERO CAUSING THE
C          MATRIX TO SPLIT INTO A DIRECT SUM OF SUBMATRICES.
C          E2(1) IS ALSO SET TO ZERO.
C
C        W CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES.
C
C        IND CONTAINS THE SUBMATRIX INDICES ASSOCIATED WITH THE
C          CORRESPONDING EIGENVALUES IN W -- 1 FOR EIGENVALUES
C          BELONGING TO THE FIRST SUBMATRIX FROM THE TOP,
C          2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC..
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C        RV1 IS A TEMPORARY STORAGE ARRAY.
C
C     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      K = 0
      TAG = 0
C
      DO 100 I = 1, N
         W(I) = D(I)
         IF (I .NE. 1) RV1(I-1) = E(I)
  100 CONTINUE
C
      E2(1) = 0.0D0
      RV1(N) = 0.0D0
C
      DO 290 L = 1, N
         J = 0
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
            IF (M .EQ. N) GO TO 120
            TST1 = DABS(W(M)) + DABS(W(M+1))
            TST2 = TST1 + DABS(RV1(M))
            IF (TST2 .EQ. TST1) GO TO 120
C     .......... GUARD AGAINST UNDERFLOWED ELEMENT OF E2 ..........
            IF (E2(M+1) .EQ. 0.0D0) GO TO 125
  110    CONTINUE
C
  120    IF (M .LE. K) GO TO 130
         IF (M .NE. N) E2(M+1) = 0.0D0
  125    K = M
         TAG = TAG + 1
  130    P = W(L)
         IF (M .EQ. L) GO TO 215
         IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         G = (W(L+1) - P) / (2.0D0 * RV1(L))
         R = PYTHAG(G,1.0D0)
         G = W(M) - P + RV1(L) / (G + DSIGN(R,G))
         S = 1.0D0
         C = 1.0D0
         P = 0.0D0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            F = S * RV1(I)
            B = C * RV1(I)
            R = PYTHAG(F,G)
            RV1(I+1) = R
            IF (R .EQ. 0.0D0) GO TO 210
            S = F / R
            C = G / R
            G = W(I+1) - P
            R = (W(I) - G) * S + 2.0D0 * C * B
            P = S * R
            W(I+1) = G + P
            G = C * R - B
  200    CONTINUE
C
         W(L) = W(L) - P
         RV1(L) = G
         RV1(M) = 0.0D0
         GO TO 105
C     .......... RECOVER FROM UNDERFLOW ..........
  210    W(I+1) = W(I+1) - P
         RV1(M) = 0.0D0
         GO TO 105
C     .......... ORDER EIGENVALUES ..........
  215    IF (L .EQ. 1) GO TO 250
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. W(I-1)) GO TO 270
            W(I) = W(I-1)
            IND(I) = IND(I-1)
  230    CONTINUE
C
  250    I = 1
  270    W(I) = P
         IND(I) = TAG
  290 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
      SUBROUTINE BINVIT(N,D,E,E2,M,W,IND,Z,IERR,RV1,RV2,RV3,RV4,RV6)
C*
C*
C*    AUTHORS -
C*       THIS ROUTINE IS A MODIFICATION OF EISPACK ROUTINE TINVIT
C*       WHICH IS A TRANSLATION OF THE INVERSE ITERATION TECHNIQUE
C*       IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON.
C*       HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
C*       THE EISPACK VERSION WAS MODIFIED BY C. MOLER AND D. SPANGLER
C*       (NRCC - 1 APR 1979).  THIS VERSION IS BY S. T. ELBERT
C*       (AMES LABORATORY - USDOE), WHO FURTHER STREAMLINED THE CODE.
C*
C*    PURPOSE -
C*       THIS ROUTINE FINDS THOSE EIGENVECTORS OF A TRIDIAGONAL
C*       SYMMETRIC MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES,
C*       USING INVERSE ITERATION.
C*
C*    METHOD -
C*       THE CALCULATIONS PROCEED AS FOLLOWS.  FIRST, THE  E2  ARRAY
C*       IS INSPECTED FOR THE PRESENCE OF A ZERO ELEMENT DEFINING A
C*       SUBMATRIX.  THE EIGENVALUES BELONGING TO THIS SUBMATRIX ARE
C*       IDENTIFIED BY THEIR COMMON SUBMATRIX INDEX IN  IND.
C*
C*       THE EIGENVECTORS OF THE SUBMATRIX ARE THEN COMPUTED BY
C*       INVERSE ITERATION.  FIRST, THE  LU  DECOMPOSITION OF THE SUB-
C*       MATRIX WITH AN APPROXIMATE EIGENVALUE SUBTRACTED FROM ITS
C*       DIAGONAL ELEMENTS IS ACHIEVED BY GAUSSIAN ELIMINATION USING
C*       PARTIAL PIVOTING.  THE MULTIPLIERS DEFINING THE LOWER
C*       TRIANGULAR MATRIX  L  ARE STORED IN THE TEMPORARY ARRAY  RV4
C*       AND THE UPPER TRIANGULAR MATRIX  U  IS STORED IN THE THREE
C*       TEMPORARY ARRAYS  RV1,  RV2, AND  RV3.  SAVING THESE
C*       QUANTITIES IN RV1,  RV2,  RV3, AND  RV4 AVOIDS REPEATING
C*       THE  LU  DECOMPOSITION IF FURTHER ITERATIONS ARE REQUIRED.
C*       AN APPROXIMATE VECTOR, STORED IN  RV6, IS COMPUTED STARTING
C*       FROM AN INITIAL VECTOR, AND THE NORM OF THE APPROXIMATE
C*       VECTOR IS COMPARED WITH A NORM OF THE SUBMATRIX TO DETERMINE
C*       WHETHER THE GROWTH IS SUFFICIENT TO ACCEPT IT AS AN
C*       EIGENVECTOR.  IF THIS VECTOR IS ACCEPTED, ITS EUCLIDIAN NORM
C*       IS MADE 1.  IF THE GROWTH IS NOT SUFFICIENT, THIS VECTOR IS
C*       USED AS THE INITIAL VECTOR IN COMPUTING THE NEXT APPROXIMATE
C*       VECTOR.  THIS ITERATION PROCESS IS REPEATED AT MOST  5  TIMES.
C*
C*       EIGENVECTORS COMPUTED IN THE ABOVE WAY CORRESPONDING TO WELL-
C*       SEPARATED EIGENVALUES OF THIS SUBMATRIX WILL BE ORTHOGONAL.
C*       HOWEVER, EIGENVECTORS CORRESPONDING TO CLOSE EIGENVALUES OF
C*       THIS SUBMATRIX MAY NOT BE SATISFACTORILY ORTHOGONAL.  HENCE,
C*       TO INSURE ORTHOGONAL EIGENVECTORS, EACH APPROXIMATE VECTOR IS
C*       MADE ORTHOGONAL TO THOSE PREVIOUSLY COMPUTED EIGENVECTORS
C*       WHOSE EIGENVALUES ARE CLOSE TO THE CURRENT EIGENVALUE, AS
C*       DETERMINED BY THE PARAMETER GRPTOL.
C*       IF THE ORTHOGONALIZATION PROCESS PRODUCES A ZERO VECTOR, A
C*       COLUMN OF THE IDENTITY MATRIX IS USED AS AN INITIAL VECTOR
C*       FOR THE NEXT ITERATION.
C*
C*       IDENTICAL EIGENVALUES ARE PERTURBED SLIGHTLY IN AN ATTEMPT TO
C*       OBTAIN INDEPENDENT EIGENVECTORS.  THESE PERTURBATIONS ARE NOT
C*       RECORDED IN THE EIGENVALUE ARRAY  W.
C*
C*       THE ABOVE STEPS ARE REPEATED ON EACH SUBMATRIX UNTIL ALL THE
C*       EIGENVECTORS ARE COMPUTED.
C*
C*
C*    ON ENTRY-
C*
C*       N      - INTEGER
C*                THE ORDER OF THE MATRIX.
C*
C*       D      - W.P. REAL (N)
C*                THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C*
C*       E      - W.P. REAL (N)
C*                CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C*                IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C*
C*       E2     - W.P. REAL (N)
C*          CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E,
C*          WITH ZEROS CORRESPONDING TO NEGLIGIBLE ELEMENTS OF E.
C*          E(I) IS CONSIDERED NEGLIGIBLE IF IT IS NOT LARGER THAN
C*          THE PRODUCT OF THE RELATIVE MACHINE PRECISION AND THE SUM
C*          OF THE MAGNITUDES OF D(I) AND D(I-1).  E2(1) MUST CONTAIN
C*          0.0 IF THE EIGENVALUES ARE IN ASCENDING ORDER, OR 2.0
C*          IF THE EIGENVALUES ARE IN DESCENDING ORDER.  IF  BISECT,
C*          TRIDIB, OR  EIMQLV  HAS BEEN USED TO FIND THE EIGENVALUES,
C*          THEIR OUTPUT E2 ARRAY IS EXACTLY WHAT IS EXPECTED HERE.
C*
C*       M      - INTEGER
C*                THE NUMBER OF SPECIFIED EIGENVALUES.
C*
C*       W      - W.P. REAL (M)
C*                THE M EIGENVALUES IN ASCENDING OR DESCENDING ORDER.
C*
C*       IND    - INTEGER (M)
C*           CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
C*           ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
C*           1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
C*           TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.
C*
C*    ON EXIT-
C*
C*       Z      - W.P. REAL (NM,M)
C*          CONTAINS THE ASSOCIATED SET OF ORTHONORMAL EIGENVECTORS.
C*          ANY VECTOR WHICH FAILS TO CONVERGE IS SET TO ZERO.
C*
C*       IERR   - INTEGER
C*                SET TO
C*                ZERO   FOR NORMAL RETURN,
C*                -R     IF THE EIGENVECTOR CORRESPONDING TO THE R-TH
C*                       EIGENVALUE FAILS TO CONVERGE IN 5 ITERATIONS.
C*
C*       RV1, RV2, RV3, RV4, AND RV6  -  W.P. REAL (N)
C*                TEMPORARY STORAGE ARRAYS.
C*
C*    COMPLEXITY -
C*       O(N*M)
C*
C*    NOTE -
C*    QUESTIONS AND COMMENTS ABOUT EISPACK SHOULD BE DIRECTED TO
C*    B. S. GARBOW, APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LAB.
C*
C***********************************************************************
C
C                    DECLARATIONS
C
      INTEGER IND(M)
      DOUBLE PRECISION D(N),E(N),E2(N),W(M),Z(N,M)
     *                ,RV1(N),RV2(N),RV3(N),RV4(N),RV6(N)
      DOUBLE PRECISION ANORM,ANORM2,EPMACH,EPS2,EPS3,EPS4,ORDER,DDOT
     *                ,TEST,U,UK,V,XU,RTPRE,RTCUR ,ABS,PARM1,PARM2
C
C-----------------------------------------------------------------------
C
      IERR = 0
      IF (M .LE. 0) RETURN
C
C     ********** EPMACH IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C
      EPMACH = 0.5D0**20
      DO 110 I=1,100
         EPMACH = EPMACH*0.5D0
         TEST = 1.D0 + EPMACH
         IF(TEST .EQ. 1.D0) GO TO 120
  110 CONTINUE
  120 CONTINUE
      EPMACH=EPMACH+EPMACH
      EPMACH=EPMACH+EPMACH
C
      ITAG = 0
      ORDER = 1.D0 - E2(1)
      IP = 1
      DO 930 IQ=1,N
         IF (IQ.NE.N .AND. E2(IQ+1) .NE. 0.D0) GO TO 930
C
C           ********** ESTABLISH AND PROCESS NEXT SUBMATRIX **********
C                           WITH INDEX RANGE  IP  TO  IQ
C
         ITAG = ITAG + 1
         IS = 0
         RTCUR=W(1)
C
         DO 920 IR = 1, M
            IF (IND(IR) .NE. ITAG) GO TO 920
               RTPRE = RTCUR
               RTCUR = W(IR)
               IPP1 = IP + 1
               IPM1 = IP - 1
               NSUB = IQ - IPM1
               IF (IS .NE. 0) GO TO 510
C
C              ********** CHECK FOR ISOLATED ROOT **********
C
                  XU = 1.D0
                  IF (IP .NE. IQ) GO TO 490
C
C                        SET ISOLATED ROOT TO UNIT VECTOR
C
                     DO 410 I=1,N
                        Z(I,IR)=0.D0
  410                CONTINUE
                     Z(IP,IR)=1.D0
                     GO TO 920
C
  490             CONTINUE
                  ANORM = DABS(D(IP))
C
                  DO 500 I = IPP1, IQ
                     ANORM = ANORM + DABS(D(I)) + DABS(E(I))
  500             CONTINUE
C
C                 ********** EPS2 IS THE CRITERION FOR GROUPING,
C                            EPS3 REPLACES ZERO PIVOTS AND EQUAL
C                            ROOTS ARE MODIFIED BY EPS3,
C                            EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW
C
                  UK = DBLE(NSUB)
                  ANORM2 = ANORM / UK
                  EPS2 = 1.D-5 * ANORM2
                  EPS3 = EPMACH * ANORM
                  EPS4 = UK * EPS3
                  UK = EPS4 / DSQRT(UK)
                  IS = IP
                  IGROUP = 0
                  GO TO 520
C
C              ********** LOOK FOR CLOSE OR COINCIDENT ROOTS **********
C
  510          CONTINUE
               IGROUP = IGROUP + 1
               IF (DABS(RTCUR-RTPRE) .GE. EPS2) IGROUP =0
               IF (ORDER * (RTCUR - RTPRE) .LE. 0.D0)
     *            RTCUR = RTPRE + ORDER * EPS3
C
C              ********** ELIMINATION WITH INTERCHANGES AND
C                         INITIALIZATION OF VECTOR **********
C
  520          CONTINUE
               V = 0.D0
               IQM1 = IQ - 1
               RV6(IP) = UK
               U = D(IP) - RTCUR
               IF(IP .EQ. IQ) GO TO 595
                  V = E(IPP1)
                  IF (IPP1 .GT. IQM1) GO TO 590
                     DO 580 I = IPP1,IQM1
                        RV6(I) = UK
                        IF (DABS(E(I)) .LT. DABS(U)) GO TO 540
C                         **** WARNING -- A DIVIDE CHECK MAY OCCUR HERE
C                          IF E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY
                           XU = U / E(I)
                           RV4(I) = XU
                           RV1(I-1) = E(I)
                           RV2(I-1) = D(I) - RTCUR
                           RV3(I-1) = E(I+1)
                           U = V - XU * RV2(I-1)
                           V = -XU * RV3(I-1)
                           GO TO 580
C
  540                   CONTINUE
                        XU = E(I) / U
                        RV4(I) = XU
                        RV1(I-1) = U
                        RV2(I-1) = V
                        RV3(I-1) = 0.D0
                        U = D(I) - RTCUR - XU * V
                        V = E(I+1)
  580                CONTINUE
  590             CONTINUE
                  RV6(IQ) = UK
                  RV3(IQM1) = 0.D0
                  IF (DABS(E(IQ)) .LT. DABS(U)) GO TO 592
                     XU = U / E(IQ)
                     RV4(IQ) = XU
                     RV1(IQM1) = E(IQ)
                     RV2(IQM1) = D(IQ) - RTCUR
                     U = V - XU * RV2(IQM1)
                     V = 0.D0
                     GO TO 595
C
  592             CONTINUE
                  XU = E(IQ) / U
                  RV4(IQ) = XU
                  RV1(IQM1) = U
                  RV2(IQM1) = V
                  U = D(IQ) - RTCUR - XU * V
  595          CONTINUE
C
               IF (DABS(U) .LT. EPS3) U = EPS3
               RV1(IQ) = U
               RV2(IQ) = 0.D0
               RV3(IQ) = 0.D0
C
C              ********** BACK SUBSTITUTION
C                         FOR I=IQ STEP -1 UNTIL IP DO -- **********
C
               DO 830 ITERS=1,5
                  RV6(IQ) = RV6(IQ)/RV1(IQ)
                  U = RV6(IQ)
                  RV6(IQM1) = (RV6(IQM1) - U * RV2(IQM1)) / RV1(IQM1)
                  V = U
                  U = RV6(IQM1)
                  ANORM = DABS(U) + DABS(V)
                  IF (IPP1 .GT. IQM1) GO TO 625
                    IQM2=IQM1-1
                    DO 620 I = IQM2,IP,-1
                      RV6(I) = (RV6(I) - U*RV2(I) - V*RV3(I)) / RV1(I)
                      V = U
                      U = RV6(I)
                      ANORM = ANORM + DABS(U)
  620               CONTINUE
  625             CONTINUE
                  IF (IGROUP .EQ. 0) GO TO 700
C
C                    ********** ORTHOGONALIZE WITH RESPECT TO PREVIOUS
C                               MEMBERS OF GROUP **********
C
                     J = IR
                     DO 680 JJ = 1,IGROUP
  630                   J = J - 1
                        IF (IND(J) .NE. ITAG) GO TO 630
                        xu = 0d0
                        nsubup = nsub+ip-1
                        do 640 i = ip, nsubup
                          xu = xu + rv6(i)*z(i,j)
  640                   continue
                        xu = -xu
                        if (xu .eq. 0d0) goto 680
                        do 650 i = ip, nsubup
                          rv6(i) = rv6(i) + xu*z(i,j)
  650                   continue
  680                CONTINUE
                     ANORM = 0.D0
                     DO 690 I = IP,IQ
                        ANORM = ANORM + DABS(RV6(I))
  690                CONTINUE
C
  700             CONTINUE
                  IF (ANORM .GE. 1.D0) GO TO 840
C
C                    ********** FORWARD SUBSTITUTION **********
C
                     IF (ANORM .NE. 0.D0) GO TO 740
                        RV6(IS) = EPS4
                        IS = IS + 1
                        IF (IS .GT. IQ) IS = IP
                        GO TO 780
C
  740                CONTINUE
                     XU = EPS4 / ANORM
C
                     DO 760 I =IP,IQ
                        RV6(I) = RV6(I) * XU
  760                CONTINUE
C
C                    ********** ELIMINATION OPERATIONS ON NEXT VECTOR
C                               ITERATE **********
C
  780                CONTINUE
                     DO 820 I = IPP1,IQ
                        U = RV6(I)
C
C                   ********** IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE
C                              WAS PERFORMED EARLIER IN THE
C                              TRIANGULARIZATION PROCESS **********
C
                        IF (RV1(I-1) .NE. E(I)) GO TO 800
                           U = RV6(I-1)
                           RV6(I-1) = RV6(I)
  800                   CONTINUE
                        RV6(I) = U - RV4(I) * RV6(I-1)
  820                CONTINUE
  830          CONTINUE
C
C                 ********** SET ERROR -- NON-CONVERGED EIGENVECTOR
C
               IERR = -IR
               DO 835 I=1,N
                  Z(I,IR)=0.D0
  835          CONTINUE
               GO TO 920
C
C              ********** NORMALIZE SO THAT SUM OF SQUARES IS
C                         1 AND EXPAND TO FULL ORDER **********
C
  840          CONTINUE
               U = 0.D0
               DO 860 I = IP,IQ
                  U = U + RV6(I)**2
  860          CONTINUE
               XU = 1.D0 / DSQRT(U)
C
               IF(IPM1.LT.1) GO TO 880
                  DO 870 I = 1, IPM1
                     Z(I,IR) = 0.D0
  870             CONTINUE
  880          CONTINUE
               DO 890 I = IP,IQ
                  Z(I,IR) = RV6(I) * XU
  890          CONTINUE
               IQP1 = IQ + 1
               IF(IQP1 .GT. N) GO TO 905
                  DO 900 I=IQP1,N
                     Z(I,IR) = 0.D0
  900             CONTINUE
  905          CONTINUE
  920    CONTINUE
         IP = IQ + 1
  930 CONTINUE
      RETURN
      END
      subroutine btribk(n,ar,ai,tau,m,zr,zi)
c
      integer i,j,k,l,m,n
      double precision ar(n,n),ai(n,n),tau(2,n),zr(n,m),zi(n,m)
      double precision h,s,si,smin
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure trbak1, num. math. 11, 181-195(1968)
c     by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine forms the eigenvectors of a complex hermitian
c     matrix by back transforming those of the corresponding
c     real symmetric tridiagonal matrix determined by  btridi.
c
c     on input
c
c        n is the order of the matrix.
c
c        ar and ai contain information about the unitary trans-
c          formations used in the reduction by  btridi  in their
c          full upper triangles except for the diagonal of ar.
c
c        tau contains further information about the transformations.
c
c        m is the number of eigenvectors to be back transformed.
c
c        zr contains the eigenvectors to be back transformed
c          in its first m columns.
c
c     on output
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the transformed eigenvectors
c          in their first m columns.
c
c     note that the last component of each returned vector
c     is real and that vector euclidean norms are preserved.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version is adapted from the eispack htridi. In contrast
c     only the upper triangle of a is used for input
c
c     -----------------------------------------------------------------
c
      if (m .eq. 0) go to 200
c     .......... transform the eigenvectors of the real symmetric
c                tridiagonal matrix to those of the hermitian
c                tridiagonal matrix. ..........

      do 50 k = 1, n
c
         do 50 j = 1, m
            zi(k,j) = -zr(k,j) * tau(2,k)
            zr(k,j) = zr(k,j) * tau(1,k)
   50 continue
c
      if (n .eq. 1) go to 200
c     .......... recover and apply the householder matrices ......... .
      do 140 i = 2, n
         l = i - 1
         h = ai(i,i)
         if (h .eq. 0.0d0) go to 140
c
         do 130 j = 1, m
            s = 0.0d0
            si = 0.0d0
c
            do 110 k = 1, l
               s = s + ar(k,i) * zr(k,j) - ai(k,i) * zi(k,j)
               si = si + ar(k,i) * zi(k,j) + ai(k,i) * zr(k,j)
  110       continue
c     .......... double divisions avoid possible underflow ..........
            s = (s / h) / h
            smin = -s
            si = - (si / h) / h
c
            do 120 k = 1, l
               zr(k,j) = zr(k,j) + smin * ar(k,i) + si * ai(k,i)
               zi(k,j) = zi(k,j) + si * ar(k,i) + s * ai(k,i)
  120       continue
c
  130    continue
c
  140 continue
c
  200 return
      end
C#elseC
C      subroutine htridi(nm,n,ar,ai,d,e,e2,tau)
Cc
C      integer i,j,k,l,n,ii,nm,jp1
C      double precision ar(nm,n),ai(nm,n),d(n),e(n),e2(n),tau(2,n)
C      double precision f,g,h,fi,gi,hh,si,scale,pythag
Cc
Cc     this subroutine is a translation of a complex analogue of
Cc     the algol procedure tred1, num. math. 11, 181-195(1968)
Cc     by martin, reinsch, and wilkinson.
Cc     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
Cc
Cc     this subroutine reduces a complex hermitian matrix
Cc     to a real symmetric tridiagonal matrix using
Cc     unitary similarity transformations.
Cc
Cc     on input
Cc
Cc        nm must be set to the row dimension of two-dimensional
Cc          array parameters as declared in the calling program
Cc          dimension statement.
Cc
Cc        n is the order of the matrix.
Cc
Cc        ar and ai contain the real and imaginary parts,
Cc          respectively, of the complex hermitian input matrix.
Cc          only the lower triangle of the matrix need be supplied.
Cc
Cc     on output
Cc
Cc        ar and ai contain information about the unitary trans-
Cc          formations used in the reduction in their full lower
Cc          triangles.  their strict upper triangles and the
Cc          diagonal of ar are unaltered.
Cc
Cc        d contains the diagonal elements of the the tridiagonal matrix
Cc
Cc        e contains the subdiagonal elements of the tridiagonal
Cc          matrix in its last n-1 positions.  e(1) is set to zero.
Cc
Cc        e2 contains the squares of the corresponding elements of e.
Cc          e2 may coincide with e if the squares are not needed.
Cc
Cc        tau contains further information about the transformations.
Cc
Cc     calls pythag for  dsqrt(a*a + b*b) .
Cc
Cc     questions and comments should be directed to burton s. garbow,
Cc     mathematics and computer science div, argonne national laboratory
Cc
Cc     this version dated august 1983.
Cc
Cc     -----------------------------------------------------------------
Cc
C      tau(1,n) = 1.0d0
C      tau(2,n) = 0.0d0
Cc
C      do 100 i = 1, n
C  100 d(i) = ar(i,i)
Cc     .......... for i=n step -1 until 1 do -- ..........
C      do 300 ii = 1, n
C         i = n + 1 - ii
C         l = i - 1
C         h = 0.0d0
C         scale = 0.0d0
C         if (l .lt. 1) go to 130
Cc     .......... scale row (algol tol then not needed) ..........
C         do 120 k = 1, l
C  120    scale = scale + dabs(ar(i,k)) + dabs(ai(i,k))
Cc
C         if (scale .ne. 0.0d0) go to 140
C         tau(1,l) = 1.0d0
C         tau(2,l) = 0.0d0
C  130    e(i) = 0.0d0
C         e2(i) = 0.0d0
C         go to 290
Cc
C  140    do 150 k = 1, l
C            ar(i,k) = ar(i,k) / scale
C            ai(i,k) = ai(i,k) / scale
C            h = h + ar(i,k) * ar(i,k) + ai(i,k) * ai(i,k)
C  150    continue
Cc
C         e2(i) = scale * scale * h
C         g = dsqrt(h)
C         e(i) = scale * g
C         f = pythag(ar(i,l),ai(i,l))
Cc     .......... form next diagonal element of matrix t ......... .
C         if (f .eq. 0.0d0) go to 160
C         tau(1,l) = (ai(i,l) * tau(2,i) - ar(i,l) * tau(1,i)) / f
C         si = (ar(i,l) * tau(2,i) + ai(i,l) * tau(1,i)) / f
C         h = h + f * g
C         g = 1.0d0 + g / f
C         ar(i,l) = g * ar(i,l)
C         ai(i,l) = g * ai(i,l)
C         if (l .eq. 1) go to 270
C         go to 170
C  160    tau(1,l) = -tau(1,i)
C         si = tau(2,i)
C         ar(i,l) = g
C  170    f = 0.0d0
Cc
C         do 240 j = 1, l
C            g = 0.0d0
C            gi = 0.0d0
Cc     .......... form element of a*u ..........
C            do 180 k = 1, j
C               g = g + ar(j,k) * ar(i,k) + ai(j,k) * ai(i,k)
C               gi = gi - ar(j,k) * ai(i,k) + ai(j,k) * ar(i,k)
C  180       continue
Cc
C            jp1 = j + 1
C            if (l .lt. jp1) go to 220
Cc
C            do 200 k = jp1, l
C               g = g + ar(k,j) * ar(i,k) - ai(k,j) * ai(i,k)
C               gi = gi - ar(k,j) * ai(i,k) - ai(k,j) * ar(i,k)
C  200       continue
Cc     .......... form element of p ..........
C  220       e(j) = g / h
C            tau(2,j) = gi / h
C            f = f + e(j) * ar(i,j) - tau(2,j) * ai(i,j)
C  240    continue
Cc
C         hh = f / (h + h)
Cc     .......... form reduced a ..........
C         do 260 j = 1, l
C            f = ar(i,j)
C            g = e(j) - hh * f
C            e(j) = g
C            fi = -ai(i,j)
C            gi = tau(2,j) - hh * fi
C            tau(2,j) = -gi
Cc
C            do 260 k = 1, j
C               ar(j,k) = ar(j,k) - f * e(k) - g * ar(i,k)
C     x                           + fi * tau(2,k) + gi * ai(i,k)
C               ai(j,k) = ai(j,k) - f * tau(2,k) - g * ai(i,k)
C     x                           - fi * e(k) - gi * ar(i,k)
C  260    continue
Cc
C  270    do 280 k = 1, l
C            ar(i,k) = scale * ar(i,k)
C            ai(i,k) = scale * ai(i,k)
C  280    continue
Cc
C         tau(2,l) = -si
C  290    hh = d(i)
C         d(i) = ar(i,i)
C         ar(i,i) = hh
C         ai(i,i) = scale * dsqrt(h)
C  300 continue
Cc
C      return
C      end
C      subroutine imtql2(nm,n,d,e,z,ierr)
Cc
C      integer i,j,k,l,m,n,ii,nm,mml,ierr
C      double precision d(n),e(n),z(nm,n)
C      double precision b,c,f,g,p,r,s,tst1,tst2,pythag
Cc
Cc     this subroutine is a translation of the algol procedure imtql2,
Cc     num. math. 12, 377-383(1968) by martin and wilkinson,
Cc     as modified in num. math. 15, 450(1970) by dubrulle.
Cc     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
Cc
Cc     this subroutine finds the eigenvalues and eigenvectors
Cc     of a symmetric tridiagonal matrix by the implicit ql method.
Cc     the eigenvectors of a full symmetric matrix can also
Cc     be found if  tred2  has been used to reduce this
Cc     full matrix to tridiagonal form.
Cc
Cc     on input
Cc
Cc        nm must be set to the row dimension of two-dimensional
Cc          array parameters as declared in the calling program
Cc          dimension statement.
Cc
Cc        n is the order of the matrix.
Cc
Cc        d contains the diagonal elements of the input matrix.
Cc
Cc        e contains the subdiagonal elements of the input matrix
Cc          in its last n-1 positions.  e(1) is arbitrary.
Cc
Cc        z contains the transformation matrix produced in the
Cc          reduction by  tred2, if performed.  if the eigenvectors
Cc          of the tridiagonal matrix are desired, z must contain
Cc          the identity matrix.
Cc
Cc      on output
Cc
Cc        d contains the eigenvalues in ascending order.  if an
Cc          error exit is made, the eigenvalues are correct but
Cc          unordered for indices 1,2,...,ierr-1.
Cc
Cc        e has been destroyed.
Cc
Cc        z contains orthonormal eigenvectors of the symmetric
Cc          tridiagonal (or full) matrix.  if an error exit is made,
Cc          z contains the eigenvectors associated with the stored
Cc          eigenvalues.
Cc
Cc        ierr is set to
Cc          zero       for normal return,
Cc          j          if the j-th eigenvalue has not been
Cc                     determined after 30 iterations.
Cc
Cc     calls pythag for  dsqrt(a*a + b*b) .
Cc
Cc     questions and comments should be directed to burton s. garbow,
Cc     mathematics and computer science div, argonne national laboratory
Cc
Cc     this version dated august 1983.
Cc
Cc     -----------------------------------------------------------------
Cc
C      ierr = 0
C      if (n .eq. 1) go to 1001
Cc
C      do 100 i = 2, n
C  100 e(i-1) = e(i)
Cc
C      e(n) = 0.0d0
Cc
C      do 240 l = 1, n
C         j = 0
Cc     .......... look for small sub-diagonal element ..........
C  105    do 110 m = l, n
C            if (m .eq. n) go to 120
C            tst1 = dabs(d(m)) + dabs(d(m+1))
C            tst2 = tst1 + dabs(e(m))
C            if (tst2 .eq. tst1) go to 120
C  110    continue
Cc
C  120    p = d(l)
C         if (m .eq. l) go to 240
C         if (j .eq. 30) go to 1000
C         j = j + 1
Cc     .......... form shift ..........
C         g = (d(l+1) - p) / (2.0d0 * e(l))
C         r = pythag(g,1.0d0)
C         g = d(m) - p + e(l) / (g + dsign(r,g))
C         s = 1.0d0
C         c = 1.0d0
C         p = 0.0d0
C         mml = m - l
Cc     .......... for i=m-1 step -1 until l do -- ..........
C         do 200 ii = 1, mml
C            i = m - ii
C            f = s * e(i)
C            b = c * e(i)
C            r = pythag(f,g)
C            e(i+1) = r
C            if (r .eq. 0.0d0) go to 210
C            s = f / r
C            c = g / r
C            g = d(i+1) - p
C            r = (d(i) - g) * s + 2.0d0 * c * b
C            p = s * r
C            d(i+1) = g + p
C            g = c * r - b
Cc     .......... form vector ..........
C            do 180 k = 1, n
C               f = z(k,i+1)
C               z(k,i+1) = s * z(k,i) + c * f
C               z(k,i) = c * z(k,i) - s * f
C  180       continue
Cc
C  200    continue
Cc
C         d(l) = d(l) - p
C         e(l) = g
C         e(m) = 0.0d0
C         go to 105
Cc     .......... recover from underflow ..........
C  210    d(i+1) = d(i+1) - p
C         e(m) = 0.0d0
C         go to 105
C  240 continue
Cc     .......... order eigenvalues and eigenvectors ..........
C      do 300 ii = 2, n
C         i = ii - 1
C         k = i
C         p = d(i)
Cc
C         do 260 j = ii, n
C            if (d(j) .ge. p) go to 260
C            k = j
C            p = d(j)
C  260    continue
Cc
C         if (k .eq. i) go to 300
C         d(k) = d(i)
C         d(i) = p
Cc
C         do 280 j = 1, n
C            p = z(j,i)
C            z(j,i) = z(j,k)
C            z(j,k) = p
C  280    continue
Cc
C  300 continue
Cc
C      go to 1001
Cc     .......... set error -- no convergence to an
Cc                eigenvalue after 30 iterations ..........
C 1000 ierr = l
C 1001 return
C      end
C      subroutine htribk(nm,n,ar,ai,tau,m,zr,zi)
Cc
C      integer i,j,k,l,m,n,nm
C      double precision ar(nm,n),ai(nm,n),tau(2,n),zr(nm,m),zi(nm,m)
C      double precision h,s,si
Cc
Cc     this subroutine is a translation of a complex analogue of
Cc     the algol procedure trbak1, num. math. 11, 181-195(1968)
Cc     by martin, reinsch, and wilkinson.
Cc     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
Cc
Cc     this subroutine forms the eigenvectors of a complex hermitian
Cc     matrix by back transforming those of the corresponding
Cc     real symmetric tridiagonal matrix determined by  htridi.
Cc
Cc     on input
Cc
Cc        nm must be set to the row dimension of two-dimensional
Cc          array parameters as declared in the calling program
Cc          dimension statement.
Cc
Cc        n is the order of the matrix.
Cc
Cc        ar and ai contain information about the unitary trans-
Cc          formations used in the reduction by  htridi  in their
Cc          full lower triangles except for the diagonal of ar.
Cc
Cc        tau contains further information about the transformations.
Cc
Cc        m is the number of eigenvectors to be back transformed.
Cc
Cc        zr contains the eigenvectors to be back transformed
Cc          in its first m columns.
Cc
Cc     on output
Cc
Cc        zr and zi contain the real and imaginary parts,
Cc          respectively, of the transformed eigenvectors
Cc          in their first m columns.
Cc
Cc     note that the last component of each returned vector
Cc     is real and that vector euclidean norms are preserved.
Cc
Cc     questions and comments should be directed to burton s. garbow,
Cc     mathematics and computer science div, argonne national laboratory
Cc
Cc     this version dated august 1983.
Cc
Cc     -----------------------------------------------------------------
Cc
C      if (m .eq. 0) go to 200
Cc     .......... transform the eigenvectors of the real symmetric
Cc                tridiagonal matrix to those of the hermitian
Cc                tridiagonal matrix. ..........
C      do 50 k = 1, n
Cc
C         do 50 j = 1, m
C            zi(k,j) = -zr(k,j) * tau(2,k)
C            zr(k,j) = zr(k,j) * tau(1,k)
C   50 continue
Cc
C      if (n .eq. 1) go to 200
Cc     .......... recover and apply the householder matrices ......... .
C      do 140 i = 2, n
C         l = i - 1
C         h = ai(i,i)
C         if (h .eq. 0.0d0) go to 140
Cc
C         do 130 j = 1, m
C            s = 0.0d0
C            si = 0.0d0
Cc
C            do 110 k = 1, l
C               s = s + ar(i,k) * zr(k,j) - ai(i,k) * zi(k,j)
C               si = si + ar(i,k) * zi(k,j) + ai(i,k) * zr(k,j)
C  110       continue
Cc     .......... double divisions avoid possible underflow ..........
C            s = (s / h) / h
C            si = (si / h) / h
Cc
C            do 120 k = 1, l
C               zr(k,j) = zr(k,j) - s * ar(i,k) - si * ai(i,k)
C               zi(k,j) = zi(k,j) - si * ar(i,k) + s * ai(i,k)
C  120       continue
Cc
C  130    continue
Cc
C  140 continue
Cc
C  200 return
C      end
C#endif
      double precision function pythag(a,b)
C
C     Finds dsqrt(a**2+b**2) without overflow or destructive underflow
C
      double precision a,b
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
C#endif
