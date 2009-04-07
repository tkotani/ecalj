      subroutine zgefa(a,lda,n,ipvt,info)
C Adapted from linpack, but uses no complex arithmetic
C     implicit none
      integer lda,n,ipvt(1),info
#if DCMPLX
c      complex*16 a(lda,1)
#else
      double precision a(2,lda,1)
#endif
c
c     zgefa factors a complex*16 matrix by gaussian elimination.
c
c     zgefa is usually called by zgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for zgeco) = (1 + 9/n)*(time for zgefa) .
c
c     on entry
c
c        a       complex*16(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that zgesl or zgedi will divide by zero
c                     if called.  use  rcond  in zgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas zaxpy,zscal,izamax
c     fortran dabs
c
c     internal variables
c
      integer izamax,j,k,kp1,l,nm1
      double precision dcabs1
#if DCMPLX
      complex*16 t,zdum,zdumr,zdumi
      double precision cabs1,dreal,dimag
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
#else
      double precision tmp,t(2)
#endif
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
#if DCMPLX
         l = izamax(n-k+1,a(k,k),1) + k - 1
#else
         l = izamax(n-k+1,a(1,k,k),1) + k - 1
#endif
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
#if DCMPLX
         if (cabs1(a(l,k)) .eq. 0.0d0) go to 40
#else
         if (dcabs1(a(1,l,k)) .eq. 0.0d0) go to 40
#endif
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
#if DCMPLX
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
#else
               tmp = a(1,l,k)
               a(1,l,k) = a(1,k,k)
               a(1,k,k) = tmp
               tmp = a(2,l,k)
               a(2,l,k) = a(2,k,k)
               a(2,k,k) = tmp
#endif
   10       continue
c
c           compute multipliers
c
#if DCMPLX
            t = -(1.0d0,0.0d0)/a(k,k)
            call zscal(n-k,t,a(k+1,k),1)
#else
            call cdiv(-1d0,0d0,a(1,k,k),a(2,k,k),t(1),t(2))
            call zscal(n-k,t,a(1,k+1,k),1)
#endif
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
#if DCMPLX
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call zaxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
#else
               t(1) = a(1,l,j)
               t(2) = a(2,l,j)
               if (l .eq. k) go to 20
                  a(1,l,j) = a(1,k,j)
                  a(2,l,j) = a(2,k,j)
                  a(1,k,j) = t(1)
                  a(2,k,j) = t(2)
   20          continue
               call zaxpy(n-k,t,a(1,k+1,k),1,a(1,k+1,j),1)
#endif
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
#if DCMPLX
      if (cabs1(a(n,n)) .eq. 0.0d0) info = n
#else
      if (dcabs1(a(1,n,n)) .eq. 0.0d0) info = n
#endif
      end
      subroutine zgedi(a,lda,n,ipvt,det,work,job)
C     implicit none
      integer lda,n,ipvt(1),job
Cr  This version uses no internal complex arithmetic.
Cr  Note: evaluation of determinant not carefully tested...
#if DCMPLX
      complex*16 a(lda,1),det(2),work(1)
#else
      double precision a(2,lda,1),work(2,1),det(2,2)
#endif
c
c     zgedi computes the determinant and inverse of a matrix
c     using the factors computed by zgeco or zgefa.
c
c     on entry
c
c        a       complex*16(lda, n)
c                the output from zgeco or zgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from zgeco or zgefa.
c
c        work    complex*16(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        a       inverse of original matrix if requested.
c                otherwise unchanged.
c
c        det     complex*16(2)
c                determinant of original matrix if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. dcabs1<(det(1)) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if zgeco has set rcond .gt. 0.0 or zgefa has set
c        info .eq. 0 .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas zaxpy,zscal,zswap
c     fortran dabs,dcmplx,mod
c
c     internal variables
c
      integer i,j,k,kb,kp1,l,nm1
      double precision ten
#if DCMPLX
      complex*16 t,zdum,zdumr,zdumi
      double precision cabs1,dreal,dimag
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
#else
      double precision t(2),dcabs1
#endif
c
c     compute determinant
c
      if (job/10 .eq. 0) go to 70
#if DCMPLX
         det(1) = (1.0d0,0.0d0)
         det(2) = (0.0d0,0.0d0)
#else
         det(1,1) = 1
         det(2,1) = 0
         det(1,2) = 0
         det(2,2) = 0
#endif
         ten = 10.0d0
         do 50 i = 1, n
#if DCMPLX
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
c        ...exit
            if (cabs1(det(1)) .eq. 0.0d0) go to 60
   10       if (cabs1(det(1)) .ge. 1.0d0) go to 20
               det(1) = dcmplx(ten,0.0d0)*det(1)
               det(2) = det(2) - (1.0d0,0.0d0)
            go to 10
   20       continue
   30       if (cabs1(det(1)) .lt. ten) go to 40
               det(1) = det(1)/dcmplx(ten,0.0d0)
               det(2) = det(2) + (1.0d0,0.0d0)
            go to 30
#else
            if (ipvt(i) .ne. i) det(1,1) = -det(1,1)
            if (ipvt(i) .ne. i) det(2,1) = -det(2,1)
            call cpy(a(1,i,i),a(2,i,i),det,det(2,1),det,det(2,1))
c        ...exit
            if (dcabs1(det) .eq. 0.0d0) go to 60
   10       if (dcabs1(det) .ge. 1.0d0) go to 20
c               det(1) = dcmplx(ten,0.0d0)*det(1)
               call cpy(ten,0d0,det(1,1),det(2,1),det(1,1),det(2,1))
               det(1,2) = det(1,2) - 1
            go to 10
   20       continue
   30       if (dcabs1(det) .lt. ten) go to 40
c complex divide ... det(1) = det(1)/dcmplx(ten,0.0d0) ...
               call cdiv(det(1,1),det(2,1),ten,0d0,det(1,1),det(2,1))
               det(1,2) = det(1,2) + 1
            go to 30
#endif
   40       continue
   50    continue
   60    continue
   70 continue
c
c     compute inverse(u)
c
      if (mod(job,10) .eq. 0) go to 150
         do 100 k = 1, n
#if DCMPLX
            a(k,k) = (1.0d0,0.0d0)/a(k,k)
            t = -a(k,k)
            call zscal(k-1,t,a(1,k),1)
#else
            call cdiv(1d0,0d0,a(1,k,k),a(2,k,k),a(1,k,k),a(2,k,k))
            t(1) = -a(1,k,k)
            t(2) = -a(2,k,k)
            call zscal(k-1,t,a(1,1,k),1)
#endif
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
#if DCMPLX
               t = a(k,j)
               a(k,j) = (0.0d0,0.0d0)
               call zaxpy(k,t,a(1,k),1,a(1,j),1)
#else
               t(1) = a(1,k,j)
               t(2) = a(2,k,j)
               a(1,k,j) = 0
               a(2,k,j) = 0
               call zaxpy(k,t,a(1,1,k),1,a(1,1,j),1)
#endif
   80       continue
   90       continue
  100    continue
c
c        form inverse(u)*inverse(l)
c
         nm1 = n - 1
         if (nm1 .lt. 1) go to 140
         do 130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do 110 i = kp1, n
#if DCMPLX
               work(i) = a(i,k)
               a(i,k) = (0.0d0,0.0d0)
#else
               work(1,i) = a(1,i,k)
               work(2,i) = a(2,i,k)
               a(1,i,k) = 0
               a(2,i,k) = 0
#endif
  110       continue
            do 120 j = kp1, n
#if DCMPLX
               t = work(j)
               call zaxpy(n,t,a(1,j),1,a(1,k),1)
#else
               t(1) = work(1,j)
               t(2) = work(2,j)
               call zaxpy(n,t,a(1,1,j),1,a(1,1,k),1)
#endif
  120       continue
            l = ipvt(k)
#if DCMPLX
c            if (l .ne. k) call zswap(n,a(1,k),1,a(1,l),1)
#else
            if (l .ne. k) call zswap(n,a(1,1,k),1,a(1,1,l),1)
#endif
  130    continue
  140    continue
  150 continue
      return
      end

