      subroutine yygefa(ar,ai,lda,n,ipvt,info)
C Adapted from linpack, but uses no complex arithmetic
C     implicit none
      integer lda,n,ipvt(1),info
      double precision ar(lda,1),ai(lda,1)
c
c     yygefa factors a complex*16 matrix by gaussian elimination.
c
c     yygefa is an adaptation of zgefa, using real arithmetic
c     on entry
c
c        ar,ai   double precision (lda, n)
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
c                     indicate that yygedi will divide by zero
c                     if called.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas-adapted yyaxpy,yscal,iyamax
c     fortran dabs
c
c     internal variables
c
      integer iyamax,j,k,kp1,l,nm1
#if DCMPLX
      complex*16 t,zdum,zdumr,zdumi
      double precision cabs1,dreal,dimag
      dreal(zdumr) = zdumr
      dimag(zdumi) = (zero,-one)*zdumi
      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
#else
      double precision tmp,t(2)
      double precision zero,one
      parameter (one=1d0, zero=0d0)
#endif
c
c     gaussian elimination with partial pivoting
c
C     call tcn('yygefa')
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = iyamax(n-k+1,ar(k,k),ai(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (dabs(ar(l,k))+dabs(ai(l,k)) .eq. zero) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               tmp = ar(l,k)
               ar(l,k) = ar(k,k)
               ar(k,k) = tmp
               tmp = ai(l,k)
               ai(l,k) = ai(k,k)
               ai(k,k) = tmp
   10       continue
c
c           compute multipliers
c
            call cdiv(-one,zero,ar(k,k),ai(k,k),t(1),t(2))
            call yscal(n-k,t(1),t(2),ar(k+1,k),ai(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t(1) = ar(l,j)
               t(2) = ai(l,j)
               if (l .eq. k) go to 20
                  ar(l,j) = ar(k,j)
                  ai(l,j) = ai(k,j)
                  ar(k,j) = t(1)
                  ai(k,j) = t(2)
   20          continue
               call yyaxpy(n-k,t(1),t(2),ar(k+1,k),ai(k+1,k),
     .           1,ar(k+1,j),ai(k+1,j),1,.true.)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (dabs(ar(n,n))+dabs(ai(n,n)) .eq. zero) info = n
C     call tcx('yygefa')
      end
      subroutine yygedi(ar,ai,lda,n,ipvt,det,wkr,wki,job)
C     implicit none
      integer lda,n,ipvt(1),job
Cr  This version is a real adaptation of zgedi
      double precision ar(lda,1),ai(lda,1),wkr(1),wki(1),det(2,2)
c
c     yygedi computes the determinant and inverse of a matrix
c     using the factors computed by yygefa.
c
c     on entry
c
c        ar,ai   double precision (lda, n)
c                the output from yygefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from yygefa.
c
c        wkr,wki double precision(n)
c                work vectors.  contents destroyed.
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
c        and if zgeco has set rcond .gt. 0.0 or yygefa has set
c        info .eq. 0 .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas-adapted yyaxpy,yscal,yswap
c
c     internal variables
c
      integer i,j,k,kb,kp1,l,nm1
      double precision zero,one,ten
      parameter (one=1d0, zero=0d0, ten=10d0)

      double precision t(2)
c
c     compute determinant
c
C     call tcn('yygedi')
      if (job/10 .eq. 0) go to 70
         det(1,1) = 1
         det(2,1) = 0
         det(1,2) = 0
         det(2,2) = 0
         do 50 i = 1, n
            if (ipvt(i) .ne. i) det(1,1) = -det(1,1)
            if (ipvt(i) .ne. i) det(2,1) = -det(2,1)
            call cpy(ar(i,i),ai(i,i),det,det(2,1),det,det(2,1))
c        ...exit
            if (dabs(det(1,1))+dabs(det(2,1)) .eq. zero) go to 60
   10       if (dabs(det(1,1))+dabs(det(2,1)) .ge. one) go to 20
c        ...   det(1) = dcmplx(ten,zero)*det(1)
               call cpy(ten,zero,det(1,1),det(2,1),det(1,1),det(2,1))
               det(1,2) = det(1,2) - 1
            go to 10
   20       continue
   30       if (dabs(det(1,1))+dabs(det(2,1)) .lt. ten) go to 40
c        ... complex divide ... det(1) = det(1)/dcmplx(ten,zero)
               call cdiv(det(1,1),det(2,1),ten,zero,det(1,1),det(2,1))
               det(1,2) = det(1,2) + 1
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
c
c     compute inverse(u)
c
      if (mod(job,10) .eq. 0) go to 150
         do 100 k = 1, n
            call cdiv(one,zero,ar(k,k),ai(k,k),ar(k,k),ai(k,k))
            t(1) = -ar(k,k)
            t(2) = -ai(k,k)
            call yscal(k-1,t(1),t(2),ar(1,k),ai(1,k),1)
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t(1) = ar(k,j)
               t(2) = ai(k,j)
               ar(k,j) = 0
               ai(k,j) = 0
               call yyaxpy(k,t(1),t(2),ar(1,k),ai(1,k),1,
     .           ar(1,j),ai(1,j),1,.true.)
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
               wkr(i) = ar(i,k)
               wki(i) = ai(i,k)
               ar(i,k) = 0
               ai(i,k) = 0
  110       continue
            do 120 j = kp1, n
               t(1) = wkr(j)
               t(2) = wki(j)
               call yyaxpy(n,t(1),t(2),ar(1,j),ai(1,j),1,
     .           ar(1,k),ai(1,k),1,.true.)
  120       continue
            l = ipvt(k)
            if (l .ne. k) then
              call dswap(n,ar(1,k),1,ar(1,l),1)
              call dswap(n,ai(1,k),1,ai(1,l),1)
            endif
  130    continue
  140    continue
  150 continue
C     call tcx('yygedi')
      end

