      subroutine dspsl(ap,n,kpvt,b)
      integer n,kpvt(1)
      double precision ap(1),b(1)
c
c     dsisl solves the double precision symmetric system
c     a * x = b
c     using the factors computed by dspfa.
c
c     on entry
c
c        ap      double precision(n*(n+1)/2)
c                the output from dspfa.
c
c        n       integer
c                the order of the matrix  a .
c
c        kpvt    integer(n)
c                the pivot vector from dspfa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero may occur if  dspco  has set rcond .eq. 0.0
c        or  dspfa  has set info .ne. 0  .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dspfa(ap,n,kpvt,info)
c           if (info .ne. 0) go to ...
c           do 10 j = 1, p
c              call dspsl(ap,n,kpvt,c(1,j))
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c     fortran iabs
c
c     internal variables.
c
      double precision ak,akm1,bk,bkm1,ddot,denom,temp
      integer ik,ikm1,ikp1,k,kk,km1k,km1km1,kp
c
c     loop backward applying the transformations and
c     d inverse to b.
c
#if CRAY
      call sspsl(ap,n,kpvt,b)
#else
      k = n
      ik = (n*(n - 1))/2
   10 if (k .eq. 0) go to 80
      kk = ik + k
      if (kpvt(k) .lt. 0) go to 40
c
c           1 x 1 pivot block.
c
      if (k .eq. 1) go to 30
      kp = kpvt(k)
      if (kp .eq. k) go to 20
c
c                 interchange.
c
      temp = b(k)
      b(k) = b(kp)
      b(kp) = temp
   20 continue
c
c              apply the transformation.
c
      call daxpy(k-1,b(k),ap(ik+1),1,b(1),1)
   30 continue
c
c           apply d inverse.
c
      b(k) = b(k)/ap(kk)
      k = k - 1
      ik = ik - k
      go to 70
   40 continue
c
c           2 x 2 pivot block.
c
      ikm1 = ik - (k - 1)
      if (k .eq. 2) go to 60
      kp = iabs(kpvt(k))
      if (kp .eq. k - 1) go to 50
c
c                 interchange.
c
      temp = b(k-1)
      b(k-1) = b(kp)
      b(kp) = temp
   50 continue
c
c              apply the transformation.
c
      call daxpy(k-2,b(k),ap(ik+1),1,b(1),1)
      call daxpy(k-2,b(k-1),ap(ikm1+1),1,b(1),1)
   60 continue
c
c           apply d inverse.
c
      km1k = ik + k - 1
      kk = ik + k
      ak = ap(kk)/ap(km1k)
      km1km1 = ikm1 + k - 1
      akm1 = ap(km1km1)/ap(km1k)
      bk = b(k)/ap(km1k)
      bkm1 = b(k-1)/ap(km1k)
      denom = ak*akm1 - 1.0d0
      b(k) = (akm1*bk - bkm1)/denom
      b(k-1) = (ak*bkm1 - bk)/denom
      k = k - 2
      ik = ik - (k + 1) - k
   70 continue
      go to 10
   80 continue
c
c     loop forward applying the transformations.
c
      k = 1
      ik = 0
   90 if (k .gt. n) go to 160
      if (kpvt(k) .lt. 0) go to 120
c
c           1 x 1 pivot block.
c
      if (k .eq. 1) go to 110
c
c              apply the transformation.
c
      b(k) = b(k) + ddot(k-1,ap(ik+1),1,b(1),1)
      kp = kpvt(k)
      if (kp .eq. k) go to 100
c
c                 interchange.
c
      temp = b(k)
      b(k) = b(kp)
      b(kp) = temp
  100 continue
  110 continue
      ik = ik + k
      k = k + 1
      go to 150
  120 continue
c
c           2 x 2 pivot block.
c
      if (k .eq. 1) go to 140
c
c              apply the transformation.
c
      b(k) = b(k) + ddot(k-1,ap(ik+1),1,b(1),1)
      ikp1 = ik + k
      b(k+1) = b(k+1) + ddot(k-1,ap(ikp1+1),1,b(1),1)
      kp = iabs(kpvt(k))
      if (kp .eq. k) go to 130
c
c                 interchange.
c
      temp = b(k)
      b(k) = b(kp)
      b(kp) = temp
  130 continue
  140 continue
      ik = ik + k + k + 1
      k = k + 2
  150 continue
      go to 90
  160 continue
      return
#endif
      end

