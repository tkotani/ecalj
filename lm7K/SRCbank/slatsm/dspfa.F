      subroutine dspfa(ap,n,kpvt,info)
      integer n,kpvt(1),info
      double precision ap(1)
c
c     dspfa factors a double precision symmetric matrix stored in
c     packed form by elimination with symmetric pivoting.
c
c     to solve  a*x = b , follow dspfa by dspsl.
c     to compute  inverse(a)*c , follow dspfa by dspsl.
c     to compute  determinant(a) , follow dspfa by dspdi.
c     to compute  inertia(a) , follow dspfa by dspdi.
c     to compute  inverse(a) , follow dspfa by dspdi.
c
c     on entry
c
c        ap      double precision (n*(n+1)/2)
c                the packed form of a symmetric matrix  a .  the
c                columns of the upper triangle are stored sequentially
c                in a one-dimensional array of length  n*(n+1)/2 .
c                see comments below for details.
c
c        n       integer
c                the order of the matrix  a .
c
c     output
c
c        ap      a block diagonal matrix and the multipliers which
c                were used to obtain it stored in packed form.
c                the factorization can be written  a = u*d*trans(u)
c                where  u  is a product of permutation and unit
c                upper triangular matrices , trans(u) is the
c                transpose of  u , and  d  is block diagonal
c                with 1 by 1 and 2 by 2 blocks.
c
c        kpvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if the k-th pivot block is singular. this is
c                     not an error condition for this subroutine,
c                     but it does indicate that dspsl or dspdi may
c                     divide by zero if called.
c
c     packed storage
c
c          the following program segment will pack the upper
c          triangle of a symmetric matrix.
c
c                k = 0
c                do 20 j = 1, n
c                   do 10 i = 1, j
c                      k = k + 1
c                      ap(k)  = a(i,j)
c             10    continue
c             20 continue
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas daxpy,dswap,idamax
c     fortran dabs,dmax1,dsqrt
c
c     internal variables
c
      double precision ak,akm1,bk,bkm1,denom,mulk,mulkm1,t
      double precision absakk,alpha,colmax,rowmax
      integer idamax,ij,ik,ikm1,im,imax,imaxp1,imim,imj,imk
      integer j,jj,jk,jkm1,jmax,jmim,k,kk,km1,km1k,km1km1,km2,kstep
      logical swap
c
c
c     initialize
c
c     alpha is used in choosing pivot block size.
#if CRAY
      call sspfa(ap,n,kpvt,info)
#else
      alpha = (1.0d0 + dsqrt(17.0d0))/8.0d0
c
      info = 0
c
c     main loop on k, which goes from n to 1.
c
      k = n
      ik = (n*(n - 1))/2
   10 continue
c
c        leave the loop if k=0 or k=1.
c
c     ...exit
      if (k .eq. 0) go to 200
      if (k .gt. 1) go to 20
      kpvt(1) = 1
      if (ap(1) .eq. 0.0d0) info = 1
c     ......exit
      go to 200
   20 continue
c
c        this section of code determines the kind of
c        elimination to be performed.  when it is completed,
c        kstep will be set to the size of the pivot block, and
c        swap will be set to .true. if an interchange is
c        required.
c
      km1 = k - 1
      kk = ik + k
      absakk = dabs(ap(kk))
c
c        determine the largest off-diagonal element in
c        column k.
c
      imax = idamax(k-1,ap(ik+1),1)
      imk = ik + imax
      colmax = dabs(ap(imk))
      if (absakk .lt. alpha*colmax) go to 30
      kstep = 1
      swap = .false.
      go to 90
   30 continue
c
c           determine the largest off-diagonal element in
c           row imax.
c
      rowmax = 0.0d0
      imaxp1 = imax + 1
      im = imax*(imax - 1)/2
      imj = im + 2*imax
      do 40 j = imaxp1, k
        rowmax = dmax1(rowmax,dabs(ap(imj)))
        imj = imj + j
   40 continue
      if (imax .eq. 1) go to 50
      jmax = idamax(imax-1,ap(im+1),1)
      jmim = jmax + im
      rowmax = dmax1(rowmax,dabs(ap(jmim)))
   50 continue
      imim = imax + im
      if (dabs(ap(imim)) .lt. alpha*rowmax) go to 60
      kstep = 1
      swap = .true.
      go to 80
   60 continue
      if (absakk .lt. alpha*colmax*(colmax/rowmax)) go to 70
      kstep = 1
      swap = .false.
      go to 80
   70 continue
      kstep = 2
      swap = imax .ne. km1
   80 continue
   90 continue
      if (dmax1(absakk,colmax) .ne. 0.0d0) go to 100
c
c           column k is zero.  set info and iterate the loop.
c
      kpvt(k) = k
      info = k
      go to 190
  100 continue
      if (kstep .eq. 2) go to 140
c
c           1 x 1 pivot block.
c
      if (.not.swap) go to 120
c
c              perform an interchange.
c
      call dswap(imax,ap(im+1),1,ap(ik+1),1)
      imj = ik + imax
      do 110 jj = imax, k
        j = k + imax - jj
        jk = ik + j
        t = ap(jk)
        ap(jk) = ap(imj)
        ap(imj) = t
        imj = imj - (j - 1)
  110 continue
  120 continue
c
c           perform the elimination.
c
      ij = ik - (k - 1)
      do 130 jj = 1, km1
        j = k - jj
        jk = ik + j
        mulk = -ap(jk)/ap(kk)
        t = mulk
        call daxpy(j,t,ap(ik+1),1,ap(ij+1),1)
        ap(jk) = mulk
        ij = ij - (j - 1)
  130 continue
c
c           set the pivot array.
c
      kpvt(k) = k
      if (swap) kpvt(k) = imax
      go to 190
  140 continue
c
c           2 x 2 pivot block.
c
      km1k = ik + k - 1
      ikm1 = ik - (k - 1)
      if (.not.swap) go to 160
c
c              perform an interchange.
c
      call dswap(imax,ap(im+1),1,ap(ikm1+1),1)
      imj = ikm1 + imax
      do 150 jj = imax, km1
        j = km1 + imax - jj
        jkm1 = ikm1 + j
        t = ap(jkm1)
        ap(jkm1) = ap(imj)
        ap(imj) = t
        imj = imj - (j - 1)
  150 continue
      t = ap(km1k)
      ap(km1k) = ap(imk)
      ap(imk) = t
  160 continue
c
c           perform the elimination.
c
      km2 = k - 2
      if (km2 .eq. 0) go to 180
      ak = ap(kk)/ap(km1k)
      km1km1 = ikm1 + k - 1
      akm1 = ap(km1km1)/ap(km1k)
      denom = 1.0d0 - ak*akm1
      ij = ik - (k - 1) - (k - 2)
      do 170 jj = 1, km2
        j = km1 - jj
        jk = ik + j
        bk = ap(jk)/ap(km1k)
        jkm1 = ikm1 + j
        bkm1 = ap(jkm1)/ap(km1k)
        mulk = (akm1*bk - bkm1)/denom
        mulkm1 = (ak*bkm1 - bk)/denom
        t = mulk
        call daxpy(j,t,ap(ik+1),1,ap(ij+1),1)
        t = mulkm1
        call daxpy(j,t,ap(ikm1+1),1,ap(ij+1),1)
        ap(jk) = mulk
        ap(jkm1) = mulkm1
        ij = ij - (j - 1)
  170 continue
  180 continue
c
c           set the pivot array.
c
      kpvt(k) = 1 - k
      if (swap) kpvt(k) = -imax
      kpvt(k-1) = kpvt(k)
  190 continue
      ik = ik - (k - 1)
      if (kstep .eq. 2) ik = ik - (k - 2)
      k = k - kstep
      go to 10
  200 continue
      return
#endif
      end

