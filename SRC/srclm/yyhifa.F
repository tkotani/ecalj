      subroutine yyhifa(ar,ai,lda,n,kpvt,info)
C- Factor an hermitian matrix with symmetric pivoting
      integer lda,n,kpvt(1),info
      double precision ar(lda,1),ai(lda,1)
c
c     to solve  a*x = b , follow yyhifa by yhisl.
c     to compute  inverse(a)*c , follow yyhifa by yhisl.
c     to compute  determinant(a) , follow yyhifa by yhidi.
c     to compute  inertia(a) , follow yyhifa by yhidi.
c     to compute  inverse(a) , follow yyhifa by yhidi.
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
     .bk(2),bkm1(2),t(2)

      double precision absakk,alpha,colmax,rowmax,dcabsr
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
      if (dcabsr(ar,ai) .eq. 0.0d0) info = 1
      go to 200
   20 continue
c
c     This section of code determines the kind of
c     elimination to be performed.  when it is completed,
c     kstep will be set to the size of the pivot block, and
c     swap will be set to .true. if an interchange is
c     required.

      km1 = k - 1
      absakk = dcabsr(ar(k,k),ai(k,k))

C ---  Determine the largest off-diagonal element in column k ---
      imax = iyamax(k-1,ar(1,k),ai(1,k),1)
      colmax = dcabsr(ar(imax,k),ai(imax,k))
      if (absakk .lt. alpha*colmax) go to 30
      kstep = 1
      swap = .false.
      go to 90
   30 continue

C ---  Determine the largest off-diagonal element in row imax. ---
      rowmax = 0.0d0
      imaxp1 = imax + 1
      do 40 j = imaxp1, k
        rowmax = dmax1(rowmax,dcabsr(ar(imax,j),ai(imax,j)))
   40 continue
      if (imax .eq. 1) go to 50
      jmax = iyamax(imax-1,ar(1,imax),ai(1,imax),1)
      rowmax=dmax1(rowmax,dcabsr(ar(jmax,imax),ai(jmax,imax)))
   50 continue
      if (dcabsr(ar(imax,imax),ai(imax,imax)) .lt. alpha*rowmax)
     .go to 60
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

C ---   Column k is zero.  set info and iterate the loop. ---
      kpvt(k) = k
      info = k
      go to 190
  100 continue
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
  110 continue
  120 continue

C ---       Perform the elimination. ---
      do 130 jj = 1, km1
        j = k - jj
C     complex divide ... mulk = -a(j,k)/a(k,k) ...
        call cdiv(-ar(j,k),-ai(j,k),
     .  ar(k,k),ai(k,k),mulk(1),mulk(2))
        call yyaxpy(j,mulk(1),-mulk(2),ar(1,k),ai(1,k),1,
     .  ar(1,j),ai(1,j),1,.true.)
        ai(j,j) = 0
        ar(j,k) = mulk(1)
        ai(j,k) = mulk(2)
  130 continue

C --- Set the pivot array. ---
      kpvt(k) = k
      if (swap) kpvt(k) = imax
      go to 190
  140 continue

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
  150 continue
      t(1) = ar(k-1,k)
      t(2) = ai(k-1,k)
      ar(k-1,k) = ar(imax,k)
      ai(k-1,k) = ai(imax,k)
      ar(imax,k) = t(1)
      ai(imax,k) = t(2)
  160 continue

C ---       Perform the elimination. ---
      km2 = k - 2
      if (km2 .eq. 0) go to 180
c ... ak = a(k,k)/a(k-1,k),  akm1 = a(k-1,k-1)/dconjg(a(k-1,k)) ...
      call cdiv(ar(k,k),ai(k,k),ar(k-1,k),ai(k-1,k),ak(1),ak(2))
      call cdiv(ar(k-1,k-1),ai(k-1,k-1),ar(k-1,k),-ai(k-1,k),
     .akm1(1),akm1(2))
c ... denom = 1.0d0 - ak*akm1 ...
      denom(1) = 1 - ak(1)*akm1(1) + ak(2)*akm1(2)
      denom(2) =   - ak(1)*akm1(2) - ak(2)*akm1(1)
      do 170 jj = 1, km2
        j = km1 - jj
c ... bk = a(j,k)/a(k-1,k),  bkm1 = a(j,k-1)/dconjg(a(k-1,k))
        call cdiv(ar(j,k),ai(j,k),
     .  ar(k-1,k),ai(k-1,k),bk(1),bk(2))
        call cdiv(ar(j,k-1),ai(j,k-1),
     .  ar(k-1,k),-ai(k-1,k),bkm1(1),bkm1(2))
c ... mulk = (akm1*bk - bkm1)/denom
        call cpy(akm1(1),akm1(2),bk(1),bk(2),mulk(1),mulk(2))
        mulk(1) = mulk(1) - bkm1(1)
        mulk(2) = mulk(2) - bkm1(2)
        call cdiv(mulk(1),mulk(2),denom(1),denom(2),
     .  mulk(1),mulk(2))
c ... mulkm1 = (ak*bkm1 - bk)/denom
        call cpy(ak(1),ak(2),bkm1(1),bkm1(2),mulkm1(1),mulkm1(2))
        mulkm1(1) = mulkm1(1) - bk(1)
        mulkm1(2) = mulkm1(2) - bk(2)
        call cdiv(mulkm1(1),mulkm1(2),denom(1),denom(2),
     .  mulkm1(1),mulkm1(2))
        call yyaxpy(j,mulk(1),-mulk(2),ar(1,k),ai(1,k),1,
     .  ar(1,j),ai(1,j),1,.true.)
        call yyaxpy(j,mulkm1(1),-mulkm1(2),ar(1,k-1),ai(1,k-1),1,
     .  ar(1,j),ai(1,j),1,.true.)
        ar(j,k) = mulk(1)
        ai(j,k) = mulk(2)
        ar(j,k-1) = mulkm1(1)
        ai(j,k-1) = mulkm1(2)
        ai(j,j) = 0
  170 continue
  180 continue

C ---    Set the pivot array. ---
      kpvt(k) = 1 - k
      if (swap) kpvt(k) = -imax
      kpvt(k-1) = kpvt(k)
  190 continue
      k = k - kstep
      go to 10
  200 continue
      end
      subroutine yyhidi(ar,ai,lda,n,kpvt,det,inert,wkr,wki,job)
C- Compute determinant, inertia and inverse using factors from yyhifa
      integer lda,n,job
      double precision ar(lda,lda),ai(lda,lda),wkr(1),wki(1)
      double precision det(2)
      integer kpvt(1),inert(3)
c
c     on entry
c
c        ar,ai   output from yyhifa
c
c        lda     integer
c                the leading dimension of the array a.
c
c        n       integer
c                the order of the matrix a.
c
c        kpvt    integer(n)
c                the pivot vector from yyhifa.
c
c        wkr,wki double precision(n)
c                work vectors.  contents destroyed.
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
c        or  yyhifa  has set  info .ne. 0 .
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
   10 continue
      if (nodet) go to 20
      det(1) = 1.0d0
      det(2) = 0.0d0
      ten = 10.0d0
   20 continue
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
   30   continue
        d = t
        t = 0.0d0
   40   continue
   50   continue
c
        if (noert) go to 60
        if (d .gt. 0.0d0) inert(1) = inert(1) + 1
        if (d .lt. 0.0d0) inert(2) = inert(2) + 1
        if (d .eq. 0.0d0) inert(3) = inert(3) + 1
   60   continue
c
        if (nodet) go to 120
        det(1) = d*det(1)
        if (det(1) .eq. 0.0d0) go to 110
   70   if (dabs(det(1)) .ge. 1.0d0) go to 80
        det(1) = ten*det(1)
        det(2) = det(2) - 1.0d0
        go to 70
   80   continue
   90   if (dabs(det(1)) .lt. ten) go to 100
        det(1) = det(1)/ten
        det(2) = det(2) + 1.0d0
        go to 90
  100   continue
  110   continue
  120   continue
  130 continue
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
      call dcopy(km1,ar(1,k),1,wkr,1)
      call dcopy(km1,ai(1,k),1,wki,1)
      do 160 j = 1, km1
        ar(j,k) = ddot(j,ar(1,j),1,wkr,1)
     .  + ddot(j,ai(1,j),1,wki,1)
        ai(j,k) = ddot(j,ar(1,j),1,wki,1)
     .  - ddot(j,ai(1,j),1,wkr,1)
        call yyaxpy(j-1,wkr(j),wki(j),ar(1,j),ai(1,j),1,
     .  ar(1,k),ai(1,k),1,.true.)
  160 continue
      ar(k,k) = ar(k,k)
     .+ddot(km1,wkr,1,ar(1,k),1)
     .+ddot(km1,wki,1,ai(1,k),1)
      ai(k,k) = 0
  170 continue
      kstep = 1
      go to 220
  180 continue

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
      call dcopy(km1,ar(1,k+1),1,wkr,1)
      call dcopy(km1,ai(1,k+1),1,wki,1)
      do 190 j = 1, km1
        ar(j,k+1) =
     .  ddot(j,ar(1,j),1,wkr,1) + 
     .  ddot(j,ai(1,j),1,wki,1)
        ai(j,k+1) =
     .  ddot(j,ar(1,j),1,wki,1) -
     .  ddot(j,ai(1,j),1,wkr,1)
        call yyaxpy(j-1,wkr(j),wki(j),ar(1,j),ai(1,j),1,
     .  ar(1,k+1),ai(1,k+1),1,.true.)
  190 continue
      ar(k+1,k+1) = ar(k+1,k+1)
     .+ddot(km1,wkr,1,ar(1,k+1),1)
     .+ddot(km1,wki,1,ai(1,k+1),1)
      ai(k+1,k+1) = 0
      ar(k,k+1) = ar(k,k+1)
     .+ddot(km1,ar(1,k),1,ar(1,k+1),1)
     .+ddot(km1,ai(1,k),1,ai(1,k+1),1)
      ai(k,k+1) = ai(k,k+1)
     .+ddot(km1,ar(1,k),1,ai(1,k+1),1)
     .-ddot(km1,ai(1,k),1,ar(1,k+1),1)
      call dcopy(km1,ar(1,k),1,wkr,1)
      call dcopy(km1,ai(1,k),1,wki,1)
      do 200 j = 1, km1
        ar(j,k) =
     .  ddot(j,ar(1,j),1,wkr,1)
     .  + ddot(j,ai(1,j),1,wki,1)
        ai(j,k) =
     .  ddot(j,ar(1,j),1,wki,1)
     .  - ddot(j,ai(1,j),1,wkr,1)
        call yyaxpy(j-1,wkr(j),wki(j),ar(1,j),ai(1,j),1,
     .  ar(1,k),ai(1,k),1,.true.)
  200 continue
      ar(k,k) = ar(k,k)
     .+ddot(km1,wkr,1,ar(1,k),1)
     .+ddot(km1,wki,1,ai(1,k),1)
      ai(k,k) = 0
  210 continue
      kstep = 2
  220 continue

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
  230 continue
      if (kstep .eq. 1) go to 240
      temp = ar(ks,k+1)
      ar(ks,k+1) = ar(k,k+1)
      ar(k,k+1) = temp
      temp = ai(ks,k+1)
      ai(ks,k+1) = ai(k,k+1)
      ai(k,k+1) = temp
  240 continue
  250 continue
      k = k + kstep
      go to 150
  260 continue
  270 continue
      end
      double precision function cdabs2(x1,x2)
      double precision x1,x2
      cdabs2 = dsqrt(x1**2 + x2**2)
      end

