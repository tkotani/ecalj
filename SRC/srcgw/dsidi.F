      subroutine dsidi(a,lda,n,kpvt,det,inert,work,job)
      integer lda,n,job
      double precision a(lda,1),work(1)
      double precision det(2)
      integer kpvt(1),inert(3)
c
c     dsidi computes the determinant, inertia and inverse
c     of a double precision symmetric matrix using the factors from
c     dsifa.
c
c     on entry
c
c        a       double precision(lda,n)
c                the output from dsifa.
c
c        lda     integer
c                the leading dimension of the array a.
c
c        n       integer
c                the order of the matrix a.
c
c        kpvt    integer(n)
c                the pivot vector from dsifa.
c
c        work    double precision(n)
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
c        and  dsico  has set rcond .eq. 0.0
c        or  dsifa  has set  info .ne. 0 .
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab
c
c     subroutines and functions
c
c     blas daxpy,dcopy,ddot,dswap
c     fortran dabs,iabs,mod
c
c     internal variables.
c
      double precision akkp1,ddot,temp
      double precision ten,d,t,ak,akp1
      integer j,jb,k,km1,ks,kstep
      logical noinv,nodet,noert
c
C#ifdefC CRAY
C      call ssidi(a,lda,n,kpvt,det,inert,work,job)
C#else
      noinv = mod(job,10) .eq. 0
      nodet = mod(job,100)/10 .eq. 0
      noert = mod(job,1000)/100 .eq. 0
c
      if (nodet .and. noert) goto 140
      if (noert) goto 10
      inert(1) = 0
      inert(2) = 0
      inert(3) = 0
   10 continue
      if (nodet) goto 20
      det(1) = 1.0d0
      det(2) = 0.0d0
      ten = 10.0d0
   20 continue
      t = 0.0d0
      do  130  k = 1, n
        d = a(k,k)
c
c           check if 1 by 1
c
        if (kpvt(k) .gt. 0) goto 50
c
c              2 by 2 block
c              use det (d  s)  =  (d/t * c - t) * t  ,  t = dabs(s)
c                      (s  c)
c              to avoid underflow/overflow troubles.
c              take two passes through scaling.  use  t  for flag.
c
        if (t .ne. 0.0d0) goto 30
        t = dabs(a(k,k+1))
        d = (d/t)*a(k+1,k+1) - t
        goto 40
   30   continue
        d = t
        t = 0.0d0
   40   continue
   50   continue
c
        if (noert) goto 60
        if (d .gt. 0.0d0) inert(1) = inert(1) + 1
        if (d .lt. 0.0d0) inert(2) = inert(2) + 1
        if (d .eq. 0.0d0) inert(3) = inert(3) + 1
   60   continue
c
        if (nodet) goto 120
        det(1) = d*det(1)
        if (det(1) .eq. 0.0d0) goto 110
   70   if (dabs(det(1)) .ge. 1.0d0) goto 80
        det(1) = ten*det(1)
        det(2) = det(2) - 1.0d0
        goto 70
   80   continue
   90   if (dabs(det(1)) .lt. ten) goto 100
        det(1) = det(1)/ten
        det(2) = det(2) + 1.0d0
        goto 90
  100   continue
  110   continue
  120   continue
  130 continue
  140 continue
c
c     compute inverse(a)
c
      if (noinv) goto 270
      k = 1
  150 if (k .gt. n) goto 260
      km1 = k - 1
      if (kpvt(k) .lt. 0) goto 180
c
c              1 by 1
c
      a(k,k) = 1.0d0/a(k,k)
      if (km1 .lt. 1) goto 170
      call dcopy(km1,a(1,k),1,work,1)
      do  160  j = 1, km1
        a(j,k) = ddot(j,a(1,j),1,work,1)
        call daxpy(j-1,work(j),a(1,j),1,a(1,k),1)
  160 continue
      a(k,k) = a(k,k) + ddot(km1,work,1,a(1,k),1)
  170 continue
      kstep = 1
      goto 220
  180 continue
c
c              2 by 2
c
      t = dabs(a(k,k+1))
      ak = a(k,k)/t
      akp1 = a(k+1,k+1)/t
      akkp1 = a(k,k+1)/t
      d = t*(ak*akp1 - 1.0d0)
      a(k,k) = akp1/d
      a(k+1,k+1) = ak/d
      a(k,k+1) = -akkp1/d
      if (km1 .lt. 1) goto 210
      call dcopy(km1,a(1,k+1),1,work,1)
      do  190  j = 1, km1
        a(j,k+1) = ddot(j,a(1,j),1,work,1)
        call daxpy(j-1,work(j),a(1,j),1,a(1,k+1),1)
  190 continue
      a(k+1,k+1) = a(k+1,k+1) + ddot(km1,work,1,a(1,k+1),1)
      a(k,k+1) = a(k,k+1) + ddot(km1,a(1,k),1,a(1,k+1),1)
      call dcopy(km1,a(1,k),1,work,1)
      do  200  j = 1, km1
        a(j,k) = ddot(j,a(1,j),1,work,1)
        call daxpy(j-1,work(j),a(1,j),1,a(1,k),1)
  200 continue
      a(k,k) = a(k,k) + ddot(km1,work,1,a(1,k),1)
  210 continue
      kstep = 2
  220 continue
c
c           swap
c
      ks = iabs(kpvt(k))
      if (ks .eq. k) goto 250
      call dswap(ks,a(1,ks),1,a(1,k),1)
      do  230  jb = ks, k
        j = k + ks - jb
        temp = a(j,k)
        a(j,k) = a(ks,j)
        a(ks,j) = temp
  230 continue
      if (kstep .eq. 1) goto 240
      temp = a(ks,k+1)
      a(ks,k+1) = a(k,k+1)
      a(k,k+1) = temp
  240 continue
  250 continue
      k = k + kstep
      goto 150
  260 continue
  270 continue
      return
C#endif
      end
