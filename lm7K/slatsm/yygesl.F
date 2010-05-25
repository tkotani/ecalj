      subroutine yygesl(ar,ai,lda,n,ipvt,br,bi,job)
C     implicit none
      integer lda,n,ipvt(1),job
      double precision ar(lda,1),ai(lda,1),br(1),bi(1)
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) goto ...
c           do  10  j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c
c     internal variables
c
      double precision t(2)
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) goto 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
      if (nm1 .lt. 1) goto 30
      do  20  k = 1, nm1
        l = ipvt(k)
        t(1) = br(l)
        t(2) = bi(l)
        if (l .eq. k) goto 10
        br(l) = br(k)
        bi(l) = bi(k)
        br(k) = t(1)
        bi(k) = t(2)
   10   continue
        call yyaxpy(n-k,t,t(2),ar(k+1,k),ai(k+1,k),1,
     .  br(k+1),bi(k+1),1,.true.)
   20 continue
   30 continue
c
c        now solve  u*x = y
c
      do  40  kb = 1, n
        k = n + 1 - kb
C ...       b(k) = b(k)/a(k,k)
        call cdiv(br(k),bi(k),ar(k,k),ai(k,k),br(k),bi(k))
        t(1) = -br(k)
        t(2) = -bi(k)
        call yyaxpy(k-1,t,t(2),ar(1,k),ai(1,k),1,
     .  br(1),bi(1),1,.true.)
   40 continue
      goto 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
      do  60  k = 1, n
        call yydot(k-1,ar(1,k),ai(1,k),1,br(1),bi(1),1,t,t(2))
C ...       b(k) = (b(k) - t)/a(k,k)
        br(k) = (br(k) - t(1))
        bi(k) = (bi(k) - t(2))
        call cdiv(br(k),bi(k),ar(k,k),ai(k,k),br(k),bi(k))
   60 continue
c
c        now solve trans(l)*x = y
c
      if (nm1 .lt. 1) goto 90
      do  80  kb = 1, nm1
        k = n - kb
C      ...  b(k) = b(k) + dot(n-k,a(k+1,k),1,b(k+1),1)
        call yydot(n-k,ar(k+1,k),ai(k+1,k),1,br(k+1),bi(k+1),1,
     .  t,t(2))
        br(k) = br(k) + t(1)
        bi(k) = bi(k) + t(2)
        l = ipvt(k)
        if (l .eq. k) goto 70
        t(1) = br(l)
        t(2) = bi(l)
        br(l) = br(k)
        bi(l) = bi(k)
        br(k) = t(1)
        bi(k) = t(2)
   70   continue
   80 continue
   90 continue
  100 continue
      end

