      module m_amix !nothing saved in this module. Bundle subroutines/functions for amix
      public:: amix 
      contains
      integer function amix(nelts,npmix,mmix,ido,beta,ipr,tm, !norm,kpvt,
     .wk,t,rmsdel)
C- Anderson mixing of a vector
C ----------------------------------------------------------------
Ci Inputs
Ci   npmix: +/- number of previous iterations to fold into mix
Ci         npmix = 0 => linear mixing (x* = x0)
Ci         For meaning of npmix<0, see Remarks
Ci   mmix: maximum number of previous iterations to fold into mix
Ci         (used for dimensioning of work array)
Ci   nelts:number of elements to mix
Ci   wk:   array of dimension (nelts,2+mmix,2) where:
Ci         wk(*,0,1) holds f(xi) (see remarks)
Ci         wk(*,i,1) holds  d(i) (see remarks) ,  i>0
Ci         wk(*,i,2) holds   xi  (see remarks) ,  i>=0
Ci   beta: new x is beta f(x*) + (1-beta) x*
Ci   tm:   upper limit to any tj:  if any tj exceeds tm, effective
Ci         npmix is decremented.
Ci   ido:  0: normal mixing; 1: calculate tj only; 2: mix with input tj
Ci         10: silly mixing (mix each element of vector independently)
Co Outputs
Co   f(x_i) => f(x_i+1); x_i => x_i+1
Co   wk(*,i,1) => wk(*,i+1,1) (thus f(x_i) => f(x_i+1))
Co   wk(*,i,2) => wk(*,i+1,2) (thus x_i => x_i+1)
Co   wk(*,0,1): f(x*)-x*; see Remarks
Co   wk(*,0,2): new x = x* + beta ( f(x*)-x* ) ; see Remarks
Co   rmsdel:rms difference between x_0 and f(x_0)
Co   amix:  returns effective npmix (see input tm)
Cr Remarks
Cr   Given a vector function f(x), where x is some
Cr   vector, we want to find x* such that f(x*) = x*.  We want to find
Cr   x* with the minimum number of computations of f(x).  Supposing
Cr   that we have x0,f(x0); x1,f(x1); x2,f(x2); ...; x_n+1,f(x_n+1).
Cr   (Usually x_j corresponds to x at the jth previous iteration.)
Cr   We take a linear combination x* of x_0, x_1, x_2, ... x_n that
Cr   minimizes <(x* - f(x*))^2>.  We then seek t_1, t_2, ... t_n in
Cr     x* = x_0 - \sum_j t_j (x_0 - x_j).                        (1)
Cr   To evaluate f(x*) we linearize d(x) = f(x)-x as
Cr     f(x*)-x*  =  d*  =  d_0  -  \sum_j t_j (d_0 - d_j)         (2)
Cr   Then \partial <(d*)^2> / \partial t_k = 0 =>
Cr     < (d_0 - \sum_j t_j (d_0 - d_j)) (d_0 - d_k) >  =  0      (3)
Cr   constitute a set of n simultaneous equations in the n unknowns t_k.
Cr   Note that d's enter into these equations, not the f's.
Cr   Given the t_k's, x* can be estimated from (2).  To dampen
Cr   instablities, a linear combination of (1) and (2) is taken as
Cr       beta f(x*)  +  (1-beta) x*  = beta d*  +  x*           (4)
Cr   beta is an input parameter, which can usually be taken to be 1.
Cr   If you want to be really timid, and constrain the program
Cr   to take a small step, you can mix  x_0 + beta d*, which as
Cr   beta -> 0 moves only a small amount away from x_0.  This feature
Cr   is set by making npmix < 0.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nelts,mmix,npmix,kpvt(mmix),ipr
      double precision norm(mmix,mmix),wk(nelts,0:mmix+1,2),t(mmix),
     .  beta,tm,rmsdel
C Local parameters
      integer i,j,nwk,inert(3),nmix,kelt,nmake,ido,iprint
      double precision ddot,det(2),sumsqr,sing,sinmax
      parameter (sinmax=100000d0)

      amix = 0
      if (nelts .eq. 0) return
C     if (ipr .ge. 20 .and. ido .ne. 2) print *
      nwk = nelts*(mmix+2)
      kelt = 0
C nmake is the number of elements mixed per matrix inversion
      nmake = nelts
      if (ido/10 .eq. 1) nmake = 1

C --  d_0 = f-x  =>  wk(*,0,1) --
      call daxpy(nelts,-1d0,wk(1,0,2),1,wk,1)

C --  Copy x_0 and d_0 to end of wk:  x*, d* constructed there --
      call dmcpy(wk,nwk,1,wk(1,mmix+1,1),nwk,1,nelts,2)

C --- Obtain the tj ---
   11 continue
      nmix = iabs(npmix)
      kelt = kelt+1

C --  Make < (d_0 - d_j) (d_0 - d_k) > and  < d_0 (d_0 - d_j) >  --
      if (ido .eq. 2) goto 40
   10 continue
      if (nmix .lt. 0) call rx( 'amix: bad nmix')


C -- Regular Anderson mixing branch --
      if (ido/10 .eq. 0) then
        sumsqr = 0
        do    i = 1, nmix
          t(i) = ddot(nelts,wk(1,0,1),1,wk(1,0,1),1) -
     .    ddot(nelts,wk(1,0,1),1,wk(1,i,1),1)
        do    j = 1, nmix
            norm(i,j) =  ddot(nelts,wk(1,0,1),1,wk(1,0,1),1)
     .      - ddot(nelts,wk(1,0,1),1,wk(1,j,1),1)
     .      - ddot(nelts,wk(1,i,1),1,wk(1,0,1),1)
     .      + ddot(nelts,wk(1,i,1),1,wk(1,j,1),1)
            sumsqr = sumsqr + norm(i,j)**2
        enddo
        enddo
        sumsqr = dsqrt(sumsqr)/(nmix+1)
C -- Silly branch --
      elseif (ido/10 .eq. 1) then
        do    i = 1, nmix
          t(i) = wk(kelt,0,1)*wk(kelt,0,1) - wk(kelt,0,1)*wk(kelt,i,1)
        do    j = 1, nmix
          norm(i,j) =  
     .          wk(kelt,0,1)*wk(kelt,0,1) - wk(kelt,0,1)*wk(kelt,j,1)
     .          - wk(kelt,i,1)*wk(kelt,0,1) + wk(kelt,i,1)*wk(kelt,j,1)
        enddo
        enddo
      endif

C --  Solve the simultaneous equations for tj --
      call dsifa(norm,mmix,nmix,kpvt,i)
      if (i .ne. 0) then
        sing = sinmax + 1
      else
        call dsisl(norm,mmix,nmix,kpvt,t)
        call dsidi(norm,mmix,nmix,kpvt,det,inert,10)
        sing = dabs(sumsqr/(det(1)*10**det(2)))
      endif

C --  Handle singular normal matrix --
      if (sing .gt. sinmax) then
        t(nmix) = 0
        nmix = nmix-1
        if (ipr .ge. 20) print 337, sinmax, nmix
  337   format(' AMIX: condition of normal eqns >',f7.0,
     .         ' Reducing nmix to',2i2)
        goto 10
      endif

C --  Reduce nmix if any t_j exceeds tm --
      do  30  j = 1, nmix
        if (dabs(t(j)) .le. dabs(tm)) goto 30
        if (ipr .ge. 20) print 338,  nmix-1, (t(i), i=1, nmix)
  338   format(' AMIX: Reducing nmix to',i3,
     .  ': t_j exceeds tm: tj=',7(f9.5,1x))
        t(nmix) = 0
        nmix = nmix-1
        goto 10
   30 continue

C --- Do mixing unless ido=1  ---
   40 continue
      amix = nmix
      if (ido .eq. 1) then
C   ... Restore f = d_0 + x  =>  wk(*,0,1)
        call daxpy(nelts,1d0,wk(1,0,2),1,wk,1)
        return
      endif

C --- Make (d,x)* = (d,x)_0 - \sum_j t_j ((d,x)_0 - (d,x)_j) ---
      do  45  j = 1, nmix
        call daxpy(nmake,-t(j),wk(kelt,0,1),1,wk(kelt,mmix+1,1),1)
        call daxpy(nmake, t(j),wk(kelt,j,1),1,wk(kelt,mmix+1,1),1)
        if (npmix .gt. 0) then
          call daxpy(nmake,-t(j),wk(kelt,0,2),1,wk(kelt,mmix+1,2),1)
          call daxpy(nmake, t(j),wk(kelt,j,2),1,wk(kelt,mmix+1,2),1)
        endif
   45 continue

C -- Do next element for silly case --
      if (ido/10 .eq. 1) then
        if (nmix .gt. 0 .and. iprint() .ge. 40)
     .  write(*,135) kelt, (t(j), j=1,nmix)
  135   format(i4,3x,'tj:',7(f8.5,2x))
        if (kelt .lt. nelts) goto 11
      endif

C --  Copy arrays to new positions --
      do  50  i = mmix, 1, -1
        call dmcpy(wk(1,i-1,1),nwk,1,wk(1,i,1),nwk,1,nelts,2)
   50 continue

C --  x* + beta d*  (or x + beta d* if npmix<0) --
      call daxpy(nelts,beta,wk(1,mmix+1,1),1,wk(1,mmix+1,2),1)

C --  Calculate rms change --
      rmsdel = dsqrt(ddot(nelts,wk,1,wk,1)/nelts)

C --- Printout ---
      if (ipr .lt. 30) goto 60
      write(*,133) nmix,mmix,nelts,beta,tm,rmsdel
  133 format(' AMIX: nmix=',i1,' mmix=',i1,'  nelts=',i6,
     .       '  beta=',f7.5,'  tm=',f8.5,'  rmsdel=',1pd9.2)
      if (nmix .gt. 0) write(*,134) (t(j), j=1,nmix)
  134 format(3x,'tj:',7(f8.5,2x))

      if (ipr .lt. 61 .and. (ipr .lt. 41 .or. nelts .gt. 100)) goto 60
      write(*,110)
      do  12  i = 1, nelts
        if (dabs(wk(i,0,1)) + dabs(wk(i,mmix+1,2)-wk(i,0,2)) .ge. 5d-7)
     .  write(*,111) i,wk(i,0,2),wk(i,0,2)+wk(i,0,1),
     .                 wk(i,0,1),wk(i,mmix+1,2)
   12 continue

C --- Restore d* and x* + beta d* from end of wk --
   60 call dmcpy(wk(1,mmix+1,1),nwk,1,wk,nwk,1,nelts,2)

  104 format(1p,4d18.11)
  111 format(i5,4f14.6)
  110 format(14x,'Old',11x,' New',9x,'Diff',10x,'Mixed')
      end

      subroutine dmcpy(a,nca,nra,b,ncb,nrb,n,m)
C- General matrix copy
C ----------------------------------------------------------------
Ci Inputs:
Ci   a,nca,nra is the source matrix and respectively the number of
Ci      elements separating columns and rows.
Ci   b,ncb,nrb is the dest. matrix and respectively the number of
Ci      elements separating columns and rows.
Ci   n,m: the number of columns and rows, respectively, to calculate
Co Outputs:
Co   result matrix stored in c
Cr Remarks:
Cr   This is a general-purpose matrix copy routine,
Cr   copying a subblock of matrix a to a subblock of matrix b.
Cr   Normally matrix nc{a,b} is the row dimension of matrix {a,b}
Cr   and nr{a,b} is 1.  Reverse nr and nc for a transposed matrix.
Cr   Arrays are locally one-dimensional so as to optimize inner loop.
Cr
Cr   Example: Set 3-by-2 block of matrix c to constant z
Cr     call dmcpy(z,0,0,c,nc,1,3,2)
Cr   Here scalar z is represented by an array of 0 dimension
Cr   Example: copy nbas-by-3 'wk' to its transpose 'bas', and the reverse:
Cr     call dmcpy(wk,nbas,1,bas,1,3,nbas,3)
Cr     call dmcpy(bas,1,3,wk,nbas,1,nbas,3)
C ----------------------------------------------------------------
C     implicit none
      integer nca,nra,ncb,nrb,n,m
      double precision a(0:*), b(0:*)
      integer i,j,ia,ib
      if (nra .eq. 1 .and. nrb .eq. 1) then
        do  11  j = 0, m-1
          ia = j*nca
          ib = j*ncb
          do  10  i = 0, n-1
            b(i+ib) = a(i+ia)
   10   continue
   11   continue
        return
      endif
      do  21  i = n-1, 0, -1
        ia = i*nra+m*nca
        ib = i*nrb+m*ncb
        do  20  j = m-1, 0, -1 
          ia = ia-nca
          ib = ib-ncb
          b(ib) = a(ia)
   20 continue
   21 continue
      end

      subroutine dsisl(a,lda,n,kpvt,b)
      integer lda,n,kpvt(1)
      double precision a(lda,1),b(1)
c
c     dsisl solves the double precision symmetric system
c     a * x = b
c     using the factors computed by dsifa.
c
c     on entry
c
c        a       double precision(lda,n)
c                the output from dsifa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        kpvt    integer(n)
c                the pivot vector from dsifa.
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
c        a division by zero may occur if  dsico  has set rcond .eq. 0.0
c        or  dsifa  has set info .ne. 0  .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dsifa(a,lda,n,kpvt,info)
c           if (info .ne. 0) goto ...
c           do  10  j = 1, p
c              call dsisl(a,lda,n,kpvt,c(1,j))
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
      integer k,kp
c
c     loop backward applying the transformations and
c     d inverse to b.
c
c#if CRAY
c      call ssisl(a,lda,n,kpvt,b)
c#else
      k = n
   10 if (k .eq. 0) goto 80
      if (kpvt(k) .lt. 0) goto 40
c
c           1 x 1 pivot block.
c
      if (k .eq. 1) goto 30
      kp = kpvt(k)
      if (kp .eq. k) goto 20
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
      call daxpy(k-1,b(k),a(1,k),1,b(1),1)
   30 continue
c
c           apply d inverse.
c
      b(k) = b(k)/a(k,k)
      k = k - 1
      goto 70
   40 continue
c
c           2 x 2 pivot block.
c
      if (k .eq. 2) goto 60
      kp = iabs(kpvt(k))
      if (kp .eq. k - 1) goto 50
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
      call daxpy(k-2,b(k),a(1,k),1,b(1),1)
      call daxpy(k-2,b(k-1),a(1,k-1),1,b(1),1)
   60 continue
c
c           apply d inverse.
c
      ak = a(k,k)/a(k-1,k)
      akm1 = a(k-1,k-1)/a(k-1,k)
      bk = b(k)/a(k-1,k)
      bkm1 = b(k-1)/a(k-1,k)
      denom = ak*akm1 - 1.0d0
      b(k) = (akm1*bk - bkm1)/denom
      b(k-1) = (ak*bkm1 - bk)/denom
      k = k - 2
   70 continue
      goto 10
   80 continue
c
c     loop forward applying the transformations.
c
      k = 1
   90 if (k .gt. n) goto 160
      if (kpvt(k) .lt. 0) goto 120
c
c           1 x 1 pivot block.
c
      if (k .eq. 1) goto 110
c
c              apply the transformation.
c
      b(k) = b(k) + ddot(k-1,a(1,k),1,b(1),1)
      kp = kpvt(k)
      if (kp .eq. k) goto 100
c
c                 interchange.
c
      temp = b(k)
      b(k) = b(kp)
      b(kp) = temp
  100 continue
  110 continue
      k = k + 1
      goto 150
  120 continue
c
c           2 x 2 pivot block.
c
      if (k .eq. 1) goto 140
c
c              apply the transformation.
c
      b(k) = b(k) + ddot(k-1,a(1,k),1,b(1),1)
      b(k+1) = b(k+1) + ddot(k-1,a(1,k+1),1,b(1),1)
      kp = iabs(kpvt(k))
      if (kp .eq. k) goto 130
c
c                 interchange.
c
      temp = b(k)
      b(k) = b(kp)
      b(kp) = temp
  130 continue
  140 continue
      k = k + 2
  150 continue
      goto 90
  160 continue
      return
c#endif
      end

      subroutine dsidi(a,lda,n,kpvt,det,inert,job)
      integer lda,n,job
      double precision a(lda,1),work(n)
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
c#if CRAY
c      call ssidi(a,lda,n,kpvt,det,inert,work,job)
c#else
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
c#endif
      end

      subroutine dsifa(a,lda,n,kpvt,info)
      integer lda,n,kpvt(1),info
      double precision a(lda,1)
c
c     dsifa factors a double precision symmetric matrix by elimination
c     with symmetric pivoting.
c
c     to solve  a*x = b , follow dsifa by dsisl.
c     to compute  inverse(a)*c , follow dsifa by dsisl.
c     to compute  determinant(a) , follow dsifa by dsidi.
c     to compute  inertia(a) , follow dsifa by dsidi.
c     to compute  inverse(a) , follow dsifa by dsidi.
c
c     on entry
c
c        a       double precision(lda,n)
c                the symmetric matrix to be factored.
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
c                     but it does indicate that dsisl or dsidi may
c                     divide by zero if called.
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
      integer imax,imaxp1,j,jj,jmax,k,km1,km2,kstep,idamax
      logical swap
c
c     initialize
c
c     alpha is used in choosing pivot block size.
c#if CRAY
c      call ssifa(a,lda,n,kpvt,info)
c#else
      alpha = (1.0d0 + dsqrt(17.0d0))/8.0d0
c
      info = 0
c
c     main loop on k, which goes from n to 1.
c
      k = n
   10 continue
c
c        leave the loop if k=0 or k=1.
c
c     ...exit
      if (k .eq. 0) goto 200
      if (k .gt. 1) goto 20
      kpvt(1) = 1
      if (a(1,1) .eq. 0.0d0) info = 1
c     ......exit
      goto 200
   20 continue
c
c        this section of code determines the kind of
c        elimination to be performed.  when it is completed,
c        kstep will be set to the size of the pivot block, and
c        swap will be set to .true. if an interchange is
c        required.
c
      km1 = k - 1
      absakk = dabs(a(k,k))
c
c        determine the largest off-diagonal element in
c        column k.
c
      imax = idamax(k-1,a(1,k),1)
      colmax = dabs(a(imax,k))
      if (absakk .lt. alpha*colmax) goto 30
      kstep = 1
      swap = .false.
      goto 90
   30 continue
c
c           determine the largest off-diagonal element in
c           row imax.
c
      rowmax = 0.0d0
      imaxp1 = imax + 1
      do  40  j = imaxp1, k
        rowmax = dmax1(rowmax,dabs(a(imax,j)))
   40 continue
      if (imax .eq. 1) goto 50
      jmax = idamax(imax-1,a(1,imax),1)
      rowmax = dmax1(rowmax,dabs(a(jmax,imax)))
   50 continue
      if (dabs(a(imax,imax)) .lt. alpha*rowmax) goto 60
      kstep = 1
      swap = .true.
      goto 80
   60 continue
      if (absakk .lt. alpha*colmax*(colmax/rowmax)) goto 70
      kstep = 1
      swap = .false.
      goto 80
   70 continue
      kstep = 2
      swap = imax .ne. km1
   80 continue
   90 continue
      if (dmax1(absakk,colmax) .ne. 0.0d0) goto 100
c
c           column k is zero.  set info and iterate the loop.
c
      kpvt(k) = k
      info = k
      goto 190
  100 continue
      if (kstep .eq. 2) goto 140
c
c           1 x 1 pivot block.
c
      if (.not.swap) goto 120
c
c              perform an interchange.
c
      call dswap(imax,a(1,imax),1,a(1,k),1)
      do  110  jj = imax, k
        j = k + imax - jj
        t = a(j,k)
        a(j,k) = a(imax,j)
        a(imax,j) = t
  110 continue
  120 continue
c
c           perform the elimination.
c
      do  130  jj = 1, km1
        j = k - jj
        mulk = -a(j,k)/a(k,k)
        t = mulk
        call daxpy(j,t,a(1,k),1,a(1,j),1)
        a(j,k) = mulk
  130 continue
c
c           set the pivot array.
c
      kpvt(k) = k
      if (swap) kpvt(k) = imax
      goto 190
  140 continue
c
c           2 x 2 pivot block.
c
      if (.not.swap) goto 160
c
c              perform an interchange.
c
      call dswap(imax,a(1,imax),1,a(1,k-1),1)
      do  150  jj = imax, km1
        j = km1 + imax - jj
        t = a(j,k-1)
        a(j,k-1) = a(imax,j)
        a(imax,j) = t
  150 continue
      t = a(k-1,k)
      a(k-1,k) = a(imax,k)
      a(imax,k) = t
  160 continue
c
c           perform the elimination.
c
      km2 = k - 2
      if (km2 .eq. 0) goto 180
      ak = a(k,k)/a(k-1,k)
      akm1 = a(k-1,k-1)/a(k-1,k)
      denom = 1.0d0 - ak*akm1
      do  170  jj = 1, km2
        j = km1 - jj
        bk = a(j,k)/a(k-1,k)
        bkm1 = a(j,k-1)/a(k-1,k)
        mulk = (akm1*bk - bkm1)/denom
        mulkm1 = (ak*bkm1 - bk)/denom
        t = mulk
        call daxpy(j,t,a(1,k),1,a(1,j),1)
        t = mulkm1
        call daxpy(j,t,a(1,k-1),1,a(1,j),1)
        a(j,k) = mulk
        a(j,k-1) = mulkm1
  170 continue
  180 continue
c
c           set the pivot array.
c
      kpvt(k) = 1 - k
      if (swap) kpvt(k) = -imax
      kpvt(k-1) = kpvt(k)
  190 continue
      k = k - kstep
      goto 10
  200 continue
      return
c#endif
      end

      end module m_amix
