module m_amix !nothing saved in this module. Bundle subroutines/functions for amix
  public:: amix
contains
  integer function amix(nelts,npmix,mmix,ido,beta,ipr,tm,  wk,t,rmsdel)
    !- Anderson mixing of a vector
    ! ----------------------------------------------------------------
    !i Inputs
    !i   npmix: +/- number of previous iterations to fold into mix
    !i         npmix = 0 => linear mixing (x* = x0)
    !i         For meaning of npmix<0, see Remarks
    !i   mmix: maximum number of previous iterations to fold into mix
    !i         (used for dimensioning of work array)
    !i   nelts:number of elements to mix
    !i   wk:   array of dimension (nelts,2+mmix,2) where:
    !i         wk(*,0,1) holds f(xi) (see remarks)
    !i         wk(*,i,1) holds  d(i) (see remarks) ,  i>0
    !i         wk(*,i,2) holds   xi  (see remarks) ,  i>=0
    !i   beta: new x is beta f(x*) + (1-beta) x*
    !i   tm:   upper limit to any tj:  if any tj exceeds tm, effective
    !i         npmix is decremented.
    !i   ido:  0: normal mixing; 1: calculate tj only; 2: mix with input tj
    !i         10: silly mixing (mix each element of vector independently)
    !o Outputs
    !o   f(x_i) => f(x_i+1); x_i => x_i+1
    !o   wk(*,i,1) => wk(*,i+1,1) (thus f(x_i) => f(x_i+1))
    !o   wk(*,i,2) => wk(*,i+1,2) (thus x_i => x_i+1)
    !o   wk(*,0,1): f(x*)-x*; see Remarks
    !o   wk(*,0,2): new x = x* + beta ( f(x*)-x* ) ; see Remarks
    !o   rmsdel:rms difference between x_0 and f(x_0)
    !o   amix:  returns effective npmix (see input tm)
    !r Remarks
    !r   Given a vector function f(x), where x is some
    !r   vector, we want to find x* such that f(x*) = x*.  We want to find
    !r   x* with the minimum number of computations of f(x).  Supposing
    !r   that we have x0,f(x0); x1,f(x1); x2,f(x2); ...; x_n+1,f(x_n+1).
    !r   (Usually x_j corresponds to x at the jth previous iteration.)
    !r   We take a linear combination x* of x_0, x_1, x_2, ... x_n that
    !r   minimizes <(x* - f(x*))^2>.  We then seek t_1, t_2, ... t_n in
    !r     x* = x_0 - \sum_j t_j (x_0 - x_j).                        (1)
    !r   To evaluate f(x*) we linearize d(x) = f(x)-x as
    !r     f(x*)-x*  =  d*  =  d_0  -  \sum_j t_j (d_0 - d_j)         (2)
    !r   Then \partial <(d*)^2> / \partial t_k = 0 =>
    !r     < (d_0 - \sum_j t_j (d_0 - d_j)) (d_0 - d_k) >  =  0      (3)
    !r   constitute a set of n simultaneous equations in the n unknowns t_k.
    !r   Note that d's enter into these equations, not the f's.
    !r   Given the t_k's, x* can be estimated from (2).  To dampen
    !r   instablities, a linear combination of (1) and (2) is taken as
    !r       beta f(x*)  +  (1-beta) x*  = beta d*  +  x*           (4)
    !r   beta is an input parameter, which can usually be taken to be 1.
    !r   If you want to be really timid, and constrain the program
    !r   to take a small step, you can mix  x_0 + beta d*, which as
    !r   beta -> 0 moves only a small amount away from x_0.  This feature
    !r   is set by making npmix < 0.
    ! ----------------------------------------------------------------
    !     implicit none
    ! Passed parameters
    integer :: nelts,mmix,npmix,kpvt(mmix),ipr
    double precision :: norm(mmix,mmix),wk(nelts,0:mmix+1,2),t(mmix), &
         beta,tm,rmsdel
    ! Local parameters
    integer :: i,j,nwk,inert(3),nmix,kelt,nmake,ido,iprint
    double precision :: ddot,det(2),sumsqr,sing,sinmax
    parameter (sinmax=100000d0)

    amix = 0
    if (nelts == 0) return
    !     if (ipr .ge. 20 .and. ido .ne. 2) print *
    nwk = nelts*(mmix+2)
    kelt = 0
    ! nmake is the number of elements mixed per matrix inversion
    nmake = nelts
    if (ido/10 == 1) nmake = 1

    ! --  d_0 = f-x  =>  wk(*,0,1) --
    call daxpy(nelts,-1d0,wk(1,0,2),1,wk,1)

    ! --  Copy x_0 and d_0 to end of wk:  x*, d* constructed there --
    call dmcpy(wk,nwk,1,wk(1,mmix+1,1),nwk,1,nelts,2)

    ! --- Obtain the tj ---
11  continue
    nmix = iabs(npmix)
    kelt = kelt+1

    ! --  Make < (d_0 - d_j) (d_0 - d_k) > and  < d_0 (d_0 - d_j) >  --
    if (ido == 2) goto 40
10  continue
    if (nmix < 0) call rx( 'amix: bad nmix')


    ! -- Regular Anderson mixing branch --
    if (ido/10 == 0) then
       sumsqr = 0
       do    i = 1, nmix
          t(i) = ddot(nelts,wk(1,0,1),1,wk(1,0,1),1) - &
               ddot(nelts,wk(1,0,1),1,wk(1,i,1),1)
          do    j = 1, nmix
             norm(i,j) =  ddot(nelts,wk(1,0,1),1,wk(1,0,1),1) &
                  - ddot(nelts,wk(1,0,1),1,wk(1,j,1),1) &
                  - ddot(nelts,wk(1,i,1),1,wk(1,0,1),1) &
                  + ddot(nelts,wk(1,i,1),1,wk(1,j,1),1)
             sumsqr = sumsqr + norm(i,j)**2
          enddo
       enddo
       sumsqr = dsqrt(sumsqr)/(nmix+1)
       ! -- Silly branch --
    elseif (ido/10 == 1) then
       do    i = 1, nmix
          t(i) = wk(kelt,0,1)*wk(kelt,0,1) - wk(kelt,0,1)*wk(kelt,i,1)
          do    j = 1, nmix
             norm(i,j) = &
                  wk(kelt,0,1)*wk(kelt,0,1) - wk(kelt,0,1)*wk(kelt,j,1) &
                  - wk(kelt,i,1)*wk(kelt,0,1) + wk(kelt,i,1)*wk(kelt,j,1)
          enddo
       enddo
    endif

    ! --  Solve the simultaneous equations for tj --
    call dsifa(norm,mmix,nmix,kpvt,i)
    if (i /= 0) then
       sing = sinmax + 1
    else
       call dsisl(norm,mmix,nmix,kpvt,t)
       call dsidi(norm,mmix,nmix,kpvt,det,inert,10)
       sing = dabs(sumsqr/(det(1)*10**det(2)))
    endif

    ! --  Handle singular normal matrix --
    if (sing > sinmax) then
       t(nmix) = 0
       nmix = nmix-1
       if (ipr >= 20) print 337, sinmax, nmix
337    format(' AMIX: condition of normal eqns >',f7.0, &
            ' Reducing nmix to',2i2)
       goto 10
    endif

    ! --  Reduce nmix if any t_j exceeds tm --
    do  30  j = 1, nmix
       if (dabs(t(j)) <= dabs(tm)) goto 30
       if (ipr >= 20) print 338,  nmix-1, (t(i), i=1, nmix)
338    format(' AMIX: Reducing nmix to',i3, &
            ': t_j exceeds tm: tj=',7(f9.5,1x))
       t(nmix) = 0
       nmix = nmix-1
       goto 10
30  enddo

    ! --- Do mixing unless ido=1  ---
40  continue
    amix = nmix
    if (ido == 1) then
       !   ... Restore f = d_0 + x  =>  wk(*,0,1)
       call daxpy(nelts,1d0,wk(1,0,2),1,wk,1)
       return
    endif

    ! --- Make (d,x)* = (d,x)_0 - \sum_j t_j ((d,x)_0 - (d,x)_j) ---
    do  45  j = 1, nmix
       call daxpy(nmake,-t(j),wk(kelt,0,1),1,wk(kelt,mmix+1,1),1)
       call daxpy(nmake, t(j),wk(kelt,j,1),1,wk(kelt,mmix+1,1),1)
       if (npmix > 0) then
          call daxpy(nmake,-t(j),wk(kelt,0,2),1,wk(kelt,mmix+1,2),1)
          call daxpy(nmake, t(j),wk(kelt,j,2),1,wk(kelt,mmix+1,2),1)
       endif
45  enddo

    ! -- Do next element for silly case --
    if (ido/10 == 1) then
       if (nmix > 0 .AND. iprint() >= 40) &
            write(*,135) kelt, (t(j), j=1,nmix)
135    format(i4,3x,'tj:',7(f8.5,2x))
       if (kelt < nelts) goto 11
    endif

    ! --  Copy arrays to new positions --
    do  50  i = mmix, 1, -1
       call dmcpy(wk(1,i-1,1),nwk,1,wk(1,i,1),nwk,1,nelts,2)
50  enddo

    ! --  x* + beta d*  (or x + beta d* if npmix<0) --
    call daxpy(nelts,beta,wk(1,mmix+1,1),1,wk(1,mmix+1,2),1)

    ! --  Calculate rms change --
    rmsdel = dsqrt(ddot(nelts,wk,1,wk,1)/nelts)

    ! --- Printout ---
    if (ipr < 30) goto 60
    write(*,133) nmix,mmix,nelts,beta,tm,rmsdel
133 format(' AMIX: nmix=',i1,' mmix=',i1,'  nelts=',i6, &
         '  beta=',f7.5,'  tm=',f8.5,'  rmsdel=',1pd9.2)
    if (nmix > 0) write(*,134) (t(j), j=1,nmix)
134 format(3x,'tj:',7(f8.5,2x))

    if (ipr < 61 .AND. (ipr < 41 .OR. nelts > 100)) goto 60
    write(*,110)
    do  12  i = 1, nelts
       if (dabs(wk(i,0,1)) + dabs(wk(i,mmix+1,2)-wk(i,0,2)) >= 5d-7) &
            write(*,111) i,wk(i,0,2),wk(i,0,2)+wk(i,0,1), &
            wk(i,0,1),wk(i,mmix+1,2)
12  enddo

    ! --- Restore d* and x* + beta d* from end of wk --
60  call dmcpy(wk(1,mmix+1,1),nwk,1,wk,nwk,1,nelts,2)

104 format(1p,4d18.11)
111 format(i5,4f14.6)
110 format(14x,'Old',11x,' New',9x,'Diff',10x,'Mixed')
  end function amix

  subroutine dmcpy(a,nca,nra,b,ncb,nrb,n,m)
    !- General matrix copy
    ! ----------------------------------------------------------------
    !i Inputs:
    !i   a,nca,nra is the source matrix and respectively the number of
    !i      elements separating columns and rows.
    !i   b,ncb,nrb is the dest. matrix and respectively the number of
    !i      elements separating columns and rows.
    !i   n,m: the number of columns and rows, respectively, to calculate
    !o Outputs:
    !o   result matrix stored in c
    !r Remarks:
    !r   This is a general-purpose matrix copy routine,
    !r   copying a subblock of matrix a to a subblock of matrix b.
    !r   Normally matrix nc{a,b} is the row dimension of matrix {a,b}
    !r   and nr{a,b} is 1.  Reverse nr and nc for a transposed matrix.
    !r   Arrays are locally one-dimensional so as to optimize inner loop.
    !r
    !r   Example: Set 3-by-2 block of matrix c to constant z
    !r     call dmcpy(z,0,0,c,nc,1,3,2)
    !r   Here scalar z is represented by an array of 0 dimension
    !r   Example: copy nbas-by-3 'wk' to its transpose 'bas', and the reverse:
    !r     call dmcpy(wk,nbas,1,bas,1,3,nbas,3)
    !r     call dmcpy(bas,1,3,wk,nbas,1,nbas,3)
    ! ----------------------------------------------------------------
    !     implicit none
    integer :: nca,nra,ncb,nrb,n,m
    double precision :: a(0:*), b(0:*)
    integer :: i,j,ia,ib
    if (nra == 1 .AND. nrb == 1) then
       do  11  j = 0, m-1
          ia = j*nca
          ib = j*ncb
          do  10  i = 0, n-1
             b(i+ib) = a(i+ia)
10        enddo
11     enddo
       return
    endif
    do  21  i = n-1, 0, -1
       ia = i*nra+m*nca
       ib = i*nrb+m*ncb
       do  20  j = m-1, 0, -1
          ia = ia-nca
          ib = ib-ncb
          b(ib) = a(ia)
20     enddo
21  enddo
  end subroutine dmcpy

  subroutine dsisl(a,lda,n,kpvt,b)
    integer :: lda,n,kpvt(1)
    double precision :: a(lda,1),b(1)

    !     dsisl solves the double precision symmetric system
    !     a * x = b
    !     using the factors computed by dsifa.

    !     on entry

    !        a       double precision(lda,n)
    !                the output from dsifa.

    !        lda     integer
    !                the leading dimension of the array  a .

    !        n       integer
    !                the order of the matrix  a .

    !        kpvt    integer(n)
    !                the pivot vector from dsifa.

    !        b       double precision(n)
    !                the right hand side vector.

    !     on return

    !        b       the solution vector  x .

    !     error condition

    !        a division by zero may occur if  dsico  has set rcond .eq. 0.0
    !        or  dsifa  has set info .ne. 0  .

    !     to compute  inverse(a) * c  where  c  is a matrix
    !     with  p  columns
    !           call dsifa(a,lda,n,kpvt,info)
    !           if (info .ne. 0) goto ...
    !           do  10  j = 1, p
    !              call dsisl(a,lda,n,kpvt,c(1,j))
    !        10 continue

    !     linpack. this version dated 08/14/78 .
    !     james bunch, univ. calif. san diego, argonne nat. lab.

    !     subroutines and functions

    !     blas daxpy,ddot
    !     fortran iabs

    !     internal variables.

    double precision :: ak,akm1,bk,bkm1,ddot,denom,temp
    integer :: k,kp

    !     loop backward applying the transformations and
    !     d inverse to b.

    ! if CRAY
    !      call ssisl(a,lda,n,kpvt,b)
    ! else
    k = n
10  if (k == 0) goto 80
    if (kpvt(k) < 0) goto 40

    !           1 x 1 pivot block.

    if (k == 1) goto 30
    kp = kpvt(k)
    if (kp == k) goto 20

    !                 interchange.

    temp = b(k)
    b(k) = b(kp)
    b(kp) = temp
20  continue

    !              apply the transformation.

    call daxpy(k-1,b(k),a(1,k),1,b(1),1)
30  continue

    !           apply d inverse.

    b(k) = b(k)/a(k,k)
    k = k - 1
    goto 70
40  continue

    !           2 x 2 pivot block.

    if (k == 2) goto 60
    kp = iabs(kpvt(k))
    if (kp == k - 1) goto 50

    !                 interchange.

    temp = b(k-1)
    b(k-1) = b(kp)
    b(kp) = temp
50  continue

    !              apply the transformation.

    call daxpy(k-2,b(k),a(1,k),1,b(1),1)
    call daxpy(k-2,b(k-1),a(1,k-1),1,b(1),1)
60  continue

    !           apply d inverse.

    ak = a(k,k)/a(k-1,k)
    akm1 = a(k-1,k-1)/a(k-1,k)
    bk = b(k)/a(k-1,k)
    bkm1 = b(k-1)/a(k-1,k)
    denom = ak*akm1 - 1.0d0
    b(k) = (akm1*bk - bkm1)/denom
    b(k-1) = (ak*bkm1 - bk)/denom
    k = k - 2
70  continue
    goto 10
80  continue

    !     loop forward applying the transformations.

    k = 1
90  if (k > n) goto 160
    if (kpvt(k) < 0) goto 120

    !           1 x 1 pivot block.

    if (k == 1) goto 110

    !              apply the transformation.

    b(k) = b(k) + ddot(k-1,a(1,k),1,b(1),1)
    kp = kpvt(k)
    if (kp == k) goto 100

    !                 interchange.

    temp = b(k)
    b(k) = b(kp)
    b(kp) = temp
100 continue
110 continue
    k = k + 1
    goto 150
120 continue

    !           2 x 2 pivot block.

    if (k == 1) goto 140

    !              apply the transformation.

    b(k) = b(k) + ddot(k-1,a(1,k),1,b(1),1)
    b(k+1) = b(k+1) + ddot(k-1,a(1,k+1),1,b(1),1)
    kp = iabs(kpvt(k))
    if (kp == k) goto 130

    !                 interchange.

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
    ! endif
  end subroutine dsisl

  subroutine dsidi(a,lda,n,kpvt,det,inert,job)
    integer :: lda,n,job
    double precision :: a(lda,1),work(n)
    double precision :: det(2)
    integer :: kpvt(1),inert(3)

    !     dsidi computes the determinant, inertia and inverse
    !     of a double precision symmetric matrix using the factors from
    !     dsifa.

    !     on entry

    !        a       double precision(lda,n)
    !                the output from dsifa.

    !        lda     integer
    !                the leading dimension of the array a.

    !        n       integer
    !                the order of the matrix a.

    !        kpvt    integer(n)
    !                the pivot vector from dsifa.

    !        work    double precision(n)
    !                work vector.  contents destroyed.

    !        job     integer
    !                job has the decimal expansion  abc  where
    !                   if  c .ne. 0, the inverse is computed,
    !                   if  b .ne. 0, the determinant is computed,
    !                   if  a .ne. 0, the inertia is computed.

    !                for example, job = 111  gives all three.

    !     on return

    !        variables not requested by job are not used.

    !        a      contains the upper triangle of the inverse of
    !               the original matrix.  the strict lower triangle
    !               is never referenced.

    !        det    double precision(2)
    !               determinant of original matrix.
    !               determinant = det(1) * 10.0**det(2)
    !               with 1.0 .le. dabs(det(1)) .lt. 10.0
    !               or det(1) = 0.0.

    !        inert  integer(3)
    !               the inertia of the original matrix.
    !               inert(1)  =  number of positive eigenvalues.
    !               inert(2)  =  number of negative eigenvalues.
    !               inert(3)  =  number of zero eigenvalues.

    !     error condition

    !        a division by zero may occur if the inverse is requested
    !        and  dsico  has set rcond .eq. 0.0
    !        or  dsifa  has set  info .ne. 0 .

    !     linpack. this version dated 08/14/78 .
    !     james bunch, univ. calif. san diego, argonne nat. lab

    !     subroutines and functions

    !     blas daxpy,dcopy,ddot,dswap
    !     fortran dabs,iabs,mod

    !     internal variables.

    double precision :: akkp1,ddot,temp
    double precision :: ten,d,t,ak,akp1
    integer :: j,jb,k,km1,ks,kstep
    logical :: noinv,nodet,noert

    ! if CRAY
    !      call ssidi(a,lda,n,kpvt,det,inert,work,job)
    ! else
    noinv = mod(job,10) .eq. 0
    nodet = mod(job,100)/10 .eq. 0
    noert = mod(job,1000)/100 .eq. 0

    if (nodet .AND. noert) goto 140
    if (noert) goto 10
    inert(1) = 0
    inert(2) = 0
    inert(3) = 0
10  continue
    if (nodet) goto 20
    det(1) = 1.0d0
    det(2) = 0.0d0
    ten = 10.0d0
20  continue
    t = 0.0d0
    do  130  k = 1, n
       d = a(k,k)

       !           check if 1 by 1

       if (kpvt(k) > 0) goto 50

       !              2 by 2 block
       !              use det (d  s)  =  (d/t * c - t) * t  ,  t = dabs(s)
       !                      (s  c)
       !              to avoid underflow/overflow troubles.
       !              take two passes through scaling.  use  t  for flag.

       if (t /= 0.0d0) goto 30
       t = dabs(a(k,k+1))
       d = (d/t)*a(k+1,k+1) - t
       goto 40
30     continue
       d = t
       t = 0.0d0
40     continue
50     continue

       if (noert) goto 60
       if (d > 0.0d0) inert(1) = inert(1) + 1
       if (d < 0.0d0) inert(2) = inert(2) + 1
       if (d == 0.0d0) inert(3) = inert(3) + 1
60     continue

       if (nodet) cycle
       det(1) = d*det(1)
       if (det(1) == 0.0d0) goto 110
70     if (dabs(det(1)) >= 1.0d0) goto 80
       det(1) = ten*det(1)
       det(2) = det(2) - 1.0d0
       goto 70
80     continue
90     if (dabs(det(1)) < ten) goto 100
       det(1) = det(1)/ten
       det(2) = det(2) + 1.0d0
       goto 90
100    continue
110    continue

130 enddo
140 continue

    !     compute inverse(a)

    if (noinv) goto 270
    k = 1
150 if (k > n) goto 260
    km1 = k - 1
    if (kpvt(k) < 0) goto 180

    !              1 by 1

    a(k,k) = 1.0d0/a(k,k)
    if (km1 < 1) goto 170
    call dcopy(km1,a(1,k),1,work,1)
    do  160  j = 1, km1
       a(j,k) = ddot(j,a(1,j),1,work,1)
       call daxpy(j-1,work(j),a(1,j),1,a(1,k),1)
160 enddo
    a(k,k) = a(k,k) + ddot(km1,work,1,a(1,k),1)
170 continue
    kstep = 1
    goto 220
180 continue

    !              2 by 2

    t = dabs(a(k,k+1))
    ak = a(k,k)/t
    akp1 = a(k+1,k+1)/t
    akkp1 = a(k,k+1)/t
    d = t*(ak*akp1 - 1.0d0)
    a(k,k) = akp1/d
    a(k+1,k+1) = ak/d
    a(k,k+1) = -akkp1/d
    if (km1 < 1) goto 210
    call dcopy(km1,a(1,k+1),1,work,1)
    do  190  j = 1, km1
       a(j,k+1) = ddot(j,a(1,j),1,work,1)
       call daxpy(j-1,work(j),a(1,j),1,a(1,k+1),1)
190 enddo
    a(k+1,k+1) = a(k+1,k+1) + ddot(km1,work,1,a(1,k+1),1)
    a(k,k+1) = a(k,k+1) + ddot(km1,a(1,k),1,a(1,k+1),1)
    call dcopy(km1,a(1,k),1,work,1)
    do  200  j = 1, km1
       a(j,k) = ddot(j,a(1,j),1,work,1)
       call daxpy(j-1,work(j),a(1,j),1,a(1,k),1)
200 enddo
    a(k,k) = a(k,k) + ddot(km1,work,1,a(1,k),1)
210 continue
    kstep = 2
220 continue

    !           swap

    ks = iabs(kpvt(k))
    if (ks == k) goto 250
    call dswap(ks,a(1,ks),1,a(1,k),1)
    do  230  jb = ks, k
       j = k + ks - jb
       temp = a(j,k)
       a(j,k) = a(ks,j)
       a(ks,j) = temp
230 enddo
    if (kstep == 1) goto 240
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
    ! endif
  end subroutine dsidi

  subroutine dsifa(a,lda,n,kpvt,info)
    integer :: lda,n,kpvt(1),info
    double precision :: a(lda,1)

    !     dsifa factors a double precision symmetric matrix by elimination
    !     with symmetric pivoting.

    !     to solve  a*x = b , follow dsifa by dsisl.
    !     to compute  inverse(a)*c , follow dsifa by dsisl.
    !     to compute  determinant(a) , follow dsifa by dsidi.
    !     to compute  inertia(a) , follow dsifa by dsidi.
    !     to compute  inverse(a) , follow dsifa by dsidi.

    !     on entry

    !        a       double precision(lda,n)
    !                the symmetric matrix to be factored.
    !                only the diagonal and upper triangle are used.

    !        lda     integer
    !                the leading dimension of the array  a .

    !        n       integer
    !                the order of the matrix  a .

    !     on return

    !        a       a block diagonal matrix and the multipliers which
    !                were used to obtain it.
    !                the factorization can be written  a = u*d*trans(u)
    !                where  u  is a product of permutation and unit
    !                upper triangular matrices , trans(u) is the
    !                transpose of  u , and  d  is block diagonal
    !                with 1 by 1 and 2 by 2 blocks.

    !        kpvt    integer(n)
    !                an integer vector of pivot indices.

    !        info    integer
    !                = 0  normal value.
    !                = k  if the k-th pivot block is singular. this is
    !                     not an error condition for this subroutine,
    !                     but it does indicate that dsisl or dsidi may
    !                     divide by zero if called.

    !     linpack. this version dated 08/14/78 .
    !     james bunch, univ. calif. san diego, argonne nat. lab.

    !     subroutines and functions

    !     blas daxpy,dswap,idamax
    !     fortran dabs,dmax1,dsqrt

    !     internal variables

    double precision :: ak,akm1,bk,bkm1,denom,mulk,mulkm1,t
    double precision :: absakk,alpha,colmax,rowmax
    integer :: imax,imaxp1,j,jj,jmax,k,km1,km2,kstep,idamax
    logical :: swap

    !     initialize

    !     alpha is used in choosing pivot block size.
    ! if CRAY
    !      call ssifa(a,lda,n,kpvt,info)
    ! else
    alpha = (1.0d0 + dsqrt(17.0d0))/8.0d0

    info = 0

    !     main loop on k, which goes from n to 1.

    k = n
10  continue

    !        leave the loop if k=0 or k=1.

    !     ...exit
    if (k == 0) goto 200
    if (k > 1) goto 20
    kpvt(1) = 1
    if (a(1,1) == 0.0d0) info = 1
    !     ......exit
    goto 200
20  continue

    !        this section of code determines the kind of
    !        elimination to be performed.  when it is completed,
    !        kstep will be set to the size of the pivot block, and
    !        swap will be set to .true. if an interchange is
    !        required.

    km1 = k - 1
    absakk = dabs(a(k,k))

    !        determine the largest off-diagonal element in
    !        column k.

    imax = idamax(k-1,a(1,k),1)
    colmax = dabs(a(imax,k))
    if (absakk < alpha*colmax) goto 30
    kstep = 1
    swap = .false.
    goto 90
30  continue

    !           determine the largest off-diagonal element in
    !           row imax.

    rowmax = 0.0d0
    imaxp1 = imax + 1
    do  40  j = imaxp1, k
       rowmax = dmax1(rowmax,dabs(a(imax,j)))
40  enddo
    if (imax == 1) goto 50
    jmax = idamax(imax-1,a(1,imax),1)
    rowmax = dmax1(rowmax,dabs(a(jmax,imax)))
50  continue
    if (dabs(a(imax,imax)) < alpha*rowmax) goto 60
    kstep = 1
    swap = .true.
    goto 80
60  continue
    if (absakk < alpha*colmax*(colmax/rowmax)) goto 70
    kstep = 1
    swap = .false.
    goto 80
70  continue
    kstep = 2
    swap = imax .ne. km1
80  continue
90  continue
    if (dmax1(absakk,colmax) /= 0.0d0) goto 100

    !           column k is zero.  set info and iterate the loop.

    kpvt(k) = k
    info = k
    goto 190
100 continue
    if (kstep == 2) goto 140

    !           1 x 1 pivot block.

    if ( .NOT. swap) goto 120

    !              perform an interchange.

    call dswap(imax,a(1,imax),1,a(1,k),1)
    do  110  jj = imax, k
       j = k + imax - jj
       t = a(j,k)
       a(j,k) = a(imax,j)
       a(imax,j) = t
110 enddo
120 continue

    !           perform the elimination.

    do  130  jj = 1, km1
       j = k - jj
       mulk = -a(j,k)/a(k,k)
       t = mulk
       call daxpy(j,t,a(1,k),1,a(1,j),1)
       a(j,k) = mulk
130 enddo

    !           set the pivot array.

    kpvt(k) = k
    if (swap) kpvt(k) = imax
    goto 190
140 continue

    !           2 x 2 pivot block.

    if ( .NOT. swap) goto 160

    !              perform an interchange.

    call dswap(imax,a(1,imax),1,a(1,k-1),1)
    do  150  jj = imax, km1
       j = km1 + imax - jj
       t = a(j,k-1)
       a(j,k-1) = a(imax,j)
       a(imax,j) = t
150 enddo
    t = a(k-1,k)
    a(k-1,k) = a(imax,k)
    a(imax,k) = t
160 continue

    !           perform the elimination.

    km2 = k - 2
    if (km2 == 0) goto 180
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
170 enddo
180 continue

    !           set the pivot array.

    kpvt(k) = 1 - k
    if (swap) kpvt(k) = -imax
    kpvt(k-1) = kpvt(k)
190 continue
    k = k - kstep
    goto 10
200 continue
    return
    ! endif
  end subroutine dsifa

end module m_amix
