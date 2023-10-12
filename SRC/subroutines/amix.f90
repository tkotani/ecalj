!> Anderson mixing of a vector
module m_amix !nothing saved in this module. Bundle linpack subroutines/functions for amix
  public:: amix
contains
  integer function amix(nelts,npmix,mmix,ido,beta,ipr,tm,  wk,t,rmsdel) 
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
    sumsqr=-1d99
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


end module m_amix
