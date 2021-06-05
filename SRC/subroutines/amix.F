      integer function amix(nelts,npmix,mmix,ido,beta,ipr,tm,norm,kpvt,
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
      integer nelts,mmix,npmix,kpvt(1),ipr
      double precision norm(mmix,mmix),wk(nelts,0:mmix+1,2),t(mmix),
     .  beta,tm,rmsdel
C Local parameters
      integer i,j,nwk,inert(3),nmix,kelt,nmake,ido,iprint,i1mach
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
        call dsidi(norm,mmix,nmix,kpvt,det,inert,kpvt,10)
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
c      call awrit6(' AMIX: nmix=%i mmix=%i  nelts=%i  beta=%1;6d  tm='//
c     .  '%1;6d  rmsdel=%1;3e',' ',80,i1mach(2),
c     .  nmix,mmix,nelts,beta,tm,rmsdel)
      write(*,133) nmix,mmix,nelts,beta,tm,rmsdel
  133 format(' AMIX: nmix=',i1,' mmix=',i1,'  nelts=',i6,
     .       '  beta=',f7.5,'  tm=',f8.5,'  rmsdel=',1pd8.2)
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

