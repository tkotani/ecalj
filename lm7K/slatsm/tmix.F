      integer function tmix(nelts,msg,nmix,mmix,ido,beta,ipr,tm,
     .  norm,kpvt,wk,t,rmsdel)
C- Accelerated mixing of a vector
C ----------------------------------------------------------------
Ci Inputs
Ci   nmix: +/- number of previous iterations to fold into mix
Ci         nmix = 0 => linear mixing (x* = x0)
Ci         For meaning of nmix<0, see Remarks
Ci   mmix: maximum number of previous iterations to fold into mix
Ci         (For dimensioning array wk)
Ci   nelts:number of elements to mix
Ci   wk:   array of dimension (nelts,2+mmix,2) where:
Ci         wk(*,0,1) holds f(xi) (see remarks)
Ci         wk(*,i,1) holds  d(i) (see remarks) ,  i>0
Ci         wk(*,i,2) holds   xi  (see remarks)
Ci   beta: new x is beta f(x*) + (1-beta) x* (linear mixing only)
Ci   tm:   upper limit to any tj:  if any tj exceeds tm, effective
Ci         nmix is decremented.
Ci   ido:  0: normal mixing; 1: calculate tj only; 2: mix with input tj
Ci         10: silly mixing (mix each element of vector independently)
Co Outputs
Co   f(x_i) => f(x_i+1); x_i => x_i+1
Co   wk(*,0,1:2): x* and f*
Co   rmsdel:rms difference between x_0 and f(x_0)
Co   tmix:  returns effective nmix (see input tm)
Cr Remarks
Cr   tmix minimizes  the dot product ( sum_j t_j (f_j - x_j) )^2
Cr   subject to the constraint sum_j t_j = 1.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nelts,mmix,nmix,kpvt(1),ipr
      double precision norm(0:mmix+1,0:mmix+1),wk(nelts,0:mmix+1,2),
     .  t(0:mmix+1),beta,tm,rmsdel
      character msg*(*)
C Local parameters
      integer i,j,nwk,inert(3),kelt,nmake,ido,iprint,is
      double precision dsqrt,ddot,det(2),sumsqr,sing,sinmax
      parameter (sinmax=1000d0)
      character s*100

C --- Setup ---
      nwk = nelts*(mmix+2)
      kelt = 0
C ... nmake is the number of elements mixed per matrix inversion
      nmake = nelts
      if (ido/10 .eq. 1) nmake = 1

C --- d_0 = f-x  =>  wk(*,0,1) ---
      call daxpy(nelts,-1d0,wk(1,0,2),1,wk,1)
C ... Copy x_0 and d_0 to end of wk:  x*, d* constructed there
      call dmcpy(wk,nwk,1,wk(1,mmix+1,1),nwk,1,nelts,2)

C --- Obtain the tj ---
   11 continue
      kelt = kelt+1
      if (ido .eq. 2) goto 40

C ... Normal equations Make < (d_j d_k) > (+ last column for const)
   10 continue
      t(0) = 1
      if (nmix .eq. 0) goto 40
      if (nmix .lt. 0) call rx('tmix: bad nmix')

C ... Regular branch
      if (ido/10 .eq. 0) then
        sumsqr = 0
        do  20  i = 0, nmix
          t(i) = 0
          norm(i,nmix+1) = 1
          norm(nmix+1,i) = 1
          do  20  j = 0, nmix
            norm(i,j) =  ddot(nelts,wk(1,i,1),1,wk(1,j,1),1)
            sumsqr = sumsqr + norm(i,j)**2
   20     continue
          sumsqr = dsqrt(sumsqr)/(nmix+1)
C ... Silly branch (not implemented)
        elseif (ido/10 .eq. 1) then
          stop 'not implemented'
          do  120  i = 1, nmix
            t(i) = wk(kelt,0,1)*wk(kelt,0,1) - wk(kelt,0,1)*wk(kelt,i,1)
            do  120  j = 1, nmix
              norm(i,j) =
     .          wk(kelt,0,1)*wk(kelt,0,1) - wk(kelt,0,1)*wk(kelt,j,1)
     .        - wk(kelt,i,1)*wk(kelt,0,1) + wk(kelt,i,1)*wk(kelt,j,1)
  120     continue
        endif
        t(nmix+1) = 1
        norm(nmix+1,nmix+1) = 0

C ...  Solve the simultaneous equations for tj
      call dsifa(norm,mmix+2,nmix+2,kpvt,i)
      if (i .ne. 0) then
        sing = sinmax + 1
      else
        call dsisl(norm,mmix+2,nmix+2,kpvt,t)
        call dsidi(norm,mmix+2,nmix+2,kpvt,det,inert,kpvt,10)
        sing = dabs(sumsqr/(det(1)*10**det(2)))
      endif

C --  Handle singular normal matrix --
      if (sing .gt. sinmax) then
        t(nmix) = 0
        nmix = nmix-1
        if (ipr .ge. 20) print 337, sinmax, nmix
  337   format(' TMIX: condition of normal eqns >',f6.0,
     .         ' Reducing nmix to',i2)
        goto 10
      endif

C --  Reduce nmix if any t_j exceeds tm --
      do  30  j = 1, nmix
        if (dabs(t(j)) .le. dabs(tm)) goto 30
        if (ipr .ge. 20) print 338,  nmix-1, (t(i), i=1, nmix)
  338   format(' TMIX: Reducing nmix to',i3,
     .         ': t_j exceeds tm: tj=',7(f8.5,2x))
        t(nmix) = 0
        nmix = nmix-1
        goto 10
   30 continue

C --- Exit before mixing if ido=1 ---
   40 continue
      tmix = nmix
      if (ido .eq. 1) then
C ...    Restore f = d_0 + x  =>  wk(*,0,1) when ido=1
        call daxpy(nelts,1d0,wk(1,0,2),1,wk,1)
        return
      endif

C --- Make (d,x)* = \sum_j t_j ((d,x)_j  ---
      call dscal(nmake,t(0),wk(kelt,mmix+1,1),1)
      call dscal(nmake,t(0),wk(kelt,mmix+1,2),1)
      do  45  j = 1, nmix
        call daxpy(nmake, t(j),wk(kelt,j,1),1,wk(kelt,mmix+1,1),1)
        call daxpy(nmake, t(j),wk(kelt,j,2),1,wk(kelt,mmix+1,2),1)
   45 continue

C ... Do next element for silly case ...
      if (ido/10 .eq. 1) then
        if (nmix .gt. 0 .and. iprint() .ge. 40)
     .  write(*,135) kelt, (t(j), j=1,nmix)
  135   format(i4,3x,'tj:',7(f8.5,2x))
        if (kelt .lt. nelts) goto 11
      endif

C --- Copy arrays to new positions ---
      do  50  i = mmix, 1, -1
   50 call dmcpy(wk(1,i-1,1),nwk,1,wk(1,i,1),nwk,1,nelts,2)

C ...  Calculate rms change
      rmsdel = dsqrt(ddot(nelts,wk,1,wk,1)/nelts)

C --  x* + beta d*  (or x + beta d* if npmix<0) --
      call daxpy(nelts,beta,wk(1,mmix+1,1),1,wk(kelt,mmix+1,2),1)

C --- Printout ---
      if (ipr .lt. 30) goto 60
      s = 'TMIX ' // msg
      call skpblb(s,80,is)
      is = is+1
      call bin2a(':  nmix=',0,0,0,1,0,80,s,is)
      call bin2a('(i5)',0,0,nmix,2,0,80,s,is)
      call bin2a('nelts=',2,0,0,1,0,80,s,is)
      call bin2a('(i5)',0,0,nelts,2,0,80,s,is)
      call bin2a('beta=',2,0,0,1,0,80,s,is)
      call bin2a('(f10.5)',0,0,beta,4,0,80,s,is)
      if (nmix .ne. 0) then
        call bin2a('tm=',2,0,0,1,0,80,s,is)
        call bin2a('(f8.2)',0,0,tm,4,0,80,s,is)
      endif
      call bin2a('rmsdel=',2,0,0,1,0,80,s,is)
      call bin2a('(1pd8.2)',0,0,rmsdel,4,0,80,s,is)
      print 147, s(1:is)
  147 format(1x,a)
      if (nmix .gt. 0) print 134, (t(j), j=0,nmix)
  134 format(3x,'tj:',7(f9.5,1x))
      if (ipr .lt. 61 .and. (ipr .lt. 41 .or. nelts .gt. 100)) goto 60
      write(*,110)
      do  12  i = 1, nelts
        if (dabs(wk(i,0,1)) + dabs(wk(i,mmix+1,2)-wk(i,0,2)) .ge. 5.d-7)
     .  write(*,111) i,wk(i,0,2),wk(i,0,2)+wk(i,0,1),
     .                 wk(i,0,1),wk(i,mmix+1,2)
   12 continue

C --  Restore  x* + beta d*  and  d*  from end of wk --
   60 call dmcpy(wk(1,mmix+1,1),nwk,1,wk,nwk,1,nelts,2)

  104 format(1p,4d18.11)
  111 format(i5,5f13.6)
  110 format(12x,'x',12X,'F(x)',9x,'Diff',9x,'Mixed')

      end


