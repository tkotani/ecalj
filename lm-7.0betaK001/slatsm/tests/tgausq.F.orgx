C ... test gausq,mklegw
      subroutine fmain
      implicit none
      integer n,n0,i,nn,i1mach,mode
      parameter (n=20,n0=500)
      double precision xp(n0),wp(n0),pwr,x2,sum,sum1,x1,slope,d1mach
      logical cmdstr,a2bin
      character*72 first,fmt*20,f2(1)
      integer iarg
      equivalence (f2,first)

C --- use this branch to make quadrature for printout ---
C      mode = 0
C      iarg = 1
C    4 if (.not. cmdstr(iarg,first)) goto 8
C      if (first(1:3) .eq. '-m=') then
C        i = 3
C        if (.not. a2bin(first,mode,2,0,' ',i,-1)) goto 8
C        iarg = iarg+1
C        goto 4
C      endif
C      i = 0
C      if (.not. a2bin(first,nn,2,0,' ',i,-1)) goto 8
C      if (.not. cmdstr(iarg+1,first)) goto 8
C      i = 0
C      if (.not. a2bin(first,x1,4,0,' ',i,-1)) goto 8
C      if (.not. cmdstr(iarg+2,first)) goto 8
C      i = 0
C      if (.not. a2bin(first,x2,4,0,' ',i,-1)) goto 8
C      call awrit2('%% rows %i cols %i',' ',80,i1mach(2),nn,2)
C      call gausq(nn,x1,x2,xp,wp,mode,000)
C
C      i = -dlog(d1mach(3))/dlog(10d0)
C      call awrit2('%x(1x,2(1pe%i.%i))',fmt,len(fmt),0,i+10,i)
C      do  400  i = 1, nn
C        write(*,fmt) xp(i),wp(i)
C  410   format(1x,2(1pe26.15))
C  400 continue
C      call cexit(-1,1)
C    8 print *, 'usage: gausq [-m=#] n x1 x2'
C      print *, '              -m=# is mode (see gausq.f)'
C      call cexit(-1,1)

      print *, 'using',n,' point quadrature'
      call gausq(n,1d0,2d0,xp,wp,0,000)
      pwr = 2.2
      sum = 0d0
      do  10  i = 1, n
   10 sum = sum + wp(i)*xp(i)**pwr
      print 333, pwr, sum, (2d0**(pwr+1)-1d0)/(pwr+1)
  333 format('check integral x**p,     p=',f6.3,2f20.15)

      pwr = 0.5
      call gausq(n,1d0,2d0,xp,wp,2,000)
      sum = 0d0
      do  20  i = 1, n
   20 sum = sum + wp(i)*xp(i)**pwr
      print 334, pwr, sum, (2d0**(pwr+3)-1d0)/(pwr+3)
  334 format('check integral x**(p+2), p=',f6.3,2f20.15)

      pwr = 1
      x2 = 1.5
      call gausq(n,1d0,60d0/x2,xp,wp,0,000)
      sum1 = 0d0
      do  30  i = 1, n
   30 sum1 = sum1 + wp(i)*xp(i)**pwr*exp(-xp(i)/x2)
      call gausq(n,1d0,x2*5,xp,wp,10,080)
      sum = 0d0
      do  35  i = 1, n
   35 sum = sum + wp(i)*xp(i)**pwr*exp(-xp(i)/x2)

      print 335, pwr, x2, sum1, sum, dlog10(abs(sum1/sum-1))
  335 format('check integral x**(p+0)exp(-x/x2), p,x2=',2f6.3,':'/
     .  '  from straight quad:  ',f20.15/
     .  '  from exp(-x/x2) quad:',f20.15,'  log10(ratio-1)',f8.3)


      x1 = 4
      pwr = 7
      x2 = 1.5
      call gausq(n,x1,x1+1.3d0*x2*dsqrt(-dlog(1d-20)),
     .  xp,wp,0,000)
      sum1 = 0d0
      do  40  i = 1, n
   40 sum1 = sum1 + wp(i)*xp(i)**pwr*exp(-(xp(i)/x2)**2)
      slope = 2*x1/(x2*x2)
      call gausq(n,x1,5/slope,xp,wp,10,000)
      sum = 0d0
      do  45  i = 1, n
   45 sum = sum + wp(i)*xp(i)**pwr*exp(-(xp(i)/x2)**2)

      print 336, pwr, x2, sum1, sum, dlog10(abs(sum1/sum-1))
  336 format('check integral x**(p+0)exp[(-x/x2)^2], p,x2=',2f6.3,':'/
     .  '  from straight quad:  ',f20.15/
     .  '  from exp(-x/x2) quad:',f20.15,'  log10(ratio-1)',f8.3)


      pwr = 1.5
      x2 = 1.5
      call gausq(n,1d0,1.3d0*x2*dsqrt(-dlog(1d-20)),
     .  xp,wp,0,000)
      sum1 = 0d0
      do  50  i = 1, n
   50 sum1 = sum1 + wp(i)*xp(i)**(pwr+2)*exp(-(xp(i)/x2)**2)
      call gausq(n,1d0,2*x2,xp,wp,12,000)
      sum = 0d0
      do  55  i = 1, n
   55 sum = sum + wp(i)*xp(i)**pwr*exp(-(xp(i)/x2)**2)

      print 337, pwr, x2, sum1, sum, dlog10(abs(sum1/sum-1))
  337 format('check integral x**(p+2)exp[(-x/x2)^2], p,x2=',2f6.3,':'/
     .  '  from straight quad:  ',f20.15/
     .  '  from exp(-x/x2) quad:',f20.15,'  log10(ratio-1)',f8.3)



      end
