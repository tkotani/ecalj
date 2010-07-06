C Tests verfc = vectorized derfc
C for quad precision, convert 'd0' to 'q0'
      subroutine fmain
      implicit none
      integer n,i,np,j
      parameter (n=500,np=60)
      double precision x(n),wk(n),f(n),fn(n),rerrmx,errmx,rmse,exact
     .  ,derfc,y0,xp(100),wp(100),pi,xx,res

      y0 = 1/sqrt(16*atan(1d0))
      do  10  i = 1, n
      x(i) = 6*(dble(i)/n) - 1
      x(i) = 10*(dble(i-1)/n)
   10 wk(i) = exp(-x(i)**2)*y0

C ... Obtain error function by numerical integration, to check derfc
C      pi = 4*atan(1d0)
C      do  40  i = 1, n, 10
C        if (x(i) .ge. 0d0 .and. x(i) .lt. 6d0) then
C          call gausq(np,x(i),x(n),xp,wp,000,000)
C          exact = 0
C          do  42  j = 1, np
C   42     exact = exact + wp(j)*exp(-xp(j)*xp(j))
C          fn(i) = 2/sqrt(pi)*exact
CC           print 339, x(i), fn(i), derfc(x(i))
CC  339      format(3f38.32)
C           print 334, x(i), derfc(x(i)), (fn(i)-derfc(x(i)))/fn(i)
C  334      format('x=',f18.8,' erfc =',f38.32,'  rel diff = ',1pe9.2)
C        endif
C   40 continue

      call cpudel(6,'start',errmx)
      call verfc(n,x,1,wk,f)
      call cpudel(6,'done',errmx)
c     pause
      rmse = 0d0
      rerrmx = 0d0
      errmx = 0d0
      print 331
  331 format('         x            f            rerr*1d16',
     .  '  aerr*1d18')
      do  20  i = 1, n
        exact = derfc(x(i))
        if (abs(f(i)/exact-1) .gt. rerrmx) rerrmx = abs(f(i)/exact-1)
        if (abs(f(i)-exact) .gt. errmx) errmx = abs(f(i)-exact)
        rmse = rmse + abs(f(i)/exact-1)**2
        print 332, i, x(i), f(i), (f(i)/exact-1), f(i)-exact
  332   format(i3,f12.6,f15.8,16p,f14.4,18p,f12.4)
   20 continue
      print 333, rerrmx,errmx,sqrt(rmse/n)
  333 format('1e10*rerrmx,errmx,rmse',10p,3f15.8)
      call cpudel(6,'exit',errmx)

   60 print *, 'enter x'
      read(*,*) xx
      if (xx .eq. 0) stop
      call verfc(1,xx,1,wk,res)
      exact = derfc(xx)
      print 332, i, xx, res, res/exact-1, res-exact
      goto 60

      end
C
C      subroutine verfc(n,x,mode,wk,f)
CC- Evaluates complement of error function for a vector of points
CC  as wk(x)*rational function, with wk(x)=exp(-x*x)/sqrt(4*pi).
CC ----------------------------------------------------------------
CCi Inputs
CCi   x,n:  vector of points, and number
CCi   mode: 0, wk is input (see Remarks); 1, wk generated internally
CCo Outputs
CCo   f = erfc(x) for n points
CCr Remarks
CCr   Relative errors in verfc (compiled on an SGI R8000) 
CCr   when compiled with IEEE quadruple-precision arithmetic,
CCr   < 1e-20 for x<2; <1e-18 for x<4; -> 2e-12 as x->10
CCr   when compiled with IEEE double-precision arithmetic
CCr   machine precision for x<5; -> 2e-12 as x->10
CC ----------------------------------------------------------------
C      implicit none
C      integer n,i,mode
C      double precision x(1),wk(1),f(1),z,w
C
C      double precision t(2,0:7),b(2,0:8),f1,f2,y0,half
C      parameter (half=.5d0)
C
C      f1(w) = wk(i)*
C     .  (((((((t(1,7)*w  + t(1,6))*w + t(1,5))*w + t(1,4))*w +
C     .         t(1,3))*w + t(1,2))*w + t(1,1))*w + t(1,0)) /
C     . ((((((((b(1,8)*w  + b(1,7))*w + b(1,6))*w + b(1,5))*w +
C     .         b(1,4))*w + b(1,3))*w + b(1,2))*w + b(1,1))*w + 1)
C      f2(w) = wk(i)*
C     .  (((((((t(2,7)*w  + t(2,6))*w + t(2,5))*w + t(2,4))*w +
C     .         t(2,3))*w + t(2,2))*w + t(2,1))*w + t(2,0)) /
C     . ((((((((b(2,8)*w  + b(2,7))*w + b(2,6))*w + b(2,5))*w +
C     .         b(2,4))*w + b(2,3))*w + b(2,2))*w + b(2,1))*w + 1)
C
CC Input data generated for x<1.3 and x>1.3 by
CC gausq -m=0 80 0 1.3 |  mc-16 -'f(f30.25,g30.20,g30.20)' . -e3 x1 \
CC 'erfc(x1)' x2 -e3 x1 x2 1 > dat.erfc1
CC gausq -m=10 80 1.3 .9 | mc-16 -'f(f30.25,g30.20,g30.20)' . -e3 x1 \
CC 'erfc(x1)' x2 -e3 x1 x2 x3'*exp(-7.5*x1)*1e7' > dat.erfc1
CC lfit-16 -pr30 '-vy0=1/sqrt(4*pi)' -y='x2/y0*exp(x1*x1)' -w 16 \
CC "z1=1 z2=x z3=z2*x z4=z3*x z5=z4*x z6=z5*x z7=z6*x z8=z7*x z9=-y*x\
CC z10=z9*x z11=z10*x z12=z11*x z13=z12*x z14=z13*x z15=z14*x z16=z15*x"\
CC dat.erfc1
CC      data t(1,0) /3.54490770181103205460333 d0 /
CC      data t(1,1) /5.35879488908299920197094 d0 /
CC      data t(1,2) /4.13171494834170603985303 d0 /
CC      data t(1,3) /1.95545653031736723219181 d0 /
CC      data t(1,4) /0.60393598762715969705170 d0 /
CC      data t(1,5) /0.12041672912546619389560 d0 /
CC      data t(1,6) /0.01431341289801243387620 d0 /
CC      data t(1,7) /0.00078269387282409683620 d0 /
CC      data b(1,1) /2.64006729549030365675172 d0 /
CC      data b(1,2) /3.14453220398292447912212 d0 /
CC      data b(1,3) /2.21203421455139074620401 d0 /
CC      data b(1,4) /1.00784627460138212754970 d0 /
CC      data b(1,5) /0.30551800710984891854570 d0 /
CC      data b(1,6) /0.06040691635637395141680 d0 /
CC      data b(1,7) /0.00715650999266887712620 d0 /
CC      data b(1,8) /0.00039135372833175242650 d0 /
CC      data t(2,0) /3.5449076592709693083486 d0 /
CC      data t(2,1) /6.5367552783003810899477 d0 /
CC      data t(2,2) /5.9459069263430207040377 d0 /
CC      data t(2,3) /3.3494256255191946772255 d0 /
CC      data t(2,4) /1.2500554382038220883901 d0 /
CC      data t(2,5) /0.3092338440247459558100 d0 /
CC      data t(2,6) /0.0474126561358463066805 d0 /
CC      data t(2,7) /0.0035794893127754589293 d0 /
CC      data b(2,1) /2.9723636627785400295532 d0 /
CC      data b(2,2) /4.0312632258660319839958 d0 /
CC      data b(2,3) /3.2735360927940406983529 d0 /
CC      data b(2,4) /1.7511343462244764429298 d0 /
CC      data b(2,5) /0.6368802721279892901290 d0 /
CC      data b(2,6) /0.1555118277486855990202 d0 /
CC      data b(2,7) /0.0237063269110756317692 d0 /
CC      data b(2,8) /0.0017897446754386780195 d0 /
C
CC gausq -m=0 80 0 1.3 | mc-16 -'f(f30.25,g30.20,g30.20)' . -e3 x1-.5\
CC 'erfc(x1)' x2 -e3 x1 x2 1 > dat.erfc1
CC gausq -m=10 80 1.3 .9 | mc-16 -'f(f30.25,g30.20,g30.20)' . -e3 x1-2\
CC 'erfc(x1)' x2 -e3 x1 x2 x3'*exp(-9.5*x1)*1' > dat.erfc1
CC lfit-16 -pr30 '-vy0=1/sqrt(4*pi)' -y='x2/y0*exp((x1+1/2)^2)' -w 16 \
CC "z1=1 z2=x z3=z2*x z4=z3*x z5=z4*x z6=z5*x z7=z6*x z8=z7*x z9=-y*x\
CC z10=z9*x z11=z10*x z12=z11*x z13=z12*x z14=z13*x z15=z14*x z16=z15*x"\
CC dat.erfc1
C      data t(1,0) /2.1825654430601881683921 d0 /
C      data t(1,1) /3.2797163457851352620353 d0 /
C      data t(1,2) /2.3678974393517268408614 d0 /
C      data t(1,3) /1.0222913982946317204515 d0 /
C      data t(1,4) /0.2817492708611548747612 d0 /
C      data t(1,5) /0.0492163291970253213966 d0 /
C      data t(1,6) /0.0050315073901668658074 d0 /
C      data t(1,7) /0.0002319885125597910477 d0 /
C      data b(1,1) /2.3353943034936909280688 d0 /
C      data b(1,2) /2.4459635806045533260353 d0 /
C      data b(1,3) /1.5026992116669133262175 d0 /
C      data b(1,4) /0.5932558960613456039575 d0 /
C      data b(1,5) /0.1544018948749476305338 d0 /
C      data b(1,6) /0.0259246506506122312604 d0 /
C      data b(1,7) /0.0025737049320207806669 d0 /
C      data b(1,8) /0.0001159960791581844571 d0 /
C
C
C      data t(2,0) /0.9053540999623491587309 d0 /
C      data t(2,1) /1.3102485359407940304963 d0 /
C      data t(2,2) /0.8466279145104747208234 d0 /
C      data t(2,3) /0.3152433877065164584097 d0 /
C      data t(2,4) /0.0729025653904144545406 d0 /
C      data t(2,5) /0.0104619982582951874111 d0 /
C      data t(2,6) /0.0008626481680894703936 d0 /
C      data t(2,7) /0.0000315486913658202140 d0 /
C      data b(2,1) /1.8653829878957091311190 d0 /
C      data b(2,2) /1.5514862329833089585936 d0 /
C      data b(2,3) /0.7521828681511442158359 d0 /
C      data b(2,4) /0.2327321308351101798032 d0 /
C      data b(2,5) /0.0471131656874722813102 d0 /
C      data b(2,6) /0.0061015346650271900230 d0 /
C      data b(2,7) /0.0004628727666611496482 d0 /
C      data b(2,8) /0.0000157743458828120915 d0 /
C
C      if (mode .eq. 1) then
C        y0 = 1/sqrt(16*atan(1d0))
C        do  5  i = 1, n
C    5   wk(i) = exp(-x(i)**2)*y0
C      endif
C
C      do  10  i = 1, n
C        z = x(i)
C        if (z .gt. 1.3d0) then
C          f(i) = f2(z-2d0)
C        elseif (z .gt. 0) then
C          f(i) = f1(z-half)
C        elseif (z .gt. -1.3d0) then
C          f(i) = 2d0 - f1(-z-half)
C        else
C          f(i) = 2d0 - f2(-z-2d0)
C        endif
C   10 continue
C      end
