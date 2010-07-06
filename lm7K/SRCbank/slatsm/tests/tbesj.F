C test dbesj0,dbesj1 by comparing to Numerical Recipes
C test dbesj by increasing accuracy in recursion
      subroutine fmain
      implicit none
      double precision x,f1,f2,dbesj0,bessj0,dbesj1,bessj1,bessj,dbesj
      print *, 'testing dbesj0 against Numerical Recipes:'
      do  10  x = 0d0, 30d0, .1d0
        f1 = dbesj0(x)
        f2 = bessj0(x)
        print *, x,f1,f1-f2
   10 continue
      print *, 'testing dbesj1 against Numerical Recipes:'
      do  20  x = 0d0, 30d0, .1d0
        f1 = dbesj1(x)
        f2 = bessj1(x)
        print *, x,f1,f1-f2
   20 continue
      print *, 'testing dbesj against greater accuracy'
      do  30  x = 0d0, 30d0, .1d0
        f1 = dbesj(4,x)
        f2 = bessj(4,x)
        print *, x,f1,f1-f2
   30 continue
      end
C ... these are taken from numerical recipes
      double precision function bessj0(x)
      implicit none
      double precision y,p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6
     .  ,s1,s2,s3,s4,s5,s6,x,ax,z,xx
      data p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,
     *  -.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-
     *  1,.1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      data r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,651619640.7d
     *  0,-11214424.18d0,77392.33017d0,-184.9052456d0/,
     *  s1,s2,s3,s4,s5,s6/57568490411.d0,1029532985.d0,
     *  9494680.718d0,59272.64853d0,267.8532712d0,1.d0/
      if(dabs(x).lt.8.d0)then
        y=x**2
        bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))
     *    /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
      else
        ax=dabs(x)
        z=8.d0/ax
        y=z**2
        xx=ax-.785398164d0
        bessj0=dsqrt(.636619772d0/ax)*(dcos(xx)*(p1+y*(p2+y*(p3+y*(p4+y
     *    *p5))))-z*dsin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      end
      double precision function bessj1(x)
      implicit none
      double precision y,p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6
     .  ,s1,s2,s3,s4,s5,s6,x,ax,z,xx
      data r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,242396853.1d0
     *  ,-2972611.439d0,15704.48260d0,-30.16036606d0/,
     *  s1,s2,s3,s4,s5,s6/144725228442.d0,2300535178.d0,
     *  18583304.74d0,99447.43394d0,376.9991397d0,1.d0/
      data p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,.2457520174d-5
     *  ,
     *  -.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0,-.2002690873d-3
     *  ,
     *  .8449199096d-5,-.88228987d-6,.105787412d-6/
      if(dabs(x).lt.8.d0)then
        y=x**2
        bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))
     *    /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
      else
        ax=dabs(x)
        z=8.d0/ax
        y=z**2
        xx=ax-2.356194491d0
        bessj1=dsqrt(.636619772d0/ax)*(dcos(xx)*(p1+y*(p2+y*(p3+y*(p4+y
     *    *p5))))-z*dsin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
     *    *dsign(1.d0,x)
      endif
      end
C      double precision function dbesj(n,x)
C      implicit none
C      integer iacc,n,j,m,jsum
C      double precision bigno,bigni,tox,x,bjm,dbesj0,bj,dbesj1,bjp,sum
C      parameter (iacc=100,bigno=1d10,bigni=1d-10)
C
CC --- Handle j=0, 1 explicitly ---
C      if (n .lt. 0) call rx('DBESJ: n lt 0')
C      if (n .eq. 0) then 
C        dbesj = dbesj0(x)
C      elseif (n .eq. 1) then 
C        dbesj = dbesj1(x)
C      endif
C
CC --- Bessel functions for larger n ---
C      tox = 2d0/x
C      if (x .gt. dble(n)) then
C        bjm = dbesj0(x)
C        bj = dbesj1(x)
C        do  11  j = 1, n-1
C          bjp = j*tox*bj-bjm
C          bjm = bj
C          bj = bjp
C   11   continue
C        dbesj = bj
C      else
C        m = 2*((n+int(dsqrt(dble(iacc*n))))/2)
C        dbesj = 0d0
C        jsum = 0
C        sum = 0d0
C        bjp = 0d0
C        bj = 1d0
C        do  12  j = m, 1,-1
C          bjm = j*tox*bj-bjp
C          bjp = bj
C          bj = bjm
C          if (dabs(bj) .gt. bigno) then
C            bj = bj*bigni
C            bjp = bjp*bigni
C            dbesj = dbesj*bigni
C            sum = sum*bigni
C          endif
C          if (jsum .ne. 0)sum = sum+bj
C          jsum = 1-jsum
C          if (j .eq. n) dbesj = bjp
C   12   continue
C        sum = 2d0*sum-bj
C        dbesj = dbesj/sum
C      endif
C      end

      double precision function bessj(n,x)
      implicit none
      integer iacc,n,j,m,jsum
      double precision bigno,bigni,tox,x,bjm,dbesj0,bj,dbesj1,bjp,sum
      parameter (iacc=150,bigno=1d12,bigni=1d-12)

C --- Handle j=0, 1 explicitly ---
      if (n .lt. 0) call rx('BESSJ: n lt 0')
      if (n .eq. 0) then 
        bessj = dbesj0(x)
      elseif (n .eq. 1) then 
        bessj = dbesj1(x)
      endif

C --- Bessel functions for larger n ---
      tox = 2d0/x
      if (x .gt. dble(n)) then
        bjm = dbesj0(x)
        bj = dbesj1(x)
        do  11  j = 1, n-1
          bjp = j*tox*bj-bjm
          bjm = bj
          bj = bjp
   11   continue
        bessj = bj
      else
        m = 2*((n+int(dsqrt(dble(iacc*n))))/2)
        bessj = 0d0
        jsum = 0
        sum = 0d0
        bjp = 0d0
        bj = 1d0
        do  12  j = m, 1,-1
          bjm = j*tox*bj-bjp
          bjp = bj
          bj = bjm
          if (dabs(bj) .gt. bigno) then
            bj = bj*bigni
            bjp = bjp*bigni
            bessj = bessj*bigni
            sum = sum*bigni
          endif
          if (jsum .ne. 0)sum = sum+bj
          jsum = 1-jsum
          if (j .eq. n) bessj = bjp
   12   continue
        sum = 2d0*sum-bj
        bessj = bessj/sum
      endif
      end

