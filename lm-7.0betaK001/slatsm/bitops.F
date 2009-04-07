C      integer function bitand(a,b)
CC- Returns the bitwise AND of two integers
C      implicit none
C      integer a,b
C      integer TF(64),i,n,i1,i2
C      if (a .lt. 0) then
C        bitand = b
C      elseif (b .lt. 0) then
C        bitand = a
C      else
C        i1 = a
C        i2 = b
C        n = 0
C        do  10  i = 1, 64
C          if (i1 .eq. 0 .or. i2 .eq. 0) goto 12
C          TF(i) = (mod(i1,2)+mod(i2,2))/2
C          i1 = i1/2
C          i2 = i2/2
C          n = n+1
C   10   continue
C   12   continue
C        bitand = 0
C        do  20  i = n, 1, -1
C   20   bitand = bitand*2 + tf(i)
C      endif
C      end
C      integer function bitor(a,b)
CC- Returns the bitwise OR of two integers
C      implicit none
C      integer a,b
C      integer TF(64),i,n,i1,i2
C      if (a .lt. 0) then
C        bitor = a
C      elseif (b .lt. 0) then
C        bitor = b
C      else
C        i1 = a
C        i2 = b
C        n = 0
C        do  10  i = 1, 64
C          if (i1 .eq. 0 .and. i2 .eq. 0) goto 12
C          TF(i) = max(mod(i1,2),mod(i2,2))
C          i1 = i1/2
C          i2 = i2/2
C          n = n+1
C   10   continue
C   12   continue
C*       print *, (tf(i), i=1,n)
C        bitor = 0
C        do  20  i = n, 1, -1
C   20   bitor = bitor*2 + tf(i)
C      endif
C      end
      integer function bitlow(a)
C- Returns an 2**n where n is the lowest nonzero bit in a number
C     implicit none
      integer a,t2n,i1
      bitlow = 0
      if (a .eq. 0) return
      i1 = a
      t2n = 1
C ... Find lowest nonzero bit of mask, corresponding ix1
   10 continue
      if (mod(i1,2) .eq. 0) then
        t2n = t2n*2
        i1 = i1/2
        goto 10
      endif
      bitlow = t2n

      end

      logical function bittst(n,bit)
C- Returns true when a bit is set in an integer
C ----------------------------------------------------------------
Ci Inputs
Ci   n: integer
Ci   bit: a bit, ie 1,2,4,8,16, etc
Co Outputs
Co   bittst: true when bit in n is set, false otherwise
C ----------------------------------------------------------------
C     implicit none
      integer n,bit
      bittst = (mod(n,bit+bit) - mod(n,bit) .eq. bit)
      end
      integer function getdig(n,i,base)
C- Extracts one digit from an integer
C ----------------------------------------------------------------
Ci Inputs
Ci   n,i,base
Co Outputs
Co   getdig = ith digit from n, base "base"; eg 4=getdig(12345,1,10)
C ----------------------------------------------------------------
C     implicit none
      integer n,i,base
      getdig = mod(n/base**i,base)
      end
      double precision function dtrunc(x,n)
C- Truncates double precision x after n digits
      integer n,is
      double precision x,s

#if QUAD
      s = qlog10(dble(abs(x)))
#else
      s = dlog10(dabs(x))
#endif
      is = s
      if (is .gt. 0) then
        is = -int(is) + n-1
      else
        is = int(-is) + n
      endif
      s = 1
      do  10  i = 1, iabs(is)
   10 s = 10*s
#if QUAD
      if (is .gt. 0) dtrunc = int(s*dble(x))/s
      if (is .lt. 0) dtrunc = int(dble(x)/s)*s
#else
      if (is .gt. 0) dtrunc = int(s*x)/s
      if (is .lt. 0) dtrunc = int(x/s)*s
#endif
      end
      double precision function dround(x,n)
C- Rounds double precision x after n digits
C     implicit none
      integer n,is,i
      double precision x,s
#if QUAD
      double precision xnint
#endif

#if QUAD
      s = qlog10(dble(abs(x)))
#else
      s = dlog10(dabs(x))
#endif
      is = s
      if (is .gt. 0) then
        is = -int(is) + n-1
      else
        is = int(-is) + n
      endif
      s = 1
      do  10  i = 1, iabs(is)
   10 s = 10*s

#if QUAD
      if (is .gt. 0) dround = xnint(s*dble(x))/s
      if (is .lt. 0) dround = xnint(dble(x)/s)*s
#else
      if (is .gt. 0) dround = dnint(s*x)/s
      if (is .lt. 0) dround = dnint(x/s)*s
#endif
      end
C ... Needed because SGI has no working nint function, quad precision
      double precision function xnint(x)
C     implicit none
      double precision x,xx

      xx = x
      xnint = nint(xx)
      end



