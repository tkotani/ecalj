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
      implicit none
      integer n,i,base
      getdig = mod(n/base**i,base)
      end
      double precision function dtrunc(x,n)
C- Truncates double precision x after n digits
      integer n,is,i
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
      do    i = 1, iabs(is)
         s = 10*s
      enddo   
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
      do    i = 1, iabs(is)
         s = 10*s
      enddo   
#if QUAD
      if (is .gt. 0) dround = xnint(s*dble(x))/s
      if (is .lt. 0) dround = xnint(dble(x)/s)*s
#else
      if (is .gt. 0) dround = dnint(s*x)/s
      if (is .lt. 0) dround = dnint(x/s)*s
#endif
      end
