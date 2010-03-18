      double precision function decrnd(x,n)
C- Truncate after n significant figures 
C ----------------------------------------------------------------
Ci Inputs
Ci   x,n
Co Outputs
Co   decrnd
Cr Remarks
Cr   
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      double precision x
      integer n
C Local parameters
      double precision scale
C This construction avoids taking integer part of negative number
      if (x .eq. 0) then
        decrnd = 0
        return
      endif
      scale = dlog10(dabs(x))
      if (scale .gt. 0) then
        scale = int(scale)
      else
        scale = -int(-scale) - 1
      endif
      scale = 10d0**(n-1-scale)
      decrnd = nint(x*scale)/scale
      end
c$$$#if TEST
c$$$      double precision decrnd,x
c$$$      x = 12345678
c$$$C     x = -x
c$$$      do  20 j = 1,7
c$$$        do  10  i = -20,20
c$$$   10   print *, decrnd(x*(10d0**i),j), decrnd(x*(10d0**i),j)/
c$$$     .      (x*(10d0**i))
c$$$        print *, '----'
c$$$   20 continue
c$$$      end
c$$$#endif

