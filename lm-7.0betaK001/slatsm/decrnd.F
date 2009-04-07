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
#if TEST
      double precision decrnd,x
      x = 12345678
C     x = -x
      do  20 j = 1,7
        do  10  i = -20,20
   10   print *, decrnd(x*(10d0**i),j), decrnd(x*(10d0**i),j)/
     .      (x*(10d0**i))
        print *, '----'
   20 continue
      end
#endif

