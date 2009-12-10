      integer function mxint(n,ivec)
C- Return the largest integer in ivec
C     implicit none
      integer n,ivec(n)
      integer imax,i
      
      mxint = 0
      if (n .le. 0) return
      imax = ivec(1)
      do  10  i = 1, n
   10 imax = max(imax,ivec(i))
      mxint = imax
      end

