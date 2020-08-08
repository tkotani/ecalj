      subroutine ishell(n,iarray)
      integer n
      integer iarray(1)
      integer lognb2,i,j,k,l,m,nn,it

      if (n .le. 1) return
      lognb2 = int(log(float(n+1))*1.4426950)
      m = n
      do  12  nn = 1, lognb2
        m = m/2
        k = n - m
        do  11  j = 1, k
          i = j
    3     continue
          l = i + m
          if (iarray(l) .lt. iarray(i)) then
            it = iarray(i)
            iarray(i) = iarray(l)
            iarray(l) = it
            i = i - m
            if (i .ge. 1) goto 3
          endif
   11   continue
   12 continue
      return
      end

