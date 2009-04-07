      integer function ichksm(n,ix)
c
c     Primitive check sum-take sum of -1**i*ix(i)

      integer ix(1)
      integer n,i
c
      ichksm = 0
      if (n .le. 0) return
      m1 = 1
      do  10  i = 1, n
        m1 = -m1
        ichksm = ichksm + m1*ix(i)
   10 continue
      end


