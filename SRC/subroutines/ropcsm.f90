      subroutine ropcsm(m,n,x,y,w,cm,sm)
C- Makes cm and sm. Must be called in sequence m=0,1,2...
      implicit none
      integer m,n,i
      double precision x(n),y(n),w(n),cm(n),sm(n)

C --- Case m=0 ---
      if (m .eq. 0) then
        do  1  i = 1, n
          cm(i) = 1d0
    1   continue
        do  2  i = 1, n
          sm(i) = 0d0
    2   continue
        return
      endif

C --- Case m=1 ---
      if (m .eq. 1) then
        do  3  i = 1, n
          cm(i) = x(i)
    3   continue
        do  4  i = 1, n
          sm(i) = y(i)
    4   continue
        return
      endif

C --- Case m ge 2 ---
      if (m .ge. 2) then
        do  5  i = 1, n
          w(i) = cm(i)
    5   continue
        do  6  i = 1, n
          cm(i) = x(i)*cm(i) - y(i)*sm(i)
    6   continue
        do  7  i = 1, n
          sm(i) = y(i)*w(i) + x(i)*sm(i)
    7   continue
        return
      endif
      end

