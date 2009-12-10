      subroutine ropcsm(m,n,x,y,w,cm,sm)
C- Makes cm and sm. Must be called in sequence m=0,1,2...
C     implicit none
      integer m,n,i
      double precision x(n),y(n),w(n),cm(n),sm(n)

C --- Case m=0 ---
      if (m .eq. 0) then
        do  1  i = 1, n
    1   cm(i) = 1d0
        do  2  i = 1, n
    2   sm(i) = 0d0
        return
      endif

C --- Case m=1 ---
      if (m .eq. 1) then
        do  3  i = 1, n
    3   cm(i) = x(i)
        do  4  i = 1, n
    4   sm(i) = y(i)
        return
      endif

C --- Case m ge 2 ---
      if (m .ge. 2) then
        do  5  i = 1, n
    5   w(i) = cm(i)
        do  6  i = 1, n
    6   cm(i) = x(i)*cm(i) - y(i)*sm(i)
        do  7  i = 1, n
    7   sm(i) = y(i)*w(i) + x(i)*sm(i)
        return
      endif
      end

