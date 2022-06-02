subroutine ropcsm(m,n,x,y,w,cm,sm)
  !- Makes cm and sm. Must be called in sequence m=0,1,2...
  implicit none
  integer :: m,n,i
  double precision :: x(n),y(n),w(n),cm(n),sm(n)

  ! --- Case m=0 ---
  if (m == 0) then
     do  1  i = 1, n
        cm(i) = 1d0
1    enddo
     do  2  i = 1, n
        sm(i) = 0d0
2    enddo
     return
  endif

  ! --- Case m=1 ---
  if (m == 1) then
     do  3  i = 1, n
        cm(i) = x(i)
3    enddo
     do  4  i = 1, n
        sm(i) = y(i)
4    enddo
     return
  endif

  ! --- Case m ge 2 ---
  if (m >= 2) then
     do  5  i = 1, n
        w(i) = cm(i)
5    enddo
     do  6  i = 1, n
        cm(i) = x(i)*cm(i) - y(i)*sm(i)
6    enddo
     do  7  i = 1, n
        sm(i) = y(i)*w(i) + x(i)*sm(i)
7    enddo
     return
  endif
end subroutine ropcsm

