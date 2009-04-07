C     Tests huntx
      subroutine fmain
C     implicit none
      integer nn
      parameter (nn=100)
      integer i,iprm(nn),low1,low2
      double precision xa(nn),wk(nn),tol,xv
      real ran1
      parameter (tol=1d-6)

      call ran1in(10)
      do  10  i= 1, nn
   10 xa(i) = 10*(2*ran1()-1)

      xa(nn-1) = xa(2) + tol/2

      call dvheap(1,nn,xa,iprm,tol,1)
      call dvprm(1,nn,xa,wk,iprm,.false.)
      call dscal(nn,-1d0,wk,1)
      call dscal(nn,-1d0,xa,1)

      print *, ' ... array is ordered in descending order'
      print *, '     compare sorted xa with permuted xa'
      low1 = 0
      low2 = 0
      print *, 'low1 low2     val  val-x(low1) val-x(low1+1)'
      do  20  xv = -10d0, 10d0, 1d0
        call huntx(wk,nn,xv,0,low1)
        call huntx(xa,nn,xv,iprm,low2)
        print 333, low1, low2,
     .    xv, wk(max(low1,1))-xv, wk(min(low1+1,nn))-xv
        if (low1 .ne. low2) stop 'oops'
        if (xv .lt. 10d0 .and. wk(max(low1,1)) .lt. xv) print *, 'oops'
        if (xv .gt.-10d0 .and. wk(min(low1+1,nn)) .gt. xv) print*,'oops'
  333   format(2i5,3f10.3)
   20 continue

      end
