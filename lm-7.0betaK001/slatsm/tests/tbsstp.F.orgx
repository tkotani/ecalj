C Test subroutine dbsstp
      subroutine fmain
      implicit none
      integer n,mx,i,ipr,ir,mi,k,ncall
C ... can make mx higher, eg, 7
      parameter(n=4,mx=4,mi=11)
      double precision y(n),dydx(n),
     .  wk(n*(3+mx)+5+1),yscal(n),yerr(n),yt(n),
     .  x,dbesj,hmin,htry,h,tol,hnext,
     .  phi(0:n),psi(0:n),hdid,xs,err
      integer nseq(0:mi)
      data nseq /11,2,3,4,6,8,12,16,24,32,48,64/
C      data nseq /11,2,4,6,8,12,16,24,32,48,64,96/

      ipr = 40
      wk(n*(3+mx)+5+1) = 1d10
      nseq(0) = 0
C ... Integration of bessel functions
      x = 5d0
      do  10  i = 0, 3
   10 y(i+1) = dbesj(i,x)
      dydx(1) = -y(2)
      dydx(2) = y(1)-y(2)
      dydx(3) = y(2)-2d0*y(3)
      dydx(4) = y(3)-3d0*y(4)
      do  11 i = 1, n
   11 yscal(i) = 1d0
      write(*,'(/1x,t8,a,t19,a,t31,a,t43,a/)')
     .  'tol','htry','hdid','hnext  calls'
      hmin = .002d0
      do  12  i = 1, 30
        htry = 1d0
        h = htry
C ...   this changes tol with increasingly strict precision
        tol = dexp(-dble(i))
C#ifdefC FIXEDPREC
C        tol = -1d-12
C#endif
        xs = x
   14   call derivs(n,x,y,dydx)
        call dbsstp(x,y,dydx,yerr,n,tol,hmin,h,yscal,mx,nseq,
     .    ir,wk,ipr)
        if (ir .eq. -1 .and. tol .lt. 0) then
          print *, '(warning): dbsstp did not meet tol'
          ir = 0
        else
          if (ir .lt. 0) stop 'we have a problem'
        endif
        if (ir .gt. 0) goto 14
        ncall = nint(wk(5))+1
        hdid = x-xs
        write(*,'(2x,e12.4,f8.2,2x,2f12.6,i4)') tol,htry,hdid,h,ncall
        do  16  k = 0, 3
   16   yt(k+1) = dbesj(k,x)
        err = 0
        do  18  k = 1, 4
   18   err = err + (y(k)-yt(k))**2
        err = dsqrt(err)
C#ifdefC FIXEDPREC
C        print 345, x, y(2),y(4), err, err/abs(tol)
C  345   format(f10.6, ' b2,b4 num int',2f10.6,' prec',f20.15,
C     .    '  prec/tol',f20.15)
C        pause
C#endif
   12 continue

      print *, 'show wk dim correctly', wk(n*(3+mx)+5+1)

      end
      
      subroutine derivs(n,x,y,dydx)
      implicit none
      integer n
      double precision x,y(n),dydx(n)
      dydx(1) = -y(2)
      dydx(2) = y(1) - (1d0/x)*y(2)
      dydx(3) = y(2) - (2d0/x)*y(3)
      dydx(4) = y(3) - (3d0/x)*y(4)
C      print 334, x, y, dydx(1)
C  334 format('derivs: x,y',f10.6,2x,5f10.6)
      end
