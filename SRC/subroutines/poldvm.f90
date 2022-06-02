!      module m_poldvm

!      interface poldvm
!        module procedure poldvm
!        subroutine poldvm(x,y,np,n,lrat,tol,lerr,yp)
!        implicit none
!        integer,intent(in):: np,n
!        logical,intent(in):: lrat
!        double precision,intent(in):: x(np),y(np),tol ! define these as double precision because of the definition in slatsm/poldvm
!        double precision, intent(out):: yp(np)
!        integer,intent(out):: lerr
!        end subroutine poldvm
!      end interface

!      contains

subroutine poldvm(x,y,np,n,lrat,tol,lerr,yp)
  !- Derivative of tabulated function by rational function interpolation
  ! ----------------------------------------------------------------
  !i Inputs
  !i    x,y,np: table of x-y pairs and number of points
  !i    n     : maximum number of points to interpolate
  !i    lrat  : T, use rational function interpolation; otherwise
  !i          : use polynomial interpolation
  !i    tol   : increase number of points in fit until n points is
  !i          : reached or when estimated value of next term < tol
  !o Outputs
  !o    yp    : interpolated value of y at each x
  !o    lerr  : (lrat only) first point which rat int had problems
  ! ----------------------------------------------------------------
  !     implicit none
  integer :: np,n
  logical :: lrat
  double precision :: dy,tol,x(np),y(np),yp(np)
  ! Local variables
  integer :: i,m,ns,ii,ip,i0,iipm,lerr!,nmax
  !      parameter (nmax=10000)
  double precision :: den,ho,hp,w,dd,t,c(n),d(n), &
       tiny,d1mach

  !      if (n .gt. nmax) call rx('poldvm: increase nmax')
  if (n+1 > np) call rx('poldvm: n+1 gt np')
  tiny = dsqrt(d1mach(1))
  if ( .NOT. lrat) tiny = 0d0
  lerr = 0

  do  10  ip = 1, np

     ! --- index ns closest table entry ---
     i0 = ip-1
     i0 = min(max(i0-n/2,0),np-n-1)
     ns = 0
     do  20  i = 1, n
        ii = i+i0
        if (ii <= ip) ns = i
        if (ii >= ip) ii = ii+1
        c(i) = (y(ii)-y(ip))/(x(ii)-x(ip))
        d(i) = (y(ii)-y(ip))/(x(ii)-x(ip)) + tiny
20   enddo

     if (ns == 0 .OR. ii > np) call rx('bug in poldvm')

     ! --- For each column of the Neville recursion update C and D ---
     yp(ip) = c(ns)
     ns = ns-1
     do  23  m = 1, n-1
        do  22  i = 1, n-m
           ii = i+i0
           iipm = i+m+i0
           if (ii >= ip) ii = ii+1
           if (iipm >= ip) iipm = iipm+1
           ho = x(ii)-x(ip)
           hp = x(iipm)-x(ip)
           w  = c(i+1)-d(i)
           if (lrat) then
              t = ho*d(i)/hp
              dd = t-c(i+1)
              ! ... Error if pole; skip this point
              if (dd == 0) then
                 yp(ip) = 0
                 if (lerr == 0) lerr = ip
                 goto 10
              endif
              dd = w/dd
              d(i) = c(i+1)*dd
              c(i) = t*dd
           else
              den = w/(ho-hp)
              d(i) = hp*den
              c(i) = ho*den
           endif
22      enddo
        ! ---   Decide which whether to accumulate correction C or D to y ---
        if (2*ns < n-m) then
           dy = c(ns+1)
        else
           dy = d(ns)
           ns = ns-1
        endif
        yp(ip) = yp(ip)+dy
        if (dabs(dy) < tol) goto 10
23   enddo

10 enddo

end subroutine poldvm
!      subroutine fmain
!C      implicit none
!      integer np,i,n,lerr
!      logical lrat
!      parameter (np=20)
!      double precision xa(np),ya(np),tol,dy,pi,fa(np),y(np)
!      pi = 4*datan(1d0)
!      tol = 1d-8
!      n = 3
!      lrat = .false.
!      print *, 'n,lrat='
!      read(*,*) n,lrat

!C ... Note rational function interpolation doesn't work well here!
!C      do 11 i = 1, np
!C      xa(i) = i*pi/np
!C      fa(i) = dcos(xa(i))
!C   11 ya(i) = dsin(xa(i))
!      do 11 i = 1, np
!      xa(i) = dble(i)/np
!      fa(i) = dexp(xa(i))
!   11 ya(i)=dexp(xa(i))
!      call poldvm(xa,ya,np,n,lrat,tol,lerr,y)
!      print *, 'poldvm returned with lerr=',lerr
!      write(*,'(t10,a1,t20,a5,t27,a12,t46,a5)')
!     .  'x','f''(x)','interpolated','error'
!      do  13  i = 1, np
!        write(*,'(1x,3f12.6,1pe15.4)') xa(i),fa(i),y(i),fa(i)-y(i)
!   13 continue
!      end

!      end module m_poldvm
