c      module m_poldvm

c      interface poldvm
c        module procedure poldvm
c        subroutine poldvm(x,y,np,n,lrat,tol,lerr,yp)
c        implicit none
c        integer,intent(in):: np,n
c        logical,intent(in):: lrat
c        double precision,intent(in):: x(np),y(np),tol ! define these as double precision because of the definition in slatsm/poldvm 
c        double precision, intent(out):: yp(np)
c        integer,intent(out):: lerr
c        end subroutine poldvm
c      end interface

c      contains

      subroutine poldvm(x,y,np,n,lrat,tol,lerr,yp)
C- Derivative of tabulated function by rational function interpolation
C ----------------------------------------------------------------
Ci Inputs
Ci    x,y,np: table of x-y pairs and number of points
Ci    n     : maximum number of points to interpolate
Ci    lrat  : T, use rational function interpolation; otherwise
Ci          : use polynomial interpolation
Ci    tol   : increase number of points in fit until n points is
Ci          : reached or when estimated value of next term < tol
Co Outputs
Co    yp    : interpolated value of y at each x
Co    lerr  : (lrat only) first point which rat int had problems
C ----------------------------------------------------------------
C     implicit none
      integer np,n
      logical lrat
      double precision dy,tol,x(np),y(np),yp(np)
C Local variables
      integer i,m,ns,ii,ip,i0,iipm,lerr!,nmax
c      parameter (nmax=10000)
      double precision den,ho,hp,w,dd,t,c(n),d(n),
     .tiny,d1mach

c      if (n .gt. nmax) call rx('poldvm: increase nmax')
      if (n+1 .gt. np) call rx('poldvm: n+1 gt np')
      tiny = dsqrt(d1mach(1))
      if (.not. lrat) tiny = 0d0
      lerr = 0

      do  10  ip = 1, np

C --- index ns closest table entry ---
        i0 = ip-1
        i0 = min(max(i0-n/2,0),np-n-1)
        ns = 0
        do  20  i = 1, n
          ii = i+i0
          if (ii .le. ip) ns = i
          if (ii .ge. ip) ii = ii+1
          c(i) = (y(ii)-y(ip))/(x(ii)-x(ip))
          d(i) = (y(ii)-y(ip))/(x(ii)-x(ip)) + tiny
   20   continue

        if (ns .eq. 0 .or. ii .gt. np) call rx('bug in poldvm')

C --- For each column of the Neville recursion update C and D ---
        yp(ip) = c(ns)
        ns = ns-1
        do  23  m = 1, n-1
          do  22  i = 1, n-m
            ii = i+i0
            iipm = i+m+i0
            if (ii .ge. ip) ii = ii+1
            if (iipm .ge. ip) iipm = iipm+1
            ho = x(ii)-x(ip)
            hp = x(iipm)-x(ip)
            w  = c(i+1)-d(i)
            if (lrat) then
              t = ho*d(i)/hp
              dd = t-c(i+1)
C ... Error if pole; skip this point
              if (dd .eq. 0) then
                yp(ip) = 0
                if (lerr .eq. 0) lerr = ip
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
   22     continue
C ---   Decide which whether to accumulate correction C or D to y ---
          if (2*ns .lt. n-m) then
            dy = c(ns+1)
          else
            dy = d(ns)
            ns = ns-1
          endif
          yp(ip) = yp(ip)+dy
          if (dabs(dy) .lt. tol) goto 10
   23   continue

   10 continue

      end subroutine poldvm
C      subroutine fmain
CC      implicit none
C      integer np,i,n,lerr
C      logical lrat
C      parameter (np=20)
C      double precision xa(np),ya(np),tol,dy,pi,fa(np),y(np)
C      pi = 4*datan(1d0)
C      tol = 1d-8
C      n = 3
C      lrat = .false.
C      print *, 'n,lrat='
C      read(*,*) n,lrat
C
CC ... Note rational function interpolation doesn't work well here!
CC      do 11 i = 1, np
CC      xa(i) = i*pi/np
CC      fa(i) = dcos(xa(i))
CC   11 ya(i) = dsin(xa(i))
C      do 11 i = 1, np
C      xa(i) = dble(i)/np
C      fa(i) = dexp(xa(i))
C   11 ya(i)=dexp(xa(i))
C      call poldvm(xa,ya,np,n,lrat,tol,lerr,y)
C      print *, 'poldvm returned with lerr=',lerr
C      write(*,'(t10,a1,t20,a5,t27,a12,t46,a5)')
C     .  'x','f''(x)','interpolated','error'
C      do  13  i = 1, np
C        write(*,'(1x,3f12.6,1pe15.4)') xa(i),fa(i),y(i),fa(i)-y(i)
C   13 continue
C      end

c      end module m_poldvm
