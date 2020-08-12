CSFCPP#define F90 1
      subroutine hrmint(xi,yi,ypi,n,nyi,nypi,x,y,yp)
C- Hermite interpolation of tabulated data at a specified point
C ----------------------------------------------------------------
Ci Inputs
Ci   xi    :set of n abscissa values, ordered by increasing value
Ci   yi    :set of n ordinate values
Ci   ypi   :set of n derivatives dyi/dx
Ci   nyi   :number of points yi used in the interpolation
Ci         := 1+order of polynomial
Ci   nypi  :number of derivatives ypi used in the interpolation
Ci   x     :point at which to evaluate function
Co Outputs
Co   y     :interpolated value of yi at x
Co   yp    :interpolated value of derivative dy/dx
Cr Remarks
Cr   This is a driver for hrminx which sorts points by increasing
Cr   distance from specified point x.  Thus points closest to x
Cr   are used for interpolation.
Cu Updates
Cu   18 Sep 00 Sorts points by increasing distance from x
Cu    9 Sep 00 Added nypi; also returns yp
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer n,nypi,nyi
      double precision x,y,yp,xi(n),yi(n),ypi(n)
C Local variables
      integer j
CSFCPP#if AUTO_ARRAY | F90
      double precision xx(n)
      integer iprm(n)
CSFCPP#else
CSFCPP      integer nmax
CSFCPP      parameter (nmax=1000)
CSFCPP      double precision xx(nmax)
CSFCPP      integer iprm(nmax)
CSFCPP      if (n .gt. nmax) call rxi('hrmint: increase nmax, need',n)
CSFCPP#endif

C     Order points by increasing distance from x
      do  j = 1, n
        xx(j) = xi(j) - x
      enddo
      call dvheap(1,n,xx,iprm,0d0,11)
      call hrminx(xi,yi,ypi,min(nyi,n),min(nypi,n),iprm,x,y,yp)

C  ... OLD ... along lines of polint ...
C      np = min(nyi,n)
C      call hunty(xi,n,np,x,jx,j)
C      jx = min(max(1,jx-np/2),n-np+1)
C
C      nmax = min(nyi,n-j+1)
C      mmax = min(nypi,n-j+1)
C      call hrminx(xi(j),yi(j),ypi(j),nmax,mmax,x,y,yp)
      end

      subroutine hrminx(xi,yi,ypi,N,M,iprm,x,y,yp)
C- Kernel called by hrmint
C ----------------------------------------------------------------
Ci Inputs
Ci   xi    :set of n abscissa values
Ci   yi    :set of n ordinate values y(xi)
Ci   ypi   :set of m derivatives (dy/dx) at xi
Ci   iprm  :if iprm(1) is nonzero, permute xi by iprm
Ci   N     :number of function values y
Ci   M     :number of function derivatives ypi; M<=N
Ci   x     :point at which to evaluate function and derivative
Co Outputs
Co   y     :interpolated value of y(x) using yi(1..n) and ypi(1..m)
Co   y     :derivative of y
Cr Remarks
Cr   The first m function derivatives are used in the interpolation
Cr   The hermite interpolating polynomial is
Cr   y(x) = sum_n=1..M L_n^M L_n^N c_n + sum_n=M+1..N L_n^M L_n^N y_n
Cr   where
Cr   c_n = (x_n - x) (a_n*y_n - y1_n) + y_n
Cr   a_n = sum_m=1..M  2/(x_n - x(m)) + sum_m=M+1..N 1/(x_n - x(m))
Cu Updates
Cu    9 Sep 00 Allow m to be different from n; also returns yp
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer M,N,iprm(N)
      double precision x,y,yp,xi(n),yi(n),ypi(n)
C Local variables
      logical lprm
      integer nn,nnx,n2
      parameter (nnx=100)
      double precision b,c
CSFCPP#if AUTO_ARRAY | F90
      double precision Pm(N),Pn(N),Pmp(N),Pnp(N),a(N)
CSFCPP#else
CSFCPP      double precision Pm(nnx),Pn(nnx),Pmp(nnx),Pnp(nnx),a(nnx)
CSFCPP      if (N .gt. nnx) call rxi('hrminx: increase nnx; need',N)
CSFCPP#endif

C ... Switch whether to permute xi
      lprm = iprm(1) .ne. 0

      call hrmins(xi,iprm,N,M,x,Pm,Pn,Pmp,Pnp,a)

      y = 0
      yp = 0
      do  nn = 1, M
        if (lprm) then
          n2 = iprm(nn)
        else
          n2 = nn
        endif
        b = a(nn)*yi(n2) - ypi(n2)
        c = (xi(n2) - x)*b + yi(n2)
        y = y + Pm(nn)*Pn(nn)*c
        yp = yp - Pm(nn)*Pn(nn)*b + (Pmp(nn)*Pn(nn)+Pm(nn)*Pnp(nn))*c
      enddo

      do  nn = M+1, N
        if (lprm) then
          n2 = iprm(nn)
        else
          n2 = nn
        endif
        y = y + Pm(nn) * Pn(nn) * yi(n2)
        yp = yp + (Pmp(nn)*Pn(nn) + Pm(nn)*Pnp(nn)) * yi(n2)
      enddo

      end

      subroutine hrmins(xi,iprm,N,M,x,Pm,Pn,Pmp,Pnp,a)
C- Setup for Hermite interpolation
C ----------------------------------------------------------------------
Ci Inputs
Ci   xi    :set of N points
Ci   iprm  :if iprm(1) is nonzero, permute xi by iprm
Ci   N     :number of points
Ci   M     :number of points for which derivative is used
Ci   x     :point at which to evaluate function and derivative
Co Outputs
Co   Pm    :Pm(n) = L_n^M  = prod_(m=1..M,m<>n) (x-xm)/(xn-xm)
Co   Pn    :Pm(n) = L_n^N  = prod_(m=1..N,m<>n) (x-xm)/(xn-xm)
Co   Pmp   :d(Pm)/dx
Co   Pnp   :d(Pn)/dx
Co   a     :a(n) = sum_(m=1..M)   2/(x(n) - x(m)) +
Co         :       sum_(m=M+1..N) 1/(x(n) - x(m))
Cr Remarks
Cr   See hrminx for usage of output for interpolation
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer N,M,iprm(N)
      double precision x,xi(n),Pm(N),Pn(N),Pmp(N),Pnp(N),a(N)
C Local variables
      logical lprm
      integer mm,nn,lzm,lzn,n2,m2
      double precision fac,fuzz,denm,denn

C ... Switch whether to permute xi
      lprm = iprm(1) .ne. 0

C ... Length scale for which x-xi is small
      fuzz = 0
      do  nn = 1, n
        if (lprm) then
          n2 = iprm(nn)
        else
          n2 = nn
        endif
        fuzz = fuzz + xi(n2)**2
      enddo
      fuzz = 1d-10*sqrt(fuzz)/n

      do  nn = 1, N
        if (lprm) then
          n2 = iprm(nn)
        else
          n2 = nn
        endif
        a(nn) = 0
        Pm(nn) = 1
        Pn(nn) = 1
        Pmp(nn) = 1
        Pnp(nn) = 1
        denm = 0
        denn = 0
        fac = 2
        lzm = 0
        lzn = 0
        do  mm = 1, N
          if (lprm) then
            m2 = iprm(mm)
          else
            m2 = mm
          endif
          if (mm .gt. M) fac = 1
          if (mm .ne. nn) then
            if (mm .le. M) then
              Pm(nn) = Pm(nn) * (x - xi(m2)) / (xi(n2) - xi(m2))
              if (abs(x - xi(m2)) .lt. fuzz) then
                Pmp(nn) = Pmp(nn)  / (xi(n2) - xi(m2))
                lzm = 1
              else
                Pmp(nn) = Pmp(nn) * (x - xi(m2)) / (xi(n2) - xi(m2))
                denm = denm + 1/(x - xi(m2))
              endif
            endif
            Pn(nn) = Pn(nn) * (x - xi(m2)) / (xi(n2) - xi(m2))
            if (abs(x - xi(m2)) .lt. fuzz) then
              Pnp(nn) = Pnp(nn) / (xi(n2) - xi(m2))
              lzn = 1
            else
              Pnp(nn) = Pnp(nn) * (x - xi(m2)) / (xi(n2) - xi(m2))
              denn = denn + 1/(x - xi(m2))
            endif
            a(nn) = a(nn) + fac/(xi(n2) - xi(m2))
          endif
        enddo
        if (lzm .eq. 0) Pmp(nn) = Pmp(nn) * denm
        if (lzn .eq. 0) Pnp(nn) = Pnp(nn) * denn
      enddo

      end
C      subroutine fmain
C      implicit none
C      integer ifi,fopng,rdm,nr,nc,np,i,npmx,nplot
C      parameter (npmx=1000)
C      double precision s(npmx),x,y,yp
C      integer nypi,nyi
C
C      write(*,'(''#  ... reading from file dat'')')
C      ifi = fopng('dat',-1,1)
C      rewind ifi
C      np = 0
C      nr = 0
C      nc = 3
C      i = rdm(ifi,1000,np,' ',s,nr,nc)
C      if (nr*nc .gt. npmx) stop 'increase npmx'
C      i = rdm(ifi,1000,npmx,' ',s,nr,nc)
C      print '(''# Read nr='',i4,''.  m,n=?'')', nr
C      nypi = nr
C      nyi = nr
C      read(*,*) nypi,nyi
C      if (nypi .gt. nyi) stop 'm must be le n'
C      nplot = 61
C      write(*,'(''% rows'',i4,'' cols 3'')') nplot
C      do  i = 1, nplot
C        x = dble(i-1)/5
C        call hrmint(s,s(1+nr),s(1+2*nr),nr,nyi,nypi,x,y,yp)
C        print *, x,y,yp
C      enddo
C
C      end

