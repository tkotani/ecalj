      subroutine politg(x,y,np,n,ip1,ip2,lopt,errmx,yp)
C- Integral by polynomial interpolation of function tabulated on mesh
C ----------------------------------------------------------------
Ci Inputs
Ci   x     :table of abscissas
Ci   y     :table of ordinates
Ci   np    :number of (x,y) pairs
Ci   n     :maximum number of points to include in poly. interpolation
Ci   ip1   :integrate points x(ip1) as low limit
Ci   ip2   :integrate points between x(ip1)..x(ip2)
Ci  lopt   :nonzero => evaluate intermediate results in quad precision
Co Outputs
Co  yp     :integrated value of y at each x(ip1..ip2)
Co  errmx  :estimate of error
Cu Updates
Cu   13 Dec 02 New arguments ip1,ip2
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer np,n,lopt,ip1,ip2
      double precision x(np),y(np),yp(ip1:np),errmx
C ... Local parameters
      integer i,nmax,ip,i0,i1
      parameter (nmax=50)
      double precision c(nmax),xmid,sum

      if (n .gt. nmax) call rx('politg:  increase nmax')
      if (n .gt. np)   call rx('politg:  n ge np')
      if (ip1 .lt. 1)  call rx('politg:  ip1 lt 1')
      if (ip2 .gt. np) call rx('politg:  ip2 gt np')
      if (n .le. 0 .or. np .lt. 1) return
      if (ip1 .gt. ip2) return

      yp(ip1) = 0
      errmx = 0
      if (ip1 .eq. ip2) return
      do  10  ip = ip1+1, ip2

C --- Get i0,i1 = left-to-leftmost and rightmost in tableau ---
      xmid = (x(ip)+x(ip-1))/2
      i0 = ip-1
      i1 = ip
      do  2  i = 1, n-1
        if (i0 .eq. 0) then
          i1 = i1+1
        elseif (i1 .eq. np) then
          i0 = i0-1
        elseif (dabs(x(i0)-xmid) .lt. dabs(x(i1+1)-xmid))  then
          i0 = i0-1
        else
          i1 = i1+1
        endif
    2  continue

       call polcof(lopt,x(ip-1),x(i0+1),y(i0+1),n,c)
       sum = c(n)/n
       errmx = errmx + sum*(x(ip)-x(ip-1))**(n+1)
       do  12  i = n-1, 1, -1
   12  sum = c(i)/i + sum*(x(ip)-x(ip-1))
       sum = sum*(x(ip)-x(ip-1))
       yp(ip) = yp(ip-1) + sum

C       print 333, (c(i), i=1,n)
C  333  format(5f15.10)
C       print *, 'sum=',sum

   10 continue

      end

