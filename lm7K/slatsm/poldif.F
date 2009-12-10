      subroutine poldif(x,y,np,n,lopt,errmx,yp)
C- Differentiate by polynomial interpolation tabulated function
C ----------------------------------------------------------------
Ci Inputs
Ci    x,y,np: table of x-y pairs and number of points
Ci    n: maximum number of points to interpolate
C     lopt nonzero => evaluate intermediate results in quad precision
Co Outputs
Co    yp: differentiated value of y at each x
Co    errmx: estimate of error
C ----------------------------------------------------------------
C     implicit none
      integer np,n,lopt
      double precision x(np),y(np),yp(np),errmx
C Local variables
      integer i,nmax,ip,i0,i1
      parameter (nmax=50)
      double precision c(nmax),xmid
C     double precision sum

      if (n .gt. nmax) call rx('poldif: increase nmax')
      if (n+1 .ge. np) call rx('poldif: n ge np')
      if (n .le. 0 .or. np .lt. 1) return
      errmx = 0
      if (np .eq. 1) return


      do  10  ip = 2, np+1

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
       i0 = min(i0,np-n)
       if (i0 .lt. 0) call rx('bug in poldif')

       call polcof(lopt,x(ip-1),x(i0+1),y(i0+1),n,c)
C       errmx = max(errmx,sum*x(ip)**(n-1))
C       sum = c(n)*(n-1)
C       do  12  i = n-1, 2, -1
C   12  sum = (i-1)*c(i) + sum*x(ip-1)
C       yp(ip-1) = sum
       yp(ip-1) = c(2)
C  ... for now
       errmx = 0

C       print *, i0,n,x(ip-1)
C       print 333, (c(i), i=1,n)
C  333  format(5f15.10)
C       print *, 'sum=',sum

   10 continue

      end

