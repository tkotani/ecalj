      double precision function dotprd(n,dx,incx,dy,incy)
C- Angle between two vectors
C ----------------------------------------------------------------------
Ci Inputs
Ci   n     : vector length
Ci   dx    : first vector
Ci   incx  : stride length of first vector
Ci   dy    : second vector
Ci   incy  : stride length of second vector
Co Outputs
Co  dotprd :angle between dx and dy, in radians
Cl Local variables
Cl         :
Cr Remarks
Cr   
Cu Updates
Cu   30 Oct 03 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      double precision dx(*),dy(*)
      integer n,incx,incy
C ... Local parameters
      double precision ddot,xx1,xx2,dsqrt

      xx1 = ddot(n,dx,incx,dx,incx)*ddot(n,dy,incy,dy,incy)
      xx2 = ddot(n,dx,incx,dy,incy)
      dotprd = dacos(min(1d0,xx2/dsqrt(xx1)))
      end
C test:
C      double precision dx(2,3),dy(1,3), dotprd,pi
C      data dx /4d0,0d0,4d0,0d0,.01d0,0d0/
C      data dy /3d0,-3d0,.01d0/
C
C      pi = 4*datan(1d0)
C
C      print *, dotprd(2,dx,2,dy,1)*180/pi
C      print *, dotprd(3,dx,2,dy,1)*180/pi
C
C      end


