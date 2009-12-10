      subroutine ume(nlme,d,v)
C- Returns Harrison's universal matrix elements
C ----------------------------------------------------------------
Ci Inputs
Ci   
Co Outputs
Co   
Cr Remarks
Cr   Converted to atomic units
Cu Updates
Cu    8 Jun 07 (MvS) Merged Klepeis's additions to TB package
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nlme
      double precision d,v(nlme)
C Local variables
      integer i
      double precision eta(10),xx
C 1981 values ..
      data eta /-1.32d0,1.42d0,2.22d0,-.63d0,0d0,0d0,0d0,0d0,0d0,0d0/
C Solid State Table ...
C      data eta /-1.4d0,1.84d0,3.24d0,-.81d0,0d0,0d0,0d0,0d0,0d0,0d0/

      call rxx(nlme .gt. 10,'UME: nlme > 10')
      call dpzero(v,nlme)
C .. eV :
C      xx = 7.62/d**2
C .. atomic units :
      xx = 2/d**2
      do  10  i = 1, nlme
        v(i) = eta(i)*xx
   10 continue

c     print *, (v(1) - 2*dsqrt(3d0)*v(2) - 3*v(3))/4
      end
