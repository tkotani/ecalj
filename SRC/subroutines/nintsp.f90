      subroutine nintsp(nx,np,nxp,x,w)
C- Points and weights for Legendre integration on a sphere surface
C ----------------------------------------------------------------------
Ci Inputs
Ci   nx    :number of points in polar angle (Legendre integration)
Ci   np    :number of points in phi (uniform points)
Ci         :np=0 => np depends on nx making dphi approx constant.
Co Outputs
Co   nxp   :total number of number of points in quadrature
Co   x     :cartesian coordinates of points on unit sphere
Co   w     :weights corresponding to x
Cr Remarks
Cr   nintsp generates a mesh of points on a unit sphere for angular
Cr   integration, using Legendre gaussian quadrature for the polar
Cr   angle, and evenly spaced points on a circle for phi integration
Cr   See fpiint as an alternative routine that uses special points.
C ----------------------------------------------------------------------
C     implicit none
      integer nx,np,nxp
      double precision x(3,1),w(1)
      integer ix,ip,i,npi
      double precision xtab(100),wx(100),sint,cost,phi,pi2

      pi2 = 8*datan(1d0)
      call mklegw(nx,xtab,wx,0)
      i = 0
      do    ip = 1, np
        do  ix = 1, nx
          i = i+1
          cost = xtab(ix)
          phi =  (pi2*ip)/np
          sint = dsqrt(1-cost**2)
          x(1,i)=dcos(phi)*sint
          x(2,i)=dsin(phi)*sint
          x(3,i)=cost
          w(i)  =wx(ix)*(pi2/np)
       enddo
       enddo
      nxp = nx*np
      if (np .gt. 0) return

C --- Let np depend on x ---
      i = 0
      do  20  ix = 1, nx
        cost = xtab(ix)
        sint = dsqrt(1-cost**2)
        npi  = max(nint(sint*2*nx),1)
        do  30  ip = 1, npi
          i = i+1
          phi =  (pi2*ip)/npi
          x(1,i)=dcos(phi)*sint
          x(2,i)=dsin(phi)*sint
          x(3,i)=cost
          w(i)  =wx(ix)*(pi2/npi)
   30   continue
   20 continue
      nxp = i
      end

