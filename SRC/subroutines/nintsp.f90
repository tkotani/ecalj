subroutine nintsp(nx,np,nxp,x,w)
  !- Points and weights for Legendre integration on a sphere surface
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nx    :number of points in polar angle (Legendre integration)
  !i   np    :number of points in phi (uniform points)
  !i         :np=0 => np depends on nx making dphi approx constant.
  !o Outputs
  !o   nxp   :total number of number of points in quadrature
  !o   x     :cartesian coordinates of points on unit sphere
  !o   w     :weights corresponding to x
  !r Remarks
  !r   nintsp generates a mesh of points on a unit sphere for angular
  !r   integration, using Legendre gaussian quadrature for the polar
  !r   angle, and evenly spaced points on a circle for phi integration
  !r   See fpiint as an alternative routine that uses special points.
  ! ----------------------------------------------------------------------
  !     implicit none
  integer :: nx,np,nxp
  double precision :: x(3,1),w(1)
  integer :: ix,ip,i,npi
  double precision :: xtab(100),wx(100),sint,cost,phi,pi2

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
  if (np > 0) return

  ! --- Let np depend on x ---
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
30   enddo
20 enddo
  nxp = i
end subroutine nintsp

