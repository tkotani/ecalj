subroutine poiss0(z,rofi,rho,nr,vhrmax, v,rhovh,vsum,nsp)!- Hartree potential for spherical rho.
  !i Inputs:
  !i   z     :nuclear charge
  !i   rho   :spherical charge density times 4*pi*r*r
  !i   nr    :number of mesh points
  !i   vhrmax:value of v(nr)
  !i   nsp   :=1 spin degenerate, =2 non-degenerate
  !i   rofi  :radial mesh points
  !o Outputs:
  !o   v     :spherical Hartree potential
  !NOTE: v do not include -2z/r
  !o   vsum  :integral over that potential which is zero at rmax.
  !o   rhovh :integral of rho*vh, where vh is the electrostatic
  !o         :potential from both electronic and nuclear charge
  !r Remarks:
  !r    Solves Poisson's equation for given spherical charge density
  !r    and a specified value vhrmax at rofi(nr).
  !r    rho =  4*pi*r*r*rhotrue  but  v = vtrue
  implicit none !double precision (a-h,o-z)
  integer :: nr,nsp
  real(8) :: rofi(nr),v(nr,nsp),rho(nr,nsp),rhovh(nsp),z,a,b
  real(8) :: rmax,r2,r3,r4,f2,f3,f4,x23,x34,cc,bb,dd,df,drdi,r,srdrdi,g,f,y2,y3,y4,ro,vnow,wgt,a2b4,vhrmax,vsum,vhat0
  real(8),parameter:: pi = 4d0*datan(1d0)
  integer :: ir,isp
  ! --- Approximate rho/r**2 by cc + bb*r + dd*r*r near zero  ---
  rmax = rofi(nr)
  r2 = rofi(2)
  r3 = rofi(3)
  r4 = rofi(4)
  f2 = 0d0
  f3 = 0d0
  f4 = 0d0
  do  75  isp = 1, nsp
     f2 = f2 + rho(2,isp)/r2**2
     f3 = f3 + rho(3,isp)/r3**2
     f4 = f4 + rho(4,isp)/r4**2
75 enddo
  x23 = (r3*r3*f2 - r2*r2*f3)/(r3 - r2)
  x34 = (r4*r4*f3 - r3*r3*f4)/(r4 - r3)
  cc = (r2*x34 - r4*x23) / (r3*(r2 - r4))
  bb = ((r2+r3)*x34 - (r3+r4)*x23) / (r3*r3*(r4-r2))
  dd = (f2 - bb*r2 - cc)/r2**2
  ! --- Numerov for inhomogeneous solution ---
  a = log(rofi(3)/rofi(2)-1d0) !rofi(i)=b(e^(a*(i-1))-1). We get a and b
  b= rofi(3)/(exp(2d0*a)-1d0)
  a2b4 = a*a/4d0
  v(1,1) = 1d0
  df = 0d0
  do  12  ir = 2, 3
     r = rofi(ir)
     drdi = a*(r + b)
     srdrdi = dsqrt(drdi)
     v(ir,1) = v(1,1) - r*r*(cc/3d0 + r*bb/6d0 + r*r*dd/10d0)
     g = v(ir,1)*r/srdrdi
     f = g*(1d0 - a2b4/12d0)
     if (ir == 2) y2 = -2d0*f2*r2*drdi*srdrdi
     if (ir == 3) y3 = -2d0*f3*r3*drdi*srdrdi
     df = f - df
12 enddo
  do 13  ir = 4, nr
     r = rofi(ir)
     drdi = a*(r + b)
     srdrdi = dsqrt(drdi)
     ro = 0d0
     do  76  isp = 1, nsp
        ro = ro + rho(ir,isp)
76   enddo
     y4 = -2d0*drdi*srdrdi*ro/r
     df = df + g*a2b4 + (y4 + 10d0*y3 + y2)/12d0
     f = f + df
     g = f/(1d0 - a2b4/12d0)
     v(ir,1) = g*srdrdi/r
     y2 = y3
     y3 = y4
13 enddo
  ! --- Add constant to get v(nr) = vhrmax ---
  vnow = v(nr,1) - 2d0*z/rmax
  do  27  ir = 1, nr
     v(ir,1) = v(ir,1) + (vhrmax - vnow)
27 enddo
  ! --- Integral vh and  rho * vh ---
  rhovh(1:nsp) = 0d0
  vsum = 0d0
  vhat0 = 0d0
  do  30  ir = 2, nr
     r = rofi(ir)
     drdi = a*(r + b)
     wgt = 2*(mod(ir+1,2) + 1)/3d0
     if (ir == nr) wgt = 1d0/3d0
     ro = 0d0
     do  31  isp = 1, nsp
        rhovh(isp)= rhovh(isp)+ wgt*drdi*rho(ir,isp)*(v(ir,1) - 2d0*z/r)
        ro = ro + rho(ir,isp)
31   enddo
     vhat0 = vhat0 + wgt*drdi*ro*(1d0/r - 1d0/rmax)
     vsum = vsum + wgt*drdi*r*r*v(ir,1)
30 enddo
  vsum  = 4d0*pi*(vsum - z*rmax*rmax)
  vhat0 = 2d0*vhat0 + 2d0*z/rmax + vhrmax
  v(1,1) = vhat0
  if(nsp == 1) return
  v(:,2)=v(:,1) ! if spin polarized ---
end subroutine poiss0
