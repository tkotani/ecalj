subroutine poinsp(z,vval,nlm,a,b,v,rofi,rho,rho0,nr,rhoves,rhves1, vnucl,vsum) !- Solves poisson Equation inside sphere
  use m_ll,only:ll
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   z     :nuclear charge
  !i   vval  :boundary conditions: potential for channel ilm = vval(ilm)
  !i   nlm   :L-cutoff for density
  !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   rofi  :radial mesh points
  !i   rho   :density, represented as sum_ilm rho(r,ilm) Y_L(ilm)
  !i   rho0  :work array, holding spherical charge density times 4*pi*r*r
  !i   nr    :number of radial mesh points
  !o Outputs
  !o   v     :electrostatic potential, satisfying boundary c. vval above
  !o         :It does not include nuclear contribution -2z/r
  !o   rhoves:electrostatic energy
  !o   rhves1:contribution to electrostatic energy from nonspherical rho
  !o   vnucl :potential at the nucleus
  !o   vsum  :integral over that potential which is zero at rmax.
  !r Remarks
  !r   Poisson's equation is d2u/dr2 = u*l*(l+1)/r**2 - 8pi*rho/r
  !u Updates
  !u   1 May 00 Adapted from nfp poinsp.f
  ! ----------------------------------------------------------------------
  implicit none
  integer :: nlm,nr
  double precision :: z,a,b,rhoves,rhves1,vnucl,vsum
  double precision :: rho(nr,nlm),vval(nlm),v(nr,nlm),rofi(nr),        vhrho(2),rho0(nr)
  integer :: nsp,ir,ilm,l,lp1
  double precision :: fpi,srfpi,y0,atepi,ea,a2b4,rpb,rmax,fllp1,df,r, drdi,srdrdi,g,x,f,y2,y3,y4,vnow,vhom,alfa,sum,wgt
  nsp = 1
  fpi = 16d0*datan(1d0)
  srfpi = dsqrt(fpi)
  y0 = 1d0/srfpi
  atepi = 2d0*fpi
  ea = dexp(a)
  a2b4 = a*a/4d0
  rpb = b
  do  5  ir = 1, nr
     rofi(ir) = rpb-b
     rpb = rpb*ea
5 enddo
  rmax = rofi(nr)
  ! --- Call poiss0 for l=0 to get good pot for small r --- ---
  do  1  ir = 1, nr
     rho0(ir) = rho(ir,1)*srfpi
1 enddo
  call poiss0(z,a,b,rofi,rho0,nr,vval(1)*y0,v,vhrho,vsum,nsp)
  rhoves = vhrho(1)
  vnucl = v(1,1)
  do  2  ir = 1, nr
     v(ir,1) = v(ir,1)*srfpi
2 enddo
  rhves1 = 0d0
  do  80  ilm = 2, nlm
     l = ll(ilm)
     lp1 = l+1
     fllp1 = l*(l+1d0)
     ! --- Numerov for inhomogeneous solution --- ---
     v(1,ilm) = 0d0
     df = 0d0
     do  12  ir = 2, 3
        r = rofi(ir)
        drdi = a*(r+b)
        srdrdi = dsqrt(drdi)
        g = (r**lp1)/srdrdi
        v(ir,ilm) = r**l
        x = fllp1*drdi*drdi/(r*r) + a2b4
        f = g*(1d0-x/12d0)
        if (ir == 2) y2 = -atepi*rho(2,ilm)*drdi*srdrdi/r
        if (ir == 3) y3 = -atepi*rho(3,ilm)*drdi*srdrdi/r
        df = f-df
12   enddo
     ir = 3
13   ir = ir+1
     r = rofi(ir)
     drdi = a*(r+b)
     srdrdi = dsqrt(drdi)
     y4 = -atepi*drdi*srdrdi*rho(ir,ilm)/r
     df = df+g*x+(y4+10d0*y3+y2)/12d0
     f = f+df
     x = fllp1*drdi*drdi/(r*r) + a2b4
     g = f/(1d0-x/12d0)
     v(ir,ilm) = g*srdrdi/r
     y2 = y3
     y3 = y4
     if (ir < nr) goto 13
     ! --- Add homogeneous solution --- ---
     vnow = v(nr,ilm)
     vhom = rmax**l
     alfa = (vval(ilm)-vnow)/vhom
     sum = 0d0
     do  10  ir = 2, nr
        r = rofi(ir)
        wgt = 2*(mod(ir+1,2)+1)
        if (ir == nr) wgt = 1d0
        v(ir,ilm) = v(ir,ilm) + alfa*(r**l)
        sum = sum+wgt*(r+b)*rho(ir,ilm)*v(ir,ilm)
10   enddo
     rhves1 = rhves1 + a*sum/3d0
     v(1,ilm) = 0d0
80 enddo
  rhoves = rhoves + rhves1
end subroutine poinsp
subroutine poiss0(z,a,b,rofi,rho,nr,vhrmax,v,rhovh,vsum,nsp)!- Hartree potential for spherical rho.
  !i Inputs:
  !i   z     :nuclear charge
  !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   rho   :spherical charge density times 4*pi*r*r
  !i   nr    :number of mesh points
  !i   vhrmax:value of v(nr)
  !i   nsp   :=1 spin degenerate, =2 non-degenerate
  !o Outputs:
  !o   rofi  :radial mesh points
  !o   v     :spherical Hartree potential
  !o   vsum  :integral over that potential which is zero at rmax.
  !o   rhovh :integral of rho*vh, where vh is the electrostatic
  !o         :potential from both electronic and nuclear charge
  !r Remarks:
  !r    Solves Poisson's equation for given spherical charge density
  !r    and a specified value vhrmax at rofi(nr).
  !r    rho =  4*pi*r*r*rhotrue  but  v = vtrue
  implicit double precision (a-h,o-z)
  integer :: nr,nsp
  double precision :: rofi(nr),v(nr,nsp),rho(nr,nsp),rhovh(nsp),pi,ea,rpb,a,b,z
  double precision :: rmax,r2,r3,r4,f2,f3,f4,x23,x34,cc,bb,dd,df,drdi,r,srdrdi,g,f,y2,y3,y4,ro,vnow,wgt,a2b4,vhrmax,vsum,vhat0
  integer :: ir,isp
  real(8):: qtotxxx
  pi = 4d0*datan(1d0)
  ea = dexp(a)
  rpb = b
  do  5  ir = 1, nr
     rofi(ir) = rpb - b
     rpb = rpb*ea
5 enddo
  rmax = rofi(nr)
  ! --- Approximate rho/r**2 by cc + bb*r + dd*r*r near zero  ---
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
  vhint=0d0
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
