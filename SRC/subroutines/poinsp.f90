subroutine poinsp(z,vval,nlm,a,b,v,rofi,rho,rho0,nr,rhoves,rhves1, &
     vnucl,vsum)
  !- Solves poisson Equation inside sphere
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
  !     implicit none
  ! ... Passed parameters
  integer :: nlm,nr
  double precision :: z,a,b,rhoves,rhves1,vnucl,vsum
  double precision :: rho(nr,nlm),vval(nlm),v(nr,nlm),rofi(nr), &
       vhrho(2),rho0(nr)
  ! ... Local parameters
  integer :: nsp,ir,ilm,l,ll,lp1
  double precision :: fpi,srfpi,y0,atepi,ea,a2b4,rpb,rmax,fllp1,df,r, &
       drdi,srdrdi,g,x,f,y2,y3,y4,vnow,vhom,alfa,sum,wgt

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

