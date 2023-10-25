subroutine radmsh(rmax,a,nr,rofi)   !- Makes shifted logarithmic mesh
  implicit none
  integer :: nr,ir
  double precision :: rmax,a,rofi(nr),b,rpb,ea
  if (mod(nr,2) == 0) call rx('radwgt: nr is even')
  b = rmax / (dexp(a*nr-a)-1d0)
  ea = dexp(a)
  rpb = b
  do    ir = 1, nr
     rofi(ir) = rpb - b
     rpb = rpb*ea
  enddo
end subroutine radmsh
subroutine radwgt(rmax,a,nr,wt)  !- Makes weights for numerical integration on shifted log mesh
  !  Integral f(r) dr = sum_i f_i * wt_i (Simpson rule)
  !  Thus sum_i wt_i = rmax
  !     implicit none
  integer :: nr,ir
  double precision :: rmax,a,wt(nr),b,xx,ea
  if (mod(nr,2) == 0) call rx('radwgt: nr is even')
  b = rmax / (dexp(a*nr-a)-1d0)
  xx = 2d0*a*b/3d0
  ea = dexp(a)
  do   ir = 1, nr
     wt(ir) = xx
     xx = xx*ea
  enddo
  do    ir = 2, nr-1, 2
     wt(ir) = 2d0*wt(ir)
  enddo
  wt(1) = wt(1)/2
  wt(nr) = wt(nr)/2
  wt(1)=0d0 !to avoid skipping ir=1 for cases. 2023-jan
end subroutine radwgt
subroutine radmwt(opt,rmax,a,nr,rofi,wt)  !- Makes mesh and weights for numerical integration on shifted log mesh
  !i   opt   :0 for uniform weight
  !i         :1 for r^2 weight
  !i   rmax  :augmentation radius, in a.u.,
  !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   nr    :number of radial mesh points
  !o Outputs
  !o   rofi  :radial mesh points for shifted log mesh
  !o   wt    :mesh weights; see Remarks
  !l Local variables
  !l         :
  !r Remarks
  !r  Integral for opt=0 : f(r) dr = sum_i f_i * wt_i (Simpson rule)
  !r  Integral for opt=1 : f(r) dr = sum_i f_i * wt_i r_i^2
  !r  Thus sum_i wt_i = rmax (opt=0), or rmax^3/3 (opt=1)
  !u Updates
  !u   22 Feb 03 First created
  implicit none
  integer :: nr,opt
  double precision :: rmax,a,rofi(nr),wt(nr)
  integer :: ir
  double precision :: b,xx,ea,rpb
  if (mod(nr,2) == 0) call rx('radwgt: nr is even')
  b = rmax / (dexp(a*nr-a)-1d0)
  xx = 2d0*a*b/3d0
  ea = dexp(a)
  rpb = b
  do   ir = 1, nr
     rofi(ir) = rpb - b
     rpb = rpb*ea
     wt(ir) = xx
     xx = xx*ea
     if (opt == 1) wt(ir) = wt(ir)*rofi(ir)**2
  enddo
  do  ir = 2, nr-1, 2
     wt(ir) = 2d0*wt(ir)
  enddo
  wt(1) = wt(1)/2
  wt(nr) = wt(nr)/2
end subroutine radmwt
subroutine radsum(nrx,nr,nlml,nsp,wt,rho,sum)   !- Numerical integration of a function on a shifted log mesh
  integer :: nrx,nr,nlml,nsp
  double precision :: wt(nr),rho(nrx,nlml,nsp),sum,ddot
  sum = ddot(nr,wt,1,rho,1)
  if (nsp == 2) sum = sum + ddot(nr,wt,1,rho(1,1,2),1)
end subroutine radsum
subroutine radgra(a,b,nr,r,v,gradv)  !- radial gradient
  !i    a,b,nr,r(nr),v(nr)
  !o Outputs:
  !o    gradv(nr)
  !r Remarks:
  !r    makes the derivative of the function v defined in a mesh of 
  !r                r(i)=B (exp A(i-1)-1)
  ! ------------------------------------------------------------------
  implicit none
  integer :: nr
  double precision :: a,b,r(nr),v(nr),gradv(nr)
  integer :: nm2,i
  ! Forward diffs for first and second point (Handbook,25.3.9 with 25.1.1)
  gradv(1) = ((6d0*v(2)+20d0/3d0*v(4)+1.2d0*v(6)) -(2.45d0*v(1)+7.5d0*v(3)+3.75d0*v(5)+1d0/6d0*v(7)))/a
  gradv(2) = ((6d0*v(3)+20d0/3d0*v(5)+1.2d0*v(7)) -(2.45d0*v(2)+7.5d0*v(4)+3.75d0*v(6)+1d0/6d0*v(8)))/a
  ! Five points formula  (25.3.6)
  nm2 = nr-2
  do i = 3, nm2
     gradv(i) = ((v(i-2)+8*v(i+1)) - (8*v(i-1)+v(i+2)))/12/a
  enddo
  ! Five points formula  (25.3.6)
  gradv(nr-1) = (-1d0/12d0*v(nr-4)+0.5d0*v(nr-3)-1.5d0*v(nr-2)+5d0/6d0*v(nr-1)+0.25d0*v(nr))/a
  gradv(nr) = (0.25d0*v(nr-4)-4d0/3d0*v(nr-3)+3d0*v(nr-2)      -4d0*v(nr-1)+25d0/12d0*v(nr))/a
  ! Three points formula  (25.3.4)
  !     gradv(nr-1)=(v(nr)-v(nr-2))/2d0/a
  !     gradv(nr)=(v(nr-2)/2d0-2d0*v(nr-1)+1.5d0*v(nr))/a
  do  i = 1, nr
     gradv(i) = gradv(i)/(r(i)+b)
  enddo
end subroutine radgra
