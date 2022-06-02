subroutine hansmd(mode,r,e,rsm,lmax,hs,dhs,ddhs,hsp,dhsp,ddhsp)
  !- Value and some derivatives of smoothed radial Hankel functions
  ! ---------------------------------------------------------------
  !i Inputs
  !i   mode  :tells hansmd what derivatives to make.
  !i         :1s digit concerns 2nd radial derivative
  !i         :0 make neither 1st or 2nd radial derivative.
  !i         :>0 make 1st and second radial derivative:
  !i         :1 ddhs = radial part of Laplacian, 1/r d^2 (r*h_l) / dr^2
  !i         :2 ddhs = 1/r d^2 (r*h_l) / dr^2  - l(l+1)/r^2 h_l
  !i         :  NB: ddhs = laplacian of 3-dimensional hs_l YL
  !i         :3 ddhs = d^2 (h_l) / dr^2
  !i         :1s digit concerns energy derivative
  !i         :0 make none of hsp,dhsp,ddhsp
  !i         :1 make all  of hsp,dhsp,ddhsp
  !i   r     :radius
  !i   e     :hankel energy
  !i   rsm   :hankel smoothing radius
  !i   lmax  :make function values for l between (0:lmax)
  !o  Outputs:
  !o   hs    : function values of the radial sm. hankel h(e,r,0:lmax)
  !o         : A solid hankel is H=h*Y_L, where Y_L are the spherical
  !o         : harmonics for unit radius (no scaling by r**l)
  !o   dhs   : radial derivative of hs
  !o   ddhs  : radial part of Laplacian of hs, i.e. 1/r d^2 (r h) /dr^2
  !o         : OR some other second derivative (see mode)
  !o   hsp   : energy derivative of hs
  !o   dhsp  : mixed energy + radial derivative of hs
  !o   ddhsp : 3-d laplacian of hsp YL
  !r Remarks
  !r  See J. Math. Phys. 39, 3393 (1998).
  !r    For radial derivative, see JMP 39, 3393, Eq. 4.7
  !r      h'  = l/r h_l - h_l+1
  !r    Second radial derivative:
  !r      h'' = l(l-1)/r^2 xi_l - (2l+1)/r xi_l+1 + xi_l+2
  !r      1/r d^2/dr^2 (r*h_l) = l(l+1)/r^2 h_l - (2l+3)/r h_l+1 + h_l+2
  !r    Energy derivative:  see JMP 39, 3393, Eq. 7.5
  !r      hp_l = r/2 h_l-1  Special case l=0: hp_0 = 1/2 h_-1 ?
  !r    Mixed energy + radial derivative:
  !r      hp'_l = h(l-1)*l/2 - h(l)*r/2
  !r    Mixed energy + kinetic energy
  !r      hp'' = -(2l+3)/2 h_l + h_l+1*r/2
  !r
  !r  Note connection with hansmr, which makes xi(l) = h(l) / r^l
  !u Updates
  !u   28 Aug 04 Also generate ddhsp
  !u   16 Jun 04 First created
  ! ---------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: mode,lmax
  double precision :: r,e,rsm
  double precision :: hs(0:lmax),dhs(0:lmax),ddhs(0:lmax)
  double precision :: hsp(0:lmax),dhsp(0:lmax),ddhsp(0:lmax)
  ! ... Local parameters
  integer :: idx,l,mode0,mode1
  double precision :: xi(-1:lmax+2) !,wk(2)

  if (lmax < 0) return
  mode0 = mod(mode,10)
  mode1 = mod(mode/10,10)
  !      call hansr(rsm,-1,lmax+2,1,lmax+2,e,r**2,1,1,idx,wk,11,xi)
  call hansr(rsm,-1,lmax+2,1,lmax+2,e,r**2,1,1,idx,11,xi)
  do  54  l = 0, lmax
     hs(l)   = xi(l)
     if (mode0 /= 0) then
        dhs(l)  = xi(l)*l/r - xi(l+1)
        if (mode0 == 1) &
             ddhs(l) = xi(l)*l*(l+1)/r**2 - (2*l+3)/r*xi(l+1) + xi(l+2)
        if (mode0 == 2) &
             ddhs(l) =                    - (2*l+3)/r*xi(l+1) + xi(l+2)
        if (mode0 == 3) &
             ddhs(l) = xi(l)*l*(l-1)/r**2 - (2*l+1)/r*xi(l+1) + xi(l+2)
     endif
     if (mode1 /= 0) then
        hsp(l)   = xi(l-1)*r/2
        dhsp(l)  = (xi(l-1)*l - xi(l)*r)/2
        ddhsp(l) = - (2*l+3)*xi(l)/2 + xi(l+1)*r/2
     endif
54 enddo
  if (mode1 /= 0) then
     hsp(0) = xi(-1)/2
  endif

end subroutine hansmd

