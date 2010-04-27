      subroutine hanszd(mode,r,e,rsm,lmax,hs,dhs,ddhs,hsp,dhsp,ddhsp)
C- Value and some derivatives of smoothed radial Hankel (e < 0) 
C  or Neumann functions (e > 0). 
C ---------------------------------------------------------------
Ci Inputs
Ci   mode  :tells hansmd what derivatives to make.
Ci         :1s digit concerns 2nd radial derivative
Ci         :0 make neither 1st or 2nd radial derivative.
Ci         :>0 make 1st and second radial derivative:
Ci         :1 ddhs = radial part of Laplacian, 1/r d^2 (r*h_l) / dr^2 
Ci         :2 ddhs = 1/r d^2 (r*h_l) / dr^2  - l(l+1)/r^2 h_l
Ci         :  NB: ddhs = laplacian of 3-dimensional hs_l YL
Ci         :3 ddhs = d^2 (h_l) / dr^2
Ci         :10s digit concerns energy derivative
Ci         :0 make none of hsp,dhsp,ddhsp
Ci         :1 and 1s digit = 0 make hsp
Ci         :  and 1s digit = 2 make all of hsp,dhsp,ddhsp
Ci   r     :radius
Ci   e     :hankel/neumann energy
Ci   rsm   :hankel/neumann smoothing radius
Ci   lmax  :make function values for l between (0:lmax)
Co  Outputs:
Co   hs    : function values of the radial sm. hankel (e < 0) 
Co         : or neumann (e > 0) function h(e,r,0:lmax)
Co         : A solid hankel/neumann is H=h*Y_L, where Y_L are the 
Co         : spherical harmonics for unit radius (no scaling by r**l)
Co   dhs   : radial derivative of hs
Co   ddhs  : radial part of Laplacian of hs, i.e. 1/r d^2 (r h) /dr^2
Co         : OR some other second derivative (see mode)
Co   hsp   : energy derivative of hs
Co   dhsp  : mixed energy + radial derivative of hs
Co   ddhsp : 3-d laplacian of hsp YL
Cr Remarks
Cr  See J. Math. Phys. 39, 3393 (1998).
Cr    For radial derivative, see JMP 39, 3393, Eq. 4.7
Cr      h'  = l/r h_l - h_l+1
Cr    Second radial derivative:
Cr      h'' = l(l-1)/r^2 xi_l - (2l+1)/r xi_l+1 + xi_l+2
Cr      1/r d^2/dr^2 (r*h_l) = l(l+1)/r^2 h_l - (2l+3)/r h_l+1 + h_l+2
Cr    Energy derivative:  see JMP 39, 3393, Eq. 7.5
Cr      hp_l = r/2 h_l-1  Special case l=0: hp_0 = 1/2 h_-1 ?
Cr    Mixed energy + radial derivative:
Cr      hp'_l = h(l-1)*l/2 - h(l)*r/2
Cr    Mixed energy + kinetic energy
Cr      hp'' = -(2l+3)/2 h_l + h_l+1*r/2
Cr
Cr  Note connection with hansmr and hansmz, which make xi(l) = h(l)/r^l
Cu Updates
Cu   29 Feb 08 (S. Lozovoi) adapted from hansmd.f to handle e > 0 
C ---------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,lmax
      double precision r,e,rsm
      double precision hs(0:lmax),dhs(0:lmax),ddhs(0:lmax)
      double precision hsp(0:lmax),dhsp(0:lmax),ddhsp(0:lmax)
C ... Local parameters
      integer l,mode0,mode1,n0
      parameter (n0 = 12)
      double precision xi(-1:n0),fi(-1:n0),rsml(-1:n0)

      if (rsm .lt. 0) call rx('hanszd: negative rsm') 
      if (lmax+2 .gt. n0) call rx('hanszd: increase parameter n0') 
      if (lmax .lt. 0) return

      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      if (mode1 .gt. 0 .and. (mode0 .eq. 1 .or. mode0 .ge. 3)) 
     .  call rx('hanszd: for energy derivatives mode0 should be '//
     .          'either 0 or 2') 

      rsml(0) = -rsm
      call hansrz(rsml,-1,lmax+2,e,r*r,1,1,11,xi,fi)

      do l = 0, lmax
        hs(l) = xi(l)
        if (mode0 .ne. 0) then
            dhs(l)  = xi(l)*l/r - xi(l+1)
          if (mode0 .eq. 1)
     .      ddhs(l) = xi(l)*l*(l+1)/r**2 - (2*l+3)/r*xi(l+1) + xi(l+2)
          if (mode0 .eq. 2)
     .      ddhs(l) =                    - (2*l+3)/r*xi(l+1) + xi(l+2)
          if (mode0 .eq. 3)
     .      ddhs(l) = xi(l)*l*(l-1)/r**2 - (2*l+1)/r*xi(l+1) + xi(l+2)
        endif
        if (mode1 .ne. 0) then
          hsp(l)   = xi(l-1)*r/2
          if (mode0 .ne. 0) then
            dhsp(l)  = (xi(l-1)*l - xi(l)*r)/2
            ddhsp(l) = - (2*l+3)*xi(l)/2 + xi(l+1)*r/2
          endif
        endif
      enddo

c The lines below are not needed since xi(-1) is now scaled inside hansrz 
c     if (mode1 .ne. 0) then
c       hsp(0) = xi(-1)/2
c     endif

      end

