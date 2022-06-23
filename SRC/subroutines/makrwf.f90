subroutine makrwf(mode,z,rmax,l,v,a,nr,rofi,pnu,nptdif,g,gp, enu,phi,dphi,phip,dphip,p)
  !- Radial wave functions and energy derivative
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :1s digit specifies boundary conditions
  !i         :0 input boundary conditions
  !i         :1 input energy
  !i         :10s digit
  !i         :1 set dphip to satisfy Wronskian condition, rather
  !i         :  than compute numerically
  !i   z     :nuclear charge
  !i   rmax  :augmentation radius, in a.u.
  !i   l     :l-quantum number
  !i   v     :spherical potential
  !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   nr    :number of radial mesh points
  !i   rofi  :radial mesh points
  !i   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
  !i          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
  !o Outputs
  !o   g     :r * wave function corresponding to b.c. pnu
  !o   gp    :r * energy derivative of g; dimensioned 8*nr
  !o   enu   :eigenvalue of g
  !o   phi   :wave function at rmax, i.e. g/rmax
  !o   dphi  :slope of wave function at rmax, i.e. (d(g/r)/dr)_rmax
  !o   phip  :energy derivative of wave function at rmax
  !o   dphip :energy derivative of slope of wave function at rmax
  !o   p     :<gp**2> (potential parameter)
  !r Remarks
  !r   This routine makes r*phi and r*phidot, where phi and phidot
  !r   are true radial wave function and energy derivatives.
  !r   phi is normalized, and p = <phidot**2>
  !u Updates
  !u   16 Aug 04 10s digit to explicitly satisfy Wronskian
  !u   22 Dec 01 Adjustments to accomodate changes in phidx
  !u   16 May 00 New routine
  ! ----------------------------------------------------------------------
  implicit none
  integer :: mode,l,nr,nptdif, konf,nn,nre,modep
  real(8):: a,rmax,z,rofi(1),v(nr,1),pnu(1:l+1),g(nr,2),gp(nr,2,4),phi,phip,dphi,dphip,p,&
       b,dnu,eb1,eb2,enu,slo(5),sum,val(5)
  real(8),parameter:: pi = 4d0*datan(1d0), tol = 1d-12
  call fsanrg(rmax,rofi(nr),rofi(nr),1d-8,'makrwf:','rmax',.true.)! rmax must match to rofi(nr)
  if (mod(mode,10) == 0) then
     b   = rmax/(dexp(a*nr-a)-1d0)
     konf = mod(pnu(l+1),10d0)
     dnu = dtan(pi*(.5d0-mod(pnu(l+1),10d0)))
     nn = konf-l-1
     val(1) = rmax
     slo(1) = dnu+1d0
     eb1 = -20d0
     eb2 =  20d0
     if (z == 1 .AND. l > 2) eb2 = 100
     enu = -0.5d0
     call rseq(eb1,eb2,enu,tol,z,l,nn,val,slo,v,g,sum,a,b,rofi,nr,  nre)
     val(1) = val(1)/dsqrt(sum)
     slo(1) = slo(1)/dsqrt(sum)
     modep = 1
  else
     modep = 2
  endif
  call phidx(modep,z,l,v,0d0,0d0,rofi,nr,nptdif,tol,enu,val,slo, &
       nn,g,gp,phi,dphi,phip,dphip,p,0d0,0d0,0d0,0d0)
end subroutine makrwf
