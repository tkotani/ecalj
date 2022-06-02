subroutine radmsh(rmax,a,nr,rofi)
  !- Makes shifted logarithmic mesh
  !     implicit none
  integer :: nr,ir
  double precision :: rmax,a,rofi(nr),b,rpb,ea
  if (mod(nr,2) == 0) call rx('radwgt: nr is even')
  b = rmax / (dexp(a*nr-a)-1d0)
  ea = dexp(a)
  rpb = b
  do  1  ir = 1, nr
     rofi(ir) = rpb - b
     rpb = rpb*ea
1 enddo
end subroutine radmsh
subroutine radwgt(rmax,a,nr,wt)
  !- Makes weights for numerical integration on shifted log mesh
  !  Integral f(r) dr = sum_i f_i * wt_i (Simpson rule)
  !  Thus sum_i wt_i = rmax
  !     implicit none
  integer :: nr,ir
  double precision :: rmax,a,wt(nr),b,xx,ea
  if (mod(nr,2) == 0) call rx('radwgt: nr is even')
  b = rmax / (dexp(a*nr-a)-1d0)
  xx = 2d0*a*b/3d0
  ea = dexp(a)
  do  1  ir = 1, nr
     wt(ir) = xx
     xx = xx*ea
1 enddo
  do  2  ir = 2, nr-1, 2
     wt(ir) = 2d0*wt(ir)
2 enddo
  wt(1) = wt(1)/2
  wt(nr) = wt(nr)/2
end subroutine radwgt
subroutine radmwt(opt,rmax,a,nr,rofi,wt)
  !- Makes mesh and weights for numerical integration on shifted log mesh
  ! ----------------------------------------------------------------------
  !i Inputs
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
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: nr,opt
  double precision :: rmax,a,rofi(nr),wt(nr)
  ! ... Local parameters
  integer :: ir
  double precision :: b,xx,ea,rpb
  if (mod(nr,2) == 0) call rx('radwgt: nr is even')
  b = rmax / (dexp(a*nr-a)-1d0)
  xx = 2d0*a*b/3d0
  ea = dexp(a)
  rpb = b
  do  1  ir = 1, nr
     rofi(ir) = rpb - b
     rpb = rpb*ea
     wt(ir) = xx
     xx = xx*ea
     if (opt == 1) wt(ir) = wt(ir)*rofi(ir)**2
1 enddo
  do  2  ir = 2, nr-1, 2
     wt(ir) = 2d0*wt(ir)
2 enddo
  wt(1) = wt(1)/2
  wt(nr) = wt(nr)/2
end subroutine radmwt

subroutine radsum(nrx,nr,nlml,nsp,wt,rho,sum)
  !- Numerical integration of a function on a shifted log mesh
  !u   19 Jun 00 added extra arguments
  !     implicit none
  integer :: nrx,nr,nlml,nsp
  double precision :: wt(nr),rho(nrx,nlml,nsp),sum,ddot

  sum = ddot(nr,wt,1,rho,1)
  if (nsp == 2) sum = sum + ddot(nr,wt,1,rho(1,1,2),1)
end subroutine radsum

subroutine radext(mode,nr,nrx,fac,a,rmax,nrbig,rbig,rofi,rwgt)
  !- Find radius, mesh suitable for extending orbitals outside MT sphere
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :1s digit
  !i         :0  nrbig,rbig are input; do not make them
  !i         :1  set rbig = smaller of   rofi(nrx)  and  fac*rmax
  !i         :10s digit
  !i         :if nonzero, make rofi and rwgt
  !i   nr    :number of radial mesh points on regular mesh
  !i   nrx   :maximum allowed number of radial mesh points
  !i   fac   :approximate factor to scale rmax, rbig ~ fac*rmax
  !i         :NB: true factor is constrained because rbig must
  !i         :conform to radial mesh specified by (rmax,a,nr)
  !i   a     :mesh points are given by
  !i         :rofi(i) = rmax [e^(a(i-1))-1] / [e^(a(nr-1))-1]
  !i   rmax  :augmentation radius, in a.u.,
  !o Outputs
  ! o  nrbig :number of points on extended mesh.
  ! o        :NB: nrbig is input if 1s digit mode=0
  ! o        :In the latter case, nrbig must be consistent with the mesh
  ! o        :points specified by (a,nr,rmax) and also rbig.
  ! o  rbig  :sphere radius of extended mesh
  ! o        :NB: rbig is input if 1s digit mode=0
  !o   rofi  :(10s digit mode > 0)
  !o         :radial mesh points: rofi(1..nrbig) will be generated
  !o         :rofi(nrbig) is rmax for extended mesh
  !o   rwgt  :(10s digit mode > 0)
  !o         :radial mesh weights: rwgt(1..nrbig) will be generated
  !o         :rwgt is actually designed for two integration radii:
  !o         :int(0,rmax) = I(1..nr) and int(rmax,rbig) = I(nr..nrbig).
  !o         :Integral int(1..nrbig) must be done in two steps, by summing
  !o         :I(1..nr) and I(nr..nrbig)
  !r Remarks
  !r
  !u Updates
  !u   24 Sep 04 First created
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: mode,nr,nrx,nrbig
  double precision :: rmax,fac,rbig,a,rofi(*),rwgt(*)
  ! ... Local parameters
  integer :: idn

  if (mod(mode,10) == 1) then
     rbig = rmax * (dexp(a*nrx-a)-1d0)/(dexp(a*nr-a)-1d0)
     !     If rbig>fac*rmax, estimate from exp((nrbig-nr)a) = fac
     if (rbig > fac*rmax) then
        idn = dlog(fac)/a
        if (mod(idn,2) == 1) idn = idn-1
        nrbig = min(nr+idn,nrx)
        rbig = rmax * (dexp(a*nrbig-a)-1d0)/(dexp(a*nr-a)-1d0)
     endif
  endif

  ! --- Points and weights on extended mesh ---
  if (mod(mode/10,10) /= 0) then
     call radmsh(rbig,a,nrbig,rofi)
     call radwgt(rbig,a,nrbig,rwgt)
     if (nr < nrbig) rwgt(nr) = rwgt(nr)/2
  endif
end subroutine radext

