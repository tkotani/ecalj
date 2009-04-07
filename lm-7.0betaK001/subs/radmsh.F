      subroutine radmsh(rmax,a,nr,rofi)
C- Makes shifted logarithmic mesh
C     implicit none
      integer nr,ir
      double precision rmax,a,rofi(nr),b,rpb,ea
      if (mod(nr,2).eq.0) call rx('radwgt: nr is even')
      b = rmax / (dexp(a*nr-a)-1d0)
      ea = dexp(a)
      rpb = b
      do  1  ir = 1, nr
        rofi(ir) = rpb - b
        rpb = rpb*ea
    1 continue
      end
      subroutine radwgt(rmax,a,nr,wt)
C- Makes weights for numerical integration on shifted log mesh
C  Integral f(r) dr = sum_i f_i * wt_i (Simpson rule)
C  Thus sum_i wt_i = rmax
C     implicit none
      integer nr,ir
      double precision rmax,a,wt(nr),b,xx,ea
      if (mod(nr,2).eq.0) call rx('radwgt: nr is even')
      b = rmax / (dexp(a*nr-a)-1d0)
      xx = 2d0*a*b/3d0
      ea = dexp(a)
      do  1  ir = 1, nr
      wt(ir) = xx
    1 xx = xx*ea
      do  2  ir = 2, nr-1, 2
    2 wt(ir) = 2d0*wt(ir)
      wt(1) = wt(1)/2
      wt(nr) = wt(nr)/2
      end
      subroutine radmwt(opt,rmax,a,nr,rofi,wt)
C- Makes mesh and weights for numerical integration on shifted log mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :0 for uniform weight
Ci         :1 for r^2 weight
Ci   rmax  :augmentation radius, in a.u.,
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Co Outputs
Co   rofi  :radial mesh points for shifted log mesh
Co   wt    :mesh weights; see Remarks
Cl Local variables
Cl         :
Cr Remarks
Cr  Integral for opt=0 : f(r) dr = sum_i f_i * wt_i (Simpson rule)
Cr  Integral for opt=1 : f(r) dr = sum_i f_i * wt_i r_i^2
Cr  Thus sum_i wt_i = rmax (opt=0), or rmax^3/3 (opt=1)
Cu Updates
Cu   22 Feb 03 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nr,opt
      double precision rmax,a,rofi(nr),wt(nr)
C ... Local parameters
      integer ir
      double precision b,xx,ea,rpb
      if (mod(nr,2).eq.0) call rx('radwgt: nr is even')
      b = rmax / (dexp(a*nr-a)-1d0)
      xx = 2d0*a*b/3d0
      ea = dexp(a)
      rpb = b
      do  1  ir = 1, nr
        rofi(ir) = rpb - b
        rpb = rpb*ea
        wt(ir) = xx
        xx = xx*ea
        if (opt .eq. 1) wt(ir) = wt(ir)*rofi(ir)**2
    1 continue
      do  2  ir = 2, nr-1, 2
    2 wt(ir) = 2d0*wt(ir)
      wt(1) = wt(1)/2
      wt(nr) = wt(nr)/2
      end

      subroutine radsum(nrx,nr,nlml,nsp,wt,rho,sum)
C- Numerical integration of a function on a shifted log mesh
Cu   19 Jun 00 added extra arguments
C     implicit none
      integer nrx,nr,nlml,nsp
      double precision wt(nr),rho(nrx,nlml,nsp),sum,ddot

      sum = ddot(nr,wt,1,rho,1)
      if (nsp .eq. 2) sum = sum + ddot(nr,wt,1,rho(1,1,2),1)
      end

      subroutine radext(mode,nr,nrx,fac,a,rmax,nrbig,rbig,rofi,rwgt)
C- Find radius, mesh suitable for extending orbitals outside MT sphere
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :0  nrbig,rbig are input; do not make them
Ci         :1  set rbig = smaller of   rofi(nrx)  and  fac*rmax
Ci         :10s digit
Ci         :if nonzero, make rofi and rwgt
Ci   nr    :number of radial mesh points on regular mesh
Ci   nrx   :maximum allowed number of radial mesh points
Ci   fac   :approximate factor to scale rmax, rbig ~ fac*rmax
Ci         :NB: true factor is constrained because rbig must
Ci         :conform to radial mesh specified by (rmax,a,nr)
Ci   a     :mesh points are given by
Ci         :rofi(i) = rmax [e^(a(i-1))-1] / [e^(a(nr-1))-1]
Ci   rmax  :augmentation radius, in a.u.,
Co Outputs
Cio  nrbig :number of points on extended mesh.
Cio        :NB: nrbig is input if 1s digit mode=0
Cio        :In the latter case, nrbig must be consistent with the mesh
Cio        :points specified by (a,nr,rmax) and also rbig.
Cio  rbig  :sphere radius of extended mesh
Cio        :NB: rbig is input if 1s digit mode=0
Co   rofi  :(10s digit mode > 0)
Co         :radial mesh points: rofi(1..nrbig) will be generated
Co         :rofi(nrbig) is rmax for extended mesh
Co   rwgt  :(10s digit mode > 0)
Co         :radial mesh weights: rwgt(1..nrbig) will be generated
Co         :rwgt is actually designed for two integration radii:
Co         :int(0,rmax) = I(1..nr) and int(rmax,rbig) = I(nr..nrbig).
Co         :Integral int(1..nrbig) must be done in two steps, by summing
Co         :I(1..nr) and I(nr..nrbig)
Cr Remarks
Cr
Cu Updates
Cu   24 Sep 04 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nr,nrx,nrbig
      double precision rmax,fac,rbig,a,rofi(*),rwgt(*)
C ... Local parameters
      integer idn
      
      if (mod(mode,10) .eq. 1) then
      rbig = rmax * (dexp(a*nrx-a)-1d0)/(dexp(a*nr-a)-1d0)
C     If rbig>fac*rmax, estimate from exp((nrbig-nr)a) = fac
      if (rbig .gt. fac*rmax) then
        idn = dlog(fac)/a
        if (mod(idn,2) .eq. 1) idn = idn-1
        nrbig = min(nr+idn,nrx)
        rbig = rmax * (dexp(a*nrbig-a)-1d0)/(dexp(a*nr-a)-1d0)
      endif
      endif

C --- Points and weights on extended mesh ---
      if (mod(mode/10,10) .ne. 0) then
        call radmsh(rbig,a,nrbig,rofi)
        call radwgt(rbig,a,nrbig,rwgt)
        if (nr .lt. nrbig) rwgt(nr) = rwgt(nr)/2
      endif
      end

