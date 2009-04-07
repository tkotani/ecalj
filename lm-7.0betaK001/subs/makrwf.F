      subroutine makrwf(mode,z,rmax,l,v,a,nr,rofi,pnu,nptdif,g,gp,
     .  enu,phi,dphi,phip,dphip,p)
C- Radial wave functions and energy derivative
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit specifies boundary conditions
Ci         :0 input boundary conditions
Ci         :1 input energy 
Ci         :10s digit 
Ci         :1 set dphip to satisfy Wronskian condition, rather
Ci         :  than compute numerically
Ci   z     :nuclear charge
Ci   rmax  :augmentation radius, in a.u.
Ci   l     :l-quantum number
Ci   v     :spherical potential
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Co Outputs
Co   g     :r * wave function corresponding to b.c. pnu
Co   gp    :r * energy derivative of g; dimensioned 8*nr
Co   enu   :eigenvalue of g
Co   phi   :wave function at rmax, i.e. g/rmax
Co   dphi  :slope of wave function at rmax, i.e. (d(g/r)/dr)_rmax
Co   phip  :energy derivative of wave function at rmax
Co   dphip :energy derivative of slope of wave function at rmax
Co   p     :<gp**2> (potential parameter)
Cr Remarks
Cr   This routine makes r*phi and r*phidot, where phi and phidot
Cr   are true radial wave function and energy derivatives.
Cr   phi is normalized, and p = <phidot**2>
Cu Updates
Cu   16 Aug 04 10s digit to explicitly satisfy Wronskian
Cu   22 Dec 01 Adjustments to accomodate changes in phidx
Cu   16 May 00 New routine
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,l,nr,nptdif
      double precision a,rmax,z,rofi(1),v(nr,1),pnu(1:l+1),
     .  g(nr,2),gp(nr,2,4),phi,phip,dphi,dphip,p
C ... Local parameters
      integer konf,nn,nre,modep
      double precision b,dnu,eb1,eb2,enu,pi,slo(5),sum,tol,val(5)

      pi = 4d0*datan(1d0)
      tol = 1d-12
C     for now, rmax must match to rofi(nr)
      call fsanrg(rmax,rofi(nr),rofi(nr),1d-8,'makrwf:','rmax',.true.)
      if (mod(mode,10) .eq. 0) then
        b   = rmax/(dexp(a*nr-a)-1d0)
        konf = mod(pnu(l+1),10d0)
        dnu = dtan(pi*(.5d0-mod(pnu(l+1),10d0)))
        nn = konf-l-1
        val(1) = rmax
        slo(1) = dnu+1d0
        eb1 = -20d0
        eb2 =  20d0
        if (z .eq. 1 .and. l .gt. 2) eb2 = 100
        enu = -0.5d0
        call rseq(eb1,eb2,enu,tol,z,l,nn,val,slo,v,g,sum,a,b,rofi,nr,
     .    nre)
        val(1) = val(1)/dsqrt(sum)
        slo(1) = slo(1)/dsqrt(sum)
C        call phidot(z,l,v,enu,a,b,rofi,nr,g,val,slo,tol,nn,
C     .    gp,phi,dphi,phip,dphip,p)
        modep = 1
      else
        modep = 2
      endif

*      print *, '!! makrwf'
C     if (mode .ge. 10) modep = modep+10

      call phidx(modep,z,l,v,0d0,0d0,rofi,nr,nptdif,tol,enu,val,slo,
     .  nn,g,gp,phi,dphi,phip,dphip,p,0d0,0d0,0d0,0d0)
C     dphip = (sloi(2)-phip)/rmax

C     call prrmsh('ul',rofi,ul,nr,nr,1+lmxa)

      end

