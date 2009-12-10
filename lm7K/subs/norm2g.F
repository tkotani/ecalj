      subroutine norm2g(a,b,kap2,g,l,liophi,nr,hcr,slo,val,rmax,rofi,
     .  phi0,phio)
C- Normalizes wave function so that phi^0(a)=1
C ----------------------------------------------------------------------
Ci Inputs:
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :                 -//-
Ci   kap2  :kinetic energy enu-vmtz
Ci   l     :angular momentum
Ci  liophi:T tabulate phi0 on the radial mesh from hcr to rmax
Ci   nr    :number of mesh points
Ci   hcr   :hard-core screening radius
Ci   slo   :slope of phi=g/r at rmax; see Remarks.
Ci   val   :value of phi=g/r at rmax; see Remarks.
Ci   rmax  :Wigner-Seitz sphere radius, in atomic units
Cio Inputs/Outputs:
Cio  g     :On input, g is phi times r
Cio        :On output, g is renormalized
Co Outputs:
Co   phio  :value of back extrapolated wavefunction at hcr
Co   phi0  :phi^0 on the radial mesh from hcr to rmax
Cr  Remarks
Cr   *phi(r) is solution to Schroedinger's equation for energy enu
Cr    in a potential v(r). 
Cr   *phi^0 is the solution to Schroedinger equation for a flat
Cr    potential vmtz with the same boundary conditions as phi.
Cr    phio is the value of phi^0 at hcr.
Cr
Cr This was adapted from the Stuttgart third-generation LMTO package.
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer l,nr
      logical liophi
      double precision a,b,kap2,g(*),hcr,slo,val,rmax,phi0(*),
     .  phio,rofi(*)
C Local variables:
      integer i,imin,iprint,lmxx
      parameter(lmxx=10)
      double precision er2,fi(0:lmxx+1),gi(0:lmxx+1),wn,wj,dphi,dn,dj,
     .                 sigl,phi,phioi,ri
C External calls:
      external bessl2,dscal

      er2 = kap2*rmax*rmax
      call bessl2(er2,0,l+1,fi(0),gi(0))
C     phi,dphi are value, slope at rmax
      phi  = val/rmax
      dphi = rmax*slo/val - 1d0
C     dj,dn are the logarithmic derivatives of bessl and hankels
      dn   = l - gi(l+1)/gi(l)*(l+l+1)
      dj   = l - fi(l+1)/fi(l)/(l+l+1)*er2
C     wj,wn are the amounts of bessl and hankel making up phi0
      wj   = (dphi-dj)*fi(l)*phi
      wn   = (dphi-dn)*gi(l)*phi
      sigl = hcr/rmax
      er2  = kap2*hcr*hcr
      call bessl2(er2,0,l+1,fi(0),gi(0))
      phio = 2d0*(wn*fi(l)*sigl**l - wj*gi(l)*sigl**(-l-1))
      if (iprint() .ge. 120) print 301, phi,phio
      call dscal(2*nr,1d0/phio,g,1)
      if (iprint() .ge. 120) print 300, kap2,val,slo,dphi,g(10),g(100)

C --- Generate phi0 on a radial mesh between hcr and rmax ---
      if (liophi) then
         phioi = 1d0/phio
         imin = int(1d0+dlog(hcr/b+1d0)/a) - 1
         do  10  i = 1, imin-1
   10    phi0(i) = 0d0
         do  20  i = imin, nr
           ri = rofi(i)
           er2 = kap2*ri*ri
           call bessl2(er2,0,l+1,fi(0),gi(0))
           phi0(i) = (2d0*(wn*fi(l)*(ri/rmax)**l -
     .                     wj*gi(l)*(ri/rmax)**(-l-1)))*phioi
   20     continue
      endif

  300 format(' NORM2G: kap2, val slo dphi g1 g2',6f10.5)
  301 format(' NORM2G: phi, phio',2f10.5)
      end

