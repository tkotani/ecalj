      subroutine rsqnri(e,l,v,z,a,rofi,i0,nr,g0,nn,val,slo,g,q)
C- Integrates a solution of the inhomogeneous Schroedinger equation
C  (H_l - e) g = g0 outward (nonrelativistic case)
C ----------------------------------------------------------------------
Ci Inputs
Ci   e     :Energy
Ci   z     :nuclear charge
Ci   l     :angular momentum
Ci   v     :spherical potential (true potential excluding nuclear part)
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci         :b is determined internally from b=rofi(nr)/(dexp(a*nr-a)-1)
Ci   rofi  :radial mesh points
Ci   i0    :*If i0<4, integrate from origin to rofi(nr)
Ci         :*Else, integrate from rofi(i0) to rofi(nr).
Ci         : NB: in this case, g(i0-2..i0) must be input!
Ci   nr    :number of radial mesh points
Ci   g0    :rhs of inhomogeneous equation, defined as
Ci         :r*(true radial function)
Ci         :OLD: r*(true radial function)/dsqrt(a*(rofi(ir)+b))
Co Outputs
Co   g     :wave function defined as  r*(real radial function)
Co         :NB: points g(1..i0-1) must already exist.
Co   nn    :number of nodes encountered
Co   val   :value of wave function at rofi(nr)
Co   slo   :slope of wave function at rofi(nr)
Co   q     :<g g>
Cb Bugs
Cb   slo calculated by from integration from orgin.  Reliability unclear
Cl Local variables
Cr Remarks
Cr   For historical reasons, this program internally integrates not
Cr   g, but g/dsqrt(a*(r+b)).  Internally g and g0 are scaled;
Cr   on exit g and g0 are restored to r*(true radial function)
Cu Updates
Cu   29 Jun 04 Adapted from old rsqnri (see old molecules program)
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer l,nn,i0,nr
      double precision e,z,val,slo,q,a
      double precision g(nr),g0(nr),v(nr),rofi(nr)
C ... Local parameters
      integer icc,ir,lp1,ncc,i1
      double precision a2b4,b2,b3,b4,blp1,blp2,fllp1,r,r2,r3,s,sum,t,
     .  vme,wgt,x,x2,x3,y,zz,b,cc(20),yold

      ncc = 15
      a2b4 = a*a/4d0
      b = rofi(nr)/(dexp(a*nr-a)-1d0)
      if (i0 .ne. 1 .and. i0 .lt. 4)
     .  call rxi('rsqnri: bad starting point',i0)

C --- For existing points, convert g to internal g/srdrdi ---
      do  ir = 1, i0-1
        r = rofi(ir)
        g(ir) = g(ir)/dsqrt(a*(r+b))
      enddo
      do  ir = 1, nr
        r = rofi(ir)
        g0(ir) = g0(ir)/dsqrt(a*(r+b))
      enddo

C --- Power series near zero ---
      lp1 = l+1
      r2 = rofi(2)
      r3 = rofi(3)
      x2 = g0(2)*dsqrt(a*(r2+b))/r2**lp1
      x3 = g0(3)*dsqrt(a*(r3+b))/r3**lp1
      blp1 = -(r2*x3-r3*x2)/(r2-r3)
C     blp20 = -(x3-x2)/(r3-r2)
      blp2 = -z*blp1/lp1
      fllp1 = l*lp1
      zz = 2d0*z
      do  icc = 1, lp1
        cc(icc) = 0d0
      enddo
      vme = v(1)-e
      cc(lp1) = 1d0
      cc(l+2) = -z*cc(lp1)/lp1
      cc(l+3) = (-zz*cc(l+2)+vme*cc(l+1)+blp1)/(4*l+6)
      cc(l+4) = (-zz*cc(l+3)+vme*cc(l+2)+blp2)/(6*l+12)
      do  icc = l+5,ncc
        cc(icc) = (-zz*cc(icc-1)+vme*cc(icc-2))/(icc*(icc-1)-fllp1)
      enddo
      q = 0d0
      nn = 0
      slo = dsqrt(a*b)*2*a*b*cc(2)
      if (i0 .lt. 4) then
        g(1) = 0d0
        do  ir = 2, 3
          r = rofi(ir)
          sum = 0d0
          do  icc = 1, ncc
            sum = r*sum + cc(ncc-icc+1)
          enddo
          g(ir) = r*sum/dsqrt(a*(r+b))
        enddo
        i1 = 4
      else
        i1 = i0
      endif

C --- Integrate to i1-1; setup b2,b3 for points i0..nr ---
      y  = 0
      b2 = 0
      b3 = 0
      do  ir = 2, i1-1
        r = rofi(ir)
        s = (a*(r+b))**2
        b4 = -g0(ir)*s
        wgt = 2*(mod(ir+1,2)+1)
        t = ((fllp1/r-zz)/r+v(ir)-e)*s + a2b4
        slo = slo + wgt*(g(ir)*t + b4)
        q = q + wgt*g(ir)*g(ir)*s
        yold = y
        y = g(ir)*(1d0-t/12)
        b2 = b3
        b3 = b4
      enddo
      x = y-yold

C --- Integrate g over points i1..nr; add to q ---
      do  ir = i1, nr
        r = rofi(ir)
        s = (a*(r+b))**2
        b4 = -g0(ir)*s
        wgt = 2*(mod(ir+1,2)+1)
        x = x + t*g(ir-1) + (b4+10*b3+b2)/12
        y = y+x
        t = ((fllp1/r-zz)/r+v(ir)-e)*s + a2b4
        g(ir) = y/(1d0-t/12)
        if (g(ir)*g(ir-1) .lt. 0d0) nn = nn+1
        slo = slo + wgt*(g(ir)*t + b4)
        q = q + wgt*g(ir)*g(ir)*s
        b2 = b3
        b3 = b4
      enddo

C --- Cleanup ---
      r = rofi(nr)
      q = (q-g(nr)*g(nr)*s)/3d0
      slo = (slo-(g(nr)*t+b4))/3d0 + dsqrt(a*b)*cc(1)
      slo = (.5d0*a*g(nr)+slo)/dsqrt(a*(r+b))/r
      val = g(nr)*dsqrt(a*(r+b))/r
      do  ir = 1, nr
        r = rofi(ir)
        g0(ir) = g0(ir)*dsqrt(a*(r+b))
        g(ir)  = g(ir)*dsqrt(a*(r+b))
      enddo
      end
      subroutine rsqnrs(e,l,v,z,rofi,i0,nr,g)
C- Small component of phi from the large component
C ----------------------------------------------------------------------
Ci Inputs
Ci   e     :Energy
Ci   z     :nuclear charge
Ci   l     :angular momentum
Ci   v     :spherical potential (true potential excluding nuclear part)
Ci   rofi  :radial mesh points
Ci   i0    :*If i0<4, integrate from origin to rofi(nr)
Ci         :*Else, integrate from rofi(i0) to rofi(nr).
Ci         : NB: in this case, g(i0-2..i0) must be input!
Ci   nr    :number of radial mesh points
Cio Inputs/Outputs
Co   g     :wave function defined as  r*(real radial function)
Co         :g(1..nr,1) is large component (input)
Co         :g(1..nr,2) is small component.
Co         :This routine generates g(i0..nr,2)
Cr Remarks
Cr   Let P = r*(large component of true radial function)
Cr       Q = r*(small component of true radial function)
Cr   (this is Skriver's notation)
Cr   Then let k = -l-1 (i.e., kappa)
Cr     Qk = (Pk' + (k/r)Pk) / D  where D=c(1-(e-V)/c^2)
Cu Updates
Cu   13 Jul 04 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer l,i0,nr
      double precision e,z
      double precision g(nr,2),v(nr),rofi(nr)
C ... Local parameters
      integer ir,nrx,nptirp
      parameter (nrx=1501,nptirp=10)
      double precision r,dg(nrx),kap,vi,tmc,c
C ... Speed of light, or infinity in nonrelativistic case
      common /cc/ c

C ... dg <- g'(1..nr,1) = P'
      call poldvm(rofi,g,nr,nptirp,.false.,1d-8,ir,dg)
      kap = -l-1

C ... g(*,2) = Q = (Pk' + (k/r)Pk) / D
      do  ir = i0, nr
        r = rofi(ir)
        vi = v(ir) - 2d0*z/r
        tmc = c - (vi-e)/c
        g(ir,2) = (dg(ir) + (1+kap)/r*g(ir,1))/tmc
      enddo

      end

