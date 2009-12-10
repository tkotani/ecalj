      subroutine gprdsr(mode,g1,g2,nrx,n1,nr,z,e1,e2,l,v,rofi,rwgt,sum)
C- Integrate inner product of two wave equations
C ----------------------------------------------------------------------
Ci Inputs
Ci  mode   :0 use both large and small components of radial w.f.
Ci         :1 use large component of radial w.f. only
Ci         :Add 10 to compute all three inner products
Ci   g1    :First  radial wave function
Ci   g2	   :Second radial wave function
Ci   nrx   :leading dimension of g1,g2
Ci   n1    :first radial mesh point for integration
Ci   nr    :last  radial mesh point for integration
Ci   z     :nuclear charge
Ci         :(not used if if large component integration only)
Ci   e1    :energy of first wave function
Ci         :(not used if if large component integration only)
Ci   e2    :energy of second wave function
Ci         :(not used if if large component integration only)
Ci   l     :l quantum number of g1,g2
Ci         :(not used if if large component integration only)
Ci   v     :spherical potential
Ci         :(not used if if large component integration only)
Ci   rofi  :radial mesh points
Ci   rwgt  :radial mesh weights
Co Outputs:
Co   sum :inner product, depending on mode
Co       :10s digit mode=0:
Co       :sum(1) = <g1 g2>; no other quantities are returned
Co       :10s digit mode=1
Co       :sum(1) = <g1 g1>; sum(2) = <g2 g2>; sum(3) = <g1 g2>
Cr Remarks:
Cu Updates
Cu   29 Jun 04 First created
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,n1,nr,nrx,l
      double precision z,e1,e2,sum(3),g1(nrx,2),g2(nrx,2),v(nr),
     .  rofi(nr),rwgt(nr)
C ... Local parameters
      integer ir
      double precision fllp1,c,r
      double precision tmcr,gf11,gf22,gf12,s1,s2,s3
C     Speed of light, or infinity in nonrelativistic case
      common /cc/ c

      fllp1 = l*(l+1)
      s1 = 0
      s2 = 0
      s3 = 0

C ... Include small component, with averaged tmcr
      if (mod(mode,10) .eq. 0) then
        do  ir = n1, nr
          r = rofi(ir)
          tmcr = r*c - (r*v(ir) - 2d0*z - r*e1)/c
          gf11 = 1d0 + fllp1/tmcr**2
          tmcr = r*c - (r*v(ir) - 2d0*z - r*e2)/c
          gf22 = 1d0 + fllp1/tmcr**2
          gf12 = (gf11 + gf22)/2
          s1 = s1 + rwgt(ir)*(gf11*g1(ir,1)*g1(ir,1)+g1(ir,2)*g1(ir,2))
          s2 = s2 + rwgt(ir)*(gf22*g2(ir,1)*g2(ir,1)+g2(ir,2)*g2(ir,2))
          s3 = s3 + rwgt(ir)*(gf12*g1(ir,1)*g2(ir,1)+g1(ir,2)*g2(ir,2))
        enddo

C ... Large component only
      else
        do  ir = n1, nr
C         r = rofi(ir)
          s1 = s1 + rwgt(ir)*g1(ir,1)*g1(ir,1)
          s2 = s2 + rwgt(ir)*g2(ir,1)*g2(ir,1)
          s3 = s3 + rwgt(ir)*g1(ir,1)*g2(ir,1)
        enddo
      endif

      if (mod(mode/10,10) .eq. 0) then
        sum(1) = s3
      else
        sum(1) = s1
        sum(2) = s2
        sum(3) = s3
      endif

      end

