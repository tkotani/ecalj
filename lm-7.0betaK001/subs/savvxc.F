      subroutine savvxc(nr,rho,rhoc,v,wt,avv)
C- Integrates density-weighted vxc to estimate avg XC-field
C ----------------------------------------------------------------------
Ci Inputs
Ci   nr    :number of radial mesh points
Ci   rho   :spherical charge density = 4 pi r**2 (true rho)
Ci   rhoc  :core charge density = 4 pi r**2 (true rhoc)
Ci   v     :spherical potential (atomsr.f)
Ci   wt    :mesh weights for radial integration
Co Outputs
Co   avv   :integral (v(i,1)-v(i,2)) * rhov / integral rhov
Cr Remarks
Cr   
Cu Updates
Cu   22 Feb 03 first created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nr
      double precision wt(nr),v(nr,2),rho(nr,2),rhoc(nr,2),avv
C ... Local parameters
      integer ir
      double precision rhov,srho,vxc

      srho = 0
      avv = 0
      do  ir = 1, nr
        rhov = (rho(ir,1)+rho(ir,2))-(rhoc(ir,1)+rhoc(ir,2))
        vxc = v(ir,1) - v(ir,2)
        srho = srho + wt(ir)*rhov
        avv = avv + wt(ir)*rhov*vxc
      enddo
      if (srho .eq. 0) then
        avv = 0
        return
      endif
      avv = avv / srho

      end

