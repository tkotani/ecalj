      subroutine dvdos(vnow,nosnow,dosnow,vhold,ztarg,dvcap,dv)
C- Estimate shift in potential shift to meet target number-of-states
C ----------------------------------------------------------------------
Ci Inputs
Ci   vnow  :current value of potential shift
Ci   nosnow:current value of charge
Ci   dosnow:current density of states, d nosnow / dv
Ci         :(used only for first iteration to estimate rfalsi step size)
Ci   vhold :vector of 12 numbers, maintained internally by rfalsi
Ci         :Initial call: vhold should be zero.
Ci   ztarg :desired charge
Ci   dvcap :maximum change in potential shift for any step
Co Outputs
Co   vnow  :updated value of estimated potential shift
Co   dv    :change in vnow this step
Cr Remarks
Cr   Routine uses regula falsi to iteratively find target nosnow.
Cu Updates
Cu   22 Sep 01 Adapted from pvpgzq
C-----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      double precision vnow,nosnow,dosnow,vhold(12),ztarg,
     .dvcap,dv
C ... Local parameters
      double precision dznow,dxmx
      integer ir

C ... First order estimate dv = (ztarg-zhave)/slope
      dznow = nosnow-ztarg
      dv = dznow/max(dosnow,1d-5)
C     Hang onto vhold(12) because rfalsi destroys it
      ir = nint(vhold(12))
C ... first iteration
      if (ir .eq. 0) then
        dxmx = dznow/max(dosnow,1d-5)
        if (abs(dxmx) .gt. dvcap) dxmx = sign(dvcap,dxmx)
        ir = 0
      else
        dxmx = dvcap
      endif
c     vhold = vnow
      call pshpr(0)
C     call rfalsi(vnow,dznow,1d-6,0d0,1d-6,dxmx,10,vhold(1),ir)
C     call rfalsi(vnow,dznow,1d-7,0d0,1d-7,dxmx,10,vhold(1),ir)
      call rfalsi(vnow,dznow,5d-8,0d0,5d-8,dxmx,10,vhold(1),ir)
C     print *, ir
      call poppr
      vhold(12) = ir
      dv = vnow - vhold(1)

      end

