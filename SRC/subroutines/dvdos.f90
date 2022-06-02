subroutine dvdos(vnow,nosnow,dosnow,vhold,ztarg,dvcap,dv)
  !- Estimate shift in potential shift to meet target number-of-states
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   vnow  :current value of potential shift
  !i   nosnow:current value of charge
  !i   dosnow:current density of states, d nosnow / dv
  !i         :(used only for first iteration to estimate rfalsi step size)
  !i   vhold :vector of 12 numbers, maintained internally by rfalsi
  !i         :Initial call: vhold should be zero.
  !i   ztarg :desired charge
  !i   dvcap :maximum change in potential shift for any step
  !o Outputs
  !o   vnow  :updated value of estimated potential shift
  !o   dv    :change in vnow this step
  !r Remarks
  !r   Routine uses regula falsi to iteratively find target nosnow.
  !u Updates
  !u   22 Sep 01 Adapted from pvpgzq
  !-----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  double precision :: vnow,nosnow,dosnow,vhold(12),ztarg, &
       dvcap,dv
  ! ... Local parameters
  double precision :: dznow,dxmx
  integer :: ir

  ! ... First order estimate dv = (ztarg-zhave)/slope
  dznow = nosnow-ztarg
  dv = dznow/max(dosnow,1d-5)
  !     Hang onto vhold(12) because rfalsi destroys it
  ir = nint(vhold(12))
  ! ... first iteration
  if (ir == 0) then
     dxmx = dznow/max(dosnow,1d-5)
     if (abs(dxmx) > dvcap) dxmx = sign(dvcap,dxmx)
     ir = 0
  else
     dxmx = dvcap
  endif
  !     vhold = vnow
  call pshpr(0)
  !     call rfalsi(vnow,dznow,1d-6,0d0,1d-6,dxmx,10,vhold(1),ir)
  !     call rfalsi(vnow,dznow,1d-7,0d0,1d-7,dxmx,10,vhold(1),ir)
  call rfalsi(vnow,dznow,5d-8,0d0,5d-8,dxmx,10,vhold(1),ir)
  !     print *, ir
  call poppr
  vhold(12) = ir
  dv = vnow - vhold(1)

end subroutine dvdos

