subroutine gtbsl2(l1,lmxh,eh,rsmh,l2)
  !- Returns the highest l with rsm,e common to those of a given l
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   l1    :current l
  !i   lmxh :basis l-cutoff
  !i   eh    :energy of smoothed Hankel
  !i   rsmh  :smoothing radius of smoothed hankel
  !o Outputs
  !o   l2    :large l for which eh and rsmh l1..l2 are in common
  !r Remarks
  !r   Routine used group functions, strux in blocks.
  !u Updates
  ! ----------------------------------------------------------------------
  !     implicit none
  integer :: l1,l2,lmxh
  double precision :: rsmh(0:lmxh),eh(0:lmxh)
  double precision :: e,rsm

  e = eh(l1)
  rsm = rsmh(l1)
  l2 = l1
10 if (l2 >= lmxh) return
  if (rsmh(l2+1) /= rsm .OR. eh(l2+1) /= e) return
  l2 = l2+1
  goto 10

end subroutine gtbsl2

