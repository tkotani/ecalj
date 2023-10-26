subroutine scglp1(mlm, kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2)  !- Makes Clebsch-Gordan coefficients coupling to l+1 for one mlm. See ropylg and m_smhankel
  use m_ll,only:ll
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mlm
  !o Outputs
  !o   kz    :ilm of z component for ll(mlm)+1; see Remarks
  !o   cz    :coefficient for z component for ll(mlm)+1
  !o   kx1   :ilm of 1st x component for ll(mlm)+1; see Remarks
  !o   kx2   :ilm of 2nd x component for ll(mlm)+1; see Remarks
  !o   cx1   :coefficient for 1st x component for ll(mlm)+1
  !o   cx2   :coefficient for 2nd x component for ll(mlm)+1
  !o   ky1   :ilm of 1st y component for ll(mlm)+1; see Remarks
  !o   ky2   :ilm of 2nd y component for ll(mlm)+1; see Remarks
  !o   cy1   :coefficient for 1st y component for ll(mlm)+1
  !o   cy2   :coefficient for 2nd y component for ll(mlm)+1
  !r Remarks
  !r   Gradients can be generated from smoothed Hankels as follows: Coefficients must still be multiplied by sqrt(3/4pi)
  implicit none
  integer :: mlm,kz,kx1,kx2,ky1,ky2,l,lav,m,mm,isg,kav,ma,mb
  real(8) :: cz,cx1,cx2,cy1,cy2,bot,top,tap,cofa,cofb
  l = ll(mlm)
  lav = l*l+l+1
  m = mlm-lav
  mm = iabs(m)
  isg = 1
  if (m < 0) isg = -1
  kav = (l+1)*(l+1)+(l+1)+1
  bot = (2*l+1)*(2*l+3)
  !     z coefficient
  kz = kav+m
  cz = dsqrt((l+mm+1)*(l-mm+1)/bot)
  !     x,y coefficients
  top = (l+1+mm)*(l+2+mm)/2d0
  tap = (l+1-mm)*(l+2-mm)/2d0
  if (mm /= 0) top = 0.5d0*top
  if (mm /= 1) tap = 0.5d0*tap
  cofa = dsqrt(top/bot)
  cofb = dsqrt(tap/bot)
  ma = isg*(mm+1)
  mb = isg*(mm-1)
  kx1 = kav+ma
  cx1 = cofa
  kx2 = kx1
  cx2 = 0d0
  if (m /= -1 .AND. m /= 0) then
     kx2 = kav+mb
     cx2 = -cofb
  endif
  ky1 = kav-ma
  cy1 = isg*cofa
  ky2 = ky1
  cy2 = 0d0
  if (m /= 0 .AND. m /= 1) then
     ky2 = kav-mb
     cy2 = isg*cofb
  endif
end subroutine scglp1
