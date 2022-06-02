subroutine config(pnu,lmax,z,konfig,lmaxc)
  !- Returns principal quantum numbers from pnu, estimating those unknown
  ! ----------------------------------------------------------------
  !i Inputs
  !i   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
  !i          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
  !i   lmax  :maximum l for which pnu is supplied
  !i   z     :nuclear charge
  !o Outputs
  !o   konfig:estimated principal quant. no. of valence state for l>lmax
  !o         :Core orbitals are specified by:
  !o         :  1, 2, ..., konf(0)-1 for s
  !o         :  2, 3, ..., konf(1)-1 for p
  !o         :  3, 4, ..., konf(2)-1 for d, and so on.
  !o   lmaxc :largest l for which there is a core state
  !r Remarks
  !u Updates
  !u   08 Feb 01 generate config for l=0..8
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  integer :: lmax,konfig(0:8),lmaxc
  double precision :: pnu(0:*),z
  ! Local parameters
  integer :: l
  integer :: LI,NA,K,RB,CS,FR,CU,AG,AU,HF
  parameter (LI=3, &
       NA=11, K=19, RB=37, CS=55, FR=87, &
       CU=29, AG=47, AU=79, HF=72)

  ! --- Calculate lmaxc ---
  lmaxc = lmax
  if (z >= LI) lmaxc = max(lmax,0)
  if (z >= NA) lmaxc = max(lmax,1)
  if (z >= CU) lmaxc = max(lmax,2)
  if (z >= HF) lmaxc = max(lmax,3)

  ! --- Estimate konfig ---
  do  10  l = 0, 8
     konfig(l) = l+1
10 enddo
  if (z >= LI) konfig(0) = 2
  if (z >= NA) konfig(0) = 3
  if (z >=  K) konfig(0) = 4
  if (z >= RB) konfig(0) = 5
  if (z >= CS) konfig(0) = 6
  if (z >= FR) konfig(0) = 7
  konfig(1) = max(konfig(0),2)

  if (z >= CU) konfig(2) = 4
  if (z >= AG) konfig(2) = 5
  if (z >= AU) konfig(2) = 6

  if (z >= HF) konfig(3) = 5

  ! --- Override konfig with given pnu ---
  do  20  l = 0, lmax
     konfig(l) = pnu(l)
20 enddo
end subroutine config

