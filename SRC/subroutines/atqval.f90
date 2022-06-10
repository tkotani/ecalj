subroutine atqval(lmxa,pnu,pnz,z,kcor,lcor,qcor,qc,qv,qsc)
  !- Return valence and core charge for one site
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   lmxa  :augmentation l-cutoff
  !i   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
  !i          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
  !i   pnz   :boundary conditions for semicore state
  !i          pnz = 0 => no semicore state
  !i   z     :nuclear charge
  !i   kcor  :(partial core occupation) p.q.n for occupation
  !i   lcor  :(partial core occupation) l quantum for occupation
  !i   qcor  :(partial core occupation) core charge
  !o Outputs
  !i   qc    :sphere core charge
  !i   qv    :nuclear charge - sphere core charge
  !i   qsc   :sphere semicore charge
  !r Remarks
  !u Updates
  !u   16 Sep 01 Added calculation of qsc.  New argument list.
  !u   30 May 00 adapted from nfp getqval
  ! ----------------------------------------------------------------------
  implicit none
  integer :: kcor,lcor,lmxa
  double precision :: qc,qcor(1),qv,qsc,z,pnu(0:lmxa),pnz(0:lmxa)
  integer :: l,konf,konfig
  double precision :: deg
  qc = 0d0
  qsc = 0d0
  do  l = 0, lmxa
     deg = 2*(2*l+1)
     if(pnz(l) /= 0 .AND. mod(int(pnz(l)),10) < int(pnu(l))) qsc = qsc + deg
     konfig = pnu(l)
     do  konf = l+1, konfig-1
        if (konf == kcor .AND. l == lcor) deg = deg + qcor(1)
        qc = qc + deg
     enddo
  enddo
  qv = z-qc
end subroutine atqval

