      subroutine atqval(lmxa,pnu,pnz,z,kcor,lcor,qcor,qc,qv,qsc)
C- Return valence and core charge for one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   lmxa  :augmentation l-cutoff
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   pnz   :boundary conditions for semicore state
Ci          pnz = 0 => no semicore state
Ci   z     :nuclear charge
Ci   kcor  :(partial core occupation) p.q.n for occupation
Ci   lcor  :(partial core occupation) l quantum for occupation
Ci   qcor  :(partial core occupation) core charge
Co Outputs
Ci   qc    :sphere core charge
Ci   qv    :nuclear charge - sphere core charge
Ci   qsc   :sphere semicore charge
Cr Remarks
Cu Updates
Cu   16 Sep 01 Added calculation of qsc.  New argument list.
Cu   30 May 00 adapted from nfp getqval
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer kcor,lcor,lmxa
      double precision qc,qcor(1),qv,qsc,z,pnu(0:lmxa),pnz(0:lmxa)
C ... Local parameters
      integer l,konf,konfig
      double precision deg

      qc = 0d0
      qsc = 0d0
      do  l = 0, lmxa
        deg = 2*(2*l+1)
        if (pnz(l) .ne. 0 .and. mod(int(pnz(l)),10) .lt. int(pnu(l)))
     .    qsc = qsc + deg
        konfig = pnu(l)
        do  konf = l+1, konfig-1
          if (konf.eq.kcor .and. l.eq.lcor) deg = deg + qcor(1)
          qc = qc + deg
        enddo
      enddo
      qv = z-qc
      end

