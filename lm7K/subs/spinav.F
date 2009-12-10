      subroutine spinav(mode,nclass,nl,nsp,pnu,qnu)
C- Averages up+down spin moments + pp's for all classes
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 average spins
Ci         :1 do not average, but exchange spins
Ci   nclass:number of inequivalent classes
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   qnu   :energy-weighted moments of the sphere charges
Co Outputs :moments are spin-averaged
Ci   pnu   :spin-averaged (mode=0) or spin-flipped (mode=1)
Ci   qnu   :spin-averaged (mode=0) or spin-flipped (mode=1)
Co   nsp   :set to 1 on output (mode=0)
Cl Local variables
Cl         :
Cr Remarks
Cr   
Cu Updates
Cu   10 Jan 06 Added mode
C ----------------------------------------------------------------------
C     implicit none
      integer mode,nclass,nl,nsp,ic
      double precision pnu(nl,nsp,nclass),qnu(3,nl,nsp,nclass)

      if (nsp .eq. 1) return

      if (mode .gt. 1) then
        call rx('spinav: bad mode')
      elseif (mode .eq. 1) then
        do  ic = 1, nclass
          call dswap(nl,pnu(1,2,ic),1,pnu(1,1,ic),1)
          call dswap(3*nl,qnu(1,1,2,ic),1,qnu(1,1,1,ic),1)
        enddo
        return
      endif

      do  ic = 1, nclass
        call daxpy(nl,1d0,pnu(1,2,ic),1,pnu(1,1,ic),1)
        call dscal(nl,.5d0,pnu(1,1,ic),1)
        call daxpy(3*nl,1d0,qnu(1,1,2,ic),1,qnu(1,1,1,ic),1)
      enddo

      do  ic = 2, nclass
        call dcopy(nl,pnu(1,1,ic),1,pnu(1,ic,1),1)
        call dcopy(3*nl,qnu(1,1,1,ic),1,qnu(1,1,ic,1),1)
      enddo
      nsp = 1
      end

