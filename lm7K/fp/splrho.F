      subroutine splrho(mode,nsp,nr,nlml,rho1,rho2,rhoc)
C- Overwrite spin pol local rho+,rho- with rho,rho+ - rho-, or reverse
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :0 input (rho+,rho-) -> (rho+ + rho-, rho+ - rho-)
Ci         :1 input (rho+ + rho-, rho+ - rho-) -> (rho+,rho-)
Ci         :10s digit
Ci         :1 suppress splitting of rho2
Ci         :2 suppress splitting of rhoc
Ci         :3 suppress both
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nr    :number of radial mesh points
Ci   nlml  :L-cutoff
Ci   rho1  :local true density, tabulated on a radial mesh
Ci   rho2  :local smoothed density, tabulated on a radial mesh
Ci   rhoc  :core density
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
      integer mode,nsp,nr,nlml
      double precision rho1(nr,nlml,nsp),rho2(nr,nlml,nsp),rhoc(nr,nsp)
      double precision fac

      if (nsp .eq. 1) return
      fac = 1
      if (mod(mode,10) .ne. 0) fac = .5d0

      call dsumdf(nr*nlml,fac,rho1,0,1,rho1(1,1,2),0,1)
      if (mod(mod(mode/10,10),2) .eq. 0)
     .  call dsumdf(nr*nlml,fac,rho2,0,1,rho2(1,1,2),0,1)
      if (mod(mod(mode/10,10)/2,2) .eq. 0)
     .  call dsumdf(nr,fac,rhoc,0,1,rhoc(1,2),0,1)
      end

      subroutine lcrho(nr,nsp,nlml1,nlml2,fac1,fac2,rho1,rho2)
C- First density is overwritten by linear combination of two densities
C ----------------------------------------------------------------------
Ci Inputs
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlml1 :number of L channels for first density
Ci   nlml2 :number of L channels for second density
Ci   fac1  :scales rho1
Ci   fac2  :scales rho2
Ci   rho2  :second density, dimensioned (nr,nlm2,nsp)
Cio Inputs/Outputs
Ci   rho1  :first density, dimensioned (nr,nlm1,nsp).  On output
Cio        :    rho1 <- fac1 * rho1 + fac2 * rho2
Cio        : fac1 scales rho1 in all nlm1 channels
Cio        : fac2 * rho2 is added into rho1 for min(nlm1,nlm2) channels
Cl Local variables
Cl         :
Cr Remarks
Cr   rho1 <- fac1 * rho1 + fac2 * rho2
Cu Updates
Cu   01 Jul 08 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nr,nsp,nlml1,nlml2
      double precision fac1,fac2
      double precision rho1(nr,nlml1,nsp),rho2(nr,nlml2,nsp)
C ... Local parameters
      integer isp,nlml

      if (nlml1 .le. 0) return
      nlml = min(nlml1,nlml2)
      do  isp = 1, nsp
        if (fac1 .ne. 1) then
          call dscal(nr*nlml1,fac1,rho1(1,1,isp),1)
        endif
        if (fac2 .ne. 0) then
          call daxpy(nr*nlml,fac2,rho2(1,1,isp),1,rho1(1,1,isp),1)
        endif
      enddo
      end

      subroutine swrho(mode,nr,nsp,nlml,nlml1,nlml2,rho1,rho2)
C- Swap two local densities, possibly with spin flip
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :1 swap spin up 1st density w/ spin down second
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlml  :number of L channels to swap
Ci   nlml1 :number of L channels for first density
Ci   nlml2 :number of L channels for second density
Ci   rho2  :second density, dimensioned (nr,nlm2,nsp)
Cio Inputs/Outputs
Cio  rho1  :first density, dimensioned (nr,nlm1,nsp).  On output
Cio  rho2  :second density, dimensioned (nr,nlm2,nsp).  On output
Cio        :densities are exchanged, and spins possibly swapped.
Cl Local variables
Cl         :
Cr Remarks
Cu Updates
Cu   19 Jul 08 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nr,nsp,nlml,nlml1,nlml2
      double precision rho1(nr,nlml1,nsp),rho2(nr,nlml2,nsp)
C ... Local parameters
      integer isp,jsp,nlmll

      nlmll = min(nlml,min(nlml1,nlml2))
      if (nlmll .le. 0) return

      do  isp = 1, nsp
        jsp = isp
        if (mode .eq. 1 .and. nsp .eq. 2) jsp = 3-isp
        call dswap(nr*nlmll,rho2(1,1,isp),1,rho1(1,1,jsp),1)
      enddo
      end

