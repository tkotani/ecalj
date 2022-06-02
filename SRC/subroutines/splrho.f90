subroutine splrho(mode,nsp,nr,nlml,rho1,rho2,rhoc)
  !- Overwrite spin pol local rho+,rho- with rho,rho+ - rho-, or reverse
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :1s digit
  !i         :0 input (rho+,rho-) -> (rho+ + rho-, rho+ - rho-)
  !i         :1 input (rho+ + rho-, rho+ - rho-) -> (rho+,rho-)
  !i         :10s digit
  !i         :1 suppress splitting of rho2
  !i         :2 suppress splitting of rhoc
  !i         :3 suppress both
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   nr    :number of radial mesh points
  !i   nlml  :L-cutoff
  !i   rho1  :local true density, tabulated on a radial mesh
  !i   rho2  :local smoothed density, tabulated on a radial mesh
  !i   rhoc  :core density
  !r Remarks
  !u Updates
  ! ----------------------------------------------------------------------
  !     implicit none
  integer :: mode,nsp,nr,nlml
  double precision :: rho1(nr,nlml,nsp),rho2(nr,nlml,nsp),rhoc(nr,nsp)
  double precision :: fac
  if (nsp == 1) return
  fac = 1
  if (mod(mode,10) /= 0) fac = .5d0
  call dsumdf(nr*nlml,fac,rho1,0,1,rho1(1,1,2),0,1)
  if (mod(mod(mode/10,10),2)  == 0) call dsumdf(nr*nlml,fac,rho2,0,1,rho2(1,1,2),0,1)
  if (mod(mod(mode/10,10)/2,2) == 0) call dsumdf(nr,fac,rhoc,0,1,rhoc(1,2),0,1)
end subroutine splrho

subroutine lcrho(nr,nsp,nlml1,nlml2,fac1,fac2,rho1,rho2)
  !- First density is overwritten by linear combination of two densities
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nr    :number of radial mesh points
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   nlml1 :number of L channels for first density
  !i   nlml2 :number of L channels for second density
  !i   fac1  :scales rho1
  !i   fac2  :scales rho2
  !i   rho2  :second density, dimensioned (nr,nlm2,nsp)
  ! o Inputs/Outputs
  !i   rho1  :first density, dimensioned (nr,nlm1,nsp).  On output
  ! o        :    rho1 <- fac1 * rho1 + fac2 * rho2
  ! o        : fac1 scales rho1 in all nlm1 channels
  ! o        : fac2 * rho2 is added into rho1 for min(nlm1,nlm2) channels
  !l Local variables
  !l         :
  !r Remarks
  !r   rho1 <- fac1 * rho1 + fac2 * rho2
  !u Updates
  !u   01 Jul 08 First created
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: nr,nsp,nlml1,nlml2
  double precision :: fac1,fac2
  double precision :: rho1(nr,nlml1,nsp),rho2(nr,nlml2,nsp)
  ! ... Local parameters
  integer :: isp,nlml

  if (nlml1 <= 0) return
  nlml = min(nlml1,nlml2)
  do  isp = 1, nsp
     if (fac1 /= 1) then
        call dscal(nr*nlml1,fac1,rho1(1,1,isp),1)
     endif
     if (fac2 /= 0) then
        call daxpy(nr*nlml,fac2,rho2(1,1,isp),1,rho1(1,1,isp),1)
     endif
  enddo
end subroutine lcrho

subroutine swrho(mode,nr,nsp,nlml,nlml1,nlml2,rho1,rho2)
  !- Swap two local densities, possibly with spin flip
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :1s digit
  !i         :1 swap spin up 1st density w/ spin down second
  !i   nr    :number of radial mesh points
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   nlml  :number of L channels to swap
  !i   nlml1 :number of L channels for first density
  !i   nlml2 :number of L channels for second density
  !i   rho2  :second density, dimensioned (nr,nlm2,nsp)
  ! o Inputs/Outputs
  ! o  rho1  :first density, dimensioned (nr,nlm1,nsp).  On output
  ! o  rho2  :second density, dimensioned (nr,nlm2,nsp).  On output
  ! o        :densities are exchanged, and spins possibly swapped.
  !l Local variables
  !l         :
  !r Remarks
  !u Updates
  !u   19 Jul 08 First created
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: mode,nr,nsp,nlml,nlml1,nlml2
  double precision :: rho1(nr,nlml1,nsp),rho2(nr,nlml2,nsp)
  ! ... Local parameters
  integer :: isp,jsp,nlmll

  nlmll = min(nlml,min(nlml1,nlml2))
  if (nlmll <= 0) return

  do  isp = 1, nsp
     jsp = isp
     if (mode == 1 .AND. nsp == 2) jsp = 3-isp
     call dswap(nr*nlmll,rho2(1,1,isp),1,rho1(1,1,jsp),1)
  enddo
end subroutine swrho

