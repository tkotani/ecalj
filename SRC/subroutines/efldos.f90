subroutine efldos(qval,nsp,emin,emax,ndos,dos,eferm,eband)
  use m_lgunit,only:stdo,stdl
  !- Linear interpolation for Fermi and band energy from dos and i-dos
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   qval  :system charge
  !i   emin  :first energy point on DOS mesh
  !i   emax  :last energy point on DOS mesh
  !i   ndos  :number of energy mesh points
  !i   dos   :dos(i,1) is the dos, dos(i,2) is the integrated dos
  !o Outputs
  !i   eferm :Fermi level
  !i   eband :sum energy bands
  !r Remarks
  !u Updates
  !u    7 Jul 00 spin polarized
  !u   30 May 00 Adapted from nfp efdos.f
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: ndos,nsp
  double precision :: eband,eferm,emax,emin,qval,dos(ndos,2,nsp)
  ! ... Local parameters
  integer :: i,i1,i2,ibot,ie,ipl,iprint,k
  double precision :: de,dferm,e1,ebot,sum,x,dossp
  dossp(k) = (dos(k,2,1)+dos(k,2,nsp))/dble(3-nsp)

  !      stdo = lgunit(1)
  !      stdl = lgunit(2)
  ipl  = iprint()

  ! ... Find interval containing fermi energy
  de = (emax-emin)/(ndos-1)
  ibot = 0
  do  ie = 1, ndos
     if (dabs(dossp(ie)) > 1d-6 .AND. ibot == 0) ibot = ie
     i2 = ie
     if (dossp(ie) > qval) goto 90
  enddo
  call rx('efldos: fermi energy lies above emax')
  eferm = 0
  eband = 0
  return

90 continue
  i1 = i2-1
  if (i1 < 0) then
     write(stdo,*)'efldos (warning): fermi energy lies below emin'
     eferm = 0
     eband = 0
     return
  endif

  ! ... Linear interpolation of i-dos to find fermi energy
  x = (qval-dossp(i1))/(dossp(i2)-dossp(i1))
  ebot = emin + (ibot-1)*de
  e1 = emin + (i1-1)*de
  eferm = e1 + de*x

  ! ... Calculate eigenvalue sum as ef*N(ef)-integral N(e)
  if (dabs(dossp(1)) >= 1d-6) call rx1( &
       'efldos: DOS nonzero at first energy point: DOS=%g',dossp(1))

  sum = 0d0
  do  i = 1, i1-1
     sum = sum + de*0.5d0*(dossp(i)+dossp(i+1))
  enddo
  sum = sum+0.5d0*(eferm-e1)*(dossp(i1)+qval)
  eband = eferm*qval - sum
  dferm = dos(i1,1,1) + x*(dos(i2,1,1)-dos(i1,1,1))
  if (nsp == 2) dferm = &
       dferm + dos(i1,1,2) + x*(dos(i2,1,2)-dos(i1,1,2))

  !      if (iprint() .ge. 20) then
  !        write(stdo,"(' efldos:  Ef = ',f16.8,' sumev = ',f16.8,'  ebot =',f16.8,
  !     .  ' DOS(Ef) = ',f16.8)")eferm,eband,ebot,dferm
  !      endif
  write(stdo,729) eferm,ebot,dferm,eband
729 format(/' efldos: fermi energy=',f10.5, &
       /' dos starts at energy  ',f11.5 &
       /' dos at fermi energy=  ',f11.5,' sumev=eband=',f11.5)

  if (ipl > 0) write(stdl,771) eferm,ebot,dferm,eband
771 format('nf dos ef',f9.5,'  ebot',f7.3,'  d(ef)',f9.2, &
       '  ebn',f11.6)

end subroutine efldos

