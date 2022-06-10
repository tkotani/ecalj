subroutine fermi2( qval,dos,ndos,emin,emax,  eferm,e1,e2,dosef)
  !- Makes fermi energy from integrated density
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   qval, number of electrons to fermi level; dos(i) integrated
  !i   density at bin i; ndos, number of bins + 1; emin, emax, energy
  !i   window.
  !o Outputs
  !o   Eferm, Fermi energy; e1, e2, confidence limits on Fermi energy
  !o   i.e., Fermi energy lies between e1 and e2.
  !o   dosef:  density of states at fermi level
  !r Remarks
  !r   emin and e1 (and emax and e2) may point to the same address.
  ! ----------------------------------------------------------------------
  implicit none
  ! Passed parameters
  integer:: ndos
  real(8):: qval,dos,emin,emax,eferm,e1,e2,dosef
  ! Local parameters
  integer:: i1,ie
  real(8):: de,q,q1,q2,d1mach
  ! External procedures
  !external:: d1mach
  DIMENSION DOS(NDOS)
  ! ccccccccccccccccc
  !      DE = (EMAX-EMIN)/(NDOS-1)
  !      do i1=1,ndos
  !        write(6,"(' e dostot =',i4,d13.6,d13.6)")i1, emin + de*(i1 - 1),dos(i1)
  !      enddo
  !      write(6,"(' qval =',d13.6)") qval
  ! ccccccccccccccccccc
  if (dos(1) > qval) print *, 'FERMI: EMIN,EMAX=', emin,emax
  if (dos(1) > qval) call rx( 'FERMI: Fermi energy lies below EMIN')
  if (dos(ndos) < qval) print *, 'FERMI: EMIN,EMAX=', emin,emax
  if (dos(ndos) < qval) then
     call rx( 'FERMI: Fermi energy lies above EMAX')
  endif
  DE = (EMAX-EMIN)/(NDOS-1)
  I1 = 1
  q = qval + d1mach(3)
  DO  1  IE = 2, NDOS
     I1 = IE
     IF ( DOS(IE) > q ) goto 2
1 enddo
2 continue
  i1 = i1 - 1
  Q1 = DOS(I1)
  Q2 = DOS(I1+1)
  ! ------------------
  e1 = emin + de*(i1 - 1)
  e2 = e1 + de
  ! rite(6,"('fermi:: i1,e1,qval,q1,q2,de',I7,5E13.5)") i1,e1,qval,q1,q2,de
  EFERM = e1 + (QVAL-Q1)/(Q2-Q1)*DE
  dosef = (q2-q1)/de
end subroutine fermi2

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
subroutine fermi_kbt( qval,dos,ndos,emin,emax, kbt, eferm_init, &
     eferm_kbt)
  implicit none
  intent(in)::   qval,dos,ndos,emin,emax, kbt, eferm_init
  intent(out)::  eferm_kbt
  !- Makes fermi energy from integrated density
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   qval, number of electrons to fermi level; dos(i) integrated
  !i   density at bin i; ndos, number of bins + 1; emin, emax, energy
  !i   window.
  !o Outputs
  !o   Eferm, Fermi energy; e1, e2, confidence limits on Fermi energy
  !o   i.e., Fermi energy lies between e1 and e2.
  !o   dosef_kbt:  density of states at fermi level
  !r Remarks
  !r   emin and e1 (and emax and e2) may point to the same address.
  ! ----------------------------------------------------------------------
  integer :: ndos
  real(8):: qval,dos,emin,emax,ddos,eferm_old,efermi_init
  real(8):: eferm_init,eferm_kbt,e0,e1,e2,dosef_kbt
  real(8):: kbt, ec
  integer:: i1,ie,ii,nii=100, ii_conv=-999, i
  real(8):: de,q,q1,q2,d1mach
  real(8):: idos, idos_old, idos_l, rydberg
  real(8):: eferm_i(4) !!Efermi(i+1),Efermi(i),Efermi(i-1),Efermi(i-2)
  ! External procedures
  external d1mach
  logical:: debug = .false.
  logical:: lessqval, moreqval
  DIMENSION dos(ndos)
  if (debug) then
     write(6,"('ndos',I8)") ndos
     write(6,"('qval,dos(1),dos(ndos)',3E13.5)") qval,dos(1),dos(ndos)
     write(6,"('emin,emax',2E13.5)") emin,emax
     write(6,"('kbt,efermi_init',2E13.5)") kbt,efermi_init
  endif
  !$$$      do ii = 2,ndos
  !$$$         write(6,*) "aaaaaaaa:",ii,dos(ii)-dos(ii-1)
  !$$$      enddo
  !$$$      call rx("end")
  if (dos(1) > qval) print *, 'FERMI: EMIN,EMAX=', emin,emax
  if (dos(1) > qval) call rx( 'FERMI: Fermi energy lies below EMIN')
  if (dos(ndos) < qval) print *, 'FERMI: EMIN,EMAX=', emin,emax
  if (dos(ndos) < qval) then
     call rx( 'FERMI: Fermi energy lies above EMAX')
  endif
  de = (emax-emin)/(ndos-1)
  q = qval + d1mach(3)
  ! rite(6,"('fermi_kbt_plot:: ## fermi level (T=0 K)',E13.5)") eferm_init
!!! self-consistently determined Fermi energy Efermi(T)
  eferm_i=-9999; eferm_i(1)=eferm_init
  lessqval = .false.
  moreqval = .false.
  write(6,*) "fermi_i: iter, EF(i) EF(i-1) EF(i-2) EF(i-3)"
  write(6,"(' fermi_i  iter =',I5,':',4E13.5)") 0,eferm_i(1:4)
  do 200 ii = 1,nii
     if (abs(eferm_i(2)-eferm_i(1)) < 1d-8 ) cycle !! check convergence
     idos   = dos(1) !!initial
     idos_l = dos(1) !!initial; left idos is necessarry for T=0 K
     do 300 ie = 2, ndos
        ddos = dos(ie)-dos(ie-1)
        e0 = emin + de*(ie - 2)
        e1 = emin + de*(ie - 1)
        e2 = e1 + de
        ec = (e1+e2)/2d0
        idos = idos + ddos*(1d0/(exp((ec-eferm_i(1))/kbt)+1))
        idos_l = idos_l + ddos*(1d0/(exp((e0-eferm_i(1))/kbt)+1))
        if (debug) then
           write(6,"('weight of fermi distribution',2E13.5)") &
                1d0/(exp((ec-eferm_old)/kbt)+1)!,1d0/(exp((e0-eferm_old)/kbt)+1)
           write(6,"('e1,e2,ec,eferm_old',4E13.5)") e1,e2,ec,eferm_old
           write(6,"('ec,dos(i1),idos',I6,5E13.5)") ec,dos(i1),idos,idos_l
        endif
300  enddo
!!! shift Efermi(i-1) <-- Efermi(i) for Efermi(i+1)
     do i = 3,1,-1
        eferm_i(i+1) = eferm_i(i)
     enddo
!!! new Fermi energy Efermi(i+1) is determined in two cases (1) and (2)
!!! (1) Fermi energy is between Efermi(i-1) and Efermi(i)
     if (lessqval .AND. moreqval .AND. ii >= 2) then
        if (idos_l <= qval) then !! Eferm(i) <= qval
           eferm_i(1)=(eferm_i(2) + eferm_i(3))/2d0
        else                !! qval < idos_l !! Eferm(i) > qval
           eferm_i(1)=(eferm_i(2) + eferm_i(4))/2d0
        endif
!!! searching Fermi energy at finite temperature
!!! (2) Fermi energy is shifted by kbt/10
     else
        if (idos_l <= qval) then
           eferm_i(1)=eferm_i(2)+kbt/real(ii) !!! E(i+1)=E(i)+kbt
           lessqval = .true.
        else                !! qval < idos_l (Efermi moves to lower energy)
           eferm_i(1)=eferm_i(2)-kbt/real(ii) !!! E(i+1)=E(i)+kbt
           moreqval = .true.
        endif
     endif
     ! ferm_i(1)=-8888
     write(6,"(' fermi_i  iter =',I5,':',4E13.5)") ii,eferm_i(1:4)
     ii_conv=ii
     continue !! unless converged
200 enddo
  if (ii_conv == nii) then
     write(6,*) "fermi_kbt:: ----------- Note -----------"
     write(6,*) "fermi_kbt:: EFermi(T) is not converged"
     write(6,*) "fermi_kbt:: ----------------------------"
  endif
  write(6,*) "--------------------------------------------------------------------"
  write(6,"(' fermi_kbt:: Efermi(T), Efermi(0) [Ry]',2f15.7)") eferm_i(1),eferm_init
  write(6,"(' fermi_kbt:: Efermi(T), Efermi(0) [eV]',2f15.7)") eferm_i(1)*rydberg(),eferm_init*rydberg()
  write(6,*) "--------------------------------------------------------------------"
!!! return
  eferm_kbt = eferm_i(1)
end subroutine fermi_kbt
