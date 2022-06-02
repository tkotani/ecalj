subroutine gtpcor(sspec,is,kcore,lcore,qcore)
  use m_struc_def  !Cgetarg
  !- Unpacks parameters related to partial core occpation
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   sspec :struct containing species-specific information
  !i   is    :species for which to unpack kcore,lcore,qcore
  !o Outputs
  !o   kcore  :p.q.n for occupation
  !o   lcore  :l quantum for occupation
  !o   qcore  :core charge and magnetic moment
  !r Remarks
  !u Updates
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: is,kcore,lcore,i_copy_size
  real(8):: qcore(2)
  type(s_spec)::sspec(*)
  ! ... Local parameters
  character(8) :: ch
  kcore = 0
  lcore = -1
  qcore(1) = 0
  qcore(2) = 0
  !      do i_spacks=is,is
  !        call spacks_copy('u',sspec(i_spacks)%coreh,is,is,ch,i_spacks)
  !      enddo
  ch=sspec(is)%coreh
  if (ch == ' ') return
  i_copy_size=size(sspec(is)%coreq)
  call dcopy(i_copy_size,sspec(is)%coreq,1,qcore,1)
  read (ch,'(i1)') kcore
  if (ch(2:2) == 's' .OR. ch(2:2) == 'S') lcore = 0
  if (ch(2:2) == 'p' .OR. ch(2:2) == 'P') lcore = 1
  if (ch(2:2) == 'd' .OR. ch(2:2) == 'D') lcore = 2
  if (ch(2:2) == 'f' .OR. ch(2:2) == 'F') lcore = 3
end subroutine gtpcor


