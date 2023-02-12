subroutine gtpcor(is,kcore,lcore,qcore)!- Unpacks parameters related to partial core occpation
  use m_struc_def 
  use m_lmfinit,only:coreq,coreh
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
  integer :: is,kcore,lcore,i_copy_size
  real(8):: qcore(2)
  character(8) :: ch
  kcore = 0
  lcore = -1
  qcore = 0d0
  ch=coreh(is)
  if (len(trim(adjustl(ch))) == 0) return
  qcore=coreq(:,is)
  read (ch,'(i1)') kcore
  if (ch(2:2) == 's' .OR. ch(2:2) == 'S') lcore = 0
  if (ch(2:2) == 'p' .OR. ch(2:2) == 'P') lcore = 1
  if (ch(2:2) == 'd' .OR. ch(2:2) == 'D') lcore = 2
  if (ch(2:2) == 'f' .OR. ch(2:2) == 'F') lcore = 3
end subroutine gtpcor


