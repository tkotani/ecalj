!> this is for the case that main routine is python (lmf.py)
subroutine sopen(output) bind(c)
  character(1):: output(*)
  character(1024):: convcchar
  open(6,file=trim(convcchar(output)))
end subroutine sopen
subroutine sclose() bind(c)
  use m_lgunit, only: m_lgunit_reset
  close(6)
  call m_lgunit_reset()
end subroutine sclose
