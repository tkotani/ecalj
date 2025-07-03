program main
  use mpi
  use m_ctrl2ctrlp,only: ConvertCtrl2ctrlpBypython
  use m_cmdpath,only: setcmdpath
  use m_args,only:    m_setargs
  use m_ext,only:     m_ext_init
  use m_lmf,only: lmf
  use m_mpi,only: setipr
  integer:: ierr,comm,procid
  call mpi_init(ierr)
  comm = MPI_COMM_WORLD
!  call MPI_Comm_set_errhandler(comm, MPI_ERRORS_RETURN, ierr)
  call mpi_comm_rank(comm, procid, ierr )
  call setipr(comm)
  call setcmdpath()
  call m_setargs()
  call m_ext_init()    ! Get sname, e.g. trim(sname)=si of ctrl.si
  if(procid==0) call ConvertCtrl2CtrlpByPython()
  call lmf()
  call mpi_finalize(ierr)
  call exit(0)
end program main
