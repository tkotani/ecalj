program main
  use mpi
  use m_ctrl2ctrlp, only: ConvertCtrl2CtrlpByPython
  use m_cmdpath,    only: setcmdpath, cmdpath
  use m_args,       only: m_setargs, argall
  use m_ext,        only: m_ext_init, sname
  use m_mpi,        only: MPI__Initialize, mpi__root, comm, setipr, mpi__rank
  use m_mlo_ovlppair, only: init_build_ovlppair, build_ovlppair_q
  integer :: ierr

  call MPI__Initialize()
  call setipr(comm)
  call setcmdpath()
  call m_setargs()
  call m_ext_init()
  if(mpi__root) then
    print *,'cmdpath:', trim(cmdpath)
    print *,'args:',    trim(argall)
    print *,'ext:',     trim(sname)
    call ConvertCtrl2CtrlpByPython()
  endif
  call mpi_barrier(comm, ierr)
  call init_build_ovlppair(comm)
  call build_ovlppair_q(spinflip=.true.)
  call rx0('Ok')
end program
