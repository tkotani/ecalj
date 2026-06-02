program main
  use mpi
  use m_cmdpath,    only: setcmdpath, cmdpath
  use m_args,       only: m_setargs, argall
  use m_ext,        only: m_ext_init, sname
  use m_mpi,        only: MPI__Initialize, MPI__Broadcast, mpi__root, comm, setipr, mpi__rank
  use m_mlo_ovlppair, only: init_build_ovlppair, build_ovlppair_q
  use m_cmdopt_registry, only: c2_sp1, c2_sp2
  integer :: ierr, isp1, isp2

  call MPI__Initialize()
  call setipr(comm)
  call setcmdpath()
  call m_setargs()
  call m_ext_init()
  call mpi_barrier(comm, ierr)
  isp1 = 2; isp2 = 1  ! default DNUP
  if(mpi__root) then
    if (c2_sp1 >= 0) isp1 = c2_sp1
    if (c2_sp2 >= 0) isp2 = c2_sp2
  endif
  call MPI__Broadcast(isp1)
  call MPI__Broadcast(isp2)
  call init_build_ovlppair(comm)
  call build_ovlppair_q(isp1, isp2)
  call rx0('Ok')
end program
