!> MPI utility routines by TK — now a thin wrapper over m_mpi
module m_MPItk
  use m_mpi, only: procid, strprocid, master, nsize, master_mpi, &
                   xmpbnd2, comm, readtk
  implicit none
  public :: m_MPItk_init, procid, strprocid, master, nsize, master_mpi, xmpbnd2, comm, readtk
  private
contains
  subroutine m_MPItk_init(commin)
    use m_mpi, only: MPI__Initialize
    integer, optional :: commin
    call MPI__Initialize(commin)
  end subroutine m_MPItk_init
end module m_MPItk
