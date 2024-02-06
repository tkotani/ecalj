module m_MPIcomm
  integer,protected:: comm,size,rank
contains
  subroutine MPIsetcomm(commin) bind(C) !this can be called from python
    implicit none
    include "mpif.h"
    integer(4) :: ierr,n,commin
    call MPI_Comm_size(commin, size, ierr)
    call MPI_Comm_rank(commin, rank, ierr)
    comm=commin
  end subroutine MPIsetcomm
end module m_MPIcomm
