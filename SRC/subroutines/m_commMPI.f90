module m_comm
  integer,protected:: comm,size,rank
contains
  subroutine setcomm(commin) bind(C)
    implicit none
    include "mpif.h"
    integer(4) :: ierr,n,commin
    call MPI_Comm_size(commin, size, ierr)
    call MPI_Comm_rank(commin, rank, ierr)
    comm=commin
  end subroutine setcomm
end module m_comm
