!>Free atom density determined by lmfa. lmf stores data into rst.*. But unchanged.
module m_fatom
  use m_MPItk,only: comm
  integer,parameter::  n0=10
  type s_spec
     ! I think lmfa detemines all the following data and write to atm.* files
     ! The data is used for lmf-MPIK and copied into rst file (unchanged).
     real(8):: qc    !core charge
     real(8):: rsmfa !rsm to fit free atom density
     real(8):: ctail !coefficients to fit of free-atom core tail by unsm. Hankel
     real(8):: etail !energy to fit of free-atom core tail
     real(8):: stc   !core kinetic energy
     integer:: nxi    ! Number of energies in fit of free-atom density tails
     real(8):: exi(n0)    ! Hankel energies for fit to c.d.; fit to free-atom density tails.
     real(8):: chfa(n0,2) ! coefficients to fit of free-atom density tails
     real(8),allocatable :: rv_a_orhoc(:) !pointer to core density
  end type s_spec
  type(s_spec),allocatable:: sspec(:) !just allocated for iors and rdovfa. Not touched.
contains
  subroutine mpibc1_s_spec(ssp)
    implicit none
    include 'mpif.h'
    type(s_spec):: ssp
    integer :: master=0,ierr
    call mpi_bcast(ssp%ctail, 1,MPI_REAL8 , master, comm,ierr)
    call mpi_bcast(ssp%etail, 1,MPI_REAL8 , master, comm,ierr)
    call mpi_bcast(ssp%stc,   1,MPI_REAL8   , master, comm,ierr)
    call mpi_bcast(ssp%nxi,   1,MPI_INTEGER , master, comm,ierr)
    call mpi_bcast(ssp%qc,    1,MPI_REAL8    , master, comm,ierr)
    call mpi_bcast(ssp%exi, size(ssp%exi), MPI_REAL8 , master, comm,ierr) 
    call mpi_bcast(ssp%chfa,size(ssp%chfa),MPI_REAL8, master, comm,ierr)
  end subroutine mpibc1_s_spec
end module m_fatom
