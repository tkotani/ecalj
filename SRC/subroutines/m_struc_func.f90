Module m_struc_func
  use m_ftox,only: ftox
  use m_lgunit,only:stml
  public  mpibc1_s_spec, mpibc1_s_site
contains
  subroutine mpibc1_s_spec(struc,label)
    use m_struc_def, only: s_spec
    use m_MPItk,only: mlog
    implicit none
    type(s_spec):: struc
    character label*(*)
!#if MPI|MPIK
    include 'mpif.h'
    integer :: numprocs, ierr
    integer :: MAX_PROCS
    parameter (MAX_PROCS = 100)
    integer :: resultlen
    character*(MPI_MAX_PROCESSOR_NAME) name
    character(10) :: shortname(0:MAX_PROCS-1)
    character(26) :: datim
    integer :: namelen(0:MAX_PROCS-1)
    character(256) :: strn
    !      integer:: lgunit
    integer :: procid,master,i_data_size
    integer:: n=0
    master=0
    call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
    call mpi_bcast(struc%z, 1,MPI_REAL8 &
         , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%rmt, 1,MPI_REAL8 &
         , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%rsmfa, 1,MPI_REAL8 &
         , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%rsma, 1,MPI_REAL8 &
         , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%rg, 1,MPI_REAL8 &
         , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%lmxa, 1,MPI_INTEGER &
         , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%lmxl, 1,MPI_INTEGER &
         , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%kmxt, 1,MPI_INTEGER &
         , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%rsmv, 1,MPI_REAL8 &
         , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%coreh, 1,MPI_REAL8 &
         , master, MPI_COMM_WORLD,ierr)
    i_data_size=size(struc%coreq) ! coreq(2)
    call mpi_bcast(struc%coreq, i_data_size,MPI_REAL8 &
         , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%a, 1,MPI_REAL8 &
         , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%nr, 1,MPI_INTEGER &
         , master, MPI_COMM_WORLD,ierr)
!    call mpi_bcast(struc%eref, 1,MPI_REAL8 &
!         , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%lfoca, 1,MPI_INTEGER &
         , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%ctail, 1,MPI_REAL8 &
         , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%etail, 1,MPI_REAL8 &
         , master, MPI_COMM_WORLD,ierr)
!    call mpi_bcast(struc%name, 1,MPI_REAL8 &
!         , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%stc, 1,MPI_REAL8 &
         , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%lmxb, 1,MPI_INTEGER &
         , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%rfoca, 1,MPI_REAL8 &
         , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%nxi, 1,MPI_INTEGER &
         , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%qc, 1,MPI_REAL8 &
         , master, MPI_COMM_WORLD,ierr)
!    call mpi_bcast(struc%eh3, 1,MPI_REAL8 &
!         , master, MPI_COMM_WORLD,ierr)
!    call mpi_bcast(struc%rs3, 1,MPI_REAL8 &
!         , master, MPI_COMM_WORLD,ierr)
!    call mpi_bcast(struc%vmtz, 1,MPI_REAL8 &
!         , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%kmxv, 1,MPI_INTEGER &
         , master, MPI_COMM_WORLD,ierr)
    i_data_size=size(struc%rcfa) ! rcfa(2)
    call mpi_bcast(struc%rcfa, i_data_size,MPI_REAL8 &
         , master, MPI_COMM_WORLD,ierr)
    i_data_size=size(struc%p) ! p(20)
    call mpi_bcast(struc%p, i_data_size,MPI_REAL8 &
         , master, MPI_COMM_WORLD,ierr)
    i_data_size=size(struc%q) ! q(20)
    call mpi_bcast(struc%q, i_data_size,MPI_REAL8 &
         , master, MPI_COMM_WORLD,ierr)
!    i_data_size=size(struc%idmod) ! idmod(10)
!    call mpi_bcast(struc%idmod, i_data_size,MPI_INTEGER &
!         , master, MPI_COMM_WORLD,ierr)
    !      i_data_size=size(struc%idxdn) ! idxdn(30)
    !      call mpi_bcast(struc%idxdn, i_data_size,MPI_INTEGER
    !     , , master, MPI_COMM_WORLD,ierr)
    i_data_size=size(struc%exi) ! exi(10)
    call mpi_bcast(struc%exi, i_data_size,MPI_REAL8 &
         , master, MPI_COMM_WORLD,ierr)
    i_data_size=size(struc%chfa) ! chfa(20)
    call mpi_bcast(struc%chfa, i_data_size,MPI_REAL8 &
         , master, MPI_COMM_WORLD,ierr)
    i_data_size=size(struc%pz) ! pz(20)
    call mpi_bcast(struc%pz, i_data_size,MPI_REAL8 &
         , master, MPI_COMM_WORLD,ierr)
!    i_data_size=size(struc%idu) ! idu(4)
!    call mpi_bcast(struc%idu, i_data_size,MPI_INTEGER &
!         , master, MPI_COMM_WORLD,ierr)
!    i_data_size=size(struc%uh) ! uh(4)
!    call mpi_bcast(struc%uh, i_data_size,MPI_REAL8 &
!         , master, MPI_COMM_WORLD,ierr)
!    i_data_size=size(struc%jh) ! jh(4)
!    call mpi_bcast(struc%jh, i_data_size,MPI_REAL8 &
!         , master, MPI_COMM_WORLD,ierr)
    if (mlog) then
       call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
       call strcop(shortname(procid),name,10,'.',ierr)
       namelen(procid) = ierr-1
       call gettime(datim)
       write(stml,ftox)'  '//datim//' Process ',procid,' of ',numprocs,' on ' &
            //shortname(procid)(1:namelen(procid))//' '//label
    endif
!#endif
  end subroutine mpibc1_s_spec
  ! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine mpibc1_s_site(struc,label)
    use m_struc_def, only: s_site
    use m_MPItk,only: mlog
    implicit none
    type(s_site):: struc
    character label*(*)
!#if MPI|MPIK
    include 'mpif.h'
    integer :: numprocs, ierr
    integer :: MAX_PROCS
    parameter (MAX_PROCS = 100)
    integer :: resultlen
    character*(MPI_MAX_PROCESSOR_NAME) name
    character(10) :: shortname(0:MAX_PROCS-1)
    character(26) :: datim
    integer :: namelen(0:MAX_PROCS-1)
    character(256) :: strn
    !      integer:: lgunit
    integer :: procid,master,i_data_size
    integer:: n=0
    master=0
    call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
    !      call mpi_bcast(struc%clabel, 1,MPI_REAL8 , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%class, 1,MPI_INTEGER , master, MPI_COMM_WORLD,ierr)
    i_data_size=size(struc%force) ! force(3)
    call mpi_bcast(struc%force, i_data_size,MPI_REAL8 , master, MPI_COMM_WORLD,ierr)
    i_data_size=size(struc%pnu) ! pnu(20)
    call mpi_bcast(struc%pnu, i_data_size,MPI_REAL8 , master, MPI_COMM_WORLD,ierr)
    i_data_size=size(struc%pos) ! pos(3)
    call mpi_bcast(struc%pos, i_data_size,MPI_REAL8 , master, MPI_COMM_WORLD,ierr)
    i_data_size=size(struc%pos0) ! pos0(3)
    call mpi_bcast(struc%pos0, i_data_size,MPI_REAL8 , master, MPI_COMM_WORLD,ierr)
    i_data_size=size(struc%pz) ! pz(20)
    call mpi_bcast(struc%pz, i_data_size,MPI_REAL8 , master, MPI_COMM_WORLD,ierr)
    i_data_size=size(struc%relax) ! relax(3)
    call mpi_bcast(struc%relax, i_data_size,MPI_INTEGER , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%spec, 1,MPI_INTEGER , master, MPI_COMM_WORLD,ierr)
    !      i_data_size=size(struc%vel) ! vel(3)
    !      call mpi_bcast(struc%vel, i_data_size,MPI_REAL8 , master, MPI_COMM_WORLD,ierr)
    !       call mpi_bcast(struc%vshft, 1,MPI_REAL8 , master, MPI_COMM_WORLD,ierr)
    if (mlog) then
       call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
       call strcop(shortname(procid),name,10,'.',ierr)
       namelen(procid) = ierr-1
       call gettime(datim)
       write(stml,ftox)'  '//datim//' Process ',procid,' of ',numprocs,' on ' &
            //shortname(procid)(1:namelen(procid))//' '//label
    endif
!#endif
  end subroutine mpibc1_s_site
end module m_struc_func
