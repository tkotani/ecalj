!>    mpi utility for wannier part R. Sakuma 2007
module rsmpi
  implicit none
  include "mpif.h"
  ! S:     FILES
  ! S:
  ! S:     nrsin: contains additional parameters
  ! S:     nrphistar: contains wavefunctions
  !     integer,parameter:: nrsin=1001,nrsj2gsmall=1002

  ! MPI parameters
  integer :: myrank_rsmpi   ! rank of the calling process
  integer :: nproc_rsmpi    ! number of processes
  integer :: io_root_rsmpi
  parameter (io_root_rsmpi = 0)
  integer :: ierror_rsmpi   ! error check

  ! ID of the current process (ex. 0000012, 0123456)
  ! used mainly for output filename (ex. VCCFP.RSMPI0000012 )
  ! if # of processes >= 10^8, increase the size of the array
  character(7) :: myrank_id_rsmpi !
  integer :: ifile_rsmpi    ! file id
  ! used for collective I/O
  integer :: bufsize_rsmpi
  parameter (bufsize_rsmpi=1024)
  character*(bufsize_rsmpi) :: buf_rsmpi

  ! measure elapsed time
  double precision :: t1,t2
contains
  !------------------------------------------------------
  subroutine RSMPI_Init()
    implicit none
    call MPI_INIT(ierror_rsmpi)
    call RSMPI_Check("MPI_INIT",ierror_rsmpi)

    t1 = MPI_WTIME()

    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc_rsmpi,ierror_rsmpi)
    call RSMPI_Check("MPI_COMM_SIZE",ierror_rsmpi)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myrank_rsmpi,ierror_rsmpi)
    call RSMPI_Check("MPI_COMM_RANK",ierror_rsmpi)
    if (Is_IO_Root_RSMPI()) then
       write(6,*) "RS: --- RSMPI_Init ---"
       write(6,*) "RS: Number of processes : ",nproc_rsmpi
    end if
    call RSMPI_Set_ID()
    !     if (Is_IO_Root_RSMPI()) then
    !      write(6,*) "RS: ID : ", myrank_id_rsmpi
    !     end if
    write(buf_rsmpi,*) "RS: RANKID = ",myrank_id_rsmpi
    call RSMPI_Write(6)
    !      write(6,*) "RS: ID : ", myrank_id_rsmpi

    !      call MPI_Barrier(MPI_COMM_WORLD,ierror_rsmpi)
    !      call RSMPI_Check("MPI_Barrier",ierror_rsmpi)

    if (Is_IO_Root_RSMPI()) then
       write(6,*) "RS: --- end of RSMPI_Init ---"
    end if
  end subroutine RSMPI_Init

  !--------------------------------------------------------
  subroutine RSMPI_Set_ID
    implicit none

    if (myrank_rsmpi/1000000 >= 10) then
       ! RS: I don't expect this will happen..
       call RSMPI_Stop("RSMPI_Set_ID: # of processes exceeds 10^6: Modify RSMPI_mod.F!")
    endif
    write(myrank_id_rsmpi,'(i7.7)') myrank_rsmpi
  end subroutine RSMPI_Set_ID
  !--------------------------------------------------------
  subroutine RSMPI_Check(func,ierror_in)
    implicit none
    character*(*) :: func
    integer,intent(in) :: ierror_in
    integer :: ireturn
    if ((ierror_in /= MPI_SUCCESS)) then
       !     if (Is_IO_Root_RSMPI()) then
       write(6,*) "RS: MPI ERROR :", func," PROCID = ",myrank_id_rsmpi," ierror =",ierror_in,".Aborted."
       ireturn = 99
       call MPI_ABORT(MPI_COMM_WORLD,ireturn,ierror_rsmpi)
    endif
  end subroutine RSMPI_Check
  !--------------------------------------------------------

  subroutine RSMPI_Stop(msg)
    implicit none
    character*(*) :: msg
    !    if (Is_IO_Root_RSMPI()) then
    write(6,*) "RS: Error: ",msg, " Aborted."
    !    endif
    call MPI_ABORT(MPI_COMM_WORLD,99,ierror_rsmpi)
  end subroutine RSMPI_Stop
  !--------------------------------------------------------
  subroutine RSMPI_Print_WTime()
    implicit none
    t2 = MPI_WTIME()
    if (Is_IO_Root_RSMPI()) then
       write(6,*) "RS: Elapsed time: ",t2-t1,"sec"
       write(6,*) "RS: Precision   : ",MPI_WTICK()
    endif
  end subroutine RSMPI_Print_WTime
  !--------------------------------------------------------
  subroutine RSMPI_Finalize()
    implicit none
    call MPI_Barrier(MPI_COMM_WORLD,ierror_rsmpi)
    call RSMPI_Check("MPI_Barrier",ierror_rsmpi)

    call RSMPI_Print_WTime()

    call MPI_FINALIZE(ierror_rsmpi)
    if (ierror_rsmpi /= MPI_SUCCESS) then
       write(6,*) "RS: MPI ERROR in RSMPI_Finalize(), ierror =", &
            ierror_rsmpi, "PROCID = ",myrank_id_rsmpi
    endif
  end subroutine RSMPI_Finalize
  !--------------------------------------------------------
  logical function Is_IO_Root_RSMPI()
    implicit none
    Is_IO_Root_RSMPI = (myrank_rsmpi .eq. io_root_rsmpi)
  end function Is_IO_Root_RSMPI
  !--------------------------------------------------------
  ! Gather information and Write

  subroutine RSMPI_Write(if_in)
    implicit none
    integer,intent(in) :: if_in
    ! function
    integer :: get_non_blank_rsmpi
    ! local
    integer :: ip !process
    integer :: tag
    integer :: status(MPI_STATUS_SIZE) ! MPI_Status
    integer :: isize
    tag = 0

    if ( .NOT. Is_IO_Root_RSMPI()) then
       call MPI_Send(buf_rsmpi,bufsize_rsmpi,MPI_CHARACTER, &
            io_root_rsmpi,tag,MPI_COMM_WORLD,ierror_rsmpi)
       call RSMPI_Check("RSMPI_Write",ierror_rsmpi)
    else
       do ip=0,nproc_rsmpi-1
          if (ip /= myrank_rsmpi) then
             call MPI_Recv(buf_rsmpi,bufsize_rsmpi,MPI_CHARACTER, &
                  ip,tag,MPI_COMM_WORLD,status,ierror_rsmpi)
             call RSMPI_Check("RSMPI_Write",ierror_rsmpi)
          endif
          isize = get_non_blank_rsmpi(bufsize_rsmpi,buf_rsmpi)
          if (isize /= 0) write(if_in,*) buf_rsmpi(1:isize)
       enddo
    endif

  end subroutine RSMPI_Write
  !--------------------------------------------------------
end module RSMPI
