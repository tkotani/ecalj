C-----------------------------------------------------------------------
C     mm.f - matrix - matrix multiply, simple self-scheduling version
C-----------------------------------------------------------------------
      subroutine fmain
      implicit none
      include "mpif.h"
      double precision , dimension(:,:) , allocatable :: a,b,c
      double precision , dimension(:) ,   allocatable :: buffer,ans
      integer myid, master, numprocs, ierr, status(MPI_STATUS_SIZE)
      integer i, j, numsent, sender
      integer anstype, row, n
      character*30 str1,str2,string
      character*26 datim
      double precision starttime, endtime
      logical cmdopt, a2bin

C Initialise command line
      call finits(2,0,0,ierr)

      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

      master = 0
C Get dimension
      j = 3
      if (cmdopt('-n=',j,0,string)) then
        if (a2bin(string,n,2,0,' ',j,72)) goto 10
      endif
      if (myid .eq. master) print *, 'Usage: mm -n=#'
      call MPI_FINALIZE(ierr)
      call fexit(1,0,' ',0)
   10 continue 
      allocate (b(1:n,1:n),stat=ierr)
      allocate (buffer(1:n),stat=ierr)
      allocate (ans(1:n),stat=ierr)
      if ( myid .eq. master ) then
        allocate (a(1:n,1:n),stat=ierr)
        allocate (c(1:n,1:n),stat=ierr)
C     master initializes and then dispatches
C     initialize a and b
        do  i = 1,n
          do  j = 1,n
            a(j,i) = i
          enddo
        enddo
        do  i = 1,n
          do  j = 1,n
            b(j,i) = i
          enddo
        enddo
        starttime = MPI_WTIME()
        numsent = 0
C     send b to each other process
        call MPI_BCAST(b, n*n, MPI_DOUBLE_PRECISION, master, 
     .                 MPI_COMM_WORLD, ierr)
C        send a row of a to each other process; tag with row number
        do  i = 1,numprocs-1
          do  j = 1,n
            buffer(j) = a(i,j)
          enddo
          call MPI_SEND(buffer, n, MPI_DOUBLE_PRECISION, i, 
     .      i, MPI_COMM_WORLD, ierr)
          numsent = numsent+1
        enddo

        do  i = 1,n
          call MPI_RECV(ans, n, MPI_DOUBLE_PRECISION, 
     .                  MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
     .                  status, ierr)
          sender     = status(MPI_SOURCE)
          anstype    = status(MPI_TAG)
          do  j = 1,n
            c(anstype,j) = ans(j)
          enddo
          if (numsent .lt. n) then
            do  j = 1,n
              buffer(j) = a(numsent+1,j)
            enddo
            call MPI_SEND(buffer, n, MPI_DOUBLE_PRECISION, 
     .                    sender, numsent+1, MPI_COMM_WORLD, ierr)
            numsent = numsent+1
          else
            call MPI_SEND(1.0, 1, MPI_DOUBLE_PRECISION, sender, 
     .                    0, MPI_COMM_WORLD, ierr)
          endif
        enddo
      else
C     slaves receive B, then compute rows of C until done message
        call MPI_BCAST(b, n*n, MPI_DOUBLE_PRECISION, 
     .                 master, MPI_COMM_WORLD, ierr)
    1   continue 
        call MPI_RECV(buffer, n, MPI_DOUBLE_PRECISION, master, 
     .                MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        if (status(MPI_TAG) .eq. 0) then
          go to 2
        else
          row = status(MPI_TAG)
          do  i = 1,n
            ans(i) = 0
            do  j = 1,n
              ans(i) = ans(i) + buffer(j)*b(j,i)
            enddo
          enddo
          call MPI_SEND(ans, n, MPI_DOUBLE_PRECISION, master, 
     .                  row, MPI_COMM_WORLD, ierr)
          go to 1
        endif
    2   continue
      endif
      if ( myid .eq. master ) then
        endtime = MPI_WTIME()
        write (*,100) endtime-starttime, MPI_WTICK()
  100   format ("Time:",f16.3,"s. Resolution ",g15.4)
      endif
      call MPI_FINALIZE(ierr)
      end

