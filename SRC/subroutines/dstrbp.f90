subroutine dstrbp(nloop,np,nblk,index)
  use m_lgunit,only:stdo
  use m_MPItk,only:master_mpi
  !- Distribute loop counter over processes
  ! ----------------------------------------------------------------------
  !i Inputs:
  !i   nloop: number of times the loop is repeated
  !i   np   : number of processors
  !i   nblk : if equal to zero then a blocking size is chosen
  !i        : if equal to -1 then do double loop balancing (see remarks)
  !o Outputs:
  !o   nblk : vector blocking size
  !o   index: list of pointers
  !r Remarks
  !r   A loop is of the form
  !r       do  i = 1, nloop, nblk
  !r   Each process will do a contiguous number of passes. Index(n)
  !r   is the counter that process n starts with. That is process n does
  !r       do i = index(n), index(n+1)-1, nblk
  !r   nblk is chosen (in the range 6-16) to give the best distribution.
  !r   The work is divided as evenly as possible over the processes.
  !r   Double loop balancing is for the case where the serial code is,
  !r       do i = 1, n
  !r       do j = i, n
  !r   In that case index is a 2D array and the loop is transformed into
  !r       do iloop = 1, index(procid,0)
  !r       i = index(procid,iloop)
  !r       do j = i, n
  !r   This distributes the work over i so that process 0 does 1 and n
  !r   process 1 does 2 and n-1 process 2 does 3 and n-2 and so on.
  !r   (Thanks to Andrew Canning)
  !u Updates
  !u   Written by Karen Johnston
  ! ----------------------------------------------------------------------
  IMPLICIT NONE
  ! Passed
  INTEGER :: nloop,np,nblk
  INTEGER, DIMENSION(0:np) :: index
  ! Local
  INTEGER :: i,step
  INTEGER, DIMENSION(0:np-1) :: xnode,inode
  integer :: iprint

  if (nblk == -1) then
     call pdstlb(nloop,np,index)
     return
  endif
  step = nblk
  ! Initialising arrays
  DO i=0,np-1
     inode(i)=0
     xnode(i)=0
  END DO
  IF (step==1) THEN
     CALL single(nloop,np,xnode,inode)
  ELSE IF (step == 0) THEN
     ! Find optimum step nblk 6-16
     CALL optimise(nloop,np,step,xnode,inode,nblk)
  ELSE
     CALL multiple(nloop,np,step,xnode,inode)
  END IF
  DO i=0,np-1
     index(i)=inode(i)+1
  END DO
  index(np)=nloop+1
  if (iprint() < 41) return
  ! Output indices
  if(master_mpi)WRITE (*,'(" DSTRBP: i   inode   xnode   block=",i5)') nblk
  DO i=0,np-1
     if(master_mpi) WRITE (*,'(3x,3i7)') i,inode(i),xnode(i)
  END DO
END SUBROUTINE dstrbp

SUBROUTINE single(nloop,np,xnode,inode)

  IMPLICIT NONE
  INTEGER :: i=0,j=0,nloop,np
  INTEGER :: min=0,rem=0
  INTEGER :: times=0,rem2=0,total=0
  INTEGER, DIMENSION(0:np-1) :: inode,xnode
  integer :: iprint

  !      if (iprint() .gt. 40) WRITE(LGUNIT(1),*) 'DSTRBP SINGLE:'
  ! Split nblock evenly with nodes
  times = nloop/np
  rem2 = MOD(nloop,np)
  call info5(41,0,0, &
       ' DSTRBP, single:  nloop = %i  np = %i  times = %i  rem2 = %i', &
       nloop,np,times,rem2,0)
  !      if (iprint() .gt. 40)
  !     .  WRITE(LGUNIT(1),*) 'nloop=',nloop,'np=',np,'times=',times
  !     .  ,'rem2=',rem2
  ! Even no. of kpts per node without remainder
  inode(0)=0
  DO i=1,np-1
     inode(i)=inode(i-1)+times
  END DO
  ! Spread over the remainder
  IF (rem2 /= 0) THEN
     DO i=1,rem2
        DO j=0,np-1
           IF (j>i) THEN
              inode(j) = inode(j)+1
           END IF
        END DO
     END DO
  END IF
  ! Number of blocks in each node
  total = 0
  DO i=0,np-2
     xnode(i) = inode(i+1)-inode(i)
     total = total + xnode(i)
  END DO
  xnode(np-1) = nloop - total
  DO i=1,np-1
     inode(i) = inode(i-1) + xnode(i-1)
  END DO

END SUBROUTINE single



SUBROUTINE optimise(nloop,np,step,xnode,inode,size)
  use m_lgunit,only:stdo

  IMPLICIT NONE
  INTEGER, PARAMETER :: stepmin=6,stepmax=16
  INTEGER :: i=0,j=0,nloop,np,step
  INTEGER :: size,rem=0
  INTEGER :: nblock=0,rem2=0,total=0
  INTEGER :: times=0,rem3=0
  INTEGER, DIMENSION(0:np-1) :: inode,xnode
  INTEGER, DIMENSION(0:nloop-1) :: iblock,xblock
  integer :: iprint

  if (iprint() > 40) WRITE(stdo,*) 'DSTRBP OPTIMISE:'
  ! Optimises the block size
  size = stepmin
  DO i=stepmin+1,stepmax
     IF ( ABS(MOD(nloop,i)) <= ABS(MOD(nloop,size))) THEN
        size = i
     END IF
  END DO
  rem = MOD(nloop,size)
  ! Split nloop into blocks
  nblock = nloop/size
  rem2 = MOD(nloop,size)
  if (iprint() > 40) then
     WRITE(stdo,*) 'size=',size,'nblock=',nblock,'rem2=',rem2
     IF (nblock < np) THEN
        WRITE(stdo,*) ""
        WRITE(stdo,*) "*** WARNING: No. blocks < no. nodes ***"
        WRITE(stdo,*) ""
     END IF
  endif
  ! Even no. of kpts per block without remainder
  iblock(0)=0
  DO i=1,nblock-1
     iblock(i)=iblock(i-1)+size
  END DO
  ! Spread over the remainder
  IF (rem2 /= 0) THEN
     DO i=0,rem2-1
        DO j=0,nblock-1
           IF (j>i) THEN
              iblock(j) = iblock(j)+1
           END IF
        END DO
     END DO
  END IF
  ! Number of nloop in each block
  DO i=0,nblock-2
     xblock(i) = iblock(i+1)-iblock(i)
     total = total + xblock(i)
  END DO
  xblock(nblock-1) = nloop - total
  DO i=1,nblock-1
     iblock(i) = iblock(i-1) + xblock(i-1)
  END DO
  ! Split blocks over nodes
  times = nblock/np
  rem3 = MOD(nblock,np)
  if (iprint() > 40) &
       WRITE(stdo,*) 'nblock=',nblock,'np=',np,'times=',times &
       ,'rem3=',rem3

  ! Even no. of blocks per node without remainder
  inode(0)=iblock(0)
  DO i=1,np-1
     inode(i)=iblock(i*times)
  END DO
  ! Spread over the remainder
  IF (rem3 /= 0) THEN
     DO i=0,rem3-2
        DO j=np-1-i,np-1
           inode(j) = inode(j)+size
        END DO
     END DO
  END IF
  if (iprint() > 40) then
     DO i=0,np-1
        WRITE(stdo,*) 'inode(',i,')=',inode(i)
     END DO
  endif
  DO i=0,np-2
     xnode(i)=inode(i+1)-inode(i)
  END DO
  xnode(np-1) = nloop-inode(np-1)

END SUBROUTINE optimise


SUBROUTINE multiple(nloop,np,step,xnode,inode)
  use m_lgunit,only:stdo

  IMPLICIT NONE
  INTEGER :: i=0,j=0,nloop,np,step
  INTEGER :: size=0,rem=0
  INTEGER :: nblock=0,rem2=0,total=0
  INTEGER :: times=0,rem3=0
  INTEGER, DIMENSION(0:np-1) :: inode,xnode
  INTEGER, DIMENSION(0:nloop-1) :: iblock,xblock
  integer :: iprint

  if (iprint() > 40) WRITE(stdo,*) 'DSTRBP MULTIPLE:'
  size = step
  rem = MOD(nloop,size)
  ! Split nloop into blocks
  nblock = nloop/size
  rem2 = MOD(nloop,size)
  if (iprint() > 40) then
     WRITE(stdo,*) 'size=',size,'nblock=',nblock,'rem2=',rem2
     IF (nblock < np) THEN
        WRITE(stdo,*) ""
        WRITE(stdo,*) "*** WARNING: No. blocks < no. nodes ***"
        WRITE(stdo,*) ""
     END IF
  endif
  ! Even no. of kpts per block without remainder
  iblock(0)=0
  DO i=1,nblock-1
     iblock(i)=iblock(i-1)+size
  END DO
  ! Spread over the remainder
  IF (rem2 /= 0) THEN
     DO i=0,rem2-1
        DO j=0,nblock-1
           IF (j>i) THEN
              iblock(j) = iblock(j)+1
           END IF
        END DO
     END DO
  END IF
  ! Number of nloop in each block
  DO i=0,nblock-2
     xblock(i) = iblock(i+1)-iblock(i)
     total = total + xblock(i)
  END DO
  xblock(nblock-1) = nloop - total
  DO i=1,nblock-1
     iblock(i) = iblock(i-1) + xblock(i-1)
  END DO
  ! Split blocks over nodes
  times = nblock/np
  rem3 = MOD(nblock,np)
  if (iprint() > 40) &
       WRITE(stdo,*) 'nblock=',nblock,'np=',np,'times=',times &
       ,'rem3=',rem3

  ! Even no. of blocks per node without remainder
  inode(0)=iblock(0)
  DO i=1,np-1
     inode(i)=iblock(i*times)
  END DO
  ! Spread over the remainder
  IF (rem3 /= 0) THEN
     DO i=0,rem3-2
        DO j=np-1-i,np-1
           inode(j) = inode(j)+size
        END DO
     END DO
  END IF
  DO i=0,np-2
     xnode(i)=inode(i+1)-inode(i)
  END DO
  xnode(np-1) = nloop-inode(np-1)

END SUBROUTINE multiple

subroutine pdstlb(nloop,np,index)
  use m_lgunit,only:stdo
  implicit none
  integer :: nloop, np, index(0:np-1,0:*)
  integer :: n, p, m, i, j, sum1, sum2, iprint
  n = nloop
  do i = 0, np-1
     index(i,0) = 0
  enddo
  p = 0
  do  i = 0, n/2-1
     j = index(p,0) + 1
     index(p,j)   = 1 + i
     index(p,j+1) = n - i
     index(p,0) = index(p,0) + 2
     p = mod(p+1,np)
  enddo
  if (mod(n,2) /= 0) then
     j = index(p,0) + 1
     index(p,j) = n/2+1
     index(p,0) = index(p,0) + 1
  endif
  sum1 = 0
  sum2 = 0
  do  i = 1, n
     sum1 = sum1 + i
  enddo
  do i = 0, np-1
     do j= 1, index(i,0)
        sum2 = sum2 + index(i,j)
     enddo
  enddo
  if (iprint() > 40) then
     write (*,1)
     do i = 0, np-1
        write (stdo,2) i,(index(i,j),j=0,index(i,0))
     enddo
  endif
  if (sum1 /= sum2) then
     if (iprint() > 0) &
          write(stdo,"('Bug in pdstlb: sum1=',i0,' sum2=',i0)") sum1,sum2
     call rx0(' ')
  endif
1 format (/' PDSTLB:'/' proc   total    loop number')
2 format (i5,3x,i5,3x,256i5)
end subroutine pdstlb
!      program test
!      implicit none
!      integer nloop,np,index(10000000),i,j
!      write(*, 10)
!   10 format ('Enter nloop,np '$)
!      read (*,*) nloop, np
!      write (*,20) nloop, np
!   20 format ('nloop=',i5,' np=',i5)
!      call pdstlb(nloop,np,index)
!      end
