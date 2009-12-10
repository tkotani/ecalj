      subroutine dstrbp(nloop,np,nblk,index)
C- Distribute loop counter over processes
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nloop: number of times the loop is repeated
Ci   np   : number of processors
Ci   nblk : if equal to zero then a blocking size is chosen
Ci        : if equal to -1 then do double loop balancing (see remarks)
Co Outputs:
Co   nblk : vector blocking size
Co   index: list of pointers
Cr Remarks
Cr   A loop is of the form
Cr       do  i = 1, nloop, nblk
Cr   Each process will do a contiguous number of passes. Index(n)
Cr   is the counter that process n starts with. That is process n does
Cr       do i = index(n), index(n+1)-1, nblk
Cr   nblk is chosen (in the range 6-16) to give the best distribution.
Cr   The work is divided as evenly as possible over the processes.
Cr   Double loop balancing is for the case where the serial code is,
Cr       do i = 1, n
Cr       do j = i, n
Cr   In that case index is a 2D array and the loop is transformed into
Cr       do iloop = 1, index(procid,0)
Cr       i = index(procid,iloop)
Cr       do j = i, n
Cr   This distributes the work over i so that process 0 does 1 and n
Cr   process 1 does 2 and n-1 process 2 does 3 and n-2 and so on.
Cr   (Thanks to Andrew Canning)
Cu Updates
Cu   Written by Karen Johnston
C ----------------------------------------------------------------------
#if MPI | MPIK
      IMPLICIT NONE
C Passed     
      INTEGER :: nloop,np,nblk
      INTEGER, DIMENSION(0:np) :: index
C Local
      INTEGER :: i,step
      INTEGER, DIMENSION(0:np-1) :: xnode,inode
      integer iprint, lgunit

      if (nblk .eq. -1) then
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
      if (iprint() .lt. 41) return
! Output indices
      WRITE (*,'(" DSTRBP: i   inode   xnode   block=",i5)') nblk
      DO i=0,np-1
        WRITE (*,'(3x,3i7)') i,inode(i),xnode(i)
      END DO

      END SUBROUTINE dstrbp
#else
      end
#endif 

#if MPI | MPIK
      SUBROUTINE single(nloop,np,xnode,inode)

      IMPLICIT NONE
      INTEGER :: i=0,j=0,nloop,np
      INTEGER :: min=0,rem=0
      INTEGER :: times=0,rem2=0,total=0
      INTEGER, DIMENSION(0:np-1) :: inode,xnode
      integer iprint,lgunit

C      if (iprint() .gt. 40) WRITE(LGUNIT(1),*) 'DSTRBP SINGLE:'
! Split nblock evenly with nodes
      times = nloop/np
      rem2 = MOD(nloop,np)
      call info5(41,0,0,
     .  ' DSTRBP, single:  nloop = %i  np = %i  times = %i  rem2 = %i',
     .  nloop,np,times,rem2,0)
C      if (iprint() .gt. 40)
C     .  WRITE(LGUNIT(1),*) 'nloop=',nloop,'np=',np,'times=',times
C     .  ,'rem2=',rem2
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

      IMPLICIT NONE
      INTEGER, PARAMETER :: stepmin=6,stepmax=16
      INTEGER :: i=0,j=0,nloop,np,step
      INTEGER :: size,rem=0
      INTEGER :: nblock=0,rem2=0,total=0
      INTEGER :: times=0,rem3=0
      INTEGER, DIMENSION(0:np-1) :: inode,xnode
      INTEGER, DIMENSION(0:nloop-1) :: iblock,xblock
      integer iprint,lgunit

      if (iprint() .gt. 40) WRITE(LGUNIT(1),*) 'DSTRBP OPTIMISE:'
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
      if (iprint() .gt. 40) then
        WRITE(LGUNIT(1),*) 'size=',size,'nblock=',nblock,'rem2=',rem2
        IF (nblock < np) THEN
          WRITE(LGUNIT(1),*) ""
          WRITE(LGUNIT(1),*) "*** WARNING: No. blocks < no. nodes ***"
          WRITE(LGUNIT(1),*) ""
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
      if (iprint() .gt. 40)
     .  WRITE(LGUNIT(1),*) 'nblock=',nblock,'np=',np,'times=',times
     .  ,'rem3=',rem3

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
      if (iprint() .gt. 40) then
        DO i=0,np-1
          WRITE(LGUNIT(1),*) 'inode(',i,')=',inode(i)
        END DO
      endif
      DO i=0,np-2
        xnode(i)=inode(i+1)-inode(i)
      END DO
      xnode(np-1) = nloop-inode(np-1)

      END SUBROUTINE optimise


      SUBROUTINE multiple(nloop,np,step,xnode,inode)

      IMPLICIT NONE
      INTEGER :: i=0,j=0,nloop,np,step
      INTEGER :: size=0,rem=0
      INTEGER :: nblock=0,rem2=0,total=0
      INTEGER :: times=0,rem3=0
      INTEGER, DIMENSION(0:np-1) :: inode,xnode
      INTEGER, DIMENSION(0:nloop-1) :: iblock,xblock
      integer iprint,lgunit

      if (iprint() .gt. 40) WRITE(LGUNIT(1),*) 'DSTRBP MULTIPLE:'
      size = step
      rem = MOD(nloop,size)
! Split nloop into blocks
      nblock = nloop/size
      rem2 = MOD(nloop,size)
      if (iprint() .gt. 40) then
        WRITE(LGUNIT(1),*) 'size=',size,'nblock=',nblock,'rem2=',rem2
        IF (nblock < np) THEN
          WRITE(LGUNIT(1),*) ""
          WRITE(LGUNIT(1),*) "*** WARNING: No. blocks < no. nodes ***"
          WRITE(LGUNIT(1),*) ""
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
      if (iprint() .gt. 40)
     .  WRITE(LGUNIT(1),*) 'nblock=',nblock,'np=',np,'times=',times
     .  ,'rem3=',rem3

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
      implicit none
      integer nloop, np, index(0:np-1,0:1)
      integer n, p, m, i, j, sum1, sum2, iprint, lgunit
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
      if (mod(n,2) .ne. 0) then
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
      if (iprint() .gt. 40) then
        write (*,1)
        do i = 0, np-1
          write (lgunit(1),2) i,(index(i,j),j=0,index(i,0))
        enddo
      endif
      if (sum1 .ne. sum2) then
        if (iprint() .gt. 0)
     .    call awrit2('Bug in pdstlb: sum1=%i sum2=%i',' ',128,
     .                lgunit(1),sum1,sum2)
        call rx0(' ')
      endif
    1 format (/' PDSTLB:'/' proc   total    loop number')
    2 format (i5,3x,i5,3x,256i5)
      end
C      program test
C      implicit none
C      integer nloop,np,index(10000000),i,j
C      write(*, 10)
C   10 format ('Enter nloop,np '$)
C      read (*,*) nloop, np
C      write (*,20) nloop, np
C   20 format ('nloop=',i5,' np=',i5)
C      call pdstlb(nloop,np,index)
C      end
#endif 

