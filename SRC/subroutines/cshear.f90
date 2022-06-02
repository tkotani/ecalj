SUBROUTINE CSHEAR(B0,B,IBTR)
  !- Tries to make microcell more compact by shearing.
  ! ----------------------------------------------------------------
  !i Inputs
  !i
  !o Outputs
  !o
  !r Remarks
  !r   Ibtr gives the transformation from bo to b.
  ! ----------------------------------------------------------------
  !      implicit none
  ! Passed parameters
  ! Local parameters
  IMPLICIT double precision (A-H,O-Z)
  integer:: ibtr,i,j
  DIMENSION B(3,3),IBTR(3,3),B0(3,3)
  !     logical print
  DO  8  I = 1, 3
     DO  7  J = 1, 3
        B(J,I) = B0(J,I)
        IBTR(J,I) = 0
7    enddo
     IBTR(I,I) = 1
8 enddo

  ! THE SHEARING TRICK DOES NOT WORK ...
  !      DO  3  I = 1, 3
  !        J = MOD(I,3)+1
  !        K = MOD(I+1,3) + 1
  !        JSHORT = 0
  !        KSHORT = 0
  !        BSHORT = B(1,I)**2 + B(2,I)**2 + B(3,I)**2 - 1.D-6
  !        DO  4  JTRY = -5, 5
  !          DO  4  KTRY = -5, 5
  !          BB1 = B(1,I) + JTRY*B0(1,J) + KTRY*B0(1,K)
  !          BB2 = B(2,I) + JTRY*B0(2,J) + KTRY*B0(2,K)
  !          BB3 = B(3,I) + JTRY*B0(3,J) + KTRY*B0(3,K)
  !          XX = BB1**2 + BB2**2 + BB3**2
  !          IF (XX .LT. BSHORT) THEN
  !            JSHORT = JTRY
  !            KSHORT = KTRY
  !            BSHORT = XX
  !          ENDIF
  !    4   CONTINUE
  !        B(1,I) = B(1,I) + JSHORT*B0(1,J) + KSHORT*B0(1,K)
  !        B(2,I) = B(2,I) + JSHORT*B0(2,J) + KSHORT*B0(2,K)
  !        B(3,I) = B(3,I) + JSHORT*B0(3,J) + KSHORT*B0(3,K)
  !        IBTR(J,I) = JSHORT
  !        IBTR(K,I) = KSHORT
  !    3 continue
  !      if (iprint() .ge. 30) then
  !        print = .false.
  !        do  1  i = 1, 3
  !    1   if (ibtr(i,i) .ne. 1) print = .true.
  !        isum = 0
  !        do  2  i = 1, 3
  !          do  2  j = i, 3
  !    2   isum = isum + ibtr(i,j)
  !        if (isum .ne. 3) print = .true.
  !        if (iprint() .gt. 100 .or. print) then
  !          print*,'CSHEAR attempting to make microcells more cube-like'
  !          print*,'B0 are old reciprocal vectors, B the new; IBTR the ',
  !     .         'transformation'
  !          WRITE(*,450)
  !          DO  45  I = 1, 3
  !   45     WRITE(*,451) (B0(J,I),J=1,3),(B(J,I),J=1,3),(IBTR(J,I),J=1,3)
  !  451     FORMAT(3F9.5,3X,3F9.5,3X,3I4)
  !  450     FORMAT(/' CSHEAR:'/14X,'B0',28X,'B',20X,'IBTR')
  !        endif
  !      endif
END SUBROUTINE CSHEAR

