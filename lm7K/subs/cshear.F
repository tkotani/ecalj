      SUBROUTINE CSHEAR(B0,B,IBTR)
C- Tries to make microcell more compact by shearing.
C ----------------------------------------------------------------
Ci Inputs
Ci
Co Outputs
Co
Cr Remarks
Cr   Ibtr gives the transformation from bo to b.
C ----------------------------------------------------------------
C      implicit none
C Passed parameters
C Local parameters
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION B(3,3),IBTR(3,3),B0(3,3)
c     logical print
      DO  8  I = 1, 3
        DO  7  J = 1, 3
          B(J,I) = B0(J,I)
          IBTR(J,I) = 0
    7   continue
        IBTR(I,I) = 1
    8 continue

c THE SHEARING TRICK DOES NOT WORK ...
c      DO  3  I = 1, 3
c        J = MOD(I,3)+1
c        K = MOD(I+1,3) + 1
c        JSHORT = 0
c        KSHORT = 0
c        BSHORT = B(1,I)**2 + B(2,I)**2 + B(3,I)**2 - 1.D-6
c        DO  4  JTRY = -5, 5
c          DO  4  KTRY = -5, 5
c          BB1 = B(1,I) + JTRY*B0(1,J) + KTRY*B0(1,K)
c          BB2 = B(2,I) + JTRY*B0(2,J) + KTRY*B0(2,K)
c          BB3 = B(3,I) + JTRY*B0(3,J) + KTRY*B0(3,K)
c          XX = BB1**2 + BB2**2 + BB3**2
c          IF (XX .LT. BSHORT) THEN
c            JSHORT = JTRY
c            KSHORT = KTRY
c            BSHORT = XX
c          ENDIF
c    4   CONTINUE
c        B(1,I) = B(1,I) + JSHORT*B0(1,J) + KSHORT*B0(1,K)
c        B(2,I) = B(2,I) + JSHORT*B0(2,J) + KSHORT*B0(2,K)
c        B(3,I) = B(3,I) + JSHORT*B0(3,J) + KSHORT*B0(3,K)
c        IBTR(J,I) = JSHORT
c        IBTR(K,I) = KSHORT
c    3 continue
c      if (iprint() .ge. 30) then
c        print = .false.
c        do  1  i = 1, 3
c    1   if (ibtr(i,i) .ne. 1) print = .true.
c        isum = 0
C        do  2  i = 1, 3
c          do  2  j = i, 3
c    2   isum = isum + ibtr(i,j)
c        if (isum .ne. 3) print = .true.
c        if (iprint() .gt. 100 .or. print) then
c          print*,'CSHEAR attempting to make microcells more cube-like'
c          print*,'B0 are old reciprocal vectors, B the new; IBTR the ',
c     .         'transformation'
c          WRITE(*,450)
c          DO  45  I = 1, 3
c   45     WRITE(*,451) (B0(J,I),J=1,3),(B(J,I),J=1,3),(IBTR(J,I),J=1,3)
c  451     FORMAT(3F9.5,3X,3F9.5,3X,3I4)
c  450     FORMAT(/' CSHEAR:'/14X,'B0',28X,'B',20X,'IBTR')
c        endif
c      endif
      END

