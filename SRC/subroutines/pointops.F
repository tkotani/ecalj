C --- The next few routines are only for the symmetry group parser ---
      SUBROUTINE TSYMOP(T,NT,G,NG)
      CHARACTER*1 T(0:NT-1)
      double precision G(9,1),H(9),HH(9)

      NG = 0
      I = 1
   90 CALL SKIPBL(T,NT,I)

      IF (I .GE. NT) RETURN
      NG = NG+1
      CALL PARSOP(T,I,G(1,NG))
      IF (T(I) .EQ. '*') THEN
        I = I+1
        CALL PARSOP(T,I,H)
        CALL DMPY(G(1,NG),3,1,H,3,1,HH,3,1,3,3,3)
        CALL DVCPY(HH,1,G(1,NG),1,9)
      ENDIF
      GOTO 90
      END
      SUBROUTINE PARSOP(T,I,A)
C- Parse a string for an operator
      double precision V(3),SP,C,S,PI2,A(3,3),ddot
      CHARACTER*1 T(0:1)
      logical parsvc
      PI2 = 8*DATAN(1.D0)
      IF (T(I) .EQ. 'R') THEN
        I = I+1
        READ(T(I),'(I1)',ERR=99) NROT
        I = I+1
        IF (.NOT. PARSVC(T,I,V)) GOTO 99
        SP = DDOT(3,V,1,V,1)
        SP = 1.D0/DSQRT(SP)
        DO  14  K = 1, 3
          V(K) = V(K)*SP
   14   CONTINUE
        C = DCOS(PI2/NROT)
        S = DSIN(PI2/NROT)
        DO  16  K = 1, 3
          DO  15  J = 1, 3
            A(J,K) = (1-C)*V(J)*V(K)
   15     CONTINUE
          A(K,K) = A(K,K) + C
   16   CONTINUE
        A(1,2) = A(1,2) + S*V(3)
        A(3,1) = A(3,1) + S*V(2)
        A(2,3) = A(2,3) + S*V(1)
        A(2,1) = A(2,1) - S*V(3)
        A(1,3) = A(1,3) - S*V(2)
        A(3,2) = A(3,2) - S*V(1)
      ELSE IF (T(I) .EQ. 'M') THEN
        I = I+1
        IF (.NOT. PARSVC(T,I,V)) GOTO 99
        SP = DDOT(3,V,1,V,1)
        DO  11  J = 1, 3
          DO  12  K = 1, 3
            A(K,J) = -2.D0*V(K)*V(J)/SP
   12     CONTINUE
          A(J,J) = A(J,J) + 1.D0
   11   CONTINUE
      ELSE IF (T(I) .EQ. 'I') THEN
        I = I+1
        CALL DVCPY(0.D0,0,A,1,9)
        CALL DVCPY(-1.D0,0,A,4,3)
      ELSE
        GOTO 99
      ENDIF
      RETURN
   99 PRINT *, 'SYMGRP: parse error at ',(T(III),III = 0,I),'  ...'
      STOP
      END
      LOGICAL FUNCTION PARSVC(T,I,V)
C- Generates the vector for a symmetry operation
      double precision V(3)
      CHARACTER*1 T(0:1)
      character*1 chr(3)
      LOGICAL CNVT
      data chr/',',',',')'/
      CALL DINIT(V,3)
      PARSVC = .TRUE.
      IF (T(I) .EQ. 'D') THEN
        V(1) = 1.D0
        V(2) = 1.D0
        V(3) = 1.D0
      ELSE IF (T(I) .EQ. '(') THEN
        I = I+1
        DO  10  M = 1, 3
          IF (.NOT. CNVT(T,V,V,V,V,4,M-1,CHR(M),I)) GOTO 99
   10   CONTINUE
      ELSE
        M = 0
        IF (T(I) .EQ. 'X') M = 1
        IF (T(I) .EQ. 'Y') M = 2
        IF (T(I) .EQ. 'Z') M = 3
        IF (M .EQ. 0) GOTO 99
        V(M) = 1.D0
      ENDIF
   90 I = I+1
      RETURN
   99 PARSVC = .FALSE.
      END

      SUBROUTINE GRPGEN(GEN,NGEN,SYMOPS,NG,NGMX)
C- Generate all symmetry operations from the generation group
C ----------------------------------------------------------------
Ci Inputs
Ci   GEN,NGEN,NGMX
Co Outputs
Co   SYMOPS,NG
Cr Remarks
Cr   This works for point groups only and is set up for integer
Cr   generators.
C ----------------------------------------------------------------
C     implicit none
      integer NGEN,NG,NGMX
      double precision GEN(9,NGEN),SYMOPS(9,NGMX),
     .                 H(9),HH(9),E(9),SIG(9)
      integer IGEN,IG,ITRY,IORD,NNOW,J,IP,I,K,N2,M1,M2,N,M
      integer iprint
      LOGICAL GRPEQL
      DATA E/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/

      if (iprint() .ge. 30) print 679
  679 FORMAT(/)

      CALL GRPCOP(E,SYMOPS)
      NG = 1
      DO  80  IGEN = 1, NGEN
        CALL GRPCOP(GEN(1,IGEN),SIG)
C --- EXTEND THE GROUP BY ALL PRODUCTS WITH SIG ---
        DO  9  IG = 1, NG
          IF (GRPEQL(SYMOPS(1,IG),SIG)) PRINT 650, IGEN,IG
  650     FORMAT(' GENERATOR',I3,'  ALREADY IN GROUP AS ELEMENT',I3)
          IF (GRPEQL(SYMOPS(1,IG),SIG)) GOTO 80
    9   CONTINUE
C --- DETERMINE ORDER ---
        CALL GRPCOP(SIG,H)
        DO  1  ITRY = 1, 100
          IORD = ITRY
          IF (GRPEQL(H,E)) GOTO 2
          CALL GRPPRD(SIG,H,H)
    1   CONTINUE
C --- PRODUCTS OF TYPE:  G1 SIG**P G2 ---
    2   NNOW = NG
        if (iprint() .ge. 40) PRINT 677, IGEN,IORD
  677   FORMAT(' ORDER OF GENERATOR',I3,'  IS',I3)
        DO  8  J = 1, NG
          CALL GRPCOP(SYMOPS(1,J),H)
          DO  10  IP = 1, IORD-1
C H = SIG**IP
            CALL GRPPRD(SIG,H,H)
            DO  11  I = 1, NG
C HH = SYMOPS_I SIG**IP
              CALL GRPPRD(SYMOPS(1,I),H,HH)
              DO  12  K = 1, NNOW
                IF ( GRPEQL(SYMOPS(1,K),HH) ) GOTO 11
   12         CONTINUE
   13         NNOW = NNOW+1
              IF (NNOW .GT. NGMX) GOTO 99
              CALL GRPCOP(HH,SYMOPS(1,NNOW))
   11       CONTINUE
   10     CONTINUE
          IF (J .EQ. 1) N2 = NNOW
    8   CONTINUE
C --- PRODUCTS WITH MORE THAN ONE SANDWICHED SIGMA-FACTOR ---
        M1 = NG+1
        M2 = NNOW
        DO  20  I = 2, 50
          DO  21  N = NG+1, N2
          DO  21  M = M1, M2
            CALL GRPPRD(SYMOPS(1,N),SYMOPS(1,M),H)
            DO  22  K = 1, NNOW
              IF (GRPEQL(SYMOPS(1,K),H)) GOTO 21
   22       CONTINUE
            NNOW = NNOW+1
            IF (NNOW .GT. NGMX) GOTO 99
            CALL GRPCOP(H,SYMOPS(1,NNOW))
   21     CONTINUE
          IF (M2 .EQ. NNOW) GOTO 25
          M1 = M2 + 1
          M2 = NNOW
   20   CONTINUE
   25   CONTINUE
        NG = NNOW
   80 CONTINUE
      if (iprint() .ge. 30) then
        print 678, NG,NGEN
  678   format(' GRPGEN:',i3,' SYMMETRY OPERATIONS FROM',I2,
     .         ' GENERATORS')
      endif
      RETURN
   99 STOP 'GRPGEN: too many elements'
      END
      SUBROUTINE GRPCOP(G,H)
C- Copy matrix
      double precision H(9),G(9)
      DO  10  I = 1, 9
        H(I) = G(I)
   10 CONTINUE
      RETURN
      END
      LOGICAL FUNCTION GRPEQL(G1,G2)
C- Checks if G1 is equal to G2
c      implicit none
      double precision g1(9),g2(9),dabs,d1mach,tiny
      intrinsic dabs
      external d1mach
      integer i

      grpeql = .false.
      tiny = 100 * d1mach(3)
      do  10  i = 1, 9
ca !!!!!!!!!!!
c       if (dabs(g1(i) - g2(i)) .gt. tiny) return
        if (dabs(g1(i) - g2(i)) .gt. 0.001) return
   10 continue
      grpeql = .true.
      return
      end
      subroutine grpprd(g1,g2,g1xg2)
C- returns the product of two point group operations
c      implicit none
      double precision g1(9),g2(9),g1xg2(9),h(9)
      integer i

      call dmpy(g1,3,1,g2,3,1,h,3,1,3,3,3)
      do  10  i = 1, 9
        g1xg2(i) = h(i)
   10 continue
      return
      end
