      subroutine fbz2(q,ginv, index_qbz,n, qbz, nqbz,
     o q_out, i_out)
      implicit none
      integer(4):: n,index_qbz(n,n,n),intq(3),i_out,nqbz
      real(8):: q(3),qbz(3,nqbz),q_out(3),ginv(3,3),qx(3)
c Fro given q, find q_out in the first BZ (FBZ)
c      print *,'fbz2:xxx1', matmul(ginv,q)
      call rangedq(matmul(ginv,q), qx)
      intq =  qx*n + 1
      i_out = index_qbz(intq(1),intq(2),intq(3))
      q_out = qbz(:,i_out)
cccccccccccccccccccccc
      call rangedq(matmul(ginv,q_out-q), qx)
      if(abs(sum(qx**2))>1d-8) then
        print *,' q  =  '
        print *,  q
        print *,' ginv*q  =  '
        print *,  matmul(ginv,q)
        print *,' intq  =  '
        print *,  intq
        print *,' qx  =  '
        print *,  qx
        print *,' q_out i_out= '
        print *,  q_out, i_out
        print *,'fbz2: qx**2>1d-8 something wrong!!!'
        stop 'fbz2: qx**2>1d-8 something wrong!!!'
      endif
      end

      SUBROUTINE TETIRR(QB,N1,N2,N3,IPQ,NQ,NTET,IDTET)
C-  Finds inequivalent tetrahedra and counts them
C ----------------------------------------------------------------------
Ci Inputs:
Ci  qb,n1,n2,n3,ipq, output from BZMESH;
Ci  nq, no. of irreducible k-points;
Co Outputs:
Co  ntet, No. of different tetrahedra
Co  idtet(1-4,i), Identifies the i'th tetrahedron in terms of the four
Co  idtet(0,i), no. of tetrahedra of the i'th kind
Cm Memory:
Cm  No large internal storage; heap not accessed.
C ----------------------------------------------------------------------

      IMPLICIT double precision (A-H,O-Z)
      DIMENSION QB(3,3),IPQ(N1,N2,N3),QB1(3,3),
     .          IBTR(3,3),KCUT(3,4,6),IMC(0:1,0:1,0:1),
     .          idtet(0:4,*),iq(4)

      CALL CCUTUP(QB,QB1,IBTR,KCUT)
      ntet = 0
C ----- START LOOPING OVER MICROCELLS ---------
      DO  20  I3 = 1, N3
      DO  20  I2 = 1, N2
      DO  20  I1 = 1, N1
C ----- SET UP IDENTIFIERS AT 8 CORNERS OF MICROCELL ------
        DO  1  K1 = 0, 1
          J1 = MOD(I1 + K1 -1,N1) + 1
        DO  1  K2 = 0, 1
          J2 = MOD(I2 + K2 -1,N2) + 1
        DO  1  K3 = 0, 1
          J3 = MOD(I3 + K3 -1,N3) + 1
          IMC(K1,K2,K3)  = IPQ(J1,J2,J3)
   1    continue
C ----- LOOP OVER TETRAHEDRA --------------
        DO  10  ITET = 1, 6
          DO  2  IC = 1, 4
            K1 = KCUT(1,IC,ITET)
            K2 = KCUT(2,IC,ITET)
            K3 = KCUT(3,IC,ITET)
            IQ(IC) = IMC(K1,K2,K3)
    2     continue

C ----- order the identifiers -----
          DO  3  J=1,3
          DO  3  I=1,4-J
            IF(IQ(I).GT.IQ(I+1))THEN
              II=IQ(I)
              IQ(I)=IQ(I+1)
              IQ(I+1)=II
            ENDIF
    3     CONTINUE
          if (ntet .eq. 0) goto 5
          do  4  n = 1, ntet
            if (       idtet(1,n) .eq. iq(1)
     .           .and. idtet(2,n) .eq. iq(2)
     .           .and. idtet(3,n) .eq. iq(3)
     .           .and. idtet(4,n) .eq. iq(4) ) then
              idtet(0,n) = idtet(0,n) + 1
              goto 10
            endif
    4     continue
    5     continue
          ntet=ntet+1
          do  6  i = 1, 4
            idtet(i,ntet) = iq(i)
    6     continue
          idtet(0,ntet)=1
   10   CONTINUE
   20 CONTINUE
      if ( iprint() .ge. 40 ) print 922, ntet, 6*n1*n2*n3
  922 format(1x,'TETIRR: ',i6,' inequivalent tetrahedra from ',i8)
      END

c$$$      SUBROUTINE BZINTS(N1,N2,N3,EP,WP,NQ,nband,NB,NSP,EMIN,EMAX,
c$$$c     .   DOS,NR,EF,JOB,NTET,IDTET)
c$$$C-  Bz integrations by linear method.
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs:
c$$$Ci  ep, energy bands;
c$$$Ci  nq, no. of irreducible k-points; nb, no. of bands; nsp, see BNDASA;
c$$$Ci  emin, emax, dos, nr : for IDOS, energy window,
c$$$Ci  emin, emax (JOB=2) tolerance in efermi
c$$$Ci  IDOS, number of points; ef, Fermi energy (job = 2); job, switch :
c$$$Ci  JOB=1: MAKES IDOS.  JOB=2:  MAKES BLOECHL-WEIGHTS.
c$$$Ci  ntet, No. of different tetrahedra
c$$$Ci  idtet(1-4,i), Identifies the i'th tetrahedron in terms of the four
c$$$Ci  irreducible k-points:
c$$$Ci  idtet(0,i), no. of tetrahedra of the i'th kind
c$$$Co Outputs:
c$$$Co  dos, Integrated Density of States (IDOS) (job = 1)
c$$$Co  wp, Bloechl quadrature weights (job = 2)
c$$$Cm Memory:
c$$$Cm  No large internal storage; heap not accessed.
c$$$C ----------------------------------------------------------------------
c$$$C      implicit none
c$$$C Passed parameters
c$$$C Local parameters
c$$$c      IMPLICIT double precision (A-H,O-Z)
c$$$c     integer fopn
c$$$c      DIMENSION EP(nb,nsp,nq),DOS(NR),EC(4),WC(4,2),WP(nband,nsp,nq),
c$$$c     .   idtet(0:4,*)
c$$$
c$$$      IF (JOB .NE. 1 .AND. JOB .NE. 2) STOP '*** BAD JOB IN BZINTS'
c$$$      IF (JOB .EQ. 1) call dinit(dos,2*nr)
c$$$      IF (JOB .EQ. 2) call dinit(wp,nband*nsp*nq)
c$$$      SEV1 = 0.D0
c$$$      SEV2 = 0.D0
c$$$      volwgt = (3.d0 - nsp) / (n1*n2*n3*6.d0)
c$$$      do  40  isp = 1, nsp
c$$$C ----- START LOOPING OVER TETRAHEDRA ---------
c$$$        DO  20  ITET = 1, NTET
c$$$         iq1=idtet(1,itet)
c$$$         iq2=idtet(2,itet)
c$$$         iq3=idtet(3,itet)
c$$$         iq4=idtet(4,itet)
c$$$         DO  20  IB = 1, nband
c$$$C ----- SET UP ENERGIES AT 4 CORNERS OF TETRAHEDRA ------
c$$$           ec(1) = ep(ib,isp,iq1)
c$$$           ec(2) = ep(ib,isp,iq2)
c$$$           ec(3) = ep(ib,isp,iq3)
c$$$           ec(4) = ep(ib,isp,iq4)
c$$$           etop = dmax1(ec(1),ec(2),ec(3),ec(4))
c$$$           ebot = dmin1(ec(1),ec(2),ec(3),ec(4))
c$$$           IF (JOB .EQ. 1) THEN
c$$$             if ( ebot .lt. emax )
c$$$     .         CALL SLINZ(volwgt*idtet(0,ITET),EC,EMIN,EMAX,DOS,NR)
c$$$           ELSE
c$$$             if ( ef .ge. ebot ) then
c$$$               CALL FSWGTS(volwgt*idtet(0,ITET),EC,EF,ETOP,WC)
c$$$               SEV1 = SEV1 + WC(1,1)*EC(1) + WC(2,1)*EC(2) +
c$$$     .                       WC(3,1)*EC(3) + WC(4,1)*EC(4)
c$$$               SEV2 = SEV2 + WC(1,2)*EC(1) + WC(2,2)*EC(2) +
c$$$     .                       WC(3,2)*EC(3) + WC(4,2)*EC(4)
c$$$               WP(ib,isp,iq1) = WP(ib,isp,iq1) + WC(1,1) + WC(1,2)
c$$$               WP(ib,isp,iq2) = WP(ib,isp,iq2) + WC(2,1) + WC(2,2)
c$$$               WP(ib,isp,iq3) = WP(ib,isp,iq3) + WC(3,1) + WC(3,2)
c$$$               WP(ib,isp,iq4) = WP(ib,isp,iq4) + WC(4,1) + WC(4,2)
c$$$             endif
c$$$           ENDIF
c$$$   20    CONTINUE
c$$$   40  continue
c$$$
c$$$      if ( iprint() .ge. 30 ) then
c$$$        IF (JOB .EQ. 2) THEN
c$$$          SUMWP = 0.D0
c$$$          do  56  isp = 1, nsp
c$$$            DO  56  IB = 1, nband
c$$$            DO  56  IQ = 1, NQ
c$$$  56      SUMWP = SUMWP + WP(ib,isp,iq)
c$$$          if (iprint() .ge. 30) print 922, ef, sumwp, sev1+sev2, sev2
c$$$c         write(fopn('LOG'),922) ef, sumwp, sev1+sev2, sev2
c$$$        ENDIF
c$$$      endif
c$$$
c$$$  922 format(1x,'BZINTS: Fermi energy:',f10.6,';',F10.6,' electrons'/
c$$$     .       9x,'Band energy:',f11.6,
c$$$     .       ', including Bloechl correction:',f10.6)
c$$$      END
c$$$
c$$$      SUBROUTINE SLINZ(VOLWGT,EC,EMIN,EMAX,DOSI,NR)
c$$$C- Adds to number-of-states for one tetrahedron
c$$$C ----------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   volwgt, weight on tetrahedron; ec energies at corners of tethdn.;
c$$$Ci   emin, emax, energy window; nr, number of bins + 1
c$$$Co Outputs
c$$$Co   dosi(k), integrated density in kth bin from tethdn.
c$$$Cr Remarks
c$$$Cr
c$$$C ----------------------------------------------------------------
c$$$C      implicit none
c$$$C Passed parameters
c$$$C Local parameters
c$$$      IMPLICIT double precision (A-H,O-Z)
c$$$      DIMENSION EC(4),DOSI(NR)
c$$$      DO  3  I = 1, 3
c$$$        DO  3  J = 1, 4-I
c$$$        IF (EC(J) .GT. EC(J+1)) THEN
c$$$        E=EC(J)
c$$$        EC(J) = EC(J+1)
c$$$        EC(J+1) = E
c$$$        ENDIF
c$$$    3 continue
c$$$      E1 = EC(1)
c$$$      E2 = EC(2)
c$$$      E3 = EC(3)
c$$$      E4 = EC(4)
c$$$      if (e4 .lt. emin) then
c$$$        i4=1
c$$$        go to 26
c$$$      endif
c$$$      DE = (EMAX-EMIN)/(NR-1)
c$$$      I01   = (E1   -EMIN)/DE + 1.9999999D0
c$$$      I02   = (E2   -EMIN)/DE + 1.9999999D0
c$$$      I03   = (E3   -EMIN)/DE + 1.9999999D0
c$$$      I04   = (E4   -EMIN)/DE + 1.9999999D0
c$$$C --------------------------------
c$$$      I1 = MAX0(I01  ,1)
c$$$      I2 = MIN0(I02  -1,NR)
c$$$      IF (I1 .LE. I2) THEN
c$$$        CC = VOLWGT/((E2-E1)*(E3-E1)*(E4-E1))
c$$$        DO  20  I = I1, I2
c$$$          X = EMIN - E1 + (I-1)*DE
c$$$          DOSI(I) = DOSI(I) + CC*X**3
c$$$   20   continue
c$$$      ENDIF
c$$$      I2 = MAX0(I02  ,1)
c$$$      I3 = MIN0(I03  -1,NR)
c$$$      IF (I2 .LE. I3) THEN
c$$$        C3 = VOLWGT*(E1+E2-E3-E4)/((E3-E1)*(E4-E1)*(E3-E2)*(E4-E2))
c$$$        C2 = VOLWGT*3.D0/((E3-E1)*(E4-E1))
c$$$        C1 = C2*(E2-E1)
c$$$        C0 = C1*(E2-E1)/3.D0
c$$$        DO  21  I = I2, I3
c$$$          X = EMIN - E2 + (I-1)*DE
c$$$          DOSI(I) = DOSI(I) + C0 + X*(C1 + X*(C2 + X*C3))
c$$$   21   continue
c$$$      ENDIF
c$$$      I3 = MAX0(I03  ,1)
c$$$      I4 = MIN0(I04  -1,NR)
c$$$      IF (I3 .LE. I4) THEN
c$$$        CC = VOLWGT/((E3-E4)*(E2-E4)*(E1-E4))
c$$$        DO  22  I = I3, I4
c$$$          X = EMIN - E4 + (I-1)*DE
c$$$          DOSI(I) = DOSI(I) + VOLWGT - CC*X**3
c$$$   22   continue
c$$$      ENDIF
c$$$      I4 = MAX0(I04  ,1)
c$$$   26 continue
c$$$      DO  25  I = I4, NR
c$$$   25 DOSI(I) = DOSI(I) + VOLWGT
c$$$      END

      SUBROUTINE FSWGTS(VOLWGT,E,EF,ETOP,W)
C- Makes weights for integration up to ef for one tetrahedron.
C ----------------------------------------------------------------
Ci Inputs
Ci   ,EF?
Co Outputs
Co
Cr Remarks
Cr   W(I,1): NORMAL WEIGHTS.  W(I,2): BLOECHL-CORRECTION.
C ----------------------------------------------------------------

      IMPLICIT double precision (A-H,O-Z)
      DIMENSION E(4),EX(4),W1(4),W2(4),EC(4),W(4,2),ISORT(4)
      VW4 = VOLWGT/4.D0
      DO  22  I = 1, 4
        W(I,1) = 0.D0
        W(I,2) = 0.D0
   22 continue
      IF (EF .GE. ETOP) THEN
        DO  23  I = 1, 4
          W(I,1) = VW4
   23   continue
        RETURN
      ENDIF
C ----- SORT ENERGIES INTO EC -----------------
      DO  2  I = 1, 4
        EX(I) = E(I)
    2 continue
      DO  3  I = 1, 4
        I00 = 1
        DO  4  J = 2, 4
          IF (EX(J) .LT. EX(I00)) I00 = J
    4   continue
        EC(I) = EX(I00)
        ISORT(I) = I00
        EX(I00) = ETOP + 1.D0
    3 continue
      E1 = EC(1)
      E2 = EC(2)
      E3 = EC(3)
      E4 = EC(4)
C ------- CASE EF BETWEEN E2,E3 -------------
      IF (E2 .LT. EF .AND. EF .LE. E3) THEN
        A31 = (EF-E1)/(E3-E1)
        A41 = (EF-E1)/(E4-E1)
        A32 = (EF-E2)/(E3-E2)
        A42 = (EF-E2)/(E4-E2)
        V1 = A31*A41
        V2 = A31*A42*(1.D0-A41)
        V3 = A42*A32*(1.D0-A31)
        W1(1) = (V1*(3-A31-A41) + V2*(2-A31-A41) + V3*(1-A31))*VW4
        W1(2) = (V1 + V2*(2-A42) + V3*(3-A32-A42))*VW4
        W1(3) = (V1*A31 + V2*A31 + V3*(A31+A32))*VW4
        W1(4) = (V1*A41 + V2*(A41+A42) + V3*A42)*VW4
        DF = ((E1+E2-E3-E4)*A32*A42 + 2*EF - E1 - E2)/((E3-E1)*(E4-E1))
        DF = 3.D0*VOLWGT*DF
C ------- CASE EF BETWEEN E1,E2 -------------------
      ELSE IF (E1 .LT. EF .AND. EF .LE. E2) THEN
        A21 = (EF-E1)/(E2-E1)
        A31 = (EF-E1)/(E3-E1)
        A41 = (EF-E1)/(E4-E1)
        XXX = A21*A31*A41*VW4
        W1(1) = XXX*(4.D0-A21-A31-A41)
        W1(2) = XXX*A21
        W1(3) = XXX*A31
        W1(4) = XXX*A41
        DF = 3.D0*VOLWGT*A31*A41/(E2-E1)
C ------- CASE EF BETWEEN E3,E4 -------------------
      ELSE IF (E3 .LT. EF .AND. EF .LE. E4) THEN
        A14 = (EF-E4)/(E1-E4)
        A24 = (EF-E4)/(E2-E4)
        A34 = (EF-E4)/(E3-E4)
        XXX = A14*A24*A34*VW4
        W1(1) = VW4 - XXX*A14
        W1(2) = VW4 - XXX*A24
        W1(3) = VW4 - XXX*A34
        W1(4) = VW4 - XXX*(4.D0-A14-A24-A34)
        DF = -3.D0*VOLWGT*A14*A24/(E3-E4)
      ENDIF
C ------ HERE FOR BLOECHL CORRECTION ---------
      DO  40  M = 1, 4
        W2(M) = 0.D0
      DO  40  N = 1, 4
        W2(M) = W2(M) + (EC(N)-EC(M))*DF*0.025D0
   40 continue
C ----------------------------------------------
      DO  35  I = 1, 4
        J = ISORT(I)
        W(J,1) = W1(I)
        W(J,2) = W2(I)
   35 continue
      RETURN
      END

      SUBROUTINE CCUTUP(B0,B,IBTR,KCUT)
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION B(3,3),IBTR(3,3),B0(3,3),KCUT0(3,4,6),KCUT(3,4,6),
     .   SHIFT(3),P(3,4)
      DATA KCUT0/
     .  0,0,0, 0,1,0, 1,1,0, 1,1,1,  0,0,0, 1,0,0, 1,1,0, 1,1,1,
     .  0,0,0, 1,0,0, 1,0,1, 1,1,1,  0,0,0, 0,1,0, 0,1,1, 1,1,1,
     .  0,0,0, 0,0,1, 0,1,1, 1,1,1,  0,0,0, 0,0,1, 1,0,1, 1,1,1 /
      DATA SHIFT/0.D0,0.D0,0.D0/
      ANRM2(X,Y,Z)=X*X*1.00001D0+Y*Y*1.00002D0+Z*Z*1.00003D0
     .  -X*0.000004D0-Y*0.000003D0-Z*0.000002D0
C ------ CALL CSHEAR TO GET MOST COMPACT CELL (doesn't work) --------
      CALL CSHEAR(B0,B,IBTR)
C ----- CHOSE CUTUP WITH SHORTEST MAX EDGE ---------
      if (iprint() .gt. 100) WRITE(*,*) 'CUTUP : '
      LZ = 0
      LXX = 0
      LYY = 0
      EDGMAX = 1.D20
      EDGMIN = 0.D0
      DO  10  LX = 0, 1
      DO  10  LY = 0, 1
        DO  12  ITET = 1, 6
        DO  12  IC = 1, 4
          CALL MXMYMZ(KCUT0(1,IC,ITET),KCUT(1,IC,ITET),LX,LY,LZ)
   12   continue
        EDMIN = 1D20
        EDMAX = 0D0
        DO  20  ITET = 1, 6
          DO  21  IC = 1, 4
            CALL GTBVEC(KCUT(1,IC,ITET),B,SHIFT,P(1,IC))
   21     continue
          DO  13  I = 1, 3
          DO  13  J = I+1, 4
            XX = ANRM2(P(1,I)-P(1,J),P(2,I)-P(2,J),P(3,I)-P(3,J))
            EDMAX = DMAX1(EDMAX,XX)
            EDMIN = DMIN1(EDMIN,XX)
   13     continue
   20   CONTINUE
        if (iprint() .gt. 100)
     .    WRITE(*,706) LX,LY,DSQRT(EDMIN),DSQRT(EDMAX)
  706   FORMAT(' LX,LY=',2I5,'   EDMIN=',F10.5,'   EDMAX=',F10.5)
        IF (EDMAX .LT. EDGMAX) THEN
          LXX = LX
          LYY = LY
          EDGMAX = EDMAX
          EDGMIN = EDMIN
        ENDIF
   10 CONTINUE
      DO  22  ITET = 1, 6
      DO  22  IC = 1, 4
        CALL MXMYMZ(KCUT0(1,IC,ITET),KCUT(1,IC,ITET),LXX,LYY,LZ)
   22 continue
      if (iprint() .gt. 100)
     .  WRITE(*,783) LXX,LYY,DSQRT(EDGMIN),DSQRT(EDGMAX)
  783 FORMAT(' LXX=',I1,'   LYY=',I1,'   EDGMIN=',F10.5,
     .  '   EDGMAX=',F10.5)

      END
      SUBROUTINE GTBVEC(K,B,SHIFT,V)
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION K(3),V(3),B(3,3),SHIFT(3)
      V(1) = SHIFT(1) + K(1)*B(1,1) + K(2)*B(1,2) + K(3)*B(1,3)
      V(2) = SHIFT(2) + K(1)*B(2,1) + K(2)*B(2,2) + K(3)*B(2,3)
      V(3) = SHIFT(3) + K(1)*B(3,1) + K(2)*B(3,2) + K(3)*B(3,3)
      RETURN
      END

      SUBROUTINE MXMYMZ(KIN,K,LX,LY,LZ)
C- Do mirrors in x,y,z if lx,ly,lz=1, respectively
C ----------------------------------------------------------------
Ci Inputs
Ci
Co Outputs
Co
Cr Remarks
Cr
C ----------------------------------------------------------------
      DIMENSION KIN(3),K(3)
      K(1) = KIN(1)
      K(2) = KIN(2)
      K(3) = KIN(3)
      IF (LX .EQ. 1) K(1) = 1-K(1)
      IF (LY .EQ. 1) K(2) = 1-K(2)
      IF (LZ .EQ. 1) K(3) = 1-K(3)
      RETURN
      END

      SUBROUTINE CSHEAR(B0,B,IBTR)
C- Tries to make microcell more compart by shearing.
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

      SUBROUTINE fermi(QVAL,DOS,NDOS,EMIN,EMAX,EFERM,e1,e2,dosef)
C- Makes fermi energy from integrated density
C ----------------------------------------------------------------------
Ci Inputs
Ci   qval, number of electrons to fermi level; dos(i) integrated
Ci   density at bin i; ndos, number of bins + 1; emin, emax, energy
Ci   window.
Co Outputs
Co   Eferm, Fermi energy; e1, e2, confidence limits on Fermi energy
Co   i.e., Fermi energy lies between e1 and e2.
Co   dosef:  density of states at fermi level
Cr Remarks
Cr   emin and e1 (and emax and e2) may point to the same address.
C ----------------------------------------------------------------------

C Passed parameters
      integer ndos
      double precision qval,dos,emin,emax,eferm,e1,e2,dosef
C Local parameters
      integer i1,ie
      double precision de,q,q1,q2,d1mach
C External procedures
      external d1mach

      DIMENSION DOS(NDOS)

ccccccccccccccccccc
c      DE = (EMAX-EMIN)/(NDOS-1)
c      do i1=1,ndos
c        write(6,"(' e dostot =',i4,d13.6,d13.6)")i1, emin + de*(i1 - 1),dos(i1)
c      enddo
c      write(6,"(' qval =',d13.6)") qval
ccccccccccccccccccccc

      if (dos(1) .gt. qval) print *, 'FERMI: EMIN,EMAX=', emin,emax
      if (dos(1) .gt. qval) stop 'FERMI: Fermi energy lies below EMIN'
      if (dos(ndos) .lt. qval) print *, 'FERMI: EMIN,EMAX=', emin,emax
      if (dos(ndos) .lt. qval) then
        stop 'FERMI: Fermi energy lies above EMAX'
      endif


      DE = (EMAX-EMIN)/(NDOS-1)
      I1 = 1
      q = qval + d1mach(3)
      DO  1  IE = 2, NDOS
c oj- changed due to bug in CRAY compiler
        I1 = IE
        IF ( DOS(IE) .gt. q ) goto 2
c       I1 = IE
    1 continue
    2 continue
      i1 = i1 - 1
      Q1 = DOS(I1)
      Q2 = DOS(I1+1)
C ------------------
      e1 = emin + de*(i1 - 1)
      e2 = e1 + de
c     e1 = emin + de*(i1 - 1)  -de
c      e2 = e1 + 3*de
c----------------------
      print *,e1,QVAL,Q1,Q2,DE
      EFERM = e1 + (QVAL-Q1)/(Q2-Q1)*DE
      dosef = (q2-q1)/de
      END


      subroutine fbzws(gu,g,q)

c takes q into the 1st BZ
c input:
c g  = base reciprocal lattice vectors
c gu = normalised base reciprocal lattice vectors,
c      (each  g is a unit vector)
c q  = the q vector
c g1 = projection of q along g1 etc.
c all vectors are in units of 2*pi/a
c output:
c q  = is now inside the 1st BZ

      implicit real*8(a-h,o-z)
      dimension gu(3,3),g(3,3),q(3)

 1111 g1         = q(1)*gu(1,1)+q(2)*gu(2,1)+q(3)*gu(3,1)
      g2         = q(1)*gu(1,2)+q(2)*gu(2,2)+q(3)*gu(3,2)
      g3         = q(1)*gu(1,3)+q(2)*gu(2,3)+q(3)*gu(3,3)
      if(dabs(g1)-1.d-5 .le. .5d0 .and.
     .   dabs(g2)-1.d-5 .le. .5d0 .and.
     .   dabs(g3)-1.d-5 .le. .5d0      ) return
      if(dabs(g1)-1.d-5 .gt. .5d0) then
cmod 93.09.30
c     ig         = g1 + dsign(4.999d0,g1)
        ig         = g1 + dsign(.4999d0,g1)
        q(1)       = q(1) - ig*g(1,1)
        q(2)       = q(2) - ig*g(2,1)
        q(3)       = q(3) - ig*g(3,1)
      endif
      g2         = q(1)*gu(1,2)+q(2)*gu(2,2)+q(3)*gu(3,2)
      if(dabs(g2)-1.d-5 .gt. .5d0) then
        ig         = g2 + dsign(.4999d0,g2)
        q(1)       = q(1) - ig*g(1,2)
        q(2)       = q(2) - ig*g(2,2)
        q(3)       = q(3) - ig*g(3,2)
      endif
      g3         = q(1)*gu(1,3)+q(2)*gu(2,3)+q(3)*gu(3,3)
      if(dabs(g3)-1.d-5 .gt. .5d0) then
        ig         = g3 + dsign(.4999d0,g3)
        q(1)       = q(1) - ig*g(1,2)
        q(2)       = q(2) - ig*g(2,2)
        q(3)       = q(3) - ig*g(3,2)
        goto 1111
      endif

      return
      end
      subroutine gentran (plat,tmax,nshell,
     o t,xt,ntsh,nt )

c generates translational lattice vectors
c may be used for both real and reciprocal space

c plat   = primitive translation lattice vectors
c tmax   = maximum length, in the same unit as plat
c nshell = minimum integer to cover tmax, from notrans.f

c t      = list of translational vectors to neighbouring cells
c          in increasing magnitude, dimension=2*(2*nshell+1)**3
c xt     = lengths, dimension (2*nshell+1)**3
c nt     = no. translational vectors
c ntsh(ishell) = no. translational vectors in a given shell
c                i.e. of the same length, no necessarily
c                related by symmetry operations

      implicit real*8(a-h,o-z)

      dimension plat(3,3),xt(*)
      dimension t(3,1),ntsh(*)

c generate neighbouring clusters
      tmax2      = tmax*tmax
      it         = 0
      do       i = -nshell,nshell
        do       j = -nshell,nshell
          do       k = -nshell,nshell
            t1         = i*plat(1,1) + j*plat(1,2) + k*plat(1,3)
            t2         = i*plat(2,1) + j*plat(2,2) + k*plat(2,3)
            t3         = i*plat(3,1) + j*plat(3,2) + k*plat(3,3)
            xl         = t1*t1 + t2*t2 + t3*t3
            if(xl .le. tmax2)then
              it         = it + 1
              t(1,it)    = t1
              t(2,it)    = t2
              t(3,it)    = t3
              xt(it)     = xl
            end if
          end do
        end do
      end do
      call dcopy(3*it,t(1,1),1,t(1,it+1),1)

c search for the shortest length
      nt         = 0
      ishell     = 0
      j          = 0
 1111 xmin       = 1.d20
      do       i = 1,it
        if(xt(i) .lt. xmin) xmin = xt(i)
      end do

c .search for other vectors of the same length
      its         = 0
      do        i = 1,it
        if(xt(i) - xmin .lt. 1.d-5)then
          its         = its + 1
          j           = j + 1
          t(1,j)      = t(1,it+i)
          t(2,j)      = t(2,it+i)
          t(3,j)      = t(3,it+i)
          xt(i)       = 1.d30
        endif
      end do

c .count the number of shells and store the number of
c  lattice vectors of a given shell i.e. of the same length
      ishell       = ishell + 1
      ntsh(ishell) = its
      nt           = nt + its
      if(j .lt. it)goto 1111

c calculate the lengths
      do         i = 1,it
        xl           = vdv(t(1,i),t(1,i),3)
        xt(i)        = dsqrt(xl)
      end do

      return
      end
c------------------------------------------------------------------
      subroutine invgrp (g,ng,
     o invg)

c index for the inverse of group rotations
c given ig, invg(ig) gives the index of the inverse of g

c g  = matrix representations of point groups
c ng = no. group elements

c invg = given ig, invg gives the index corresponding to
c        the inverse of g

      implicit real*8(a-h,o-z)

      dimension g(3*3,ng)
      dimension invg(ng)
      dimension ginv(3*3),w1(9)
      data tol /1.d-5/

      do       i = 1,ng

c invert g
        call dinv33(g(1,i),0,ginv,det)

c find g = ginv
        do       j = 1,ng
          call vminv   (g(1,j),ginv,9,w1)
          diff       = vdv (w1,w1,9)
          if (diff .lt. tol) goto 1111
        end do
 1111   invg(i)    = j
        if (j .gt. ng) stop 'invgrp: cannot identify inverse of g'
      end do

      return
      end
c---------------------------------------------------------------------
      integer function invrot (irot,invg,ngrp)

c 92.02.26
c given a rotation labelled by irot, invrot gives the label for
c the inverse of the rotation

c irot = label for R
c invg(R) = R^(-1)
c ngrp = no. rotation matrices

      implicit real*8 (a-h,o-z)
      dimension invg(ngrp)

      invrot     = invg(irot)

      return
      end
c------------------------------------------------------------------
      subroutine nkstar (qibz,qbz,grp,ginv,
     w iw1,
     d nqibz,nqbz,ngrp,
     o nstar,irotk)
c 92.02.22
c generates the no. stars of k
c i.e. the no. times k appears in the FBZ
c qibz  = k { IBZ
c qbz   = k { FBZ
c grp   = rotation matrices
c nqibz = no. k { IBZ
c nqbz  = no. k { FBZ
c ngrp  = no. rotation matrices
c nstar(k) = no. times k appears in the FBZ
c irotk(k{IBZ,R) = index to k{FBZ
      implicit real*8 (a-h,o-z)
      dimension qibz(3,nqibz),qbz(3,nqbz),grp(3,3,ngrp),ginv(3,3)
      dimension nstar(nqibz),irotk(nqibz,ngrp)
      integer(4):: verbose
      real(8):: diff(3),diff2(3),tol=1d-6,gg(3,3)
c initialise index
c      call iinit   (irotk,nqibz*ngrp)
c      call iinit   (iw1,nqbz)
c      call iinit   (nstar,nqibz)

      do kp = 1,nqbz
        if(verbose()>104) write(6,"(' === kp=',i3,' qbz=',3f8.3)")kp,qbz(:,kp)
      enddo

      irotk=0
      nstar=0
      if(verbose()>104) print *,' nkstar:'
      do kp = 1,nqbz
        if(verbose()>104) write(6,"(' ===== kp=',i3,' qbz=',3f8.3)")kp,qbz(:,kp)
        do k  = 1,nqibz
          do ir = 1,ngrp
c        if(verbose()>104) print *,' grp=',ir !;      print *, grp(:,ir)
            diff = matmul(grp(:,:,ir),qibz(:,k)) - qbz(:,kp)
            call rangedq(matmul(ginv,diff), diff2)
            if(verbose()>104)
     &      write(6,"(' matmul(ginv,diff)=',3f8.3,' ',3f8.3)") diff, matmul(ginv,diff)
            if(sum(abs(diff2))< tol) then
              irotk(k,ir)= kp
              nstar(k)   = nstar(k) + 1
              goto 1022
            endif
          enddo
        enddo
        stop 'nkstar: can not find irotk'
 1022   continue
      enddo
c check that the sum of stars equal to the no. k{FBZ
      nsum = ivsum (nstar,nqibz)
      if (nsum .ne. nqbz) then
        print *,' nums nqbz=',nsum,nqbz
        do k  = 1,nqibz
          do ir = 1,ngrp
            if(irotk(k,ir)/=0) write(6,"(' k ir irotk=',3i6)") k,ir, irotk(k,ir)
          enddo
        enddo
        stop 'nkstar: wrong no. stars'
      endif
c write k { IBZ and no. stars to file KPNT
      ifkpnt     = ifile('KPNT')
      if (ifkpnt .gt. 0) then
        write (ifkpnt,*) 'irreducible k-points and no. stars'
        write (ifkpnt,*) 'k, k-vector, nstar '
        do       k = 1,nqibz
          write(ifkpnt,"(1x,i5,3f8.5,i3)")k,qibz(1,k),qibz(2,k),qibz(3,k),nstar(k)
        end do
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      stop 'test end'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end








c------------------------------------------------------------------
      subroutine nkstar_old (qibz,qbz,grp,ginv,
     w iw1,
     d nqibz,nqbz,ngrp,
     o nstar,irotk)

c 92.02.22
c generates the no. stars of k
c i.e. the no. times k appears in the FBZ

c qibz  = k { IBZ
c qbz   = k { FBZ
c grp   = rotation matrices
c nqibz = no. k { IBZ
c nqbz  = no. k { FBZ
c ngrp  = no. rotation matrices

c nstar(k) = no. times k appears in the FBZ
c irotk(k{IBZ,R) = index to k{FBZ

      implicit real*8 (a-h,o-z)
      dimension qibz(3,nqibz),qbz(3,nqbz),grp(9,ngrp),ginv(9),
     w          iw1(nqbz)
      dimension nstar(nqibz),irotk(nqibz,ngrp)
      dimension w1(3),w2(3),w3(3)
      data      tol/1.d-5/

c initialise index
      call iinit   (irotk,nqibz*ngrp)
      call iinit   (iw1,nqbz)
      call iinit   (nstar,nqibz)

c loop over  k { IBZ
      do       k = 1,nqibz
ccccccccccccccccccccccccccccccccccccccc
c         print *,' k =',k,' qibz=',qibz(:,k)
ccccccccccccccccccccccccccccccccccccccc

c loop over rotations
        do      ir = 1,ngrp

cccccccccccccccccccccccc
c         print *,' grp='
c         print *, grp(:,ir)
cccccccccccccccccccccccc

c rotate k { IBZ
          call mv      (grp(1,ir),qibz(1,k),
     d              3,3,3,
     o              w1)

c identify the rotated k with k { FBZ
          kp         = 0

c> check if k{FBZ has not been identified
 1111     kp         = kp + 1
          if (iw1(kp) .eq. 0) then
cccccccccccccccccccccccccccccccccc
c         print *,' k ir kp=',k,ir,kp,' qibz=',qibz(:,k)
c         print *,' grp='
c         print *, grp(:,ir)
c         print *,' ginv '
c         print *, ginv
c         print *,' qbz '
c         print *, qbz(:,kp)
cccccccccccccccccccccccccccccccccc
c> check if the rotated k is connected to k{FBZ by a G
            call vminv   (qbz(1,kp),w1,3,w2)
            call mv      (ginv,w2,
     d              3,3,3, 
     o              w3)
            xsign      = .001d0*dsign(1.d0,w3(1))
            ysign      = .001d0*dsign(1.d0,w3(2))
            zsign      = .001d0*dsign(1.d0,w3(3))
            w2(1)      = w3(1) - dble(idint((w3(1)+xsign)))
            w2(2)      = w3(2) - dble(idint((w3(2)+ysign)))
            w2(3)      = w3(3) - dble(idint((w3(3)+zsign)))

            if ( vdv(w2,w2,3) .lt. tol ) then
              nstar(k)   = nstar(k) + 1
              irotk(k,ir)= kp
              iw1(kp)    = 1
            endif
          endif
cccccccccccccccccccccccccccccccccccc
c      print *,' qbz='
c      print *, qbz(:,kp)
c      print *,' w3='
c      print *, w3
c      print *,' w2='
c      print *, w2
cccccccccccccccccccccccccccccccccccc
          if (kp .lt. nqbz) goto 1111

c end of rotation-loop
cccccccccccccccccccccccccccc
c      print *,' nstar(k) =',nstar(k)
cccccccccccccccccccccccccccc
        end do

c end of k { IBZ -loop
      end do

c check that the sum of stars equal to the no. k{FBZ
      nsum       = ivsum (nstar,nqibz)
      if (nsum .ne. nqbz) stop 'nkstar: wrong no. stars'

c write k { IBZ and no. stars to file KPNT
      ifkpnt     = ifile('KPNT')
      if (ifkpnt .gt. 0) then
        write (ifkpnt,*) 'irreducible k-points and no. stars'
        write (ifkpnt,*) 'k, k-vector, nstar '
        do       k = 1,nqibz
          write (ifkpnt,6000) k,qibz(1,k),qibz(2,k),qibz(3,k),nstar(k)
        end do
      endif

c formats
 6000 format (1x,i5,3f8.5,i3)

      return
      end

