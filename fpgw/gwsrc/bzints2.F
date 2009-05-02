c taken from the gwx1/lib/tetra.for
      SUBROUTINE BZINTS2x(volwgt,EP,WP,NQ,nband,NB,NSP,EMIN,EMAX,
     .   DOS,NR,EF,JOB,NTET,IDTET)
C-  Bz integrations by linear method.
C ----------------------------------------------------------------------
Ci Inputs:
Ci  ep, energy bands;
Ci  nq, no. of irreducible k-points; nb, no. of bands; nsp, see BNDASA;
Ci  emin, emax, dos, nr : for IDOS, energy window,
Ci  emin, emax (JOB=2) tolerance in efermi
Ci  IDOS, number of points; ef, Fermi energy (job = 2); job, switch :
Ci  JOB=1: MAKES IDOS.  JOB=2:  MAKES BLOECHL-WEIGHTS.
Ci  ntet, No. of different tetrahedra
Ci  idtet(1-4,i), Identifies the i'th tetrahedron in terms of the four
Ci  irreducible k-points:
Ci  idtet(0,i), no. of tetrahedra of the i'th kind
Co Outputs:
Co  dos, Integrated Density of States (IDOS) (job = 1)
Co  wp, Bloechl quadrature weights (job = 2)
Cm Memory:
Cm  No large internal storage; heap not accessed.
C ----------------------------------------------------------------------
c      implicit none
C Passed parameters
C Local parameters
      IMPLICIT double precision (A-H,O-Z)
      integer fopn
      DIMENSION EP(nb,nsp,nq),DOS(NR),EC(4),WC(4,2),WP(nband,nsp,nq),
     .   idtet(0:4,*)
ccccccccccccccccccccccccccccccccccc
c        do  is  = 1,nsp !Readin eband
c        do  iqi = 1,nq
c        do   it = 1,nb
c          write(6,'("is q ib eband =",3i5,d13.6)')is,iqi,it,ep(it,is,iqi)
c        enddo
c        enddo
c        enddo
c        stop 'test end xxxxxxxxxx'
cccccccccccccccccccccccccccccccccccc


      IF (JOB .NE. 1 .AND. JOB .NE. 2) STOP '*** BAD JOB IN BZINTS2x'
c      IF (JOB .EQ. 1) call dinit(dos,2*nr)
c takao
      IF (JOB .EQ. 1) call dinit(dos,nr)
      IF (JOB .EQ. 2) call dinit(wp,nband*nsp*nq)
      SEV1 = 0.D0
      SEV2 = 0.D0
c      volwgt = (3.d0 - nsp) / (n1*n2*n3*6.d0)
      do  40  isp = 1, nsp
C ----- START LOOPING OVER TETRAHEDRA ---------
        DO  20  ITET = 1, NTET
          iq1=idtet(1,itet)
          iq2=idtet(2,itet)
          iq3=idtet(3,itet)
          iq4=idtet(4,itet)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      if(iq1==iq2.and.iq2==iq3.and.iq3==iq4) then
c        write(6,"( 'itet= ',i5,' iq= ',4i5,' vol=',i4)")
c     &  itet,iq1,iq2,iq3,iq4,idtet(0,ITET)
c      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        DO  20  IB = 1, nb !nband
C ----- SET UP ENERGIES AT 4 CORNERS OF TETRAHEDRA ------
          ec(1) = ep(ib,isp,iq1)
          ec(2) = ep(ib,isp,iq2)
          ec(3) = ep(ib,isp,iq3)
          ec(4) = ep(ib,isp,iq4)
cccccccccccccccccccccc
c       write(6,"('bzint2x ib E=',i4,4d13.6)") ib,EC(1:4)
cccccccccccccccccccccc
          etop = dmax1(ec(1),ec(2),ec(3),ec(4))
          ebot = dmin1(ec(1),ec(2),ec(3),ec(4))
          IF (JOB .EQ. 1) THEN
            if ( ebot .lt. emax )
c     .         CALL SLINZ(volwgt*idtet(0,ITET),EC,EMIN,EMAX,DOS,NR)
c takao
     .         CALL SLINZ2(volwgt*idtet(0,ITET),EC,EMIN,EMAX,DOS,NR)
cccccccccccccccccccccccccccccccccccc
c      if(abs(emax-emin).lt.1.0d-6) then
c      if(abs(dos(46)-dos(47))>8*1d-4) then
c        write(6,*)dos(46),dos(47)
c        write(6,*)isp,itet,iq1,iq2,iq3,iq4,ec
c      stop
c      endif
c      endif
cccccccccccccccccccccccccccccccccccc

          ELSE
            if ( ef .ge. ebot ) then
              CALL FSWGTS(volwgt*idtet(0,ITET),EC,EF,ETOP,WC)
              SEV1 = SEV1 + WC(1,1)*EC(1) + WC(2,1)*EC(2) +
     .                       WC(3,1)*EC(3) + WC(4,1)*EC(4)
              SEV2 = SEV2 + WC(1,2)*EC(1) + WC(2,2)*EC(2) +
     .                       WC(3,2)*EC(3) + WC(4,2)*EC(4)
              WP(ib,isp,iq1) = WP(ib,isp,iq1) + WC(1,1) + WC(1,2)
              WP(ib,isp,iq2) = WP(ib,isp,iq2) + WC(2,1) + WC(2,2)
              WP(ib,isp,iq3) = WP(ib,isp,iq3) + WC(3,1) + WC(3,2)
              WP(ib,isp,iq4) = WP(ib,isp,iq4) + WC(4,1) + WC(4,2)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c takao for test
c              print 123, itet,
c     &                   WC(1,1)+WC(2,1)+WC(3,1)+WC(4,1),
c     &                   WC(1,2)+WC(2,2)+WC(3,2)+WC(4,2),
c     &                   WC(1,1)+WC(2,1)+WC(3,1)+WC(4,1)
c     &                 + WC(1,2)+WC(2,2)+WC(3,2)+WC(4,2),
c     &                   volwgt*idtet(0,ITET)
c  123         format(' itet=',i4,' wc1 sum=',d23.16,' wc2 sum=',d23.16,
c     &               ' wc1+wc2 sum=',d23.16,
c     &               ' volwgt*idtet(0)=',d23.16)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            endif
          ENDIF
   20   CONTINUE
   40 continue
ccccccccccccccccccccccccccccc
c      stop 'takao test'
ccccccccccccccccccccccccccccc

c$$$      if ( iprint() .ge. 30 ) then
c$$$        IF (JOB .EQ. 2) THEN
c$$$          SUMWP = 0.D0
c$$$          do  56  isp = 1, nsp
c$$$            DO  56  IB = 1, nband
c$$$            DO  56  IQ = 1, NQ
c$$$  56      SUMWP = SUMWP + WP(ib,isp,iq)
c$$$          if (iprint() .ge. 30) print 922, ef, sumwp, sev1+sev2, sev2
c$$$c          kkk=fopn('LOG')
c$$$c          write(kkk,922) ef, sumwp, sev1+sev2, sev2
c$$$        ENDIF
c$$$      endif

  922 format(1x,'BZINTS2x: Fermi energy:',f10.6,';',F20.16,' electrons'/
     .       9x,'Band energy:',f11.6,
     .       ', including Bloechl correction:',f10.6)
      END
c==========================================================================



      SUBROUTINE SLINZ2(VOLWGT,EC,EMIN,EMAX,DOSI,NR)
C- Adds to number-of-states for one tetrahedron
C ----------------------------------------------------------------
Ci Inputs
Ci   volwgt, weight on tetrahedron; ec energies at corners of tethdn.;
Ci   emin, emax, energy window; nr, number of bins + 1
Co Outputs
Co   dosi(k), integrated density in kth bin from tethdn.
Cr Remarks
Cr
C ----------------------------------------------------------------
C      implicit none
C Passed parameters
C Local parameters
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION EC(4),DOSI(NR)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      print *,' slinz2: volwgt=',volwgt
c      print *, d1mach(3)
c      volsum=0.0d0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DO  3  I = 1, 3
      DO  3  J = 1, 4-I
        IF (EC(J) .GT. EC(J+1)) THEN
          E=EC(J)
          EC(J) = EC(J+1)
          EC(J+1) = E
        ENDIF
    3 continue
      E1 = EC(1)
      E2 = EC(2)
      E3 = EC(3)
      E4 = EC(4)
      if (e4 .lt. emin) then
        i4=1
        go to 26
      endif
      DE = (EMAX-EMIN)/(NR-1)

      d2 = 2.0d0*( 1.0d0-d1mach(3) )



cx takao ---------------
cx This is a correction in order to get the very accurate Fermi energy.
cx
      E1x=(E1-EMIN)/DE
      if( e1x .ge. dble(nr+10) ) e1x = nr + 10
      if( e1x. le. -10.0d0)      e1x =  - 10
      I01 = e1x + d2

      E2x=(E2-EMIN)/DE
      if( e2x .ge. dble(nr+10) ) e2x = nr + 10
      if( e2x. le. -10.0d0)      e2x =  - 10
      I02 = e2x + d2

      E3x=(E3-EMIN)/DE
      if( e3x .ge. dble(nr+10) ) e3x = nr + 10
      if( e3x. le. -10.0d0)      e3x =  - 10
      I03 = e3x + d2

      E4x=(E4-EMIN)/DE
      if( e4x .ge. dble(nr+10) ) e4x = nr + 10
      if( e4x. le. -10.0d0)      e4x =  - 10
      I04 = e4x + d2
cx this was sometimes poor and cause a problem when you want to get very accurate Ef.
c      write(6,"('EEEE=',4d13.6)") E1, E2, E3, E4
c      print *, I01, I02, I03, I04
c      I01   = (E1   -EMIN)/DE + 1.9999999D0
c      I02   = (E2   -EMIN)/DE + 1.9999999D0
c      I03   = (E3   -EMIN)/DE + 1.9999999D0
c      I04   = (E4   -EMIN)/DE + 1.9999999D0


C --------------------------------
      I1 = MAX0(I01  ,1)
      I2 = MIN0(I02-1, NR)
      IF (I1 .LE. I2) THEN
        CC = VOLWGT/((E2-E1)*(E3-E1)*(E4-E1))
        DO  20  I = I1, I2
          X = EMIN - E1 + (I-1)*DE
          DOSI(I) = DOSI(I) + CC*X**3
   20   continue
      ENDIF

      I2 = MAX0(I02  ,1)
      I3 = MIN0(I03-1,NR)
      IF (I2 .LE. I3) THEN
        C3 = VOLWGT*(E1+E2-E3-E4)/((E3-E1)*(E4-E1)*(E3-E2)*(E4-E2))
        C2 = VOLWGT*3.D0/((E3-E1)*(E4-E1))
        C1 = C2*(E2-E1)
        C0 = C1*(E2-E1)/3.D0
        DO  21  I = I2, I3
          X = EMIN - E2 + (I-1)*DE
          DOSI(I) = DOSI(I) + C0 + X*(C1 + X*(C2 + X*C3))
   21   continue
      ENDIF

      I3 = MAX0(I03  ,1)
      I4 = MIN0(I04  -1,NR)
      IF (I3 .LE. I4) THEN
        CC = VOLWGT/((E3-E4)*(E2-E4)*(E1-E4))
        DO  22  I = I3, I4
          X = EMIN - E4 + (I-1)*DE
          DOSI(I) = DOSI(I) + VOLWGT - CC*X**3
   22   continue
      ENDIF

      I4 = MAX0(I04  ,1)
   26 continue
      DO  25  I = I4, NR
        DOSI(I) = DOSI(I) + VOLWGT
   25 continue
      END


c -------------------------------------
c modified a bit from LMTO3 /ferdi/lmto/bndasa2.f/efrang
c e1=e2 is just middle of HOMO and LUMO for insulator
      subroutine efrang3(nsp,nkp,nbmax,nband,zval,e1,e2,eband,
     .                  ebbot,ebtop,ebbott,ebtopp,elo,ehi,bandgap)
C- Find range of Fermi energy.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nsp,nkp,nbmax,nband,eband
Ci   zval no. of valence electrons
Co Outputs
Co   e1,e2: e1 < ef < e2
Co   elo, ehi:  lowest and highest band found
Cr Remarks
Cr    For an even no. of electrons ef is between the bottom of the
Cr    zval/2+1'th band and the top of the zval/2 'th band. If the
Cr    former is higher that the latter we have an insulator and e1=e2.
Cr    For an odd no. of electrons ef is between the bottom and top
Cr    of the (zval+1)/2 'th band.
Cr    For spin pol case: bottom of the zval+1'th band < ef <
Cr                       top of the zval'th band.
Cr                       If the bottom is higher that the top then
Cr                       we have an insulator and e1=e2.
Cr    If nsp is zero, the upper and lower band limits are returned
C ----------------------------------------------------------------------

C Passed parameters
      integer nsp,nkp,nbmax,nband
      double precision zval,e1,e2,eband(nbmax,nsp,nkp),ebbot(nband,nsp),
     .                 ebtop(nband,nsp),ebbott(*),ebtopp(*),elo,ehi
C Local parameters
      double precision e,d1mach
      integer i,j,ikp,isp,iba,nval,nbbot,nbtop
      logical ::  zvalisinteger
C External calls:
      external d1mach
      real(8):: bandgap
      elo = d1mach(2)
      ehi = -elo
      nval = zval + 1d-4

      if ( nval .eq. zval ) then
        zvalisinteger=.true.
      else
        zvalisinteger=.false.
        write(*,*) 'efrang3, comment:  zval is not an integer'
      endif

      if(nsp==1) goto 1
      if(nsp==2) goto 2

C --- nsp = 1 ---
    1 continue
c      if( (nval/2) * 2 .eq. nval ) then
cmod Apr.21.1999
      if ( zvalisinteger ) then
        if (nval.eq.0) then
          nbtop = 1
          nbbot = 1
        else if ( (nval/2) * 2 .eq. nval ) then  ! even integer
          nbtop = nval / 2
          nbbot = nbtop + 1
        else
          nbtop = nval / 2 + 1    ! odd integer or others
          nbbot = nbtop
        endif
      else    ! zval is not integer
        nbbot = nval / 2  +1
        nbtop = nbbot
      endif
c e1 is for bottom of conduction.
c e2 is for top of valence.  for insulator.
c e1<e2 for metal
      e1 = eband(nbbot,1,1)
      e2 = eband(nbtop,1,1)
      elo = eband(1,1,1)
      ehi = eband(nband,1,1)
      do  20  ikp = 2 , nkp
        if( eband(1,1,ikp) .lt. elo ) elo = eband(1,1,ikp)
        if( eband(nband,1,ikp) .gt. ehi ) ehi = eband(nband,1,ikp)
        if( eband(nbbot,1,ikp) .lt. e1 ) e1 = eband(nbbot,1,ikp)
        if( eband(nbtop,1,ikp) .gt. e2 ) e2 = eband(nbtop,1,ikp)
   20 continue
cxxx
c      print *,' xxxxxx efrang: e1 e2=',e1,e2
c      if( e1 .gt. e2 ) e1 = e2
c if e1>e2 . This is insulator case.
      if( e1 > e2 ) then
        bandgap = e1-e2
        em = 0.5d0*(e1+e2)
        e1 = em
        e2 = em
      else
        bandgap=0d0
      endif
      return
cxxx

C --- nsp = 2 ---
    2 continue
      if (zvalisinteger) then
        do  30  isp = 1, nsp
        do  30  iba = 1, nband
          ebbot(iba,isp) = eband(iba,isp,1)
          ebtop(iba,isp) = eband(iba,isp,1)
   30   continue
      else
        do  31  isp = 1, nsp
        do  31  iba = 1, nband
          ebbot(iba,isp) = minval(eband(iba,isp,1:nkp))
          ebtop(iba,isp) = maxval(eband(iba,isp,1:nkp))
   31   continue
      endif

      do  40  ikp = 1, nkp
      do  40  isp = 1, nsp
      do  40  iba = 1, nband
        if( eband(1,isp,ikp) .lt. elo ) elo = eband(1,isp,ikp)
        if( eband(nband,isp,ikp) .gt. ehi ) ehi = eband(nband,isp,ikp)
        if( eband(iba,isp,ikp) .lt. ebbot(iba,isp) )
     .                ebbot(iba,isp) =  eband(iba,isp,ikp)
        if( eband(iba,isp,ikp) .gt. ebtop(iba,isp) )
     .                ebtop(iba,isp) =  eband(iba,isp,ikp)
   40 continue
      do  50  i = 1, nband*nsp - 1
      do  50  j = 1, nband*nsp - i
        if( ebbott(j) .gt. ebbott(j+1) ) then
          e = ebbott(j)
          ebbott(j) = ebbott(j+1)
          ebbott(j+1) = e
        endif
        if( ebtopp(j) .gt. ebtopp(j+1) ) then
          e = ebtopp(j)
          ebtopp(j) = ebtopp(j+1)
          ebtopp(j+1) = e
        endif
   50 continue

      if (zvalisinteger) then
        e1 = ebbott(nval+1)
        if( nval+1 .gt. nband*nsp ) e1 = ebbott(nval)
        e2 = ebtopp(nval)
      else
        e1 = ebbott(nval+1)
        e2 = ebtopp(nval+1)
      endif

cxxx
c      if( e1 .gt. e2 ) e1 = e2
      if( e1 > e2 ) then !insulator case
        bandgap = e1-e2
        em = 0.5d0*(e1+e2)
        e1 = em
        e2 = em
      else
        bandgap=0d0
      endif
cxxx
      end

c--------------------------------------------------
      subroutine getvaln(
     i             konfig,z,nl,natom,iclass,nclass, 
     i             valn)
C - Get valn
Co valn    = number of valence electron.
c
      implicit none
      integer(4):: nclass,natom,nl,ia,ic,l
      real(8)   :: valn,ef, z(nclass)
      integer(4):: iclass(natom),konfig(0:nl-1,nclass)
      valn    = 0d0
      do ia   = 1,natom
        ic    = iclass(ia)
        valn  = valn + z(ic)
        do    l = 0,nl-1
          valn  = valn - (konfig(l,ic)-l-1) *( 2*l +1)*2
        end do
      end do
      print *,' getvaln: valn=',valn
      end









