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
Cstop2rx 2013.08.09 kino        stop 'fbz2: qx**2>1d-8 something wrong!!!'
        call rx( 'fbz2: qx**2>1d-8 something wrong!!!')
      endif
      end

      SUBROUTINE TETIRR_notused(QB,N1,N2,N3,IPQ,NQ,NTET,IDTET)
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
      implicit integer(i-n)
      DIMENSION QB(3,3),IPQ(N1,N2,N3),QB1(3,3),
     .          IBTR(3,3),KCUT(3,4,6),IMC(0:1,0:1,0:1),
     .          idtet(0:4,*),iq(4)

      CALL CCUTUP(QB,QB1,IBTR,KCUT)
      ntet = 0
C ----- START LOOPING OVER MICROCELLS ---------
      DO  202  I3 = 1, N3
      DO  201  I2 = 1, N2
      DO  20  I1 = 1, N1
C ----- SET UP IDENTIFIERS AT 8 CORNERS OF MICROCELL ------
        DO   K1 = 0, 1
          J1 = MOD(I1 + K1 -1,N1) + 1
        DO   K2 = 0, 1
          J2 = MOD(I2 + K2 -1,N2) + 1
        DO   K3 = 0, 1
          J3 = MOD(I3 + K3 -1,N3) + 1
          IMC(K1,K2,K3)  = IPQ(J1,J2,J3)
        enddo
        enddo
        enddo
C ----- LOOP OVER TETRAHEDRA --------------
        DO  10  ITET = 1, 6
          DO  2  IC = 1, 4
            K1 = KCUT(1,IC,ITET)
            K2 = KCUT(2,IC,ITET)
            K3 = KCUT(3,IC,ITET)
            IQ(IC) = IMC(K1,K2,K3)
    2     continue

C ----- order the identifiers -----
          DO   J=1,3
          DO   I=1,4-J
           IF(IQ(I).GT.IQ(I+1))THEN
              II=IQ(I)
              IQ(I)=IQ(I+1)
              IQ(I+1)=II
           ENDIF
          enddo
          enddo
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
          write(6,*)' nnnnn i1i2i3=',i1i2i3,idtet(:,ntet)
   10   CONTINUE
   20 CONTINUE
 201  CONTINUE
 202  CONTINUE
      print 922, ntet, 6*n1*n2*n3
  922 format(1x,'TETIRR: ',i6,' inequivalent tetrahedra from ',i8)
      END

      subroutine fermi2(
     i     qval,dos,ndos,emin,emax,
     o     eferm,e1,e2,dosef)
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
      implicit none
C Passed parameters
      integer:: ndos
      real(8):: qval,dos,emin,emax,eferm,e1,e2,dosef
C Local parameters
      integer:: i1,ie
      real(8):: de,q,q1,q2,d1mach
C External procedures
      external:: d1mach

      DIMENSION DOS(NDOS)

ccccccccccccccccccc
c      DE = (EMAX-EMIN)/(NDOS-1)
c      do i1=1,ndos
c        write(6,"(' e dostot =',i4,d13.6,d13.6)")i1, emin + de*(i1 - 1),dos(i1)
c      enddo
c      write(6,"(' qval =',d13.6)") qval
ccccccccccccccccccccc

      if (dos(1) .gt. qval) print *, 'FERMI: EMIN,EMAX=', emin,emax
      if (dos(1) .gt. qval) call rx( 'FERMI: Fermi energy lies below EMIN')
      if (dos(ndos) .lt. qval) print *, 'FERMI: EMIN,EMAX=', emin,emax
      if (dos(ndos) .lt. qval) then
        call rx( 'FERMI: Fermi energy lies above EMAX')
      endif

      DE = (EMAX-EMIN)/(NDOS-1)
      I1 = 1
      q = qval + d1mach(3)
      DO  1  IE = 2, NDOS
        I1 = IE
        IF ( DOS(IE) .gt. q ) goto 2
    1 continue
    2 continue
      i1 = i1 - 1
      Q1 = DOS(I1)
      Q2 = DOS(I1+1)
C ------------------
      e1 = emin + de*(i1 - 1)
      e2 = e1 + de
      !write(6,"('fermi:: i1,e1,qval,q1,q2,de',I7,5E13.5)") i1,e1,qval,q1,q2,de
      EFERM = e1 + (QVAL-Q1)/(Q2-Q1)*DE
      dosef = (q2-q1)/de
      end subroutine fermi2
c-------------------------------------------------------------
      subroutine fermi_kbt(
     i               qval,dos,ndos,emin,emax, kbt, eferm_init,
     o               eferm_kbt)
      implicit none
      intent(in)::   qval,dos,ndos,emin,emax, kbt, eferm_init
      intent(out)::  eferm_kbt
C- Makes fermi energy from integrated density
C ----------------------------------------------------------------------
Ci Inputs
Ci   qval, number of electrons to fermi level; dos(i) integrated
Ci   density at bin i; ndos, number of bins + 1; emin, emax, energy
Ci   window.
Co Outputs
Co   Eferm, Fermi energy; e1, e2, confidence limits on Fermi energy
Co   i.e., Fermi energy lies between e1 and e2.
Co   dosef_kbt:  density of states at fermi level
Cr Remarks
Cr   emin and e1 (and emax and e2) may point to the same address.
C ----------------------------------------------------------------------
      integer ndos
      real(8):: qval,dos,emin,emax,ddos,eferm_old,efermi_init
      real(8):: eferm_init,eferm_kbt,e0,e1,e2,dosef_kbt
      real(8):: kbt, ec
      integer:: i1,ie,ii,nii=100, ii_conv=-999, i
      real(8):: de,q,q1,q2,d1mach
      real(8):: idos, idos_old, idos_l, rydberg
      real(8):: eferm_i(4) !!Efermi(i+1),Efermi(i),Efermi(i-1),Efermi(i-2)
C External procedures
      external d1mach
      logical:: debug = .false.
      logical:: lessqval, moreqval

      DIMENSION dos(ndos)

      if (debug) then
         write(6,"('ndos',I8)") ndos
         write(6,"('qval,dos(1),dos(ndos)',3E13.5)") qval,dos(1),dos(ndos)
         write(6,"('emin,emax',2E13.5)") emin,emax
         write(6,"('kbt,efermi_init',2E13.5)") kbt,efermi_init
      endif
c$$$      do ii = 2,ndos
c$$$         write(6,*) "aaaaaaaa:",ii,dos(ii)-dos(ii-1)
c$$$      enddo
c$$$      call rx("end")      
      if (dos(1) .gt. qval) print *, 'FERMI: EMIN,EMAX=', emin,emax
      if (dos(1) .gt. qval) call rx( 'FERMI: Fermi energy lies below EMIN')
      if (dos(ndos) .lt. qval) print *, 'FERMI: EMIN,EMAX=', emin,emax
      if (dos(ndos) .lt. qval) then
        call rx( 'FERMI: Fermi energy lies above EMAX')
      endif

      de = (emax-emin)/(ndos-1)
      q = qval + d1mach(3)
      !write(6,"('fermi_kbt_plot:: ## fermi level (T=0 K)',E13.5)") eferm_init
      !!! self-consistently determined Fermi energy Efermi(T)
      eferm_i=-9999; eferm_i(1)=eferm_init

      lessqval = .false.
      moreqval = .false.
      write(6,*) "fermi_i: iter, EF(i) EF(i-1) EF(i-2) EF(i-3)"
      write(6,"(' fermi_i  iter =',I5,':',4E13.5)") 0,eferm_i(1:4)
      do 200 ii = 1,nii
         if (abs(eferm_i(2)-eferm_i(1)) < 1d-8 ) cycle !! check convergence
         idos   = dos(1) !!initial
         idos_l = dos(1) !!initial; left idos is necessarry for T=0 K
         do 300 ie = 2, ndos
            ddos = dos(ie)-dos(ie-1)
            e0 = emin + de*(ie - 2)
            e1 = emin + de*(ie - 1)
            e2 = e1 + de
            ec = (e1+e2)/2d0
            idos = idos + ddos*(1d0/(exp((ec-eferm_i(1))/kbt)+1))
            idos_l = idos_l + ddos*(1d0/(exp((e0-eferm_i(1))/kbt)+1))
            if (debug) then
               write(6,"('weight of fermi distribution',2E13.5)") 
     &              1d0/(exp((ec-eferm_old)/kbt)+1)!,1d0/(exp((e0-eferm_old)/kbt)+1)
               write(6,"('e1,e2,ec,eferm_old',4E13.5)") e1,e2,ec,eferm_old
               write(6,"('ec,dos(i1),idos',I6,5E13.5)") ec,dos(i1),idos,idos_l
            endif
300      continue

         !!! shift Efermi(i-1) <-- Efermi(i) for Efermi(i+1)
         do i = 3,1,-1
            eferm_i(i+1) = eferm_i(i)
         enddo
         !!! new Fermi energy Efermi(i+1) is determined in two cases (1) and (2)
         !!! (1) Fermi energy is between Efermi(i-1) and Efermi(i)
         if (lessqval .and. moreqval .and. ii >= 2) then
            if (idos_l <= qval) then !! Eferm(i) <= qval
               eferm_i(1)=(eferm_i(2) + eferm_i(3))/2d0
            else                !! qval < idos_l !! Eferm(i) > qval
               eferm_i(1)=(eferm_i(2) + eferm_i(4))/2d0
            endif
         !!! searching Fermi energy at finite temperature
         !!! (2) Fermi energy is shifted by kbt/10
         else
            if (idos_l <= qval) then
               eferm_i(1)=eferm_i(2)+kbt/real(ii) !!! E(i+1)=E(i)+kbt
               lessqval = .true.
            else                !! qval < idos_l (Efermi moves to lower energy)
               eferm_i(1)=eferm_i(2)-kbt/real(ii) !!! E(i+1)=E(i)+kbt
               moreqval = .true.
            endif
         endif
         !eferm_i(1)=-8888
         write(6,"(' fermi_i  iter =',I5,':',4E13.5)") ii,eferm_i(1:4)
         ii_conv=ii
200   continue !! unless converged 
      if (ii_conv == nii) then
         write(6,*) "fermi_kbt:: ----------- Note -----------"
         write(6,*) "fermi_kbt:: EFermi(T) is not converged"
         write(6,*) "fermi_kbt:: ----------------------------"
      endif
      write(6,*) "--------------------------------------------------------------------"
      write(6,"(' fermi_kbt:: Efermi(T), Efermi(0) [Ry]',2f15.7)") eferm_i(1),eferm_init
      write(6,"(' fermi_kbt:: Efermi(T), Efermi(0) [eV]',2f15.7)") eferm_i(1)*rydberg(),eferm_init*rydberg()
      write(6,*) "--------------------------------------------------------------------"
      !!! return
      eferm_kbt = eferm_i(1)
      end subroutine fermi_kbt
c-------------------------------------------------------------
ccc old version
c$$$
c$$$      if (debug) then
c$$$         write(6,"('ndos',I8)") ndos
c$$$         write(6,"('qval,dos(1),dos(ndos)',3E13.5)") qval,dos(1),dos(ndos)
c$$$         write(6,"('emin,emax',2E13.5)") emin,emax
c$$$         write(6,"('kbt,efermi_init',2E13.5)") kbt,efermi_init
c$$$      endif
c$$$c$$$      do ii = 2,ndos
c$$$c$$$         write(6,*) "aaaaaaaa:",ii,dos(ii)-dos(ii-1)
c$$$c$$$      enddo
c$$$c$$$      call rx("end")      
c$$$      if (dos(1) .gt. qval) print *, 'FERMI: EMIN,EMAX=', emin,emax
c$$$      if (dos(1) .gt. qval) call rx( 'FERMI: Fermi energy lies below EMIN')
c$$$      if (dos(ndos) .lt. qval) print *, 'FERMI: EMIN,EMAX=', emin,emax
c$$$      if (dos(ndos) .lt. qval) then
c$$$        call rx( 'FERMI: Fermi energy lies above EMAX')
c$$$      endif
c$$$
c$$$      de = (emax-emin)/(ndos-1)
c$$$      q = qval + d1mach(3)
c$$$      eferm  = eferm_init
c$$$      !write(6,"('fermi_kbt_plot:: ## fermi level (T=0 K)',E13.5)") eferm_init
c$$$
c$$$      !!! self-consistently determined Fermi energy Efermi(T)
c$$$      do 200 ii = 1,nii
c$$$         if (abs(eferm-eferm_old) < 1d-8 ) cycle !! check convergence
c$$$         idos   = dos(1)        !!initial
c$$$         idos_l = dos(1)        !!initial; left idos is necessarry for T=0 K
c$$$         eferm_old = eferm
c$$$         do 300 ie = 2, ndos
c$$$            if ( idos_l .gt. q ) cycle
c$$$            idos_old = idos
c$$$            ddos = dos(ie)-dos(ie-1)
c$$$            e0 = emin + de*(ie - 2)
c$$$            e1 = emin + de*(ie - 1)
c$$$            e2 = e1 + de
c$$$            ec = (e1+e2)/2d0
c$$$            idos = idos + ddos*(1d0/(exp((ec-eferm_old)/kbt)+1))
c$$$            idos_l = idos_l + ddos*(1d0/(exp((e0-eferm_old)/kbt)+1))
c$$$            i1 = ie
c$$$            if (debug) then
c$$$               write(6,"('weight of fermi distribution',2E13.5)") 
c$$$     &              1d0/(exp((ec-eferm_old)/kbt)+1)!,1d0/(exp((e0-eferm_old)/kbt)+1)
c$$$               write(6,"('e1,e2,ec,eferm_old',4E13.5)") e1,e2,ec,eferm_old
c$$$               write(6,"('i1,ec,dos(i1),idos',I6,5E13.5)") i1,ec,dos(i1),idos,idos_l
c$$$            endif
c$$$300      continue
c$$$         q1 = dos(i1)
c$$$         q2 = dos(i1+1)
c$$$         eferm = e1 + (qval-q1)/(q2-q1)*de
c$$$         !dosef_kbt = (q2-q1)/de
c$$$         !write(6,"('fermi_kbt: i1,e1,qval,q1,q2,de',I7,5E13.5)") i1,e1,qval,q1,q2,de
c$$$         write(6,"('fermi_kbt: iter,efermi,efermi_old',I5,2E13.5)") ii,eferm,eferm_old
c$$$         ii_conv=ii
c$$$200   continue !! unless converged 
c$$$      if (ii_conv == nii+1) then
c$$$         write(6,*) "fermi_kbt:: ----------- Note -----------"
c$$$         write(6,*) "fermi_kbt:: EFermi(T) is not converged"
c$$$         write(6,*) "fermi_kbt:: ----------------------------"
c$$$      endif
c$$$      write(6,*) "----------------------------------------------------------"
c$$$      write(6,"('fermi_kbt:: Efermi(T), Efermi(0)',2E17.8)") eferm,eferm_init
c$$$      !call rx("STOP routine")
c$$$      write(6,"('emin,emax',2E13.5)") emin,emax
c-------------------------------------------------------------

c-------------------------------------------------------------
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
      implicit integer(i-n)
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
      implicit integer(i-n)
      dimension g(3*3,ng)
      dimension invg(ng)
      dimension ginv(3*3),w1(9)
      data tol /1.d-5/
      do       i = 1,ng
        call minv33(g(1,i),ginv)
!! find g = ginv
        do       j = 1,ng
c           call vminv   (g(1,j),ginv,9,w1)
          w1=g(:,j)-ginv
          diff       = sum(w1**2) !vdv (w1,w1,9)
          if (diff .lt. tol) goto 1111
        end do
 1111   invg(i)    = j
        if (j .gt. ng) call rx( 'invgrp: cannot identify inverse of g')
      end do
      return
      end

c$$$      integer function invrot (irot,invg,ngrp)
c$$$c irot = label for R
c$$$c invg(R) = R^(-1)
c$$$c ngrp = no. rotation matrices
c$$$      implicit real*8 (a-h,o-z)
c$$$      integer:: irot,ngrp, invg(ngrp)
c$$$      invrot     = invg(irot)
c$$$      return
c$$$      end
