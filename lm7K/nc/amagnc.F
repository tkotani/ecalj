      subroutine amagnc(nbas,nl,ipc,rhos,nrhos,qnu,eula,neul,mode,
     .  amag,aamom)
C- Printout magnetic moments in unit cell
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   rhos  :spin density matrix (for mode=0)
Ci         :rhos should be hermitian in spin space, but may not be owing
Ci         :to energy integration errors in the complex plane.
Ci         :amagnc uses a symmetrized form to minimize errors.
Ci   qnu   :moments (for mode=1)
Ci   eula  :Euler angles for noncollinear spins
Ci   neul  :1 if Euler angles are l- and m-independent,
Ci         :nl if Euler are l-dependent and m-independent
Ci         :nl**2 if Euler are l- and m-dependent
Ci   mode  :0 use spin density matrix to make moments along T
Ci         :1 use qnus to make moments along qnu
Co Outputs
Co   amag(1..3): net system magnetic moment
Co   aamom :local magnetic moments
Cr Remarks
Cr   Definition of rho in terms of M: (standard definition of sigma)
Cr      rho = M . sigma/2
Cr   Pauli matrices sigma:
Cr
Cr              (0  1)             (0 -i)           (1  0)
Cr     sigmax = (    )    sigmay = (    )  sigmaz = (    )
Cr              (1  0)             (i  0)           (0 -1)
Cr   Given rho, M can be obtain from:
Cr     M_x =  2 Re(rho21) = Re (rho12+rho21)
Cr     M_y =  2 Im(rho21) = Im (rho21-rho12)
Cr     M_z =  (rho11)-(rho22)
Cr   Second (symmetrized) form is used because for numerical reasons,
Cr   rhos may not be properly hermitian, e.g. when rhos is generated
Cr   by a Green's function technique.
Cu Updates
Cu   21 Apr 04 Revised to properly accomodate m-dependent Euler angles
Cu   17 Feb 03 Revised amagnc; cleaner and bug fixes.
C ----------------------------------------------------------------------
C     implicit none
C Passed variables
      integer nbas,neul,nl,nrhos,ipc(*),mode
      double precision eula(nbas,neul,3),rhos(2,0:2,nrhos,2,2,1),
     .  qnu(3,nl,2,1),amag(3),aamom(nbas)
C Local variables
      logical lrhol
      integer i,ib,ic,lp1,lgunit,ipr,k,stdo,ilm,l,m
      double precision alphan,betan,gamman,arhol,arhom,sarho,
     .  rotg(3,3),amom,amlm(3),amgm(3),saml(3),samg(3),amgl(3),aml(3)
      integer PRT1,PRT2
      parameter (PRT1=30,PRT2=40)

C     character*1 cmode(2)
C     data cmode /'o','i'/

C --- Setup ---
      call getpr(ipr)
      stdo = lgunit(1)
      call isanrg(mode,0,1,' amagnc','mode',.true.)
      amag(1) = 0
      amag(2) = 0
      amag(3) = 0
      sarho = 0
      lrhol = nrhos .eq. nl .or. mode .eq. 1
      if (nrhos .ne. nl .and. nrhos .ne. nl**2 .and. mode .eq. 0)
     .  call rxi('amagnc: bad nrhos',nrhos)

      if (ipr .ge. 20) then
        if (mode.eq.0) write(stdo,332) 'density matrix'
        if (mode.eq.1) write(stdo,332) 'Euler angles and sphere charges'
        if (lrhol .and. ipr .ge. PRT1) write(stdo,335) '   l'
        if (.not. lrhol .and. ipr .ge. PRT1) write(stdo,335) ' ilm'
      endif

  332 format(/' AMAGNC: magnetic moments from ',a,':')
  335 format(16x,'global coordinates',15x,'local coordinates'/
     .        3x,'ib',a,'    Mx        My        Mz',
     .       10x,           'Mx        My        Mz        |M|')

C --- Get moments for each spin ---
      do  ib = 1, nbas
        call dpzero(saml,3)
        call dpzero(samg,3)
        ic = ipc(ib)
        aamom(ib) = 0
        ilm = 0
        do  l = 0, nl-1
        lp1 = l+1
        call dpzero(amgl,3)
        call dpzero(aml,3)
        arhol = 0
        do  m = -l, l
          ilm = ilm+1
          k = 1
          if (neul .eq. nl) k = lp1
          if (neul .eq. nl**2) k = ilm
          alphan = eula(ib,k,1)
          betan  = eula(ib,k,2)
          gamman = eula(ib,k,3)

          if (mode .eq. 0) then
C       ... Rotation matrix for T in local coordinates
C           This is M in local coordinates, symmetrized form
            k = lp1
            if (.not. lrhol) k = ilm
            amlm(1) = rhos(1,0,k,1,2,ic) + rhos(1,0,k,2,1,ic)
            amlm(2) = rhos(2,0,k,2,1,ic) - rhos(2,0,k,1,2,ic)
            amlm(3) = rhos(1,0,k,1,1,ic) - rhos(1,0,k,2,2,ic)
            if (lrhol) call dscal(3,1/dble(2*l+1),amlm,1)
          else
            amlm(1) = 0
            amlm(2) = 0
            amlm(3) = qnu(1,lp1,1,ic) - qnu(1,lp1,2,ic)
            call dscal(3,1/dble(2*l+1),amlm,1)
          endif
          arhom = dsqrt(amlm(1)**2 + amlm(2)**2 + amlm(3)**2)

C         Rotation matrix for input Euler angles
          call eua2rm(alphan,betan,gamman,rotg)
C         print 336, 'rotg',((rotg(i,k),k=1,3),i=1,3)
C 336     format(' rot matrix:',a,/(3f15.10))

C         Rotate amlm to global coordinates (amgm)
          call dgemm('T','N',3,1,3,1d0,rotg,3,amlm,3,0d0,amgm,3)
C         call dgemm('N','N',3,1,3,1d0,rotg,3,amlm,3,0d0,amgm,3)

          call daxpy(3,1d0,amgm,1,amgl,1)
          call daxpy(3,1d0,amlm,1,aml,1)
          arhol = arhol + arhom
          aamom(ib) = aamom(ib) + arhom

          if (.not. lrhol .and. ipr.ge.PRT2)
     .      write(stdo,333) ib,ilm,amgm,amlm,arhom
  333     format(i5,i4,3f10.6,2x,3f10.6,f10.6)
C         End of loop over m
          enddo

          call daxpy(3,1d0,aml,1,saml,1)
          call daxpy(3,1d0,amgl,1,samg,1)
          if (lrhol .and. ipr.ge.PRT2)
     .      write(stdo,333) ib,l,amgl,aml,arhol
        enddo
        aamom(ib) = dabs(aamom(ib))
        amom = dsqrt(samg(1)**2 + samg(2)**2 + samg(3)**2)
        sarho = sarho + dabs(amom)
        if (ipr.ge.PRT1) write(stdo,334) ib,samg,saml,amom
  334   format(i5,'   *',3f10.6,2x,3f10.6,f10.6)
        if (ipr.ge.PRT2) write(stdo,'(1x)')
        amag(1) = amag(1) + samg(1)/nbas
        amag(2) = amag(2) + samg(2)/nbas
        amag(3) = amag(3) + samg(3)/nbas
      enddo
      amom = dsqrt(amag(1)**2 + amag(2)**2 + amag(3)**2)
      call info5(20,0,0,' <Mx>=%,6;6d  <My>=%,6;6d  <Mz>=%,6;6d'//
     .  '  <M>=%,6;6d  <|M|>=%,6;6d',amag(1),amag(2),amag(3),amom,
     .  sarho/nbas)
      call awrit5(' <Mx>=%,6;6d  <My>=%,6;6d  <Mz>=%,6;6d'//
     .  '  <M>=%,6;6d  <|M|>=%,6;6d',
     .  ' ',80,lgunit(2),amag(1),amag(2),amag(3),amom,sarho/nbas)

      end
      subroutine amagn2(mode,nl,nlo,nclass,nbas,ipc,eula,neul,bxc,qnu,
     .  nrhos,rhos)
C- When rhos++ and rhos-- is only l-dependent, approximately adjust
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 use spin density matrix to make moments m3 along zhat
Ci         :1 use qnus to make moments along qnu
Ci   mode  :use qnu for rhos (debugging only)
Ci   nl    :(global maximum l) + 1
Ci   nbas  :size of basis
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   eula  :Euler angles for noncollinear spins
Ci   neul  :1 if Euler angles are l-independent, nl otherwise
Ci   qnu   :energy-weighted moments of the sphere charges
Ci   nrhos :number of channels which spin density-matrix is stored
Cio Inputs/Outputs
Ci   rhos  :spin density-matrix.  rhos++ and rhos-- are modified.
Cl Local variables
Cl         :
Cr Remarks
Cr   If because rhos++ and rhos-- were averaged over m (which can happen
Cr   if makwts accumulates doswts by l only; also for the ASA moments
Cr   they are always averaged over m), no account is
Cr   taken of their m-dependent orientation.  A better average to
Cr   rhos++ and rhos-- would be to make
Cr     Mbar = 1/(2l+1) sum_m Rm^-1 Mm(loc)
Cr   where Mm(loc) = m-dependent local magnetization in its local
Cr   coordinate system and Rm is the rotation matrix to loc. coords.
Cr   This is what would obtain if the entire rhos were averaged over m.
Cr   Finally, from Mbar we obtain the z-part of Mm in its local
Cr   coordinate system as Mmz(loc) = Rm Mbar.  We leave the x and y
Cr   parts untouched as they were explicitly resolved.
Cr   This is of course an approximation.  Better to allow the entire
Cr   rhos to be accumulated by m, which is accomplished by resolving
Cr   doswts by m (nlo = nl**2).
Cr
Cr   Correcting the ASA moments qnu+-qnu-.  A similar principle applies.
Cr   qnu+-qnu- has less information than the full rhos, because only
Cr   the z component is retained.  We construct qnu+-qnu- by requiring
Cr   that the rotation of qnu+-qnu- to the B-field axis matches
Cr   the rhos++ - rhos-- in this same reference.
Cr   If qnu+ - qnu- = (rho++ - rho--) in the B-field reference frame,
Cr   we assign qnu+ - qnu-  = (rho++ - rho--)/R_33 in the coordinates
Cr   of the orbital, where R is the rotation from the orbital frame
Cr   to the Bxc frame: R=Rb*Rm+ (see below).
Cr   When the orbital axis (spec by rm) and and Bxc align then
Cr   (Rb*Rm+)_33=1 and the correction vanishes.  This should always
Cr   occur in the rigid-spin approximation, when the each orbital
Cr   is collinear with Bxc.
Cu Updates
Cu   09 Apr 04  First created.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nl,nlo,nclass,nbas,neul,nrhos,ipc(nbas)
      double precision eula(nbas,neul,3),rhos(2,0:2,nrhos,2,2,nclass)
      double precision qnu(3,nl,2,1),bxc(3,*)
C ... Local parameters
      integer ilm,l,j,m,k,i,ic,ib,iclbsj
      double precision mbar(3,0:2),amg(3,0:2),ammi(3),ddot,rotm(3,3),
     .  alpha,beta,gamma,ql,amli(0:2),qp,qm,eulat(neul,3),ph,th,
     .  hatbxc(3),pi,rotb(3,3),rotmb(3,3),rotmb3,rhom(3),qli(0:2),fac

      pi = 4*datan(1d0)
      do  ic = 1, nclass
      ib = iclbsj(ic,ipc,-nbas,1)
      if (ib .gt. 0) then
      call dmscop(eulat,1,eula,nbas,ib,ib,1,neul*3,1,1,1d0)
      call dcopy(3,bxc(1,ic),1,hatbxc,1)
      ilm = 0
      do  l = 0, nl-1

C       Retain starting ilm for this l
        j = ilm

C   --- Adjust rhos when ++,-- part averaged over m ---
        if (nrhos .gt. nlo) then

C   ... Determine average M for this l
        call dpzero(mbar,3*3)
        do  m = -l, l
          ilm = ilm+1
C         Euler angles for this channel
          k = 1
          if (neul .eq. nl) k = l+1
          if (neul .eq. nl*nl) k = ilm
          alpha = eulat(k,1)
          beta  = eulat(k,2)
          gamma = eulat(k,3)
          call eua2rm(alpha,beta,gamma,rotm)
          do  i = 0, 2
            if (mode .eq. 0) then
              ammi(1) = rhos(1,i,ilm,1,2,ic) + rhos(1,i,ilm,2,1,ic)
              ammi(2) = rhos(2,i,ilm,2,1,ic) - rhos(2,i,ilm,1,2,ic)
              ammi(3) = rhos(1,i,ilm,1,1,ic) - rhos(1,i,ilm,2,2,ic)
            elseif (mode .eq. 1) then
              ammi(1) = 0
              ammi(2) = 0
              ammi(3) = (qnu(i+1,l+1,1,ic) - qnu(i+1,l+1,2,ic))/(2*l+1)
            endif
            call dgemm('T','N',3,1,3,1d0,rotm,3,ammi,3,0d0,amg(1,i),3)
            call daxpy(3,1d0/(2*l+1),amg(1,i),1,mbar(1,i),1)
          enddo
        enddo

C   ... Distribute z part over rhom
        ilm = j
        do  m = -l, l
          ilm = ilm+1
C         Euler angles for this channel
          k = 1
          if (neul .eq. nl) k = l+1
          if (neul .eq. nl*nl) k = ilm
          alpha = eulat(k,1)
          beta  = eulat(k,2)
          gamma = eulat(k,3)
          call eua2rm(alpha,beta,gamma,rotm)
          do  i = 0, 2
            if (mode .eq. 0) then
              ql  = rhos(1,i,ilm,1,1,ic) + rhos(1,i,ilm,2,2,ic)
C             rhom(3) = rhos(1,i,ilm,1,1,ic) - rhos(1,i,ilm,2,2,ic)
            elseif (mode .eq. 1) then
              ql  = (qnu(i+1,l+1,1,ic) + qnu(i+1,l+1,2,ic))/(2*l+1)
C             rhom(3) = (qnu(i+1,l+1,1,ic) - qnu(i+1,l+1,2,ic))/(2*l+1)
            endif
C           rhoz_m(loc) = rotm * M
            rhom(3) = ddot(3,rotm(3,1),3,mbar(1,i),1)
            qp = (ql+rhom(3))/2
            qm = (ql-rhom(3))/2
C            print 533, i,ilm,rhos(1,i,ilm,1,1,ic)-qp,
C     .                       rhos(1,i,ilm,2,2,ic)-qm
C  533       format(2i3,2f15.10)
            rhos(1,i,ilm,1,1,ic) = qp
            rhos(1,i,ilm,2,2,ic) = qm
          enddo
        enddo
        endif

C   --- Generate qnu+-qnu- from projection of rhos onto Bxc-field ---
C   ... Projection of m-averaged rhos onto Bxc-field
C       The coordinate system for rhos for this orbital is specified by
C       Euler angles (which generate Rm).  Then rhos(loc)=Rm*rhos(glob).
C       Let Rb = rotation of of Bxc to zhat.  (NB: |Bxc|=1).
C       Bxc.rhos = Rb*Bxc . Rb*rhos = zhat . (Rb*Rm+) rhos(loc)
C       The z component of (Rb*(Rm+)) (rho++ - rho--) is qnu+ - qnu-
C       in the Bxc coordinate system.  Because qnu+ - qnu- is stored
C       in the local coordinate system of Euler angles, we compute
C       qnu+ - qnu-  = m-averaged rotation of qnu+ - qnu- in the
C       Bxc frame: qnu+ - qnu- (loc-frame) = R qnu+ - qnu- (Bxc-frame)
C       where  R = 1/(2l+1) sum m (Rb*(Rm+))^-1 = 1/(2l+1) sum m Rm*Rb+
C       Branch has no effect if Bxc || all orbitals in this l channel
        if (neul .ne. 1 .or. nrhos .gt. nlo) then
        ilm = j
        call dpzero(amli,3)
        call dpzero(qli,3)
        rotmb3 = 0
C       Determine average M for this l, neglecting rhox and rhoy
        do  m = -l, l
          ilm = ilm+1
C         Euler angles for magnetic moment this channel
          k = 1
          if (neul .eq. nl) k = l+1
          if (neul .eq. nl*nl) k = ilm
          alpha = eula(ib,k,1)
          beta  = eula(ib,k,2)
          gamma = eula(ib,k,3)
          call eua2rm(alpha,beta,gamma,rotm)
C         Rotation matrix Rb rotates Bxc to zhat
          if (hatbxc(1)**2+hatbxc(2)**2 .ne. 0) then
            ph = datan2(hatbxc(2),hatbxc(1))
            th = datan2(dsqrt(hatbxc(1)**2+hatbxc(2)**2),hatbxc(3))
          else
            ph = 0
            th = 0
          endif
          call rotma(ph+pi/2,pi/2,th,rotb)
C         Rb * Rm+
          call dgemm('N','T',3,3,3,1d0,rotb,3,rotm,3,0d0,rotmb,3)
C         ammi: j=3 component of rhom = qnu+-qnu- in Bxc coord sys.
          if (nrhos .eq. nl*nl) then
            k = ilm
            fac = 1
          else
            fac = 1/dble(2*l+1)
            k = l+1
          endif
          do  i = 0, 2
            rhom(1) = (rhos(1,i,k,1,2,ic) + rhos(1,i,k,2,1,ic))*fac
            rhom(2) = (rhos(2,i,k,2,1,ic) - rhos(2,i,k,1,2,ic))*fac
            rhom(3) = (rhos(1,i,k,1,1,ic) - rhos(1,i,k,2,2,ic))*fac
            call dgemm('N','N',3,1,3,1d0,rotmb,3,rhom,3,0d0,ammi,3)
C           m-averaged rho+ + rho-, rho+ - rho- in Bxc frame
            qli(i) = qli(i) +(rhos(1,i,k,1,1,ic)+rhos(1,i,k,2,2,ic))*fac
            amli(i) = amli(i) + ammi(3)
          enddo
C         Accumulate sum_m (Rb * Rm+) zhat to rotate amli->loc frame
          rotmb3 = rotmb3 + rotmb(3,3)/(2*l+1)
        enddo
C       Make (qp-qm)_i = amli_i/R_33.  Then
C       sum Rb Rm+ (qnu+ - qnu-) = m-averaged rhos++-rhos-- in Bxc frame
        do  i = 0, 2
          amli(i) = amli(i)/rotmb3
          qp = (qli(i)+amli(i))/2
          qm = (qli(i)-amli(i))/2
C          if (i .eq. 0) print533,i,l,ic,
C     .      qnu(i+1,l+1,1,ic)+qnu(i+1,l+1,2,ic)-qp-qm,
C     .      qnu(i+1,l+1,1,ic)-qnu(i+1,l+1,2,ic)-(qp-qm),
C     .      amli(i)*rotmb3,amli(0),qp,qm
C  533     format(3i3,6f13.8)
          qnu(i+1,l+1,1,ic) = qp
          qnu(i+1,l+1,2,ic) = qm
        enddo
      endif
      enddo
      endif
      enddo

      end

