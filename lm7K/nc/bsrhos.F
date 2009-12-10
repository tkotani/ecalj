      subroutine bsrhos(nbas,nl,ipc,rhos,nrhos,qnu,pp,sop,eula,neul,
     .  bxc,bsite,nbf,lihdim,indxsh,mode,bsigr)
C- Double-counting term <B.sigma.rho>
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   rhos  :spin density matrix (for mode=0), by class
Ci         :rhos should be hermitian in spin space, but may not be owing
Ci         :to energy integration errors in the complex plane.
Ci         :bsrhos uses a symmetrized form to minimize errors.
Ci   qnu   :moments (for mode=1)
Ci   eula  :Euler angles for noncollinear spins
Ci   neul  :1 if Euler angles are l- and m-independent,
Ci         :nl if Euler are l-dependent and m-independent
Ci         :nl**2 if Euler are l- and m-dependent
Ci   bsite :magnetic field by site
Ci   nbf   :1 if bsite is l- and m-independent,
Ci         :nl if bsite is l-dependent and m-independent
Ci         :nl**2 if bsite is l- and m-dependent
Ci         :99 : flag that there exists no external field
Ci   mode  :0, use spin density matrix to make moments along T
Ci         :1, use qnus to make moments along qnu
Co Outputs
Cr Remarks
Cr   Definition of rho in terms of M: (standard definition of sigma)
Cr      rho = M . sigma/2
Cr   Pauli matrices sigma: (NB: normal def is twice this)
Cr
Cr              (0  1)             (0 -i)           (1  0)
Cr     sigmax = (    )    sigmay = (    )  sigmaz = (    )
Cr              (1  0)             (i  0)           (0 -1)
Cr
Cr   Given rho, M can be obtain from:
Cr     M_x =  2 Re(rho21) = Re (rho12+rho21)
Cr     M_y =  2 Im(rho21) = Im (rho21-rho12)
Cr            M_z =  (rho11)-(rho22)
Cr
Cr   Second (symmetrized) form is used because for numerical reasons,
Cr   rhos may not be quite hermitian; e.g. when rhos is generated
Cr   by a Green's function technique.
Cr
Cr   Double counting term is
Cr     Tr <(B.sigma)(rho)>
Cr
Cr   Input B.sigma is
Cr                1 (Bz   Bx-iBy)
Cr    B.sigma =   - (           )
Cr                2 (Bx+iBy  -Bz)
Cr
Cr   Then (CHECK wrong factors of 2 in both b and sigma)
Cr      Tr <(B.sigma)(rho)>
Cr      = 1/2 (Bz(rho11-rho22) + (Bx-iBy) rho21 + (Bx+iBy) rho12)
Cr      = Bz(rho11-rho22)/2 + (Bx-iBy)(Mx+iMy)/4 + (Bx+iBy)(Mx-iMy)/4
Cr      = Bz(rho11-rho22)/2 + Bx Mx / 2 + By My /2
Cr      = 1/2 (B . M)
Cr   This formula can be computed either with the moments qnu
Cr   or from the spin-density matrix.
Cu Updates
Cu   07 Apr 04 First created
C ----------------------------------------------------------------------
C     implicit none
C Passed variables
      integer nbas,nrhos,nl,nsp,ipc(nbas),mode,lihdim,indxsh(lihdim),
     .  neul,nbf
      parameter (nsp=2)
      double precision rhos(2,3,nrhos,2,2,1),eula(nbas,neul,3),bxc(3,*),
     .  qnu(3,0:nl-1,2,*),pp(6,0:nl-1,nsp,*),bsite(nbas,nbf,3),bsigr(2),
     .  sop(0:nl-1,nsp,nsp,9,*)
C Local variables
      logical lrhol,lbf
      integer i,j,ib,ic,lgunit,ipr,k,stdo,lmr,l,m,ilm,
     .  kb,kb2
      integer PRT1,PRT2
      double precision ammi(3),amli(3,3),qli(3),qm(3),ab,ddot,
     .  rotm(3,3),rotb(3,3),rotmb(3,3),alpha,beta,gamma,facB,
     .  bhat(3),sam(3),bxtbxc,pnu(nl,nsp),qnum(3,nsp),bdotm,bdotma,
     .  bdotml,ph,th,rhom(3),amloc(3,3),pi,bxcpp(3,2)
      double precision fac1,fac2
      character strn0*12,strn1*25,strn2*4,strn3*10
C     double precision bloc(3),hatm(3)
C     logical pass1,lprt
C     facB covers convention for B:
C     facB= 1 => +B induces positive -M
C     facB=-1 => +B induces positive +M
      parameter (facB=-1d0)
      parameter (PRT1=30,PRT2=50)

C --- Setup ---
      call dpzero(pnu,nl*nsp)
      pi = 4*datan(1d0)
      call getpr(ipr)
      lrhol = nrhos .eq. nl .or. mode .eq. 1
      lbf = nbf .ne. 99 .and. nbf .gt. 0
      stdo = lgunit(1)
      call isanrg(mode,0,1,' bsrhos','mode',.true.)
      kb2 = 1
      if (lbf) kb2 = 2
      bsigr(1) = 0
      bsigr(2) = 0

C --- Loop over bfield type: kb=1 for bxc, kb=2 for bext ---
      do  kb = 1, kb2

      if (ipr .ge. PRT2) then
        strn0 = 'B(exch-corr)'
        if (kb .eq. 2) strn0 = 'B (external)'

        if (mode.eq.0) strn1 = 'density matrix.'
        if (mode.eq.1) strn1 = 'sphere multipole moments.'

        strn2 =  ' '
        if (neul .eq. 1)     strn2 = '  l '
        if (neul .eq. nl)    strn2 = '  l '
        if (neul .eq. nl**2) strn2 = ' ilm'
        if (strn2 .eq. ' ') call rxi('BSRHOS: bad neul',neul)

        strn3 = ' '
        if (kb .eq. 2) strn3 = 'Bxc^.Bext^'
        write(stdo,333) strn0, strn1, strn2, strn3
  333   format(/' BSRHOS: double counting terms for ',a,
     .    '.  Spin density M'/
     .    ' (loc coordinates B || z) taken from ',a/
     .    '   ib',a,5x,'Mx',10x,'My',10x,'Mz',11x,'<M.B>',4x,a)
      endif

C --- Loop over sites ---
      lmr = 0
      do  ib = 1, nbas
        bdotma = 0
        call dpzero(sam,3)
        ic = ipc(ib)
C       Assume for now mag B = Bxc
        call dcopy(3,bxc(1,ic),1,bhat,1)

        ilm = 0
        do  l = 0, nl-1
          bdotml = 0
          k = l+1
          if (indxsh(lmr+1) .gt. lihdim) then
            lmr = lmr + 2*l+1
            ilm = ilm + 2*l+1
            goto 2
          endif

C     ... Matrix elements of Bxc for this channel
C         (If kb=2, overwrite later with M.E. to bext)
          do  i = 1, 3
            bxcpp(i,1) = sop(l,1,1,i+6,ic)
            bxcpp(i,2) = sop(l,2,2,i+6,ic)
          enddo

C     ... Make B.M for this channel
          call dpzero(qli,3)
          call dpzero(amli,3*3)
          do   m = -l, l
            ilm = ilm + 1
            lmr = lmr+1
            bxtbxc = 0

C           kb=2 : use Bext in place of Bxc
            if (kb .eq. 2) then
              k = 1
              if (nbf .eq. nl) k = l+1
              if (nbf .eq. nl*nl) k = ilm
              call dcopy(3,bxc(1,ic),1,bhat,1)
              bxtbxc = dsqrt(
     .          ddot(3,bsite(ib,k,1),nbas*nbf,bsite(ib,k,1),nbas*nbf))
              if (bxtbxc .ne. 0) then
                bxtbxc = ddot(3,bsite(ib,k,1),nbas*nbf,bhat,1)/bxtbxc
              endif
              bhat(1) = bsite(ib,k,1)
              bhat(2) = bsite(ib,k,2)
              bhat(3) = bsite(ib,k,3)
              ab = dsqrt(ddot(3,bhat,1,bhat,1))
C             bhat should be normalized; but it doesn't matter
C             since the bhat is only used for orientation
C             call dscal(3,1/ab,bhat,1)
C             Normalization defined s.t. d.c. = facB/2 Bext.M
              do  i = 1, 3
                bxcpp(i,1) = facB*ab/2*sop(l,1,1,i+3,ic)
                bxcpp(i,2) = facB*ab/2*sop(l,2,2,i+3,ic)
              enddo
            endif

C           Euler angles for magnetic moment this channel
            k = 1
            if (neul .eq. nl) k = l+1
            if (neul .eq. nl*nl) k = ilm
            alpha = eula(ib,k,1)
            beta  = eula(ib,k,2)
            gamma = eula(ib,k,3)

C       ... Moments of M in coordinate system where bhat.zhat=1
C           In this coordinate system, double-counting terms
C           involve only z component of M.
C           This block makes amloc(i,j) and amli(i,j) where
C           i=1,2,3 for 0th,1st,2nd moment; j=Cartesian coordinate.
C
C           The coordinate system for M for this orbital is specified by
C           Euler angles (which generate Rm).  Then M(loc)=Rm*M(global).
C           Let Rb = rotation of of Bxc to zhat.  (NB: |Bxc|=1).
C           Bxc.M = Rb*Bxc . Rb*M = zhat . (Rb*(Rm+)) M(loc)
C           M in global coordinates = Rm+ M(loc)
C           M in the 'Bxc' coordinate system is (Rb*(Rm+)) M(loc)
            call eua2rm(alpha,beta,gamma,rotm)
C           Rotation matrix Rb that rotates Bxc to zhat
            if (bhat(1)**2+bhat(2)**2 .ne. 0) then
              ph = datan2(bhat(2),bhat(1))
              th = datan2(dsqrt(bhat(1)**2+bhat(2)**2),bhat(3))
            else
              ph = 0
              th = 0
            endif

            call rotma(ph+pi/2,pi/2,th,rotb)
C           Check that B is properly rotated to zhat
C           call dgemm('N','N',3,1,3,1d0,rotm,3,bhat,3,0d0,bloc,3)
C           print *, bhat
C           print *, bloc
C           Rb * Rm+
            call dgemm('N','T',3,3,3,1d0,rotb,3,rotm,3,0d0,rotmb,3)

C           Alternative <B.M> when M spec'd by Euler angles
C           Rotate B to coordinate system of moment (Euler angles)
C           call eua2rm(alpha,beta,gamma,rotm)
C           hatm(1) = 0
C           hatm(2) = 0
C           hatm(3) = 1
C           call dgemm('N','N',3,1,3,1d0,rotm,3,bhat,3,0d0,bloc,3)
C           dotmb = ddot(3,bloc,1,hatm,1)

C           Alternative <B.M> when M spec'd by Euler angles
C           Same product in global coordinates
C           hatm(1) = dcos(alpha)*dsin(beta)
C           hatm(2) = dsin(alpha)*dsin(beta)
C           hatm(3) = dcos(beta)
C           Alternatively,
C           call eua2rm(alpha,beta,gamma,rotm)
C           hatm(1) = rotm(3,1)
C           hatm(2) = rotm(3,2)
C           hatm(3) = rotm(3,3)
C           dotmb = ddot(3,bhat,1,hatm,1)

            do  i = 1, 3

C             Case M(loc) specified by rhos
              if (mode .eq. 0) then
                k = l+1
                if (.not. lrhol) k = ilm
                rhom(1) = rhos(1,i,k,1,2,ic) + rhos(1,i,k,2,1,ic)
                rhom(2) = rhos(2,i,k,2,1,ic) - rhos(2,i,k,1,2,ic)
                rhom(3) = rhos(1,i,k,1,1,ic) - rhos(1,i,k,2,2,ic)
                qm(i)   = rhos(1,i,k,1,1,ic) + rhos(1,i,k,2,2,ic)
                if (lrhol) then
                  call dscal(3,1/dble(2*l+1),rhom,1)
                  qm(i) = qm(i)/(2*l+1)
                endif
C             Case M(loc) specified by qnu
              elseif (mode .eq. 1) then
                rhom(1) = 0
                rhom(2) = 0
                rhom(3) = (qnu(i,l,1,ic) - qnu(i,l,2,ic))/(2*l+1)
                qm(i)   = (qnu(i,l,1,ic) + qnu(i,l,2,ic))/(2*l+1)
              endif

              call dgemm('N','N',3,1,3,1d0,rotmb,3,rhom,3,0d0,ammi,3)
              qli(i) = qli(i) + qm(i)
              do  j = 1, 3
                amloc(i,j) = ammi(j)
                amli(i,j)  = amli(i,j) + ammi(j)
              enddo

            enddo

C       ... Compute double-counting, m-dependent Euler angles
            if (neul .eq. nl**2) then
              do  i = 1, 3
                qnum(i,1) = (qm(i) + amloc(i,3))/2
                qnum(i,2) = (qm(i) - amloc(i,3))/2
              enddo
              call asbxcm(pp(1,l,1,ic),pp(1,l,2,ic),
     .          qnum(1,1),qnum(1,2),bxcpp,fac1,fac2)
              bdotm = fac1-fac2
C             print *, '!!'
C             bsigr(kb) = bsigr(kb) + bdotm
              bdotml = bdotml + bdotm
              if (ipr .ge. PRT2) then
                if (bxtbxc .ne. 0)
     .            write(stdo,'(i5,i3,3f12.6,f14.7,f10.5)')
     .            ib,ilm,(amloc(1,j),j=1,3),fac1-fac2,bxtbxc
                if (bxtbxc .eq. 0) write(stdo,'(i5,i3,3f12.6,2f14.7)')
     .            ib,ilm,(amloc(1,j),j=1,3),fac1-fac2
              endif
            endif

C         End loop over m
          enddo

C     ... Double-counting, m-independent case
          do  i = 1, 3
            qnum(i,1) = (qli(i) + amli(i,3))/2
            qnum(i,2) = (qli(i) - amli(i,3))/2
          enddo
          call asbxcm(pp(1,l,1,ic),pp(1,l,2,ic),
     .      qnum(1,1),qnum(1,2),bxcpp,fac1,fac2)
          bdotm = fac1-fac2
C         Sanity check that m-summed d.c. = l- contr
C         Only satisfied when B or rho is m-independent
          if (neul .eq. nl**2) then
            if (nbf .le. nl) call
     .        fsanrg(bdotm,bdotml,bdotml,1d-9,'bsrhos','d.c.',.true.)
            bdotm = bdotml
          endif

          bsigr(kb) = bsigr(kb) + bdotm
          bdotma = bdotma + bdotm
          if (ipr .ge. PRT2 .and. neul .le. nl) then
            write(stdo,'(i5,i3,3f12.6,2f14.7)')
     .        ib,l,(amli(1,j),j=1,3),bdotm
          endif
          if (ipr .ge. PRT2 .and. neul .gt. nl) then
            write(stdo,'(''  l ='',i3,3f12.6,2f14.7)')
     .           l,(amli(1,j),j=1,3),bdotm
          endif

C         Add to net magnetization
          do  j = 1, 3
            sam(j) = sam(j) + amli(1,j)
          enddo

C          if (ipr .gt. PRT2 .or. .true.) then
C            write(stdo,'(''  l isp'',19x,''qnu'',37x,''qloc'')')
C            do  k = 1, 2
C              write(stdo,'(2i3,1x,3f13.7:1x,3f13.7)')
C     .          l,k,(qnu(i,l,k,ic),i=1,3),(qnum(i,k),i=1,3)
C            enddo
C            print *, '---------'
C          endif

C       End loop over l
        enddo
    2   continue

      if (ipr .ge. PRT2) then
        write(stdo,'(4x,''atom'',3f12.6,f14.7)')(amli(1,j),j=1,3),bdotma
      endif

C     End loop over sites
      enddo

      if (ipr .ge. PRT2) then
        write(stdo,'(5x,''sum'',3f12.6,f14.7)') sam,bsigr(kb)
      endif

C     End loop over field type
      enddo

      if (ipr .ge. PRT1) then
        write(stdo,345) (bsigr(kb), kb=1,kb2)
  345   format(' BSRHOS: <Bxc.sigma.rho>=',f12.7:
     .                '  <Bext.sigma.rho>=',f12.7)
      endif

      end
      subroutine asbxcm(pp1,pp2,qnu1,qnu2,bpp,fac1,fac2)
C- ASA Double-counting term <B.M> for one channel
C ----------------------------------------------------------------------
Ci Inputs
Ci   pp1   :potential parameters for spin 1
Ci   pp2   :potential parameters for spin 2
Ci   qnu1  :sphere moments for spin 1
Ci   qnu2  :sphere moments for spin 2
Ci   bpp   :matrix elements of orbitals and B-field (soprm)
Ci         : bxc(1:3) = <phi|B|phi> <phi|B|phidot> <phidot|B|phidot>
Co Outputs
Co   fac1  :m0<phi|B|phi> + m1<phi|B|dot> + m2<dot|B|dot> spin 1
Co         :where m0=q0-p*q2, m1=q1, m2=q2
Co   fac2  :m0<phi|B|phi> + m1<phi|B|dot> + m2<dot|B|dot> spin 2
Cl Local variables
Cl         :
Cr Remarks
Cr   pp are used to extract the true <phi|phi> portion qnu
Cu Updates
Cu   07 Apr 04 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      double precision pp1(5),pp2(5),qnu1(3),qnu2(3),bpp(3,2)
      double precision fac1,fac2
C ... Local parameters
      double precision ovl,pgam(2)
      double precision oalpha,ql(2)

C     p^gamma = amount of phidot**2 in qnu(1)
      ovl = oalpha(pp1(1),pp1(2),pp1(3)**2,pp1(6),pp1(5))
      pgam(1) = pp1(4)-ovl**2
      ovl = oalpha(pp2(1),pp2(2),pp2(3)**2,pp2(6),pp2(5))
      pgam(2) = pp2(4)-ovl**2
      ql(1) = qnu1(1) - pgam(1)*qnu1(3)
      ql(2) = qnu2(1) - pgam(2)*qnu2(3)

      fac1 = ql(1)*bpp(1,1) + qnu1(2)*bpp(2,1) + qnu1(3)*bpp(3,1)
      fac2 = ql(2)*bpp(1,2) + qnu2(2)*bpp(2,2) + qnu2(3)*bpp(3,2)

      end

