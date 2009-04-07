      subroutine magtrq(nbas,nl,nclass,ipc,sdmod,sdprm,ssite,sspec,
     .  dclabl,pp,rhos,nrhos,ehf,eula,neul,frc,aamom)
C- Magnetic torque and update of the Euler angles
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   nclass:number of inequivalent classes
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   sdmod :specifies kind of spin statics or dynamics; see Remarks
Ci         :1s digit
Ci           0: Euler angles updated by zeroing out off-diagonal
Ci              parts of the spin-density matrix
Ci           1: Euler angles updated by force * sdprm(2)
Ci           2: Simple spin dynamics w/; see remarks
Ci           3: Spin dynamics with friction terms
Ci         :10's digit affects mixing of eula (not used here)
Ci         :100's digit
Ci            1: Nose thermostat (see Remarks)
Ci         :1000's digit
Ci            1: Reserved by lm to prevent updating of atomic pp's
Ci   sdprm :parameters governing dynamics (see Remarks)
Ci         :(1) scale factor amplifying magnetic forces
Ci         :(2) time step tau for Landau dynamics
Ci         :(3) reference energy etot0 for Nose thermostat
Ci              (not used unless fscal is zero)
Ci         :(4) maximum allowed change in angle (not used unless nonneg)
Ci         :(5) etol: set tau=0 this iter if etot-ehf>etol
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec class relax
Ci     Stored:    *
Ci     Passed to: *
Ci   dclabl:class name, packed as a real number
Ci   pp    :potential parameters (atomsr.f)
Ci   rhos  : spin density matrix
Ci         : Note: rhos should be hermitian in spin space, but may not
Ci         : be owing to energy integration errors in the complex
Ci         : plane. Use a symmetrized form to minimize errors.
Ci   neul  :1 if Euler angles are l- and m-independent,
Ci         :nl if Euler are l-dependent and m-independent
Ci         :nl**2 if Euler are l- and m-dependent
Ci   ehf
Co Outputs
Co   eula  :updated, if sdmod<>3
Co   frc   :if sdmod=2
Co   aamom :local magnetic moments
Cl Local variables
Cl   ifrlx :suppress updating any Euler angles for any ifrlx nonzero
Cl   ila   :1 if Euler angles are l-dependent, otherwise 0
Cr Remarks
Cr   Spin statics: (sdmod=1):  Direction vector e (=z in loc. coord)
Cr     updated by fscal * J, where J = "force"
Cr   Spin dynamics: (sdmod=2):  Direction vector rotated by
Cr     z = tau*gyro/M (z cross J) + Nose thermostat
Cr     tau, M, gyro:  time step, moment amplitude, gyromagnetic ratio.
Cr     Nose thermostat accounts for coupling to the lattice, and adds
Cr     'statics' rotation, whose magnitude magtrq estimates to change
Cr     the internal energy to the target value etot0.  It estimates
Cr     the rotation angle from dE = 1/2 sum f_j th_j, and sets each
Cr     th_j in proportion to the f_j.  The th_j estimated by the above
Cr     is additionally scaled by sdprm(1).
Cr   Formulas for moments:
Cr            M_x =  2 Re(rho+-) = (rho+-)+(rho-+)
Cr            M_y =  2 Im(rho+-) = ((rho+-)-(rho-+))/i
Cr            M_z =  (rho++)-(rho--)
Cr   Second (symmetrized) form is used because for numerical reasons,
Cr   rhos may not be quite hermitian; e.g. when rhos is generated
Cr   by a Green's function technique.
Cb Bugs
Cr   No attempt to symmetrize updated Euler angles by class.
Cu Updates
Cu   21 Apr 04 Additions for an m-dependent spin-density matrix
Cu   04 Sep 03 Fix improper weighting of spin density matrix for
Cu             sdmod=0
Cu   27 Jan 03 Some alterations : sdmod is now consistent with docs
Cu   22 May 02 Bug fix to handle high and neglected orbitals
C ----------------------------------------------------------------------
C     implicit none
C Passed variables
      integer nbas,nl,nrhos,neul,nclass,ipc(*),sdmod
      double precision pp(6,nl,2,nclass),ssite(1),sspec(1),
     .  eula(nbas,neul,3),frc(nbas,3),aamom(nbas),ehf,
     .  sdprm(5),rhos(2,0:2,nrhos,2,2,nclass),dclabl(nclass)
C local variables
      logical nose,ltmp,lrhol
      integer ipr,ib,ic,m,lp1,lgunit,idm,ila,ipass,ifrlx(4),lsdmod,ml,
     .  mmax,stdo,is,ilm,kpass
      integer n0,nkap0
      parameter (nkap0=3,n0=10)
      integer idxdn(n0,nkap0)
      double precision amag(3),sumq,wt,swt,etot0,diffe,difold,htau,
     .  qm(64),amxm(64),amym(64),amzm(64),q,amx(0:2),amy(0:2),amz(0:2),
     .  fscal,rhox,rhoy,rhoz,coff0,coff1,vx,vy,fx,fy,sumJ,sumJ2,
     .  rotg(3,3)
      character*8 clabl
      character outs*80
      character hmod2(3)*20
      save difold
      data hmod2 /' ',', l-dep',', lm-dep'/ difold /0d0/

C     CALL PSHPR(55)

C --- Setup ---
      stdo = lgunit(1)
      lrhol = nrhos .eq. nl
      call getpr(ipr)
      if (neul .eq. nl*nl .and. lrhol) then
        call info0(10,0,0,
     .    'MAGTRQ: rhos is not m-dependent; nothing calculated')
        return
      endif
C ... Printout system moments, get total moment of unit cell
      call pshpr(min(ipr-30,29))
      call amagnc(nbas,nl,ipc,rhos,nrhos,0,eula,neul,0,amag,aamom)
      call poppr
      call query('sdmod',2,sdmod)
      lsdmod = iabs(sdmod)
C ... Spin statics or dynamics mode
      idm = mod(lsdmod,10)
      if (iabs(idm) .le. 1) call query('fscal',4,sdprm(1))
      if (idm .eq. 2) call query('tau',4,sdprm(2))
      call getpr(ipr)

      fscal = sdprm(1)
C ... ila=2,1,0 if Euler angles lm-, l-, or just site- dependent
      ila = 0
      if (neul .eq. nl) ila = 1
      if (neul .eq. nl**2) ila = 2
      if (nrhos .lt. neul)
     .  call rx('magtrq: Euler angles m-dependent but rhos is not')

C ... nose = .true. if Nose thermostat
      nose = mod(lsdmod/100,10) .ne. 0 .and. idm .eq. 2
      if (nose) fscal = 0
      etot0 = sdprm(3)

C ... Check for consistency in modes
      if (ila .gt. 0 .and. idm .eq. 3)
     .  call rxi('no spin-dependent angles in mode',idm)
      ipass = 0

C ... save tau into backup
      htau = sdprm(2)

C ... Printout
      outs = ' '
      if (idm .eq. 0) call info2(10,0,0,' MAGTRQ:  '//
     .  'mode %i (align P with density-matrix)  fscal %d',sdmod,fscal)
      if (idm .eq. 1) then
        if (sdprm(4) .eq. 0) call info2(10,0,0,' MAGTRQ:  mode %i '//
     .    '(follow forces'//hmod2(ila+1)//'%a)  fscal %d',sdmod,fscal)
        if (sdprm(4) .ne. 0) call info2(10,0,0,' MAGTRQ:  mode %i '//
     .    '(follow forces'//hmod2(ila+1)//
     .    '%a)  fscal %d max d(theta)=%d',fscal,sdprm(4))
      endif
      if (ipr .ge. 10 .and. idm .eq. 2) then
        call awrit1(' MAGTRQ:  mode %i (simple Landau dynamics'//
     .    hmod2(ila+1),outs,80,0,sdmod)
        if (nose) call awrit1('%a Nose etot0 %d',outs,80,0,etot0)
        call awrit1('%a)  tau %d',outs,80,-lgunit(1),sdprm(2))
      endif
      if (idm .eq. 3) call info2(10,0,0,' MAGTRQ:  mode %i '//
     .  '(%?;n;l-dep ;%j;spin dynamics)',sdmod,ila)

C ... Check that etot0-ehf<tol
      if (idm .eq. 2 .and. sdprm(2) .gt. 0 .and. sdprm(5) .gt. 0
     .  .and. dabs(etot0-ehf) .gt. sdprm(5)) then
        if (ipr .ge. 10)
     .    call awrit2(' set tau=0 this iter:  etot0-ehf=%d exceeds '//
     .    'etol=%d',outs,80,-lgunit(1),dabs(etot0-ehf),sdprm(5))
        sdprm(2) = 0
      endif

      if (ipr.ge.30 .and. ipr.lt.45 .and. idm.ne.3) write(stdo,335)
C  335 format(10x,'local phi,d(theta)',14x,'New Euler angles')
  335 format(11x,'local phi,d(theta)',6x,'alpha     beta',10x,
     .  'changed to')


C --- Entry point for second pass ---
    5 continue
      ipass = ipass+1
      sumJ  = 0
      sumJ2 = 0
      do  10  ib = 1, nbas
      call upack('site spec class relax',ssite,ib,is,ic,ifrlx,0)
      call r8tos8(dclabl(ic),clabl)
      call upack2('spec idxdn',sspec,is,idxdn)
      sumq = 0
      swt  = 0
      rhox = 0
      rhoy = 0
      rhoz = 0
      fx   = 0
      fy   = 0
      if (ipr .ge. 45) write(stdo,333)
  333 format(/' Class mom l      Q       rhox      rhoy      rhoz',
     .  '        dC        Vx        Vy')

C --- Accumulate force for this atom from each l and moment ---
      do  11  kpass = 1, 2
      if (ipr.ge.45 .and. idm.ne.3 .and. ila.gt.0 .and. kpass.eq.2)
     .  write(stdo,335)
      do  12  m = 0, 2
      do  14  lp1 = 1, nl
      if (idxdn(lp1,1) .le. 2) then

C     coff0 = (C+)(sqrdel-/sqrdel+) - (C-)(sqrdel+/sqrdel-)
      coff0      = pp(2,lp1,1,ic)*(pp(3,lp1,2,ic)/pp(3,lp1,1,ic)) -
     .             pp(2,lp1,2,ic)*(pp(3,lp1,1,ic)/pp(3,lp1,2,ic))
      coff1      = (pp(3,lp1,2,ic)/pp(3,lp1,1,ic)) -
     .             (pp(3,lp1,1,ic)/pp(3,lp1,2,ic))
      q = 0
      amz(m) = 0
      amy(m) = 0
      amx(m) = 0
      mmax = 0
      ilm = lp1-1
      if (.not. lrhol) then
        mmax = lp1-1
        ilm = (lp1-1)**2
      endif
      do  15  ml = -mmax, mmax
        ilm = ilm+1
        qm(ilm)   = rhos(1,m,ilm,1,1,ic) + rhos(1,m,ilm,2,2,ic)
        amzm(ilm) = rhos(1,m,ilm,1,1,ic) - rhos(1,m,ilm,2,2,ic)
        amym(ilm) = rhos(2,m,ilm,2,1,ic) - rhos(2,m,ilm,1,2,ic)
        amxm(ilm) = rhos(1,m,ilm,1,2,ic) + rhos(1,m,ilm,2,1,ic)
C       See note about bug, above
C       amym(m,ilm) = -amym(m,ilm)
C       These are the m-summed versions
        q = q + qm(ilm)
        amz(m) = amz(m) + amzm(ilm)
        amy(m) = amy(m) + amym(ilm)
        amx(m) = amx(m) + amxm(ilm)
   15 continue
      
C     Torques on the Euler angles, res. by lm, l, or just site
      mmax = 0
      ilm = lp1-1
      if (ila .eq. 2) then
        mmax = lp1-1
        ilm = (lp1-1)**2
      endif
      do  16  ml = -mmax, mmax
        ilm = ilm+1
        if (ila .eq. 2) then
          amz(m) = amzm(ilm)
          amy(m) = amym(ilm)
          amx(m) = amxm(ilm)
        endif

C   ... Forces, in local coordinate system
        vx = coff0*amx(m)
        vy = coff0*amy(m)

C   ... zeroth moment: make sum_l fx,fy,rhox,rhoy,rhoz
        if (m .eq. 0) then
C         l- or lm- dependent torque
          if (ila .ge. 1 .and. kpass .eq. 2) then
            call xxmgtq(lp1-1,idm,vx,vy,amx(0),amy(0),
     .        amz(0),sdprm,fscal,clabl,ifrlx(lp1),
     .        eula(ib,ilm,1),eula(ib,ilm,2),eula(ib,ilm,3),ipr)
          endif
          if (kpass .eq. 1) then
            wt = 1
            swt = 1
            sumq = sumq + q
            rhox = rhox + wt*amx(0)
            rhoy = rhoy + wt*amy(0)
            rhoz = rhoz + wt*amz(0)
            fx   = fx   + vx
            fy   = fy   + vy
            sumJ = sumJ + dsqrt(vx**2+vy**2)
            sumJ2= sumJ2+       vx**2+vy**2
          endif
        endif
        outs = clabl
        if (m .ne. 0 .or. lp1 .ne. 1) outs = ' '
        ltmp = ipr.ge.50 .or. ipr.ge.45 .and. m.eq.0
        if (ltmp .and. kpass.eq.1) write(stdo,334) outs(1:4), m,
     .    lp1-1, q,amx(m),amy(m),amz(m), coff0, vx, vy
  334   format(2x,a4,2i3,7f10.6)
        if (ltmp .and. kpass.eq.1 .and. lp1.eq.nl .and. ml.eq.mmax) then
          write(stdo,336)
     .      sumq, rhoz/swt, rhox/swt, rhoy/swt, fx, fy
  336     format('      sphere',4f10.6,10x,2f10.6)
        endif
   16 continue
      endif
   14 continue
   12 continue
   11 continue
      rhox = rhox / swt
      rhoy = rhoy / swt
      rhoz = rhoz / swt

C --- Update Euler angles for this atom depending on mode ---
      if (ipr .ge. 45 .and. idm .ne. 3 .and. ila .eq. 0) write(stdo,335)

      if (ila .eq. 0 .and. idm .ne. 3) then
        call xxmgtq(-1,idm,fx,fy,rhox,rhoy,rhoz,sdprm,fscal,clabl,
     .    ifrlx(1),eula(ib,1,1),eula(ib,1,2),eula(ib,1,3),ipr)

C ... Spin-dynamics mode
      elseif (ila .eq. 0 .and. idm .eq. 3) then
C   ... Rotation matrix corresponding to current Euler angles
        call eua2rm(eula(ib,1,1),eula(ib,1,2),eula(ib,1,3),rotg)
C   ... Forces in glob. xyz coords: 1/2 reproduces true dE/d(phi,th)
C       Different signs because  fx <-> phi-hat but fy <-> -theta-hat
        do  18  m = 1, 3
   18   frc(ib,m) = rotg(1,m)*fx/2 - rotg(2,m)*fy/2
        write(stdo,357) ib,fx,fy,(frc(ib,m),m=1,3)
  357   format(i4,2f12.6,1x,3f12.6)
      elseif (idm .eq. 3) then
        call rx('MAGTRQ: no SD for l-dependent Euler now')
      endif
   10 continue

c      call prmx('f',frc,nbas,nbas,3)
      if (idm .eq. 3) return

C --- Nose thermostat: make second pass relaxing along J ---
      if (nose .and. ipass .eq. 1 .and. idm .gt. 1) then
        diffe = etot0-ehf
C ...   Rotate each angle in proportion to force.
C       fscal from estimated desired change in energy:
        fscal = -2*diffe/sumJ2
        idm = 1
        if (ipr .gt. 10) then
          write(stdo,'(1x)')
          if (difold .eq. 0)
     .      call awrit3(' Nose etot0-ehf=%d  sumJ,J2=%d %d',
     .      ' ',80,lgunit(1),diffe,sumJ,sumJ2)
          if (difold .ne. 0) call awrit4(
     .      ' Nose etot0-ehf=%,6;6d  old=%,6;6d  sum J,J2=%d %d',
     .      ' ',80,lgunit(1),diffe,difold,sumJ,sumJ2)
          call awrit3(' est fscal = %d  scaled by %d  fscal = %d',
     .      ' ',80,lgunit(1),fscal,sdprm(1),fscal*sdprm(1))
        endif
        fscal = fscal*sdprm(1)
        difold = diffe
        goto 5
      endif
      sdprm(2) = htau

C      print *, '!!'
C      CALL POPPR
C     call wkchk('exit magtrq')
C     call rx('done')
      end
      subroutine xxmgtq(l,idm,fx,fy,rhox,rhoy,rhoz,sdprm,fscal,clabl,
     .  ifrlx,alpha,beta,gamma,ipr)
C- Update Euler angles from forces according to mode
C ----------------------------------------------------------------------
Ci Inputs
Ci   l     :l quantum number (used for printout only)
Ci   idm   :mode:
Ci           0: Euler angles corresponding to zero of off-diagonal rho
Ci           1: Euler angles updated by force
Ci           2: Simple spin dynamics
Ci           3: Spin dynamics with friction terms
Ci   fx    :force in x, local coordinate
Ci   fy    :force in y, local coordinate
Ci   rhox  :x-component of s.d. matrix, local coordinates
Ci   rhoy  :y-component of s.d. matrix, local coordinates
Ci   rhoz  :z-component of s.d. matrix, local coordinates
Ci   sdprm :spin-dynamics parameters
Ci   fscal :scale factor amplifying rotation matrix
Ci   clabl :class name
Ci   ifrlx :0 suppress updating alpha,beta,gamma
Ci         :1 update alpha,beta,gamma
Ci   ipr   :print verbosity
Cio Inputs/Outputs
Cio   alpha :On input, 1st Euler angle defining local coordinate
Cio         :On output, this Euler angle is updated according to mode
Cio   beta  :On input, 2nd Euler angle defining local coordinate
Cio         :On output, this Euler angle is updated according to mode
Cio   gamma :On input, 3rd Euler angle defining local coordinate
Cio         :On output, this Euler angle is updated according to mode
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   27 Jan 03
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer l,idm,ipr,ifrlx
      character*(*) clabl
      double precision fx,fy,rhox,rhoy,rhoz,sdprm(4),fscal,
     .  alpha,beta,gamma
C ... Local parameters
      integer lgunit
      double precision alphan,betan,gamman,ph,th,thx,
     .  rotg(3,3),rotm(3,3),rotl(3,3),pi,gyro,tau
      character outs*80,c*4

C --- Setup ---
      pi = 4*datan(1d0)
      alphan = alpha
      betan  = beta
      gamman = gamma
      tau    = sdprm(2)
      gyro   = 2
      ph     = 0
      th     = 0

C --- (idm=0) Align p, T ---
C ... th,ph in local coordinates corresponding to rhox,rhoy=0
C     Check sign of rhoz to allow local mom<0 restore to equil.
      if (idm .eq. 0) then
        if (rhox**2+rhoy**2 .ne. 0) then

C         old
C         ph = datan2(rhox,-rhoy)
C         th = datan(dsqrt(rhox**2+rhoy**2)/abs(rhoz))

C         Revised 19 Apr 04
          ph = datan2(rhoy,rhox)
          if (rhoz .gt. 0) then
            th = datan2(dsqrt(rhox**2+rhoy**2),rhoz)
          else
            th = datan2(-dsqrt(rhox**2+rhoy**2),-rhoz)
          endif
C         Debugging: confirm that rotl * rho points along z
C         call rotma(ph+pi/2,pi/2,th,rotl)
C         call rm2eua(rotl,alphan,betan,gamman)
C         rotg(1,1) = rhox
C         rotg(2,1) = rhoy
C         rotg(3,1) = rhoz
C         call dgemm('N','N',3,1,3,1d0,rotl,3,rotg,3,0d0,rotg(1,2),3)
C         print '(3f12.6)', rotg(1:3,2)
C         goto 10

        endif

C --- (idm=1) Move e along J = force ---
C ... th,ph in local coordinates corresponding to (fx,fy)
      elseif (idm .eq. 1) then

C       old
C        if (fx .ne. 0 .or. fy .ne. 0) ph = datan2(fx,-fy)
C        th = -datan(dsqrt(fx**2+fy**2))

C       Debugging: compare to new
C        call rotma(ph,pi/2,th,rotl)
C        call rm2eua(rotl,alphan,betan,gamman)
C        print *, alphan,betan

C       Revised 19 Apr 04
        if (fx .ne. 0 .or. fy .ne. 0) ph = datan2(fy,fx)
        th = -datan(dsqrt(fx**2+fy**2))
C       Debugging: compare to old
C       call rotma(ph+pi/2,pi/2,th,rotl)
C       call rm2eua(rotl,alphan,betan,gamman)
C       print *, alphan,betan
C       pause
C
      endif

C ... Finish idm=0,1 : make local rotation matrix
      if (idm .eq. 0 .or. idm .eq. 1) then
C       Upper bound on allowed angle change
        thx = th
        c = ' '
        if (sdprm(4) .gt. 0 .and. fscal .ne. 0) then
          thx = -min(-th,sdprm(4)/dabs(fscal))
          if (thx .ne. th) c = ' (*)'
        endif
C   ... matrix for rotation by fscal*th about ph in xy plane
C       use this for old conventions
C       call rotma(ph,pi/2,fscal*thx,rotl)
C       Revised 19 Apr 04
        ph = ph+pi/2
        call rotma(ph,pi/2,fscal*thx,rotl)
C       print *, 'ph,th=',sngl(ph),sngl(th)
C       call prmx('rotl',rotl,3,3,3)
C       call rm2eua(rotl,alphan,betan,gamman)
C       print *, alphan,betan
C

C --- (idm=2) Move angles along tgm * z x J + fscal * J ---
      elseif (idm .eq. 2) then
        call rxi('update this SD branch, mode',idm)

C  ...  Landau term, store for now in rotg (NB: sign not checked)
        ph = ph+pi/2
        call rotma(ph,pi/2,tau*gyro/rhoz*th,rotg)
C  ...  Friction term, store temporarily rotm
        call rotma(ph-pi/2,pi/2,fscal*th,rotm)
C  ...  Net rotation is accumulation of Landau + friction terms
        call dmpy(rotm,3,1,rotg,3,1,rotl,3,1,3,3,3)

      else
        call rxi('MAGTRQ: no mode %i implemented',idm)
      endif

C --- New Euler angles from local rotation matrix rotl ---
C     The local coordinate system is specified by Euler angles which
C     generate Rg below.  In the local coordinate system the input
C     density points along z^: thus Rg*M^_in = z^ => M^_in(loc) = z^.
C     In these local coordinates, M^_out(loc) points in some other
C     direction specified by Rl, which is a matrix rotating M^_out(loc)
C     to z^.  Thus Rl M^_out(loc) = z^.  This block constructs a new
C     Rm that rotatates M^_out to z^.
C       M = Rg+ M(loc)  for any M
C       M^_out = Rg+ M^_out(loc) =  Rg+ Rl+ z^
C     Then
C       (Rl Rg) M^_out = z^
C     so the coordinate system that rotates M^_out to z^ is
C       Rm = (Rl Rg)
C     Rg = rotation matrix corresponding to current Euler angles
      call eua2rm(alphan,betan,gamman,rotg)
C     call eua2rm(0d0,0d0,0d0,rotg)
C ... Rm = rotation matrix to rotate M^_out to z^
      call dmpy(rotl,3,1,rotg,3,1,rotm,3,1,3,3,3)
C ... Euler angles for Rm
      call rm2eua(rotm,alphan,betan,gamman)
      call euasmo(.8d0,alpha,alphan)
      call euasmo(.8d0,gamma,gamman)
C      outs = ' '//clabl
C      if (l .gt. 0) outs = ' '
C      if (ifrlx .eq. 0) c = ' (f)'
C      call awrit5('%a%6p%;12,8D%;12,8D  %;12,8D%;12,8D%;12,8D'//c,outs,
C     .  80,0,ph,th,alphan,betan,gamman)
C      if (l .ge. 0) call awrit1('%a l=%i',outs,80,0,l)
C      if (ipr .ge. 30) call awrit0('%a',outs,-80,-lgunit(1))
C      call awrit0('%5o%0pEuler%a',outs,-80,-lgunit(2))
C      if (ipr .ge. 40) call awrit3('%21pold Euler  %;12,8D%;12,8D'//
C     .  '%;12,8D',outs,80,lgunit(1),alpha,beta,gamma)
      outs = '  '//clabl
      if (l .gt. 0) outs = ' '
      if (ifrlx .eq. 0) c = ' (f)'
      call awrit6('%a%7p%;10,6D%;12,8D%;12,6D%;10,6D%;12,6D%;10,6D'//
     .  c,outs,80,0,ph,th,alpha,beta,alphan,betan)
      if (l .ge. 0) call awrit1('%a l=%i',outs,80,0,l)
      if (ipr .ge. 30) call awrit0('%a',outs,-80,-lgunit(1))
      call awrit0('%5o%0pEuler%a',outs,-80,-lgunit(2))

C ... Update Euler angles if ifrlx is nonzero
      if (ifrlx .eq. 0) return
      alpha = alphan
      beta  = betan
C     gamma = gamman
      end

      subroutine euasmo(tol,aold,anew)
C- Rotate new angle by +/- 2 pi to make it as close as possible to old.
C     implicit none
      double precision aold,anew,tol,pi2

      pi2 = 8*datan(1d0)
      if (dabs((aold-anew)/pi2) .lt. tol) return
      if (anew .gt. aold) then
        anew = anew - pi2
      else
        anew = anew + pi2
      endif
      end

