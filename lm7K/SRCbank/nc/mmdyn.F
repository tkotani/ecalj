      subroutine mmdyn(smove,nbas,xsi,eula,qspirl,sdir,ebil,ibil,passf,
     .  iax,nttab,wkint,frc,amag,aamom,etot,tnow,nvario,ifs)
C- Micromagnetics dynamics
C ----------------------------------------------------------------------
Ci Inputs
Ci   smove :struct pertaining to dynamics
Ci     Elts read: nmodt modt ct kt ts tstot tsequ gyro ts0 tol prmint
Ci     Passed to: mmstp
Ci   nbas  :size of basis
Ci   frc   :corresponding magnetic forces (generated internally for
Ci         :empirical hamiltonian).
Ci   xsi   :starting global daemon friction parameters
Ci   qspirl:Parameters defining spin-spiral (not used at present)
Ci   sdir  :work array of dimension 3*nbas+3 holding dir vecs of eula
Ci   ebil  :coefficients to pair (heisenberg) hamiltonian
Ci   ibil  :table of which rule was used to make ebil
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   nttab :total number of pairs in neighbor and iax (pairc.f)
Ci   wkint :work array holding internal integration parameters.
Ci         :Bulirsch-Stoer integration: dim (neq*(3+mx)+5)) (see dbsstp)
Ci   aamom :absolute magnitude of magnetic moments on each atom
Ci   nvario:number of variables to include in save file
Ci   ifs   :save file logical unit
Ci   etot  :total energy (passed when forces computed externally)
Cio Inputs/Outputs
Cio   eula :On input, Euler angles for current time. If forces are
Cio        :calculated externally, they should correspond to these eula.
Cio        :On output, Euler angles for a new time
cio   passf:parameter characterizing present stage of integration
Cio        :On input,
Cio        :  -1, magnetic forces calculated internally.
Cio        :      mmdyn only returns when dynamics are completed
Cio        :      or aborted.
Cio        : >=0, magnetic forces calculated externally
Cio        :   0, first call to mmdyn: some initialization performed
Cio        :   1, start of a new time step
Cio        :  >1, continuation of current time step
Cio        :On output,
Cio        :  -1, dynamics are completed
Cio        :   1, start of a new time step
Cio        :  >1, continuation of current time step
Cio   tnow :current time since start of simulation
Cio        :This is the integration variable the integrator uses
Cio        :as it proceeds through the integration step.
Cio        :Euler angles should correspond to tnow
Co Outputs
Co   amag  :system magnetic moment
Co   etot  :total energy computed when forces calculated internally
Cl Local variables
Cl   ir    :parameter characterizing present stage of integration
Cl         :It is set by the numerical integrator.
Cl         :   0, start of a new integration step
Cl         :  >0, continuation of current integration step
Cl   ts    :suggest time interval for numerical integration to next point.
Cl         :ts may change after the integrator has completed a time step.
Cr Remarks
Cr   mmham integrates the Landau-Lifshitz equations of motion
Cr   with global daemons.  It operates in one of two modes:
Cr
Cr 1.  An empirical (heisenberg) hamiltonian H is supplied in the
Cr     through array ebil.  The coefficients J are normalized in
Cr     a way that the total energy is
Cr         E = 1/2 sum_ij J_ij e_i e_j
Cr     where e_j are the (unit) direction vectors.
Cr     In this mode mmdyn computes the forces internally and
Cr     integrates the equations of motion over the specified
Cr     interval (see struc smove).
Cr
Cr 2.  The magnetic forces are supplied by the calling program.
Cr     In this mode mmdyn carries out one step of the numerical
Cr     integration of equations of motion.  Given an time tnow and
Cr     Euler angles and forces corresponding to that time, mmdyn
Cr     estimates the Euler angles for a new time tnow by integrating
Cr     the equations of motion.  Caller generates the forces
Cr     corresponding to the updated tnow; this cycle is repeated
Cr     until the total elaspsed time is reached.
Cf Files:
Cf  *Writes time,energy,moment to
Cf   save file on completion of an integration step.
Cf  *Writes euler angles to file eula-sv at the start of 
Cf   an integration step.
Cu Updates
Cu   29 May 03 Change save file output
Cu   17 Aug 01 adapted from routine mm_dyn from ASA version 5.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,nvario,nttab,niax,passf,ifs,ibil(nttab),ir
      parameter (niax=10)
      integer iax(niax,nttab)
      double precision xsi(3),eula(nbas,3),tnow,
     .  sdir(nbas,3),qspirl(4),etot,smove(1),ebil(nttab),wkint(1),
     .  frc(nbas,3),amag(3),aamom(nbas)
C ... Local variables
      integer iprint,i,it,lgunit,fopna,j,ipr,onrul,opcorr,mxprul,
     .  ifi,neq,modt(3),iprt,nmodt
C ... For micromagnetics hamiltonian
      double precision tmom,kt,ts,ttot,tequ,ct(3),etscl,
     .  pi,tpi,vvars(4),ehf
      equivalence (tmom,vvars(2)),(ehf,vvars(3))
      character*80 outs
C ... heap:
      integer w(1)
      common /w/ w

C --- Setup repeated every call ---
      call upack('move nmodt modt ct',smove,nmodt,modt,ct,0,0)
      call upack('move kt ts tstot tsequ',smove,kt,ts,ttot,tequ,0)
      ts = min(ts,ttot-tnow)
      neq = nbas*2
      pi = 4*datan(1d0)
      vvars(1) = tnow
C ... Add nmodt for xsidot if thermostat term
      if (nmodt .gt. 0) then
        neq = neq+nmodt
        do  5  j = 1, nmodt
    5   eula(nbas+j,2) = xsi(j)
      endif
C ... Determine max pair rule for pair correlation function
      if (passf .lt. 0) then
        mxprul = 0
        do  6  j = 1, nttab
    6   mxprul = max(mxprul,ibil(j))
        call defi(onrul,-(mxprul+1))
        call defr(opcorr,(mxprul+1))
        do  7  j = 1, nttab
    7   w(onrul+ibil(j)) = w(onrul+ibil(j)) + 1
      endif

C --- Pick up where we left off if passf>0 ---
      call pshpr(iprint())
      if (passf .gt. 0) goto 16
C     We don't to rewind save file each iteration; so just do it once.
      call poseof(ifs)

C ... Printout
      call togpr
      if (iprint() .ge. 10) then
        call awrit7(' mmdyn:  simulation to'
     .  //' t=%;4d  kT=%;4d%?#n#  modet =%n:1i  ct =%n:1;4d#%2j#',
     .  outs,80,-lgunit(2),ttot,kt,nmodt,nmodt,modt,nmodt,ct)
        print *
        call awrit0('%a',outs,80,-lgunit(1))
      endif
      call togpr

C --- Re-entry for new time step ---
   10 continue
      if (iprint() .ge. 20)
     .  call awrit5('%x mmdyn:  new iteration:  tnow=%;4d  t-step=%;4d'
     .  //'%?#n#  xsi =%n:1;4,4d#%j#',outs,80,-lgunit(2),tnow,ts,
     .  nmodt,nmodt,eula(nbas+1,2))
      if (iprint() .ge. 20) then
        print *
        call awrit0('%a',outs,80,-lgunit(1))
      endif

C --- Re-entry for continued integration of this time step ---
      ir = 0
   20 continue
        call setpr(iprt(min(ir+2,3)))

C   --- Convert Euler angles into vectors ---
        do  110  i = 1, nbas
          sdir(i,1) = dcos(eula(i,1))*dsin(eula(i,2))
          sdir(i,2) = dsin(eula(i,1))*dsin(eula(i,2))
          sdir(i,3) = dcos(eula(i,2))
  110   continue
C   ... Rationalize Euler angles from vector
C        do  111  i = 1, nbas
C          eula(i,1) = datan2(sdir(i,2),sdir(i,1))
C          eula(i,2) = dacos(sdir(i,3))
C          eula(i,3) = 0
C  111   continue
        if (nmodt .gt. 0) then
          do  112  j = 1, nmodt
  112     sdir(nbas+j,3) = eula(nbas+j,2)
        endif

C   --- Calculate the forces (empirical hamiltonian only) ---
        if (passf .lt. 0) then
          etot = 0
          call dpzero(frc,nbas*3+3)
          call dpzero(w(opcorr),mxprul+1)
          call mmpair(nbas,nttab,iax,ebil,ibil,aamom,sdir,
     .      etot,amag,w(opcorr),frc)
          do  119  i = 1, mxprul+1
            if (w(onrul+i-1) .ne. 0) then
              call dpscop(w(opcorr),w(opcorr),1,i,i,1d0/w(onrul+i-1))
            endif
  119     continue
        endif
        tmom = dsqrt(amag(1)**2 + amag(2)**2 + amag(3)**2)

C   --- Printout ---
        ipr = iprint()
        if (ipr .gt. 30) then
          j = lgunit(2)
          print *
          call awrit2('%% Euler angles rows %i cols %i mode 0',
     .      ' ',80,j,nbas,3)
          call awrit3('%x mmdyn : starting Euler angles'//
     .      '%?#n#  xsi =%n:1;4,4d#%j#:',outs,80,-lgunit(1),nmodt,
     .      nmodt,sdir(nbas+1,3))
          if (ipr .gt. 30) print 344
  344     format('  ib   alpha    beta',7x,'sx',7x,
     .      'sy',7x,'sz',8x,'fx',7x,'fy',7x,'fz')
          do  118  i = 1, nbas
            write(j,346) (eula(i,it), it=1,2)
            if (ipr .gt. 30) print 345, i, (eula(i,it), it=1,2),
     .        (sdir(i,it), it=1,3), (frc(i,it), it=1,3)
  118     continue
          print *
  345     format(i4,2f9.5,1x,3f9.5,1x,3f9.5)
  346     format(2f12.6,' 0')
        endif

C   --- Case start of a new iteration  ---
        if (ir .eq. 0) then
C     ... Printout
          etscl = 1000/dble(nbas)
          etscl = 1
C         call u_sdyn(smove,xx,xx,tnow,xx,xx,xx)
          ehf = etot*etscl
          vvars(1) = tnow
C         call iosave('mxxf','time,mmom,ehf',vvars,-ifs,nvario)
          write(ifs,366) vvars(1),vvars(2),vvars(3)
  366     format(f14.4,f12.6,f12.6)
          if (ipr .ge. 10) then
            call awrit4(
     .        ' time = %;4g  mag =%3:1;4d |mag| = %;4d '//
     .        'etot = %;6,6d',outs,80,-lgunit(1),tnow,amag,tmom,etot*
     .        etscl)
            call awrit0('%a',outs,-80,-lgunit(2))
            if (passf .lt. 0) call
     .      awrit2('<s.s> %n:1;4d',' ',200,lgunit(2),mxprul+1,w(opcorr))
          endif
C     ... Save euler angles to save file
          ifi = fopna('eula-sv',-1,0)
          rewind ifi
          call ioeula(nbas,1,eula,1,sdir(nbas+1,3),-ifi)
          call fclose(ifi)
C     ... Exit when time exceeds ttot
          call query('ttot',4,ttot)
          if (tnow .gt. ttot+1d-8 .or. dabs(ts) .lt. 1d-8) then
            call poppr
            passf = -1
            return
          endif
        endif

C   --- Project out components of force along spin direction ---
c       call mm_prj(nbas,sdir,iprint(),w(ofrc),etot,amag,tmom)

C   --- Do another piece of this integration step ---
        call mmstp(nbas,neq,smove,sdir,frc,aamom,ir,tnow,ts,wkint,eula)
        call pack1('move ts',smove,ts)
C   ... Exit to caller when passf>=0
        if (passf .ge. 0) then
          do  15  j = 1, nmodt
   15     xsi(j) = eula(nbas+j,2)
          do  17  i = 1, nbas
            sdir(i,1) = dcos(eula(i,1))*dsin(eula(i,2))
            sdir(i,2) = dsin(eula(i,1))*dsin(eula(i,2))
            sdir(i,3) = dcos(eula(i,2))
   17     continue
          passf = ir+1
          call poppr
          return
        endif
C   ... Re-entry point for passf>0 (forces generated externally)
   16   continue
        if (ir .gt. 0) goto 20

C   --- End of time step:  update ts for new step ---
        ts = min(ts,ttot-tnow)
        call pack1('move ts',smove,ts)

C   ... Rationalize Euler angles
        continue
        tpi = 2*pi
        do  121  i = 1, nbas
          if (eula(i,2) .lt. 0d0) eula(i,2) = eula(i,2) + tpi
          if (eula(i,2) .gt. pi) then
            eula(i,2) = tpi - eula(i,2)
            eula(i,1) = eula(i,1) - pi
          endif
          if (eula(i,1) .gt. pi)  eula(i,1) = eula(i,1) - tpi
          if (eula(i,1) .lt. -pi) eula(i,1) = eula(i,1) + tpi
  121   continue
      goto 10

      end

      subroutine mmstp(nbas,neq,smove,sdir,frc,aamom,ir,tnow,ts,wkint,
     .  eula)
C- Integrate Euler angles part of one step, given magnetic forces
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   neq   :nbas*2 + nxsi where nxsi is the number of global daemons
Ci   smove :struct pertaining to dynamics (subroutine umove)
Ci     Elts read: kt ts gyro nmodt modt ct prmint ts0 tol
Ci   sdir  :unit magnetic direction vectors
Ci   aamom :local magnetic moments scaling the Heisenberg hamiltonian
Ci   ir    :internal integration parameter keeping track of current
Ci         :stage of integration step
Ci   wkint  work array holding internal integration parameters.
Ci         Bulirsch-Stoer integration: dim (neq*(3+mx)+5)) (dbsstp)
Cio Inputs/Outputs
Cio  ts    :On input,  suggested time step
Cio        :On output, time step actually taken
Cio  tnow  :On input,  time corresponding to eula,frc
Cio        :On output, integrator will set to a new value for which it
Cio        :           wants an updated eula,frc
Cio        :This is the independent variable the integrator assigns
Cio        :as it proceeds through the integration step.
Cio        :eula and frc should correspond to tnow
Cio  eula  :On input,  current values of Euler angles
Cio        :On output, updated values by integrator
Cio  frc   :On input,  magnetic forces
Cio        :On output, e-dot wrt phi,theta
Co Outputs
Cb Bugs:
Cb   xsi tacked on to end of eula; destroying any eula(1..3,3)
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,neq,ir
      double precision smove(1),gyro,ts,wkint(1),tnow,
     .  eula(nbas,3),sdir(nbas,3),frc(nbas,3),aamom(nbas)
C ... Local parameters
      integer oyerr,ofeu,nseq(0:15),mode,modt(3),mx,iprint,ode,nmodt,
     .  igets
      double precision ts0,tol,kt,ct(3),prmint(20)
      equivalence (ts0,prmint(3)),(tol,prmint(4))
C heap:
      integer w(1)
      common /w/ w

C ... Unpack dynamics and bs parameters
C     call u_sdyn(smove,!kt,!ts,!tnow,!tstot,!tsequ,!gyro)
      call upack('move kt gyro',smove,kt,gyro,0,0,0)
C     nseq(0) = 0
C     call u_bsi(smove,mode!,nmodt!,modt!,ct!,ts0!,tol!,xx,mx,nseq)
      call upack('move nmodt modt ct',smove,nmodt,modt,ct,0,0)
      mode = igets('move prmint:1.1',smove)

C --- Make e-dot, overwrite into frc  ---
      call defrr(ode,3*nbas+nmodt)
      call mmdife(nbas,neq,gyro,aamom,nmodt,modt,ct,kt,sdir,frc,w(ode))
c     print *, 'comment  out overwrite frc'
      call dcopy(3*nbas+nmodt,w(ode),1,frc,1)
      call rlse(ode)

C --- e-dot wrt phi,theta ---
      call defrr(ofeu, neq)
c     call fx2fth(nbas,nmodt,modt,eula,sdir,nbas,frc,w(ofeu))
      call dx2dth(nbas,nmodt,modt,eula,sdir,nbas,frc,w(ofeu))
c     stop

C --- Bulirsch-Stoer integration step ---
      if (mode .eq. 1) then
C   ... Unpack BS parameters
        call upack1('move prmint',smove,prmint)
        mx = prmint(5)
        call discop(prmint(6),nseq,12,1,1,0)
        call defrr(oyerr, neq)
C        print *, ir,tnow
C        print *, eula(1,1),eula(1,2),eula(1,3)
C        print *, eula(2,1),eula(2,2),eula(2,3)
        call dbsstp(tnow,eula,w(ofeu),w(oyerr),neq,tol,ts0,ts,0d0,mx,
     .    nseq,ir,wkint,iprint())
C        print *, ir,tnow
C        print *, eula(1,1),eula(1,2),eula(1,3)
C        print *, eula(2,1),eula(2,2),eula(2,3)
        call rlse(ofeu)
      endif

      end

      subroutine dx2dth(nbas,nmodt,modt,eula,sdir,ndim,dxyz,dpol)
C- Transform d(x,y,z) into d(phi,th)
C  Also copies nmodt elements to end of dpol
C  phi = atan(y/x) => phidot = d/de atan (y + e f_y)/(x + e f_x)
C                            = (x*f_y - y*f_x) / (x*x + y*y)
C  th  = acos(z)  => thdot = d/de acos(z+ e f_z) = -f_z/sqrt(1-z*z)??
C  (can't remember how I got dpol(2) ... connected w/ cons. of length?)
C     implicit none
      integer nbas,i,ndim,iprint,m,nmodt,modt(3),i1mach
      double precision eula(nbas,3),dpol(ndim,2),dxyz(ndim,3),
     .  sdir(nbas,3)
      double precision st,ftop,xx,p,pi,pi8
C     double precision yy,zz,cp,ct,t,sp,fac
      character*80 outs

      pi = 4*datan(1d0)
      pi8 = 8*pi

C     call dscal(nbas*3,1d3,dxyz,1)

      do  10  i = 1, nbas

C       dpol(i,2) =  cp*ct*dxyz(i,1) + sp*ct*dxyz(i,2) - st*dxyz(i,3)

C   ... Get the sign of theta directly from the Euler angle
C       At st identically zero, dpol is indeterminant
        p = sdir(i,1)**2 + sdir(i,2)**2 + 1d-20
        st = dsqrt(p)*(1-2*mod(int((eula(i,2)+pi8)/pi),2)) + 1d-20
        dpol(i,1) = sdir(i,1)*dxyz(i,2)/p
     .            - sdir(i,2)*dxyz(i,1)/p
        dpol(i,2) = (sdir(i,1)*sdir(i,3)*dxyz(i,1) +
     .               sdir(i,2)*sdir(i,3)*dxyz(i,2) -
     .               (1-sdir(i,3)**2)*dxyz(i,3))/st

C   ... test numerically
C        fac = 1d-6
C        cp = dcos(eula(i,1))
C        ct = dcos(eula(i,2))
C        sp = dsin(eula(i,1))
C        st = dsin(eula(i,2))
C        xx = cp*st + fac*dxyz(i,1)
C        yy = sp*st + fac*dxyz(i,2)
C        zz = ct    + fac*dxyz(i,3)
C        p = sqrt(xx**2 + yy**2 + zz**2)
C        xx = xx/p
C        yy = yy/p
C        zz = zz/p
C
C        p = datan2(yy,xx)
C        t = dacos(zz)
C        if (eula(i,2) .gt. pi) then
C          t = pi+pi-t
C          p = p - pi
C        endif
C        if (dabs(p-eula(i,1)) .gt. 6) p = p - pi-pi
C        if ((p-eula(i,1)) .lt. -6) p = p +pi+pi
C
C        print 369, i,
C     .    dpol(i,1),(p-eula(i,1))/fac,
C     .    dpol(i,2),(t-eula(i,2))/fac,
C     .    st,t,eula(i,2)
C  369   format(i4,2f10.6,2x,2f10.6,3x,5f8.4)

   10 continue

C     pause

      do  11  i = 1, nmodt
   11 dpol(nbas+i,2) = dxyz(nbas+i,3)

C --- Printout ---
      if (iprint() .gt. 40) then
        print 387
        ftop = 0
        do  20  i = 1, nbas
          xx = dsqrt(dxyz(i,1)**2 + dxyz(i,2)**2 + dxyz(i,3)**2)
          ftop = dmax1(ftop,xx)
          print 345, i, eula(i,1),eula(i,2), (dpol(i,m),m=1,2), xx,
     .      (dxyz(i,m), m=1,3)
  345     format(i4,4f10.6,f8.4,1x,3f8.4)
   20   continue
        print 390, ftop
      endif
      if (iprint() .lt. 40) return
      outs = ' dx2dth:'
      if (nmodt .ne. 0) call awrit4('%a  xsi=%n:1;4d  xsidot=%n:1;4d',
     .  outs,len(outs),0,nmodt,sdir(nbas+1,3),nmodt,dpol(nbas+1,2))
      call awrit1('%a  fmax=%;4d',outs,len(outs),0,ftop)
      if (iprint() .ge. 40) call awrit0('%a',outs,-80,-i1mach(2))

  387 format('  ib    phi',6x,'theta      gph       gth      abs',7x,
     .  'gx      gy      gz')
  390 format(37x,'gmax = ',f8.4)

      end

      subroutine mmdife(nbas,neq,gyro,aamom,nmodt,modt,ct,kt,e,f,de)
C- Make de = gyro s x f + modt-function of T
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   neq   :not used!
Ci   gyro  :gyromagnetic ratio
Ci   aamom :local magnetic moments
Ci   nmodt :number of thermodynamics modes
Ci   modt  :thermodynamic modes:
Ci         :0, no temperature term
Ci         :one's digit:  1, a = (ey,ex,ez)
Ci         :              2, a = (ez,ey,ex)
Ci         :              3, a = (ex,ez,ey)
Ci         :              4, a = (ey,ez,ex)
Ci         :ten's digit:  1, h(xi) = xi
Ci         :              3, h(xi) = xi**3
Ci   ct    :thermodynamic mode coefficients
Ci   kt    :temperature
Ci   e     :direction vectors
Ci   f     :forces
Co Outputs
Co   de    :differential changes to direction vectors
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,nmodt,modt(3),neq
      double precision gyro,kt,e(nbas,3),f(nbas,3),de(nbas,3),
     .  aamom(nbas),ct(3)
C ... Local parameters
      double precision xsi,hxsi,xsidot,a(3),ae(3),eae(3),
     .  edotf,edota,adotf,xx
      integer i,ipwr,im,modj,ikt

C     call prmx('e',e,nbas,nbas,3)
C     call prmx('f',f,nbas,nbas,3)

C --- de = -gyro/aamom e x f ---
      do  5  i = 1, nbas
        xx = -gyro/aamom(i)
        a(1) = xx * (e(i,2)*f(i,3) - e(i,3)*f(i,2))
        a(2) = xx * (e(i,3)*f(i,1) - e(i,1)*f(i,3))
        a(3) = xx * (e(i,1)*f(i,2) - e(i,2)*f(i,1))
        de(i,1) = a(1)
        de(i,2) = a(2)
        de(i,3) = a(3)
    5 continue

C     call prmx('e-dot-0',de,nbas,nbas,3)

C --- Make thermostat terms ---
      do  6  im = 1, nmodt
      modj = mod(modt(im),10)
      if (modj .gt. 4) goto 999
      ipwr = mod(modt(im)/10,10)
      xsi = e(nbas+im,3)
      if (ipwr .ne. 1 .and. ipwr .ne. 3) goto 999
      hxsi = xsi**ipwr
      xsidot = 0

C --- Thermostat, modes 1,2,3,4 ---
      ikt = 1
      if (modj .eq. 4) ikt = 0
      do  10  i = 1, nbas
        if (modj .eq. 1) then
          a(1) = e(i,2)
          a(2) = e(i,1)
          a(3) = e(i,3)
        elseif (modj .eq. 2) then
          a(1) = e(i,3)
          a(2) = e(i,2)
          a(3) = e(i,1)
        elseif (modj .eq. 3) then
          a(1) = e(i,1)
          a(2) = e(i,3)
          a(3) = e(i,2)
        elseif (modj .eq. 4) then
          a(1) = e(i,2)
          a(2) = e(i,3)
          a(3) = e(i,1)
        endif

        ae(1) = a(2)*e(i,3) - a(3)*e(i,2)
        ae(2) = a(3)*e(i,1) - a(1)*e(i,3)
        ae(3) = a(1)*e(i,2) - a(2)*e(i,1)
        eae(1) = e(i,2)*ae(3) - e(i,3)*ae(2)
        eae(2) = e(i,3)*ae(1) - e(i,1)*ae(3)
        eae(3) = e(i,1)*ae(2) - e(i,2)*ae(1)
        de(i,1) = de(i,1) - hxsi*eae(1)
        de(i,2) = de(i,2) - hxsi*eae(2)
        de(i,3) = de(i,3) - hxsi*eae(3)
        edotf = e(i,1)*f(i,1) + e(i,2)*f(i,2) + e(i,3)*f(i,3)
        edota = e(i,1)*a(1) + e(i,2)*a(2) + e(i,3)*a(3)
        adotf = a(1)*f(i,1) + a(2)*f(i,2) + a(3)*f(i,3)
        xsidot = xsidot + (adotf - edotf*edota) - kt*(ikt - 3*edota)
   10 continue

      xsidot = ct(im)/nbas*xsidot
      de(nbas+im,3) = xsidot

C ... End of this thermostat term
    6 continue

c     call prmx('e x a x e',de,nbas,nbas,3)
c     call prmx('e-dot',de,nbas,nbas+1,3)
      return

C --- Error exit ---
  999 call fexit(-1,111,'  Exit -1 MMDIFE:  modet=%i not implemented',
     .  modt(im))

      end

C      subroutine mm_prj(nbas,e,ipr,f,etot,amag,tmom)
CC- Project out component of force along spin direction
C      implicit none
C      integer nbas
C      double precision etot,f(nbas,3),e(nbas,3),amag(3),tmom
C      integer ipr,i1mach,lgunit,ib,m
C      double precision ftop,xx,r2,fs
C      character*100 outs
C
CC --- Project out component along e ---
C      do  10  ib = 1, nbas
C        r2 = e(ib,1)**2 + e(ib,2)**2 + e(ib,3)**2
C        fs = e(ib,1)*f(ib,1) + e(ib,2)*f(ib,2) + e(ib,3)*f(ib,3)
C        f(ib,1) = f(ib,1) - fs*e(ib,1)
C        f(ib,2) = f(ib,2) - fs*e(ib,2)
C        f(ib,3) = f(ib,3) - fs*e(ib,3)
C        fs = e(ib,1)*f(ib,1) + e(ib,2)*f(ib,2) + e(ib,3)*f(ib,3)
C   10 continue
C
CC --- Printout ---
C      if (ipr .ge. 30) then
C        print *
C        xx = dsqrt(amag(1)**2 + amag(2)**2 + amag(3)**2)
C        call awrit3(' mm_prj: etot = %;6,6d  <mag> = %;5,5d '//
C     .    '<x,y,z> =%3:1;5,5d',outs,80,-lgunit(1),etot,xx,amag)
C        call awrit0('%a',outs,-80,-lgunit(2))
C      endif
C      if (ipr .gt. 40) then
C        print 389
C        ftop = 0
C        do  20  ib = 1, nbas
C          xx = dsqrt(f(ib,1)**2 + f(ib,2)**2 + f(ib,3)**2)
C          ftop = dmax1(ftop,xx)
C          write(i1mach(2),388)
C     .      ib, (f(ib,m),m=1,3),xx, (e(ib,m),m=1,3)
C   20   continue
C        write(i1mach(2),390) ftop
C      endif
C
C  388 format(i4,3f10.6,f8.4,3f10.6)
C  389 format('  ib',10x,'total force',10x,'abs val',14x,'dir')
C  390 format(18x,'mm_prj:  fmax = ',f8.4)
C      end

C      subroutine fx2fth(nbas,nmodt,modt,eula,sdir,ndim,fxyz,fpol)
CC- Transform f_xyz  into f_phi,th
C      implicit none
C      integer nbas,i,ndim,iprint,m,nmodt,modt(3),i1mach
C      double precision eula(nbas,3),fpol(ndim,2),fxyz(ndim,3),
C     .  sdir(nbas,1)
C      double precision st,ftop,xx,p,pi,pi8
C      double precision yy,zz,cp,ct,t,sp,x,y,z,cpe,spe,xe,ye,dephi,cte
C     .  ,ste,ze,deth
C      character*80 outs
C
C      pi = 4*datan(1d0)
C      pi8 = 8*pi
C
C      do  10  i = 1, nbas
C
CC   ... Get the sign of theta directly from the Euler angle
CC       At st identically zero, f is indeterminant
C        p = sdir(i,1)**2 + sdir(i,2)**2 + 1d-20
C        st = dsqrt(p)*(1-2*mod(int((eula(i,2)+pi8)/pi),2)) + 1d-20
C
C        fpol(i,1) = sdir(i,1)*fxyz(i,2) - sdir(i,2)*fxyz(i,1)
C
C        cp = dcos(eula(i,1))
C        sp = dsin(eula(i,1))
C        ct = dcos(eula(i,2))
C        fpol(i,2) = (fxyz(i,1)*cp*ct + fxyz(i,2)*sp*ct - fxyz(i,3)*st)
CC        fpol(i,2) = (fxyz(i,1)*sdir(i,1)*sdir(i,3)/st +
CC     .               fxyz(i,2)*sdir(i,2)*sdir(i,3)/st - fxyz(i,3)*st)
C
CC ...   test numerically ...
CC        cp = dcos(eula(i,1))
CC        sp = dsin(eula(i,1))
CC        ct = dcos(eula(i,2))
CC        x  = cp*st
CC        y  = sp*st
CC        z  = ct
CC        cpe = dcos(eula(i,1)+1d-6)
CC        spe = dsin(eula(i,1)+1d-6)
CC        xe = cpe*st
CC        ye = spe*st
CC        dephi = fxyz(i,1)*(xe-x) + fxyz(i,2)*(ye-y)
CC
CC        cte = dcos(eula(i,2)+1d-6)
CC        ste = dsin(eula(i,2)+1d-6)
CC        xe = cp*ste
CC        ye = sp*ste
CC        ze = cte
CC        deth= fxyz(i,1)*(xe-x) + fxyz(i,2)*(ye-y) + fxyz(i,3)*(ze-z)
CC        print 369, fpol(i,1), dephi/1d-6, fpol(i,2), deth/1d-6
CC  369   format('xx',4f10.6,3x,5f8.4)
C
C   10 continue
C
CC --- Printout ---
C      if (iprint() .gt. 40) then
C        print 387
C        do  20  i = 1, nbas
C          xx = dsqrt(fxyz(i,1)**2 + fxyz(i,2)**2 + fxyz(i,3)**2)
C          ftop = dmax1(ftop,xx)
C          print 345, i, eula(i,1),eula(i,2), (fpol(i,m),m=1,2), xx,
C     .      (fxyz(i,m), m=1,3)
C  345     format(i4,4f10.6,f8.4,1x,3f8.4)
C   20   continue
C        print 390, ftop
C      endif
C      if (iprint() .lt. 40) return
C      outs = ' dx2dth:'
C      if (nmodt .ne. 0) call awrit4('%a  xsi=%n:1;4d  xsidot=%n:1;4d',
C     .  outs,len(outs),0,nmodt,sdir(nbas+1,3),nmodt,fpol(nbas+1,2))
C      call awrit1('%a  fmax=%;4d',outs,len(outs),0,ftop)
C      if (iprint() .ge. 40) call awrit0('%a',outs,-80,-i1mach(2))
C
C  387 format('  ib    phi',6x,'theta      fph       fth      abs',7x,
C     .  'fx      fy      fz')
C  390 format(37x,'fmax = ',f8.4)
C
C      stop
C      end

