c#define ONE_ATOM_UNIT 1

      subroutine atomsc(lgdd,nl,nsp,lmax,z,rhozbk,kcor,lcor,qcor,rmax,a,
     .nr,rofi,ec,ev,pnu,qnu,idmod,v,dv,rhoin,rho,rhoc,nrmix,qc,sumec,
     .sumtc,sumev,ekin,utot,rhoeps,etot,amgm,rhrmx,vrmax,qtot,exrmax,
     .job,niter,lfrz)

      use m_globalvariables
C- Makes an atomic sphere self-consistent and get atomic charges
C ----------------------------------------------------------------
Ci Inputs
Ci   lgdd  :T  add q2 phi phidd into the density
Ci         :F  add q2 <phi phidd> phi phi into the density
Ci         :NB: both produce the same integrated density.
Ci         :Which one should be used depends on the context.
Ci   nl    :leading dimension of pnu,qnu.
Ci         :Also, total charge cutoff; see Remarks
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   lmax  :maximum l for this site (but see Remarks)
Ci   z     :nuclear charge
Ci   rhozbk:constant nuclear background density (jellium) added to z
Ci   kcor  :(partial core occupation) p.q.n for occupation
Ci   lcor  :(partial core occupation) l quantum for occupation
Ci   qcor  :(partial core occupation) core charge and moment
Ci   rmax  :potential, density calculated to rmax
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :                 -//-
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,,
Ci         :pnu = .5 - atan(dnu)/pi + (princ.quant.number).
Ci   qnu   :energy moments of charge (see Remarks)
Ci   idmod :0,1 or 2, specifing how the enu is set for an l-channel
CCi   pz   :pnu for second panel (npan=2)
CCi   qz   :qnu for second panel (npan=2)
CCi   idmoz:idmod for second panel; specifies how semicore state
CCi        :is treated
CCi   npan :>1 for two-panel calculations
Ci   v     :spherical potential (job='pot', else input v not used)
Ci   dv    :constant potential shift added to sphere.
Ci         :Not used now.
Ci   rhoin:input density (job='rho', else input rhoin not used)
Ci         :used internally as a work array
Ci   nrmix :number of prior densities to mix to accelerate convergence
Ci         :to self-consistency (Anderson mixing); see Remarks.
Ci         :nrmix<0, linear mixing, with mixing beta |nrmix/100|
Ci   job   :(see Remarks)
Ci         :job='pot': start with potential
Ci         :job='rho': start with rhoin.
Ci         :job='gue': start with guessed charge density.
Ci   niter :number of iterations to attempt convergence to
Ci          self-consistency (see Remarks)
Ci   lfrz  :0 for soft core, >0 for frozen core
Cio Input/Outputs
Cio  ec    :core eigenvalues.  On input, these are guessed values.
Cio        :if input ec(1) = 0, atomsc makes an internal choice fo ec.
Cio  ev    :valence eigenvalues  On input, these are guessed values.
Cio        :if input ec(1) = 0, atomsc makes an internal choice for ev.
Co Outputs
Co   rofi  :dimensioned (nr,2).
Co         :(*,1) radial mesh points
Co         :(*,2) weights
Co   rho   :spherical charge density = 4 pi r**2 rhotrue
Co   rhoc  :core charge density (unchanged if lfrz=1)
Co   qc:   :core electronic charge
Co   sumec :core single particle energy
Co   sumtc :core kinetic energy (unchanged if lfrz=1)
Co   ekin  :total kinetic energy
Co   utot  :electrostatic energy
Co   rhoeps:exchange-correlation energy
Co   etot  :sphere total energy
Co   amgm  :difference between spin up and spin down charge
Co   rhrmx :true density at rmax
Co   vrmax :true potential at rmax
Co   qtot  :net charge in sphere
Co   exrmax:exchange-correlation energy at rmax
Cr Remarks
Cr   Boundary conditions pnu, moments qnu, and the electrostatic
Cr   potential at rmax are all that is required to uniquely determine
Cr   a self-consistent spherical charge and potential.  'job'
Cr   determines how the starting potential is made, but the final
Cr   potential  and density should be independent of the initial choice.
Cr
Cr   atomsc uses the boundary condition that the potential at rmax
Cr   is zero, regardless of the total net charge inside the sphere.
Cr   See subroutine madpot for discussion of how this choice affects
Cr   the Madelung energy.
Cr
Cr   Sphere total energy is sum of K.E., Hartree energy, XC energy:
Cr      etot = ekin + utot + rhoeps
Cr   The kinetic energy is computed via double-counting terms
Cr     ekin = sumev + sumec + dsumec - rhovh - rhomu
Cb Bugs
Cb   Total energy terms need to be cleaned up and simplified.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      logical lfrz,lgdd
      character job*3
      integer nr,nsp,nl,nrmix,niter,kcor,lcor,ncmx,nvmx,lmax,
     .idmod(0:nl-1)
      parameter (ncmx=50, nvmx=20)
      double precision ec(ncmx),ev(nvmx),rofi(nr,2),v(nr,nsp),
     .rho(nr,nsp),rhoc(nr,nsp),rhoin(nr,nsp),pnu(nl,2),
     .qnu(3,nl,nsp),z,rmax,a,qc,vrmax(2),exrmax(2),dv,rhrmx,
     .rhozbk,qcor(2),amgm,ekin,etot,qtot,rhoeps,sumec,sumev,sumtc,
     .utot
C     double precision z(nl,2),qz(3,nl,3)
C     integer idmoz(nl)
C  Local parameters
      logical last,cmdopt,ltmp
      integer iprint,ir,isp,iter,jmix,lgunit,lrel,nglob,ipr1,l,ii,stdo
      character strn*10
      double precision b,ddot,decay,dl,dold,drdi,drho,dsumec,ea,fac,pi,
     .rho0t,rhomu,rhovh,rmsdel,ro,sum,tolrsq,vhrmax,vnucl,vrhoc,vsum,
     .zvnucl,rvh(2),rho0(2),reps(2),rmu(2),sec(2),stc(2),sev(2),dasum
C for Anderson mixing:
      integer:: amix , nmix
      real(8) ,allocatable :: rho_rv(:)

      double precision norm(10,10),awk(10,2),beta,beta1
      parameter (beta = 0.3d0)
c Heap allocation
c      integer w(1)
c      common /w/ w
C tolch is tolerance for change in the charge density, tolv for rseq
      double precision tolch,tl
      parameter (tolch=5d-5,tolrsq=1d-12)
      external lgunit

      if (lmax .ge. nl) call rx('atomsc:  lmax too large')
Changenglob      lrel = nglob('lrel')
      lrel = globalvariables%lrel
      stdo = lgunit(1)
      pi = 4d0*datan(1d0)
      b = rmax/(dexp(a*(nr-1)) - 1)
      vhrmax = 0d0
      nmix = min(max(nrmix,0),10)
      sec(2) = 0d0
      sec(1) = 0d0
      stc(2) = 0d0
      stc(1) = 0d0
      allocate(rho_rv(nr*nsp*2*(nmix+2)))
Ctakao_ZeroClear_NotRequiered        if (nr*nsp*2*(nmix+2)<0) rho_rv(:)=0.0d0


C --- Core charge, radial mesh points and weights ---
      if (kcor .ne. 0) then
        if (qcor(1) .ne. 0 .or. qcor(2) .ne. 0) then
          call info5(30,0,0,'%9fAdd core hole:  kcor=%i  lcor=%i'//
     .    '  qcor=%d  amom=%d',kcor,lcor,qcor,qcor(2),0)
        endif
      endif
      call getqvc(nsp,nl,lmax,z,pnu,qnu,0,0,kcor,lcor,qcor,qc,qtot,amgm,
     .ec,ev)
C --- Guesses for core and valence eigenvalues ---
      if (ec(1) .eq. 0)
     .call getqvc(nsp,nl,lmax,z,pnu,qnu,ncmx,nvmx,kcor,lcor,qcor,qc,
     .qtot,amgm,ec,ev)
      call radmsh(rmax,a,nr,rofi)
      call radwgt(rmax,a,nr,rofi(1,2))

C --- Moments printout ---
      if (iprint() .ge. 20) then
        ltmp=dasum(lmax+1,qnu(2,1,1),3)+dasum(lmax+1,qnu(3,1,1),3).eq.0
        if (ltmp) then
          call awrit5('  Pl= %n:-1,1;5#8d%?!n==2! spn 2  %n:-1,1;5#8d',
     .    ' ',180,stdo,lmax+1,pnu,nsp,lmax+1,pnu(1,2))
          call dcopy(lmax+1,qnu,3,awk,1)
          call dcopy(lmax+1,qnu(1,1,nsp),3,awk(1,2),1)
          call awrit5('  Ql= %n:-1,1;5#8d%?!n==2! spn 2  %n:-1,1;5#8d',
     .    ' ',180,stdo,lmax+1,awk,nsp,lmax+1,awk(1,2))
        elseif (iprint() .ge. 30) then
          write(stdo,891)
  891     format(
     .    '   l',8x,'pl',11x,'q0',11x,'q1',11x,'q2',5x,' id ',6x,'dl')
          do  11  isp = 1, nsp
            do  11  l = 0, lmax
              dl = dtan(pi*(.5d0 - pnu(l+1,isp)))
              if (dabs(dl) .gt. 9999) dl = 0
              write(stdo,100) l,pnu(l+1,isp),(qnu(ii,l+1,isp),ii=1,3),
     .        idmod(l),dl
  100         format(i4,4f13.7,i4,f13.7)
   11     continue
        endif
      endif

C --- Initial charge density ---
      if (job .eq. 'pot') then
        call newrho(z,lrel,lgdd,nl,1,lmax,a,b,nr,rofi,v,rhoin,rhoc,kcor,
     .  lcor,qcor,pnu,qnu,sec,stc,sev,ec,ev,tolrsq,nsp,lfrz,000)
        if (niter .eq. 0) then
          call dcopy(nr*nsp,rhoin,1,rho,1)
Cki#error, have return with len(w_varlist)>0 at line 189
          if (allocated(rho_rv)) deallocate(rho_rv)
          return

        endif
      else if (job .eq. 'gue') then
        decay = 1d0+z/10d0
        decay = dmin1(decay,5d0)
        decay = 5
C         write(6,335) decay
C  335   format(/' initialize exponential density with decay',f7.3)
        sum = 0d0
        do  4  ir = 1, nr
          ro = dexp(-decay*rofi(ir,1))*rofi(ir,1)**2
          rhoin(ir,1) = ro
          sum = sum+ro*a*(rofi(ir,1)+b)
    4   continue
        fac = z/(sum*nsp)
        do  5  ir = 1, nr
          rhoin(ir,1) = rhoin(ir,1)*fac
    5   continue
        if (nsp .eq. 2) then
          do  6  ir = 1, nr
            rhoin(ir,2) = rhoin(ir,1)
    6     continue
        endif
      else if (job .ne. 'rho') then
        call rx('atomsc: job not pot|rho|gue')
      endif
      call dpscop ( rhoin , rho_rv , nr * nsp , 1 , 1 + nr * nsp 
     .* ( nmix + 2 ) , 1d0 )


C --- Start self-consistency loop ---
      drho = 100d0
      last = .false.
      if (iprint() .ge. 30) write(stdo,341)
      jmix = 0
      dold = 1
      beta1 = beta
      if (nrmix .lt. 0) beta1 = dble(-nrmix)/100
      do  35  iter = 1, niter
        tl = tolrsq
C       Hartree potential
        call addzbk(rofi,nr,nsp,rhoin,rhozbk,-1d0)
        call poiss0(z,a,b,rofi,rhoin,nr,vhrmax,v,rvh,vsum,nsp)
C       call prrmsh('vh atom ',rofi,v,nr,nr,nsp)

        call addzbk(rofi,nr,nsp,rhoin,rhozbk,1d0)
        vnucl = v(1,1)
C       Exchange-correlation potential
        if (last .and. iprint() .ge. 50) call pshpr(80)
c        print *,'xxx goto vxc0sp'
        call vxc0sp(a,b,rofi,rhoin,nr,v,rho0,reps,rmu,nsp,exrmax)
        if (last .and. iprint() .ge. 50) call poppr
C       Get rhrmx, exrmax
        fac = 4*pi*rofi(nr,1)**2
        rhrmx = rhoin(nr,1)/fac
        if (nsp .eq. 2) then
          exrmax(2) = exrmax(1)
          rhrmx = rhrmx + rhoin(nr,2)/fac
        endif
        ipr1 = 0
        if (last .and. iprint().ge.40) ipr1 = 1
        if (last .and. iprint().gt.40) ipr1 = 2
        call newrho(z,lrel,lgdd,nl,1,lmax,a,b,nr,rofi,v,rho,rhoc,
     .  kcor,lcor,qcor,pnu,qnu,sec,stc,sev,ec,ev,tl,nsp,lfrz,ipr1)
        drho = 0
        sum = 0
        vrhoc = 0
        rho0t = 0
        do  40  isp = 1, nsp
          rho0t = rho0t + rho0(isp)
          sum = sum + ddot(nr,rofi(1,2),1,rho(1,isp),1)
          do  42  ir = 1, nr
            drdi = a*(rofi(ir,1) + b)
            drho = drho + rofi(ir,2)/drdi*dabs(rho(ir,isp)-rhoin(ir,isp))
   42     continue
          do  41  ir = 2, nr
            ea = (v(ir,isp)-2*z/rofi(ir,1))
            vrhoc = vrhoc + rofi(ir,2)*ea*rhoc(ir,isp)
   41     continue
   40   continue
C       call prrmsh('rho for atom ',rofi,rho,nr,nr,nsp)
        call dcopy ( nr * nsp , rho , 1 , rho_rv , 1 )

        jmix = amix ( nr * nsp , min ( jmix , nmix ) , nmix , 0 , beta1 
     .  , iprint ( ) - 70 , .9d0 , norm , awk ( 1 , 2 ) , rho_rv , 
     .  awk , rmsdel )

        call dpscop ( rho_rv , rhoin , nr * nsp , 1 + nr * nsp * ( 
     .  nmix + 2 ) , 1 , 1d0 )

C
        if (last) goto 90
        if (iprint() .ge. 41 .or. iprint() .ge. 30 .and.
     .  (drho .lt. tolch .or. iter .eq. niter-1 .or. iter .eq. 1))
     .  write(stdo,340) iter,sum,drho,vnucl,rho0t,vsum,beta1
  340   format(i5,f12.6,1p,e12.3,0p,f14.4,e14.4,f14.4,f7.2)
  341   format(/'  iter     qint',9x,'drho',10x,'vh0',10x,'rho0',
     .  10x,'vsum',5x,'beta')
        last = (drho .lt. tolch .or. iter .eq. niter-1)
        if (iprint() .gt. 100)  call query(' ',-1,0)
        jmix = jmix+1
C       Beta for next iteration
        beta1 = min(max((1-drho/dold)/beta1,beta1-.2d0,beta),
     .  1d0,beta1+.2d0)
        if (nmix .gt. 0 .and. drho .lt. 1) beta1 = 1
        if (nrmix .lt. 0) beta1 = dble(-nrmix)/100
        dold = drho
   35 continue
   90 continue
      if (allocated(rho_rv)) deallocate(rho_rv)

C     End of iteration loop
      if (iprint() .ge. 30) write(stdo,'(1x)')

C --- Collect terms for total energy ---
      rhoeps = 0
      rhomu  = 0
      sumev  = 0
      if (.not. lfrz) then
        sumec = 0
        sumtc = 0
      endif
      rhovh  = 0
      do  80  isp = 1, nsp
        if (nsp .eq. 2 .and. iprint() .gt. 30)
     .  write(stdo,230) isp,v(nr,isp)-2*z/rmax,
     .  sev(isp),sec(isp),rvh(isp),reps(isp),rmu(isp)
  230   format(' Spin',i2,':',
     .  /' vrmax=',f12.5,'    sumev= ',f12.5,'    sumec=',f12.5,
     .  /' rhovh=',f12.5,'    rhoeps=',f12.5,'    rhomu=',f12.5)
        rhoeps = rhoeps + reps(isp)
        rhomu = rhomu + rmu(isp)
        sumev = sumev + sev(isp)
        if (.not. lfrz) then
          sumec = sumec + sec(isp)
          sumtc = sumtc + stc(isp)
        endif
        rhovh = rhovh + rvh(isp)
  80  continue
      zvnucl = -z*vnucl
      utot = .5d0*(rhovh + zvnucl)
C ... Correction to core eigenvalues if sumec not obtained from this V
      dsumec = vrhoc - (sumec-sumtc)
      ekin = sumev + sumec + dsumec - rhovh - rhomu
      etot = ekin + utot + rhoeps
      if (iprint() .ge. 40) write(stdo,139) sumev,sumec,vnucl,rhovh,
     .zvnucl,utot,rhomu,rhoeps,dsumec,ekin,sumtc,etot
  139 format(/' sumev=',f13.6,'    sumec =',f13.6,'   vnucl =',f13.6
     ./' rhovh=',f13.6,'    zvnucl=',f13.6,'   utot  =',f13.6
     ./' rhomu=',f13.6,'    rhoeps=',f13.6,'   dsumec=',f13.6
     ./' ekin= ',f13.6,'    tcore =',f13.6,'   etot  =',f13.6)
      vrmax(1) = -2*z/rmax
      do  55  isp = 1, nsp
        vrmax(1) = vrmax(1) + v(nr,isp)/nsp
   55 continue
      vrmax(2) = 0
      if (nsp .eq. 2) vrmax(2) = v(nr,1)-v(nr,2)

C --- write out rho if requested ---
      if (cmdopt('--dumprho',8,0,strn)) then
        call prrmsh('rho for atom ',rofi,rho,nr,nr,nsp)
        call prrmsh('pot for atom ',rofi,v,nr,nr,nsp)
        call prrmsh('rhoc for atom ',rofi,rhoc,nr,nr,nsp)
        allocate(rho_rv(abs(-nr*nl*nsp)))
        if (-nr*nl*nsp<0) rho_rv(:)=0.0d0

        call newrho ( z , lrel , lgdd , nl , nl , lmax , a , b , nr , 
     .  rofi , v , rho_rv , rhoc , kcor , lcor , qcor , pnu , qnu 
     .  , sec , stc , sev , ec , ev , tl , nsp , lfrz , 000 )

        call prrmsh ( 'rhol for atom ' , rofi , rho_rv , nr , nr , 
     .  nl * nsp )

        if (allocated(rho_rv)) deallocate(rho_rv)

      endif

      end subroutine atomsc

      subroutine addzbk(rofi,nr,nsp,rho,rhozbk,scale)

C     implicit none
      integer nr,nsp
      double precision rofi(1),rho(nr,1),rhozbk,scale
      integer ir,isp
      double precision s

      if (rhozbk .eq. 0) return
      s = 16*datan(1d0)*scale*rhozbk
      do  10  isp = 1, nsp
        do  10  ir = 2, nr
          rho(ir,isp) = rho(ir,isp) + s*rofi(ir)*rofi(ir)
   10 continue
      end subroutine addzbk

      subroutine poiss0(z,a,b,rofi,rho,nr,vhrmax,v,rhovh,vsum,nsp)

C- Hartree potential for spherical rho.
C  ---------------------------------------------------------------------
Ci Inputs:
Ci   z     :nuclear charge
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   rho   :spherical charge density times 4*pi*r*r
Ci   nr    :number of mesh points
Ci   vhrmax:value of v(nr)
Ci   nsp   :=1 spin degenerate, =2 non-degenerate
Co Outputs:
Co   rofi  :radial mesh points
Co   v     :spherical Hartree potential
Co   vsum  :integral over that potential which is zero at rmax.
Co   rhovh :integral of rho*vh, where vh is the electrostatic
Co         :potential from both electronic and nuclear charge
Cr Remarks:
Cr    Solves Poisson's equation for given spherical charge density
Cr    and a specified value vhrmax at rofi(nr).
Cr    rho =  4*pi*r*r*rhotrue  but  v = vtrue
C ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer nr,nsp
      double precision rofi(nr),v(nr,nsp),rho(nr,nsp),rhovh(nsp)
      double precision pi,ea,rpb,a,b,z
      double precision rmax,r2,r3,r4,f2,f3,f4,x23,x34,cc,bb,dd,df,drdi,
     .r,srdrdi,g,f,y2,y3,y4,ro,vnow,wgt,a2b4,vhrmax,vsum,vhat0
      integer ir,isp

      pi = 4d0*datan(1d0)
      ea = dexp(a)
      rpb = b
      do  5  ir = 1, nr
        rofi(ir) = rpb - b
        rpb = rpb*ea
    5 continue
      rmax = rofi(nr)

C --- Approximate rho/r**2 by cc + bb*r + dd*r*r near zero  ---
      r2 = rofi(2)
      r3 = rofi(3)
      r4 = rofi(4)
      f2 = 0d0
      f3 = 0d0
      f4 = 0d0
      do  75  isp = 1, nsp
        f2 = f2 + rho(2,isp)/r2**2
        f3 = f3 + rho(3,isp)/r3**2
        f4 = f4 + rho(4,isp)/r4**2
   75 continue
      x23 = (r3*r3*f2 - r2*r2*f3)/(r3 - r2)
      x34 = (r4*r4*f3 - r3*r3*f4)/(r4 - r3)
      cc = (r2*x34 - r4*x23) / (r3*(r2 - r4))
      bb = ((r2+r3)*x34 - (r3+r4)*x23) / (r3*r3*(r4-r2))
      dd = (f2 - bb*r2 - cc)/r2**2

C --- Numerov for inhomogeneous solution ---
      a2b4 = a*a/4d0
      v(1,1) = 1d0
      df = 0d0
      do  12  ir = 2, 3
        r = rofi(ir)
        drdi = a*(r + b)
        srdrdi = dsqrt(drdi)
        v(ir,1) = v(1,1) - r*r*(cc/3d0 + r*bb/6d0 + r*r*dd/10d0)
        g = v(ir,1)*r/srdrdi
        f = g*(1d0 - a2b4/12d0)
        if (ir .eq. 2) y2 = -2d0*f2*r2*drdi*srdrdi
        if (ir .eq. 3) y3 = -2d0*f3*r3*drdi*srdrdi
        df = f - df
   12 continue
      do 13  ir = 4, nr
        r = rofi(ir)
        drdi = a*(r + b)
        srdrdi = dsqrt(drdi)
        ro = 0d0
        do  76  isp = 1, nsp
          ro = ro + rho(ir,isp)
   76   continue
        y4 = -2d0*drdi*srdrdi*ro/r
        df = df + g*a2b4 + (y4 + 10d0*y3 + y2)/12d0
        f = f + df
        g = f/(1d0 - a2b4/12d0)
        v(ir,1) = g*srdrdi/r
        y2 = y3
        y3 = y4
   13 continue

C --- Add constant to get v(nr) = vhrmax ---
      vnow = v(nr,1) - 2d0*z/rmax
      do  27  ir = 1, nr
        v(ir,1) = v(ir,1) + (vhrmax - vnow)
   27 continue

C --- Integral vh and  rho * vh ---
      rhovh(1) = 0d0
      rhovh(nsp) = 0d0
      vsum = 0d0
      vhat0 = 0d0
      do  30  ir = 2, nr
        r = rofi(ir)
        drdi = a*(r + b)
        wgt = 2*(mod(ir+1,2) + 1)/3d0
        if (ir .eq. nr) wgt = 1d0/3d0
        ro = 0d0
        do  31  isp = 1, nsp
          rhovh(isp) =  rhovh(isp)
     .    + wgt*drdi*rho(ir,isp)*(v(ir,1) - 2d0*z/r)
          ro = ro + rho(ir,isp)
   31   continue
        vhat0 = vhat0 + wgt*drdi*ro*(1d0/r - 1d0/rmax)
        vsum = vsum + wgt*drdi*r*r*(v(ir,1) - vhrmax)
   30 continue
      vsum = 4d0*pi*(vsum - z*rmax*rmax)
      vhat0 = 2d0*vhat0 + 2d0*z/rmax + vhrmax
      v(1,1) = vhat0

C --- Copy to second spin channel if spin polarized ---
      if (nsp .eq. 1) return
      do  84  ir = 1, nr
        v(ir,2) = v(ir,1)
   84 continue
      end subroutine poiss0


C      subroutine fctp0(l,rofi,v,z,nr,nctp0,xrim,xmin,nsave)
CC- Initialize things for FCTP, which finds classical turning point
CC  ---------------------------------------------------
CCi Inputs
CCi   l     :angular quantum number
CCi   rofi  :radial mesh points
CCi   v     :spherical potential (atomsr.f)
CCi   z     :nuclear charge
CCi   nr    :number of radial mesh points
CCo Outputs:
CCo   nctp0:
CCo   xrim:
CCo   xmin:
CCo   nsave: Estimate for nctp
CC  ---------------------------------------------------
C      implicit none
C      integer l,nctp0,nr,nsave
C      double precision xmin,xrim,z,v(nr),rofi(nr)
CC Local
C      integer ir
C      double precision fllp1,r,x,xlast
C      fllp1 = l*(l + 1)
C      r = rofi(10)
C      x = fllp1/r/r - 2d0*z/r + v(10)
C      ir = 10
C   80 ir = ir + 1
C      xlast = x
C      r = rofi(ir)
C      x = fllp1/r/r - 2d0*z/r + v(ir)
C      if (x .le. xlast .and. ir .lt. nr) goto 80
C      nctp0 = ir - 1
C      xmin = xlast
C      r = rofi(nr)
C      xrim = fllp1/r/r - 2d0*z/r + v(nr)
C      if (xmin .ge. xrim-3d0) nctp0 = nr
C      if (xmin .ge. xrim-3d0) xmin = xrim
C      nsave = (nctp0 + nr)/2
C      end
C      subroutine fctp(e,nctp,nctp0,xrim,xmin,nsave,l,rofi,v,z,nr,a,b)
CC- Finds classical turning point for wave function
CC ----------------------------------------------------------------------
CCi Inputs
CCi   e     :energy
CCi   nctp0 :
CCi   xrim  :
CCi   xmin  :
CCi   nsave: Estimate for nctp
CCi   l     :l quantum number
CCi   rofi  :radial mesh points
CCi   v     :spherical potential
CCi   z     :nuclear charge
CCi   nr    :number of radial mesh points
CCi   a     :mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
CCi   b     :mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
CCo Outputs:
CCo   nctp:  rofi(nctp) is classical turning point
CCo   nsave: New estimate for NCTP, next iteration
CCr Remarks
CC  ---------------------------------------------------
C      implicit none
C      integer l,nctp,nctp0,nr,nsave
C      double precision a,b,e,xmin,xrim,z,v(nr),rofi(nr)
CC Local
C      integer irep,k,n1,n2,nlast,ntry
C      double precision dfdr,dvdr,fllp1,fntry,fofr,r,rtry,vme
C
C      fllp1 = l*(l+1)
C      if (nctp0 .eq. nr) then
C        nctp = nr
C        return
C      endif
C      if (e .gt. xrim) then
C        nctp = nr
C        return
C        endif
C      if (e .lt. xmin) then
C        nctp = 2
C        return
C      endif
C      n1 = nctp0
C      n2 = nr
C      nctp = nsave
C      nlast = -10
C      do  20  irep = 1, 20
C      if (nctp.gt.n2 .or. nctp.lt.n1) nctp=(n1+n2+1)/2
C        r = rofi(nctp)
C        vme = (v(nctp)-e)
C        if (nctp .gt. nr) call rx('fctp0: nctp gt nr')
C        k = min(nctp,nr-1)
C        dvdr = (v(k+1)-v(k-1))/(2.d0*a*(r+b))
CC       dvdr = (v(nctp+1)-v(nctp-1))/(2.d0*a*(r+b))
C        fofr = fllp1/r/r - 2.d0*z/r + vme
C        dfdr = -2.d0*fllp1/r/r/r+2.d0*z/r/r + dvdr
C        rtry = r - fofr/dfdr
C        rtry = dmax1(rtry,rofi(2))
C        fntry = dlog(rtry/b+1.d0)/a + 1.d0
C        ntry = fntry + .5d0
CC|      write(6,810) irep,n1,nctp,n2,fntry,ntry
CC|810   format(i6,'   n1,nctp,n2=',3i5,'   ntry=',f8.3,i6)
C        if (nlast .eq. nctp) goto 98
C        if (fofr .gt. 0.d0) n2=nctp
C        if (fofr .lt. 0.d0) n1=nctp
C        nlast = nctp
C        nctp = ntry
C   20 continue
C   98 if (nctp .eq. nctp0+1) nctp = 2
C      nsave = nctp
C      end

      subroutine fctp0(l,nr,rofi,v,z,nctp0)

C- Initialize things for FCTP, which finds classical turning point
C ----------------------------------------------------------------------
Ci Inputs:
Ci   l     :angular momentum
Ci   rofi  :radial mesh points
Ci   v     :spherical potential
Ci   z     :nuclear charge
Ci   nr    :number of mesh points
Co Outputs:
Co   nctp0 :minimum of effective potential
Co          or nr if v(nr) > vmin + 3 Ry
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer l,nr,nctp0
      double precision z,v(nr),rofi(nr)
C Local variables:
      integer ir,irmin
      double precision veff,zz,fllp1,vi,vim1
      parameter (irmin=11)
C Statement functions:
      veff(ir)=fllp1/(rofi(ir)*rofi(ir))-zz/rofi(ir)+v(ir)

      zz = z+z
      fllp1 = l*(l+1)

      ir = irmin
      vi = veff(ir)
      vim1 = veff(ir-1)
   10 if (vi.le.vim1 .and. ir.lt.nr) then
        ir = ir+1
        vim1 = vi
        vi = veff(ir)
        goto 10
      endif

      nctp0 = ir-1
      if (veff(nctp0).ge.veff(nr)-3.d0) nctp0 = nr

      end subroutine fctp0

      subroutine fctp(a,b,e,l,nctp0,nr,rofi,v,z,nctp)

C- Finds classical turning point
C  ---------------------------------------------------------------------
Ci Inputs:
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :                 -//-
Ci   e     :energy
Ci   l     :angular momentum
Ci   nctp0 :minimum of effective potential (fctp0)
Ci   nr    :number of mesh points
Ci   rofi  :radial mesh points
Ci   v     :spherical potential
Ci   z     :nuclear charge
Ci Inputs/Outputs:
Cio  nctp  :rofi(nctp) is classical turning point
Cu Updates
Cu   13 Jun 00 Added safety checks to guard against jumps in potential
C  ---------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer nctp,nctp0,l,nr
      double precision e,rofi(nr),v(nr),z,a,b
C Local variables:
      integer ir,irep,n1,n2,nlast,ntry
      double precision r,veff,fllp1,fofr,dfdr,rtry,zz
C Intrinsic functions:
      intrinsic dlog,dmax1,min0
C Statement functions:
      veff(ir)=fllp1/(rofi(ir)*rofi(ir))-zz/rofi(ir)+v(ir)

      zz=z+z
      fllp1 = l*(l+1)

      if (nctp0.eq.nr .or. e.gt.veff(nr)) then
        nctp = nr
      elseif (e.lt.veff(nctp0)) then
        nctp = 2
      else
        n1 = nctp0
        n2 = nr-1
        nlast = -10
        do  10  irep = 1, 20
          if (nctp.gt.n2 .or. nctp.lt.n1) nctp = (n1 + n2)/2
          r = rofi(nctp)
          fofr =  veff(nctp)-e
          dfdr =-(fllp1+fllp1)/r/r/r + zz/r/r +
     .    (v(nctp+1) - v(nctp-1))/(2.d0*a*(r+b))
          rtry = dmax1(r-fofr/dfdr,rofi(2))
          ntry = dlog(rtry/b+1.d0)/a + 1.5d0

C         If there was a large change, check for safety
          if (nlast .eq. nctp .and. iabs(ntry-nctp) .gt. 10) then
            if (fofr .gt. 0) then
              do  ntry = nctp-1, nctp0+1, -1
                fofr =  veff(ntry)-e
                if (fofr .lt. 0) goto 20
                nctp = ntry
                n2 = ntry
              enddo
   20         continue
            else
              do  ntry = nctp+1, nr-1
                fofr =  veff(ntry)-e
                if (fofr .gt. 0) goto 21
                nctp = ntry
                n1 = ntry
              enddo
   21         continue
            endif
          endif

C         Exit point
          if (nlast .eq. nctp) then
            if (nctp .eq. nctp0+1) nctp = 2
            return
          endif
          if (fofr .gt. 0.d0)  then
            n2 = nctp
          else
            n1 = nctp
          endif
          nlast = nctp
          nctp = min0(ntry,nr-1)
   10   continue
        if (nctp .eq. nctp0+1) nctp = 2
      endif
      end subroutine fctp


      subroutine newrho(z,lrel,lgdd,nl,nlr,lmax,a,b,nr,rofi,v,rho,rhoc,
     .kcor,lcor,qcor,pnu,qnu,sumec,sumtc,sumev,ec,ev,tol,nsp,lfrz,ipr)

C- Makes spherical charge density for a spherical potential.
C  ---------------------------------------------------
Ci Inputs:
Ci   z     :nuclear charge
Ci   lrel  :0 for nonrelativistic, 1 for relativistic
Ci   lgdd  :T q2 is coefficient to phidot**2 + phi*phidotdot
Ci         :F q2 is coefficient to phidot**2 - p phi**2
Ci         :Both produce the same integrated density; see Remarks.
Ci         :lgdd=F follows Stuttgart conventions.
Ci   nl    :(global maximum l) + 1
Ci   nlr   :second dimension of rho:
Cr         :1 if spherical part of rho is to be generated
Ci         :nl if generated rho is to be decomposed by l
Ci         :In the latter case, the core is not calculated
Ci   lmax  :maximum l for a given site
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Ci   v     :spherical potential (atomsr.f)
Ci   kcor  :(partial core occupation) p.q.n for occupation
Ci   lcor  :(partial core occupation) l quantum for occupation
Ci   qcor  :(partial core occupation) core charge and moment
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,,
Ci         :pnu = .5 - atan(dnu)/pi + (princ.quant.number).
Ci   qnu   :energy moments of charge (see Remarks)
Ci   tol   :precision to which wave functions are integrated.
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   g     :normalized wave function times r
Ci   gp    :energy derivatives of g
Ci   lfrz  :T, do not make core rho
Ci   ipr   :0 no printout
Ci         :1 summary printout of core
Ci         :2 detailed printout of core
Co Outputs:
Co   rho   :spherical charge density times 4*pi*r*r
Co   rhoc  :core density times 4*pi*r*r
Co   sumec :sum of core eigenvalues
Co   sumtc :core kinetic energy = sumec - v*rhoc
Co   sumev :sum of valence eigenvalues
Co   ec    :core eigenvalues
Co   ev    :valence eigenvalues
Cr Remarks:
Cr   rho is determined by boundary conditions pnu and moments qnu.
Cr   For rmax>10 sets phidot, phidotdot are made zero.
Cr
Cr   Recall that val,slo correspond to u = r*phi, so
Cr     rmax * slo / val = D_u = D_phi + 1.
Cr     For val=rmax  slo = D + 1
Cr
Cr   Switch lgdd concerns the contribution of <phidot|phidot> to the
Cr   sphere charge.  The zeroth moment may be defined as the amount of
Cr   charge in channel l, i.e.
Cr     q^0 = (amount of <phi|phi>) +  p* (amount of <phidot|phidot>)
Cr   where
Cr     p = <phidot|phidot>
Cr   Then q^0 is not the amount of phi*phi inside the sphere, but rather
Cr   but instead, it is
Cr     (q_0 - p q_2)
Cr   since q2 is the amount of <phidot|phidot>.  The charge density is
Cr     rho= (q_0 - p q_2) phi*phi + 2 q_1 phi*phidot + q_2 phidot*phidot
Cr   This is the Stuttgart convention (lgdd=F)
Cr
Cr   Methfessel convention:  to avoid explicit dependence of rho on
Cr   p, he approximated p*phi*phi with -phi*phidotdot (they have the
Cr   same integrated charge).  Then
Cr     rho= q_0 phi*phi + 2 q_1 phi*phidot +
Cr          q_2 (phidot*phidot + phi*phidotdot)
C  ---------------------------------------------------
C     implicit none
C Passed parameters
      logical lgdd,lfrz
      integer nl,nlr,lmax,nr,nsp,lrel,kcor,lcor,ipr
      double precision z,a,b,tol,sumev(nsp),sumec(nsp),sumtc(nsp),ec(1),
     .ev(1),qcor(2),v(nr,nsp),rofi(nr),rho(nr,nlr,nsp),rhoc(nr,nsp),
     .qnu(3,nl,nsp),pnu(nl,nsp)
C Local parameters
      logical free
      integer konfig(0:10),l,isp,ir,ival,nn,nre,jr,lmaxc,lr,k,nrmx
      parameter (nrmx=1501)
      double precision rocrit,pi,eb1,eb2,q0,q1,q2,rmax,eval,dl,val(5),
     .slo(5),sum,ro,phi,dphi,phip,dphip,p,fllp1,r,tmc,gfac,c,q00
      double precision g(2*nrmx),gp(2*nrmx*4)
C     Speed of light, or infinity in nonrelativistic case
      common /cc/ c

      call setcc(lrel)

      if (nr .gt. nrmx) call rxi(' newrho: increase nrx, need',nr)
      lr = 1
      rocrit = 0.002d0/4
      pi = 4d0*datan(1d0)
      eb1 = -50d0
      eb2 =  50d0
      rmax = rofi(nr)
      free = (rmax .gt. 9.99d0)
      call config(pnu,lmax,z,konfig,lmaxc)
      if (kcor .gt. 0) then
        lmaxc = max(lmaxc,lcor)
        konfig(lcor) = max(konfig(lcor),kcor+1)
      endif

C --- Calculate core density ---
      if (nlr .eq. 1) then
        if (.not. lfrz) then
          call dpzero(rhoc,nr*nsp)
          call rhocor(0,z,lmaxc,nsp,konfig,a,b,nr,rofi,v,g,
     .    kcor,lcor,qcor,tol,ec,sumec,sumtc,rhoc,0d0,ipr)
        endif
        call dcopy(nr*nsp,rhoc,1,rho,1)
      endif

C --- Loop over valence states ---
      ival = 0
      do  20  isp = 1, nsp
        sumev(isp) = 0d0
        do  20  l = 0, lmax

          if (nlr .eq. nl) lr = l+1
          q0 = max(qnu(1,l+1,isp),0d0)
          q1 = qnu(2,l+1,isp)
          q2 = qnu(3,l+1,isp)
          if (q0 .lt. 1d-6) goto 20
          nn = int(pnu(l+1,isp)) - l - 1
          ival = ival+1
          eval = ev(ival)

          val(1) = rmax
          dl = dtan(pi*(0.5d0 - pnu(l+1,isp)))
          slo(1) = dl + 1

          if (free) val(1) = 1d-30
          if (free) slo(1) = -val(1)
          call rseq(eb1,eb2,eval,tol,z,l,nn,val,slo,v(1,isp),g,
     .    sum,a,b,rofi,nr,nre)
          ev(ival) = eval
          sumev(isp) = sumev(isp) + eval*q0 + q1
          ro = g(nr)**2
          if (.not.free .and. ro .lt. rocrit) write(*,766) l,nn,nre,ro
  766     format(' NEWRHO (warning): PHP,PHPP set to zero,l,nn,nre,rho=',
     .    3i5,2f8.4)
          if (free .or. ro .lt. rocrit) then
            call dpzero(gp,8*nr)
            p = 0
          else
            val(1) = val(1)/dsqrt(sum)
            slo(1) = slo(1)/dsqrt(sum)
C          call phidot(z,l,v(1,isp),eval,a,b,rofi,nr,g,val,slo,tol,
C     .                nn,gp,phi,dphi,phip,dphip,p)
            call phidx(1,z,l,v(1,isp),0d0,0d0,rofi,nr,2,tol,eval,val,
     .      slo,nn,g,gp,phi,dphi,phip,dphip,p,0d0,0d0,0d0,0d0)
          endif
          fllp1 = l*(l+1)
C  ...  Case add q2 phi phidd rho
          if (lgdd) then
            k = 2*nr
            do  21  ir = 2, nre
              jr = ir + nr
              r = rofi(ir)
              tmc = c - (v(ir,isp) - 2d0*z/r - eval)/c
              gfac = 1d0 + fllp1/(tmc*r)**2
              rho(ir,lr,isp) =  rho(ir,lr,isp) +
     .        q0*(gfac*g(ir)**2 + g(jr)**2) +
     .        2*q1*(gfac*g(ir)*gp(ir) + g(jr)*gp(jr)) +
     .        q2*(gfac*(gp(ir)**2 + g(ir)*gp(ir+k)) +
     .        gp(jr)**2 + g(jr)*gp(jr+k))
  21        continue
C  ...  Case add -p q2 phi phi into rho
          else
            q00 = q0-p*q2
            do  22  ir = 2, nre
              jr = ir + nr
              r = rofi(ir)
              tmc = c - (v(ir,isp) - 2d0*z/r - eval)/c
              gfac = 1d0 + fllp1/(tmc*r)**2
              rho(ir,lr,isp) = rho(ir,lr,isp) +
     .        q00*(gfac*g(ir)**2 + g(jr)**2) +
     .        2*q1*(gfac*g(ir)*gp(ir) + g(jr)*gp(jr)) +
     .        q2*(gfac*gp(ir)**2 + gp(jr)**2)
   22       continue
          endif
  20  continue
      end subroutine newrho

      subroutine rhocor(isw,z,lmax,nsp,konfig,a,b,nr,rofi,v,g,kcor,lcor,
     .qcor,tol,ec,sumec,sumtc,rho,gcore,ipr)

C- Generates the (spherical) charge density from the core states
C ----------------------------------------------------------------------
Ci Inputs:
Ci   isw   :1 return the core wave functions in gcore
Ci   z     :nuclear charge
Ci   lmax  :maximum l
Ci   nsp   :=1 for non-spin-polarized =2 for spin-polarized calculations
Ci   konfig:core configuration. Core orbitals are specified by:
Ci         :  1, 2, ..., konf(0)-1 for s
Ci         :  2, 3, ..., konf(1)-1 for p
Ci         :  3, 4, ..., konf(2)-1 for d, and so on.
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Ci   v     :spherical potential (electronic contribution)
Ci   g     :work array holding normalized wave function times r
Ci   kcor  :(partial core occupation) p.q.n for occupation
Ci   lcor  :(partial core occupation) l quantum for occupation
Ci   qcor  :(partial core occupation) core charge and moment
Ci   tol   :wave function tolerance
Ci   ipr   :0 no printout
Ci         :1 summary printout
Ci         :2 detailed printout
Cio Inputs/Outputs:
Cio  ec    :guessed core eigenvalues (input)
Cio         core eigenvalues (output)
Co Outputs:
Co   rho   :spherical charge density times 4*pi*r*r
Co         :the spherical charge density of the core states is added
Co   sumec :sum of core eigenvalues
Co   sumtc :core kinetic energy
Co   gcore :(isw=0) not used
Co         :(isw=1) core wave functions
Cl Local variables
Cl   deg   :orbital occupation number
Cr Remarks:
Cu Updates
Cu   19 Apr 01 core wave functions may be saved in gcore
C  ---------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer isw,ipr,lmax,konfig(0:lmax),nr,nsp,kcor,lcor
      double precision a,b,z,g(nr,2),rho(nr,nsp),rofi(nr),tol,qcor(2),
     .sumec(nsp),sumtc(nsp),v(nr,nsp),gcore(nr,2,*),ec(*)
C ... Local parameters
      character*1 pqn(9),ang(6)
      integer icore,isp,konf,l,ll,nodes,nre,stdo,lgunit
      double precision deg,dlml,e1,e2,ecor0,ecore,qcore,rhorim,rmax,
     .rorim,slo,sum,tcore,val,vrho,kappa2
C ... Heap
      data pqn /'1','2','3','4','5','6','7','8','9'/
      data ang /'s','p','d','f','g','h'/

      stdo = lgunit(1)
      rmax  = rofi(nr)
      b     = rmax/(dexp(a*nr-a)-1.d0)
      e1    = -2.5d0*z*z - 5d0
      e2    = 20d0
      icore = 0
      do  80  isp = 1, nsp
        sumec(isp) = 0d0
        sumtc(isp) = 0d0
        qcore = 0d0
        rhorim = 0d0
        if (ipr .ge. 2) write(stdo,757)
        if(ipr.ge.2 .and. nsp.eq.2)
     .  write(stdo,'('' spin'',i2,'':'')') isp
  757   format(/' state  chg          ecor0',10x,'ecore',10x,'tcore',
     .  4x,'nre',2x,'rho(rmax)')
        do  10  l = 0, lmax
          do  10  konf = l, konfig(l)-2
            deg = (2*(2*l+1))/nsp
            if (konf+1.eq.kcor .and. l.eq.lcor) then
              deg = deg + (qcor(1)+(3-2*isp)*qcor(2))/nsp
            endif
            icore = icore+1
            nodes = konf - l
            ecor0 = ec(icore)
            val = 1d-30
            slo = -val
            call rseq(e1,e2,ecor0,tol,z,l,nodes,val,slo,v(1,isp),
     .      g,sum,a,b,rofi,nr,nre)
            ecore = ecor0
C     ... Correct core energy by using hankel bc's
            kappa2 = ecor0 - v(nr,isp) + 2*z/rmax
            if (nre .eq. nr .and. kappa2 .lt. 0d0) then
              dlml = -1d0 - dsqrt(-kappa2)*rmax
              do  31  ll = 1, l
                dlml = -kappa2*rmax*rmax/dlml - (2*ll+1)
   31         continue
              slo = val*(dlml+l+1)/rmax
              call rseq(e1,e2,ecore,tol,z,l,nodes,val,slo,v(1,isp),
     .        g,sum,a,b,rofi,nr,nre)
            endif
            ec(icore) = ecore
            if (isw .eq. 1) call dcopy(2*nr,g,1,gcore(1,1,icore),1)

C     --- Add to rho, make integral v*rho ---
            call xyrhsr(ecore,l,z,a,b,nr,nre,g,rofi,v(1,isp),rho(1,isp),
     .      deg,vrho,rorim)
            rhorim = rhorim + rorim
            tcore = ecore - vrho
            qcore = qcore + deg
            sumec(isp) = sumec(isp) + deg*ecore
            sumtc(isp) = sumtc(isp) + deg*tcore
            if (ipr.ge.2) write(stdo,758)
     .      pqn(konf+1),ang(l+1),deg,ecor0,ecore,tcore,nre,rorim
  758       format(1x,2a1,f8.2,3f15.6,i7,f9.5)
   10   continue
        if (ipr.gt.0) write(stdo,230) qcore,sumec(isp),sumtc(isp),rhorim
  230   format(' sum q=',f5.2,'  sum ec=',f12.5,'  sum tc=',f12.5,
     .  '  rho(rmax)',f8.5)
   80 continue
      end subroutine rhocor


      subroutine xyrhsr(ecore,l,z,a,b,nr,nre,g,rofi,v,rho,deg,vrho,
     .rhormx)

C- Make density and integrate potential*density for one core state
C ----------------------------------------------------------------------
Ci Inputs
Ci   ecore :core eigenvalue
Ci   l     :l quantum number
Ci   z     :nuclear charge
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   nre   :Make density to the smaller of nr and nre
Ci   g     :normalized wave function times r
Ci   rofi  :radial mesh points
Ci   v     :electronic part of spherical potential
Ci   deg   :occupation number
Co Outputs
Co   rho   :contribution core density from this state from 1..min(nr,nre)
Co   vrho  :integral
Co   rhormx:contribution to true density at nr from this state
Cr Remarks
Cr   xyrhsr makes the density at points 1..min(nr,nre), and the integral
Cr   of rho*density from 1..nre.  Thus:
Cr     if nre .eq. nr, the density and integral are the same
Cr     if nre .lt. nr  g**2 is negligible at nr (deep core state)
Cr     if nre .gt. nr  (large sphere) contribution to rho*v is included
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nr,nre,l
      double precision a,b,z,ecore,g(nr,2),rofi(1),v(nr),rho(nr),rhormx
C ... Local parameters
      integer nrmx,ir
      double precision fpi,fllp1,vrho,r,wgt,tmc,c,gfac,rhoir,deg,rmax
C     Speed of light, or infinity in nonrelativistic case
      common /cc/ c

      fpi = 16d0*datan(1d0)
      fllp1 = l*(l+1)
      vrho = 0
      nrmx = min0(nr,nre)
C ... Make rho, and integrate vrho for points 1..nrmx
      do  11  ir = 2, nrmx
        r = rofi(ir)
        wgt = 2*(mod(ir+1,2)+1)
        if (ir .eq. nre) wgt = 1
        tmc = c - (v(ir) - 2*z/r - ecore)/c
        gfac = 1 + fllp1/(tmc*r)**2
        rhoir = gfac*g(ir,1)**2 + g(ir,2)**2
        vrho = vrho + wgt*rhoir * (v(ir)-2*z/r) * (r+b)
        rho(ir) = rho(ir) + deg*rhoir
   11 continue
      rmax = rofi(nr)
      rhormx = 0
C     nfp has the following line:
C     if (nre .ge. nr) rhormx = deg*g(nr,1)**2/(rmax*rmax*fpi)
      if (nre .ge. nr) rhormx = deg*rhoir/(fpi*rmax**2)

C ... Integrate rho*v from nrmx+1 .. nre
      do  12  ir = nrmx+1, nre
        r = rofi(ir)
        wgt = 2*(mod(ir+1,2)+1)
        if (ir .eq. nre) wgt = 1
        tmc = c - (v(ir) - 2*z/r - ecore)/c
        gfac = 1 + fllp1/(tmc*r)**2
        rhoir = gfac*g(ir,1)**2 + g(ir,2)**2
        vrho = vrho + wgt*rhoir * (v(ir)-2*z/r) * (r+b)
   12 continue
      vrho = vrho*a/3

      end subroutine xyrhsr


C      subroutine phidot(z,l,v,e,a,b,rofi,nr,g,val,slo,tol,nn,
C     .  gp,phi,dphi,phip,dphip,p)
CC- Generate Phidot,Phidotdot for a prescribed energy
CC ----------------------------------------------------------------
CCi Inputs:
CCi   z:     nuclear charge
CCi   l:     l quantum number for this g
CCi   v:     potential
CCi   a,b:   defines shifted logarithmic mesh (rmesh.f)
CCi   rofi:  list of points; must be consistent with a,b
CCi   nr:    number of mesh points
CCi   e:     Energy
CCi   val:   Value of r * wave function at sphere radius (rseq)
CCi   slo:   Derivative of r * wave function at sphere radius (rseq)
CCi   g:     Wave function times r normalized so that int (g*g) dr = 1
CCi   tol:   precision to which wave function is integrated
CCi   nn:    number of nodes
CCo Outputs:
CCo   gp:    first four energy derivatives to G
CCo   phi:   wave function at rmax, i.e. g/rmax
CCo   dphi:  slope of wave function at rmax, i.e. d(g/rmax)/dr
CCo   phip:  energy derivative of wave function at rmax
CCo   dphip: energy derivative of slope of wave function at rmax
CCo   p:     <gp**2> (potential parameter)
CCr Remarks:
CCr   This version makes energy derivatives by numerical differentiation
CCr   of wave function phi, and has the same calling sequence as the
CCr   analytic version phidot.  The only difference is that this routine
CCr   returns four derivatives of phi, whereas the analytic version only
CCr   returns two and is applicable only to the nonrelativistic case.
CC ----------------------------------------------------------------
C      implicit none
CC passed parameters
C      integer l,nr,nn
C      double precision z,e,a,b,val,slo,phi,dphi,phip,dphip,p,tol
C      double precision v(1),rofi(1)
C      double precision g(2*nr),gp(2*nr,4)
CC local variables
C      integer nre,i,iprint,nptdif
C      double precision rmax,eb1,eb2,dele,ddde,sum1,
C     .                 vali(5),sloi(5),ei(4),de1,de2,del1,del2
C      parameter (nptdif = 2)
C
C      rmax = rofi(nr)
C      eb1 = -50d0
C      eb2 = 20d0
C
Cc      dele = tol**.2D0
C      dele = .002d0
C      if (tol .gt. 1d-9 .and. iprint() .ge. 0)
C     .  print *, 'phidot: tol too high for reliable num. diff'
C
C      ddde = -rmax/g(nr)**2
C      ei(1) = 1
C      ei(2) = -1
C      ei(3) = 1.5d0
C      ei(4) = -1.5d0
C      do  10  i = 1, nptdif
C        sloi(i) = slo + dele*ei(i)*ddde*val/rmax
C        ei(i) = e + dele*ei(i)
C        call rseq(eb1,eb2,ei(i),tol,z,l,nn,val,sloi(i),v,gp(1,i),
C     .            sum1,a,b,rofi,nr,nre)
C        vali(i) = val/dsqrt(sum1)
C        sloi(i) = sloi(i)/dsqrt(sum1)
C   10 continue
C      de1  = (ei(1) - ei(2))/2
C      del1 = (ei(1) + ei(2))/2 - e
C      de2  = (ei(3) - ei(4))/2
C      del2 = (ei(3) + ei(4))/2 - e
CC     Energy derivatives of value and slope
C      call dfphi(de1,del1,de2,del2,1,val,vali,nptdif.eq.4)
C      call dfphi(de1,del1,de2,del2,1,slo,sloi,nptdif.eq.4)
C
C      call dfphi(de1,del1,de2,del2,2*nr,g,gp,nptdif.eq.4)
C      call gintsr(gp,gp,a,b,nr,z,e,l,v,rofi,p)
C
CC ... Get phi,dphi from val,slo = (r*phi),(r*phi)' at rmax
C      phi = val/rmax
C      dphi = (slo - phi)/rmax
C      phip = vali(1)/rmax
C      dphip = (sloi(1) - phip)/rmax
C
C      if (iprint() .ge. 111) print 749, phi,dphi,phip,dphip
C  749 format(' PHIDOT:  phi,phip,phip,dphip=',4f12.6)
C      end
      subroutine gintsr(g1,g2,a,b,nr,z,e,l,v,rofi,sum)

C- Integrate inner product of two wave equations
C ----------------------------------------------------------------
Ci   g1,g2 :First and second radial wave functions
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   z     :nuclear charge
Ci   e     :energy
Ci   l     :l quantum number of g1,g2
Ci   v     :spherical potential
Ci   rofi  :radial mesh points
Co Outputs:
Co   sum   :inner product
Cr Remarks:
Cr   Uses Simpson's rule
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nr,l
      double precision a,b,z,e,sum,g1(nr,2),g2(nr,2),v(nr),rofi(nr)
C ... Local parameters
      integer i,ir
      double precision fllp1,c,r
      double precision tmc,fi
C     Speed of light, or infinity in nonrelativistic case
      common /cc/ c
      tmc(i,r) = c - (v(i) - 2d0*z/r - e)/c
      fi(i,r) = (r+b)*(g1(i,1)*g2(i,1)*(1 + fllp1/(tmc(i,r)*r)**2)
     .+ g1(i,2)*g2(i,2))

      fllp1 = l*(l+1)
      sum = 0d0
      do  10  ir = 2, nr-1, 2
        sum = sum + fi(ir,rofi(ir))
   10 continue
      sum = 2*sum
      do  11  ir = 3, nr-2, 2
        sum = sum + fi(ir,rofi(ir))
   11 continue
      sum = (2*sum + fi(nr,rofi(nr)))*a/3

      end subroutine gintsr

      subroutine gintsl(g1,g2,a,b,nr,rofi,sum)

C- Integrate inner product of two wave equations, large component only
C ----------------------------------------------------------------------
Ci Inputs
Ci   g1    :first wave function
Ci   g2    :second wave function
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Co Outputs
Co   sum:   inner product
Cr Remarks
Cr   Like gintsr, but uses large component only (corresponds to c->infty)
Cu Updates
Cu   20 Apr 01
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nr
      double precision a,b,g1(nr),g2(nr),rofi(nr)
C ... Local parameters
      integer ir
      double precision sum,r

      sum = 0
      do  10  ir = 2, nr-1,2
        r = rofi(ir)
        sum = sum+(r+b)*(g1(ir)*g2(ir))
   10 continue
      sum = sum+sum
      do  11  ir = 3, nr-2, 2
        r = rofi(ir)
        sum = sum+(r+b)*(g1(ir)*g2(ir))
   11 continue
      sum = sum+sum
      r = rofi(nr)
      sum = sum+(r+b)*(g1(nr)*g2(nr))
      sum = sum*a/3
      end subroutine gintsl

      subroutine asprjq(mode,clabl,nl,nsp,eula,neul,pnu,qnu,
     .pnuloc,qnuloc,bhat,amom)

C- Find average magnetization axis and project ASA moments onto it
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit affects bhat
Ci         :0 use bhat as input
Ci         :1 Print out mag. dir. from eula and qnu; make amom
Ci         :2 make bhat from eula and qnu; make amom
Ci         :3 combination of 1+2
Ci         :10s digit affects qnuloc,pnuloc
Ci         :0 do nothing to qnuloc,pnuloc
Ci         :1 copy pnu->pnuloc and qnu->qnuloc
Ci         :2 copy and rotate to bhat coordinates
Ci         :10s digit affects B,qnuloc when magnetization < 0
Ci         :1 scale B by -1 when
Ci   mode  :0 do nothing; just return
Ci   clabl :class name (for printout only)
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   eula  :Euler angles for noncollinear spins
Ci   neul  :1, nl, or nl**2 if Euler angles are: l-independent,
Ci         :l-dependent, or lm-dependent, respectively
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   qnu   :energy-weighted moments of the sphere charges
Cio Inputs/Outputs
Cio  bhat  :direction vector for B-field; see 1s digit mode
Co Outputs
Co   amom  :projection of magnetization onto bhat
Co   pnuloc:projection of pnu onto bhat
Co   qnuloc:projection of qnu onto bhat
Cr Remarks
Cr   If asprjq computes the B-field, it is taken to be parallel to
Cr   the average magnetization, which is computed by summing the
Cr   (vector) magnetization over all orbitals.
Cu Updates
Cu   06 Apr 04 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nl,nsp,neul
      double precision pnu(nl,nsp),qnu(3,nl,nsp),eula(neul,3),
     .amom,bhat(3),pnuloc(nl,nsp),qnuloc(3,nl,nsp)
      character clabl*8
C ... Local parameters
      logical lwrite
      integer i,k,l,m,ilm,stdo,lgunit,ipr,PRT1,PRT2,mode0,mode1,mode2
      double precision ql,al,alpha,beta,gamma,hatm(3),ploc,qloc(3),
     .dploc,aloc(3),dotp,ddot,dsqrt,rotm(3,3),sal
      character strn*5
      parameter (PRT1=40,PRT2=50)

      stdo = lgunit(1)
      call getpr(ipr)

C --- Construct the average magnetization (direction, amplitude) ---
      mode0 = mod(mode,10)
      mode2 = mod(mode/100,10)
      if (mode0 .ne. 0 .and. nsp .eq. 2) then
        lwrite = .true.
        if (mod(mode0,2) .eq. 0 .or. ipr .lt. PRT2) lwrite = .false.
        if (neul .le. 1) lwrite = .false.
        strn = '    l'
        if (neul .eq. nl*nl) strn = '  ilm'
        if (lwrite .and. mode0 .ge. 2) write(stdo,345) strn
        if (lwrite .and. mode0 .lt. 2) write(stdo,345) strn,'mhat.bxc'
  345   format(a,'     ql       mom',9x,'alpha     beta      gamma',
     .  19x,'mhat':17x,a)

        ilm = 0
        call dpzero(aloc,3)
        sal = 0
        do  l = 0, nl-1
          do   m = -l, l
            ilm = ilm+1
            if (neul .eq. nl) then
              alpha = eula(l+1,1)
              beta  = eula(l+1,2)
              gamma = eula(l+1,3)
            elseif (neul .eq. nl*nl) then
              alpha = eula(ilm,1)
              beta  = eula(ilm,2)
              gamma = eula(ilm,3)
            elseif (neul .eq. 1) then
              alpha = eula(1,1)
              beta  = eula(1,2)
              gamma = eula(1,3)
            else
              call rxi('atscpp: bad value neul=',neul)
            endif

C         Charge, magnetic moments, quantization axis for these angles
            ql = qnu(1,l+1,1)+qnu(1,l+1,2)
            al = qnu(1,l+1,1)-qnu(1,l+1,2)
            sal = sal + al/(2*l+1)
C         The local magnetization points along V = zhat(loc).
C         In global coordinates V = rotm^-1 zhat(loc) because
C         rotm*V = zhat(loc) when V points along M as
C         described in eua2rm.  V is then
C         hatm(1) = dcos(alpha)*dsin(beta)
C         hatm(2) = dsin(alpha)*dsin(beta)
C         hatm(3) = dcos(beta)
            call eua2rm(alpha,beta,gamma,rotm)
            hatm(1) = rotm(3,1)
            hatm(2) = rotm(3,2)
            hatm(3) = rotm(3,3)

C         Add to total magnetization
            call daxpy(3,al/(2*l+1),hatm,1,aloc,1)

C         Printout
            if (neul .eq. nl**2 .and. lwrite .and. mode0 .ge. 2) then
              write(stdo,'(i5,2f10.6,3x,3f10.6,3x,3f10.6,3x,f10.6)')
     .        ilm, ql/(2*l+1), al/(2*l+1), alpha, beta, gamma, hatm
            elseif (neul .eq. nl**2 .and. lwrite .and. mode0 .lt. 2) then
              write(stdo,'(i5,2f10.6,3x,3f10.6,3x,3f10.6,3x,f10.6)')
     .        ilm, ql/(2*l+1), al/(2*l+1), alpha, beta, gamma, hatm,
     .        ddot(3,hatm,1,bhat,1)
            elseif (neul .eq. nl .and. lwrite .and. mode0 .ge. 2) then
              write(stdo,'(i5,2f10.6,3x,3f10.6,3x,3f10.6,3x,f10.6)')
     .        l, ql, al, alpha, beta, gamma, hatm
              lwrite = .false.
            elseif (neul .eq. nl .and. lwrite .and. mode0 .lt. 2) then
              write(stdo,'(i5,2f10.6,3x,3f10.6,3x,3f10.6,3x,f10.6)')
     .        l, ql, al, alpha, beta, gamma, hatm, ddot(3,hatm,1,bhat,1)
              lwrite = .false.
            endif

          enddo
          lwrite = .true.
          if (mod(mode0,2) .eq. 0 .or. ipr .lt. PRT2) lwrite = .false.
        enddo
        amom = dsqrt(ddot(3,aloc,1,aloc,1))

        if (amom .ne. 0) then
C       Assign bhat to point along magnetization direction
          if (mode0 .ge. 2) call dpcopy(aloc,bhat,1,3,1/amom)

C       If sum moments < 0, optionally reverse B
          if (sal .lt. 0 .and. mode2 .eq. 1 .and. mode0 .ge. 2) then
            if (mode0 .ge. 2) call dpcopy(aloc,bhat,1,3,-1/amom)
          endif

C       Printout of bhat, moment, angle
          call info5(PRT1,0,0,' ATOM='//clabl//
     .    '%a  bhat=%3;10,6D  |M|=%;6,6d  bhat.M/|M|=%;6d',bhat,
     .    amom,ddot(3,aloc,1,bhat,1)/amom,0,0)

        else
          bhat(1) = 0
          bhat(2) = 0
          bhat(3) = 1
        endif

C     End of block to contruct bhat
      endif

C --- Rotate the qnu along projection bhat ---
      mode1 = mod(mode/10,10)
      if (mode1 .eq. 0) return

      if (nsp .eq. 1 .or. mode1 .eq. 1 .or.
     .ddot(3,bhat,1,bhat,1) .eq. 0d0) then
        call dcopy(3*nl*nsp,qnu,1,qnuloc,1)
        call dcopy(1*nl*nsp,pnu,1,pnuloc,1)
        return
      endif

      call dpzero(qnuloc,3*nl*nsp)
      ilm = 0
      do  l = 0, nl-1
        call dpzero(aloc,3)
        dploc = 0
        do   m = -l, l
          ilm = ilm+1
          if (neul .eq. nl) then
            alpha = eula(l+1,1)
            beta  = eula(l+1,2)
C           gamma = eula(l+1,3)
          elseif (neul .eq. nl*nl) then
            alpha = eula(ilm,1)
            beta  = eula(ilm,2)
C           gamma = eula(ilm,3)
          elseif (neul .eq. 1) then
            alpha = eula(1,1)
            beta  = eula(1,2)
C           gamma = eula(1,3)
          else
            call rxi('atscpp: bad value neul=',neul)
          endif

C         Charge, magnetic moments, quantization axis for these angles
          hatm(1) = dcos(alpha)*dsin(beta)
          hatm(2) = dsin(alpha)*dsin(beta)
          hatm(3) = dcos(beta)
          dotp = hatm(1)*bhat(1) + hatm(2)*bhat(2) + hatm(3)*bhat(3)

          do  i = 1, 3
            aloc(i) = aloc(i) + (qnu(i,l+1,1)-qnu(i,l+1,2))/(2*l+1)*dotp
          enddo
          dploc = dploc + (pnu(l+1,1)-pnu(l+1,2))/(2*l+1)*dotp

        enddo

        do  i = 1, 3
          qloc(i) = qnu(i,l+1,1) + qnu(i,l+1,2)
C         aloc(i) = qnu(i,l+1,1) - qnu(i,l+1,2)
          qnuloc(i,l+1,1) = (qloc(i) + aloc(i))/2
          qnuloc(i,l+1,2) = (qloc(i) - aloc(i))/2
        enddo
        ploc = pnu(l+1,1) + pnu(l+1,2)
        pnuloc(l+1,1) = (ploc + dploc)/2
        pnuloc(l+1,2) = (ploc - dploc)/2

      enddo

      if (ipr .ge. PRT2) then
        write(stdo,'(''  l isp'',19x,''qnu'',37x,''qloc'')')
        do  k = 1, 2
          do  l = 0, nl-1
            write(stdo,'(2i3,1x,3f13.7:1x,3f13.7)')
     .      l,k,(qnu(i,l+1,k),i=1,3),(qnuloc(i,l+1,k),i=1,3)
          enddo
        enddo
      endif

      end subroutine asprjq


