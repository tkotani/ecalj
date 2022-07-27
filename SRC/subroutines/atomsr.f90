subroutine atomsc(lgdd,nl,nsp,lmax,z,rhozbk,kcor,lcor,qcor,rmax,a, &
     nr,rofi,ec,ev,pnu,qnu,idmod,v,dv,rhoin,rho,rhoc,nrmix,qc,sumec, &
     sumtc,sumev,ekin,utot,rhoeps,etot,amgm,rhrmx,vrmax,qtot,exrmax, &
     job,niter,lfrz,plplus,qlplus,nmcore,qelectron,vsum)
  use m_lmfinit,only: lrel
  use m_getqvc
  use m_lgunit,only: stdo
  use m_amix,only: amix
  !- Makes an atomic sphere self-consistent and get atomic charges
  ! ----------------------------------------------------------------
  !i Inputs
  !i   lgdd  :T  add q2 phi phidd into the density
  !i         :F  add q2 <phi phidd> phi phi into the density
  !i         :NB: both produce the same integrated density.
  !i         :Which one should be used depends on the context.
  !i   nl    :leading dimension of pnu,qnu.
  !i         :Also, total charge cutoff; see Remarks
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   lmax  :maximum l for this site (but see Remarks)
  !i   z     :nuclear charge
  !i   rhozbk:constant nuclear background density (jellium) added to z
  !i   kcor  :(partial core occupation) p.q.n for occupation
  !i   lcor  :(partial core occupation) l quantum for occupation
  !i   qcor  :(partial core occupation) core charge and moment
  !i   rmax  :potential, density calculated to rmax
  !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   b     :                 -//-
  !i   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,,
  !i         :pnu = .5 - atan(dnu)/pi + (princ.quant.number).
  !i   qnu   :energy moments of charge (see Remarks)
  !i   idmod :0,1 or 2, specifing how the enu is set for an l-channel
  ! i   pz   :pnu for second panel (npan=2)
  ! i   qz   :qnu for second panel (npan=2)
  ! i   idmoz:idmod for second panel; specifies how semicore state
  ! i        :is treated
  ! i   npan :>1 for two-panel calculations
  !i   v     :spherical potential (job='pot', else input v not used)
  !i   dv    :constant potential shift added to sphere.
  !i         :Not used now.
  !i   rhoin:input density (job='rho', else input rhoin not used)
  !i         :used internally as a work array
  !i   nrmix :number of prior densities to mix to accelerate convergence
  !i         :to self-consistency (Anderson mixing); see Remarks.
  !i         :nrmix<0, linear mixing, with mixing beta |nrmix/100|
  !i   job   :(see Remarks)
  !i         :job='pot': start with potential
  !i         :job='rho': start with rhoin.
  !i         :job='gue': start with guessed charge density.
  !i   niter :number of iterations to attempt convergence to
  !i          self-consistency (see Remarks)
  !i   lfrz  :0 for soft core, >0 for frozen core
  ! o Input/Outputs
  ! o  ec    :core eigenvalues.  On input, these are guessed values.
  ! o        :if input ec(1) = 0, atomsc makes an internal choice fo ec.
  ! o  ev    :valence eigenvalues  On input, these are guessed values.
  ! o        :if input ec(1) = 0, atomsc makes an internal choice for ev.
  !o Outputs
  !o   rofi  :dimensioned (nr,2).
  !o         :(*,1) radial mesh points
  !o         :(*,2) weights
  !o   rho   :spherical charge density = 4 pi r**2 rhotrue
  !o   rhoc  :core charge density (unchanged if lfrz=1)
  !o   qc:   :core electronic charge
  !o   sumec :core single particle energy
  !o   sumtc :core kinetic energy (unchanged if lfrz=1)
  !o   ekin  :total kinetic energy
  !o   utot  :electrostatic energy
  !o   rhoeps:exchange-correlation energy
  !o   etot  :sphere total energy
  !o   amgm  :difference between spin up and spin down charge
  !o   rhrmx :true density at rmax
  !o   vrmax :true potential at rmax
  !o   qtot  :net charge in sphere
  !o   exrmax:exchange-correlation energy at rmax
  !r Remarks
  !r   Boundary conditions pnu, moments qnu, and the electrostatic
  !r   potential at rmax are all that is required to uniquely determine
  !r   a self-consistent spherical charge and potential.  'job'
  !r   determines how the starting potential is made, but the final
  !r   potential  and density should be independent of the initial choice.
  !r
  !r   atomsc uses the boundary condition that the potential at rmax
  !r   is zero, regardless of the total net charge inside the sphere.
  !r   See subroutine madpot for discussion of how this choice affects
  !r   the Madelung energy.
  !r
  !r   Sphere total energy is sum of K.E., Hartree energy, XC energy:
  !r      etot = ekin + utot + rhoeps
  !r   The kinetic energy is computed via double-counting terms
  !r     ekin = sumev + sumec + dsumec - rhovh - rhomu
  !b Bugs
  !b   Total energy terms need to be cleaned up and simplified.
  ! ----------------------------------------------------------------
  use m_ftox
  use m_vxcatom,only: vxc0sp
  implicit none
  logical :: lfrz,lgdd
  character job*3
  integer :: nr,nsp,nl,nrmix,niter,kcor,lcor,ncmx,nvmx,lmax, &
       idmod(0:nl-1)
  parameter (ncmx=50, nvmx=20)
  double precision :: ec(ncmx),ev(nvmx),rofi(nr,2),v(nr,nsp), &
       rho(nr,nsp),rhoc(nr,nsp),rhoin(nr,nsp),pnu(nl,2), &
       qnu(3,nl,nsp),z,rmax,a,qc,vrmax(2),exrmax(2),dv,rhrmx, &
       rhozbk,qcor(2),amgm,ekin,etot,qtot,rhoeps,sumec,sumev,sumtc, &
       utot
  logical :: last,cmdopt,ltmp
  integer :: iprint,ir,isp,iter,jmix,nglob,ipr1,l,ii
  character strn*10
  double precision :: b,ddot,decay,dl,dold,drdi,drho,dsumec,ea,fac,pi, &
       rho0t,rhomu,rhovh,rmsdel,ro,sum,tolrsq,vhrmax,vnucl,vrhoc,vsum, &
       zvnucl,rvh(2),rho0(2),reps(2),rmu(2),sec(2),stc(2),sev(2),dasum
  ! for Anderson mixing:
  integer:: nmix !amix , 
  real(8) ,allocatable :: rho_rv(:)
  real(8):: plplus(0:lmax,nsp),qlplus(0:lmax,nsp),qelectron
  double precision :: norm(10,10),awk(10,2),beta,beta1
  parameter (beta = 0.3d0)
  ! tolch is tolerance for change in the charge density, tolv for rseq
  double precision :: tolch,tl
  parameter (tolch=5d-5,tolrsq=1d-12)
!  external lgunit
  integer:: nmcore
  print *,'atomsc nmcore=',nmcore
  if (lmax >= nl) call rx('atomsc:  lmax too large')
  !      lrel = globalvariables%lrel
  pi = 4d0*datan(1d0)
  b = rmax/(dexp(a*(nr-1)) - 1)
  vhrmax = 0d0
  nmix = min(max(nrmix,0),10)
  sec(2) = 0d0
  sec(1) = 0d0
  stc(2) = 0d0
  stc(1) = 0d0
  allocate(rho_rv(nr*nsp*2*(nmix+2)))
  ! --- Core charge, radial mesh points and weights ---
  if (kcor /= 0) then
     if (qcor(1) /= 0 .OR. qcor(2) /= 0) then
        write(stdo,ftox)'Add core hole:  kcor=',kcor,' lcor=',lcor,' qcor=',ftof(qcor),'amom=',ftof(qcor(2))
     endif
  endif
  call getqvc(nsp,nl,lmax,z,pnu,qnu,0,0,kcor,lcor,qcor,qc,qtot,amgm)
  ! --- Guesses for core and valence eigenvalues ---
  ! akao probably not needed here
  if (ec(1) == 0) &
       call getqvc(nsp,nl,lmax,z,pnu,qnu,ncmx,nvmx,kcor,lcor,qcor,qc, &
       qtot,amgm,ec,ev) !only for core number check since we set ev and ev below.
  ! ill be replace by simple code.
  call radmsh(rmax,a,nr,rofi)
  call radwgt(rmax,a,nr,rofi(1,2))
  !! initialize ev here(overide getqvc now).feb2011
  ev = -0.5d0
  ec = -5.0d0
  !! takao2012jun. The followings are confusing when we use PZ and P together; because of historical reason, ql is used for PZ and so on.Need to improved...
  !$$$C --- Moments printout ---
  !$$$      if (iprint() .ge. 20) then
  !$$$        ltmp=dasum(lmax+1,qnu(2,1,1),3)+dasum(lmax+1,qnu(3,1,1),3).eq.0
  !$$$        if (ltmp) then
  !$$$          call awrit5('  Pl= %n:-1,1;5#8d%?!n==2! spn 2  %n:-1,1;5#8d',
  !$$$     .    ' ',180,stdo,lmax+1,pnu,nsp,lmax+1,pnu(1,2))
  !$$$          call dcopy(lmax+1,qnu,3,awk,1)
  !$$$          call dcopy(lmax+1,qnu(1,1,nsp),3,awk(1,2),1)
  !$$$          call awrit5('  Ql= %n:-1,1;5#8d%?!n==2! spn 2  %n:-1,1;5#8d',
  !$$$     .    ' ',180,stdo,lmax+1,awk,nsp,lmax+1,awk(1,2))
  !$$$        elseif (iprint() .ge. 30) then
  !$$$          write(stdo,891)
  !$$$  891     format(
  !$$$     .    '   l',8x,'pl',11x,'q0',11x,'q1',11x,'q2',5x,' id ',6x,'dl')
  !$$$          do  11  isp = 1, nsp
  !$$$            do  11  l = 0, lmax
  !$$$              dl = dtan(pi*(.5d0 - pnu(l+1,isp)))
  !$$$              if (dabs(dl) .gt. 9999) dl = 0
  !$$$              write(stdo,100) l,pnu(l+1,isp),(qnu(ii,l+1,isp),ii=1,3),
  !$$$     .        idmod(l),dl
  !$$$  100         format(i4,4f13.7,i4,f13.7)
  !$$$   11     continue
  !$$$        endif
  !$$$      endif

  ! --- Initial charge density ---
  if (job == 'pot') then
     call newrho(z,lrel,lgdd,nl,1,lmax,a,b,nr,rofi,v,rhoin,rhoc,kcor, &
          lcor,qcor,pnu,qnu,sec,stc,sev,ec,ev,tolrsq,nsp,lfrz,000,plplus,qlplus,nmcore)
     if (niter == 0) then
        call dcopy(nr*nsp,rhoin,1,rho,1)
        ! i#error, have return with len(w_varlist)>0 at line 189
        if (allocated(rho_rv)) deallocate(rho_rv)
        return

     endif
  else if (job == 'gue') then
     decay = 1d0+z/10d0
     decay = dmin1(decay,5d0)
     decay = 5
     !         write(6,335) decay
     !  335   format(/' initialize exponential density with decay',f7.3)
     sum = 0d0
     do  4  ir = 1, nr
        ro = dexp(-decay*rofi(ir,1))*rofi(ir,1)**2
        rhoin(ir,1) = ro
        sum = sum+ro*a*(rofi(ir,1)+b)
4    enddo
     fac = z/(sum*nsp)
     do  5  ir = 1, nr
        rhoin(ir,1) = rhoin(ir,1)*fac
5    enddo
     if (nsp == 2) then
        do  6  ir = 1, nr
           rhoin(ir,2) = rhoin(ir,1)
6       enddo
     endif
  else if (job /= 'rho') then
     call rx('atomsc: job not pot|rho|gue')
  endif
!  call dpscop ( rhoin , rho_rv , nr * nsp , 1 , 1 + nr * nsp* ( nmix + 2 ) , 1d0 )
  rho_rv(1+nr*nsp*(nmix+2):nr*nsp+nr*nsp*(nmix+2))= reshape(rhoin(1:nr,1:nsp),[nr*nsp])
  ! --- Start self-consistency loop ---
  drho = 100d0
  last = .false.
  if (iprint() >= 30) write(stdo,341)
  jmix = 0
  dold = 1
  beta1 = beta
  vrhoc=-1d99
  if (nrmix < 0) beta1 = dble(-nrmix)/100
  do  35  iter = 1, niter
     tl = tolrsq
     call addzbk(rofi,nr,nsp,rhoin,rhozbk,-1d0)
     if(abs(rmax-rofi(nr,1))>1d-6) call rx('atomsr.F:something wrong. abs(rmax-rofi(nr,1))>1d-3')
     vhrmax=2d0*(qelectron-z)/rmax
     call poiss0(z,a,b,rofi,rhoin,nr,vhrmax,v,rvh,vsum,nsp) !  Hartree potential
     vsum = vsum + 4d0*pi*(z-qelectron)*rmax**2
     call addzbk(rofi,nr,nsp,rhoin,rhozbk,1d0)
     vnucl = v(1,1)  
     if (last .AND. iprint() >= 50) call pshpr(80)
     call vxc0sp(a,b,rofi,rhoin,nr,v,rho0,reps,rmu,nsp,exrmax(1))! Exchange-correlation potential
     if (last .AND. iprint() >= 50) call poppr
     fac = 4*pi*rofi(nr,1)**2
     rhrmx = rhoin(nr,1)/fac
     if (nsp == 2) then !       Get rhrmx, exrmax
        exrmax(2) = exrmax(1)
        rhrmx = rhrmx + rhoin(nr,2)/fac
     endif
     ipr1 = 0
     if (last .AND. iprint() >= 40) ipr1 = 1
     if (last .AND. iprint() > 40) ipr1 = 2
     if( .NOT. last) call pshpr(15) !low print()
     call newrho(z,lrel,lgdd,nl,1,lmax,a,b,nr,rofi,v,rho,rhoc, &
          kcor,lcor,qcor,pnu,qnu,sec,stc,sev,ec,ev,tl,nsp,lfrz,ipr1,plplus,qlplus,nmcore)
     if( .NOT. last) call poppr !set back to original print()
     drho = 0d0
     sum = 0d0
     vrhoc = 0d0
     rho0t = 0d0
     do  40  isp = 1, nsp
        rho0t = rho0t + rho0(isp)
        sum = sum + ddot(nr,rofi(1,2),1,rho(1,isp),1)
        do  42  ir = 1, nr
           drdi = a*(rofi(ir,1) + b)
           drho = drho + rofi(ir,2)/drdi*dabs(rho(ir,isp)-rhoin(ir,isp))
42      enddo
        do  41  ir = 2, nr
           ea = (v(ir,isp)-2*z/rofi(ir,1))
           vrhoc = vrhoc + rofi(ir,2)*ea*rhoc(ir,isp)
41      enddo
40   enddo
!     call dcopy ( nr * nsp , rho , 1 , rho_rv , 1 )
     rho_rv(1:nr*nsp)=reshape(rho(1:nr,1:nsp),[nr*nsp])
     jmix = amix ( nr * nsp , min ( jmix , nmix ) , nmix , 0 , beta1 &
          , iprint ( ) - 70 , .9d0 ,  rho_rv , awk , rmsdel )
     !call dpscop ( rho_rv , rhoin , nr * nsp , 1 + nr * nsp * (nmix + 2 ) , 1 , 1d0 )
     rhoin(1:nr,1:nsp) = reshape(rho_rv(1+nr*nsp*(nmix+2):nr*nsp+nr*nsp*(nmix+2)),[nr,nsp])
     if (last) goto 90
     if (iprint() >= 41 .OR. iprint() >= 30 .AND. &
          (drho < tolch .OR. iter == niter-1 .OR. iter == 1)) &
          write(stdo,340) iter,sum,drho,vnucl,rho0t,vsum,beta1
340  format(i5,f12.6,1p,e12.3,0p,f14.4,e14.4,f14.4,f7.2)
341  format(/'  iter     qint',9x,'drho',10x,'vh0',10x,'rho0',10x,'vsum',5x,'beta')
     last = (drho .lt. tolch .or. iter .eq. niter-1)
     jmix = jmix+1
     !       Beta for next iteration
     beta1 = min(max((1-drho/dold)/beta1,beta1-.2d0,beta), 1d0,beta1+.2d0)
     if (nmix > 0 .AND. drho < 1) beta1 = 1
     if (nrmix < 0) beta1 = dble(-nrmix)/100
     dold = drho
35 enddo
90 continue
  if (allocated(rho_rv)) deallocate(rho_rv)
  if (iprint() >= 30) write(stdo,'(1x)')
  ! --- Collect terms for total energy ---
  rhoeps = 0d0
  rhomu  = 0d0
  sumev  = 0d0
  if ( .NOT. lfrz) then
     sumec = 0d0
     sumtc = 0d0
  endif
  rhovh  = 0d0
  do  80  isp = 1, nsp
     if (nsp == 2 .AND. iprint() > 30) &
          write(stdo,230) isp,v(nr,isp)-2*z/rmax, &
          sev(isp),sec(isp),rvh(isp),reps(isp),rmu(isp)
230  format(' Spin',i2,':', &
          /' vrmax=',f12.5,'    sumev= ',f12.5,'    sumec=',f12.5, &
          /' rhovh=',f12.5,'    rhoeps=',f12.5,'    rhomu=',f12.5)
     rhoeps = rhoeps + reps(isp)
     rhomu = rhomu + rmu(isp)
     sumev = sumev + sev(isp)
     if ( .NOT. lfrz) then
        sumec = sumec + sec(isp)
        sumtc = sumtc + stc(isp)
     endif
     rhovh = rhovh + rvh(isp)
80 enddo
  zvnucl = -z*vnucl
  utot = .5d0*(rhovh + zvnucl)
  ! ... Correction to core eigenvalues if sumec not obtained from this V
  dsumec = vrhoc - (sumec-sumtc)
  ekin = sumev + sumec + dsumec - rhovh - rhomu
  etot = ekin + utot + rhoeps
  if (iprint() >= 40) write(stdo,139) sumev,sumec,vnucl,rhovh, &
       zvnucl,utot,rhomu,rhoeps,dsumec,ekin,sumtc,etot
139 format(/' sumev=',f13.6,'    sumec =',f13.6,'   vnucl =',f13.6 &
       /' rhovh=',f13.6,'    zvnucl=',f13.6,'   utot  =',f13.6 &
       /' rhomu=',f13.6,'    rhoeps=',f13.6,'   dsumec=',f13.6 &
       /' ekin= ',f13.6,'    tcore =',f13.6,'   etot  =',f13.6)
  vrmax(1) = -2*z/rmax
  do  55  isp = 1, nsp
     vrmax(1) = vrmax(1) + v(nr,isp)/nsp
55 enddo
  vrmax(2) = 0d0
  if (nsp == 2) vrmax(2) = v(nr,1)-v(nr,2)
end subroutine atomsc

subroutine addzbk(rofi,nr,nsp,rho,rhozbk,scale)
  !     implicit none
  integer :: nr,nsp
  double precision :: rofi(nr),rho(nr,*),rhozbk,scale
  integer :: ir,isp
  double precision :: s
  if (rhozbk == 0) return
  s = 16d0*datan(1d0)*scale*rhozbk
  do isp = 1, nsp
     rho(:,isp) = rho(:,isp) + s*rofi(:)**2
  enddo
end subroutine addzbk

subroutine poiss0(z,a,b,rofi,rho,nr,vhrmax,v,rhovh,vsum,nsp)
  !- Hartree potential for spherical rho.
  !  ---------------------------------------------------------------------
  !i Inputs:
  !i   z     :nuclear charge
  !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   rho   :spherical charge density times 4*pi*r*r
  !i   nr    :number of mesh points
  !i   vhrmax:value of v(nr)
  !i   nsp   :=1 spin degenerate, =2 non-degenerate
  !o Outputs:
  !o   rofi  :radial mesh points
  !o   v     :spherical Hartree potential
  !o   vsum  :integral over that potential which is zero at rmax.
  !o   rhovh :integral of rho*vh, where vh is the electrostatic
  !o         :potential from both electronic and nuclear charge
  !r Remarks:
  !r    Solves Poisson's equation for given spherical charge density
  !r    and a specified value vhrmax at rofi(nr).
  !r    rho =  4*pi*r*r*rhotrue  but  v = vtrue
  ! ----------------------------------------------------------------------
  implicit double precision (a-h,o-z)
  integer :: nr,nsp
  double precision :: rofi(nr),v(nr,nsp),rho(nr,nsp),rhovh(nsp)
  double precision :: pi,ea,rpb,a,b,z
  double precision :: rmax,r2,r3,r4,f2,f3,f4,x23,x34,cc,bb,dd,df,drdi, &
       r,srdrdi,g,f,y2,y3,y4,ro,vnow,wgt,a2b4,vhrmax,vsum,vhat0
  integer :: ir,isp
  real(8):: qtotxxx
  pi = 4d0*datan(1d0)
  ea = dexp(a)
  rpb = b
  do  5  ir = 1, nr
     rofi(ir) = rpb - b
     rpb = rpb*ea
5 enddo
  rmax = rofi(nr)
  ! --- Approximate rho/r**2 by cc + bb*r + dd*r*r near zero  ---
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
75 enddo
  x23 = (r3*r3*f2 - r2*r2*f3)/(r3 - r2)
  x34 = (r4*r4*f3 - r3*r3*f4)/(r4 - r3)
  cc = (r2*x34 - r4*x23) / (r3*(r2 - r4))
  bb = ((r2+r3)*x34 - (r3+r4)*x23) / (r3*r3*(r4-r2))
  dd = (f2 - bb*r2 - cc)/r2**2

  ! --- Numerov for inhomogeneous solution ---
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
     if (ir == 2) y2 = -2d0*f2*r2*drdi*srdrdi
     if (ir == 3) y3 = -2d0*f3*r3*drdi*srdrdi
     df = f - df
12 enddo
  do 13  ir = 4, nr
     r = rofi(ir)
     drdi = a*(r + b)
     srdrdi = dsqrt(drdi)
     ro = 0d0
     do  76  isp = 1, nsp
        ro = ro + rho(ir,isp)
76   enddo
     y4 = -2d0*drdi*srdrdi*ro/r
     df = df + g*a2b4 + (y4 + 10d0*y3 + y2)/12d0
     f = f + df
     g = f/(1d0 - a2b4/12d0)
     v(ir,1) = g*srdrdi/r
     y2 = y3
     y3 = y4
13 enddo
  ! --- Add constant to get v(nr) = vhrmax ---
  vnow = v(nr,1) - 2d0*z/rmax
  do  27  ir = 1, nr
     v(ir,1) = v(ir,1) + (vhrmax - vnow)
27 enddo
  ! --- Integral vh and  rho * vh ---
  vhint=0d0
  rhovh(1) = 0d0
  rhovh(nsp) = 0d0
  vsum = 0d0
  vhat0 = 0d0
  do  30  ir = 2, nr
     r = rofi(ir)
     drdi = a*(r + b)
     wgt = 2*(mod(ir+1,2) + 1)/3d0
     if (ir == nr) wgt = 1d0/3d0
     ro = 0d0
     do  31  isp = 1, nsp
        rhovh(isp) =  rhovh(isp) &
             + wgt*drdi*rho(ir,isp)*(v(ir,1) - 2d0*z/r)
        ro = ro + rho(ir,isp)
31   enddo
     vhat0 = vhat0 + wgt*drdi*ro*(1d0/r - 1d0/rmax)
     vsum = vsum + wgt*drdi*r*r*v(ir,1)
30 enddo
  !      vsum = 4d0*pi*(vsum - z*rmax*rmax)
  vsum = 4d0*pi*(vsum - z*rmax*rmax)
  vhat0 = 2d0*vhat0 + 2d0*z/rmax + vhrmax
  v(1,1) = vhat0
  ! --- Copy to second spin channel if spin polarized ---
  if (nsp == 1) return
  do  84  ir = 1, nr
     v(ir,2) = v(ir,1)
84 enddo
end subroutine poiss0

subroutine fctp0(l,nr,rofi,v,z,nctp0)
  !- Initialize things for FCTP, which finds classical turning point
  ! ----------------------------------------------------------------------
  !i Inputs:
  !i   l     :angular momentum
  !i   rofi  :radial mesh points
  !i   v     :spherical potential
  !i   z     :nuclear charge
  !i   nr    :number of mesh points
  !o Outputs:
  !o   nctp0 :minimum of effective potential
  !o          or nr if v(nr) > vmin + 3 Ry
  ! ----------------------------------------------------------------------
  !     implicit none
  ! Passed variables:
  integer :: l,nr,nctp0
  double precision :: z,v(nr),rofi(nr)
  ! Local variables:
  integer :: ir,irmin
  double precision :: veff,zz,fllp1,vi,vim1
  parameter (irmin=11)
  ! Statement functions:
  veff(ir)=fllp1/(rofi(ir)*rofi(ir))-zz/rofi(ir)+v(ir)
  zz = z+z
  fllp1 = l*(l+1)

  ir = irmin
  vi = veff(ir)
  vim1 = veff(ir-1)
10 if (vi <= vim1 .AND. ir < nr) then
     ir = ir+1
     vim1 = vi
     vi = veff(ir)
     goto 10
  endif
  nctp0 = ir-1
  if (veff(nctp0) >= veff(nr)-3.d0) nctp0 = nr
end subroutine fctp0

subroutine fctp(a,b,e,l,nctp0,nr,rofi,v,z,nctp)
  !- Finds classical turning point
  !  ---------------------------------------------------------------------
  !i Inputs:
  !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   b     :                 -//-
  !i   e     :energy
  !i   l     :angular momentum
  !i   nctp0 :minimum of effective potential (fctp0)
  !i   nr    :number of mesh points
  !i   rofi  :radial mesh points
  !i   v     :spherical potential
  !i   z     :nuclear charge
  !i Inputs/Outputs:
  ! o  nctp  :rofi(nctp) is classical turning point
  !u Updates
  !u   13 Jun 00 Added safety checks to guard against jumps in potential
  !  ---------------------------------------------------------------------
  !     implicit none
  ! Passed variables:
  integer :: nctp,nctp0,l,nr
  double precision :: e,rofi(nr),v(nr),z,a,b
  ! Local variables:
  integer :: ir,irep,n1,n2,nlast,ntry
  double precision :: r,veff,fllp1,fofr,dfdr,rtry,zz
  ! Intrinsic functions:
  intrinsic dlog,dmax1,min0
  ! Statement functions:
  veff(ir)=fllp1/(rofi(ir)*rofi(ir))-zz/rofi(ir)+v(ir)
  zz=z+z
  fllp1 = l*(l+1)
  if (nctp0 == nr .OR. e > veff(nr)) then
     nctp = nr
  elseif (e < veff(nctp0)) then
     nctp = 2
  else
     n1 = nctp0
     n2 = nr-1
     nlast = -10
     do  10  irep = 1, 20
        if (nctp > n2 .OR. nctp < n1) nctp = (n1 + n2)/2
        r = rofi(nctp)
        fofr =  veff(nctp)-e
        dfdr =-(fllp1+fllp1)/r/r/r + zz/r/r + &
             (v(nctp+1) - v(nctp-1))/(2.d0*a*(r+b))
        rtry = dmax1(r-fofr/dfdr,rofi(2))
        ntry = dlog(rtry/b+1.d0)/a + 1.5d0

        !         If there was a large change, check for safety
        if (nlast == nctp .AND. iabs(ntry-nctp) > 10) then
           if (fofr > 0) then
              do  ntry = nctp-1, nctp0+1, -1
                 fofr =  veff(ntry)-e
                 if (fofr < 0) goto 20
                 nctp = ntry
                 n2 = ntry
              enddo
20            continue
           else
              do  ntry = nctp+1, nr-1
                 fofr =  veff(ntry)-e
                 if (fofr > 0) goto 21
                 nctp = ntry
                 n1 = ntry
              enddo
21            continue
           endif
        endif
        !         Exit point
        if (nlast == nctp) then
           if (nctp == nctp0+1) nctp = 2
           return
        endif
        if (fofr > 0.d0)  then
           n2 = nctp
        else
           n1 = nctp
        endif
        nlast = nctp
        nctp = min0(ntry,nr-1)
10   enddo
     if (nctp == nctp0+1) nctp = 2
  endif
end subroutine fctp


subroutine newrho(z,lrel,lgdd,nl,nlr,lmax,a,b,nr,rofi,v,rho,rhoc, &
     kcor,lcor,qcor,pnu,qnu,sumec,sumtc,sumev,ec,ev,tol,nsp,lfrz,ipr,plplus,qlplus,nmcore)
  !! ev is dummy now. feb2010 takao
  !- Makes spherical charge density for a spherical potential.
  !  ---------------------------------------------------
  !i Inputs:
  !i   z     :nuclear charge
  !i   lrel  :0 for nonrelativistic, 1 for relativistic
  !i   lgdd  :T q2 is coefficient to phidot**2 + phi*phidotdot
  !i         :F q2 is coefficient to phidot**2 - p phi**2
  !i         :Both produce the same integrated density; see Remarks.
  !i         :lgdd=F follows Stuttgart conventions.
  !i   nl    :(global maximum l) + 1
  !i   nlr   :second dimension of rho:
  !r         :1 if spherical part of rho is to be generated
  !i         :nl if generated rho is to be decomposed by l
  !i         :In the latter case, the core is not calculated
  !i   lmax  :maximum l for a given site
  !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   nr    :number of radial mesh points
  !i   rofi  :radial mesh points
  !i   v     :spherical potential (atomsr.f)
  !i   kcor  :(partial core occupation) p.q.n for occupation
  !i   lcor  :(partial core occupation) l quantum for occupation
  !i   qcor  :(partial core occupation) core charge and moment
  !i   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,,
  !i         :pnu = .5 - atan(dnu)/pi + (princ.quant.number).
  !i   qnu   :energy moments of charge (see Remarks)
  !i   tol   :precision to which wave functions are integrated.
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   g     :normalized wave function times r
  !i   gp    :energy derivatives of g
  !i   lfrz  :T, do not make core rho
  !i   ipr   :0 no printout
  !i         :1 summary printout of core
  !i         :2 detailed printout of core
  !o Outputs:
  !o   rho   :spherical charge density times 4*pi*r*r
  !o   rhoc  :core density times 4*pi*r*r
  !o   sumec :sum of core eigenvalues
  !o   sumtc :core kinetic energy = sumec - v*rhoc
  !o   sumev :sum of valence eigenvalues
  !o   ec    :core eigenvalues
  !o   ev    :valence eigenvalues !not output now
  !r Remarks:
  !r   rho is determined by boundary conditions pnu and moments qnu.
  !r   For rmax>10 sets phidot, phidotdot are made zero.
  !r
  !r   Recall that val,slo correspond to u = r*phi, so
  !r     rmax * slo / val = D_u = D_phi + 1.
  !r     For val=rmax  slo = D + 1
  !r
  !r   Switch lgdd concerns the contribution of <phidot|phidot> to the
  !r   sphere charge.  The zeroth moment may be defined as the amount of
  !r   charge in channel l, i.e.
  !r     q^0 = (amount of <phi|phi>) +  p* (amount of <phidot|phidot>)
  !r   where
  !r     p = <phidot|phidot>
  !r   Then q^0 is not the amount of phi*phi inside the sphere, but rather
  !r   but instead, it is
  !r     (q_0 - p q_2)
  !r   since q2 is the amount of <phidot|phidot>.  The charge density is
  !r     rho= (q_0 - p q_2) phi*phi + 2 q_1 phi*phidot + q_2 phidot*phidot
  !r   This is the Stuttgart convention (lgdd=F)
  !r
  !r   Methfessel convention:  to avoid explicit dependence of rho on
  !r   p, he approximated p*phi*phi with -phi*phidotdot (they have the
  !r   same integrated charge).  Then
  !r     rho= q_0 phi*phi + 2 q_1 phi*phidot +
  !r          q_2 (phidot*phidot + phi*phidotdot)
  !  ---------------------------------------------------
  use m_rhocor,only:rhocor
  use m_phidx,only: phidx
  implicit none
  ! Passed parameters
  logical :: lgdd,lfrz
  integer :: nl,nlr,lmax,nr,nsp,lrel,kcor,lcor,ipr,iz
  double precision :: z,a,b,tol,sumev(nsp),sumec(nsp),sumtc(nsp),ec(*), &
       ev(*),qcor(2),v(nr,nsp),rofi(nr),rho(nr,nlr,nsp),rhoc(nr,nsp), &
       qnu(3,nl,nsp),pnu(nl,nsp)
  ! Local parameters
  logical :: free
  integer :: konfig(0:10),l,isp,ir,ival,nn,nre,jr,lmaxc,lr,k,nrmx
  parameter (nrmx=1501)
  double precision :: rocrit,pi,eb1,eb2,q0,q1,q2,rmax,eval,dl,val(5), &
       slo(5),sum,ro,phi,dphi,phip,dphip,p,fllp1,r,tmc,gfac,c,q00
  double precision :: g(2*nrmx),gp(2*nrmx*4)
  !     Speed of light, or infinity in nonrelativistic case
  common /cc/ c
  real(8):: plplus(0:lmax,nsp),qlplus(0:lmax,nsp)
  integer::nmcore !jun2012
  call setcc(lrel)
  !       call tcn('newrho')

  if (nr > nrmx) call rxi(' newrho: increase nrx, need',nr)
  lr = 1
  rocrit = 0.002d0/4
  pi = 4d0*datan(1d0)
  eb1 = -50d0
  eb2 =  50d0
  rmax = rofi(nr)
  free = (rmax .gt. 9.99d0)
  call config(pnu,lmax,z,konfig,lmaxc)
  if (kcor > 0) then
     lmaxc = max(lmaxc,lcor)
     konfig(lcor) = max(konfig(lcor),kcor+1)
  endif

  ! --- Calculate core density ---
  if (nlr == 1) then
     if ( .NOT. lfrz) then
        call dpzero(rhoc,nr*nsp)
        !          print *,'goto rhocor nmcore=',nmcore
        call rhocor(0,z,lmaxc,nsp,konfig,a,b,nr,rofi,v,g, &
             kcor,lcor,qcor,tol,ec,sumec,sumtc,rhoc,nmcore=nmcore,ipr=ipr)
     endif
     call dcopy(nr*nsp,rhoc,1,rho,1)
  endif

  ! --- Loop over valence states ---
  !      ival = 0
  !          if(iz==0) then
  !          ival = ival+1
  !          eval = ev(ival)
  !          endif
  eval=-0.5d0 !initial condition. the same as getqvc
  do  202  isp = 1, nsp
     sumev(isp) = 0d0
     do  201  l = 0, lmax
        do  20  iz=0,1  !takao feb2011
           if (nlr == nl) lr = l+1
           q0 = max(qnu(1,l+1,isp),0d0)
           q1 = qnu(2,l+1,isp)
           q2 = qnu(3,l+1,isp)
           if (q0 < 1d-6) goto 20
           nn = int(pnu(l+1,isp)) - l - 1
           val(1) = rmax
           dl = dtan(pi*(0.5d0 - pnu(l+1,isp)))
           slo(1) = dl + 1
           !! Setting when iz=1 override setting when iz=1
           if(iz==1)then
              if(plplus(l,isp)>0d0 .AND. qlplus(l,isp)>0d0) then
                 !!  note this is a case of Qv(lower principle quantum number).
                 !!  Search "=== Charge for l" in freeat.F, and NOTE: above it.
                 !             print *,'vvvvvv', plplus(l,isp)
                 q0=qlplus(l,isp)
                 q1=0d0
                 q2=0d0
                 nn=int(plplus(l,isp)) - l - 1
                 dl = dtan(pi*(0.5d0 - plplus(l,isp)))
                 slo(1) = dl + 1
              else
                 cycle
              endif
           endif

           ! ccccccccccccccccccccc
           !          print *,'tttt:l nn iz q0 dl =',l,nn,iz, q0,dl
           ! cccccccccccccccccccccc

           if (free) val(1) = 1d-30
           if (free) slo(1) = -val(1)
           call rseq(eb1,eb2,eval,tol,z,l,nn,val,slo,v(1,isp),g, &
                sum,a,b,rofi,nr,nre)
           !          ev(ival) = eval !takao think this is not necessary. ev is a dummy
           ! eb2010 just for initial condition.
           sumev(isp) = sumev(isp) + eval*q0 + q1
           ro = g(nr)**2
           if ( .NOT. free .AND. ro < rocrit) write(*,766) l,nn,nre,ro
766        format(' NEWRHO (warning): PHP,PHPP set to zero,l,nn,nre,rho=', &
                3i5,2f8.4)
           if (free .OR. ro < rocrit) then
              call dpzero(gp,8*nr)
              p = 0
           else
              val(1) = val(1)/dsqrt(sum)
              slo(1) = slo(1)/dsqrt(sum)
              !          call phidot(z,l,v(1,isp),eval,a,b,rofi,nr,g,val,slo,tol,
              !     .                nn,gp,phi,dphi,phip,dphip,p)
              call phidx(1,z,l,v(1,isp),0d0,0d0,rofi,nr,2,tol,eval,val, &
                   slo,nn,g,gp,phi,dphi,phip,dphip,p,0d0,[0d0],0d0,[0d0])
           endif
           fllp1 = l*(l+1)
           !  ...  Case add q2 phi phidd rho
           if (lgdd) then
              k = 2*nr
              do  21  ir = 2, nre
                 jr = ir + nr
                 r = rofi(ir)
                 tmc = c - (v(ir,isp) - 2d0*z/r - eval)/c
                 gfac = 1d0 + fllp1/(tmc*r)**2
                 rho(ir,lr,isp) =  rho(ir,lr,isp) + &
                      q0*(gfac*g(ir)**2 + g(jr)**2) + &
                      2*q1*(gfac*g(ir)*gp(ir) + g(jr)*gp(jr)) + &
                      q2*(gfac*(gp(ir)**2 + g(ir)*gp(ir+k)) + &
                      gp(jr)**2 + g(jr)*gp(jr+k))
21            enddo
              !  ...  Case add -p q2 phi phi into rho
           else
              q00 = q0-p*q2
              do  22  ir = 2, nre
                 jr = ir + nr
                 r = rofi(ir)
                 tmc = c - (v(ir,isp) - 2d0*z/r - eval)/c
                 gfac = 1d0 + fllp1/(tmc*r)**2
                 rho(ir,lr,isp) = rho(ir,lr,isp) + &
                      q00*(gfac*g(ir)**2 + g(jr)**2) + &
                      2*q1*(gfac*g(ir)*gp(ir) + g(jr)*gp(jr)) + &
                      q2*(gfac*gp(ir)**2 + gp(jr)**2)
22            enddo
           endif
20      enddo
201  enddo
202 enddo
  !      call tcx('newrho')
end subroutine newrho





subroutine gintsl(g1,g2,a,b,nr,rofi,sum)
  !- Integrate inner product of two wave equations, large component only
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   g1    :first wave function
  !i   g2    :second wave function
  !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   nr    :number of radial mesh points
  !i   rofi  :radial mesh points
  !o Outputs
  !o   sum:   inner product
  !r Remarks
  !r   Like gintsr, but uses large component only (corresponds to c->infty)
  !u Updates
  !u   20 Apr 01
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: nr
  double precision :: a,b,g1(nr),g2(nr),rofi(nr)
  ! ... Local parameters
  integer :: ir
  double precision :: sum,r
  sum = 0
  do  10  ir = 2, nr-1,2
     r = rofi(ir)
     sum = sum+(r+b)*(g1(ir)*g2(ir))
10 enddo
  sum = sum+sum
  do  11  ir = 3, nr-2, 2
     r = rofi(ir)
     sum = sum+(r+b)*(g1(ir)*g2(ir))
11 enddo
  sum = sum+sum
  r = rofi(nr)
  sum = sum+(r+b)*(g1(nr)*g2(nr))
  sum = sum*a/3
end subroutine gintsl

