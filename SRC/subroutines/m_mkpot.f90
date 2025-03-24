!>Get one-particle potential. See http://dx.doi.org/10.7566/JPSJ.84.034702
module m_mkpot !How to learng this? Instead of reading all source, understand I/O.
  use m_lgunit,only:stdo,stdl
  use m_lmfinit,only: z_i=>z,lmxa_i=>lmxa,lmxb_i=>lmxb,kmxt_i=>kmxt,nlmxlx
  use m_lmfinit,only: nbas,qbg=>zbak,ham_frzwf,lmaxu,nsp,nlibu,n0,nppn,lfrce, nchan=>pot_nlma, nvl=>pot_nlml
  use m_struc_def,only: s_rv1,s_cv1,s_sblock,s_rv4,s_cv5
  public:: m_mkpot_init, m_mkpot_energyterms, m_mkpot_novxc, m_mkpot_deallocate 
  ! Potential terms, call m_mkpot_init. Generated at mkpot-locpot-augmat-gaugm
  type(s_sblock),allocatable,protected,public  :: ohsozz(:,:),ohsopm(:,:) !SOC matrix
  type(s_cv5),allocatable,protected,public  :: oppi(:,:) !pi-integral Eq.(C.6)  in the JPSJ.84.034702
  type(s_rv4),allocatable,protected,public  :: otau(:,:) !tau            (C.5) 
  type(s_rv4),allocatable,protected,public  :: osig(:,:) !sigma          (C.4) 
  complex(8),allocatable,protected ,public  :: osmpot(:,:,:,:)!0th component of Eq.(34)
  
  real(8),allocatable,protected,public:: fes1_rv(:), fes2_rv(:) !force terms
  real(8),allocatable,protected,public:: hab_rv(:,:,:), sab_rv(:,:,:,:,:), qmom(:,:),vesrmt(:)
  real(8),protected,public:: qval,vconst,qsc
  real(8),allocatable,protected,public:: phzdphz(:,:,:,:) !val and slo at Rmt for local orbitals.
  
  ! Energy terms by call m_mkpot_energyterms
  real(8),protected,public:: utot,rhoexc,xcore,valvef,amom, valves,cpnves,rhovxc
  ! NoVxc terms  by call m_mkpot_novxc
  type(s_cv5),allocatable,protected,public  :: oppix(:,:)    !pi-integral without xc term
  complex(8),allocatable,protected ,public  :: spotx(:,:,:,:)!0th component of Eq.(34) without xc term
  private
  !! nov2021 dipole contribution added  (not working...)! oppixd,spotxd: dipole part to oppix,spotx
  !      type(s_cv1),allocatable,protected,public  :: oppixd(:,:,:)
  !      complex(8),allocatable,protected ,public  :: spotxd(:,:,:,:,:)
contains
  subroutine m_mkpot_novxc() ! outputs are oppix and spotx (for no vxc terms).
    use m_supot,only: n1,n2,n3
    use m_density,only: osmrho, orhoat !main input density
    logical:: novxc_
    real(8),allocatable :: fes1_rvx(:)
    type(s_sblock),allocatable:: ohsozzx(:,:),ohsopmx(:,:) !dummy for SOC
    type(s_rv4),allocatable   :: otaux(:,:) !dummy
    type(s_rv4),allocatable   :: osigx(:,:) !dummy
    write(stdo,"(a)")' m_mkpot_novxc: Making one-particle potential without XC part ...'
    allocate( vesrmt(nbas))
    allocate( qmom(nlmxlx,nbas)) !rhomom
    allocate( hab_rv(3,3,n0*nsp*nbas))
    allocate( sab_rv(3,3,n0,nsp,nbas))
    allocate( phzdphz(nppn,n0,nsp,nbas))
    allocate( fes1_rvx(3*nbas))
    allocate( spotx(n1,n2,n3,nsp),source=(0d0,0d0)) !smooth potential without XC
    allocate( ohsozzx(3,nbas), ohsopmx(3,nbas)) !dummy
    allocate( osigx(3,nbas), otaux(3,nbas), oppix(3,nbas))
    call dfaugm(osigx, otaux, oppix, ohsozzx,ohsopmx) !for sig,tau,ppi without XC(LDA)
    call mkpot(1, osmrho,orhoat, spotx,osigx,otaux,oppix,fes1_rvx,ohsozzx,ohsopmx, novxc_) !obtain osigx,otaux,oppix,smpotx  (without XC)
    !                                                                         When novxc_ exists, we exclud XC(LDA) part. We only need spotx and oppix
    deallocate(phzdphz,vesrmt,qmom,hab_rv,sab_rv,fes1_rvx)
  end subroutine m_mkpot_novxc
  subroutine m_mkpot_init()
    use m_supot,only: n1,n2,n3
    use m_density,only: osmrho, orhoat !main input to represent electron density
    use m_lmfinit,only: lso,nbas, nlibu,lmaxu,lldau,nsp,lat_alat,lxcf,lpzex
    use m_struc_def,only: s_rv1,s_sblock
    integer:: i,is,ib,kmax,lmxa,lmxh,nelt2,nlma,nlmh,iprint
    call tcn('m_mkpot_init')
    if(iprint()>=10) write(stdo,"(a)")'m_mkpot_init: Making one-particle potential ...'
    allocate( vesrmt(nbas))
    allocate( osmpot(n1,n2,n3,nsp)) 
    allocate( qmom(nlmxlx,nbas))
    allocate( hab_rv(3,3,n0*nsp*nbas))
    allocate( sab_rv(3,3,n0,nsp,nbas))
    allocate( phzdphz(nppn,n0,nsp,nbas))
    allocate( fes1_rv(3*nbas))
    allocate( osig(3,nbas), otau(3,nbas), oppi(3,nbas))
    allocate( ohsozz(3,nbas), ohsopm(3,nbas))
    call dfaugm(osig,otau,oppi,ohsozz,ohsopm) !allocate pointer array for sig,tau,ppi integrals
    call mkpot(1,osmrho,orhoat,osmpot,osig , otau , oppi, fes1_rv, ohsozz,ohsopm)
    call tcx('m_mkpot_init')
  end subroutine m_mkpot_init
  subroutine m_mkpot_energyterms(smrho_out,orhoat_out) 
    use m_MPItk,only: master_mpi
    use m_struc_def
    type(s_rv1):: orhoat_out(:,:)
    complex(8) :: smrho_out(:,:,:,:)
    call tcn('m_mkpot_energyterms')
    if(master_mpi) write(stdo,"('m_mkpot_energyterms')")
    if(allocated(fes2_rv)) deallocate(fes2_rv)
    allocate(fes2_rv(3*nbas))
    call mkpot(0, smrho_out,orhoat_out, osmpot,osig,otau,oppi,fes2_rv,ohsozz,ohsopm) !job=0 is for no augmentation term
    call tcx('m_mkpot_energyterms')
  end subroutine m_mkpot_energyterms
  subroutine mkpot(job,smrho,orhoat, smpot,osig,otau,oppi,fes,ohsozz,ohsopm, novxc_)!- Make the potential from the density (smrho, orhoat) !dipole_) 
    ! job=0 => not make core and augmentation matrices
    ! job=1 => make core and augmentation matrices     ! xxxxxx problematic option dipole_ removed. (for <i|{\bf r}|j> matrix for novxc)
    use m_lmfinit,only:lso,nbas,ispec,nlibu,lmaxu,lldau,nsp,alat=>lat_alat,lxcf,lpzex,nlmxlx
    use m_lattic,only: plat=>lat_plat, vol=>lat_vol,rv_a_opos
    use m_supot,only: n1,n2,n3
    use m_MPItk,only: master_mpi
    use m_struc_def
    use m_bstrux,only: m_bstrux_init
    use m_elocp,only: elocp
    use m_locpot,only: Locpot
    use m_smvxcm,only: smvxcm
    use m_smves,only: smves
    use m_ftox
    use m_rhomom,only: rhomom
    ! for job=0
    !o         utot   = total electrostatic energy
    !o         valves = valence rho * estat potential
    !o         cpnves = core+nuc * estat potential
    !o         rhoexc = rho * exc
    !o         rhovxc = rho * vxc
    !o         xcore  = rhoc * total potential
    !o         valvef = smrhov * vsm + sum_ib valvef_ib
    !o           valvef_ib = rhov * vtrue - smrhov * vsm)_ib
    !o         amom   = system magnetic moment

    !! documents below are under construction.
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   lfrce :nonzero =>  contribution to forces
    !i   lcplxp=1 only ::0 if ppi is real; 1 if ppi is complex
    !i   n1,n2,n3 dimensions of smrho for smooth crystal density
    !i   smrho :smooth crystal density, on a uniform mesh
    !i   orhoat:local atomic densities (true and smooth parts)
    !i   qbg   :homogeneous background charge
    
    !i ... The following are LDA+U inputs
    !i   vorb  :orbital dependent potential
    !i   nlibu :number of U blocks  used to dimension vorb
    !i   lmaxu :max l for U blocks  used to dimension vorb
    !i   lldau :lldau(ib)=0 => no U on this site otherwise
    !i         :U on site ib with dmat beginning at dmats(*,lldau(ib))
    
    !o Outputs:
    !o   smpot :smooth potential on a uniform mesh:
    !o         :Ves~ + vxc = Ves(n0 + compensating gaussians) + vxc
    !o   qmom  :multipole moments of valence sphere densities
    !o   vconst:additional constant potential term
    !o   vesrmt:electrostatic potential at rmt
    !o   qval  :total valence charge, including semicore states
    !o   qsc   :total charge from semicore states (local orbitals)
    !o   osig,otau,oppis: augmentation matrices
    !o   ppn   :nmto-like potential parameters
    !o   hab,vab,sab: augmentation matrices for local parts of the ham.
    !o         :See Remarks in augmat.f for their generation and doc.
    !o   vval  :coffs to YL expansion of es potential at MT boundary
    !o   gpot0 :integrals of gaussians times electrostatic potential
    !o   fes   :contribution to the force from electrostatic + xc potential
    !o   sham->eterms various integrals for the total energy are stored:
    !l Local variables
    !l   rvmusm:int rhosm * vxc(rhosm+smcor1) where smcor1 is portion
    !l         :of smooth core density treated nonperturbatively.

    ! xxx   rvepsm:int rhosm * exc(rhosm+smcor1) where smcor1 is portion
    ! xxx         :of smooth core density treated nonperturbatively.
    ! !Cl   rvepsv:integral of valence density times exc(valence density) !removed now
    ! !Cl   rvvxcv:integral of valence density times vxc(valence density) !removed now
    !l   rhvsm :integral n0~ phi0~
    !l   sgp0  :compensating gaussians * sm-Ves = int (n0~-n0) phi0~
    !l   valfsm:rhvsm - sgp0 + rvmusm + fcvxc0 = n0 (phi0~ + Vxc(n0))
    !l         :NB sgp0 associated with the local parts, as is done in
    !l         :making the ppi matrix elements.
    !l   valftr:local contribution to density * veff
    !l         := sgp0 - fcvxca + sum_ib valvfa_ib
    !l   valvfa:generated by locpot, called valvef there.  Sum_ib of:
    !l         :vefv1 - vefv2
    !l         := int[rho1*(v1-2*Z/r) - (rho2*v2)] - (n0~-n0)*gpotb - focvxc
    !l   lso   :nonzero => spin orbit coupling
    
    !r Remarks
    !r *The total density is a sum of three terms,
    !r
    !r    n0(mesh) + sum_RL (n_RL(r) - n0_RL(r))
    !r
    !r  The first term is the smooth density on a mesh of points; the
    !r  second is the true density and is defined on a radial mesh for each
    !r  sphere; the last is the 1-center expansion of the smooth density on
    !r  the radial mesh.  (Note: because of l-truncation, n0_R(r) is not
    !r  identical to the one-center expansion of n0(mesh).  The sum of the
    !r  three terms converges rapidly with l because errors in n_R(r) are
    !r  mostly canceled by errors in n0_R(r).)
    !r
    !r *Computation of the electrostatic energy:
    !r  We add and subtract a set of compensating gaussian orbitals
    !r
    !r    n0 + sum_RL Q_RL g_RL + sum_RL (n_RL(r) - n0_RL(r) - Q_RL g_RL)
    !r
    !r  which render the integral of the local part (the last 3 terms)
    !r  zero in each RL channel.  The g_RL must be localized enough that
    !r  their spillout beyond the MT radius is negligible.
    !r
    !r  We define
    !r
    !r    n0~ = n0 + compensating gaussians
    !r
    !r  In the interstitial, the electrostatic potential of n0~ is the true
    !r  estat potential.  The potential of n0 is called phi0 and the
    !r  potential of n0~ is called phi0~.  The total electrostatic energy
    !r  is computed as
    !r    the electrostatic energy of  n0~  +
    !r    the electrostatic energy of (neutral) local parts
    !r
    implicit none
    type(s_rv1) :: orhoat(3,nbas)
    type(s_cv5) :: oppi(*)
    type(s_sblock) :: ohsozz(*),ohsopm(*)
    type(s_rv4) :: otau(*), osig(*)
    complex(8):: smrho(n1,n2,n3,nsp),smpot(n1,n2,n3,nsp)
    logical:: cmdopt0,novxc,secondcall=.false.     !      integer,optional:: dipole_
    logical,optional:: novxc_
    integer:: job,i1,i2,i3,i,iprint,isw,isum, ifi,isp,j,k,ix,ismpot(4) !dipole,
    real(8):: hpot0_rv(nbas), dq,cpnvsa,qsmc,smq,smag,sum2,rhoex,rhoec,rhvsm,sqloc,sqlocc,saloc,uat,usm,valfsm, &
         valvfa,vvesat,vsum,zsum,rvvxcv(nsp),rvvxc(nsp),rvmusm(nsp),rmusm(nsp), rvepsm(nsp),vxcavg(nsp),repat(nsp),&
         repatx(nsp),repatc(nsp),rmuat(nsp),repsm(nsp),repsmx(nsp),repsmc(nsp),rhobg,gpot0(nlmxlx,nbas),vab_rv(3,3,n0*nsp*nbas),&
         vval(nlmxlx,nbas),fes(3,nbas),rhvsm0
    real(8),parameter:: minimumrho=1d-14,pi=4d0*datan(1d0),tpi=2d0*pi
    real(8),external:: rydberg
    character(80) :: outs
    character strn*120
    call tcn('mkpot')
    Printsmoothbackgroundcharge: if (qbg /= 0) then !
      rhobg = (3d0/4d0/pi*vol)**(1d0/3d0)
      if(master_mpi)write(stdo,ftox)' Energy for background charge q=',ftod(qbg),'radius r=',rhobg,'E=9/5*q*q/r=',1.8d0*qbg*qbg/rhobg
    endif Printsmoothbackgroundcharge
    SmoothPart: block
      integer:: ife
      call rhomom(orhoat, qmom,vsum) ! Multipole moments qmom is calculated
      call smves(qmom,gpot0,vval,hpot0_rv,smrho,smpot,vconst,smq,qsmc,fes,rhvsm0,rhvsm,zsum,vesrmt,qbg)!0th comp. of Estatic potential Ves and Ees
      if(master_mpi) then !.and.cmdopt0('--espot')) then !writeout electrostatic potential 
         open(newunit=ife,file='estaticpot.dat')
         !write(ife) n1,n2,n3,alat,plat
         ismpot = shape(smpot)
         iy=1
         do ix=1, ismpot(3),4
         do iz=1, ismpot(3)
            write(ife,'(i5,2e16.8," !eV ")') iz,ix,iy,(smpot(1,1,ix,1)+vconst)*rydberg()
            write(ife,'(i5,2e16.8," !eV ")') iz,ix,iy,(smpot(1,1,ix,1)+vconst)*rydberg()
            write(ife,'(i5,2e16.8," !eV ")') iz,ix,iy,(smpot(1,1,ix,1)+vconst)*rydberg()
            write(ife,'(i5,2e16.8," !eV ")') iz,ix,iy,(smpot(1,1,ix,1)+vconst)*rydberg()
         enddo   
            write(ife,*)
         enddo   
         close(ife)
      endif   
      smag = merge(2d0*dreal(sum(smrho(:,:,:,1)))*vol/(n1*n2*n3) - smq,0d0,nsp==2) !mag mom  ! 2*nup-(nup-ndn) = nup-ndn
      novxc= present(novxc_)
      repsm=0d0;  repsmx=0d0;  repsmc=0d0;   rmusm=0d0;  rvmusm=0d0
      ADDsmoothExchangeCorrelationPotential: if(.not.novxc) then
        block
          complex(8):: smvxc_zv(n1*n2*n3*nsp),smvx_zv(n1*n2*n3*nsp), smvc_zv(n1*n2*n3*nsp),smexc_zv(n1*n2*n3)
          real(8):: fxc_rv(3,nbas)
          smvxc_zv=0d0; smvx_zv=0d0; smvc_zv=0d0; smexc_zv=0d0; fxc_rv=0d0 !We use n0+n^sH_a to obtain smpot.
          call smvxcm(lfrce, smrho,smpot,smvxc_zv,smvx_zv,smvc_zv, smexc_zv,repsm,repsmx,repsmc,rmusm,rvmusm,rvepsm, fxc_rv )
          if( lfrce /= 0 ) fes = fes+fxc_rv
        endblock
      endif ADDsmoothExchangeCorrelationPotential
    endblock SmoothPart
    StructureConstantWhenExtendedLO: block
      call elocp() ! set ehl and rsml for extendet local orbitals
      if(sum(lpzex)/=0) call m_bstrux_init()!computes structure constant (C_akL Eq.(38) in /JPSJ.84.034702) when we have extended local orbital.
    endblock StructureConstantWhenExtendedLO
    MTpartsIntegrals: block
      call locpot(job,novxc,orhoat,qmom,vval,gpot0, & !Make local potential at atomic sites and augmentation matrices 
           osig,otau,oppi,ohsozz,ohsopm, phzdphz,hab_rv,vab_rv,sab_rv,  &
           vvesat,repat,repatx,repatc,rmuat, valvfa,xcore, sqloc,sqlocc,saloc,qval,qsc )
      valfsm = rhvsm0 + sum(rvmusm) - vconst*qbg ! 0th comp. of rho_val*Veff= rho0*Ves +rho0*Vxc -vconst*qbg 
      valvef = valfsm + valvfa                   ! Veff*n_val= veff0*rho0_val + veff1*rho1_val-veff2*rho2_val
      usm    = 0.5d0*rhvsm   ! 0th comp. of Eq.(27). rhvsm= \int 0thEes*(n0+n^c_sH +gaussians)
      uat    = 0.5d0*vvesat  ! vvesat= \int 1stEes*(n1+n^c)+\int (1stEes-zcontribution)*z  - \int 2ndEes*(n2+n_sH+gaussians)
      utot = usm + uat !Ees total electro static energy. Eq.(27)
      rhoexc = sum(repsm) + sum(repat) ! Exc=\int rho*exc 
      rhoex  = sum(repsmx)+ sum(repatx)! Ex 
      rhoec  = sum(repsmc)+ sum(repatc)! Ec
      rhovxc = sum(rmusm) + sum(rmuat) ! \int rho*Vxc
      dq = smq+sqloc + qsmc+sqlocc + qbg -zsum !smooth part + local part + smoothcore + core local + qbackground -Z
      amom = smag+saloc !magnetic moment
      if(dabs(dq)>1d-3) write(stdo,"(a,f13.6)")' (warning) system not neutral, dq=',dq
    endblock MTpartsIntegrals
    PrintEnergies: if(iprint() >= 30) then
       write(stdo,"('mkpot:',/'   Energy terms(Ry):',7x,'smooth',11x,'local',11x,'total')")
       write(stdo,680) &
            'rhoval*veff ',valfsm,valvfa,valvef, & !\int rho Veff
            'Eestatic    ',usm,uat,utot, & ! electrostatic energy
            'rho*exc     ',sum(repsm),sum(repat),rhoexc, &
            'rho*vxc     ',sum(rmusm),sum(rmuat),rhovxc, &
            'valence chg ',smq,sqloc,smq+sqloc !valence electron density, smooth part + local part
       if(nsp == 2) write (stdo,680) 'valence mag ',smag,saloc,amom
       write(stdo,680) 'core charge ',qsmc,sqlocc,qsmc+sqlocc
       write(stdo,"('   Charges:  valence',f12.5,'   cores',f12.5,'   nucleii',f12.5/'   hom background',f12.5, &
         '   deviation from neutrality: ',f12.5)") smq+sqloc,qsmc+sqlocc,-zsum,qbg,dq
680    format(3x,a,4x,3f17.6)
    endif PrintEnergies
    Printlogfile: if(iprint()>0) then
       write (stdl,"('fp qvl',f11.6,'  sm',f11.6,'  loc',f11.6,'  qbg',f11.6,' dQ',f11.6)") smq+sqloc,smq,sqloc,qbg,dq
       if (nsp == 2) write (stdl,"('fp mag ',f11.5,'  sm ',f11.5,'  loc ',f11.5)") smag+saloc,smag,saloc
       write (stdl,"('fp pot  rvxc',f18.7,'  rexc ',f18.7,'  rves ',f16.7)") rhovxc,rhoexc,utot
       write (stdl,"('fp pot  rex ',f18.7,'  rec ',f19.7)") rhoex,rhoec
    endif Printlogfile
    call tcx('mkpot')
  end subroutine mkpot
  subroutine dfaugm(osig, otau, oppi, ohsozz,ohsopm )!Allocate pointer arrays
    use m_struc_def,only:s_rv1,s_cv1,s_sblock
    use m_lmfinit,only: lso,nkaphh,nsp,nbas,ispec
    !o  osig,otau,oppi,ohsozz  :memory allocated (ohsoz is for SOC)
    !r Remarks
    !r   sig and tau are l diagonal, ppi is full matrix
    !r   Thus integral (P~_kL P~_k'L' - P_kL P_k'L') is diagonal in LL', sig(nf1,nf2,0..lmax) with lmax the l-cutoff
    !r   For sig(Pkl,Pkl), nf1=nf2==1+kmax; lmax=lmxa
    !r   For sig(Hsm,Pkl), nf1=nkaph and nf2=1+kmax; lmax=lmxh
    !r   For sig(Hsm,Hsm), nf1=nf2=nkaph; lmax = lmxh
    !
    !ohsozz(2,ib): Hsm*Pkl for Lz(diag)
    !ohsopm(2,ib): Hsm*Pkl soffd(:,1) = <H|L-(isp=1,isp=2)|P>,
    !               soffd(:,2)= <H|L+(isp=2,isp=1)|P>
    !       dagger(soffd(:,2))= <P|L-(isp=1,isp=2)|H>
    implicit none
    type(s_cv5) :: oppi(3,nbas) !poiner array for pi integral
    type(s_sblock):: ohsozz(3,nbas),ohsopm(3,nbas)
    type(s_rv4) :: otau(3,nbas)
    type(s_rv4)::  osig(3,nbas)
    integer :: ib,is,kmax,lmxa,lmxh,nelt1,nlma,nlmh,nelt,nkaph
    logical:: cmdopt0
    do  ib = 1, nbas
       is = ispec(ib)
       nkaph=nkaphh(is)
       lmxa=lmxa_i(is) !max l of augmenation
       lmxh=lmxb_i(is) !max l of head 
       kmax=kmxt_i(is) !0:kmax for radial index, ! nkaph:number of orbital types for a given L quantum no. in basis
       nlma = (lmxa+1)**2
       nlmh = (lmxh+1)**2
       if (lmxa == -1) cycle
       if (lmxh > lmxa) call rx('dfaugm: lmxh > lmxa unexpected')
       allocate(osig(1,ib)%v(0:kmax,0:kmax,0:lmxa,nsp) )    ! Pkl*Pkl
       allocate(otau(1,ib)%v(0:kmax,0:kmax,0:lmxa,nsp) )
       allocate(oppi(1,ib)%cv(0:kmax,0:kmax,nlma,nlma,nsp))
       allocate(osig(2,ib)%v(nkaph,0:kmax,0:lmxh,nsp))      ! Hsm*Pkl  
       allocate(otau(2,ib)%v(nkaph,0:kmax,0:lmxh,nsp))
       allocate(oppi(2,ib)%cv(nkaph,0:kmax,nlmh,nlma,nsp)) 
       allocate(osig(3,ib)%v(nkaph,nkaph,0:lmxh,nsp))       ! Hsm*Hsm
       allocate(otau(3,ib)%v(nkaph,nkaph,0:lmxh,nsp)) 
       allocate(oppi(3,ib)%cv(nkaph,nkaph,nlmh,nlmh,nsp))
       if(lso/=0 .OR. cmdopt0('--socmatrix')) then !spin-orbit copling matrix elements
          nelt = (kmax+1)*(kmax+1)*nlma*nlma !ohsopm (L- and L+) is irrelevant for lso=2
          allocate(ohsozz(1,ib)%sdiag(0:kmax,0:kmax,nlma,nlma,nsp),ohsopm(1,ib)%soffd(0:kmax,0:kmax,nlma,nlma,nsp)) ! Pkl*Pkl zz and pm component
          nelt = nkaph*(kmax+1)*nlmh*nlma 
          allocate(ohsozz(2,ib)%sdiag(nkaph,0:kmax,nlmh,nlma,nsp), ohsopm(2,ib)%soffd(nkaph,0:kmax,nlmh,nlma,nsp))! Hsm*Pkl
          nelt = nkaph*nkaph*nlmh*nlmh 
          allocate(ohsozz(3,ib)%sdiag(nkaph,nkaph,nlmh,nlmh,nsp), ohsopm(3,ib)%soffd(nkaph,nkaph,nlmh,nlmh,nsp))! Hsm*Hsm
       endif
    enddo
  end subroutine dfaugm
  subroutine m_mkpot_deallocate()
    if(allocated(vesrmt)) deallocate(vesrmt,fes1_rv,phzdphz,sab_rv,hab_rv,qmom,oppi,otau,osig,osmpot,ohsozz,ohsopm)
  end subroutine m_mkpot_deallocate
end module m_mkpot
