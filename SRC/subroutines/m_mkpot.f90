module m_mkpot ! Potential terms. See http://dx.doi.org/10.7566/JPSJ.84.034702
  use m_lmfinit,only: nbas,stdo,qbg=>zbak,ham_frzwf,lmaxu,nsp,nlibu,n0,nppn &
       ,lfrce,stdl, nchan=>pot_nlma, nvl=>pot_nlml,nkaph
  use m_struc_def,only: s_rv1,s_cv1,s_sblock,s_rv4,s_cv5
  public:: m_mkpot_init, m_mkpot_energyterms, m_mkpot_novxc, m_mkpot_deallocate !,m_Mkpot_novxc_dipole
  ! Potential terms, call m_mkpot_init. Generated at mkpot-locpot-augmat-gaugm
  type(s_sblock),allocatable,protected,public  :: ohsozz(:,:),ohsopm(:,:) !SOC matrix
  type(s_cv5),allocatable,protected,public  :: oppi(:,:) !pi-integral Eq.(C.6) 
  type(s_rv4),allocatable,protected,public  :: otau(:,:) !tau            (C.5) 
  type(s_rv4),allocatable,protected,public  :: osig(:,:) !sigma          (C.4) 
  complex(8),allocatable,protected ,public  :: osmpot(:,:,:,:)!0th component of Eq.(34)
  real(8),allocatable,protected,public:: fes1_rv(:), fes2_rv(:) !force terms
  real(8),allocatable,protected,public:: hab_rv(:,:,:), sab_rv(:,:,:), qmom(:),vesrmt(:)
  real(8),protected,public:: qval,vconst,qsc
  real(8),allocatable,protected,public:: phzdphz(:,:,:,:) !val and slo at Rmt for local orbitals.
  ! Energy terms by call m_mkpot_energyterms
  real(8),protected,public:: utot,rhoexc,xcore,valvef,amom, valves,cpnves,rhovxc
  ! NoVxc terms  by call m_mkpot_novxc
  type(s_cv5),allocatable,protected,public  :: oppix(:,:) !pi-integral without xc term
  complex(8),allocatable,protected ,public  :: spotx(:,:,:,:)!0th component of Eq.(34) without xc term
  private
  !! nov2021 dipole contribution added  (not working...)! oppixd,spotxd: dipole part to oppix,spotx
  !      type(s_cv1),allocatable,protected,public  :: oppixd(:,:,:)
  !      complex(8),allocatable,protected ,public  :: spotxd(:,:,:,:,:)
contains
  subroutine m_mkpot_novxc() ! outputs are oppix and spotx (for no vxc terms).
    use m_supot,only: k1,k2,k3
    use m_density,only: osmrho, orhoat !main input density
    logical:: novxc_
    real(8),allocatable :: fes1_rvx(:)! gpot0(:),vval(:),vab_rv(:,:,:) !dummy
    type(s_sblock),allocatable:: ohsozzx(:,:),ohsopmx(:,:) !dummy for SOC
    type(s_rv4),allocatable   :: otaux(:,:) !dummy
    type(s_rv4),allocatable   :: osigx(:,:) !dummy
    write(stdo,"(a)")' m_mkpot_novxc: Making one-particle potential without XC part ...'
    allocate( vesrmt(nbas))
    allocate( qmom(nvl)) !rhomom
    allocate( hab_rv(3,3,n0*nsp*nbas))
    allocate( sab_rv(3,3,n0*nsp*nbas))
    allocate( phzdphz(nppn,n0,nsp,nbas))
    allocate( fes1_rvx(3*nbas))
    allocate( spotx(k1,k2,k3,nsp)) !smooth potential without XC
    spotx=0d0
    allocate( ohsozzx(3,nbas), ohsopmx(3,nbas)) !dummy
    allocate( osigx(3,nbas), otaux(3,nbas), oppix(3,nbas))
    call dfaugm(osigx, otaux, oppix, ohsozzx,ohsopmx) !for sig,tau,ppi without XC(LDA)
    call mkpot(1, osmrho,orhoat, spotx,osigx,otaux,oppix,fes1_rvx,ohsozzx,ohsopmx, novxc_) !obtain osigx,otaux,oppix,smpotx  (without XC)
    !                                                                         When novxc_ exists, we exclud XC(LDA) part. We only need spotx and oppix
    deallocate(phzdphz,vesrmt,qmom,hab_rv,sab_rv,fes1_rvx)
  end subroutine m_mkpot_novxc
  subroutine m_mkpot_init()
    use m_supot,only: k1,k2,k3
    use m_density,only: osmrho, orhoat !main input density
    use m_lmfinit,only: lso,nbas,sspec=>v_sspec, nlibu,lmaxu,lldau,nsp,lat_alat,lxcf,lpzex
    use m_struc_def,only: s_rv1,s_sblock
    integer:: i,is,ib,kmax,lmxa,lmxh,nelt2,nlma,nlmh,iprint
    call tcn('m_mkpot_init')
    if(iprint()>=10) write(stdo,"(a)")' m_mkpot_init: Making one-particle potential ...'
    allocate( vesrmt(nbas))
    allocate( osmpot(k1,k2,k3,nsp)) 
    allocate( qmom(nvl))
    allocate( hab_rv(3,3,n0*nsp*nbas))
    allocate( sab_rv(3,3,n0*nsp*nbas))
    allocate( phzdphz(nppn,n0,nsp,nbas))
    allocate( fes1_rv(3*nbas))
    allocate( osig(3,nbas), otau(3,nbas), oppi(3,nbas))
    allocate( ohsozz(3,nbas), ohsopm(3,nbas))
    call dfaugm(osig,otau,oppi,ohsozz,ohsopm) !allocation for sig,tau,ppi integrals
    call mkpot(1,osmrho,orhoat,osmpot,osig , otau , oppi, fes1_rv, ohsozz,ohsopm)
    call tcx('m_mkpot_init')
  end subroutine m_mkpot_init
  subroutine m_mkpot_energyterms(smrho_out,orhoat_out) 
    use m_MPItk,only: master_mpi
    use m_struc_def
    type(s_rv1):: orhoat_out(:,:)
    complex(8) :: smrho_out(:)
    call tcn('m_mkpot_energyterms')
    if(master_mpi) write(stdo,*)
    if(master_mpi) write(stdo,"(' m_mkpot_energyterms')")
    if(allocated(fes2_rv)) deallocate(fes2_rv)
    allocate(fes2_rv(3*nbas))
    call mkpot(0, smrho_out,orhoat_out, osmpot,osig,otau,oppi,fes2_rv,ohsozz,ohsopm) !job=0 is for no augmentation term
    call tcx('m_mkpot_energyterms')
  end subroutine m_mkpot_energyterms
  !- Make the potential from the density (smrho, orhoat)
  subroutine mkpot(job,smrho,orhoat, smpot,osig,otau,oppi,fes,ohsozz,ohsopm, novxc_) !dipole_)
    ! job=0 => not make core and augmentation matrices
    ! job=1 => make core and augmentation matrices
    ! xxxxxx problematic option dipole_ removed. (for <i|{\bf r}|j> matrix for novxc)
    use m_lmfinit,only:lso,nbas,ispec,sspec=>v_sspec,nlibu,lmaxu,lldau,nsp,lat_alat,lxcf,lpzex
    use m_lattic,only: lat_plat,lat_qlat, lat_vol,rv_a_opos
    use m_supot,only: k1,k2,k3,lat_nabc
    use m_MPItk,only: master_mpi
    use m_struc_def
    use m_bstrux,only: m_bstrux_init
    use m_elocp,only: elocp
    use m_locpot,only: Locpot
    use m_smvxcm,only: smvxcm
    use m_smves,only: smves
    use m_ftox
    ! other other output in module variables
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
    !i   k1,k2,k3 dimensions of smrho for smooth crystal density
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
    integer :: job,i1,i2,i3
    type(s_rv1) :: orhoat(*)
    type(s_cv5) :: oppi(*)
    type(s_sblock) :: ohsozz(*),ohsopm(*)
    type(s_rv4) :: otau(*)
    type(s_rv4) :: osig(*)
    real(8)::  fes(3,nbas)
    complex(8):: smrho(k1,k2,k3,nsp),smpot(k1,k2,k3,nsp)
    logical:: cmdopt0,novxc
    logical,optional:: novxc_
    character(80) :: outs
    integer :: i,ipr,iprint,n1,n2,n3,ngabc(3),lxcfun,isw,isum
    real(8):: hpot0_rv(nbas), dq,cpnvsa, & 
         qsmc,smq,smag,sum2,rhoex,rhoec,rhvsm,sgp0, &
         sqloc,sqlocc,saloc,uat,usm,valfsm,valftr, &
         valvfa,vvesat, & !rvepva,rvexva,rvecva,rvvxva,rvepsa,rvvxca,
         vol,vsum,zsum,zvnsm, & !rvepsv(nsp),rvexv(nsp),rvecv(nsp), &
         rvvxcv(nsp),rvvxc(nsp),  & !,rveps(nsp)
         rvmusm(nsp),rmusm(nsp), rvepsm(nsp), &
         vxcavg(nsp),repat(nsp),repatx(nsp),repatc(nsp), & !,fcvxca(nsp),fcvxc0(nsp)
         rmuat(nsp),repsm(nsp),repsmx(nsp),repsmc(nsp),rhobg
    equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
    complex(8),allocatable:: smpotbk(:,:,:),smpotbkx(:,:,:)
    real(8),parameter:: minimumrho=1d-14
    real(8):: plat(3,3),alat
    integer:: ifi,isp
    character strn*120
    logical:: secondcall=.false.     !      integer,optional:: dipole_
    integer:: j,k !dipole,
    real(8),parameter::  pi=4d0*datan(1d0),tpi=2d0*pi
    real(8):: gpot0(nvl), vab_rv(3,3,n0*nsp*nbas),vval(nchan)
    call tcn('mkpot')
    if(cmdopt0('--density') .AND. master_mpi .AND. secondcall) then ! new density mode
       plat =lat_plat
       alat =lat_alat
       open(newunit=ifi,file='smrho.xsf')
       do isp = 1, nsp
          write(ifi,'("CRYSTAL")')
          write(ifi,'("PRIMVEC")')
          write(ifi,'(3f10.5)') ((plat(i1,i2)*alat*0.529177208,i1=1,3) &
               ,i2=1,3)
          write(ifi,'("PRIMCOORD")')
          write(ifi,'(2i5)') nbas,1
          do i = 1, nbas
             write(ifi,'(i4,2x,3f10.5)')sspec(ispec(i))%z,(rv_a_opos(i2,i)*alat*0.529177208,i2=1,3)
          enddo
          write(ifi,'("BEGIN_BLOCK_DATAGRID_3D")')
          write(ifi,'("charge_density_spin_",i1)') isp
          write(ifi,'("BEGIN_DATAGRID_3D_isp_",i1)') isp
          write(ifi,'(3i4)') k1,k2,k3
          write(ifi,'(3f10.5)') 0.,0.,0.
          write(ifi,'(3f10.5)') ((plat(i1,i2)*alat*0.529177208,i1=1,3),i2=1,3)
          write(ifi,'(8e14.6)') (((dble(smrho(i1,i2,i3,isp)),i1=1,k1),i2=1,k2),i3=1,k3)
          write(ifi,'("END_DATAGRID_3D_isp_",i1)') isp
          write(ifi,'("END_BLOCK_DATAGRID_3D")')
       enddo
       close(ifi)
    else
       secondcall=.true.
    endif
    ipr = iprint()
    lxcfun = lxcf
    ngabc=lat_nabc
    vol=lat_vol
    if (qbg /= 0) then !Printout for smooth background charge 
       rhobg = (3d0/4d0/pi*vol)**(1d0/3d0)
       if(master_mpi) write(stdo,ftox)' Energy for background charge', &
            ' q=',ftod(qbg),'radius r=',rhobg,'E=9/5*q*q/r=',1.8d0*qbg*qbg/rhobg
    endif
    call rhomom(orhoat, qmom,vsum) !multipole moments
    call smves(qmom,gpot0,vval,hpot0_rv,sgp0,smrho,smpot,vconst,smq,qsmc,fes,rhvsm,zvnsm,zsum,vesrmt,qbg )  ! Smooth electrostatic potential
    smag = 0
    if(nsp == 2) smag = 2d0*dreal(sum(smrho(:,:,:,1)))*vol/(n1*n2*n3) - smq !spin part.
    if( .NOT. present(novxc_)) then ! Add smooth exchange-correlation potential 
       novxc=.false.
       block
         complex(8):: smvxc_zv(k1*k2*k3*nsp),smvx_zv(k1*k2*k3*nsp), smvc_zv(k1*k2*k3*nsp),smexc_zv(k1*k2*k3)
         real(8):: fxc_rv(3,nbas)
         smvxc_zv=0d0; smvx_zv=0d0; smvc_zv=0d0; smexc_zv=0d0; fxc_rv=0d0
         call smvxcm(lfrce, smrho,smpot,smvxc_zv,smvx_zv,smvc_zv, smexc_zv,repsm,repsmx,repsmc,rmusm,rvmusm,rvepsm, fxc_rv )
         if( lfrce /= 0 ) fes = fes+fxc_rv
       endblock
    else
       novxc=.true.
       repsm=0d0;  repsmx=0d0;  repsmc=0d0;   rmusm=0d0;  rvmusm=0d0;
    endif
    !! Add dipole contribution (x,y,z) to smpot (we only need <i|x|j> and so on, but because of technical reason,
    !! Calculate <i|x|j> by subtraction. nov2021 (this is stupid implementation--- See MLWF paper.
    !$$$      dipole=0
    !$$$      if(present(dipole_)) then
    !$$$         if(dipole_/=0) dipole=dipole_
    !$$$         write(6,*)' mkpot dipole=', present(dipole_),dipole_,dipole,k1,k2,k3
    !$$$      endif
    !$$$      if(dipole/=0) then
    !$$$         do isp=1,nsp;  do i=1,k1;    do j=1,k2;  do k=1,k3
    !$$$            smpot(i,j,k,isp)=smpot(i,j,k,isp)+ lat_alat*sum(plat(dipole,:)*[(i-1)/dble(k1),(j-1)/dble(k2),(k-1)/dble(k3)])
    !$$$         enddo;   enddo;    enddo;    enddo
    !$$$      endif
    call elocp() ! set ehl and rsml for extendet local orbitals
    if(sum(lpzex)/=0) call m_bstrux_init()!computes structure constant (C_akL Eq.(38) in /JPSJ.84.034702) when we have extended local orbital.
    call locpot(job,novxc,orhoat,qmom,vval,gpot0, & !,idipole ) !Make local potential at atomic sites and augmentation matrices 
         osig,otau,oppi,ohsozz,ohsopm, phzdphz,hab_rv,vab_rv,sab_rv,  &
         vvesat,cpnvsa, repat,repatx,repatc,rmuat, valvfa,xcore, sqloc,sqlocc,saloc,qval,qsc )
    if(cmdopt0('--density') .AND. master_mpi .AND. secondcall) return
    ! Integral of valence density times estatic potential
    ! Associate term (n0~-n0) Ves(n0~) with local part because of the ppi matrix elements
    ! Also add fcvxc0(1) to smooth part because rvmusm+fcvxc0 is perturbative approximation for rvmusm when cores are not treated perturbatively.
    valves = rhvsm + vvesat ! ... Valence density times VEelectroStatic
    valfsm = rhvsm + sum(rvmusm) - sgp0 - vconst*qbg !rho*Ves +rho*Vxc - Qmom*Ves -vconst*qbg
    valftr = valvfa + sgp0    ! atomic rho*veff + Qmom*Ves
    valvef = valfsm + valftr
    cpnves = zvnsm + cpnvsa! ... Integral of core+nucleus times Ves(estatic potential)
    rhoexc = sum(repsm) + sum(repat) ! Exc=\int rho*exc 
    rhoex  = sum(repsmx)+ sum(repatx)! Ex 
    rhoec  = sum(repsmc)+ sum(repatc)! Ec
    rhovxc = sum(rmusm) + sum(rmuat) ! \int rho*Vxc
    usm = 0.5d0*(rhvsm+zvnsm)
    uat = 0.5d0*(vvesat+cpnvsa)
    utot = usm + uat !total electro static energy
    dq = smq+sqloc + qsmc+sqlocc + qbg -zsum !smooth part + local part + smoothcore + core local + qbackground -Z
    amom = smag+saloc !magnetic moment
    if(ipr >= 30) then
       write(stdo,"('  mkpot:',/'   Energy terms:',11x,'smooth',11x,'local',11x,'total')")
       write(stdo,680) 'rhoval*veff ',valfsm,valftr,valvef, & !\int rho Veff
            'rhoval*ves ',rhvsm,vvesat,valves, & !\int rho Ves
            'psnuc*ves  ',zvnsm,cpnvsa,cpnves, & !\int rho(Z+core) Ves
            'utot       ',usm,uat,utot, & !total electrostatic energy
            'rho*exc    ',sum(repsm),sum(repat),rhoexc, &
            'rho*vxc    ',sum(rmusm),sum(rmuat),rhovxc, &
            'valence chg',smq,sqloc,smq+sqloc !valence electron density, smooth part + local part
       if (nsp == 2) write (stdo,680) 'valence mag',smag,saloc,amom
       write (stdo,680) 'core charge',qsmc,sqlocc,qsmc+sqlocc
       write (stdo,670) smq+sqloc,qsmc+sqlocc,-zsum,qbg,dq
670    format('   Charges:  valence',f12.5,'   cores',f12.5,'   nucleii',f12.5/'   hom background',f12.5, &
         '   deviation from neutrality: ',f12.5)
680    format(3x,a,4x,3f17.6)
    endif
    if(ipr>0) then
       write (stdl,"('fp qvl',f11.6,'  sm',f11.6,'  loc',f11.6,'  qbg',f11.6,' dQ',f11.6)") smq+sqloc,smq,sqloc,qbg,dq
       if (nsp == 2) write (stdl,"('fp mag ',f11.5,'  sm ',f11.5,'  loc ',f11.5)") smag+saloc,smag,saloc
       write (stdl,"('fp pot  rvxc',f18.7,'  rexc ',f18.7,'  rves ',f16.7)") rhovxc,rhoexc,utot
       write (stdl,"('fp pot  rex ',f18.7,'  rec ',f19.7)") rhoex,rhoec
    endif
    if(dabs(dq)>1d-3) write(stdo,"(a,f13.6)")' (warning) system not neutral, dq=',dq
    call tcx('mkpot')
  end subroutine mkpot
  subroutine dfaugm(osig, otau, oppi, ohsozz,ohsopm )
    use m_struc_def,only:s_rv1,s_cv1,s_sblock
    use m_lmfinit,only: lso,nkaph,nsp,nbas,ispec,sspec=>v_sspec
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
    type(s_cv5) :: oppi(3,nbas)
    type(s_sblock):: ohsozz(3,nbas),ohsopm(3,nbas)
    type(s_rv4) :: otau(3,nbas)
    type(s_rv4)::  osig(3,nbas)
    integer :: ib,is,kmax,lmxa,lmxh,nelt1,nlma,nlmh,nelt
    logical:: cmdopt0
    do  ib = 1, nbas
       is = ispec(ib)
       lmxa=sspec(is)%lmxa !max l of augmenation
       lmxh=sspec(is)%lmxb !max l of head 
       kmax=sspec(is)%kmxt !0:kmax for radial index, ! nkaph:number of orbital types for a given L quantum no. in basis
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
          allocate(ohsozz(1,ib)%sdiag(nelt,nsp),ohsopm(1,ib)%soffd(nelt,nsp)) ! Pkl*Pkl zz and pm component
          nelt = nkaph*(kmax+1)*nlmh*nlma 
          allocate(ohsozz(2,ib)%sdiag(nelt,nsp), ohsopm(2,ib)%soffd(nelt,nsp))! Hsm*Pkl
          nelt = nkaph*nkaph*nlmh*nlmh 
          allocate(ohsozz(3,ib)%sdiag(nelt,nsp), ohsopm(3,ib)%soffd(nelt,nsp))! Hsm*Hsm
       endif
    enddo
  end subroutine dfaugm
  subroutine m_mkpot_deallocate()
    if (allocated(vesrmt)) then
       deallocate(vesrmt,fes1_rv,phzdphz,sab_rv,hab_rv,qmom,oppi,otau,osig,osmpot,ohsozz,ohsopm)
    endif
  end subroutine m_mkpot_deallocate
end module m_mkpot

  !$$$      subroutine m_mkpot_novxc_dipole()
  !$$$! output: oppixd, spotxd
  !$$$! Assuming osigx,otaux are allocated already by calling m_mkpot_novxc in advance.
  !$$$! nov2021:   dipole option is added whether we calculate <i|{\bf r}|j> matrix (by differece).
  !$$$      use m_supot,only: k1,k2,k3
  !$$$      use m_density,only: osmrho, orhoat !main input density
  !$$$      use m_lmfinit,only: nkaph,nsp,nbas,ssite=>v_ssite,sspec=>v_sspec
  !$$$      integer:: i,lfrzw,ib,is,kmax,lmxa,lmxh,nlma,nlmh,iidipole,ilfzw
  !$$$      logical:: novxc_
  !$$$      write(stdo,"(a)")' m_mkpot_novxc_dipole: Making one-particle potential without XC part ...'
  !$$$      allocate( vesrmt(nbas))
  !$$$      allocate(  qmom(nvl)) !rhomom
  !$$$      allocate( phzdphz(nppn,n0,nsp,nbas))
  !$$$      allocate(  hab_rv(nab*n0*nsp*nbas))
  !$$$      allocate(  vab_rv(nab*n0*nsp*nbas))
  !$$$      allocate(  sab_rv(nab*n0*nsp*nbas))
  !$$$      allocate(  gpot0(nvl))
  !$$$      allocate(  vval(nchan))
  !$$$      allocate(  fes1_rv(3*nbas))
  !$$$      allocate( ohsozzx(3,nbas), ohsopmx(3,nbas)) !dummy
  !$$$      lfrzw = 0
  !$$$      if(ham_frzwf) lfrzw = 1   !freeze all augmentation wave
  !$$$      ilfzw = 1 + 10*lfrzw !+ 100    ! Adding 100 means excluding XC(LDA) part. nolxc=T
  !$$$      allocate( spotxd(k1,k2,k3,nsp,3), oppixd(3,nbas,3)) !smooth potential without XC
  !$$$      spotxd=0d0
  !$$$      do  ib = 1, nbas
  !$$$         is = int(ssite(ib)%spec)
  !$$$         lmxa=sspec(is)%lmxa
  !$$$         lmxh=sspec(is)%lmxb
  !$$$         kmax=sspec(is)%kmxt
  !$$$         nlma = (lmxa+1)**2
  !$$$         nlmh = (lmxh+1)**2
  !$$$         do iidipole=1,3        !see dfaugm
  !$$$            allocate(oppixd(1,ib,iidipole)%cv((kmax+1)*(kmax+1)*nlma*nlma*nsp))
  !$$$            allocate(oppixd(3,ib,iidipole)%cv( nkaph*nkaph*nlmh*nlmh*nsp))
  !$$$            allocate(oppixd(2,ib,iidipole)%cv( nkaph*(kmax+1)*nlmh*nlma*nsp))
  !$$$         enddo
  !$$$      enddo
  !$$$      do iidipole=1,3
  !$$$         call mkpot(ilfzw,osmrho,orhoat,
  !$$$     o        spotxd(:,:,:,:,iidipole),osigx,otaux,oppixd(:,:,iidipole),fes1_rv,ohsozzx,ohsopmx,
  !$$$     &        novxc_,dipole_=iidipole) !!when novxc_ exists, we exclud XC(LDA) part.
  !$$$      enddo
  !$$$      deallocate(vesrmt,qmom,phzdphz,hab_rv,vab_rv,sab_rv,gpot0,vval,fes1_rv,ohsozzx,ohsopmx)
  !$$$      end subroutine
