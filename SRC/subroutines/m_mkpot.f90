module m_mkpot ! http://dx.doi.org/10.7566/JPSJ.84.034702
  use m_lmfinit,only: nbas,stdo,qbg=>zbak,ham_frzwf,lmaxu,nsp,nlibu,n0,nab,nppn &
       ,lfrce=>ctrl_lfrce,stdl, nchan=>pot_nlma, nvl=>pot_nlml,nkaph
  use m_struc_def,only: s_rv1,s_cv1,s_sblock
!  use m_MPItk,only: master_mpi,mlog
  !! output ------------------------------------------------------------------
  public::  m_Mkpot_init, m_Mkpot_deallocate, m_Mkpot_energyterms, m_Mkpot_novxc!,m_Mkpot_novxc_dipole
  type(s_sblock),allocatable,protected,public  :: ohsozz(:,:),ohsopm(:,:) !SOC matrix
  ! Generated at mkpot-locpot-augmat-gaugm
  type(s_cv1),allocatable,protected,public  :: sv_p_oppi(:,:) !pi-integral Eq.(C.6)
  type(s_rv1),allocatable,protected,public  :: sv_p_otau(:,:) !tau            (C.5)
  type(s_rv1),allocatable,protected,public  :: sv_p_osig(:,:) !sigma          (C.4)
  complex(8),allocatable,protected ,public  :: osmpot(:,:,:,:)!0th component of Eq.(34)

  type(s_cv1),allocatable,protected,public  :: sv_p_oppix(:,:) !pi-integral without xc term
  complex(8),allocatable,protected ,public  :: spotx(:,:,:,:)

  real(8),protected,public:: utot,rhoexc,xcore,valvef,amom,  valves,cpnves,rhovxc !energy terms
  real(8),allocatable,protected,public:: fes1_rv(:), fes2_rv(:) !force terms
  real(8),allocatable,protected,public:: hab_rv(:), sab_rv(:), ppnl_rv(:,:,:,:), qmom(:),vesrmt(:)
  real(8),protected,public:: qval,vconst,qsc

  !! nov2021 dipole contribution added
  !! sv_p_oppixd: add dipole part to sv_p_oppix
  !! spotxd:      add dipole part to spotx
  !      type(s_cv1),allocatable,protected,public  :: sv_p_oppixd(:,:,:)
  !      complex(8),allocatable,protected ,public  :: spotxd(:,:,:,:,:)

  private
  real(8),allocatable,protected,private :: gpot0(:),vval(:),vab_rv(:) !dummy
  type(s_sblock),allocatable,private  :: ohsozzx(:,:),ohsopmx(:,:) !dummy for SOC
  type(s_rv1),allocatable,protected,private  :: sv_p_otaux(:,:) !dummy
  type(s_rv1),allocatable,protected,private  :: sv_p_osigx(:,:) !dummy
contains
  subroutine m_mkpot_novxc()
    ! output: sv_p_oppix, spotx
    use m_supot,only: k1,k2,k3
    use m_density,only: osmrho, orhoat !main input density
    integer:: lfrzw,ilfzw
    logical:: novxc_
    write(stdo,"(a)")' m_mkpot_novxc: Making one-particle potential without XC part ...'
    allocate( vesrmt(nbas))
    allocate(  qmom(nvl)) !rhomom
    allocate( ppnl_rv(nppn,n0,nsp,nbas))
    allocate(  hab_rv(nab*n0*nsp*nbas))
    allocate(  vab_rv(nab*n0*nsp*nbas))
    allocate(  sab_rv(nab*n0*nsp*nbas))
    allocate(  gpot0(nvl))
    allocate(  vval(nchan))
    allocate(  fes1_rv(3*nbas))
    allocate( ohsozzx(3,nbas), ohsopmx(3,nbas)) !dummy
    allocate( spotx(k1,k2,k3,nsp)) !smooth potential without XC
    spotx=0d0
    allocate( sv_p_osigx(3,nbas), sv_p_otaux(3,nbas), sv_p_oppix(3,nbas))
    call dfaugm(sv_p_osigx, sv_p_otaux, sv_p_oppix, ohsozzx,ohsopmx)!for sig,tau,ppi without XC(LDA)
    !     We obtain sv_p_osigx, sv_p_otaux, sv_p_oppix, smpotx  without XC
    call mkpot(1, osmrho , orhoat , &
         spotx,sv_p_osigx,sv_p_otaux,sv_p_oppix,  fes1_rv,ohsozzx,ohsopmx, &
         novxc_) !when novxc_ exists, we exclud XC(LDA) part.
    deallocate(vesrmt,qmom,ppnl_rv,hab_rv,vab_rv,sab_rv,gpot0,vval,fes1_rv,ohsozzx,ohsopmx)
  end subroutine m_mkpot_novxc
  !!
  !$$$      subroutine m_mkpot_novxc_dipole()
  !$$$! output: sv_p_oppixd, spotxd
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
  !$$$      allocate( ppnl_rv(nppn,n0,nsp,nbas))
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
  !$$$      allocate( spotxd(k1,k2,k3,nsp,3), sv_p_oppixd(3,nbas,3)) !smooth potential without XC
  !$$$      spotxd=0d0
  !$$$      do  ib = 1, nbas
  !$$$         is = int(ssite(ib)%spec)
  !$$$         lmxa=sspec(is)%lmxa
  !$$$         lmxh=sspec(is)%lmxb
  !$$$         kmax=sspec(is)%kmxt
  !$$$         nlma = (lmxa+1)**2
  !$$$         nlmh = (lmxh+1)**2
  !$$$         do iidipole=1,3        !see dfaugm
  !$$$            allocate(sv_p_oppixd(1,ib,iidipole)%cv((kmax+1)*(kmax+1)*nlma*nlma*nsp))
  !$$$            allocate(sv_p_oppixd(3,ib,iidipole)%cv( nkaph*nkaph*nlmh*nlmh*nsp))
  !$$$            allocate(sv_p_oppixd(2,ib,iidipole)%cv( nkaph*(kmax+1)*nlmh*nlma*nsp))
  !$$$         enddo
  !$$$      enddo
  !$$$      do iidipole=1,3
  !$$$         call mkpot(ilfzw,osmrho,orhoat,
  !$$$     o        spotxd(:,:,:,:,iidipole),sv_p_osigx,sv_p_otaux,sv_p_oppixd(:,:,iidipole),fes1_rv,ohsozzx,ohsopmx,
  !$$$     &        novxc_,dipole_=iidipole) !!when novxc_ exists, we exclud XC(LDA) part.
  !$$$      enddo
  !$$$      deallocate(vesrmt,qmom,ppnl_rv,hab_rv,vab_rv,sab_rv,gpot0,vval,fes1_rv,ohsozzx,ohsopmx)
  !$$$      end subroutine

  ! ssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine m_mkpot_init(llmfgw)
    use m_supot,only: k1,k2,k3
    use m_density,only: osmrho, orhoat !main input density
    use m_lmfinit,only: lso,nbas, ssite=>v_ssite, sspec=>v_sspec, nlibu,lmaxu,lldau,nsp,lat_alat,lxcf,lpzex
    integer:: i,lfrzw,is,ib,kmax,lmxa,lmxh,nelt2,nlma,nlmh,iprint!,lfrce
    real(8) ::qbz
    logical :: llmfgw,fsmode,fullmesh,cmdopt0
    call tcn('m_mkpot_init')
    if(iprint()>=10) write(stdo,"(a)")' m_mkpot_init: Making one-particle potential ...'
    lfrzw = 0
    if(ham_frzwf) lfrzw = 1   !freeze all augmentation wave
    !! Make the potential and total energy terms for given density (smrho,rhoat,qbg)  ---
    !! mkpot calls locpot. and locpot calls augmat. augmat calculates sig,tau,ppi.
    i = 1 + 10*lfrzw
!    if(llmfgw) i = i + 10000 !GW driver mode
    !! Arrays used in the generation of the potential ---
    allocate( vesrmt(nbas))
    allocate( osmpot(k1,k2,k3,nsp)) !smooth potential without XC
    allocate( qmom(nvl))
    allocate( gpot0(nvl))
    allocate( vval(nchan))
    allocate(  hab_rv(nab*n0*nsp*nbas))
    allocate(  vab_rv(nab*n0*nsp*nbas))
    allocate(  sab_rv(nab*n0*nsp*nbas))
    allocate( ppnl_rv(nppn,n0,nsp,nbas))
    allocate( fes1_rv(3*nbas))
    allocate( sv_p_osig(3,nbas), sv_p_otau(3,nbas), sv_p_oppi(3,nbas))
    allocate( ohsozz(3,nbas), ohsopm(3,nbas))
    call dfaugm(sv_p_osig,sv_p_otau,sv_p_oppi,ohsozz,ohsopm) !allocation for sig,tau,ppi integrals
    call mkpot (1, osmrho , orhoat, &
         osmpot, sv_p_osig , sv_p_otau , sv_p_oppi, fes1_rv, ohsozz,ohsopm)
    call tcx('m_mkpot_init')
  end subroutine m_mkpot_init
  ! ssssssssssssssssssssssssssssssssssssssssssssssssssssss
  
  subroutine m_mkpot_energyterms(smrho_out,orhoat_out)
    use m_struc_def
    type(s_rv1):: orhoat_out(:,:)
    complex(8) :: smrho_out(:)
    call tcn('m_mkpot_energyterms')
    if(allocated(fes2_rv)) deallocate(fes2_rv)
    allocate( fes2_rv(3*nbas))
    call mkpot(0, smrho_out , orhoat_out , &
         osmpot, sv_p_osig , sv_p_otau , sv_p_oppi, fes2_rv, ohsozz,ohsopm)
    call tcx('m_mkpot_energyterms')
  end subroutine m_mkpot_energyterms
   
  subroutine mkpot(job, smrho, sv_p_orhoat, &
       smpot, sv_p_osig, sv_p_otau, sv_p_oppi, fes, ohsozz,ohsopm, novxc_) !dipole_)
    ! xxx problematic option dipole_ removed. (for <i|{\bf r}|j> matrix for novxc)
    use m_lmfinit,only:lso,nbas,ssite=>v_ssite,sspec=>v_sspec,nlibu,lmaxu,lldau,nsp,lat_alat,lxcf,lpzex
    use m_lattic,only: lat_plat,lat_qlat, lat_vol,rv_a_opos
    use m_supot,only: k1,k2,k3,lat_nabc
    use m_MPItk,only: master_mpi
    use m_struc_def
    use m_ldau,only: vorb !input. 'U-V_LDA(couter term)' of LDA+U
    use m_bstrux,only: m_bstrux_init
    use m_elocp,only: elocp
    use m_locpot,only: Locpot
    use m_smvxcm,only: smvxcm
    use m_ftox
    ! ote other output in module area
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
    !- Make the potential from the density.
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nbas  :size of basis
    !i   ssite :struct containing site-specific information
    !i     Passed to: rhomom smves smvxcm locpot
    !i   sspec :struct containing species-specific information
    !i     Elts read:
    !i     Passed to: rhomom smves smvxcm locpot
    !i   lfrce :nonzero =>  contribution to forces
    !i   lcplxp=1 only ::0 if ppi is real; 1 if ppi is complex
    !i   k1,k2,k3 dimensions of smrho for smooth crystal density
    !i   smrho :smooth crystal density, on a uniform mesh
    !i   orhoat:local atomic densities (true and smooth parts)
    !i   qbg   :homogeneous background charge
    !i   job   :1s digit
    !i         : 0 stops after energy terms
    !i         : 1 makes potpars also
    !i         :10s digit
    !i         :1 suppress updating potential used to make potpars
    !i         :100s digit
    !i         :1 exclude exchange-correlation potential
    ! xxx         :1000s digit
    ! xxx         :1 Make rveps and rvvxc
    !i         :10000s digit
    !i         :1 write sphere density for each site ib to file rhoMT.ib
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
    ! Cl   rvepsv:integral of valence density times exc(valence density) !removed now
    ! Cl   rvvxcv:integral of valence density times vxc(valence density) !removed now
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
    !r  The first term is computed in subroutine smves;
    !r  the second term is computed in subroutine locpot.
    !u Updates
    !u   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
    !u   27 Apr 05 Added LDA+U stuff (Lambrecht)
    !u   24 Dec 04 Changes for full L.S coupling
    !u    1 Sep 04 Adapted mkpot to handle complex ppi; fold so into ppi
    !u   15 Jul 04 (Chantis) Matrix elements Lz.Sz for spin-orbit coupling
    !u   14 Jan 02 rvexv and rvecv (T. Miyake)
    !u   17 Sep 01 Returns qsc.  Altered argument list.
    !u   24 Aug 01 Extended to local orbitals.  Altered argument list.
    !u             Local potentials now have correct boundary conditions
    !u   15 Aug 01 Generates rvepsv and rvvxcv.  Changed call to locpot.
    !u   20 Apr 01 Generates vesrmt
    !u   18 Apr 01 Added ability to exclude exchange-correlation potential
    !u   20 Feb 01 Added ppn to potential parameters generated
    !u   15 Jun 00 spin polarized
    !u   22 Apr 00 Adapted from nfp mk_potential.f
    ! ------------------------------------------------------------
    ! to do:
    ! 1. check that rhov*vxc approximately the same in locpot as smvxc
    ! 2. ditto for rhov*exc; need finish making rhov*exc in smvxc
    ! 3. ? See about understanding and changing foxexc both in locpot
    !    and smvxc.  The total energy seems to come out just right,
    !    but focexc isn't real correction to rhov*exc.  Does it matter?
    !    Maybe can't use lfoca=2 is all.
    ! 4. enable 1000s digit job to pass through
    ! 5. change locpot so that rvvxcv etc is make only if appropriate
    !    1000s job bit set.
    ! 6. Make <vxc_nn> and compare to integrals made here.
    !    If problems, Maybe even look at sm part only, or local part only.
    implicit none
    integer :: job,i1,i2,i3,ispec
    type(s_rv1) :: sv_p_orhoat(*)
    type(s_cv1) :: sv_p_oppi(*)
    type(s_sblock) :: ohsozz(*),ohsopm(*)
    type(s_rv1) :: sv_p_otau(*)
    type(s_rv1) :: sv_p_osig(*)
    real(8)::  fes(3,nbas)
    complex(8):: smrho(k1,k2,k3,nsp),smpot(k1,k2,k3,nsp)
    logical:: cmdopt0,novxc
    logical,optional:: novxc_
    character(80) :: outs
    integer :: i,lgunit,nglob,ipl,ipr,iprint,n1,n2,n3,ngabc(3),lxcfun,isw,isum
    real(8) ,allocatable :: fxc_rv(:,:)
    real(8) ,allocatable :: hpot0_rv(:)
    complex(8) ,allocatable :: smvxc_zv(:)
    complex(8) ,allocatable :: smvx_zv(:)
    complex(8) ,allocatable :: smvc_zv(:)
    complex(8) ,allocatable :: smexc_zv(:)
    double precision :: dq,focexc,cpnvsa, &
         focvxc,qsmc,smq,smag,sum2,rhoex,rhoec,rhvsm,sgp0, &
         sqloc,sqlocc,saloc,uat,usm,valfsm,valftr, &
         rvepva,rvexva,rvecva,rvvxva,rvepsa,rvvxca,valvfa,vvesat, &
         vol,vsum,zsum,zvnsm,rvepsv(2),rvexv(2),rvecv(2), &
         rvvxcv(2),fcexc0(2),rveps(2),rvvxc(2),fcex0(2),fcec0(2), &
         fcexca(2),fcexa(2),fceca(2),rvmusm(2),rvepsm(2),rmusm(2), &
         vxcavg(2),fcvxc0(2),fcvxca(2),repat(2),repatx(2),repatc(2), &
         rmuat(2),repsm(2),repsmx(2),repsmc(2),rhobg
    equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
    integer:: itest,wdummy,nnn
    complex(8),allocatable:: smpotbk(:,:,:),smpotbkx(:,:,:)
    real(8),parameter:: minimumrho=1d-14
    real(8)::sss,smmin(2),srshift
    real(8):: plat(3,3),alat
    integer:: ifi,isp
    character strn*120
    logical:: secondcall=.false.
    !      integer,optional:: dipole_
    integer:: j,k !dipole,
    real(8),parameter::  pi=4d0*datan(1d0),tpi=2d0*pi
    call tcn('mkpot')
    !! new density mode
    if(cmdopt0('--density') .AND. master_mpi .AND. secondcall) then
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
             ispec=ssite(i)%spec
             write(ifi,'(i4,2x,3f10.5)') sspec(ispec)%z,(rv_a_opos(i2,i)*alat*0.529177208,i2=1,3)
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
    ipl = ipr
    lxcfun = lxcf
    ngabc=lat_nabc
    vol=lat_vol
    ! --- Printout for smooth background charge ---
    if (qbg /= 0) then
       rhobg = (3d0/4d0/pi*vol)**(1d0/3d0)
       write(stdo,ftox)' Energy for background charge', &
            ' q=',ftod(qbg),'radius r=',rhobg,'E=9/5*q*q/r=',1.8d0*qbg*qbg/rhobg
    endif
    ! --- Smooth electrostatic potential ---
    call rhomom (sv_p_orhoat,  qmom,vsum)
    allocate(hpot0_rv(nbas))
    !      i = 1
    !      if (cmdopt0('--oldvc')) i = 0
    call smves(nbas , ssite , sspec ,  k1 , k2 , k3 , &
         qmom , gpot0 , vval , hpot0_rv , sgp0 , smrho , smpot , vconst &
         , smq , qsmc , fes , rhvsm , zvnsm , zsum , vesrmt , qbg )
    smag = 0
    if (nsp == 2) then
       call mshint(vol,1,n1,n2,n3,k1,k2,k3,smrho,smag,sum2)
       smag = 2*smag - smq
    endif
    if (allocated(hpot0_rv)) deallocate(hpot0_rv)
    ! --- Add smooth exchange-correlation potential ---
    if( .NOT. present(novxc_)) then
       novxc=.false.
       allocate(smvxc_zv(k1*k2*k3*nsp),smvx_zv(k1*k2*k3*nsp), &
            smvc_zv(k1*k2*k3*nsp),smexc_zv(k1*k2*k3), fxc_rv(3,nbas))
       smvxc_zv(:)=0d0
       smvx_zv(:) =0d0
       smvc_zv(:) =0d0
       smexc_zv(:)=0d0
       fxc_rv(:,:)  =0d0
       !   ... Smooth exchange-correlation potential
       smvxc_zv=0d0
       smexc_zv=0d0
       call smvxcm ( ssite , sspec , nbas , lfrce , k1 , k2 ,&
       k3 , smrho , smpot , smvxc_zv , smvx_zv , smvc_zv , &
            smexc_zv , repsm , repsmx , repsmc , rmusm , rvmusm , rvepsm &
            , fcexc0 , fcex0 , fcec0 , fcvxc0 , fxc_rv )
       if ( lfrce /= 0 ) fes = fes+fxc_rv
       if (allocated(fxc_rv)) deallocate(fxc_rv)
       if (allocated(smexc_zv)) deallocate(smexc_zv)
       if (allocated(smvc_zv)) deallocate(smvc_zv)
       if (allocated(smvx_zv)) deallocate(smvx_zv)
       if (allocated(smvxc_zv)) deallocate(smvxc_zv)
    else
       novxc=.true.
       repsm=0d0
       repsmx=0d0
       repsmc=0d0
       rmusm=0d0
       rvmusm=0d0
       fcexc0=0d0
       fcex0=0d0
       fcec0=0d0
       fcvxc0=0d0
    endif
    !! Add dipole contribution (x,y,z) to smpot (we only need <i|x|j> and so on, but because of technical reason,
    !! Calculate <i|x|j> by subtraction. nov2021 (this is stupid implementation--- See MLWF paper.
    !$$$      dipole=0
    !$$$      if(present(dipole_)) then
    !$$$         if(dipole_/=0) dipole=dipole_
    !$$$         write(6,*)' mkpot dipole=', present(dipole_),dipole_,dipole,k1,k2,k3
    !$$$      endif
    !$$$      if(dipole/=0) then
    !$$$         do isp=1,nsp
    !$$$         do i=1,k1
    !$$$         do j=1,k2
    !$$$         do k=1,k3
    !$$$            smpot(i,j,k,isp)=smpot(i,j,k,isp)+ lat_alat*sum(plat(dipole,:)
    !$$$     &         * [(i-1)/dble(k1),(j-1)/dble(k2),(k-1)/dble(k3)])
    !$$$         enddo
    !$$$         enddo
    !$$$         enddo
    !$$$         enddo
    !$$$      endif
    ! --- Make parameters for extended local orbitals ---
    call elocp() ! set ehl and rsml for local orbitals
    if(sum(lpzex)/=0) call m_bstrux_init()!computes structure constant (C_akL Eq.(38) in /JPSJ.84.034702)
    ! when we have extended local orbital.

    ! --- Make local potential at atomic sites and augmentation matrices ---
    rhobg=qbg/vol
    call locpot ( &
         sv_p_orhoat , qmom , vval , gpot0 , sv_p_osig , sv_p_otau , &
         sv_p_oppi,ohsozz,ohsopm, ppnl_rv , hab_rv , vab_rv , sab_rv , vvesat , cpnvsa , repat , &
         repatx , repatc , rmuat , rvepva , rvexva , rvecva , rvvxva , &
         rvepsa , rvvxca , valvfa , xcore , fcexca , fcexa , fceca , fcvxca &
         , sqloc , sqlocc , saloc , qval , qsc , job , rhobg , &
         nlibu, lmaxu, vorb, lldau,  novxc)!,dipole)
    if(cmdopt0('--density') .AND. master_mpi .AND. secondcall) return

    ! ... Combine spin-up and spin-down integrals
    repsm(1)  = repsm(1) + repsm(2)
    repsmx(1) = repsmx(1)+ repsmx(2)
    repsmc(1) = repsmc(1)+ repsmc(2)
    rmusm(1)  = rmusm(1) + rmusm(2)
    rvmusm(1) = rvmusm(1) + rvmusm(2)
    fcexc0(1) = fcexc0(1) + fcexc0(2)
    fcex0(1)  = fcex0(1) + fcex0(2)
    fcec0(1)  = fcec0(1) + fcec0(2)
    fcvxc0(1) = fcvxc0(1) + fcvxc0(2)
    repat(1)  = repat(1) + repat(2)
    repatx(1) = repatx(1)+ repatx(2)
    repatc(1) = repatc(1)+ repatc(2)
    rmuat(1)  = rmuat(1) + rmuat(2)
    fcexca(1) = fcexca(1) + fcexca(2)
    fcexa(1)  = fcexa(1) + fcexa(2)
    fceca(1)  = fceca(1) + fceca(2)
    fcvxca(1) = fcvxca(1) + fcvxca(2)
    ! ... Integral of valence density times estatic potential
    valves = rhvsm + vvesat

    ! ... Valence density times veff.
    !    *Associate term (n0~-n0) Ves(n0~) with local part
    !     because of the ppi matrix elements
    !    *Also add fcvxc0(1) to smooth part because
    !     rvmusm+fcvxc0 is perturbative approximation for rvmusm
    !     when cores are not treated perturbatively.
    valfsm = rhvsm + rvmusm(1) - sgp0 - vconst*qbg
    valftr = valvfa + sgp0
    valfsm = valfsm + fcvxc0(1)
    valftr = valftr - fcvxca(1)
    valvef = valfsm + valftr
    ! ... Integral of core+nucleus times estatic potential
    cpnves = zvnsm + cpnvsa
    ! ... Total xc energy and potential integral
    focexc = fcexc0(1) - fcexca(1)
    !     focex  = fcex0(1)  - fcexa(1)
    !     focec  = fcec0(1)  - fceca(1)
    focvxc = fcvxc0(1) - fcvxca(1)
    if (ipr >= 30 .AND. dabs(focexc) > 1d-6) write (stdo,850) focexc,focvxc
850 format(' foca xc integrals for spillout charge:',2f12.6)
    repsm(1) = repsm(1) + fcexc0(1)
    repsmx(1)= repsmx(1)+ fcex0(1)
    repsmc(1)= repsmc(1)+ fcec0(1)
    repat(1) = repat(1) - fcexca(1)
    repatx(1)= repatx(1)- fcexa(1)
    repatc(1)= repatc(1)- fceca(1)
    rmusm(1) = rmusm(1) + fcvxc0(1) + fcexc0(1)
    rmuat(1) = rmuat(1) - fcvxca(1) - fcexca(1)
    rhoexc = repsm(1) + repat(1)
    rhoex  = repsmx(1)+ repatx(1)
    rhoec  = repsmc(1)+ repatc(1)
    rhovxc = rmusm(1) + rmuat(1)
    ! ... Total electrostatic energy
    usm = 0.5d0*(rhvsm+zvnsm)
    uat = 0.5d0*(vvesat+cpnvsa)
    utot = usm + uat
    dq = smq+sqloc + qsmc+sqlocc + qbg -zsum
    amom = smag+saloc

    ! --- Printout ---
    if (ipr >= 20) write(stdo,'(1x)')
    if (ipr >= 30) then
       write (stdo,681)
       write (stdo,680) 'rhoval*vef ',valfsm,valftr,valvef, &
            'rhoval*ves ',rhvsm,vvesat,valves, &
            'psnuc*ves  ',zvnsm,cpnvsa,cpnves, &
            'utot       ',usm,uat,utot, &
            'rho*exc    ',repsm(1),repat(1),rhoexc, &
            'rho*vxc    ',rmusm(1),rmuat(1),rhovxc, &
            'valence chg',smq,sqloc,smq+sqloc
       if (nsp == 2) &
            write (stdo,680) 'valence mag',smag,saloc,amom
       write (stdo,680) 'core charge',qsmc,sqlocc,qsmc+sqlocc
       write (stdo,670) smq+sqloc,qsmc+sqlocc,-zsum,qbg,dq
    endif
680 format(3x,a,4x,3f17.6)
681 format(' Energy terms:',13x,'smooth',11x,'local',11x,'total')
670 format(/' Charges:  valence',f12.5,'   cores',f12.5, &
         '   nucleii',f12.5/'    hom background',f12.5, &
         '   deviation from neutrality: ',f12.5)
    if (ipl >= 1) then
       write (stdl,710) smq+sqloc,smq,sqloc,qbg,dq
710    format('fp qvl',f11.6,'  sm',f11.6,'  loc',f11.6, &
            '  qbg',f11.6,' dQ',f11.6)
       if (nsp == 2) write (stdl,711) smag+saloc,smag,saloc
711    format('fp mag ',f11.5,'  sm ',f11.5,'  loc ',f11.5)
       write (stdl,720) rhovxc,rhoexc,utot
720    format('fp pot  rvxc',f18.7,'  rexc ',f18.7,'  rves ',f16.7)
       write (stdl,721) rhoex,rhoec
721    format('fp pot  rex ',f18.7,'  rec ',f19.7)
    endif
    if (dabs(dq)>1d-3) write(stdo,"(a,f13.6)")' (warning) system not neutral, dq=',dq
    call tcx('mkpot')
  end subroutine mkpot

  subroutine dfaugm ( sv_p_osig, sv_p_otau , sv_p_oppi, ohsozz,ohsopm )
    use m_struc_def,only:s_rv1,s_cv1,s_sblock
    use m_lmfinit,only: lso,nkaph,nsp,nbas,ssite=>v_ssite,sspec=>v_sspec
    !- Allocate augmentation matrices sigma,tau,pi for all atoms
    ! ----------------------------------------------------------------------
    !o Outputs
    !o   osig  :memory allocated
    !o   otau  :memory allocated
    !o   oppi  :memory allocated
    !r Remarks
    !r   Pointers are specified as osig(itype,ibas) where
    !r     type=1: case Pkl*Pkl
    !r     type=2: case Pkl*Hsm
    !r     type=3: case Hsm*Hsm
    !r   sig and tau are l diagonal, ppi is full matrix
    !r   Thus integral (P~_kL P~_k'L' - P_kL P_k'L') is diagonal in LL',
    !r       sig(nf1,nf2,0..lmax) with lmax the l-cutoff
    !r   For sig(Pkl,Pkl), nf1=nf2==1+kmax; lmax=lmxa
    !r   For sig(Hsm,Pkl), nf1=nkaph and nf2=1+kmax; lmax=lmxh
    !r   For sig(Hsm,Hsm), nf1=nf2=nkaph; lmax = lmxh
    !l   nkaph :number of orbital types for a given L quantum no. in basis
    ! ----------------------------------------------------------------------
    implicit none
    type(s_cv1) :: sv_p_oppi(3,nbas)
    type(s_sblock):: ohsozz(3,nbas),ohsopm(3,nbas)
    type(s_rv1) :: sv_p_otau(3,nbas), sv_p_osig(3,nbas)
    integer :: ib,igetss,is,kmax,lmxa,lmxh,nelt1,nelt2,nglob,nlma,nlmh,nelt !,nso
    logical:: cmdopt0
    do  ib = 1, nbas
       is = ssite(ib)%spec
       lmxa=sspec(is)%lmxa
       lmxh=sspec(is)%lmxb
       kmax=sspec(is)%kmxt
       nlma = (lmxa+1)**2
       nlmh = (lmxh+1)**2
       if (lmxa == -1) cycle
       !   ... Case Pkl*Pkl
       nelt1 = (kmax+1)*(kmax+1)*(lmxa+1)*nsp
       nelt2 = (kmax+1)*(kmax+1)*nlma*nlma*nsp!*nspc!*nso
       allocate(sv_p_osig(1,ib)%v(abs(nelt1)))
       allocate(sv_p_otau(1,ib)%v(abs(nelt1)))
       allocate(sv_p_oppi(1,ib)%cv(nelt2))
       !   ... Case Hsm*Hsm
       nelt1 = nkaph*nkaph*(lmxh+1)*nsp
       nelt2 = nkaph*nkaph*nlmh*nlmh*nsp!*nspc !*nso
       allocate(sv_p_osig(3,ib)%v(abs(nelt1)))
       allocate(sv_p_otau(3,ib)%v(abs(nelt1)))
       allocate(sv_p_oppi(3,ib)%cv(nelt2))
       !   ...  Case Hsm*Pkl
       if (lmxh > lmxa) call rx('dfaugm: lmxh > lmxa unexpected')
       nelt1 = nkaph*(kmax+1)*(lmxh+1)*nsp
       nelt2 = nkaph*(kmax+1)*nlmh*nlma*nsp!*nspc !*nso
       allocate(sv_p_osig(2,ib)%v(abs(nelt1)))
       allocate(sv_p_otau(2,ib)%v(abs(nelt1)))
       allocate(sv_p_oppi(2,ib)%cv(nelt2))
       if(lso/=0 .OR. cmdopt0('--socmatrix')) then !spin-orbit copling matrix elements
          nelt = (kmax+1)*(kmax+1)*nlma*nlma !ohsopm (L- and L+) is irrelevant for lso=2
          allocate(ohsozz(1,ib)%sdiag(nelt,nsp), ohsopm(1,ib)%soffd(nelt,nsp)) ! Pkl*Pkl zz and pm component
          nelt = nkaph*nkaph*nlmh*nlmh !*nsp
          allocate(ohsozz(3,ib)%sdiag(nelt,nsp), ohsopm(3,ib)%soffd(nelt,nsp)) ! Hsm*Hsm
          nelt = nkaph*(kmax+1)*nlmh*nlma !*nsp
          allocate(ohsozz(2,ib)%sdiag(nelt,nsp), ohsopm(2,ib)%soffd(nelt,nsp)) !
          !     ohsozz(2,ib): Hsm*Pkl for Lz(diag)
          !     ohsopm(2,ib): Hsm*Pkl soffd(:,1) = <H|L-(isp=1,isp=2)|P>, soffd(:,2)     = <H|L+(isp=2,isp=1)|P>
          !                                                            dagger(soffd(:,2))= <P|L-(isp=1,isp=2)|H>
       endif
    enddo
  end subroutine dfaugm
  !!------------------------------------------------------
  
  subroutine m_mkpot_deallocate()
    if (allocated(vesrmt)) then
       deallocate(vesrmt,fes1_rv,ppnl_rv,sab_rv,vab_rv,hab_rv,vval,gpot0,qmom, &
            sv_p_oppi,sv_p_otau,sv_p_osig,osmpot,ohsozz,ohsopm)
    endif
  end subroutine m_mkpot_deallocate
  
end module m_mkpot
