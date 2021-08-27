      module m_mkpot ! http://dx.doi.org/10.7566/JPSJ.84.034702
      use m_lmfinit,only: nbas,stdo,qbg=>zbak,ham_frzwf,lmaxu,nsp,nlibu,n0,nab,nppn
     & ,lfrce=>ctrl_lfrce,stdl, nchan=>pot_nlma, nvl=>pot_nlml,nkaph
      use m_struc_def,only: s_rv1,s_cv1
      use m_MPItk,only: master_mpi,mlog
!! output ------------------------------------------------------------------
      public::  m_Mkpot_init, m_Mkpot_deallocate, m_Mkpot_energyterms, m_Mkpot_novxc
      type(s_cv1),allocatable,protected,public  :: ohsozz(:,:),ohsopm(:,:) !SOC matrix 
                              ! Generated at mkpot-locpot-augmat-gaugm
      type(s_cv1),allocatable,protected,public  :: sv_p_oppi(:,:) !pi-integral Eq.(C.6)
      type(s_rv1),allocatable,protected,public  :: sv_p_otau(:,:) !tau            (C.5)
      type(s_rv1),allocatable,protected,public  :: sv_p_osig(:,:) !sigma          (C.4)
      complex(8),allocatable,protected ,public  :: osmpot(:,:,:,:)!0th component of Eq.(34)
!
      type(s_cv1),allocatable,protected,public  :: sv_p_oppix(:,:) !pi-integral without xc term
      type(s_rv1),allocatable,protected,public  :: sv_p_otaux(:,:) !tau
      type(s_rv1),allocatable,protected,public  :: sv_p_osigx(:,:)
      complex(8),allocatable,protected ,public  :: spotx(:,:,:,:)
!
      real(8),protected,public:: utot,rhoexc,xcore,valvef,amom,  valves,cpnves,rhovxc !energy terms
      real(8),allocatable,protected,public:: fes1_rv(:), fes2_rv(:) !force terms
      real(8),allocatable,protected,public:: hab_rv(:), sab_rv(:), ppnl_rv(:,:,:,:), qmom(:),vesrmt(:)
      real(8),protected,public:: qval,vconst,qsc
!!---------------------------------------------------------------------------      
      private
      real(8),allocatable,protected,private :: gpot0(:),vval(:),vab_rv(:)
      contains
      subroutine m_mkpot_novxc()
      use m_supot,only: k1,k2,k3
      use m_density,only: osmrho, orhoat !main input density
      integer:: i,lfrzw
      type(s_cv1),allocatable  :: ohsozzx(:,:),ohsopmx(:,:) 
      write(stdo,"(a)")' m_mkpot_novxc: Making one-particle potential without XC part ...'
      allocate( sv_p_osigx(3,nbas), sv_p_otaux(3,nbas), sv_p_oppix(3,nbas))
      allocate( ohsozzx(3,nbas), ohsopmx(3,nbas)) !dummy
      allocate( spotx(k1,k2,k3,nsp)) !smooth potential without XC
      call dfaugm(sv_p_osigx, sv_p_otaux, sv_p_oppix, ohsozzx,ohsopmx)!for sig,tau,ppi without XC(LDA)
      allocate( vesrmt(nbas))
      allocate(  qmom(nvl)) !rhomom
      allocate( ppnl_rv(nppn,n0,nsp,nbas))
      allocate(  hab_rv(nab*n0*nsp*nbas))
      allocate(  vab_rv(nab*n0*nsp*nbas))
      allocate(  sab_rv(nab*n0*nsp*nbas))
      allocate(  gpot0(nvl))
      allocate(  vval(nchan))
      allocate(  fes1_rv(3*nbas))
      spotx=0d0
      lfrzw = 0
      if(ham_frzwf) lfrzw = 1   !freeze all augmentation wave
      i = 1 + 10*lfrzw + 100    ! Adding 100 means excluding XC(LDA) part. nolxc=T
!     We obtain sv_p_osigx, sv_p_otaux, sv_p_oppix, smpotx  without XC
      call mkpot (i, osmrho , orhoat , 
     o     spotx, sv_p_osigx , sv_p_otaux , sv_p_oppix, fes1_rv, ohsozzx,ohsopmx)
      deallocate(vesrmt,qmom,ppnl_rv,hab_rv,vab_rv,sab_rv,gpot0,vval,fes1_rv)
      end subroutine 
!!ssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine m_mkpot_init(llmfgw)
      use m_supot,only: k1,k2,k3
      use m_density,only: osmrho, orhoat !main input density
      use m_lmfinit,only: lso,nbas, ssite=>v_ssite, sspec=>v_sspec, nlibu,lmaxu,lldau,nsp,lat_alat,lxcf,lpzex 
      integer:: i,lfrzw,is,ib,kmax,lmxa,lmxh,nelt2,nlma,nlmh!,lfrce
      real(8) ::qbz
      logical :: llmfgw,fsmode,fullmesh,cmdopt0
c      lfrce=ctrl_lfrce
c      fullmesh = cmdopt0('--fullmesh') !full mesh stop just after do 2010 iq loop.
c      if(fsmode) fullmesh=.true.
c      if(cmdopt0('--band').or.fullmesh) lfrce = 0
      lfrzw = 0
      if(ham_frzwf) lfrzw = 1   !freeze all augmentation wave
!! Make the potential and total energy terms for given density (smrho,rhoat,qbg)  ---
!! mkpot calls locpot. and locpot calls augmat. augmat calculates sig,tau,ppi.
      i = 1 + 10*lfrzw   
      if (cmdopt0('--wrhomt')) then
         i = i + 10000 !write rhomt mode
      else if (cmdopt0('--wpotmt')) then
         i = i + 20000 !write potmt mode
      endif
      if(llmfgw) i = i + 10000 !GW driver mode
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
      call togpr()
      allocate( sv_p_osig(3,nbas), sv_p_otau(3,nbas), sv_p_oppi(3,nbas))
      allocate( ohsozz(3,nbas), ohsopm(3,nbas))
      call dfaugm(sv_p_osig,sv_p_otau,sv_p_oppi,ohsozz,ohsopm) !allocation for sig,tau,ppi integrals
      call mkpot (i, osmrho , orhoat, 
     o     osmpot, sv_p_osig , sv_p_otau , sv_p_oppi, fes1_rv, ohsozz,ohsopm)
      call togpr()
      end subroutine
!!ssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine m_mkpot_energyterms(smrho_out,orhoat_out)
      use m_struc_def
      type(s_rv1):: orhoat_out(:,:)
      complex(8) :: smrho_out(:)
      if(allocated(fes2_rv)) deallocate(fes2_rv)
      allocate( fes2_rv(3*nbas))
         call mkpot (0, smrho_out , orhoat_out , 
     o     osmpot, sv_p_osig , sv_p_otau , sv_p_oppi, fes2_rv, ohsozz,ohsopm)
      end subroutine
!!------------------------------------------------------
      subroutine m_mkpot_deallocate()
      if (allocated(vesrmt)) then
         deallocate(vesrmt,fes1_rv,ppnl_rv,sab_rv,vab_rv,hab_rv,vval,gpot0,qmom,
     &    sv_p_oppi,sv_p_otau,sv_p_osig,osmpot,ohsozz,ohsopm)
      endif
      end
!!------------------------------------------------------
      subroutine mkpot (
     i    job, smrho, sv_p_orhoat,  
     o    smpot, sv_p_osig, sv_p_otau, sv_p_oppi, fes, ohsozz,ohsopm ) 
      use m_lmfinit,only: lso,nbas, ssite=>v_ssite, sspec=>v_sspec, nlibu,lmaxu,lldau,nsp,lat_alat,lxcf,lpzex 
      use m_lattic,only: lat_plat,lat_qlat, lat_vol
      use m_supot,only: k1,k2,k3,lat_nabc
      use m_MPItk,only: master_mpi
      use m_struc_def  
      use m_ldau,only: vorb !input. 'U-V_LDA(couter term)' of LDA+U 
      use m_bstrux,only: m_bstrux_init
      use m_elocp,only: elocp
!Note other output in module area
c for job=0      
Co         utot   = total electrostatic energy
Co         valves = valence rho * estat potential
Co         cpnves = core+nuc * estat potential
Co         rhoexc = rho * exc
Co         rhovxc = rho * vxc
Co         xcore  = rhoc * total potential
Co         valvef = smrhov * vsm + sum_ib valvef_ib
Co           valvef_ib = rhov * vtrue - smrhov * vsm)_ib
Co         amom   = system magnetic moment
!      
!! documents below are under construction.
C- Make the potential from the density.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ssite :struct containing site-specific information
Ci     Passed to: rhomom smves smvxcm locpot
Ci   sspec :struct containing species-specific information
Ci     Elts read:
Ci     Passed to: rhomom smves smvxcm locpot
Ci   lfrce :nonzero =>  contribution to forces
Ci   lcplxp=1 only ::0 if ppi is real; 1 if ppi is complex
Ci   k1,k2,k3 dimensions of smrho for smooth crystal density
Ci   smrho :smooth crystal density, on a uniform mesh
Ci   orhoat:local atomic densities (true and smooth parts)
Ci   qbg   :homogeneous background charge
Ci   job   :1s digit
Ci         : 0 stops after energy terms
Ci         : 1 makes potpars also
Ci         :10s digit
Ci         :1 suppress updating potential used to make potpars
Ci         :100s digit
Ci         :1 exclude exchange-correlation potential
Cixxx         :1000s digit
Cixxx         :1 Make rveps and rvvxc
Ci         :10000s digit
Ci         :1 write sphere density for each site ib to file rhoMT.ib
Ci ... The following are LDA+U inputs
Ci   vorb  :orbital dependent potential
Ci   nlibu :number of U blocks  used to dimension vorb
Ci   lmaxu :max l for U blocks  used to dimension vorb
Ci   lldau :lldau(ib)=0 => no U on this site otherwise
Ci         :U on site ib with dmat beginning at dmats(*,lldau(ib))
Co Outputs:
Co   smpot :smooth potential on a uniform mesh:
Co         :Ves~ + vxc = Ves(n0 + compensating gaussians) + vxc
Co   qmom  :multipole moments of valence sphere densities
Co   vconst:additional constant potential term
Co   vesrmt:electrostatic potential at rmt
Co   qval  :total valence charge, including semicore states
Co   qsc   :total charge from semicore states (local orbitals)
Co   osig,otau,oppis: augmentation matrices
Co   ppn   :nmto-like potential parameters
Co   hab,vab,sab: augmentation matrices for local parts of the ham.
Co         :See Remarks in augmat.f for their generation and doc.
Co   vval  :coffs to YL expansion of es potential at MT boundary
Co   gpot0 :integrals of gaussians times electrostatic potential
Co   fes   :contribution to the force from electrostatic + xc potential
Co   sham->eterms various integrals for the total energy are stored:
Cl Local variables
Cl   rvmusm:int rhosm * vxc(rhosm+smcor1) where smcor1 is portion
Cl         :of smooth core density treated nonperturbatively.
Clxxx   rvepsm:int rhosm * exc(rhosm+smcor1) where smcor1 is portion
Clxxx         :of smooth core density treated nonperturbatively.
ccCl   rvepsv:integral of valence density times exc(valence density) !removed now
ccCl   rvvxcv:integral of valence density times vxc(valence density) !removed now
Cl   rhvsm :integral n0~ phi0~
Cl   sgp0  :compensating gaussians * sm-Ves = int (n0~-n0) phi0~
Cl   valfsm:rhvsm - sgp0 + rvmusm + fcvxc0 = n0 (phi0~ + Vxc(n0))
Cl         :NB sgp0 associated with the local parts, as is done in
Cl         :making the ppi matrix elements.
Cl   valftr:local contribution to density * veff
Cl         := sgp0 - fcvxca + sum_ib valvfa_ib
Cl   valvfa:generated by locpot, called valvef there.  Sum_ib of:
Cl         :vefv1 - vefv2
Cl         := int[rho1*(v1-2*Z/r) - (rho2*v2)] - (n0~-n0)*gpotb - focvxc
Cl   lso   :nonzero => spin orbit coupling
Cr Remarks
Cr *The total density is a sum of three terms,
Cr
Cr    n0(mesh) + sum_RL (n_RL(r) - n0_RL(r))
Cr
Cr  The first term is the smooth density on a mesh of points; the
Cr  second is the true density and is defined on a radial mesh for each
Cr  sphere; the last is the 1-center expansion of the smooth density on
Cr  the radial mesh.  (Note: because of l-truncation, n0_R(r) is not
Cr  identical to the one-center expansion of n0(mesh).  The sum of the
Cr  three terms converges rapidly with l because errors in n_R(r) are
Cr  mostly canceled by errors in n0_R(r).)
Cr
Cr *Computation of the electrostatic energy:
Cr  We add and subtract a set of compensating gaussian orbitals
Cr
Cr    n0 + sum_RL Q_RL g_RL + sum_RL (n_RL(r) - n0_RL(r) - Q_RL g_RL)
Cr
Cr  which render the integral of the local part (the last 3 terms)
Cr  zero in each RL channel.  The g_RL must be localized enough that
Cr  their spillout beyond the MT radius is negligible.
Cr
Cr  We define
Cr
Cr    n0~ = n0 + compensating gaussians
Cr
Cr  In the interstitial, the electrostatic potential of n0~ is the true
Cr  estat potential.  The potential of n0 is called phi0 and the
Cr  potential of n0~ is called phi0~.  The total electrostatic energy
Cr  is computed as
Cr    the electrostatic energy of  n0~  +
Cr    the electrostatic energy of (neutral) local parts
Cr
Cr  The first term is computed in subroutine smves;
Cr  the second term is computed in subroutine locpot.
Cu Updates
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   27 Apr 05 Added LDA+U stuff (Lambrecht)
Cu   24 Dec 04 Changes for full L.S coupling
Cu    1 Sep 04 Adapted mkpot to handle complex ppi; fold so into ppi
Cu   15 Jul 04 (Chantis) Matrix elements Lz.Sz for spin-orbit coupling
Cu   14 Jan 02 rvexv and rvecv (T. Miyake)
Cu   17 Sep 01 Returns qsc.  Altered argument list.
Cu   24 Aug 01 Extended to local orbitals.  Altered argument list.
Cu             Local potentials now have correct boundary conditions
Cu   15 Aug 01 Generates rvepsv and rvvxcv.  Changed call to locpot.
Cu   20 Apr 01 Generates vesrmt
Cu   18 Apr 01 Added ability to exclude exchange-correlation potential
Cu   20 Feb 01 Added ppn to potential parameters generated
Cu   15 Jun 00 spin polarized
Cu   22 Apr 00 Adapted from nfp mk_potential.f
C ------------------------------------------------------------
C to do:
C 1. check that rhov*vxc approximately the same in locpot as smvxc
C 2. ditto for rhov*exc; need finish making rhov*exc in smvxc
C 3. ? See about understanding and changing foxexc both in locpot
C    and smvxc.  The total energy seems to come out just right,
C    but focexc isn't real correction to rhov*exc.  Does it matter?
C    Maybe can't use lfoca=2 is all.
C 4. enable 1000s digit job to pass through
C 5. change locpot so that rvvxcv etc is make only if appropriate
C    1000s job bit set.
C 6. Make <vxc_nn> and compare to integrals made here.
C    If problems, Maybe even look at sm part only, or local part only.
      implicit none
      integer job,i1,i2,i3,ispec
      type(s_rv1) :: sv_p_orhoat(*)
      type(s_cv1) :: sv_p_oppi(*),ohsozz(*),ohsopm(*)
      type(s_rv1) :: sv_p_otau(*)
      type(s_rv1) :: sv_p_osig(*)
      real(8)::  fes(3,nbas) 
      complex(8):: smrho(k1,k2,k3,nsp),smpot(k1,k2,k3,nsp)
      logical cmdopt0
      character*80 outs
      integer i,lgunit,nglob,ipl,ipr,iprint,n1,n2,n3,ngabc(3),lxcfun,isw,isum
      real(8) ,allocatable :: fxc_rv(:,:)
      real(8) ,allocatable :: hpot0_rv(:)
      complex(8) ,allocatable :: smvxc_zv(:)
      complex(8) ,allocatable :: smvx_zv(:)
      complex(8) ,allocatable :: smvc_zv(:)
      complex(8) ,allocatable :: smexc_zv(:)
      double precision dq,focexc,cpnvsa, 
     .focvxc,qsmc,smq,smag,sum2,rhoex,rhoec,rhvsm,sgp0,
     .sqloc,sqlocc,saloc,uat,usm,valfsm,valftr, 
     .rvepva,rvexva,rvecva,rvvxva,rvepsa,rvvxca,valvfa,vvesat, 
     .vol,vsum,zsum,zvnsm,rvepsv(2),rvexv(2),rvecv(2), 
     .rvvxcv(2),fcexc0(2),rveps(2),rvvxc(2),fcex0(2),fcec0(2),
     .fcexca(2),fcexa(2),fceca(2),rvmusm(2),rvepsm(2),rmusm(2),
     .vxcavg(2),fcvxc0(2),fcvxca(2),repat(2),repatx(2),repatc(2),
     .rmuat(2),repsm(2),repsmx(2),repsmc(2),rhobg,pi
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      integer:: itest,wdummy,nnn
      complex(8),allocatable:: smpotbk(:,:,:),smpotbkx(:,:,:)
      real(8),parameter:: minimumrho=1d-14
      real(8)::sss,smmin(2),srshift
      real(8):: plat(3,3),alat
      integer:: ifi,isp
      character strn*120
      logical:: secondcall=.false.
      
      call tcn('mkpot')
!! new density mode
      if(cmdopt0('--density').and.master_mpi.and.secondcall) then
        plat =lat_plat
        alat =lat_alat
        open(newunit=ifi,file='smrho.xsf')
        do isp = 1, nsp
          write(ifi,'("CRYSTAL")')
          write(ifi,'("PRIMVEC")')
          write(ifi,'(3f10.5)') ((plat(i1,i2)*alat*0.529177208,i1=1,3)
     .     ,i2=1,3)
          write(ifi,'("PRIMCOORD")')
          write(ifi,'(2i5)') nbas,1
          do i = 1, nbas
            ispec=ssite(i)%spec
            write(ifi,'(i4,2x,3f10.5)') int(sspec(ispec)%z),(ssite(i)%pos(i2)*alat*0.529177208,i2=1,3)
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
c      if (lso/=0 .and. lcplxp==0) call rx('mkpot: incompatible lso,lcplxp')
c      if (isum(nbas,lldau,1)/=0 .and. lcplxp==0) call rx('mkpot: incompatible ldau,lcplxp')
      lxcfun = lxcf 
      pi   = 4d0*datan(1d0)
      ngabc=lat_nabc
      vol=lat_vol
C$$$ enforce_positive_smrho() is now used only in smvxc2 Nov30 2010.
C$$$      if(enforce_positive_smrho()) then
C$$$!!== negative smrho check==
C$$$!! We need full mesh method(all G numner is needed) or so,
C$$$!! to assure positive definite, xcpbe seems to interpret abs(rho) when rho is negative.
C$$$!!  ---> see xcpbe.F90 (especially, invcb).
C$$$       do isp=1,nsp
C$$$        nnn=0
C$$$        smmin(isp)=0d0
C$$$        do i1=1,n1
C$$$        do i2=1,n2
C$$$        do i3=1,n3
C$$$c          print *,i1,i2,i3,isp
C$$$c          print *,i1,i2,i3,isp,smrho(i1,i2,i3,isp)
C$$$          sss=dreal(smrho(i1,i2,i3,isp))
C$$$          if(sss<0d0) then
C$$$            nnn=nnn+1
C$$$            if(sss<smmin(isp)) then
C$$$              smmin(isp)=sss
C$$$            endif
C$$$          endif
C$$$        enddo
C$$$        enddo
C$$$        enddo
C$$$        if(nnn>0) then
C$$$          write(6,*) 'mkpot negative smrho; isp,number,min(smrho)='
C$$$     &    ,isp,nnn,smmin(isp)
C$$$        else
C$$$          write(6,*) 'all smrho is positive for isp=',isp
C$$$        endif
C$$$       enddo
C$$$       if(sum(smmin(1:nsp))/=0d0) then
C$$$         srshift = minimumrho + maxval(abs(smmin(1:nsp)))
C$$$         smrho(:,:,:,1:nsp) = smrho(:,:,:,1:nsp) + srshift
C$$$         print *,'enforce positive smrho, to which we add srshift=',srshift
C$$$       endif
C$$$      endif

C --- Printout for smooth background charge ---
      if (qbg .ne. 0) then
        rhobg = (3d0/4d0/pi*vol)**(1d0/3d0)
        call info5(20,1,0,' Energy for background charge'//
     .  ' q=%d, radius r=%;3d :  E = 9/5*q*q/r = %;4d', qbg,rhobg,1.8d0*qbg*qbg/rhobg,0,0)
      endif

C --- Smooth electrostatic potential ---
      call rhomom (sv_p_orhoat,  qmom,vsum)
      vconst = -vsum/vol
C      if (ipr .ge. 40) write (stdo,334) vsum,vconst
C  334 format(' vsum=',f14.6,'   vconst=',f14.6)
      allocate(hpot0_rv(nbas))

      i = 1
      if (cmdopt0('--oldvc')) i = 0
      call smves ( i , nbas , ssite , sspec ,  k1 , k2 , k3 ,  !slat ,
     .qmom , gpot0 , vval , hpot0_rv , sgp0 , smrho , smpot , vconst 
     ., smq , qsmc , fes , rhvsm , zvnsm , zsum , vesrmt , qbg )
      smag = 0
      if (nsp .eq. 2) then
        call mshint(vol,1,n1,n2,n3,k1,k2,k3,smrho,smag,sum2)
        smag = 2*smag - smq
      endif
      if (allocated(hpot0_rv)) deallocate(hpot0_rv)

C --- Add smooth exchange-correlation potential ---
      if (mod(job,1000) .lt. 100) then
        allocate(smvxc_zv(k1*k2*k3*nsp),smvx_zv(k1*k2*k3*nsp),
     &  smvc_zv(k1*k2*k3*nsp),smexc_zv(k1*k2*k3), fxc_rv(3,nbas))
        smvxc_zv(:)=0d0
        smvx_zv(:) =0d0
        smvc_zv(:) =0d0
        smexc_zv(:)=0d0
        fxc_rv(:,:)  =0d0
C   ... Smooth exchange-correlation potential
        smvxc_zv=0d0
        smexc_zv=0d0
        call smvxcm ( ssite , sspec , nbas , lfrce , k1 , k2 , !slat , 
     .  k3 , smrho , smpot , smvxc_zv , smvx_zv , smvc_zv , 
     .  smexc_zv , repsm , repsmx , repsmc , rmusm , rvmusm , rvepsm 
     .  , fcexc0 , fcex0 , fcec0 , fcvxc0 , fxc_rv )
        if ( lfrce .ne. 0 ) fes = fes+fxc_rv
        if (allocated(fxc_rv)) deallocate(fxc_rv)
        if (allocated(smexc_zv)) deallocate(smexc_zv)
        if (allocated(smvc_zv)) deallocate(smvc_zv)
        if (allocated(smvx_zv)) deallocate(smvx_zv)
        if (allocated(smvxc_zv)) deallocate(smvxc_zv)
      else
        call dpzero(repsm,2)
        call dpzero(repsmx,2)
        call dpzero(repsmc,2)
        call dpzero(rmusm,2)
        call dpzero(rvmusm,2)
        call dpzero(fcexc0,2)
        call dpzero(fcex0,2)
        call dpzero(fcec0,2)
        call dpzero(fcvxc0,2)
      endif

C --- Make parameters for extended local orbitals ---
      call elocp() !nbas,nsp,ssite,sspec,job) !set ehl and rsml for local orbitals
      if(sum(lpzex)/=0) call m_bstrux_init() !We make structure constant (C_akL Eq.(38) in /JPSJ.84.034702)
                                             !when we have extended local orbital.

C --- Make local potential at atomic sites and augmentation matrices ---
      rhobg=qbg/vol
      call locpot ( !nbas , nsp , lso ,  ssite , sspec !, slat lcplxp ,
     .    sv_p_orhoat , qmom , vval , gpot0 , sv_p_osig , sv_p_otau , 
     .    sv_p_oppi,ohsozz,ohsopm, ppnl_rv , hab_rv , vab_rv , sab_rv , vvesat , cpnvsa , repat , 
     .    repatx , repatc , rmuat , rvepva , rvexva , rvecva , rvvxva , 
     .    rvepsa , rvvxca , valvfa , xcore , fcexca , fcexa , fceca , fcvxca 
     .    , sqloc , sqlocc , saloc , qval , qsc , job , rhobg , 
     .    nlibu , lmaxu , vorb , lldau )
!!
      if(cmdopt0('--density').and.master_mpi.and.secondcall) then
        return
      endif

C ... Combine spin-up and spin-down integrals
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

C ... Integral of valence density times estatic potential
      valves = rhvsm + vvesat

C ... Valence density times veff.
C    *Associate term (n0~-n0) Ves(n0~) with local part
C     because of the ppi matrix elements
C    *Also add fcvxc0(1) to smooth part because
C     rvmusm+fcvxc0 is perturbative approximation for rvmusm
C     when cores are not treated perturbatively.
      valfsm = rhvsm + rvmusm(1) - sgp0 - vconst*qbg
      valftr = valvfa + sgp0
      valfsm = valfsm + fcvxc0(1)
      valftr = valftr - fcvxca(1)
      valvef = valfsm + valftr
C ... Integral of core+nucleus times estatic potential
      cpnves = zvnsm + cpnvsa
C ... Total xc energy and potential integral
      focexc = fcexc0(1) - fcexca(1)
C     focex  = fcex0(1)  - fcexa(1)
C     focec  = fcec0(1)  - fceca(1)
      focvxc = fcvxc0(1) - fcvxca(1)
      if (ipr.ge.30 .and. dabs(focexc).gt.1d-6)
     .write (stdo,850) focexc,focvxc
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

C ... Total electrostatic energy
      usm = 0.5d0*(rhvsm+zvnsm)
      uat = 0.5d0*(vvesat+cpnvsa)
      utot = usm + uat
      dq = smq+sqloc + qsmc+sqlocc + qbg -zsum
      amom = smag+saloc

C --- Printout ---
      if (ipr .ge. 20) write(stdo,'(1x)')
      if (ipr .ge. 30) then
        write (stdo,681)
        write (stdo,680) 'rhoval*vef ',valfsm,valftr,valvef,
     .  'rhoval*ves ',rhvsm,vvesat,valves,
     .  'psnuc*ves  ',zvnsm,cpnvsa,cpnves,
     .  'utot       ',usm,uat,utot,
     .  'rho*exc    ',repsm(1),repat(1),rhoexc,
     .  'rho*vxc    ',rmusm(1),rmuat(1),rhovxc,
     .  'valence chg',smq,sqloc,smq+sqloc
        if (nsp .eq. 2)
     .  write (stdo,680) 'valence mag',smag,saloc,amom
        write (stdo,680) 'core charge',qsmc,sqlocc,qsmc+sqlocc
        write (stdo,670) smq+sqloc,qsmc+sqlocc,-zsum,qbg,dq
      endif
  680 format(3x,a,4x,3f17.6)
  681 format(' Energy terms:',13x,'smooth',11x,'local',11x,'total')
  670 format(/' Charges:  valence',f12.5,'   cores',f12.5,
     .'   nucleii',f12.5/'    hom background',f12.5,
     .'   deviation from neutrality: ',f12.5)
      if (ipl .ge. 1) then
        write (stdl,710) smq+sqloc,smq,sqloc,qbg,dq
  710   format('fp qvl',f11.6,'  sm',f11.6,'  loc',f11.6,
     .  '  qbg',f11.6,' dQ',f11.6)
        if (nsp .eq. 2) write (stdl,711) smag+saloc,smag,saloc
  711   format('fp mag ',f11.5,'  sm ',f11.5,'  loc ',f11.5)
        write (stdl,720) rhovxc,rhoexc,utot
  720   format('fp pot  rvxc',f18.7,'  rexc ',f18.7,'  rves ',f16.7)
        write (stdl,721) rhoex,rhoec
  721   format('fp pot  rex ',f18.7,'  rec ',f19.7)
      endif
      if (dabs(dq)>1d-3) write(stdo,"(a,f13.6)")
     .     ' (warning) system not neutral, dq=',dq
      call tcx('mkpot')
      end subroutine mkpot
     
      subroutine dfaugm ( sv_p_osig, sv_p_otau , sv_p_oppi, ohsozz,ohsopm )
      use m_struc_def,only:s_rv1,s_cv1
      use m_lmfinit,only: lso,nkaph,nsp,nbas,ssite=>v_ssite,sspec=>v_sspec 
C- Allocate augmentation matrices sigma,tau,pi for all atoms
C ----------------------------------------------------------------------
Co Outputs
Co   osig  :memory allocated
Co   otau  :memory allocated
Co   oppi  :memory allocated
Cr Remarks
Cr   Pointers are specified as osig(itype,ibas) where
Cr     type=1: case Pkl*Pkl
Cr     type=2: case Pkl*Hsm
Cr     type=3: case Hsm*Hsm
Cr   sig and tau are l diagonal, ppi is full matrix
Cr   Thus integral (P~_kL P~_k'L' - P_kL P_k'L') is diagonal in LL',
Cr       sig(nf1,nf2,0..lmax) with lmax the l-cutoff
Cr   For sig(Pkl,Pkl), nf1=nf2==1+kmax; lmax=lmxa
Cr   For sig(Hsm,Pkl), nf1=nkaph and nf2=1+kmax; lmax=lmxh
Cr   For sig(Hsm,Hsm), nf1=nf2=nkaph; lmax = lmxh
Cl   nkaph :number of orbital types for a given L quantum no. in basis
C ----------------------------------------------------------------------
      implicit none
      type(s_cv1) :: sv_p_oppi(3,nbas),ohsozz(3,nbas),ohsopm(3,nbas)
      type(s_rv1) :: sv_p_otau(3,nbas), sv_p_osig(3,nbas)
      integer ib,igetss,is,kmax,lmxa,lmxh,nelt1,nelt2,nglob,nlma,nlmh!,nso
      do  ib = 1, nbas
        is = int(ssite(ib)%spec)
        lmxa=sspec(is)%lmxa
        lmxh=sspec(is)%lmxb
        kmax=sspec(is)%kmxt
        nlma = (lmxa+1)**2
        nlmh = (lmxh+1)**2
        if (lmxa .eq. -1) cycle
C   ... Case Pkl*Pkl
        nelt1 = (kmax+1)*(kmax+1)*(lmxa+1)*nsp
        nelt2 = (kmax+1)*(kmax+1)*nlma*nlma*nsp!*nspc!*nso
        allocate(sv_p_osig(1,ib)%v(abs(nelt1)))
        allocate(sv_p_otau(1,ib)%v(abs(nelt1)))
        allocate(sv_p_oppi(1,ib)%cv(nelt2))
C   ... Case Hsm*Hsm
        nelt1 = nkaph*nkaph*(lmxh+1)*nsp
        nelt2 = nkaph*nkaph*nlmh*nlmh*nsp!*nspc !*nso
        allocate(sv_p_osig(3,ib)%v(abs(nelt1)))
        allocate(sv_p_otau(3,ib)%v(abs(nelt1)))
        allocate(sv_p_oppi(3,ib)%cv(nelt2))
C   ...  Case Hsm*Pkl
        if (lmxh.gt.lmxa) call rx('dfaugm: lmxh.gt.lmxa unexpected')
        nelt1 = nkaph*(kmax+1)*(lmxh+1)*nsp
        nelt2 = nkaph*(kmax+1)*nlmh*nlma*nsp!*nspc !*nso
        allocate(sv_p_osig(2,ib)%v(abs(nelt1)))
        allocate(sv_p_otau(2,ib)%v(abs(nelt1)))
        allocate(sv_p_oppi(2,ib)%cv(nelt2))
!        
        if(lso/=0) then         !spin-orbit copling matrix elements
           nelt2 = (kmax+1)*(kmax+1)*nlma*nlma*nsp
           allocate(ohsozz(1,ib)%cv(nelt2), ohsopm(1,ib)%cv(nelt2)) ! Pkl*Pkl zz and pm component
           nelt2 = nkaph*nkaph*nlmh*nlmh*nsp
           allocate(ohsozz(3,ib)%cv(nelt2), ohsopm(3,ib)%cv(nelt2)) ! Hsm*Hsm
           nelt2 = nkaph*(kmax+1)*nlmh*nlma*nsp
           allocate(ohsozz(2,ib)%cv(nelt2), ohsopm(2,ib)%cv(nelt2)) ! Hsm*Pkl
        endif
      enddo
      end subroutine dfaugm
      end module m_mkpot
