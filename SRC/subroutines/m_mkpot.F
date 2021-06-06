      module m_mkpot
      use m_ldau,only: vorb
      use m_supot,only: k1,k2,k3
      use m_lmfinit,only: nbas,stdo,qbg=>zbak,ham_frzwf,lmaxu,nsp,nlibu,n0,nab,nppn,ctrl_lfrce
      use m_suham,only: nchan=>pot_nlma, nvl=>pot_nlml
      use m_MPItk,only: master_mpi,mlog
      use m_struc_def,only: s_rv1
      use m_density,only: osmrho, orhoat
      
      type(s_rv1),allocatable,protected  :: sv_p_oppi(:,:)
      type(s_rv1),allocatable,protected  :: sv_p_otau(:,:)
      type(s_rv1),allocatable,protected  :: sv_p_osig(:,:)
      complex(8),allocatable,protected   :: osmpot(:,:,:,:)
      real(8):: qval,protected
      real(8),allocatable,protected :: fes1_rv(:), gpot0_rv(:), hab_rv(:),ppnl_rv(:,:,:,:), qmom_rv(:)
      real(8),allocatable,protected :: sab_rv(:),  vab_rv(:),vval_rv(:),sav_rv(:)
      real(8),protected:: qsc,vconst 
      real(8),allocatable,protected :: vrmt(:)
      real(8) ,allocatable :: fes2_rv(:)
      
      integer,private:: lfrce
      contains
      subroutine m_mkpot_deallocate()
      if (allocated(vrmt)) then
         deallocate(vrmt,fes1_rv,ppnl_rv,sab_rv,vab_rv,hab_rv,vval_rv,gpot0_rv,qmom_rv,
     &    sv_p_oppi,sv_p_otau,sv_p_osig,osmpot)
      endif
      end
      
      subroutine m_mkpot_init(llmfgw)
      integer:: i,lfrzw,lfrce
      real(8) ::qbz
      
      logical :: llmfgw,fsmode,fullmesh,cmdopt0
c      qbg = ctrl_zbak(1) !homogenious background charge
!
      lfrce=ctrl_lfrce
      fullmesh = cmdopt0('--fullmesh') !full mesh stop just after do 2010 iq loop.
      if(fsmode) fullmesh=.true.
      if(cmdopt0('--band').or.fullmesh) lfrce = 0
!      
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
      allocate( vrmt(nbas))
      allocate( osmpot(k1,k2,k3,nsp)) !smooth potential without XC
      allocate( qmom_rv(nvl))
      allocate( gpot0_rv(nvl))
      allocate( vval_rv(nchan))
      allocate(  hab_rv(nab*n0*nsp*nbas))
      allocate(  vab_rv(nab*n0*nsp*nbas))
      allocate(  sab_rv(nab*n0*nsp*nbas))
      allocate( ppnl_rv(nppn,n0,nsp,nbas))
      allocate( fes1_rv(3*nbas))
c      allocate( vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu))
      call togpr()
      allocate( sv_p_osig(3,nbas), sv_p_otau(3,nbas), sv_p_oppi(3,nbas))
      call dfaugm(sv_p_osig, sv_p_otau, sv_p_oppi ) !allocation for sig,tau,ppi integrals
      
      call mkpot (lfrce,
     i    osmrho, orhoat , qbg  
     o  , osmpot, qmom_rv, vconst, vrmt, sv_p_osig, sv_p_otau, sv_p_oppi 
     o  , ppnl_rv , hab_rv , vab_rv , sab_rv , qval , qsc , gpot0_rv 
     o  , vval_rv , fes1_rv , i , vorb )
      call togpr()
      end subroutine
      
!!------------------------------------------------------
      subroutine m_mkpot_energyget(smrho_out,orhoat_out)
      use m_struc_def
      type(s_rv1):: orhoat_out(:,:)
      complex(8) :: smrho_out(:)
c      lfrzw = 0
c      if(ham_frzwf) lfrzw = 1   !freeze all augmentation wave
!! Arrays used in the generation of the potential ---
c      allocate( vrmt(nbas))
c      allocate( osmpot(k1,k2,k3,nsp)) !smooth potential without XC
c      allocate( qmom_rv(nvl))
c      allocate( gpot0_rv(nvl))
c      allocate( vval_rv(nchan))
c      allocate(  hab_rv(nab*n0*nsp*nbas))
c      allocate(  vab_rv(nab*n0*nsp*nbas))
c      allocate(  sab_rv(nab*n0*nsp*nbas))
c     allocate( ppnl_rv(nppn,n0,nsp,nbas))
      if(allocated(fes2_rv)) deallocate(fes2_rv)
      allocate( fes2_rv(3*nbas))
      call mkpot (lfrce,
     i    smrho_out, orhoat_out , qbg  
     o  , osmpot, qmom_rv, vconst, vrmt, sv_p_osig, sv_p_otau, sv_p_oppi 
     o  , ppnl_rv , hab_rv , vab_rv , sab_rv , qval , qsc , gpot0_rv 
     o  , vval_rv , fes2_rv , 0 , vorb ) !job=0
      end subroutine
      
      end module
