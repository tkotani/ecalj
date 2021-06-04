      module m_lmfgw
      use m_density,only: zv_a_osmrho=>osmrho,sv_a_oorhat=>orhoat
      use m_lmfinit,only:nbas,stdo,qbg=>zbak,ham_frzwf,lmaxu,nsp,nlibu,n0,nab,nppn
      use m_supot,only: k1,k2,k3
      use m_struc_def,only: s_rv1
      use m_MPItk,only: master_mpi,mlog
      use m_suham,only: nchan=>pot_nlma, nvl=>pot_nlml
      
      integer,protected:: jobgw=-999
      type(s_rv1),allocatable,protected  :: sv_p_oppix(:,:)
      type(s_rv1),allocatable,protected  :: sv_p_otaux(:,:)
      type(s_rv1),allocatable,protected  :: sv_p_osigx(:,:)
      complex(8),allocatable,protected   :: spotx(:,:,:,:)
      
      contains
      subroutine m_lmfgw_init()
      integer:: i,lfrzw
      real(8) ::qbz,vconst,vrmt(nbas),qval,qsc!,qbg
      complex(8),allocatable:: vorb(:,:,:,:)
      real(8) ,allocatable :: fes1_rv(:), gpot0_rv(:), hab_rv(:),ppnl_rv(:,:,:,:), qmom_rv(:)
      real(8) ,allocatable :: sab_rv(:), vab_rv(:),vval_rv(:),sav_rv(:)
c      qbg = ctrl_zbak(1) !homogenious background charge
      lfrzw = 0
      if(ham_frzwf) lfrzw = 1  !freeze all augmentation wave
      if(master_mpi) then
         write(stdo,*) ' lmfgw: input one of the following jobs:'
         write(stdo,*) '   0 : init mode; creates files SYMOPS, LATTC, CLASS, NLAindx, ldima'
         write(stdo,*) '   1 : GW setup mode'
         write(stdo,*) ' job 0 or 1?'
         read (5,*) jobgw
      endif
      call mpibc1(jobgw,1,2, mlog,'bndfp','jobgw')
      if(.not.( jobgw==0.or.jobgw==1)) call rx(' lmf-MPIK is only for jobgw=0,1 now')
      if(jobgw/=0) then
         write(stdo,"(a)")' Making one-particle potential without XC part ...'
         allocate( sv_p_osigx(3,nbas), sv_p_otaux(3,nbas), sv_p_oppix(3,nbas))
         call dfaugm(sv_p_osigx, sv_p_otaux, sv_p_oppix) !for sig,tau,ppi without XC(LDA)
         allocate( spotx(k1,k2,k3,nsp)) !smooth potential without XC
         allocate( vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu))
         allocate(  qmom_rv(nvl))
         allocate( gpot0_rv(nvl))
         allocate( vval_rv(nchan))
         allocate(  hab_rv(nab*n0*nsp*nbas))
         allocate(  vab_rv(nab*n0*nsp*nbas))
         allocate(  sab_rv(nab*n0*nsp*nbas))
         allocate( ppnl_rv(nppn,n0,nsp,nbas))
         allocate( fes1_rv(3*nbas))
         spotx=0d0
         i = 1 + 10*lfrzw + 100 ! Adding 100 means excluding XC(LDA) part. nolxc=T
!     We obtain  sv_p_osigx , sv_p_otaux , sv_p_oppix, smpotx  without XC by the 'call mkpot'.
!     Other output are overwritten by the 2nd call of mkpot 
         call mkpot (0,         !nbas , ssite , sspec ,   0 , lcplxp , k1 , k2 , k3 ,
     i        zv_a_osmrho , sv_a_oorhat , qbg  
     o        , spotx, qmom_rv , vconst , vrmt , sv_p_osigx , sv_p_otaux , sv_p_oppix 
     o        , ppnl_rv , hab_rv , vab_rv , sab_rv , qval , qsc , gpot0_rv 
     o        , vval_rv , fes1_rv , i , vorb ) 
         deallocate(vorb,qmom_rv,gpot0_rv,vval_rv,hab_rv,vab_rv,sab_rv,ppnl_rv,fes1_rv)
      endif
      end subroutine
      end module
