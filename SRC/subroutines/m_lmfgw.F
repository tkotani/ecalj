      module m_lmfgw
      public:: m_lmfgw_init
      private
      contains
      subroutine m_lmfgw_init() !writing files for GW calculaiton
      use m_lmfinit,only:nbas,stdo,qbg=>zbak,ham_frzwf,lmaxu,nsp,nlibu,n0,nab,nppn,sspec=>v_sspec,ssite=>v_ssite,
     & nkap0,alat=>lat_alat,ham_ldham,iprmb
      use m_supot,only: k1,k2,k3
      use m_MPItk,only: master_mpi,mlog
      use m_lmfinit,only: nchan=>pot_nlma, nvl=>pot_nlml
      use m_lattic,only: qlat=>lat_qlat, plat=>lat_plat
      use m_suham,only: ndham=>ham_ndham
      integer:: i,lfrzw,i1,i_copy_size,ib,ifildima,ificlass,ipqn,is,isp,konf,l,lmaxa,ifilnlaindx,
     & lmxax,ifinlaindx,ndima,npqn,mxint,nat
      real(8) ::qbz,vconst,vrmt(nbas),qval,qsc
      complex(8),allocatable:: vorb(:,:,:,:)
      real(8) ,allocatable :: fes1_rv(:), gpot0_rv(:), hab_rv(:),ppnl_rv(:,:,:,:), qmom_rv(:)
      real(8) ,allocatable :: sab_rv(:), vab_rv(:),vval_rv(:),sav_rv(:),bas(:,:)
      integer,allocatable:: lmxa(:),konft(:,:,:),ipc(:),ipcx(:)
c      if(.not.master_mpi) return
      end subroutine
      end module m_lmfgw
