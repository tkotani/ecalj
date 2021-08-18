      module m_supot
      use m_struc_def,only: s_rv1

      complex(8) ,allocatable,protected ::  zv_a_obgv (:)
      integer ,   allocatable,protected ::  iv_a_oips0 (:)
      real(8) ,   allocatable,protected ::  rv_a_ogv (:,:)
      integer ,   allocatable,protected ::  iv_a_okv (:,:)
      real(8),protected:: lat_gmax
      integer,protected:: lat_ng,lat_nabc(3), k1 , k2 , k3 
      contains
      
      subroutine m_supot_init()  
      use m_lattic,only: rv_a_odlv,rv_a_oqlv,lat_plat,rv_a_opos
      use m_mksym,only:   rv_a_osymgr,rv_a_oag 
      use m_lmfinit,only : lcd4,ctrl_nbas,ctrl_nspin,lat_alat,ftmesh,lat_gmaxin
      use m_lattic,only: lat_vol, lat_awald
      use m_lattic,only: lat_nkd, lat_nkq
      use m_mksym,only:  lat_nsgrp
C- Initialization for G vectors bgv,ips0,gv,kv !See gvlst2 and sgvsym 
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode=0 for FP clculaiton (=1 make Madelung matrix for monopoles (ASA))
C ----------------------------------------------------------------------
      implicit none
      integer,parameter:: mode=0
      integer:: nbas ,  nsp , nkd , nkq , igets , ngabc(3) 
     .        , n1 , n2 , n3 , ngmx , ng , ngrp , iprint
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision awald,alat,vol,plat(3,3),gmax,xx
      integer ::iwdummy
      real(8):: wdummy(3)=0d0
      call info0(50,1,0,' supot : potential setup ... ')
      nbas=ctrl_nbas
      nsp=ctrl_nspin
      alat=lat_alat
      vol=lat_vol
      awald=lat_awald
      nkd=lat_nkd
      nkq=lat_nkq
C --- Setup for FT charge density, potential representation ---
      gmax = lat_gmaxin
      ngabc= ftmesh    
      if (lcd4) then
        alat = lat_alat
        plat = lat_plat
C   ... Generate energy cutoff gmax or n1..n3
        call mshsiz(alat,plat,0,gmax,ngabc,ngmx)
        call fftz30(n1,n2,n3,k1,k2,k3)
C   ... Make list of lattice vectors within cutoff
        allocate(rv_a_ogv(ngmx,3))
        allocate(iv_a_okv(ngmx,3))
        call gvlst2(alat, plat, wdummy, n1,n2,n3, 0d0,gmax,0,8, ngmx, ng, iv_a_okv, rv_a_ogv, xx, xx)
        if (ng /= ngmx) then
          print *,' gmax,ng ngmx=',gmax,ng,ngmx
          call rx('supot: bug in gvlst2')
        endif
        lat_ng = ng
        lat_gmax = gmax
        lat_nabc = ngabc
        k1=lat_nabc(1)
        k2=lat_nabc(2)
        k3=lat_nabc(3)
        allocate(iv_a_oips0(ng))
        allocate(zv_a_obgv(ng))
        iv_a_oips0(:)=0.0d0
        zv_a_obgv(:)=0.0d0
        ngrp = lat_nsgrp
        call sgvsym ( ngrp , rv_a_osymgr , rv_a_oag , ng , rv_a_ogv , iv_a_oips0 , zv_a_obgv )
      endif
      end subroutine m_supot_init
      end module m_supot
