      module m_supot
      use m_struc_def,only: s_rv1
      type(s_rv1), allocatable :: sv_a_oorhat(:)
      complex(8) , allocatable,protected ::  zv_a_obgv (:)
      integer , allocatable,protected ::  iv_a_oips0 (:)
      real(8) , allocatable,protected ::  rv_a_ogv (:)
      integer , allocatable,protected ::  iv_a_okv (:)
      complex(8) , allocatable ::  zv_a_osmrho (:)
      complex(8) , allocatable ::  zv_a_osmpot (:)
      real(8) , allocatable ::  rv_a_omad (:)
      contains
      subroutine supot(mode)    !,spot) !,sctrl,slat
      use m_lattic,only: rv_a_odlv,rv_a_oqlv,lat_plat
      use m_mksym,only:    rv_a_osymgr,rv_a_oag 
      use m_lmfinit,only : lcd4,ctrl_nbas,ctrl_nbasp,ctrl_nspin
      use m_lmfinit,only: lat_alat,lat_vol,lat_awald,lat_nkd,lat_gmax,lat_ng,
     &     lat_nabc,lat_nkq,rv_a_opos
      use m_mksym,only:  lat_nsgrp
C- Initialization for potential
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1 make Madelung matrix for monopoles (ASA)
Ci   sctrl :struct containing parameters governing program flow
Ci   slat  :struct containing information about the lattice
Ci     Stored:    ng ogv okv gmax nabc oips0 obgv
Ci   spot  :struct containing information about the potential
Ci     Stored:    omad osmrho osmpot oorhat
Co Outputs
Cr Remarks
Cr   Permanent arrays are created to handle the potential, depending
Cr   on the mode
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode
c      type(s_ctrl)::sctrl
c      type(s_lat)::slat
c      type(s_pot)::spot
C ... Local parameters
      integer:: nbas , nbasp , nbaspp , nsp , nkd , nkq , igets , ngabc(3) 
     .             , n1 , n2 , n3 , k1 , k2 , k3 , ngmx , ng , ngrp , iprint
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision awald,alat,vol,plat(3,3),gmax,xx
      integer ::iwdummy
      real(8):: wdummy(3)=0d0
      call info0(50,1,0,' supot : potential setup ... ')
      nbas=ctrl_nbas
      nbasp=ctrl_nbasp
      nsp=ctrl_nspin
      alat=lat_alat
      vol=lat_vol
      awald=lat_awald
      nkd=lat_nkd
      nkq=lat_nkq
      nbaspp = 2*nbasp - nbas

C --- Monopole Madelung matrix ---
      if (mod(mode,10) .eq. 1) then
        allocate(rv_a_omad(nbaspp**2))
        rv_a_omad(:)=0.0d0
        call madmat ( nbaspp , rv_a_opos , awald , alat , vol , rv_a_odlv
     .      , nkd , rv_a_oqlv , nkq , rv_a_omad )
      endif
C --- Setup for FT charge density, potential representation ---
c      if ( iand(4,int(ctrl_lcd)) .ne. 0 ) then
      if (lcd4) then
        alat = lat_alat
        plat = lat_plat
        ngabc= lat_nabc
        gmax = lat_gmax
C   ... Generate energy cutoff gmax or n1..n3
c        call pshpr(iprint()-0)
       call mshsiz(alat,plat,0,gmax,ngabc,ngmx)
c       print *,' xxxx mshsiz=',ngmx
        call fftz30(n1,n2,n3,k1,k2,k3)
C   ... Make list of lattice vectors within cutoff
        allocate(rv_a_ogv(abs(ngmx*3)))
        if (allocated(iv_a_okv)) deallocate(iv_a_okv)
        allocate(iv_a_okv(abs(ngmx*3)))
c        call gvlist ( alat , plat , iwdummy , n1 , n2 , n3 , gmax , 8, ngmx , ng , lat_iv_a_okv , lat_rv_a_ogv , xx , xx )
        call gvlst2(alat, plat, wdummy, n1,n2,n3, 0d0,gmax,0,8, ngmx, ng, iv_a_okv, rv_a_ogv, xx, xx)
        if (ng .ne. ngmx) then
          print *,' gmax,ng ngmx=',gmax,ng,ngmx
          call rx('supot: bug in gvlst2')
        endif
        lat_ng = ng
        lat_gmax = gmax
        lat_nabc = ngabc !it was call ii8copy(i_copy_size=3,ngabc,1,lat_nabc,1)
        if (allocated(zv_a_osmrho)) deallocate(zv_a_osmrho)
        allocate(zv_a_osmrho(k1*k2*k3*nsp))
        zv_a_osmrho(:)=0d0
        if (allocated(zv_a_osmpot)) deallocate(zv_a_osmpot)
        allocate(zv_a_osmpot(k1*k2*k3*nsp))
        zv_a_osmpot(:)=0.0d0
        if (allocated(sv_a_oorhat)) deallocate(sv_a_oorhat)
        allocate(sv_a_oorhat(3*nbaspp))
c        do i=1,3*nbaspp
c        enddo
C   ... Setup for symmetrization
        allocate(iv_a_oips0(ng))
        iv_a_oips0(:)=0.0d0
        allocate(zv_a_obgv(ng))
        zv_a_obgv(:)=0.0d0
        ngrp=lat_nsgrp
        call sgvsym ( ngrp , rv_a_osymgr , rv_a_oag , ng , rv_a_ogv , iv_a_oips0
     .       , zv_a_obgv )
      endif
      end subroutine supot
      end module m_supot
