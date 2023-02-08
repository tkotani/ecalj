module m_supot
  use m_struc_def,only: s_rv1

  complex(8) ,allocatable,protected ::  zv_a_obgv (:)
  integer,   allocatable,protected ::  iv_a_oips0 (:)
  real(8),   allocatable,protected ::  rv_a_ogv (:,:)
  integer,   allocatable,protected ::  iv_a_okv (:,:)
  real(8),protected:: lat_gmax
  integer,protected:: lat_nabc(3),lat_ng,k1 , k2 , k3
  integer,protected,target,private::  ngabc(3)
  integer,protected,pointer:: n1,n2,n3
contains

  subroutine m_supot_init()
    use m_lattic,only: rv_a_odlv,rv_a_oqlv,lat_plat,rv_a_opos,lat_qlat
    use m_mksym,only:   rv_a_osymgr,rv_a_oag
    use m_lmfinit,only : lcd4,nsp,lat_alat,ftmesh,lat_gmaxin,stdo
    use m_lattic,only: lat_vol, lat_awald
    use m_lattic,only: lat_nkd, lat_nkq
    use m_mksym,only:  lat_nsgrp
    use m_shortn3,only: shortn3_initialize,shortn3,nout,nlatout
    use m_ftox
    !- Initialization for G vectors bgv,ips0,gv,kv !See gvlst2 and sgvsym
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   mode=0 for FP clculaiton (=1 make Madelung matrix for monopoles (ASA))
    ! ----------------------------------------------------------------------
    implicit none
    integer,parameter:: mode=0
    integer:: nkd , nkq , ngmx , ng , ngrp , iprint
    double precision :: awald,alat,vol,plat(3,3),gmax,xx
    real(8):: qpg(3),gg,qlat1(3,3),gmax2,gs(3),qlat(3,3),tpiba,  tol=1d-8
    integer:: ig,j1,j2,j3
    call tcn('m_supot_init')
    n1=>ngabc(1)
    n2=>ngabc(2)
    n3=>ngabc(3)
    ngabc=ftmesh !initial condition for mshsiz
!    nbas=ctrl_nbas
!    nsp=ctrl_nspin
    alat=lat_alat
    vol=lat_vol
    awald=lat_awald
    nkd=lat_nkd
    nkq=lat_nkq
    ! --- Setup for FFT charge density, potential representation ---
    gmax = lat_gmaxin
    if (lcd4) then
       alat = lat_alat
       plat = lat_plat
       qlat = lat_qlat
       call mshsiz(alat,plat,gmax,ngabc,ngmx) !return n1 n2 n3 (=ngabc) satisfying gmax
       !write(stdo,ftox)' 000 gmax ngmx=',ftof(gmax),ngmx
       call fftz30(n1,n2,n3,k1,k2,k3)
       !   ... Make list of lattice vectors within cutoff
       gvblock: block
         real(8):: ogv(ngmx,3)
         integer:: okv(ngmx,3)
         call gvlst2(alat, plat, [0d0,0d0,0d0], n1,n2,n3, 0d0,gmax,0,8+1000, ngmx, ng, okv, ogv, xx)  !+1000 for symmetry cheker for sgvsym
         !write(stdo,ftox)' supot: gmax ng ngmx=',ftof(gmax),ng,ngmx
         ngmx = ng
         allocate(rv_a_ogv(ng,3))
         allocate(iv_a_okv(ng,3))
         rv_a_ogv(1:ng,1:3)=ogv(1:ng,1:3)
         iv_a_okv(1:ng,1:3)=okv(1:ng,1:3)
         !print *,'ogv(1:ng,1:3)',ogv(1,1:3),okv(1,1:3)
       endblock gvblock
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
    call tcx('m_supot_init')
  end subroutine m_supot_init
end module m_supot
