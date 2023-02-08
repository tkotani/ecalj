module m_mkqp
  integer,allocatable,protected :: iv_a_oidtet(:), iv_a_oipq(:), iv_a_ostar(:)
  real(8),allocatable,protected :: rv_p_oqp (:,:), rv_a_owtkp(:)
  integer,protected:: bz_nabc(3),bz_ntet, bz_nkp
contains
  subroutine m_mkqp_init()
    use m_lmfinit,only: bz_lshft,bz_tetrahedron,bz_lmet,ldos,bz_nabcin
    use m_lattic,only:  lat_plat
    use m_mksym,only:   rv_a_osymgr,lat_npgrp
    use m_tetirr,only: tetirr
    !! Set up k-points and related quantities for BZ integration
    !! ----------------------------------------------------------------------
    !! gettet: = bz_lmet/=0 .or. ldos/=0 given at setup
    !!            T read or generate tetrahedra corners, if
    !!             tetrahedron integration set
    !!   lnoirr: =F (T suppress generation of inequivalent tetrahedra)
    !!   lreduc: xxx 0 do not save array ipq
    !!         : =1 save array ipq
    !!         : xxx -1 ignore symmetry operations, make qp for full BZ.
    !!   lgstar: nozero, generate igstar according to bzmesh, which see
    !!         : xxx 0 igstar is not made
    !!         : xxx 2 igstar contains inverse mapping of ipq
    !!         : =-2 igstar contains group ops rotating irreducible
    !!         :    to to full BZ.
    !! Outputs protected
    !! idtet, iqp, qp, star, wtkp : for the Brillouin Zone
    implicit none
    logical :: gettet
    integer::  lgstar=-2, i_copy_size,i_data_size !lreduc=1,
    logical:: lgors,ltet,lnoirr=.false.,llshft(3),lipq !lsx,
    integer:: mxkp , nfilqp , nkp , nkxyz(3) , npgrp &
         , lshft(3) , lpbc , ntet , i , iprint , igets
    integer,allocatable :: iv_a_owk(:)
    integer,allocatable :: iv_a_tmp(:)
    double precision :: plat(3,3),qlat(3,3),vol
    character outs*80
    integer ::iwdummy
    call tcn('m_mkqp_init')
    gettet = bz_lmet/=0 .or. ldos/=0
    ntet = 0
    nkxyz=bz_nabcin
    lshft=bz_lshft
    plat =lat_plat
!    nsgrp=lat_nsgrp
    npgrp=lat_npgrp
    lpbc = 0
    ltet = gettet .and. bz_tetrahedron
    call dinv33(plat,1,qlat,vol)
    ! ... Make the qp list from bzmesh
    mxkp = nkxyz(1)*nkxyz(2)*nkxyz(3)
    allocate(iv_a_ostar(0:mxkp))
    allocate(rv_a_owtkp(mxkp))
    allocate(rv_p_oqp(3,mxkp))
    allocate(iv_a_oipq(6*mxkp))
    iv_a_ostar=0
    rv_a_owtkp=0d0
    llshft = [(lshft(i)/=0,i=1,3)]
    iv_a_ostar(0)=lgstar
    call bzmesh ( plat , qlat , nkxyz ( 1 ) , nkxyz ( 2 ) , nkxyz &
         ( 3 ) , llshft , rv_a_osymgr, npgrp , iv_a_oipq , rv_p_oqp , &
         rv_a_owtkp , nkp , mxkp , iv_a_ostar , lpbc )
    deallocate(rv_a_owtkp,      rv_p_oqp)
    allocate  (rv_a_owtkp(nkp), rv_p_oqp(3,nkp))
    rv_a_owtkp=0d0
    call pshpr(0)
    iv_a_ostar(0)=lgstar
    call bzmesh ( plat , qlat , nkxyz ( 1 ) , nkxyz ( 2 ) , nkxyz &
         ( 3 ) , llshft , rv_a_osymgr, npgrp , iv_a_oipq , rv_p_oqp , &
         rv_a_owtkp , nkp , mxkp , iv_a_ostar , lpbc )
    call poppr
    if (ltet) then ! ... Generate inequivalent tetrahedra
       allocate(iv_a_tmp(mxkp*30))
       call tetirr(qlat, nkxyz(1), nkxyz(2), nkxyz(3), iv_a_oipq,ntet,iv_a_tmp)
       allocate(iv_a_oidtet(5*ntet))
       iv_a_oidtet(:)=iv_a_tmp(1:5*ntet)
       deallocate(iv_a_tmp)
    endif
    bz_nkp   = nkp
    bz_nabc = nkxyz
    bz_ntet  = ntet
    call tcx('m_mkqp_init')
  end subroutine m_mkqp_init
end module m_mkqp

