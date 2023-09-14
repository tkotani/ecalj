module m_mkqp ! Set up k-points and related quantities for BZ integration
  integer,allocatable,protected :: iv_a_oidtet(:), iv_a_oipq(:) 
  real(8),allocatable,protected :: rv_p_oqp (:,:), rv_a_owtkp(:)
  integer,protected:: bz_nabc(3),bz_ntet, bz_nkp
contains
  subroutine m_mkqp_init() ! Set up k-points and related quantities for BZ integration
    use m_lmfinit,only: lshft=>bz_lshft, bz_tetrahedron,bz_lmet,ldos,nkk=>bz_nabcin
    use m_lattic,only:  plat=>lat_plat
    use m_mksym,only:   rv_a_osymgr, npgrp=>lat_npgrp
    use m_tetirr,only: tetirr
!    use m_bzmesh,only: bzmesh
    implicit none
    logical:: ltet,llshft(3)
    integer:: mxkp, nkp, ntet , i , iprint,ifac(3)
    integer,allocatable :: iv_a_owk(:),iv_a_tmp(:)
    real(8),allocatable:: rv_a_owtkp_temp(:),rv_p_oqp_temp(:,:)
    character outs*80
    real(8):: qb(3,3)
    call tcn('m_mkqp_init')
    ltet   = (bz_lmet/=0 .or. ldos/=0) .and. bz_tetrahedron !=T means reading or generating tetrahedra corners, if tetrahedron integration set.
    mxkp   = nkk(1)*nkk(2)*nkk(3)
    llshft = [(lshft(i)/=0,i=1,3)]
    allocate(rv_a_owtkp_temp(mxkp), rv_p_oqp_temp(3,mxkp), iv_a_oipq(6*mxkp))
    call bzmesh(plat,qb,ifac,nkk(1),nkk(2),nkk(3),llshft,rv_a_osymgr, npgrp,iv_a_oipq, rv_p_oqp_temp,rv_a_owtkp_temp,nkp,mxkp)
    allocate(rv_a_owtkp, source=rv_a_owtkp_temp(1:nkp))
    allocate(rv_p_oqp,   source=rv_p_oqp_temp(1:3,1:nkp))
    ntet=0
    if(ltet) then ! ... Generate inequivalent tetrahedra
       allocate(iv_a_tmp(mxkp*30))
       call tetirr(qb, nkk(1), nkk(2), nkk(3), iv_a_oipq,ntet,iv_a_tmp)
       allocate(iv_a_oidtet, source=iv_a_tmp(1:5*ntet))
    endif
    bz_nkp  = nkp
    bz_nabc = nkk
    bz_ntet = ntet
    call tcx('m_mkqp_init')
  end subroutine m_mkqp_init
end module m_mkqp
