module m_mkqp ! Set up k-points and related quantities for BZ integration
  integer,allocatable,protected :: iv_a_oidtet(:), iv_a_oipq(:) 
  real(8),allocatable,protected :: rv_p_oqp (:,:), rv_a_owtkp(:)
  integer,protected:: bz_nabc(3),bz_ntet, bz_nkp
contains
  subroutine m_mkqp_init() ! Set up k-points and related quantities for BZ integration
    use m_lmfinit,only: lshft=>bz_lshft, bz_tetrahedron,bz_lmet,ldos,nkxyz=>bz_nabcin
    use m_lattic,only:  plat=>lat_plat!,qlat=>lat_qlat
    use m_mksym,only:   rv_a_osymgr, npgrp=>lat_npgrp
    use m_tetirr,only: tetirr
    implicit none
    logical:: ltet,llshft(3)
    integer:: mxkp, nkp, ntet , i , iprint
    integer,allocatable :: iv_a_owk(:),iv_a_tmp(:)
    character outs*80
    real(8):: qlat(3,3),vol
    call tcn('m_mkqp_init')
    ltet   = (bz_lmet/=0 .or. ldos/=0) .and. bz_tetrahedron !=T means reading or generating tetrahedra corners, if tetrahedron integration set.
    mxkp   = nkxyz(1)*nkxyz(2)*nkxyz(3)
    llshft = [(lshft(i)/=0,i=1,3)]
    allocate(rv_a_owtkp(mxkp), rv_p_oqp(3,mxkp), iv_a_oipq(6*mxkp))
    call dinv33(plat,1,qlat,vol)
    call bzmesh(plat,qlat,nkxyz(1),nkxyz(2),nkxyz(3),llshft,rv_a_osymgr, npgrp,iv_a_oipq,rv_p_oqp,rv_a_owtkp,nkp,mxkp)
    deallocate(rv_a_owtkp,      rv_p_oqp)
    allocate  (rv_a_owtkp(nkp), rv_p_oqp(3,nkp))
    call pshpr(0)
    call bzmesh(plat,qlat,nkxyz(1),nkxyz(2),nkxyz(3),llshft,rv_a_osymgr, npgrp,iv_a_oipq,rv_p_oqp,rv_a_owtkp,nkp,mxkp)
    call poppr
    ntet=0
    if(ltet) then ! ... Generate inequivalent tetrahedra
       allocate(iv_a_tmp(mxkp*30))
       call tetirr(qlat, nkxyz(1), nkxyz(2), nkxyz(3), iv_a_oipq,ntet,iv_a_tmp)
       allocate(iv_a_oidtet(5*ntet))
       iv_a_oidtet(:)=iv_a_tmp(1:5*ntet)
       deallocate(iv_a_tmp)
    endif
    bz_nkp  = nkp
    bz_nabc = nkxyz
    bz_ntet = ntet
    call tcx('m_mkqp_init')
  end subroutine m_mkqp_init
end module m_mkqp
