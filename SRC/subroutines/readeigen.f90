!> Return eigenvalus and eigenfunctions for given q and isp.
! -----------------------------------------------------------
! We can get eigenfunctions for Wannier, as well. See hmagnon.F
! note: we have to call init_foobar to call readeval, readcphi, readgeig.
! ----------------
module m_readeigen
  use m_ftox
  use m_mpiio,only: openm,readm,closem
  use m_lgunit,only:stdo
  use m_iqindx_qtt,only: Iqindx2_, Init_iqindx_qtt
  use m_hamindex,only:   ngpmx, nqtt, nqi, qtt,iqimap, iqmap,igmap,shtvg,qlat,symops,ngrp
  use m_hamindex,only:   plat,invgx, miat,tiat,dlmm,shtvg,symops,lmxax,nbas
  use m_read_bzdata,only: ginv
  use m_struct_from_lmf,only: nsp=>nspin, mrecb,mrece,mrecg,nband,nspc; use m_gw_product_basis,only: ndima,ndimanspc,nspx
  use m_keyvalue,only: getkeyvalue
  use m_GWinput, only: gwinput_init, gwinput_loaded, tg_KeepQG => KeepQG
  use m_keep_wfs,only: keep_wfs_init, update_keep_geig, update_keep_cphi, set_geig_from_keep, set_cphi_from_keep
  use m_sharedmem, only: shm_init, shm_rank, shm_alloc_c8_4d, shm_barrier
  use m_mpi, only: ipr, MPI__AllreduceAND, MPI__zBcast => MPI__zBcast_h, get_mpi_master
#ifdef __GPU
  use m_blas, only : zmm => zmm_d
#else
  use m_blas, only : zmm => zmm_h
#endif
  use,intrinsic :: ieee_arithmetic
  !! qtt(1:3, nqtt)  :q-vector in full BZ (no symmetry) in QGpsi, QGcou
  !! qtti(1:3,nqi)   :eivenvalues, eigenvectors are calculated only for irr=1 in QGpsi (See lqg4gw).
  implicit none
  public:: Init_readeigen,Init_readeigen2, Lowesteval, Readeval,Readgeigf,Readcphif
  public:: readgeigf_mpi, readcphif_mpi
  public:: Readcphifq
  integer, allocatable, public, protected :: ngp(:)
  private
  integer:: norbtx,imx,ifcphim,ifgeigm,nqixx  !ifgeigW,ifcphiW,
  real(8),private:: leval, quu(3)
  logical,private:: init=.true.,init2=.true.,keepeig
  logical,private:: debug=.false.
  character(8),external :: xt
  real(8),allocatable,private:: evud(:,:,:)
  ! geig/cphi: node-shared via MPI shared memory on both CPU and GPU builds.
  ! GPU build transfers one q-point slice to device in readgeigf_mpi/readcphif_mpi
  ! rather than keeping the full arrays in VRAM.
  complex(8),pointer,private:: geig(:,:,:,:)=>null(),cphi(:,:,:,:)=>null()
  integer,save,private:: geig_shm_id=-1, cphi_shm_id=-1
  integer,allocatable,private:: ngvecp(:,:,:), ngvecprev(:,:,:,:)
  integer,allocatable,private:: l_tbl(:),k_tbl(:),ibas_tbl(:),offset_tbl(:),offset_rev_tbl(:,:,:)
  logical,private:: keepqg
contains
  function readcphifq() result(qu)
    real(8):: qu(3)                  ! I think qu=q now.
    qu=quu
  end function readcphifq
  pure function readeval(q,isp) result(ev) ! Return ev(1:nband) for given q(1:3) and isp
    intent(in)  ::       q,isp
    integer :: isp
    real(8) :: q(3), ev(nband)
    integer:: iq,iqindx,i
    real(8):: qu(3)
    call iqindx2_(q, iq, qu) !qu is used q. q-qu is a G vector.
    ev(1:nband) = evud(1:nband,iqimap(iq),isp) !iqimap is given in suham.F/gen_hamindex
    !if(debug) then
    !   if(ipr) write(6,*)'iq iqimap(iq)=',iq,iqimap(iq)
    !   if(ipr) write(6,"('iq iqimap(iq) q=',2i8,3f13.5)")iq,iqimap(iq),q
    !   if(ipr) write(6,"(9f9.4)")ev(1:9)
    !endif
  end function readeval
  !> Return ev(1:nband) for given q(1:3) and isp
  function readgeigf(q,isp) result(geigen)
    real(8),intent(in):: q(3)
    integer,intent(in):: isp
    real(8):: qu(3)
    complex(8):: geigen(ngpmx*nspc,nband)
    geigen=0d0 !2024-5-17 for ifort NaN initialization
    call readgeig(q,isp,qu,geigen)
  end function readgeigf
  subroutine readgeig(q,isp, qu,geigen)
    use m_ftox
    use m_rotwave,only: Rotipw
    use m_ftox
    implicit none
    !logical,optional,intent(in):: fpmt
    real(8),intent(in) :: q(3)
    integer,intent(in) :: isp
    real(8),intent(out) :: qu(3)
    complex(8), intent(out) :: geigen(ngpmx*nspc,nband)
    integer:: iq,iqindx,ikpisp,napw,iqq,nnn(3),ig,igg,ig2,iqi,igxt,i,ioff,ispc
    real(8)   :: ddd(3),platt(3,3),qpg(3),qpgr(3),qtarget(3),qout(3),qin(3)
    complex(8):: geigenr(ngpmx*nspc,nband),img=(0d0,1d0),img2pi
    integer :: ifiqg
    img2pi=2d0*4d0*datan(1d0)*img
    platt=transpose(plat) !this is inverse of qlat
    if(init2) call rx( 'readgeig: modele is not initialized yet')
    call iqindx2_(q, iq, qu) !qu is used q. q-qu is a G vector.
    if(debug.and.ipr) write(6,*)' readgeig:xxx iq=',iq
    iqq=iqmap(iq)
    iqi=iqimap(iq)
    igg=igmap(iq)
    qtarget=qtt(:,iq) ! iqq is mapped to qtarget=qu=qtt(:,iq)
    !!  qtt(iqq) is rotated to qtt(iq) by sympos(  ,igg).
    if(ngp(iq)==0) return
    if(ngp(iq)/=ngp(iqq)) then
       if(ipr) write(6,*)' ddddd readgeig: iq iqq igg=',iq,iqq,igg,q,qu
       if(ipr) write(6,*)' ddddd qtarget=',qtarget,' ddddd q (iqq)=',qtt(:,iqq)
       if(ipr) write(6,"(a,3i5,3f10.4,2i5)")' ngp(iq) ngp(iqq)=',iq,iqq,igg,q,ngp(iq),ngp(iqq)
       call rx( 'readgeig:x2 ngp(iq)/=ngp(iqq)')
    endif
    if(keepeig) then
       geigenr(1:ngpmx*nspc,1:nband) = geig(1:ngpmx*nspc,1:nband,iqi,isp)
    else
       ikpisp= isp + nsp*(iqi-1)
       i=readm(ifgeigm,rec=ikpisp, data=geigenr(1:ngpmx*nspc,1:nband) )
!       open(newunit=ifgeig, file='GEIG'//trim(xt(iqi))//trim(xt(isp)),form='unformatted')
!       read(ifgeig) geigenr(1:ngpmx*nspc,1:nband) 
!       close(ifgeig)
    endif
    !!   qinput: qtt(:,iqq)  ---> qtarget: qtt(:,iq) ( G-vector difference from symops*qtt(:,iqq) )
    igxt=1 !not timereversal
    if(.not.keepqg) then
      allocate( ngvecp(3,ngpmx,iqq:iqq))
      allocate( ngvecprev(-imx:imx,-imx:imx,-imx:imx,iq:iq))
      BLOCK
        integer:: ngvecp_tmp(3,ngpmx)
        open(newunit=ifiqg, file='__QGpsi_rec',form='unformatted', access='direct', recl=4*(3*ngpmx+(2*imx+1)**3), status='old')
        read(ifiqg, rec=iq)  ngvecp_tmp(1:3,1:ngpmx),ngvecprev(-imx:imx,-imx:imx,-imx:imx,iq)
        read(ifiqg, rec=iqq) ngvecp(1:3, 1:ngp(iqq),iqq)
        close(ifiqg)
      END BLOCK
    endif
    do ispc=1,nspc
       ioff=(ispc-1)*ngpmx
       call rotipw(qtt(:,iqq),qtt(:,iq),ngp(iqq),nband, &
            platt,qlat,symops(:,:,igg),ngvecp(:,:,iqq),ngvecprev(:,:,:,iq),shtvg(:,igg),igxt,imx, &
            geigenr(ioff+1:ioff+ngp(iqq),1:nband), geigen(ioff+1:ioff+ngp(iq),1:nband))
    enddo
    if(.not.keepqg) deallocate(ngvecp,ngvecprev)
  end subroutine readgeig
  function readcphif(q,isp) result(cphif)
    integer,intent(in):: isp
    real(8),intent(in):: q(3)
    real(8) :: qu(3)
    complex(8):: cphif(ndima*nspc,nband)
    call readcphi(q,isp, qu, cphif)
    quu=qu
  end function readcphif
  subroutine readcphi(q,isp,  qu,cphif)!, fpmt)
    use m_rotwave,only: Rotmto
    implicit none
    !logical,optional,intent(in):: fpmt
    !!-- return mto part of eigenfunction for given q(1:3) and isp
    real(8), intent(in) :: q(3)
    integer, intent(in)  :: isp
    real(8), intent(out)  :: qu(3)
    complex(8), intent(out)  :: cphif(ndima*nspc,nband)
    integer:: iq,iqindx,ikpisp,iqq,iorb,ibaso,ibas,k,l,ini1,ini2,iend1,iend2, igg,ig,iqi,i,igxt,ioff,ispc,ix
    real(8)   :: qrot(3) ,qout(3)
    complex(8):: phase,cphifr(ndima*nspc,nband),phaseatom !takao 1->*->nband
    complex(8),parameter:: img=(0d0,1d0) ! MIZUHO-IR
    complex(8):: img2pi = 2d0*4d0*datan(1d0)*img ! MIZUHO-IR
    if(init2) call rx( 'readcphi: modele is not initialized yet')
    call iqindx2_(q, iq, qu) !for given q, get iq. qu is used q. q-qu= G vectors. qu=qtt(:,iq)
    igg=igmap(iq)  ! qtt(:,iq)= matmul(sympos(  ,igg),qtt(:,iqq))
    iqq=iqmap(iq)  ! mapped from qtt(:,iqq) to qtt(:,iq);
    ! qtt(:,iq)=matmul(sym(igg),qtt(:,iqq))+some Gvector(see iqindx2 above)
    iqi=iqimap(iq) ! iqi is index for irr.=1 (cphi calculated. See qg4gw and sugw.F).
    ! qtt(:,iqq) = qtti(:,iqi) is satisfied.
    ! we have eigenfunctions calculated only for qtti(:,iqi).
    if(keepeig) then
       cphifr(1:ndima*nspc,1:nband) = cphi(1:ndima*nspc,1:nband,iqi,isp)
    else 
       ikpisp= isp + nsp*(iqi-1)
       i=readm(ifcphim,rec=ikpisp, data=cphifr(1:ndima*nspc,1:nband)) ! , rec=ikpisp
!     open(newunit=ifcphi, file='CPHI'//trim(xt(iqi))//trim(xt(isp)),form='unformatted')
!       read(ifcphi) cphifr(1:ndima*nspc,1:nband) 
!       close(ifcphi)
    endif
    if(debug) write(6,"('readcphi:: xxx sum of cphifr=',3i4,4d23.16)")ndimanspc,ndimanspc,norbtx, &
         sum(cphifr(1:ndimanspc,1:nband)),sum(abs(cphifr(1:ndimanspc,1:nband)))
    igxt=1 !not timereversal (for future)
    do ispc=1,nspc
       ioff=ndima*(ispc-1)
       call rotmto(qtt(:,iqq),ndima,nband,norbtx,ibas_tbl,l_tbl,k_tbl,offset_tbl,offset_rev_tbl, &
            maxval(ibas_tbl),maxval(l_tbl),maxval(k_tbl), &
            symops(1,1,igg),shtvg(:,igg),dlmm(:,:,:,igg),lmxax,miat(:,igg),tiat(:,:,igg),igxt,nbas, &
            cphifr(ioff+1:ioff+ndima,:), cphif(ioff+1:ioff+ndima,:))
     enddo
  end subroutine readcphi

  function readgeigf_mpi(q, isp, mpi_mode, comm) result(geigen)
    real(8), intent(in) :: q(3)
    integer, intent(in) :: isp
    logical, intent(in), optional :: mpi_mode
    integer, intent(in), optional :: comm
    complex(8) :: geigen(ngpmx*nspc,nband), geigenr(ngpmx*nspc,nband)
    integer :: i, iq, ikpisp, iqq, igg, iqi, igxt, ioff, ispc, ifiqg
    real(8) :: platt(3,3), qtarget(3), qu(3)
    logical :: has_geig, mpi_master, mpi_mode_in
    integer, save :: iqq_prev = -99999, iq_prev = -99999
#ifdef __GPU
    attributes(device) :: geigen
#endif
    mpi_master = .true.
    mpi_mode_in = .false.
    if(present(mpi_mode)) mpi_mode_in = mpi_mode
    if(mpi_mode_in .and. .not. present(comm)) call rx( 'readgeigf_mpi: mpi_mode requires comm')
    if(mpi_mode_in) mpi_master = get_mpi_master(comm)
    !$acc kernels
    geigen(:,:) = (0d0, 0d0) !2024-5-17 for ifort NaN initialization
    !$acc end kernels
    platt = transpose(plat) !this is inverse of qlat
    if(init2) call rx( 'readgeig: modele is not initialized yet')
    call iqindx2_(q, iq, qu) !qu is used q. q-qu is a G vector.
    quu = qu
    if(debug) write(stdo,*)' readgeig:xxx iq=',iq
    iqq=iqmap(iq)
    iqi=iqimap(iq)
    igg=igmap(iq)
    qtarget=qtt(:,iq) ! iqq is mapped to qtarget=qu=qtt(:,iq)
    if(ngp(iq)==0) return
    if(ngp(iq)/=ngp(iqq)) then
       if(ipr) write(stdo,*)' ddddd readgeig: iq iqq igg=',iq,iqq,igg,q,qu
       if(ipr) write(stdo,*)' ddddd qtarget=',qtarget,' ddddd q (iqq)=',qtt(:,iqq)
       if(ipr) write(stdo,"(a,3i5,3f10.4,2i5)")' ngp(iq) ngp(iqq)=',iq,iqq,igg,q,ngp(iq),ngp(iqq)
       call rx( 'readgeig:x2 ngp(iq)/=ngp(iqq)')
    endif
    !$acc enter data create(geigenr)
    if(keepeig) then
      geigenr(1:ngpmx*nspc,1:nband) = geig(1:ngpmx*nspc,1:nband,iqi,isp)  ! host SHM slice
      !$acc update device(geigenr)
    else
      !$acc host_data use_device(geigenr)
      has_geig = set_geig_from_keep(iqi,isp,geigenr) !set geigenr if it is stored
      !$acc end host_data
      if(mpi_mode_in) call MPI__AllreduceAND(has_geig, communicator=comm)
      if(.not.has_geig) then
        ikpisp= isp + nsp*(iqi-1)
        if(mpi_master) i=readm(ifgeigm,rec=ikpisp, data=geigenr(1:ngpmx*nspc,1:nband))
        if(mpi_mode_in) call MPI__zBcast(geigenr, ngpmx*nspc*nband, communicator=comm)
        !$acc update device(geigenr)
        !$acc host_data use_device(geigenr)
        call update_keep_geig(iqi,isp,geigenr)
        !$acc end host_data
      endif
    endif
    igxt=1 !not timereversal
    if(.not.keepqg .and. (iqq_prev /= iqq .or. iq_prev /= iq)) then
      ReadQGpsi: BLOCK
        integer:: ngvecp_tmp(3,ngpmx)
        open(newunit=ifiqg, file='__QGpsi_rec',form='unformatted', access='direct', recl=4*(3*ngpmx+(2*imx+1)**3), status='old')
        if(iqq_prev /= iqq)  then
          if(allocated(ngvecp)) then
            !$acc exit data delete(ngvecp)
            deallocate(ngvecp)
          endif
          allocate(ngvecp(3,ngpmx,iqq:iqq))
          read(ifiqg, rec=iqq) ngvecp(1:3, 1:ngp(iqq),iqq)
          !$acc enter data copyin(ngvecp)
          iqq_prev = iqq
        endif
        if(iq_prev /= iq)  then
          if(allocated(ngvecprev)) then
            !$acc exit data delete(ngvecprev)
            deallocate(ngvecprev)
          endif
          allocate(ngvecprev(-imx:imx,-imx:imx,-imx:imx,iq:iq))
          read(ifiqg, rec=iq) ngvecp_tmp(1:3,1:ngpmx),ngvecprev(-imx:imx,-imx:imx,-imx:imx,iq)
          !$acc enter data copyin(ngvecprev)
          iq_prev = iq
        endif
        close(ifiqg)
      END BLOCK ReadQGpsi
    endif
    do ispc=1,nspc
      ioff=(ispc-1)*ngpmx
      rotipw: block
        complex(8), parameter :: img=(0d0,1d0), img2pi = 2d0*4d0*datan(1d0)*img
        integer :: ig, ig2, iband, nnn(3)
        complex(8) :: cphase
        real(8) :: qpg(3), qpgr(3), qin(3)
        qin(:) = qtt(:,iqq)
        !$acc data copyin(shtvg(1:3,igg), qin, qtarget, qlat, symops(1:3,1:3,igg), platt) &
        !$acc      present(ngvecp, ngvecprev, geigenr)
        !$acc parallel
        !$acc loop gang independent private(nnn, qpgr, qpg)
        do ig = 1,ngp(iqq)
          !$acc loop vector
          do i = 1, 3
            qpg(i) = qin(i) + sum(qlat(i,:)*ngvecp(:,ig,iqq))
          enddo
          !$acc loop vector
          do i = 1, 3
            qpgr(i) = sum(symops(i,:,igg)*qpg(:))
          enddo
          if(igxt==-1) qpgr(:) = -qpgr(:)                               ! xxxxxxxx need to check!
          !$acc loop vector
          do i = 1, 3
            nnn(i) = nint(sum(platt(i,:)*(qpgr(:)-qtarget(:))))
          enddo
          ig2 = ngvecprev(nnn(1),nnn(2),nnn(3),iq)   ! index for G
          cphase = exp(-img2pi*sum(qpgr(:)*shtvg(:,igg)))
          !$acc loop vector
          do iband = 1, nband
            geigen(ioff+ig2,iband) = geigenr(ioff+ig,iband)*cphase
          enddo
        enddo
        !$acc end parallel
        !$acc end data
      endblock rotipw
    enddo
    !$acc exit data delete(geigenr)
     if(debug) then
       if(any(ieee_is_nan(dble(geigen)))) write(stdo,ftox) "xxx NaN in Real geig"
       if(any(ieee_is_nan(imag(geigen)))) write(stdo,ftox) "xxx NaN in Imag "
     endif
  end function readgeigf_mpi

  function readcphif_mpi(q, isp, mpi_mode, comm) result(cphif)
    real(8), intent(in) :: q(3)
    integer, intent(in) :: isp
    logical, intent(in), optional :: mpi_mode
    integer, intent(in), optional :: comm
    complex(8) :: cphif(ndima*nspc,nband), cphifr(ndima*nspc,nband)
    integer:: i, iq, ikpisp, iqq, igg, iqi, igxt, ioff, ispc
    real(8) ::  qu(3)
    logical :: has_cphi, mpi_master, mpi_mode_in
#ifdef __GPU
    attributes(device) :: cphif
#endif
    mpi_master = .true.
    mpi_mode_in = .false.
    if(present(mpi_mode)) mpi_mode_in = mpi_mode
    if(mpi_mode_in .and. .not. present(comm)) call rx( 'readcphif_mpi: mpi_mode requires comm')
    if(mpi_mode_in) mpi_master = get_mpi_master(comm)
    if(init2) call rx( 'readcphi: modele is not initialized yet')
    call iqindx2_(q, iq, qu) !for given q, get iq. qu is used q. q-qu= G vectors. qu=qtt(:,iq)
    igg=igmap(iq)  ! qtt(:,iq)= matmul(sympos(  ,igg),qtt(:,iqq))
    iqq=iqmap(iq)  ! mapped from qtt(:,iqq) to qtt(:,iq);
    ! qtt(:,iq)=matmul(sym(igg),qtt(:,iqq))+some Gvector(see iqindx2 above)
    iqi=iqimap(iq) ! iqi is index for irr.=1 (cphi calculated. See qg4gw and sugw.F).
    ! qtt(:,iqq) = qtti(:,iqi) is satisfied.
    ! we have eigenfunctions calculated only for qtti(:,iqi).
    quu(:) = qu(:)
    !$acc enter data create(cphifr)
    if(keepeig) then
      cphifr(1:ndima*nspc,1:nband) = cphi(1:ndima*nspc,1:nband,iqi,isp)  ! host SHM slice
      !$acc update device(cphifr)
    else
      !$acc host_data use_device(cphifr)
      has_cphi = set_cphi_from_keep(iqi,isp,cphifr)
      !$acc end host_data
      if(mpi_mode_in) call MPI__AllreduceAND(has_cphi, communicator=comm)
      if(.not.has_cphi)then
        ikpisp= isp + nsp*(iqi-1)
        if(mpi_master) i=readm(ifcphim,rec=ikpisp, data=cphifr(1:ndima*nspc,1:nband))
        if(mpi_mode_in) call MPI__zBcast(cphifr, ndima*nspc*nband, communicator=comm)
        !$acc update device(cphifr)
        !$acc host_data use_device(cphifr)
        call update_keep_cphi(iqi,isp,cphifr)
        !$acc end host_data
      endif
    endif
    if(debug) write(stdo,"('readcphi:: xxx sum of cphifr=',3i4,4d23.16)")ndimanspc,ndimanspc,norbtx, &
         sum(cphifr(1:ndimanspc,1:nband)),sum(abs(cphifr(1:ndimanspc,1:nband)))
    call flush(6)
    igxt=1 !not timereversal (for future)
    do ispc=1,nspc
       ioff=ndima*(ispc-1)
       rotmto: block
         real(8) :: qrot(3), qin(3)
         complex(8) :: phase(nbas), dlmm_tmp(-lmxax:lmxax,-lmxax:lmxax, 0:lmxax)
         complex(8), parameter :: img=(0d0,1d0), img2pi = 2d0*4d0*datan(1d0)*img
         integer :: iorb, ibas, l, k, ini1, ini2, ierr
#ifdef __GPU
        attributes(device) :: dlmm_tmp
#endif
         dlmm_tmp(:,:,:) = dlmm(:,:,:,igg) !copy to device
         qin(:) = qtt(:,iqq)
         qrot = matmul(symops(:,:,igg),qin)
         if(igxt==-1) qrot=-qrot !july2012takao
         phase = [(exp(-img2pi*sum(qrot*tiat(:,ibas, igg))),ibas=1,nbas)]
         !$acc host_data use_device(cphifr)
         do iorb=1, norbtx
           ibas = ibas_tbl(iorb)
           l = l_tbl(iorb)
           k = k_tbl(iorb)
           ini1 = offset_tbl(iorb)+1
           ini2 = offset_rev_tbl(miat(ibas,igg),l,k)+1
           ierr = zmm(dlmm_tmp(-l,-l,l), cphifr(ini1+ioff,1), cphif(ini2+ioff,1), m=(2*l+1), n=nband, k=(2*l+1), &
                        alpha=cmplx(phase(ibas), kind=8), lda=(2*lmxax+1), ldb=ndima*nspc, ldc=ndima*nspc)
         enddo
         !$acc end host_data
       endblock rotmto
    enddo
    !$acc exit data delete(cphifr)
     if(debug) then
       if(any(ieee_is_nan(dble(cphif)))) write(stdo,ftox) "xxx NaN in Real cphi"
       if(any(ieee_is_nan(imag(cphif)))) write(stdo,ftox) "xxx NaN in Imag cphi"
     endif
    if(debug) write(stdo,*) 'end of readcphif_d'; call flush(6)
  end function readcphif_mpi

  subroutine init_readeigen() ! initialization. Save QpGpsi EVU EVD to arrays.--
    integer:: iq,is,ifiqg,nnnn,ikp,isx,ik,ib,verbose
    integer:: ifev,nband_ev, nqi_, nsp_ev ,ngpmx_ ,nqtt_,nspc_
    real(8):: QpGcut_psi
    real(8),allocatable:: qtt_(:,:),qtti_(:,:)
    if(ipr) write(stdo,ftox) 'init_readeigen:'
    if(nsp<0 .OR. nsp>2) call rx( 'init_reaeigen:nsp wrong')
    !write(*,*)'nqi=',nqi!,nqtt
    call init_iqindx_qtt()
    open(newunit=ifiqg ,file='__QGpsi',form='unformatted')
    read(ifiqg) nqtt_ , ngpmx_, QpGcut_psi, nnnn,nqi_ ,imx
    if(ipr) write(6,*)'read(ifiqg)', nqtt , ngpmx_, QpGcut_psi, nnnn,nqi
    if(nqi  /=  nqi_) call rx( 'init_readeigen:nqi/=nqi_ 11111')
    if(nqtt/=  nqtt_) call rx( 'init_readeigen:nqtt/=nqtt_ 11111')
    if(ngpmx_/=ngpmx) call rx('ngpmx error: 1111111 readeigen')
    allocate( qtt_(3,nqtt),ngp(nqtt) )
    call gwinput_init()
    if (gwinput_loaded) then
       keepqg = tg_KeepQG
    else
       call rx('m_GWinput: legacy GWinput reader is disabled. GWinput.toml is required.')
!       call getkeyvalue("GWinput","KeepQG",keepqg,default=.true.)
    endif
    if((.not.keepqg).and.ipr) write(6,*) 'keepQG = .false. in readeigen'
    if(keepqg) then
      allocate( ngvecp(3,ngpmx,nqtt))
      allocate( ngvecprev(-imx:imx,-imx:imx,-imx:imx,nqtt) )
    endif
    do ikp = 1,nqtt
       read (ifiqg) qtt_(:,ikp), ngp(ikp)
       if(keepqg) then
         read (ifiqg) ngvecp(1:3, 1:ngp(ikp),ikp),ngvecprev(-imx:imx,-imx:imx,-imx:imx,ikp)
       else
         read (ifiqg)
       endif
    enddo
    if(keepqg) then
      !$acc enter data copyin(ngvecp, ngvecprev)
    endif
    close(ifiqg)
    deallocate(qtt_)
    open(newunit=ifev,file='__EValue',form='unformatted')
    read(ifev) nband_ev, nqi_, nsp_ev, nspc_
    if(ipr) write(stdo,ftox)'Read EValue: nband nqi nsp nspc nspx', nband, nqi, nsp,nspc,nspx
    if(nband_ev/=nband) call rx( 'init_readeigen:nband_ev/=nband')
    if(nsp_ev  /= nsp)  call rx( 'init_readeigen:nsp_ev/=nsp')
    if(nqi     /= nqi_) call rx( 'init_readeigen:nqi/=nqi_')
    if(nspc    /= nspc_)call rx( 'init_readeigen:nspc/=nspc_')
    nqixx=nqi
    allocate(evud(nband,nqi_,nsp),qtti_(3,nqi_))
    read(ifev) qtti_(1:3,1:nqi_)
    read(ifev) evud(1:nband, 1:nqi, 1:nspx )
    close(ifev)
    if(debug) then
       do is= 1,nspx
          do ik= 1,nqi
             do ib= 1,nband
               if(evud(ib,ik,is)<1d10.and.ipr) & !Set huge number for padding in sugw.f90
                    write(6,"('ib ik e=',3i5,f13.5,2x,3f9.4)") ib,ik,is,evud(ib,ik,is), qtti_(1:3,ik)
             enddo
          enddo
       enddo
       write(6,*)'init_readeigen:end'
    endif
    leval= minval(evud)
    init=.false.
  end subroutine init_readeigen
  real(8) function lowesteval()
    lowesteval=leval
  end function lowesteval
  
  subroutine init_readeigen2()    ! this should be called after init_readgeigen
    implicit none
    integer:: iq,is,ifiqg,ikp, isx,ikpisp,verbose,ifoc, i1,i2,i3,i4,i5,iorb,iorbold,i
    logical :: keepeigen
    character(8) :: xt
    call readmnla_cphi()
    keepeig = keepeigen()
    init2=.false.
    if(Keepeig       .and.ipr) write(6,*)' KeepEigen=T; readin geig and cphi into m_readeigen'
    if(( .NOT. Keepeig).and.ipr) write(6,*)' KeepEigen=F; not keep geig and cphi in m_readeigen'
    i=openm(newunit=ifcphim,file='__CPHI',recl=mrecb) ! Obata moved openm here, bug was 'openm after return
    i=openm(newunit=ifgeigm,file='__GEIG',recl=mrecg) ! in the case of keepeig=F ' fix at 2024-10-15
    if (mrecb /= ndima*nspc*nband*16) then
      if(ipr) write(stdo,'(a,2i0)') ' readeigen: __CPHI mrecb vs ndima*nspc*nband*16 = ', mrecb, ndima*nspc*nband*16
      call rx('readeigen: __CPHI record size inconsistent with GW product basis ndima.'// &
              ' Re-run sugw and regenerate PB.toml to match current ctrlg parameters.')
    endif
    if( .NOT. Keepeig) call keep_wfs_init() ! allocate for keep wfs
    if( .NOT. keepeig) return
    ! CPU and GPU: allocate geig/cphi in node-shared memory; only the node-root rank
    ! reads from disk; a barrier publishes the data to all ranks on the node.
    ! GPU build transfers one q-point slice to device per readgeigf_mpi/readcphif_mpi
    ! call instead of keeping the full arrays in VRAM.
    call shm_init()
    call shm_alloc_c8_4d(geig, ngpmx*nspc, nband, nqi, nspx, geig_shm_id)
    call shm_alloc_c8_4d(cphi, ndima*nspc, nband, nqi, nspx, cphi_shm_id)
    if(shm_rank()==0) then
       do ikp= 1,nqi
          do is= 1,nspx
             ikpisp= is + nsp*(ikp-1)
             i=readm(ifcphim,rec=ikpisp, data=cphi(1:ndima*nspc,1:nband,ikp,is))
             if(ngpmx/=0) i=readm(ifgeigm,rec=ikpisp,data=geig(1:ngpmx*nspc,1:nband,ikp,is))
          enddo
       enddo
    endif
    call shm_barrier(geig_shm_id)
    call shm_barrier(cphi_shm_id)
    if(keepeig)i=closem(ifcphim)
    if(keepeig)i=closem(ifgeigm)
  end subroutine init_readeigen2
  subroutine readmnla_cphi()
    !! === readin @MNLA_CPHI for rotation of MTO part of eigenfunction cphi ===
    implicit none
    integer:: iq,is,ifiqg,ikp, isx,ikpisp,verbose,ifoc, i1,i2,i3,i4,i5,iorb,iorbold
    open(newunit=ifoc,file='@MNLA_CPHI')
    read(ifoc,*)
    norbtx=0
    do
       read(ifoc,*,end=106) i1,i2,i3,i4,i5,iorb
       if(iorb>norbtx) norbtx=iorb
    enddo
106 continue
    if(ipr) write(6,*) ' end of readin @MNLA_CPHI: norbtx=',norbtx
    rewind ifoc
    read(ifoc,*)
    allocate(l_tbl(norbtx),k_tbl(norbtx),ibas_tbl(norbtx),offset_tbl(norbtx))
    iorbold=0
    do
       read(ifoc,*,end=107)i1,i2,i3,i4,i5,iorb
       if(iorbold/=iorb) then
          k_tbl(iorb)=i2
          l_tbl(iorb)=i3
          ibas_tbl(iorb)=i4
          offset_tbl(iorb)=i5-1
          iorbold=iorb
       endif
    enddo
107 continue
    allocate(offset_rev_tbl(maxval(ibas_tbl),0:maxval(l_tbl),maxval(k_tbl)))
    offset_rev_tbl=-99999
    do iorb=1,norbtx
       offset_rev_tbl(ibas_tbl(iorb),l_tbl(iorb),k_tbl(iorb))= offset_tbl(iorb)
    enddo
    close(ifoc)
  end subroutine readmnla_cphi
  !reaadcmlo moved to m_mlo_wfs
  !readcphiw/readgeigw moved to m_wan_wfs
end module m_readeigen
