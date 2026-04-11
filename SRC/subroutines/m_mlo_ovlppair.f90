module m_mlo_ovlppair
  use m_zmel,     only: build_zmel, zmel, Mptauof_zmel, set_m2e_prod_basis
  use m_mlo_scrw, only: nnmlo_init => nnwf_init, nnmlo => nnwf, nnmlo_mask => nnwf_mask
  use m_mlo_ham,  only: read_ham_rs, nmlo => ndimMTO
  use m_mpi,      only: MPI__root, MPI__rank, MPI__size, MPI__reduceSum, comm
  use m_lgunit,   only: m_lgunit_init, stdo
  use m_genallcf_v3, only: Genallcf_v3, nspin
  use m_read_bzdata, only: Read_BZDATA, nqbz, qbz, nqibz, qibz, nq0i, q0i
  use m_readgwinput, only: ReadGWinputKeys
  use m_readqg,      only: Readngmx2
  use m_hamindex,    only: Readhamindex
  use m_readeigen,   only: Init_readeigen, Init_readeigen2
  use m_itq,         only: Setitq_mlo
  use m_readVcoud,   only: Readvcoud, ngb, ReleaseZcousq
  use m_blas,        only: zmm => zmm_h, zmv => zmv_h, m_op_C, m_op_T
  use m_lattic,   only: plat=>lat_plat
  use m_mpiio, only: openm, readm, writem, openedm, closem

  use m_mlo_ham, only: ib_tableI
  use m_mlo_scrw, only: pair_site
  use m_ftox

  use m_lmfinit,only: m_lmfinit_init
  use m_lattic,only: m_lattic_init
  use m_mksym,only: m_mksym_init
  use m_mpitk, only: m_mpitk_init
  implicit none
  !! ovlppair_t and geometric data alongside it. Populated by read_ovlppair_t.
  complex(8), allocatable :: ovlppair_t(:,:,:)   !! (npairmx_f, nnmlo, nnmlo) real-space overlap pairs
  integer :: npairmx_f = 0, nbas_f = 0
  real(8) :: plat_f(3,3)
  integer, allocatable :: npair_f(:,:)            !! (nbas_f, nbas_f)
  integer, allocatable :: nqwgt_f(:,:,:)          !! (npairmx_f, nbas_f, nbas_f)
  integer, allocatable :: nlat_f(:,:,:,:)         !! (3, npairmx_f, nbas_f, nbas_f)
  public :: init_ovlppair, build_ovlppair_q, get_ovlppair_q
contains
  subroutine init_ovlppair(comm_in)
    integer, intent(in) :: comm_in
    call M_lgunit_init()
    call m_mpitk_init(comm_in)
    call m_lmfinit_init('ovlppair', int(comm_in)) ! Read ctrlp into module m_lmfinit.
    call m_lattic_init()                        ! lattice setup (for ewald sum)
    call m_mksym_init()                         ! symmetry go into m_lattic and m_mksym
    call Genallcf_v3(incwfx=0) !Basic data. incwfin= 0 takes 'ForX0 for core' in GWinput
    call Read_BZDATA()         !Readin BZDATA. See m_read_bzdata in gwsrc/rwbzdata.f
    call ReadGWinputKeys()     !Readin dataset in GWinput
    call Readngmx2()           !Get ngpmx and ngcmx in m_readqg
    !! Get space-group transformation information. See header of mptaouof.
    !! But we only use symops=E in hx0fp0 mode. c.f. hsfp0.sc
    call Mptauof_zmel(symops=reshape([1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0],[3,3]),ng=1)
    !! Rdpp gives ppbrd: radial integrals and cgr = rotated cg coeffecients. --> call Rdpp(ngrpx,symope) is moved to Mptauof_zmel \in m_zmel
    call Readhamindex()
    call Init_readeigen() ! Initialization of readEigen !readin m_hamindex
    call Init_readeigen2()
    call read_ham_rs()
    call nnmlo_init(nnwf_size_reduction=.true.)
    call Setitq_mlo(nmlo) ! Set itq in m_zmel
    ! call set_qibz(plat, qbz, nqbz, symops, ngrp)
  end subroutine init_ovlppair

  ! subroutine build_ovlppair_t(spinflip)
  !   logical, intent(in) :: spinflip
  !   integer :: isp, isp_k, isp_kq, iq, ik, ipr, npr, istat
  !   real(8) :: q(3)
  !   complex(8), allocatable :: ovlppair_q(:,:,:), zmel_nnmlo(:,:)
  !
  !   IspLoop: do isp = 1, nspin
  !     isp_k = isp
  !     isp_kq = isp
  !     if(spinflip) isp_kq = 3-isp
  !     allocate(ovlppair_q(nnmlo,nnmlo,nqibz), source = (0d0,0d0))
  !     IqLoop: do iq = 1, nqibz
  !       if(mod(iq-1,mpi__size)/=mpi__rank) cycle IqLoop
  !       q = qibz(:,iq)
  !       call Readvcoud(q, iq, NoVcou=.false.) !Readin vcousq,zcousq ngb ngc for the Coulomb matrix
  !       npr = ngb
  !       call set_m2e_prod_basis(npr=npr)
  !       call ReleaseZcousq()
  !       allocate(zmel_nnmlo(npr,nnmlo))
  !       do ik=1, nqbz
  !         write(stdo,ftox) 'yyy', iq, ik
  !         call build_zmel(q=q+qbz(:,ik), kvec=q, irot=1, rkvec=q, ns1=1, ns2=nmlo, ispm=isp_k, &
  !                         nqini=1, nqmax=nmlo, ispq=isp_kq, nctot=0, ncc=0, zmelconjg=.true., &
  !                         is_m_basis=.false., mpi_mode=.false., mlo_mode=.true.)
  !         do ipr=1, npr
  !           zmel_nnmlo(ipr,1:nnmlo) = pack(reshape(zmel(ipr,:,:), shape=[nmlo*nmlo]), mask=nnmlo_mask)
  !         enddo
  !         istat = zmm(zmel_nnmlo, zmel_nnmlo, ovlppair_q(:,:,iq), nnmlo, nnmlo, npr, beta=(1d0,0d0), opA=m_op_C)
  !       enddo
  !       deallocate(zmel_nnmlo)
  !     enddo IqLoop
  !
  !     FourieTransform:block
  !       use m_gennlat,  only: m_gennlat_init, npairmx, npair, nlat, nqwgt
  !       use m_keyvalue, only: getkeyvalue
  !       use m_lmfinit,  only: nbas
  !       use m_mlo_ham,  only: ib_tableI, ib_tableM
  !       use m_mlo_scrw, only: pair_site
  !       use m_setqibz_lmfham,only: qibz,irotq,irotg,ndiff,iqbzrep,qbzii,igiqibz,nqibz,iqii,wiqibz,ngx,igx
  !       integer :: nnn(3), ii, jj, ib1, ib2, np, ibt1, ibt2, it, iqbz, i, j, ifile
  !       integer, allocatable :: jdims(:), idims(:)
  !       complex(8), allocatable :: phases(:)
  !       complex(8), parameter :: img=(0d0,1d0)
  !       real(8), parameter ::pi=4d0*atan(1d0)
  !       real(8) :: qp(3)
  !       call getkeyvalue("GWinput", "n1n2n3", nnn,3)
  !       call m_gennlat_init(nnn) !for interpolation of Hamiltonian
  !       if(mpi__root) write(stdo,ftox) 'info about FT', nnn
  !       if(mpi__root) write(stdo,ftox) 'info about FT', npairmx, npair(:,:)
  !       allocate(ovlppair_t(npairmx,nnmlo,nnmlo), phases(npairmx))
  !       ovlppair_t(:,:,:) = (0d0, 0d0)
  !       do iqbz = 1, nqibz
  !         if(mod(iqbz-1,mpi__size)/=mpi__rank) cycle
  !         qp = qbz(:,iqbz)
  !         iqibz = irotq(iqbz)
  !         !! 対称性どーやるんや。。
  !         call rotmatMTO(igg=irotg(iqbz),q=qibz(:,iqibz),qtarget=qp+matmul(qlat,ndiff(:,iqibz)),ndimh=nMTO,rotmat=rotmat)
  !         forall(i=1:ndimMTO,j=1:ndimMTO) rotmatt(i,j)=rotmat(ix(i),ix(j))
  !         do ibt1 = 1, size(ib_tableI)
  !           do ibt2 = 1, size(ib_tableI)
  !             ib1 = ib_tableI(ibt1)
  !             ib2 = ib_tableI(ibt2)
  !             np = npair(ib1,ib2)
  !             idims = pack([(i,i=1,nnmlo)], mask=(pair_site(:,1)==ib1))
  !             jdims = pack([(j,j=1,nnmlo)], mask=(pair_site(:,2)==ib2))
  !             phases(1:np) = [(1d0/dble(nqbz)*exp(img*2d0*pi* sum(qp*(matmul(plat,nlat(:,it,ib1,ib2))))),it=1,np)]
  !             do concurrent(it=1:np, ii=1:size(idims), jj=1:size(jdims))
  !               ovlppair_t(it,idims(ii),jdims(jj)) = ovlppair_t(it,idims(ii),jdims(jj)) + ovlppair_q(idims(ii),jdims(jj),iqbz)*phases(it)
  !             enddo
  !           enddo
  !         enddo
  !       enddo
  !       call MPI__ReduceSum(0, ovlppair_t, size(ovlppair_t))
  !       WriteOvlpPairT:if(mpi__root) then
  !         open(newunit=ifile, file=ovlppair_t_fname(isp,logical(spinflip)), form='unformatted')
  !         write(ifile) npairmx, nbas
  !         write(ifile) plat(1:3,1:3), npair(1:nbas,1:nbas), nlat(1:3,1:npairmx,1:nbas,1:nbas), nqwgt(1:npairmx,1:nbas,1:nbas)
  !         write(ifile) ovlppair_t(:,:,:)
  !         close(ifile)
  !         do ibt2 = 1, size(ib_tableI)
  !           do ibt1 = 1, size(ib_tableI)
  !             ib1 = ib_tableI(ibt1)
  !             ib2 = ib_tableI(ibt2)
  !             np = npair(ib1,ib2)
  !             idims = pack([(i,i=1,nnmlo)], mask=(pair_site(:,1)==ib1))
  !             jdims = pack([(j,j=1,nnmlo)], mask=(pair_site(:,2)==ib2))
  !             phases(1:np) = 1d0/dble(nqwgt(1:np,ib1,ib2))
  !             do jj = 1, size(jdims)
  !               do ii = 1, size(idims)
  !                 write(stdo,ftox) 'nwf nwf index, sum ovlppair:', idims(ii), jdims(jj), &
  !                                  sum([(ovlppair_t(it,idims(ii),jdims(jj))*phases(it),it=1,np)])
  !               enddo
  !             enddo
  !           enddo
  !         enddo
  !       endif WriteOvlpPairT
  !     endblock FourieTransform
  !     deallocate(ovlppair_q)
  !   enddo IspLoop
  ! end subroutine build_ovlppair_t

  subroutine read_ovlppair_t(isp, spinflip)
    integer, intent(in) :: isp
    logical, intent(in) :: spinflip
    integer :: ifile
    open(newunit=ifile, file=ovlppair_t_fname(isp,spinflip), form='unformatted')
    if(allocated(ovlppair_t)) deallocate(ovlppair_t)
    if(allocated(npair_f))    deallocate(npair_f)
    if(allocated(nlat_f))     deallocate(nlat_f)
    if(allocated(nqwgt_f))    deallocate(nqwgt_f)
    read(ifile) npairmx_f, nbas_f
    allocate(ovlppair_t(npairmx_f, nnmlo, nnmlo))
    allocate(npair_f(nbas_f,nbas_f), nlat_f(3,npairmx_f,nbas_f,nbas_f), nqwgt_f(npairmx_f,nbas_f,nbas_f))
    read(ifile) plat_f(1:3,1:3), npair_f, nlat_f, nqwgt_f
    read(ifile) ovlppair_t(:,:,:)
    close(ifile)
  end subroutine read_ovlppair_t

  function get_ovlppair_q_ft(q, isp) result(ovlppair_q)
    real(8), intent(in) :: q(3)
    integer, intent(in) :: isp
    complex(8) :: ovlppair_q(nnmlo, nnmlo)
    integer :: ibt1, ibt2, ib1, ib2, it, i, j, jj, istat
    integer, allocatable :: jdims(:), idims(:)
    complex(8), allocatable :: phases(:)
    complex(8), parameter :: img=(0d0,1d0)
    real(8), parameter ::pi=4d0*atan(1d0)
    allocate(phases(npairmx_f))
    ovlppair_q(:,:) = (0d0, 0d0)
    do ibt2 = 1, size(ib_tableI)
      do ibt1 = 1, size(ib_tableI)
        ib2 = ib_tableI(ibt2)
        ib1 = ib_tableI(ibt1)
        idims = pack([(i,i=1,nnmlo)], mask=(pair_site(:,1)==ib1))
        jdims = pack([(j,j=1,nnmlo)], mask=(pair_site(:,2)==ib2))
        do it = 1, npair_f(ib1,ib2)
          phases(it) = 1d0/dble(nqwgt_f(it,ib1,ib2))*exp(-img*2d0*pi*sum(q*matmul(plat_f,nlat_f(:,it,ib1,ib2))))
        enddo
        do jj = 1, size(jdims)
          istat = zmv(ovlppair_t(1,idims(1),jdims(jj)), phases, ovlppair_q(idims(1),jdims(jj)), &
                      n=size(idims), m=npair_f(ib1,ib2), opA=m_op_T, lda=npairmx_f)
        enddo
      enddo
    enddo
  end function get_ovlppair_q_ft

  subroutine build_ovlppair_q(spinflip)
    use m_read_bzdata, only: nq0i, q0i
    logical, intent(in) :: spinflip
    integer :: isp, isp_k, isp_kq, iq, ik, ipr, npr, istat, mrecl, iq0i, ngc
    real(8) :: q(3), quu(3)
    complex(8), allocatable :: ovlppair_q(:,:), zmel_nnmlo(:,:)
    integer :: ifile, ifile_info
    mrecl = nnmlo*nnmlo*16
    if(mpi__root) then
      open(newunit=ifile_info, file='__MLOOvlpPairQ.info', form='unformatted', action='write')
      write(ifile_info) mrecl, nq0i
      write(ifile_info) q0i
    endif
    IspLoop: do isp = 1, nspin
      isp_k = isp
      isp_kq = isp
      if(spinflip) isp_kq = 3-isp
      if(openedm(ifile)) istat = closem(ifile)
      istat = openm(newunit=ifile, file=ovlppair_q_fname(isp,logical(spinflip)), recl=mrecl)
      IqLoop: do iq = nqibz+1, nqibz + nq0i
        iq0i = iq - nqibz
        if(mod(iq, mpi__size)/=mpi__rank) cycle IqLoop
        q = q0i(:,iq0i)
        call Readqg0('QGcou', q, quu, ngc) ! ngc: the number of IPW for the interaction matrix (in QGcou),
        call Readvcoud(q, iq, NoVcou=.false.) !Readin vcousq,zcousq ngb ngc for the Coulomb matrix
        npr = ngb
        call set_m2e_prod_basis(npr=npr) !set npr for zmel, E basis is not used for q=0
        call ReleaseZcousq()
        allocate(zmel_nnmlo(npr,nnmlo))
        allocate(ovlppair_q(nnmlo,nnmlo), source = (0d0,0d0))
        IkLoop: do ik=1, nqbz
          call build_zmel(q=q+qbz(:,ik), kvec=q, irot=1, rkvec=q, ns1=1, ns2=nmlo, ispm=isp_k, &
                          nqini=1, nqmax=nmlo, ispq=isp_kq, nctot=0, ncc=0, zmelconjg=.true., &
                          is_m_basis=.true., mpi_mode=.false., mlo_mode=.true.)
          do ipr=1, npr
            zmel_nnmlo(ipr,1:nnmlo) = pack(reshape(zmel(ipr,:,:), shape=[nmlo*nmlo]), mask=nnmlo_mask)
          enddo
          istat = zmm(zmel_nnmlo, zmel_nnmlo, ovlppair_q(:,:), nnmlo, nnmlo, npr, beta=(1d0,0d0), opA=m_op_C)
        enddo IkLoop
        ovlppair_q(:,:) = ovlppair_q(:,:)/dble(nqbz)
        istat = writem(ifile, rec=iq0i, data=ovlppair_q(:,:))
        write(stdo,ftox) 'q, npr, sum S^2(q)',q, npr, sum(ovlppair_q)
        deallocate(ovlppair_q, zmel_nnmlo)
      enddo IqLoop
      istat = closem(ifile)
    enddo IspLoop
  end subroutine build_ovlppair_q

  function get_ovlppair_q(q, isp, spinflip) result(ovlppair_q)
    real(8), intent(in) :: q(3)
    integer, intent(in) :: isp
    logical, intent(in) :: spinflip
    complex(8) :: ovlppair_q(nnmlo, nnmlo)
    integer, save :: ifile = -1, isp_prev = -1
    logical, save :: spinflip_prev = .false.
    integer, save :: mrecl_f, nq0i_f
    real(8), save, allocatable :: q0i_f(:,:)
    logical :: opened
    integer :: iq, istat, ifile_info
    if(isp /= isp_prev .or. spinflip .neqv. spinflip_prev) then
      open(newunit=ifile_info, file='__MLOOvlpPairQ.info', form='unformatted', action='read')
      read(ifile_info) mrecl_f, nq0i_f
      if(allocated(q0i_f)) deallocate(q0i_f)
      allocate(q0i_f(3,nq0i_f))
      read(ifile_info) q0i_f
      close(ifile_info)

      inquire(unit=ifile, opened=opened)
      if(opened) close(ifile)
      open(newunit=ifile, file=ovlppair_q_fname(isp, spinflip), &
           form='unformatted', access='direct', recl=mrecl_f, action='read')
      isp_prev = isp
      spinflip_prev = spinflip
    endif
    do iq = 1, nq0i_f
      if(all(abs(q0i_f(:,iq) - q) < 1d-8)) then
        read(ifile, rec=iq) ovlppair_q(:,:)
        return
      endif
    enddo
    call rx('get_ovlppair_q: q not found in q0i')
  end function get_ovlppair_q

  pure function ovlppair_t_fname(isp, spinflip) result(fname)
    integer, intent(in) :: isp
    logical, intent(in) :: spinflip
    character(:), allocatable :: fname
    if(isp==1 .and. .not.spinflip) fname = '__MLOOvlpPairT.UP'
    if(isp==2 .and. .not.spinflip) fname = '__MLOOvlpPairT.DN'
    if(isp==1 .and. spinflip)      fname = '__MLOOvlpPairT.UPDN'
    if(isp==2 .and. spinflip)      fname = '__MLOOvlpPairT.DNUP'
  end function ovlppair_t_fname

  pure function ovlppair_q_fname(isp, spinflip) result(fname)
    integer, intent(in) :: isp
    logical, intent(in) :: spinflip
    character(:), allocatable :: fname
    if(isp==1 .and. .not.spinflip) fname = '__MLOOvlpPairQ.UP'
    if(isp==2 .and. .not.spinflip) fname = '__MLOOvlpPairQ.DN'
    if(isp==1 .and. spinflip)      fname = '__MLOOvlpPairQ.UPDN'
    if(isp==2 .and. spinflip)      fname = '__MLOOvlpPairQ.DNUP'
  end function ovlppair_q_fname
end module m_mlo_ovlppair
