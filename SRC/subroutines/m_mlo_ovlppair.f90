module m_mlo_ovlppair
  use m_zmel,     only: build_zmel, zmel, Mptauof_zmel, set_nbb_zmel
  use m_mlo_scrw, only: nnmlo_init => nnwf_init, nnmlo => nnwf, nnmlo2_mask => nnwf2_mask, trace
  use m_mlo_ham,  only: read_ham_rs, nmlo => ndimMTO
  use m_mpi,      only: MPI__root, MPI__rank, MPI__size, MPI__reduceSum, comm
  use m_lgunit,   only: m_lgunit_init, stdo
  use m_genallcf_v3, only: Genallcf_v3, nspin
  use m_read_bzdata, only: Read_BZDATA, nqbz, qbz, nqibz, qibz, nq0i, q0i, wbz, irk, ngrp
  use m_readgwinput, only: ReadGWinputKeys
  use m_readqg,      only: Readngmx2
  use m_hamindex,    only: Readhamindex
  use m_readeigen,   only: Init_readeigen, Init_readeigen2
  use m_itq,         only: Setitq_mlo
  use m_blas,        only: zmm => zmm_h, zmv => zmv_h, m_op_C, m_op_T
  use m_lattic,   only: plat=>lat_plat
  use m_mpiio, only: openm, readm, writem, openedm, closem

    use m_rdpp, only: nbloch
    use m_readVcoud, only: ReadVcoud, ngc, npr=>ngb
    use m_read_ppovl, only: ppovlp, getppx2, ngcread
    use m_lapack, only: zminv => zminv_h, zhev => zhev_h

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
  public :: init_build_ovlppair
  public :: build_ovlppair_q, get_ovlppair_q
  ! public :: build_ovlppair_t
contains
  subroutine init_build_ovlppair(comm_in)
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
  end subroutine init_build_ovlppair
  !
  ! subroutine build_ovlppair_t(spinflip)
  !   logical, intent(in) :: spinflip
  !   integer :: isp, kx, ipr, istat, mrecl, iq0i, i, isp_k, isp_kq, irot, kr, iq
  !   real(8) :: q(3), quu(3), evl(nnmlo), qibz_k(3), qbz_kr(3)
  !   complex(8), allocatable :: ovlppair_q(:,:), pbmovlp_inv(:,:), oinv_zmel(:,:,:), ovlppair4(:,:,:,:)
  !   complex(8), allocatable :: zmelk(:,:,:)
  !   integer :: ifile, ifile_info
  !
  !   mrecl = nnmlo*nnmlo*16
  !   if(mpi__root) then
  !     open(newunit=ifile_info, file='__MLOOvlpPairQ.info', form='unformatted', action='write')
  !     write(ifile_info) mrecl, nq0i
  !     write(ifile_info) q0i
  !   endif
  !
  !   IspLoop: do isp = 1, nspin
  !     isp_k = isp
  !     isp_kq = isp
  !     if(spinflip) isp_kq = 3-isp
  !
  !     if(openedm(ifile)) istat = closem(ifile)
  !     istat = openm(newunit=ifile, file=ovlppair_q_fname(isp,logical(spinflip)), recl=mrecl)
  !
  !     IqLoop: do iq = 1, nqibz
  !       if(mod(iq-1,mpi__size)/=mpi__rank) cycle IqLoop
  !       q = qibz(:,iq)
  !       call ReadVcoud(q, iq, NoVcou=.false.) !Readin vcousq,zcousq ngb ngc for the Coulomb matrix
  !       call set_nbb_zmel(npr)                !set npr for zmel, E basis is not used for q=0
  !
  !       call getppx2(q, get_ppovlp=.true.)     !set ppovlp
  !       allocate(pbmovlp_inv(npr,npr), source = (0d0,0d0))
  !       forall(i=1:nbloch) pbmovlp_inv(i,i) = (1d0,0d0)
  !       pbmovlp_inv(nbloch+1:npr,nbloch+1:npr) = ppovlp(1:ngc,1:ngc)
  !       istat = zminv(pbmovlp_inv, n=npr)
  !       if(ngc /= ngcread) call rx('ngc /= ngcread')
  !
  !       allocate(zmelk(npr,nmlo,nmlo), source = (0d0, 0d0))
  !       KxLoop: do kx=1, nqibz
  !         IrotLoop: do irot=1, ngrp
  !           kr = irk(kx, irot)
  !           if(kr==0) cycle IrotLoop
  !           qibz_k = qibz(:,kx)
  !           qbz_kr = qbz (:,kr)
  !           call build_zmel(q=q, kvec=qibz_k, irot=irot, rkvec=qbz_kr, ns1=1, ns2=nmlo, ispm=isp_k, &
  !                           nqini=1, nqmax=nmlo, ispq=isp_kq, nctot=0, ncc=0, zmelconjg=.false., &
  !                           is_m_basis=.true., mpi_mode=.false., mlo_mode=.true.)
  !           zmelk(:,:,:) = zmelk(:,:,:) + zmel(:,:,:)
  !         enddo IrotLoop
  !       enddo KxLoop
  !       zmelk = zmelk/dble(nqbz)
  !
  !       allocate(ovlppair4(nmlo,nmlo,nmlo,nmlo), oinv_zmel(npr,nmlo,nmlo))
  !       allocate(ovlppair_q(nnmlo,nnmlo), source=(0d0,0d0))
  !       istat = zmm(pbmovlp_inv, zmelk, oinv_zmel, npr, nmlo**2, npr)
  !       istat = zmm(zmelk, oinv_zmel, ovlppair4, nmlo**2, nmlo**2, npr, opA=m_op_C)
  !       ovlppair_q(:,:) = reshape(pack(reshape(reshape(ovlppair4, shape(ovlppair4), order=[2,1,4,3]), &
  !                                   [size(ovlppair4)]), mask=nnmlo2_mask), shape = [nnmlo, nnmlo])
  !       istat = writem(ifile, rec=iq, data=ovlppair_q(:,:))
  !       write(stdo,ftox) 'q, npr, sum S^2(q)',q, npr, sum(ovlppair_q), trace(ovlppair_q)
  !       deallocate(ovlppair_q)
  !     enddo IqLoop

      ! FourieTransform:block
      !   use m_gennlat,  only: m_gennlat_init, npairmx, npair, nlat, nqwgt
      !   use m_keyvalue, only: getkeyvalue
      !   use m_lmfinit,  only: nbas
      !   use m_mlo_ham,  only: ib_tableI, ib_tableM
      !   use m_mlo_scrw, only: pair_site
      !   use m_setqibz_lmfham,only: qibz,irotq,irotg,ndiff,iqbzrep,qbzii,igiqibz,nqibz,iqii,wiqibz,ngx,igx
      !   integer :: nnn(3), ii, jj, ib1, ib2, np, ibt1, ibt2, it, iqbz, i, j, ifile
      !   integer, allocatable :: jdims(:), idims(:)
      !   complex(8), allocatable :: phases(:)
      !   complex(8), parameter :: img=(0d0,1d0)
      !   real(8), parameter ::pi=4d0*atan(1d0)
      !   real(8) :: qp(3)
      !   call getkeyvalue("GWinput", "n1n2n3", nnn,3)
      !   call m_gennlat_init(nnn) !for interpolation of Hamiltonian
      !   if(mpi__root) write(stdo,ftox) 'info about FT', nnn
      !   if(mpi__root) write(stdo,ftox) 'info about FT', npairmx, npair(:,:)
      !   allocate(ovlppair_t(npairmx,nnmlo,nnmlo), phases(npairmx))
      !   ovlppair_t(:,:,:) = (0d0, 0d0)
      !   do iqbz = 1, nqibz
      !     if(mod(iqbz-1,mpi__size)/=mpi__rank) cycle
      !     qp = qbz(:,iqbz)
      !     iqibz = irotq(iqbz)
      !     !! 対称性どーやるんや。。
      !     call rotmatMTO(igg=irotg(iqbz),q=qibz(:,iqibz),qtarget=qp+matmul(qlat,ndiff(:,iqibz)),ndimh=nMTO,rotmat=rotmat)
      !     forall(i=1:ndimMTO,j=1:ndimMTO) rotmatt(i,j)=rotmat(ix(i),ix(j))
      !     do ibt1 = 1, size(ib_tableI)
      !       do ibt2 = 1, size(ib_tableI)
      !         ib1 = ib_tableI(ibt1)
      !         ib2 = ib_tableI(ibt2)
      !         np = npair(ib1,ib2)
      !         idims = pack([(i,i=1,nnmlo)], mask=(pair_site(:,1)==ib1))
      !         jdims = pack([(j,j=1,nnmlo)], mask=(pair_site(:,2)==ib2))
      !         phases(1:np) = [(1d0/dble(nqbz)*exp(img*2d0*pi* sum(qp*(matmul(plat,nlat(:,it,ib1,ib2))))),it=1,np)]
      !         do concurrent(it=1:np, ii=1:size(idims), jj=1:size(jdims))
      !           ovlppair_t(it,idims(ii),jdims(jj)) = ovlppair_t(it,idims(ii),jdims(jj)) + ovlppair_q(idims(ii),jdims(jj),iqbz)*phases(it)
      !         enddo
      !       enddo
      !     enddo
      !   enddo
      !   call MPI__ReduceSum(0, ovlppair_t, size(ovlppair_t))
      !   WriteOvlpPairT:if(mpi__root) then
      !     open(newunit=ifile, file=ovlppair_t_fname(isp,logical(spinflip)), form='unformatted')
      !     write(ifile) npairmx, nbas
      !     write(ifile) plat(1:3,1:3), npair(1:nbas,1:nbas), nlat(1:3,1:npairmx,1:nbas,1:nbas), nqwgt(1:npairmx,1:nbas,1:nbas)
      !     write(ifile) ovlppair_t(:,:,:)
      !     close(ifile)
      !     do ibt2 = 1, size(ib_tableI)
      !       do ibt1 = 1, size(ib_tableI)
      !         ib1 = ib_tableI(ibt1)
      !         ib2 = ib_tableI(ibt2)
      !         np = npair(ib1,ib2)
      !         idims = pack([(i,i=1,nnmlo)], mask=(pair_site(:,1)==ib1))
      !         jdims = pack([(j,j=1,nnmlo)], mask=(pair_site(:,2)==ib2))
      !         phases(1:np) = 1d0/dble(nqwgt(1:np,ib1,ib2))
      !         do jj = 1, size(jdims)
      !           do ii = 1, size(idims)
      !             write(stdo,ftox) 'nwf nwf index, sum ovlppair:', idims(ii), jdims(jj), &
      !                              sum([(ovlppair_t(it,idims(ii),jdims(jj))*phases(it),it=1,np)])
      !           enddo
      !         enddo
      !       enddo
      !     enddo
      !   endif WriteOvlpPairT
      ! endblock FourieTransform
  !   enddo IspLoop
  ! end subroutine build_ovlppair_t

  ! subroutine read_ovlppair_t(isp, spinflip)
  !   integer, intent(in) :: isp
  !   logical, intent(in) :: spinflip
  !   integer :: ifile
  !   open(newunit=ifile, file=ovlppair_t_fname(isp,spinflip), form='unformatted')
  !   if(allocated(ovlppair_t)) deallocate(ovlppair_t)
  !   if(allocated(npair_f))    deallocate(npair_f)
  !   if(allocated(nlat_f))     deallocate(nlat_f)
  !   if(allocated(nqwgt_f))    deallocate(nqwgt_f)
  !   read(ifile) npairmx_f, nbas_f
  !   allocate(ovlppair_t(npairmx_f, nnmlo, nnmlo))
  !   allocate(npair_f(nbas_f,nbas_f), nlat_f(3,npairmx_f,nbas_f,nbas_f), nqwgt_f(npairmx_f,nbas_f,nbas_f))
  !   read(ifile) plat_f(1:3,1:3), npair_f, nlat_f, nqwgt_f
  !   read(ifile) ovlppair_t(:,:,:)
  !   close(ifile)
  ! end subroutine read_ovlppair_t
  !
  ! function get_ovlppair_q_ft(q, isp) result(ovlppair_q)
  !   real(8), intent(in) :: q(3)
  !   integer, intent(in) :: isp
  !   complex(8) :: ovlppair_q(nnmlo, nnmlo)
  !   integer :: ibt1, ibt2, ib1, ib2, it, i, j, jj, istat
  !   integer, allocatable :: jdims(:), idims(:)
  !   complex(8), allocatable :: phases(:)
  !   complex(8), parameter :: img=(0d0,1d0)
  !   real(8), parameter ::pi=4d0*atan(1d0)
  !   allocate(phases(npairmx_f))
  !   ovlppair_q(:,:) = (0d0, 0d0)
  !   do ibt2 = 1, size(ib_tableI)
  !     do ibt1 = 1, size(ib_tableI)
  !       ib2 = ib_tableI(ibt2)
  !       ib1 = ib_tableI(ibt1)
  !       idims = pack([(i,i=1,nnmlo)], mask=(pair_site(:,1)==ib1))
  !       jdims = pack([(j,j=1,nnmlo)], mask=(pair_site(:,2)==ib2))
  !       do it = 1, npair_f(ib1,ib2)
  !         phases(it) = 1d0/dble(nqwgt_f(it,ib1,ib2))*exp(-img*2d0*pi*sum(q*matmul(plat_f,nlat_f(:,it,ib1,ib2))))
  !       enddo
  !       do jj = 1, size(jdims)
  !         istat = zmv(ovlppair_t(1,idims(1),jdims(jj)), phases, ovlppair_q(idims(1),jdims(jj)), &
  !                     n=size(idims), m=npair_f(ib1,ib2), opA=m_op_T, lda=npairmx_f)
  !       enddo
  !     enddo
  !   enddo
  ! end function get_ovlppair_q_ft

  ! pure function ovlppair_t_fname(isp, spinflip) result(fname)
  !   integer, intent(in) :: isp
  !   logical, intent(in) :: spinflip
  !   character(:), allocatable :: fname
  !   if(isp==1 .and. .not.spinflip) fname = '__MLOOvlpPairT.UP'
  !   if(isp==2 .and. .not.spinflip) fname = '__MLOOvlpPairT.DN'
  !   if(isp==1 .and. spinflip)      fname = '__MLOOvlpPairT.UPDN'
  !   if(isp==2 .and. spinflip)      fname = '__MLOOvlpPairT.DNUP'
  ! end function ovlppair_t_fname

  subroutine build_ovlppair_q(spinflip)
   use,intrinsic :: ieee_arithmetic
    logical, intent(in) :: spinflip
    integer :: isp, kx, ipr, istat, mrecl, iq0i, i, ispin1, ispin2, irot, kr
    real(8) :: q(3), quu(3), evl(nnmlo), qibz_k(3), qbz_kr(3)
    complex(8), allocatable :: ovlppair_q(:,:), pbmovlp_inv(:,:), oinv_zmel(:,:,:), ovlppair4(:,:,:,:)
    complex(8), allocatable :: zmelk(:,:,:)
    integer :: ifile, ifile_info
    mrecl = nnmlo*nnmlo*16
    mrecl = 16*nmlo**4
    if(mpi__root) then
      open(newunit=ifile_info, file='__MLOOvlpPairQ.info', form='unformatted', action='write')
      write(ifile_info) mrecl, nq0i
      write(ifile_info) q0i
      ! write(stdo,ftox) mrecl, nq0i
      ! write(stdo,ftox) q0i
    endif
    IspLoop: do isp = 1, nspin
      ispin1 = isp
      ispin2 = isp
      if(spinflip) ispin2 = 3-isp !oppsite spin

      if(openedm(ifile)) istat = closem(ifile)
      istat = openm(newunit=ifile, file=ovlppair_q_fname(isp,logical(spinflip)), recl=mrecl)

      IqLoop: do iq0i = 1, nq0i
        if(mod(iq0i-1, mpi__size)/=mpi__rank) cycle IqLoop
        q = q0i(:,iq0i)
        call ReadVcoud(q, iq0i, NoVcou=.true.) !set ngc, ngb (=npr) used in zmel
        call set_nbb_zmel(npr)                 !set npr for zmel, E basis is not used for q=0
        call getppx2(q, get_ppovlp=.true.)     !set ppovlp
        if(ngc /= ngcread) call rx('ngc /= ngcread')

        write(stdo,ftox) 'nbloch, ngc, npr', nbloch, ngc, npr, nbloch + ngc
        ! block
        !   complex(8), allocatable :: ppovlp_check(:,:)
        !   real(8), allocatable :: evl_check(:)
        !   allocate(ppovlp_check, source=ppovlp)
        !   allocate(evl_check(ngc))
        !   istat = zhev(ppovlp_check, n=ngc, evl=evl_check)
        !   write(*,*) 'ppovlp max eigenvalue:', maxval(evl_check)
        !   write(*,*) 'ppovlp min eigenvalue:', minval(evl_check)
        ! endblock

        allocate(pbmovlp_inv(npr,npr), source=(0d0,0d0))
        forall(i=1:nbloch) pbmovlp_inv(i,i) = (1d0,0d0)
        pbmovlp_inv(nbloch+1:npr,nbloch+1:npr) = ppovlp(1:ngc,1:ngc)
        istat = zminv(pbmovlp_inv, n=npr)

        allocate(zmelk(npr,nmlo,nmlo), source = (0d0, 0d0))
        IkLoop: do kx=1, nqbz
          call build_zmel(q=qbz(:,kx), kvec=q, irot=1, rkvec=q, ns1=1, ns2=nmlo, ispm=ispin1, &
                          nqini=1, nqmax=nmlo, ispq=ispin2, nctot=0, ncc=0, zmelconjg=.false., &
                          is_m_basis=.true., mpi_mode=.false., mlo_mode=.true.)
          ! call build_zmel(q=q+qbz(:,kx), kvec=q, irot=1, rkvec=q, ns1=1, ns2=nmlo, ispm=ispin1, &
          !                 nqini=1, nqmax=nmlo, ispq=ispin2, nctot=0, ncc=0, zmelconjg=.false., &
          !                 is_m_basis=.true., mpi_mode=.false., mlo_mode=.true.)
          if(any(ieee_is_nan(dble(zmel)))) write(stdo,ftox) "NaN in Real zmel", iq0i, kx, q, npr, nmlo
          if(any(ieee_is_nan(imag(zmel)))) write(stdo,ftox) "NaN in Imag zmel", iq0i, kx, q, npr, nmlo
          zmelk(:,:,:) = zmelk(:,:,:) + zmel(:,:,:)
          ! istat = zmm(pbmovlp_inv, zmel, oinv_zmel, npr, nmlo**2, npr)
          ! istat = zmm(zmel, oinv_zmel, ovlppair4, nmlo**2, nmlo**2, npr, opA=m_op_C)
          ! ovlppair_q = ovlppair_q + wbz(kx)*reshape(pack(reshape(reshape(ovlppair4, shape(ovlppair4), order=[1,2,3,4]), &
          !                                                [size(ovlppair4)]), mask=nnmlo2_mask), shape = [nnmlo, nnmlo])
        enddo IkLoop
        zmelk = zmelk/dble(nqbz)
        allocate(ovlppair4(nmlo,nmlo,nmlo,nmlo), oinv_zmel(npr,nmlo,nmlo))
        allocate(ovlppair_q(nnmlo,nnmlo), source=(0d0,0d0))

        istat = zmm(pbmovlp_inv, zmelk, oinv_zmel, npr, nmlo**2, npr)
        istat = zmm(zmelk, oinv_zmel, ovlppair4, nmlo**2, nmlo**2, npr, opA=m_op_C)

        ! istat = zmm(zmelk, zmelk, ovlppair4, nmlo**2, nmlo**2, npr, opA=m_op_C)

        ! ovlppair_q = reshape(pack(reshape(reshape(ovlppair4, shape(ovlppair4), order=[1,2,3,4]), &
        !                      [size(ovlppair4)]), mask=nnmlo2_mask), shape = [nnmlo, nnmlo])
        ! istat = writem(ifile, rec=iq0i, data=ovlppair_q(:,:))
        ovlppair4 = reshape(ovlppair4, shape(ovlppair4), order=[2,1,4,3])
        istat = writem(ifile, rec=iq0i, data=ovlppair4(:,:,:,:))
        block
          complex(8) :: ovlpp_trace
          integer :: i, j
          ovlpp_trace = 0d0
          do i =1,  nmlo
            do j =1, nmlo
             ovlpp_trace = ovlpp_trace + ovlppair4(i,j,i,j)
            enddo
          enddo
          write(stdo,ftox) 'ovlpp_trace:', q, ovlpp_trace
       endblock
       block
         real(8) :: evl(nmlo**2)
        ! ovlppair4 = (ovlppair4 + conjg(transpose(ovlppair4)))*0.5d0
        istat = zhev(ovlppair4, n=nmlo**2, evl=evl)
        write(*,*) 'max eigenvalue:', maxval(evl)
        write(*,*) 'min eigenvalue:', minval(evl)
        write(*,*) 'condition number:', maxval(evl)/minval(evl)
        write(*,*)  evl
      endblock
        deallocate(zmelk, ovlppair_q, ovlppair4, oinv_zmel, pbmovlp_inv)
      enddo IqLoop
      istat = closem(ifile)
    enddo IspLoop
  end subroutine build_ovlppair_q

  function get_ovlppair_q(q, isp, spinflip, minus_q) result(ovlppair_q)
    real(8), intent(in) :: q(3)
    integer, intent(in) :: isp
    logical, intent(in) :: spinflip
    logical, intent(in), optional :: minus_q
    complex(8) :: ovlppair_q(nnmlo, nnmlo)
    complex(8) :: ovlppair4(nmlo,nmlo,nmlo,nmlo)
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
        ! read(ifile, rec=iq) ovlppair_q(:,:)
        read(ifile, rec=iq) ovlppair4
        if(present(minus_q)) then
          if(minus_q) ovlppair4 = conjg(reshape(ovlppair4, shape=[nmlo,nmlo,nmlo,nmlo], order=[2,1,4,3]))
        endif
        ovlppair_q = reshape(pack(reshape(reshape(ovlppair4, shape(ovlppair4), order=[1,2,4,3]), &
                             [size(ovlppair4)]), mask=nnmlo2_mask), shape = [nnmlo, nnmlo])
        return
      endif
    enddo
    call rx('get_ovlppair_q: q not found in q0i')
  end function get_ovlppair_q

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
