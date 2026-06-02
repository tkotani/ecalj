module m_mlo_ham
  use m_HamPMT, only: plat, npair, nlat, nqwgt,nspc
  use m_lgunit, only: stdo
  use m_mpi, only: ipr
  use m_blas, only: zmm => zmm_h, zmv => zmv_h, m_op_T
  use m_lapack, only: zhgv => zhgv_h, zsv => zsv_h
  use m_ftox, only: ftox
  use m_cmdopt_registry, only: c0_socmatrix
  implicit none
  public :: read_ham_rs, calc_ham_eigen
  integer, protected, target :: ndimMTO, npairmx, nspx, nsite
  integer, allocatable, protected:: ib_tableM(:), l_tableM(:), k_tableM(:), ib_tableI(:)
  complex(8),allocatable, protected:: ovlmr(:,:,:,:), hammr(:,:,:,:), hammhsor(:,:,:,:) !npairmx, ndimMTO, ndimMTO, nspx order
  logical, protected :: socmatrix = .false.
contains
  subroutine read_ham_rs()! read RealSpace MTO Hamiltonian
    use m_cmdopt_registry, only: c0_socmatrix
    integer:: ifihmto, i
    socmatrix = c0_socmatrix
    open(newunit=ifihmto,file='HamRsMLO',form='unformatted', action='read')
    read(ifihmto) ndimMTO,npairmx,nspx !    allocate(ix(ndimMTO))
    if(ipr) write(stdo,ftox)'MTOHamiltonian: ndimMTO,npairmx,nspx=',ndimMTO,npairmx,nspx
    allocate(hammr(npairmx,ndimMTO,ndimMTO,nspx))
    allocate(ovlmr(npairmx,ndimMTO,ndimMTO,nspx))
    read(ifihmto) hammr
    if(socmatrix) then
      allocate(hammhsor(npairmx,ndimMTO,ndimMTO,3))
      read(ifihmto) hammhsor
    endif
    read(ifihmto) ovlmr
    allocate(ib_tableM(1:ndimMTO),k_tableM(1:ndimMTO),l_tableM(1:ndimMTO))
    read(ifihmto) ib_tableM(1:ndimMTO),k_tableM(1:ndimMTO),l_tableM(1:ndimMTO)
    close(ifihmto)
    if(ipr) write(stdo,*)'OK: Read HamRsMLO file! Use i-ioffib for setting <Worb>'
    ib_tableI = pack(ib_tableM(1:ndimMTO), [(all(ib_tableM(:i-1)/=ib_tableM(i)), i=1,ndimMTO)])
    if(ipr) write(stdo,ftox) 'Atomic sites in the primitive cell for MLO Hamiltonian: ', ib_tableI
    nsite = size(ib_tableI)
  end subroutine read_ham_rs

  subroutine calc_ham_eigen(isp, q, ev, evec, ovlp_evec, dual_evec)
    integer, intent(in) :: isp
    real(8), intent(in) :: q(3)
    real(8), intent(out) :: ev(:)
    complex(8), optional, intent(out) :: evec(:,:), ovlp_evec(:,:), dual_evec(:,:)
    complex(8), allocatable :: ovlm(:,:,:,:)
    complex(8), allocatable, target :: hamm(:,:,:,:)
    complex(8), allocatable :: ovlm_buf(:,:)
    complex(8), pointer :: hamm2d(:,:)
    complex(8), allocatable :: hammhso(:,:,:)
    real(8), parameter :: pi=4d0*atan(1d0)
    complex(8), parameter :: img=(0d0,1d0)
    complex(8) :: phases(npairmx)
    integer, allocatable :: idims(:), jdims(:)
    integer :: i, j, ib1, ib2, it, jj, ibt1, ibt2, istat, io, js, jsp_src, nspinor

    nspinor = merge(2, 1, socmatrix)
    allocate(ovlm(ndimMTO,nspinor,ndimMTO,nspinor), source=(0d0,0d0))
    allocate(hamm(ndimMTO,nspinor,ndimMTO,nspinor), source=(0d0,0d0))
    if(socmatrix) allocate(hammhso(ndimMTO,ndimMTO,3) ,source=(0d0,0d0))

    FourierTransform: do ibt2 = 1, size(ib_tableI)
      ib2 = ib_tableI(ibt2)
      jdims = pack([(j,j=1,ndimMTO)], mask=(ib_tableM(:)==ib2))
      do ibt1 = 1, size(ib_tableI)
        ib1 = ib_tableI(ibt1)
        idims = pack([(i,i=1,ndimMTO)], mask=(ib_tableM(:)==ib1))
        do it = 1, npair(ib1,ib2)
          phases(it) = 1d0/dble(nqwgt(it,ib1,ib2))*exp(-img*2d0*pi*sum(q*matmul(plat,nlat(:,it,ib1,ib2))))
        enddo
        do js = 1, nspinor
          jsp_src = merge(isp, js, .not. socmatrix)
          do jj = 1, size(jdims)
            istat = zmv(hammr(1,idims(1),jdims(jj),jsp_src), phases, hamm(idims(1),js,jdims(jj),js), n=size(idims), m=npair(ib1,ib2), opA=m_op_T, lda=npairmx)
            istat = zmv(ovlmr(1,idims(1),jdims(jj),jsp_src), phases, ovlm(idims(1),js,jdims(jj),js), n=size(idims), m=npair(ib1,ib2), opA=m_op_T, lda=npairmx)
          enddo
        enddo
        if(socmatrix) then
          do io = 1, 3
            do jj = 1, size(jdims)
              istat = zmv(hammhsor(1,idims(1),jdims(jj),io), phases, hammhso(idims(1),jdims(jj),io), n=size(idims), m=npair(ib1,ib2), opA=m_op_T, lda=npairmx)
            enddo
          enddo
        endif
      enddo
    enddo FourierTransform

    if(socmatrix) then !nspinor == 2
      hamm(:,1,:,1) = hamm(:,1,:,1) + hammhso(:,:,1)
      hamm(:,2,:,2) = hamm(:,2,:,2) + hammhso(:,:,2)
      hamm(:,1,:,2) = hammhso(:,:,3)
      hamm(:,2,:,1) = dconjg(transpose(hammhso(:,:,3)))
      if(present(ovlp_evec) .or. present(dual_evec)) then
        allocate(ovlm_buf(2*ndimMTO, 2*ndimMTO))
        ovlm_buf = reshape(ovlm, [2*ndimMTO, 2*ndimMTO])
      endif
      istat = zhgv(hamm, ovlm, n=2*ndimMTO, evl=ev)
      if(present(evec) .or. present(ovlp_evec) .or. present(dual_evec)) then
        hamm2d(1:2*ndimMTO, 1:2*ndimMTO) => hamm
        if(present(evec))      evec = hamm2d
        if(present(ovlp_evec)) istat = zmm(ovlm_buf, hamm2d, ovlp_evec, m=2*ndimMTO, n=2*ndimMTO, k=2*ndimMTO)
        if(present(dual_evec)) then
          dual_evec = hamm2d
          istat = zsv(ovlm_buf, dual_evec, n=2*ndimMTO, nrhs=2*ndimMTO)
        endif
      endif
    else !nspinor == 1
      if(present(ovlp_evec) .or. present(dual_evec)) then
        allocate(ovlm_buf(ndimMTO, ndimMTO))
        ovlm_buf = ovlm(:,1,:,1)
      endif
      istat = zhgv(hamm, ovlm, n=ndimMTO, evl=ev)
      if(present(evec))      evec(:,:) = hamm(:,1,:,1)
      if(present(ovlp_evec)) istat = zmm(ovlm_buf, hamm, ovlp_evec, m=ndimMTO, n=ndimMTO, k=ndimMTO)
      if(present(dual_evec)) then
        dual_evec(:,:) = hamm(:,1,:,1)
        istat = zsv(ovlm_buf, dual_evec, n=ndimMTO, nrhs=ndimMTO)
      endif
    endif
  end subroutine calc_ham_eigen
end module m_mlo_ham
