module m_mlo_ham
  use m_HamPMT, only: plat, npair, nlat, nqwgt,nspc
  use m_lgunit, only: stdo
  use m_mpi, only: ipr
  use m_blas, only: zmm => zmm_h, zmv => zmv_h, m_op_T
  use m_lapack, only: zhgv => zhgv_h, zsv => zsv_h
  use m_ftox, only: ftox
  implicit none
  public :: read_ham_rs, calc_ham_eigen
  integer, protected, target :: ndimMTO, npairmx, nspx, nsite
  integer, allocatable, protected:: ib_tableM(:), l_tableM(:), k_tableM(:), ib_tableI(:)
  complex(8),allocatable, protected:: ovlmr(:,:,:,:), hammr(:,:,:,:),hammhsor(:,:,:,:) !npairmx, ndimMTO, ndimMTO, nspx order
  logical, protected :: socmatrix = .false.
contains
  subroutine read_ham_rs()! read RealSpace MTO Hamiltonian
    integer:: ifihmto, i
    logical:: cmdopt0
    socmatrix=cmdopt0('--socmatrix')
    open(newunit=ifihmto,file='HamRsMLO',form='unformatted')
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

  subroutine calc_ham_eigen(nspi,nspe, q, ev, evec, ovlp_evec, dual_evec)
    real(8), intent(in) :: q(3)
    integer, intent(in) :: nspi, nspe
    real(8), intent(out) :: ev(:,:) !MLO eigenvalue
    complex(8), optional, intent(out) :: evec(:,:), ovlp_evec(:,:), dual_evec(:,:) !MLO wavefunction
    complex(8), allocatable :: ovlm(:,:,:,:), hamm(:,:,:,:)
    complex(8) :: ovlm_buf(ndimMTO,ndimMTO)
    complex(8),allocatable:: hammhso(:,:,:)
    real(8), parameter :: oveps=1d-15, pi=4d0*atan(1d0)
    complex(8), parameter :: img=(0d0,1d0)
    complex(8) :: phases(npairmx)
    integer, allocatable :: idims(:), jdims(:)
    integer ::i, j, nev, ib1, ib2, it, ii, jj, ibt1, ibt2, istat, io, isp, ns, sp
    ! complex(8) :: phase
    ! ovlm = 0d0
    ! hamm = 0d0
    ! FourierTransormationFROMrealspcaeTOqspace: do i=1,ndimMTO !MLO Hamiltonian
    !   do   j=1,ndimMTO
    !     ib1 = ib_tableM(i) !atomic-site index in the primitive cell
    !     ib2 = ib_tableM(j)
    !     do it =1,npair(ib1,ib2)
    !       phase=1d0/dble(nqwgt(it,ib1,ib2))*exp(-img*2d0*pi* sum(q*matmul(plat,nlat(:,it,ib1,ib2))))
    !       ! hamm(i,j)= hamm(i,j)+ hammr(i,j,it,isp)*phase !MLO Hamiltonian at qp
    !       ! ovlm(i,j)= ovlm(i,j)+ ovlmr(i,j,it,isp)*phase
    !       hamm(i,j)= hamm(i,j)+ hammr(it,i,j,isp)*phase !MLO Hamiltonian at qp
    !       ovlm(i,j)= ovlm(i,j)+ ovlmr(it,i,j,isp)*phase
    !     enddo
    !   enddo
    ! enddo FourierTransormationFROMrealspcaeTOqspace
    ns = merge(2, 1, socmatrix)
    allocate(ovlm(ndimMTO,ns,ndimMTO,ns), hamm(ndimMTO,ns,ndimMTO,ns))
    if(socmatrix) allocate(hammhso(ndimMTO,ndimMTO,3) ,source=(0d0,0d0))
    ovlm=0d0
    hamm=0d0
    do isp=nspi,nspe
      sp = merge(isp, 1, socmatrix) !spin block index into hamm/ovlm
      if(.not.socmatrix) then !reset per-spin: each isp is diagonalized independently
        hamm(:,1,:,1)=0d0; ovlm(:,1,:,1)=0d0
      endif
      FourierTransormationFROMrealspcaeTOqspace:do ibt2 = 1, size(ib_tableI)
        ib2 = ib_tableI(ibt2)
        jdims = pack([(j,j=1,ndimMTO)], mask=(ib_tableM(:)==ib2))
        do ibt1 = 1, size(ib_tableI)
          ib1 = ib_tableI(ibt1)
          idims = pack([(i,i=1,ndimMTO)], mask=(ib_tableM(:)==ib1))
          do it =1, npair(ib1,ib2)
            phases(it) = 1d0/dble(nqwgt(it,ib1,ib2))*exp(-img*2d0*pi* sum(q*matmul(plat,nlat(:,it,ib1,ib2))))
          enddo
          do jj = 1, size(jdims)
            istat = zmv(hammr(1,idims(1),jdims(jj),isp), phases, hamm(idims(1),sp,jdims(jj),sp), n=size(idims), m=npair(ib1,ib2), opA=m_op_T, lda=npairmx)
            istat = zmv(ovlmr(1,idims(1),jdims(jj),isp), phases, ovlm(idims(1),sp,jdims(jj),sp), n=size(idims), m=npair(ib1,ib2), opA=m_op_T, lda=npairmx)
          enddo
          if(socmatrix.and.isp==nspx) then
          do jj = 1, size(jdims)
          do io=1,3
            istat = zmv(hammhsor(1,idims(1),jdims(jj),io), phases,hammhso(idims(1),jdims(jj),io),n=size(idims),m=npair(ib1,ib2),opA=m_op_T, lda=npairmx)
          enddo
          enddo
          endif
        enddo
      enddo FourierTransormationFROMrealspcaeTOqspace
      if(present(ovlp_evec) .or. present(dual_evec)) ovlm_buf(:,:) = ovlm(:,1,:,1) !keep ovlm
      if(socmatrix.and.isp==nspx) then
        hamm(:,1,:,1)= hamm(:,1,:,1)+ hammhso(:,:,1)
        hamm(:,2,:,2)= hamm(:,2,:,2)+ hammhso(:,:,2)
        hamm(:,1,:,2)= hammhso(:,:,3) !see m_bandcal
        hamm(:,2,:,1)= dconjg(transpose(hammhso(:,:,3)))
        istat = zhgv(hamm, ovlm, n=2*ndimMTO, evl=ev(:,isp-nspi+1)) !in-place hamm -> evec
        if(present(evec))      call rx('present(evec) not yet')
        if(present(ovlp_evec)) call rx('present(ovlp_evec) not yet')
        if(present(dual_evec)) call rx('present(dual_evec) not yet')
        return
      elseif(socmatrix) then
        cycle !accumulate; diagonalize at isp==nspx
      else
        istat = zhgv(hamm, ovlm, n=ndimMTO, evl=ev(:,isp-nspi+1)) !in-place hamm -> evec
      endif
      if(present(evec)) evec(:,:) = hamm(:,1,:,1)
      if(present(ovlp_evec)) istat = zmm(ovlm_buf, hamm, ovlp_evec, m=ndimMTO, n=ndimMTO, k=ndimMTO)
      if(present(dual_evec)) then
        dual_evec(:,:) = hamm(:,1,:,1)
        istat = zsv(ovlm_buf, dual_evec, n=ndimMTO, nrhs=ndimMTO)
      endif
    enddo
  end subroutine calc_ham_eigen
end module m_mlo_ham
