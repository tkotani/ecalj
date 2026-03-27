module m_mlo_uovlp !u overlap
  use m_blas, only:  zmv => zmv_h, m_op_T
  use m_mlo_ham, only: nwf => ndimMTO, ib_tableM, ib_tableI, nsite
  use m_mlo_scrw, only: nnwf, nnwf_mask
  use m_lgunit, only: stdo
  use m_ftox, only: ftox
  use m_mpi, only: ipr
  implicit none
  public :: read_uovlpt, calc_uovlpq
  private
  complex(8), allocatable :: uovlpt(:,:,:,:)
  real(8) :: plat(3,3)
  integer, allocatable :: g0list(:,:), npair(:,:), nqwgt(:,:,:), nlat(:,:,:,:)
  integer :: npairmx, ng0, nbas
contains
  subroutine read_uovlpt(isp, spin_flip)
    integer, intent(in) :: isp
    logical, intent(in), optional :: spin_flip
    logical :: spin_flip_in
    integer:: ifile
    if(isp==1 .and. spin_flip) open(newunit=ifile,file='__UOVLPT.UPDN',form='unformatted')
    if(isp==2 .and. spin_flip) open(newunit=ifile,file='__UOVLPT.DNUP',form='unformatted')
    if(isp==1 .and. (.not. spin_flip)) open(newunit=ifile,file='__UOVLPT.UP',form='unformatted') !not implemented
    if(isp==2 .and. (.not. spin_flip)) open(newunit=ifile,file='__UOVLPT.DN',form='unformatted') !not implemented
    read(ifile) npairmx, nbas, ng0
    allocate(npair(nbas,nbas), nlat(3,npairmx,nbas,nbas), nqwgt(npairmx,nbas,nbas))
    allocate(g0list(3,ng0), source = 0)
    read(ifile) plat(:,:), npair(:,:), nlat(:,:,:,:), nqwgt(:,:,:)
    allocate(uovlpt(npairmx,nwf,nwf,ng0))
    read(ifile) uovlpt(:,:,:,:)
    close(ifile)
  end subroutine
  subroutine calc_uovlpq(q, uovlpq, g0)
    real(8), intent(in) :: q(3)
    integer, intent(in), optional :: g0(3)
    integer :: ig0, ibt1, ibt2, ib1, ib2, it, i, j, jj, istat
    integer, allocatable :: jdims(:), idims(:)
    complex(8), intent(out) :: uovlpq(nnwf)
    complex(8) :: uovlpq_mat(nwf,nwf)
    complex(8), allocatable :: phases(:)
    complex(8), parameter :: img=(0d0,1d0)
    real(8), parameter ::pi=4d0*atan(1d0)
    logical, allocatable :: mask(:)
    allocate(phases(npairmx))
    if(present(g0)) then
      if(any(g0(1:3) /= 0)) call rx('non zero g0 is not implemented')
    endif
    ig0 = 1
    do ibt2 = 1, nsite
      ib2 = ib_tableI(ibt2)
      jdims = pack([(j,j=1,nwf)], mask=(ib_tableM(:)==ib2))
      do ibt1 = 1, nsite
        ib1 = ib_tableI(ibt1)
        idims = pack([(i,i=1,nwf)], mask=(ib_tableM(:)==ib1))
        do it =1, npair(ib1,ib2)
          phases(it) = 1d0/dble(nqwgt(it,ib1,ib2))*exp(-img*2d0*pi* sum(q*matmul(plat,nlat(:,it,ib1,ib2))))
        enddo
        do jj = 1, size(jdims)
          istat = zmv(uovlpt(1,idims(1),jdims(jj),ig0), phases, uovlpq_mat(idims(1),jdims(jj)), &
                      n=size(idims), m=npair(ib1,ib2), opA=m_op_T, lda=npairmx)
        enddo
      enddo
    enddo
    uovlpq(1:nnwf) = pack(reshape(uovlpq_mat, shape=[nwf*nwf]), mask = nnwf_mask)
  end subroutine
end module
