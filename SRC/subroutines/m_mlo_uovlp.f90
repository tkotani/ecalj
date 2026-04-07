module m_mlo_uovlp !u overlap
! use m_blas, only:  zmv => zmv_h, m_op_T  ! uovlpt disabled
  use m_mlo_ham, only: nwf => ndimMTO, ib_tableM, ib_tableI, nsite
  use m_mlo_scrw, only: nnwf, nnwf_mask
  use m_lgunit, only: stdo
  use m_ftox, only: ftox
  use m_mpi, only: ipr
  implicit none
! public :: read_uovlpt, calc_formfactor, ng0, g0list  ! uovlpt is theoretically incorrect; disabled
  public :: read_formfactor, get_formfactor, q0i_q, nq0i_q
  private
! --- uovlpt: theoretically incorrect; disabled ---
! complex(8), allocatable :: uovlpt(:,:,:)  ! (npairmx, nwf, nwf): one G at a time
! real(8) :: plat(3,3), qlat(3,3)
! real(8), allocatable :: g0list(:,:)  ! (3, ng0): G vectors in qbz units; g0list(:,1)=0 (G=0)
! integer, allocatable :: npair(:,:), nqwgt(:,:,:), nlat(:,:,:,:)
! integer :: npairmx, ng0, nbas
! integer :: ifile_uovlp_dat = -1  ! direct-access file unit for uovlpt data
! integer :: ig0_loaded = 0        ! which ig0 is currently in uovlpt (0 = none loaded)
! -------------------------------------------------
  complex(8), allocatable :: formfactor_q(:,:)    ! (nwf, nwf): one q at a time
  real(8), allocatable :: q0i_q(:,:)        ! (3, nq0i_q): q-points from file
  integer :: nq0i_q = 0
  integer :: ifile_formfactor_dat = -1          ! direct-access file unit for formfactor_q data
  integer :: iq_loaded = 0                  ! which iq is currently in formfactor_q (0 = none)
contains
! --- uovlpt: theoretically incorrect; disabled ---
! subroutine read_uovlpt(isp, spin_flip)
!   integer, intent(in) :: isp
!   logical, intent(in), optional :: spin_flip
!   integer :: ihdr
!   character(len=20) :: datfile
!   ! Close any previously open data file from a prior spin call
!   if(ifile_uovlp_dat >= 0) close(ifile_uovlp_dat)
!   ig0_loaded = 0
!   open(newunit=ihdr, file='__UovlpTG.info', form='unformatted')
!   read(ihdr) npairmx, nbas, ng0
!   allocate(npair(nbas,nbas), nlat(3,npairmx,nbas,nbas), nqwgt(npairmx,nbas,nbas))
!   read(ihdr) plat(:,:), qlat(:,:), npair(:,:), nlat(:,:,:,:), nqwgt(:,:,:)
!   allocate(g0list(3,ng0), source=0d0)
!   read(ihdr) g0list(:,:)
!   close(ihdr)
!   ! Spin-dependent stream data file
!   if(isp==1 .and. spin_flip) then
!     datfile = '__UovlpTG.UPDN'
!   elseif(isp==2 .and. spin_flip) then
!     datfile = '__UovlpTG.DNUP'
!   elseif(isp==1) then
!     datfile = '__UovlpTG.UP'   ! not implemented
!   else
!     datfile = '__UovlpTG.DN'   ! not implemented
!   endif
!   open(newunit=ifile_uovlp_dat, file=trim(datfile), form='unformatted', access='direct', recl=npairmx*nwf*nwf*16)
!   if(.not. allocated(uovlpt)) allocate(uovlpt(npairmx,nwf,nwf))
! end subroutine
!
! subroutine load_uovlpt_ig0(Gvec_in)
!   real(8), intent(in) :: Gvec_in(3)
!   integer :: ig0
!   do ig0 = 1, ng0
!     if(all(abs(g0list(:,ig0) - Gvec_in) < 1d-8)) then
!       if(ig0 == ig0_loaded) return
!       write(stdo,ftox) 'read ig0:', ig0, g0list(:,ig0)
!       read(ifile_uovlp_dat, rec=ig0) uovlpt(:,:,:)
!       ig0_loaded = ig0
!       return
!     endif
!   enddo
!   call rx('Gvec not found in load_uovlpt_ig0')
! end subroutine
!
! subroutine calc_formfactor(qg, formfactor)
!   real(8), intent(in) :: qg(3) !q+G
!   real(8) :: qvec(3), Gvec(3), dist, min_dist
!   integer :: ibt1, ibt2, ib1, ib2, it, i, j, jj, istat, ig0
!   integer, allocatable :: jdims(:), idims(:)
!   complex(8), intent(out) :: formfactor(nnwf)
!   complex(8) :: formfactor_mat(nwf,nwf)
!   complex(8), allocatable :: phases(:)
!   complex(8), parameter :: img=(0d0,1d0)
!   real(8), parameter ::pi=4d0*atan(1d0)
!   allocate(phases(npairmx))
!   ! Find G in g0list that minimizes |q| = |Q - G|
!   min_dist = huge(1d0)
!   do ig0 = 1, ng0
!     dist = sum((qg - g0list(:,ig0))**2)
!     if(dist < min_dist) then
!       min_dist = dist
!       Gvec = g0list(:,ig0)
!     endif
!   enddo
!   qvec = qg - Gvec
!     ! write(stdo,ftox) 'q+G', qg
!     ! write(stdo,ftox) 'gvec', Gvec
!     ! write(stdo,ftox) 'qvec', qvec
!     ! write(stdo,ftox) 'plat', plat
!     ! write(stdo,ftox) 'qlat', qlat
!   call load_uovlpt_ig0(Gvec)
!   do ibt2 = 1, nsite
!     ib2 = ib_tableI(ibt2)
!     jdims = pack([(j,j=1,nwf)], mask=(ib_tableM(:)==ib2))
!     do ibt1 = 1, nsite
!       ib1 = ib_tableI(ibt1)
!       idims = pack([(i,i=1,nwf)], mask=(ib_tableM(:)==ib1))
!       do it =1, npair(ib1,ib2)
!         ! Phase: exp(-i*q*R); G*R=0 for integer lattice R so G does not enter the phase
!         phases(it) = 1d0/dble(nqwgt(it,ib1,ib2))*exp(-img*2d0*pi* sum(qvec*matmul(plat,nlat(:,it,ib1,ib2))))
!       enddo
!       do jj = 1, size(jdims)
!         istat = zmv(uovlpt(1,idims(1),jdims(jj)), phases, formfactor_mat(idims(1),jdims(jj)), &
!                     n=size(idims), m=npair(ib1,ib2), opA=m_op_T, lda=npairmx)
!       enddo
!     enddo
!   enddo
!   write(stdo,ftox)  'TrTuovlpxuvolp',ig0, qg, sqrt(sum(qg**2)), sum(abs(formfactor_mat(:,:))**2)
!   formfactor(1:nnwf) = pack(reshape(formfactor_mat, shape=[nwf*nwf]), mask = nnwf_mask)
! end subroutine
! -------------------------------------------------

  subroutine read_formfactor(isp, spin_flip)
    integer, intent(in) :: isp
    logical, intent(in) :: spin_flip
    integer :: ihdr, nwf_file, nqbz_file, nspin_file
    character(20) :: datfile
    if(ifile_formfactor_dat >= 0) close(ifile_formfactor_dat)
    iq_loaded = 0
    ! Read info file (sequential): header + q-point list
    open(newunit=ihdr, file='__MLOFormFactor.info', form='unformatted')
    read(ihdr) nwf_file, nqbz_file, nspin_file, nq0i_q
    if(allocated(q0i_q))  deallocate(q0i_q)
    if(allocated(formfactor_q)) deallocate(formfactor_q)
    allocate(q0i_q(3, nq0i_q))
    allocate(formfactor_q(nwf_file, nwf_file))
    read(ihdr) q0i_q(:,:)
    close(ihdr)
    ! Open spin-dependent data file (direct-access): one record per q-point
    if(isp==1 .and. (.not.spin_flip)) datfile = '__MLOFormFactor.UP'
    if(isp==2 .and. (.not.spin_flip)) datfile = '__MLOFormFactor.DN'
    if(isp==1 .and. spin_flip) datfile = '__MLOFormFactor.UPDN'
    if(isp==2 .and. spin_flip) datfile = '__MLOFormFactor.DNUP'
    open(newunit=ifile_formfactor_dat, file=trim(datfile), form='unformatted', access='direct', recl=nwf_file*nwf_file*16)
    write(stdo,ftox) 'read formfactor_q: isp, nwf, nq0i_q =', isp, nwf_file, nq0i_q
  end subroutine

  subroutine get_formfactor(q, formfactor)
    real(8), intent(in) :: q(3)
    complex(8), intent(out) :: formfactor(nnwf)
    integer :: iq
    do iq = 1, nq0i_q
      if(all(abs(q0i_q(:,iq) - q) < 1d-8)) then
        if(iq /= iq_loaded) then
          read(ifile_formfactor_dat, rec=iq) formfactor_q(:,:)
          iq_loaded = iq
        endif
        formfactor(1:nnwf) = pack(reshape(formfactor_q(:,:), shape=[nwf*nwf]), mask = nnwf_mask)
        return
      endif
    enddo
    call rx('q not found in get_formfactor')
  end subroutine
end module
