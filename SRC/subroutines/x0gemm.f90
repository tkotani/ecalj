!>Accumulating rxcq
subroutine x0gemm(rcxq, npr, nwhis, npm, ns1, ns2, iw_lo_in, iw_hi_in)
  use m_kind, only: kp => kindrcxq
  use m_mpi, only: comm_b => comm_b_xq, mpi__rank_b => mpi__rank_b_xq, mpi__size_b => mpi__size_b_xq
  use m_x0kf, only: icounkmink, icounkmaxk, iwini, iwend, itc, itpc, jpmc, icouini, whwc
  use m_blas, only: m_op_c, m_op_n, m_op_t
  use m_zmel, only: zmel
  use m_lgunit, only: stdo
  use m_ftox
#if defined(__MP) && defined(__GPU)
  use m_blas, only: gemm => cmm_d
#elif defined(__MP)
  use m_blas, only: gemm => cmm_h
#elif defined(__GPU)
  use m_blas, only: gemm => zmm_d
#else
  use m_blas, only: gemm => zmm_h
#endif
#ifdef __GPU
  use openacc
  use cudafor
#endif
  use, intrinsic :: ieee_arithmetic
  !$ use omp_lib
  implicit none
  integer, intent(in) :: npr, nwhis, npm, ns1, ns2
  integer, intent(in) :: iw_lo_in, iw_hi_in
  complex(kind=kp), intent(inout) :: rcxq(npr,npr,iw_lo_in:iw_hi_in)
  integer :: icoun, igb1, igb2, iw, jpm, iw_pos, it, itp, ittp, nttp_max, ierr
  integer :: iw_lo, iw_hi
  integer :: pos_lo(2), pos_hi(2)
  integer, allocatable :: nttp(:,:),  itw(:,:,:), itpw(:,:,:)
  complex(kind=kp), allocatable :: zw(:,:), wzw(:,:)
  complex(kind=kp), parameter :: CONE = (1_kp, 0_kp)
  real(8), allocatable :: whw(:,:,:)
  logical :: debug = .false.
#ifdef __GPU
  attributes(device) :: rcxq, zw, wzw
#endif
  iw_lo = iw_lo_in
  iw_hi = iw_hi_in

  ! Owned positive-iw_pos ranges derived from flat iw_lo:iw_hi.
  ! jpm=1: flat iw = +iw_pos → owned when iw_lo <= iw_pos <= iw_hi (positive part)
  ! jpm=2: flat iw = -iw_pos → owned when iw_lo <= -iw_pos <= iw_hi (negative part)
  pos_lo(1) = max(iw_lo, 1);      pos_hi(1) = min(iw_hi, nwhis)
  pos_lo(2) = max(1, -iw_hi);     pos_hi(2) = min(nwhis, -iw_lo)

  allocate(nttp(nwhis,npm), source = 0)
  do icoun = icounkmink, icounkmaxk
    jpm = jpmc(icoun)
    do iw = max(iwini(icoun), pos_lo(jpm)), min(iwend(icoun), pos_hi(jpm))
      nttp(iw,jpm) = nttp(iw,jpm) + 1
    enddo
  enddo

  nttp_max = maxval(nttp(1:nwhis,1:npm))
  if(debug) write(stdo, ftox)'nttp_max = ', nttp_max
  allocate (itw(nttp_max,nwhis,npm), source = 0)
  allocate (itpw(nttp_max,nwhis,npm), source = 0)
  allocate (whw(nttp_max,nwhis,npm), source = 0d0)

  nttp(1:nwhis,1:npm) = 0
  do icoun = icounkmink, icounkmaxk
    jpm = jpmc(icoun)
    it  = itc (icoun)
    itp = itpc(icoun)
    if(it  < ns1 .or. it > ns2) cycle
    do iw = max(iwini(icoun), pos_lo(jpm)), min(iwend(icoun), pos_hi(jpm))
      nttp(iw,jpm) = nttp(iw,jpm) + 1
      ittp = nttp(iw,jpm)
      itw(ittp,iw,jpm) = it
      itpw(ittp,iw,jpm) = itp
      whw(ittp,iw,jpm) = whwc(iw-iwini(icoun)+icouini(icoun))
    enddo
  enddo

  allocate(zw(nttp_max,npr), wzw(nttp_max,npr))
  !$acc data copyin(whw, itw, itpw, zmel)
  do iw = iw_lo, iw_hi
    if (iw == 0) cycle
    if (iw > 0) then; jpm = 1; iw_pos = iw
    else;             jpm = 2; iw_pos = -iw
    endif
    if (nttp(iw_pos,jpm) < 1) cycle
    !$acc kernels loop independent collapse(2)
    do ittp = 1, nttp(iw_pos,jpm)
      do igb1 = 1, npr
        it  = itw(ittp,iw_pos,jpm); itp = itpw(ittp,iw_pos,jpm)
        zw(ittp,igb1) = cmplx(zmel(igb1,it,itp),kind=kp)
      enddo
    enddo
    !$acc end kernels
    !$acc kernels loop independent collapse(2)
    do igb2 = 1, npr
      do ittp = 1, nttp(iw_pos,jpm)
        wzw(ittp,igb2) = cmplx(zw(ittp,igb2)*whw(ittp,iw_pos,jpm),kind=kp)
      enddo
    enddo
    !$acc end kernels
    ierr = gemm(zw, wzw, rcxq(1,1,iw), npr, npr, nttp(iw_pos,jpm), &
            &  opA = m_op_C, beta = CONE, ldA = nttp_max, ldB = nttp_max)
  enddo
  !$acc end data

  deallocate(itw, itpw, whw, wzw, zw, nttp)

end subroutine x0gemm
