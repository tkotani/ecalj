!>Accumulating rxcq
subroutine x0gemm(rcxq, npr, ipr_col, npr_col, nwhis, npm, ns1, ns2) 
  use m_kind, only: kp => kindrcxq
  use m_mpi, only: comm_b, mpi__rank_b, mpi__size_b
  use m_x0kf, only: icounkmink, icounkmaxk, iwini, iwend, itc, itpc, jpmc, icouini, whwc
  use m_blas, only: m_op_c, m_op_n, m_op_t
#ifdef __MP
  use m_blas, only: gemm => cmm
#else
  use m_blas, only: gemm => zmm
#endif
  use m_zmel, only: zmel
  use m_lgunit, only: stdo
  use m_ftox
#ifdef __GPU
  use openacc
  use cudafor
#endif
  use, intrinsic :: ieee_arithmetic
  !$ use omp_lib
  implicit none
  integer, intent(in) :: npr, ipr_col, npr_col, nwhis, npm, ns1, ns2
  complex(kind=kp), intent(inout) :: rcxq(npr,npr_col,nwhis,npm) !accumulating to rcxq
  integer :: icoun, igb1, igb2, iw, jpm, it, itp, ittp, nttp_max, ierr
  integer, allocatable :: nttp(:,:),  itw(:,:,:), itpw(:,:,:)
  complex(kind=kp), allocatable :: zw(:,:), wzw(:,:)
  complex(kind=kp), parameter :: CONE = (1_kp, 0_kp)
  real(8), allocatable :: whw(:,:,:)
  logical :: debug = .false.
#ifdef __GPU
  attributes(device) :: zw, wzw
#endif
  allocate(nttp(nwhis,npm), source = 0)
  do icoun = icounkmink, icounkmaxk
    jpm = jpmc(icoun)
    do iw = iwini(icoun), iwend(icoun)
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
    do iw = iwini(icoun), iwend(icoun)
      nttp(iw,jpm) = nttp(iw,jpm) + 1
      ittp = nttp(iw,jpm)
      itw(ittp,iw,jpm) = it
      itpw(ittp,iw,jpm) = itp
      whw(ittp,iw,jpm) = whwc(iw-iwini(icoun)+icouini(icoun))
    enddo
  enddo

  allocate(zw(nttp_max,npr), wzw(nttp_max,npr_col))
  !$acc host_data use_device(rcxq)
  !$acc data copyin(whw, itw, itpw, zmel)
  do jpm = 1, npm
    do iw = 1, nwhis
      if (nttp(iw,jpm) < 1) cycle
      !$acc kernels loop independent collapse(2)
      do ittp = 1, nttp(iw,jpm)
        do igb1 = 1, npr
          it  = itw(ittp,iw,jpm); itp = itpw(ittp,iw,jpm)
          zw(ittp,igb1) = cmplx(zmel(igb1,it,itp),kind=kp)
        enddo
      enddo
      !$acc end kernels
      !$acc kernels loop independent collapse(2)
      do ittp = 1, nttp(iw,jpm)
        do igb2 = 1, npr_col
          it  = itw(ittp,iw,jpm); itp = itpw(ittp,iw,jpm)
          wzw(ittp,igb2) = cmplx(zmel(igb2+ipr_col-1,it,itp)*whw(ittp,iw,jpm),kind=kp)
        enddo
      enddo
      !$acc end kernels
      ierr = gemm(zw, wzw, rcxq(1,1,iw,jpm), npr, npr_col, nttp(iw,jpm), &
              &  opA = m_op_C, beta = CONE , ldA = nttp_max, ldB = nttp_max)
    enddo
  enddo
  !$acc end data
  !$acc end host_data
  deallocate(itw, itpw, whw, wzw, zw, nttp)

end subroutine x0gemm
