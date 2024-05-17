!>Accumulating rxcq
subroutine x0gpu(rcxq, npr, ipr_col, npr_col, nwhis, npm) 
  use m_kind,only: kindrcxq
  use m_mpi, only: comm_b, mpi__rank_b, mpi__size_b
  use m_x0kf,only: icounkmink, icounkmaxk, iwini, iwend, itc, itpc, jpmc, icouini, whwc
  use m_blas, only: m_op_c, m_op_n, m_op_t, zmm, zmm_batch !CPU or GPU versions specifed by macro __GPU
  use m_zmel, only: zmel
#ifdef __GPU
  use openacc
  use cudafor
#endif
  !$ use omp_lib
  implicit none
  integer, intent(in) :: npr, ipr_col, npr_col, nwhis, npm
  complex(kindrcxq), intent(inout) :: rcxq(npr,npr_col,nwhis,npm) !accumulating to rcxq
  integer :: icoun, igb1, igb2, iw, jpm, it, itp, ittp, nttp_max, ierr
  integer, allocatable :: nttp(:,:),  itw(:,:,:), itpw(:,:,:)
  complex(8), allocatable :: zw(:,:), wzw(:,:)
  real(8), allocatable :: whw(:,:,:)
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
  allocate (itw(nttp_max,nwhis,npm), source = 0)
  allocate (itpw(nttp_max,nwhis,npm), source = 0)
  allocate (whw(nttp_max,nwhis,npm), source = 0d0)

  nttp(1:nwhis,1:npm) = 0
  do icoun = icounkmink, icounkmaxk
    jpm = jpmc(icoun)
    it  = itc (icoun)
    itp = itpc(icoun)
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
    !!$omp parallel do private(iw, ittp, igb1, it, itp, zw, wzw) shared(rcxq) schedule(dynamic,10) 
    do iw = 1, nwhis
      if (nttp(iw,jpm) < 1) cycle
      !$acc kernels loop independent collapse(2)
      do ittp = 1, nttp(iw,jpm)
        do igb1 = 1, npr
          it  = itw(ittp,iw,jpm); itp = itpw(ittp,iw,jpm)
          zw(ittp,igb1) = zmel(igb1,it,itp)
        enddo
      enddo
      !$acc end kernels

      !$acc kernels loop independent collapse(2)
      do ittp = 1, nttp(iw,jpm)
        do igb2 = 1, npr_col
          it  = itw(ittp,iw,jpm); itp = itpw(ittp,iw,jpm)
          wzw(ittp,igb2) = zmel(igb2+ipr_col-1,it,itp)*whw(ittp,iw,jpm)
        enddo
      enddo
      !$acc end kernels

      ierr = zmm(zw, wzw, rcxq(1,1,iw,jpm), npr, npr_col, nttp(iw,jpm), &
              &  opA = m_op_C, beta = (1d0,0d0), ldA = nttp_max, ldB = nttp_max)
    enddo
    !!$omp end parallel do
  enddo
  !$acc end data
  !$acc end host_data
  deallocate(itw, itpw, whw, wzw, zw, nttp)

end subroutine x0gpu
