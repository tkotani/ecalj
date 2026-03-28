module m_lapack
#ifdef __GPU
  use cublas_v2
  use cudafor
  use cusolverdn
#endif
  use m_blas
  implicit none
  public :: zhgv_h, zhev_h, zminv_h, zsv_h, zgev_h
#ifdef __GPU
  public :: zhgv_d, zhev_d, zminv_d, zsv_d, zgev_d, cusolver_finalize
#endif
  private
#ifdef __GPU
  integer :: cuda_runtime_version
  type(cusolverDnHandle), value :: cusolver_handle
  type(cusolverDnParams), value :: cusolver_params
  logical, save :: set_cusolver_handle = .false.
#endif
contains
  integer function zminv_h(a, n, lda) result(istat)
    complex(8) :: a(*)
    integer, intent(in) :: n
    integer, optional :: lda
    integer :: lda_in
    integer, allocatable :: ipvt(:)
    complex(8),allocatable:: work(:)
    integer :: lwork, original_num_threads
    complex(8) :: wkopt
    lda_in = n; if(present(lda)) lda_in = lda
#ifdef __MKL_ZMINV_SEQUENTIAL
    original_num_threads = mkl_get_max_threads()
    call mkl_set_num_threads(1)
#endif
    allocate(ipvt(n))
    call zgetrf(n,n,a,lda_in,ipvt,istat)
    lwork = -1
    call zgetri(n, a, lda_in, ipvt, wkopt, lwork, istat)
    lwork = int(dble(wkopt))
    allocate(work(lwork))
    call zgetri(n, a, lda_in, ipvt, work, lwork, istat)
#ifdef __MKL_ZMINV_SEQUENTIAL
    call mkl_set_num_threads(original_num_threads)
#endif
    deallocate(work,ipvt)
  end function zminv_h
  integer function zsv_h(a, b, n, nrhs, lda, ldb) result(istat)
    integer, intent(in) :: n, nrhs
    complex(8) :: a(*), b(*)
    integer, allocatable :: ipiv(:)
    integer, intent(in), optional :: lda, ldb
    integer :: lda_in, ldb_in
    allocate(ipiv(n))
    lda_in = n; ldb_in = n
    if(present(lda)) lda_in = lda
    if(present(ldb)) ldb_in = ldb
    call zgetrf(n, n, a, lda_in, ipiv, istat)
    if(istat /= 0) return
    call zgetrs('N', n, nrhs, a, lda_in, ipiv, b, ldb_in, istat)
    deallocate(ipiv)
  end function zsv_h
  integer function zgev_h(a, n, evl, evl_vec, evr_vec, lda) result(istat)
    integer, intent(in) :: n
    complex(8) :: a(*)
    complex(8), intent(out) :: evl(n)
    complex(8), intent(out), optional :: evl_vec(*), evr_vec(*)
    integer, intent(in), optional :: lda
    integer :: lda_in, ldvl, ldvr, lwork, ilo, ihi
    complex(8), allocatable :: work(:), vl_loc(:), vr_loc(:)
    real(8),    allocatable :: rwork(:), scale(:), rconde(:), rcondv(:)
    real(8) :: abnrm
    character :: jobvl, jobvr
    lda_in = n; if(present(lda)) lda_in = lda
    jobvl = 'N'; ldvl = 1
    jobvr = 'N'; ldvr = 1
    if(present(evl_vec)) then; jobvl = 'V'; ldvl = lda_in; endif
    if(present(evr_vec)) then; jobvr = 'V'; ldvr = lda_in; endif
    allocate(vl_loc(ldvl*n), vr_loc(ldvr*n), rwork(2*n), work(1), scale(n), rconde(n), rcondv(n))
    lwork = -1
    call zgeevx('B', jobvl, jobvr, 'N', n, a, lda_in, evl, &
                vl_loc, ldvl, vr_loc, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, &
                work, lwork, rwork, istat)
    lwork = int(dble(work(1))); deallocate(work); allocate(work(lwork))
    call zgeevx('B', jobvl, jobvr, 'N', n, a, lda_in, evl, &
                vl_loc, ldvl, vr_loc, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, &
                work, lwork, rwork, istat)
    if(present(evl_vec)) evl_vec(1:ldvl*n) = vl_loc
    if(present(evr_vec)) evr_vec(1:ldvr*n) = vr_loc
    deallocate(work, rwork, vl_loc, vr_loc, scale, rconde, rcondv)
  end function zgev_h
  integer function zhev_h(A, n, evl, il, iu, lda) result(istat)
  ! Solving the standard eigenvalue problem Az = lambda z, where A is a Hermitian matrix
  ! Eigenvalues are stored in evl, eigenvectors are stored in A
  !!! range: 1<=il<=iu<=n
    integer, intent(in) :: n
    real(8), intent(out) :: evl(n)
    complex(8) :: A(*)
    integer, intent(in), optional :: lda, il, iu
    integer :: lda_in, il_in, iu_in
    complex(8), allocatable :: work(:), z(:)
    real(8), allocatable :: rwork(:)
    integer, allocatable :: isuppz(:), iwork(:)
    integer :: m, lwork, lrwork, liwork, info
    real(8) :: abstol, vl, vu, dlamch
    lda_in = n; if(present(lda)) lda_in = lda
    il_in = 1; iu_in = n
    if(present(il)) il_in = il
    if(present(iu)) iu_in = iu
    vl = 0d0; vu = 0d0; abstol = 2d0*dlamch('S')
    allocate(z(lda_in*n), isuppz(2*n))
    lwork = -1; lrwork = -1; liwork = -1
    allocate(work(1), rwork(1), iwork(1))
    call zheevr('V', 'I', 'U', n, a, lda_in, vl, vu, il_in, iu_in, abstol, m, evl, z, lda_in, isuppz, &
                work, lwork, rwork, lrwork, iwork, liwork, info)
    lwork = int(dble(work(1))); lrwork = int(rwork(1)); liwork = iwork(1)
    deallocate(work, rwork, iwork)
    allocate(work(lwork), rwork(lrwork), iwork(liwork))
    call zheevr('V', 'I', 'U', n, a, lda_in, vl, vu, il_in, iu_in, abstol, m, evl, z, lda_in, isuppz, &
                work, lwork, rwork, lrwork, iwork, liwork, info)
    istat = info
    a(1:lda_in*n) = z(1:lda_in*n)
    deallocate(work, rwork, iwork, z, isuppz)
  end function zhev_h
  integer function zhgv_h(A, B, n, evl, il, iu, lda, ldb) result(istat)
  ! Solving the generalized eigenvalue problem Az = lambda Bz, where A, B are Hermitian matrixes, z is eigenfunction
  ! Eigenvalues are stores in evl, eigenvectors are stored in A
  !!! that range is 1<=IL <= IU <= N
    integer, intent(in) :: n !size of matrix
    real(8), intent(out) :: evl(n) !eigenvalues
    complex(8) :: A(*), B(*)
    integer, intent(in), optional :: lda, ldb, il, iu
    integer :: lda_in, ldb_in, il_in, iu_in
    integer, parameter :: nb = 64
    complex(8), allocatable:: work(:), z(:)
    complex(8) :: dummy(1)
    integer, allocatable:: ifail(:), iwork(:)
    real(8), allocatable:: rwork(:)
    integer :: m, lwork, info
    real(8) :: dlamch, abstol, vl = 0d0, vu = 0d0
    lda_in = n; ldb_in = n
    if(present(lda)) lda_in = lda
    if(present(ldb)) ldb_in = ldb
    il_in = 1; iu_in = n
    if(present(il)) il_in = il
    if(present(iu)) iu_in = iu
    abstol = 2d0*dlamch('S')
    allocate(z(lda_in*n), source = (0d0, 0d0))
    allocate(work(1))
    allocate(rwork(7*n), ifail(n), iwork(5*n))
    lwork = -1
    call zhegvx( 1, 'V', 'I', 'U', n, a, lda_in, b, ldb_in, &
       vl, vu, il_in, iu_in, abstol, m, evl, z, lda_in, &
       work, lwork, rwork, iwork, ifail, info )
    lwork = max((nb+1)*n, nint(dble(dummy(1))))
    deallocate(work)
    allocate(work(lwork))
    call zhegvx( 1, 'V', 'I', 'U', n, a, lda_in, b, ldb_in, &
       vl, vu, il_in, iu_in, abstol, m, evl, z, lda_in, &
       work, lwork, rwork, iwork, ifail, info)
    istat = info
    a(1:lda_in*n) = z(1:lda_in*n)
    deallocate(work,rwork,iwork,ifail)
  end function zhgv_h
#ifdef __GPU
  !cusolverDnXtrtri used in zminv has internal compiler bug before cuda 12.5
  !https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/index.html
  integer function zminv_d(a, n, lda) result(istat)
    complex(8), device :: a(*)
    integer, intent(in) :: n
    integer, optional :: lda
    integer :: lda_in
    complex(8), allocatable, device :: awork(:)
    integer(4), device :: devinfo
    integer(8), device :: ipiv(n)
    integer(8)         :: ipiv_cpu(n)
    integer(1), allocatable, device :: buffer_d(:)
    integer(1), allocatable         :: buffer_h(:)
    integer(8) :: n_8, lda_8, lbuffer_d, lbuffer_h
    integer(4) :: i, j
    lda_in = n; if(present(lda)) lda_in = lda

    n_8 = int(n,8)
    lda_8 = int(lda_in,8)
    istat = cublas_init()
    istat = cusolver_init()
    ! LU factorization PA = LU in-place on a, L: Lower triangular with unit diagonal components, U: Upper triangular with
    ! non-unit-diagonal,  P: permutation of rows stored in ipiv
    ! A^-1 = U^-1 * L^-1 * P
    istat = cusolverDnXgetrf_buffersize(cusolver_handle, cusolver_params, n_8, n_8, cudaDataType(CUDA_C_64F), a, lda_8, &
                                        cudaDataType(CUDA_C_64F), lbuffer_d, lbuffer_h)
    allocate(buffer_d(lbuffer_d), buffer_h(lbuffer_h))
    istat = cusolverDnXgetrf(cusolver_handle, cusolver_params, n_8, n_8, cudaDataType(CUDA_C_64F), a, lda_8, &
                             ipiv, cudaDataType(CUDA_C_64F), buffer_d, lbuffer_d, buffer_h, lbuffer_h, devinfo)
    deallocate(buffer_d, buffer_h)
    allocate(awork(lda_in*n))
    !$acc kernels
    awork(1:lda_in*n) = a(1:lda_in*n) !awork -> Used as a upper triangular matrix with diagonal part
    !$acc end kernels
    !$acc kernels loop collapse(2) independent
    do j = 1, n
      do i = 1, n
        if(i < j) a(lda_in*(j-1)+i) = (0d0, 0d0) !Used as lower triangular matrix, diagonal components will be replaced as 1 in the following function
      enddo
    enddo
    !$acc end kernels
    !Get inverse of L = L^-1
    istat = cusolverDnXtrtri_buffersize(cusolver_handle, CUBLAS_FILL_MODE_LOWER, CUBLAS_DIAG_UNIT, n_8,  &
                                        cudaDataType(CUDA_C_64F), a, lda_8, lbuffer_d, lbuffer_h)
    if(cuda_runtime_version <= 12040)  then !prescription for bug of cusolverDnXtrtri_buffersize
      lbuffer_d = lbuffer_d*16; lbuffer_h = lbuffer_h*16
    endif
    allocate(buffer_d(lbuffer_d), buffer_h(lbuffer_h))
    istat = cusolverDnXtrtri(cusolver_handle, CUBLAS_FILL_MODE_LOWER, CUBLAS_DIAG_UNIT, n_8, cudaDataType(CUDA_C_64F), a, lda_8, &
                             buffer_d, lbuffer_d, buffer_h, lbuffer_h, devinfo)
    deallocate(buffer_d, buffer_h)
    !Get inverse of U = U^-1
    istat = cusolverDnXtrtri_buffersize(cusolver_handle, CUBLAS_FILL_MODE_UPPER, CUBLAS_DIAG_NON_UNIT, n_8, &
                                        cudaDataType(CUDA_C_64F), awork, lda_8, lbuffer_d, lbuffer_h)
    if(cuda_runtime_version <= 12040)  then !prescription for bug of cusolverDnXtrtri_buffersize
      lbuffer_d = lbuffer_d*16; lbuffer_h = lbuffer_h*16
    endif
    allocate(buffer_d(lbuffer_d), buffer_h(lbuffer_h))
    istat = cusolverDnXtrtri(cusolver_handle, CUBLAS_FILL_MODE_UPPER, CUBLAS_DIAG_NON_UNIT, n_8, cudaDataType(CUDA_C_64F), &
                             awork, lda_8, buffer_d, lbuffer_d, buffer_h, lbuffer_h, devinfo)
    deallocate(buffer_d, buffer_h)
    ! Get U^-1 * L^-1 = awork * a
    istat = cublasZtrmm(cublas_handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, n, n, &
                        (1d0, 0d0), awork, lda_in, a, lda_in, a, lda_in) !result in a (in-place)
    ipiv_cpu = ipiv !copy to CPU from GPU
    do i = n, 1, -1
      if(ipiv_cpu(i) > i) istat = cublasZswap(cublas_handle, n, a(lda_in*(i-1)+1), 1, a(lda_in*(ipiv_cpu(i)-1)+1), 1)
    enddo
    deallocate(awork)
  end function zminv_d
  integer function zsv_d(a, b, n, nrhs, lda, ldb) result(istat)
    integer, intent(in) :: n, nrhs
    integer, intent(in), optional :: lda, ldb
    complex(8), device :: a(*), b(*)
    integer(8), device :: ipiv(n)
    integer(8) :: n_8, lda_8, ldb_8, nrhs_8
    integer(1), allocatable, device :: buffer_d(:)
    integer(1), allocatable         :: buffer_h(:)
    integer(8) :: lbuffer_d, lbuffer_h
    integer(4), device :: devinfo
    integer :: lda_in, ldb_in
    istat = cusolver_init()
    lda_in = n; if(present(lda)) lda_in = lda
    ldb_in = n; if(present(ldb)) ldb_in = ldb

    n_8 = int(n,8); lda_8 = int(lda_in,8); ldb_8 = int(ldb_in,8); nrhs_8 = int(nrhs,8)
    istat = cusolverDnXgetrf_buffersize(cusolver_handle, cusolver_params, n_8, n_8, cudaDataType(CUDA_C_64F), a, lda_8, &
                                        cudaDataType(CUDA_C_64F), lbuffer_d, lbuffer_h)
    allocate(buffer_d(lbuffer_d), buffer_h(lbuffer_h))
    istat = cusolverDnXgetrf(cusolver_handle, cusolver_params, n_8, n_8, cudaDataType(CUDA_C_64F), a, lda_8, &
                             ipiv, cudaDataType(CUDA_C_64F), buffer_d, lbuffer_d, buffer_h, lbuffer_h, devinfo)
    deallocate(buffer_d, buffer_h)
    istat = cusolverDnXgetrs(cusolver_handle, cusolver_params, CUBLAS_OP_N, n_8, nrhs_8, cudaDataType(CUDA_C_64F), a, lda_8, &
                             ipiv, cudaDataType(CUDA_C_64F), b, ldb_8, devinfo)
  end function zsv_d
  integer function zgev_d(a, n, evl, evl_vec, evr_vec, lda) result(istat)
  ! CPU fallback: cusolverDnZgeev was removed from cuSOLVER (CUDA >= 13)
  ! Device arrays are copied to host, solved via LAPACK zgeevx, then copied back.
    integer, intent(in) :: n
    complex(8), device :: a(*)
    complex(8), device :: evl(n)
    complex(8), device, optional :: evl_vec(*), evr_vec(*)
    integer, intent(in), optional :: lda
    integer :: lda_in, ldvl, ldvr, lwork, ilo, ihi
    complex(8), allocatable :: a_h(:), evl_h(:), work(:), vl_loc(:), vr_loc(:)
    real(8),    allocatable :: rwork(:), scale_(:), rconde(:), rcondv(:)
    real(8) :: abnrm
    character :: jobvl, jobvr
    lda_in = n; if(present(lda)) lda_in = lda
    jobvl = 'N'; ldvl = 1
    jobvr = 'N'; ldvr = 1
    if(present(evl_vec)) then; jobvl = 'V'; ldvl = lda_in; endif
    if(present(evr_vec)) then; jobvr = 'V'; ldvr = lda_in; endif
    allocate(a_h(lda_in*n), evl_h(n))
    a_h(1:lda_in*n) = a(1:lda_in*n)  ! device -> host
    allocate(vl_loc(ldvl*n), vr_loc(ldvr*n), rwork(2*n), work(1), scale_(n), rconde(n), rcondv(n))
    lwork = -1
    call zgeevx('B', jobvl, jobvr, 'N', n, a_h, lda_in, evl_h, &
                vl_loc, ldvl, vr_loc, ldvr, ilo, ihi, scale_, abnrm, rconde, rcondv, &
                work, lwork, rwork, istat)
    lwork = int(dble(work(1))); deallocate(work); allocate(work(lwork))
    call zgeevx('B', jobvl, jobvr, 'N', n, a_h, lda_in, evl_h, &
                vl_loc, ldvl, vr_loc, ldvr, ilo, ihi, scale_, abnrm, rconde, rcondv, &
                work, lwork, rwork, istat)
    a(1:lda_in*n) = a_h(1:lda_in*n)  ! host -> device
    evl(1:n) = evl_h(1:n)            ! host -> device
    if(present(evl_vec)) evl_vec(1:ldvl*n) = vl_loc
    if(present(evr_vec)) evr_vec(1:ldvr*n) = vr_loc
    deallocate(a_h, evl_h, work, rwork, vl_loc, vr_loc, scale_, rconde, rcondv)
  end function zgev_d
  integer function zhev_d(A, n, evl, il, iu, lda) result(istat) !Not tested
  ! Solving the standard eigenvalue problem Az = lambda z, where A is a Hermitian matrix (GPU version)
  ! Eigenvalues are stored in evl, eigenvectors are stored in A
  !!! range: 1<=il<=iu<=n
    integer, intent(in) :: n
    real(8), intent(out), device :: evl(n)
    complex(8), device :: A(*)
    integer, intent(in), optional :: lda, il, iu
    integer :: lda_in, il_in, iu_in
    real(8) :: vl, vu
    integer, device :: devInfo
    complex(8), allocatable, device :: work(:)
    integer :: m, lwork
    lda_in = n; if(present(lda)) lda_in = lda
    il_in = 1; iu_in = n
    if(present(il)) il_in = il
    if(present(iu)) iu_in = iu
    vl = 0d0; vu = 0d0
    istat = cusolver_init()
    istat = cusolverDnZheevdx_bufferSize(cusolver_handle, CUSOLVER_EIG_MODE_VECTOR, &
                                          CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_UPPER, &
                                          n, a, lda_in, vl, vu, il_in, iu_in, m, evl, lwork)
    allocate(work(lwork))
    istat = cusolverDnZheevdx(cusolver_handle, CUSOLVER_EIG_MODE_VECTOR, &
                               CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_UPPER, &
                               n, a, lda_in, vl, vu, il_in, iu_in, m, evl, work, lwork, devInfo)
    deallocate(work)
  end function zhev_d
  integer function zhgv_d(A, B, n, evl, il, iu, lda, ldb) result(istat)
    integer, intent(in) :: n !size of matrix
    real(8), intent(out), device :: evl(n) !eigenvalues
    complex(8), device :: A(*), B(*)
    integer, intent(in), optional :: lda, ldb, il, iu
    integer :: lda_in, ldb_in, il_in, iu_in
    real(8):: vu = 0d0, vl = 0d0
    integer, device :: devInfo
    complex(8), allocatable, device :: work(:)
    integer :: m, lwork
    lda_in = n; ldb_in = n
    if(present(lda)) lda_in = lda
    if(present(ldb)) ldb_in = ldb
    il_in = 1; iu_in = n
    if(present(il)) il_in = il
    if(present(iu)) iu_in = iu
    istat = cusolver_init()
    istat = cusolverDnZhegvdx_bufferSize(cusolver_handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, &
                                          CUSOLVER_EIG_RANGE_I,  CUBLAS_FILL_MODE_UPPER, &
                                          n, a, lda_in, b, ldb_in, vl, vu, il_in, iu_in, m, evl, lwork)
    allocate(work(lwork))
    istat = cusolverDnZhegvdx(cusolver_handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, &
                                          CUSOLVER_EIG_RANGE_I,  CUBLAS_FILL_MODE_UPPER, &
                                          n, a, lda_in, b, ldb_in, vl, vu, il_in, iu_in, m, evl, work, lwork, devInfo)
    deallocate(work)
  end function zhgv_d
  integer function cusolver_init() result(istat)
    istat = 0
    if(.not.set_cusolver_handle) then 
      istat = cusolverDnCreate(cusolver_handle)
      if(istat /= CUSOLVER_STATUS_SUCCESS) then
        print *, 'Error in cusolverDnCreate'
      endif
      istat = cusolverDnCreateParams(cusolver_params)
      istat = cudaRuntimeGetversion(cuda_runtime_version)
      set_cusolver_handle = .true.
    endif
  end function cusolver_init
  integer function cusolver_finalize() result(istat)
    istat = 0
    if(set_cusolver_handle) then
      istat = cusolverDnDestroy(cusolver_handle)
      if(istat /= CUSOLVER_STATUS_SUCCESS) then
        print *, 'Error in cusolverDnDestroy'
      endif
      istat = cusolverDnDestroyParams(cusolver_params)
      if(istat /= CUSOLVER_STATUS_SUCCESS) then
        print *, 'Error in cusolverDnDestroyParams'
      endif
      set_cusolver_handle = .false.
    endif
  end function cusolver_finalize
#endif
end module m_lapack
