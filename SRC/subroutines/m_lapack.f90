module m_lapack
#ifdef __GPU
  use cublas_v2
  use cudafor
  use cusolverdn
#endif
  use m_blas
  implicit none
  include "mpif.h"
  integer :: cuda_runtime_version
  public :: zminv, zminv_batch
  public :: zhgv
#ifdef __GPU
  type(cusolverDnHandle), value :: cusolver_handle 
  type(cusolverDnParams), value :: cusolver_params
  logical, save :: set_cusolver_init = .false.
#endif
  private
  integer, parameter :: kp = 8
contains
#ifdef __GPU
  integer function cusolver_init() result(istat)
    istat = 0
    if(.not.set_cusolver_init) then 
      istat = cusolverDnCreate(cusolver_handle)
      if(istat /= CUSOLVER_STATUS_SUCCESS) then
        print *, 'Error in cusolverDnCreate'
      endif
      istat = cusolverDnCreateParams(cusolver_params)
      istat = cudaRuntimeGetversion(cuda_runtime_version)
      set_cusolver_init = .true.
    endif
  end function cusolver_init
#endif
  !cusolverDnXtrtri used in zminv has internal compiler bug before cuda 12.5
  !https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/index.html
  function zminv(a, n, lda) result(istat)
    implicit none
    complex(8) :: a(*)
    integer, intent(in) :: n
    integer, optional :: lda
    integer :: lda_in, istat
#ifdef __GPU
    attributes(device) :: a
#endif
    lda_in = n; if(present(lda)) lda_in = lda
#ifdef __GPU
    GPU: block
      complex(8), device :: awork(lda_in*n)
      integer(4), device :: devinfo
      integer(8), device :: ipiv(n)
      integer(8)         :: ipiv_cpu(n)
      integer(1), allocatable, device :: buffer_d(:)
      integer(1), allocatable         :: buffer_h(:)
      integer(8) :: n_8, lda_8, lbuffer_d, lbuffer_h
      integer(4) :: i, j

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
    end block GPU
#else
     CPU: block
       integer, allocatable :: ipvt(:)
       complex(8),allocatable:: work(:)
       integer :: lwork
       lwork = 64*n
       allocate(work(lwork))
       allocate(ipvt(n))
       call zgetrf(n,n,a,lda_in,ipvt,istat)
       call zgetri(n,  a,lda_in,ipvt,work,lwork,istat)
       deallocate(work,ipvt)
     end block CPU
#endif
  end function zminv

  function zminv_batch(a, n, lda, nbatch) result(istat)
    use, intrinsic :: iso_c_binding
    implicit none
    complex(kind=kp) :: a(*)
    integer, intent(in) :: n
    integer, optional :: lda, nbatch
    integer :: lda_in, nbatch_in, ibatch, istat
    integer, allocatable :: ipvt(:,:), info(:)
#ifdef __GPU
    attributes(device) :: a, ipvt, info
#endif
    lda_in = n; nbatch_in = 1
    if(present(lda)) lda_in = lda
    if(present(nbatch)) nbatch_in = nbatch
    allocate (ipvt(n,nbatch_in), info(nbatch_in))
#ifdef __GPU
    GPU: block
      complex(kind=kp), allocatable, device :: c(:)
      type(c_devptr), allocatable, device :: a_dptr(:), c_dptr(:)
      integer :: array_size
      array_size = lda_in*lda_in*nbatch_in
      allocate (c(array_size))
      allocate (a_dptr(nbatch_in), c_dptr(nbatch_in))
      istat = cudadevicesynchronize() 
      do ibatch = 1, nbatch_in
        a_dptr(ibatch) = c_devloc(a(lda_in*lda_in*(ibatch-1)+1))
        c_dptr(ibatch) = c_devloc(c(lda_in*lda_in*(ibatch-1)+1))
      enddo
      istat = cublas_init()
      istat = cublasZgetrfBatched(cublas_handle, n, a_dptr, lda_in, ipvt, info, nbatch_in)
      istat = cublasZgetriBatched(cublas_handle, n, a_dptr, lda_in, ipvt, c_dptr, lda_in, info, nbatch_in)
      !$acc kernels
         a(1:array_size) = c(1:array_size)
      !$acc end kernels
      deallocate(a_dptr, c_dptr, c)
    end block GPU
#else
    CPU: block
      complex(kind=kp),allocatable:: work(:)
      integer :: lwork 
      lwork = 64*n
      allocate(work(lwork))
      do ibatch = 1, nbatch_in
        call zgetrf(n,n,a(lda_in*lda_in*(ibatch-1)+1),lda_in,ipvt(1,ibatch),info(ibatch))
        call zgetri(n,  a(lda_in*lda_in*(ibatch-1)+1),lda_in,ipvt(1,ibatch),work,lwork,info(ibatch))
      enddo
      istat = info(1)
    end block CPU
#endif
    deallocate(ipvt, info)
  end function zminv_batch
  function zhgv(A, B, n, evl, il, iu, lda, ldb) result(istat)
  ! Solving the generalized eigenvalue problem Az = lambda Bz, where A, B are Hermitian matrixes, z is eigenfunction
  ! Eigenvalues are stores in evl, eigenvectors are stored in A
  !!! that range is 1<=IL <= IU <= N
    implicit none
    integer, intent(in) :: n !size of matrix
    real(8), intent(out) :: evl(n) !eigenvalues
    complex(8) :: A(*), B(*)
    integer, intent(in), optional :: lda, ldb, il, iu
    integer :: lda_in, ldb_in, il_in, iu_in, istat
#ifdef __GPU
    attributes(device) :: a, b, evl
#endif

    lda_in = n; ldb_in = n
    if(present(lda)) lda_in = lda
    if(present(ldb)) ldb_in = ldb
    il_in = 1; iu_in = n
    if(present(il)) il_in = il
    if(present(iu)) iu_in = iu

#ifdef __GPU
    GPU: block
      real(8):: vu = 0d0, vl = 0d0
      integer, device :: devInfo
      complex(8), allocatable, device :: work(:)
      integer :: m, lwork
      istat = cusolver_init()
      istat = cusolverDnZhegvdx_bufferSize(cusolver_handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, &
                                            CUSOLVER_EIG_RANGE_I,  CUBLAS_FILL_MODE_UPPER, &
                                            n, a, lda_in, b, ldb_in, vl, vu, il_in, iu_in, m, evl, lwork)
      allocate(work(lwork))
      istat = cusolverDnZhegvdx(cusolver_handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, &
                                            CUSOLVER_EIG_RANGE_I,  CUBLAS_FILL_MODE_UPPER, &
                                            n, a, lda_in, b, ldb_in, vl, vu, il_in, iu_in, m, evl, work, lwork, devInfo)
      deallocate(work)
    endblock GPU
#else
    CPU: block
  !U: Upper traingal of A and B are stored
  !V: compute eigenvalues and eigenvectors
  !I: the IL-th through IU-th eigenvalues will be found
      real(8) :: dlamch, abstol, vl = 0d0, vu = 0d0
      integer, parameter :: nb = 64
      complex(8), allocatable:: work(:), z(:)
      complex(8) :: dummy(1)
      integer, allocatable:: ifail(:), iwork(:)
      real(8), allocatable:: rwork(:)
      integer :: lwork, info
      integer :: m !number of eigenvalues found
      abstol = 2d0*dlamch('S')
      allocate(z(lda_in*n), source = (0d0, 0d0))
      allocate(work(1))
      allocate(rwork(7*n), ifail(n), iwork(5*n))
      lwork = -1
      call zhegvx( 1, 'V', 'I', 'U', n, a, lda_in, b, ldb_in, &
         vl, vu, il_in, iu_in, abstol, m, evl, z, lda_in, &
         work, lwork, rwork, iwork, ifail, info )
      lwork = max((nb+1)*n, nint(real(dummy(1))))
      deallocate(work)
      allocate(work(lwork))
      call zhegvx( 1, 'V', 'I', 'U', n, a, lda_in, b, ldb_in, &
         vl, vu, il_in, iu_in, abstol, m, evl, z, lda_in, &
         work, lwork, rwork, iwork, ifail, info)
      istat = info
      a(1:lda_in*n) = z(1:lda_in*n)
      deallocate(work,rwork,iwork,ifail)
    endblock CPU
#endif
  end function
end module m_lapack
