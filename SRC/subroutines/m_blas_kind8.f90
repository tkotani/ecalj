submodule(m_blas) m_blas_kind8
  !$use omp_lib
  implicit none
  integer, parameter :: kp = 8
contains
  module function zmm(a, b, c, m, n, k, opa, opb, alpha, beta, lda, ldb, ldc) result(istat)
    complex(kind=kp) :: a(*), b(*), c(*)
    integer, intent(in) :: m, n, k
    character, intent(in), optional :: opa, opb
    complex(kind=kp), intent(in), optional :: alpha, beta
    integer, optional :: lda, ldb, ldc

    complex(kind=kp) :: alpha_in, beta_in
    integer :: lda_in, ldb_in, ldc_in, istat
    character :: opa_in, opb_in
#ifdef __GPU
    attributes(device) :: a, b, c
#endif
    if (m < 1 .or. n < 1 .or. k < 1) return

    alpha_in = (1_kp, 0_kp); beta_in = (0_kp, 0_kp)
    if(present(alpha)) alpha_in = alpha
    if(present(beta)) beta_in = beta

    opa_in = m_op_n; opb_in = m_op_n
    if(present(opa)) opa_in = opa
    if(present(opb)) opb_in = opb

    !opa(a) = m x k, opb(b) = k x n, c = m x n
    lda_in = m; ldb_in = k; ldc_in = m
    if(opa_in == m_op_t .or. opa_in == m_op_c) lda_in = k !a = k x m
    if(opb_in == m_op_t .or. opb_in == m_op_c) ldb_in = n !b = n x k

    if(present(lda)) lda_in = lda
    if(present(ldb)) ldb_in = ldb
    if(present(ldc)) ldc_in = ldc

#ifdef __GPU
    cublas_gemm3m: block 
      integer :: opa_in_cublas, opb_in_cublas
      istat = cublas_init()
      opa_in_cublas = get_m_op_cublas(opa_in)
      opb_in_cublas = get_m_op_cublas(opb_in)
      istat = cublaszgemm3m(cublas_handle, opa_in_cublas, opb_in_cublas,  m, n, k, &
                          & alpha_in, a, lda_in , b, ldb_in, beta_in, c, ldc_in)
    end block cublas_gemm3m
#else
    call zgemm3m(opa_in, opb_in, m, n, k, alpha_in, a, lda_in, b, ldb_in, beta_in, c, ldc_in)
    istat = 0
#endif
  end function zmm

  module function zmm_batch(a, b, c, m, n, k, nbatch, opa, opb, alpha, beta, lda, ldb, ldc, samea, sameb, comm) result(istat)
    complex(kind=kp) :: a(*), b(*), c(*)
    integer, intent(in) :: m, n, k, nbatch
    character, intent(in), optional :: opa, opb
    complex(kind=kp), intent(in), optional :: alpha, beta
    integer, optional :: lda, ldb, ldc
    logical, optional :: samea, sameb
    integer, optional :: comm

    integer(8) :: stridea, strideb, stridec
    complex(kind=kp) :: alpha_in, beta_in
    integer :: lda_in, ldb_in, ldc_in, istat
    character :: opa_in, opb_in
#ifdef __GPU
    attributes(device) :: a, b, c
#endif

    if (nbatch < 1) return
    if (m < 1 .or. n < 1 .or. k < 1) return

    alpha_in = (1_kp, 0_kp); beta_in = (0_kp, 0_kp)
    if(present(alpha)) alpha_in = alpha
    if(present(beta)) beta_in = beta

    opa_in = m_op_n; opb_in = m_op_n
    if(present(opa)) opa_in = opa
    if(present(opb)) opb_in = opb

    !opa(a) = m x k, opb(b) = k x n, c = m x n
    lda_in = m; ldb_in = k; ldc_in = m
    if(opa_in == m_op_t .or. opa_in == m_op_c) lda_in = k !a = k x m
    if(opb_in == m_op_t .or. opb_in == m_op_c) ldb_in = n !b = n x k

    if(present(lda)) lda_in = lda
    if(present(ldb)) ldb_in = ldb
    if(present(ldc)) ldc_in = ldc

    stridea = int(lda_in*k,8); strideb = int(ldb_in*n,8); stridec = int(ldc_in*n,8)
    if (opa_in == m_op_t .or. opa_in == m_op_c) stridea = int(lda_in*m,8)
    if (opb_in == m_op_t .or. opb_in == m_op_c) strideb = int(ldb_in*k,8)

    if(present(samea)) then
      if(samea) stridea = 0_8
    endif
    if(present(sameb)) then
      if(sameb) strideb = 0_8
    endif

#ifdef __GPU
    cublas_gemmstridedbatched: block
      integer :: opa_in_cublas, opb_in_cublas
      istat = cublas_init()
      opa_in_cublas = get_m_op_cublas(opa_in)
      opb_in_cublas = get_m_op_cublas(opb_in)
      istat = cublaszgemmstridedbatched(cublas_handle, opa_in_cublas, opb_in_cublas,  m, n, k,  &
                 &  alpha_in, a, lda_in, stridea, b, ldb_in, strideb, beta_in, c, ldc_in, stridec, nbatch)
    end block cublas_gemmstridedbatched
#else
    blas_gemmbatch: block
      integer :: i
      integer :: ini_batch, end_batch, mpi_size, mpi_rank, ierr, nbatch_irank, irank
      ! MPI parallelization version, but, it is not efficient for small matrix
      if(present(comm)) then
        call mpi_comm_size(comm, mpi_size, ierr)
        call mpi_comm_rank(comm, mpi_rank, ierr)
        call int_split(nbatch, mpi_size, mpi_rank, ini_batch, end_batch, nbatch_irank)
        do i = ini_batch, end_batch
          call zgemm3m(opa_in, opb_in, m, n, k, alpha_in, a(stridea*(i-1)+1), lda_in, &
                     & b(strideb*(i-1)+1), ldb_in, beta_in, c(stridec*(i-1)+1), ldc_in)
        enddo
        do irank = 0, mpi_size - 1 
          call int_split(nbatch, mpi_size, irank, ini_batch, end_batch, nbatch_irank)
          call mpi_bcast(c(stridec*(ini_batch-1)+1), stridec*nbatch_irank, mpi_complex16, irank, comm, ierr)
        enddo
      else
        do i = 1, nbatch
          call zgemm3m(opa_in, opb_in, m, n, k, alpha_in, a(stridea*(i-1)+1), lda_in, &
                     & b(strideb*(i-1)+1), ldb_in, beta_in, c(stridec*(i-1)+1), ldc_in)
        enddo
      endif
      istat = nbatch
    end block blas_gemmbatch
#endif
  end function zmm_batch

  !cusolverDnXtrtri used in zminv has internal compiler bug before cuda 12.5
  !https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/index.html
  module function zminv(a, n, lda) result(istat)
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

  module function zminv_batch(a, n, lda, nbatch) result(istat)
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
    end block CPU
#endif
    deallocate(ipvt, info)
  end function zminv_batch
end submodule m_blas_kind8
