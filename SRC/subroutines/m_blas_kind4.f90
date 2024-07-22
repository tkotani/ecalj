submodule(m_blas) m_blas_kind4
  !$use omp_lib
  use cublas_v2, m_type =>CUDA_C_32F, compute_type => CUBLAS_COMPUTE_32F_FAST_TF32, algo => cublas_gemm_default
  implicit none
  integer, parameter :: kp = 4
contains
  module function cmm(a, b, c, m, n, k, opa, opb, alpha, beta, lda, ldb, ldc) result(istat)
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
    cublas_gemm: block 
      integer :: opa_in_cublas, opb_in_cublas
      istat = cublas_init()
      opa_in_cublas = get_m_op_cublas(opa_in)
      opb_in_cublas = get_m_op_cublas(opb_in)
      ! istat = cublascgemm3m(cublas_handle, opa_in_cublas, opb_in_cublas,  m, n, k, &
      !                     & alpha_in, a, lda_in , b, ldb_in, beta_in, c, ldc_in)
      istat = cublasGemmEX(cublas_handle, opa_in_cublas, opb_in_cublas, m, n, k,  &
                           alpha_in, a, m_type, lda_in, b, m_type, ldb_in, beta_in, c, m_type, ldc_in,&
                           compute_type, algo)
    endblock cublas_gemm
#else
    call cgemm3m(opa_in, opb_in, m, n, k, alpha_in, a, lda_in, b, ldb_in, beta_in, c, ldc_in)
    istat = 0
#endif
  end function cmm

  module function cmm_batch(a, b, c, m, n, k, nbatch, opa, opb, alpha, beta, lda, ldb, ldc, samea, sameb, comm) result(istat)
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
      ! istat = cublascgemmstridedbatched(cublas_handle, opa_in_cublas, opb_in_cublas,  m, n, k,  &
      !            &  alpha_in, a, lda_in, stridea, b, ldb_in, strideb, beta_in, c, ldc_in, stridec, nbatch)
      istat = cublasGemmStridedBatchedEX(cublas_handle, opa_in_cublas, opb_in_cublas, m, n, k, &
                                         alpha_in, a, m_type, lda_in, stridea, b, m_type, ldb_in, strideb, beta_in, &
                                         c, m_type, ldc_in, stridec, nbatch, compute_type, algo)
    endblock cublas_gemmstridedbatched
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
          call cgemm3m(opa_in, opb_in, m, n, k, alpha_in, a(stridea*(i-1)+1), lda_in, &
                     & b(strideb*(i-1)+1), ldb_in, beta_in, c(stridec*(i-1)+1), ldc_in)
        enddo
        do irank = 0, mpi_size - 1 
          call int_split(nbatch, mpi_size, irank, ini_batch, end_batch, nbatch_irank)
          call mpi_bcast(c(stridec*(ini_batch-1)+1), stridec*nbatch_irank, mpi_complex8, irank, comm, ierr)
        enddo
      else
        do i = 1, nbatch
          call cgemm3m(opa_in, opb_in, m, n, k, alpha_in, a(stridea*(i-1)+1), lda_in, &
                     & b(strideb*(i-1)+1), ldb_in, beta_in, c(stridec*(i-1)+1), ldc_in)
        enddo
      endif
      istat = nbatch
    end block blas_gemmbatch
#endif
  end function cmm_batch
end submodule m_blas_kind4
