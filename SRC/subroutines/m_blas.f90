module m_blas !wrapper for BLAS and cuBLAS
  !$use omp_lib
#ifdef __GPU
  use cublas_v2
  use cudafor
#endif
  implicit none
  include "mpif.h"
  public :: int_split
  public :: m_op_n, m_op_t, m_op_c
  public :: cmm_h, cmm_batch_h, zmm_h, zmm_batch_h, dmm_h, dmv_h, zmv_h
#ifdef __GPU
  public :: cmm_d, cmm_batch_d, zmm_d, zmm_batch_d, dmm_d, dmv_d, zmv_d
  public :: cublas_init, cublas_handle
  type(cublashandle), value :: cublas_handle
  logical, save :: set_cublas_handle = .false.
#endif
  character, parameter :: m_op_n = 'N', m_op_t = 'T', m_op_c = 'C'
contains
  integer function cmm_h(a, b, c, m, n, k, opa, opb, alpha, beta, lda, ldb, ldc) result(istat)
    complex(4) :: a(*), b(*), c(*)
    integer, intent(in) :: m, n, k
    character, intent(in), optional :: opa, opb
    complex(4), intent(in), optional :: alpha, beta
    integer, optional :: lda, ldb, ldc
    complex(4) :: alpha_in, beta_in
    integer :: lda_in, ldb_in, ldc_in
    character :: opa_in, opb_in
    if (m < 1 .or. n < 1 .or. k < 1) return
    alpha_in = (1.0, 0.0); beta_in = (0.0, 0.0)
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
    call cgemm3m(opa_in, opb_in, m, n, k, alpha_in, a, lda_in, b, ldb_in, beta_in, c, ldc_in)
    istat = 0
  end function cmm_h
  integer function cmm_batch_h(a, b, c, m, n, k, nbatch, opa, opb, alpha, beta, lda, ldb, ldc, samea, sameb, comm) result(istat)
    complex(4) :: a(*), b(*), c(*)
    integer, intent(in) :: m, n, k, nbatch
    character, intent(in), optional :: opa, opb
    complex(4), intent(in), optional :: alpha, beta
    integer, optional :: lda, ldb, ldc
    logical, optional :: samea, sameb
    integer, optional :: comm
    integer(8) :: stridea, strideb, stridec
    complex(4) :: alpha_in, beta_in
    integer :: lda_in, ldb_in, ldc_in
    character :: opa_in, opb_in
    integer :: i
    integer :: ini_batch, end_batch, mpi_size, mpi_rank, ierr, nbatch_irank, irank
    if (nbatch < 1) return
    if (m < 1 .or. n < 1 .or. k < 1) return
    alpha_in = (1.0, 0.0); beta_in = (0.0, 0.0)
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
  end function cmm_batch_h
  integer function dmv_h(a, x, y, m, n, opa, alpha, beta, lda, incx, incy) result(istat)
    implicit none
    real(8) :: a(*), x(*), y(*)
    integer, intent(in) :: m, n
    character, intent(in), optional :: opa
    real(8), intent(in), optional :: alpha, beta
    integer, optional :: lda, incx, incy
    real(8) :: alpha_in, beta_in
    integer :: lda_in, incx_in, incy_in
    character :: opa_in
    if (m < 1 .or. n < 1) return
    alpha_in = 1d0; beta_in = 0d0
    if(present(alpha)) alpha_in = alpha
    if(present(beta)) beta_in = beta
    opa_in = m_op_n
    if(present(opa)) opa_in = opa
    lda_in = m; incx_in = 1; incy_in = 1
    if(present(lda)) lda_in = lda
    if(present(incx)) incx_in = incx
    if(present(incy)) incy_in = incy
    call dgemv(opa_in, m, n, alpha_in, a, lda_in, x, incx_in, beta_in, y, incy_in)
    istat = 0
  end function dmv_h
  integer function zmv_h(a, x, y, m, n, opa, alpha, beta, lda, incx, incy) result(istat)
    implicit none
    complex(8) :: a(*), x(*), y(*)
    integer, intent(in) :: m, n
    character, intent(in), optional :: opa
    complex(8), intent(in), optional :: alpha, beta
    integer, optional :: lda, incx, incy
    complex(8) :: alpha_in, beta_in
    integer :: lda_in, incx_in, incy_in
    character :: opa_in
    if (m < 1 .or. n < 1) return
    alpha_in = (1d0, 0d0); beta_in = (0d0, 0d0)
    if(present(alpha)) alpha_in = alpha
    if(present(beta)) beta_in = beta
    opa_in = m_op_n
    if(present(opa)) opa_in = opa
    lda_in = m; incx_in = 1; incy_in = 1
    if(present(lda)) lda_in = lda
    if(present(incx)) incx_in = incx
    if(present(incy)) incy_in = incy
    call zgemv(opa_in, m, n, alpha_in, a, lda_in, x, incx_in, beta_in, y, incy_in)
    istat = 0
  end function zmv_h
  integer function dmm_h(a, b, c, m, n, k, opa, opb, alpha, beta, lda, ldb, ldc) result(istat)
    real(8) :: a(*), b(*), c(*)
    integer, intent(in) :: m, n, k
    character, intent(in), optional :: opa, opb
    real(8), intent(in), optional :: alpha, beta
    integer, optional :: lda, ldb, ldc
    real(8) :: alpha_in, beta_in
    integer :: lda_in, ldb_in, ldc_in
    character :: opa_in, opb_in
    if (m < 1 .or. n < 1 .or. k < 1) return
    alpha_in = 1d0; beta_in = 0d0
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
    call dgemm(opa_in, opb_in, m, n, k, alpha_in, a, lda_in, b, ldb_in, beta_in, c, ldc_in)
    istat = 0
  end function dmm_h
  integer function zmm_h(a, b, c, m, n, k, opa, opb, alpha, beta, lda, ldb, ldc) result(istat)
    complex(8) :: a(*), b(*), c(*)
    integer, intent(in) :: m, n, k
    character, intent(in), optional :: opa, opb
    complex(8), intent(in), optional :: alpha, beta
    integer, optional :: lda, ldb, ldc
    complex(8) :: alpha_in, beta_in
    integer :: lda_in, ldb_in, ldc_in
    character :: opa_in, opb_in
    if (m < 1 .or. n < 1 .or. k < 1) return
    alpha_in = (1d0, 0d0); beta_in = (0d0, 0d0)
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
    call zgemm3m(opa_in, opb_in, m, n, k, alpha_in, a, lda_in, b, ldb_in, beta_in, c, ldc_in)
    istat = 0
  end function zmm_h
  integer function zmm_batch_h(a, b, c, m, n, k, nbatch, opa, opb, alpha, beta, lda, ldb, ldc, samea, sameb, comm) result(istat)
    complex(8) :: a(*), b(*), c(*)
    integer, intent(in) :: m, n, k, nbatch
    character, intent(in), optional :: opa, opb
    complex(8), intent(in), optional :: alpha, beta
    integer, optional :: lda, ldb, ldc
    logical, optional :: samea, sameb
    integer, optional :: comm
    integer(8) :: stridea, strideb, stridec
    complex(8) :: alpha_in, beta_in
    integer :: lda_in, ldb_in, ldc_in
    character :: opa_in, opb_in
    integer :: i
    integer :: ini_batch, end_batch, mpi_size, mpi_rank, ierr, nbatch_irank, irank
    if (nbatch < 1) return
    if (m < 1 .or. n < 1 .or. k < 1) return
    alpha_in = (1d0, 0d0); beta_in = (0d0, 0d0)
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
  end function zmm_batch_h
#ifdef __GPU
  integer function cmm_d(a, b, c, m, n, k, opa, opb, alpha, beta, lda, ldb, ldc) result(istat)
    use cublas_v2, m_type =>CUDA_C_32F, compute_type => CUBLAS_COMPUTE_32F_FAST_TF32, algo => cublas_gemm_default
    complex(4), device :: a(*), b(*), c(*)
    integer, intent(in) :: m, n, k
    character, intent(in), optional :: opa, opb
    complex(4), intent(in), optional :: alpha, beta
    integer, optional :: lda, ldb, ldc
    complex(4) :: alpha_in, beta_in
    integer :: lda_in, ldb_in, ldc_in
    character :: opa_in, opb_in
    integer :: opa_in_cublas, opb_in_cublas
    if (m < 1 .or. n < 1 .or. k < 1) return
    alpha_in = (1.0, 0.0); beta_in = (0.0, 0.0)
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
    istat = cublas_init()
    opa_in_cublas = get_m_op_cublas(opa_in)
    opb_in_cublas = get_m_op_cublas(opb_in)
    istat = cublasGemmEX(cublas_handle, opa_in_cublas, opb_in_cublas, m, n, k,  &
                         alpha_in, a, m_type, lda_in, b, m_type, ldb_in, beta_in, c, m_type, ldc_in,&
                         compute_type, algo)
  end function cmm_d
  integer function cmm_batch_d(a, b, c, m, n, k, nbatch, opa, opb, alpha, beta, lda, ldb, ldc, samea, sameb, comm) result(istat)
    use cublas_v2, m_type =>CUDA_C_32F, compute_type => CUBLAS_COMPUTE_32F_FAST_TF32, algo => cublas_gemm_default
    complex(4), device :: a(*), b(*), c(*)
    integer, intent(in) :: m, n, k, nbatch
    character, intent(in), optional :: opa, opb
    complex(4), intent(in), optional :: alpha, beta
    integer, optional :: lda, ldb, ldc
    logical, optional :: samea, sameb
    integer, optional :: comm
    integer(8) :: stridea, strideb, stridec
    complex(4) :: alpha_in, beta_in
    integer :: lda_in, ldb_in, ldc_in
    character :: opa_in, opb_in
    integer :: opa_in_cublas, opb_in_cublas
    if (nbatch < 1) return
    if (m < 1 .or. n < 1 .or. k < 1) return
    alpha_in = (1.0, 0.0); beta_in = (0.0, 0.0)
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
    istat = cublas_init()
    opa_in_cublas = get_m_op_cublas(opa_in)
    opb_in_cublas = get_m_op_cublas(opb_in)
    istat = cublasGemmStridedBatchedEX(cublas_handle, opa_in_cublas, opb_in_cublas, m, n, k, &
                                       alpha_in, a, m_type, lda_in, stridea, b, m_type, ldb_in, strideb, beta_in, &
                                       c, m_type, ldc_in, stridec, nbatch, compute_type, algo)
  end function cmm_batch_d
  integer function dmv_d(a, x, y, m, n, opa, alpha, beta, lda, incx, incy) result(istat)
    implicit none
    real(8), device :: a(*), x(*), y(*)
    integer, intent(in) :: m, n
    !caution: size of matrix a is m x n (not size of op(A))
    character, intent(in), optional :: opa
    real(8), intent(in), optional :: alpha, beta
    integer, optional :: lda, incx, incy
    real(8) :: alpha_in, beta_in
    integer :: lda_in, incx_in, incy_in
    character :: opa_in
    integer :: opa_in_cublas, opb_in_cublas
    if (m < 1 .or. n < 1) return
    alpha_in = 1d0; beta_in = 0d0
    if(present(alpha)) alpha_in = alpha
    if(present(beta)) beta_in = beta
    opa_in = m_op_n
    if(present(opa)) opa_in = opa
    lda_in = m; incx_in = 1; incy_in = 1
    if(present(lda)) lda_in = lda
    if(present(incx)) incx_in = incx
    if(present(incy)) incy_in = incy
    istat = cublas_init()
    opa_in_cublas = get_m_op_cublas(opa_in)
    istat = cublasdgemv(cublas_handle, opa_in_cublas,  m, n,  &
                      & alpha_in, a, lda_in , x, incx_in, beta_in, y, incy_in)
  end function dmv_d
  integer function zmv_d(a, x, y, m, n, opa, alpha, beta, lda, incx, incy) result(istat)
    implicit none
    complex(8), device :: a(*), x(*), y(*)
    integer, intent(in) :: m, n
    character, intent(in), optional :: opa
    complex(8), intent(in), optional :: alpha, beta
    integer, optional :: lda, incx, incy
    complex(8) :: alpha_in, beta_in
    integer :: lda_in, incx_in, incy_in
    character :: opa_in
    integer :: opa_in_cublas, opb_in_cublas
    if (m < 1 .or. n < 1) return
    alpha_in = (1d0, 0d0); beta_in = (0d0, 0d0)
    if(present(alpha)) alpha_in = alpha
    if(present(beta)) beta_in = beta
    opa_in = m_op_n
    if(present(opa)) opa_in = opa
    lda_in = m; incx_in = 1; incy_in = 1
    if(present(lda)) lda_in = lda
    if(present(incx)) incx_in = incx
    if(present(incy)) incy_in = incy
    istat = cublas_init()
    opa_in_cublas = get_m_op_cublas(opa_in)
    istat = cublaszgemv(cublas_handle, opa_in_cublas,  m, n,  &
                      & alpha_in, a, lda_in , x, incx_in, beta_in, y, incy_in)
  end function zmv_d
  integer function dmm_d(a, b, c, m, n, k, opa, opb, alpha, beta, lda, ldb, ldc) result(istat)
    real(8), device :: a(*), b(*), c(*)
    integer, intent(in) :: m, n, k
    character, intent(in), optional :: opa, opb
    real(8), intent(in), optional :: alpha, beta
    integer, optional :: lda, ldb, ldc

    real(8) :: alpha_in, beta_in
    integer :: lda_in, ldb_in, ldc_in
    character :: opa_in, opb_in
    integer :: opa_in_cublas, opb_in_cublas
    if (m < 1 .or. n < 1 .or. k < 1) return
    alpha_in = 1d0; beta_in = 0d0
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
    istat = cublas_init()
    opa_in_cublas = get_m_op_cublas(opa_in)
    opb_in_cublas = get_m_op_cublas(opb_in)
    istat = cublasdgemm(cublas_handle, opa_in_cublas, opb_in_cublas,  m, n, k, &
                      & alpha_in, a, lda_in , b, ldb_in, beta_in, c, ldc_in)
  end function dmm_d
  integer function zmm_d(a, b, c, m, n, k, opa, opb, alpha, beta, lda, ldb, ldc) result(istat)
    complex(8), device :: a(*), b(*), c(*)
    integer, intent(in) :: m, n, k
    character, intent(in), optional :: opa, opb
    complex(8), intent(in), optional :: alpha, beta
    integer, optional :: lda, ldb, ldc
    complex(8) :: alpha_in, beta_in
    integer :: lda_in, ldb_in, ldc_in
    character :: opa_in, opb_in
    integer :: opa_in_cublas, opb_in_cublas
    if (m < 1 .or. n < 1 .or. k < 1) return
    alpha_in = (1d0, 0d0); beta_in = (0d0, 0d0)
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
    istat = cublas_init()
    opa_in_cublas = get_m_op_cublas(opa_in)
    opb_in_cublas = get_m_op_cublas(opb_in)
    istat = cublaszgemm3m(cublas_handle, opa_in_cublas, opb_in_cublas,  m, n, k, &
                        & alpha_in, a, lda_in , b, ldb_in, beta_in, c, ldc_in)
  end function zmm_d
  integer function zmm_batch_d(a, b, c, m, n, k, nbatch, opa, opb, alpha, beta, lda, ldb, ldc, samea, sameb, comm) result(istat)
    complex(8), device :: a(*), b(*), c(*)
    integer, intent(in) :: m, n, k, nbatch
    character, intent(in), optional :: opa, opb
    complex(8), intent(in), optional :: alpha, beta
    integer, optional :: lda, ldb, ldc
    logical, optional :: samea, sameb
    integer, optional :: comm
    integer(8) :: stridea, strideb, stridec
    complex(8) :: alpha_in, beta_in
    integer :: lda_in, ldb_in, ldc_in
    character :: opa_in, opb_in
    integer :: opa_in_cublas, opb_in_cublas
    if (nbatch < 1) return
    if (m < 1 .or. n < 1 .or. k < 1) return
    alpha_in = (1d0, 0d0); beta_in = (0d0, 0d0)
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
    istat = cublas_init()
    opa_in_cublas = get_m_op_cublas(opa_in)
    opb_in_cublas = get_m_op_cublas(opb_in)
    istat = cublaszgemmstridedbatched(cublas_handle, opa_in_cublas, opb_in_cublas,  m, n, k,  &
               &  alpha_in, a, lda_in, stridea, b, ldb_in, strideb, beta_in, c, ldc_in, stridec, nbatch)
  end function zmm_batch_d
#endif

#ifdef __GPU
  integer function cublas_init() result(istat)
    istat = 0
    if(.not.set_cublas_handle) then 
      istat = cublascreate(cublas_handle)
      set_cublas_handle = .true.
    endif
  end function cublas_init
  integer function get_m_op_cublas(m_op_blas) result(m_op_cublas)
    character, intent(in) :: m_op_blas
    select case (m_op_blas)
      case(m_op_c) ; m_op_cublas = cublas_op_c
      case(m_op_t) ; m_op_cublas = cublas_op_t
      case default ; m_op_cublas = cublas_op_n
    end select
  end  function get_m_op_cublas
#endif
  subroutine int_split(ndata, nsplit, irank, iini, iend, n)
    integer, intent(in) :: ndata, nsplit, irank
    integer, intent(out) :: iini, iend, n
    n = (ndata + irank)/nsplit
    iini = (ndata/nsplit)*irank + max(irank + mod(ndata, nsplit) - nsplit, 0) + 1  
    iend = iini + n - 1
  end subroutine
end module m_blas
