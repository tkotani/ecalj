module m_blas !wrapper for BLAS and cuBLAS
  !$use omp_lib
  use m_gemmul8, only: gemmul8_init, gemmul8_pays, num_moduli_d, num_moduli_z, num_moduli_c, fastmode_gemmul8
#ifdef __GPU
  use cublas_v2
  use cudafor
#endif
  implicit none
  public :: int_split, la_cache_reset
  public :: m_op_n, m_op_t, m_op_c
  public :: cmm_h, cmm_batch_h, zmm_h, zmm_batch_h, dmm_h, dmv_h, zmv_h, zvv_h
#ifdef __GPU
  public :: cmm_d, cmm_batch_d, zmm_d, zmm_batch_d, dmm_d, dmv_d, zmv_d, zvv_d, cmm_h16_d, sigma_fp16
  public :: cublas_init, cublas_handle, cublas_finalize, cublas_set_stream
  type(cublashandle), target :: cublas_handle
  logical, save :: set_cublas_handle = .false.
  integer(cuda_stream_kind), save :: blas_stream = 0   ! stream of all device products (cublas_set_stream)
#endif
  character, parameter :: m_op_n = 'N', m_op_t = 'T', m_op_c = 'C'
  integer, parameter :: BACKEND_BLAS = 0 !BLAS/cuBLAS
  integer, parameter :: BACKEND_GEMMUL8 = 1
  integer, parameter :: BACKEND_AUTO = 2
  integer, parameter :: BACKEND_BLAS_FP32 = 3 !cuBLAS with FP32 arithmetic also at level tf32 (single precision only)
  integer, parameter :: BACKEND_SIGMA = 4 !a product of Sigma_c: under --sigma_tf32 the rows of level tf32 and TF32
  public :: BACKEND_BLAS_FP32, BACKEND_SIGMA
contains
  integer function cmm_h(a, b, c, m, n, k, opa, opb, alpha, beta, lda, ldb, ldc, policy, key) result(istat)
    complex(4) :: a(*), b(*), c(*)
    integer, intent(in) :: m, n, k
    character, intent(in), optional :: opa, opb
    complex(4), intent(in), optional :: alpha, beta
    integer, intent(in), optional :: lda, ldb, ldc, policy !policy is dummy
    integer, intent(in), optional :: key   ! accepted for the device version's interface; unused on the host
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
    include "mpif.h"   ! NOT `use mpi`: nvfortran cannot resolve the generic
                      ! mpi_bcast with a complex(4) buffer + integer(8) count;
                      ! the legacy include keeps the call at implicit interface.
    complex(4) :: a(*), b(*), c(*)
    integer, intent(in) :: m, n, k, nbatch
    character, intent(in), optional :: opa, opb
    complex(4), intent(in), optional :: alpha, beta
    integer, intent(in), optional :: lda, ldb, ldc
    logical, intent(in), optional :: samea, sameb
    integer, intent(in), optional :: comm
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
    real(8) :: a(*), x(*), y(*)
    integer, intent(in) :: m, n
    character, intent(in), optional :: opa
    real(8), intent(in), optional :: alpha, beta
    integer, intent(in), optional :: lda, incx, incy
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
  integer function zvv_h(x, y, n, res, incx, incy) result(istat)
    implicit none
#ifdef __GPU
    complex(8),external :: zdotc
#endif
    complex(8) :: x(*), y(*)
    complex(8) :: res
    integer, intent(in) :: n
    integer, intent(in), optional :: incx, incy
    integer :: incx_in, incy_in
    if (n < 1) return
    incx_in = 1; incy_in = 1
    if(present(incx)) incx_in = incx
    if(present(incy)) incy_in = incy
#ifdef __GPU
    res = zdotc   (n, x, incx_in, y, incy_in)
#else
    call zdotc(res,n, x, incx_in, y, incy_in) !because of https://gitlab.com/QEF/q-e/-/wikis/Support/zdotc-crash
#endif
    istat = 0
  end function zvv_h
  integer function zmv_h(a, x, y, m, n, opa, alpha, beta, lda, incx, incy) result(istat)
    complex(8) :: a(*), x(*), y(*)
    integer, intent(in) :: m, n
    character, intent(in), optional :: opa
    complex(8), intent(in), optional :: alpha, beta
    integer, intent(in), optional :: lda, incx, incy
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
    integer, intent(in), optional :: lda, ldb, ldc
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
  integer function zmm_h(a, b, c, m, n, k, opa, opb, alpha, beta, lda, ldb, ldc, policy, key) result(istat)
    complex(8) :: a(*), b(*), c(*)
    integer, intent(in) :: m, n, k
    character, intent(in), optional :: opa, opb
    complex(8), intent(in), optional :: alpha, beta
    integer, intent(in), optional :: lda, ldb, ldc, policy !policy is dummy
    integer, intent(in), optional :: key   ! accepted for the device version's interface; unused on the host
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
    include "mpif.h"   ! NOT `use mpi`: see cmm_batch_h above.
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
  integer function cmm_d(a, b, c, m, n, k, opa, opb, alpha, beta, lda, ldb, ldc, policy, key) result(istat)
    !> C = alpha op(A) op(B) + beta C, complex single precision on the device.  The backend (cuBLAS, the real-SGEMM
    !> route or GEMMul8) comes from m_linalg_policy; policy=BACKEND_BLAS / BACKEND_GEMMUL8 forces one.
    !> key >= 0 (optional): same key = same A until la_cache_reset, so backends may keep their form of A.
    use cublas_v2, m_type =>CUDA_C_32F, algo => cublas_gemm_default
    use m_linalg_policy, only: la_backend, la_moduli, la_level, la_sigma_tf32, OP_CGEMM, BK_CUBLAS, BK_REALSGEMM, BK_GEMMUL8, &
                               BK_REALHGEMM
    use m_la_realsgemm, only: realsgemm_c, realhgemm_c
    complex(4), device, target :: a(*), b(*), c(*)
    integer, intent(in) :: m, n, k
    character, intent(in), optional :: opa, opb
    complex(4), intent(in), optional :: alpha, beta
    integer, intent(in), optional :: lda, ldb, ldc, policy, key
    complex(4) :: alpha_in, beta_in
    integer :: lda_in, ldb_in, ldc_in, policy_in, key_in, bk, ctype
    character :: opa_in, opb_in
    integer :: opa_in_cublas, opb_in_cublas
    logical :: issigma
    istat = 0
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
    policy_in = BACKEND_AUTO
    if(present(policy)) policy_in = policy
    key_in = -1
    if(present(key)) key_in = key
    istat = cublas_init()
    opa_in_cublas = get_m_op_cublas(opa_in)
    opb_in_cublas = get_m_op_cublas(opb_in)
    ctype = merge(CUBLAS_COMPUTE_32F_FAST_TF32, CUBLAS_COMPUTE_32F, la_level() == 'tf32')  ! level tf32: TF32, else FP32
    issigma = .false.
    select case (policy_in)
    case (BACKEND_BLAS);      bk = BK_CUBLAS
    case (BACKEND_BLAS_FP32); bk = BK_CUBLAS; ctype = CUBLAS_COMPUTE_32F
    case (BACKEND_SIGMA)
      issigma = la_sigma_tf32()
      bk = la_backend(OP_CGEMM, m, n, k, sigma=issigma)
      if (issigma) ctype = CUBLAS_COMPUTE_32F_FAST_TF32
    case (BACKEND_GEMMUL8);   bk = BK_GEMMUL8
    case default;             bk = la_backend(OP_CGEMM, m, n, k)
    end select
    if (bk == BK_REALSGEMM) then
      istat = realsgemm_c(cublas_handle, opa_in, opb_in, m, n, k, alpha_in, a, lda_in, b, ldb_in, beta_in, c, ldc_in, &
                          ctype, key_in)
      if (istat /= -1) return
      bk = BK_CUBLAS                                 ! the route needs opB = N
    endif
    if (bk == BK_REALHGEMM) then
      istat = realhgemm_c(cublas_handle, opa_in, opb_in, m, n, k, alpha_in, a, lda_in, b, ldb_in, beta_in, c, ldc_in, key_in)
      if (istat /= -1) return
      bk = BK_CUBLAS                                 ! opB /= N or B not contiguous: cuBLAS
    endif
    if (bk == BK_GEMMUL8 .and. .not. gemmul8_pays(m, n, k)) bk = BK_CUBLAS   ! too small for the split
    if (bk == BK_GEMMUL8) then
#ifdef __GEMMUL8
      block
        use m_gemmul8, only: gemmul8_handle, gemmul8_cgemm, gemmul8_set_stream
        use iso_c_binding
        integer :: nm
        istat = gemmul8_init()
        call gemmul8_set_stream(gemmul8_handle, blas_stream)
        nm = la_moduli(OP_CGEMM, m, n, k, sigma=issigma)
        if (nm <= 0) nm = num_moduli_c
        call gemmul8_cgemm(gemmul8_handle, opa_in_cublas, opb_in_cublas, m, n, k, alpha_in, &
                           c_loc(a), lda_in, c_loc(b), ldb_in, beta_in, c_loc(c), ldc_in, nm, fastmode_gemmul8, key_in)
      endblock
      return
#else
      call rx0('Error: the linalg policy chose gemmul8, but the GEMMul8 library is not linked.')
#endif
    endif
    istat = cublasGemmEX(cublas_handle, opa_in_cublas, opb_in_cublas, m, n, k,  &
                         alpha_in, a, m_type, lda_in, b, m_type, ldb_in, beta_in, c, m_type, ldc_in, ctype, algo)
  end function cmm_d
  logical function sigma_fp16(m, n, k)
    !> Does the Sigma_c product of this size run on the FP16 route (realhgemm, rows of level tf32 under --sigma_tf32)?
    !> Then the caller may hand B over in FP16 itself (cmm_h16_d).
    use m_linalg_policy, only: la_backend, la_sigma_tf32, OP_CGEMM, BK_REALHGEMM
    integer, intent(in) :: m, n, k
    sigma_fp16 = .false.
    if (.not. la_sigma_tf32()) return
    sigma_fp16 = la_backend(OP_CGEMM, m, n, k, sigma=.true.) == BK_REALHGEMM
  end function sigma_fp16
  integer function cmm_h16_d(a, bh, fb, c, m, n, k, opa, beta, key) result(istat)
    !> C = op(A) B + beta C with B given in FP16 by the caller: bh = the (2k x n) real form of B (complex B read as real,
    !> column by column) times fb, a power of 2 with |bh| < 2^14.  A is kept under the key in FP16 (m_la_realsgemm).
    !> Only when sigma_fp16(m,n,k); opB = N, ldb = k, ldc = m.
    use m_la_realsgemm, only: realhgemm_bh_c
    complex(4), device :: a(*), c(*)
    real(2), device :: bh(*)
    real(4), intent(in) :: fb
    integer, intent(in) :: m, n, k
    character, intent(in) :: opa
    complex(4), intent(in), optional :: beta
    integer, intent(in), optional :: key
    complex(4) :: beta_in
    integer :: lda_in, key_in
    beta_in = (0.0, 0.0)
    if (present(beta)) beta_in = beta
    key_in = -1
    if (present(key)) key_in = key
    lda_in = m
    if (opa == m_op_t .or. opa == m_op_c) lda_in = k
    istat = cublas_init()
    istat = realhgemm_bh_c(cublas_handle, opa, m, n, k, (1.0, 0.0), a, lda_in, bh, fb, beta_in, c, m, key_in)
  end function cmm_h16_d
  integer function cmm_batch_d(a, b, c, m, n, k, nbatch, opa, opb, alpha, beta, lda, ldb, ldc, samea, sameb, comm) result(istat)
    use cublas_v2, m_type =>CUDA_C_32F, algo => cublas_gemm_default
    use m_linalg_policy, only: la_level
    complex(4), device :: a(*), b(*), c(*)
    integer, intent(in) :: m, n, k, nbatch
    character, intent(in), optional :: opa, opb
    complex(4), intent(in), optional :: alpha, beta
    integer, intent(in), optional :: lda, ldb, ldc
    logical, intent(in), optional :: samea, sameb
    integer, intent(in), optional :: comm
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
                                       c, m_type, ldc_in, stridec, nbatch, &
                                       merge(CUBLAS_COMPUTE_32F_FAST_TF32, CUBLAS_COMPUTE_32F, la_level() == 'tf32'), algo)
  end function cmm_batch_d
  integer function zvv_d(x, y, n, res, incx, incy) result(istat)
    implicit none
    complex(8), device :: x(*), y(*)
    complex(8) :: res
    integer, intent(in) :: n
    integer, optional :: incx, incy
    integer :: incx_in, incy_in
    if (n < 1) return
    incx_in = 1; incy_in = 1
    if(present(incx)) incx_in = incx
    if(present(incy)) incy_in = incy
    istat = cublas_init()
    istat = cublaszdotc(cublas_handle, n, x, incx_in, y, incy_in, res)
  end function zvv_d
  integer function dmv_d(a, x, y, m, n, opa, alpha, beta, lda, incx, incy) result(istat)
    real(8), device :: a(*), x(*), y(*)
    integer, intent(in) :: m, n
    !caution: size of matrix a is m x n (not size of op(A))
    character, intent(in), optional :: opa
    real(8), intent(in), optional :: alpha, beta
    integer, intent(in), optional :: lda, incx, incy
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
    complex(8), device :: a(*), x(*), y(*)
    integer, intent(in) :: m, n
    character, intent(in), optional :: opa
    complex(8), intent(in), optional :: alpha, beta
    integer, intent(in), optional :: lda, incx, incy
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
  integer function dmm_d(a, b, c, m, n, k, opa, opb, alpha, beta, lda, ldb, ldc, policy, key) result(istat)
    !> Real double precision on the device; backend cuBLAS or GEMMul8 from m_linalg_policy (see cmm_d).
    use m_linalg_policy, only: la_backend, la_moduli, OP_DGEMM, BK_CUBLAS, BK_GEMMUL8
    real(8), device, target :: a(*), b(*), c(*)
    integer, intent(in) :: m, n, k
    character, intent(in), optional :: opa, opb
    real(8), intent(in), optional :: alpha, beta
    integer, intent(in), optional :: lda, ldb, ldc, policy, key
    real(8) :: alpha_in, beta_in
    integer :: lda_in, ldb_in, ldc_in, policy_in, key_in, bk
    character :: opa_in, opb_in
    integer :: opa_in_cublas, opb_in_cublas
    istat = 0
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
    policy_in = BACKEND_AUTO
    if(present(policy)) policy_in = policy
    key_in = -1
    if(present(key)) key_in = key
    istat = cublas_init()
    opa_in_cublas = get_m_op_cublas(opa_in)
    opb_in_cublas = get_m_op_cublas(opb_in)
    select case (policy_in)
    case (BACKEND_BLAS);    bk = BK_CUBLAS
    case (BACKEND_GEMMUL8); bk = BK_GEMMUL8
    case default;           bk = la_backend(OP_DGEMM, m, n, k)
    end select
    if (bk == BK_GEMMUL8 .and. .not. gemmul8_pays(m, n, k)) bk = BK_CUBLAS   ! too small for the split
    if (bk == BK_GEMMUL8) then
#ifdef __GEMMUL8
      block
        use m_gemmul8, only: gemmul8_handle, gemmul8_dgemm, gemmul8_set_stream
        use iso_c_binding
        integer :: nm
        istat = gemmul8_init()
        call gemmul8_set_stream(gemmul8_handle, blas_stream)
        nm = la_moduli(OP_DGEMM, m, n, k)
        if (nm <= 0) nm = num_moduli_d
        call gemmul8_dgemm(gemmul8_handle, opa_in_cublas, opb_in_cublas, m, n, k, alpha_in, &
                           c_loc(a), lda_in, c_loc(b), ldb_in, beta_in, c_loc(c), ldc_in, nm, fastmode_gemmul8, key_in)
      endblock
      return
#else
      call rx0('Error: the linalg policy chose gemmul8, but the GEMMul8 library is not linked.')
#endif
    endif
    istat = cublasdgemm(cublas_handle, opa_in_cublas, opb_in_cublas,  m, n, k, &
                      & alpha_in, a, lda_in , b, ldb_in, beta_in, c, ldc_in)
  end function dmm_d
  integer function zmm_d(a, b, c, m, n, k, opa, opb, alpha, beta, lda, ldb, ldc, policy, key) result(istat)
    !> Complex double precision on the device; backend cuBLAS (zgemm3m) or GEMMul8 from m_linalg_policy (see cmm_d).
    use m_linalg_policy, only: la_backend, la_moduli, OP_ZGEMM, BK_CUBLAS, BK_GEMMUL8
    complex(8), device, target :: a(*), b(*), c(*)
    integer, intent(in) :: m, n, k
    character, intent(in), optional :: opa, opb
    complex(8), intent(in), optional :: alpha, beta
    integer, intent(in), optional :: lda, ldb, ldc, policy, key
    complex(8) :: alpha_in, beta_in
    integer :: lda_in, ldb_in, ldc_in, policy_in, key_in, bk
    character :: opa_in, opb_in
    integer :: opa_in_cublas, opb_in_cublas
    istat = 0
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
    policy_in = BACKEND_AUTO
    if(present(policy)) policy_in = policy
    key_in = -1
    if(present(key)) key_in = key
    istat = cublas_init()
    opa_in_cublas = get_m_op_cublas(opa_in)
    opb_in_cublas = get_m_op_cublas(opb_in)
    select case (policy_in)
    case (BACKEND_BLAS);    bk = BK_CUBLAS
    case (BACKEND_GEMMUL8); bk = BK_GEMMUL8
    case default;           bk = la_backend(OP_ZGEMM, m, n, k)
    end select
    if (bk == BK_GEMMUL8 .and. .not. gemmul8_pays(m, n, k)) bk = BK_CUBLAS   ! too small for the split
    if (bk == BK_GEMMUL8) then
#ifdef __GEMMUL8
      block
        use m_gemmul8, only: gemmul8_handle, gemmul8_zgemm, gemmul8_set_stream
        use iso_c_binding
        integer :: nm
        istat = gemmul8_init()
        call gemmul8_set_stream(gemmul8_handle, blas_stream)
        nm = la_moduli(OP_ZGEMM, m, n, k)
        if (nm <= 0) nm = num_moduli_z
        call gemmul8_zgemm(gemmul8_handle, opa_in_cublas, opb_in_cublas, m, n, k, alpha_in, &
                           c_loc(a), lda_in, c_loc(b), ldb_in, beta_in, c_loc(c), ldc_in, nm, fastmode_gemmul8, key_in)
      endblock
      return
#else
      call rx0('Error: the linalg policy chose gemmul8, but the GEMMul8 library is not linked.')
#endif
    endif
    istat = cublaszgemm3m(cublas_handle, opa_in_cublas, opb_in_cublas,  m, n, k, &
                        & alpha_in, a, lda_in , b, ldb_in, beta_in, c, ldc_in)
  end function zmm_d
  integer function zmm_batch_d(a, b, c, m, n, k, nbatch, opa, opb, alpha, beta, lda, ldb, ldc, samea, sameb, comm) result(istat)
    complex(8), device :: a(*), b(*), c(*)
    integer, intent(in) :: m, n, k, nbatch
    character, intent(in), optional :: opa, opb
    complex(8), intent(in), optional :: alpha, beta
    integer, intent(in), optional :: lda, ldb, ldc
    logical, intent(in), optional :: samea, sameb
    integer, intent(in), optional :: comm
    integer(8) :: stridea, strideb, stridec
    complex(8) :: alpha_in, beta_in
    integer :: lda_in, ldb_in, ldc_in
    character :: opa_in, opb_in
    integer :: opa_in_cublas, opb_in_cublas
    integer :: i
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
  subroutine cublas_set_stream(stream)
    !> From now on every device product of this module runs on this CUDA stream: cuBLAS, the realsgemm route
    !> (it takes the stream of the handle) and GEMMul8 (set before each call).  0 = the default stream.
    !> Sigma_c (m_sxcf_sc) puts them on the stream of OpenACC queue 1 with its async(1) kernels.
    use cudafor
    implicit none
    integer(cuda_stream_kind), intent(in) :: stream
    integer :: istat
    istat = cublas_init()
    istat = cublasSetStream(cublas_handle, stream)
    blas_stream = stream
  end subroutine
  integer function cublas_finalize() result(istat)
    istat = 0
    if(set_cublas_handle) then
        istat = cublasdestroy(cublas_handle)
        set_cublas_handle = .false.
    endif
  end function cublas_finalize
  integer function get_m_op_cublas(m_op_blas) result(m_op_cublas)
    character, intent(in) :: m_op_blas
    select case (m_op_blas)
      case(m_op_c) ; m_op_cublas = cublas_op_c
      case(m_op_t) ; m_op_cublas = cublas_op_t
      case default ; m_op_cublas = cublas_op_n
    end select
  end  function get_m_op_cublas
#endif
  subroutine la_cache_reset()
    !> Forget the matrices kept under keys (A' of realsgemm, the split A of GEMMul8).  Call when the matrices
    !> behind the keys change, e.g. at the start of each q point in hgw.  No-op without GPU.
#ifdef __GPU
    use m_la_realsgemm, only: realsgemm_reset
#ifdef __GEMMUL8
    use m_gemmul8, only: gemmul8_cache_reset
#endif
    call realsgemm_reset()
#ifdef __GEMMUL8
    call gemmul8_cache_reset()
#endif
#endif
  end subroutine la_cache_reset
  subroutine int_split(ndata, nsplit, irank, iini, iend, n, start_index)
    integer, intent(in) :: ndata, nsplit, irank
    integer, intent(out) :: iini, iend, n
    integer, optional, intent(in) :: start_index
    n = (ndata + irank)/nsplit
    iini = (ndata/nsplit)*irank + max(irank + mod(ndata, nsplit) - nsplit, 0) + 1  
    iend = iini + n - 1
    if(present(start_index)) then
      iini = iini + start_index - 1
      iend = iend + start_index - 1
    endif
  end subroutine
end module m_blas
