!> Backend "realsgemm" of m_linalg_policy: complex single-precision C = alpha op(A) B + beta C as ONE real SGEMM
!> on the same memory.  B (k x n complex, opB = N) is read as a real (2k x n) matrix and C (m x n) is written as a
!> real (2m x n) one; only op(A) is copied, into A' (2k x 2m real): with a = alpha op(A)(j,l),
!>   column 2j-1 of A' holds (Re a, -Im a) at rows (2l-1, 2l),  column 2j holds (Im a, Re a).
!> Then A'^T B is C with Re and Im interleaved; the same 8mnk flops in FP32.  cuBLAS cgemm runs at about 31 TFLOPS
!> on RTX 5090 against 50-54 for this route (m=k=1053, n >= 1000; TOOLS/ozbench, 2026-09-27).
!> A non-real beta is applied to C first (one pass over C).  opB /= N returns -1 and the caller uses cuBLAS.
!> key >= 0: the caller promises that the same key means the same A, op and alpha until realsgemm_reset; A' is then
!> kept in a pool (ECALJ_LA_CACHE_GB, default 4, and at most 1/4 of the free device memory) and reused.
!> Without a key A' is built for every call (4km floats written); below n = 256 that costs more than the faster GEMM
!> saves (break-even n ~ 160 on RTX 5090), so such calls also return -1 and go to cuBLAS.
module m_la_realsgemm
#ifdef __GPU
  use cudafor
  use cublas_v2
  implicit none
  private
  public :: realsgemm_c, realsgemm_reset
  real(4), device, allocatable, save :: ap(:)        ! A' of a call without key
  real(4), device, allocatable, save :: pool(:)      ! A' of keyed calls, back to back
  integer(8), save :: pool_cap = -1, pool_used = 0
  integer, parameter :: maxslot = 4096
  integer, parameter :: nminkey = 256                ! smallest n worth building A' for without a key
  integer, save :: nslot = 0
  integer, save :: skey(maxslot), sm(maxslot), sk(maxslot)
  character, save :: sop(maxslot)
  complex(4), save :: salpha(maxslot)
  integer(8), save :: soff(maxslot)
contains
  integer function realsgemm_c(handle, opa, opb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, ctype, key) result(istat)
    type(cublasHandle) :: handle
    character, intent(in) :: opa, opb
    integer, intent(in) :: m, n, k, lda, ldb, ldc, ctype, key
    complex(4), intent(in) :: alpha, beta
    complex(4), device :: a(*), b(*), c(*)
    real(4) :: one_r, beta_r
    integer(8) :: off
    integer :: is
    istat = -1
    if (opb /= 'N' .and. opb /= 'n') return
    if (4_8*k*m >= huge(1)) return
    if (key < 0 .and. n < nminkey) return             ! building A' for one small product does not pay
    if (aimag(beta) /= 0.0) then                    ! complex beta: C := beta C first, then add with beta 1
      call scale_c(c, m, n, ldc, beta)
      beta_r = 1.0
    else
      beta_r = real(beta)
    endif
    one_r = 1.0
    is = 0
    if (key >= 0) is = findslot(key, opa, m, k, alpha)
    if (is > 0) then                                 ! A' of this key is in the pool
      off = soff(is)
      istat = cublasGemmEx(handle, CUBLAS_OP_T, CUBLAS_OP_N, 2*m, n, 2*k, one_r, pool(off+1:off+4_8*k*m), CUDA_R_32F, 2*k, &
                           b, CUDA_R_32F, 2*ldb, beta_r, c, CUDA_R_32F, 2*ldc, ctype, CUBLAS_GEMM_DEFAULT)
      return
    endif
    if (key >= 0) is = newslot(key, opa, m, k, alpha)
    if (is > 0) then
      off = soff(is)
      call make_ap(pool(off+1:off+4_8*k*m), opa, m, k, alpha, a, lda)
      istat = cublasGemmEx(handle, CUBLAS_OP_T, CUBLAS_OP_N, 2*m, n, 2*k, one_r, pool(off+1:off+4_8*k*m), CUDA_R_32F, 2*k, &
                           b, CUDA_R_32F, 2*ldb, beta_r, c, CUDA_R_32F, 2*ldc, ctype, CUBLAS_GEMM_DEFAULT)
      return
    endif
    if (allocated(ap)) then
      if (size(ap, kind=8) < 4_8*k*m) deallocate(ap)
    endif
    if (.not. allocated(ap)) allocate(ap(4_8*k*m))
    call make_ap(ap, opa, m, k, alpha, a, lda)
    istat = cublasGemmEx(handle, CUBLAS_OP_T, CUBLAS_OP_N, 2*m, n, 2*k, one_r, ap, CUDA_R_32F, 2*k, &
                         b, CUDA_R_32F, 2*ldb, beta_r, c, CUDA_R_32F, 2*ldc, ctype, CUBLAS_GEMM_DEFAULT)
  end function realsgemm_c

  subroutine make_ap(w, opa, m, k, alpha, a, lda)
    !> A' (2k x 2m) from alpha op(A); op(A)(j,l) = A(j,l) (N), A(l,j) (T), conjg(A(l,j)) (C).
    real(4), device :: w(*)
    character, intent(in) :: opa
    integer, intent(in) :: m, k, lda
    complex(4), intent(in) :: alpha
    complex(4), device :: a(*)
    integer :: j, l, iop
    complex(4) :: v
    logical :: unit
    iop = 0                                          ! 0: N, 1: T, 2: C
    if (opa == 'T' .or. opa == 't') iop = 1
    if (opa == 'C' .or. opa == 'c') iop = 2
    unit = alpha == (1.0, 0.0)
    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, m
      do l = 1, k
        if (iop == 0) then
          v = a(j + (l-1)*lda)
        else
          v = a(l + (j-1)*lda)
        endif
        if (iop == 2) v = conjg(v)
        if (.not. unit) v = alpha*v
        w(2*l-1 + (2*j-2)*2*k) =  real(v)
        w(2*l   + (2*j-2)*2*k) = -aimag(v)
        w(2*l-1 + (2*j-1)*2*k) =  aimag(v)
        w(2*l   + (2*j-1)*2*k) =  real(v)
      enddo
    enddo
  end subroutine make_ap

  subroutine scale_c(c, m, n, ldc, beta)
    complex(4), device :: c(*)
    integer, intent(in) :: m, n, ldc
    complex(4), intent(in) :: beta
    integer :: i, j
    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, n
      do i = 1, m
        c(i + (j-1)*ldc) = beta*c(i + (j-1)*ldc)
      enddo
    enddo
  end subroutine scale_c

  integer function findslot(key, opa, m, k, alpha) result(is)
    integer, intent(in) :: key, m, k
    character, intent(in) :: opa
    complex(4), intent(in) :: alpha
    integer :: i
    is = 0
    do i = 1, nslot
      if (skey(i) == key .and. sop(i) == opa .and. sm(i) == m .and. sk(i) == k .and. salpha(i) == alpha) then
        is = i
        return
      endif
    enddo
  end function findslot

  integer function newslot(key, opa, m, k, alpha) result(is)
    !> A new pool slot for this key, or 0 when the pool is full (the call then runs without keeping A').
    integer, intent(in) :: key, m, k
    character, intent(in) :: opa
    complex(4), intent(in) :: alpha
    integer(8) :: need
    is = 0
    if (pool_cap < 0) call init_pool()
    need = 4_8*k*m
    if (nslot >= maxslot .or. pool_used + need > pool_cap) return
    nslot = nslot + 1
    is = nslot
    skey(is) = key; sop(is) = opa; sm(is) = m; sk(is) = k; salpha(is) = alpha
    soff(is) = pool_used
    pool_used = pool_used + need
  end function newslot

  subroutine init_pool()
    character(32) :: cv
    integer :: st, ios, istat
    real(8) :: gb
    integer(8) :: free, total
    gb = 4d0
    call get_environment_variable('ECALJ_LA_CACHE_GB', cv, status=st)
    if (st == 0) then
      read(cv,*,iostat=ios) gb
      if (ios /= 0) gb = 4d0
    endif
    istat = cudaMemGetInfo(free, total)
    pool_cap = int(min(gb*1d9, 0.25d0*real(free,8))/4d0, 8)      ! in reals
    if (pool_cap < 1) then
      pool_cap = 0
      return
    endif
    allocate(pool(pool_cap), stat=istat)
    if (istat /= 0) pool_cap = 0
  end subroutine init_pool

  subroutine realsgemm_reset()
    !> Forget the kept A' (call when the matrices behind the keys change, e.g. at the next q point).
    nslot = 0
    pool_used = 0
  end subroutine realsgemm_reset
#endif
end module m_la_realsgemm
