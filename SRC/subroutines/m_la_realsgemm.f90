!> Backend "realsgemm" of m_linalg_policy: complex single-precision C = alpha op(A) B + beta C as ONE real SGEMM
!> on the same memory.  B (k x n complex, opB = N) is read as a real (2k x n) matrix and C (m x n) is written as a
!> real (2m x n) one; only op(A) is copied, with a = alpha op(A)(j,l):
!>   opA = C, T: A' (2k x 2m): column 2j-1 holds (Re a, -Im a) at rows (2l-1, 2l), column 2j holds (Im a, Re a);
!>               C = A'^T B.
!>   opA = N:    A'' (2m x 2k): column 2l-1 holds (Re a, Im a) at rows (2j-1, 2j) (column l of A read as real),
!>               column 2l holds (-Im a, Re a) (i times it); C = A'' B.
!> Both copies read and write A along its columns (coalesced).  The same 8mnk flops in FP32; on GPUs where cuBLAS
!> cgemm runs well below the SGEMM rate this route is faster (2026-09-27, RTX 5090, m=k=1053: cgemm about 31 TFLOPS,
!> this route 50-54).  linalgtune measures both (time and error) per shape.
!> A non-real beta is applied to C first (one pass over C).  opB /= N returns -1 and the caller uses cuBLAS.
!> The copies run on the stream of the cuBLAS handle, as the SGEMM does (m_blas cublas_set_stream).  Scratch arrays
!> and pool slots are reused without a device sync: calls must stay on one stream, with a device sync at each stream
!> switch (m_sxcf_sc sigma_stream_begin/end).
!> key >= 0: the caller promises that the same key means the same matrix A until realsgemm_reset; A' is kept per
!> (key, op, m, k, alpha) in a pool (ECALJ_LA_CACHE_GB, default 4, and at most 1/4 of the free device memory; GEMMul8
!> has its own pool of the same budget) and reused.
!> Without a key A' is built for every call (4km floats written, coalesced); that pays back from n of a few tens,
!> so products without a key and n < nminkey return -1 and go to cuBLAS.
!>
!> Backend "realhgemm": the same route with A' and B in FP16, FP32 accumulation and C (tensor cores).  FP16 has the
!> 10-bit mantissa of TF32, so the input rounding is that of TF32 (BF16 would give 8 times the error; 2026-09-27,
!> RTX 5090: 1.75-1.9 times as fast as TF32).  FP16
!> overflows at 65504, so A' and B are scaled by powers of 2 to a largest element in [2^13, 2^14): elements down to
!> 2^-27 of the largest keep the full relative precision, smaller ones lose it (their absolute error stays below
!> 2^-38 of the largest, far under the rounding of the large ones).  A': the scale is found on the
!> host when A' is built (once per key).  B: its largest |element| by cublasIsamax, then the scale and the alpha of
!> the GEMM on the device (cuBLAS in device pointer mode), so the host never waits (Sigma_c runs asynchronously).
!> Needs a contiguous B (ldb = k) besides opB = N.  The FP16 pool of kept A' has half the budget of the FP32 one.
module m_la_realsgemm
#ifdef __GPU
  use cudafor
  use cublas_v2
  use iso_c_binding
  implicit none
  private
  public :: realsgemm_c, realhgemm_c, realhgemm_bh_c, realsgemm_reset
  real(4), device, allocatable, save :: ap(:)        ! A' of a call without key
  real(4), device, allocatable, save :: pool(:)      ! A' of keyed calls, back to back
  real(2), device, allocatable, save :: aph(:), poolh(:), bh(:)   ! realhgemm: the same in FP16, and B in FP16
  real(4), device, allocatable, save :: abd(:)       ! realhgemm: alpha and beta of the GEMM, set on the device
  integer, device, allocatable, save :: imaxd(:)     ! realhgemm: position of the largest |B| (cublasIsamax)
  integer, parameter :: ip32 = 1, ip16 = 2          ! the two pools of kept A'
  integer(8), save :: pool_cap(2) = -1, pool_used(2) = 0   ! in elements
  integer, parameter :: maxslot = 4096
  integer, parameter :: nminkey = 16                 ! smallest n worth building A' for without a key
  integer, parameter :: ehalf = 14                   ! FP16 data scaled to a largest element in [2^13, 2^14)
  integer, save :: nslot(2) = 0
  integer, save :: skey(maxslot,2), sm(maxslot,2), sk(maxslot,2)
  character, save :: sop(maxslot,2)
  complex(4), save :: salpha(maxslot,2)
  integer(8), save :: soff(maxslot,2)
  real(4), save :: sfa(maxslot)                      ! the scale of the FP16 A' of a slot
contains
  integer function realsgemm_c(handle, opa, opb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, ctype, key) result(istat)
    type(cublasHandle) :: handle
    character, intent(in) :: opa, opb
    integer, intent(in) :: m, n, k, lda, ldb, ldc, ctype, key
    complex(4), intent(in) :: alpha, beta
    complex(4), device :: a(*), b(*), c(*)
    real(4) :: one_r, beta_r
    integer(8) :: off
    integer :: is, opap, ldap, ist
    integer(cuda_stream_kind) :: st
    istat = -1
    if (opb /= 'N' .and. opb /= 'n') return
    if (4_8*k*m >= huge(1)) return                    ! make_ap indexes A' with default integers
    if (key < 0 .and. n < nminkey) return             ! building A' for one small product does not pay
    ist = cublasGetStream(handle, st)
    if (aimag(beta) /= 0.0) then                    ! complex beta: C := beta C first, then add with beta 1
      call scale_c(c, m, n, ldc, beta, st)
      beta_r = 1.0
    else
      beta_r = real(beta)
    endif
    one_r = 1.0
    if (opa == 'N' .or. opa == 'n') then               ! A'' (2m x 2k), C = A'' B
      opap = CUBLAS_OP_N; ldap = 2*m
    else                                               ! A' (2k x 2m), C = A'^T B
      opap = CUBLAS_OP_T; ldap = 2*k
    endif
    is = 0
    if (key >= 0) is = findslot(key, opa, m, k, alpha, ip32)
    if (is > 0) then                                 ! A' of this key is in the pool
      off = soff(is,ip32)
      istat = cublasGemmEx(handle, opap, CUBLAS_OP_N, 2*m, n, 2*k, one_r, pool(off+1:off+4_8*k*m), CUDA_R_32F, ldap, &
                           b, CUDA_R_32F, 2*ldb, beta_r, c, CUDA_R_32F, 2*ldc, ctype, CUBLAS_GEMM_DEFAULT)
      return
    endif
    if (key >= 0) is = newslot(key, opa, m, k, alpha, ip32)
    if (is > 0) then
      off = soff(is,ip32)
      call make_ap(pool(off+1:off+4_8*k*m), opa, m, k, alpha, a, lda, st)
      istat = cublasGemmEx(handle, opap, CUBLAS_OP_N, 2*m, n, 2*k, one_r, pool(off+1:off+4_8*k*m), CUDA_R_32F, ldap, &
                           b, CUDA_R_32F, 2*ldb, beta_r, c, CUDA_R_32F, 2*ldc, ctype, CUBLAS_GEMM_DEFAULT)
      return
    endif
    if (allocated(ap)) then
      if (size(ap, kind=8) < 4_8*k*m) deallocate(ap)
    endif
    if (.not. allocated(ap)) allocate(ap(4_8*k*m))
    call make_ap(ap, opa, m, k, alpha, a, lda, st)
    istat = cublasGemmEx(handle, opap, CUBLAS_OP_N, 2*m, n, 2*k, one_r, ap, CUDA_R_32F, ldap, &
                         b, CUDA_R_32F, 2*ldb, beta_r, c, CUDA_R_32F, 2*ldc, ctype, CUBLAS_GEMM_DEFAULT)
  end function realsgemm_c

  integer function realhgemm_c(handle, opa, opb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, key) result(istat)
    !> realsgemm in FP16 (see the head of this module).  -1: not done here, the caller uses cuBLAS.
    type(cublasHandle) :: handle
    character, intent(in) :: opa, opb
    integer, intent(in) :: m, n, k, lda, ldb, ldc, key
    complex(4), intent(in) :: alpha, beta
    complex(4), device, target :: a(*), b(*), c(*)
    real(4) :: beta_r, fa
    real(4), device, pointer :: br(:)
    integer(8) :: off
    integer :: is, opap, ldap, ist, nb
    integer(cuda_stream_kind) :: st
    istat = -1
    if (opb /= 'N' .and. opb /= 'n') return
    if (ldb /= k) return                               ! B is scanned and converted as one vector
    if (4_8*k*m >= huge(1) .or. 2_8*k*n >= huge(1)) return
    if (key < 0 .and. n < nminkey) return
    ist = cublasGetStream(handle, st)
    if (aimag(beta) /= 0.0) then
      call scale_c(c, m, n, ldc, beta, st)
      beta_r = 1.0
    else
      beta_r = real(beta)
    endif
    if (opa == 'N' .or. opa == 'n') then
      opap = CUBLAS_OP_N; ldap = 2*m
    else
      opap = CUBLAS_OP_T; ldap = 2*k
    endif
    is = 0
    if (key >= 0) then
      is = findslot(key, opa, m, k, alpha, ip16)
      if (is == 0) then
        is = newslot(key, opa, m, k, alpha, ip16)
        if (is > 0) call make_aph(poolh(soff(is,ip16)+1:soff(is,ip16)+4_8*k*m), opa, m, k, alpha, a, lda, st, sfa(is))
      endif
    endif
    if (is == 0) then                                  ! no key, or the pool is full
      if (allocated(aph)) then
        if (size(aph, kind=8) < 4_8*k*m) deallocate(aph)
      endif
      if (.not. allocated(aph)) allocate(aph(4_8*k*m))
      call make_aph(aph, opa, m, k, alpha, a, lda, st, fa)
    else
      fa = sfa(is)
    endif
    nb = 2*k*n
    if (allocated(bh)) then
      if (size(bh) < nb) deallocate(bh)
    endif
    if (.not. allocated(bh)) allocate(bh(nb))
    if (.not. allocated(abd)) allocate(abd(2), imaxd(1))
    call c_f_pointer(c_devloc(b), br, [nb])
    ist = cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_DEVICE)
    ist = cublasIsamax(handle, nb, br, 1, imaxd(1))
    call make_bh(br, nb, bh, imaxd, abd, fa, beta_r, st)
    if (is > 0) then
      off = soff(is,ip16)
      istat = cublasGemmEx(handle, opap, CUBLAS_OP_N, 2*m, n, 2*k, abd(1), poolh(off+1:off+4_8*k*m), CUDA_R_16F, ldap, &
                           bh, CUDA_R_16F, 2*k, abd(2), c, CUDA_R_32F, 2*ldc, CUBLAS_COMPUTE_32F, CUBLAS_GEMM_DEFAULT)
    else
      istat = cublasGemmEx(handle, opap, CUBLAS_OP_N, 2*m, n, 2*k, abd(1), aph, CUDA_R_16F, ldap, &
                           bh, CUDA_R_16F, 2*k, abd(2), c, CUDA_R_32F, 2*ldc, CUBLAS_COMPUTE_32F, CUBLAS_GEMM_DEFAULT)
    endif
    ist = cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST)
  end function realhgemm_c

  integer function realhgemm_bh_c(handle, opa, m, n, k, alpha, a, lda, bh, fb, beta, c, ldc, key) result(istat)
    !> realhgemm with B given by the caller already in FP16: bh = the (2k x n) real form of B times fb (a power of 2 the
    !> caller chose from a bound of |B|, so that |bh| < 2^14).  No scan of B and no conversion here; alpha of the GEMM
    !> 1/(fa fb) is known on the host.  Sigma_c writes bh in its weighting kernel (m_sxcf_sc).
    type(cublasHandle) :: handle
    character, intent(in) :: opa
    integer, intent(in) :: m, n, k, lda, ldc, key
    complex(4), intent(in) :: alpha, beta
    complex(4), device :: a(*), c(*)
    real(2), device :: bh(*)
    real(4), intent(in) :: fb
    real(4) :: beta_r, fa, alpha_r
    integer(8) :: off
    integer :: is, opap, ldap, ist
    integer(cuda_stream_kind) :: st
    istat = -1
    if (4_8*k*m >= huge(1)) return
    ist = cublasGetStream(handle, st)
    if (aimag(beta) /= 0.0) then
      call scale_c(c, m, n, ldc, beta, st)
      beta_r = 1.0
    else
      beta_r = real(beta)
    endif
    if (opa == 'N' .or. opa == 'n') then
      opap = CUBLAS_OP_N; ldap = 2*m
    else
      opap = CUBLAS_OP_T; ldap = 2*k
    endif
    is = 0
    if (key >= 0) then
      is = findslot(key, opa, m, k, alpha, ip16)
      if (is == 0) then
        is = newslot(key, opa, m, k, alpha, ip16)
        if (is > 0) call make_aph(poolh(soff(is,ip16)+1:soff(is,ip16)+4_8*k*m), opa, m, k, alpha, a, lda, st, sfa(is))
      endif
    endif
    if (is == 0) then
      if (allocated(aph)) then
        if (size(aph, kind=8) < 4_8*k*m) deallocate(aph)
      endif
      if (.not. allocated(aph)) allocate(aph(4_8*k*m))
      call make_aph(aph, opa, m, k, alpha, a, lda, st, fa)
    else
      fa = sfa(is)
    endif
    alpha_r = 1.0/(fa*fb)
    if (is > 0) then
      off = soff(is,ip16)
      istat = cublasGemmEx(handle, opap, CUBLAS_OP_N, 2*m, n, 2*k, alpha_r, poolh(off+1:off+4_8*k*m), CUDA_R_16F, ldap, &
                           bh, CUDA_R_16F, 2*k, beta_r, c, CUDA_R_32F, 2*ldc, CUBLAS_COMPUTE_32F, CUBLAS_GEMM_DEFAULT)
    else
      istat = cublasGemmEx(handle, opap, CUBLAS_OP_N, 2*m, n, 2*k, alpha_r, aph, CUDA_R_16F, ldap, &
                           bh, CUDA_R_16F, 2*k, beta_r, c, CUDA_R_32F, 2*ldc, CUBLAS_COMPUTE_32F, CUBLAS_GEMM_DEFAULT)
    endif
  end function realhgemm_bh_c

  subroutine make_ap(w, opa, m, k, alpha, a, lda, st)
    !> From alpha op(A): A' (2k x 2m) for opA = T, C (op(A)(j,l) = A(l,j), conjg(A(l,j))), A'' (2m x 2k) for opA = N.
    !> The inner loop runs down a column of A in both cases (coalesced reads and writes).
    real(4), device :: w(*)
    character, intent(in) :: opa
    integer, intent(in) :: m, k, lda
    complex(4), intent(in) :: alpha
    complex(4), device :: a(*)
    integer(cuda_stream_kind), intent(in) :: st
    integer :: j, l, iop
    complex(4) :: v
    logical :: unit
    iop = 0                                          ! 0: N, 1: T, 2: C
    if (opa == 'T' .or. opa == 't') iop = 1
    if (opa == 'C' .or. opa == 'c') iop = 2
    unit = alpha == (1.0, 0.0)
    if (iop == 0) then
      !$cuf kernel do(2) <<<*,*,stream=st>>>
      do l = 1, k
        do j = 1, m
          v = a(j + (l-1)*lda)
          if (.not. unit) v = alpha*v
          w(2*j-1 + (2*l-2)*2*m) =  real(v)
          w(2*j   + (2*l-2)*2*m) =  aimag(v)
          w(2*j-1 + (2*l-1)*2*m) = -aimag(v)
          w(2*j   + (2*l-1)*2*m) =  real(v)
        enddo
      enddo
      return
    endif
    !$cuf kernel do(2) <<<*,*,stream=st>>>
    do j = 1, m
      do l = 1, k
        v = a(l + (j-1)*lda)
        if (iop == 2) v = conjg(v)
        if (.not. unit) v = alpha*v
        w(2*l-1 + (2*j-2)*2*k) =  real(v)
        w(2*l   + (2*j-2)*2*k) = -aimag(v)
        w(2*l-1 + (2*j-1)*2*k) =  aimag(v)
        w(2*l   + (2*j-1)*2*k) =  real(v)
      enddo
    enddo
  end subroutine make_ap

  subroutine make_aph(w, opa, m, k, alpha, a, lda, st, fa)
    !> make_ap in FP16: the same layout, times fa = 2^(ehalf - exponent(largest element)).  The largest element is a
    !> reduction returned to the host (one host wait per A' built: once per key, at every call without a key).
    real(2), device :: w(*)
    character, intent(in) :: opa
    integer, intent(in) :: m, k, lda
    complex(4), intent(in) :: alpha
    complex(4), device :: a(*)
    integer(cuda_stream_kind), intent(in) :: st
    real(4), intent(out) :: fa
    integer :: j, l, iop
    complex(4) :: v
    real(4) :: amax
    logical :: unit
    iop = 0
    if (opa == 'T' .or. opa == 't') iop = 1
    if (opa == 'C' .or. opa == 'c') iop = 2
    unit = alpha == (1.0, 0.0)
    amax = 0.0
    if (iop == 0) then
      !$cuf kernel do(2) <<<*,*,stream=st>>>
      do l = 1, k
        do j = 1, m
          v = a(j + (l-1)*lda)
          if (.not. unit) v = alpha*v
          amax = max(amax, abs(real(v)), abs(aimag(v)))
        enddo
      enddo
    else
      !$cuf kernel do(2) <<<*,*,stream=st>>>
      do j = 1, m
        do l = 1, k
          v = a(l + (j-1)*lda)
          if (.not. unit) v = alpha*v
          amax = max(amax, abs(real(v)), abs(aimag(v)))
        enddo
      enddo
    endif
    fa = 1.0
    if (amax > 0.0) fa = scale(1.0, ehalf - exponent(amax))
    if (iop == 0) then
      !$cuf kernel do(2) <<<*,*,stream=st>>>
      do l = 1, k
        do j = 1, m
          v = a(j + (l-1)*lda)
          if (.not. unit) v = alpha*v
          v = fa*v
          w(2*j-1 + (2*l-2)*2*m) = real( real(v), kind=2)
          w(2*j   + (2*l-2)*2*m) = real( aimag(v), kind=2)
          w(2*j-1 + (2*l-1)*2*m) = real(-aimag(v), kind=2)
          w(2*j   + (2*l-1)*2*m) = real( real(v), kind=2)
        enddo
      enddo
      return
    endif
    !$cuf kernel do(2) <<<*,*,stream=st>>>
    do j = 1, m
      do l = 1, k
        v = a(l + (j-1)*lda)
        if (iop == 2) v = conjg(v)
        if (.not. unit) v = alpha*v
        v = fa*v
        w(2*l-1 + (2*j-2)*2*k) = real( real(v), kind=2)
        w(2*l   + (2*j-2)*2*k) = real(-aimag(v), kind=2)
        w(2*l-1 + (2*j-1)*2*k) = real( aimag(v), kind=2)
        w(2*l   + (2*j-1)*2*k) = real( real(v), kind=2)
      enddo
    enddo
  end subroutine make_aph

  subroutine make_bh(br, nb, bh, imaxd, abd, fa, beta_r, st)
    !> B (read as nb reals) to FP16 times fb = 2^(ehalf - exponent(|B(imaxd)|)), and the alpha 1/(fa fb) and beta of
    !> the GEMM, all on the device: nothing comes back to the host.
    integer, intent(in) :: nb
    real(4), device :: br(nb)
    real(2), device :: bh(nb)
    integer, device :: imaxd(1)
    real(4), device :: abd(2)
    real(4), intent(in) :: fa, beta_r
    integer(cuda_stream_kind), intent(in) :: st
    integer :: i
    real(4) :: bmax, fb
    !$cuf kernel do <<<*,*,stream=st>>>
    do i = 1, nb
      bmax = abs(br(imaxd(1)))
      fb = 1.0
      if (bmax > 0.0) fb = scale(1.0, ehalf - exponent(bmax))
      bh(i) = real(fb*br(i), kind=2)
      if (i == 1) then
        abd(1) = 1.0/(fa*fb)
        abd(2) = beta_r
      endif
    enddo
  end subroutine make_bh

  subroutine scale_c(c, m, n, ldc, beta, st)
    complex(4), device :: c(*)
    integer, intent(in) :: m, n, ldc
    complex(4), intent(in) :: beta
    integer(cuda_stream_kind), intent(in) :: st
    integer :: i, j
    !$cuf kernel do(2) <<<*,*,stream=st>>>
    do j = 1, n
      do i = 1, m
        c(i + (j-1)*ldc) = beta*c(i + (j-1)*ldc)
      enddo
    enddo
  end subroutine scale_c

  integer function findslot(key, opa, m, k, alpha, ip) result(is)
    integer, intent(in) :: key, m, k, ip
    character, intent(in) :: opa
    complex(4), intent(in) :: alpha
    integer :: i
    is = 0
    do i = 1, nslot(ip)
      if (skey(i,ip) == key .and. sop(i,ip) == opa .and. sm(i,ip) == m .and. sk(i,ip) == k .and. salpha(i,ip) == alpha) then
        is = i
        return
      endif
    enddo
  end function findslot

  integer function newslot(key, opa, m, k, alpha, ip) result(is)
    !> A new slot in pool ip for this key, or 0 when the pool is full (the call then runs without keeping A').
    integer, intent(in) :: key, m, k, ip
    character, intent(in) :: opa
    complex(4), intent(in) :: alpha
    integer(8) :: need
    is = 0
    if (pool_cap(ip) < 0) call init_pool(ip)
    need = 4_8*k*m
    if (nslot(ip) >= maxslot .or. pool_used(ip) + need > pool_cap(ip)) return
    nslot(ip) = nslot(ip) + 1
    is = nslot(ip)
    skey(is,ip) = key; sop(is,ip) = opa; sm(is,ip) = m; sk(is,ip) = k; salpha(is,ip) = alpha
    soff(is,ip) = pool_used(ip)
    pool_used(ip) = pool_used(ip) + need
  end function newslot

  subroutine init_pool(ip)
    !> Pool ip32: ECALJ_LA_CACHE_GB (default 4) GB, at most 1/4 of the free device memory.  ip16: half of that.
    integer, intent(in) :: ip
    character(32) :: cv
    integer :: st, ios, istat
    real(8) :: gb, bytes
    integer(8) :: free, total
    gb = 4d0
    call get_environment_variable('ECALJ_LA_CACHE_GB', cv, status=st)
    if (st == 0) then
      read(cv,*,iostat=ios) gb
      if (ios /= 0) gb = 4d0
    endif
    istat = cudaMemGetInfo(free, total)
    bytes = min(gb*1d9, 0.25d0*real(free,8))
    if (ip == ip16) bytes = 0.5d0*bytes
    pool_cap(ip) = int(bytes/merge(4d0, 2d0, ip == ip32), 8)      ! in elements
    if (pool_cap(ip) < 1) then
      pool_cap(ip) = 0
      return
    endif
    if (ip == ip32) allocate(pool(pool_cap(ip)), stat=istat)
    if (ip == ip16) allocate(poolh(pool_cap(ip)), stat=istat)
    if (istat /= 0) pool_cap(ip) = 0
  end subroutine init_pool

  subroutine realsgemm_reset()
    !> Forget the kept A' of both routes (call when the matrices behind the keys change, e.g. at the next q point).
    nslot = 0
    pool_used = 0
  end subroutine realsgemm_reset
#endif
end module m_la_realsgemm
