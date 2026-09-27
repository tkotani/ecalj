module m_gemmul8
  use iso_c_binding
  implicit none
  ! Moduli and scaling mode of the Ozaki-II emulation (INT8).  Measured on RTX 5090 for the Sigma_c
  ! shape m=k=1053, n=49928 (TOOLS/ozbench, 2026-09-27; cuBLAS FP32 cgemm 31 TFLOPS, err 1.3e-6):
  !   cgemm 7 fast 43 TFLOPS err 4e-7, 8 fast 38 TFLOPS err 4e-8;  zgemm 14 fast 21 TFLOPS err 4e-15.
  ! Which products go to GEMMul8 is decided by m_linalg_policy (small products run slower emulated).
  ! ECALJ_GEMMUL8_MODULI_C / _Z / _D and ECALJ_GEMMUL8_FAST (0/1) override these for tests.
  integer :: num_moduli_d = 14, num_moduli_z = 14, num_moduli_c = 7, fastmode_gemmul8 = 1
  type(c_ptr) :: gemmul8_handle
#ifdef __GEMMUL8
  interface
    subroutine gemmul8_init_handle(handle) bind(C, name="gemmul8_init_handle_")
      use iso_c_binding
      type(c_ptr) :: handle
    end subroutine
    subroutine gemmul8_finalize_handle(handle) bind(C, name="gemmul8_finalize_handle_")
      use iso_c_binding
      type(c_ptr), value :: handle
    end subroutine
    subroutine gemmul8_zgemm(handle, transa, transb, m, n, k, alpha, devA, lda, devB, ldb, beta, devC, ldc, &
                             num_moduli, fastmode, key) bind(C, name="gemmul8_zgemm_")
      use iso_c_binding
      type(c_ptr), value :: handle
      integer, value :: transa, transb, m, n, k, lda, ldb, ldc
      complex(8), value :: alpha, beta
      type(c_ptr), value :: devA, devB, devC
      integer, value :: num_moduli
      integer, value :: fastmode, key   ! key >= 0: A is kept and reused under this key (gemmul8_cache_reset)
    endsubroutine
    subroutine gemmul8_dgemm(handle, transa, transb, m, n, k, alpha, devA, lda, devB, ldb, beta, devC, ldc, &
                             num_moduli, fastmode, key) bind(C, name="gemmul8_dgemm_")
      use iso_c_binding
      type(c_ptr), value :: handle
      integer, value :: transa, transb, m, n, k, lda, ldb, ldc
      real(8), value :: alpha, beta
      type(c_ptr), value :: devA, devB, devC
      integer, value :: num_moduli
      integer, value :: fastmode, key   ! key >= 0: A is kept and reused under this key (gemmul8_cache_reset)
    endsubroutine
    subroutine gemmul8_cgemm(handle, transa, transb, m, n, k, alpha, devA, lda, devB, ldb, beta, devC, ldc, &
                             num_moduli, fastmode, key) bind(C, name="gemmul8_cgemm_")
      use iso_c_binding
      type(c_ptr), value :: handle
      integer, value :: transa, transb, m, n, k, lda, ldb, ldc
      complex(4), value :: alpha, beta
      type(c_ptr), value :: devA, devB, devC
      integer, value :: num_moduli
      integer, value :: fastmode, key   ! key >= 0: A is kept and reused under this key (gemmul8_cache_reset)
    endsubroutine
    subroutine gemmul8_cache_reset() bind(C, name="gemmul8_cache_reset_")
    end subroutine
    subroutine gemmul8_set_stream(handle, stream) bind(C, name="gemmul8_set_stream_")
      use iso_c_binding
      type(c_ptr), value :: handle
      integer(c_intptr_t), value :: stream   ! a CUDA stream (cuda_stream_kind)
    end subroutine
  endinterface
#endif
contains
  logical function gemmul8_pays(m, n, k)
    !> GEMMul8 splits A and B into moduli and runs one INT8 GEMM per modulus: a fixed cost per call that a small
    !> product does not pay back.  The 7x7 x nband rotations in readeigen (thousands per k point) cost hgw 110 s with
    !> GEMMul8 against cuBLAS (LiTi2O4 6^3, 2026-09-27).  Below this size m_blas uses cuBLAS whatever the table says.
    integer, intent(in) :: m, n, k
    gemmul8_pays = min(m, n, k) >= 64 .and. real(m,8)*real(n,8)*real(k,8) >= 1d8
  end function gemmul8_pays
  integer function gemmul8_init() result(istat)
    !> Create the GEMMul8 handle at the first product routed to it (m_linalg_policy decides which ones).
    use m_lgunit, only: stdo
    use m_ftox, only: ftox
    use m_mpi,only: ipr
    logical, save :: is_gemmul8_inited = .false.
    istat = 0
    if(is_gemmul8_inited) return
    is_gemmul8_inited = .true.
    call envint('ECALJ_GEMMUL8_MODULI_C', num_moduli_c)
    call envint('ECALJ_GEMMUL8_MODULI_Z', num_moduli_z)
    call envint('ECALJ_GEMMUL8_MODULI_D', num_moduli_d)
    call envint('ECALJ_GEMMUL8_FAST',     fastmode_gemmul8)
#ifdef __GEMMUL8
    call gemmul8_init_handle(gemmul8_handle)
    if(ipr) write(stdo,ftox) 'gemmul8: default moduli c z d =', &
         num_moduli_c, num_moduli_z, num_moduli_d, 'fastmode=', fastmode_gemmul8
#else
    call rx0('Error: the policy chose gemmul8 but the GEMMul8 library is not linked (InstallAll.py --gemmul8).')
#endif
  contains
    subroutine envint(name, val)
      character(*), intent(in) :: name
      integer, intent(inout) :: val
      character(32) :: cv
      integer :: st, ios, v
      call get_environment_variable(name, cv, status=st)
      if (st /= 0) return
      read(cv, *, iostat=ios) v
      if (ios == 0) val = v
    end subroutine envint
  end function gemmul8_init
end module m_gemmul8
