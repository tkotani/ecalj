module m_gemmul8
  use iso_c_binding
  use m_cmdopt_registry, only: c0_use_gemmul8
  implicit none
  ! Moduli and scaling mode of the Ozaki-II emulation (INT8).  Measured on RTX 5090 for the Sigma_c
  ! shape m=k=1053, n=49928 (TOOLS/ozbench, 2026-09-27; cuBLAS FP32 cgemm 31 TFLOPS, err 1.3e-6):
  !   cgemm 7 fast 43 TFLOPS err 4e-7, 8 fast 38 TFLOPS err 4e-8;  zgemm 14 fast 21 TFLOPS err 4e-15.
  ! Smaller products (n of a few hundred, or m=n=158) run slower than cuBLAS, hence gemmul8_worth.
  ! ECALJ_GEMMUL8_MODULI_C / _Z / _D and ECALJ_GEMMUL8_FAST (0/1) override these for tests.
  integer :: num_moduli_d = 14, num_moduli_z = 14, num_moduli_c = 7, fastmode_gemmul8 = 1
  logical :: use_gemmul8 = .false.
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
                             num_moduli, fastmode, enable_skip_A, enable_skip_B, skip_scalA, skip_scalB) bind(C, name="gemmul8_zgemm_")
      use iso_c_binding
      type(c_ptr), value :: handle
      integer, value :: transa, transb, m, n, k, lda, ldb, ldc
      complex(8), value :: alpha, beta
      type(c_ptr), value :: devA, devB, devC
      integer, value :: num_moduli
      integer, value :: fastmode, enable_skip_A, enable_skip_B, skip_scalA, skip_scalB
    endsubroutine
    subroutine gemmul8_dgemm(handle, transa, transb, m, n, k, alpha, devA, lda, devB, ldb, beta, devC, ldc, &
                             num_moduli, fastmode, enable_skip_A, enable_skip_B, skip_scalA, skip_scalB) bind(C, name="gemmul8_dgemm_")
      use iso_c_binding
      type(c_ptr), value :: handle
      integer, value :: transa, transb, m, n, k, lda, ldb, ldc
      real(8), value :: alpha, beta
      type(c_ptr), value :: devA, devB, devC
      integer, value :: num_moduli
      integer, value :: fastmode, enable_skip_A, enable_skip_B, skip_scalA, skip_scalB
    endsubroutine
    subroutine gemmul8_cgemm(handle, transa, transb, m, n, k, alpha, devA, lda, devB, ldb, beta, devC, ldc, &
                             num_moduli, fastmode, enable_skip_A, enable_skip_B, skip_scalA, skip_scalB) bind(C, name="gemmul8_cgemm_")
      use iso_c_binding
      type(c_ptr), value :: handle
      integer, value :: transa, transb, m, n, k, lda, ldb, ldc
      complex(4), value :: alpha, beta
      type(c_ptr), value :: devA, devB, devC
      integer, value :: num_moduli
      integer, value :: fastmode, enable_skip_A, enable_skip_B, skip_scalA, skip_scalB
    endsubroutine
  endinterface
#endif
contains
  integer function gemmul8_init() result(istat)
    use m_lgunit, only: stdo
    use m_ftox, only: ftox
    use m_mpi,only: ipr
    logical, save :: is_gemmul8_inited = .false.
    if(is_gemmul8_inited) return
    use_gemmul8 = c0_use_gemmul8
    is_gemmul8_inited = .true.
    call envint('ECALJ_GEMMUL8_MODULI_C', num_moduli_c)
    call envint('ECALJ_GEMMUL8_MODULI_Z', num_moduli_z)
    call envint('ECALJ_GEMMUL8_MODULI_D', num_moduli_d)
    call envint('ECALJ_GEMMUL8_FAST',     fastmode_gemmul8)
#ifdef __GEMMUL8
    if(use_gemmul8) call gemmul8_init_handle(gemmul8_handle)
#endif
    if(use_gemmul8 .and. ipr) write(stdo,ftox) 'Using gemmul8 for large GPU matrix products: moduli c z d =', &
         num_moduli_c, num_moduli_z, num_moduli_d, 'fastmode=', fastmode_gemmul8
#ifndef __GEMMUL8
    if(use_gemmul8) call rx0('Error: gemmul8 library is not linked.')
#endif
    istat = 0
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
  logical function gemmul8_worth(m, n, k)
    !> Emulation pays off only for large products: all three sizes >= 1000 and m*n*k >= 1e10.
    integer, intent(in) :: m, n, k
    gemmul8_worth = min(m, n, k) >= 1000 .and. real(m,8)*real(n,8)*real(k,8) >= 1d10
  end function gemmul8_worth
end module m_gemmul8
