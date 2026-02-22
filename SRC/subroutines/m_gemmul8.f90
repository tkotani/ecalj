module m_gemmul8
  use iso_c_binding
  implicit none
  integer, parameter :: num_moduli_d = 15, num_moduli_z = 15, num_moduli_c = 7
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
    logical :: cmdopt0
    if(is_gemmul8_inited) return
    use_gemmul8 = cmdopt0('--use_gemmul8')
    is_gemmul8_inited = .true.
#ifdef __GEMMUL8
    if(use_gemmul8) call gemmul8_init_handle(gemmul8_handle)
#endif
    if(use_gemmul8 .and. ipr) write(stdo,ftox), 'Using gemmul8 for GPU matrix multiplication'
#ifndef __GEMMUL8
    if(use_gemmul8) call rx0('Error: gemmul8 library is not linked.')
#endif
    istat = 0
  end function gemmul8_init
end module m_gemmul8
