module m_blas !wrapper for BLAS and cuBLAS
  !$use omp_lib
#ifdef __GPU
  use cublas_v2
  use cudafor
#endif
  implicit none
  include "mpif.h"
  public :: cmm, cmm_batch, zmm, zmm_batch, m_op_n, m_op_t, m_op_c
  public :: int_split
#ifdef __GPU
  type(cublashandle), value :: cublas_handle
  logical, save :: set_cublas_handle = .false.
#endif
  character, parameter :: m_op_n = 'N', m_op_t = 'T', m_op_c = 'C'
  interface 
    module function cmm(a, b, c, m, n, k, opa, opb, alpha, beta, lda, ldb, ldc) result(istat)
      complex(kind=4) :: a(*), b(*), c(*)
      integer, intent(in) :: m, n, k
      character, intent(in), optional :: opa, opb
      complex(kind=4), intent(in), optional :: alpha, beta
      integer, optional :: lda, ldb, ldc
      integer :: istat
#ifdef __GPU
      attributes(device) :: a, b, c
#endif
    end function
    module function cmm_batch(a, b, c, m, n, k, nbatch, opa, opb, alpha, beta, lda, ldb, ldc, samea, sameb, comm) result(istat)
      complex(kind=4) :: a(*), b(*), c(*)
      integer, intent(in) :: m, n, k, nbatch
      character, intent(in), optional :: opa, opb
      complex(kind=4), intent(in), optional :: alpha, beta
      integer, optional :: lda, ldb, ldc
      logical, optional :: samea, sameb
      integer, optional :: comm
      integer :: istat
#ifdef __GPU
      attributes(device) :: a, b, c
#endif
    end function
    module function zmm(a, b, c, m, n, k, opa, opb, alpha, beta, lda, ldb, ldc) result(istat)
      complex(kind=8) :: a(*), b(*), c(*)
      integer, intent(in) :: m, n, k
      character, intent(in), optional :: opa, opb
      complex(kind=8), intent(in), optional :: alpha, beta
      integer, optional :: lda, ldb, ldc
      integer :: istat
#ifdef __GPU
      attributes(device) :: a, b, c
#endif
    end function
    module function zmm_batch(a, b, c, m, n, k, nbatch, opa, opb, alpha, beta, lda, ldb, ldc, samea, sameb, comm) result(istat)
      complex(kind=8) :: a(*), b(*), c(*)
      integer, intent(in) :: m, n, k, nbatch
      character, intent(in), optional :: opa, opb
      complex(kind=8), intent(in), optional :: alpha, beta
      integer, optional :: lda, ldb, ldc
      logical, optional :: samea, sameb
      integer, optional :: comm
      integer :: istat
#ifdef __GPU
      attributes(device) :: a, b, c
#endif
    end function
  end interface
  private
contains
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