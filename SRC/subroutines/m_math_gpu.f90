module m_math_gpu
  use cublas_v2
  implicit none
  public :: zmm, mm_op_n, mm_op_t, mm_op_c
  integer, parameter :: mm_op_n = cublas_op_n 
  integer, parameter :: mm_op_t = cublas_op_t 
  integer, parameter :: mm_op_c = cublas_op_c
  private
  type(cublashandle), value :: handle
  logical, save :: tfirst = .true.

  contains
  function zmm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) result(istat)
    integer :: istat
    integer, value :: transa, transb
    integer, value :: m, n, k
    integer, intent(in), value :: lda, ldb, ldc
    complex(8) :: alpha, beta
    complex(8), dimension(lda,*), device :: a
    complex(8), dimension(ldb,*), device :: b
    complex(8), dimension(ldc,*), device :: c
    if(tfirst) then
      istat = cublascreate(handle)
      tfirst = .false.
    endif
    istat = cublaszgemm(handle, transa, transb,  m, n, k, alpha, a, lda , b, ldb, beta, c, ldc)
  end function zmm
end module m_math_gpu
