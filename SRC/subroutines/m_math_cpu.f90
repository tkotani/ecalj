module m_math_cpu
  !$use omp_lib
  implicit none
  public :: zmm, zmm_sb, mm_op_n, mm_op_t, mm_op_c
  character, parameter :: mm_op_n = 'N'
  character, parameter :: mm_op_t = 'T'
  character, parameter :: mm_op_c = 'C'
  private

  contains
  function zmm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) result(istat)
    integer :: istat
    character, value :: transa, transb
    integer, value :: m, n, k
    integer, intent(in), value :: lda, ldb, ldc
    complex(8) :: alpha, beta
    complex(8), dimension(lda,*) :: a
    complex(8), dimension(ldb,*) :: b
    complex(8), dimension(ldc,*) :: c
    call zgemm3m(transa, transb,  m, n, k, alpha, a, lda , b, ldb, beta, c, ldc)
    istat = 0
  end function zmm

  function cmm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) result(istat)
    integer :: istat
    character, value :: transa, transb
    integer, value :: m, n, k
    integer, intent(in), value :: lda, ldb, ldc
    complex(4) :: alpha, beta
    complex(4), dimension(lda,*) :: a
    complex(4), dimension(ldb,*) :: b
    complex(4), dimension(ldc,*) :: c
    call cgemm3m(transa, transb,  m, n, k, alpha, a, lda , b, ldb, beta, c, ldc)
    istat = 0
  end function cmm

  function zmm_sb(transa, transb, m, n, k, alpha, a, lda, stridea, b, ldb, strideb, beta, c, ldc, stridec, nbatch) &
     & result(istat)
    integer :: istat
    character, value :: transa, transb
    integer, value :: m, n, k, nbatch
    integer, intent(in), value :: lda, ldb, ldc
    integer(8) :: stridea, strideb, stridec
    complex(8) :: alpha, beta
    complex(8), dimension(lda,*) :: a
    complex(8), dimension(ldb,*) :: b
    complex(8), dimension(ldc,*) :: c
    integer :: i 
    !!$omp parallel do shared(a, b, c)
    do i=1, nbatch
      call zgemm3m(transa, transb, m, n, k, alpha, a(stridea*(i-1)+1,1), lda, &
                 & b(strideb*(i-1)+1,1), ldb, beta, c(stridec*(i-1)+1,1), ldc)
    enddo
    !!$omp end parallel do
    istat = 0
  end function zmm_sb
end module m_math_cpu
