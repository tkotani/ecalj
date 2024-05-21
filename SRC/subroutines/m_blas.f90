module m_blas !wrapper for BLAS and cuBLAS
  !$use omp_lib
#ifdef __GPU
  use cublas_v2
  use cudafor
#endif
  implicit none
  include "mpif.h"
  public :: zmm, zmm_batch, m_op_n, m_op_t, m_op_c
  public :: int_split
#ifdef __GPU
  type(cublashandle), value :: cublas_handle
  logical, save :: set_cublas_handle = .false.
#endif
  character, parameter :: m_op_n = 'N', m_op_t = 'T', m_op_c = 'C'
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

  function zmm(a, b, c, m, n, k, opa, opb, alpha, beta, lda, ldb, ldc) result(istat)
    complex(8) :: a(*), b(*), c(*)
    integer, intent(in) :: m, n, k
    character, intent(in), optional :: opa, opb
    complex(8), intent(in), optional :: alpha, beta
    integer, optional :: lda, ldb, ldc

    complex(8) :: alpha_in, beta_in
    integer :: lda_in, ldb_in, ldc_in, istat
    character :: opa_in, opb_in
#ifdef __GPU
    attributes(device) :: a, b, c
#endif
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

#ifdef __GPU
    cublas_zgemm3m: block 
      integer :: opa_in_cublas, opb_in_cublas
      istat = cublas_init()
      opa_in_cublas = get_m_op_cublas(opa_in)
      opb_in_cublas = get_m_op_cublas(opb_in)
      istat = cublaszgemm3m(cublas_handle, opa_in_cublas, opb_in_cublas,  m, n, k, &
                          & alpha_in, a, lda_in , b, ldb_in, beta_in, c, ldc_in)
    end block cublas_zgemm3m
#else
    call zgemm3m(opa_in, opb_in, m, n, k, alpha_in, a, lda_in, b, ldb_in, beta_in, c, ldc_in)
    istat = 0
#endif
  end function zmm

  function zmm_batch(a, b, c, m, n, k, nbatch, opa, opb, alpha, beta, lda, ldb, ldc, samea, sameb, comm) result(istat)
    complex(8) :: a(*), b(*), c(*)
    integer, intent(in) :: m, n, k, nbatch
    character, intent(in), optional :: opa, opb
    complex(8), intent(in), optional :: alpha, beta
    integer, optional :: lda, ldb, ldc
    logical, optional :: samea, sameb
    integer, optional :: comm

    integer(8) :: stridea, strideb, stridec
    complex(8) :: alpha_in, beta_in
    integer :: lda_in, ldb_in, ldc_in, istat
    character :: opa_in, opb_in
#ifdef __GPU
    attributes(device) :: a, b, c
#endif

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

#ifdef __GPU
    cublas_zgemmstridedbatched: block
      integer :: opa_in_cublas, opb_in_cublas
      istat = cublas_init()
      opa_in_cublas = get_m_op_cublas(opa_in)
      opb_in_cublas = get_m_op_cublas(opb_in)
      istat = cublaszgemmstridedbatched(cublas_handle, opa_in_cublas, opb_in_cublas,  m, n, k,  &
                 &  alpha_in, a, lda_in, stridea, b, ldb_in, strideb, beta_in, c, ldc_in, stridec, nbatch)
    end block cublas_zgemmstridedbatched
#else
    blas_zgemmbatch: block
      integer :: i
      integer :: ini_batch, end_batch, mpi_size, mpi_rank, ierr, nbatch_irank, irank
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
    end block blas_zgemmbatch
#endif
  end function zmm_batch

  subroutine int_split(ndata, nsplit, irank, iini, iend, n)
    integer, intent(in) :: ndata, nsplit, irank
    integer, intent(out) :: iini, iend, n
    n = (ndata + irank)/nsplit
    iini = (ndata/nsplit)*irank + max(irank + mod(ndata, nsplit) - nsplit, 0) + 1  
    iend = iini + n - 1
  end subroutine
end module m_blas
