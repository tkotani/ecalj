module m_lapack
#ifdef __GPU
  use cublas_v2
  use cudafor
  use cusolverdn
#endif
  use m_blas
  use m_lgunit,      only: stdo
  implicit none
  public :: zhgv_h, zhgv_lindep_h, zhev_h, zminv_h, zsv_h, zgev_h, zggv_h
#ifdef __GPU
  public :: zhgv_d, zhev_d, zminv_d, zsv_d, cusolver_finalize !, zgev_d
#endif
  private
#ifdef __GPU
  integer :: cuda_runtime_version
  type(cusolverDnHandle), value :: cusolver_handle
  type(cusolverDnParams), value :: cusolver_params
  logical, save :: set_cusolver_handle = .false.
#endif
contains
  integer function zminv_h(a, n, lda) result(istat)
    complex(8) :: a(*)
    integer, intent(in) :: n
    integer, optional :: lda
    integer :: lda_in
    integer, allocatable :: ipvt(:)
    complex(8),allocatable:: work(:)
    integer :: lwork, original_num_threads
    complex(8) :: wkopt
    lda_in = n; if(present(lda)) lda_in = lda
#ifdef __MKL_ZMINV_SEQUENTIAL
    original_num_threads = mkl_get_max_threads()
    call mkl_set_num_threads(1)
#endif
    allocate(ipvt(n))
    call zgetrf(n,n,a,lda_in,ipvt,istat)
    lwork = -1
    call zgetri(n, a, lda_in, ipvt, wkopt, lwork, istat)
    lwork = int(dble(wkopt))
    allocate(work(lwork))
    call zgetri(n, a, lda_in, ipvt, work, lwork, istat)
#ifdef __MKL_ZMINV_SEQUENTIAL
    call mkl_set_num_threads(original_num_threads)
#endif
    deallocate(work,ipvt)
  end function zminv_h
  integer function zsv_h(a, b, n, nrhs, lda, ldb) result(istat)
    integer, intent(in) :: n, nrhs
    complex(8) :: a(*), b(*)
    integer, allocatable :: ipiv(:)
    integer, intent(in), optional :: lda, ldb
    integer :: lda_in, ldb_in
    allocate(ipiv(n))
    lda_in = n; ldb_in = n
    if(present(lda)) lda_in = lda
    if(present(ldb)) ldb_in = ldb
    call zgetrf(n, n, a, lda_in, ipiv, istat)
    if(istat /= 0) return
    call zgetrs('N', n, nrhs, a, lda_in, ipiv, b, ldb_in, istat)
    deallocate(ipiv)
  end function zsv_h
  integer function zgev_h(a, n, evl, evl_vec, evr_vec, lda) result(istat)
    integer, intent(in) :: n
    complex(8) :: a(*)
    complex(8), intent(out) :: evl(n)
    complex(8), intent(out), optional :: evl_vec(*), evr_vec(*)
    integer, intent(in), optional :: lda
    integer :: lda_in, ldvl, ldvr, lwork, ilo, ihi
    complex(8), allocatable :: work(:), vl_loc(:), vr_loc(:)
    real(8),    allocatable :: rwork(:), scale(:), rconde(:), rcondv(:)
    real(8) :: abnrm
    character :: jobvl, jobvr
    lda_in = n; if(present(lda)) lda_in = lda
    jobvl = 'N'; ldvl = 1
    jobvr = 'N'; ldvr = 1
    if(present(evl_vec)) then; jobvl = 'V'; ldvl = lda_in; endif
    if(present(evr_vec)) then; jobvr = 'V'; ldvr = lda_in; endif
    allocate(vl_loc(ldvl*n), vr_loc(ldvr*n), rwork(2*n), work(1), scale(n), rconde(n), rcondv(n))
    lwork = -1
    call zgeevx('B', jobvl, jobvr, 'N', n, a, lda_in, evl, &
                vl_loc, ldvl, vr_loc, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, &
                work, lwork, rwork, istat)
    lwork = int(dble(work(1))); deallocate(work); allocate(work(lwork))
    call zgeevx('B', jobvl, jobvr, 'N', n, a, lda_in, evl, &
                vl_loc, ldvl, vr_loc, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, &
                work, lwork, rwork, istat)
    if(present(evl_vec)) evl_vec(1:ldvl*n) = vl_loc
    if(present(evr_vec)) evr_vec(1:ldvr*n) = vr_loc
    deallocate(work, rwork, vl_loc, vr_loc, scale, rconde, rcondv)
  end function zgev_h
  integer function zggv_h(a, b, n, evl, evl_vec, evr_vec, lda, ldb, alpha, beta) result(istat)
  ! Solving the generalized eigenvalue problem Az = lambda Bz for general matrices
  ! Eigenvalues lambda(i) = alpha(i)/beta(i) are stored in evl
  ! Left eigenvectors stored in evl_vec (optional), right eigenvectors in evr_vec (optional)
  ! alpha, beta are available as optional outputs (e.g. to detect infinite eigenvalues where beta=0)
    integer, intent(in) :: n
    complex(8) :: a(*), b(*)
    complex(8), intent(out) :: evl(n)
    complex(8), intent(out), optional :: evl_vec(*), evr_vec(*), alpha(n), beta(n)
    integer, intent(in), optional :: lda, ldb
    integer :: lda_in, ldb_in, ldvl, ldvr, lwork, i
    complex(8), allocatable :: work(:), vl_loc(:), vr_loc(:), alpha_loc(:), beta_loc(:)
    real(8), allocatable :: rwork(:)
    character :: jobvl, jobvr
    lda_in = n; if(present(lda)) lda_in = lda
    ldb_in = n; if(present(ldb)) ldb_in = ldb
    jobvl = 'N'; ldvl = 1
    jobvr = 'N'; ldvr = 1
    if(present(evl_vec)) then; jobvl = 'V'; ldvl = lda_in; endif
    if(present(evr_vec)) then; jobvr = 'V'; ldvr = lda_in; endif
    allocate(vl_loc(ldvl*n), vr_loc(ldvr*n), rwork(8*n), work(1), alpha_loc(n), beta_loc(n))
    lwork = -1
    call zggev(jobvl, jobvr, n, a, lda_in, b, ldb_in, alpha_loc, beta_loc, &
               vl_loc, ldvl, vr_loc, ldvr, work, lwork, rwork, istat)
    lwork = int(dble(work(1))); deallocate(work); allocate(work(lwork))
    call zggev(jobvl, jobvr, n, a, lda_in, b, ldb_in, alpha_loc, beta_loc, &
               vl_loc, ldvl, vr_loc, ldvr, work, lwork, rwork, istat)
    do i = 1, n
      if(abs(beta_loc(i)) == 0d0) then
        evl(i) = cmplx(huge(1d0), 0d0, 8)
      else
        evl(i) = alpha_loc(i) / beta_loc(i)
      endif
    enddo
    if(present(evl_vec)) evl_vec(1:ldvl*n) = vl_loc
    if(present(evr_vec)) evr_vec(1:ldvr*n) = vr_loc
    if(present(alpha)) alpha(1:n) = alpha_loc
    if(present(beta))  beta(1:n)  = beta_loc
    deallocate(work, rwork, vl_loc, vr_loc, alpha_loc, beta_loc)
  end function zggv_h
  integer function zhev_h(A, n, evl, il, iu, lda) result(istat)
  ! Solving the standard eigenvalue problem Az = lambda z, where A is a Hermitian matrix
  ! Eigenvalues are stored in evl, eigenvectors are stored in A
  !!! range: 1<=il<=iu<=n
    integer, intent(in) :: n
    real(8), intent(out) :: evl(n)
    complex(8) :: A(*)
    integer, intent(in), optional :: lda, il, iu
    integer :: lda_in, il_in, iu_in
    complex(8), allocatable :: work(:), z(:)
    real(8), allocatable :: rwork(:)
    integer, allocatable :: isuppz(:), iwork(:)
    integer :: m, lwork, lrwork, liwork, info
    real(8) :: abstol, vl, vu, dlamch
    lda_in = n; if(present(lda)) lda_in = lda
    il_in = 1; iu_in = n
    if(present(il)) il_in = il
    if(present(iu)) iu_in = iu
    vl = 0d0; vu = 0d0; abstol = 2d0*dlamch('S')
    allocate(z(lda_in*n), isuppz(2*n))
    lwork = -1; lrwork = -1; liwork = -1
    allocate(work(1), rwork(1), iwork(1))
    call zheevr('V', 'I', 'U', n, a, lda_in, vl, vu, il_in, iu_in, abstol, m, evl, z, lda_in, isuppz, &
                work, lwork, rwork, lrwork, iwork, liwork, info)
    lwork = int(dble(work(1))); lrwork = int(rwork(1)); liwork = iwork(1)
    deallocate(work, rwork, iwork)
    allocate(work(lwork), rwork(lrwork), iwork(liwork))
    call zheevr('V', 'I', 'U', n, a, lda_in, vl, vu, il_in, iu_in, abstol, m, evl, z, lda_in, isuppz, &
                work, lwork, rwork, lrwork, iwork, liwork, info)
    istat = info
    a(1:lda_in*n) = z(1:lda_in*n)
    deallocate(work, rwork, iwork, z, isuppz)
  end function zhev_h
  integer function zhgv_h(A, B, n, evl, il, iu, lda, ldb, keep_ab, evec) result(istat)
  ! Solving the generalized eigenvalue problem Az = lambda Bz, where A, B are Hermitian matrixes, z is eigenfunction
  ! Eigenvalues are stores in evl, eigenvectors are stored in A (or evec if present)
  ! keep_ab=.true.: A and B are not overwritten (internal copies used for LAPACK call)
  ! evec present: eigenvectors are stored in evec instead of A
  ! Use keep_ab=.true. with evec to make the call fully non-destructive
  !!! that range is 1<=IL <= IU <= N
    integer, intent(in) :: n !size of matrix
    real(8), intent(out) :: evl(n) !eigenvalues
    complex(8) :: A(*), B(*)
    integer, intent(in), optional :: lda, ldb, il, iu
    logical, intent(in), optional :: keep_ab
    complex(8), intent(out), optional :: evec(*)
    integer :: lda_in, ldb_in, il_in, iu_in
    integer, parameter :: nb = 64
    complex(8), allocatable:: work(:), z(:), a_work(:), b_work(:)
    integer, allocatable:: ifail(:), iwork(:)
    real(8), allocatable:: rwork(:)
    integer :: m, lwork, info
    real(8) :: dlamch, abstol, vl = 0d0, vu = 0d0
    logical :: keep_ab_in
    lda_in = n; ldb_in = n
    if(present(lda)) lda_in = lda
    if(present(ldb)) ldb_in = ldb
    il_in = 1; iu_in = n
    if(present(il)) il_in = il
    if(present(iu)) iu_in = iu
    keep_ab_in = .false.; if(present(keep_ab)) keep_ab_in = keep_ab
    abstol = 2d0*dlamch('S')
    allocate(z(lda_in*n), source = (0d0, 0d0))
    allocate(work(1))
    allocate(rwork(7*n), ifail(n), iwork(5*n))
    if(keep_ab_in) then
      allocate(a_work(lda_in*n), b_work(ldb_in*n))
      a_work(1:lda_in*n) = a(1:lda_in*n)
      b_work(1:ldb_in*n) = b(1:ldb_in*n)
      lwork = -1
      call zhegvx( 1, 'V', 'I', 'U', n, a_work, lda_in, b_work, ldb_in, &
         vl, vu, il_in, iu_in, abstol, m, evl, z, lda_in, &
         work, lwork, rwork, iwork, ifail, info )
      lwork = max((nb+1)*n, nint(dble(work(1))))
      deallocate(work); allocate(work(lwork))
      call zhegvx( 1, 'V', 'I', 'U', n, a_work, lda_in, b_work, ldb_in, &
         vl, vu, il_in, iu_in, abstol, m, evl, z, lda_in, &
         work, lwork, rwork, iwork, ifail, info)
      deallocate(a_work, b_work)
    else
      lwork = -1
      call zhegvx( 1, 'V', 'I', 'U', n, a, lda_in, b, ldb_in, &
         vl, vu, il_in, iu_in, abstol, m, evl, z, lda_in, &
         work, lwork, rwork, iwork, ifail, info )
      lwork = max((nb+1)*n, nint(dble(work(1))))
      deallocate(work); allocate(work(lwork))
      call zhegvx( 1, 'V', 'I', 'U', n, a, lda_in, b, ldb_in, &
         vl, vu, il_in, iu_in, abstol, m, evl, z, lda_in, &
         work, lwork, rwork, iwork, ifail, info)
    endif
    istat = info
    if(info < 0) then
      write(*,'(a,i0,a)') 'zhgv_h: zhegvx error: illegal value in argument ', -info, '. Aborting.'
      ! stop
    elseif(info > n) then
      write(*,'(a,i0,a)') 'zhgv_h: zhegvx error: B is not positive definite (Cholesky failed at minor ', info-n, '). Aborting.'
      ! stop
    elseif(info > 0) then
      write(*,'(a,i0,a)') 'zhgv_h: zhegvx warning: ', info, ' eigenvector(s) failed to converge (see ifail).'
    endif
    if(present(evec)) then
      evec(1:lda_in*n) = z(1:lda_in*n)
    else
      a(1:lda_in*n) = z(1:lda_in*n)
    endif
    deallocate(work,rwork,iwork,ifail,z)
  end function zhgv_h
  integer function zhgv_lindep_h(A, B, n, evl, il, iu, lda, ldb, evec, thr_lo, nev, nkeep) result(istat)
  ! Solving the generalized eigenvalue problem Az = lambda Bz for Hermitian A, B
  ! Handles near-singular/indefinite B via canonical orthogonalization:
  !   1. Diagonalize B: eigenvectors V, eigenvalues b_evl (ascending)
  !   2. Discard eigenvectors with b_evl <= thr_lo  (linear dependencies)
  !   3. Form reduced basis X(:,j) = V(:,j) / sqrt(b_evl(j))  for kept vectors
  !   4. Solve reduced standard eigenvalue problem (X^H A X) y = lambda y
  !   5. Back-transform eigenvectors: z = X * y
  ! A and B are never overwritten (copies used internally)
  ! evec (optional out): eigenvectors; if absent, stored in A
  ! thr_lo (optional in): threshold for linear dependence (default 1e-1)
  ! nev (optional out): number of eigenvalues found (may be < iu-il+1 after reduction)
  ! nkeep (optional in): if present, keep only the nkeep largest B eigenvalues.
  !   When both thr_lo and nkeep are given, the stricter criterion (fewer vectors) is applied.
    integer, intent(in) :: n
    real(8), intent(out) :: evl(n)
    complex(8) :: A(*), B(*)
    integer, intent(in), optional :: lda, ldb, il, iu, nkeep
    complex(8), intent(out), optional :: evec(*)
    real(8), intent(in), optional :: thr_lo
    integer, intent(out), optional :: nev
    integer :: lda_in, ldb_in, il_in, iu_in, il_lo, iu_lo, n_lo, n_ev, i, j, i0
    real(8) :: thr_lo_in
    complex(8), allocatable :: b_evec(:), x(:), tmp(:), a_lo(:), z_out(:)
    real(8), allocatable :: b_evl(:), evl_lo(:)
    lda_in = n; ldb_in = n
    if(present(lda)) lda_in = lda
    if(present(ldb)) ldb_in = ldb
    il_in = 1; iu_in = n
    if(present(il)) il_in = il
    if(present(iu)) iu_in = iu
    thr_lo_in = 1d-8; if(present(thr_lo)) thr_lo_in = thr_lo
    ! Step 1: Diagonalize B (copy so B is not modified)
    allocate(b_evec(ldb_in*n), b_evl(n))
    b_evec(1:ldb_in*n) = b(1:ldb_in*n)
    istat = zhev_h(b_evec, n, b_evl, lda=ldb_in)
    if(istat /= 0) then
      write(*,'(a,i0)') 'zhgv_lindep_h: diagonalization of B failed, info=', istat
      return
    endif
    ! Step 2: Count linearly independent vectors (b_evl in ascending order)
    ! 2a: threshold criterion
    i0 = n + 1
    do i = 1, n
      if(b_evl(i) > thr_lo_in) then; i0 = i; exit; endif
    enddo
    ! 2b: nkeep criterion — keep at most nkeep largest eigenvalues
    if(present(nkeep)) i0 = max(i0, n - nkeep + 1)
    n_lo = n - i0 + 1
    if(n_lo == 0) then
      write(stdo,'(a,es10.3)') 'zhgv_lindep_h: no B eigenvalue above thr_lo=', thr_lo_in
      istat = n + 1
      if(present(nev)) nev = 0
      deallocate(b_evec, b_evl)
      return
    endif
    ! if(n_lo < n) write(stdo,'(a,i0,a,i0,a,es10.3)') &
    !   'zhgv_lindep_h: ', n-n_lo, ' linear dependence(s) removed, n_lo=', n_lo, ', thr=', thr_lo_in
    ! Step 3: Form X(n x n_lo): column j = b_evec_col(i0+j-1) / sqrt(b_evl(i0+j-1))
    allocate(x(n*n_lo))
    do j = 1, n_lo
      x((j-1)*n+1:(j-1)*n+n) = b_evec(ldb_in*(i0+j-2)+1:ldb_in*(i0+j-2)+n) / sqrt(b_evl(i0+j-1))
    enddo
    deallocate(b_evec, b_evl)
    ! Step 4: Form a_lo(n_lo x n_lo) = X^H * A * X
    ! tmp(n x n_lo) = A * X  (A is Hermitian, upper triangle used)
    allocate(tmp(n*n_lo))
    call zhemm('L', 'U', n, n_lo, (1d0,0d0), a(1), lda_in, x(1), n, (0d0,0d0), tmp(1), n)
    allocate(a_lo(n_lo*n_lo))
    call zgemm('C', 'N', n_lo, n_lo, n, (1d0,0d0), x(1), n, tmp(1), n, (0d0,0d0), a_lo(1), n_lo)
    deallocate(tmp)
    ! Step 5: Solve reduced standard eigenvalue problem, mapping il/iu to reduced size
    il_lo = il_in
    iu_lo = min(iu_in, n_lo)
    if(il_lo > iu_lo) then
      if(present(nev)) nev = 0
      istat = 0; deallocate(x, a_lo); return
    endif
    n_ev = iu_lo - il_lo + 1
    allocate(evl_lo(n_lo))
    istat = zhev_h(a_lo, n_lo, evl_lo, il=il_lo, iu=iu_lo)
    if(istat /= 0) then
      write(stdo,'(a,i0)') 'zhgv_lindep_h: diagonalization of reduced A failed, info=', istat
      deallocate(x, a_lo, evl_lo); return
    endif
    evl(1:n_ev) = evl_lo(1:n_ev)
    if(present(nev)) nev = n_ev
    deallocate(evl_lo)
    ! Step 6: Back-transform: z_out(n x n_ev) = X(n x n_lo) * a_lo_evec(n_lo x n_ev)
    ! After zhev_h, a_lo columns 1..n_ev hold the eigenvectors (leading dim n_lo)
    allocate(z_out(n*n_ev))
    call zgemm('N', 'N', n, n_ev, n_lo, (1d0,0d0), x(1), n, a_lo(1), n_lo, (0d0,0d0), z_out(1), n)
    deallocate(x, a_lo)
    ! Step 7: Output eigenvectors to evec or A (with proper lda stride)
    if(present(evec)) then
      do j = 1, n_ev
        evec(lda_in*(j-1)+1:lda_in*(j-1)+n) = z_out((j-1)*n+1:(j-1)*n+n)
      enddo
    else
      do j = 1, n_ev
        a(lda_in*(j-1)+1:lda_in*(j-1)+n) = z_out((j-1)*n+1:(j-1)*n+n)
      enddo
    endif
    deallocate(z_out)
  end function zhgv_lindep_h
#ifdef __GPU
  !cusolverDnXtrtri used in zminv has internal compiler bug before cuda 12.5
  !https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/index.html
  integer function zminv_d(a, n, lda) result(istat)
    complex(8), device :: a(*)
    integer, intent(in) :: n
    integer, optional :: lda
    integer :: lda_in
    complex(8), allocatable, device :: awork(:)
    integer(4), device :: devinfo
    integer(8), device :: ipiv(n)
    integer(8)         :: ipiv_cpu(n)
    integer(1), allocatable, device :: buffer_d(:)
    integer(1), allocatable         :: buffer_h(:)
    integer(8) :: n_8, lda_8, lbuffer_d, lbuffer_h
    integer(4) :: i, j
    lda_in = n; if(present(lda)) lda_in = lda

    n_8 = int(n,8)
    lda_8 = int(lda_in,8)
    istat = cublas_init()
    istat = cusolver_init()
    ! LU factorization PA = LU in-place on a, L: Lower triangular with unit diagonal components, U: Upper triangular with
    ! non-unit-diagonal,  P: permutation of rows stored in ipiv
    ! A^-1 = U^-1 * L^-1 * P
    istat = cusolverDnXgetrf_buffersize(cusolver_handle, cusolver_params, n_8, n_8, cudaDataType(CUDA_C_64F), a, lda_8, &
                                        cudaDataType(CUDA_C_64F), lbuffer_d, lbuffer_h)
    allocate(buffer_d(lbuffer_d), buffer_h(lbuffer_h))
    istat = cusolverDnXgetrf(cusolver_handle, cusolver_params, n_8, n_8, cudaDataType(CUDA_C_64F), a, lda_8, &
                             ipiv, cudaDataType(CUDA_C_64F), buffer_d, lbuffer_d, buffer_h, lbuffer_h, devinfo)
    deallocate(buffer_d, buffer_h)
    allocate(awork(lda_in*n))
    !$acc kernels
    awork(1:lda_in*n) = a(1:lda_in*n) !awork -> Used as a upper triangular matrix with diagonal part
    !$acc end kernels
    !$acc kernels loop collapse(2) independent
    do j = 1, n
      do i = 1, n
        if(i < j) a(lda_in*(j-1)+i) = (0d0, 0d0) !Used as lower triangular matrix, diagonal components will be replaced as 1 in the following function
      enddo
    enddo
    !$acc end kernels
    !Get inverse of L = L^-1
    istat = cusolverDnXtrtri_buffersize(cusolver_handle, CUBLAS_FILL_MODE_LOWER, CUBLAS_DIAG_UNIT, n_8,  &
                                        cudaDataType(CUDA_C_64F), a, lda_8, lbuffer_d, lbuffer_h)
    if(cuda_runtime_version <= 12040)  then !prescription for bug of cusolverDnXtrtri_buffersize
      lbuffer_d = lbuffer_d*16; lbuffer_h = lbuffer_h*16
    endif
    allocate(buffer_d(lbuffer_d), buffer_h(lbuffer_h))
    istat = cusolverDnXtrtri(cusolver_handle, CUBLAS_FILL_MODE_LOWER, CUBLAS_DIAG_UNIT, n_8, cudaDataType(CUDA_C_64F), a, lda_8, &
                             buffer_d, lbuffer_d, buffer_h, lbuffer_h, devinfo)
    deallocate(buffer_d, buffer_h)
    !Get inverse of U = U^-1
    istat = cusolverDnXtrtri_buffersize(cusolver_handle, CUBLAS_FILL_MODE_UPPER, CUBLAS_DIAG_NON_UNIT, n_8, &
                                        cudaDataType(CUDA_C_64F), awork, lda_8, lbuffer_d, lbuffer_h)
    if(cuda_runtime_version <= 12040)  then !prescription for bug of cusolverDnXtrtri_buffersize
      lbuffer_d = lbuffer_d*16; lbuffer_h = lbuffer_h*16
    endif
    allocate(buffer_d(lbuffer_d), buffer_h(lbuffer_h))
    istat = cusolverDnXtrtri(cusolver_handle, CUBLAS_FILL_MODE_UPPER, CUBLAS_DIAG_NON_UNIT, n_8, cudaDataType(CUDA_C_64F), &
                             awork, lda_8, buffer_d, lbuffer_d, buffer_h, lbuffer_h, devinfo)
    deallocate(buffer_d, buffer_h)
    ! Get U^-1 * L^-1 = awork * a
    istat = cublasZtrmm(cublas_handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, n, n, &
                        (1d0, 0d0), awork, lda_in, a, lda_in, a, lda_in) !result in a (in-place)
    ipiv_cpu = ipiv !copy to CPU from GPU
    do i = n, 1, -1
      if(ipiv_cpu(i) > i) istat = cublasZswap(cublas_handle, n, a(lda_in*(i-1)+1), 1, a(lda_in*(ipiv_cpu(i)-1)+1), 1)
    enddo
    deallocate(awork)
  end function zminv_d
  integer function zsv_d(a, b, n, nrhs, lda, ldb) result(istat)
    integer, intent(in) :: n, nrhs
    integer, intent(in), optional :: lda, ldb
    complex(8), device :: a(*), b(*)
    integer(8), device :: ipiv(n)
    integer(8) :: n_8, lda_8, ldb_8, nrhs_8
    integer(1), allocatable, device :: buffer_d(:)
    integer(1), allocatable         :: buffer_h(:)
    integer(8) :: lbuffer_d, lbuffer_h
    integer(4), device :: devinfo
    integer :: lda_in, ldb_in
    istat = cusolver_init()
    lda_in = n; if(present(lda)) lda_in = lda
    ldb_in = n; if(present(ldb)) ldb_in = ldb

    n_8 = int(n,8); lda_8 = int(lda_in,8); ldb_8 = int(ldb_in,8); nrhs_8 = int(nrhs,8)
    istat = cusolverDnXgetrf_buffersize(cusolver_handle, cusolver_params, n_8, n_8, cudaDataType(CUDA_C_64F), a, lda_8, &
                                        cudaDataType(CUDA_C_64F), lbuffer_d, lbuffer_h)
    allocate(buffer_d(lbuffer_d), buffer_h(lbuffer_h))
    istat = cusolverDnXgetrf(cusolver_handle, cusolver_params, n_8, n_8, cudaDataType(CUDA_C_64F), a, lda_8, &
                             ipiv, cudaDataType(CUDA_C_64F), buffer_d, lbuffer_d, buffer_h, lbuffer_h, devinfo)
    deallocate(buffer_d, buffer_h)
    istat = cusolverDnXgetrs(cusolver_handle, cusolver_params, CUBLAS_OP_N, n_8, nrhs_8, cudaDataType(CUDA_C_64F), a, lda_8, &
                             ipiv, cudaDataType(CUDA_C_64F), b, ldb_8, devinfo)
  end function zsv_d
  ! integer function zgev_d(a, n, evl, evl_vec, evr_vec, lda) result(istat) !Not tested
  !   integer, intent(in) :: n
  !   complex(8), device :: a(*)
  !   complex(8), device :: evl(n)
  !   complex(8), device, optional :: evl_vec(*), evr_vec(*)
  !   integer, intent(in), optional :: lda
  !   integer :: lda_in, ldvl, ldvr, lwork
  !   complex(8), allocatable, device :: work(:), vl_loc(:), vr_loc(:)
  !   integer(4), device :: devinfo
  !   integer(4) :: jobvl_mode, jobvr_mode
  !   lda_in = n; if(present(lda)) lda_in = lda
  !   jobvl_mode = CUSOLVER_EIG_MODE_NOVECTOR; ldvl = 1
  !   jobvr_mode = CUSOLVER_EIG_MODE_NOVECTOR; ldvr = 1
  !   if(present(evl_vec)) then; jobvl_mode = CUSOLVER_EIG_MODE_VECTOR; ldvl = lda_in; endif
  !   if(present(evr_vec)) then; jobvr_mode = CUSOLVER_EIG_MODE_VECTOR; ldvr = lda_in; endif
  !   istat = cusolver_init()
  !   allocate(vl_loc(ldvl*n), vr_loc(ldvr*n))
  !   istat = cusolverDnZgeev_bufferSize(cusolver_handle, jobvl_mode, jobvr_mode, n, a, lda_in, &
  !                                       evl, vl_loc, ldvl, vr_loc, ldvr, lwork)
  !   allocate(work(lwork))
  !   istat = cusolverDnZgeev(cusolver_handle, jobvl_mode, jobvr_mode, n, a, lda_in, &
  !                            evl, vl_loc, ldvl, vr_loc, ldvr, work, lwork, devinfo)
  !   if(present(evl_vec)) evl_vec(1:ldvl*n) = vl_loc
  !   if(present(evr_vec)) evr_vec(1:ldvr*n) = vr_loc
  !   deallocate(work, vl_loc, vr_loc)
  ! end function zgev_d
  integer function zhev_d(A, n, evl, il, iu, lda) result(istat) !Not tested
  ! Solving the standard eigenvalue problem Az = lambda z, where A is a Hermitian matrix (GPU version)
  ! Eigenvalues are stored in evl, eigenvectors are stored in A
  !!! range: 1<=il<=iu<=n
    integer, intent(in) :: n
    real(8), intent(out), device :: evl(n)
    complex(8), device :: A(*)
    integer, intent(in), optional :: lda, il, iu
    integer :: lda_in, il_in, iu_in
    real(8) :: vl, vu
    integer, device :: devInfo
    complex(8), allocatable, device :: work(:)
    integer :: m, lwork
    lda_in = n; if(present(lda)) lda_in = lda
    il_in = 1; iu_in = n
    if(present(il)) il_in = il
    if(present(iu)) iu_in = iu
    vl = 0d0; vu = 0d0
    istat = cusolver_init()
    istat = cusolverDnZheevdx_bufferSize(cusolver_handle, CUSOLVER_EIG_MODE_VECTOR, &
                                          CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_UPPER, &
                                          n, a, lda_in, vl, vu, il_in, iu_in, m, evl, lwork)
    allocate(work(lwork))
    istat = cusolverDnZheevdx(cusolver_handle, CUSOLVER_EIG_MODE_VECTOR, &
                               CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_UPPER, &
                               n, a, lda_in, vl, vu, il_in, iu_in, m, evl, work, lwork, devInfo)
    deallocate(work)
  end function zhev_d
  integer function zhgv_d(A, B, n, evl, il, iu, lda, ldb) result(istat)
    integer, intent(in) :: n !size of matrix
    real(8), intent(out), device :: evl(n) !eigenvalues
    complex(8), device :: A(*), B(*)
    integer, intent(in), optional :: lda, ldb, il, iu
    integer :: lda_in, ldb_in, il_in, iu_in
    real(8):: vu = 0d0, vl = 0d0
    integer, device :: devInfo
    complex(8), allocatable, device :: work(:)
    integer :: m, lwork
    lda_in = n; ldb_in = n
    if(present(lda)) lda_in = lda
    if(present(ldb)) ldb_in = ldb
    il_in = 1; iu_in = n
    if(present(il)) il_in = il
    if(present(iu)) iu_in = iu
    istat = cusolver_init()
    istat = cusolverDnZhegvdx_bufferSize(cusolver_handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, &
                                          CUSOLVER_EIG_RANGE_I,  CUBLAS_FILL_MODE_UPPER, &
                                          n, a, lda_in, b, ldb_in, vl, vu, il_in, iu_in, m, evl, lwork)
    allocate(work(lwork))
    istat = cusolverDnZhegvdx(cusolver_handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, &
                                          CUSOLVER_EIG_RANGE_I,  CUBLAS_FILL_MODE_UPPER, &
                                          n, a, lda_in, b, ldb_in, vl, vu, il_in, iu_in, m, evl, work, lwork, devInfo)
    deallocate(work)
  end function zhgv_d
  integer function cusolver_init() result(istat)
    istat = 0
    if(.not.set_cusolver_handle) then 
      istat = cusolverDnCreate(cusolver_handle)
      if(istat /= CUSOLVER_STATUS_SUCCESS) then
        print *, 'Error in cusolverDnCreate'
      endif
      istat = cusolverDnCreateParams(cusolver_params)
      istat = cudaRuntimeGetversion(cuda_runtime_version)
      set_cusolver_handle = .true.
    endif
  end function cusolver_init
  integer function cusolver_finalize() result(istat)
    istat = 0
    if(set_cusolver_handle) then
      istat = cusolverDnDestroy(cusolver_handle)
      if(istat /= CUSOLVER_STATUS_SUCCESS) then
        print *, 'Error in cusolverDnDestroy'
      endif
      istat = cusolverDnDestroyParams(cusolver_params)
      if(istat /= CUSOLVER_STATUS_SUCCESS) then
        print *, 'Error in cusolverDnDestroyParams'
      endif
      set_cusolver_handle = .false.
    endif
  end function cusolver_finalize
#endif
end module m_lapack
