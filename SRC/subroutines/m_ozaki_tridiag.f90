module m_ozaki_tridiag
  !! Eigensolver: CPU tridiag + Ozaki GEMM back-transform
  !! Supports both epsovl>0 (eigensolve reduction) and epsovl=0 (Cholesky reduction)
#ifdef __GPU
  use m_blas, only: zmm => zmm_d, m_op_C
  use cudafor
#endif
  implicit none
  public :: ozaki_tridiag_zhegv
#ifdef __GPU
contains
  subroutine ozaki_tridiag_zhegv(n, nev, H_h, lda, S_h, evals, evecs_d, epsovl_in, info)
    integer, intent(in) :: n, nev, lda
    real(8), intent(in) :: epsovl_in
    complex(8), intent(in) :: H_h(lda, n), S_h(lda, n)
    real(8), intent(out) :: evals(nev)
    complex(8), device, intent(out) :: evecs_d(n, nev)
    integer, intent(out) :: info
    complex(8), allocatable :: Hp(:,:), zz(:,:), evecs_hp(:,:)
    real(8), allocatable :: eo(:), rwork(:), evals_hp(:)
    complex(8), allocatable :: work(:)
    integer, allocatable :: iwork(:), ifail(:)
    integer :: info_l, istat_g, nm, nev_hp, nevout
    real(8) :: vl, vu, abstol
    info = 0
    ! === Step 1: Reduce generalized → standard eigenvalue problem ===
    call reduce_generalized(n, H_h, lda, S_h, epsovl_in, nm, Hp, zz, info_l)
    if(info_l /= 0) then; info = info_l; return; endif
    ! Hp is nm×nm standard Hermitian, zz is n×nm transformation
    ! === Step 2: Solve standard eigenvalue problem on CPU ===
    nev_hp = min(nev, nm)
    allocate(evals_hp(nev_hp), evecs_hp(nm, nev_hp))
    allocate(work(nm*nm), rwork(7*nm), iwork(5*nm), ifail(nm))
    abstol = 1d-12
    call zheevx('V','I','U', nm, Hp, nm, vl, vu, 1, nev_hp, abstol, nevout, &
         evals_hp, evecs_hp, nm, work, nm*nm, rwork, iwork, ifail, info_l)
    evals(1:nev_hp) = evals_hp(1:nev_hp)
    if(nev_hp < nev) evals(nev_hp+1:nev) = 1d99
    deallocate(work, rwork, iwork, ifail, evals_hp, Hp)
    ! === Step 3: Back-transform on GPU via Ozaki GEMM: evec = zz * evecs_hp ===
    block
      complex(8), device, allocatable :: zz_d(:,:), Y_d(:,:)
      allocate(zz_d(n, nm), Y_d(nm, nev_hp))
      zz_d = zz
      Y_d = evecs_hp(1:nm, 1:nev_hp)
      istat_g = zmm(zz_d, Y_d, evecs_d, m=n, n=nev_hp, k=nm)
      deallocate(zz_d, Y_d)
    endblock
    deallocate(zz, evecs_hp)
  end subroutine

  subroutine reduce_generalized(n, H_h, lda, S_h, epsovl, nm, Hp, zz, info)
    !! Reduce H*x=λ*S*x to standard H'*y=λ*y
    !! epsovl=0: Cholesky, epsovl>0: overlap eigensolve with cutoff
    !! Output: Hp(nm,nm), zz(n,nm) such that H' = zz^H*H*zz, x = zz*y
    integer, intent(in) :: n, lda
    complex(8), intent(in) :: H_h(lda, n), S_h(lda, n)
    real(8), intent(in) :: epsovl
    integer, intent(out) :: nm, info
    complex(8), allocatable, intent(out) :: Hp(:,:), zz(:,:)
    complex(8), allocatable :: tmp(:,:), ovl(:,:)
    real(8), allocatable :: eo(:)
    integer :: ix, ni, jj, info_l
    info = 0
    if(epsovl < 1d-14) then
      ! === Cholesky path: S = L*L^H, zz = L^{-H} ===
      allocate(ovl(n,n))
      ovl = S_h(1:n, 1:n)
      call zpotrf('L', n, ovl, n, info_l)
      if(info_l /= 0) then; info = -1; return; endif
      call ztrtri('L', 'N', n, ovl, n, info_l)
      do jj = 2, n; ovl(1:jj-1, jj) = (0d0,0d0); enddo
      ! zz = (L^{-1})^H = L^{-H} (upper triangular stored as full)
      nm = n
      allocate(zz(n, n))
      zz = conjg(transpose(ovl))  ! L^{-H}
      deallocate(ovl)
    else
      ! === Eigensolve path: S = U*D*U^H, zz = U*D^{-1/2} with cutoff ===
      allocate(ovl(n,n), eo(n))
      ovl = S_h(1:n, 1:n)
      block
        complex(8), allocatable :: work_s(:)
        real(8), allocatable :: rwork_s(:)
        allocate(work_s(n*n), rwork_s(3*n))
        call zheev('V', 'U', n, ovl, n, eo, work_s, n*n, rwork_s, info_l)
        deallocate(work_s, rwork_s)
      endblock
      ni = 1
      do ix = 1, n
        if(eo(ix) > epsovl) then; ni = ix; exit; endif
      enddo
      nm = n - ni + 1
      allocate(zz(n, nm))
      zz = ovl(:, ni:n)
      do ix = ni, n
        zz(:, ix-ni+1) = zz(:, ix-ni+1) / sqrt(eo(ix))
      enddo
      deallocate(ovl, eo)
    endif
    ! H' = zz^H * H * zz
    allocate(tmp(nm, n), Hp(nm, nm))
    call zgemm('C','N', nm, n, n, (1d0,0d0), zz, n, H_h, lda, (0d0,0d0), tmp, nm)
    call zgemm('N','N', nm, nm, n, (1d0,0d0), tmp, nm, zz, n, (0d0,0d0), Hp, nm)
    deallocate(tmp)
  end subroutine
#endif
end module m_ozaki_tridiag
