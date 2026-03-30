module m_chefsi
  !! Chebyshev-filtered subspace iteration eigensolver
  !! All heavy GEMM via Ozaki (zmm_d). FP64 only for small Rayleigh-Ritz and Cholesky.
#ifdef __GPU
  use m_blas, only: zmm => zmm_d, m_op_C
  use cudafor
#endif
  implicit none
  public :: chefsi_zhegv
contains

#ifdef __GPU
  subroutine chefsi_zhegv(n, nev, H_h, lda, S_h, evals, evecs_d, max_iter, deg, tol, info)
    !! Solve generalized Hermitian eigenvalue problem H*x = λ*S*x
    !! via Cholesky reduction + CheFSI. Returns lowest nev eigenvalues.
    !! H_h, S_h: host arrays (lda × n). evecs_d: device array (n × nev).
    integer, intent(in) :: n, nev, lda, max_iter, deg
    complex(8), intent(in) :: H_h(lda, n), S_h(lda, n)
    real(8), intent(out) :: evals(nev)
    complex(8), device, intent(out) :: evecs_d(n, nev)
    real(8), intent(in) :: tol
    integer, intent(out) :: info
    ! Local
    complex(8), allocatable :: Linv(:,:), Hp(:,:), G(:,:), M(:,:), work(:), Xtmp(:,:)
    real(8), allocatable :: evals_all(:), rwork(:), rr(:,:)
    complex(8), device, allocatable :: Hp_d(:,:), X_d(:,:), HX_d(:,:), Y_d(:,:)
    complex(8), device, allocatable :: Yk_d(:,:), Ykm1_d(:,:), Ykm2_d(:,:), HY_d(:,:)
    complex(8), device, allocatable :: G_d(:,:), Rinv_d(:,:), Xnew_d(:,:)
    real(8) :: lam_upper, lam_max, e_cf, c_cf, rnorm
    integer :: nev_w, nbuf, iter, k, jj, info_l, istat
    logical :: converged

    info = 0
    nbuf = max(20, nev/3)
    nev_w = min(nev + nbuf, n)

    ! === Cholesky reduction: H' = L^{-1} * H * L^{-H} ===
    allocate(Linv(n,n), Hp(n,n))
    Linv(1:n,1:n) = S_h(1:n,1:n)
    call zpotrf('L', n, Linv, n, info_l)
    if(info_l /= 0) then; info = -1; return; endif
    call ztrtri('L', 'N', n, Linv, n, info_l)
    do jj = 2, n; Linv(1:jj-1, jj) = (0d0,0d0); enddo
    call zgemm('N','N',n,n,n,(1d0,0d0),Linv,n,H_h,lda,(0d0,0d0),Hp,n)
    allocate(Xtmp(n,n))
    call zgemm('N','C',n,n,n,(1d0,0d0),Hp,n,Linv,n,(0d0,0d0),Xtmp,n)
    Hp = Xtmp

    ! === Spectral bounds ===
    allocate(evals_all(n), work(n*n), rwork(3*n))
    Xtmp = Hp
    call zheev('N','U',n,Xtmp,n,evals_all,work,n*n,rwork,info_l)
    lam_upper = evals_all(nev_w) + 0.1d0*(evals_all(nev_w)-evals_all(1))
    lam_max = evals_all(n)
    deallocate(evals_all, work, rwork, Xtmp)

    ! === Copy H' to GPU ===
    allocate(Hp_d(n,n)); Hp_d = Hp

    ! === Random initial vectors ===
    allocate(X_d(n,nev_w), HX_d(n,nev_w))
    allocate(Xtmp(n,nev_w), rr(n,nev_w))
    call random_number(rr); Xtmp = rr
    X_d = Xtmp; deallocate(Xtmp, rr)

    ! === Allocate work arrays ===
    allocate(G_d(nev_w,nev_w), G(nev_w,nev_w), Rinv_d(nev_w,nev_w), Xnew_d(n,nev_w))
    allocate(Y_d(n,nev_w), Yk_d(n,nev_w), Ykm1_d(n,nev_w), Ykm2_d(n,nev_w), HY_d(n,nev_w))
    allocate(evals_all(nev_w))

    ! === Orthogonalize + initial Rayleigh-Ritz ===
    call ozaki_chol_orth(n, nev_w, X_d, G_d, G, Rinv_d, Xnew_d)
    istat = zmm(Hp_d, X_d, HX_d, m=n, n=nev_w, k=n)
    call ozaki_rayleigh_ritz(n, nev_w, X_d, HX_d, evals_all, G_d, G, Rinv_d, Xnew_d)

    ! === Main CheFSI loop ===
    do iter = 1, max_iter
      e_cf = (lam_max - lam_upper) / 2d0
      c_cf = (lam_max + lam_upper) / 2d0
      ! Chebyshev filter: Y = T_deg(scaled H) * X
      call ozaki_chebyshev_filter(n, nev_w, Hp_d, X_d, Y_d, Yk_d, Ykm1_d, Ykm2_d, HY_d, deg, c_cf, e_cf)
      ! Orthogonalize (twice)
      call ozaki_chol_orth(n, nev_w, Y_d, G_d, G, Rinv_d, Xnew_d)
      call ozaki_chol_orth(n, nev_w, Y_d, G_d, G, Rinv_d, Xnew_d)
      X_d = Y_d
      ! Rayleigh-Ritz
      istat = zmm(Hp_d, X_d, HX_d, m=n, n=nev_w, k=n)
      call ozaki_rayleigh_ritz(n, nev_w, X_d, HX_d, evals_all, G_d, G, Rinv_d, Xnew_d)
      ! Update bounds
      lam_upper = evals_all(nev_w) + 0.05d0*abs(evals_all(nev_w)-evals_all(1))
      ! Check convergence (first nev only)
      rnorm = ozaki_residual(n, nev_w, nev, HX_d, X_d, evals_all)
      if(rnorm < tol) exit
    enddo
    evals = evals_all(1:nev)

    ! === Back-transform eigenvectors: x = L^{-H} * y ===
    allocate(Xtmp(n,nev))
    Xtmp = X_d(1:n, 1:nev)  ! D2H
    block
      complex(8), device, allocatable :: Linv_d(:,:), evec_tmp_d(:,:)
      allocate(Linv_d(n,n), evec_tmp_d(n,nev))
      Linv_d = Linv
      ! evecs = (L^{-1})^H * y = Linv^H * y
      istat = zmm(Linv_d, X_d, evec_tmp_d, m=n, n=nev, k=n, opA=m_op_C)
      evecs_d = evec_tmp_d
      deallocate(Linv_d, evec_tmp_d)
    endblock

    deallocate(Linv, Hp, Hp_d, X_d, HX_d, Y_d, Yk_d, Ykm1_d, Ykm2_d, HY_d)
    deallocate(G_d, G, Rinv_d, Xnew_d, evals_all, Xtmp)
  end subroutine

  subroutine ozaki_chol_orth(n, p, X_d, G_d, G_h, Rinv_d, Xnew_d)
    integer, intent(in) :: n, p
    complex(8), device, intent(inout) :: X_d(n,p)
    complex(8), device :: G_d(p,p), Rinv_d(p,p), Xnew_d(n,p)
    complex(8) :: G_h(p,p)
    integer :: istat, info_l, jj
    istat = zmm(X_d, X_d, G_d, m=p, n=p, k=n, opA=m_op_C)
    G_h = G_d
    call zpotrf('U', p, G_h, p, info_l)
    if(info_l /= 0) return
    call ztrtri('U', 'N', p, G_h, p, info_l)
    do jj = 1, p-1; G_h(jj+1:p, jj) = (0d0,0d0); enddo
    Rinv_d = G_h
    istat = zmm(X_d, Rinv_d, Xnew_d, m=n, n=p, k=p)
    X_d = Xnew_d
  end subroutine

  subroutine ozaki_rayleigh_ritz(n, p, X_d, HX_d, evals, G_d, G_h, C_d, Xnew_d)
    integer, intent(in) :: n, p
    complex(8), device, intent(inout) :: X_d(n,p), HX_d(n,p)
    real(8), intent(out) :: evals(p)
    complex(8), device :: G_d(p,p), C_d(p,p), Xnew_d(n,p)
    complex(8) :: G_h(p,p)
    complex(8), allocatable :: work(:)
    real(8), allocatable :: rwork(:)
    integer :: istat, info_l
    istat = zmm(X_d, HX_d, G_d, m=p, n=p, k=n, opA=m_op_C)
    G_h = G_d
    allocate(work(p*p), rwork(3*p))
    call zheev('V', 'U', p, G_h, p, evals, work, p*p, rwork, info_l)
    deallocate(work, rwork)
    C_d = G_h
    istat = zmm(X_d, C_d, Xnew_d, m=n, n=p, k=p)
    X_d = Xnew_d
    istat = zmm(HX_d, C_d, Xnew_d, m=n, n=p, k=p)
    HX_d = Xnew_d
  end subroutine

  subroutine ozaki_chebyshev_filter(n, p, H_d, X_d, Y_d, Yk_d, Ykm1_d, Ykm2_d, HY_d, deg, c, e)
    integer, intent(in) :: n, p, deg
    complex(8), device, intent(in) :: H_d(n,n), X_d(n,p)
    complex(8), device, intent(out) :: Y_d(n,p)
    complex(8), device :: Yk_d(n,p), Ykm1_d(n,p), Ykm2_d(n,p), HY_d(n,p)
    real(8), intent(in) :: c, e
    complex(8), allocatable :: tmp1(:,:), tmp2(:,:), tmp3(:,:)
    integer :: k, istat
    ! Y0 = X
    Ykm1_d = X_d
    ! Y1 = (H*X - c*X) / e
    istat = zmm(H_d, Ykm1_d, HY_d, m=n, n=p, k=n)
    allocate(tmp1(n,p), tmp2(n,p))
    tmp1 = HY_d; tmp2 = Ykm1_d
    tmp1 = (tmp1 - c*tmp2) / e
    Yk_d = tmp1
    ! Three-term recurrence
    do k = 2, deg
      Ykm2_d = Ykm1_d
      Ykm1_d = Yk_d
      istat = zmm(H_d, Ykm1_d, HY_d, m=n, n=p, k=n)
      tmp1 = HY_d; tmp2 = Ykm1_d; tmp3 = Ykm2_d
      tmp1 = 2d0*(tmp1 - c*tmp2)/e - tmp3
      Yk_d = tmp1
    enddo
    Y_d = Yk_d
    deallocate(tmp1, tmp2)
    if(allocated(tmp3)) deallocate(tmp3)
  end subroutine

  real(8) function ozaki_residual(n, nev_w, nev, HX_d, X_d, lambda)
    integer, intent(in) :: n, nev_w, nev
    complex(8), device, intent(in) :: HX_d(n,nev_w), X_d(n,nev_w)
    real(8), intent(in) :: lambda(nev_w)
    complex(8), allocatable :: HX_h(:,:), X_h(:,:)
    real(8) :: rmax, ri
    integer :: i
    allocate(HX_h(n,nev_w), X_h(n,nev_w))
    HX_h = HX_d; X_h = X_d
    rmax = 0d0
    do i = 1, nev
      ri = sqrt(sum(abs(HX_h(:,i) - lambda(i)*X_h(:,i))**2)) / max(abs(lambda(i)), 1d-30)
      rmax = max(rmax, ri)
    enddo
    deallocate(HX_h, X_h)
    ozaki_residual = rmax
  end function

#endif
end module m_chefsi
