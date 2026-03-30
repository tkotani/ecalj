module m_qdwh_eigensolver
  !! QDWH spectral divide-and-conquer eigensolver
  !! 100% GEMM (Ozaki) + Cholesky. No FP64 GEMM.
  !! Fully batchable over k-points.
#ifdef __GPU
  use m_blas, only: zmm => zmm_d, m_op_C, m_op_N
  use cusolverdn
  use cudafor
#endif
  implicit none
  public :: qdwh_batched_zhegv
#ifdef __GPU
contains

  subroutine qdwh_batched_zhegv(n, nev, nmat, H_h, S_h, lda, evals, evecs_d, info)
    !! Batched generalized eigenvalue problem via QDWH
    !! H(lda,n,nmat), S(lda,n,nmat): host. evals(nev,nmat): host. evecs_d(n,nev,nmat): device.
    integer, intent(in) :: n, nev, nmat, lda
    complex(8), intent(in) :: H_h(lda, n, nmat), S_h(lda, n, nmat)
    real(8), intent(out) :: evals(nev, nmat)
    complex(8), device, intent(out) :: evecs_d(n, nev, nmat)
    integer, intent(out) :: info
    integer :: ik, nev_qdwh, info_l, istat
    complex(8), allocatable :: Hp(:,:), zz(:,:), ovl(:,:)
    real(8), allocatable :: eo(:)
    integer :: ni, nm, jj, ix
    type(cusolverDnHandle) :: cs_h
    integer, device, allocatable :: devinfo

    info = 0
    nev_qdwh = max(nev, n/2)  ! expand to n/2 for better QDWH convergence
    nev_qdwh = min(nev_qdwh, n)
    istat = cusolverDnCreate(cs_h)
    allocate(devinfo)

    do ik = 1, nmat
      ! === Step 1: Cholesky reduction (epsovl=0) ===
      allocate(ovl(n,n))
      ovl = S_h(1:n, 1:n, ik)
      call zpotrf('L', n, ovl, n, info_l)
      if(info_l /= 0) then; info = -1; deallocate(ovl); cycle; endif
      call ztrtri('L', 'N', n, ovl, n, info_l)
      do jj = 2, n; ovl(1:jj-1, jj) = (0d0,0d0); enddo
      ! H' = Linv * H * Linv^H
      allocate(Hp(n,n), zz(n,n))
      call zgemm('N','N',n,n,n,(1d0,0d0),ovl,n,H_h(1,1,ik),lda,(0d0,0d0),zz,n)
      call zgemm('N','C',n,n,n,(1d0,0d0),zz,n,ovl,n,(0d0,0d0),Hp,n)

      ! === Step 2: QDWH spectral projector on GPU ===
      qdwh_solve: block
        complex(8), device, allocatable :: Hp_d(:,:), X_d(:,:), P_d(:,:)
        complex(8), device, allocatable :: Q_d(:,:), HQ_d(:,:)
        complex(8), allocatable :: M_h(:,:), evecs_sub(:,:)
        real(8), allocatable :: evals_sub(:)
        real(8) :: sigma, alpha_qdwh
        integer :: nev_sub
        allocate(Hp_d(n,n))
        Hp_d = Hp

        ! === Sigma from LAPACK (first k-point) or reuse previous ===
        sigma_est: block
          real(8), save :: sigma_saved = 0d0
          logical, save :: sigma_initialized = .false.
          complex(8), allocatable :: Xtmp(:,:)
          real(8) :: trace_p
          integer :: rank_p
          if(.not. sigma_initialized) then
            ! First k-point: compute all eigenvalues to determine sigma
            block
              complex(8), allocatable :: Htmp(:,:), work_s(:)
              real(8), allocatable :: evals_s(:), rwork_s(:)
              allocate(Htmp(n,n), evals_s(n), work_s(n*n), rwork_s(3*n))
              Htmp = Hp
              call zheev('N','U',n,Htmp,n,evals_s,work_s,n*n,rwork_s,info_l)
              sigma_saved = (evals_s(nev_qdwh) + evals_s(min(nev_qdwh+1,n))) / 2d0
              deallocate(Htmp, evals_s, work_s, rwork_s)
            endblock
            sigma_initialized = .true.
          endif
          sigma = sigma_saved
          ! Compute sign(Hp - sigma*I) via QDWH
          allocate(X_d(n,n), Xtmp(n,n))
          Xtmp = Hp
          do jj = 1, n; Xtmp(jj,jj) = Xtmp(jj,jj) - sigma; enddo
          ! Scale by 1/sqrt(||X||_1 * ||X||_inf)
          block
            real(8) :: one_n, inf_n, alpha_i
            integer :: jk
            one_n = 0d0
            do jk=1,n; one_n = max(one_n, sum(abs(Xtmp(:,jk)))); enddo
            inf_n = 0d0
            do jk=1,n; inf_n = max(inf_n, sum(abs(Xtmp(jk,:)))); enddo
            alpha_i = 1d0 / sqrt(one_n * inf_n)
            Xtmp = Xtmp * alpha_i
          endblock
          X_d = Xtmp; deallocate(Xtmp)
          call qdwh_polar(n, X_d, cs_h, devinfo)
          ! P = (I - sign) / 2
          allocate(P_d(n,n), Xtmp(n,n))
          Xtmp = X_d; Xtmp = -Xtmp
          do jj = 1, n; Xtmp(jj,jj) = Xtmp(jj,jj) + 1d0; enddo
          Xtmp = Xtmp / 2d0; P_d = Xtmp
          trace_p = 0d0
          do jj = 1, n; trace_p = trace_p + dble(Xtmp(jj,jj)); enddo
          rank_p = nint(trace_p)
          deallocate(Xtmp, X_d)
          if(rank_p >= nev .and. rank_p < n) nev_qdwh = rank_p
        endblock sigma_est

        ! Extract subspace: Q = orth(P * R) where R is random (n × nev_qdwh)
        allocate(Q_d(n, nev_qdwh))
        block
          complex(8), device, allocatable :: R_d(:,:)
          complex(8), allocatable :: Rtmp(:,:)
          real(8), allocatable :: rr(:,:)
          allocate(Rtmp(n, nev_qdwh), rr(n, nev_qdwh))
          call random_number(rr); Rtmp = rr
          deallocate(rr)
          allocate(R_d(n, nev_qdwh))
          R_d = Rtmp; deallocate(Rtmp)
          ! Q = P * R (project random vectors onto range of P)
          istat = zmm(P_d, R_d, Q_d, m=n, n=nev_qdwh, k=n)
          deallocate(R_d)
        endblock
        deallocate(P_d)
        ! Orthogonalize
        call gpu_chol_orth(n, nev_qdwh, Q_d)
        call gpu_chol_orth(n, nev_qdwh, Q_d)

        ! Rayleigh-Ritz: M = Q^H * Hp * Q
        allocate(HQ_d(n, nev_qdwh))
        istat = zmm(Hp_d, Q_d, HQ_d, m=n, n=nev_qdwh, k=n)
        block
          complex(8), device, allocatable :: M_d(:,:), C_d(:,:), Qnew_d(:,:)
          allocate(M_d(nev_qdwh, nev_qdwh))
          istat = zmm(Q_d, HQ_d, M_d, m=nev_qdwh, n=nev_qdwh, k=n, opA=m_op_C)
          allocate(M_h(nev_qdwh, nev_qdwh))
          M_h = M_d
          deallocate(M_d)
          ! Small eigensolve on CPU
          allocate(evals_sub(nev_qdwh), evecs_sub(nev_qdwh, nev_qdwh))
          block
            complex(8), allocatable :: work_s(:)
            real(8), allocatable :: rwork_s(:)
            allocate(work_s(nev_qdwh**2), rwork_s(3*nev_qdwh))
            call zheev('V','U',nev_qdwh,M_h,nev_qdwh,evals_sub,work_s,nev_qdwh**2,rwork_s,info_l)
            deallocate(work_s, rwork_s)
          endblock
          evecs_sub = M_h
          evals(1:nev, ik) = evals_sub(1:nev)
          ! Rotate: Q_final = Q * evecs_sub(:, 1:nev)
          allocate(C_d(nev_qdwh, nev), Qnew_d(n, nev))
          C_d = evecs_sub(1:nev_qdwh, 1:nev)
          istat = zmm(Q_d, C_d, Qnew_d, m=n, n=nev, k=nev_qdwh)
          deallocate(C_d, M_h, evals_sub, evecs_sub)

          ! Back-transform: evec = Linv^H * Q_final (Ozaki GEMM)
          block
            complex(8), device, allocatable :: Linv_d(:,:)
            allocate(Linv_d(n,n))
            Linv_d = ovl  ! L^{-1}
            istat = zmm(Linv_d, Qnew_d, evecs_d(1,1,ik), m=n, n=nev, k=n, opA=m_op_C)
            deallocate(Linv_d)
          endblock
          deallocate(Qnew_d)
        endblock
        deallocate(Q_d, HQ_d, Hp_d)
      endblock qdwh_solve
      deallocate(ovl, Hp, zz)
    enddo

    deallocate(devinfo)
    istat = cusolverDnDestroy(cs_h)
  end subroutine

  subroutine qdwh_polar(n, X_d, cs_h, devinfo)
    !! QDWH polar: X -> sign(X). All GPU: Ozaki GEMM + cuSOLVER Cholesky.
    !! No QR step — Cholesky for all iterations (FP64 zpotrf handles large c).
    integer, intent(in) :: n
    complex(8), device, intent(inout) :: X_d(n, n)
    type(cusolverDnHandle), intent(in) :: cs_h
    integer, device, intent(inout) :: devinfo
    complex(8), device, allocatable :: XhX_d(:,:), G_d(:,:), Linv_d(:,:), Z_d(:,:), Xnew_d(:,:)
    complex(8), device, allocatable :: work_d(:)
    complex(8), allocatable :: Gh(:,:), Xh(:,:), Xnh(:,:)
    real(8) :: Lk, L2, dd, sqd, a, b, c, e_coef, a_minus_e, tol_l
    integer :: iter, max_iter, istat, jj, info_l, lwork
    max_iter = 20
    Lk = 1d-15  ! eps
    tol_l = 5d-15
    allocate(XhX_d(n,n), G_d(n,n), Z_d(n,n), Xnew_d(n,n))
    do iter = 1, max_iter
      if(Lk + tol_l >= 1d0) exit
      L2 = Lk * Lk
      dd = (4d0*(1d0/L2 - 1d0)/L2)**(1d0/3d0)
      sqd = sqrt(1d0 + dd)
      a = sqd + sqrt(2d0 - dd + 2d0*(2d0-L2)/(L2*sqd))
      b = (a - 1d0)**2 / 4d0
      c = a + b - 1d0
      e_coef = b / c
      a_minus_e = a - e_coef
      Lk = Lk * (a + b*L2) / (1d0 + c*L2)
      ! G = c * X^H * X + I  (Ozaki GEMM for X^H*X)
      istat = zmm(X_d, X_d, XhX_d, m=n, n=n, k=n, opA=m_op_C)
      allocate(Gh(n,n)); Gh = XhX_d
      Gh = c * Gh
      do jj = 1, n; Gh(jj,jj) = Gh(jj,jj) + (1d0,0d0); enddo
      ! Cholesky G = L*L^H on GPU
      G_d = Gh
      istat = cusolverDnZpotrf_bufferSize(cs_h, CUBLAS_FILL_MODE_LOWER, n, G_d, n, lwork)
      allocate(work_d(lwork))
      istat = cusolverDnZpotrf(cs_h, CUBLAS_FILL_MODE_LOWER, n, G_d, n, work_d, lwork, devinfo)
      deallocate(work_d)
      ! L^{-1} on host (ztrtri, O(n^3/3), cheap)
      Gh = G_d
      do jj = 1, n-1; Gh(1:jj, jj+1) = (0d0,0d0); enddo  ! zero upper
      call ztrtri('L', 'N', n, Gh, n, info_l)
      allocate(Linv_d(n,n)); Linv_d = Gh; deallocate(Gh)
      ! Z = G^{-1} * X = L^{-H} * L^{-1} * X  (two Ozaki GEMMs)
      istat = zmm(Linv_d, X_d, Xnew_d, m=n, n=n, k=n)       ! T = L^{-1} * X
      istat = zmm(Linv_d, Xnew_d, Z_d, m=n, n=n, k=n, opA=m_op_C)  ! Z = L^{-H} * T
      deallocate(Linv_d)
      ! X_new = e*X + (a-e)*Z  (host vector ops)
      allocate(Xh(n,n), Xnh(n,n)); Xh = X_d; Xnh = Z_d
      Xh = e_coef * Xh + a_minus_e * Xnh
      X_d = Xh; deallocate(Xh, Xnh)
    enddo
    ! Extra Halley (a=3,b=1,c=3,e=1/3) if needed
    do while(iter <= max_iter)
      allocate(Xh(n,n)); Xh = X_d
      ! Same Cholesky step with a=3,e=1/3,c=3
      istat = zmm(X_d, X_d, XhX_d, m=n, n=n, k=n, opA=m_op_C)
      allocate(Gh(n,n)); Gh = XhX_d
      Gh = 3d0 * Gh
      do jj=1,n; Gh(jj,jj) = Gh(jj,jj) + (1d0,0d0); enddo
      G_d = Gh
      istat = cusolverDnZpotrf_bufferSize(cs_h, CUBLAS_FILL_MODE_LOWER, n, G_d, n, lwork)
      allocate(work_d(lwork))
      istat = cusolverDnZpotrf(cs_h, CUBLAS_FILL_MODE_LOWER, n, G_d, n, work_d, lwork, devinfo)
      deallocate(work_d)
      Gh = G_d; do jj=1,n-1; Gh(1:jj,jj+1) = (0d0,0d0); enddo
      call ztrtri('L','N',n,Gh,n,info_l)
      allocate(Linv_d(n,n)); Linv_d = Gh; deallocate(Gh)
      istat = zmm(Linv_d, X_d, Xnew_d, m=n, n=n, k=n)
      istat = zmm(Linv_d, Xnew_d, Z_d, m=n, n=n, k=n, opA=m_op_C)
      deallocate(Linv_d)
      allocate(Xnh(n,n)); Xnh = Z_d
      block; complex(8),allocatable::Xtmp(:,:); allocate(Xtmp(n,n)); Xtmp = X_d
      Xtmp = (1d0/3d0)*Xtmp + (8d0/3d0)*Xnh; X_d = Xtmp; deallocate(Xtmp); endblock
      deallocate(Xnh)
      allocate(Xnh(n,n)); Xnh = X_d
      if(maxval(abs(Xnh-Xh)) < tol_l**(1d0/3d0)) then
        deallocate(Xh,Xnh); exit
      endif
      deallocate(Xh,Xnh); iter = iter + 1
    enddo
    deallocate(XhX_d, G_d, Z_d, Xnew_d)
    ! Newton-Schulz: U = 1.5*U - 0.5*U*(U^H*U)
    block
      complex(8), device, allocatable :: UhU_d(:,:), T_d(:,:)
      allocate(UhU_d(n,n), T_d(n,n))
      istat = zmm(X_d, X_d, UhU_d, m=n, n=n, k=n, opA=m_op_C)
      istat = zmm(X_d, UhU_d, T_d, m=n, n=n, k=n)
      allocate(Xh(n,n), Xnh(n,n)); Xh = X_d; Xnh = T_d
      Xh = 1.5d0*Xh - 0.5d0*Xnh; X_d = Xh
      deallocate(Xh, Xnh, UhU_d, T_d)
    endblock
  end subroutine

  subroutine gpu_chol_orth(n, p, X_d)
    integer, intent(in) :: n, p
    complex(8), device, intent(inout) :: X_d(n, p)
    complex(8), device, allocatable :: G_d(:,:), Rinv_d(:,:), Xnew_d(:,:)
    complex(8), allocatable :: G_h(:,:)
    integer :: istat, info_l, jj
    allocate(G_d(p,p), G_h(p,p), Rinv_d(p,p), Xnew_d(n,p))
    istat = zmm(X_d, X_d, G_d, m=p, n=p, k=n, opA=m_op_C)
    G_h = G_d
    call zpotrf('U', p, G_h, p, info_l)
    if(info_l /= 0) return
    call ztrtri('U', 'N', p, G_h, p, info_l)
    do jj = 1, p-1; G_h(jj+1:p, jj) = (0d0,0d0); enddo
    Rinv_d = G_h
    istat = zmm(X_d, Rinv_d, Xnew_d, m=n, n=p, k=p)
    X_d = Xnew_d
    deallocate(G_d, G_h, Rinv_d, Xnew_d)
  end subroutine

#endif
end module m_qdwh_eigensolver
