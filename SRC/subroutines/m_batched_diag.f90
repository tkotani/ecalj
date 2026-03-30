module m_batched_diag
  !! All-GPU batched eigenvalue solver: Cholesky + cuSOLVER streams + Ozaki GEMM
  !! Solves H*x = λ*S*x for nmat independent k-points with epsovl=0 (Cholesky).
#ifdef __GPU
  use m_blas, only: zmm => zmm_d, m_op_C
  use cusolverdn
  use cudafor
#endif
  implicit none
  public :: batched_diag_gpu
#ifdef __GPU
contains

  subroutine batched_diag_gpu(nd, nmat, nbandmx, nev_list, &
       hamm, ovlm, evals_out, evecs_out, numprocs, ngpu_ranks, info)
    !! Batched generalized eigenvalue: H*x = λ*S*x via Cholesky + Zheevd on GPU streams
    !!
    !! Input:
    !!   nd          : matrix dimension (padded, same for all k-points)
    !!   nmat        : number of k-points (matrices)
    !!   nbandmx     : leading dimension of hamm/ovlm
    !!   nev_list(nmat): number of eigenvalues requested per k-point
    !!   hamm(nbandmx, nbandmx, nmat): Hamiltonians (host)
    !!   ovlm(nbandmx, nbandmx, nmat): Overlap matrices (host)
    !!   numprocs    : total MPI ranks (for NSTR computation)
    !!   ngpu_ranks  : number of GPU ranks
    !!
    !! Output:
    !!   evals_out(nbandmx, nmat): eigenvalues (host), padded with 1d99
    !!   evecs_out(nd, maxval(nev_list), nmat): eigenvectors (host)
    !!   info        : 0 on success
    integer, intent(in) :: nd, nmat, nbandmx, numprocs, ngpu_ranks
    integer, intent(in) :: nev_list(nmat)
    complex(8), intent(in) :: hamm(nbandmx, nbandmx, nmat)
    complex(8), intent(in) :: ovlm(nbandmx, nbandmx, nmat)
    real(8), intent(out) :: evals_out(nbandmx, nmat)
    complex(8), intent(out) :: evecs_out(nd, maxval(nev_list), nmat)
    integer, intent(out) :: info

    integer :: NSTR, is, jd, jj, nev, lw_pf, lw_ev, info_l, istat
    integer(kind=8), allocatable :: strms(:)
    type(cusolverDnHandle), allocatable :: cs_str(:)
    integer, device, allocatable :: devinfo
    complex(8), device, allocatable :: omat_all(:,:,:), Linv_all(:,:,:)
    real(8), device, allocatable :: eo_all(:,:)
    complex(8), device, allocatable :: S_d(:,:), work_d(:)
    complex(8), allocatable :: Lh(:,:)
    real(8), allocatable :: eo_h(:,:)

    info = 0
    NSTR = min(max(1, numprocs / max(1, ngpu_ranks)), nmat)

    ! Create streams and cuSOLVER handles
    allocate(strms(NSTR), cs_str(NSTR), devinfo)
    do is = 1, NSTR
      istat = cudaStreamCreate(strms(is))
      istat = cusolverDnCreate(cs_str(is))
      istat = cusolverDnSetStream(cs_str(is), strms(is))
    enddo

    allocate(omat_all(nd, nd, nmat), Linv_all(nd, nd, nmat), eo_all(nd, nmat))

    ! === Phase A1: Cholesky zpotrf on GPU streams ===
    allocate(S_d(nd, nd))
    istat = cusolverDnZpotrf_bufferSize(cs_str(1), CUBLAS_FILL_MODE_LOWER, nd, S_d, nd, lw_pf)
    deallocate(S_d)
    allocate(work_d(lw_pf * NSTR))
    do jd = 1, nmat
      is = mod(jd-1, NSTR) + 1
      allocate(S_d(nd, nd)); S_d = ovlm(1:nd, 1:nd, jd)
      omat_all(:,:,jd) = S_d; deallocate(S_d)
      istat = cusolverDnZpotrf(cs_str(is), CUBLAS_FILL_MODE_LOWER, nd, &
           omat_all(1,1,jd), nd, work_d((is-1)*lw_pf+1), lw_pf, devinfo)
    enddo
    do is = 1, NSTR; istat = cudaStreamSynchronize(strms(is)); enddo
    deallocate(work_d)

    ! === Phase A2: L^{-1} via cuBLAS ztrsm on streams ===
    block
      use cublas_v2
      type(cublasHandle), allocatable :: cb_str(:)
      complex(8), allocatable :: Ih(:,:)
      allocate(cb_str(NSTR))
      do is = 1, NSTR
        istat = cublasCreate(cb_str(is))
        istat = cublasSetStream(cb_str(is), strms(is))
      enddo
      allocate(Ih(nd,nd)); Ih = (0d0,0d0)
      do jj = 1, nd; Ih(jj,jj) = (1d0,0d0); enddo
      do jd = 1, nmat
        Linv_all(:,:,jd) = omat_all(:,:,jd)
        omat_all(:,:,jd) = Ih
      enddo
      do jd = 1, nmat
        is = mod(jd-1, NSTR) + 1
        istat = cublasZtrsm_v2(cb_str(is), CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER, &
             CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, nd, nd, (1d0,0d0), &
             Linv_all(1,1,jd), nd, omat_all(1,1,jd), nd)
      enddo
      do is = 1, NSTR; istat = cudaStreamSynchronize(strms(is)); enddo
      ! Swap: Linv_all = L^{-1} (in omat_all), omat_all freed for H'
      block
        complex(8), device, allocatable :: tmp(:,:,:)
        allocate(tmp(nd,nd,nmat)); tmp = omat_all; Linv_all = tmp; deallocate(tmp)
      endblock
      do is = 1, NSTR; istat = cublasDestroy(cb_str(is)); enddo
      deallocate(cb_str, Ih)
    endblock

    ! === Phase A3: H' = Linv * H * Linv^H via Ozaki GEMM ===
    do jd = 1, nmat
      block
        complex(8), device, allocatable :: Li_d(:,:), H_d(:,:), T_d(:,:), Hp_d(:,:)
        allocate(Li_d(nd,nd), H_d(nd,nd), T_d(nd,nd), Hp_d(nd,nd))
        Li_d = Linv_all(:,:,jd)
        H_d = hamm(1:nd, 1:nd, jd)
        istat = zmm(Li_d, H_d, T_d, m=nd, n=nd, k=nd)
        istat = zmm(T_d, Li_d, Hp_d, m=nd, n=nd, k=nd, opB=m_op_C)
        omat_all(:,:,jd) = Hp_d
        deallocate(Li_d, H_d, T_d, Hp_d)
      endblock
    enddo

    ! === Phase B: Zheevd on GPU streams ===
    block
      complex(8), device, allocatable :: dd(:,:)
      real(8), device, allocatable :: de(:)
      allocate(dd(nd,nd), de(nd))
      istat = cusolverDnZheevd_bufferSize(cs_str(1), CUSOLVER_EIG_MODE_VECTOR, &
           CUBLAS_FILL_MODE_UPPER, nd, dd, nd, de, lw_ev)
      deallocate(dd, de)
    endblock
    allocate(work_d(lw_ev * NSTR))
    do jd = 1, nmat
      is = mod(jd-1, NSTR) + 1
      istat = cusolverDnZheevd(cs_str(is), CUSOLVER_EIG_MODE_VECTOR, &
           CUBLAS_FILL_MODE_UPPER, nd, omat_all(1,1,jd), nd, eo_all(1,jd), &
           work_d((is-1)*lw_ev+1), lw_ev, devinfo)
    enddo
    do is = 1, NSTR; istat = cudaStreamSynchronize(strms(is)); enddo
    deallocate(work_d)

    ! === Phase C: D2H eigenvalues + Ozaki back-transform ===
    allocate(eo_h(nd, nmat)); eo_h = eo_all; deallocate(eo_all)
    do jd = 1, nmat
      nev = nev_list(jd)
      evals_out(1:nev, jd) = eo_h(1:nev, jd)
      evals_out(nev+1:nbandmx, jd) = 1d99
      block
        complex(8), device, allocatable :: Li_d(:,:), Y_d(:,:), Z_d(:,:)
        complex(8), allocatable :: Ztmp(:,:)
        allocate(Li_d(nd,nd), Y_d(nd,nev), Z_d(nd,nev))
        Li_d = Linv_all(:,:,jd)
        Y_d = omat_all(1:nd, 1:nev, jd)
        istat = zmm(Li_d, Y_d, Z_d, m=nd, n=nev, k=nd, opA=m_op_C)
        allocate(Ztmp(nd,nev)); Ztmp = Z_d
        evecs_out(1:nd, 1:nev, jd) = Ztmp
        deallocate(Li_d, Y_d, Z_d, Ztmp)
      endblock
    enddo
    deallocate(eo_h, omat_all, Linv_all)

    ! Cleanup
    do is = 1, NSTR
      istat = cusolverDnDestroy(cs_str(is))
      istat = cudaStreamDestroy(strms(is))
    enddo
    deallocate(strms, cs_str, devinfo)
  end subroutine
#endif
end module m_batched_diag
