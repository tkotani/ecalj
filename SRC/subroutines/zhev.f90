module m_zhev
  public zhev_tk4, zhev_gpu_cleanup
  private
contains
  subroutine zhev_gpu_cleanup()
    ! Release GPU handles and workspace after k-loop to free memory
#ifdef __GPU
    use cusolverdn
    use cublas_v2
    use cudafor
    use m_zhev_gpu_handles
#endif
    implicit none
#ifdef __GPU
    integer :: istat
    if(zhev_gpu_handles_init) then
      istat = cudaDeviceSynchronize()
      istat = cusolverDnDestroy(zhev_cusolver_handle)
      istat = cublasDestroy(zhev_cublas_handle)
      zhev_gpu_handles_init = .false.
    endif
#endif
  end subroutine
  subroutine zhev_tk4(n,h,s,nmx,nev, e,z, epsovl)
#ifdef __GPU
    use m_gpu, only: use_gpu
#endif
    implicit none
    integer :: n,nev,nmx,ltime,ngv,ncut
    !      logical ipr
    complex(8) :: h(n,n),s(n,n),z(n,nmx)
    complex(8),allocatable:: work(:)
    integer:: i,j,ni=999999,ik,ik2,k
    real(8):: epsovl,epsx,e(n),emx,eee,fac,vldummy,vudummy,abstol
    real(8),allocatable:: rwork(:)
    integer:: ier,lwork,ix,ifi,ifig,nmx0,ifail(n),nevx
    integer,allocatable:: iwork(:)
    character*1:: jobz
    logical ::nexist
    integer,save:: lworksave=0

    complex(8),allocatable::ii(:,:)
    real(8):: eo(n)
    complex(8),allocatable:: omat(:,:),wk11(:), &
         zz(:,:),hh(:,:),hhm(:,:),znm(:,:)
    integer:: nevl,nm,nmout,nevout
    logical:: debug=.false.

    if(allocated(omat)) deallocate(omat)
    allocate(omat(n,n))
    omat = s !reserved
    !
    if(epsovl< 1d-14) then
       call zhev_tk2(n,h,omat,nmx,nev, e,z)
       return
    endif
    call tcn('zhev_tk4')
#ifdef __GPU
    if(use_gpu) then
    gpudiag: block
      use cusolverdn
      use cublas_v2
      use cudafor
      use m_zhev_gpu_handles
      complex(8), device, allocatable :: omat_d(:,:), h_d(:,:), zz_d(:,:), hhm_d(:,:), hh_d(:,:), z_d(:,:)
      real(8), device, allocatable :: eo_d(:), e_d(:)
      complex(8), device, allocatable :: work_d(:)
      complex(8), allocatable :: zz_h(:,:)
      integer, device, allocatable :: devinfo
      integer :: istat2, lwork2, m_out
      if(.not. zhev_gpu_handles_init) then
        istat2 = cusolverDnCreate(zhev_cusolver_handle)
        istat2 = cublasCreate(zhev_cublas_handle)
        zhev_gpu_handles_init = .true.
      endif
      allocate(devinfo)
      ! Step 1: Diag overlap on GPU
      allocate(omat_d(n,n), eo_d(n))
      omat_d = omat
      istat2 = cusolverDnZheevdx_bufferSize(zhev_cusolver_handle, CUSOLVER_EIG_MODE_VECTOR, &
           CUSOLVER_EIG_RANGE_ALL, CUBLAS_FILL_MODE_UPPER, n, omat_d, n, &
           0d0, 0d0, 1, n, m_out, eo_d, lwork2)
      allocate(work_d(lwork2))
      istat2 = cusolverDnZheevdx(zhev_cusolver_handle, CUSOLVER_EIG_MODE_VECTOR, &
           CUSOLVER_EIG_RANGE_ALL, CUBLAS_FILL_MODE_UPPER, n, omat_d, n, &
           0d0, 0d0, 1, n, m_out, eo_d, work_d, lwork2, devinfo)
      deallocate(work_d)
      eo = eo_d
      do ix = 1, n
        if(eo(ix) > epsovl) then; ni = ix; exit; endif
      enddo
      nm = n - ni + 1
      nevl = nm
      ! Step 2: Build projection zz (host, then copy to GPU)
      allocate(zz_h(n, nm))
      omat = omat_d
      do ix = ni, n
        zz_h(:, ix-ni+1) = omat(:, ix) / sqrt(eo(ix))
      enddo
      allocate(zz_d(n, nm))
      zz_d = zz_h
      deallocate(zz_h, omat_d, eo_d)
      ! Step 3: Project H: hh = zz^H * H * zz
      allocate(h_d(n,n), hhm_d(nm,n), hh_d(nm,nm))
      h_d = h
      istat2 = cublasZgemm_v2(zhev_cublas_handle, CUBLAS_OP_C, CUBLAS_OP_N, nm, n, n, &
           (1d0,0d0), zz_d, n, h_d, n, (0d0,0d0), hhm_d, nm)
      istat2 = cublasZgemm_v2(zhev_cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N, nm, nm, n, &
           (1d0,0d0), hhm_d, nm, zz_d, n, (0d0,0d0), hh_d, nm)
      deallocate(hhm_d, h_d)
      ! Step 4: Diag reduced H on GPU
      if(nmx==0) then
        nev = nm
      else
        nev = min(nmx, nm)
      endif
      allocate(e_d(nm))
      istat2 = cusolverDnZheevdx_bufferSize(zhev_cusolver_handle, CUSOLVER_EIG_MODE_VECTOR, &
           CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_UPPER, nm, hh_d, nm, &
           0d0, 0d0, 1, nev, m_out, e_d, lwork2)
      allocate(work_d(lwork2))
      istat2 = cusolverDnZheevdx(zhev_cusolver_handle, CUSOLVER_EIG_MODE_VECTOR, &
           CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_UPPER, nm, hh_d, nm, &
           0d0, 0d0, 1, nev, m_out, e_d, work_d, lwork2, devinfo)
      deallocate(work_d)
      e(1:nev) = e_d(1:nev)
      nev = m_out
      deallocate(e_d)
      ! Step 5: Back-transform z = zz * hh_d(:,1:nev)
      z = 1d99
      allocate(z_d(n, nev))
      istat2 = cublasZgemm_v2(zhev_cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N, n, nev, nm, &
           (1d0,0d0), zz_d, n, hh_d, nm, (0d0,0d0), z_d, n)
      z(1:n, 1:nev) = z_d
      deallocate(zz_d, hh_d, z_d)
      deallocate(devinfo)
      istat2 = cudaDeviceSynchronize()
    endblock gpudiag
    else
#endif
    !! ====== CPU path: LAPACK ======
    !! ... eigenvalue of ovarlap matrix
    jobz = 'V'
    lwork = n*n
    allocate(wk11(lwork),rwork(max(1,3*n-2)))
    call zheev(jobz,'U',n,omat,n,eo,wk11,lwork,rwork,ier)
    deallocate(wk11,rwork)
    if(debug) then
       write(6,*)'zhev_tk4: ovlmat='
       do ix=0,n,5
          write(6,"(5(i5,d10.2))") (i,eo(i),i=ix+1,min(ix+5,n))
       enddo
    endif
    do ix= 1,n
       if(eo(ix)>epsovl) then
          ni = ix               ! take i = ni...n
          exit
       endif
    enddo
    nm = n-ni+1            ! this is the dimension.
    nevl=nm
    allocate(zz(n,nm))     ! zz is the projection matrix
    do ix=ni,n
       zz(:,ix-ni+1) = omat(:,ix)/sqrt(eo(ix))
    enddo
    if(debug) then
       write(6,*)' reduced by OVEPS: n--> nm=',nm
    endif
    !! Hamiltonian  <zz|H|zz>
    allocate(hh(nm,nm),hhm(nm,n))
    call zgemm('C','N',nm,n,n,(1d0,0d0),zz,n,h,n,(0d0,0d0),hhm,nm)
    call zgemm('N','N',nm,nm,n,(1d0,0d0),hhm,nm,zz,n,(0d0,0d0),hh,nm)
    deallocate(hhm)
    if(nmx==0) then
       jobz='N'
       nev=nm
    else
       jobz = 'V'
       nev = min(nmx,nm)
    endif
    abstol= 1d-10
    lwork = max(1,2*nm,lworksave)
    allocate(work(lwork),rwork(7*nm),iwork(5*nm),znm(nm,max(1,nev)))
    call zheevx(jobz,'I','U',nm,hh,nm,vldummy,vudummy,1,nev,abstol,nevout,e,znm,nm,work,lwork,rwork,iwork,ifail,ier)
    lworksave= WORK(1)
    call rxx(nev/=nevout,'zhev_tk4: nev /=nevout something wrong. ')
    call rxx(ier.ne.0, 'zhev_tk4: zheev for hh cause error.')
    deallocate(work,iwork,rwork)
    z=1d99
    call zgemm('N','N',n, min(nmx,nm),nm,(1d0,0d0),zz,n,znm,nm,(0d0,0d0),z,n)
#ifdef __GPU
    endif
#endif
    if(allocated(znm)) deallocate(znm)
    if(allocated(zz)) deallocate(zz)
    if( .FALSE. ) then !! === diagonalize === (this part is in zhev_tk2), Kept here for debug purpose
       if(nmx==0) then
          jobz='N'
          nev=n
       else
          jobz = 'V'
          nev=nmx
       endif
       abstol=1d-10            ! OK?
       lwork=max(1,2*n,lworksave) !OK? efficient?
       allocate(work(lwork),rwork(7*n),iwork(5*n))
       call zhegvx(1,jobz,'I','U',n,h,n,s,n,vldummy,vudummy,1,nev,abstol,nevx,e,z,n,work,lwork,rwork,iwork,ifail,ier)
       lworksave= WORK(1)      !this is optimum lwork
       call rxx(nev/=nevx,'zhev_tk4: nev /=nevx something wrong. ')
       call rxx(ier.ne.0, 'zhev_tk4: zhegvx cannot find all eigen.')
       deallocate(work,iwork,rwork)
    endif
    if(nmx/=0) then
       phaselock: block ! a phaselock for continuity of evec as for ham and ovl. Phase of evec is paralell to (1,1,1,...1)*ovl^-1 
       integer::iev !z-->evec
       complex(8)::z0(n)
       complex(8),parameter::img=(0d0,1d0)
       z0=[(1d0/(1d0+0.01d0*i),i=1,n)] !1,1,1,... may cause z//z0 because of some symmetry ; I am afraid that sum=0 causing error.
       forall(iev=1:nev) z(:,iev)=z(:,iev)*exp(-img*dimag(log(sum(z0*z(:,iev)))))
       endblock phaselock
    endif
    !
    call tcx('zhev_tk4')
  end subroutine zhev_tk4
  subroutine zhev_tk2(n,h,s,nmx,nev, e,z)
    !!== Eigenvalues and/or some eigenvectors of a Hermitian matrix (weighted for first nlmto basis).==
    !! ----------------------------------------------------------------
    !! Inputs:
    ! c   nlmto:dimension of MTO space of 1:nlmto i respected when diagonalization.
    !!     n:    dimension of h
    !!   h,n:  hermitian matrix, dimensioned h(n,n)
    !!   s:    hermitian overlap matrix,
    !!   nmx:  requested number of eigenvectors to be found (and eigenvalues). If nmx>n, nmx is taken to be n.
    !!         if nmx=0, nev=n (see NOTE below).
    !!   ipr :print switch
    !!   ifig, savez,getz:dummy
    !! Outputs:
    !!   e:    eigenvalues
    !!   nev:  number of eigenvectors (=nmx) or (=n if nmx=0)
    !!   z:    eigenvectors (1..nev)  (declared as z(n,*)
    !!   h and s are destroyed on exit.
    !!   july2012takao
    !! NOTE: this can be called in the loop of ikp,isp loop. Then data are appended to a file ifig.
    !!   If nmx==0, all eigenvalues are returned but without eigenfunctions.
    !! -----------------------------------------------------------------------
    implicit none
    integer :: n,nev,nmx,ltime,ngv,ncut
    !      logical ipr
    complex(8) :: h(n,n),s(n,n),z(n,*)
    complex(8),allocatable:: work(:)
    integer:: i,j,ni=999999,ik,ik2,k
    real(8):: epsovl,epsx,e(n),emx,eee,fac,vldummy,vudummy,abstol
    real(8),allocatable:: rwork(:)
    integer:: ier,lwork,ix,ifi,ifig,nmx0,ifail(n),nevx
    integer,allocatable:: iwork(:)
    character*1:: jobz
    logical ::nexist
    integer,save:: lworksave=0
    call tcn('zhev_tk2')
    !! === diagonalize ===
    if(nmx==0) then
       jobz='N'
       nev=n
    else
       jobz = 'V'
       nev=nmx
    endif
    abstol=1d-10 ! OK?
    lwork=max(1,2*n,lworksave) !OK? efficient?
    allocate(work(lwork),rwork(7*n),iwork(5*n))
    call zhegvx(1,jobz,'I','U',n,h,n,s,n,vldummy,vudummy,1,nev,abstol,nevx,e,z,n,work,lwork,rwork,iwork,ifail,ier)
    lworksave= WORK(1)  !this is optimum lwork
    !      print *,'nev nevx n=',nev,nevx,n
    call rxx(nev/=nevx,'zhev_tk2: nev /=nevx something wrong. ')
    call rxx(ier.ne.0, 'zhev_tk2: zhegvx cannot find all eigen.')
    deallocate(work,iwork,rwork)
    call tcx('zhev_tk2')
  end subroutine zhev_tk2
  subroutine zhevx(n,lh,h,s,lov,lx,nmx,emx,nev,wk,linv,e,lz,z)! Eigenvalues and/or some eigenvectors of a Hermitian matrix
    !i Inputs:
    !i   n:    order of h and s
    !i   lh:   leading dimension of h and s
    !i   h:    hermitian matrix, dimensioned h(n,n)
    !i   s:    hermitian overlap matrix, (used only if lov is true)
    !i   nmx:  maximum number of eigenvectors to be found
    !i   emx:  eigenvalue limit for eigenvectors to be found
    !i         (not used if LAPACK zhegv is invoked)
    !i   wk:   work array of length at least 11n
    !i         NB: If LAPACK version is used, and eigenvectors are sought
    !i         wk should be dimensioned (n*nmx*2)
    !i   lov:  0 no overlap matrix
    !i         1 overlap matrix, return evecs of nonorthogonal H
    !i   lx:   if T, calls routines to exploit unit stride lengths (risc)
    !i         Not used if LAPACK zhegv is invoked.
    !i   linv: if T, using inverse iteration
    !i         Not used if LAPACK zhegv is invoked.
    !i   lz:   leading dimension of z
    !o Outputs:
    !o   e:    eigenvalues
    !o   nev:  number of eigenvectors found
    !o   z:    eigenvectors (1..nev)  (declared as z(n,*)
    !o   s:    has been decomposed into and LL+ decomposition.
    !o         You can call zhev2 to scale a vector by L
    !r Remarks:
    !r   z must be at least of dimension z(n,n), even though nev<n.
    !r   h and s are destroyed on exit.
    !r   Aborts on exit
    !p Procedures used:
    !p   (lapack)  zhegv
    !p   (eispack) htribk, htridx, imtql2, tqlrat
    !u Updates
    !u   24 Feb 07 Bug fix when nmx=0
    !u   17 May 03 Adapted from zhev, intended to supersede zhev.
    !u   14 Aug 02 Added zheev when lov is F; new zhev2.
    !u   21 Jan 02 Added code to invoke LAPACK zhegv in place of diagno
    ! ----------------------------------------------------------------
    implicit none
    logical :: linv,lx
    integer :: lov,lh,lz,nn
    integer :: n,nev,nmx
    double precision :: h(n,n),s(*),e(n),wk(*),emx
    integer :: ier,lwork
    character jobz
    complex(8)::z(n,*)
    call tcn('zhevx')
    if(lh/=n) call rx('zhevx: we assume lh=n now')
    if (nmx <= 0) then
       jobz = 'N'
       lwork = 4*n
       if (lov > 0) then
          call zhegv(1,jobz,'U',n,h,lh,s,lh,e,wk(1+3*n),lwork,wk(1),ier)
          call rxx(ier.ne.0,'zhevx: zhegv cannot find all evals')
       else
          call zheev(jobz,'U',n,h,lh,e,wk(1+3*n),lwork,wk(1),ier)
          call rxx(ier.ne.0,'zhevx: zheev cannot find all evals')
       endif
       nev = 0
    else
       jobz = 'V'
       lwork = n*nmx
       if (lov > 0) then
          call zhegv(1,jobz,'U',n,h,lh,s,lh,e,z,lwork,wk(1),ier)
          call rxx(ier.ne.0,'zhevx: zhegv cannot find all evals')
       else
          call zheev(jobz,'U',n,h,lh,e,z,lwork,wk(1),ier)
          call rxx(ier.ne.0,'zhevx: zheev cannot find all evals')
       endif
       nn=min(n,nmx)
       z(1:n,1:nn)=h(1:n,1:nn)
       !        call zmcpy('N',h,lh,1,z,lz,1,n,min(n,nmx))
       !       call zprm('evecs',2,z,lz,n,nmx)
       nev = min(n,nmx)
    endif
100 call tcx('zhevx')
  end subroutine zhevx
  subroutine zhev(n,h,s,lov,lx,nmx,emx,nev,wk,linv,ltime,e,z)
    !- Eigenvalues and/or some eigenvectors of a Hermitian matrix
    ! ----------------------------------------------------------------
    !i Inputs:
    !i   n:    dimension of h
    !i   h,n:  hermitian matrix, dimensioned h(n,n)
    !i   s:    hermitian overlap matrix, (used only if lov is true)
    !i   nmx:  maximum number of eigenvectors to be found
    !i   emx:  eigenvalue limit for eigenvectors to be found
    !i   wk:   work array of length at least 11n
    !i   lov:  if T, non-orthogonal
    !i   lx:   if T, calls routines to exploit unit stride lengths (risc)
    !i         Not used if LAPACK zhegv is invoked.
    !i   linv: if T, using inverse iteration
    !i         Not used if LAPACK zhegv is invoked.
    !o Outputs:
    !o   e:    eigenvalues
    !o   nev:  number of eigenvectors found
    !o   z:    eigenvectors (1..nev)  (declared as z(n,*)
    !o   s:    has been decomposed into and LL+ decomposition.
    !o         You can call zhev2 to scale a vector by L
    !r Remarks:
    !r   z must be at least of dimension z(n,n), even though nev<n.
    !r   h and s are destroyed on exit.
    !r   Aborts on exit
    !p Procedures used:
    !p   (lapack)  zhegv
    !p   (eispack) htribk, htridx, imtql2, tqlrat
    !u Updates
    !u   14 Aug 02 Added zheev when lov is F; new zhev2.
    !u   21 Jan 02 Added code to invoke LAPACK zhegv in place of diagno
    ! ----------------------------------------------------------------
    implicit none
    logical :: linv,lx
    integer :: n,nev,nmx,ltime
    double precision :: h(*),s(*),e(n),wk(*),z(*),emx
    logical :: lov
    integer :: ier,lwork
    character jobz
    call tcn('zhev')
    if (nmx <= 0) then
       jobz = 'N'
       lwork = 4*n
       if (lov) then
          call zhegv(1,jobz,'U',n,h,n,s,n,e,wk(1+3*n),lwork,wk(1),ier)
       else
          call zheev(jobz,'U',n,h,n,e,wk(1+3*n),lwork,wk(1),ier)
       endif
       nev = 0
    else
       jobz = 'V'
       lwork = n*min(n,nmx)
       if (lov) then
          call zhegv(1,jobz,'U',n,h,n,s,n,e,z,lwork,wk(1),ier)
       else
          call zheev(jobz,'U',n,h,n,e,z,lwork,wk(1),ier)
       endif
       call zcopy(n*min(n,nmx),h,1,z,1)
       nev = min(n,nmx)
    endif
    call rxx(ier.ne.0,'zhev: zhegv cannot find all evals')
100 call tcx('zhev')
  end subroutine zhev
end module m_zhev
